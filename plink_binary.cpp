#include "plink_binary.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using std::string;
using std::stringstream;
using std::vector;

// final '1' for snp major mode; only supporting this at present
#define MAGIC_LEN 3
static int magic_number[MAGIC_LEN] = { 108, 27, 1 };

plink_binary::plink_binary(void)
{}

plink_binary::plink_binary(string dataset) {
    // single argument: open as read (default)
    plink_binary::open(dataset);
}

plink_binary::~plink_binary() {}

void plink_binary::close() {
    if (open_for_write) {
        if (individuals.size() == 0) {
            throw gftools::malformed_data("No individuals to write");
        }
        write_fam(individuals);
        if (snps.size() == 0) {
            throw gftools::malformed_data("No SNPs to write");
        }
        write_bim(snps);

        if (bed_file->is_open()) {
            bed_file->close();
        }
        delete bed_file;
    }
    else {
        if (is_mem_mapped) {
            munmap(fmap, flen);
            ::close(fd);
        }
        else {
            if (bed_file->is_open()) {
                bed_file->close();
            }
            delete bed_file;
        }
    }
    individuals.resize(0);
    snps.resize(0);
    snp_index.clear();
}

void plink_binary::open(string dataset) {
    init(dataset, 0);
}

void plink_binary::open(string dataset, bool mode) {
    init(dataset, mode);
}

void plink_binary::init(string dataset, bool mode) {
    this->dataset = dataset;
    // default "NN" for no call
    missing_genotype = 'N';
    quell_mem_mapping = false;
    if (mode) {
        open_for_write = 1;
        open_bed_write(dataset + ".bed");
    }
    else {
        open_for_write = 0;
        read_bim(snps);
        if (snps.size() == 0) {
            throw gftools::malformed_data("No SNPs read");
        }
        read_fam(individuals);
        if (individuals.size() == 0) {
            throw gftools::malformed_data("No individuals read");
        }
        bytes_per_snp = (3 + individuals.size()) / 4;
        open_bed_read(dataset + ".bed", quell_mem_mapping);
    }
}

bool plink_binary::is_empty(std::ifstream &ifstream) {
    return ifstream.peek() == std::ifstream::traits_type::eof();
}

bool plink_binary::next_snp(gftools::snp &s, vector<string> &genotypes) {
    vector<int> gt_int;
    if (snp_ptr >= snps.size()) {
        return false;
    }

    get_snp(3 + snp_ptr * bytes_per_snp, gt_int);
    genotypes_itoa(snps[snp_ptr], gt_int, genotypes);
    s = snps[snp_ptr++];
    return true;
}

void plink_binary::read_snp(std::string snp, vector<string> &genotypes) {
    vector<int> gt_int;
    int index = snp_index[snp];
    read_snp(3 + index * bytes_per_snp, gt_int);
    genotypes_itoa(snps[index], gt_int, genotypes);
    snp_ptr = index + 1;
}

void plink_binary::genotypes_itoa(gftools::snp s, vector<int> g_num, vector<string> &g_str) {
    g_str.resize(0);

    for (unsigned int i = 0; i < g_num.size(); i++) {
        string genotype;
        string missing;
        missing.append(&missing_genotype);
        missing.append(&missing_genotype);
        switch (g_num[i]) {
            case 0: genotype = missing; break;
            case 1: genotype = s.allele_a + s.allele_a; break;
            case 2: genotype = s.allele_a + s.allele_b; break;
            case 3: genotype = s.allele_b + s.allele_b; break;
            default: throw gftools::malformed_data("integer genotypes must be 0..3");
        }
        g_str.push_back(genotype);
    }
}

void plink_binary::genotypes_atoi(gftools::snp &snp, const vector<string> g_str, vector<int> &g_num) {
    std::vector<string> alleles = collate_alleles(g_str);
    string a = alleles[0];
    string b = alleles[1];
    string missing = string(1, missing_genotype);
    string no_call = string(2, missing_genotype);

    std::map<string, int> lookup;
    if ((a == missing && b != missing) || (a != missing && b == missing)) {
        throw gftools::malformed_data("Missing a call for only one allele");
    }
    else {
        // Strand is being normalised here, for heterozygotes
        lookup[no_call] = 0;
        lookup[snp.allele_a + snp.allele_a] = 1;
        lookup[snp.allele_a + snp.allele_b] = 2;
        lookup[snp.allele_b + snp.allele_a] = 2;
        lookup[snp.allele_b + snp.allele_b] = 3;
    }

    for (vector<string>::const_iterator i = g_str.begin(); i < g_str.end(); i++) {
        string call = *i;
        if (snp.is_known() && lookup.count(call) == 0) {
            throw gftools::malformed_data("Unexpected call for SNP " +
                                          snp.name + " " +
                                          snp.allele_a + snp.allele_b +
                                          ": " + call);
        }
        g_num.push_back(lookup[call]);
    }
}

vector<string> plink_binary::collate_alleles(const vector<string> &g_str) {
    std::set<string> allele_a;
    std::set<string> allele_b;
    std::set<string> allele_ab;
    string missing = string(1, missing_genotype);

    vector<string>::const_iterator gi;
    for (gi = g_str.begin(); gi < g_str.end(); gi++) {
        string call = *gi;

        if (call.length() == 2) {
            string a = call.substr(0, 1);
            string b = call.substr(1, 2);

            if (a != missing) {
                allele_a.insert(a);
                allele_ab.insert(a);
            }
            if (b != missing) {
                allele_b.insert(b);
                allele_ab.insert(b);
            }
        }
        else {
            throw gftools::malformed_data("Unrecognised genotype " + call);
        }
    }

    if (allele_ab.size() > 3) {
        std::stringstream ss;
        ss << "[";
        for (gi = g_str.begin(); gi < g_str.end(); gi++) {
            if (gi != g_str.begin()) {
                ss << ", ";
            }
            ss << *gi;
        }
        ss << "]";

        throw gftools::malformed_data("Collated genotype was not bi-allelic: " +
                                      ss.str());
    }

    vector<string> alleles;
    std::set<string>::iterator it;
    if (allele_a.empty()) {
        alleles.push_back(missing);
    }
    else {
        for (it = allele_a.begin(); it != allele_a.end(); it++) {
            alleles.push_back(*it);
        }
    }
    if (allele_b.empty()) {
        alleles.push_back(missing);
    }
    else {
        for (it = allele_b.begin(); it != allele_b.end(); it++) {
            alleles.push_back(*it);
        }
    }

    return alleles;
}

void plink_binary::read_snp(int snp_index, vector<int> &genotypes) {
    get_snp(3 + snp_index * bytes_per_snp, genotypes);
}

gftools::snp plink_binary::from_bim(string record) {
    stringstream ss(record);
    gftools::snp s;
    ss >> s.chromosome;
    ss >> s.name;
    ss >> s.genetic_position;
    ss >> s.physical_position;
    ss >> s.allele_a;
    ss >> s.allele_b;
    return s;
}

string plink_binary::to_bim(gftools::snp s) {
    stringstream record;
    record << (s.chromosome.length() ? s.chromosome : "0");
    record << "\t" << s.name;
    record << "\t" << (s.genetic_position ? s.genetic_position : 0);
    record << "\t" << (s.physical_position ? s.physical_position : 0);
    record << "\t" << s.allele_a;
    record << "\t" << s.allele_b;
    return record.str();
}

gftools::individual plink_binary::from_fam(string record) {
    stringstream ss(record, stringstream::in | stringstream::out);
    gftools::individual ind;
    ss >> ind.family;
    ss >> ind.name;
    ss >> ind.father;
    ss >> ind.mother;
    ss >> ind.sex;
    ss >> ind.phenotype;
    return ind;
}

string plink_binary::to_fam(gftools::individual ind) {
    stringstream record;
    record << (ind.family.length() ? ind.family : ind.name);
    record << "\t" << ind.name;
    record << "\t" << (ind.father.length() ? ind.father : "-9");
    record << "\t" << (ind.mother.length() ? ind.mother : "-9");
    record << "\t" << (ind.sex.length() ? ind.sex : "-9");
    record << "\t" << (ind.phenotype.length() ? ind.phenotype : "-9");
    return record.str();
}

void plink_binary::write_bim(vector<gftools::snp> snps) {
    std::ofstream file;
    string fn = dataset + ".bim";
    file.open(fn.c_str(), std::fstream::out);

    for (unsigned int i = 0; i < snps.size(); i++) {
        file << to_bim(snps[i]) << std::endl;
    }

    file.close();
}

void plink_binary::read_bim(vector<gftools::snp> &snps) {
    std::ifstream file;
    string fn = dataset + ".bim";
    file.open(fn.c_str());

    if (!file) {
        throw gftools::malformed_data("Failed to open BIM file for " + dataset +
                                      ": " + error_message());
    }
    else {
        if (is_empty(file)) {
            throw gftools::malformed_data("Empty BIM file for " + dataset);
        }
    }

    string str;
    int i = 0;

    while (getline(file, str)) {
        gftools::snp s = from_bim(str);
        snps.push_back(s);
        snp_index[s.name] = i++;
    }

    file.close();
}

void plink_binary::write_fam(vector<gftools::individual> individuals) {
    std::ofstream file;
    string fn = dataset + ".fam";
    file.open(fn.c_str(), std::fstream::out);

    for (unsigned int i = 0; i < individuals.size(); i++) {
        file << to_fam(individuals[i]) << std::endl;
    }

    file.close();
}

void plink_binary::read_fam(vector<gftools::individual> &individuals) {
    std::ifstream file;
    string fn = dataset + ".fam";
    file.open(fn.c_str());
    if (!file.is_open()) {
        throw gftools::malformed_data("Missing FAM file for " + dataset);
    }
    string s;
    int i = 0;

    while (getline(file, s)) {
        gftools::individual ind = from_fam(s);
        individuals.push_back(ind);
        i++;
    }

    file.close();
}

void plink_binary::open_bed_write(string filename) {
    bed_file = new std::fstream();
    bed_file->open(filename.c_str(), std::fstream::out | std::fstream::binary);

    // in theory can use mmap here but you'd have to declare the
    // snps (or at least the snp count) before opening. as most
    // files will likely be opened for reading, (yet) supported
    is_mem_mapped = 0;

    write_bed_header();
}

void plink_binary::open_bed_read(string filename, bool quell_mem_mapping) {
    if (quell_mem_mapping) {
        is_mem_mapped = 0;
        bed_file = new std::fstream();
        bed_file->open(filename.c_str(), std::ios::in);
        if (!bed_file) {
            throw gftools::malformed_data("Failed to open BED file for " +
                                          dataset + ": " + error_message());
        }
    }
    else {
        struct stat result;

        if (stat(filename.c_str(), &result) == -1) {
            throw gftools::malformed_data("Failed to stat BED file: " +
                                          error_message());
        }

        flen = result.st_size;
        fd = ::open(filename.c_str(), O_RDONLY);
        if (fd == -1) {
            throw gftools::malformed_data("Failed to open BED file for " +
                                          dataset + ": " + error_message());
        }

        fmap = (char *) mmap(0, flen, PROT_READ, MAP_PRIVATE, fd, 0);
        if (fmap == MAP_FAILED) {
            throw gftools::malformed_data("Failed to memory map BED file for " +
                                          dataset + ": " + error_message());
        }
        else {
            is_mem_mapped = 1;
        }
    }

    read_bed_header();
    snp_ptr = 0;

}

void plink_binary::write_bed_header() {
    char buffer[MAGIC_LEN];

    for (int i = 0; i < MAGIC_LEN; i++) {
        buffer[i] = magic_number[i];
    }

    if (is_mem_mapped) {
        memcpy(fmap, buffer, MAGIC_LEN);
    }
    else {
        bed_file->write(buffer, MAGIC_LEN);
    }
}

void plink_binary::read_bed_header() {
    char buffer[MAGIC_LEN + 1];

    if (is_mem_mapped) {
        memcpy(buffer, fmap, MAGIC_LEN + 1);
    }
    else {
        bed_file->seekg(0, std::ios_base::beg);
        bed_file->read(buffer, MAGIC_LEN + 1);
    }

    for (int i = 0; i < MAGIC_LEN; i++) {
        if (buffer[i] != magic_number[i]) {
            throw gftools::malformed_data("Corrupt or incompatible BED file?");
        }
    }
}

void plink_binary::extract_bed(size_t pos, size_t len, char *buffer) {
    if (is_mem_mapped) {
        memcpy(buffer, fmap + pos, len);
    }
    else {
        bed_file->seekg(pos, std::ios_base::beg);
        bed_file->read(buffer, len);
    }
}

void plink_binary::get_snp(size_t pos, vector<int> &genotypes) {
    char *buffer = (char *)malloc(bytes_per_snp);

    extract_bed(pos, bytes_per_snp, buffer);
    genotypes.resize(0);
    uncompress_calls(buffer, individuals.size(), genotypes);
    free(buffer);
}

void plink_binary::write_snp(gftools::snp s, vector<int> genotypes) {
    if (genotypes.size() == 0) {
        throw gftools::malformed_data("No genotypes defined");
    }
    if (individuals.size() == 0) {
        throw gftools::malformed_data("No individuals defined");
    }
    if (genotypes.size() != individuals.size()) {
        stringstream ss;
        ss << "Incorrect individual count: ";
        ss << genotypes.size();
        ss << " genotypes for SNP ";
        ss << s.name;
        ss << " whereas ";
        ss << individuals.size();
        ss << " individuals defined";
        throw gftools::malformed_data(ss.str());
    }
    bed_write(s, genotypes);
}

void plink_binary::write_snp(gftools::snp s, vector<string> genotypes) {
    vector<int> g_num;
    genotypes_atoi(s, genotypes, g_num);
    write_snp(s, g_num);
}

void plink_binary::bed_write(gftools::snp s, vector<int> genotypes) {
    int len = (3 + individuals.size()) / 4;
    char buffer[len];
    compress_calls(buffer, genotypes);
    snps.push_back(s);
    bed_file->write(buffer, len);
}

// Encode/decode from plink encoding.
// In this class we store genotypes as int: 0 (no call); 1 (AA); 2 (AB); 3 (BB)
// Plink encodes in two bits as:            1 (no call); 0 (AA); 2 (AB); 3 (BB)
// Note that the bytes in a bed file are written in reverse order:
// in the plink documentation a missing genotype is binary 10;
// 01 (as above) after reversal.
void plink_binary::uncompress_calls(char *buffer, size_t len, vector<int> &calls) {
    char c = buffer[0];
    for (unsigned int i = 0; i < len; i++) {
        if (!(i % 4)) {
            c = buffer[i / 4];             // if div by 4 => get byte
        }
        char d = c & 3;
        c >>= 2;
        switch(d) {
            case 1: calls.push_back(0); break;
            case 2: calls.push_back(2); break;
            case 3: calls.push_back(3); break;
            case 0: calls.push_back(1); break;
        }
    }
}

void plink_binary::compress_calls(char *buffer, vector<int> calls) {
    int pos = 0;
    unsigned char c = 0;
    for (unsigned int i = 0; i < calls.size(); i++) {
        c >>= 2;
        switch(calls[i]) {
            case 1: c |= (0 << 6); break;
            case 2: c |= (2 << 6); break;
            case 3: c |= (3 << 6); break;
            // treat anything else as a no call;
            // may want to have a strict option here?
            default: c |= (1 << 6);
        }
        // write every fourth i
        if ((i % 4) == 3) {
            buffer[pos++] = c;
            c = 0;
        }
    }
    // if not a multiple of four, need to right shift stored data
    // to low end of byte
    if (calls.size() % 4) {
        int unfilled = calls.size();
        while (unfilled++ % 4) {
            c >>= 2;
        }
        buffer[pos] = c;
    }
}

string plink_binary::error_message() {
    char *msg  = strerror(errno);

    return msg ? string(msg) : "unknown error";
}