#include "plink_binary.h"
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// final '1' for snp major mode; only supporting this at present
#define MAGIC_LEN 3
static int magic_number[MAGIC_LEN] = { 108, 27, 1 };

plink_binary::plink_binary(void)
{
}

plink_binary::plink_binary(std::string dataset)
{
    // single argument: open as read (default)
    plink_binary::open(dataset);
}

plink_binary::~plink_binary() {
}

void plink_binary::close()
{
    if (open_for_write) {
        if (individuals.size() == 0) throw gftools::malformed_data("No individuals to write");
        write_fam(individuals);
        if (snps.size() == 0) throw gftools::malformed_data("No SNPs to write");
        write_bim(snps);
        bed_file->close();
        delete bed_file;
    } else {
        if (is_mem_mapped) {
            munmap(fmap, flen);
            ::close(fd);
        } else {
            bed_file->close();
            delete bed_file;
        }
    }
    individuals.resize(0);
    snps.resize(0);
    snp_index.clear();
}

void plink_binary::open(std::string dataset)
{
    init(dataset, 0);
}

void plink_binary::open(std::string dataset, bool mode)
{
    init(dataset, mode);
}

void plink_binary::init(std::string dataset, bool mode)
{
    this->dataset = dataset;
    // default "NN" for no call
    missing_genotype = 'N';
    quell_mem_mapping = false;
    if (mode) {
        open_for_write = 1;
        open_bed_write(dataset + ".bed");
    } else {
        open_for_write = 0;
        read_bim(snps);
        if (snps.size() == 0) throw gftools::malformed_data("No SNPs read");
        read_fam(individuals);
        if (individuals.size() == 0) throw gftools::malformed_data("No individuals read");
        bytes_per_snp = (3 + individuals.size()) / 4;
        open_bed_read(dataset + ".bed", quell_mem_mapping);
    }
}

bool plink_binary::next_snp(gftools::snp &s, std::vector<std::string> &genotypes)
{
    std::vector<int> gt_int;
    if (snp_ptr >= snps.size())
        return false;
    get_snp(3 + snp_ptr * bytes_per_snp, gt_int);
    genotypes_itoa(snps[snp_ptr], gt_int, genotypes);
    s = snps[snp_ptr++];
    return true;
}

void plink_binary::read_snp(std::string snp, std::vector<std::string> &genotypes)
{
    std::vector<int> gt_int;
    int index = snp_index[snp];
    read_snp(3 + index * bytes_per_snp, gt_int);
    genotypes_itoa(snps[index], gt_int, genotypes);
    snp_ptr = index + 1;
}

void plink_binary::genotypes_itoa(gftools::snp s, std::vector<int> g_num, std::vector<std::string> &g_str)
{
    g_str.resize(0);

    for (unsigned int i = 0; i < g_num.size(); i++) {
        std::string genotype;
        std::string missing;
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

void plink_binary::genotypes_atoi(gftools::snp &s, std::vector<std::string> g_str, std::vector<int> &g_num)
{
    // scan genotypes to get the allele codes required
    std::set<std::string> alleles;
    std::string missing;
    missing.append(&missing_genotype);
    missing.append(&missing_genotype);
    if (! s.allele_a.length() || ! s.allele_b.length()) {
        for (unsigned int i = 0; i < g_str.size(); i++) {
            if (g_str[i].substr(0, 2) != missing) {
                // better way of doing this?
                if (alleles.find(g_str[i].substr(0, 1)) == alleles.end())
                    alleles.insert(g_str[i].substr(0, 1));
                if (alleles.find(g_str[i].substr(1, 1)) == alleles.end())
                    alleles.insert(g_str[i].substr(1, 1));
                if (g_str[i].length() != 2)
                    throw gftools::malformed_data("Unrecognised genotype " + g_str[i]);
            }
        }
    }

    std::vector<std::string> labels;   // should be 0, 1 or 2 allele codes
    std::set<std::string>::iterator it;
    for (it = alleles.begin(); it != alleles.end(); it++)
        labels.push_back(*it);

    if (labels.size() == 3) {
        std::stringstream ss;
        ss << "More than two alleles found for SNP " + s.name;
        for (unsigned int i = 0; i < labels.size(); i++)
            ss << " " << labels[i];
        throw gftools::malformed_data(ss.str());
    }
    if (labels.size() == 0) {
        // no genotypes - return all no-calls
        for (unsigned int i = 0; i < g_str.size(); i++)
            g_num.push_back(0);
        s.allele_a = s.allele_b = "0";
        return;
    }
    if (labels.size() == 1)
        labels.push_back("0");

    s.allele_a = labels[0];
    s.allele_b = labels[1];

    std::map<std::string, int> lookup;

    lookup[missing] = 0;
    lookup[labels[0] + labels[0]] = 1;
    lookup[labels[0] + labels[1]] = 2;
    lookup[labels[1] + labels[0]] = 2;
    lookup[labels[1] + labels[1]] = 3;
    for (unsigned int i = 0; i < g_str.size(); i++)
        g_num.push_back(lookup[g_str[i]]);
}

void plink_binary::read_snp(int snp_index, std::vector<int> &genotypes)
{
    get_snp(3 + snp_index * bytes_per_snp, genotypes);
}

gftools::snp plink_binary::from_bim(std::string record) {
    std::stringstream ss(record);
    gftools::snp s;
    ss >> s.chromosome;
    ss >> s.name;
    ss >> s.genetic_position;
    ss >> s.physical_position;
    ss >> s.allele_a;
    ss >> s.allele_b;
    return s;
}

std::string plink_binary::to_bim(gftools::snp s) {
    std::stringstream record;
    record << (s.chromosome.length() ? s.chromosome : "0");
    record << "\t" << s.name;
    record << "\t" << (s.genetic_position ? s.genetic_position : 0);
    record << "\t" << (s.physical_position ? s.physical_position : 0);
    record << "\t" << s.allele_a;
    record << "\t" << s.allele_b;
    return record.str();
}

gftools::individual plink_binary::from_fam(std::string record) {
    std::stringstream ss(record, std::stringstream::in | std::stringstream::out);
    gftools::individual ind;
    ss >> ind.family;
    ss >> ind.name;
    ss >> ind.father;
    ss >> ind.mother;
    ss >> ind.sex;
    ss >> ind.phenotype;
    return ind;
}

std::string plink_binary::to_fam(gftools::individual ind) {
    std::stringstream record;
    record << (ind.family.length() ? ind.family : ind.name);
    record << "\t" << ind.name;
    record << "\t" << (ind.father.length() ? ind.father : "-9");
    record << "\t" << (ind.mother.length() ? ind.mother : "-9");
    record << "\t" << (ind.sex.length() ? ind.sex : "-9");
    record << "\t" << (ind.phenotype.length() ? ind.phenotype : "-9");
    return record.str();
}

void plink_binary::write_bim(std::vector<gftools::snp> snps)
{
    std::ofstream file;
    std::string fn = dataset + ".bim";
    file.open(fn.c_str(), std::fstream::out);

    for (unsigned int i = 0; i < snps.size(); i++)
        file << to_bim(snps[i]) << std::endl;

    file.close();
}

void plink_binary::read_bim(std::vector<gftools::snp> &snps)
{
    std::ifstream file;
    std::string fn = dataset + ".bim";
    file.open(fn.c_str());
    if (! file.is_open())
        throw gftools::malformed_data("Missing bim file for " + dataset);
    std::string str;
    int i = 0;

    while (getline(file, str)) {
        gftools::snp s = from_bim(str);
        snps.push_back(s);
        snp_index[s.name] = i++;
    }

    file.close();
}

void plink_binary::write_fam(std::vector<gftools::individual> individuals)
{
    std::ofstream file;
    std::string fn = dataset + ".fam";
    file.open(fn.c_str(), std::fstream::out);

    for (unsigned int i = 0; i < individuals.size(); i++)
        file << to_fam(individuals[i]) << std::endl;

    file.close();
}

void plink_binary::read_fam(std::vector<gftools::individual> &individuals)
{
    std::ifstream file;
    std::string fn = dataset + ".fam";
    file.open(fn.c_str());
    if (! file.is_open())
        throw gftools::malformed_data("Missing fam file for " + dataset);
    std::string s;
    int i = 0;

    while (getline(file, s)) {
        gftools::individual ind = from_fam(s);
        individuals.push_back(ind);
        i++;
    }

    file.close();
}

void plink_binary::open_bed_write(std::string filename)
{
    bed_file = new std::fstream();
    bed_file->open(filename.c_str(), std::fstream::out | std::fstream::binary);

    // in theory can use mmap here but you'd have to declare the
    // snps (or at least the snp count) before opening. as most
    // files will likely be opened for reading, (yet) supported
    is_mem_mapped = 0;

    write_bed_header();
}

void plink_binary::open_bed_read(std::string filename, bool quell_mem_mapping)
{
    if ((fd = ::open(filename.c_str(), O_RDONLY)) == -1)
        throw gftools::malformed_data("Missing bed file for " + dataset);
    fd = ::open(filename.c_str(), O_RDONLY);

    struct stat buf;
    if (stat(filename.c_str(), &buf) == -1)
        throw gftools::malformed_data("Missing bed file");
    flen = buf.st_size;

    if (quell_mem_mapping || (fmap = (char *)mmap(0, flen, PROT_READ, MAP_PRIVATE, fd, 0)) == MAP_FAILED) {
        // error mem mapping - use a stream
        is_mem_mapped = 0;
        bed_file = new std::fstream();
        bed_file->open(filename.c_str());
    } else {
        is_mem_mapped = 1;
    }
    read_bed_header();
    snp_ptr = 0;
}

void plink_binary::write_bed_header()
{
    char buffer[MAGIC_LEN];

    for (int i = 0; i < MAGIC_LEN; i++) 
        buffer[i] = magic_number[i];

    if (is_mem_mapped)
        memcpy(fmap, buffer, MAGIC_LEN);
    else
        bed_file->write(buffer, MAGIC_LEN);
}

void plink_binary::read_bed_header()
{
    char buffer[MAGIC_LEN + 1];

    if (is_mem_mapped) {
        memcpy(buffer, fmap, MAGIC_LEN + 1);
    }
    else {
        bed_file->seekg(0, std::ios_base::beg);
        bed_file->read(buffer, MAGIC_LEN + 1);
    }

    for (int i = 0; i < MAGIC_LEN; i++)
        if (buffer[i] != magic_number[i])
            throw gftools::malformed_data("Corrupt or incompatible bed file?");
}

void plink_binary::extract_bed(size_t pos, size_t len, char *buffer)
{
    if (is_mem_mapped) {
        memcpy(buffer, fmap + pos, len);
    } else {
        bed_file->seekg(pos, std::ios_base::beg);
        bed_file->read(buffer, len);
    }
}

void plink_binary::get_snp(size_t pos, std::vector<int> &genotypes)
{
    char *buffer = (char *)malloc(bytes_per_snp);

    extract_bed(pos, bytes_per_snp, buffer);
    genotypes.resize(0);
    uncompress_calls(buffer, individuals.size(), genotypes);
    free(buffer);
}

void plink_binary::write_snp(gftools::snp s, std::vector<int> genotypes) {
    if (genotypes.size() == 0)
        throw gftools::malformed_data("No genotypes defined");
    if (individuals.size() == 0)
        throw gftools::malformed_data("No individuals defined");
    if (genotypes.size() != individuals.size()) {
        std::stringstream ss;
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

void plink_binary::write_snp(gftools::snp s, std::vector<std::string> genotypes)
{
    std::vector<int> g_num;
    genotypes_atoi(s, genotypes, g_num);
    write_snp(s, g_num);
}

void plink_binary::bed_write(gftools::snp s, std::vector<int> genotypes)
{
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
void plink_binary::uncompress_calls(char *buffer, size_t len, std::vector<int> &calls)
{
    char c = buffer[0];
    for (unsigned int i = 0; i < len; i++) {
        if (! (i % 4)) c = buffer[i / 4];  // if div by 4 => get byte
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

void plink_binary::compress_calls(char *buffer, std::vector<int> calls)
{
    int pos = 0;
    unsigned char c = 0;
    for (unsigned int i = 0; i < calls.size(); i++) {
        c >>= 2;
        switch(calls[i]) {
            case 1:  c |= (0 << 6); break;
            case 2:  c |= (2 << 6); break;
            case 3:  c |= (3 << 6); break;
            // treat anything else as a no call;
            // may want to have a strict option here?
            default: c |= (1 << 6);
        }
        // write every forth i
        if ((i % 4) == 3) {
            buffer[pos++] = c;
            c = 0;
        }
    }
    // if not a multiple of four, need to right shift stored data
    // to low end of byte
    if (calls.size() % 4) {
        int unfilled = calls.size();
        while (unfilled++ % 4)
            c >>= 2;
        buffer[pos] = c;
    }
}

