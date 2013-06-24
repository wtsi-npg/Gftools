/*
 * Copyright (c) 2011-2012 Genome Research Ltd. All rights reserved.
 *
 * This file is part of Gftools.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "utilities.h"
#include "plink_binary.h"

using std::fstream;
using std::ifstream;
using std::ofstream;
using std::string;
using std::stringstream;
using std::vector;

using gftools::individual;
using gftools::snp;
using gftools::error_message;
using gftools::at_eof;

// final '1' for snp major mode; only supporting this at present
#define MAGIC_LEN 3
static int magic_number[MAGIC_LEN] = { 108, 27, 1 };

const char DEFAULT_MISSING_ALLELE = 'N';

plink_binary::plink_binary(void) {
    // default "NN" for no call
    missing_genotype = DEFAULT_MISSING_ALLELE;
}

plink_binary::plink_binary(string dataset) {
    // default "NN" for no call
    missing_genotype = DEFAULT_MISSING_ALLELE;
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
    quell_mem_mapping = false;
    if (mode) {
        open_for_write = 1;
        open_bed_write(dataset + ".bed");
    }
    else {
        open_for_write = 0;
        read_bim(snps);
        if (snps.empty()) {
            throw gftools::malformed_data("No SNPs read");
        }
        read_fam(individuals);
        if (individuals.empty()) {
            throw gftools::malformed_data("No individuals read");
        }
        bytes_per_snp = (3 + individuals.size()) / 4;
        open_bed_read(dataset + ".bed", quell_mem_mapping);
    }
}

bool plink_binary::next_snp(snp &snp, vector<string> &genotypes) {
    vector<int> gt_int;
    if (snp_ptr >= snps.size()) {
        return false;
    }

    get_snp(3 + snp_ptr * bytes_per_snp, gt_int);
    genotypes_itoa(snps[snp_ptr], gt_int, genotypes);
    snp = snps[snp_ptr++];
    return true;
}

void plink_binary::read_snp(string snp, vector<string> &genotypes) {
    vector<int> gt_int;
    int index = snp_index[snp];
    read_snp(3 + index * bytes_per_snp, gt_int);
    genotypes_itoa(snps[index], gt_int, genotypes);
    snp_ptr = index + 1;
}

void plink_binary::genotypes_itoa(snp snp, const vector<int> g_num, vector<string> &g_str) {
    g_str.resize(0);

    for (vector<int>::const_iterator i = g_num.begin(); i < g_num.end(); i++) {
        string genotype;
        string missing = string(2, missing_genotype);
        switch (*i) {
            case 0: genotype = missing; break;
            case 1: genotype = snp.allele_a + snp.allele_a; break;
            case 2: genotype = snp.allele_a + snp.allele_b; break;
            case 3: genotype = snp.allele_b + snp.allele_b; break;
            default: throw gftools::malformed_data("integer genotypes must be 0..3");
        }
        g_str.push_back(genotype);
    }
}

void plink_binary::genotypes_atoi(gftools::snp &snp, const vector<string> g_str, vector<int> &g_num) {
    vector<string> alleles = collate_alleles(g_str);
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
    vector<string> alleles;

    if (!g_str.empty()) {
        string allele_a = "0";
        string allele_b = "0";
        std::set<string> alleles_ab;
        string missing = string(1, missing_genotype);

        vector<string>::const_iterator gi;
        for (gi = g_str.begin(); gi < g_str.end(); gi++) {
            string call = *gi;

            if (call.length() == 2) {
                string a = call.substr(0, 1);
                string b = call.substr(1, 2);

                if (a != missing) {
                    if (allele_a == "0") {
                        allele_a = a;
                    }
                    else if (allele_a != a && allele_a != allele_b && a != b) {
                        throw gftools::malformed_data("Ambiguous allele A in " +
                                                      call + " of " +
                                                      call_str(g_str));
                    }
                    alleles_ab.insert(a);
                }

                if (b != missing) {
                    if (allele_b == "0") {
                        allele_b = b;
                    }
                    else if (allele_b != b && allele_a != allele_b && a != b) {
                        throw gftools::malformed_data("Ambiguous allele B in " +
                                                      call + " of " +
                                                      call_str(g_str));
                    }
                    alleles_ab.insert(b);
                }
            }
            else {
                throw gftools::malformed_data("Unrecognised genotype " + call);
            }
        }

        if (alleles_ab.size() > 2) {
             std::stringstream ss;
             ss << "[";

             std::set<string>::const_iterator ai;
             for (ai = alleles_ab.begin(); ai != alleles_ab.end(); ai++) {
               if (ai != alleles_ab.begin()) {
                 ss << ", ";
               }
               ss << *ai;
             }
             ss << "]";

            throw gftools::malformed_data("Collated genotype was not bi-allelic: " +
                                          ss.str());
        }

        alleles.push_back(allele_a);

        // This seems arbitrary. If you have only homozygous calls,
        // how do you know whether this is allele A or B, without
        // reference to the annotation? I'm copying the implementation
        // in the original code, but I
        if (alleles_ab.size() == 1) {
            alleles.push_back("0");
        }
        else {
            alleles.push_back(allele_b);
        }
    }

    return alleles;
}

void plink_binary::read_snp(int snp_index, vector<int> &genotypes) {
    get_snp(3 + snp_index * bytes_per_snp, genotypes);
}

snp plink_binary::from_bim(string record) {
    stringstream ss(record);
    snp snp;
    ss >> snp.chromosome;
    ss >> snp.name;
    ss >> snp.genetic_position;
    ss >> snp.physical_position;
    ss >> snp.allele_a;
    ss >> snp.allele_b;
    return snp;
}

string plink_binary::to_bim(snp snp) {
    stringstream record;
    record << (snp.chromosome.length() ? snp.chromosome : "0");
    record << "\t" << snp.name;
    record << "\t" << (snp.genetic_position ? snp.genetic_position : 0);
    record << "\t" << (snp.physical_position ? snp.physical_position : 0);
    record << "\t" << snp.allele_a;
    record << "\t" << snp.allele_b;
    return record.str();
}

individual plink_binary::from_fam(string record) {
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

string plink_binary::to_fam(individual ind) {
    stringstream record;
    record << (ind.family.length() ? ind.family : ind.name);
    record << "\t" << ind.name;
    record << "\t" << (ind.father.length() ? ind.father : "-9");
    record << "\t" << (ind.mother.length() ? ind.mother : "-9");
    record << "\t" << (ind.sex.length() ? ind.sex : "-9");
    record << "\t" << (ind.phenotype.length() ? ind.phenotype : "-9");
    return record.str();
}

void plink_binary::write_bim(const vector<snp> snps) {
    ofstream file;
    string fn = dataset + ".bim";
    file.open(fn.c_str(), fstream::out);

    for (unsigned int i = 0; i < snps.size(); i++) {
        file << to_bim(snps[i]) << std::endl;
    }

    file.close();
}

void plink_binary::read_bim(vector<snp> &snps) {
    ifstream file;
    string fn = dataset + ".bim";
    file.open(fn.c_str());

    if (!file) {
        throw gftools::malformed_data("Failed to open BIM file for " + dataset +
                                      ": " + error_message());
    }
    else {
        if (at_eof(file)) {
            throw gftools::malformed_data("Empty BIM file for " + dataset);
        }
    }

    string str;
    int i = 0;

    while (getline(file, str)) {
        snp snp = from_bim(str);
        snps.push_back(snp);
        snp_index[snp.name] = i++;
    }

    file.close();
}

void plink_binary::write_fam(const vector<individual> individuals) {
    ofstream file;
    string fn = dataset + ".fam";
    file.open(fn.c_str(), fstream::out);

    for (unsigned int i = 0; i < individuals.size(); i++) {
        file << to_fam(individuals[i]) << std::endl;
    }

    file.close();
}

void plink_binary::read_fam(vector<individual> &individuals) {
    ifstream file;
    string fn = dataset + ".fam";
    file.open(fn.c_str());
    if (!file.is_open()) {
        throw gftools::malformed_data("Missing FAM file for " + dataset);
    }

    string str;
    int i = 0;
    while (getline(file, str)) {
        individual ind = from_fam(str);
        individuals.push_back(ind);
        i++;
    }

    file.close();
}

void plink_binary::open_bed_write(string filename) {
    bed_file = new fstream();
    bed_file->open(filename.c_str(), fstream::out | fstream::binary);

    // in theory can use mmap here but you'd have to declare the
    // snps (or at least the snp count) before opening. as most
    // files will likely be opened for reading, (yet) supported
    is_mem_mapped = 0;

    write_bed_header();
}

void plink_binary::open_bed_read(string filename, bool quell_mem_mapping) {
    if (quell_mem_mapping) {
        is_mem_mapped = 0;
        bed_file = new fstream();
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

void plink_binary::write_snp(snp snp, const vector<int> genotypes) {
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
        ss << snp.name;
        ss << " whereas ";
        ss << individuals.size();
        ss << " individuals defined.";
        throw gftools::malformed_data(ss.str());
    }
    bed_write(snp, genotypes);
}

void plink_binary::write_snp(snp snp, vector<string> genotypes) {
    vector<int> g_num;
    genotypes_atoi(snp, genotypes, g_num);
    write_snp(snp, g_num);
}

void plink_binary::bed_write(snp snp, vector<int> genotypes) {
    int len = (3 + individuals.size()) / 4;
    char buffer[len];
    compress_calls(buffer, genotypes);
    snps.push_back(snp);
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

string plink_binary::call_str(const vector<string> &g_str) {
    std::stringstream ss;
    ss << "[";

    vector<string>::const_iterator gi;
    for (gi = g_str.begin(); gi < g_str.end(); gi++) {
        if (gi != g_str.begin()) {
            ss << ", ";
        }
        ss << *gi;
    }
    ss << "]";

    return ss.str();
}
