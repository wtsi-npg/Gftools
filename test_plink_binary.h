#ifndef TEST_PLINK_BINARY_H
#define TEST_PLINK_BINARY_H

#include <cstdio>
#include <map>
#include <fstream>
#include <sstream>

#include <cxxtest/TestSuite.h>
#include "plink_binary.h"

using std::ifstream;
using std::string;
using std::vector;

using gftools::individual;
using gftools::snp;

class ReadTest : public CxxTest::TestSuite {

    vector<string> expected_a;
    vector<string> expected_b;
    vector<string> expected_snp;
    vector<string> expected_ind;
    std::map<string, vector<string> > expected_gen;

public:
    void setUp() {
        expected_a = vector<string>();
        string ea[4] = {"0", "A", "A", "A"};
        for (unsigned int i = 0; i < 4; i++) {
            expected_a.push_back(string(ea[i]));
        }

        string eb[4] = {"0", "0", "C", "G"};
        for (unsigned int i = 0; i < 4; i++) {
            expected_b.push_back(string(eb[i]));
        }

        string esn[4] = {"rs1000", "rs1001", "rs1002", "rs1003"};
        string geno[][4] = {{"00", "00", "00", "00"},
                            {"AA", "AA", "AA", "AA"},
                            {"AC", "CC", "AC", "CC"},
                            {"GG", "GG", "AG", "GG"}};

        expected_gen = std::map<string, vector<string> >();
        for (unsigned int i = 0; i < 4; i++) {
            string name = string(esn[i]);

            vector<string> genotypes = vector<string>();
            for (int j = 0; j < 4; j++) {
                genotypes.push_back(string(geno[i][j]));
            }

            expected_snp.push_back(name);
            expected_gen[name] = genotypes;
        }

        string eid[] = {"sample_000", "sample_001", "sample_002", "sample_003"};
        for (unsigned int i = 0; i < 4; i++) {
            expected_ind.push_back(string(eid[i]));
        }
    }

    void test_open_read() {
        // Implicit open for reading
        plink_binary pb = plink_binary("data");
        pb.missing_genotype = '0';

        TS_ASSERT_EQUALS(4, pb.individuals.size());
        for (int i = 0; i < 4; i++) {
            TS_ASSERT_EQUALS(expected_ind[i], pb.individuals[i].name);
        }

        vector<string> genotypes;
        snp snp;
        for (int i = 0; i < 4; i++) {
            TS_ASSERT(pb.next_snp(snp, genotypes));
            TS_ASSERT_EQUALS(expected_a[i], snp.allele_a);
            TS_ASSERT_EQUALS(expected_b[i], snp.allele_b);

            vector<string> expected_genotypes = expected_gen[snp.name];
            TS_ASSERT_EQUALS(genotypes.size(), expected_genotypes.size());

            for (int j = 0; j < 4; j++) {
                TS_ASSERT_EQUALS(expected_genotypes[j], genotypes[j]);
            }
        }
        TS_ASSERT(!pb.next_snp(snp, genotypes));

        pb.close();
    }

    void test_open_empty() {
        // All these Plink files are present, but empty
        plink_binary pb = plink_binary();
        TS_ASSERT_THROWS_ANYTHING(pb.open("empty"));
    }

    void test_open_missing() {
        // This file does not exist
        plink_binary pb = plink_binary();
        TS_ASSERT_THROWS_ANYTHING(pb.open("no such dataset"));
    }

    void test_open_write() {
        char *tmpname = NULL;
        tmpname = tmpnam(NULL);
        if (!tmpname) {
            TS_FAIL("Failed to create a temporary file name");
        }
        else {
            plink_binary pbi = plink_binary("data");
            plink_binary pbo = plink_binary();

            string tmpfile = string(tmpname);
            pbo.open(tmpfile, true);

            vector<string> genotypes_in;
            snp snp_in;
            pbo.individuals = pbi.individuals;
            pbo.snp_index = pbi.snp_index;

            // Write a new BED dataset copy
            while (pbi.next_snp(snp_in, genotypes_in)) {
                pbo.write_snp(snp_in, genotypes_in);
            }
            pbi.close();
            pbo.close();

            pbi = plink_binary("data");
            pbo = plink_binary(tmpfile);

            // Compare individuals
            vector<individual>::iterator ii;
            vector<individual>::iterator io = pbo.individuals.begin();
            for (ii = pbi.individuals.begin(); ii < pbi.individuals.end(); ii++) {
                TS_ASSERT_EQUALS((*ii).name, (*io).name);
                io++;
            }

            vector<string> genotypes_out;
            snp snp_out;

            while (pbi.next_snp(snp_in, genotypes_in)) {
                TS_ASSERT(pbo.next_snp(snp_out, genotypes_out));

                // Compare SNPs
                TS_ASSERT_EQUALS(snp_in.name, snp_out.name);
                TS_ASSERT_EQUALS(snp_in.chromosome, snp_out.chromosome);
                TS_ASSERT_EQUALS(snp_in.allele_a, snp_out.allele_a);
                TS_ASSERT_EQUALS(snp_in.allele_b, snp_out.allele_b);

                // Compare genotypes
                vector<string>::iterator gi;
                vector<string>::iterator go = genotypes_out.begin();
                for (gi = genotypes_in.begin(); gi < genotypes_in.end(); gi++) {
                    TS_ASSERT_EQUALS(*gi, *go);
                    go++;
                }

                TS_ASSERT(go == genotypes_out.end());
            }

            vector<string> suffixes = vector<string>(3);
            suffixes[0] = ".bed";
            suffixes[1] = ".bim";
            suffixes[2] = ".fam";

            for (vector<string>::iterator fi = suffixes.begin(); fi < suffixes.end(); fi++) {
                string fn = tmpfile + *fi;
                if (remove(fn.c_str()) != 0) {
                    TS_FAIL("Failed to remove tmp file " + fn);
                }
            }
        }
    }

    void test_write_snp() {
        snp snp;
        snp.name = "test_snp";
        snp.allele_a = "A";
        snp.allele_b = "C";

        vector<individual> inds = vector<individual>();
        for (int i = 0; i < 4; i++) {
            individual ind;
            std::stringstream name;
            name << "sample_00" << i;
            ind.name = name.str();
            inds.push_back(ind);
        }

        vector<string> bad_gt = vector<string>();
        bad_gt.push_back("GG");
        bad_gt.push_back("AG");
        bad_gt.push_back("GG");
        bad_gt.push_back("GA");

        vector<string> good_gt = vector<string>();
        good_gt.push_back("AA");
        good_gt.push_back("AC");
        good_gt.push_back("CA");
        good_gt.push_back("NN");

        char *tmpname = NULL;
        tmpname = tmpnam(NULL);
        if (!tmpname) {
            TS_FAIL("Failed to create a temporary file name");
        }
        else {
            string tmpfile = string(tmpname);
            plink_binary pb = plink_binary();
            pb.open(tmpfile, true);
            pb.individuals = inds;

            // Should fail to write all bad genotypes
            vector<string>::iterator i;
            for (i = bad_gt.begin(); i < bad_gt.end(); i++) {
                vector<string> tmp = vector<string>(1, *i);
                TS_ASSERT_THROWS_ANYTHING(pb.write_snp(snp, tmp));
            }
            pb.write_snp(snp, good_gt);
            pb.close();

            // Check written genotypes
            pb = plink_binary(tmpfile);
            vector<string> genotypes;
            TS_ASSERT(pb.next_snp(snp, genotypes));

            vector<string> expected_gt = vector<string>();
            expected_gt.push_back("AA");
            expected_gt.push_back("AC");
            expected_gt.push_back("AC"); // Note, was CA
            expected_gt.push_back("NN");

            vector<string>::iterator ie;
            vector<string>::iterator io = genotypes.begin();
            for (ie = expected_gt.begin(); ie < expected_gt.end(); ie++) {
                TS_ASSERT_EQUALS(*ie, *io);
                io++;
            }
            TS_ASSERT(io == genotypes.end());
            pb.close();

            // Check the SNP annotation
            string bim = tmpfile + ".bim";
            ifstream bim_file;
            bim_file.open(bim.c_str(), std::fstream::in);
            TS_ASSERT(bim_file.is_open());

            string line;
            std::getline(bim_file, line);
            TS_ASSERT_EQUALS("0	test_snp	0	0	A	C", line);
            bim_file.close();
        }
    }

};

#endif
