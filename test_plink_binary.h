#ifndef TEST_PLINK_BINARY_H
#define TEST_PLINK_BINARY_H

#include <stdio.h>
#include <cxxtest/TestSuite.h>
#include "plink_binary.h"

using std::string;
using std::vector;

using gftools::individual;
using gftools::snp;

class ReadTest : public CxxTest::TestSuite {
public:
    void testOpen() {
        // Implicit open for reading
        plink_binary pb = plink_binary("data");
        TS_ASSERT_EQUALS(4, pb.individuals.size());
        pb.close();
    }

    void testOpenEmpty() {
        // All these Plink files are present, but empty
        plink_binary pb = plink_binary();
        TS_ASSERT_THROWS_ANYTHING(pb.open("empty"));
    }

    void testOpenMissing() {
        // This file does not exist
        plink_binary pb = plink_binary();
        TS_ASSERT_THROWS_ANYTHING(pb.open("no such dataset"));
    }

    void testOpenWrite() {
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

            // Write a new BED dataset
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

};

#endif