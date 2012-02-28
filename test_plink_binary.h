#ifndef TEST_PLINK_BINARY_H
#define TEST_PLINK_BINARY_H

#include <cxxtest/TestSuite.h>
#include "plink_binary.h"

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
    
//    void testOpenWrite() {
//        plink_binary pbi = plink_binary("data");
//        plink_binary pbo = plink_binary();
//        pbo.open("test_write", true);
//        
//        std::vector<std::string> genotypes;
//        gftools::snp snp;
//        pbo.individuals = pbi.individuals;
//        pbo.snp_index = pbi.snp_index;
//        
//        while (pbi.next_snp(snp, genotypes)) {
//            pbo.write_snp(snp, genotypes);
//        }
//        pbo.close();
//    }

};

#endif