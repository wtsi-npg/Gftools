#include "plink_utils.h"
#include "plink_binary.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>

void shift_right(char *in, size_t len, int count);

bool plink_utils::qmerge(std::vector<std::string> inputs, std::string output) {
    if (inputs.size() < 2)
        return true;

    std::vector<plink_binary> pb;
    for (size_t i = 0; i < inputs.size(); i++) {
        plink_binary p(inputs[i]);
        pb.push_back(p);
    }
    plink_binary merged;
    merged.open(output, 1);

    bool err = qmerge(pb, merged);

    merged.close();
    for (size_t i = 0; i < pb.size(); i++)
        pb[i].close();

    return err;
}

bool plink_utils::qmerge(std::vector<plink_binary> pb, plink_binary &merged) {
    if (pb.size() < 2)
        return true;

    bool err = false;

    /*
     * compare snps/samples in first two files to check whether
     * we have SNPs or samples in common, but not both
     * (this method is not intended to cover such cases)
     */
    bool inds_match = ind_cmp(pb[0].individuals, pb[1].individuals);
    bool snps_match = snp_cmp(pb[0].snps, pb[1].snps);
    bool merge_snps = false;

    if (snps_match && inds_match) {
        throw gftools::malformed_data(pb[0].dataset + " and " + pb[1].dataset + " are identical");
        err = true;
    }
    if (!(snps_match || inds_match)) {
        throw gftools::malformed_data("First two data sets have neither SNPs or samples in common");
        err = true;
    }
    if (inds_match)
        merge_snps = true;

    for (size_t i = 2; i < pb.size(); i++) {
        if (merge_snps) {
            if (!ind_cmp(pb[0].individuals, pb[i].individuals)) {
                throw gftools::malformed_data(pb[i].dataset + ": mismatching individuals");
                err = true;
            }
        }
        else {
            if (!snp_cmp(pb[0].snps, pb[i].snps)) {
                throw gftools::malformed_data(pb[i].dataset + ": mismatching SNPs");
                err = true;
            }
        }
    }

    if (!err) {
        if (merge_snps) {
            plink_utils::qmerge_snps(pb, merged);
        }
        else {
            // plink_utils::qmerge_samples(pb, merged);
        }
    }

    return err;
}

void plink_utils::qmerge_snps(std::vector<plink_binary> pb, plink_binary &merged) {
    for (size_t j = 0; j < pb[0].individuals.size(); j++) {
        merged.individuals.push_back(pb[0].individuals[j]);
    }
    for (size_t i = 0; i < pb.size(); i++) {
        for (size_t j = 0; j < pb[i].snps.size(); j++) {
            merged.snps.push_back(pb[i].snps[j]);
        }
    }
    size_t len = pb[0].bytes_per_snp;
    char *buffer = (char *)malloc(len);
    for (size_t i = 0; i < pb.size(); i++) {
        for (size_t j = 0; j < pb[i].snps.size(); j++) {
            pb[i].extract_bed(3 + j * len, len, buffer);
            merged.bed_file->write(buffer, len);
        }
    }
    free(buffer);

    return;
}

bool plink_utils::ind_cmp(std::vector<gftools::individual> ind1, std::vector<gftools::individual> ind2) {
    if (ind1.size() != ind2.size())
        return false;

    for(size_t i = 0; i < ind1.size(); i++) {
        if (ind1[i].name.compare(ind2[i].name))
            return false;

        if (ind1[i].family.compare(ind2[i].family))
            return false;
    }
    return true;
}

bool plink_utils::snp_cmp(std::vector<gftools::snp> snp1, std::vector<gftools::snp> snp2) {
    if (snp1.size() != snp2.size())
        return false;

    for(size_t i = 0; i < snp1.size(); i++) {
        if (snp1[i].name.compare(snp2[i].name))
            return false;
    }
    return true;
}