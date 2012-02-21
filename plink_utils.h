#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include "plink_binary.h"

class plink_utils {
private:
    bool snp_cmp(std::vector<gftools::snp> snp1, std::vector<gftools::snp> snp2);

    bool ind_cmp(std::vector<gftools::individual> ind1, std::vector<gftools::individual> ind2);

    void qmerge_snps(std::vector<plink_binary> pb, plink_binary &merged);

    // void qmerge_samples(std::vector<plink_binary> pb, plink_binary &merged);
public:
    bool qmerge(std::vector<plink_binary> inputs, plink_binary &output);

    bool qmerge(std::vector<std::string> inputs, std::string output);

};

#endif // UTILS_H