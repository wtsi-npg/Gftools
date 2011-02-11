#ifndef _SNP_H_
#define _SNP_H_

namespace gftools
{
    class snp {
        public:
            int genetic_position, physical_position;
            std::string name;
            std::string chromosome;
            std::string allele_a, allele_b;
        snp(std::string name) {
            this->name = name;
            genetic_position = 0;
            physical_position = 0;
        }
        snp() {
            genetic_position = 0;
            physical_position = 0;
        }
    };
}

#endif // _SNP_H_
