#ifndef SNP_H
#define SNP_H

namespace gftools {
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

#endif // SNP_H