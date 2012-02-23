#ifndef SNP_H
#define SNP_H

namespace gftools {
    /** A SNP.
     *
     * A named, bi-allelic SNP with chromsomal location.
     */
    class snp {
public:
        /// The position of the SNP
        int genetic_position, physical_position;
        /// The SNP name.
        std::string name;
        /// The name of the chromosome on which the SNP is located.
        std::string chromosome;
        // The alleles of the SNP
        std::string allele_a, allele_b;

        /** Creates a new, named SNP and genetic and physical positions = 0.
         *
         * @param name The SNP name.
         */
        snp(std::string name) {
            this->name = name;
            genetic_position = 0;
            physical_position = 0;
        }

        /** Creates a new, unnamed SNP and genetic and physical positions = 0.
         */
        snp() {
            genetic_position = 0;
            physical_position = 0;
        }

    };

}

#endif // SNP_H