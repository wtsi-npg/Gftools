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

#ifndef GFTOOLS_SNP_H
#define GFTOOLS_SNP_H

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

        /** Returns true if both alleles of the SNP have been assigned.
         */
        bool is_known() {
            return (!allele_a.empty() || !allele_b.empty());
        }

    };

}

#endif // GFTOOLS_SNP_H