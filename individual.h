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

#ifndef GFTOOLS_INDIVIDUAL_H
#define GFTOOLS_INDIVIDUAL_H

namespace gftools {
    /** An individual whose genotype has been determined.
     */
    class individual {
public:
        std::string family;
        /// An unique identifier for this individual.
        std::string name;
        /// The identifier of the father individual.
        std::string father;
        /// The identifier of the mother individual.
        std::string mother;
        std::string sex;
        std::string phenotype;
        individual() {};
        individual(std::string family, std::string name,
                   std::string father, std::string mother,
                   std::string sex, std::string phenotype) {
            this->family = family;
            this->name = name;
            this->father = father;
            this->mother = mother;
            this->sex = sex;
            this->phenotype = phenotype;
        }

    };

}

#endif // GFTOOLS_INDIVIDUAL_H