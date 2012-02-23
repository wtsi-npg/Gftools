#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

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

#endif // INDIVIDUAL_H