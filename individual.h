#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H

namespace gftools {
    class individual {
public:
        std::string family;
        std::string name;
        std::string father;
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