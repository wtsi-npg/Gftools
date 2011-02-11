%module plink_binary

%include "std_vector.i"
%include "std_string.i"

%{
#include "individual.h"
#include "snp.h"
#include "plink_binary.h"
%}

%include "individual.h"
%include "snp.h"
%include "plink_binary.h"

namespace std {
    %template(vectorstr) std::vector<string>;
    %template(vectori) std::vector<int>;
    %template(vectorind) std::vector<gftools::individual>;
    %template(vectorsnp) std::vector<gftools::snp>;
}
