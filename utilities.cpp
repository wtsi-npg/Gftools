
#include <cerrno>
#include <cstring>
#include <fstream>
#include <iostream>

#include "utilities.h"

using std::ifstream;
using std::string;

namespace gftools {

    const string error_message() {
        char *msg = strerror(errno);
        return msg ? string(msg) : "unknown error";
    }

    bool at_eof(ifstream &ifstream) {
        return ifstream.peek() == ifstream::traits_type::eof();
    }
}
