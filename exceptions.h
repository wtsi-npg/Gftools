#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <stdexcept>

namespace gftools {
    /** An exception raised for all types of data error.
     */
    class malformed_data : public std::runtime_error {
public:
        malformed_data()
            : runtime_error("Malformed data") {}
        malformed_data(std::string info)
            : runtime_error(info) {}
    };

}

#endif // EXCEPTIONS_H