#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

#include <stdexcept>

namespace gftools
{
    class malformed_data : public std::runtime_error
    {
        public:
        malformed_data()
            : runtime_error("Malformed data") {}
        malformed_data(std::string info)
            : runtime_error(info) {}
    };
}

#endif // _EXCEPTIONS_H_
