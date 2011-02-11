#include <stdexcept>

namespace gftools
{
    class MalformedData : public std::runtime_error
    {
        public:
        MalformedData()
            : runtime_error("Malformed data") {}
        MalformedData(std::string info)
            : runtime_error(info) {}
    };
}
