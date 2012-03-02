#ifndef GFTOOLS_UTILITIES_H
#define GFTOOLS_UTILITIES_H

#include <fstream>

namespace gftools {

    /** Returns the current C error message.
     *
     * @returns The error message or 'unknown error',
     */
    const std::string error_message();

    /** Returns true if the next element in the stream is eof.
     */
    bool at_eof(std::ifstream &ifstream);
}

#endif // GFTOOLS_UTILITIES_H
