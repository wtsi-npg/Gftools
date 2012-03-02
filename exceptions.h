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

#ifndef GFTOOLS_EXCEPTIONS_H
#define GFTOOLS_EXCEPTIONS_H

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

#endif // GFTOOLS_EXCEPTIONS_H