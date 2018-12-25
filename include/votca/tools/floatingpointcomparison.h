/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#ifndef _VOTCA_TOOLS_FLOATINGPOINTCOMPARISON_H
#define _VOTCA_TOOLS_FLOATINGPOINTCOMPARISON_H

#include <algorithm>
#include <cmath>

/**
 * \brief Provides a means for comparing floating point numbers
 *
 * Implements relative method - do not use for comparing with zero
 * use this most of the time, tolerance needs to be meaningful in your context
 *
 * Function taken from https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
 * user ShitalShal
 */
namespace votca {
namespace tools {

template<typename T>
static bool isApproximatelyEqual(T a, T b, T tolerance )
{
    T diff = std::abs(a - b);
    if (diff <= tolerance)
        return true;

    if (diff < std::max(std::abs(a), std::abs(b)) * tolerance)
        return true;
    return false;
}

}
}
#endif // _VOTCA_TOOLS_FLOATINGPOINTCOMPARISON_H
