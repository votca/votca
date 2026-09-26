/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

#ifndef VOTCA_TOOLS_FLOATINGPOINTCOMPARISON_H
#define VOTCA_TOOLS_FLOATINGPOINTCOMPARISON_H

// Standard includes
#include <algorithm>
#include <cmath>

/**
 * \brief Provides a means for comparing floating point numbers
 *
 * Implements relative method - do not use for comparing with zero
 * use this most of the time, tolerance needs to be meaningful in your context
 *
 * Function taken from
 * https://stackoverflow.com/questions/17333/what-is-the-most-effective-way-for-float-and-double-comparison
 * user ShitalShal
 */

// Standard includes
#include <cstdlib>

namespace votca {
namespace tools {

// NOT `static`. A `static` function template in a header has internal
// linkage, so every translation unit that includes the header gets its
// own private copy -- and in every TU that includes the header without
// calling the function, that copy is dead code. Clang reports it under
// -Wunused-template, which newer releases enable as part of the general
// warning set this project already builds with, and -DENABLE_WERROR=ON
// turns into a build failure. The first casualty is csg's
// lammpsdatareader.cc, whose #include of this header is left over: it
// calls nothing from it.
//
// Dropping `static` is not a workaround for the warning, it is the
// correct declaration: a function template needs neither `static` nor
// `inline` to be defined in a header. Implicit instantiations are
// already exempt from the one-definition rule, so several TUs
// instantiating the same specialisation share one copy rather than
// colliding.
template <typename T>
bool isApproximatelyEqual(T a, T b, T tolerance) {
  T diff = std::abs(a - b);
  if (diff <= tolerance) {
    return true;
  }

  if (diff < std::max(std::abs(a), std::abs(b)) * tolerance) {
    return true;
  }
  return false;
}

}  // namespace tools
}  // namespace votca
#endif  // VOTCA_TOOLS_FLOATINGPOINTCOMPARISON_H
