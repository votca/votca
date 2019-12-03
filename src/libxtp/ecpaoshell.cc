/*
 *            Copyright 2009-2019 The VOTCA Development Team
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
#include "votca/xtp/ecpaoshell.h"
#include "votca/xtp/ecpaobasis.h"
#include <votca/xtp/aomatrix.h>

namespace votca {
namespace xtp {

std::ostream& operator<<(std::ostream& out, const ECPAOShell& shell) {
  out << "AtomIndex:" << shell.getAtomIndex();
  out << " Shelltype:" << shell.getType() << " L:" << shell.getL()
      << " NonLocal:" << shell.isNonLocal() << " Func:" << shell.getNumFunc()
      << "\n";
  for (const auto& gaussian : shell) {
    out << " Gaussian Decay: " << gaussian.getDecay();
    out << " Power: " << gaussian.getPower();
    out << " Contractions: " << gaussian.getContraction() << "\n";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
