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

// Local VOTCA includes
#include "votca/xtp/ecpaoshell.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/checkpointtable.h"
#include "votca/xtp/ecpaobasis.h"

namespace votca {
namespace xtp {

void ECPAOGaussianPrimitive::WriteData(data& d, const ECPAOShell& shell) const {
  d.atomid = shell.getAtomIndex();
  d.l = static_cast<Index>(shell.getL());
  d.startindex = shell.getStartIndex();
  d.decay = getDecay();
  d.contraction = getContraction();
  d.power = _power;
  d.x = shell.getPos().x();
  d.y = shell.getPos().y();
  d.z = shell.getPos().z();
  d.lmax = static_cast<Index>(shell.getLmaxElement());
}

void ECPAOGaussianPrimitive::SetupCptTable(CptTable& table) {

  table.addCol<Index>("atomidx", HOFFSET(data, atomid));
  table.addCol<Index>("L", HOFFSET(data, l));
  table.addCol<Index>("startidx", HOFFSET(data, startindex));
  table.addCol<Index>("power", HOFFSET(data, power));
  table.addCol<double>("decay", HOFFSET(data, decay));
  table.addCol<double>("contr", HOFFSET(data, contraction));
  table.addCol<double>("pos.x", HOFFSET(data, x));
  table.addCol<double>("pos.y", HOFFSET(data, y));
  table.addCol<double>("pos.z", HOFFSET(data, z));
  table.addCol<Index>("lmax", HOFFSET(data, lmax));
}

std::ostream& operator<<(std::ostream& out, const ECPAOShell& shell) {
  out << "AtomIndex:" << shell.getAtomIndex();
  out << " Shelltype:" << xtp::EnumToString(shell.getL())
      << " L:" << Index(shell.getL()) << " NonLocal:" << shell.isNonLocal()
      << " Func:" << shell.getNumFunc() << "\n";
  for (const auto& gaussian : shell) {
    out << " Gaussian Decay: " << gaussian.getDecay();
    out << " Power: " << gaussian.getPower();
    out << " Contractions: " << gaussian.getContraction() << "\n";
  }
  return out;
}

}  // namespace xtp
}  // namespace votca
