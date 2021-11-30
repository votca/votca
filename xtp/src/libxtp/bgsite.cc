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

#include "votca/xtp/bgsite.h"

namespace votca {
namespace xtp {

BGSite::BGSite(const PolarSite& pol) {
  id_ = pol.getId();
  position_ = pol.getPos();
  element_ = pol.getElement();

  // Multipole data
  mp_.rank = pol.getRank();
  mp_.charge = pol.getCharge();
  mp_.dipole = pol.getStaticDipole();
  // the 1/3 saves factors in further calculations.The quadrupole in the BG site
  // should by multiplied by 3 to obtain the true quadrupole
  mp_.quadrupole = (1.0 / 3.0) * pol.CalculateCartesianMultipole();
  induced_dipole_ = pol.getInducedDipole();
  polarization_ = pol.getpolarization();
}

}  // namespace xtp
}  // namespace votca