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

#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/staticregion.h>

namespace votca {
namespace xtp {

void PolarRegion::Initialize(const tools::Property& prop) {

  std::string filename =
      prop.ifExistsReturnElseThrowRuntimeError<std::string>("options");
  tools::Property polar_xml;
  tools::load_property_from_xml(polar_xml, filename);
  _max_iter = polar_xml.ifExistsReturnElseReturnDefault("max_iter", _max_iter);
  _deltaE = polar_xml.ifExistsReturnElseReturnDefault("tolerance", _deltaE);
  _exp_damp = polar_xml.ifExistsReturnElseReturnDefault("exp_damp", _exp_damp);
  _induce_intra_mol = polar_xml.ifExistsReturnElseReturnDefault(
      "induce_intra_molecule", _induce_intra_mol);

  return;
}

bool PolarRegion::Converged() const { return false; }

void PolarRegion::Evaluate() {}

void PolarRegion::ResetRegion() { return; }
void PolarRegion::InteractwithQMRegion(QMRegion& region) { return; }
void PolarRegion::InteractwithPolarRegion(PolarRegion& region) { return; }
void PolarRegion::InteractwithStaticRegion(StaticRegion& region) { return; }

}  // namespace xtp
}  // namespace votca
