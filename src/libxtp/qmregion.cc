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

void QMRegion::Initialize(const tools::Property& prop) { return; }

bool QMRegion::Converged() const { return true; }

void QMRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {
  return;
}

void QMRegion::ResetRegion() { return; }
void QMRegion::InteractwithQMRegion(QMRegion& region) { return; }
void QMRegion::InteractwithPolarRegion(PolarRegion& region) { return; }
void QMRegion::InteractwithStaticRegion(StaticRegion& region) { return; }

void QMRegion::WritePDB(csg::PDBWriter& writer) const {
  writer.WriteContainer(_orb.QMAtoms());
}

void QMRegion::WriteToCpt(CheckpointWriter& w) const {
  w(_id, "id");
  w(identify(), "type");
  CheckpointWriter v = w.openChild("orbitals");
  _orb.WriteToCpt(v);
}

void QMRegion::ReadFromCpt(CheckpointReader& r) {
  r(_id, "id");
  CheckpointReader rr = r.openChild("orbitals");
  _orb.ReadFromCpt(rr);
}

}  // namespace xtp
}  // namespace votca
