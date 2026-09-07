/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#include "votca/xtp/ewaldshapecorrection.h"

namespace votca {
namespace xtp {

namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

EwaldShapeCorrection::EwaldShapeCorrection(double volume,
                                           const EwaldRegistry& registry,
                                           EwaldShape shape)
    : volume_(volume), registry_(registry), shape_(shape) {}

Eigen::Vector3d EwaldShapeCorrection::TotalDipoleMoment(
    EwaldChargeState source_state) const {
  Eigen::Vector3d M = Eigen::Vector3d::Zero();
  for (Index id : registry_.AllIds()) {
    if (!registry_.Has(id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(id, source_state);
    for (const PolarSite& site : segment) {
      M += site.getCharge() * site.getPos();
      M += site.getStaticDipole();
      M += site.getInducedDipole();
    }
  }
  return M;
}

template <enum Estatic CE>
void EwaldShapeCorrection::AddFieldAt(PolarSite& target,
                                      EwaldChargeState source_state) const {
  const Eigen::Vector3d M = TotalDipoleMoment(source_state);

  Eigen::Vector3d field;
  if (shape_ == EwaldShape::Cube) {
    field = -(4.0 * kPi / (3.0 * volume_)) * M;
  } else {
    field = Eigen::Vector3d::Zero();
    field.z() = -(4.0 * kPi / volume_) * M.z();
  }

  if (CE == Estatic::noE_V) {
    target.V_noE() += field;
  } else {
    target.V() += field;
  }
}

template void EwaldShapeCorrection::AddFieldAt<Estatic::V>(
    PolarSite&, EwaldChargeState) const;
template void EwaldShapeCorrection::AddFieldAt<Estatic::noE_V>(
    PolarSite&, EwaldChargeState) const;

}  // namespace xtp
}  // namespace votca
