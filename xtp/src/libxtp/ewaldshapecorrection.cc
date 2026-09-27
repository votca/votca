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

// Standard includes
#include <cmath>

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

void EwaldShapeCorrection::Accumulate(Moments& m, const PolarSite& site,
                                      const Eigen::Vector3d& pos) {
  // Rank 1: charge and dipole. See CalcStaticEnergyBetween for why the
  // sites' intrinsic quadrupoles are not folded in here, and for what
  // that costs relative to legacy.
  const double q = site.getCharge();
  const Eigen::Vector3d mu = site.getStaticDipole();
  m.q0 += q;
  m.q1 += q * pos;
  m.q1 += mu;
  m.q2 += 0.5 * q * pos * pos.transpose();
  m.q2 += mu * pos.transpose();
}

EwaldShapeCorrection::Moments EwaldShapeCorrection::BackgroundMoments(
    EwaldChargeState source_state,
    const std::vector<const PolarSite*>& exclusions) const {
  Moments m;
  for (Index id : registry_.AllIds()) {
    if (!registry_.Has(id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(id, source_state);
    for (const PolarSite& site : segment) {
      // Identity by address, exactly as the reciprocal-space cross
      // energy does it, so a site held out there is held out here too.
      bool excluded = false;
      for (const PolarSite* skip : exclusions) {
        if (skip == &site) {
          excluded = true;
          break;
        }
      }
      if (excluded) {
        continue;
      }
      Accumulate(m, site, site.getPos());
    }
  }
  return m;
}

double EwaldShapeCorrection::CalcStaticEnergyBetween(
    const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>>& foreground,
    const std::vector<const PolarSite*>& background_exclusions,
    EwaldChargeState source_state) const {
  // See this method's own declaration for the formula, for why all three
  // moments are needed, and for the permanent-only convention.

  // The foreground's moments come from the sites themselves -- so in the
  // job's OWN charge state -- but at the positions supplied, which are
  // the positions those sites actually occupy.
  Moments fg;
  for (const auto& entry : foreground) {
    Accumulate(fg, *entry.first, entry.second);
  }

  const Moments bg = BackgroundMoments(source_state, background_exclusions);

  if (shape_ == EwaldShape::Cube) {
    const double bracket =
        fg.q0 * bg.q2.trace() + bg.q0 * fg.q2.trace() - fg.q1.dot(bg.q1);
    return -(4.0 * kPi / (3.0 * volume_)) * bracket;
  }
  const double bracket =
      fg.q0 * bg.q2(2, 2) + bg.q0 * fg.q2(2, 2) - fg.q1.z() * bg.q1.z();
  return -(4.0 * kPi / volume_) * bracket;
}

EwaldShapeCorrection::Moments EwaldShapeCorrection::BackgroundInducedMoments(
    EwaldChargeState source_state,
    const std::vector<const PolarSite*>& exclusions) const {
  Moments m;
  for (Index id : registry_.AllIds()) {
    if (!registry_.Has(id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(id, source_state);
    for (const PolarSite& site : segment) {
      bool excluded = false;
      for (const PolarSite* skip : exclusions) {
        if (skip == &site) {
          excluded = true;
          break;
        }
      }
      if (excluded) {
        continue;
      }
      // Induced dipole only: no charge, hence no q0 and no 0.5*q*r*r^T.
      const Eigen::Vector3d mu = site.getInducedDipole();
      m.q1 += mu;
      m.q2 += mu * site.getPos().transpose();
    }
  }
  return m;
}

double EwaldShapeCorrection::CalcInducedSourceEnergyBetween(
    const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>>& foreground,
    const std::vector<const PolarSite*>& background_exclusions,
    EwaldChargeState source_state) const {
  // See this method's own declaration.
  Moments fg;
  for (const auto& entry : foreground) {
    Accumulate(fg, *entry.first, entry.second);
  }

  const Moments bg =
      BackgroundInducedMoments(source_state, background_exclusions);

  // bg.q0 is zero by construction, so its term is dropped rather than
  // written out and multiplied by zero.
  if (shape_ == EwaldShape::Cube) {
    const double bracket = fg.q0 * bg.q2.trace() - fg.q1.dot(bg.q1);
    return -(4.0 * kPi / (3.0 * volume_)) * bracket;
  }
  const double bracket = fg.q0 * bg.q2(2, 2) - fg.q1.z() * bg.q1.z();
  return -(4.0 * kPi / volume_) * bracket;
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
