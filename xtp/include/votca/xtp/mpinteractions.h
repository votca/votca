
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

#pragma once
#ifndef VOTCA_XTP_MPINTERACTIONS_H
#define VOTCA_XTP_MPINTERACTIONS_H

#include <algorithm>
#include <vector>
// Local VOTCA includes
#include "votca/xtp/bgsite.h"
#include "votca/xtp/interactiontensor.h"
#include "votca/xtp/multipole.h"

namespace votca {
namespace xtp {

template <Screening screen, Index maxMPRank>
class MPField {
 public:
  Eigen::Vector3d fieldAtBy(BGSite site1, BGSite site2) {
    Eigen::Vector3d dr = site1.getPos() - site2.getPos();
    return field(site2.getMP(), dr);
  };

 private:
  InteractionTensor<screen, maxMPRank + 1> T_;

  Eigen::Vector3d field(const Multipole& mp, const Eigen::Vector3d& dr);
};

template <Screening screen, Index maxMPRank>
class MPEnergy {
 public:
  double energy(BGSite site1, BGSite site2) {
    Eigen::Vector3d dr = site1.getPos() - site2.getPos();
    return energy(site1.getMP(), site2.getMP(), dr);
  };

 private:
  InteractionTensor<screen, 2 * maxMPRank> T_;

  double energy(const Multipole& mp1, const Multipole& mp2,
                const Eigen::Vector3d& dr);
};

template <Screening screen, Index maxMPRank>
Eigen::Vector3d MPField<screen, maxMPRank>::field(const Multipole& mp,
                                                  const Eigen::Vector3d& dr) {
  T_.computeTensors(dr);
  Eigen::Vector3d res = Eigen::Vector3d::Zero();
  res -= T_.interact(mp.charge);
  if constexpr (maxMPRank > 0) {
    res += T_.interact(mp.dipole);
  }
  if constexpr (maxMPRank > 1) {
    res -= (1.0 / 3.0) * T_.interact(mp.quadrupole);
  }
  return -res;
}

template <Screening screen, Index maxMPRank>
double MPEnergy<screen, maxMPRank>::energy(const Multipole& mp1,
                                           const Multipole& mp2,
                                           const Eigen::Vector3d& dr) {
  T_.computeTensors(dr);
  double res = 0.0;
  res += T_.interact(mp1.charge, mp2.charge);
  if constexpr (maxMPRank > 0) {
    res += -T_.interact(mp1.charge, mp2.dipole) +
           T_.interact(mp2.charge, mp1.dipole) -
           T_.interact(mp1.dipole, mp2.dipole);
  }
  if constexpr (maxMPRank > 1) {
    res += (1.0 / 3.0) * T_.interact(mp1.charge, mp2.quadrupole) +
           (1.0 / 3.0) * T_.interact(mp2.charge, mp1.quadrupole) + 
           (1.0 / 3.0) * T_.interact(mp1.dipole, mp2.quadrupole) -
           (1.0 / 3.0) * T_.interact(mp2.dipole, mp1.quadrupole) +
           (1.0 / 9.0) * T_.interact(mp1.quadrupole, mp2.quadrupole);
  }
  return res;
}

}  // namespace xtp
}  // namespace votca

#endif