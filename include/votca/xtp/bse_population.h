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

#pragma once
#ifndef VOTCA_XTP_BSE_POPULATIONS_H
#define VOTCA_XTP_BSE_POPULATIONS_H

#include <votca/xtp/checkpoint.h>
#include <votca/xtp/eigen.h>

/**
 * \brief Small container to keep occupation of BSE states for each atom
 *
 *
 */

namespace votca {
namespace xtp {
struct BSE_Population {
  Eigen::VectorXd H;
  Eigen::VectorXd E;
  double Gs = 0;

  void Initialize(long size) {
    H = Eigen::VectorXd::Zero(size);
    E = Eigen::VectorXd::Zero(size);
    Gs = 0;
  }

  void WriteToCpt(CheckpointWriter& w) const {
    w(H, "holeinfo");
    w(E, "electroninfo");
    w(Gs, "groundstate");
  }

  void ReadFromCpt(CheckpointReader& r) {
    r(H, "holeinfo");
    r(E, "electroninfo");
    r(Gs, "groundstate");
  }

  friend std::ostream& operator<<(std::ostream& out,
                                  const BSE_Population& pop) {
    if (pop.H.size() < 1) {
      return out;
    }
    Eigen::VectorXd diff = pop.H - pop.E;
    out << "GroundstateCharge:" << pop.Gs << "\n";
    out << "Index hole electron dQ Qeff\n";
    for (long i = 0; i < pop.H.size(); ++i) {
      out << i << " " << pop.H(i) << " " << pop.E(i) << " " << diff(i) << " "
          << diff(i) + pop.Gs << "\n";
    }
    return out;
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_POPULATIONS_H
