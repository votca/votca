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
#ifndef VOTCA_XTP_FORCES_H
#define VOTCA_XTP_FORCES_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "gwbseengine.h"
#include "logger.h"
#include "qmatom.h"
#include "segment.h"

namespace votca {
namespace xtp {

class StateTracker;

class Forces {
 public:
  Forces(GWBSEEngine& gwbse_engine, const StateTracker& tracker)
      : gwbse_engine_(gwbse_engine), tracker_(tracker){};

  void Initialize(tools::Property& options);
  void Calculate(const Orbitals& orbitals);

  void setLog(Logger* pLog) { pLog_ = pLog; }

  const Eigen::MatrixX3d& GetForces() const { return forces_; };
  void Report() const;

 private:
  Eigen::Vector3d NumForceForward(Orbitals orbitals, Index atom_index);
  Eigen::Vector3d NumForceCentral(Orbitals orbitals, Index atom_index);
  void RemoveTotalForce();

  double displacement_;
  std::string force_method_;

  GWBSEEngine& gwbse_engine_;
  const StateTracker& tracker_;
  bool remove_total_force_ = true;

  Eigen::MatrixX3d forces_;
  Logger* pLog_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_FORCES_H
