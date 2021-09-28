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
#ifndef VOTCA_XTP_GEOMETRY_OPTIMIZATION_H
#define VOTCA_XTP_GEOMETRY_OPTIMIZATION_H

// Standard includes
#include <cstdio>

// Local VOTCA includes
#include "bfgs_trm.h"
#include "energy_costfunction.h"
#include "gwbseengine.h"
#include "logger.h"
#include "qmatom.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

class GeometryOptimization {
 public:
  GeometryOptimization(GWBSEEngine& gwbse_engine, Orbitals& orbitals)
      : gwbse_engine_(gwbse_engine),
        orbitals_(orbitals){

        };

  void Initialize(tools::Property& options);

  void setLog(Logger* pLog) { pLog_ = pLog; }

  void Evaluate();

 private:
  static void Report(const BFGSTRM& bfgstrm, const Forces& forces,
                     Logger& pLog);
  static void WriteTrajectory(const std::string& filename, QMMolecule& atoms,
                              const BFGSTRM& bfgstrm);

  QMState opt_state_;
  std::string optimizer_;
  std::string trajfile_;
  GWBSEEngine& gwbse_engine_;
  Orbitals& orbitals_;

  Energy_costfunction::conv_paras conv_;
  Index max_iteration_;
  double trust_radius_;

  tools::Property statetracker_options_;
  tools::Property force_options_;

  Logger* pLog_;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_GEOMETRY_OPTIMIZATION_H
