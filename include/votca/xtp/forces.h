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
#ifndef VOTCA_XTP_FORCES_H
#define VOTCA_XTP_FORCES_H

#include <stdio.h>
#include <votca/xtp/gwbseengine.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmatom.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/statefilter.h>

namespace votca {
namespace xtp {

class Forces {
 public:
  Forces(GWBSEEngine& gwbse_engine, const Statefilter& filter)
      : _gwbse_engine(gwbse_engine),
        _filter(filter),
        _remove_total_force(false){};

  void Initialize(tools::Property& options);
  void Calculate(const Orbitals& orbitals);

  void setLog(Logger* pLog) { _pLog = pLog; }

  const Eigen::MatrixX3d& GetForces() const { return _forces; };
  void Report() const;

 private:
  Eigen::Vector3d NumForceForward(Orbitals orbitals, int atom_index);
  Eigen::Vector3d NumForceCentral(Orbitals orbitals, int atom_index);
  void RemoveTotalForce();

  double _displacement;
  std::string _force_method;

  GWBSEEngine& _gwbse_engine;
  const Statefilter& _filter;
  bool _remove_total_force;

  Eigen::MatrixX3d _forces;
  Logger* _pLog;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_FORCES_H
