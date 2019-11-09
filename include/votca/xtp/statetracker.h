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
#ifndef VOTCA_XTP_STATETRACKER_H
#define VOTCA_XTP_STATETRACKER_H

#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/populationanalysis.h>
#include <votca/xtp/qmfragment.h>
#include <votca/xtp/qmstate.h>

namespace votca {
namespace xtp {
/**
 *  \brief  Tracks from a spectrum of states the state, which fullfills certain
 * criteria
 *
 *
 */

class StateTracker {

 public:
  void Initialize(tools::Property& options);
  void setLogger(Logger* log) { _log = log; }
  void setInitialState(const QMState& state) { _statehist.push_back(state); }
  void PrintInfo() const;
  QMState InitialState() const { return _statehist[0]; }
  QMState CalcStateAndUpdate(const Orbitals& orbitals);
  QMState CalcState(const Orbitals& orbitals) const;

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

 private:
  std::vector<Index> OscTracker(const Orbitals& orbitals) const;
  std::vector<Index> LocTracker(const Orbitals& orbitals) const;
  std::vector<Index> DeltaQTracker(const Orbitals& orbitals) const;
  std::vector<Index> OverlapTracker(const Orbitals& orbitals) const;

  Eigen::VectorXd CalculateOverlap(const Orbitals& orbitals) const;

  void UpdateLastCoeff(const Orbitals& orbitals);
  Eigen::MatrixXd CalcOrthoCoeffs(const Orbitals& orbitals) const;

  std::vector<Index> CollapseResults(
      std::vector<std::vector<Index> >& results) const;
  std::vector<Index> ComparePairofVectors(std::vector<Index>& vec1,
                                          std::vector<Index>& vec2) const;

  Logger* _log;

  std::vector<QMState> _statehist;

  bool _use_osctracker = false;
  double _oscthreshold = 0.0;

  bool _use_overlaptracker = false;
  Eigen::VectorXd _laststatecoeff;
  double _overlapthreshold = 0.0;

  bool _use_localizationtracker = false;
  QMFragment<double> _fragment_loc;

  bool _use_dQtracker = false;
  QMFragment<double> _fragment_dQ;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATETRACKER_H
