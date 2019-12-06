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

#include <memory>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/statefilter_base.h>

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
  void Initialize(const tools::Property& options);
  void setLogger(Logger* log) { _log = log; }
  void setInitialState(const QMState& state) { _statehist.push_back(state); }
  void PrintInfo() const;
  QMState InitialState() const { return _statehist[0]; }
  QMState CalcStateAndUpdate(const Orbitals& orbitals);
  QMState CalcState(const Orbitals& orbitals) const;

  void WriteToCpt(CheckpointWriter& w) const;

  void ReadFromCpt(CheckpointReader& r);

 private:
  void UpdateLastCoeff(const Orbitals& orbitals);
  std::vector<Index> CollapseResults(
      std::vector<std::vector<Index> >& results) const;
  std::vector<Index> ComparePairofVectors(std::vector<Index>& vec1,
                                          std::vector<Index>& vec2) const;

  Logger* _log;

  std::vector<QMState> _statehist;
  std::vector<std::unique_ptr<StateFilter_base> > _filters;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATETRACKER_H
