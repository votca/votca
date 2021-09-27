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
#ifndef VOTCA_XTP_STATETRACKER_H
#define VOTCA_XTP_STATETRACKER_H

// Standard includes
#include <memory>

// Local VOTCA includes
#include "logger.h"
#include "orbitals.h"
#include "qmstate.h"
#include "statefilter_base.h"

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
  void setLogger(Logger* log) { log_ = log; }
  void setInitialState(const QMState& state) { statehist_.push_back(state); }
  void PrintInfo() const;
  QMState InitialState() const { return statehist_[0]; }
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

  Logger* log_;

  std::vector<QMState> statehist_;
  std::vector<std::unique_ptr<StateFilter_base> > filters_;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATETRACKER_H
