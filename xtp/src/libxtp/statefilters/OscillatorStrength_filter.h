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
#ifndef VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H
#define VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H

// Local VOTCA includes
#include "votca/xtp/statefilter_base.h"

namespace votca {
namespace xtp {

/**
    \brief OscillatorStrength_filter
    tracks states according to their f-value
    only works for singlets


 */

class OscillatorStrength_filter : public StateFilter_base {
 public:
  std::string Identify() const final { return "oscillatorstrength"; }

  void Initialize(const tools::Property& options) final;

  void Info(Logger& log) const final;

  bool NeedsInitialState() const final { return false; }

  void UpdateHist(const Orbitals& orb, QMState state) final;

  std::vector<Index> CalcIndeces(const Orbitals& orb,
                                 QMStateType type) const final;

  void WriteToCpt(CheckpointWriter& w) final;

  void ReadFromCpt(CheckpointReader& r) final;

 private:
  double threshold_ = 0.0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H
