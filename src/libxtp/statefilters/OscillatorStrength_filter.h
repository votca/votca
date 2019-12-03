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
#ifndef __VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H
#define __VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H

#include "statefilter_base.h"

namespace votca {
namespace xtp {

/**
    \brief OscillatorStrength_filter
    tracks states according to their f-value
    only works for singlets


 */

class OscillatorStrength_filter : public StateFilter_base {
 public:
  std::string Identify() const final { return "Oscillatorstrength"; }

  void Initialize(const tools::Property& options) final;

  void Info(Logger& log) const final;

  void UpdateHist(const Orbitals& orb) final;

  std::vector<Index> CalcIndeces(const Orbitals& orb) const final;

  void WriteToCpt(CheckpointWriter& w) final;

  void ReadFromCpt(CheckpointReader& r) final;

 private:
  double _threshold = 0.0;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_OSCILLATORSTRENGTH_FILTER_H */
