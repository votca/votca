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

#include "OscillatorStrength_filter.h"
namespace votca {
namespace xtp {

void OscillatorStrength_filter::Initialize(const tools::Property& options) {
  _threshold = options.ifExistsReturnElseThrowRuntimeError<double>(
      "oscillator_strength");
}

void OscillatorStrength_filter::Info(Logger& log) const {
  XTP_LOG(Log::error, log)
      << "Using oscillator strength tracker with threshold " << _threshold
      << std::flush;
}

void OscillatorStrength_filter::UpdateHist(const Orbitals&, QMState) { return; }

std::vector<Index> OscillatorStrength_filter::CalcIndeces(const Orbitals& orb,
                                                          QMStateType) const {
  Eigen::VectorXd oscs = orb.Oscillatorstrengths();
  std::vector<Index> indexes;
  for (Index i = 0; i < oscs.size(); i++) {
    if (oscs[i] > _threshold) {
      indexes.push_back(i);
    }
  }
  return indexes;
}

void OscillatorStrength_filter::WriteToCpt(CheckpointWriter& w) {

  w(_threshold, "threshold");
}

void OscillatorStrength_filter::ReadFromCpt(CheckpointReader& r) {
  r(_threshold, "threshold");
}

}  // namespace xtp
}  // namespace votca