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

#include "Density_filter.h"
namespace votca {
namespace xtp {

void Density_filter::Initialize(const tools::Property& options) {
  _threshold = options.ifExistsReturnElseThrowRuntimeError<double>(".");
}

void Density_filter::Info(Logger& log) const {
  if (_threshold == 0.0) {
    XTP_LOG(Log::error, log)
        << "Using density filter with no threshold " << std::flush;
  } else {
    XTP_LOG(Log::error, log)
        << "Using density filter with threshold " << _threshold << std::flush;
  }
}

void Density_filter::UpdateHist(const Orbitals& orb, QMState state) {
  _laststate_dmat = orb.DensityMatrixFull(state);
}

Eigen::VectorXd Density_filter::CalculateDNorm(const Orbitals& orb,
                                               QMStateType type) const {
  Index nostates = orb.NumberofStates(type);
  Eigen::VectorXd norm = Eigen::VectorXd::Zero(nostates);
  for (Index i = 0; i < nostates; i++) {
    QMState state(type, i, false);
    norm(i) = (orb.DensityMatrixFull(state) - _laststate_dmat).norm();
  }
  double lastnorm = _laststate_dmat.norm();
  return norm / lastnorm;
}

std::vector<Index> Density_filter::CalcIndeces(const Orbitals& orb,
                                               QMStateType type) const {
  Eigen::VectorXd Overlap = CalculateDNorm(orb, type);
  Index offset = 0;
  if (type == QMStateType::DQPstate) {
    offset = orb.getGWAmin();
  }
  return ReduceAndSortIndecesDown(Overlap, offset, _threshold);
}

void Density_filter::WriteToCpt(CheckpointWriter& w) {
  w(_laststate_dmat, "laststatedmat");
  w(_threshold, "threshold");
}

void Density_filter::ReadFromCpt(CheckpointReader& r) {
  r(_laststate_dmat, "laststatedmat");
  r(_threshold, "threshold");
}

}  // namespace xtp
}  // namespace votca