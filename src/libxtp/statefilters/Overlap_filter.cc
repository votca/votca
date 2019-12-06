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

#include "Overlap_filter.h"
namespace votca {
namespace xtp {

void Overlap_filter::Initialize(const tools::Property& options) {
  _threshold = options.ifExistsReturnElseThrowRuntimeError<double>(".");
}

void Overlap_filter::Info(Logger& log) const {
  if (_threshold == 0.0) {
    XTP_LOG(Log::error, log)
        << "Using overlap filter with no threshold " << std::flush;
  } else {
    XTP_LOG(Log::error, log)
        << "Using overlap filter with threshold " << _threshold << std::flush;
  }
}

Eigen::VectorXd Overlap_filter::CalculateOverlap(const Orbitals& orb,
                                                 QMStateType type) const {
  Eigen::MatrixXd coeffs = CalcOrthoCoeffs(orb, type);
  Eigen::VectorXd overlap = (coeffs * _laststatecoeff).cwiseAbs2();
  return overlap;
}

Eigen::MatrixXd Overlap_filter::CalcOrthoCoeffs(const Orbitals& orb,
                                                QMStateType type) const {
  Eigen::MatrixXd coeffs;
  if (type.isSingleParticleState()) {
    if (type == QMStateType::DQPstate) {
      coeffs = orb.CalculateQParticleAORepresentation();
    } else {
      coeffs = orb.MOs().eigenvectors();
    }
  } else {
    throw std::runtime_error(
        "Overlap for excitons is implemented within TDA. Use a different "
        "filter");
  }
  return coeffs;
}

void Overlap_filter::UpdateHist(const Orbitals& orb, QMState state) {
  Eigen::MatrixXd ortho_coeffs = CalcOrthoCoeffs(orb, state.Type());
  Index offset = 0;
  if (state.Type() == QMStateType::DQPstate) {
    offset = orb.getGWAmin();
  }
  _laststatecoeff = ortho_coeffs.col(state.StateIdx() - offset);
}

std::vector<Index> Overlap_filter::CalcIndeces(const Orbitals& orb,
                                               QMStateType type) const {
  Index offset = 0;
  if (type.isGWState()) {
    offset = orb.getGWAmin();
  }
  Eigen::VectorXd Overlap = CalculateOverlap(orb, type);
  return ReduceAndSortIndecesUp(Overlap, offset, _threshold);
}

void Overlap_filter::WriteToCpt(CheckpointWriter& w) {
  w(_laststatecoeff, "laststatecoeff");
  w(_threshold, "threshold");
}

void Overlap_filter::ReadFromCpt(CheckpointReader& r) {
  r(_laststatecoeff, "laststatecoeff");
  r(_threshold, "threshold");
}

}  // namespace xtp
}  // namespace votca