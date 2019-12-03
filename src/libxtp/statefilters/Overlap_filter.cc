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
  _threshold = options.ifExistsReturnElseThrowRuntimeError<double>("overlap");
}

void Overlap_filter::Info(Logger& log) const {
  if (_threshold == 0.0) {
    XTP_LOG(Log::error, *_log)
        << "Using overlap filer with no threshold " << flush;
  } else {
    XTP_LOG(Log::error, *_log)
        << "Using overlap filer with threshold " << _threshold << flush;
  }
}

Eigen::VectorXd Overlap_filter::CalculateOverlap(
    const Orbitals& orbitals) const {
  Eigen::MatrixXd coeffs = CalcOrthoCoeffs(orbitals);
  Eigen::VectorXd overlap = (coeffs * _laststatecoeff).cwiseAbs2();
  return overlap;
}

Eigen::MatrixXd Overlap_filter::CalcOrthoCoeffs(
    const Orbitals& orbitals) const {
  QMStateType type = _statehist[0].Type();
  Eigen::MatrixXd coeffs;
  if (type.isSingleParticleState()) {
    if (type == QMStateType::DQPstate) {
      coeffs = orbitals.CalculateQParticleAORepresentation();
    } else {
      coeffs = orbitals.MOs().eigenvectors();
    }
  } else {
    throw std::runtime_error(
        "Overlap for excitons is implemented within TDA. Use a different "
        "filter");
  }
  return coeffs;
}

void Overlap_filter::UpdateHist(const Orbitals&) {
  Eigen::MatrixXd ortho_coeffs = CalcOrthoCoeffs(orbitals);
  Index offset = 0;
  if (_statehist[0].Type() == QMStateType::DQPstate) {
    offset = orbitals.getGWAmin();
  }
  _laststatecoeff = ortho_coeffs.col(_statehist.back().StateIdx() - offset);
}

std::vector<Index> Overlap_filter::CalcIndeces(const Orbitals& orb) const {

  std::vector<Index> indexes;
  if (_statehist.size() <= 1) {
    indexes = std::vector<Index>{_statehist[0].StateIdx()};
    return indexes;
  }

  Eigen::VectorXd Overlap = CalculateOverlap(orbitals);
  Index validelements = Index(Overlap.size());
  for (Index i = 0; i < Index(Overlap.size()); i++) {
    if (Overlap(i) < _threshold) {
      validelements--;
    }
  }

  std::vector<Index> index = std::vector<Index>(Overlap.size());
  std::iota(index.begin(), index.end(), 0);
  std::stable_sort(index.begin(), index.end(), [&Overlap](Index i1, Index i2) {
    return Overlap[i1] > Overlap[i2];
  });

  Index offset = 0;
  if (_statehist[0].Type().isGWState()) {
    offset = orbitals.getGWAmin();
  }

  for (Index i : index) {
    if (Index(indexes.size()) == validelements) {
      break;
    }
    indexes.push_back(i + offset);
  }
  return indexes;
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