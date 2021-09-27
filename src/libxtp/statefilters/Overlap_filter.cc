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

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"

// Local private VOTCA includes
#include "Overlap_filter.h"

namespace votca {
namespace xtp {

void Overlap_filter::Initialize(const tools::Property& options) {
  threshold_ = options.get(".").as<double>();
}

void Overlap_filter::Info(Logger& log) const {
  if (threshold_ == 0.0) {
    XTP_LOG(Log::error, log)
        << "Using overlap filter with no threshold " << std::flush;
  } else {
    XTP_LOG(Log::error, log)
        << "Using overlap filter with threshold " << threshold_ << std::flush;
  }
}

Eigen::VectorXd Overlap_filter::CalculateOverlap(const Orbitals& orb,
                                                 QMStateType type) const {
  AOOverlap S_ao;
  S_ao.Fill(orb.SetupDftBasis());

  Eigen::MatrixXd coeffs = CalcAOCoeffs(orb, type);
  if (type.isSingleParticleState()) {
    return (coeffs.transpose() * S_ao.Matrix() * laststatecoeff_).cwiseAbs();
  } else {

    Index basis = orb.getBasisSetSize();
    Eigen::VectorXd overlap = Eigen::VectorXd::Zero(coeffs.cols());
#pragma omp parallel for schedule(dynamic)
    for (Index i = 0; i < coeffs.cols(); i++) {
      {
        Eigen::VectorXd state = coeffs.col(i).head(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> mat(state.data(), basis, basis);
        Eigen::VectorXd laststate = laststatecoeff_.head(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> lastmat(laststate.data(), basis,
                                                  basis);

        overlap(i) = (mat * S_ao.Matrix() * lastmat.transpose())
                         .cwiseProduct(S_ao.Matrix())
                         .sum();
      }
      if (!orb.getTDAApprox()) {
        Eigen::VectorXd state = coeffs.col(i).tail(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> mat(state.data(), basis, basis);
        Eigen::VectorXd laststate = laststatecoeff_.tail(basis * basis);
        Eigen::Map<const Eigen::MatrixXd> lastmat(laststate.data(), basis,
                                                  basis);

        overlap(i) -= (mat * S_ao.Matrix() * lastmat.transpose())
                          .cwiseProduct(S_ao.Matrix())
                          .sum();
      }
    }
    return overlap.cwiseAbs();
  }
}

Eigen::MatrixXd Overlap_filter::CalcExcitonAORepresentation(
    const Orbitals& orb, QMStateType type) const {
  Eigen::MatrixXd coeffs;
  Index nostates = orb.NumberofStates(type);
  Index bse_cmax = orb.getBSEcmax();
  Index bse_cmin = orb.getBSEcmin();
  Index bse_vmax = orb.getBSEvmax();
  Index bse_vmin = orb.getBSEvmin();
  Index bse_vtotal = bse_vmax - bse_vmin + 1;
  Index bse_ctotal = bse_cmax - bse_cmin + 1;
  Index basis = orb.getBasisSetSize();
  Index bse_size_ao = basis * basis;
  auto occlevels = orb.MOs().eigenvectors().middleCols(bse_vmin, bse_vtotal);
  auto virtlevels = orb.MOs().eigenvectors().middleCols(bse_cmin, bse_ctotal);

  if (orb.getTDAApprox()) {
    coeffs.resize(bse_size_ao, nostates);
  } else {
    coeffs.resize(2 * bse_size_ao, nostates);
  }
#pragma omp parallel for schedule(dynamic)
  for (Index i = 0; i < nostates; i++) {
    {
      Eigen::VectorXd exciton;
      if (type == QMStateType::Singlet) {
        exciton = orb.BSESinglets().eigenvectors().col(i);
      } else {
        exciton = orb.BSETriplets().eigenvectors().col(i);
      }
      Eigen::Map<const Eigen::MatrixXd> mat(exciton.data(), bse_ctotal,
                                            bse_vtotal);
      const Eigen::MatrixXd aomatrix =
          occlevels * mat.transpose() * virtlevels.transpose();
      coeffs.col(i).head(bse_size_ao) =
          Eigen::Map<const Eigen::VectorXd>(aomatrix.data(), bse_size_ao);
    }
    if (!orb.getTDAApprox()) {
      Eigen::VectorXd exciton;
      if (type == QMStateType::Singlet) {
        exciton = orb.BSESinglets().eigenvectors2().col(i);
      } else {
        exciton = orb.BSETriplets().eigenvectors2().col(i);
      }
      Eigen::Map<const Eigen::MatrixXd> mat(exciton.data(), bse_ctotal,
                                            bse_vtotal);
      const Eigen::MatrixXd aomatrix =
          occlevels * mat.transpose() * virtlevels.transpose();
      coeffs.col(i).tail(bse_size_ao) =
          Eigen::Map<const Eigen::VectorXd>(aomatrix.data(), bse_size_ao);
    }
  }
  return coeffs;
}

Eigen::MatrixXd Overlap_filter::CalcAOCoeffs(const Orbitals& orb,
                                             QMStateType type) const {
  Eigen::MatrixXd coeffs;
  if (type.isSingleParticleState()) {
    if (type == QMStateType::DQPstate) {
      coeffs = orb.CalculateQParticleAORepresentation();
    } else {
      coeffs = orb.MOs().eigenvectors();
    }
  } else {
    coeffs = CalcExcitonAORepresentation(orb, type);
  }
  return coeffs;
}

void Overlap_filter::UpdateHist(const Orbitals& orb, QMState state) {
  Eigen::MatrixXd aocoeffs = CalcAOCoeffs(orb, state.Type());
  Index offset = 0;
  if (state.Type() == QMStateType::DQPstate) {
    offset = orb.getGWAmin();
  }
  laststatecoeff_ = aocoeffs.col(state.StateIdx() - offset);
}

std::vector<Index> Overlap_filter::CalcIndeces(const Orbitals& orb,
                                               QMStateType type) const {
  Index offset = 0;
  if (type.isGWState()) {
    offset = orb.getGWAmin();
  }
  Eigen::VectorXd Overlap = CalculateOverlap(orb, type);
  return ReduceAndSortIndecesUp(Overlap, offset, threshold_);
}

void Overlap_filter::WriteToCpt(CheckpointWriter& w) {
  w(laststatecoeff_, "laststatecoeff");
  w(threshold_, "threshold");
}

void Overlap_filter::ReadFromCpt(CheckpointReader& r) {
  r(laststatecoeff_, "laststatecoeff");
  r(threshold_, "threshold");
}

}  // namespace xtp
}  // namespace votca
