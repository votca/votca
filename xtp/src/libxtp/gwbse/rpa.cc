/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#include "votca/xtp/rpa.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/screening_kernels.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace votca {
namespace xtp {

void RPA::UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                                 const Eigen::VectorXd& gwaenergies,
                                 Index qpmin) {
  Index rpatotal = rpamax_ - rpamin_ + 1;
  energies_ = dftenergies.segment(rpamin_, rpatotal);
  Index gwsize = Index(gwaenergies.size());

  energies_.segment(qpmin - rpamin_, gwsize) = gwaenergies;

  ShiftUncorrectedEnergies(dftenergies, qpmin, gwsize);
}

// Shifts the energies of levels in the RPA but outside the QP window, see
// OutOfWindowShift.
void RPA::ShiftUncorrectedEnergies(const Eigen::VectorXd& dftenergies,
                                   Index qpmin, Index gwsize) {
  out_of_window_.Apply(energies_, dftenergies, rpamin_, rpamax_, homo_, qpmin,
                       qpmin + gwsize - 1);
}

void RPA::VisitHoleVirtualRows(Index m_level,
                               const WeightedGram::RowsVisitor& use) const {
  const Index n_unocc = rpamax_ - homo_;
  use(Mmn_[m_level].bottomRows(n_unocc));
}

Eigen::MatrixXd RPA::ResponseSum(const WeightedGram::WeightsFn& weights) const {
  const Index size = Mmn_.auxsize();
  const Index n_occ = homo_ + 1 - rpamin_;
  // holes left out of the response (TCMatrix_gwbse::ExcludedHoles)
  const Index skip = std::min(Mmn_.ExcludedHoles(), n_occ);
  OptionalTiming timing(timings_, "screening: epsilon (RPA)");
  return WeightedGram::Compute(
      n_occ - skip, size,
      [&weights, skip](Index m, Eigen::VectorXd& w) { weights(m + skip, w); },
      [this, skip](Index m, const WeightedGram::RowsVisitor& use) {
        VisitHoleVirtualRows(m + skip, use);
      });
}

template <bool imag>
Eigen::MatrixXd RPA::calculate_epsilon(double frequency) const {
  const Index n_unocc = rpamax_ - homo_;
  const double freq2 = frequency * frequency;
  const double eta2 = eta_ * eta_;
  Eigen::MatrixXd result =
      ResponseSum([&](Index m_level, Eigen::VectorXd& denom) {
        const Eigen::ArrayXd deltaE =
            energies_.tail(n_unocc).array() - energies_(m_level);
        if (imag) {
          denom = 4 * deltaE / (deltaE.square() + freq2);
        } else {
          Eigen::ArrayXd deltEf = deltaE - frequency;
          Eigen::ArrayXd sum = deltEf / (deltEf.square() + eta2);
          deltEf = deltaE + frequency;
          sum += deltEf / (deltEf.square() + eta2);
          denom = 2 * sum;
        }
      });
  result.diagonal().array() += 1.0;
  if (Mmn_.HasAdditionalResponse()) {
    result += Mmn_.AdditionalResponse(imag, frequency);
  }
  if (Mmn_.AdditionalPoleCount() > 0) {
    // the poles as transitions, with the weights above
    const Eigen::ArrayXd& W = Mmn_.AdditionalPoleEnergies().array();
    Eigen::ArrayXd w;
    if (imag) {
      w = 4 * W / (W.square() + freq2);
    } else {
      w = 2 * ((W - frequency) / ((W - frequency).square() + eta2) +
               (W + frequency) / ((W + frequency).square() + eta2));
    }
    const Eigen::MatrixXd C = Mmn_.AdditionalPoleVectors();
    result.noalias() += C * w.matrix().asDiagonal() * C.transpose();
  }
  return result;
}

template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency) const;
template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency) const;

Eigen::MatrixXd RPA::calculate_epsilon_r(std::complex<double> frequency) const {
  if (Mmn_.HasAdditionalResponse()) {
    throw std::runtime_error(
        "RPA: an additional response (inactive_screening with the "
        "plasmon-pole model) is only available at the plasmon-pole fitting "
        "points, not for complex frequencies (sigma_integrator cda).");
  }
  const Index n_unocc = rpamax_ - homo_;
  const double sigma_1 = std::pow(frequency.imag() + eta_, 2);
  const double sigma_2 = std::pow(frequency.imag() - eta_, 2);
  Eigen::MatrixXd result =
      ResponseSum([&](Index m_level, Eigen::VectorXd& chi) {
        const Eigen::ArrayXd deltaE =
            energies_.tail(n_unocc).array() - energies_(m_level);
        const Eigen::ArrayXd deltaEm = frequency.real() - deltaE;
        const Eigen::ArrayXd deltaEp = frequency.real() + deltaE;
        // the factor -2 of the response is folded into the weights
        chi = -2.0 * (deltaEm * (deltaEm.abs2() + sigma_1).inverse() -
                      deltaEp * (deltaEp.abs2() + sigma_2).inverse());
      });
  if (Mmn_.AdditionalPoleCount() > 0) {
    const Eigen::ArrayXd W = Mmn_.AdditionalPoleEnergies().array();
    const Eigen::ArrayXd deltaEm = frequency.real() - W;
    const Eigen::ArrayXd deltaEp = frequency.real() + W;
    const Eigen::VectorXd chi =
        -2.0 * (deltaEm * (deltaEm.square() + sigma_1).inverse() -
                deltaEp * (deltaEp.square() + sigma_2).inverse());
    const Eigen::MatrixXd C = Mmn_.AdditionalPoleVectors();
    result.noalias() += C * chi.asDiagonal() * C.transpose();
  }
  result.diagonal().array() += 1.0;
  return result;
}

RPA::rpa_eigensolution RPA::Diagonalize_H2p() const {
  // transitions of M, then the additional poles (TCMatrix_gwbse)
  Eigen::VectorXd AmB = Calculate_H2p_AmB();
  Eigen::MatrixXd ApB = Calculate_H2p_ApB();
  const Index rpasize = AmB.size();

  RPA::rpa_eigensolution sol;
  sol.ERPA_correlation = -0.25 * (ApB.trace() + AmB.sum());

  // C = AmB^1/2 * ApB * AmB^1/2
  Eigen::MatrixXd& C = ApB;
  C.applyOnTheLeft(AmB.cwiseSqrt().asDiagonal());
  C.applyOnTheRight(AmB.cwiseSqrt().asDiagonal());

  const SymmetricEigenSystem es = Diagonalize_H2p_C(C);

  // Do not remove this line! It has to be there for MKL to not crash
  sol.omega = Eigen::VectorXd::Zero(es.values.size());
  sol.omega = es.values.cwiseSqrt();
  sol.ERPA_correlation += 0.5 * sol.omega.sum();

  {
    std::ostringstream oss;
    oss << TimeStamp() << " RKS H2p poles (first 20, eV): ";
    const Index nprint = std::min<Index>(20, sol.omega.size());
    for (Index i = 0; i < nprint; ++i) {
      oss << std::fixed << std::setprecision(6)
          << tools::conv::hrt2ev * sol.omega(i);
      if (i + 1 < nprint) {
        oss << ", ";
      }
    }
    XTP_LOG(Log::error, log_) << oss.str() << std::flush;
  }

  XTP_LOG(Log::info, log_) << TimeStamp()
                           << " Lowest neutral excitation energy (eV): "
                           << tools::conv::hrt2ev * sol.omega.minCoeff()
                           << std::flush;

  // RPA correlation energy calculated from Eq.9 of J. Chem. Phys. 132, 234114
  // (2010)
  XTP_LOG(Log::error, log_)
      << TimeStamp()
      << " RPA correlation energy (Hartree): " << sol.ERPA_correlation
      << std::flush;

  sol.XpY = Eigen::MatrixXd(rpasize, rpasize);

  Eigen::VectorXd AmB_sqrt = AmB.cwiseSqrt();
  Eigen::VectorXd Omega_sqrt_inv = sol.omega.cwiseSqrt().cwiseInverse();
  for (int s = 0; s < rpasize; s++) {
    sol.XpY.col(s) =
        Omega_sqrt_inv(s) * AmB_sqrt.cwiseProduct(es.vectors.col(s));
  }

  return sol;
}

Eigen::VectorXd RPA::Calculate_H2p_AmB() const {
  const Index lumo = homo_ + 1;
  const Index n_occ = lumo - rpamin_;
  const Index n_unocc = rpamax_ - lumo + 1;
  // excluded holes (TCMatrix_gwbse::ExcludedHoles) have no transitions
  const Index skip = std::min(Mmn_.ExcludedHoles(), n_occ);
  const Index rpasize = (n_occ - skip) * n_unocc;
  const Index npoles = Mmn_.AdditionalPoleCount();
  vc2index vc = vc2index(skip, 0, n_unocc);
  Eigen::VectorXd AmB = Eigen::VectorXd::Zero(rpasize + npoles);
  for (Index v = skip; v < n_occ; v++) {
    Index i = vc.I(v, 0);
    AmB.segment(i, n_unocc) =
        energies_.segment(n_occ, n_unocc).array() - energies_(v);
  }
  AmB.tail(npoles) = Mmn_.AdditionalPoleEnergies();
  return AmB;
}

Eigen::MatrixXd RPA::Calculate_H2p_ApB() const {
  const Index lumo = homo_ + 1;
  const Index n_occ = lumo - rpamin_;
  const Index n_unocc = rpamax_ - lumo + 1;
  const Index skip = std::min(Mmn_.ExcludedHoles(), n_occ);
  const Index rpasize = (n_occ - skip) * n_unocc;
  const Index npoles = Mmn_.AdditionalPoleCount();
  vc2index vc = vc2index(skip, 0, n_unocc);
  Eigen::MatrixXd ApB =
      Eigen::MatrixXd::Zero(rpasize + npoles, rpasize + npoles);
#pragma omp parallel for schedule(guided)
  for (Index v2 = skip; v2 < n_occ; v2++) {
    Index i2 = vc.I(v2, 0);
    const Eigen::MatrixXd Mmn_v2T =
        Mmn_[v2].middleRows(n_occ, n_unocc).transpose();
    for (Index v1 = v2; v1 < n_occ; v1++) {
      Index i1 = vc.I(v1, 0);
      // Multiply with factor 2 to sum over both (identical) spin states
      ApB.block(i1, i2, n_unocc, n_unocc) =
          2 * 2 * Mmn_[v1].middleRows(n_occ, n_unocc) * Mmn_v2T;
    }
  }
  if (npoles > 0) {
    // the poles couple like transitions with integrals c_k (same factor 4)
    const Eigen::MatrixXd C = Mmn_.AdditionalPoleVectors();
    ApB.bottomRightCorner(npoles, npoles) = 4 * C.transpose() * C;
#pragma omp parallel for schedule(dynamic)
    for (Index v = skip; v < n_occ; v++) {
      ApB.block(rpasize, vc.I(v, 0), npoles, n_unocc) =
          4 * C.transpose() * Mmn_[v].middleRows(n_occ, n_unocc).transpose();
    }
  }
  ApB.diagonal() += Calculate_H2p_AmB();
  return ApB;
}

SymmetricEigenSystem RPA::Diagonalize_H2p_C(const Eigen::MatrixXd& C) const {
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Diagonalizing two-particle Hamiltonian "
      << std::flush;
  const SymmetricEigenSystem es = SymmetricEigen(C);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Diagonalization done " << std::flush;
  double minCoeff = es.values.minCoeff();
  if (minCoeff <= 0.0) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Detected non-positive eigenvalue: " << minCoeff
        << std::flush;
    throw std::runtime_error("Detected non-positive eigenvalue.");
  }
  return es;
}

}  // namespace xtp
}  // namespace votca
