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
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/screening_kernels.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"
#include <iomanip>
#include <sstream>

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
  const Index n_unocc = rpamax_ - homo_;
  OptionalTiming timing(timings_, "screening: epsilon (RPA)");
  if (OpenMP_CUDA::UsingGPUs() == 0) {
    return WeightedGram::Compute(
        n_occ, size, weights,
        [this](Index m, const WeightedGram::RowsVisitor& use) {
          VisitHoleVirtualRows(m, use);
        });
  }
  // GPU path: per-thread products M^T diag(w) M on the devices
  OpenMP_CUDA transform;
  transform.createTemporaries(n_unocc, size);
#pragma omp parallel
  {
    Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
    for (Index m_level = 0; m_level < n_occ; m_level++) {
      Eigen::VectorXd w;
      weights(m_level, w);
      VisitHoleVirtualRows(m_level,
                           [&](const Eigen::Ref<const Eigen::MatrixXd>& rows) {
                             transform.PushMatrix(rows, threadid);
                           });
      transform.A_TDA(w, threadid);
    }
  }
  return transform.getReductionVar();
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
  return result;
}

template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency) const;
template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency) const;

Eigen::MatrixXd RPA::calculate_epsilon_r(std::complex<double> frequency) const {
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
  result.diagonal().array() += 1.0;
  return result;
}

RPA::rpa_eigensolution RPA::Diagonalize_H2p() const {
  const Index lumo = homo_ + 1;
  const Index n_occ = lumo - rpamin_;
  const Index n_unocc = rpamax_ - lumo + 1;
  const Index rpasize = n_occ * n_unocc;

  Eigen::VectorXd AmB = Calculate_H2p_AmB();
  Eigen::MatrixXd ApB = Calculate_H2p_ApB();

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
  const Index rpasize = n_occ * n_unocc;
  vc2index vc = vc2index(0, 0, n_unocc);
  Eigen::VectorXd AmB = Eigen::VectorXd::Zero(rpasize);
  for (Index v = 0; v < n_occ; v++) {
    Index i = vc.I(v, 0);
    AmB.segment(i, n_unocc) =
        energies_.segment(n_occ, n_unocc).array() - energies_(v);
  }
  return AmB;
}

Eigen::MatrixXd RPA::Calculate_H2p_ApB() const {
  const Index lumo = homo_ + 1;
  const Index n_occ = lumo - rpamin_;
  const Index n_unocc = rpamax_ - lumo + 1;
  const Index rpasize = n_occ * n_unocc;
  vc2index vc = vc2index(0, 0, n_unocc);
  Eigen::MatrixXd ApB = Eigen::MatrixXd::Zero(rpasize, rpasize);
#pragma omp parallel for schedule(guided)
  for (Index v2 = 0; v2 < n_occ; v2++) {
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
