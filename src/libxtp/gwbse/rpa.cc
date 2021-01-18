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
#include "votca/xtp/rpa.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/threecenter.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

void RPA::UpdateRPAInputEnergies(const Eigen::VectorXd& dftenergies,
                                 const Eigen::VectorXd& gwaenergies,
                                 Index qpmin) {
  Index rpatotal = _rpamax - _rpamin + 1;
  _energies = dftenergies.segment(_rpamin, rpatotal);
  Index gwsize = Index(gwaenergies.size());

  _energies.segment(qpmin - _rpamin, gwsize) = gwaenergies;

  ShiftUncorrectedEnergies(dftenergies, qpmin, gwsize);
}

// Shifts energies of levels that are not QP corrected but
// used in the RPA:
// between rpamin and qpmin: by maximum abs of explicit QP corrections
//                           from qpmin to HOMO
// between qpmax and rpamax: by maximum abs of explicit QP corrections
//                           from LUMO to qpmax
void RPA::ShiftUncorrectedEnergies(const Eigen::VectorXd& dftenergies,
                                   Index qpmin, Index gwsize) {

  Index lumo = _homo + 1;
  Index qpmax = qpmin + gwsize - 1;

  // get max abs QP corrections for occupied/virtual levels
  double max_correction_occ = getMaxCorrection(dftenergies, qpmin, _homo);
  double max_correction_virt = getMaxCorrection(dftenergies, lumo, qpmax);

  // shift energies
  Index offset_qp_virt = qpmax + 1 - _rpamin;
  _energies.segment(0, qpmin - _rpamin).array() -= max_correction_occ;
  _energies.segment(offset_qp_virt, _rpamax - qpmax).array() +=
      max_correction_virt;
}

double RPA::getMaxCorrection(const Eigen::VectorXd& dftenergies, Index min,
                             Index max) {

  Index range = max - min + 1;
  Eigen::VectorXd corrections = _energies.segment(min - _rpamin, range) -
                                dftenergies.segment(min - _rpamin, range);

  return (corrections.cwiseAbs()).maxCoeff();
}

template <bool imag>
Eigen::MatrixXd RPA::calculate_epsilon(double frequency) const {
  const Index size = _Mmn.auxsize();

  const Index lumo = _homo + 1;
  const Index n_occ = lumo - _rpamin;
  const Index n_unocc = _rpamax - lumo + 1;
  const double freq2 = frequency * frequency;
  const double eta2 = _eta * _eta;

  OpenMP_CUDA transform;
  transform.createTemporaries(n_unocc, size);

#pragma omp parallel for schedule(dynamic)
  for (Index m_level = 0; m_level < n_occ; m_level++) {
    const double qp_energy_m = _energies(m_level);

    const Eigen::MatrixXd Mmn_RPA = _Mmn[m_level].bottomRows(n_unocc);

    const Eigen::ArrayXd deltaE = _energies.tail(n_unocc).array() - qp_energy_m;
    Eigen::VectorXd denom;
    if (imag) {
      denom = 4 * deltaE / (deltaE.square() + freq2);
    } else {
      Eigen::ArrayXd deltEf = deltaE - frequency;
      Eigen::ArrayXd sum = deltEf / (deltEf.square() + eta2);
      deltEf = deltaE + frequency;
      sum += deltEf / (deltEf.square() + eta2);
      denom = 2 * sum;
    }

    transform.A_TDA(Mmn_RPA, denom);
  }
  Eigen::MatrixXd result = transform.A_TDA_result();
  result.diagonal().array() += 1.0;
  return result;
}

template Eigen::MatrixXd RPA::calculate_epsilon<true>(double frequency) const;
template Eigen::MatrixXd RPA::calculate_epsilon<false>(double frequency) const;

Eigen::MatrixXd RPA::calculate_epsilon_r(std::complex<double> frequency) const {

  const Index size = _Mmn.auxsize();

  const Index lumo = _homo + 1;
  const Index n_occ = lumo - _rpamin;
  const Index n_unocc = _rpamax - lumo + 1;
  OpenMP_CUDA transform;
  transform.createTemporaries(n_unocc, size);

#pragma omp parallel for schedule(dynamic)
  for (Index m_level = 0; m_level < n_occ; m_level++) {

    const double qp_energy_m = _energies(m_level);
    const Eigen::MatrixXd Mmn_RPA = _Mmn[m_level].bottomRows(n_unocc);
    const Eigen::ArrayXd deltaE = _energies.tail(n_unocc).array() - qp_energy_m;

    Eigen::ArrayXd deltaEm = frequency.real() - deltaE;
    Eigen::ArrayXd deltaEp = frequency.real() + deltaE;

    double sigma_1 = std::pow(frequency.imag() + _eta, 2);
    double sigma_2 = std::pow(frequency.imag() - _eta, 2);

    Eigen::VectorXd chi =
        deltaEm * (deltaEm.cwiseAbs2() + sigma_1).cwiseInverse() -
        deltaEp * (deltaEp.cwiseAbs2() + sigma_2).cwiseInverse();
    transform.A_TDA(Mmn_RPA, chi);
  }
  Eigen::MatrixXd result = -2 * transform.A_TDA_result();
  result.diagonal().array() += 1.0;
  return result;
}

RPA::rpa_eigensolution RPA::Diagonalize_H2p() const {
  const Index lumo = _homo + 1;
  const Index n_occ = lumo - _rpamin;
  const Index n_unocc = _rpamax - lumo + 1;
  const Index rpasize = n_occ * n_unocc;

  Eigen::VectorXd AmB = Calculate_H2p_AmB();
  Eigen::MatrixXd ApB = Calculate_H2p_ApB();

  RPA::rpa_eigensolution sol;
  sol.ERPA_correlation = -0.25 * (ApB.trace() + AmB.sum());

  // C = AmB^1/2 * ApB * AmB^1/2
  Eigen::MatrixXd& C = ApB;
  C.applyOnTheLeft(AmB.cwiseSqrt().asDiagonal());
  C.applyOnTheRight(AmB.cwiseSqrt().asDiagonal());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es = Diagonalize_H2p_C(C);

  // Do not remove this line! It has to be there for MKL to not crash
  sol.omega = Eigen::VectorXd::Zero(es.eigenvalues().size());
  sol.omega = es.eigenvalues().cwiseSqrt();
  sol.ERPA_correlation += 0.5 * sol.omega.sum();

  XTP_LOG(Log::info, _log) << TimeStamp()
                           << " Lowest neutral excitation energy (eV): "
                           << tools::conv::hrt2ev * sol.omega.minCoeff()
                           << std::flush;

  // RPA correlation energy calculated from Eq.9 of J. Chem. Phys. 132, 234114
  // (2010)
  XTP_LOG(Log::error, _log)
      << TimeStamp()
      << " RPA correlation energy (Hartree): " << sol.ERPA_correlation
      << std::flush;

  sol.XpY = Eigen::MatrixXd(rpasize, rpasize);

  Eigen::VectorXd AmB_sqrt = AmB.cwiseSqrt();
  Eigen::VectorXd Omega_sqrt_inv = sol.omega.cwiseSqrt().cwiseInverse();
  for (int s = 0; s < rpasize; s++) {
    sol.XpY.col(s) =
        Omega_sqrt_inv(s) * AmB_sqrt.cwiseProduct(es.eigenvectors().col(s));
  }

  return sol;
}

Eigen::VectorXd RPA::Calculate_H2p_AmB() const {
  const Index lumo = _homo + 1;
  const Index n_occ = lumo - _rpamin;
  const Index n_unocc = _rpamax - lumo + 1;
  const Index rpasize = n_occ * n_unocc;
  vc2index vc = vc2index(0, 0, n_unocc);
  Eigen::VectorXd AmB = Eigen::VectorXd::Zero(rpasize);
  for (Index v = 0; v < n_occ; v++) {
    Index i = vc.I(v, 0);
    AmB.segment(i, n_unocc) =
        _energies.segment(n_occ, n_unocc).array() - _energies(v);
  }
  return AmB;
}

Eigen::MatrixXd RPA::Calculate_H2p_ApB() const {
  const Index lumo = _homo + 1;
  const Index n_occ = lumo - _rpamin;
  const Index n_unocc = _rpamax - lumo + 1;
  const Index rpasize = n_occ * n_unocc;
  const Index auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, n_unocc);
  Eigen::MatrixXd ApB = Eigen::MatrixXd::Zero(rpasize, rpasize);
#pragma omp parallel for schedule(guided)
  for (Index v2 = 0; v2 < n_occ; v2++) {
    Index i2 = vc.I(v2, 0);
    const Eigen::MatrixXd Mmn_v2T =
        _Mmn[v2].block(n_occ, 0, n_unocc, auxsize).transpose();
    for (Index v1 = v2; v1 < n_occ; v1++) {
      Index i1 = vc.I(v1, 0);
      // Multiply with factor 2 to sum over both (identical) spin states
      ApB.block(i1, i2, n_unocc, n_unocc) =
          2 * 2 * _Mmn[v1].block(n_occ, 0, n_unocc, auxsize) * Mmn_v2T;
    }
  }
  ApB.diagonal() += Calculate_H2p_AmB();
  return ApB;
}

Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> RPA::Diagonalize_H2p_C(
    const Eigen::MatrixXd& C) const {
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Diagonalizing two-particle Hamiltonian "
      << std::flush;
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);  // Uses lower triangle
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Diagonalization done " << std::flush;
  double minCoeff = es.eigenvalues().minCoeff();
  if (minCoeff <= 0.0) {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Detected non-positive eigenvalue: " << minCoeff
        << std::flush;
    throw std::runtime_error("Detected non-positive eigenvalue.");
  }
  return es;
}

}  // namespace xtp
}  // namespace votca
