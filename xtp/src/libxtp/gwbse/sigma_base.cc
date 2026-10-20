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

// Standard includes
#include <cmath>

// Third party includes
#include <boost/math/constants/constants.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/sigma_base.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

Eigen::MatrixXd Sigma_base::PairContraction(
    Index nrows,
    const std::function<void(Index, Eigen::MatrixXd&)>& makeX) const {
  const Index qpmin = opt_.qpmin - opt_.rpamin;
  const Index auxsize = Mmn_.auxsize();
  // X of a block of levels kept at once, and the chunk of all M_n gathered
  // for one product; the block size sets how often the integrals are read.
  const double kBlockBytes = block_bytes_;
  const double kChunkBytes = chunk_bytes_;
  const Index block = std::max<Index>(
      1, std::min<Index>(qptotal_, Index(kBlockBytes / (8.0 * double(nrows) *
                                                        double(auxsize)))));
  const Index chunk = std::max<Index>(
      1, std::min<Index>(auxsize, Index(kChunkBytes / (8.0 * double(nrows) *
                                                       double(qptotal_)))));
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
  std::vector<Eigen::MatrixXd> X(static_cast<std::size_t>(block));
  // chunk buffers kept across chunks (no reallocation per chunk)
  Eigen::MatrixXd Mc;
  Eigen::MatrixXd Xc;
  for (Index m0 = 0; m0 < qptotal_; m0 += block) {
    const Index nb = std::min(block, qptotal_ - m0);
#pragma omp parallel for schedule(dynamic)
    for (Index j = 0; j < nb; ++j) {
      makeX(m0 + j, X[std::size_t(j)]);
    }
    for (Index s0 = 0; s0 < auxsize; s0 += chunk) {
      const Index ns = std::min(chunk, auxsize - s0);
      Mc.resize(nrows * ns, qptotal_);
      Xc.resize(nrows * ns, nb);
#pragma omp parallel for schedule(static)
      for (Index n = 0; n < qptotal_; ++n) {
        Eigen::Map<Eigen::MatrixXd>(Mc.col(n).data(), nrows, ns) =
            Mmn_[n + qpmin].block(0, s0, nrows, ns);
      }
      for (Index j = 0; j < nb; ++j) {
        Eigen::Map<Eigen::MatrixXd>(Xc.col(j).data(), nrows, ns) =
            X[std::size_t(j)].middleCols(s0, ns);
      }
      A.middleRows(m0, nb).noalias() += Xc.transpose() * Mc;
    }
  }
  return A;
}

Sigma_base::DiagPoleForm Sigma_base::PoleForm(Index) const {
  throw std::runtime_error(
      "The pole form of Sigma_c is only available for sigma_integrator ppm "
      "and exact");
}

Sigma_base::ScreeningData Sigma_base::ScreeningDump(
    const Eigen::VectorXd&) const {
  throw std::runtime_error(
      "Screening data are only available for sigma_integrator ppm and exact");
}

void Sigma_base::FrameSums(const Eigen::MatrixXd& a,
                           const Eigen::VectorXd& omegas,
                           const Eigen::VectorXd& nu, ScreeningData& data) {
  const Index S = omegas.size();
  // columns: Omega^-1, Omega^0, Omega^1, then Omega/(nu^2+Omega^2)
  Eigen::MatrixXd g(S, 3 + nu.size());
  for (Index s = 0; s < S; ++s) {
    const double w = omegas(s);
    g(s, 0) = 1.0 / w;
    g(s, 1) = 1.0;
    g(s, 2) = w;
    for (Index j = 0; j < nu.size(); ++j) {
      g(s, 3 + j) = -w / (nu(j) * nu(j) + w * w);
    }
  }
  const Eigen::MatrixXd sums = a * g;
  data.frame_moments = sums.leftCols(3);
  data.frame_eps_inv_imag = sums.rightCols(nu.size());
}

Eigen::MatrixXd Sigma_base::CalcExchangeMatrix() const {
  const Index occlevel = opt_.homo - opt_.rpamin + 1;
  const Index qpmin = opt_.qpmin - opt_.rpamin;
  // Exchange is with the bare v, which is M M^T only as long as the
  // auxiliary frame is orthogonal. Once the index is dressed for an
  // environment (TCMatrix_gwbse::DressAuxIndex) it takes the bare
  // interaction in the current frame as an explicit kernel.
  const bool bare = Mmn_.AuxFrameIsOrthogonal();
  Eigen::MatrixXd v;
  if (!bare) {
    v = Mmn_.ToCurrentAuxFrame(
        Eigen::MatrixXd::Identity(Mmn_.auxsize(), Mmn_.auxsize()));
  }
  // Sigma_x(m,n) = - sum_{i occ} sum_P M_m(i,P) [v] M_n(i,P)
  const Eigen::MatrixXd A =
      PairContraction(occlevel, [&](Index m, Eigen::MatrixXd& X) {
        if (bare) {
          X = Mmn_[m + qpmin].topRows(occlevel);
        } else {
          X = Mmn_[m + qpmin].topRows(occlevel) * v;
        }
      });
  return -0.5 * (A + A.transpose());
}

Eigen::MatrixXd Sigma_base::CalcHartreeMatrix(const Eigen::MatrixXd& R) const {
  const Index occlevel = opt_.homo - opt_.rpamin + 1;
  const Index qpmin = opt_.qpmin - opt_.rpamin;
  // d_P = sum_{i occ} M_ii,P (the density in the auxiliary basis)
  Eigen::VectorXd d = Eigen::VectorXd::Zero(Mmn_.auxsize());
  for (Index i = 0; i < occlevel; ++i) {
    d += Mmn_[i].row(i).transpose();
  }
  // Kernel 1 (bare v) or 1 + R (u = v + v_reac), given in the fill-time
  // frame and brought into the current one (see CalcExchangeMatrix). On
  // dressed integrals F F^T = 1 + R, so the screened kernel is close to the
  // identity there; it is still taken explicitly, which keeps this right
  // in any frame.
  if (R.size() > 0) {
    if (R.rows() != Mmn_.auxsize() || R.cols() != Mmn_.auxsize()) {
      throw std::runtime_error(
          "Sigma_base::CalcHartreeMatrix: reaction field is " +
          std::to_string(R.rows()) + "x" + std::to_string(R.cols()) +
          ", auxiliary basis has " + std::to_string(Mmn_.auxsize()) +
          " functions.");
    }
    d = Mmn_.ToCurrentAuxFrame(
            Eigen::MatrixXd::Identity(Mmn_.auxsize(), Mmn_.auxsize()) + R) *
        d;
  } else if (!Mmn_.AuxFrameIsOrthogonal()) {
    d = Mmn_.ToCurrentAuxFrame(
            Eigen::MatrixXd::Identity(Mmn_.auxsize(), Mmn_.auxsize())) *
        d;
  }
  Eigen::MatrixXd J(qptotal_, qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index m = 0; m < qptotal_; ++m) {
    J.col(m) = 2.0 * Mmn_[m + qpmin].middleRows(qpmin, qptotal_) * d;
  }
  return 0.5 * (J + J.transpose());
}

Eigen::MatrixXd Sigma_base::CalcReactionFieldMatrix(
    const Eigen::MatrixXd& R) const {
  const Eigen::MatrixXd Rc = Mmn_.ToCurrentAuxFrame(R);
  const Index occlevel = opt_.homo - opt_.rpamin + 1;
  const Index qpmin = opt_.qpmin - opt_.rpamin;
  // Row l of X_m is s_l M_m(l,:) R: the occupation sign folded in once.
  const Eigen::MatrixXd A =
      PairContraction(Mmn_.nsize(), [&](Index m, Eigen::MatrixXd& X) {
        X = Mmn_[m + qpmin] * Rc;
        X.topRows(occlevel) *= -1.0;
      });
  return 0.25 * (A + A.transpose());
}

Eigen::VectorXd Sigma_base::CalcCorrelationDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::VectorXd result = Eigen::VectorXd::Zero(qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level = 0; gw_level < qptotal_; gw_level++) {
    result(gw_level) =
        CalcCorrelationDiagElement(gw_level, frequencies[gw_level]);
  }
  return result;
}

Eigen::MatrixXd Sigma_base::CalcCorrelationSRG(
    const Eigen::VectorXd& frequencies, double s) const {
  if (!HasPoleForm()) {
    throw std::runtime_error(
        "SRG-regularized QSGW needs a self-energy with a pole form (exact or "
        "ppm)");
  }
  if (!(s > 0.0)) {
    throw std::runtime_error("SRG flow parameter s must be positive");
  }
  if (frequencies.size() != qptotal_) {
    throw std::runtime_error("CalcCorrelationSRG: one frequency per level");
  }
  std::vector<DiagPoleForm> forms;
  forms.reserve(std::size_t(qptotal_));
  for (Index n = 0; n < qptotal_; ++n) {
    forms.push_back(PoleForm(n));
  }
  const DiagPoleForm& first = forms.front();
  const ExactPoles& poles = first.poles;
  const Index B = poles.nlevels();
  const Index S = poles.nmodes();
  const Index P = poles.npoles();
  // pole positions and factors, p = m + B s
  Eigen::ArrayXd z(P);
  Eigen::ArrayXd c = Eigen::ArrayXd::Ones(P);
  for (Index mode = 0; mode < S; ++mode) {
    for (Index m = 0; m < B; ++m) {
      z(m + B * mode) = poles.energies(m) +
                        ((m < poles.n_occ) ? -1.0 : 1.0) * poles.omegas(mode);
      if (first.mode_factors != nullptr) {
        c(m + B * mode) = first.mode_factors[mode];
      }
    }
  }
  Eigen::MatrixXd result(qptotal_, qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index n1 = 0; n1 < qptotal_; ++n1) {
    const Eigen::Map<const Eigen::ArrayXd> r1(forms[std::size_t(n1)].residues,
                                              P);
    const Eigen::ArrayXd d1 = frequencies(n1) - z;
    const Eigen::ArrayXd cr1 = c * r1;
    for (Index n2 = 0; n2 <= n1; ++n2) {
      const Eigen::Map<const Eigen::ArrayXd> r2(forms[std::size_t(n2)].residues,
                                                P);
      const Eigen::ArrayXd d2 = frequencies(n2) - z;
      const Eigen::ArrayXd D = d1.square() + d2.square();
      // (1 - exp(-s D)) / D, with its limit s at D -> 0
      const Eigen::ArrayXd damp = (s * D < 1e-12)
                                      .select(Eigen::ArrayXd::Constant(P, s),
                                              -(-s * D).unaryExpr([](double x) {
                                                return std::expm1(x);
                                              }) / D.max(1e-300));
      const double value =
          first.prefactor * (cr1 * r2 * (d1 + d2) * damp).sum();
      result(n1, n2) = value;
      result(n2, n1) = value;
    }
  }
  return result;
}

Eigen::MatrixXd Sigma_base::CalcCorrelationOffDiag(
    const Eigen::VectorXd& frequencies) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(qptotal_, qptotal_);
#pragma omp parallel for schedule(dynamic)
  for (Index gw_level1 = 0; gw_level1 < qptotal_; gw_level1++) {
    for (Index gw_level2 = gw_level1 + 1; gw_level2 < qptotal_; gw_level2++) {
      double sigma_c = CalcCorrelationOffDiagElement(
          gw_level1, gw_level2, frequencies[gw_level1], frequencies[gw_level2]);
      result(gw_level2, gw_level1) = sigma_c;
    }
  }
  result = result.selfadjointView<Eigen::Lower>();
  return result;
}

}  // namespace xtp
}  // namespace votca
