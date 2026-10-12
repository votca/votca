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
