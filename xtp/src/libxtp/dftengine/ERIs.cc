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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include <algorithm>
#include <cmath>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/ERIs.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/symmetric_matrix.h"
namespace votca {
namespace xtp {

void ERIs::Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis,
                      double pair_threshold) {
  threecenter_.Fill(auxbasis, dftbasis, pair_threshold);
  return;
}

void ERIs::Initialize_4c(const AOBasis& dftbasis) {

  basis_ = dftbasis.GenerateLibintBasis();
  shellpairs_ = dftbasis.ComputeShellPairs();
  starts_ = dftbasis.getMapToBasisFunctions();
  maxnprim_ = dftbasis.getMaxNprim();
  maxL_ = dftbasis.getMaxL();

  shellpairdata_ = ComputeShellPairData(basis_, shellpairs_);

  schwarzscreen_ = ComputeSchwarzShells(dftbasis);
  return;
}

std::vector<std::vector<libint2::ShellPair>> ERIs::ComputeShellPairData(
    const std::vector<libint2::Shell>& basis,
    const std::vector<std::vector<Index>>& shellpairs) const {
  std::vector<std::vector<libint2::ShellPair>> shellpairdata(basis.size());
  const double ln_max_engine_precision =
      std::log(std::numeric_limits<double>::epsilon() * 1e-10);

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0; s1 < Index(shellpairs.size()); s1++) {
    for (Index s2 : shellpairs[s1]) {
      shellpairdata[s1].emplace_back(
          libint2::ShellPair(basis[s1], basis[s2], ln_max_engine_precision));
    }
  }
  return shellpairdata;
}

Eigen::MatrixXd ERIs::ComputeShellBlockNorm(const Eigen::MatrixXd& dmat) const {
  Eigen::MatrixXd result =
      Eigen::MatrixXd::Zero(starts_.size(), starts_.size());
#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0l; s1 < Index(basis_.size()); ++s1) {
    Index bf1 = starts_[s1];
    Index n1 = basis_[s1].size();
    for (Index s2 = 0l; s2 <= s1; ++s2) {
      Index bf2 = starts_[s2];
      Index n2 = basis_[s2].size();

      result(s2, s1) = dmat.block(bf2, bf1, n2, n1).cwiseAbs().maxCoeff();
    }
  }
  return result.selfadjointView<Eigen::Upper>();
}

// J[D]_mu nu = sum_P B_P,mu nu c_P with c_P = sum_kl B_P,kl D_kl: two
// matrix-vector products with the stored pair x aux tensor.
Eigen::MatrixXd ERIs::CalculateERIs_3c(const Eigen::MatrixXd& DMAT) const {
  assert(threecenter_.size() > 0 &&
         "Please call Initialize before running this");
  const Eigen::MatrixXd& data = threecenter_.Data();
  // lower triangle of DMAT, as before
  const Eigen::VectorXd dpacked = threecenter_.PackWeighted(DMAT);
  Eigen::VectorXd c(data.cols());
#pragma omp parallel for schedule(static)
  for (Index P = 0; P < data.cols(); ++P) {
    c(P) = data.col(P).dot(dpacked);
  }
  Eigen::VectorXd jpacked(data.rows());
  const Index nblocks = std::max<Index>(1, OPENMP::getMaxThreads()) * 8;
  const Index blocksize = (data.rows() + nblocks - 1) / nblocks;
#pragma omp parallel for schedule(static)
  for (Index blk = 0; blk < nblocks; ++blk) {
    const Index row = blk * blocksize;
    const Index len = std::min(blocksize, data.rows() - row);
    if (len > 0) {
      jpacked.segment(row, len).noalias() = data.middleRows(row, len) * c;
    }
  }
  return threecenter_.Unpack(jpacked);
}

// K[D] = -sum_P B_P D B_P is linear in D. Writing D = X X^T - Y Y^T from
// its eigendecomposition turns every term into products with N x rank(D)
// factors -- the work of the occupied-MO route -- instead of two N^3
// products per aux function. Any symmetric D is handled exactly, including
// indefinite ones (incremental differences); the gain comes from the low
// rank of the densities the SCF produces: a density mixed from idempotent
// ones has at most a few times as many nonzero eigenvalues as there are
// occupied orbitals.
Eigen::MatrixXd ERIs::CalculateEXX_dmat(const Eigen::MatrixXd& DMAT) const {
  assert(threecenter_.size() > 0 &&
         "Please call Initialize before running this");
  Eigen::MatrixXd EXX = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      0.5 * (DMAT + DMAT.transpose()));
  const Eigen::VectorXd& lambda = es.eigenvalues();
  if (lambda.size() == 0) {
    return EXX;
  }
  // Eigenvalues below this are rounding noise of a lower-rank matrix.
  const double cutoff = kDensityRankCutoff * lambda.cwiseAbs().maxCoeff();
  std::vector<Index> positive;
  std::vector<Index> negative;
  for (Index k = 0; k < lambda.size(); ++k) {
    if (lambda(k) > cutoff) {
      positive.push_back(k);
    } else if (lambda(k) < -cutoff) {
      negative.push_back(k);
    }
  }
  const Index npos = Index(positive.size());
  const Index nneg = Index(negative.size());
  if (npos + nneg == 0) {
    return EXX;
  }
  // columns: sqrt(|lambda_k|) u_k, positive eigenvalues first
  Eigen::MatrixXd factors(DMAT.rows(), npos + nneg);
  for (Index c = 0; c < npos; ++c) {
    const Index k = positive[c];
    factors.col(c) = std::sqrt(lambda(k)) * es.eigenvectors().col(k);
  }
  for (Index c = 0; c < nneg; ++c) {
    const Index k = negative[c];
    factors.col(npos + c) = std::sqrt(-lambda(k)) * es.eigenvectors().col(k);
  }
  return ExchangeFromFactors(factors, npos);
}

Eigen::MatrixXd ERIs::CalculateEXX_mos(const Eigen::MatrixXd& occMos) const {
  assert(threecenter_.size() > 0 &&
         "Please call Initialize before running this");
  // D = 2 C C^T
  return ExchangeFromFactors(std::sqrt(2.0) * occMos, occMos.cols());
}

// Each aux function contributes (B_P X)(B_P X)^T. Per thread, B_P is unpacked
// into one reused buffer, B_P X is formed with a plain matrix product, and
// the rows of several aux functions are stacked so that the N x N
// accumulator is touched by one rank update per batch instead of once per
// aux function.
Eigen::MatrixXd ERIs::ExchangeFromFactors(const Eigen::MatrixXd& factors,
                                          Index npos) const {
  const Index n = factors.rows();
  const Index nneg = factors.cols() - npos;
  Eigen::MatrixXd EXX = Eigen::MatrixXd::Zero(n, n);
  if (factors.cols() == 0) {
    return EXX;
  }
  const Index batch =
      std::max<Index>(1, kExchangeBatchRows / std::max(npos, nneg));
  const Eigen::MatrixXd xt = factors.leftCols(npos).transpose();
  const Eigen::MatrixXd yt = factors.rightCols(nneg).transpose();

#pragma omp parallel
  {
    Eigen::MatrixXd local = Eigen::MatrixXd::Zero(n, n);
    // entries of screened pairs are never written and stay zero
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd tx(batch * npos, n);
    Eigen::MatrixXd ty(batch * nneg, n);
    Index filled = 0;
    auto flush = [&]() {
      if (filled == 0) {
        return;
      }
      local.selfadjointView<Eigen::Lower>().rankUpdate(
          tx.topRows(filled * npos).transpose(), -1.0);
      if (nneg > 0) {
        local.selfadjointView<Eigen::Lower>().rankUpdate(
            ty.topRows(filled * nneg).transpose(), 1.0);
      }
      filled = 0;
    };
#pragma omp for schedule(dynamic)
    for (Index i = 0; i < threecenter_.size(); i++) {
      threecenter_.FillFullMatrix(i, b);
      if (npos > 0) {
        tx.middleRows(filled * npos, npos).noalias() = xt * b;
      }
      if (nneg > 0) {
        ty.middleRows(filled * nneg, nneg).noalias() = yt * b;
      }
      if (++filled == batch) {
        flush();
      }
    }
    flush();
#pragma omp critical
    { EXX.triangularView<Eigen::Lower>() += local; }
  }
  return EXX.selfadjointView<Eigen::Lower>();
}

}  // namespace xtp
}  // namespace votca
