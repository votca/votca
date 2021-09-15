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

// Local VOTCA includes
#include "votca/xtp/ERIs.h"
#include "votca/xtp/aobasis.h"
#include "votca/xtp/symmetric_matrix.h"
namespace votca {
namespace xtp {

void ERIs::Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis) {
  threecenter_.Fill(auxbasis, dftbasis);
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

Eigen::MatrixXd ERIs::CalculateERIs_3c(const Eigen::MatrixXd& DMAT) const {
  assert(threecenter_.size() > 0 &&
         "Please call Initialize before running this");
  Eigen::MatrixXd ERIs2 = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
  Symmetric_Matrix dmat_sym = Symmetric_Matrix(DMAT);
#pragma omp parallel for schedule(guided) reduction(+ : ERIs2)
  for (Index i = 0; i < threecenter_.size(); i++) {
    const Symmetric_Matrix& threecenter = threecenter_[i];
    // Trace over prod::DMAT,I(l)=componentwise product over
    const double factor = threecenter.TraceofProd(dmat_sym);
    Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Upper> m =
        ERIs2.selfadjointView<Eigen::Upper>();
    threecenter.AddtoEigenUpperMatrix(m, factor);
  }

  return ERIs2.selfadjointView<Eigen::Upper>();
}

Eigen::MatrixXd ERIs::CalculateEXX_dmat(const Eigen::MatrixXd& DMAT) const {
  assert(threecenter_.size() > 0 &&
         "Please call Initialize before running this");
  Eigen::MatrixXd EXX = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

#pragma omp parallel for schedule(guided) reduction(+ : EXX)
  for (Index i = 0; i < threecenter_.size(); i++) {
    const Eigen::MatrixXd threecenter = threecenter_[i].UpperMatrix();
    EXX -= threecenter.selfadjointView<Eigen::Upper>() * DMAT *
           threecenter.selfadjointView<Eigen::Upper>();
  }
  return EXX;
}

Eigen::MatrixXd ERIs::CalculateEXX_mos(const Eigen::MatrixXd& occMos) const {
  assert(threecenter_.size() > 0 &&
         "Please call Initialize before running this");
  Eigen::MatrixXd EXX = Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows());

#pragma omp parallel for schedule(guided) reduction(+ : EXX)
  for (Index i = 0; i < threecenter_.size(); i++) {
    const Eigen::MatrixXd TCxMOs_T =
        occMos.transpose() *
        threecenter_[i].UpperMatrix().selfadjointView<Eigen::Upper>();
    EXX -= TCxMOs_T.transpose() * TCxMOs_T;
  }
  return 2 * EXX;
}

}  // namespace xtp
}  // namespace votca
