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
#include "votca/xtp/make_libint_work.h"
#include "votca/xtp/symmetric_matrix.h"
#include <libint2.hpp>
namespace votca {
namespace xtp {

void ERIs::Initialize(const AOBasis& dftbasis, const AOBasis& auxbasis) {
  _threecenter.Fill(auxbasis, dftbasis);
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

Eigen::MatrixXd ERIs::ComputeSchwarzShells(const AOBasis& basis) const {

  Index noshells = basis.getNumofShells();

  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(noshells, noshells);
  Index nthreads = OPENMP::getMaxThreads();
  std::vector<libint2::Engine> engines(nthreads);
  double epsilon = 0.0;
  engines[0] = libint2::Engine(libint2::Operator::coulomb, basis.getMaxNprim(),
                               static_cast<int>(basis.getMaxL()), 0, epsilon);

  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }

  std::vector<libint2::Shell> shells = basis.GenerateLibintBasis();
  using MatrixLibInt =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

#pragma omp parallel for schedule(dynamic)
  for (Index s1 = 0l; s1 < basis.getNumofShells(); ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const libint2::Engine::target_ptr_vec& buf = engine.results();
    Index n1 = shells[s1].size();

    for (Index s2 = 0l; s2 <= s1; ++s2) {
      Index n2 = shells[s2].size();
      Index n12 = n1 * n2;

      engines[thread_id]
          .compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
              shells[s1], shells[s2], shells[s1], shells[s2]);

      Eigen::Map<const MatrixLibInt> buf_mat(buf[0], n12, n12);

      result(s2, s1) = std::sqrt(buf_mat.cwiseAbs().maxCoeff());
    }
  }
  return result.selfadjointView<Eigen::Upper>();
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

template <bool with_exchange>
std::array<Eigen::MatrixXd, 2> ERIs::Compute4c(const Eigen::MatrixXd& dmat,
                                               double error) const {

  Index nthreads = OPENMP::getMaxThreads();
  std::array<Eigen::MatrixXd, 2> result;
  result[0] = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  if (with_exchange) {
    result[1] = Eigen::MatrixXd::Zero(dmat.rows(), dmat.cols());
  }

  Eigen::MatrixXd& hartree = result[0];
  Eigen::MatrixXd& exchange = result[1];
  Eigen::MatrixXd dnorm_block = ComputeShellBlockNorm(dmat);
  double fock_precision = error;
  // engine precision controls primitive truncation, assume worst-case scenario
  // (all primitive combinations add up constructively)
  Index max_nprim4 = maxnprim_ * maxnprim_ * maxnprim_ * maxnprim_;
  double engine_precision = std::min(fock_precision / dnorm_block.maxCoeff(),
                                     std::numeric_limits<double>::epsilon()) /
                            double(max_nprim4);
  std::cout << engine_precision << std::endl;
  std::vector<libint2::Engine> engines(nthreads);
  engines[0] = libint2::Engine(libint2::Operator::coulomb, int(maxnprim_),
                               int(maxL_), 0);
  engines[0].set_precision(engine_precision);  // shellset-dependent precision
                                               // control will likely break
                                               // positive definiteness
                                               // stick with this simple recipe
  for (Index i = 1; i < nthreads; ++i) {
    engines[i] = engines[0];
  }
  Index nshells = basis_.size();

  //#pragma omp parallel for schedule(dynamic)reduction(+ : hartree)reduction(+
  //: exchange)
  for (Index s1 = 0; s1 < nshells; ++s1) {
    Index thread_id = OPENMP::getThreadId();
    libint2::Engine& engine = engines[thread_id];
    const auto& buf = engine.results();
    Index start_1 = starts_[s1];
    const libint2::Shell& shell1 = basis_[s1];
    Index n1 = shell1.size();

    auto sp12_iter = shellpairdata_[s1].begin();
    for (Index s2 : shellpairs_[s1]) {
      Index start_2 = starts_[s2];
      const libint2::Shell& shell2 = basis_[s2];
      Index n2 = shell2.size();
      double dnorm_12 = dnorm_block(s1, s2);
      const libint2::ShellPair* sp12 = &(*sp12_iter);
      ++sp12_iter;

      for (Index s3 = 0; s3 <= s1; ++s3) {

        Index start_3 = starts_[s3];
        const libint2::Shell& shell3 = basis_[s3];
        Index n3 = shell3.size();
        auto sp34_iter = shellpairdata_[s3].begin();
        double dnorm_123 = std::max(dnorm_block(s1, s3),
                                    std::max(dnorm_block(s2, s3), dnorm_12));
        Index s4max = (s1 == s3) ? s2 : s3;
        for (Index s4 : shellpairs_[s3]) {
          if (s4 > s4max) {
            break;
          }  // for each s3, s4 are stored in monotonically increasing
             // order

          const libint2::ShellPair* sp34 = &(*sp34_iter);
          // must update the iter even if going to skip s4
          ++sp34_iter;
          double dnorm_1234 =
              std::max(dnorm_block(s1, s4),
                       std::max(dnorm_block(s2, s4),
                                std::max(dnorm_block(s3, s4), dnorm_123)));

          if (dnorm_1234 * schwarzscreen_(s1, s2) * schwarzscreen_(s3, s4) <
              fock_precision) {
            continue;
          }

          const libint2::Shell& shell4 = basis_[s4];
          engine
              .compute2<libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
                  shell1, shell2, shell3, shell4, sp12, sp34);
          const auto* buf_1234 = buf[0];
          if (buf_1234 == nullptr) {
            continue;  // if all integrals screened out, skip to next quartet
          }
          Index start_4 = starts_[s4];
          Index n4 = shell4.size();
          Index s12_deg = (s1 == s2) ? 1 : 2;
          Index s34_deg = (s3 == s4) ? 1 : 2;
          Index s12_34_deg = (s1 == s3) ? (s2 == s4 ? 1 : 2) : 2;
          Index s1234_deg = s12_deg * s34_deg * s12_34_deg;

          for (Index f1 = 0, f1234 = 0; f1 != n1; ++f1) {
            const Index bf1 = f1 + start_1;
            for (Index f2 = 0; f2 != n2; ++f2) {
              const Index bf2 = f2 + start_2;
              for (Index f3 = 0; f3 != n3; ++f3) {
                const Index bf3 = f3 + start_3;
                for (Index f4 = 0; f4 != n4; ++f4, ++f1234) {
                  const Index bf4 = f4 + start_4;

                  const double value = buf_1234[f1234];

                  const double value_scal_by_deg = value * double(s1234_deg);

                  hartree(bf1, bf2) += dmat(bf3, bf4) * value_scal_by_deg;
                  hartree(bf3, bf4) += dmat(bf1, bf2) * value_scal_by_deg;
                  if (with_exchange) {
                    exchange(bf1, bf3) -= dmat(bf2, bf4) * value_scal_by_deg;
                    exchange(bf2, bf4) -= dmat(bf1, bf3) * value_scal_by_deg;
                    exchange(bf1, bf4) -= dmat(bf2, bf3) * value_scal_by_deg;
                    exchange(bf2, bf3) -= dmat(bf1, bf4) * value_scal_by_deg;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  std::array<Eigen::MatrixXd, 2> result2;
  result2[0] = 0.25 * (result[0] + result[0].transpose());
  if (with_exchange) {
    result2[1] = 0.125 * (result[1] + result[1].transpose());
  }
  return result2;
}

template std::array<Eigen::MatrixXd, 2> ERIs::Compute4c<true>(
    const Eigen::MatrixXd& dmat, double error) const;
template std::array<Eigen::MatrixXd, 2> ERIs::Compute4c<false>(
    const Eigen::MatrixXd& dmat, double error) const;

Eigen::MatrixXd ERIs::CalculateERIs_3c(const Eigen::MatrixXd& DMAT) const {
  Eigen::MatrixXd ERIs2 = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());
  Symmetric_Matrix dmat_sym = Symmetric_Matrix(DMAT);
#pragma omp parallel for schedule(guided) reduction(+ : ERIs2)
  for (Index i = 0; i < _threecenter.size(); i++) {
    const Symmetric_Matrix& threecenter = _threecenter[i];
    // Trace over prod::DMAT,I(l)=componentwise product over
    const double factor = threecenter.TraceofProd(dmat_sym);
    Eigen::SelfAdjointView<Eigen::MatrixXd, Eigen::Upper> m =
        ERIs2.selfadjointView<Eigen::Upper>();
    threecenter.AddtoEigenUpperMatrix(m, factor);
  }

  return ERIs2.selfadjointView<Eigen::Upper>();
}

Eigen::MatrixXd ERIs::CalculateEXX_dmat(const Eigen::MatrixXd& DMAT) const {

  Eigen::MatrixXd EXX = Eigen::MatrixXd::Zero(DMAT.rows(), DMAT.cols());

#pragma omp parallel for schedule(guided) reduction(+ : EXX)
  for (Index i = 0; i < _threecenter.size(); i++) {
    const Eigen::MatrixXd threecenter = _threecenter[i].UpperMatrix();
    EXX += threecenter.selfadjointView<Eigen::Upper>() * DMAT *
           threecenter.selfadjointView<Eigen::Upper>();
  }
  return EXX;
}

Eigen::MatrixXd ERIs::CalculateEXX_mos(const Eigen::MatrixXd& occMos) const {

  Eigen::MatrixXd EXX = Eigen::MatrixXd::Zero(occMos.rows(), occMos.rows());

#pragma omp parallel for schedule(guided) reduction(+ : EXX)
  for (Index i = 0; i < _threecenter.size(); i++) {
    const Eigen::MatrixXd TCxMOs_T =
        occMos.transpose() *
        _threecenter[i].UpperMatrix().selfadjointView<Eigen::Upper>();
    EXX += TCxMOs_T.transpose() * TCxMOs_T;
  }
  return 2 * EXX;
}

}  // namespace xtp
}  // namespace votca
