/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE screening_kernels_test

// Standard includes
#include <iostream>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/accelerate_lapack.h"
#include "votca/xtp/screening_kernels.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(screening_kernels_test)

namespace {
// blocks of different heights, weights of both signs and exact zeros
struct Blocks {
  std::vector<Eigen::MatrixXd> R;
  std::vector<Eigen::VectorXd> w;
  Index ncols;

  Blocks(Index nblocks, Index ncols_in) : ncols(ncols_in) {
    for (Index m = 0; m < nblocks; ++m) {
      const Index rows = 3 + (5 * m) % 11;
      R.push_back(Eigen::MatrixXd::Random(rows, ncols));
      Eigen::VectorXd wm = Eigen::VectorXd::Random(rows);
      wm(0) = 0.0;
      w.push_back(wm);
    }
  }

  Eigen::MatrixXd Reference() const {
    Eigen::MatrixXd S = Eigen::MatrixXd::Zero(ncols, ncols);
    for (std::size_t m = 0; m < R.size(); ++m) {
      S += R[m].transpose() * w[m].asDiagonal() * R[m];
    }
    return S;
  }

  Eigen::MatrixXd Kernel() const {
    return WeightedGram::Compute(
        Index(R.size()), ncols,
        [this](Index m, Eigen::VectorXd& wm) { wm = w[std::size_t(m)]; },
        [this](Index m, const WeightedGram::RowsVisitor& use) {
          // a block of a larger matrix, as for the virtual rows of a slice
          Eigen::MatrixXd padded(R[std::size_t(m)].rows() + 2, ncols);
          padded.bottomRows(R[std::size_t(m)].rows()) = R[std::size_t(m)];
          use(padded.bottomRows(R[std::size_t(m)].rows()));
        });
  }
};
}  // namespace

BOOST_AUTO_TEST_CASE(weighted_gram_equals_direct_sum) {
  for (Index ncols : {1, 7, 300, 700}) {
    Blocks b(13, ncols);
    const Eigen::MatrixXd ref = b.Reference();
    const Eigen::MatrixXd S = b.Kernel();
    BOOST_CHECK_LE((S - ref).cwiseAbs().maxCoeff(),
                   1e-12 * (1.0 + ref.cwiseAbs().maxCoeff()));
    BOOST_CHECK_EQUAL((S - S.transpose()).cwiseAbs().maxCoeff(), 0.0);
  }
}

BOOST_AUTO_TEST_CASE(weighted_gram_any_batch_size) {
  Blocks b(9, 120);
  const Eigen::MatrixXd ref = b.Reference();
  // budgets below one row, about one block, and everything at once
  for (double bytes : {1.0, 8.0 * 120 * 10, 1e9}) {
    WeightedGram::setBatchMemory(bytes);
    const Eigen::MatrixXd S = b.Kernel();
    BOOST_CHECK_LE((S - ref).cwiseAbs().maxCoeff(), 1e-12);
  }
  WeightedGram::setBatchMemory(5e8);
}

BOOST_AUTO_TEST_CASE(weighted_gram_no_blocks) {
  const Eigen::MatrixXd S = WeightedGram::Compute(
      0, 4, [](Index, Eigen::VectorXd&) {},
      [](Index, const WeightedGram::RowsVisitor&) {});
  BOOST_CHECK_EQUAL(S.rows(), 4);
  BOOST_CHECK_EQUAL(S.cwiseAbs().maxCoeff(), 0.0);
}

BOOST_AUTO_TEST_CASE(accelerate_lapack_directly) {
  // printed so that a test log shows which path the helpers below took
  std::cout << "Accelerate LAPACK called directly: "
            << (accelerate::LapackAvailable() ? "yes" : "no") << std::endl;
  if (!accelerate::LapackAvailable()) {
    BOOST_CHECK_EQUAL(accelerate::dsyevd('V', 'L', 1, nullptr, 1, nullptr),
                      accelerate::kAccelerateUnavailable);
    return;
  }
  // 2 x 2 symmetric matrix with eigenvalues 1 and 3
  double a[4] = {2.0, 1.0, 1.0, 2.0};
  double w[2] = {0.0, 0.0};
  BOOST_CHECK_EQUAL(accelerate::dsyevd('V', 'L', 2, a, 2, w), 0);
  BOOST_CHECK_CLOSE(w[0], 1.0, 1e-10);
  BOOST_CHECK_CLOSE(w[1], 3.0, 1e-10);
  // Cholesky of a non-positive-definite matrix reports it
  double b[4] = {1.0, 2.0, 2.0, 1.0};
  BOOST_CHECK_GT(accelerate::dpotrf('L', 2, b, 2), 0);
}

BOOST_AUTO_TEST_CASE(symmetric_eigen_and_spd_inverse) {
  const Index n = 60;
  const Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, n);
  const Eigen::MatrixXd A = X.transpose() * X + Eigen::MatrixXd::Identity(n, n);
  const SymmetricEigenSystem es = SymmetricEigen(A);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ref(A);
  BOOST_CHECK_LE((es.values - ref.eigenvalues()).cwiseAbs().maxCoeff(), 1e-10);
  BOOST_CHECK_LE((A * es.vectors - es.vectors * es.values.asDiagonal()).norm(),
                 1e-9);
  BOOST_CHECK_LE(
      (es.vectors.transpose() * es.vectors - Eigen::MatrixXd::Identity(n, n))
          .norm(),
      1e-11);

  const Eigen::MatrixXd inv = InverseSPD(A);
  BOOST_CHECK_LE((A * inv - Eigen::MatrixXd::Identity(n, n)).norm(), 1e-10);
  BOOST_CHECK_EQUAL((inv - inv.transpose()).cwiseAbs().maxCoeff(), 0.0);

  // not positive definite: falls back to LU
  Eigen::MatrixXd B = A;
  B(0, 0) = -B(0, 0);
  const Eigen::MatrixXd invB = InverseSPD(B);
  BOOST_CHECK_LE((B * invB - Eigen::MatrixXd::Identity(n, n)).norm(), 1e-8);
}

BOOST_AUTO_TEST_CASE(eigenvalues_only_and_cholesky) {
  const Index n = 50;
  const Eigen::MatrixXd X = Eigen::MatrixXd::Random(n, n);
  const Eigen::MatrixXd A = X.transpose() * X + Eigen::MatrixXd::Identity(n, n);
  const SymmetricEigenSystem full = SymmetricEigen(A);
  const SymmetricEigenSystem values = SymmetricEigen(A, false);
  BOOST_CHECK_EQUAL(values.vectors.size(), 0);
  BOOST_CHECK_LE((full.values - values.values).cwiseAbs().maxCoeff(), 1e-10);

  const CholeskyFactor chol(A);
  BOOST_CHECK(chol.ok());
  const Eigen::MatrixXd& L = chol.matrixL();
  BOOST_CHECK_LE((L * L.transpose() - A).cwiseAbs().maxCoeff(), 1e-10);
  BOOST_CHECK_EQUAL(L.triangularView<Eigen::StrictlyUpper>()
                        .toDenseMatrix()
                        .cwiseAbs()
                        .maxCoeff(),
                    0.0);
  const Eigen::MatrixXd B = Eigen::MatrixXd::Random(n, 3);
  BOOST_CHECK_LE((A * chol.Solve(B) - B).cwiseAbs().maxCoeff(), 1e-9);
  BOOST_CHECK_LE((L * chol.SolveL(B) - B).cwiseAbs().maxCoeff(), 1e-10);

  Eigen::MatrixXd C = A;
  C(0, 0) = -1.0;
  BOOST_CHECK(!CholeskyFactor(C).ok());
}

BOOST_AUTO_TEST_SUITE_END()
