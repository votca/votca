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

#define BOOST_TEST_MODULE bse_fullsolver_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/bse_fullsolver.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(bse_fullsolver_test)

namespace {
// dense stand-in for a matrix-free BSE operator
struct DenseOp {
  Eigen::MatrixXd M;
  Index rows() const { return M.rows(); }
  Eigen::VectorXd diagonal() const { return M.diagonal(); }
  Eigen::MatrixXd matmul(const Eigen::MatrixXd& x) const { return M * x; }
};

// A, B shaped like a BSE: diagonal pair energies, weak couplings
std::pair<DenseOp, DenseOp> MakeAB(Index n, unsigned seed) {
  std::srand(seed);
  Eigen::MatrixXd K = 0.02 * Eigen::MatrixXd::Random(n, n);
  K = 0.5 * (K + K.transpose()).eval();
  Eigen::MatrixXd Bm = 0.015 * Eigen::MatrixXd::Random(n, n);
  Bm = 0.5 * (Bm + Bm.transpose()).eval();
  Eigen::MatrixXd Am = K;
  for (Index i = 0; i < n; ++i) {
    Am(i, i) += 0.3 + 0.8 * double(i) / double(n);
  }
  return {DenseOp{Am}, DenseOp{Bm}};
}

Eigen::VectorXd DenseRoots(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> amb(A - B);
  const Eigen::MatrixXd root = amb.operatorSqrt();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(root * (A + B) * root);
  return es.eigenvalues().cwiseSqrt();
}
}  // namespace

BOOST_AUTO_TEST_CASE(roots_and_vectors_match_dense) {
  const Index n = 300;
  const Index nroots = 8;
  auto [A, B] = MakeAB(n, 7);
  Logger log;
  FullBSEDavidson::Options opt;
  opt.tolerance = 1e-7;
  opt.max_iterations = 100;
  opt.max_subspace = 4 * nroots;  // forces restarts
  FullBSEDavidson solver(log, opt);
  votca::tools::EigenSystem es = solver.Solve(A, B, nroots);

  const Eigen::VectorXd ref = DenseRoots(A.M, B.M).head(nroots);
  BOOST_CHECK_SMALL((es.eigenvalues() - ref).cwiseAbs().maxCoeff(), 1e-9);
  const Eigen::MatrixXd& X = es.eigenvectors();
  const Eigen::MatrixXd& Y = es.eigenvectors2();
  for (Index i = 0; i < nroots; ++i) {
    const double w = es.eigenvalues()(i);
    // A X + B Y = w X, B X + A Y = -w Y
    BOOST_CHECK_SMALL((A.M * X.col(i) + B.M * Y.col(i) - w * X.col(i)).norm(),
                      1e-6);
    BOOST_CHECK_SMALL((B.M * X.col(i) + A.M * Y.col(i) + w * Y.col(i)).norm(),
                      1e-6);
    BOOST_CHECK_SMALL(X.col(i).squaredNorm() - Y.col(i).squaredNorm() - 1.0,
                      1e-10);
  }
}

BOOST_AUTO_TEST_CASE(unstable_reference_throws) {
  const Index n = 40;
  auto [A, B] = MakeAB(n, 3);
  B.M = A.M;  // A - B = 0: not positive definite
  Logger log;
  FullBSEDavidson solver(log, FullBSEDavidson::Options());
  BOOST_CHECK_THROW(solver.Solve(A, B, 4), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
