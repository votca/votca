/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE anderson_test

// Standard includes
#include <cmath>
#include <iostream>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/anderson_mixing.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(anderson_test)

BOOST_AUTO_TEST_CASE(coeffs_test) {

  Anderson mixing_;
  mixing_.Configure(3, 0.7);
  Eigen::VectorXd in1 = Eigen::VectorXd::Zero(7);
  in1 << -0.580533, -0.535803, -0.476481, -0.380558, 0.0969526, 0.133036,
      0.164243;
  mixing_.UpdateInput(in1);

  Eigen::VectorXd out1 = Eigen::VectorXd::Zero(7);
  out1 << -0.604342, -0.548675, -0.488088, -0.385654, 0.106193, 0.139172,
      0.170433;
  mixing_.UpdateOutput(out1);

  Eigen::VectorXd mixed = mixing_.MixHistory();

  Eigen::VectorXd ref1 = Eigen::VectorXd::Zero(7);
  ref1 << -0.597199, -0.544813, -0.484606, -0.384126, 0.103421, 0.137331,
      0.168576;

  bool check_linear = mixed.isApprox(ref1, 0.00001);
  if (!check_linear) {
    std::cout << "Ref:" << ref1 << std::endl;
    std::cout << "Linear:" << mixed << std::endl;
  }

  BOOST_CHECK_EQUAL(check_linear, 1);

  mixing_.UpdateInput(mixed);

  Eigen::VectorXd out2 = Eigen::VectorXd::Zero(7);
  out2 << -0.605576, -0.549458, -0.488876, -0.385821, 0.106788, 0.139509,
      0.170768;
  mixing_.UpdateOutput(out2);

  mixed = mixing_.MixHistory();

  Eigen::VectorXd ref2 = Eigen::VectorXd::Zero(7);
  ref2 << -0.606303, -0.549862, -0.489247, -0.385968, 0.10708, 0.139698,
      0.170959;

  bool check_nonlinear_order2 = mixed.isApprox(ref2, 0.00001);
  if (!check_nonlinear_order2) {
    std::cout << "Ref:" << ref2 << std::endl;
    std::cout << "Nonlinear 2nd order:" << mixed << std::endl;
  }

  BOOST_CHECK_EQUAL(check_nonlinear_order2, 1);
  mixing_.UpdateInput(mixed);

  Eigen::VectorXd out3 = Eigen::VectorXd::Zero(7);
  out3 << -0.606242, -0.549887, -0.489296, -0.385898, 0.107162, 0.139718,
      0.170977;
  mixing_.UpdateOutput(out3);

  mixed = mixing_.MixHistory();
  Eigen::VectorXd ref3 = Eigen::VectorXd::Zero(7);
  ref3 << -0.606242, -0.549888, -0.489296, -0.385897, 0.107163, 0.139718,
      0.170977;
  bool check_nonlinear_order3 = mixed.isApprox(ref3, 0.00001);
  if (!check_nonlinear_order3) {
    std::cout << "Ref:" << ref3 << std::endl;
    std::cout << "Nonlinear 3rd order:" << mixed << std::endl;
  }

  BOOST_CHECK_EQUAL(check_nonlinear_order3, 1);
}

namespace {
// g(x) = M x + b with a symmetric M of spectral radius 0.9
struct AffineMap {
  Eigen::MatrixXd M;
  Eigen::VectorXd b;
  explicit AffineMap(votca::Index n) {
    Eigen::MatrixXd Q =
        Eigen::HouseholderQR<Eigen::MatrixXd>(
            Eigen::MatrixXd::NullaryExpr(n, n,
                                         [](votca::Index i, votca::Index j) {
                                           return std::sin(1.3 * double(i + 1) +
                                                           0.7 * double(j * j));
                                         }))
            .householderQ();
    Eigen::VectorXd lambda(n);
    for (votca::Index i = 0; i < n; ++i) {
      lambda(i) = -0.9 + 1.8 * double(i) / double(n - 1);
    }
    M = Q * lambda.asDiagonal() * Q.transpose();
    b = Eigen::VectorXd::LinSpaced(n, -1.0, 1.0);
  }
  Eigen::VectorXd operator()(const Eigen::VectorXd& x) const {
    return M * x + b;
  }
  Eigen::VectorXd FixedPoint() const {
    const votca::Index n = b.size();
    return (Eigen::MatrixXd::Identity(n, n) - M).lu().solve(b);
  }
};

// iterations until |g(x) - x| < tol (or max_iter)
votca::Index Iterate(const AffineMap& g, votca::Index order, double alpha,
                     double tol, votca::Index max_iter, Eigen::VectorXd& x) {
  Anderson mixing;
  mixing.Configure(order, alpha);
  x = Eigen::VectorXd::Zero(g.b.size());
  for (votca::Index it = 1; it <= max_iter; ++it) {
    mixing.UpdateInput(x);
    const Eigen::VectorXd gx = g(x);
    if ((gx - x).cwiseAbs().maxCoeff() < tol) {
      return it;
    }
    mixing.UpdateOutput(gx);
    x = mixing.MixHistory();
  }
  return max_iter + 1;
}
}  // namespace

// On an affine map Anderson with a history at least the dimension is a
// Krylov method: it reaches the fixed point in about n steps, where plain
// mixing needs hundreds (spectral radius 0.9).
BOOST_AUTO_TEST_CASE(affine_map_converges_in_dimension_steps) {
  const votca::Index n = 6;
  AffineMap g(n);
  Eigen::VectorXd x;
  const votca::Index it_anderson = Iterate(g, n, 1.0, 1e-10, 100, x);
  BOOST_CHECK_LE(it_anderson, n + 3);
  BOOST_CHECK_SMALL((x - g.FixedPoint()).cwiseAbs().maxCoeff(), 1e-9);

  const votca::Index it_linear = Iterate(g, 0, 0.5, 1e-10, 1000, x);
  BOOST_CHECK_GT(it_linear, 10 * it_anderson);

  // a short history still converges, more slowly
  const votca::Index it_short = Iterate(g, 2, 0.7, 1e-10, 500, x);
  BOOST_CHECK_LE(it_short, 500);
  BOOST_CHECK_SMALL((x - g.FixedPoint()).cwiseAbs().maxCoeff(), 1e-9);
}

// Only the last order + 1 input/output pairs are used.
BOOST_AUTO_TEST_CASE(history_is_bounded) {
  const votca::Index n = 5;
  const votca::Index order = 2;
  AffineMap g(n);
  std::vector<Eigen::VectorXd> in;
  std::vector<Eigen::VectorXd> out;
  Anderson long_run;
  long_run.Configure(order, 0.7);
  Eigen::VectorXd x = Eigen::VectorXd::Constant(n, 0.3);
  Eigen::VectorXd mixed;
  for (votca::Index it = 0; it < 8; ++it) {
    in.push_back(x);
    out.push_back(g(x));
    long_run.UpdateInput(in.back());
    long_run.UpdateOutput(out.back());
    mixed = long_run.MixHistory();
    x = mixed;
  }
  Anderson fresh;
  fresh.Configure(order, 0.7);
  for (std::size_t k = in.size() - std::size_t(order + 1); k < in.size(); ++k) {
    fresh.UpdateInput(in[k]);
    fresh.UpdateOutput(out[k]);
  }
  BOOST_CHECK_SMALL((fresh.MixHistory() - mixed).cwiseAbs().maxCoeff(), 1e-14);
}

// A stagnating history (the same input/output pair twice, as when evGW
// stalls) makes the coefficient problem singular: the mix must stay finite
// and must not increase the residual of the latest pair.
BOOST_AUTO_TEST_CASE(repeated_history_stays_finite) {
  const votca::Index n = 4;
  AffineMap g(n);
  Anderson mixing;
  mixing.Configure(5, 1.0);
  const Eigen::VectorXd x0 = Eigen::VectorXd::Constant(n, 0.1);
  const Eigen::VectorXd x1 = Eigen::VectorXd::Constant(n, -0.2);
  for (const Eigen::VectorXd* x : {&x0, &x1, &x1, &x1}) {
    mixing.UpdateInput(*x);
    mixing.UpdateOutput(g(*x));
  }
  const Eigen::VectorXd mixed = mixing.MixHistory();
  BOOST_CHECK(mixed.allFinite());
  // with alpha = 1 on an affine map, the mixed output is g of the mixed
  // input; its residual is the least-squares one, at most the latest
  BOOST_CHECK_LE(
      (g(mixed) - mixed).norm(),
      (g(x1) - x1).norm() * (1.0 + 1e-12) + (g.M.norm() + 1.0) * 1e-12);
}

BOOST_AUTO_TEST_SUITE_END()
