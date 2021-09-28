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

#define BOOST_TEST_MODULE eigen_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/eigen.h"

using namespace votca;

BOOST_AUTO_TEST_SUITE(eigen_test)

BOOST_AUTO_TEST_CASE(symmetric_test) {

  Eigen::MatrixXd notsym = Eigen::MatrixXd::Random(6, 6);

  Eigen::MatrixXd sym = notsym + notsym.transpose();

  bool symmetry = true;
  for (votca::Index i = 0; i < sym.rows(); ++i) {
    for (votca::Index j = 0; j <= i; ++j) {
      if (std::abs(sym(i, j) - sym(j, i)) > 1e-9) {
        symmetry = false;
        break;
      }
    }
  }

  BOOST_CHECK_EQUAL(symmetry, 1);
}

BOOST_AUTO_TEST_CASE(inverse_test) {
  Eigen::MatrixXd m = Eigen::MatrixXd::Zero(17, 17);
  m << -0.0107438, 0.177662, -0.000893645, -0.000893645, -0.000893645, 0.345174,
      -0.000591276, -0.000591276, -0.000591276, 0.136613, 0.110948, 0.140527,
      0.115047, 0.140527, 0.115047, 0.140527, 0.115047, 2.183e-15, 1.39936e-14,
      -0.2034, 0.00238481, 0.201015, -7.8295e-14, -0.174952, 0.00205127,
      0.172901, -5.75953e-13, -7.55155e-13, -0.234829, -0.324284, 0.232076,
      0.320482, 0.00275331, 0.00380215, 1.92294e-15, 5.30409e-14, -0.114679,
      0.233489, -0.11881, -1.08691e-13, -0.0986399, 0.200833, -0.102193,
      1.10298e-12, 1.56415e-12, -0.132399, -0.182835, -0.137168, -0.189421,
      0.269568, 0.372256, 9.81255e-05, -0.00180964, -0.165679, -0.165679,
      -0.165679, -0.00572861, -0.144589, -0.144589, -0.144589, -0.283887,
      -0.391777, 0.0942211, 0.13068, 0.0942211, 0.13068, 0.0942211, 0.13068,
      0.0257368, -0.0923997, 0.00531379, 0.00531379, 0.00531379, -2.45842,
      0.0248521, 0.0248521, 0.0248521, -0.0235266, 0.870391, -0.021766,
      0.938627, -0.021766, 0.938627, -0.021766, 0.938627, -6.05072e-15,
      -3.50553e-14, 0.0320188, -0.217226, 0.185207, 2.84606e-13, 0.164774,
      -1.11788, 0.953107, -3.43579e-13, -7.51681e-12, -0.0107087, -0.224765,
      -0.0619426, -1.30011, 0.0726513, 1.52488, -1.85212e-14, 3.71612e-14,
      0.232345, -0.0884434, -0.143902, 1.65531e-12, 1.19569, -0.455144,
      -0.740541, -5.57708e-13, -1.28855e-11, -0.0777079, -1.63101, 0.0481279,
      1.01016, 0.0295799, 0.620853, -0.00067642, 0.0030588, 0.164587, 0.164587,
      0.164587, 0.0706733, 0.85369, 0.85369, 0.85369, -0.0833538, -1.77428,
      0.0291688, 0.555286, 0.0291688, 0.555286, 0.0291688, 0.555286,
      2.67521e-15, 1.49908e-14, -0.499711, 0.0302687, 0.469443, -8.25069e-14,
      0.794429, -0.0481205, -0.746308, -1.68595e-12, 3.88068e-13, -0.631236,
      0.170551, 0.593, -0.16022, 0.0382355, -0.0103307, -2.92301e-15,
      -4.20358e-14, 0.253557, -0.559541, 0.305984, 2.74225e-14, -0.403099,
      0.889545, -0.486446, -4.63627e-12, 1.28551e-12, 0.320293, -0.0865389,
      0.386519, -0.104432, -0.706813, 0.190971, 0.000103179, -0.000739357,
      -0.399322, -0.399322, -0.399322, -0.00417793, 0.631921, 0.631921,
      0.631921, -0.743087, 0.196419, 0.252682, -0.0664109, 0.252682, -0.0664109,
      0.252682, -0.0664109, 0.00149008, -0.210613, 0.0151116, 0.0151116,
      0.0151116, -0.131439, -0.0188722, -0.0188722, -0.0188722, 0.641757,
      -0.258506, 0.662618, -0.303913, 0.662618, -0.303913, 0.662618, -0.303913,
      -1.01256e-14, 8.31141e-14, 0.154221, -0.726494, 0.572273, 3.21687e-13,
      -0.180079, 0.848303, -0.668224, -6.90137e-12, 1.08126e-11, -0.16221,
      0.269927, -0.601918, 1.00163, 0.764127, -1.27156, -1.4385e-14,
      4.61332e-14, 0.749843, -0.241362, -0.508481, 3.69788e-13, -0.875567,
      0.281831, 0.593736, -1.00864e-11, 1.59982e-11, -0.788686, 1.31242,
      0.534821, -0.889976, 0.253865, -0.422448, 8.45112e-05, -0.00427837,
      -0.539015, -0.539015, -0.539015, -0.00364986, 0.62842, 0.62842, 0.62842,
      0.873707, -1.42717, -0.270316, 0.466699, -0.270316, 0.466699, -0.270316,
      0.466699, 0.220745, 1.96298, 0.000217492, 0.000217492, 0.000217492,
      -3.94464, -0.000657288, -0.000657288, -0.000657288, 0.170958, 0.775841,
      0.170979, 0.77469, 0.170979, 0.77469, 0.170979, 0.77469, -1.00858,
      0.364513, 2.52534e-06, 2.52534e-06, 2.52534e-06, -0.192712, -5.74555e-07,
      -5.74555e-07, -5.74555e-07, -0.0243332, 0.0351697, -0.0243273, 0.0351687,
      -0.0243273, 0.0351687, -0.0243273, 0.0351687;

  Eigen::FullPivLU<Eigen::MatrixXd> dec(m);

  bool check_inv =
      (m * dec.inverse()).isApprox(Eigen::MatrixXd::Identity(17, 17), 0.01);

  if (!check_inv) {
    std::cout << dec.isInvertible() << " " << dec.rank() << " " << dec.rcond()
              << std::endl;
    std::cout << m * dec.solve(Eigen::MatrixXd::Identity(17, 17)) << std::endl;
  }
  BOOST_CHECK_EQUAL(check_inv, 1);
}

BOOST_AUTO_TEST_CASE(PCG_test) {

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3, 3);
  A << 1, 3, -2, 3, 5, 6, 2, 4, 3;
  Eigen::MatrixXd A_sym = A * A.transpose();
  Eigen::VectorXd x_ref = Eigen::VectorXd::Zero(3);
  x_ref << -15, 8, 2;
  Eigen::VectorXd b = A_sym * x_ref;

  Eigen::LLT<Eigen::MatrixXd> lltOfA(A_sym);
  Eigen::VectorXd x_chol = lltOfA.solve(b);
  bool check_chol = x_chol.isApprox(x_ref, 1e-7);
  BOOST_CHECK_EQUAL(check_chol, 1);
  if (!check_chol) {
    std::cout << "result" << std::endl;
    std::cout << x_chol.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << x_ref.transpose() << std::endl;
  }

  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg1;
  cg1.compute(A_sym);
  Eigen::VectorXd x_cg1 = cg1.solve(b);

  bool check_cg1 = x_cg1.isApprox(x_ref, 1e-7);
  BOOST_CHECK_EQUAL(check_cg1, 1);
  if (!check_cg1) {
    std::cout << "result" << std::endl;
    std::cout << x_cg1.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << x_ref.transpose() << std::endl;
  }

  Eigen::ConjugateGradient<Eigen::MatrixXd, Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      cg2;
  cg2.compute(A_sym);
  Eigen::VectorXd x_cg2 = cg2.solveWithGuess(b, x_ref);

  bool check_cg2 = x_cg2.isApprox(x_ref, 1e-7);
  BOOST_CHECK_EQUAL(check_cg2, 1);
  if (!check_cg2) {
    std::cout << "result" << std::endl;
    std::cout << x_cg2.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << x_ref.transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
