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

#define BOOST_TEST_MODULE gaussian_quadratures_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/quadrature_factory.h"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace votca;

// defines a Gaussian as integration test
// result should by ~sqrt(pi)
class FunctionEvaluation {
 public:
  FunctionEvaluation(){};

  double operator()(Index j, double point, bool symmetry) const {
    double factor = 0.0;
    // this is only here to staisfy the unused variable warning
    if (j < 100) {
      if (symmetry) {
        factor = 2.0;
      } else {
        factor = 1.0;
      }
    }

    return factor * exp(-std::pow(point, 2));
  }
};

BOOST_AUTO_TEST_SUITE(gaussian_quadratures_test)

BOOST_AUTO_TEST_CASE(gauss_legendre) {

  QuadratureFactory::RegisterAll();
  std::unique_ptr<GaussianQuadratureBase> _gq =
      std::unique_ptr<GaussianQuadratureBase>(Quadratures().Create("legendre"));

  std::vector<int> orders{8, 10, 12, 14, 16, 18, 20, 40, 100};
  FunctionEvaluation f = FunctionEvaluation();

  Eigen::VectorXd integrals(9);
  for (Index i = 0; i < 9; i++) {
    _gq->configure(orders[i]);
    integrals(i) = _gq->Integrate(f);
  }

  Eigen::VectorXd integrals_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/gaussian_quadratures/gauss_legendre.mm");

  bool check_integral = integrals.isApprox(integrals_ref, 1e-10);
  if (!check_integral) {
    std::cout << "Gauss-Legendre" << std::endl;
    std::cout << integrals << std::endl;
    std::cout << "Gauss-Legendre ref" << std::endl;
    std::cout << integrals_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(check_integral, true);
}

BOOST_AUTO_TEST_CASE(gauss_laguerre) {

  QuadratureFactory::RegisterAll();
  std::unique_ptr<GaussianQuadratureBase> _gq =
      std::unique_ptr<GaussianQuadratureBase>(Quadratures().Create("laguerre"));
  std::vector<int> orders{8, 10, 12, 14, 16, 18, 20, 40, 100};
  FunctionEvaluation f = FunctionEvaluation();

  Eigen::VectorXd integrals(9);
  for (Index i = 0; i < 9; i++) {
    _gq->configure(orders[i]);
    integrals(i) = _gq->Integrate(f);
  }

  Eigen::VectorXd integrals_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/gaussian_quadratures/gauss_laguerre.mm");

  bool check_integral = integrals.isApprox(integrals_ref, 1e-10);
  if (!check_integral) {
    std::cout << "Gauss-Laguerre" << std::endl;
    std::cout << integrals << std::endl;
    std::cout << "Gauss-Laguerre ref" << std::endl;
    std::cout << integrals_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(check_integral, true);
}

BOOST_AUTO_TEST_CASE(gauss_hermite) {

  QuadratureFactory::RegisterAll();
  std::unique_ptr<GaussianQuadratureBase> _gq =
      std::unique_ptr<GaussianQuadratureBase>(Quadratures().Create("hermite"));
  std::vector<int> orders{8, 10, 12, 14, 16, 18, 20, 40, 100};
  FunctionEvaluation f = FunctionEvaluation();

  Eigen::VectorXd integrals(9);
  for (Index i = 0; i < 9; i++) {
    _gq->configure(orders[i]);
    integrals(i) = _gq->Integrate(f);
  }

  Eigen::VectorXd integrals_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/gaussian_quadratures/gauss_hermite.mm");

  bool check_integral = integrals.isApprox(integrals_ref, 1e-10);
  if (!check_integral) {
    std::cout << "Gauss-Hermite" << std::endl;
    std::cout << integrals << std::endl;
    std::cout << "Gauss-Hermite ref" << std::endl;
    std::cout << integrals_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(check_integral, true);
}

BOOST_AUTO_TEST_SUITE_END()
