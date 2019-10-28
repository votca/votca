/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE bfgs_test
#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/adiis_costfunction.h>
#include <votca/xtp/bfgs-trm.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/optimiser_costfunction.h>
using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(bfgs_test)

BOOST_AUTO_TEST_CASE(parabola_test) {
  class parabola : public Optimiser_costfunction {

    double EvaluateCost(const Eigen::VectorXd& parameters) override {
      Eigen::VectorXd value = parameters;
      value(0) -= 2;
      double cost = value.cwiseAbs2().sum();
      return cost;
    }

    Eigen::VectorXd EvaluateGradient(
        const Eigen::VectorXd& parameters) override {
      Eigen::VectorXd gradient = 2 * parameters;
      gradient(0) -= 4;
      return gradient;
    }

    bool Converged(const Eigen::VectorXd&, double,
                   const Eigen::VectorXd& gradient) override {
      if (gradient.cwiseAbs().maxCoeff() < 1e-8) {
        return true;
      } else {
        return false;
      }
    }

    Index NumParameters() const override { return 5; }
  };

  parabola p5;
  BFGSTRM bfgstrm(p5);
  bfgstrm.setNumofIterations(100);
  bfgstrm.setTrustRadius(0.1);
  bfgstrm.Optimize(5 * Eigen::VectorXd::Ones(5));

  Eigen::VectorXd ref = Eigen::VectorXd::Zero(5);
  ref(0) = 2;
  bool equal = bfgstrm.getParameters().isApprox(ref, 0.00001);
  if (!equal) {
    cout << "minimum found:" << endl;
    cout << bfgstrm.getParameters() << endl;
    cout << "minimum ref:" << endl;
    cout << ref << endl;
  } else {
    cout << bfgstrm.getIteration() << endl;
  }
  BOOST_CHECK_EQUAL(equal, 1);
}

BOOST_AUTO_TEST_CASE(booth_test) {
  class booth : public Optimiser_costfunction {

    double EvaluateCost(const Eigen::VectorXd& parameters) override {
      double x = parameters[0];
      double y = parameters[1];

      return (x + 2 * y - 7) * (x + 2 * y - 7) +
             (2 * x + y - 5) * (2 * x + y - 5);
    }

    Eigen::VectorXd EvaluateGradient(
        const Eigen::VectorXd& parameters) override {
      double x = parameters[0];
      double y = parameters[1];
      Eigen::VectorXd gradient = Eigen::VectorXd::Zero(2);
      gradient[0] = 2 * (5 * x + 4 * y - 17);
      gradient[1] = 2 * (4 * x + 5 * y - 19);
      return gradient;
    }

    bool Converged(const Eigen::VectorXd&, double,
                   const Eigen::VectorXd& gradient) override {
      if (gradient.cwiseAbs().maxCoeff() < 1e-8) {
        return true;
      } else {
        return false;
      }
    }

    Index NumParameters() const override { return 2; }
  };

  booth p2;
  BFGSTRM bfgstrm(p2);
  bfgstrm.setNumofIterations(1000);
  bfgstrm.setTrustRadius(0.1);
  bfgstrm.Optimize(5 * Eigen::VectorXd::Ones(2));

  Eigen::VectorXd ref = Eigen::VectorXd::Zero(2);
  ref << 1, 3;
  bool equal = bfgstrm.getParameters().isApprox(ref, 0.00001);
  if (!equal) {
    cout << "minimum found:" << endl;
    cout << bfgstrm.getParameters() << endl;
    cout << "minimum ref:" << endl;
    cout << ref << endl;
  } else {
    cout << bfgstrm.getIteration() << endl;
  }
  BOOST_CHECK_EQUAL(equal, 1);
}

BOOST_AUTO_TEST_CASE(adiis_test) {
  Index size = 5;

  Eigen::VectorXd DiF = Eigen::VectorXd::Zero(size);
  DiF << 0.679243, 0.562675, 0.39399, -0.0258519, 0;
  Eigen::MatrixXd DiFj = Eigen::MatrixXd::Zero(size, size);
  DiFj << 0.613998, 0.192684, 0.0326914, -0.193661, 0, 0.192371, 0.0653754,
      0.0131825, -0.0631915, 0, 0.032739, 0.0131883, 0.0038493, -0.0110991, 0,
      -0.192873, -0.0631203, -0.0111063, 0.0633221, 0, 0, 0, 0, 0, 0;

  Logger log;

  ADIIS_costfunction a_cost = ADIIS_costfunction(DiF, DiFj);
  BFGSTRM optimizer = BFGSTRM(a_cost);
  optimizer.setNumofIterations(1000);
  optimizer.setTrustRadius(0.01);
  optimizer.setLog(&log);
  // Starting point: equal weights on all matrices
  Eigen::VectorXd coeffs = Eigen::VectorXd::Constant(size, 1.0 / size);
  optimizer.Optimize(coeffs);
  bool success = optimizer.Success();
  coeffs = optimizer.getParameters().cwiseAbs2();
  double xnorm = coeffs.sum();
  coeffs /= xnorm;

  BOOST_CHECK_EQUAL(success, true);
  Eigen::VectorXd ref = Eigen::VectorXd::Zero(size);
  ref << 0, 0, 0, 0.40826075912352, 0.59173924087648;
  bool equal = coeffs.isApprox(ref, 0.00001);
  if (!equal) {
    cout << "minimum found:" << endl;
    cout << coeffs << endl;
    cout << "minimum ref:" << endl;
    cout << ref << endl;
  } else {
    cout << optimizer.getIteration() << endl;
  }
  BOOST_CHECK_EQUAL(equal, 1);
}

BOOST_AUTO_TEST_SUITE_END()
