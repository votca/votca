/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE akimaspline_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/tools/akimaspline.h>

using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(akimaspline_test)

BOOST_AUTO_TEST_CASE(interpolate_test) {

  int size = 80;
  Eigen::VectorXd x = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd y = Eigen::VectorXd::Zero(size);
  for (int i = 0; i < size; ++i) {
    x(i) = 0.25 * i;
    y(i) = std::sin(x(i));
  }
  AkimaSpline cspline;
  cspline.setBCInt(0);
  cspline.Interpolate(x, y);
  Eigen::VectorXd Fref = Eigen::VectorXd::Zero(size);
  Eigen::VectorXd F2ref = Eigen::VectorXd::Zero(size);
  Fref << 0.0, 0.247404, 0.479426, 0.681639, 0.841471, 0.948985, 0.997495,
      0.983986, 0.909297, 0.778073, 0.598472, 0.381661, 0.14112, -0.108195,
      -0.350783, -0.571561, -0.756802, -0.894989, -0.97753, -0.999293,
      -0.958924, -0.858934, -0.70554, -0.508279, -0.279415, -0.0331792, 0.21512,
      0.450044, 0.656987, 0.823081, 0.938, 0.994599, 0.989358, 0.922604,
      0.798487, 0.624724, 0.412118, 0.173889, -0.0751511, -0.319519, -0.544021,
      -0.734698, -0.879696, -0.969998, -0.99999, -0.967808, -0.875452,
      -0.728665, -0.536573, -0.311119, -0.0663219, 0.182599, 0.420167, 0.631611,
      0.803784, 0.925982, 0.990607, 0.993641, 0.934895, 0.818022, 0.650288,
      0.442122, 0.206467, -0.0420244, -0.287903, -0.515882, -0.711785,
      -0.863433, -0.961397, -0.999586, -0.975626, -0.891006, -0.750987,
      -0.564276, -0.342481, -0.0993915, 0.149877, 0.389827, 0.60554, 0.783603;

  F2ref << -2.16652e-19, -0.248695, -0.481928, -0.685196, -0.845863, -0.953937,
      -1.0027, -0.989121, -0.914043, -0.782134, -0.601596, -0.383653, -0.141857,
      0.10876, 0.352614, 0.574544, 0.760752, 0.89966, 0.982632, 1.00451,
      0.963929, 0.863417, 0.709223, 0.510932, 0.280874, 0.0333524, -0.216243,
      -0.452393, -0.660415, -0.827377, -0.942896, -0.99979, -0.994522,
      -0.927419, -0.802655, -0.627984, -0.414269, -0.174797, 0.0755433,
      0.321187, 0.54686, 0.738533, 0.884287, 0.97506, 1.00521, 0.972859,
      0.880021, 0.732468, 0.539373, 0.312743, 0.066668, -0.183552, -0.42236,
      -0.634907, -0.807979, -0.930815, -0.995777, -0.998827, -0.939774,
      -0.822291, -0.653682, -0.44443, -0.207545, 0.0422437, 0.289406, 0.518574,
      0.7155, 0.86794, 0.966415, 1.0048, 0.980712, 0.895677, 0.754829, 0.567512,
      0.34318, 0.103971, -0.165813, -0.335308, -0.819762, 0;
  Eigen::VectorXd F = cspline.getSplineF();
  Eigen::VectorXd F2 = cspline.getSplineF2();

  Eigen::VectorXd Xref = Eigen::VectorXd::Zero(size);
  Xref << 0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25,
      3.5, 3.75, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7,
      7.25, 7.5, 7.75, 8, 8.25, 8.5, 8.75, 9, 9.25, 9.5, 9.75, 10, 10.25, 10.5,
      10.75, 11, 11.25, 11.5, 11.75, 12, 12.25, 12.5, 12.75, 13, 13.25, 13.5,
      13.75, 14, 14.25, 14.5, 14.75, 15, 15.25, 15.5, 15.75, 16, 16.25, 16.5,
      16.75, 17, 17.25, 17.5, 17.75, 18, 18.25, 18.5, 18.75, 19, 19.25, 19.5,
      19.75;

  bool equalX = Xref.isApprox(cspline.getX(), 1e-5);
  if (!equalX) {
    std::cout << "result X" << std::endl;
    std::cout << cspline.getX().transpose() << std::endl;
    std::cout << "ref X" << std::endl;
    std::cout << Xref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equalX, true);

  bool equal1 = Fref.isApprox(F, 1e-5);
  if (!equal1) {
    std::cout << "result F" << std::endl;
    std::cout << F.transpose() << std::endl;
    std::cout << "ref F" << std::endl;
    std::cout << Fref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal1, true);

  bool equal2 = F2ref.isApprox(F2, 1e-5);

  if (!equal2) {
    std::cout << "result F2" << std::endl;
    std::cout << F2.transpose() << std::endl;
    std::cout << "ref F2" << std::endl;
    std::cout << F2ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal2, true);

  Eigen::VectorXd rs = Eigen::VectorXd::Zero(10);
  rs << 0.45, 0.47, 0.8, 0.75, 0.6, 0.4, 0.9, 0.55;
  Eigen::VectorXd values_ref = Eigen::VectorXd::Zero(10);
  values_ref << 0.434362, 0.452449, 0.717761, 0.681639, 0.564937, 0.388734,
      0.783279, 0.52316, 0, 0;
  Eigen::VectorXd derivatives_ref = Eigen::VectorXd::Zero(10);
  derivatives_ref << 0.906562, 0.902216, 0.698382, 0.747323, 0.818045, 0.918909,
      0.615268, 0.854075, 1.02038, 1.02038;
  Eigen::VectorXd values = cspline.Calculate(rs);
  Eigen::VectorXd derivatives = cspline.CalculateDerivative(rs);

  bool equal_val = values_ref.isApprox(values, 1e-5);

  if (!equal_val) {
    std::cout << "result value" << std::endl;
    std::cout << values.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << values_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_val, true);

  bool equal_derivative = derivatives_ref.isApprox(derivatives, 1e-5);

  if (!equal_derivative) {
    std::cout << "result value" << std::endl;
    std::cout << derivatives.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << derivatives_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_derivative, true);
}

BOOST_AUTO_TEST_SUITE_END()
