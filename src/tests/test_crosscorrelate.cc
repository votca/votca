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

#define BOOST_TEST_MODULE crosscorrelate_test
#include "../../include/votca/tools/crosscorrelate.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <iostream>
using namespace votca;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(crosscorrelate_test)

BOOST_AUTO_TEST_CASE(crosscorrelate) {

  DataCollection<double> d;
  DataCollection<double>::array* x = d.CreateArray("first");
  x->resize(50);
  for (Index i = 0; i < Index(x->size()); i++) {
    (*x)[i] = std::sin(double(i) * 0.5) + std::cos(double(i) * 1.5);
  }

  DataCollection<double>::selection s;
  s.push_back(x);

  CrossCorrelate cor;
  cor.AutoCorrelate(s);
  std::vector<double>& corr = cor.getData();
  Eigen::Map<Eigen::VectorXd> m(corr.data(), corr.size());
  Eigen::VectorXd corr_ref = Eigen::VectorXd::Zero(50);
  corr_ref << 1, 0.467627, -0.228828, -0.0618146, 0.270531, -0.24973, -0.956717,
      -0.674499, 0.117233, 0.172174, -0.246548, 0.0477983, 0.84253, 0.842251,
      0.0484615, -0.245153, 0.170656, 0.11432, -0.672333, -0.951362, -0.250719,
      0.264426, -0.0608222, -0.22103, 0.466958, 0.989178, 0.466958, -0.22103,
      -0.0608222, 0.264426, -0.250719, -0.951362, -0.672333, 0.11432, 0.170656,
      -0.245153, 0.0484615, 0.842251, 0.84253, 0.0477983, -0.246548, 0.172174,
      0.117233, -0.674499, -0.956717, -0.24973, 0.270531, -0.0618146, -0.228828,
      0.467627;
  bool equal_val = corr_ref.isApprox(m, 1e-5);

  if (!equal_val) {
    std::cout << "result value" << std::endl;
    std::cout << m.transpose() << std::endl;
    std::cout << "ref value" << std::endl;
    std::cout << corr_ref.transpose() << std::endl;
  }
  BOOST_CHECK_EQUAL(equal_val, true);
}

BOOST_AUTO_TEST_SUITE_END()
