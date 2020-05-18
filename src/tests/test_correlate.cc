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

#define BOOST_TEST_MODULE correlate_test
#include "../../include/votca/tools/correlate.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <iostream>
using namespace votca;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(correlate_test)

BOOST_AUTO_TEST_CASE(correlation_test) {

  DataCollection<double> d;
  DataCollection<double>::array* x = d.CreateArray("first");
  x->resize(50);
  for (Index i = 0; i < Index(x->size()); i++) {
    (*x)[i] = double(i) * 0.5;
  }

  DataCollection<double>::array* y = d.CreateArray("second");
  *y = *x;

  DataCollection<double>::array* z = d.CreateArray("third");
  *z = *x;
  std::reverse(z->begin(), z->end());
  DataCollection<double>::selection s;
  s.push_back(x);
  s.push_back(y);
  s.push_back(z);
  Correlate cor;
  cor.CalcCorrelations(s);
  BOOST_CHECK_CLOSE(cor.getData()[0], 1, 1e-9);
  BOOST_CHECK_CLOSE(cor.getData()[1], -1, 1e-9);
}

BOOST_AUTO_TEST_SUITE_END()
