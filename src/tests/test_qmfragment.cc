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

#define BOOST_TEST_MODULE qmfragment_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/qmfragment.h>

using namespace std;
using namespace votca::xtp;
BOOST_AUTO_TEST_SUITE(qmfragment_test)

BOOST_AUTO_TEST_CASE(stringprocessing) {
  std::string index = "1 2 3 5...9";
  QMFragment<double> frag("check", 0, index);
  std::vector<int> ref = {1, 2, 3, 5, 6, 7, 8, 9};

  bool checked = (frag.getIndices() == ref);
  if (!checked) {
    std::cout << "result" << std::endl;
    for (int e : frag.getIndices()) {
      std::cout << e;
    }
    std::cout << std::endl;
    std::cout << "ref" << std::endl;
    for (int e : ref) {
      std::cout << e;
    }
  }
}

BOOST_AUTO_TEST_CASE(readinandwritinghdf5) {}

BOOST_AUTO_TEST_SUITE_END()
