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
  QMFragment<double> frag(0, index);
  std::vector<long> ref = {1, 2, 3, 5, 6, 7, 8, 9};

  bool checked = (frag.getIndices() == ref);
  if (!checked) {
    std::cout << "result" << std::endl;
    for (long e : frag.getIndices()) {
      std::cout << e;
    }
    std::cout << std::endl;
    std::cout << "ref" << std::endl;
    for (long e : ref) {
      std::cout << e;
    }
  }
}

BOOST_AUTO_TEST_CASE(readinandwritinghdf5) {
  std::string index = "1 2 3 5...9";
  QMFragment<double> frag(0, index);

  frag.value() = -0.5;

  {
    CheckpointFile f("QMFragment_test.hdf5");
    CheckpointWriter w = f.getWriter();
    frag.WriteToCpt(w);
  }

  CheckpointFile f("QMFragment_test.hdf5");
  CheckpointReader r = f.getReader();
  QMFragment<double> frag2(r);
  BOOST_CHECK_EQUAL(frag2.getId(), frag.getId());
  BOOST_CHECK_CLOSE(frag.value(), frag2.value(), 1e-5);
  for (long i = 0; i < frag.size(); i++) {
    BOOST_CHECK_EQUAL(frag.getIndices()[i], frag2.getIndices()[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END()
