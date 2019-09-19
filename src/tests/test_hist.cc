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

#define BOOST_TEST_MODULE hist_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <votca/xtp/hist.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(hist_test)

BOOST_AUTO_TEST_CASE(readwrite_hdf5) {

  hist<double> h;

  BOOST_CHECK_EQUAL(h.filled(), false);
  h.push_back(1.0);
  BOOST_CHECK_EQUAL(h.filled(), true);
  BOOST_CHECK_CLOSE(h.getDiff(), 1.0, 1e-9);
  h.push_back(1.5);
  BOOST_CHECK_CLOSE(h.getDiff(), 0.5, 1e-9);
  BOOST_CHECK_EQUAL(h.filled(), true);

  {
    CheckpointFile f("hist_test.hdf5");
    CheckpointWriter w = f.getWriter();
    h.WriteToCpt(w);
  }

  hist<double> h2;
  CheckpointFile f("hist_test.hdf5");
  CheckpointReader r = f.getReader();
  h2.ReadFromCpt(r);

  BOOST_CHECK_CLOSE(h.getDiff(), h2.getDiff(), 1e-9);
  BOOST_CHECK_CLOSE(h.back(), h2.back(), 1e-9);
  BOOST_CHECK_EQUAL(h.filled(), h2.filled());

  Energy_terms t;
  t.E_indu_indu() = 1.0;
  hist<Energy_terms> a;
  a.push_back(t);
  CheckpointFile ff("hist_test2.hdf5");
  CheckpointWriter ww = ff.getWriter();
  a.WriteToCpt(ww);

  hist<Energy_terms> aa;
  CheckpointReader rr = ff.getReader();
  aa.ReadFromCpt(rr);

  bool is_close = aa.getDiff().data().isApprox(a.getDiff().data(), 1e-9);
  BOOST_CHECK_EQUAL(is_close, true);
}

BOOST_AUTO_TEST_SUITE_END()
