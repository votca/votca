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

#define BOOST_TEST_MODULE segment_test
#include <boost/test/unit_test.hpp>
#include <vector>
#include <votca/xtp/checkpoint.h>
#include <votca/xtp/segment.h>

using namespace votca::tools;
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(segment_test)

BOOST_AUTO_TEST_CASE(getElementtest) {
  Segment seg("one", 1);
  Atom atm1(1, "C", Eigen::Vector3d::Zero());
  Atom atm2(2, "H", Eigen::Vector3d::UnitX());
  Atom atm3(3, "N", Eigen::Vector3d::UnitY());
  Atom atm4(4, "Au", Eigen::Vector3d::UnitZ());
  Atom atm5(5, "C", 2 * Eigen::Vector3d::UnitZ());

  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm3);
  seg.push_back(atm5);
  std::vector<std::string> unique_ele = seg.FindUniqueElements();

  std::vector<std::string> unique_ele_ref = {"C", "H", "N", "Au"};

  BOOST_CHECK_EQUAL(unique_ele.size(), unique_ele_ref.size());

  for (Index i = 0; i < unique_ele.size(); i++) {
    BOOST_CHECK_EQUAL(unique_ele[i], unique_ele_ref[i]);
  }
}

BOOST_AUTO_TEST_CASE(writehdf5) {
  Segment seg("one", 1);
  Atom atm1(1, "C", Eigen::Vector3d::Zero());
  Atom atm2(2, "H", Eigen::Vector3d::UnitX());
  Atom atm3(3, "N", Eigen::Vector3d::UnitY());
  Atom atm4(4, "Au", Eigen::Vector3d::UnitZ());
  Atom atm5(5, "C", 2 * Eigen::Vector3d::UnitZ());

  seg.push_back(atm1);
  seg.push_back(atm2);
  seg.push_back(atm2);
  seg.push_back(atm3);
  seg.push_back(atm4);
  seg.push_back(atm3);
  seg.push_back(atm5);
  {
    CheckpointFile f("segment_test.hdf5");
    CheckpointWriter w = f.getWriter();
    seg.WriteToCpt(w);
  }

  CheckpointFile f("segment_test.hdf5");
  CheckpointReader r = f.getReader();
  Segment seg2(r);
  for (Index i = 0; i < seg2.size(); i++) {
    const Atom& a1 = seg[i];
    const Atom& a2 = seg2[i];
    BOOST_CHECK_EQUAL(a1.getName(), a2.getName());
    bool pos_equal = a1.getPos().isApprox(a2.getPos(), 1e-9);
    BOOST_CHECK_EQUAL(pos_equal, true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
