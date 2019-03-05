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

#define BOOST_TEST_MODULE topology_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <votca/xtp/atom.h>
#include <votca/xtp/qmpair.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/topology.h>

#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

#include <votca/csg/boundarycondition.h>

using namespace votca::tools;

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(topology_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Topology top; }

BOOST_AUTO_TEST_CASE(box_test) {

  Eigen::Matrix3d box = 2 * Eigen::Matrix3d::Identity();
  Topology top;
  top.setBox(box);

  auto vol = top.BoxVolume();
  BOOST_CHECK_CLOSE(vol, 8, 0.0001);
  auto box2 = top.getBox();

  BOOST_CHECK_EQUAL(box.isApprox(box2, 1e-6), true);
}

BOOST_AUTO_TEST_CASE(simple_test) {

  Topology top;
  top.setStep(1);
  BOOST_CHECK_EQUAL(top.getStep(), 1);
  top.setTime(1.21);
  BOOST_CHECK_CLOSE(top.getTime(), 1.21, 0.0001);
}

BOOST_AUTO_TEST_SUITE_END()
