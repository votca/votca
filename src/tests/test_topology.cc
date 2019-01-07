/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <votca/xtp/qmpair.h>
#include <votca/xtp/topology.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/atom.h>

#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>

#include <votca/csg/boundarycondition.h>

using namespace votca::tools;

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(topology_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Topology top; }

BOOST_AUTO_TEST_CASE(box_test) { 

  // Box takes a vector
  double x1 = 2.0;
  double y1 = 0.0;
  double z1 = 0.0;

  double x2 = 0.0;
  double y2 = 2.0;
  double z2 = 0.0;

  double x3 = 0.0;
  double y3 = 0.0;
  double z3 = 2.0;

  vec v1(x1,y1,z1);
  vec v2(x2,y2,z2);
  vec v3(x3,y3,z3);

  matrix box(v1,v2,v3);

  Topology top;
  top.setBox(box);

  auto vol = top.BoxVolume();
  BOOST_CHECK_CLOSE( vol, 8, 0.0001 );
  auto box2 = top.getBox();
  
  auto v1_2 = box2.getCol(0);
  auto v2_2 = box2.getCol(1);
  auto v3_2 = box2.getCol(2);
  BOOST_CHECK_EQUAL(v1,v1_2);
  BOOST_CHECK_EQUAL(v2,v2_2);
  BOOST_CHECK_EQUAL(v3,v3_2);
}

BOOST_AUTO_TEST_CASE(simple_test){

  Topology top;
  top.setStep(1);
  BOOST_CHECK_EQUAL(top.getStep(),1);
  top.setTime(1.21);
  BOOST_CHECK_EQUAL(static_cast<int>(top.getTime()*100),121);
  top.setDatabaseId(3);
  BOOST_CHECK_EQUAL(top.getDatabaseId(),3);

}

BOOST_AUTO_TEST_SUITE_END()
