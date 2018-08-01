/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE basebead_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/basebead.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/topology.h>
#include <votca/csg/molecule.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

BOOST_AUTO_TEST_SUITE(basebead_test)

BOOST_AUTO_TEST_CASE(test_basebead_constructor) {
  BaseBead basebead;
}

BOOST_AUTO_TEST_CASE(test_basebead_getters_setters) {

  BaseBead basebead; 
  BOOST_CHECK_EQUAL(_round(basebead.getMass(),3),_round(0.0,3));
  BOOST_CHECK_EQUAL(_round(basebead.getQ(),3),_round(0.0,3));
  BOOST_CHECK(!basebead.HasPos(());


	Topology top;
	string bead_type_name = "C1";
	BeadType * b_type = top.GetOrCreateBeadType(bead_type_name);

	int symmetry = 1;
	string name = "dummy";
	int resnr = 0;
	double mass = 1.21;
	double charge = -0.87;

	Bead * b = top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

	BOOST_CHECK_EQUAL(round_(b->getM(),3),round_(mass,3));
	BOOST_CHECK_EQUAL(round_(b->getQ(),3),round_(charge,3));
	BOOST_CHECK_EQUAL(b->getId(),0);
	BOOST_CHECK_EQUAL(b->getName(),name);
	BOOST_CHECK_EQUAL(b->getResnr(),resnr);
	BOOST_CHECK_EQUAL(b->getSymmetry(),symmetry);

}


BOOST_AUTO_TEST_SUITE_END()
