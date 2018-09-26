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

#define BOOST_TEST_MODULE elements_test
#include <boost/test/unit_test.hpp>
#include <exception>
#include <cmath>
#include <votca/tools/elements.h>

using namespace std;
using namespace votca::tools;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

BOOST_AUTO_TEST_SUITE(elements_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Elements ele; }

BOOST_AUTO_TEST_CASE(accessors_test) {
  Elements ele;
  BOOST_CHECK_EQUAL(ele.getVdWChelpG("H"), 1.45);
  BOOST_CHECK_THROW(ele.getVdWChelpG("Blah"), invalid_argument);

  BOOST_CHECK_EQUAL(ele.getMass("K"), 39.098);
  BOOST_CHECK_EQUAL(ele.getEleNum("Li"), 3);

  BOOST_CHECK_EQUAL(ele.getEleName(17), "Cl");
  BOOST_CHECK_EQUAL(ele.getEleShort("MAGNESIUM"), "Mg");
  BOOST_CHECK_EQUAL(ele.getEleFull("Ge"), "GERMANIUM");
  BOOST_CHECK_EQUAL(ele.getVdWMK("F"), 1.35);
  BOOST_CHECK_THROW(ele.getVdWMK("Pb"), invalid_argument);
  BOOST_CHECK_EQUAL(round_(ele.getCovRad("Cl","ang"),3),1.02);
  BOOST_CHECK_EQUAL(round_(ele.getCovRad("Cl","nm"),3),0.102);
  BOOST_CHECK_THROW(round_(ele.getCovRad("Cl","Blah"),3),invalid_argument);
  
  BOOST_CHECK_EQUAL(ele.getPolarizability("F"), 0.440e-3);
  BOOST_CHECK_THROW(ele.getPolarizability("Pb"), invalid_argument);

	BOOST_CHECK(ele.isMassAssociatedWithElement(12.01,0.01));
	BOOST_CHECK(!ele.isMassAssociatedWithElement(12.51,0.01));
	BOOST_CHECK_EQUAL("C",ele.getEleShortClosestInMass(12.01,0.01));
  
}

BOOST_AUTO_TEST_SUITE_END()
