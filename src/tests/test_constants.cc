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

#define BOOST_TEST_MODULE constants_test
#include <boost/test/unit_test.hpp>
#include <votca/tools/constants.h>

using namespace votca::tools;
using namespace votca::tools::conv;
using namespace votca::tools::topology_constants;

BOOST_AUTO_TEST_SUITE(constants_test)

BOOST_AUTO_TEST_CASE(constants_test1) {
  BOOST_CHECK_CLOSE(Pi, boost::math::constants::pi<double>(), 5E-5);
  BOOST_CHECK_CLOSE(kB, 8.617332478E-5, 5E-5);
  BOOST_CHECK_CLOSE(hbar, 6.5821192815E-16, 5E-5);
  BOOST_CHECK_CLOSE(bohr2nm, 0.052917721092, 5E-5);
  BOOST_CHECK_CLOSE(nm2bohr, 18.897259886, 5E-5);
  BOOST_CHECK_CLOSE(ang2bohr, 1.8897259886, 5E-5);
  BOOST_CHECK_CLOSE(bohr2ang, 1.0 / 1.8897259886, 5E-5);
  BOOST_CHECK_CLOSE(nm2ang, 10.0, 5E-5);
  BOOST_CHECK_CLOSE(ang2nm, 0.1, 5E-5);
  BOOST_CHECK_CLOSE(hrt2ev, 27.21138602, 5E-5);
  BOOST_CHECK_CLOSE(ev2hrt, 1.0 / 27.21138602, 5E-5);
  BOOST_CHECK_CLOSE(ev2kj_per_mol, 96.485, 5E-5);
  BOOST_CHECK_EQUAL(unassigned_element, "unassigned");
  BOOST_CHECK_EQUAL(unassigned_bead_type, "unassigned");
  BOOST_CHECK_EQUAL(unassigned_residue_type, "unassigned");
  BOOST_CHECK_EQUAL(unassigned_molecule_type, "unassigned");
  BOOST_CHECK_EQUAL(unassigned_residue_id, -1);
  BOOST_CHECK_EQUAL(unassigned_molecule_id, -1);
}

BOOST_AUTO_TEST_SUITE_END()
