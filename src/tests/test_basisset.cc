/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE basisset_test

// Standard includes
#include <fstream>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aoshell.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(basisset_test)

BOOST_AUTO_TEST_CASE(FreeFunctions_test) {

  BOOST_CHECK_EQUAL(EnumToString(L::H), "H");
  BOOST_CHECK_EQUAL(votca::Index(StringToEnum("F")), votca::Index(L::F));

  BOOST_CHECK_EQUAL(OffsetFuncShell(L::S), 0);

  BOOST_CHECK_EQUAL(OffsetFuncShell(L::D), 4);

  BOOST_CHECK_EQUAL(NumFuncShell(L::I), 13);

  BOOST_CHECK_EQUAL(NumFuncShell_cartesian(L::G), 15);

  BOOST_CHECK_EQUAL(OffsetFuncShell_cartesian(L::P), 1);

  BOOST_CHECK_EQUAL(CheckShellType("FG"), true);
  BOOST_CHECK_EQUAL(CheckShellType("S"), true);
  BOOST_CHECK_EQUAL(CheckShellType("PD"), true);
  BOOST_CHECK_EQUAL(CheckShellType("s"), false);
  BOOST_CHECK_EQUAL(CheckShellType("apdf"), false);
  BOOST_CHECK_EQUAL(CheckShellType("PDS"), false);
  BOOST_CHECK_EQUAL(CheckShellType("PDG"), false);
  BOOST_CHECK_EQUAL(CheckShellType("SDFGI"), false);
}

BOOST_AUTO_TEST_SUITE_END()
