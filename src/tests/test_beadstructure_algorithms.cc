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

#define BOOST_TEST_MODULE beadstructurealgorithms_test
#include <boost/test/unit_test.hpp>

#include "../../include/votca/csg/beadstructurealgorithms.h"
#include <string>
#include <votca/csg/basebead.h>
#include <votca/csg/beadstructure.h>

using namespace std;
using namespace votca::csg;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

class TestBead : public BaseBead {
 public:
  TestBead() : BaseBead(){};
};

BOOST_AUTO_TEST_SUITE(beadstructurealgorithms_test)

BOOST_AUTO_TEST_CASE(test_beadstructure_breakIntoStructures) {

  // Beads for bead structure 1
  // Make a methane molecule
  //
  //     H
  //     |
  // H - C - H
  //     |
  //     H
  //
  TestBead testbead1;
  testbead1.setName("Hydrogen");
  testbead1.setId(1);

  TestBead testbead2;
  testbead2.setName("Carbon");
  testbead2.setId(2);

  TestBead testbead3;
  testbead3.setName("Hydrogen");
  testbead3.setId(3);

  TestBead testbead4;
  testbead4.setName("Hydrogen");
  testbead4.setId(4);

  TestBead testbead5;
  testbead5.setName("Hydrogen");
  testbead5.setId(5);

  // Make a Water molecule
  //
  // H - O - H
  //

  TestBead testbead6;
  testbead6.setName("Hydrogen");
  testbead6.setId(6);

  TestBead testbead7;
  testbead7.setName("Oxygen");
  testbead7.setId(7);

  TestBead testbead8;
  testbead8.setName("Hydrogen");
  testbead8.setId(8);

  // Methane
  BeadStructure beadstructure1;
  beadstructure1.AddBead(&testbead1);
  beadstructure1.AddBead(&testbead2);
  beadstructure1.AddBead(&testbead3);
  beadstructure1.AddBead(&testbead4);
  beadstructure1.AddBead(&testbead5);

  // Water
  BeadStructure beadstructure2;
  beadstructure2.AddBead(&testbead6);
  beadstructure2.AddBead(&testbead7);
  beadstructure2.AddBead(&testbead8);

  // Methane and Water
  BeadStructure beadstructure;
  beadstructure.AddBead(&testbead1);
  beadstructure.AddBead(&testbead2);
  beadstructure.AddBead(&testbead3);
  beadstructure.AddBead(&testbead4);
  beadstructure.AddBead(&testbead5);
  beadstructure.AddBead(&testbead6);
  beadstructure.AddBead(&testbead7);
  beadstructure.AddBead(&testbead8);

  // Connect beads
  // Methane
  beadstructure1.ConnectBeads(1, 2);
  beadstructure1.ConnectBeads(3, 2);
  beadstructure1.ConnectBeads(4, 2);
  beadstructure1.ConnectBeads(5, 2);

  // Water
  beadstructure2.ConnectBeads(6, 7);
  beadstructure2.ConnectBeads(7, 8);

  // Methane and Water
  beadstructure.ConnectBeads(1, 2);
  beadstructure.ConnectBeads(3, 2);
  beadstructure.ConnectBeads(4, 2);
  beadstructure.ConnectBeads(5, 2);
  beadstructure.ConnectBeads(6, 7);
  beadstructure.ConnectBeads(7, 8);
  vector<BeadStructure> structures = breakInToStructures(beadstructure);

  bool structure1_found = false;
  bool structure2_found = false;
  for (auto structure : structures) {
    if (structure.isStructureEquivalent(beadstructure1)) {
      structure1_found = true;
    }
    if (structure.isStructureEquivalent(beadstructure2)) {
      structure2_found = true;
    }
  }
  BOOST_CHECK(structure1_found);
  BOOST_CHECK(structure2_found);

  // Adding another water
  //
  // H - O - H
  //

  TestBead testbead9;
  testbead9.setName("Hydrogen");
  testbead9.setId(9);

  TestBead testbead11;
  testbead11.setName("Hydrogen");
  testbead11.setId(11);

  TestBead testbead10;
  testbead10.setName("Oxygen");
  testbead10.setId(10);

  // Adding the water
  beadstructure.AddBead(&testbead9);
  beadstructure.AddBead(&testbead10);
  beadstructure.AddBead(&testbead11);

  beadstructure.ConnectBeads(9, 10);
  beadstructure.ConnectBeads(11, 10);

  structures = breakInToStructures(beadstructure);

  structure1_found = false;
  structure2_found = false;
  int structure2_count = 0;
  for (auto structure : structures) {
    if (structure.isStructureEquivalent(beadstructure1)) {
      structure1_found = true;
    }
    if (structure.isStructureEquivalent(beadstructure2)) {
      structure2_found = true;
      ++structure2_count;
    }
  }
  BOOST_CHECK(structure1_found);
  BOOST_CHECK(structure2_found);
  BOOST_CHECK_EQUAL(structure2_count, 2);
}

BOOST_AUTO_TEST_SUITE_END();
