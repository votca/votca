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

#include "../../include/votca/csg/basebead.h"
#include "../../include/votca/csg/beadstructure.h"
#include "../../include/votca/csg/beadstructurealgorithms.h"  // IWYU pragma: keep

using namespace std;
using namespace votca::csg;

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

  // Adding a Helium
  //
  // He
  //

  TestBead testbead12;
  testbead12.setName("Helium");
  testbead12.setId(12);

  // Methane
  BeadStructure<BaseBead> beadstructure_methane;
  beadstructure_methane.AddBead(&testbead1);
  beadstructure_methane.AddBead(&testbead2);
  beadstructure_methane.AddBead(&testbead3);
  beadstructure_methane.AddBead(&testbead4);
  beadstructure_methane.AddBead(&testbead5);

  // Water
  BeadStructure<BaseBead> beadstructure_water;
  beadstructure_water.AddBead(&testbead6);
  beadstructure_water.AddBead(&testbead7);
  beadstructure_water.AddBead(&testbead8);

  // Helium
  BeadStructure<BaseBead> beadstructure_helium;
  beadstructure_helium.AddBead(&testbead12);

  // Methane and Water and Helium
  BeadStructure<BaseBead> beadstructure;
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
  beadstructure_methane.ConnectBeads(1, 2);
  beadstructure_methane.ConnectBeads(3, 2);
  beadstructure_methane.ConnectBeads(4, 2);
  beadstructure_methane.ConnectBeads(5, 2);

  // Water
  beadstructure_water.ConnectBeads(6, 7);
  beadstructure_water.ConnectBeads(7, 8);

  // Methane and Water
  beadstructure.ConnectBeads(1, 2);
  beadstructure.ConnectBeads(3, 2);
  beadstructure.ConnectBeads(4, 2);
  beadstructure.ConnectBeads(5, 2);
  beadstructure.ConnectBeads(6, 7);
  beadstructure.ConnectBeads(7, 8);
  cout << "Calling break into structures " << endl;
  vector<BeadStructure<BaseBead>> structures =
      breakIntoStructures(beadstructure);

  BOOST_CHECK_EQUAL(structures.size(), 2);
  cout << "Bead Count 1 " << structures.at(0).BeadCount() << endl;
  cout << "Bead Count 2 " << structures.at(1).BeadCount() << endl;
  bool methane_found = false;
  bool water_found = false;
  for (auto structure : structures) {
    if (structure.isStructureEquivalent(beadstructure_methane)) {
      methane_found = true;
    }
    if (structure.isStructureEquivalent(beadstructure_water)) {
      water_found = true;
    }
  }
  BOOST_CHECK(methane_found);
  BOOST_CHECK(water_found);

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

  // Adding a Helium
  TestBead testbead13;
  testbead13.setName("Helium");
  testbead13.setId(13);

  beadstructure.AddBead(&testbead13);

  structures = breakIntoStructures(beadstructure);

  BOOST_CHECK_EQUAL(structures.size(), 4);

  methane_found = false;
  water_found = false;
  bool structure3_found = false;
  int structure2_count = 0;
  for (auto structure : structures) {
    if (structure.isStructureEquivalent(beadstructure_methane)) {
      methane_found = true;
    }
    if (structure.isStructureEquivalent(beadstructure_water)) {
      water_found = true;
      ++structure2_count;
    }
    if (structure.isStructureEquivalent(beadstructure_helium)) {
      structure3_found = true;
    }
  }
  BOOST_CHECK(methane_found);
  BOOST_CHECK(water_found);
  BOOST_CHECK(structure3_found);
  BOOST_CHECK_EQUAL(structure2_count, 2);
}

BOOST_AUTO_TEST_SUITE_END()
