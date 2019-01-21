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

#define BOOST_TEST_MODULE beadmotifalgorithms_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/beadmotifalgorithms.h>

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

BOOST_AUTO_TEST_SUITE(beadmotif_algorithms_test)

BOOST_AUTO_TEST_CASE(test_breakintomotifs) {

  BeadStructure beadstructure1;

  // Beads for bead structure 1
  // Make a methane molecule
  //
  //     H
  //     |
  // H - C - H
  //     |
  //     H
  //
  // Should return this as type single_structure

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

  beadstructure1.AddBead(&testbead1);
  beadstructure1.AddBead(&testbead2);
  beadstructure1.AddBead(&testbead3);
  beadstructure1.AddBead(&testbead4);
  beadstructure1.AddBead(&testbead5);

  beadstructure1.ConnectBeads(1, 2);
  beadstructure1.ConnectBeads(3, 2);
  beadstructure1.ConnectBeads(4, 2);
  beadstructure1.ConnectBeads(5, 2);
  // Make a Water molecule
  //
  // H - O - H
  //

  // Should turn into a motif of type edge

  TestBead testbead6;
  testbead6.setName("Hydrogen");
  testbead6.setId(6);

  TestBead testbead7;
  testbead7.setName("Oxygen");
  testbead7.setId(7);

  TestBead testbead8;
  testbead8.setName("Hydrogen");
  testbead8.setId(8);

  beadstructure1.AddBead(&testbead6);
  beadstructure1.AddBead(&testbead7);
  beadstructure1.AddBead(&testbead8);

  // Connect beads
  beadstructure1.ConnectBeads(6, 7);
  beadstructure1.ConnectBeads(7, 8);

  // Make a carbon ring
  //
  // C9 - C10
  // |    |
  // C11- C12
  //
  // Should return as type loop
  TestBead testbead9;
  testbead9.setName("Carbon");
  testbead9.setId(9);

  TestBead testbead10;
  testbead10.setName("Carbon");
  testbead10.setId(10);

  TestBead testbead11;
  testbead11.setName("Carbon");
  testbead11.setId(11);

  TestBead testbead12;
  testbead12.setName("Carbon");
  testbead12.setId(12);

  beadstructure1.AddBead(&testbead9);
  beadstructure1.AddBead(&testbead10);
  beadstructure1.AddBead(&testbead11);
  beadstructure1.AddBead(&testbead12);

  beadstructure1.ConnectBeads(9, 10);
  beadstructure1.ConnectBeads(9, 11);
  beadstructure1.ConnectBeads(10, 12);
  beadstructure1.ConnectBeads(11, 12);

  // Make a fused ring
  //
  // C13- C14- C15
  // |    |    |
  // C16- C17- C18
  //
  // Should return type fused ring

  vector<TestBead> fused_ring;
  for (int index = 0; index < 6; ++index) {
    int id = index + 13;
    TestBead temp;
    temp.setName("Carbon");
    temp.setId(id);
    fused_ring.push_back(temp);
  }
  for (int index = 0; index < 6; ++index) {
    beadstructure1.AddBead(&fused_ring.at(index));
  }

  beadstructure1.ConnectBeads(13, 14);
  beadstructure1.ConnectBeads(14, 15);
  beadstructure1.ConnectBeads(15, 18);
  beadstructure1.ConnectBeads(13, 16);
  beadstructure1.ConnectBeads(14, 17);
  beadstructure1.ConnectBeads(16, 17);
  beadstructure1.ConnectBeads(17, 18);

  // Make a single
  //
  // He19

  TestBead testbead19;
  testbead19.setName("Helium");
  testbead19.setId(19);
  beadstructure1.AddBead(&testbead19);

  list<BeadMotif> bead_motifs;
  bead_motifs = breakIntoMotifs<list<BeadMotif>>(beadstructure1);

  BOOST_CHECK_EQUAL(bead_motifs.size(), 5);

  bool found_type_line = false;
  bool found_type_loop = false;
  bool found_type_single_structure = false;
  bool found_type_single = false;
  bool found_type_fused_ring = false;

  for (BeadMotif& motif : bead_motifs) {
    if (motif.getType() == BeadMotif::MotifType::single_bead) {
      found_type_single = true;
    } else if (motif.getType() == BeadMotif::MotifType::line) {
      found_type_line = true;
    } else if (motif.getType() == BeadMotif::MotifType::loop) {
      found_type_loop = true;
    } else if (motif.getType() == BeadMotif::MotifType::fused_ring) {
      found_type_fused_ring = true;
    } else if (motif.getType() == BeadMotif::MotifType::single_structure) {
      found_type_single_structure = true;
    }
  }

  BOOST_CHECK(found_type_line);
  BOOST_CHECK(found_type_loop);
  BOOST_CHECK(found_type_single_structure);
  BOOST_CHECK(found_type_single);
  BOOST_CHECK(found_type_fused_ring);
}
/*
  BOOST_AUTO_TEST_CASE(test_breakintosimplemotifs){
  unordered_map<int, BeadMotif> breakIntoSimpleMotifs(BeadMotif beadmotif);

  }*/
BOOST_AUTO_TEST_SUITE_END()
