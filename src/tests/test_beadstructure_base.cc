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

#define BOOST_TEST_MODULE beadstructure_test
#include <boost/test/unit_test.hpp>

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

BOOST_AUTO_TEST_SUITE(beadstructure_test)

BOOST_AUTO_TEST_CASE(test_beadstructure_constructor) {
  BeadStructure beadstructure;
}

BOOST_AUTO_TEST_CASE(test_beadstructure_beadcount) {
  BeadStructure beadstructure;
  BOOST_CHECK_EQUAL(beadstructure.BeadCount(), 0);
}

BOOST_AUTO_TEST_CASE(test_beadstructure_add_and_getbead) {
  BeadStructure beadstructure;
  TestBead testbead;
  testbead.setId(2);
  beadstructure.AddBead(&testbead);
  BOOST_CHECK_EQUAL(beadstructure.BeadCount(), 1);
  auto testbead2 = beadstructure.getBead(2);
  BOOST_CHECK_EQUAL(testbead2->getId(), testbead.getId());
}

BOOST_AUTO_TEST_CASE(test_beadstructure_ConnectBeads) {
  BeadStructure beadstructure;
  TestBead testbead1;
  testbead1.setId(1);
  testbead1.setName("Carbon");
  TestBead testbead2;
  testbead2.setId(2);
  testbead2.setName("Carbon");
  beadstructure.AddBead(&testbead1);
  beadstructure.AddBead(&testbead2);
  beadstructure.ConnectBeads(1, 2);
}

BOOST_AUTO_TEST_CASE(test_beadstructure_isSingleStructure) {
  BeadStructure beadstructure;

  TestBead testbead1;
  testbead1.setName("Carbon");
  testbead1.setId(1);

  TestBead testbead2;
  testbead2.setName("Carbon");
  testbead2.setId(2);

  TestBead testbead3;
  testbead3.setName("Oxygen");
  testbead3.setId(3);

  TestBead testbead4;
  testbead4.setName("Hydrogen");
  testbead4.setId(4);

  TestBead testbead5;
  testbead5.setName("Hydrogen");
  testbead5.setId(5);

  beadstructure.AddBead(&testbead1);
  beadstructure.AddBead(&testbead2);
  BOOST_CHECK(!beadstructure.isSingleStructure());

  // C - C
  beadstructure.ConnectBeads(1, 2);
  BOOST_CHECK(beadstructure.isSingleStructure());

  // C - C  O
  beadstructure.AddBead(&testbead3);
  BOOST_CHECK(!beadstructure.isSingleStructure());

  // C - C - O
  beadstructure.ConnectBeads(2, 3);
  BOOST_CHECK(beadstructure.isSingleStructure());

  // C - C - O  H - H
  beadstructure.AddBead(&testbead4);
  beadstructure.AddBead(&testbead5);
  beadstructure.ConnectBeads(4, 5);
  BOOST_CHECK(!beadstructure.isSingleStructure());
}

BOOST_AUTO_TEST_CASE(test_beadstructure_isStructureEquivalent) {
  BeadStructure beadstructure1;
  BeadStructure beadstructure2;

  // Beads for bead structure 1
  TestBead testbead1;
  testbead1.setName("Carbon");
  testbead1.setId(1);

  TestBead testbead2;
  testbead2.setName("Carbon");
  testbead2.setId(2);

  TestBead testbead3;
  testbead3.setName("Oxygen");
  testbead3.setId(3);

  TestBead testbead4;
  testbead4.setName("Hydrogen");
  testbead4.setId(4);

  TestBead testbead5;
  testbead5.setName("Hydrogen");
  testbead5.setId(5);

  // Beads for bead structure 2
  TestBead testbead6;
  testbead6.setName("Carbon");
  testbead6.setId(6);

  TestBead testbead7;
  testbead7.setName("Carbon");
  testbead7.setId(7);

  TestBead testbead8;
  testbead8.setName("Oxygen");
  testbead8.setId(8);

  TestBead testbead9;
  testbead9.setName("Hydrogen");
  testbead9.setId(9);

  TestBead testbead10;
  testbead10.setName("Hydrogen");
  testbead10.setId(10);

  BOOST_CHECK(beadstructure1.isStructureEquivalent(beadstructure2));
  beadstructure1.AddBead(&testbead1);
  BOOST_CHECK(!beadstructure1.isStructureEquivalent(beadstructure2));
  beadstructure2.AddBead(&testbead6);
  BOOST_CHECK(beadstructure1.isStructureEquivalent(beadstructure2));

  beadstructure1.AddBead(&testbead2);
  beadstructure2.AddBead(&testbead7);

  beadstructure1.ConnectBeads(1, 2);
  BOOST_CHECK(!beadstructure1.isStructureEquivalent(beadstructure2));
  beadstructure2.ConnectBeads(6, 7);
  BOOST_CHECK(beadstructure1.isStructureEquivalent(beadstructure2));

  beadstructure1.AddBead(&testbead3);
  beadstructure1.AddBead(&testbead4);
  beadstructure1.AddBead(&testbead5);
  beadstructure1.ConnectBeads(2, 3);
  beadstructure1.ConnectBeads(4, 5);

  beadstructure2.AddBead(&testbead10);
  beadstructure2.AddBead(&testbead8);
  beadstructure2.AddBead(&testbead9);
  beadstructure2.ConnectBeads(7, 8);
  beadstructure2.ConnectBeads(9, 10);
  BOOST_CHECK(beadstructure1.isStructureEquivalent(beadstructure2));
}

BOOST_AUTO_TEST_CASE(test_beadstructure_getNeighBeads) {
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

  beadstructure1.AddBead(&testbead1);
  beadstructure1.AddBead(&testbead2);
  beadstructure1.AddBead(&testbead3);
  beadstructure1.AddBead(&testbead4);
  beadstructure1.AddBead(&testbead5);
  beadstructure1.AddBead(&testbead6);
  beadstructure1.AddBead(&testbead7);
  beadstructure1.AddBead(&testbead8);

  // At this point non of the beads are connected so should return a vector of
  // size 0
  auto v1 = beadstructure1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v1.size(), 0);
  auto v2 = beadstructure1.getNeighBeads(2);
  BOOST_CHECK_EQUAL(v2.size(), 0);
  auto v3 = beadstructure1.getNeighBeads(3);
  BOOST_CHECK_EQUAL(v3.size(), 0);
  auto v4 = beadstructure1.getNeighBeads(4);
  BOOST_CHECK_EQUAL(v4.size(), 0);
  auto v5 = beadstructure1.getNeighBeads(5);
  BOOST_CHECK_EQUAL(v5.size(), 0);
  auto v6 = beadstructure1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v6.size(), 0);
  auto v7 = beadstructure1.getNeighBeads(7);
  BOOST_CHECK_EQUAL(v7.size(), 0);
  auto v8 = beadstructure1.getNeighBeads(8);
  BOOST_CHECK_EQUAL(v8.size(), 0);

  // Connect beads
  beadstructure1.ConnectBeads(1, 2);
  beadstructure1.ConnectBeads(3, 2);
  beadstructure1.ConnectBeads(4, 2);
  beadstructure1.ConnectBeads(5, 2);
  beadstructure1.ConnectBeads(6, 7);
  beadstructure1.ConnectBeads(7, 8);

  v1 = beadstructure1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v1.size(), 1);
  v2 = beadstructure1.getNeighBeads(2);
  BOOST_CHECK_EQUAL(v2.size(), 4);
  v3 = beadstructure1.getNeighBeads(3);
  BOOST_CHECK_EQUAL(v3.size(), 1);
  v4 = beadstructure1.getNeighBeads(4);
  BOOST_CHECK_EQUAL(v4.size(), 1);
  v5 = beadstructure1.getNeighBeads(5);
  BOOST_CHECK_EQUAL(v5.size(), 1);
  v6 = beadstructure1.getNeighBeads(1);
  BOOST_CHECK_EQUAL(v6.size(), 1);
  v7 = beadstructure1.getNeighBeads(7);
  BOOST_CHECK_EQUAL(v7.size(), 2);
  v8 = beadstructure1.getNeighBeads(8);
  BOOST_CHECK_EQUAL(v8.size(), 1);
}

BOOST_AUTO_TEST_CASE(test_beadstructure_catchError) {

  {
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

    TestBead testbead6;
    testbead6.setName("Hydrogen");
    testbead6.setId(5);

    BeadStructure beadstructure;
    beadstructure.AddBead(&testbead1);
    beadstructure.AddBead(&testbead2);
    beadstructure.AddBead(&testbead3);
    beadstructure.AddBead(&testbead4);
    beadstructure.AddBead(&testbead5);

    BOOST_CHECK_THROW(beadstructure.AddBead(&testbead6), invalid_argument);
    BOOST_CHECK_THROW(beadstructure.ConnectBeads(0, 1), invalid_argument);
    BOOST_CHECK_THROW(beadstructure.ConnectBeads(5, 6), invalid_argument);
    BOOST_CHECK_THROW(beadstructure.ConnectBeads(1, 1), invalid_argument);
  }
}
BOOST_AUTO_TEST_SUITE_END()
