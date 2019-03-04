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

#define BOOST_TEST_MODULE beadmotif_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/basebead.h>
#include <votca/csg/beadmotif.h>

using namespace std;
using namespace votca::csg;

class TestBead : public BaseBead {
 public:
  TestBead() : BaseBead(){};
};

BOOST_AUTO_TEST_SUITE(beadmotif_test)

BOOST_AUTO_TEST_CASE(test_beadmotif_constructor) { BeadMotif beadmotif; }

BOOST_AUTO_TEST_CASE(test_beadmotif_beadcount) {
  BeadMotif beadmotif;
  BOOST_CHECK_EQUAL(beadmotif.BeadCount(), 0);
}

BOOST_AUTO_TEST_CASE(test_beadmotif_getType) {
  BeadMotif beadmotif;
  auto      type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::empty);

  TestBead testbead;
  testbead.setId(2);
  testbead.setName("Helium");
  beadmotif.AddBead(&testbead);
  BOOST_CHECK_EQUAL(beadmotif.BeadCount(), 1);
  type = beadmotif.getType();

  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::single_bead);

  TestBead testbead2;
  testbead2.setId(3);
  testbead2.setName("Helium");
  beadmotif.AddBead(&testbead2);
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::multiple_structures);
}

BOOST_AUTO_TEST_CASE(test_beadmotif_getType2) {
  BeadMotif beadmotif;
  auto      type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::empty);
  TestBead testbead1;
  testbead1.setId(1);
  testbead1.setName("Carbon");
  beadmotif.AddBead(&testbead1);

  // C1
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::single_bead);

  // C1   C2
  TestBead testbead2;
  testbead2.setId(2);
  testbead2.setName("Carbon");
  beadmotif.AddBead(&testbead2);
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::multiple_structures);

  // C1 - C2
  beadmotif.ConnectBeads(1, 2);
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::line);

  TestBead testbead3;
  testbead3.setId(3);
  testbead3.setName("Carbon");
  beadmotif.AddBead(&testbead3);
  beadmotif.ConnectBeads(2, 3);

  // C1 - C2 - C3
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::line);

  beadmotif.ConnectBeads(3, 1);

  // C1 - C2
  //   \  :
  //     C3
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::loop);

  TestBead testbead4;
  testbead4.setId(4);
  testbead4.setName("Carbon");
  beadmotif.AddBead(&testbead4);

  beadmotif.ConnectBeads(1, 4);

  // C4
  // |
  // C1 - C2
  //   \  :
  //     C3
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::single_structure);

  TestBead testbead5;
  testbead5.setId(5);
  testbead5.setName("Carbon");
  beadmotif.AddBead(&testbead5);

  beadmotif.ConnectBeads(5, 4);
  beadmotif.ConnectBeads(5, 2);

  // C4 - C5
  // |    |
  // C1 - C2
  //   \  :
  //     C3
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::fused_ring);

  TestBead testbead6;
  testbead6.setId(6);
  testbead6.setName("Carbon");
  beadmotif.AddBead(&testbead6);

  TestBead testbead7;
  testbead7.setId(7);
  testbead7.setName("Carbon");
  beadmotif.AddBead(&testbead7);

  beadmotif.ConnectBeads(5, 6);
  beadmotif.ConnectBeads(6, 7);
  beadmotif.ConnectBeads(7, 5);

  // C4 - C5 - C6
  // |    |  \ |
  // C1 - C2   C7
  //   \  :
  //     C3
  // This is not a fused ring
  type = beadmotif.getType();
  BOOST_CHECK_EQUAL(type, BeadMotif::MotifType::single_structure);
}

BOOST_AUTO_TEST_SUITE_END()
