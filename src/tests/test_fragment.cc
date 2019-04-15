/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE fragment_test
#include <boost/test/unit_test.hpp>
#include <votca/tools/matrix.h>
#include <votca/tools/vec.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/fragment.h>

using namespace votca::xtp;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(fragment_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Fragment frag(1, "frag1"); }

BOOST_AUTO_TEST_CASE(getters_test) {
  Fragment frag(3, "frag3");
  BOOST_CHECK_EQUAL(frag.getId(), 3);
  BOOST_CHECK_EQUAL(frag.getName(), "frag3");
}

BOOST_AUTO_TEST_CASE(position_test) {

  BOOST_TEST_MESSAGE("Testing: setPos & getPos");
  {
    vec v;
    v.setX(0.0);
    v.setY(0.0);
    v.setZ(0.0);
    Fragment frag(3, "frag3");
    frag.setPos(v);
    auto v2 = frag.getPos();
    BOOST_CHECK_EQUAL(v2.x(), 0.0);
    BOOST_CHECK_EQUAL(v2.y(), 0.0);
    BOOST_CHECK_EQUAL(v2.z(), 0.0);
  }

  BOOST_TEST_MESSAGE("Testing: getCoMD");
  {
    vec v;
    v.setX(0.0);
    v.setX(1.0);
    v.setY(0.0);
    v.setZ(0.0);
    Fragment frag(3, "frag3");
    frag.setPos(v);
    auto v2 = frag.getCoMD();
    cerr << v2.x() << endl;
    BOOST_CHECK_EQUAL(v2.x(), 1.0);
    cerr << v2.y() << endl;
    BOOST_CHECK_EQUAL(v2.y(), 0.0);
    cerr << v2.z() << endl;
    BOOST_CHECK_EQUAL(v2.z(), 0.0);
  }

  BOOST_TEST_MESSAGE("Testing: getCoQM");
  {
    // The qmpos is not used to calculate the positions of the fragment
    // when calcPos is called with the MD tag which is also the default
    vec qmpos;
    qmpos.setX(2.0);
    qmpos.setY(2.0);
    qmpos.setZ(2.0);

    vec pos;
    pos.setX(3.0);
    pos.setY(3.0);
    pos.setZ(3.0);

    // If set to false cannot calculate the QM position of fragment when
    // the atom is added
    bool hasQM = true;
    Atom* atm =
        new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
    atm->setPos(pos);

    Fragment frag(3, "frag3");
    frag.AddAtom(atm);
    // Default is MD
    frag.calcPos();
    auto v2 = frag.getPos();
    BOOST_CHECK_EQUAL(v2.x(), 3.0);
    BOOST_CHECK_EQUAL(v2.y(), 3.0);
    BOOST_CHECK_EQUAL(v2.z(), 3.0);

    frag.calcPos("QM");
    v2 = frag.getCoQM();
    BOOST_CHECK_EQUAL(v2.x(), 2.0);
    BOOST_CHECK_EQUAL(v2.y(), 2.0);
    BOOST_CHECK_EQUAL(v2.z(), 2.0);
    delete atm;
  }
}

BOOST_AUTO_TEST_CASE(translate_test) {

  BOOST_TEST_MESSAGE("Testing: TranslateBy");
  vec v;
  v.setX(0.0);
  v.setY(0.0);
  v.setZ(0.0);
  Fragment frag(3, "frag3");
  frag.setPos(v);
  vec shift;
  shift.setX(1.0);
  shift.setY(-2.0);
  shift.setZ(4.0);
  frag.TranslateBy(shift);
  auto v2 = frag.getPos();
  BOOST_CHECK_EQUAL(v2.x(), 1.0);
  BOOST_CHECK_EQUAL(v2.y(), -2.0);
  BOOST_CHECK_EQUAL(v2.z(), 4.0);
}

BOOST_AUTO_TEST_CASE(ptr_set_test) {

  BOOST_TEST_MESSAGE("Testing: setTopology & setMolecule & setSegment");
  Fragment frag(3, "frag3");
  frag.setTopology(nullptr);
  frag.setMolecule(nullptr);
  frag.setSegment(nullptr);
}

BOOST_AUTO_TEST_CASE(rotate_pos_test) {

  BOOST_TEST_MESSAGE("Testing: Rotate");
  // Setting values for a rotation matrix
  // 90 degree rotation
  //
  // 1  0  0
  // 0  0 -1
  // 0  1  0
  vec row1;
  row1.setX(1);
  row1.setY(0);
  row1.setZ(0);
  vec row2;
  row2.setX(0);
  row2.setY(0);
  row2.setZ(-1);
  vec row3;
  row3.setX(0);
  row3.setY(1);
  row3.setZ(0);

  // Rotation center
  vec center;
  center.setX(0.0);
  center.setY(0.0);
  center.setZ(0.0);

  matrix rot_mat(row1, row2, row3);

  // The qmpos is not used to calculate the positions of the fragment
  // when calcPos is called with the MD tag which is also the default
  vec qmpos;
  qmpos.setX(2.0);
  qmpos.setY(2.0);
  qmpos.setZ(2.0);

  vec pos;
  pos.setX(3.0);
  pos.setY(3.0);
  pos.setZ(3.0);

  // If set to false cannot calculate the QM position of fragment when
  // the atom is added
  bool hasQM = true;
  Atom* atm = new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
  atm->setPos(pos);

  Fragment frag(3, "frag3");
  frag.AddAtom(atm);
  // Default is MD
  frag.calcPos();
  auto v2 = frag.getPos();
  BOOST_CHECK_EQUAL(v2.x(), 3.0);
  BOOST_CHECK_EQUAL(v2.y(), 3.0);
  BOOST_CHECK_EQUAL(v2.z(), 3.0);

  frag.calcPos("QM");
  v2 = frag.getCoQM();
  BOOST_CHECK_EQUAL(v2.x(), 2.0);
  BOOST_CHECK_EQUAL(v2.y(), 2.0);
  BOOST_CHECK_EQUAL(v2.z(), 2.0);

  frag.Rotate(rot_mat, center);
  // Demonstrate that the rotate function does not change the md position
  frag.calcPos();
  v2 = frag.getPos();
  BOOST_CHECK_EQUAL(v2.x(), 3.0);
  BOOST_CHECK_EQUAL(v2.y(), 3.0);
  BOOST_CHECK_EQUAL(v2.z(), 3.0);

  // Rotate does change the qm position
  frag.calcPos("QM");
  v2 = frag.getCoQM();
  BOOST_CHECK_EQUAL(v2.x(), 2.0);
  BOOST_CHECK_EQUAL(v2.y(), 2.0);
  BOOST_CHECK_EQUAL(v2.z(), -2.0);
  delete atm;
}

BOOST_AUTO_TEST_SUITE_END()
