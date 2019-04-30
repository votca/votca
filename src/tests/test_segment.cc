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

#define BOOST_TEST_MODULE segment_test
#include <boost/test/unit_test.hpp>
#include <vector>
#include <votca/tools/vec.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/segment.h>

using namespace votca::tools;
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(segment_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Segment seg(1, "seg1"); }

BOOST_AUTO_TEST_CASE(simple_getters_setters_test) {
  Segment seg(3, "seg2");
  BOOST_CHECK_EQUAL(seg.getId(), 3);
  BOOST_CHECK_EQUAL(seg.getName(), "seg2");
  vec v;
  v.setX(1.1);
  v.setY(2.2);
  v.setZ(3.3);
  seg.setPos(v);
  auto v2 = seg.getPos();
  BOOST_CHECK_EQUAL(v2.getX(), v.getX());
  BOOST_CHECK_EQUAL(v2.getY(), v.getY());
  BOOST_CHECK_EQUAL(v2.getZ(), v.getZ());
  seg.setOcc(1.1, 1);
  BOOST_CHECK_EQUAL(seg.getOcc(1), 1.1);
}

BOOST_AUTO_TEST_CASE(add_atom_test) {

  bool hasQM = true;
  vec qmpos;
  qmpos.setX(2.0);
  qmpos.setY(2.0);
  qmpos.setZ(2.0);

  vec pos;
  pos.setX(3.0);
  pos.setY(3.0);
  pos.setZ(3.0);

  Atom* atm = new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
  atm->setPos(pos);
  Segment seg(4, "seg4");
  seg.AddAtom(atm);
  vector<Atom*> v_atoms = seg.Atoms();
  BOOST_CHECK_EQUAL(v_atoms.size(), 1);
  delete atm;
}

BOOST_AUTO_TEST_CASE(calc_pos_test) {

  bool hasQM = true;
  vec qmpos;
  qmpos.setX(2.0);
  qmpos.setY(2.0);
  qmpos.setZ(2.0);

  vec pos;
  pos.setX(3.0);
  pos.setY(3.0);
  pos.setZ(3.0);

  Atom* atm = new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
  atm->setPos(pos);
  Segment seg(4, "seg4");
  seg.AddAtom(atm);
  vector<Atom*> v_atoms = seg.Atoms();
  BOOST_CHECK_EQUAL(v_atoms.size(), 1);
  // Calculates the MD position of the segment (center of mass)
  seg.calcPos();
  auto v_seg = seg.getPos();
  BOOST_CHECK_EQUAL(v_seg.getX(), 3.0);
  BOOST_CHECK_EQUAL(v_seg.getY(), 3.0);
  BOOST_CHECK_EQUAL(v_seg.getZ(), 3.0);
  delete atm;
}

BOOST_AUTO_TEST_SUITE_END()
