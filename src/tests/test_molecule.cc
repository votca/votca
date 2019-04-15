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

#define BOOST_TEST_MODULE molecule_test
#include <boost/test/unit_test.hpp>
#include <vector>
#include <votca/tools/vec.h>
#include <votca/xtp/atom.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/molecule.h>
#include <votca/xtp/segment.h>

using namespace votca::tools;
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(molecule_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  Molecule mol(1, "molecule");
  Molecule mol2;
}

BOOST_AUTO_TEST_CASE(simple_getters_setters_test) {

  Molecule mol(1, "molecule");
  BOOST_CHECK_EQUAL(mol.getId(), 1);
  BOOST_CHECK_EQUAL(mol.getName(), "molecule");
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

  Atom *atm = new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
  atm->setPos(pos);

  Molecule mol(2, "molecule");
  mol.AddAtom(atm);
  vector<Atom *> v_atoms = mol.Atoms();
  BOOST_CHECK_EQUAL(v_atoms.size(), 1);
  delete atm;
}

BOOST_AUTO_TEST_CASE(add_segment_test) {

  bool hasQM = true;
  vec qmpos;
  qmpos.setX(2.0);
  qmpos.setY(2.0);
  qmpos.setZ(2.0);

  vec pos;
  pos.setX(3.0);
  pos.setY(3.0);
  pos.setZ(3.0);

  Atom *atm = new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
  atm->setPos(pos);
  Segment seg(4, "seg4");
  seg.AddAtom(atm);
  vector<Atom *> v_atoms = seg.Atoms();
  BOOST_CHECK_EQUAL(v_atoms.size(), 1);
  // Calculates the MD position of the segment (center of mass)
  seg.calcPos();

  Molecule mol(3, "molecule");
  mol.AddSegment(&seg);

  vector<Segment *> v_seg = mol.Segments();
  BOOST_CHECK_EQUAL(v_seg.size(), 1);
  auto pos_seg = v_seg.at(0)->getPos();
  BOOST_CHECK_EQUAL(pos_seg.getX(), 3.0);
  BOOST_CHECK_EQUAL(pos_seg.getY(), 3.0);
  BOOST_CHECK_EQUAL(pos_seg.getZ(), 3.0);
  delete atm;
}

BOOST_AUTO_TEST_CASE(add_fragment_test) {
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
  Atom *atm = new Atom(nullptr, "res1", 1, "CSP", 2, hasQM, 3, qmpos, "C", 1.0);
  atm->setPos(pos);

  Fragment frag(3, "frag3");
  frag.AddAtom(atm);
  // Default is MD
  frag.calcPos();
  auto v2 = frag.getPos();
  BOOST_CHECK_EQUAL(v2.x(), 3.0);
  BOOST_CHECK_EQUAL(v2.y(), 3.0);
  BOOST_CHECK_EQUAL(v2.z(), 3.0);

  Molecule mol(4, "molecule");
  mol.AddFragment(&frag);
  vector<Fragment *> v_frag = mol.Fragments();
  BOOST_CHECK_EQUAL(v_frag.size(), 1);
  v2 = v_frag.at(0)->getPos();
  BOOST_CHECK_EQUAL(v2.x(), 3.0);
  BOOST_CHECK_EQUAL(v2.y(), 3.0);
  BOOST_CHECK_EQUAL(v2.z(), 3.0);
  delete atm;
}

BOOST_AUTO_TEST_SUITE_END()
