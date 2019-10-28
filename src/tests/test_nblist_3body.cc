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

#define BOOST_TEST_MODULE nblist_3body_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>
#include <votca/csg/bead.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/nblist_3body.h>
#include <votca/csg/topology.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(nblist_3body_test)

BOOST_AUTO_TEST_CASE(test_nblist_3body_constructor) { NBList_3Body nb; }

BOOST_AUTO_TEST_CASE(test_nblist_3body_generate_list) {
  NBList_3Body *nb;
  nb = new NBList_3Body();

  nb->setCutoff(2.0);

  Topology top;

  Eigen::Matrix3d m = 5 * Eigen::Matrix3d::Identity();

  top.setBox(m);

  Eigen::Vector3d pos = Eigen::Vector3d::Zero();

  Molecule *mol;
  mol = top.CreateMolecule("UNKNOWN");

  string bead_type_name = "CG";
  top.RegisterBeadType(bead_type_name);

  votca::votca::Index symmetry = 1;
  string name = "dummy1";
  votca::Index resnr = 0;
  double mass = 1.0;
  double charge = -1.0;
  Bead *b;
  b = top.CreateBead(symmetry, name, bead_type_name, resnr, mass, charge);
  b->setPos(pos);
  mol->AddBead(b, bead_type_name);
  b->setMolecule(mol);

  symmetry = 1;
  name = "dummy2";
  resnr = 0;
  mass = 2.0;
  charge = -2.0;
  b = top.CreateBead(symmetry, name, bead_type_name, resnr, mass, charge);
  mol->AddBead(b, bead_type_name);
  b->setMolecule(mol);
  pos = Eigen::Vector3d::UnitX();
  b->setPos(pos);

  symmetry = 1;
  name = "dummy3";
  resnr = 0;
  mass = 3.0;
  charge = -3.0;
  b = top.CreateBead(symmetry, name, bead_type_name, resnr, mass, charge);
  mol->AddBead(b, bead_type_name);
  b->setMolecule(mol);
  pos[0] = 1.0;
  pos[1] = 1.0;
  pos[2] = 0.0;
  b->setPos(pos);

  BeadList beads;
  beads.Generate(top, "CG");

  nb->Generate(beads, true);

  BOOST_CHECK_EQUAL(nb->size(), 3);

  NBList_3Body::iterator triple_iter;
  triple_iter = nb->begin();
  BOOST_CHECK_EQUAL((*triple_iter)->bead1()->getId(), 0);
  BOOST_CHECK_EQUAL((*triple_iter)->bead2()->getId(), 1);
  BOOST_CHECK_EQUAL((*triple_iter)->bead3()->getId(), 2);
  BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.0, 1e-4);
  BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.414214, 1e-4);
  BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.0, 1e-4);

  ++triple_iter;

  BOOST_CHECK_EQUAL((*triple_iter)->bead1()->getId(), 1);
  BOOST_CHECK_EQUAL((*triple_iter)->bead2()->getId(), 0);
  BOOST_CHECK_EQUAL((*triple_iter)->bead3()->getId(), 2);
  BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.0, 1e-4);
  BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.0, 1e-4);
  BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.414214, 1e-4);

  ++triple_iter;

  BOOST_CHECK_EQUAL((*triple_iter)->bead1()->getId(), 2);
  BOOST_CHECK_EQUAL((*triple_iter)->bead2()->getId(), 0);
  BOOST_CHECK_EQUAL((*triple_iter)->bead3()->getId(), 1);
  BOOST_CHECK_CLOSE((*triple_iter)->dist12(), 1.414214, 1e-4);
  BOOST_CHECK_CLOSE((*triple_iter)->dist13(), 1.0, 1e-4);
  BOOST_CHECK_CLOSE((*triple_iter)->dist23(), 1.0, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
