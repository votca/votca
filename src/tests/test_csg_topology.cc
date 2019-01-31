/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE csg_topology_test
#include <boost/test/floating_point_comparison.hpp>
#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../../include/votca/csg/topology.h"
#include "../../include/votca/csg/beadtype.h"



using namespace std;
using namespace votca::tools;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(csg_topology_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Topology top; }

BOOST_AUTO_TEST_CASE(box_test) { 

  // Box takes a vector
  double x1 = 2.0;
  double y1 = 0.0;
  double z1 = 0.0;

  double x2 = 0.0;
  double y2 = 2.0;
  double z2 = 0.0;

  double x3 = 0.0;
  double y3 = 0.0;
  double z3 = 2.0;

  vec v1(x1,y1,z1);
  vec v2(x2,y2,z2);
  vec v3(x3,y3,z3);

  matrix box(v1,v2,v3);

  Topology top;
  top.setBox(box);

  auto vol = top.BoxVolume();
  BOOST_CHECK_CLOSE(vol,8,1e-5);
  auto box2 = top.getBox();
  
  BOOST_CHECK_EQUAL(box2.isClose(box,1e-5),true);
}

BOOST_AUTO_TEST_CASE(simple_test){

  Topology top;
  top.setStep(1);
  BOOST_CHECK_EQUAL(top.getStep(),1);
  top.setTime(1.21);
  BOOST_CHECK_CLOSE(top.getTime(),1.21,1e-5);

}

BOOST_AUTO_TEST_CASE(create_bead_type){
  Topology top;
  string bead_type_name = "type1";
  std::weak_ptr<BeadType> weak_type = top.GetOrCreateBeadType(bead_type_name);

  auto bead_type_ptr=weak_type.lock();
  BOOST_CHECK(bead_type_ptr->getName()==bead_type_name);

  top.Cleanup();
}
/**
 * Test is creating a bead with the topology object and then ensuring that the
 * bead has the correct properties.
 **/
BOOST_AUTO_TEST_CASE(create_bead) {
  Topology top;
  // 1 - for spherical bead
  byte_t symmetry = 1;
  string bead_name = "bead_test";

  string bead_type_name = "type1";
  weak_ptr<BeadType> weak_type = top.GetOrCreateBeadType(bead_type_name);

  int residue_number = 1;
  double mass = 1.1;
  double charge = 0.3;

  auto bead_ptr = top.CreateBead(symmetry,
      bead_name,weak_type,residue_number,mass,charge);

  BOOST_CHECK_CLOSE(bead_ptr->getQ(),0.3,1e-5);
  BOOST_CHECK_CLOSE(bead_ptr->getMass(),1.1,1e-5);
  BOOST_CHECK_EQUAL(bead_ptr->getResnr(),residue_number);
  BOOST_CHECK_EQUAL(bead_ptr->getSymmetry(),symmetry);
  BOOST_CHECK(bead_ptr->getName() == bead_name);

  std::weak_ptr<BeadType> weak_type2 = bead_ptr->getType();
  auto bead_type_ptr = weak_type.lock();
  auto bead_type_ptr2=weak_type2.lock();
  BOOST_CHECK(bead_type_ptr2->getName() == bead_type_ptr->getName());
  BOOST_CHECK_EQUAL(top.BeadCount(),1);

  top.Cleanup();

}

/**
 * This test ensures that the interactions are stored correctly. Three beads
 * are created and two interactions are created connecting the three beads. The
 * interactions are then checked to be sure they hold the right information.
 **/
BOOST_AUTO_TEST_CASE(add_bonded_interation_test){
  Topology top;
  // 1 - for spherical bead
  byte_t symmetry = 1;

  string bead_type_name = "type1";
  std::weak_ptr<BeadType> weak_type = top.GetOrCreateBeadType(bead_type_name);

  int residue_number = 1;
  double mass = 1.1;
  double charge = 0.3;

  // Create 3 beads
  string bead_name = "bead_test";
  auto bead_ptr = top.CreateBead(symmetry,
      bead_name,weak_type,residue_number,mass,charge);
  bead_ptr->setId(0);

  string bead_name2 = "bead_test2";
  auto bead_ptr2 = top.CreateBead(symmetry,
      bead_name2,weak_type,residue_number,mass,charge);
  bead_ptr2->setId(1);

  string bead_name3 = "bead_test3";
  auto bead_ptr3 = top.CreateBead(symmetry,
      bead_name3,weak_type,residue_number,mass,charge);
  bead_ptr3->setId(2);

  BOOST_CHECK_EQUAL(top.BeadCount(),3);

  // Create two bonded interactions
  string interaction_group = "bond";
  auto bond1 = new IBond(0,1);
  bond1->setGroup(interaction_group);
  auto bond2 = new IBond(1,2);
  bond2->setGroup(interaction_group);

  top.AddBondedInteraction(bond1);
  top.AddBondedInteraction(bond2);

  auto interaction_container = top.BondedInteractions();
  BOOST_CHECK_EQUAL(interaction_container.size(),2);

  cout << "interaction name " << interaction_container.at(0)->getName() << endl;
  cout << "interaction name " << interaction_container.at(1)->getName() << endl;
  BOOST_CHECK_EQUAL(interaction_container.at(0)->getBeadId(0),0);
  BOOST_CHECK_EQUAL(interaction_container.at(0)->getBeadId(1),1);
  BOOST_CHECK_EQUAL(interaction_container.at(1)->getBeadId(0),1);
  BOOST_CHECK_EQUAL(interaction_container.at(1)->getBeadId(1),2);

  top.Cleanup();
}

BOOST_AUTO_TEST_SUITE_END()
