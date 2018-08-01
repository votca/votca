/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE interaction_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/interaction.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/topology.h>
#include <votca/csg/molecule.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(interaction_test)

BOOST_AUTO_TEST_CASE(test_interaction_constructor) {
  IBond bond1(1,2);
  IAngle angle1(1,2,3);
  IDihedral dihedral(1,2,3,4);
}

BOOST_AUTO_TEST_CASE(test_interaction_setters_getters){

  IBond bond1(1,2);
  bond1.setGroup("large");
  bond1.setGroupId(1);
  bond1.setIndex(1);
  bond1.setMolecule(1);
  
  string name = bond1.getName();
  cout << name << endl;
  bool correctName = name.compare("molecule 1:large 1:index 1")==0;
  BOOST_CHECK(correctName);
  BOOST_CHECK_EQUAL(bond1.getGroupId(),1);  
  BOOST_CHECK_EQUAL(bond1.getIndex(),1);  
  BOOST_CHECK_EQUAL(bond1.getMolecule(),1);  
  BOOST_CHECK_EQUAL(bond1.BeadCount(),2);  
  BOOST_CHECK_EQUAL(bond1.getBeadId(0),1);  
  BOOST_CHECK_EQUAL(bond1.getBeadId(1),2);  
  string groupName = bond1.getGroup();
  correctName = groupName.compare("large")==0;
  
  BOOST_CHECK(correctName);

  IAngle angle1(1,2,3);
  angle1.setGroup("medium");
  angle1.setGroupId(1);
  angle1.setIndex(1);
  angle1.setMolecule(1);

  name = angle1.getName();
  cout << name << endl;
  correctName = name.compare("molecule 1:medium 1:index 1")==0;
  BOOST_CHECK(correctName);
  BOOST_CHECK_EQUAL(angle1.getGroupId(),1);  
  BOOST_CHECK_EQUAL(angle1.getIndex(),1);  
  BOOST_CHECK_EQUAL(angle1.getMolecule(),1);  
  BOOST_CHECK_EQUAL(angle1.BeadCount(),3);  
  BOOST_CHECK_EQUAL(angle1.getBeadId(0),1);  
  BOOST_CHECK_EQUAL(angle1.getBeadId(1),2);  
  BOOST_CHECK_EQUAL(angle1.getBeadId(2),3);  
  groupName = angle1.getGroup();
  correctName = groupName.compare("medium")==0;
 
  IDihedral dihedral1(1,2,3,4);
  dihedral1.setGroup("small");
  dihedral1.setGroupId(1);
  dihedral1.setIndex(1);
  dihedral1.setMolecule(1);

  name = dihedral1.getName();
  cout << name << endl;
  correctName = name.compare("molecule 1:small 1:index 1")==0;
  BOOST_CHECK(correctName);
  BOOST_CHECK_EQUAL(dihedral1.getGroupId(),1);  
  BOOST_CHECK_EQUAL(dihedral1.getIndex(),1);  
  BOOST_CHECK_EQUAL(dihedral1.getMolecule(),1);  
  BOOST_CHECK_EQUAL(dihedral1.BeadCount(),4);  
  BOOST_CHECK_EQUAL(dihedral1.getBeadId(0),1);  
  BOOST_CHECK_EQUAL(dihedral1.getBeadId(1),2);  
  BOOST_CHECK_EQUAL(dihedral1.getBeadId(2),3);  
  BOOST_CHECK_EQUAL(dihedral1.getBeadId(3),4);  
  groupName = dihedral1.getGroup();
  correctName = groupName.compare("small")==0;
 
}

BOOST_AUTO_TEST_SUITE_END()
