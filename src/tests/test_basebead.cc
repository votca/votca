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

#define BOOST_TEST_MODULE basebead_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/basebead.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/topology.h>
#include <votca/csg/molecule.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

BOOST_AUTO_TEST_SUITE(basebead_test)

BOOST_AUTO_TEST_CASE(test_basebead_constructor) {
  BaseBead basebead;
}

BOOST_AUTO_TEST_CASE(test_basebead_getters_setters) {

  BaseBead basebead; 
  BOOST_CHECK_EQUAL(round_(basebead.getMass(),3),round_(0.0,3));
  BOOST_CHECK_EQUAL(round_(basebead.getQ(),3),round_(0.0,3));
  BOOST_CHECK(!basebead.HasPos());

  basebead.setId(0);
  BOOST_CHECK_EQUAL(basebead.getId(),0);
  
  basebead.setName("Bead1");
  string name = "Bead1";
  BOOST_CHECK(name == basebead.getName());

  basebead.setMass(1.0);
  BOOST_CHECK_EQUAL(round_(basebead.getMass(),3),round_(1.0,3));
  
  basebead.setQ(2.0);
  BOOST_CHECK_EQUAL(round_(basebead.getQ(),3),round_(2.0,3));
 
  vec xyz(-1.3,2.9,9.2);
  basebead.setPos(xyz);
  BOOST_CHECK(basebead.HasPos());
  auto xyz2 = basebead.getPos();
  BOOST_CHECK_EQUAL(round_(xyz2.x(),3),round_(-1.3,3));
  BOOST_CHECK_EQUAL(round_(xyz2.y(),3),round_(2.9,3));
  BOOST_CHECK_EQUAL(round_(xyz2.z(),3),round_(9.2,3));
  
  auto xyz3 = basebead.Pos();
  BOOST_CHECK_EQUAL(round_(xyz3.x(),3),round_(-1.3,3));
  BOOST_CHECK_EQUAL(round_(xyz3.y(),3),round_(2.9,3));
  BOOST_CHECK_EQUAL(round_(xyz3.z(),3),round_(9.2,3));

  Topology top;
  auto mol = top.CreateMolecule("Molecule1");
  basebead.setMolecule(mol);
  auto mol2 = basebead.getMolecule();
  bool molecules_equal = mol2->getName()=="Molecule1";
  BOOST_CHECK(molecules_equal);

}
BOOST_AUTO_TEST_SUITE_END()
