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

#define BOOST_TEST_MODULE triplelist_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/beadtriple.h>
#include <votca/csg/triplelist.h>
#include <votca/csg/topology.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(triplelist_test)

BOOST_AUTO_TEST_CASE(triplelist_constructor) {
    TripleList triplelist;
}

BOOST_AUTO_TEST_CASE(triplelist_add_triple) {
    TripleList triplelist;
    
    Topology top;

    string bead_type_name = "CG";
    BeadType * b_type = top.GetOrCreateBeadType(bead_type_name);

    int symmetry = 1;
    string name = "dummy1";
    int resnr = 0;
    double mass = 1.0;
    double charge = -1.0;

    top.CreateBead(symmetry,name,b_type,resnr,mass,charge);
    
    symmetry = 1;
    name = "dummy2";
    resnr = 0;
    mass = 2.0;
    charge = -2.0;
    
    top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

    symmetry = 1;
    name = "dummy3";
    resnr = 0;
    mass = 3.0;
    charge = -3.0;
    
    top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

    vec dist12(0.1,0.2,0.3);
    vec dist13(0.2,0.4,0.3);
    vec dist23(0.1,0.2,0.0);

    BeadTriple *testtriple(top.getBead(0),top.getBead(1),top.getBead(2),dist12,dist13,dist23);
    
    triplelist.AddTriple(testtriple);
    
    BeadTriple *triplefront, tripleback;
    
    *triplefront = triplelist.front();
    BOOST_CHECK(triplefront->bead1->mass,1.0);
    BOOST_CHECK(triplefront->bead2->mass,2.0);
    BOOST_CHECK(triplefront->bead3->mass,3.0);
    BOOST_CHECK(triplefront->bead1->name,"dummy1");
    BOOST_CHECK(triplefront->bead2->name,"dummy2");
    BOOST_CHECK(triplefront->bead3->name,"dummy3");
    
}

BOOST_AUTO_TEST_SUITE_END()
