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

#define BOOST_TEST_MODULE nblist_3body_test
#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>
#include <votca/csg/bead.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/beadlist.h>
#include <votca/csg/nblistgrid_3body.h>
#include <votca/csg/topology.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(nblist_3body_test)

BOOST_AUTO_TEST_CASE(test_nblist_3body_constructor) {
    NBList_3Body nb;
}

BOOST_AUTO_TEST_CASE(test_nblist_3body_generate_list) {
    NBList_3Body *nb;
    nb = new NBListGrid_3Body();
    
    nb->setCutoff(2.0);    
    
    Topology top;
    
    matrix m;
    m.ZeroMatrix();
    m[0][0] = 5.0;
    m[1][1] = 5.0;
    m[2][2] = 5.0;    

    top.setBox(m);    

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
    
    vec pos;

    Bead *b=top.getBead(0);
    pos[0]=0.0;
    pos[1]=0.0;
    pos[2]=0.0;
    b->setPos(pos);
    
    b=top.getBead(1);
    pos[0]=1.0;
    pos[1]=0.0;
    pos[2]=0.0;
    b->setPos(pos);

    b = top.getBead(2);
    pos[0]=1.0;
    pos[1]=1.0;
    pos[2]=0.0;
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
