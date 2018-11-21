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

#define BOOST_TEST_MODULE beadtriple_test
#include <boost/test/unit_test.hpp>

#include <tuple>
#include <string>
#include <votca/csg/bead.h>
#include <votca/csg/beadtype.h>
#include <votca/csg/beadtriple.h>
#include <votca/csg/topology.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(beadtriple_test)

BOOST_AUTO_TEST_CASE(test_beadtriple_constructor) {

    Topology top;

    string bead_type_name = "C1";
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
    mass = 1.0;
    charge = -1.0;
    
    top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

    symmetry = 1;
    name = "dummy3";
    resnr = 0;
    mass = 1.0;
    charge = -1.0;
    
    top.CreateBead(symmetry,name,b_type,resnr,mass,charge);

    vec dist12(0.1,0.2,0.3);
    vec dist13(0.2,0.4,0.3);
    vec dist23(0.1,0.2,0.0);

    BeadTriple testtriple(top.getBead(0),top.getBead(1),top.getBead(2),dist12,dist13,dist23);

    double d12ref=0.3741657;
    double d13ref=0.5385165;
    double d23ref=0.2236068;
    
    d12 = testtriple.dist12();
    d13 = testtriple.dist12();
    d23 = testtriple.dist12();    

    bool equald12=d12ref.isApprox(d12,1e-5);
    if(!equald12){
        std::cout<<"result d12"<<std::endl;  
        std::cout<<d12<<std::endl;
        std::cout<<"ref d12"<<std::endl;
        std::cout<<d12ref<<std::endl;
    }
    BOOST_CHECK_EQUAL(equald12, true);

    bool equald13=d12ref.isApprox(d13,1e-5);
    if(!equald13){
        std::cout<<"result d13"<<std::endl;
        std::cout<<d13<<std::endl;
        std::cout<<"ref d13"<<std::endl;
        std::cout<<d13ref<<std::endl;
    }
    BOOST_CHECK_EQUAL(equald13, true);

    bool equald23=d23ref.isApprox(d23,1e-5);
    if(!equald23){
        std::cout<<"result d23"<<std::endl;
        std::cout<<d23<<std::endl;
        std::cout<<"ref d23"<<std::endl;
        std::cout<<d23ref<<std::endl;
    }
    BOOST_CHECK_EQUAL(equald23, true);
}

BOOST_AUTO_TEST_SUITE_END()
