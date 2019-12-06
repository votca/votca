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

#define BOOST_TEST_MODULE pdbreader_test
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include "../../include/votca/csg/openbox.h"
#include "../../include/votca/csg/triclinicbox.h"
#include "../../include/votca/csg/orthorhombicbox.h"

#include <vector>
#include <memory>
#include <iostream>

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(boundarycondition_test)

BOOST_AUTO_TEST_CASE(test_boundarycondition_initiatialization) {
    vector<unique_ptr<BoundaryCondition>> boundaries;
    
    boundaries.push_back(unique_ptr<OpenBox>(new OpenBox));   
    boundaries.push_back(unique_ptr<TriclinicBox>(new TriclinicBox));   
    boundaries.push_back(unique_ptr<OrthorhombicBox>(new OrthorhombicBox));   

    BOOST_CHECK_EQUAL(boundaries.at(0)->getBoxType(),BoundaryCondition::eBoxtype::typeOpen);
    BOOST_CHECK_EQUAL(boundaries.at(1)->getBoxType(),BoundaryCondition::eBoxtype::typeTriclinic);
    BOOST_CHECK_EQUAL(boundaries.at(2)->getBoxType(),BoundaryCondition::eBoxtype::typeOrthorhombic);
}

BOOST_AUTO_TEST_CASE(test_boundarycondition_boxvolume) {
    vector<unique_ptr<BoundaryCondition>> boundaries;
    
    boundaries.push_back(unique_ptr<OpenBox>(new OpenBox));   
    boundaries.push_back(unique_ptr<TriclinicBox>(new TriclinicBox));   
    boundaries.push_back(unique_ptr<OrthorhombicBox>(new OrthorhombicBox));   

    Eigen::Matrix3d box;
    box << 0.0, 0.0, 0.0, 
           0.0, 0.0, 0.0,
           0.0, 0.0, 0.0;

    boundaries.at(0)->setBox(box);
    boundaries.at(1)->setBox(box);
    boundaries.at(2)->setBox(box);

    BOOST_CHECK_EQUAL(boundaries.at(0)->BoxVolume(),0.0);
    BOOST_CHECK_EQUAL(boundaries.at(1)->BoxVolume(),0.0);
    BOOST_CHECK_EQUAL(boundaries.at(2)->BoxVolume(),0.0);

    box(0,0) = 1.0;
    box(1,1) = 2.0;
    box(2,2) = 1.0;

    boundaries.at(0)->setBox(box);
    boundaries.at(1)->setBox(box);
    boundaries.at(2)->setBox(box);
    
    BOOST_CHECK_EQUAL(boundaries.at(0)->BoxVolume(),2.0);
    BOOST_CHECK_EQUAL(boundaries.at(1)->BoxVolume(),2.0);
    BOOST_CHECK_EQUAL(boundaries.at(2)->BoxVolume(),2.0);
}

BOOST_AUTO_TEST_CASE(test_boundarycondition_clone) {
    vector<unique_ptr<BoundaryCondition>> boundaries;
    
    boundaries.push_back(unique_ptr<OpenBox>(new OpenBox));   
    boundaries.push_back(unique_ptr<TriclinicBox>(new TriclinicBox));   
    boundaries.push_back(unique_ptr<OrthorhombicBox>(new OrthorhombicBox));   

    Eigen::Matrix3d box;
    box << 0.0, 0.0, 0.0, 
           0.0, 0.0, 0.0,
           0.0, 0.0, 0.0;

    boundaries.at(0)->setBox(box);
    boundaries.at(1)->setBox(box);
    boundaries.at(2)->setBox(box);

    vector<unique_ptr<BoundaryCondition>> boundaries_copy;

    boundaries_copy.push_back(boundaries.at(0)->Clone());
    boundaries_copy.push_back(boundaries.at(1)->Clone());
    boundaries_copy.push_back(boundaries.at(2)->Clone());

    BOOST_CHECK_EQUAL(boundaries_copy.at(0)->getBoxType(),BoundaryCondition::eBoxtype::typeOpen);
    BOOST_CHECK_EQUAL(boundaries_copy.at(1)->getBoxType(),BoundaryCondition::eBoxtype::typeTriclinic);
    BOOST_CHECK_EQUAL(boundaries_copy.at(2)->getBoxType(),BoundaryCondition::eBoxtype::typeOrthorhombic);

    box(0,0) = 1.0;
    box(1,1) = 2.0;
    box(2,2) = 1.0;

    boundaries_copy.at(0)->setBox(box);
    boundaries_copy.at(1)->setBox(box);
    boundaries_copy.at(2)->setBox(box);
    
    BOOST_CHECK_EQUAL(boundaries_copy.at(0)->BoxVolume(),2.0);
    BOOST_CHECK_EQUAL(boundaries_copy.at(1)->BoxVolume(),2.0);
    BOOST_CHECK_EQUAL(boundaries_copy.at(2)->BoxVolume(),2.0);

    /* Ensure that the original boundaries were not altered */
    BOOST_CHECK_EQUAL(boundaries.at(0)->BoxVolume(),0.0);
    BOOST_CHECK_EQUAL(boundaries.at(1)->BoxVolume(),0.0);
    BOOST_CHECK_EQUAL(boundaries.at(2)->BoxVolume(),0.0);

}
BOOST_AUTO_TEST_SUITE_END()
