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

#define BOOST_TEST_MODULE bondedstatistics_test
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <map>
#include <string>
#include "../csg_boltzmann/bondedstatistics.h"

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

BOOST_AUTO_TEST_SUITE(bondedstatistics_test)

BOOST_AUTO_TEST_CASE(test_bondedstatistics_constructor) { 
  BondedStatistics bonded_statistics; 
}

BOOST_AUTO_TEST_CASE(test_bondedstatistics_begin) {
  Topology top;                                                                  
  // Create two bonded interactions                                              
  string interaction_group = "covalent_bond";         
  string interaction_group_compare = ":covalent_bond";  
  auto bond1 = new IBond(0,1);                                                   
  bond1->setGroup(interaction_group);                                            
  auto bond2 = new IBond(1,2);                                                   
  bond2->setGroup(interaction_group);                                            

  top.AddBondedInteraction(bond1);                                               
  top.AddBondedInteraction(bond2);                                               

  BondedStatistics bonded_statistics;
  bonded_statistics.BeginCG(&top,nullptr);

  DataCollection<double>& data_collection = bonded_statistics.BondedValues();
  vector<DataCollection<double>::array *>& vector_of_arrays = data_collection.Data(); 
  BOOST_CHECK_EQUAL(vector_of_arrays.size(),2);
  BOOST_CHECK_EQUAL(vector_of_arrays.at(0)->getName(),interaction_group_compare);
  BOOST_CHECK_EQUAL(vector_of_arrays.at(1)->getName(),interaction_group_compare);
  // The arrays do not store any numbers at this point
  BOOST_CHECK_EQUAL(vector_of_arrays.at(0)->size(),0);
  BOOST_CHECK_EQUAL(vector_of_arrays.at(1)->size(),0);
  top.Cleanup(); 
}

BOOST_AUTO_TEST_CASE(test_evalconfiguration_begin) {
  Topology top;     

  // Setup topology class
  { 
    // Set the system size
    double x1 = 10.0;                                                               
    double y1 = 0.0;                                                               
    double z1 = 0.0;                                                               

    double x2 = 0.0;                                                               
    double y2 = 10.0;                                                               
    double z2 = 0.0;                                                               

    double x3 = 0.0;                                                               
    double y3 = 0.0;                                                               
    double z3 = 10.0;                                                               

    vec v1(x1,y1,z1);                                                              
    vec v2(x2,y2,z2);                                                              
    vec v3(x3,y3,z3);                                                              

    matrix box(v1,v2,v3);         
    top.setBox(box);

    // Create three beads
    byte_t symmetry = 1;                                                           

    string bead_type_name = "type1";                                               
    auto bead_type_ptr = top.GetOrCreateBeadType(bead_type_name);                  

    int residue_number = 1;                                                        
    double mass = 1.1;                                                             
    double charge = 0.3;                                                           

    // Create 3 beads                                                              
    string bead_name = "bead_test";                                                
    vec pos_bead1(5.0,3.0,5.0);
    auto bead_ptr = top.CreateBead(symmetry,                                       
        bead_name,bead_type_ptr,residue_number,mass,charge);                       
    bead_ptr->setId(0);                                                            
    bead_ptr->setPos(pos_bead1);

    string bead_name2 = "bead_test2";                                              
    vec pos_bead2(5.0,4.0,5.0);
    auto bead_ptr2 = top.CreateBead(symmetry,                                      
        bead_name2,bead_type_ptr,residue_number,mass,charge);                      
    bead_ptr2->setId(1);                                                           
    bead_ptr2->setPos(pos_bead2);

    string bead_name3 = "bead_test3";                                              
    vec pos_bead3(5.0,6.0,5.0);
    auto bead_ptr3 = top.CreateBead(symmetry,                                      
        bead_name3,bead_type_ptr,residue_number,mass,charge);                      
    bead_ptr3->setId(2);                                                           
    bead_ptr3->setPos(pos_bead3);

    // Create two bonded interactions                                              
    string interaction_group = "covalent_bond";         
    auto bond1 = new IBond(0,1);                                                   
    bond1->setGroup(interaction_group);                                            
    auto bond2 = new IBond(1,2);                                                   
    bond2->setGroup(interaction_group);                                            

    top.AddBondedInteraction(bond1);                                               
    top.AddBondedInteraction(bond2);                                               
  }

  BondedStatistics bonded_statistics;
  bonded_statistics.BeginCG(&top,nullptr);

  // Calling EvalConfiguration on IBond structures will store the distances
  // between the beads in the BondedStatitics class
  bonded_statistics.EvalConfiguration(&top,nullptr);

  DataCollection<double>& data_collection = bonded_statistics.BondedValues();
  vector<DataCollection<double>::array *>& vector_of_arrays = data_collection.Data(); 
  
  // Distance between bead 0 and bead 1 
  BOOST_CHECK_EQUAL(vector_of_arrays.at(0)->at(0),1.0);
  // Distance between bead 1 and bead 2
  BOOST_CHECK_EQUAL(vector_of_arrays.at(1)->at(0),2.0);

  top.Cleanup(); 
}
BOOST_AUTO_TEST_SUITE_END()
