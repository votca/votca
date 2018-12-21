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

BOOST_AUTO_TEST_SUITE_END()
