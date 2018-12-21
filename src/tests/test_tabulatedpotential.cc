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

#define BOOST_TEST_MODULE tabulatedpotential_test
#include <boost/test/unit_test.hpp>

#include <map>
#include <cstdlib>
#include <vector>
#include <string>
#include "../csg_boltzmann/tabulatedpotential.h"

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

double randomDouble(double min_val, double max_val){
  double f = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  return min_val+f*(max_val-min_val); 
}


// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

BOOST_AUTO_TEST_SUITE(tabulatedpotential_test)

BOOST_AUTO_TEST_CASE(test_tabulatedpotential_constructor) { 
  TabulatedPotential tablulatedpotential; 
}

BOOST_AUTO_TEST_CASE(test_register){
  TabulatedPotential tablulatedpotential; 
  map<string,AnalysisTool *> commands;
  tablulatedpotential.Register(commands);
}

BOOST_AUTO_TEST_CASE(test_command){

  Topology top;
  BondedStatistics bonded_statistics;
  // Setup BondedStatistics Object
  {

    // Set the system size                                                       
    double x1 = 20.0;                                                            
    double y1 = 0.0;                                                             
    double z1 = 0.0;                                                             

    double x2 = 0.0;                                                             
    double y2 = 20.0;                                                            
    double z2 = 0.0;                                                             

    double x3 = 0.0;                                                             
    double y3 = 0.0;                                                             
    double z3 = 20.0;                                                            

    vec v1(x1,y1,z1);                                                            
    vec v2(x2,y2,z2);                                                            
    vec v3(x3,y3,z3);                                                            

    matrix box(v1,v2,v3);                                                        
    top.setBox(box);                                                             

    // Create three beads                                                        
    byte_t symmetry = 1;                                                         

    string bead_type_name = "H2";                                             
    auto bead_type_ptr = top.GetOrCreateBeadType(bead_type_name);                

    double mass = 0.9;                                                           
    double charge = 0.0;                                                         

    // Create a bunch of H2 molecules where the distance between them is 1.2
    string interaction_group = "covalent_bond";                                  
    double bond_length = 1.2;
    int number_of_H2 = 10;
    for(int index=0; index<number_of_H2;++index){
      int residue_number = index+1;

      string bead_name = to_string(index*2)+"H";     
      double x = randomDouble(0.0,x1-bond_length*2);
      double y = randomDouble(0.0,y2-bond_length*2);
      double z = randomDouble(0.0,z3-bond_length*2);
      vec bead_pos(x,y,z);  
      auto bead_ptr = top.CreateBead(symmetry,                                     
          bead_name,bead_type_ptr,residue_number,mass,charge);
      bead_ptr->setId(index*2);
      bead_ptr->setPos(bead_pos);

      string bead_name2 = to_string(index*2+1)+"H";     
      x += bond_length*(0.5+randomDouble(0.0,1.0));
      y += bond_length*(0.5+randomDouble(0.0,1.0));
      z += bond_length*(0.5+randomDouble(0.0,1.0));
      vec bead_pos2(x,y,z);  
      bead_ptr = top.CreateBead(symmetry,                                     
          bead_name2,bead_type_ptr,residue_number,mass,charge);
      bead_ptr->setId(index*2+1);
      bead_ptr->setPos(bead_pos2);

      auto bond = new IBond(index*2,index*2+1);                                                 
      bond->setGroup(interaction_group);                                          
      top.AddBondedInteraction(bond);                                             
       
    }
    bonded_statistics.BeginCG(&top,nullptr);
    bonded_statistics.EvalConfiguration(&top,nullptr);
  }

  TabulatedPotential tabulatedpotential; 
  map<string,AnalysisTool *> commands;
  tabulatedpotential.Register(commands);

  vector<string> arguments;
  arguments.push_back("set");
  arguments.push_back("smooth_pdf");
  arguments.push_back("2");

  string command = "hist";
  tabulatedpotential.Command(bonded_statistics,command,arguments);
  top.Cleanup();
}

BOOST_AUTO_TEST_SUITE_END()
