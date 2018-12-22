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
  double random_num = static_cast<double>(rand())/static_cast<double>(RAND_MAX);
  return min_val+random_num*(max_val-min_val); 
}

vector<double> getColumnFromFile(string file_name, int column){
  vector<double> data;
  ifstream file;
  file.open(file_name);
  string line;
  if(file.is_open()){
    while( getline(file,line)){
      string word;
      istringstream ss(line);
      for(int i=0;i<column;++i){
        ss >> word;
      }
      data.push_back(stod(word));
    }
    file.close();
  }
  return data;
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
  TabulatedPotential tabulatedpotential;
  BOOST_CHECK_EQUAL(tabulatedpotential.getTemperature(),300);
  auto iterations_to_smooth = tabulatedpotential.getSmoothIterations();
  BOOST_CHECK_EQUAL(iterations_to_smooth.first,0);
  BOOST_CHECK_EQUAL(iterations_to_smooth.second,0);
}

BOOST_AUTO_TEST_CASE(test_register){
  TabulatedPotential tablulatedpotential; 
  map<string,AnalysisTool *> commands;
  tablulatedpotential.Register(commands);
}

BOOST_AUTO_TEST_CASE(test_command){

  Topology top;
  BondedStatistics bonded_statistics;
  string interaction_group = "interaction";                                  
  string interaction_group_name = ":interaction";
  vector<string> interactions;
  interactions.push_back("file_interactions.txt");
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

    // Create a bunch of H2 molecules, each bead is considered a molecule they
    // are placed on a regular grid so the table properties can be compared 
    // consistently
    
    int number_of_H2 = 0;
    int residue_number = 0;
    for(double x=2.0; x<(x1-2.0); x+=4.0){
      for(double y=2.0; y<(y2-2.0); y+=3.0){
        for(double z=2.0;z<(z3-2.0);z+=4.0){
          residue_number++;

          string bead_name =to_string(number_of_H2)+"_H2";
          vec bead_pos(x,y,z);  
          auto bead_ptr = top.CreateBead(symmetry,                                     
              bead_name,bead_type_ptr,residue_number,mass,charge);
          bead_ptr->setId(number_of_H2);
          bead_ptr->setPos(bead_pos);
          number_of_H2++;
        }
      }
    }

    cout << "Number of H2 " << number_of_H2 << endl;
    // Every molecule interacts with every other molecule
    for(int index=0; index<number_of_H2;++index){
      for(int index2=index+1;index2<number_of_H2;++index2){
        auto bond = new IBond(index,index2);                                                 
        bond->setGroup(interaction_group+to_string(index)+"_"+to_string(index2));                                          
        top.AddBondedInteraction(bond);                                             
        interactions.push_back(interaction_group_name+to_string(index)+"_"+to_string(index2));
      }
    }

    bonded_statistics.BeginCG(&top,nullptr);
    bonded_statistics.EvalConfiguration(&top,nullptr);
  }// End of setup

  DataCollection<double> & bonded_values = bonded_statistics.BondedValues();
  cout << "bonded_values after pulling out of statistics " << bonded_values.size() << endl;
  size_t number_of_interactions = bonded_values.size();

  TabulatedPotential tabulatedpotential; 
  map<string,AnalysisTool *> commands;
  tabulatedpotential.Register(commands);

  // Test 1
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    string command = "tab";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);
    vector<double> column3 = getColumnFromFile(interactions.at(0),3);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),3.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),5.86);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),8.73);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),11.59);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),14.46);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),5.16);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),1.84);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),0.00);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),0.36);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),5.70);

    BOOST_CHECK_EQUAL(round_(column3.at(0),2),1.16);
    BOOST_CHECK_EQUAL(round_(column3.at(1),2),0.90);
    BOOST_CHECK_EQUAL(round_(column3.at(2),2),0.26);
    BOOST_CHECK_EQUAL(round_(column3.at(3),2),-1.00);
    BOOST_CHECK_EQUAL(round_(column3.at(4),2),-1.87);

  } // End of Test 1


  // Test 2
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","smooth_pdf","2"};
    string command = "tab";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);
    vector<double> column3 = getColumnFromFile(interactions.at(0),3);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),3.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),5.86);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),8.73);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),11.59);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),14.46);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),0.79);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),0.24);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),0.00);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),0.14);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),0.77);

    BOOST_CHECK_EQUAL(round_(column3.at(0),2),0.19);
    BOOST_CHECK_EQUAL(round_(column3.at(1),2),0.14);
    BOOST_CHECK_EQUAL(round_(column3.at(2),2),0.02);
    BOOST_CHECK_EQUAL(round_(column3.at(3),2),-0.13);
    BOOST_CHECK_EQUAL(round_(column3.at(4),2),-0.22);

  } // End of Test 2

  // Test 3
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","smooth_pdf","2"};
    vector<string> arguments3{"set","smooth_pot","1"};
    string command = "tab";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,arguments3);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);
    vector<double> column3 = getColumnFromFile(interactions.at(0),3);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),3.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),5.86);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),8.73);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),11.59);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),14.46);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),0.58);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),0.36);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),0.26);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),0.33);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),0.59);

    BOOST_CHECK_EQUAL(round_(column3.at(0),2),0.08);
    BOOST_CHECK_EQUAL(round_(column3.at(1),2),0.06);
    BOOST_CHECK_EQUAL(round_(column3.at(2),2),0.01);
    BOOST_CHECK_EQUAL(round_(column3.at(3),2),-0.06);
    BOOST_CHECK_EQUAL(round_(column3.at(4),2),-0.09);

  } // End of Test 3

  // Test 4
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","periodic","1"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),3.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),5.86);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),8.73);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),11.59);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),14.46);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),0.03);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),0.06);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),0.12);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),0.11);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),0.03);

  } // End of Test 4

  // Test 5
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","periodic","1"};
    vector<string> arguments3{"set","normalize","0"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,arguments3);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),3.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),5.86);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),8.73);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),11.59);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),14.46);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),404.0);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),848.0);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),1772.0);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),1536.0);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),404.0);

  } // End of Test 5

  // Test 6
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","periodic","1"};
    vector<string> arguments3{"set","extend","0"};
    vector<string> arguments4{"set","auto","0"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,arguments3);
    tabulatedpotential.Command(bonded_statistics,command,arguments4);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),0.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),0.25);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),0.50);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),0.75);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),1.00);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),1508.0);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),1164.0);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),436.0);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),1452.0);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),1508.0);

  } // End of Test 6

  // Test 7
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","scale","bond"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),0.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),0.25);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),0.50);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),0.75);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),1.00);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),19372.0);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),18624.0);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),1744.0);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),2581.33);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),19372.0);

  } // End of Test 7

  // Test 8
  {
    vector<string> arguments;
    // Set the table properties so that only 11 points are used, this way we do
    // not need to compare as many
    arguments.push_back("set");
    arguments.push_back("n");
    arguments.push_back("5");

    vector<string> arguments2{"set","scale","angle"};
    string command = "hist";
    tabulatedpotential.Command(bonded_statistics,command,arguments);
    tabulatedpotential.Command(bonded_statistics,command,arguments2);
    tabulatedpotential.Command(bonded_statistics,command,interactions);

    vector<double> column1 = getColumnFromFile(interactions.at(0),1);
    vector<double> column2 = getColumnFromFile(interactions.at(0),2);

    BOOST_CHECK_EQUAL(round_(column1.at(0),2),0.00);
    BOOST_CHECK_EQUAL(round_(column1.at(1),2),0.25);
    BOOST_CHECK_EQUAL(round_(column1.at(2),2),0.50);
    BOOST_CHECK_EQUAL(round_(column1.at(3),2),0.75);
    BOOST_CHECK_EQUAL(round_(column1.at(4),2),1.00);

    BOOST_CHECK_EQUAL(round_(column2.at(0),2),5593.78);
    BOOST_CHECK_EQUAL(round_(column2.at(1),2),4704.86);
    BOOST_CHECK_EQUAL(round_(column2.at(2),2),909.42);
    BOOST_CHECK_EQUAL(round_(column2.at(3),2),2130.16);
    BOOST_CHECK_EQUAL(round_(column2.at(4),2),5593.78);

  } // End of Test 8


  top.Cleanup();
}

BOOST_AUTO_TEST_SUITE_END()
