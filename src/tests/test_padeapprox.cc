/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE padeapprox_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <fstream>
#include <votca/xtp/padeapprox.h>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <votca/xtp/orbitals.h>
using namespace votca::xtp;
using namespace std; 

BOOST_AUTO_TEST_SUITE(padeapprox_test)

BOOST_AUTO_TEST_CASE(padeapprox_test){

    
  std::cout << "Started Pade Test" << std::endl;
 
  std::vector<std::complex<double>> grid;
  
  std::complex<double> im(0,1);
  
  int num_points=200;

  PadeApprox pade_1;
  pade_1.initialize(num_points);
  PadeApprox pade_2;
  pade_2.initialize(num_points);
  PadeApprox pade_3;
  pade_3.initialize(num_points);
  PadeApprox pade_4;
  pade_4.initialize(num_points);
  PadeApprox pade_5;
  pade_5.initialize(num_points);
  PadeApprox pade_6;
  pade_6.initialize(num_points);
  PadeApprox pade_7;
  pade_7.initialize(num_points);
  PadeApprox pade_8;
  pade_8.initialize(num_points);
  PadeApprox pade_9;
  pade_9.initialize(num_points);
  
  for(int i=0;i<num_points;i++){
      grid.push_back(i);
  }
  
  for (int i = 0; i < grid.size(); i++) {
    std::cout << "Gridpoint: " << grid.at(i) << std::endl;
  }
  
  std::vector<Eigen::Matrix3d> val;
  
  for(int i=1;i<num_points+1;i++){
      val.push_back(Eigen::Matrix3d::Random());
  }
  
  for (int i = 0; i < val.size(); i++) {
    std::cout << "Values: " << val.at(i)(0,0) << std::endl;
  }
  
  std::cout << "Init done"<< std::endl;
  for (int j = 0; j < grid.size(); j++) {
    //std::cout << "loop "<<j<< std::endl; 
    pade_1.addPoint(grid.at(j), val.at(j)(0,0));
    pade_2.addPoint(grid.at(j), val.at(j)(0,1));
    pade_3.addPoint(grid.at(j), val.at(j)(0,2));
    pade_4.addPoint(grid.at(j), val.at(j)(1,0));
    pade_5.addPoint(grid.at(j), val.at(j)(1,1));
    pade_6.addPoint(grid.at(j), val.at(j)(1,2));
    pade_7.addPoint(grid.at(j), val.at(j)(2,0));
    pade_8.addPoint(grid.at(j), val.at(j)(2,1));
    pade_9.addPoint(grid.at(j), val.at(j)(2,2));
  }
  std::vector<Eigen::Matrix3cd> results;
  for (std::complex<double> w:grid){
    results.push_back(Eigen::Matrix3d::Zero());
    results[results.size()-1](0,0)=pade_1.evaluatePoint(w);
    results[results.size()-1](0,1)=pade_2.evaluatePoint(w);
    results[results.size()-1](0,2)=pade_3.evaluatePoint(w);
    results[results.size()-1](1,0)=pade_4.evaluatePoint(w);
    results[results.size()-1](1,1)=pade_5.evaluatePoint(w);
    results[results.size()-1](1,2)=pade_6.evaluatePoint(w);
    results[results.size()-1](2,0)=pade_7.evaluatePoint(w);
    results[results.size()-1](2,1)=pade_8.evaluatePoint(w);
    results[results.size()-1](2,2)=pade_9.evaluatePoint(w);
  }
  
  for(int i=0;i<num_points;i++){
      std::cout<<"at point "<<grid[i]<<" diff="<<results[i](0,0)-val[i](0,0)<<std::endl;
  }
  bool test1 = results[0].isApprox(val.at(0),1e-9);
  bool test2 = results[2].isApprox(val.at(2),1e-9);
  bool test3 = results[4].isApprox(val.at(4),1e-9);
  bool test4 = results[6].isApprox(val.at(6),1e-9);
  bool test5 = results[results.size()-1].isApprox(val.at(num_points-1),1e-9);
  BOOST_CHECK_EQUAL(test1, true);
  BOOST_CHECK_EQUAL(test2, true);
  BOOST_CHECK_EQUAL(test3, true);
  BOOST_CHECK_EQUAL(test4, true);
  BOOST_CHECK_EQUAL(test5, true);
}BOOST_AUTO_TEST_SUITE_END()