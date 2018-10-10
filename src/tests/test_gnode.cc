/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE gnode_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/gnode.h>
#include <votca/xtp/glink.h>
#include <iostream>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

BOOST_AUTO_TEST_SUITE(gnode_test)

BOOST_AUTO_TEST_CASE(deterministic_test) {
  votca::xtp::GNode g;
  g.events=std::vector<votca::xtp::GLink>(6);
  g.events[0].destination=0;
  g.events[0].rate =10;
  g.events[1].destination=1;  
  g.events[1].rate =20;
  g.events[2].destination=2;
  g.events[2].rate =15;
  g.events[3].destination=3;
  g.events[3].rate =18;
  g.events[4].destination=4;
  g.events[4].rate =12;
  g.events[5].destination=5;
  g.events[5].rate =25;
  g.escape_rate=100;  
  g.MakeHuffTree();
  std::cout<<g.findHoppingDestination(0.55)->destination<<std::endl;
  std::cout<<g.findHoppingDestination(0.85)->destination<<std::endl;
  std::cout<<g.findHoppingDestination(0.25)->destination<<std::endl;
  std::cout<<g.findHoppingDestination(0.15)->destination<<std::endl;
  std::cout<<g.findHoppingDestination(0.35)->destination<<std::endl;
  std::cout<<g.findHoppingDestination(0.65)->destination<<std::endl;
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.55)->destination,0);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.85)->destination,1);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.25)->destination,2);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.15)->destination,3);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.35)->destination,4);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.65)->destination,5);
}

BOOST_AUTO_TEST_CASE(statistical_test) {
      srand( (unsigned)time( NULL ) );
  votca::xtp::GNode g;
  g.events=std::vector<votca::xtp::GLink>(6);
  g.events[0].destination=0;
  g.events[0].rate =10;
  g.events[1].destination=1;
  g.events[1].rate =20;
  g.events[2].destination=2;
  g.events[2].rate =15;
  g.events[3].destination=3;
  g.events[3].rate =22;
  g.events[4].destination=4;
  g.events[4].rate =8;
  g.events[5].destination=5;
  g.events[5].rate =25;
  g.MakeHuffTree();  
  vector<int> count(6);
  srand( (unsigned)time( NULL ) );
  double d;

    for (int i=0;i<6;i++){
        count[i]=0;
    }
    int ind;
    long max=1e6;
    for (long l=0;l<max;l++){
       d=(double)rand()/RAND_MAX;
       votca::xtp::GLink * L=g.findHoppingDestination(d);
       ind=L->destination;
       count[ind]++;
    }
  BOOST_CHECK_EQUAL(abs(count[0]-100000)<5000,true);
  BOOST_CHECK_EQUAL(abs(count[1]-200000)<5000,true);
  BOOST_CHECK_EQUAL(abs(count[2]-150000)<5000,true);
  BOOST_CHECK_EQUAL(abs(count[3]-220000)<5000,true);
  BOOST_CHECK_EQUAL(abs(count[4]-80000)<5000,true);
  BOOST_CHECK_EQUAL(abs(count[5]-250000)<5000,true);
}
BOOST_AUTO_TEST_SUITE_END()
