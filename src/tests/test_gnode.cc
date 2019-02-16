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
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <votca/xtp/glink.h>
#include <votca/xtp/gnode.h>

BOOST_AUTO_TEST_SUITE(gnode_test)

BOOST_AUTO_TEST_CASE(chosen_id_test) {
  votca::xtp::GNode g;
  g.events = std::vector<votca::xtp::GLink>(6);
  g.events[0].destination = 0;
  g.events[0].rate = 10;
  g.events[1].destination = 1;
  g.events[1].rate = 20;
  g.events[2].destination = 2;
  g.events[2].rate = 15;
  g.events[3].destination = 3;
  g.events[3].rate = 18;
  g.events[4].destination = 4;
  g.events[4].rate = 12;
  g.events[5].destination = 5;
  g.events[5].rate = 25;
  g.escape_rate = 100;
  g.MakeHuffTree();
  std::cout << g.findHoppingDestination(0.55)->destination << std::endl;
  std::cout << g.findHoppingDestination(0.85)->destination << std::endl;
  std::cout << g.findHoppingDestination(0.25)->destination << std::endl;
  std::cout << g.findHoppingDestination(0.15)->destination << std::endl;
  std::cout << g.findHoppingDestination(0.35)->destination << std::endl;
  std::cout << g.findHoppingDestination(0.65)->destination << std::endl;
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.55)->destination, 0);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.85)->destination, 1);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.25)->destination, 2);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.15)->destination, 3);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.35)->destination, 4);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.65)->destination, 5);
}

BOOST_AUTO_TEST_CASE(count_test) {
  votca::xtp::GNode g;
  g.events = std::vector<votca::xtp::GLink>(11);
  g.events[0].destination = 0;
  g.events[0].rate = 15;
  g.events[1].destination = 1;
  g.events[1].rate = 9;
  g.events[2].destination = 2;
  g.events[2].rate = 11;
  g.events[3].destination = 3;
  g.events[3].rate = 8;
  g.events[4].destination = 4;
  g.events[4].rate = 12;
  g.events[5].destination = 5;
  g.events[5].rate = 7;
  g.events[6].destination = 6;
  g.events[6].rate = 13;
  g.events[7].destination = 7;
  g.events[7].rate = 6;
  g.events[8].destination = 8;
  g.events[8].rate = 14;
  g.events[9].destination = 9;
  g.events[9].rate = 5;
  g.events[10].destination = 10;
  g.events[10].rate = 100;
  g.MakeHuffTree();
  vector<int> count(11);
  double d = 0;
  for (int i = 0; i < 11; i++) {
    count[i] = 0;
  }
  int ind;
  while (d < 1) {
    votca::xtp::GLink* L = g.findHoppingDestination(d);
    ind = L->destination;
    count[ind]++;
    d += 0.000001;
  }
  std::cout << count[0] << endl;
  std::cout << count[1] << endl;
  std::cout << count[2] << endl;
  std::cout << count[3] << endl;
  std::cout << count[4] << endl;
  std::cout << count[5] << endl;
  std::cout << count[6] << endl;
  std::cout << count[7] << endl;
  std::cout << count[8] << endl;
  std::cout << count[9] << endl;
  std::cout << count[10] << endl;
  BOOST_CHECK_EQUAL(count[0], 75000);
  BOOST_CHECK_EQUAL(count[1], 45000);
  BOOST_CHECK_EQUAL(count[2], 55000);
  BOOST_CHECK_EQUAL(count[3], 40000);
  BOOST_CHECK_EQUAL(count[4], 60000);
  BOOST_CHECK_EQUAL(count[5], 35000);
  BOOST_CHECK_EQUAL(count[6], 65000);
  BOOST_CHECK_EQUAL(count[7], 30000);
  BOOST_CHECK_EQUAL(count[8], 70000);
  BOOST_CHECK_EQUAL(count[9], 25001);
  BOOST_CHECK_EQUAL(count[10], 499999);
}
BOOST_AUTO_TEST_SUITE_END()
