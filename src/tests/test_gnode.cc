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

#define BOOST_TEST_MODULE gnode_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <votca/xtp/glink.h>
#include <votca/xtp/gnode.h>

using namespace std;
using namespace votca::xtp;
BOOST_AUTO_TEST_SUITE(gnode_test)

BOOST_AUTO_TEST_CASE(chosen_id_test) {

  QMStateType electron = QMStateType::Electron;

  vector<GNode> dests;
  for (int i = 0; i < 6; i++) {
    Segment seg("one", i);
    dests.push_back(GNode(seg, electron, true));
  }
  Segment seg("one", 6);
  GNode g(seg, electron, true);
  g.AddEvent(&dests[0], Eigen::Vector3d::Zero(), 10);
  g.AddEvent(&dests[1], Eigen::Vector3d::Zero(), 20);
  g.AddEvent(&dests[2], Eigen::Vector3d::Zero(), 15);
  g.AddEvent(&dests[3], Eigen::Vector3d::Zero(), 18);
  g.AddEvent(&dests[4], Eigen::Vector3d::Zero(), 12);
  g.AddEvent(&dests[5], Eigen::Vector3d::Zero(), 25);
  g.InitEscapeRate();
  g.MakeHuffTree();
  std::cout << g.findHoppingDestination(0.55)->getDestination()->getId()
            << std::endl;
  std::cout << g.findHoppingDestination(0.85)->getDestination()->getId()
            << std::endl;
  std::cout << g.findHoppingDestination(0.25)->getDestination()->getId()
            << std::endl;
  std::cout << g.findHoppingDestination(0.15)->getDestination()->getId()
            << std::endl;
  std::cout << g.findHoppingDestination(0.35)->getDestination()->getId()
            << std::endl;
  std::cout << g.findHoppingDestination(0.65)->getDestination()->getId()
            << std::endl;
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.55)->getDestination()->getId(),
                    0);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.85)->getDestination()->getId(),
                    1);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.25)->getDestination()->getId(),
                    2);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.15)->getDestination()->getId(),
                    3);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.35)->getDestination()->getId(),
                    4);
  BOOST_CHECK_EQUAL(g.findHoppingDestination(0.65)->getDestination()->getId(),
                    5);
}

BOOST_AUTO_TEST_CASE(count_test) {
  QMStateType electron = QMStateType::Electron;

  vector<GNode> dests;
  for (int i = 0; i < 11; i++) {
    Segment seg("one", i);
    dests.push_back(GNode(seg, electron, true));
  }
  Segment seg("one", 12);
  GNode g(seg, electron, true);

  g.AddEvent(&dests[0], Eigen::Vector3d::Zero(), 15);
  g.AddEvent(&dests[1], Eigen::Vector3d::Zero(), 9);
  g.AddEvent(&dests[2], Eigen::Vector3d::Zero(), 11);
  g.AddEvent(&dests[3], Eigen::Vector3d::Zero(), 8);
  g.AddEvent(&dests[4], Eigen::Vector3d::Zero(), 12);
  g.AddEvent(&dests[5], Eigen::Vector3d::Zero(), 7);
  g.AddEvent(&dests[6], Eigen::Vector3d::Zero(), 13);
  g.AddEvent(&dests[7], Eigen::Vector3d::Zero(), 6);
  g.AddEvent(&dests[8], Eigen::Vector3d::Zero(), 14);
  g.AddEvent(&dests[9], Eigen::Vector3d::Zero(), 5);
  g.AddEvent(&dests[10], Eigen::Vector3d::Zero(), 100);

  g.InitEscapeRate();
  g.MakeHuffTree();
  vector<int> count(11, 0);
  double d = 0;
  while (d < 1) {
    GLink* L = g.findHoppingDestination(d);
    int ind = L->getDestination()->getId();
    count[ind]++;
    d += 0.000001;
  }

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
