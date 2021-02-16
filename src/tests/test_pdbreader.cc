/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <cmath>
#include <fstream>
#include <string>

// Third party includes
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/csg/bead.h"
#include "votca/csg/topologyreader.h"

using namespace std;
using namespace votca::csg;

BOOST_AUTO_TEST_SUITE(pdbreader_test)

BOOST_AUTO_TEST_CASE(test_topologyreader) {

  votca::tools::Elements ele;

  Topology top;
  TopologyReader::RegisterPlugins();
  string str = std::string(CSG_TEST_DATA_FOLDER) + "/pdbreader/Molecule1.pdb";
  auto reader = std::unique_ptr<TopologyReader>(TopReaderFactory().Create(str));
  reader->ReadTopology(str, top);
  BOOST_CHECK_EQUAL(reader != nullptr, true);
  BOOST_CHECK_EQUAL(top.BeadCount(), 10);

  vector<votca::Index> resnr = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1};
  vector<string> bd_name = {"C", "H", "H", "H", "H", "C", "H", "H", "H", "H"};
  vector<double> x = {-0.5249, -0.6202, -0.539,  -0.4682, -0.4724,
                      -0.2248, -0.1518, -0.3153, -0.2442, -0.1880};
  vector<double> y = {0.1055, 0.1521, 0.0026, 0.1124, 0.1550,
                      0.1671, 0.2451, 0.1999, 0.1430, 0.0804};
  vector<double> z = {-0.000, -0.0141, 0.0255, -0.0904, 0.079,
                      -0.000, 0.0051,  0.0467, -0.1024, 0.0507};
  Bead *bd;
  Eigen::Vector3d v;
  for (votca::Index i = 0; i < 10; i++) {
    bd = top.getBead(i);
    BOOST_CHECK_EQUAL(bd->getId(), i);
    BOOST_CHECK_EQUAL(bd->getResnr(), resnr.at(i));
    BOOST_CHECK_EQUAL(bd->getName(), bd_name.at(i));
    // BOOST_CHECK_EQUAL(bd->getMass(),ele.getMass(bd->getName()));
    v = bd->getPos();
    BOOST_CHECK_CLOSE(bd->getQ(), 0, 1e-5);
    BOOST_CHECK_CLOSE(v.x(), x.at(i), 1e-5);
    BOOST_CHECK_CLOSE(v.y(), y.at(i), 1e-5);
    BOOST_CHECK_CLOSE(v.z(), z.at(i), 1e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
