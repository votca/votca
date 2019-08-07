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

#define BOOST_TEST_MODULE grid_test
#include <boost/test/unit_test.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <votca/xtp/grid.h>

using namespace votca::xtp;
BOOST_AUTO_TEST_SUITE(grid_test)

BOOST_AUTO_TEST_CASE(setupgrid_test) {
  QMMolecule qm = QMMolecule("bla", 1);
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  qm.push_back(QMAtom(0, "C", pos));

  Grid grid;
  grid.setupCHELPGGrid(qm);
  // BOOST_CHECK_EQUAL(grid.size(),2910);Thischeck is extremly sensitive to
  // numerical precision e.g. 1e-18 so it is not a good test

  Eigen::Vector3d start(-5.1022601692, -1.1338355932, -0.56691779658);
  Eigen::Vector3d end(5.10226, 1.13384, 0.566918);
  bool start_check = start.isApprox(grid.getGridPositions()[0], 0.00001);
  bool end_check = end.isApprox(grid.getGridPositions().back(), 0.00001);

  if (!start_check) {
    std::cout << "ref" << std::endl;
    std::cout << start << std::endl;
    std::cout << "result" << std::endl;
    std::cout << grid.getGridPositions()[0] << std::endl;
  }
  BOOST_CHECK_EQUAL(start_check, true);
  if (!end_check) {
    std::cout << "ref" << std::endl;
    std::cout << end << std::endl;
    std::cout << "result" << std::endl;
    std::cout << grid.getGridPositions().back() << std::endl;
  }
  BOOST_CHECK_EQUAL(end_check, true);
}
BOOST_AUTO_TEST_SUITE_END()
