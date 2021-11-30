/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
#include <libint2/initialize.h>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE background_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/background.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(background_test)

BOOST_AUTO_TEST_CASE(all_test) {
  libint2::initialize();

  // Create a fake polar background
  std::vector<PolarSegment> backgroundSegments;
  PolarSegment seg0("segment0", 0);
  PolarSegment seg1("segment1", 1);
  seg0.push_back(PolarSite(0, "H", Eigen::Vector3d(1, 1, 1)));
  seg0.push_back(PolarSite(1, "H", Eigen::Vector3d(1, 1, 1.5)));
  // seg0 has pos 1,1,1.25
  seg1.push_back(PolarSite(0, "H", Eigen::Vector3d(2.5, 2.5, 2.5)));
  seg1.push_back(PolarSite(1, "H", Eigen::Vector3d(2.5, 2.5, 3.0)));
  // seg1 has pos 2.5, 2.5, 2.75
  backgroundSegments.push_back(seg0);
  backgroundSegments.push_back(seg1);

  // Other necessary input
  Logger log;
  EwaldOptions options;
  options.alpha = 1;
  options.k_cutoff = 1;
  options.r_cutoff = 10;
  options.shape = Shape::cube;
  options.sharpness = 0.39;

  // Create background
  Background bg(log, 5 * Eigen::Matrix3d::Identity(), options,
                backgroundSegments);

  BGNbList nbList =  bg.getNbList();
  // Check the first two neighbours of segment 0
  BOOST_CHECK(nbList.getNeighboursOf(0)[0].getId() == 1);
  BOOST_CHECK(nbList.getNeighboursOf(0)[1].getId() == 1);
  BOOST_CHECK_CLOSE(nbList.getNeighboursOf(0)[0].getDist(), 2.5981, 1e-3);
  BOOST_CHECK_CLOSE(nbList.getNeighboursOf(0)[1].getDist(), 4.0927, 1e-3);

  //
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
