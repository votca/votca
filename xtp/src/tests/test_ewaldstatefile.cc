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

#define BOOST_TEST_MODULE ewald_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/background.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/orbitals.h"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using votca::Index;
BOOST_AUTO_TEST_SUITE(ewald_calculations)

BOOST_AUTO_TEST_CASE(writereadcycle_statefile) {
  libint2::initialize();

  // Create a fake polar background
  std::vector<PolarSegment> backgroundSegments;
  PolarSegment seg0("segment0", 0);
  PolarSegment seg1("segment1", 1);
  seg0.push_back(PolarSite(0, "H", Eigen::Vector3d(1, 1, 1)));
  seg0.push_back(PolarSite(1, "H", Eigen::Vector3d(1, 1, 1.5)));
  seg1.push_back(PolarSite(0, "H", Eigen::Vector3d(3, 3, 3)));
  seg1.push_back(PolarSite(1, "H", Eigen::Vector3d(2.7, 3, 3.7)));
  backgroundSegments.push_back(seg0);
  backgroundSegments.push_back(seg1);

  // Other necessary input
  Logger log;
  EwaldOptions options;
  options.alpha = 1;
  options.k_cutoff = 1;
  options.r_cutoff = 1;
  options.shape = Shape::cube;
  options.sharpness = 0.39;

  // Create background
  Background bg(log, 5 * Eigen::Matrix3d::Identity(), options,
                backgroundSegments);

  // Write to state file
  bg.writeToStateFile("testFile.hdf5");

  // Read From state file
  Background bg2(log);
  bg2.readFromStateFile("testFile.hdf5");
  // Compare
  bool background_is_equal = (bg == bg2);
  BOOST_CHECK(background_is_equal);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(energy_calc_NaCl){
  // Index cry_l = 16;
  // Index N = cry_l * cry_l * cry_l;
  // double l = (double)cry_l / 2.0;
  // std::vector<PolarSegment> backgroundSegments;
  // for (Index i = 0; i < N; ++i) {
  //   PolarSegment seg("seg", i);
  //   // note element type doesn't matter for the calcul
  //   seg.push_back(PolarSite(0, "H", ));

  //   Index ix = i % cry_l;
  //   Index iy = (i % (cry_l * cry_l)) / cry_l;
  //   Index iz = i / (cry_l * cry_l);

  //   xyz[i][0] = (double)(ix)*0.5;
  //   xyz[i][1] = (double)(iy)*0.5;
  //   xyz[i][2] = (double)(iz)*0.5;

  //   charges[i] = (((ix + iy + iz) % 2) ? 1.0 : -1.0);
  // }
  // EwaldOptions options;
  // options.alpha = 4;
  // options.r_cutoff = 1.9;
  // options.k_cutoff = 30;
  // options.shape = Shape::cube;
  // options.sharpness = 0.39;
  // Eigen::Matrix3d cellMatrix;
  // cellMatrix << l, 0, 0, 0, l, 0, 0, 0, l;

  // Logger log;
  // Background bg(log, cellMatrix, options, backgroundSegments);

  // double energy = bg.interactionEnergy();




}

BOOST_AUTO_TEST_CASE(energy_eval_with_induced_dipoles){
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(field_eval_with_induced_dipoles){
BOOST_CHECK(true);
}

BOOST_AUTO_TEST_CASE(background_polarization){
  BOOST_CHECK(true);
}

BOOST_AUTO_TEST_SUITE_END()
