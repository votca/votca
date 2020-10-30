/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#include <memory>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE lammpdatatopologyreaderwriter_test
#include <boost/test/unit_test.hpp>

#include "../../include/votca/csg/bead.h"
#include "../../include/votca/csg/orthorhombicbox.h"
#include "../../include/votca/csg/topologyreader.h"
#include "../../include/votca/csg/trajectoryreader.h"
#include "../../include/votca/csg/trajectorywriter.h"
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <votca/tools/elements.h>
#include <votca/tools/types.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(lammpsdatareaderwriter_test)

/**
 * \brief Test the topology reader
 *
 * This test is designed to test the topology reader this is done by
 * creating a small lammps data file. The file is then read in with the
 * topology reader and the values in the top object are then examined
 * to ensure they represent the values from the file.
 */
BOOST_AUTO_TEST_CASE(test_topologyreader) {

  string lammpsdatafilename =
      std::string(CSG_TEST_DATA_FOLDER) + "/lammpsdatareader/test_polymer.data";

  Topology top;

  TopologyReader::RegisterPlugins();
  std::unique_ptr<TopologyReader> lammpsDataReader =
      std::unique_ptr<TopologyReader>(
          TopReaderFactory().Create(lammpsdatafilename));
  lammpsDataReader->ReadTopology(lammpsdatafilename, top);

  BOOST_CHECK_EQUAL(top.BeadCount(), 100);
  Eigen::Vector3d first_bead_correct_pos(62.806, 52.5127, 49.8873);
  Bead *firstBead = top.getBead(0);
  auto first_bead_pos = firstBead->getPos();
  BOOST_CHECK(first_bead_correct_pos.isApprox(first_bead_pos, 1e-3));

  Eigen::Vector3d last_bead_correct_pos(102.78495, 78.0388, 59.9629);
  Bead *lastBead = top.getBead(99);
  auto last_bead_pos = lastBead->getPos();
  BOOST_CHECK(last_bead_correct_pos.isApprox(last_bead_pos, 1e-3));

  auto mol = top.getMolecule(0);

  BOOST_CHECK_EQUAL(mol->getName(), "N100");
  BOOST_CHECK_EQUAL(mol->BeadCount(), 100);

  auto interaction_cont = top.BondedInteractions();
  votca::Index numBondInter = 99;
  votca::Index numAngleInter = 98;
  votca::Index numDihedralInter = 97;
  votca::Index totalInter = numBondInter + numAngleInter + numDihedralInter;
  BOOST_CHECK_EQUAL(interaction_cont.size(), totalInter);
}

/**
 * \brief Test the trajectory reader
 *
 * This test is designed to test the trajectory reader this is done by
 * creating a small lammps data file. A topology object is created and
 * it is filled by first reading information from a single data file
 * using the topology reader. A trajectory reader is then used to read
 * a second data file. This leads to updating the positions of the atoms
 */

BOOST_AUTO_TEST_CASE(test_trajectoryreader) {

  string lammpsdatafilename = std::string(CSG_TEST_DATA_FOLDER) +
                              "/lammpsdatareader/test_polymer3.data";
  Topology top;

  TopologyReader::RegisterPlugins();
  std::unique_ptr<TopologyReader> lammpsDataReader =
      std::unique_ptr<TopologyReader>(
          TopReaderFactory().Create(lammpsdatafilename));
  lammpsDataReader->ReadTopology(lammpsdatafilename, top);

  string lammpsdatafilename2 = std::string(CSG_TEST_DATA_FOLDER) +
                               "/lammpsdatareader/test_polymer4.data";

  TrajectoryReader::RegisterPlugins();
  std::unique_ptr<TrajectoryReader> lammpsDataReaderTrj =
      std::unique_ptr<TrajectoryReader>(
          TrjReaderFactory().Create(lammpsdatafilename2));

  lammpsDataReaderTrj->Open(lammpsdatafilename2);
  lammpsDataReaderTrj->FirstFrame(top);
  lammpsDataReaderTrj->Close();

  BOOST_CHECK_EQUAL(top.BeadCount(), 100);

  Eigen::Vector3d first_bead_correct_pos(65.7991, 51.04235, 58.480193);
  Bead *firstBead = top.getBead(0);
  auto first_bead_pos = firstBead->getPos();

  cout << first_bead_correct_pos << endl;
  cout << first_bead_pos << endl;

  BOOST_CHECK(first_bead_correct_pos.isApprox(first_bead_pos, 1e-3));
  Eigen::Vector3d last_bead_correct_pos(108.431, 83.94695, 68.5254);
  Bead *lastBead = top.getBead(99);
  auto last_bead_pos = lastBead->getPos();

  cout << last_bead_correct_pos << endl;
  cout << last_bead_pos << endl;
  BOOST_CHECK(last_bead_correct_pos.isApprox(last_bead_pos, 1e-3));

  auto mol = top.getMolecule(0);
  BOOST_CHECK_EQUAL(mol->BeadCount(), 100);

  auto interaction_cont = top.BondedInteractions();
  votca::Index numBondInter = 99;
  votca::Index numAngleInter = 98;
  votca::Index numDihedralInter = 97;
  votca::Index totalInter = numBondInter + numAngleInter + numDihedralInter;
  BOOST_CHECK_EQUAL(interaction_cont.size(), totalInter);
}

BOOST_AUTO_TEST_SUITE_END()
