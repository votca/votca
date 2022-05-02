
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

#define BOOST_TEST_MODULE lammpdatatopologyreaderwriter_test

// Standard includes
#include <cmath>
#include <cstdio>
#include <fstream>
#include <memory>
#include <string>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/elements.h>
#include <votca/tools/types.h>

// Local VOTCA includes
#include "votca/csg/bead.h"
#include "votca/csg/orthorhombicbox.h"
#include "votca/csg/topologyreader.h"
#include "votca/csg/trajectoryreader.h"
#include "votca/csg/trajectorywriter.h"

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
      TopReaderFactory().Create(lammpsdatafilename);
  lammpsDataReader->ReadTopology(lammpsdatafilename, top);

  BOOST_CHECK_EQUAL(top.BeadCount(), 100);
  Eigen::Vector3d first_bead_correct_pos(6.2806, 5.25127, 4.98873);
  Bead *firstBead = top.getBead(0);
  auto first_bead_pos = firstBead->getPos();
  BOOST_CHECK(first_bead_correct_pos.isApprox(first_bead_pos, 1e-3));

  Eigen::Vector3d last_bead_correct_pos(10.278495, 7.80388, 5.99629);
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
      TopReaderFactory().Create(lammpsdatafilename);
  lammpsDataReader->ReadTopology(lammpsdatafilename, top);

  string lammpsdatafilename2 = std::string(CSG_TEST_DATA_FOLDER) +
                               "/lammpsdatareader/test_polymer4.data";

  TrajectoryReader::RegisterPlugins();
  std::unique_ptr<TrajectoryReader> lammpsDataReaderTrj =
      TrjReaderFactory().Create(lammpsdatafilename2);

  lammpsDataReaderTrj->Open(lammpsdatafilename2);
  lammpsDataReaderTrj->FirstFrame(top);
  lammpsDataReaderTrj->Close();

  BOOST_CHECK_EQUAL(top.BeadCount(), 100);

  Eigen::Vector3d first_bead_correct_pos(6.57991, 5.104235, 5.8480193);
  Bead *firstBead = top.getBead(0);
  auto first_bead_pos = firstBead->getPos();

  cout << first_bead_correct_pos << endl;
  cout << first_bead_pos << endl;

  BOOST_CHECK(first_bead_correct_pos.isApprox(first_bead_pos, 1e-3));
  Eigen::Vector3d last_bead_correct_pos(10.8431, 8.394695, 6.85254);
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

BOOST_AUTO_TEST_CASE(test_molecules) {

  string lammpsdatafilename = std::string(CSG_TEST_DATA_FOLDER) +
                              "/lammpsdatareader/test_bondedmolecules.data";
  Topology top;

  TopologyReader::RegisterPlugins();
  std::unique_ptr<TopologyReader> lammpsDataReader =
      TopReaderFactory().Create(lammpsdatafilename);
  BOOST_REQUIRE_THROW(lammpsDataReader->ReadTopology(lammpsdatafilename, top),
                      std::runtime_error);
  string lammpsdatafilename2 = std::string(CSG_TEST_DATA_FOLDER) +
                               "/lammpsdatareader/test_twomolecules.data";

  Topology top2;
  std::unique_ptr<TopologyReader> lammpsDataReader2 =
      TopReaderFactory().Create(lammpsdatafilename2);

  lammpsDataReader2->ReadTopology(lammpsdatafilename2, top2);
  std::vector<std::string> refnames{"H2O1-1", "H2O1-0"};
  // the first molecule has order OHH while the second has HHO so the second is
  // lexicographically first
  votca::Index i = 0;
  for (const auto &mol : top2.Molecules()) {
    BOOST_CHECK_EQUAL(mol.getName(), refnames[i]);
    i++;
  }
}

BOOST_AUTO_TEST_SUITE_END()
