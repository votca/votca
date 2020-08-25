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
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE statetracker_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/statetracker.h"

using namespace votca::xtp;
using namespace std;
BOOST_AUTO_TEST_SUITE(statetracker_test)
BOOST_AUTO_TEST_CASE(osc) {

  Orbitals orb;
  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/statetracker/molecule.xyz");
  orb.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                      "/statetracker/3-21G.xml");
  Logger log;
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(4);
  orb.setBSEindices(0, 16);
  orb.setTDAApprox(true);

  Eigen::MatrixXd& MOs = orb.MOs().eigenvectors();

  orb.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/statetracker/MOs.mm");

  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;
  orb.BSESinglets().eigenvalues() = se_ref;

  // reference coefficients
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/statetracker/spsi_ref.mm");

  orb.BSESinglets().eigenvectors() = spsi_ref;
  orb.CalcCoupledTransition_Dipoles();

  {
    votca::tools::Property prop;
    prop.LoadFromXML(std::string(XTP_TEST_DATA_FOLDER) +
                     "/statetracker/statetracker.xml");
    StateTracker tracker;
    tracker.setLogger(&log);
    QMState s("s1");
    tracker.setInitialState(s);
    tracker.Initialize(prop.get("statetracker"));
    QMState newstate = tracker.CalcState(orb);
    BOOST_CHECK_EQUAL(newstate.Type().ToString(), "s");
    BOOST_CHECK_EQUAL(newstate.StateIdx(), 1);
  }
}

BOOST_AUTO_TEST_CASE(readwrite_hdf5) {
  votca::tools::Property prop;
  prop.LoadFromXML(std::string(XTP_TEST_DATA_FOLDER) +
                   "/statetracker/statetracker2.xml");
  StateTracker tracker;
  QMState s("s1");
  tracker.setInitialState(s);
  tracker.Initialize(prop.get("statetracker"));

  {
    CheckpointFile f("statetracker_test.hdf5");
    CheckpointWriter w = f.getWriter();
    tracker.WriteToCpt(w);
  }
  StateTracker tracker2;
  CheckpointFile f("statetracker_test.hdf5");
  CheckpointReader r = f.getReader();
  tracker2.ReadFromCpt(r);

  BOOST_CHECK_EQUAL(tracker.InitialState().ToString(),
                    tracker2.InitialState().ToString());
}

BOOST_AUTO_TEST_SUITE_END()
