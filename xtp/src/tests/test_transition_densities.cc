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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE diabatization_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/transition_densities.h"
#include <votca/tools/types.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(transition_densities_test)

BOOST_AUTO_TEST_CASE(matrix_test) {
  libint2::initialize();

  // populate orbitals object
  Orbitals dimer;
  dimer.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                               "/diabatization/dimer.xyz");
  dimer.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                      "/diabatization/def2-svp.xml");
  dimer.SetupAuxBasis(std::string(XTP_TEST_DATA_FOLDER) +
                      "/diabatization/aux-def2-svp.xml");
  dimer.setNumberOfAlphaElectrons(44);
  dimer.setNumberOfOccupiedLevels(44);

  Eigen::VectorXd ref_MOvals = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/diabatization/dimer_MOvals.mm");
  Eigen::MatrixXd ref_MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/diabatization/dimer_MOs.mm");
  dimer.MOs().eigenvalues() = ref_MOvals;
  dimer.MOs().eigenvectors() = ref_MOs;

  dimer.setGWindices(0, 130);
  dimer.setBSEindices(0, 130);
  dimer.setTDAApprox(false);

  Eigen::VectorXd ref_singletvals =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/diabatization/dimer_singletE.mm");

  Eigen::MatrixXd ref_spsi = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/diabatization/dimer_singlet_spsi.mm");

  Eigen::MatrixXd ref_spsi2 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/diabatization/dimer_singlet_spsi2.mm");

  dimer.BSESinglets().eigenvalues() = ref_singletvals;
  dimer.BSESinglets().eigenvectors() = ref_spsi;
  dimer.BSESinglets().eigenvectors2() = ref_spsi2;

  // set logger
  Logger log;
  log.setReportLevel(Log::error);
  log.setMultithreading(true);
  log.setCommonPreface("\n... ...");

  TransitionDensities tdmat(dimer, dimer, &log);
  tdmat.configure();
  QMStateType qmtype;
  qmtype.FromString("singlet");
  QMState state1 = QMState(qmtype, 0, false);
  QMState state2 = QMState(qmtype, 1, false);

  Eigen::MatrixXd tmat = tdmat.Matrix(state1, state2);

  Eigen::MatrixXd tmat_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/transition_densities/matrix_ref.mm");

  bool check_tmat = tmat.isApprox(tmat_ref, 0.00001);
  BOOST_CHECK_EQUAL(check_tmat, 1);
  if (!check_tmat) {
    std::cout << "ref" << std::endl;
    std::cout << tmat_ref << std::endl;
    std::cout << "result" << std::endl;
    std::cout << tmat << std::endl;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
