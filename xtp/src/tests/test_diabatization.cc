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

#define BOOST_TEST_MODULE diabatization_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/eigenio_matrixmarket.h"
#include "votca/xtp/erdiabatization.h"
#include "votca/xtp/logger.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(erdiabatization_test)

BOOST_AUTO_TEST_CASE(coupling_test) {
  libint2::initialize();

  // populate orbitals object
  Orbitals dimer;
  dimer.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                               "/erdiabatization/dimer.xyz");
  dimer.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                      "/erdiabatization/def2-svp.xml");
  dimer.SetupAuxBasis(std::string(XTP_TEST_DATA_FOLDER) +
                      "/erdiabatization/aux-def2-svp.xml");
  dimer.setNumberOfAlphaElectrons(44);
  dimer.setNumberOfOccupiedLevels(44);

  Eigen::VectorXd ref_MOvals = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/erdiabatization/dimer_MOvals.mm");
  Eigen::MatrixXd ref_MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/erdiabatization/dimer_MOs.mm");
  dimer.MOs().eigenvalues() = ref_MOvals;
  dimer.MOs().eigenvectors() = ref_MOs;

  dimer.setGWindices(0, 130);
  dimer.setBSEindices(0, 130);
  dimer.setTDAApprox(false);

  Eigen::VectorXd ref_singletvals =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/erdiabatization/dimer_singletE.mm");

  Eigen::MatrixXd ref_spsi = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/erdiabatization/dimer_singlet_spsi.mm");

  Eigen::MatrixXd ref_spsi2 = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/erdiabatization/dimer_singlet_spsi2.mm");

  dimer.BSESinglets().eigenvalues() = ref_singletvals;
  dimer.BSESinglets().eigenvectors() = ref_spsi;
  dimer.BSESinglets().eigenvectors2() = ref_spsi2;

  // set logger
  Logger log;
  log.setReportLevel(Log::error);
  log.setMultithreading(true);
  log.setCommonPreface("\n... ...");
  ERDiabatization ERDiabatization(dimer, dimer, &log, 1, 2, "singlet", true);
  ERDiabatization.configure();
  ERDiabatization.setUpMatrices();

  // Calculate angle
  double angle = ERDiabatization.Calculate_angle();
  double angle_ref = 0.71542472271498847;
  BOOST_CHECK_CLOSE(angle_ref, angle, 1e-4);

  // diabatic Hamiltonian
  Eigen::MatrixXd diabatic_H = ERDiabatization.Calculate_diabatic_H(angle);

  Eigen::MatrixXd diabatic_H_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/erdiabatization/Hdiab_ref.mm");

  BOOST_CHECK_CLOSE(diabatic_H_ref(0, 0), diabatic_H(0, 0), 1e-6);
  BOOST_CHECK_CLOSE(diabatic_H_ref(0, 1), diabatic_H(0, 1), 1e-6);
  BOOST_CHECK_CLOSE(diabatic_H_ref(1, 0), diabatic_H(1, 0), 1e-6);
  BOOST_CHECK_CLOSE(diabatic_H_ref(1, 1), diabatic_H(1, 1), 1e-6);

  libint2::finalize();
}
BOOST_AUTO_TEST_SUITE_END()
