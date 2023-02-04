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
#include "votca/xtp/fcddiabatization.h"
#include "votca/xtp/gmhdiabatization.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmfragment.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(diabatization_test)

BOOST_AUTO_TEST_CASE(ER_coupling_test) {
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
  ERDiabatization ERDiabatization(dimer, dimer, &log, 1, 2, "singlet", true);
  ERDiabatization.configure();
  ERDiabatization.setUpMatrices();

  // Calculate angle
  double angle = ERDiabatization.Calculate_angle();
  double angle_ref = 0.7063121313715891;
  BOOST_CHECK_CLOSE_FRACTION(angle_ref, angle, 5e-3);

  // diabatic Hamiltonian
  std::pair<Eigen::VectorXd, Eigen::MatrixXd> rotated_H =
      ERDiabatization.Calculate_diabatic_H(angle);

  double E1_ref = 0.21995615824390968;
  double E2_ref = 0.22090555616542826;
  BOOST_CHECK_CLOSE(E1_ref, rotated_H.first(0), 1e-4);
  BOOST_CHECK_CLOSE(E2_ref, rotated_H.first(1), 1e-4);

  Eigen::MatrixXd& diabatic_H = rotated_H.second;

  Eigen::MatrixXd diabatic_H_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/diabatization/Hdiab_ref.mm");

  BOOST_CHECK_CLOSE_FRACTION(diabatic_H_ref(0, 0), diabatic_H(0, 0), 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(diabatic_H_ref(0, 1), diabatic_H(0, 1), 1e-3);
  BOOST_CHECK_CLOSE_FRACTION(diabatic_H_ref(1, 0), diabatic_H(1, 0), 1e-3);
  BOOST_CHECK_CLOSE_FRACTION(diabatic_H_ref(1, 1), diabatic_H(1, 1), 1e-4);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(GMH_coupling_test) {
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

  GMHDiabatization GMHDiabatization(dimer, dimer, &log, 1, 2, "singlet");
  GMHDiabatization.configure();
  std::pair<double, double> coupling = GMHDiabatization.calculate_coupling();

  double J_ref = 0.011865823910700513;
  double J_proj_ref = 0.0009197395063185544;

  double J = coupling.first * votca::tools::conv::hrt2ev;
  double J_proj = coupling.second * votca::tools::conv::hrt2ev;

  BOOST_CHECK_CLOSE(J_ref, J, 1e-6);
  BOOST_CHECK_CLOSE(J_proj_ref, J_proj, 1e-6);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(FCD_coupling_test) {
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

  std::vector<QMFragment<BSE_Population> > fragments(2);

  fragments[0] = QMFragment<BSE_Population>(0, "0:8");
  fragments[1] = QMFragment<BSE_Population>(1, "9:17");

  FCDDiabatization FCDDiabatization(dimer, dimer, &log, 1, 2, "singlet",
                                    fragments);
  FCDDiabatization.configure();
  double coupling = FCDDiabatization.calculate_coupling();

  double J_ref = 0.00071879817182406039;
  double J = coupling * votca::tools::conv::hrt2ev;
  BOOST_CHECK_CLOSE_FRACTION(J_ref, J, 1e-5);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
