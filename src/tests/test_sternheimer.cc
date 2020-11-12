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
#include <libint2/initialize.h>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE sternheimer_test
#include <boost/test/unit_test.hpp>
#include <fstream>
#include <iostream>
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/sternheimer.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sternheimer_test)

BOOST_AUTO_TEST_CASE(sternheimer_polar) {
  libint2::initialize();

  std::cout << "Started Sternheimer polar test" << std::endl;

  Eigen::VectorXd mo_eigenvalues = Eigen::VectorXd::Zero(17);
  mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
      0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
      1.05054, 1.13372, 1.13372, 1.13372, 1.72964;
  Eigen::MatrixXd mo_eigenvectors =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gw/mo_eigenvectors.mm");

  Eigen::MatrixXd vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/vxc.mm");

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  BasisSet auxbasis;
  auxbasis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.setAuxbasisName(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");

  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  orbitals.MOs().eigenvectors() = mo_eigenvectors;
  orbitals.MOs().eigenvalues() = mo_eigenvalues;
  orbitals.setXCFunctionalName("XC_LDA_X XC_LDA_C_VWN");

  orbitals.setNumberOfAlphaElectrons(10);

  Logger log;

  Sternheimer sternheimer(orbitals, &log);

  Sternheimer::options_sternheimer options;

  options.do_precalc_fxc = true;

  sternheimer.configurate(options);

  sternheimer.setUpMatrices();

  std::cout << "Test Set Up Complete" << std::endl;

  std::vector<Eigen::Matrix3cd> polar = sternheimer.Polarisability();

  Eigen::Vector3d ref1;
  ref1 << -2.60731, -2.61376, -2.61588;
  Eigen::Vector3d ref2;
  ref2 << -3.33651, -3.33445, -3.33956;

  bool check_0 = polar[0].diagonal().real().isApprox(ref1, 1e-4);
  bool check_20 = polar[20].diagonal().real().isApprox(ref2, 1e-4);

  if (!check_0) {
    std::cout << "*** 0 ***" << std::endl;
    std::cout << "res" << std::endl;
    std::cout << polar[0].diagonal().real() << std::endl;
    std::cout << "ref " << std::endl;
    std::cout << ref1 << std::endl;
  }
  BOOST_CHECK_EQUAL(check_0, true);

  if (!check_20) {
    std::cout << "*** 20 ***" << std::endl;
    std::cout << "res" << std::endl;
    std::cout << polar[20].diagonal().real() << std::endl;
    std::cout << "ref " << std::endl;
    std::cout << ref2 << std::endl;
  }
  BOOST_CHECK_EQUAL(check_20, true);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(sternheimer_energy_gradient) {

  std::cout << "Started Sternheimer energy gradient test" << std::endl;

  Eigen::VectorXd mo_eigenvalues = Eigen::VectorXd::Zero(17);
  mo_eigenvalues << -10.6784, -0.746424, -0.394948, -0.394948, -0.394948,
      0.165212, 0.227713, 0.227713, 0.227713, 0.763971, 0.763971, 0.763971,
      1.05054, 1.13372, 1.13372, 1.13372, 1.72964;
  Eigen::MatrixXd mo_eigenvectors =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/gw/mo_eigenvectors.mm");

  Eigen::MatrixXd vxc = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/gw/vxc.mm");

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/gw/molecule.xyz");

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");

  BasisSet auxbasis;
  auxbasis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");
  orbitals.setAuxbasisName(std::string(XTP_TEST_DATA_FOLDER) + "/gw/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  orbitals.MOs().eigenvectors() = mo_eigenvectors;
  orbitals.MOs().eigenvalues() = mo_eigenvalues;
  orbitals.setXCFunctionalName("XC_LDA_X XC_LDA_C_VWN");
  orbitals.setNumberOfAlphaElectrons(10);

  Logger log;

  Sternheimer sternheimer(orbitals, &log);

  Sternheimer::options_sternheimer options;

  options.do_precalc_fxc = true;
  options.max_mixing_history = 20;

  sternheimer.configurate(options);

  sternheimer.setUpMatrices();

  std::cout << "Test Set Up Complete" << std::endl;

  std::vector<Eigen::Vector3cd> EPC = sternheimer.EnergyGradient();
  std::vector<Eigen::Vector3cd> EPCmo_0 = sternheimer.MOEnergyGradient(0, 0);
  std::vector<Eigen::Vector3cd> EPCmo_4 = sternheimer.MOEnergyGradient(4, 4);

  std::cout << "EPC[0]" << EPC[0].real() << std::endl;
  std::cout << "EPC[4]" << EPC[4].real() << std::endl;
  std::cout << "EPCmo0[0]" << EPCmo_0[0].real() << std::endl;
  std::cout << "EPCmo0[4]" << EPCmo_0[4].real() << std::endl;
  std::cout << "EPCmo4[0]" << EPCmo_4[0].real() << std::endl;
  std::cout << "EPCmo4[4]" << EPCmo_4[4].real() << std::endl;

  Eigen::Vector3d ref1;
  ref1 << -4.02724e+54, 0.838474, 0.233327;
  Eigen::Vector3d ref2;
  ref2 << -0.196696, 0.535469, 0.43545;
  Eigen::Vector3d ref3;
  ref3 << -1.29072e-14, 23.3186, -1.3997;
  Eigen::Vector3d ref4;
  ref4 << 0.136145, -4.68696e+36, 0.136145;
  Eigen::Vector3d ref5;
  ref5 << 7.06401, 1.50099, 0.307372;
  Eigen::Vector3d ref6;
  ref6 << -2.68194e+123, -5.08698e+114, 3.06678e+120;

  bool check_1 = EPC[0].real().isApprox(ref1, 1e-4);
  bool check_2 = EPC[4].real().isApprox(ref2, 1e-4);
  bool check_3 = EPCmo_0[0].real().isApprox(ref1, 1e-4);
  bool check_4 = EPCmo_0[4].real().isApprox(ref2, 1e-4);
  bool check_5 = EPCmo_4[0].real().isApprox(ref1, 1e-4);
  bool check_6 = EPCmo_4[4].real().isApprox(ref2, 1e-4);

  BOOST_CHECK_EQUAL(check_1, true);
  BOOST_CHECK_EQUAL(check_2, true);
  BOOST_CHECK_EQUAL(check_3, true);
  BOOST_CHECK_EQUAL(check_4, true);
  BOOST_CHECK_EQUAL(check_5, true);
  BOOST_CHECK_EQUAL(check_6, true);
}

BOOST_AUTO_TEST_SUITE_END()
