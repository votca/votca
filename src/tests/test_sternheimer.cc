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

#define BOOST_TEST_MODULE sternheimer_test
#include <boost/test/unit_test.hpp>
#include <eigen3/Eigen/src/Core/Matrix.h>
#include <fstream>
#include <iostream>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/sternheimer.h>
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sternheimer_test)

BOOST_AUTO_TEST_CASE(sternheimer_full) {

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
  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  orbitals.MOs().eigenvectors() = mo_eigenvectors;
  Logger log;

  Sternheimer sternheimer(orbitals, &log);

  Sternheimer::options_sternheimer options;

  options.do_precalc_fxc=true;

  sternheimer.configurate(options);

  std::vector<Eigen::Matrix3cd> polar = sternheimer.Polarisability();
  std::vector<Eigen::Vector3cd> EPC = sternheimer.EnergyGradient();
  std::vector<Eigen::Vector3cd> EPC = sternheimer.MOEnergyGradient(n, n);



}

BOOST_AUTO_TEST_SUITE_END()
