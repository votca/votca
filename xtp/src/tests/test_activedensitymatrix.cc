/*
 * Copyright 2009-2022 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE activedensitymatrix_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/activedensitymatrix.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(activedensitymatrix_test)
BOOST_AUTO_TEST_CASE(activematrix_test) {

  libint2::initialize();
  Orbitals orbitals_;
  Orbitals orb2;
  orbitals_.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                   "/activedensitymatrix/ch3oh.xyz");
  orbitals_.setNumberOfOccupiedLevels(9);
  orbitals_.setNumberOfAlphaElectrons(9);

  orbitals_.SetupDftBasis(std::string(XTP_TEST_DATA_FOLDER) +
                          "/activedensitymatrix/def2-tzvp.xml");
  Eigen::MatrixXd LMOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/activedensitymatrix/LMOs.mm");

  orbitals_.setLMOs(LMOs);
  orb2.setLMOs(LMOs);

  Logger log;
  std::vector<Index> activeatoms = {{1, 5}};

  ActiveDensityMatrix DMAT_A(orbitals_, activeatoms);

  Eigen::MatrixXd Dmat = DMAT_A.compute_Dmat_A();
  Eigen::MatrixXd DmatA = Dmat;

  Eigen::MatrixXd test_DmatA = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/activedensitymatrix/ch3oh.mm");
  bool checkDmatA = DmatA.isApprox(test_DmatA, 2e-6);
  BOOST_CHECK_EQUAL(checkDmatA, 1);

  libint2::finalize();
}
BOOST_AUTO_TEST_SUITE_END()