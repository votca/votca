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

#define BOOST_TEST_MODULE pmdecomposition_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/pmdecomposition.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(pmdecomposition_test)
BOOST_AUTO_TEST_CASE(decomposedorbitals_test) {

  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/pmdecomposition/benzene.xyz");
  orbitals.setBasisSetSize(222);
  orbitals.setNumberOfOccupiedLevels(21);

  orbitals.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                           "/pmdecomposition/def2-tzvp.xml");

  orbitals.MOs().eigenvectors() =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/pmdecomposition/orbitalsMOs_ref.mm");

  orbitals.setNumberOfAlphaElectrons(21);

  Logger log;

  PMDecomposition pmd(log);
  pmd.computePMD(orbitals);

  Eigen::MatrixXd LMOs = orbitals.getPMLocalizedOrbitals();
  Eigen::MatrixXd test_MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/pmdecomposition/test_pm_orbitals.mm");

  bool checkMOs = LMOs.isApprox(test_MOs, 1e-3);
  BOOST_CHECK_EQUAL(checkMOs, 1);

  libint2::finalize();
}
BOOST_AUTO_TEST_SUITE_END()