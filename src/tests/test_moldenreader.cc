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

#define BOOST_TEST_MODULE moldenreader_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/moldenreader.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(moldenreader_test)

BOOST_AUTO_TEST_CASE(moldenreader_test) {

  Eigen::MatrixXd coeffs_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/molden/orbitalsMOs_ref.mm");

  Orbitals orbitals_ref;
  orbitals_ref.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                      "/molden/benzene.xyz");

  Logger log;
  MoldenReader molden(log);
  molden.setBasissetInfo(
      std::string(XTP_TEST_DATA_FOLDER) + "/molden/def2-tzvp.xml",
      "aux-def2-tzvp");
  Orbitals orbitals;
  molden.parseMoldenFile(
      std::string(XTP_TEST_DATA_FOLDER) + "/molden/benzene.molden.input",
      orbitals);

  // Check if MO's are read correctly
  BOOST_CHECK(orbitals.MOs().eigenvectors().isApprox(coeffs_ref, 1e-5));

  // Check if atoms are read correctly
  BOOST_CHECK(orbitals.QMAtoms().size() == orbitals_ref.QMAtoms().size());
  for (int i = 0; i < orbitals.QMAtoms().size(); i++) {
    BOOST_CHECK(orbitals.QMAtoms()[i].getPos().isApprox(
        orbitals_ref.QMAtoms()[i].getPos(), 1e-3));
  }
}
}