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

#define BOOST_TEST_MODULE orca_test

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

  Orbitals orbitalsReference;
  orbitalsReference.ReadFromCpt(std::string(XTP_TEST_DATA_FOLDER) +
                                "/molden/benzene.orb");

  Logger log;
  MoldenReader molden(log);
  molden.setBasissetInfo("def2-tzvp", "aux-def2-tzvp");
  Orbitals orbitals;
  molden.parseMoldenFile(
      std::string(XTP_TEST_DATA_FOLDER) + "/molden/benzene.molden.input",
      orbitals);

  BOOST_CHECK(orbitalsReference.MOs().eigenvectors().isApprox(
      orbitals.MOs().eigenvectors(), 1e-5));

  BOOST_CHECK(orbitals.QMAtoms().size() == orbitalsReference.QMAtoms().size());

  for (int i = 0; i < orbitals.QMAtoms().size(); i++) {
    BOOST_CHECK(orbitals.QMAtoms()[i].getPos().isApprox(
        orbitalsReference.QMAtoms()[i].getPos(), 1e-3));
  }
}
}