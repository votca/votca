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

#define BOOST_TEST_MODULE aomatrix3ddipole_test

// Third party includes
#include <boost/test/unit_test.hpp>

#include <votca/tools/eigenio_matrixmarket.h>
// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/orbitals.h"
#include <libint2/initialize.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(aomatrix3ddipole_test)

BOOST_AUTO_TEST_CASE(aomatrices3ddipole_test) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/aomatrix3d/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/aomatrix3d/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  AO3ddipole dip;
  dip.Fill(aobasis);

  for (unsigned i= 0; i < 3; i++){
    std::cout << "Matrix " << i << std::endl;
    std::cout << dip.Matrix()[i] << std::endl;
  }


  // Just here so we can generate output with --output-on-failure
  BOOST_CHECK_EQUAL(false, true);

  libint2::finalize();
}


BOOST_AUTO_TEST_SUITE_END()
