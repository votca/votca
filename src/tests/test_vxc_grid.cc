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

#define BOOST_TEST_MODULE vxc_grid_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(vxc_grid_test)

AOBasis CreateBasis(const QMMolecule& mol) {

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/vxc_grid/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

BOOST_AUTO_TEST_CASE(vxc_grid_build) {
  libint2::initialize();
  QMMolecule mol("none", 0);

  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/vxc_grid/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);

  BOOST_CHECK_EQUAL(grid.getGridSize(), grid.getGridpoints().size());
  BOOST_CHECK_EQUAL(grid.getGridSize(), 53404);
  BOOST_CHECK_EQUAL(grid.getBoxesSize(), 51);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
