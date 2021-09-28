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

#define BOOST_TEST_MODULE regular_grid_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/regular_grid.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(regular_grid_test)

AOBasis CreateBasis(const QMMolecule& mol) {

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/regular_grid/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

BOOST_AUTO_TEST_CASE(regular_grid_build) {
  libint2::initialize();
  QMMolecule mol("none", 0);

  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/regular_grid/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Regular_Grid grid;
  Eigen::Array<votca::Index, 3, 1> steps(5, 5, 5);
  Eigen::Array3d padding(1.0, 1.0, 1.0);
  grid.GridSetup(steps, padding, mol, aobasis);

  BOOST_CHECK_EQUAL(grid.getGridSize(), 125);
  BOOST_CHECK_EQUAL(grid.getBoxesSize(), 1);

  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[0].x(), -2.1888606344960548, 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[1].z(), -1.0944303172480274, 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[2].z(), 0, 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[30].y(), -1.0944303172480274, 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[45].z(), -2.18886063, 1e-5);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(regular_grid_build_large) {
  libint2::initialize();
  QMMolecule mol("none", 0);

  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/regular_grid/molecule.xyz");
  auto extend = mol.CalcSpatialMinMax();
  AOBasis aobasis = CreateBasis(mol);

  Regular_Grid grid;
  Eigen::Array<votca::Index, 3, 1> steps(30, 30, 30);
  Eigen::Array3d padding(1.0, 1.0, 1.0);
  grid.GridSetup(steps, padding, mol, aobasis);
  auto max = extend.second.array() + padding + 1e-9;
  auto min = extend.first.array() - padding - 1e-9;
  for (auto& a : grid) {
    for (auto& point : a.getGridPoints()) {
      bool check_larger = (point.array() > max).any();
      BOOST_CHECK_EQUAL(check_larger, false);
      if (check_larger) {
        std::cout << point.transpose() << std::endl;
        std::cout << max.transpose() << std::endl;
      }
      bool check_smaller = (point.array() < min).any();
      BOOST_CHECK_EQUAL(check_smaller, false);
      if (check_smaller) {
        std::cout << point.transpose() << std::endl;
        std::cout << min.transpose() << std::endl;
      }
    }
  }

  BOOST_CHECK_EQUAL(grid.getGridSize(), 27000);
  BOOST_CHECK_EQUAL(grid.getBoxesSize(), 54);

  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[0].x(), min.x(), 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[0].y(), min.y(), 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[0].z(), min.z(), 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[30].y(), -2.0379047286687406, 1e-5);
  BOOST_CHECK_CLOSE(grid[0].getGridPoints()[45].z(), 0.075477952913657109,
                    1e-5);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
