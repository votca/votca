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

#define BOOST_TEST_MODULE orbreorder_test

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>
#include <votca/tools/filesystem.h>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(orbreorder_test)

BOOST_AUTO_TEST_CASE(orbreorder_test) {
  libint2::initialize();
  // clang-format off
  std::array<Index,49> multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            -1,1,1,1,1,1,-1, //f 
            -1,-1,1,1,1,1,1,-1,-1, //g
            -1,-1,-1,1,1,1,1,1,-1,-1,-1, //h
            -1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1 //i
            };
  std::array<Index, 49> reorderList={
            0, //s
            1,-1,0, //p
            0,1,-1,2,-2, //d
            0,1,-1,2,-2,3,-3, //f 
            0,1,-1,2,-2,3,-3,4,-4, //g
            0,1,-1,2,-2,3,-3,4,-4,5,-5, //h
            0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6 //i
            };
  // clang-format on

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/orbreorder/methane.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/orbreorder/def2-tzvp.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd moldenCoeffs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbreorder/MOLDEN_order.mm");

  Eigen::MatrixXd votcaCoeffs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbreorder/VOTCA_order.mm");

  // Convert moldenCoeffs to votcaCoeffs
  OrbReorder reorder(reorderList, multipliers);
  reorder.reorderOrbitals(moldenCoeffs, aobasis);

  std::array<votca::Index, 49> votcaOrder_old = {
      0,                                           // s
      0, -1, 1,                                    // p
      0, -1, 1, -2, 2,                             // d
      0, -1, 1, -2, 2, -3, 3,                      // f
      0, -1, 1, -2, 2, -3, 3, -4, 4,               // g
      0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5,        // h
      0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5, -6, 6  // i
  };

  std::array<votca::Index, 49> multiplier2;
  multiplier2.fill(1);
  OrbReorder reorder2(votcaOrder_old, multiplier2);

  reorder2.reorderOrbitals(votcaCoeffs, aobasis);

  bool check = moldenCoeffs.isApprox(votcaCoeffs, 1e-5);
  BOOST_CHECK(check);
  if (!check) {
    std::cout << "votcaCoeffs" << std::endl;
    std::cout << votcaCoeffs << std::endl;
    std::cout << "moldenCoeffs" << std::endl;
    std::cout << moldenCoeffs << std::endl;
  }

  libint2::finalize();
}
}