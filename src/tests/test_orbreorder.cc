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
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(orbreorder_test)

BOOST_AUTO_TEST_CASE(orbreorder_test) {

  // clang-format off
  std::array<Index,25> _multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,-1,-1, //f 
            1,1,1,1,1,-1,-1,-1,-1 //g
            };

  OrbTranspositions _transpositions { 
    std::vector<std::array<Index, 2>> {}, //s
    std::vector<std::array<Index, 2>> {   //p
      {0, 2}
    }, 
    std::vector<std::array<Index, 2>> {   //d
      {1, 2},
      {3, 4}
      }, 
    std::vector<std::array<Index, 2>> {   //f
      {1, 2},  
      {3, 4},
      {5, 6}
    }, 
    std::vector<std::array<Index, 2>> {   //g
      {1, 2},
      {3, 4},
      {5, 6},
      {7, 8}
    }
  };
  //clang-format on


  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/orbreorder/methane.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/orbreorder/def2-tzvp.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());


  Eigen::MatrixXd moldenCoeffs =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/orbreorder/MOLDEN_order.mm");
  //moldenCoeffs.transposeInPlace();
  std::cout << moldenCoeffs << std::endl;
  std::cout << std::endl;

  Eigen::MatrixXd votcaCoeffs =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/orbreorder/VOTCA_order.mm");
  //votcaCoeffs.transposeInPlace();
  
  std::cout << votcaCoeffs << std::endl;
  std::cout << std::endl;
 
  // Convert moldenCoeffs to votcaCoeffs
  OrbReorder reorder(_transpositions, _multipliers);
  reorder.reorderOrbitals(moldenCoeffs, aobasis);

  std::cout << moldenCoeffs << std::endl;
  std::cout << std::endl;

  for(int i = 0; i < moldenCoeffs.cols(); i++){
    for(int j = 0; j < moldenCoeffs.cols(); j++){
      std::cout << moldenCoeffs(j,i) << " " << votcaCoeffs(j,i) << std::endl;
      BOOST_CHECK(std::abs(moldenCoeffs(j,i) - votcaCoeffs(j,i)) < 1e-5);
      
    }
  }
}
}


