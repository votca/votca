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

#define BOOST_TEST_MODULE vxc_potential_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(vxc_potential_test)

AOBasis CreateBasis(const QMMolecule& mol) {
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

Eigen::MatrixXd DMat() {
  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/dmat.mm");
  return dmat;
}

BOOST_AUTO_TEST_CASE(vxc_test) {

  QMMolecule mol("none", 0);

  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER)+ "/vxc_potential/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = DMat();
  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);
  Vxc_Potential<Vxc_Grid> num(grid);
  num.setXCfunctional("XC_GGA_X_PBE XC_GGA_C_PBE");

  BOOST_CHECK_EQUAL(grid.getGridSize(), grid.getGridpoints().size());
  BOOST_CHECK_EQUAL(grid.getGridSize(), 53404);
  BOOST_CHECK_EQUAL(grid.getBoxesSize(), 51);

  BOOST_CHECK_CLOSE(num.getExactExchange("XC_GGA_X_PBE XC_GGA_C_PBE"), 0.0,
                    1e-5);
  BOOST_CHECK_CLOSE(num.getExactExchange("XC_HYB_GGA_XC_PBEH"), 0.25, 1e-5);
  Mat_p_Energy e_vxc = num.IntegrateVXC(dmat);


  Eigen::MatrixXd vxc_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
    std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/vxc_ref.mm");
    
  bool check_vxc = e_vxc.matrix().isApprox(vxc_ref, 0.0001);

  BOOST_CHECK_CLOSE(e_vxc.energy(), -4.6303432151572643, 1e-5);
  if (!check_vxc) {
    std::cout << "ref" << std::endl;
    std::cout << vxc_ref << std::endl;
    std::cout << "calc" << std::endl;
    std::cout << e_vxc.matrix() << std::endl;
  }
  BOOST_CHECK_EQUAL(check_vxc, 1);
}

BOOST_AUTO_TEST_SUITE_END()
