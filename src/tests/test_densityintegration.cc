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

#define BOOST_TEST_MODULE density_integration_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/density_integration.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/tools/eigenio_matrixmarket.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(density_integration_test)

AOBasis CreateBasis(const QMMolecule& mol) {

  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/densityintegration/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

BOOST_AUTO_TEST_CASE(density_test) {

  QMMolecule mol("none", 0);

  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/densityintegration/molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/densityintegration/dmat.mm");

  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);
  DensityIntegration<Vxc_Grid> num(grid);

  double ntot = num.IntegrateDensity(dmat);
  BOOST_CHECK_CLOSE(ntot, 8.000000, 1e-5);

  Eigen::Vector3d pos = {3, 3, 3};

  BOOST_CHECK_CLOSE(num.IntegratePotential(pos), -1.543242, 1e-4);

  Eigen::Vector3d field = num.IntegrateField(pos);
  Eigen::Vector3d field_ref = {0.172802, 0.172802, 0.172802};
  bool field_check = field.isApprox(field_ref, 1e-5);
  if (!field_check) {
    std::cout << "field" << std::endl;
    std::cout << field.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << field_ref.transpose() << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(gyration_test) {

  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            1.000000     1.000000     1.000000" << endl;
  xyzfile << " H            1.629118     1.629118     1.629118" << endl;
  xyzfile << " H           0.370882    0.370882     1.629118" << endl;
  xyzfile << " H            1.629118    0.370882    0.370882" << endl;
  xyzfile << " H           0.370882     1.629118   0.370882" << endl;
  xyzfile.close();

  QMMolecule mol("none", 0);

  mol.LoadFromFile("molecule.xyz");
  AOBasis aobasis = CreateBasis(mol);

  Eigen::MatrixXd dmat = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/densityintegration/dmat.mm");
  Vxc_Grid grid;
  grid.GridSetup("medium", mol, aobasis);
  DensityIntegration<Vxc_Grid> num(grid);

  Gyrationtensor tensor = num.IntegrateGyrationTensor(dmat);
  BOOST_CHECK_CLOSE(tensor.mass, 8.0000005, 1e-5);

  Eigen::Vector3d dip_ref = Eigen::Vector3d::Zero();
  dip_ref << 1.88973, 1.88973, 1.88973;
  bool centroid_check = dip_ref.isApprox(tensor.centroid, 1e-5);
  BOOST_CHECK_EQUAL(centroid_check, true);
  if (!centroid_check) {
    std::cout << "centroid" << std::endl;
    std::cout << tensor.centroid.transpose() << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << dip_ref.transpose() << std::endl;
  }
  Eigen::Matrix3d gyro_ref = Eigen::Matrix3d::Zero();
  gyro_ref << 0.596158, 2.85288e-12, 2.86873e-12, 2.85289e-12, 0.596158,
      2.87163e-12, 2.86874e-12, 2.87161e-12, 0.596158;
  bool gyro_check = gyro_ref.isApprox(tensor.gyration, 1e-5);
  BOOST_CHECK_EQUAL(gyro_check, true);
  if (!gyro_check) {
    std::cout << "gyro" << std::endl;
    std::cout << tensor.gyration << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << gyro_ref << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
