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

#define BOOST_TEST_MODULE orbitals_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/orbitals.h"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using votca::Index;
BOOST_AUTO_TEST_SUITE(orbitals_test)
BOOST_AUTO_TEST_CASE(readxyztest) {
  std::ofstream xyzfile("molecule.xyz");
  xyzfile << " C 0.0 3 1" << std::endl;
  xyzfile << " methane" << std::endl;
  xyzfile << " C            .000000     .000000     .000000" << std::endl;
  xyzfile << " H            .629118     .629118     .629118" << std::endl;
  xyzfile << " H           -.629118    -.629118     .629118" << std::endl;
  xyzfile << " H            .629118    -.629118    -.629118" << std::endl;
  xyzfile << " H           -.629118     .629118    -.629118" << std::endl;
  xyzfile.close();

  Orbitals orb;

  BOOST_CHECK_THROW(orb.QMAtoms().LoadFromFile("molecule.xyz"),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(sortEnergies) {

  Orbitals orb;
  Eigen::VectorXd Energies = Eigen::VectorXd::LinSpaced(10, -5, 5);
  Eigen::VectorXd switched = Energies;
  switched(3) = Energies(5);
  switched(5) = Energies(3);
  orb.MOs().eigenvalues() = switched;
  orb.MOs().eigenvectors() = Eigen::MatrixXd::Zero(10, 10);
  orb.OrderMOsbyEnergy();
  bool issorted = Energies.isApprox(orb.MOs().eigenvalues(), 0.001);
  if (!issorted) {
    std::cout << "before" << std::endl;
    std::cout << Energies << std::endl;
    std::cout << "after" << std::endl;
    std::cout << orb.MOs().eigenvalues() << std::endl;
  }
  BOOST_CHECK_EQUAL(issorted, true);
}

BOOST_AUTO_TEST_CASE(densmat_test) {

  Orbitals orb;
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(4);
  orb.setBSEindices(0, 9);
  orb.setNumberOfAlphaElectrons(5);
  orb.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);
  orb.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/MOs.mm");

  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/orbitals/molecule.xyz");
  QMState s = QMState("n");
  Eigen::MatrixXd dmat_gs = orb.DensityMatrixFull(s);

  Eigen::MatrixXd dmat_gs_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/dmat_gs_ref.mm");

  bool check_dmat_gs = dmat_gs.isApprox(dmat_gs_ref, 0.0001);
  if (!check_dmat_gs) {
    std::cout << "Result gs" << std::endl;
    std::cout << dmat_gs << std::endl;
    std::cout << "Ref" << std::endl;
    std::cout << dmat_gs_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(check_dmat_gs, 1);
  orb.setTDAApprox(false);
  orb.BSESinglets().eigenvectors() =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/BSE_vectors1.mm");

  orb.BSESinglets().eigenvectors2() =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/BSE_vectors2.mm");

  QMState s1 = QMState("s1");
  Eigen::MatrixXd dmat_s1 = orb.DensityMatrixFull(s1);

  Eigen::MatrixXd dmat_s1_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/dmat_s1_ref.mm");

  bool check_dmat_s1 = dmat_s1.isApprox(dmat_s1_ref, 0.0001);
  if (!check_dmat_s1) {
    std::cout << "Result s1" << std::endl;
    std::cout << dmat_s1 << std::endl;
    std::cout << "Ref" << std::endl;
    std::cout << dmat_s1_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(check_dmat_s1, 1);

  QMState n2s1 = QMState("n2s1");
  Eigen::MatrixXd dmat_n2s1 = orb.DensityMatrixFull(n2s1);
  Eigen::MatrixXd dmat_n2s1_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/dmat_n2s1_ref.mm");

  bool check_dmat_n2s1 = dmat_n2s1.isApprox(dmat_n2s1_ref, 0.0001);
  if (!check_dmat_n2s1) {
    std::cout << "Result n2s1" << std::endl;
    std::cout << dmat_n2s1 << std::endl;
    std::cout << "Ref" << std::endl;
    std::cout << dmat_n2s1_ref << std::endl;
  }

  BOOST_CHECK_EQUAL(check_dmat_n2s1, 1);
}

BOOST_AUTO_TEST_CASE(dipole_test) {

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/orbitals/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/3-21G.xml");
  orbitals.setDFTbasisName("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
  orbitals.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/MOs2.mm");

  orbitals.setBSEindices(0, 16);
  orbitals.setTDAApprox(true);

  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/spsi_ref.mm");

  orbitals.BSESinglets().eigenvectors() = spsi_ref;
  QMState state_trans = QMState("n2s1");

  Eigen::Vector3d res_trans = orbitals.CalcElDipole(state_trans);
  Eigen::Vector3d ref_trans = Eigen::Vector3d::Zero();
  ref_trans << 0.118565, 0.0444239, -0.0505149;

  bool check_trans = ref_trans.isApprox(res_trans, 0.0001);
  if (!check_trans) {
    std::cout << "Result transition dipole" << std::endl;
    std::cout << res_trans << std::endl;
    std::cout << "Ref transition dipole" << std::endl;
    std::cout << ref_trans << std::endl;
  }
  BOOST_CHECK_EQUAL(check_trans, 1);

  QMState state_s1 = QMState("s1");
  Eigen::Vector3d res_s1 = orbitals.CalcElDipole(state_s1);
  Eigen::Vector3d ref_s1 = Eigen::Vector3d::Zero();
  ref_s1 << -0.15153501734, -0.42406579479, 0.033954362839;
  bool check_s1 = ref_s1.isApprox(res_s1, 0.0001);
  if (!check_s1) {
    std::cout << "Result s1 dipole" << std::endl;
    std::cout << res_s1 << std::endl;
    std::cout << "Ref s1 dipole" << std::endl;
    std::cout << ref_s1 << std::endl;
  }
  BOOST_CHECK_EQUAL(check_s1, 1);
}

BOOST_AUTO_TEST_CASE(osc_strength) {
  Orbitals orb;
  orb.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                             "/orbitals/molecule.xyz");
  orb.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                      "/orbitals/3-21G.xml");

  QMState s("s1");
  orb.setBasisSetSize(17);
  orb.setNumberOfOccupiedLevels(4);
  orb.setBSEindices(0, 16);
  orb.setTDAApprox(true);

  Eigen::MatrixXd& MOs = orb.MOs().eigenvectors();
  orb.MOs().eigenvalues() = Eigen::VectorXd::Ones(17);
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/MOs3.mm");
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;
  orb.BSESinglets().eigenvalues() = se_ref;
  // reference coefficients
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/orbitals/spsi_ref2.mm");

  orb.BSESinglets().eigenvectors() = spsi_ref;
  orb.CalcCoupledTransition_Dipoles();

  std::vector<Eigen::Vector3d> dipoles = orb.TransitionDipoles();
  std::vector<Eigen::Vector3d> dipoles_ref;
  dipoles_ref.push_back(Eigen::Vector3d(0.110512, 0.048776, -0.0515914));
  dipoles_ref.push_back(Eigen::Vector3d(-0.13408, 0.0969472, 0.0261392));
  dipoles_ref.push_back(Eigen::Vector3d(0.0586073, 0.121606, -0.0606862));

  for (Index i = 0; i < 3; i++) {
    bool check = dipoles[i].isApprox(dipoles_ref[i], 1e-5);
    BOOST_CHECK_EQUAL(check, true);
    if (!check) {
      std::cout << "ref" << i << std::endl;
      std::cout << dipoles_ref[i].transpose() << std::endl;
      std::cout << "result" << i << std::endl;
      std::cout << dipoles[i].transpose() << std::endl;
    }
  }

  Eigen::VectorXd oscs = orb.Oscillatorstrengths();
  Eigen::VectorXd oscs_ref = Eigen::VectorXd::Zero(3);
  oscs_ref << 0.001236, 0.00201008, 0.00156925;

  bool check_oscs = oscs.isApprox(oscs_ref, 1e-5);
  BOOST_CHECK_EQUAL(check_oscs, true);
  if (!check_oscs) {
    std::cout << "result" << std::endl;
    std::cout << oscs << std::endl;
    std::cout << "ref" << std::endl;
    std::cout << oscs_ref << std::endl;
  }
}

BOOST_AUTO_TEST_SUITE_END()
