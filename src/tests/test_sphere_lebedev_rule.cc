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

#define BOOST_TEST_MODULE sphere_lebedev_rule_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/sphere_lebedev_rule.h"
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(sphere_lebedev_rule_test)

BOOST_AUTO_TEST_CASE(medium_test) {

  QMMolecule mol("noname", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/sphere_lebedev_rule/CH4.xyz");

  LebedevGrid spheregrid;

  auto grid = spheregrid.CalculateSphericalGrids(mol, "medium");

  auto Hgrid = grid.at("H");
  auto Cgrid = grid.at("C");

  Eigen::VectorXd C_phi_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/C_phi_ref_medium.mm");
  C_phi_ref *= votca::tools::conv::Pi / 180.0;

  Eigen::VectorXd C_theta_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/C_theta_ref_medium.mm");

  Eigen::VectorXd C_weight_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/C_weight_ref_medium.mm");
  BOOST_CHECK_EQUAL(Cgrid.phi.size(), C_phi_ref.size());
  BOOST_CHECK_EQUAL(Cgrid.theta.size(), C_theta_ref.size());
  BOOST_CHECK_EQUAL(Cgrid.weight.size(), C_weight_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.phi.size(), C_phi_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.theta.size(), C_theta_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.weight.size(), C_weight_ref.size());

  bool Cphi = C_phi_ref.isApprox(Cgrid.phi, 0.001);
  bool Ctheta = C_theta_ref.isApprox(Cgrid.theta, 0.001);
  if (!Cphi || !Ctheta) {
    std::cout << "phi_ref : Phi_comp | theta_ref : theta_comp" << std::endl;
    for (Index i = 0; i < C_phi_ref.size(); i++) {
      std::cout << Cgrid.phi[i] << ":" << C_phi_ref[i] << " | "
                << Cgrid.theta[i] << ":" << C_theta_ref[i] << std::endl;
    }
  }
  bool Cweight = C_weight_ref.isApprox(Cgrid.weight, 0.0001);
  BOOST_CHECK_EQUAL(Cphi, true);
  BOOST_CHECK_EQUAL(Ctheta, true);
  BOOST_CHECK_EQUAL(Cweight, true);

  bool Hphi = C_phi_ref.isApprox(Hgrid.phi, 0.001);
  bool Htheta = C_theta_ref.isApprox(Hgrid.theta, 0.001);
  bool Hweight = C_weight_ref.isApprox(Hgrid.weight, 0.0001);
  BOOST_CHECK_EQUAL(Hphi, true);
  BOOST_CHECK_EQUAL(Htheta, true);
  BOOST_CHECK_EQUAL(Hweight, true);
}

BOOST_AUTO_TEST_CASE(fine_test) {

  QMMolecule mol("noname", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/sphere_lebedev_rule/molecule.xyz");
  LebedevGrid spheregrid;

  auto grid = spheregrid.CalculateSphericalGrids(mol, "fine");

  auto Hgrid = grid.at("H");
  auto Gegrid = grid.at("Ge");

  Eigen::VectorXd Ge_phi_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/Ge_phi_ref_fine.mm");

  Eigen::VectorXd Ge_theta_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/Ge_theta_ref_fine.mm");

  Eigen::VectorXd Ge_weight_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/sphere_lebedev_rule/Ge_weight_ref_fine.mm");

  Eigen::VectorXd H_phi_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/H_phi_ref_fine.mm");

  Eigen::VectorXd H_theta_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/H_theta_ref_fine.mm");

  Eigen::VectorXd H_weight_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/H_weight_ref_fine.mm");

  BOOST_CHECK_EQUAL(Gegrid.phi.size(), Ge_phi_ref.size());
  BOOST_CHECK_EQUAL(Gegrid.theta.size(), Ge_theta_ref.size());
  BOOST_CHECK_EQUAL(Gegrid.weight.size(), Ge_weight_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.phi.size(), H_phi_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.theta.size(), H_theta_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.weight.size(), H_weight_ref.size());

  bool Gephi = Ge_phi_ref.isApprox(Gegrid.phi, 0.001);
  bool Getheta = Ge_theta_ref.isApprox(Gegrid.theta, 0.001);
  if (!Gephi || !Getheta) {
    std::cout << "phi_ref : Phi_comp | theta_ref : theta_comp" << std::endl;
    for (Index i = 0; i < Ge_phi_ref.size(); i++) {
      std::cout << Gegrid.phi[i] << ":" << Ge_phi_ref[i] << " | "
                << Gegrid.theta[i] << ":" << Ge_theta_ref[i] << std::endl;
    }
  }

  bool Geweight = Ge_weight_ref.isApprox(Gegrid.weight, 0.0001);
  BOOST_CHECK_EQUAL(Gephi, true);
  BOOST_CHECK_EQUAL(Getheta, true);
  BOOST_CHECK_EQUAL(Geweight, true);

  bool Hphi = H_phi_ref.isApprox(Hgrid.phi, 0.001);
  bool Htheta = H_theta_ref.isApprox(Hgrid.theta, 0.001);
  bool Hweight = H_weight_ref.isApprox(Hgrid.weight, 0.0001);
  BOOST_CHECK_EQUAL(Hphi, true);
  BOOST_CHECK_EQUAL(Htheta, true);
  BOOST_CHECK_EQUAL(Hweight, true);
}
BOOST_AUTO_TEST_CASE(element_not_implemented) {

  QMMolecule mol("noname", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/sphere_lebedev_rule/hg.xyz");

  LebedevGrid spheregrid;

  BOOST_REQUIRE_THROW(spheregrid.CalculateSphericalGrids(mol, "xfine"),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(xfine_test) {

  QMMolecule mol("noname", 0);
  mol.LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/sphere_lebedev_rule/molecule.xyz");

  LebedevGrid spheregrid;

  auto grid = spheregrid.CalculateSphericalGrids(mol, "xfine");

  auto Hgrid = grid.at("H");
  auto Gegrid = grid.at("Ge");

  Eigen::VectorXd Ge_phi_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/Ge_phi_ref_xfine.mm");

  Eigen::VectorXd Ge_theta_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/Ge_theta_ref_xfine.mm");

  Eigen::VectorXd Ge_weight_ref =
      votca::tools::EigenIO_MatrixMarket::ReadVector(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/sphere_lebedev_rule/Ge_weight_ref_xfine.mm");

  Eigen::VectorXd H_phi_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/H_phi_ref_xfine.mm");

  Eigen::VectorXd H_theta_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/H_theta_ref_xfine.mm");
  Eigen::VectorXd H_weight_ref = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) +
      "/sphere_lebedev_rule/H_weight_ref_xfine.mm");
  BOOST_CHECK_EQUAL(Gegrid.phi.size(), Ge_phi_ref.size());
  BOOST_CHECK_EQUAL(Gegrid.theta.size(), Ge_theta_ref.size());
  BOOST_CHECK_EQUAL(Gegrid.weight.size(), Ge_weight_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.phi.size(), H_phi_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.theta.size(), H_theta_ref.size());
  BOOST_CHECK_EQUAL(Hgrid.weight.size(), H_weight_ref.size());

  bool Gephi = Ge_phi_ref.isApprox(Gegrid.phi, 0.001);
  bool Getheta = Ge_theta_ref.isApprox(Gegrid.theta, 0.001);
  if (!Gephi || !Getheta) {
    std::cout << "phi_ref : Phi_comp | theta_ref : theta_comp" << std::endl;
    for (Index i = 0; i < Ge_phi_ref.size(); i++) {
      std::cout << Gegrid.phi[i] << ":" << Ge_phi_ref[i] << " | "
                << Gegrid.theta[i] << ":" << Ge_theta_ref[i] << std::endl;
    }
  }

  bool Geweight = Ge_weight_ref.isApprox(Gegrid.weight, 0.0001);
  BOOST_CHECK_EQUAL(Gephi, true);
  BOOST_CHECK_EQUAL(Getheta, true);
  BOOST_CHECK_EQUAL(Geweight, true);

  bool Hphi = H_phi_ref.isApprox(Hgrid.phi, 0.001);
  bool Htheta = H_theta_ref.isApprox(Hgrid.theta, 0.001);
  bool Hweight = H_weight_ref.isApprox(Hgrid.weight, 0.0001);
  BOOST_CHECK_EQUAL(Hphi, true);
  BOOST_CHECK_EQUAL(Htheta, true);
  BOOST_CHECK_EQUAL(Hweight, true);
}

BOOST_AUTO_TEST_SUITE_END()
