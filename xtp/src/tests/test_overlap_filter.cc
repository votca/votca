/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE overlap_filter_test

// Standard includes
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include <votca/xtp/filterfactory.h>

// VOTCA includes
#include <libint2/initialize.h>
#include <votca/tools/eigenio_matrixmarket.h>
using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(overlap_filter_test)

BOOST_AUTO_TEST_CASE(coeffs_test) {

  libint2::initialize();
  FilterFactory factory;
  std::unique_ptr<StateFilter_base> rho_f = factory.Create("overlap");

  std::ofstream opt("overlap_filter.xml");
  opt << "<overlap>0.0045</overlap>" << std::endl;
  opt.close();
  votca::tools::Property prop;
  prop.LoadFromXML("overlap_filter.xml");
  rho_f->Initialize(prop.get("overlap"));

  Orbitals A;
  A.setDFTbasisName(std::string(XTP_TEST_DATA_FOLDER) +
                    "/overlap_filter/3-21G.xml");
  A.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                           "/overlap_filter/molecule.xyz");
  A.setBasisSetSize(17);
  A.setNumberOfAlphaElectrons(5);
  A.setNumberOfOccupiedLevels(5);
  A.MOs().eigenvalues() = Eigen::VectorXd::Zero(17);
  A.MOs().eigenvalues() << -19.8117, -6.22408, -6.14094, -6.14094, -6.14094,
      -3.72889, -3.72889, -3.72889, -3.64731, -3.09048, -3.09048, -3.09048,
      -2.63214, -2.08206, -2.08206, -2.08206, -2.03268;

  A.QPpertEnergies() = Eigen::VectorXd::Zero(17);
  A.QPpertEnergies() << -10.189, -1.01045, -0.620145, -0.620146, -0.620148,
      0.261183, 0.348342, 0.348343, 0.348342, 0.920829, 0.920829, 0.920829,
      1.18002, 1.27325, 1.27325, 1.27325, 1.96983;

  A.MOs().eigenvectors() = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/overlap_filter/MOs_A.mm");

  A.setBSEindices(0, 16);
  A.setTDAApprox(true);
  A.setGWindices(0, 16);
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/overlap_filter/spsi_ref.mm");
  A.BSESinglets().eigenvectors() = spsi_ref;

  // reference energy
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;

  A.BSESinglets().eigenvalues() = se_ref;
  A.CalcCoupledTransition_Dipoles();

  BOOST_CHECK_EQUAL(rho_f->NeedsInitialState(), true);

  rho_f->UpdateHist(A, QMState("s1"));

  std::vector<votca::Index> ref = {0};
  std::vector<votca::Index> results =
      rho_f->CalcIndeces(A, QMStateType::Singlet);

  BOOST_CHECK_EQUAL(results.size(), ref.size());
  for (votca::Index i = 0; i < votca::Index(ref.size()); i++) {
    BOOST_CHECK_EQUAL(ref[i], results[i]);
  }

  rho_f->UpdateHist(A, QMState("s2"));
  std::vector<votca::Index> results2 =
      rho_f->CalcIndeces(A, QMStateType::Singlet);
  std::vector<votca::Index> ref2 = {1};
  for (votca::Index i = 0; i < votca::Index(ref2.size()); i++) {
    BOOST_CHECK_EQUAL(ref2[i], results2[i]);
  }
  BOOST_CHECK_EQUAL(results2.size(), ref2.size());
  std::unique_ptr<StateFilter_base> rho_f2 = factory.Create("overlap");

  rho_f2->Initialize(prop.get("overlap"));

  rho_f2->UpdateHist(A, QMState("pqp8"));
  std::vector<votca::Index> results3 =
      rho_f2->CalcIndeces(A, QMStateType::PQPstate);

  std::vector<votca::Index> ref3 = {8};
  BOOST_CHECK_EQUAL(results3.size(), ref3.size());
  for (votca::Index i = 0; i < votca::Index(ref3.size()); i++) {
    BOOST_CHECK_EQUAL(ref3[i], results3[i]);
  }

  // reference energy
  Eigen::VectorXd se_ref_btda = Eigen::VectorXd::Zero(3);
  se_ref_btda << 0.0887758, 0.0887758, 0.0887758;

  // reference coefficients
  Eigen::MatrixXd spsi_ref_btda =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/overlap_filter/spsi_ref_btda.mm");

  // // reference coefficients AR
  Eigen::MatrixXd spsi_ref_btda_AR =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) +
          "/overlap_filter/spsi_ref_btda_AR.mm");

  A.BSESinglets().eigenvectors() = spsi_ref_btda;
  A.BSESinglets().eigenvectors2() = spsi_ref_btda_AR;

  A.BSESinglets().eigenvalues() = se_ref_btda;
  A.setTDAApprox(false);
  A.CalcCoupledTransition_Dipoles();

  std::unique_ptr<StateFilter_base> rho_f3 = factory.Create("overlap");

  rho_f3->Initialize(prop.get("overlap"));

  rho_f3->UpdateHist(A, QMState("s1"));

  std::vector<votca::Index> ref_btda = {0};
  std::vector<votca::Index> results_btda =
      rho_f3->CalcIndeces(A, QMStateType::Singlet);

  BOOST_CHECK_EQUAL(results_btda.size(), ref_btda.size());
  for (votca::Index i = 0; i < votca::Index(ref_btda.size()); i++) {
    BOOST_CHECK_EQUAL(ref_btda[i], results_btda[i]);
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
