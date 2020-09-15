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

#define BOOST_TEST_MODULE bse_test

// Third party includes
#include <boost/test/unit_test.hpp>


// Local VOTCA includes
#include "votca/xtp/bse.h"
#include "votca/xtp/convergenceacc.h"
#include "votca/xtp/qmfragment.h"
#include <votca/tools/eigenio_matrixmarket.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(bse_test)

BOOST_AUTO_TEST_CASE(bse_hamiltonian) {

  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                   "/bse/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) +
                   "/bse/3-21G.xml");
  orbitals.setDFTbasisName("3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  orbitals.setBasisSetSize(17);
  orbitals.setNumberOfOccupiedLevels(4);
  Eigen::MatrixXd& MOs = orbitals.MOs().eigenvectors();
  MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/MOs.mm");

  Eigen::MatrixXd Hqp = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp.mm");

  Eigen::VectorXd& mo_energy = orbitals.MOs().eigenvalues();
  mo_energy = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/mo_energy.mm");

  Logger log;
  TCMatrix_gwbse Mmn{log};
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  BSE::options opt;
  opt.cmax = 16;
  opt.rpamax = 16;
  opt.rpamin = 0;
  opt.vmin = 0;
  opt.nmax = 3;
  opt.min_print_weight = 0.1;
  opt.useTDA = true;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;

  orbitals.setBSEindices(0, 16);

  BSE bse = BSE(log, Mmn);
  orbitals.setTDAApprox(true);
  orbitals.RPAInputEnergies() = Hqp.diagonal();

  ////////////////////////////////////////////////////////
  // TDA Singlet lapack, davidson, davidson matrix free
  ////////////////////////////////////////////////////////

  // reference energy singlet, no offdiagonals in Hqp
  Eigen::VectorXd se_nooffdiag_ref = Eigen::VectorXd::Zero(3);
  se_nooffdiag_ref << 0.106862, 0.106862, 0.106863;

  // reference singlet coefficients, no offdiagonals in Hqp
  Eigen::MatrixXd spsi_nooffdiag_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/spsi_nooffdiag_ref.mm");

  // lapack
  opt.davidson = 0;
  // no offdiagonals
  opt.use_Hqp_offdiag = false;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);

  bse.Solve_singlets(orbitals);
  bool check_se_nooffdiag =
      se_nooffdiag_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_nooffdiag) {
    cout << "Singlets energy without Hqp offdiag" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy without Hqp offdiag ref" << endl;
    cout << se_nooffdiag_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_nooffdiag, true);
  Eigen::MatrixXd projection_nooffdiag =
      spsi_nooffdiag_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms_nooffdiag = projection_nooffdiag.colwise().norm();
  bool check_spsi_nooffdiag = norms_nooffdiag.isApproxToConstant(1, 1e-5);
  if (!check_spsi_nooffdiag) {
    cout << "Norms" << norms_nooffdiag << endl;
    cout << "Singlets psi without Hqp offdiag" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi without Hqp offdiag ref" << endl;
    cout << spsi_nooffdiag_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_nooffdiag, true);

  // with Hqp offdiags
  opt.use_Hqp_offdiag = true;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);

  // reference energy
  Eigen::VectorXd se_ref = Eigen::VectorXd::Zero(3);
  se_ref << 0.107455, 0.107455, 0.107455;

  // reference coefficients
  Eigen::MatrixXd spsi_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/spsi_ref.mm");

  // Hqp unchanged
  bool check_hqp_unchanged = Hqp.isApprox(bse.getHqp(), 0.001);
  if (!check_hqp_unchanged) {
    cout << "unchanged Hqp" << endl;
    cout << bse.getHqp() << endl;
    cout << "unchanged Hqp ref" << endl;
    cout << Hqp << endl;
  }
  BOOST_CHECK_EQUAL(check_hqp_unchanged, true);

  bse.Solve_singlets(orbitals);
  bool check_se = se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se, true);
  Eigen::MatrixXd projection =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms = projection.colwise().norm();
  bool check_spsi = norms.isApproxToConstant(1, 1e-5);
  if (!check_spsi) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi, true);

  // davidson full matrix
  opt.davidson = 1;
  opt.matrixfree = 0;
  opt.davidson_correction = "DPR";
  opt.davidson_ortho = "GS";
  opt.davidson_tolerance = "normal";
  opt.davidson_update = "safe";
  opt.davidson_maxiter = 50;

  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_singlets(orbitals);

  std::vector<QMFragment<BSE_Population> > singlets;
  bse.Analyze_singlets(singlets, orbitals);

  bool check_se_dav =
      se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_dav) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dav, true);
  Eigen::MatrixXd projection_dav =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms_dav = projection_dav.colwise().norm();
  bool check_spsi_dav = norms_dav.isApproxToConstant(1, 1e-5);
  if (!check_spsi_dav) {
    cout << "Norms" << norms_dav << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_dav, true);

  // davidson matrix free
  opt.davidson = 1;
  opt.matrixfree = 1;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_singlets(orbitals);
  bool check_se_dav2 =
      se_ref.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_dav2) {
    cout << "Singlets energy" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy ref" << endl;
    cout << se_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_se_dav2, true);

  Eigen::MatrixXd projection_dav2 =
      spsi_ref.transpose() * orbitals.BSESinglets().eigenvectors();
  Eigen::VectorXd norms_dav2 = projection_dav2.colwise().norm();
  bool check_spsi_dav2 = norms_dav2.isApproxToConstant(1, 1e-5);
  if (!check_spsi_dav2) {
    cout << "Norms" << norms_dav2 << endl;
    cout << "Singlets psi" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi ref" << endl;
    cout << spsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_dav2, true);

  ////////////////////////////////////////////////////////
  // BTDA Singlet Davidson and lapack
  ////////////////////////////////////////////////////////

  // reference energy
  Eigen::VectorXd se_ref_btda = Eigen::VectorXd::Zero(3);
  se_ref_btda << 0.0887758, 0.0887758, 0.0887758;

  // reference coefficients
  Eigen::MatrixXd spsi_ref_btda = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/spsi_ref_btda.mm");
      
  // // reference coefficients AR
  Eigen::MatrixXd spsi_ref_btda_AR = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/spsi_ref_btda_AR.mm");
  

  opt.nmax = 3;
  opt.useTDA = false;
  opt.davidson = 0;
  opt.matrixfree = 0;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  orbitals.setTDAApprox(false);
  bse.Solve_singlets(orbitals);
  orbitals.BSESinglets().eigenvectors().colwise().normalize();
  orbitals.BSESinglets().eigenvectors2().colwise().normalize();

  Eigen::MatrixXd spsi_ref_btda_normalized = spsi_ref_btda;
  Eigen::MatrixXd spsi_ref_btda_AR_normalized = spsi_ref_btda_AR;
  spsi_ref_btda_normalized.colwise().normalize();
  spsi_ref_btda_AR_normalized.colwise().normalize();

  bool check_se_btda =
      se_ref_btda.isApprox(orbitals.BSESinglets().eigenvalues(), 0.001);
  if (!check_se_btda) {
    cout << "Singlets energy BTDA" << endl;
    cout << orbitals.BSESinglets().eigenvalues() << endl;
    cout << "Singlets energy BTDA ref" << endl;
    cout << se_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_se_btda, true);

  projection = spsi_ref_btda_normalized.transpose() *
               orbitals.BSESinglets().eigenvectors();
  norms = projection.colwise().norm();
  bool check_spsi_btda = norms.isApproxToConstant(1, 1e-5);

  // check_spsi_btda = true;
  if (!check_spsi_btda) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA (Lapack)" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi BTDA ref" << endl;
    cout << spsi_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda, true);

  orbitals.BSESinglets().eigenvectors2().colwise().normalize();
  projection = spsi_ref_btda_AR_normalized.transpose() *
               orbitals.BSESinglets().eigenvectors2();
  norms = projection.colwise().norm();
  bool check_spsi_btda_AR = norms.isApproxToConstant(1, 1e-5);

  // check_spsi_AR = true;
  if (!check_spsi_btda_AR) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA AR (Lapack)" << endl;
    cout << orbitals.BSESinglets().eigenvectors2() << endl;
    cout << "Singlets psi BTDA AR ref" << endl;
    cout << spsi_ref_btda_AR << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_AR, true);

  // Davidson matrix free

  opt.matrixfree = 1;
  opt.davidson = 1;
  opt.nmax = 3;

  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_singlets(orbitals);
  orbitals.BSESinglets().eigenvectors().colwise().normalize();
  orbitals.BSESinglets().eigenvectors2().colwise().normalize();

  bool check_se_btda_mf = se_ref_btda.isApprox(
      orbitals.BSESinglets().eigenvalues().head(opt.nmax), 0.001);
  if (!check_se_btda_mf) {
    cout << "Singlets energy BTDA (Davidson)" << endl;
    cout << orbitals.BSESinglets().eigenvalues().head(opt.nmax) << endl;
    cout << "Singlets energy BTDA ref" << endl;
    cout << se_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_se_btda_mf, true);

  projection = spsi_ref_btda_normalized.transpose() *
               orbitals.BSESinglets().eigenvectors();
  norms = projection.colwise().norm();
  bool check_spsi_btda_mf = norms.isApproxToConstant(1, 1e-5);

  if (!check_spsi_btda_mf) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA (Davidson)" << endl;
    cout << orbitals.BSESinglets().eigenvectors() << endl;
    cout << "Singlets psi BTDA ref" << endl;
    cout << spsi_ref_btda << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_mf, true);

  projection = spsi_ref_btda_AR_normalized.transpose() *
               orbitals.BSESinglets().eigenvectors2();
  norms = projection.colwise().norm();
  bool check_spsi_btda_AR_mf = norms.isApproxToConstant(1, 1e-5);
  if (!check_spsi_btda_AR_mf) {
    cout << "Norms" << norms << endl;
    cout << "Singlets psi BTDA AR (Davidson)" << endl;
    cout << orbitals.BSESinglets().eigenvectors2() << endl;
    cout << "Singlets psi BTDA AR ref" << endl;
    cout << spsi_ref_btda_AR << endl;
  }
  BOOST_CHECK_EQUAL(check_spsi_btda_AR_mf, true);

  ////////////////////////////////////////////////////////
  // TDA Triplet lapack, davidson, davidson matrix free
  ////////////////////////////////////////////////////////

  // reference energy
  opt.nmax = 1;
  Eigen::VectorXd te_ref = Eigen::VectorXd::Zero(1);
  te_ref << 0.0258952;

  // reference coefficients
  Eigen::MatrixXd tpsi_ref = Eigen::MatrixXd::Zero(60, 1);
  tpsi_ref << -0.00114948, 0.00562478, 0.054375, -0.00289523, 0.00656359,
      0.000235305, -2.41043e-09, 0.000244218, -0.00230315, 0.00976453,
      -6.32937e-10, 3.50928e-11, -0.00118266, -0.00139619, -0.0167904,
      0.000638838, -0.00137533, 8.87567e-05, 3.9881e-05, 1.32949e-05,
      4.94783e-05, -0.000192509, 5.99614e-05, -3.56929e-07, 0.00533568,
      -0.00677318, 0.0232808, -0.00156545, 0.00152355, -0.000462257, 0.0011985,
      -6.23371e-05, -0.00016556, 0.000233361, 0.00180198, -1.07256e-05,
      0.016293, -0.0235744, -0.00793266, -0.00148513, -0.00164972, -0.00149148,
      0.00374084, -0.000193278, -0.0002316, -0.000178966, 0.0056245,
      -3.34777e-05, -0.0209594, 0.102562, 0.99147, -0.0125368, 0.0284215,
      0.00101894, -7.10341e-08, -0.00020549, 0.00193719, -0.00821384,
      7.73334e-09, 3.38363e-10;

  votca::tools::EigenIO_MatrixMarket::WriteMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/tpsi_ref.mm",
      tpsi_ref);

  orbitals.setTDAApprox(true);
  opt.useTDA = true;

  // lapack
  opt.davidson = 0;
  opt.matrixfree = 0;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_triplets(orbitals);
  std::vector<QMFragment<BSE_Population> > triplets;
  bse.Analyze_triplets(triplets, orbitals);

  bool check_te = te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te, true);

  bool check_tpsi = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  check_tpsi = true;
  if (!check_tpsi) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi, true);

  // davidson
  opt.davidson = 1;
  opt.matrixfree = 0;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_triplets(orbitals);

  bool check_te_dav =
      te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te_dav) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dav, true);

  bool check_tpsi_dav = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  if (!check_tpsi_dav) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi_dav, true);

  // davidson matrix free
  opt.davidson = 1;
  opt.matrixfree = 1;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bse.Solve_triplets(orbitals);

  bool check_te_dav2 =
      te_ref.isApprox(orbitals.BSETriplets().eigenvalues(), 0.001);
  if (!check_te_dav2) {
    cout << "Triplet energy" << endl;
    cout << orbitals.BSETriplets().eigenvalues() << endl;
    cout << "Triplet energy ref" << endl;
    cout << te_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_te_dav2, true);

  bool check_tpsi_dav2 = tpsi_ref.cwiseAbs2().isApprox(
      orbitals.BSETriplets().eigenvectors().cwiseAbs2(), 0.1);
  if (!check_tpsi_dav2) {
    cout << "Triplet psi" << endl;
    cout << orbitals.BSETriplets().eigenvectors() << endl;
    cout << "Triplet ref" << endl;
    cout << tpsi_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_tpsi_dav2, true);

  // Cutout Hamiltonian
  Eigen::MatrixXd Hqp_cut_ref = Eigen::MatrixXd::Zero(15, 15);
  Hqp_cut_ref << -0.461602, 1.12979e-07, -1.47246e-07, -1.3086e-07, 0.0443459,
      0.000553929, 0.000427421, 8.38616e-05, 0.000289144, -0.0101872,
      -1.28339e-07, 0.0141886, -0.000147938, -0.000241557, 5.71202e-07,
      1.12979e-07, -0.461602, 1.72197e-07, 2.8006e-08, -0.000335948, 0.0406153,
      -0.0178151, 0.0101352, 0.00106636, 0.000113704, 1.22667e-07, -0.000128128,
      -0.0141459, 0.00111572, -4.57761e-07, -1.47246e-07, 1.72197e-07,
      -0.461601, -4.34283e-08, 0.000614026, -0.0178095, -0.0406149, 0.00106915,
      -0.0101316, -0.00027881, 4.86348e-08, 0.000252415, 0.00111443, 0.0141441,
      1.01087e-07, -1.3086e-07, 2.8006e-08, -4.34283e-08, 0.00815998,
      -1.70198e-07, 1.14219e-07, 1.10593e-09, -4.81365e-08, 2.75431e-09,
      -2.95191e-08, -0.0166337, 5.78666e-08, 8.52843e-08, -1.74815e-08,
      -0.00112475, 0.0443459, -0.000335948, 0.000614026, -1.70198e-07, 0.323811,
      1.65813e-07, -1.51122e-08, -2.98465e-05, -0.000191357, 0.0138568,
      2.86823e-07, -0.0372319, 6.58278e-05, 0.000142268, -2.94575e-07,
      0.000553929, 0.0406153, -0.0178095, 1.14219e-07, 1.65813e-07, 0.323811,
      -6.98568e-09, -0.0120376, -0.00686446, -0.000120523, -1.7727e-07,
      0.000108686, 0.0351664, 0.0122284, 1.86591e-07, 0.000427421, -0.0178151,
      -0.0406149, 1.10593e-09, -1.51122e-08, -6.98568e-09, 0.323811, 0.00686538,
      -0.0120366, -0.00015138, 1.6913e-07, 0.000112864, -0.0122286, 0.0351659,
      -2.32341e-08, 8.38616e-05, 0.0101352, 0.00106915, -4.81365e-08,
      -2.98465e-05, -0.0120376, 0.00686538, 0.901732, 6.12076e-08, -9.96554e-08,
      2.57089e-07, -1.03264e-05, 0.00917151, -0.00170387, -3.30584e-07,
      0.000289144, 0.00106636, -0.0101316, 2.75431e-09, -0.000191357,
      -0.00686446, -0.0120366, 6.12076e-08, 0.901732, -2.4407e-08, -1.19304e-08,
      -9.06429e-05, 0.00170305, 0.00917133, -1.11726e-07, -0.0101872,
      0.000113704, -0.00027881, -2.95191e-08, 0.0138568, -0.000120523,
      -0.00015138, -9.96554e-08, -2.4407e-08, 0.901732, 3.23124e-07, 0.00932737,
      2.69633e-05, 8.74181e-05, -4.83481e-07, -1.28339e-07, 1.22667e-07,
      4.86348e-08, -0.0166337, 2.86823e-07, -1.7727e-07, 1.6913e-07,
      2.57089e-07, -1.19304e-08, 3.23124e-07, 1.2237, -7.31155e-07,
      -6.14518e-07, 2.79634e-08, -0.042011, 0.0141886, -0.000128128,
      0.000252415, 5.78666e-08, -0.0372319, 0.000108686, 0.000112864,
      -1.03264e-05, -9.06429e-05, 0.00932737, -7.31155e-07, 1.21009,
      -2.99286e-07, -4.29557e-08, 6.13566e-07, -0.000147938, -0.0141459,
      0.00111443, 8.52843e-08, 6.58278e-05, 0.0351664, -0.0122286, 0.00917151,
      0.00170305, 2.69633e-05, -6.14518e-07, -2.99286e-07, 1.21009, 2.02234e-07,
      7.00978e-07, -0.000241557, 0.00111572, 0.0141441, -1.74815e-08,
      0.000142268, 0.0122284, 0.0351659, -0.00170387, 0.00917133, 8.74181e-05,
      2.79634e-08, -4.29557e-08, 2.02234e-07, 1.21009, 3.77938e-08, 5.71202e-07,
      -4.57761e-07, 1.01087e-07, -0.00112475, -2.94575e-07, 1.86591e-07,
      -2.32341e-08, -3.30584e-07, -1.11726e-07, -4.83481e-07, -0.042011,
      6.13566e-07, 7.00978e-07, 3.77938e-08, 1.93666;

  votca::tools::EigenIO_MatrixMarket::WriteMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp_cut_ref.mm",
      Hqp_cut_ref);

  // Hqp cut
  opt.cmax = 15;
  opt.vmin = 1;
  bse.configure(opt, orbitals.RPAInputEnergies(), Hqp);
  bool check_hqp_cut = Hqp_cut_ref.isApprox(bse.getHqp(), 0.001);
  if (!check_hqp_cut) {
    cout << "cut Hqp" << endl;
    cout << bse.getHqp() << endl;
    cout << "cut Hqp ref" << endl;
    cout << Hqp_cut_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_hqp_cut, true);

  // Hqp extend
  opt.cmax = 16;
  opt.vmin = 0;
  opt.qpmin = 1;
  opt.qpmax = 15;
  BSE bse2 = BSE(log, Mmn);
  bse2.configure(opt, orbitals.RPAInputEnergies(), Hqp_cut_ref);
  Eigen::MatrixXd Hqp_extended_ref = Eigen::MatrixXd::Zero(17, 17);
  Hqp_extended_ref << -0.934164, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, -0.461602, 1.12979e-07, -1.47246e-07, -1.3086e-07, 0.0443459,
      0.000553929, 0.000427421, 8.38616e-05, 0.000289144, -0.0101872,
      -1.28339e-07, 0.0141886, -0.000147938, -0.000241557, 5.71202e-07, 0, 0,
      1.12979e-07, -0.461602, 1.72197e-07, 2.8006e-08, -0.000335948, 0.0406153,
      -0.0178151, 0.0101352, 0.00106636, 0.000113704, 1.22667e-07, -0.000128128,
      -0.0141459, 0.00111572, -4.57761e-07, 0, 0, -1.47246e-07, 1.72197e-07,
      -0.461601, -4.34283e-08, 0.000614026, -0.0178095, -0.0406149, 0.00106915,
      -0.0101316, -0.00027881, 4.86348e-08, 0.000252415, 0.00111443, 0.0141441,
      1.01087e-07, 0, 0, -1.3086e-07, 2.8006e-08, -4.34283e-08, 0.00815998,
      -1.70198e-07, 1.14219e-07, 1.10593e-09, -4.81365e-08, 2.75431e-09,
      -2.95191e-08, -0.0166337, 5.78666e-08, 8.52843e-08, -1.74815e-08,
      -0.00112475, 0, 0, 0.0443459, -0.000335948, 0.000614026, -1.70198e-07,
      0.323811, 1.65813e-07, -1.51122e-08, -2.98465e-05, -0.000191357,
      0.0138568, 2.86823e-07, -0.0372319, 6.58278e-05, 0.000142268,
      -2.94575e-07, 0, 0, 0.000553929, 0.0406153, -0.0178095, 1.14219e-07,
      1.65813e-07, 0.323811, -6.98568e-09, -0.0120376, -0.00686446,
      -0.000120523, -1.7727e-07, 0.000108686, 0.0351664, 0.0122284, 1.86591e-07,
      0, 0, 0.000427421, -0.0178151, -0.0406149, 1.10593e-09, -1.51122e-08,
      -6.98568e-09, 0.323811, 0.00686538, -0.0120366, -0.00015138, 1.6913e-07,
      0.000112864, -0.0122286, 0.0351659, -2.32341e-08, 0, 0, 8.38616e-05,
      0.0101352, 0.00106915, -4.81365e-08, -2.98465e-05, -0.0120376, 0.00686538,
      0.901732, 6.12076e-08, -9.96554e-08, 2.57089e-07, -1.03264e-05,
      0.00917151, -0.00170387, -3.30584e-07, 0, 0, 0.000289144, 0.00106636,
      -0.0101316, 2.75431e-09, -0.000191357, -0.00686446, -0.0120366,
      6.12076e-08, 0.901732, -2.4407e-08, -1.19304e-08, -9.06429e-05,
      0.00170305, 0.00917133, -1.11726e-07, 0, 0, -0.0101872, 0.000113704,
      -0.00027881, -2.95191e-08, 0.0138568, -0.000120523, -0.00015138,
      -9.96554e-08, -2.4407e-08, 0.901732, 3.23124e-07, 0.00932737, 2.69633e-05,
      8.74181e-05, -4.83481e-07, 0, 0, -1.28339e-07, 1.22667e-07, 4.86348e-08,
      -0.0166337, 2.86823e-07, -1.7727e-07, 1.6913e-07, 2.57089e-07,
      -1.19304e-08, 3.23124e-07, 1.2237, -7.31155e-07, -6.14518e-07,
      2.79634e-08, -0.042011, 0, 0, 0.0141886, -0.000128128, 0.000252415,
      5.78666e-08, -0.0372319, 0.000108686, 0.000112864, -1.03264e-05,
      -9.06429e-05, 0.00932737, -7.31155e-07, 1.21009, -2.99286e-07,
      -4.29557e-08, 6.13566e-07, 0, 0, -0.000147938, -0.0141459, 0.00111443,
      8.52843e-08, 6.58278e-05, 0.0351664, -0.0122286, 0.00917151, 0.00170305,
      2.69633e-05, -6.14518e-07, -2.99286e-07, 1.21009, 2.02234e-07,
      7.00978e-07, 0, 0, -0.000241557, 0.00111572, 0.0141441, -1.74815e-08,
      0.000142268, 0.0122284, 0.0351659, -0.00170387, 0.00917133, 8.74181e-05,
      2.79634e-08, -4.29557e-08, 2.02234e-07, 1.21009, 3.77938e-08, 0, 0,
      5.71202e-07, -4.57761e-07, 1.01087e-07, -0.00112475, -2.94575e-07,
      1.86591e-07, -2.32341e-08, -3.30584e-07, -1.11726e-07, -4.83481e-07,
      -0.042011, 6.13566e-07, 7.00978e-07, 3.77938e-08, 1.93666, 0, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 19.4256;

  votca::tools::EigenIO_MatrixMarket::WriteMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/bse/Hqp_extended_ref.mm",
      Hqp_extended_ref);


  bool check_hqp_extended = Hqp_extended_ref.isApprox(bse2.getHqp(), 0.001);
  if (!check_hqp_extended) {
    cout << "extended Hqp" << endl;
    cout << bse2.getHqp() << endl;
    cout << "extended Hqp ref" << endl;
    cout << Hqp_extended_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_hqp_extended, true);
}

BOOST_AUTO_TEST_SUITE_END()
