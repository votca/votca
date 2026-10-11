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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE rpa_test

// Third party includes
#include "boost/test/unit_test.hpp"

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace votca;
using namespace std;

BOOST_AUTO_TEST_SUITE(rpa_test)

BOOST_AUTO_TEST_CASE(rpa_calcenergies) {

  Logger log;
  TCMatrix_gwbse Mmn;
  Eigen::VectorXd eigenvals;
  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 9);
  Eigen::VectorXd dftenergies = Eigen::VectorXd::Zero(10);
  dftenergies << -0.5, -0.4, -0.3, -0.2, -0.2, -0.1, 0, 0.1, 0.2, 0.3;
  Eigen::VectorXd gwenergies = Eigen::VectorXd::Zero(7);
  gwenergies << -0.15, -0.05, 0.05, 0.15, 0.45, 0.55, 0.65;
  votca::Index qpmin = 1;
  rpa.UpdateRPAInputEnergies(dftenergies, gwenergies, qpmin);
  Eigen::VectorXd rpaenergies = rpa.getRPAInputEnergies();
  Eigen::VectorXd rpaenergies_ref = Eigen::VectorXd::Zero(10);
  rpaenergies_ref << -0.85, -0.15, -0.05, 0.05, 0.15, 0.45, 0.55, 0.65, 0.75,
      0.85;
  bool e_check = rpaenergies_ref.isApprox(rpaenergies, 0.0001);

  if (!e_check) {
    cout << "energy" << endl;
    cout << rpaenergies << endl;
    cout << "energy_ref" << endl;
    cout << rpaenergies_ref << endl;
  }
  BOOST_CHECK_EQUAL(e_check, true);
}

BOOST_AUTO_TEST_CASE(rpa_full) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/rpa/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/rpa/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::VectorXd eigenvals = votca::tools::EigenIO_MatrixMarket::ReadVector(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/eigenvals.mm");

  Eigen::MatrixXd eigenvectors = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/eigenvectors.mm");
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, eigenvectors);

  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(eigenvals);
  Eigen::MatrixXd e_i = rpa.calculate_epsilon_i(0.5);

  Eigen::MatrixXd i_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/i_ref.mm");
  bool i_check = i_ref.isApprox(e_i, 0.0001);

  if (!i_check) {
    cout << "Epsilon_i" << endl;
    cout << e_i << endl;
    cout << "Epsilon_i_ref" << endl;
    cout << i_ref << endl;
  }
  BOOST_CHECK_EQUAL(i_check, 1);

  Eigen::MatrixXd e_r = rpa.calculate_epsilon_r(0.0);

  Eigen::MatrixXd r_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/rpa/r_ref.mm");
  bool r_check = r_ref.isApprox(e_r, 0.0001);

  if (!r_check) {
    cout << "Epsilon_r" << endl;
    cout << e_r << endl;
    cout << "Epsilon_r_ref" << endl;
    cout << r_ref << endl;
  }

  BOOST_CHECK_EQUAL(r_check, 1);

  Eigen::MatrixXd e_r_complex =
      rpa.calculate_epsilon_r(std::complex<double>(0.5, 0.5));

  Eigen::MatrixXd r_complex_ref =
      votca::tools::EigenIO_MatrixMarket::ReadMatrix(
          std::string(XTP_TEST_DATA_FOLDER) + "/rpa/r_complex_ref.mm");
  bool r_complex_check = r_complex_ref.isApprox(e_r_complex, 0.0001);

  if (!r_complex_check) {
    cout << "Epsilon_r_complex" << endl;
    cout << e_r_complex << endl;
    cout << "Epsilon_r_compelx_ref" << endl;
    cout << r_complex_ref << endl;
  }

  BOOST_CHECK_EQUAL(r_complex_check, 1);

  libint2::finalize();
}

// With rpamin > 0 the RPA energies are stored from rpamin on, while qpmin and
// homo are absolute level numbers.
BOOST_AUTO_TEST_CASE(rpa_calcenergies_with_rpamin) {
  Logger log;
  TCMatrix_gwbse Mmn;
  RPA rpa(log, Mmn);
  rpa.configure(4, 2, 9);  // homo 4, levels 2..9
  Eigen::VectorXd dftenergies = Eigen::VectorXd::Zero(10);
  dftenergies << -0.5, -0.4, -0.3, -0.2, -0.2, -0.1, 0, 0.1, 0.2, 0.3;
  // QP levels 3..7: occupied corrections -0.15, -0.1; virtual 0.15, 0.2, 0.15
  Eigen::VectorXd gwenergies = Eigen::VectorXd::Zero(5);
  gwenergies << -0.35, -0.3, 0.05, 0.2, 0.25;
  rpa.UpdateRPAInputEnergies(dftenergies, gwenergies, 3);
  Eigen::VectorXd rpaenergies_ref = Eigen::VectorXd::Zero(8);
  // level 2 shifted by -0.15, levels 8, 9 by +0.2
  rpaenergies_ref << -0.45, -0.35, -0.3, 0.05, 0.2, 0.25, 0.4, 0.5;
  BOOST_CHECK_SMALL(
      (rpa.getRPAInputEnergies() - rpaenergies_ref).cwiseAbs().maxCoeff(),
      1e-12);
}

// Levels outside the GW window: nearest takes the boundary correction,
// linear fits the corrections next to the boundary and extrapolates over at
// most the fitted span.
BOOST_AUTO_TEST_CASE(out_of_window_shift_modes) {
  // levels 0..11, homo 5, GW window 2..9, RPA 0..11
  Eigen::VectorXd dft = Eigen::VectorXd::LinSpaced(12, -1.1, 1.1);
  const Index homo = 5;
  const Index qpmin = 2;
  const Index qpmax = 9;
  // corrections linear in the DFT energy, separately occupied and virtual
  auto occ = [](double e) { return -0.1 + 0.2 * (e + 0.7); };
  auto virt = [](double e) { return 0.2 + 0.3 * (e - 0.7); };
  Eigen::VectorXd gw(qpmax - qpmin + 1);
  for (Index l = qpmin; l <= qpmax; ++l) {
    gw(l - qpmin) = dft(l) + ((l <= homo) ? occ(dft(l)) : virt(dft(l)));
  }
  auto run = [&](OutOfWindowShift::Mode mode, double fraction) {
    Logger log;
    TCMatrix_gwbse Mmn;
    RPA rpa(log, Mmn);
    rpa.configure(homo, 0, 11);
    OutOfWindowShift shift;
    shift.mode = mode;
    shift.fit_fraction = fraction;
    rpa.setOutOfWindowShift(shift);
    rpa.UpdateRPAInputEnergies(dft, gw, qpmin);
    return Eigen::VectorXd(rpa.getRPAInputEnergies() - dft);
  };
  const double h = dft(1) - dft(0);

  const Eigen::VectorXd nearest = run(OutOfWindowShift::Mode::Nearest, 0.3);
  BOOST_CHECK_CLOSE(nearest(0), occ(dft(2)), 1e-10);
  BOOST_CHECK_CLOSE(nearest(1), occ(dft(2)), 1e-10);
  BOOST_CHECK_CLOSE(nearest(10), virt(dft(9)), 1e-10);
  BOOST_CHECK_CLOSE(nearest(11), virt(dft(9)), 1e-10);

  // fit over 2 levels per side (span h): level 1 and 10 on the line, levels
  // 0 and 11 held at one span beyond the boundary
  const Eigen::VectorXd linear = run(OutOfWindowShift::Mode::Linear, 0.3);
  BOOST_CHECK_CLOSE(linear(1), occ(dft(1)), 1e-8);
  BOOST_CHECK_CLOSE(linear(0), occ(dft(2) - h), 1e-8);
  BOOST_CHECK_CLOSE(linear(10), virt(dft(10)), 1e-8);
  BOOST_CHECK_CLOSE(linear(11), virt(dft(9) + h), 1e-8);
  // fit over all window levels (span 3h): both levels on the line
  const Eigen::VectorXd linear_all = run(OutOfWindowShift::Mode::Linear, 1.0);
  BOOST_CHECK_CLOSE(linear_all(0), occ(dft(0)), 1e-8);
  BOOST_CHECK_CLOSE(linear_all(11), virt(dft(11)), 1e-8);

  // max: largest |correction| per side, downwards / upwards
  const Eigen::VectorXd max = run(OutOfWindowShift::Mode::Max, 0.3);
  BOOST_CHECK_CLOSE(max(0), -std::abs(occ(dft(2))), 1e-10);
  BOOST_CHECK_CLOSE(max(11), virt(dft(9)), 1e-10);
  // window levels untouched
  for (Index l = qpmin; l <= qpmax; ++l) {
    BOOST_CHECK_CLOSE(linear(l), gw(l - qpmin) - dft(l), 1e-10);
  }
}

BOOST_AUTO_TEST_SUITE_END()
