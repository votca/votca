/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
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
#include "votca/xtp/sigma_base.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE sigma_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/ppm.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigmafactory.h"
#include "votca/xtp/threecenter.h"
#include "xtp_libint2.h"
using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sigma_test)

BOOST_AUTO_TEST_CASE(sigma_full) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/sigma_ppm/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/MOs.mm");

  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(mo_energy);

  std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("ppm", Mmn, rpa);

  Sigma_base::options opt;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.rpamin = 0;
  opt.rpamax = 16;
  opt.eta = 1e-3;
  sigma->configure(opt);

  Eigen::MatrixXd x = sigma->CalcExchangeMatrix();

  Eigen::MatrixXd x_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/x_ref.mm");

  bool check_x = x_ref.isApprox(x, 1e-5);
  if (!check_x) {
    cout << "Sigma X" << endl;
    cout << x << endl;
    cout << "Sigma X ref" << endl;
    cout << x_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_x, true);

  sigma->PrepareScreening();
  Eigen::MatrixXd c = sigma->CalcCorrelationOffDiag(mo_energy);
  c.diagonal() = sigma->CalcCorrelationDiag(mo_energy);

  Eigen::MatrixXd c_ref = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/c_ref.mm");

  bool check_c_diag = c.diagonal().isApprox(c_ref.diagonal(), 1e-5);
  if (!check_c_diag) {
    cout << "Sigma C" << endl;
    cout << c.diagonal() << endl;
    cout << "Sigma C ref" << endl;
    cout << c_ref.diagonal() << endl;
  }
  BOOST_CHECK_EQUAL(check_c_diag, true);

  bool check_c = c.isApprox(c_ref, 1e-5);
  if (!check_c) {
    cout << "Sigma C" << endl;
    cout << c << endl;
    cout << "Sigma C ref" << endl;
    cout << c_ref << endl;
  }
  BOOST_CHECK_EQUAL(check_c, true);
  libint2::finalize();
}

// Ignoring core levels in the RPA (rpamin > 0) must give the same self-energy
// as including them with vanishing three-centre integrals: the levels below
// rpamin then neither screen nor enter the sum over states. This checks that
// the occupied/virtual split of the PPM self-energy uses the level numbering
// of its rpamin-based arrays.
BOOST_AUTO_TEST_CASE(sigma_rpamin_equals_decoupled_core) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/sigma_ppm/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  const Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/MOs.mm");
  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  Logger log;
  const votca::Index homo = 4;
  const votca::Index core = 1;  // levels below this are ignored

  // reference: all levels, the core level decoupled
  TCMatrix_gwbse Mmn_all;
  Mmn_all.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn_all.Fill(aobasis, aobasis, MOs);
  for (votca::Index m = 0; m < Mmn_all.msize(); ++m) {
    Mmn_all[m].topRows(core).setZero();
  }
  for (votca::Index m = 0; m < core; ++m) {
    Mmn_all[m].setZero();
  }
  RPA rpa_all(log, Mmn_all);
  rpa_all.configure(homo, 0, 16);
  rpa_all.setRPAInputEnergies(mo_energy);

  // the same with the core level left out of the RPA
  TCMatrix_gwbse Mmn_core;
  Mmn_core.Initialize(aobasis.AOBasisSize(), core, 16, core, 16);
  Mmn_core.Fill(aobasis, aobasis, MOs);
  RPA rpa_core(log, Mmn_core);
  rpa_core.configure(homo, core, 16);
  rpa_core.setRPAInputEnergies(mo_energy.segment(core, 17 - core));

  Sigma_base::options opt;
  opt.homo = homo;
  opt.qpmin = core;
  opt.qpmax = 16;
  opt.rpamax = 16;
  opt.eta = 1e-3;

  std::unique_ptr<Sigma_base> sigma_all =
      SigmaFactory().Create("ppm", Mmn_all, rpa_all);
  opt.rpamin = 0;
  sigma_all->configure(opt);
  std::unique_ptr<Sigma_base> sigma_core =
      SigmaFactory().Create("ppm", Mmn_core, rpa_core);
  opt.rpamin = core;
  sigma_core->configure(opt);

  sigma_all->PrepareScreening();
  sigma_core->PrepareScreening();
  const Eigen::VectorXd freq = mo_energy.segment(core, 17 - core);
  const Eigen::VectorXd c_all = sigma_all->CalcCorrelationDiag(freq);
  const Eigen::VectorXd c_core = sigma_core->CalcCorrelationDiag(freq);
  BOOST_CHECK_GT(c_all.cwiseAbs().maxCoeff(), 1e-3);
  BOOST_CHECK_SMALL((c_all - c_core).cwiseAbs().maxCoeff(), 1e-10);
  const Eigen::MatrixXd off_all = sigma_all->CalcCorrelationOffDiag(freq);
  const Eigen::MatrixXd off_core = sigma_core->CalcCorrelationOffDiag(freq);
  BOOST_CHECK_SMALL((off_all - off_core).cwiseAbs().maxCoeff(), 1e-10);
  libint2::finalize();
}

// The blocked matrix-product evaluation of the exchange and of the
// off-diagonal correlation must not depend on the block and chunk sizes,
// and must equal the element-wise formula.
BOOST_AUTO_TEST_CASE(pair_contraction_blocking) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/sigma_ppm/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  const Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_ppm/MOs.mm");
  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << -0.612601, -0.341755, -0.341755, -0.341755, 0.137304, 0.16678,
      0.16678, 0.16678, 0.671592, 0.671592, 0.671592, 0.974255, 1.01205,
      1.01205, 1.01205, 1.64823, 19.4429;
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);
  RPA rpa(log, Mmn);
  rpa.configure(4, 0, 16);
  rpa.setRPAInputEnergies(mo_energy);
  std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("ppm", Mmn, rpa);
  Sigma_base::options opt;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.rpamin = 0;
  opt.rpamax = 16;
  opt.eta = 1e-3;
  sigma->configure(opt);

  const Eigen::MatrixXd x_default = sigma->CalcExchangeMatrix();
  sigma->PrepareScreening();
  const Eigen::MatrixXd c_default = sigma->CalcCorrelationOffDiag(mo_energy);
  // a few levels per block and a few auxiliary functions per chunk
  sigma->setPairContractionMemory(8.0 * 17 * 17 * 3, 8.0 * 17 * 17 * 2);
  const Eigen::MatrixXd x_small = sigma->CalcExchangeMatrix();
  const Eigen::MatrixXd c_small = sigma->CalcCorrelationOffDiag(mo_energy);
  BOOST_CHECK_SMALL((x_small - x_default).cwiseAbs().maxCoeff(), 1e-12);
  BOOST_CHECK_SMALL((c_small - c_default).cwiseAbs().maxCoeff(), 1e-12);
  double maxdiff = 0.0;
  for (votca::Index i = 0; i < 17; ++i) {
    BOOST_CHECK_EQUAL(c_small(i, i), 0.0);
    for (votca::Index j = i + 1; j < 17; ++j) {
      const double element = sigma->CalcCorrelationOffDiagElement(
          i, j, mo_energy(i), mo_energy(j));
      maxdiff = std::max(maxdiff, std::abs(element - c_small(i, j)));
      maxdiff = std::max(maxdiff, std::abs(element - c_small(j, i)));
    }
  }
  BOOST_CHECK_SMALL(maxdiff, 1e-12);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
