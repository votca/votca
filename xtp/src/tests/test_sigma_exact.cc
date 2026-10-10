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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE sigma_test

// Standard includes
#include <algorithm>
#include <cmath>
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/aobasis.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/sigmafactory.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(sigma_test)

BOOST_AUTO_TEST_CASE(sigma_full) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/sigma_exact/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/3-21G.xml");

  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());

  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << 0.0468207, 0.0907801, 0.0907801, 0.104563, 0.592491, 0.663355,
      0.663355, 0.768373, 1.69292, 1.97724, 1.97724, 2.50877, 2.98732, 3.4418,
      3.4418, 4.81084, 17.1838;

  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/MOs.mm");

  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);

  RPA rpa(log, Mmn);
  rpa.setRPAInputEnergies(mo_energy);
  rpa.configure(4, 0, 16);
  std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("exact", Mmn, rpa);

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
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/x_ref.mm");

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
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/c_ref.mm");

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

// The batched evaluation (far poles by Chebyshev interpolation) and the
// blocked off-diagonal products agree with term-by-term sums.
BOOST_AUTO_TEST_CASE(sigma_exact_batched_equals_direct) {
  libint2::initialize();
  Orbitals orbitals;
  orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                  "/sigma_exact/molecule.xyz");
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, orbitals.QMAtoms());
  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  mo_energy << 0.0468207, 0.0907801, 0.0907801, 0.104563, 0.592491, 0.663355,
      0.663355, 0.768373, 1.69292, 1.97724, 1.97724, 2.50877, 2.98732, 3.4418,
      3.4418, 4.81084, 17.1838;
  Eigen::MatrixXd MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/MOs.mm");
  Logger log;
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(aobasis.AOBasisSize(), 0, 16, 0, 16);
  Mmn.Fill(aobasis, aobasis, MOs);
  RPA rpa(log, Mmn);
  rpa.setRPAInputEnergies(mo_energy);
  rpa.configure(4, 0, 16);
  std::unique_ptr<Sigma_base> sigma = SigmaFactory().Create("exact", Mmn, rpa);
  Sigma_base::options opt;
  opt.homo = 4;
  opt.qpmin = 0;
  opt.qpmax = 16;
  opt.rpamin = 0;
  opt.rpamax = 16;
  opt.eta = 1e-3;
  sigma->configure(opt);
  sigma->PrepareScreening();

  // grids like the QP scans: 1001 points of 0.001 Ha around each level, and
  // a wide one of 3001 points over 6 Ha
  for (votca::Index level = 0; level < 17; ++level) {
    for (double halfwidth : {0.5, 3.0, 0.75}) {
      // the last like the coarse shell scan of the QP search (63 points)
      const votca::Index n =
          (halfwidth == 0.75) ? 63 : ((halfwidth < 1.0) ? 1001 : 3001);
      const Eigen::VectorXd w = Eigen::VectorXd::LinSpaced(
          n, mo_energy(level) - halfwidth, mo_energy(level) + halfwidth);
      const Eigen::VectorXd batched =
          sigma->CalcCorrelationDiagElements(level, w);
      // reference: the same sums term by term in chunks too small for the
      // Chebyshev split; and the single-frequency sums, which round
      // t = w - z differently: near a pole (|t| ~ eta) that alone moves a
      // term by up to r^2 eps |w| / eta^2
      double cheb = 0.0;    // batched - chunked
      double single = 0.0;  // batched - single
      double scale = 0.0;
      for (votca::Index i0 = 0; i0 < n; i0 += 16) {
        const votca::Index len = std::min<votca::Index>(16, n - i0);
        const Eigen::VectorXd chunk =
            sigma->CalcCorrelationDiagElements(level, w.segment(i0, len));
        for (votca::Index i = 0; i < len; ++i) {
          const double direct =
              sigma->CalcCorrelationDiagElement(level, w(i0 + i));
          cheb = std::max(cheb, std::abs(batched(i0 + i) - chunk(i)));
          single = std::max(single, std::abs(batched(i0 + i) - direct));
          scale = std::max(scale, std::abs(direct));
        }
      }
      BOOST_CHECK_SMALL(cheb, 1e-13 * std::max(1.0, scale));
      BOOST_CHECK_SMALL(single, 1e-11 * std::max(1.0, scale));
    }
  }
  // derivative against a central difference
  for (votca::Index level : {3, 4, 9}) {
    const double w0 = mo_energy(level) + 0.0123;
    const double h = 1e-6;
    const double fd = (sigma->CalcCorrelationDiagElement(level, w0 + h) -
                       sigma->CalcCorrelationDiagElement(level, w0 - h)) /
                      (2 * h);
    BOOST_CHECK_CLOSE(sigma->CalcCorrelationDiagElementDerivative(level, w0),
                      fd, 1e-4);
  }
  // off-diagonal: blocked products against the pair loop of the base class
  const Eigen::VectorXd w = mo_energy.array() + 0.01;
  const Eigen::MatrixXd blocked = sigma->CalcCorrelationOffDiag(w);
  const Eigen::MatrixXd pairs = sigma->Sigma_base::CalcCorrelationOffDiag(w);
  BOOST_CHECK_SMALL((blocked - pairs).cwiseAbs().maxCoeff(),
                    1e-12 * std::max(1.0, pairs.cwiseAbs().maxCoeff()));
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
