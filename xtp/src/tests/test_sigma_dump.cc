/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
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

#define BOOST_TEST_MODULE sigma_dump_test

// Standard includes
#include <cmath>
#include <vector>

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
using votca::Index;

namespace {

// methane, 3-21G (the sigma_exact test data); levels 0..16, homo 4. The
// QP window starts above the lowest level, so the offsets between QP and
// RPA numbering are exercised.
struct Setup {
  Orbitals orbitals;
  AOBasis aobasis;
  Eigen::VectorXd mo_energy = Eigen::VectorXd::Zero(17);
  Eigen::MatrixXd MOs;
  Setup() {
    orbitals.QMAtoms().LoadFromFile(std::string(XTP_TEST_DATA_FOLDER) +
                                    "/sigma_exact/molecule.xyz");
    BasisSet basis;
    basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/3-21G.xml");
    aobasis.Fill(basis, orbitals.QMAtoms());
    mo_energy << -0.6, -0.3, -0.29, -0.28, -0.27, 0.06, 0.09, 0.09, 0.10, 0.59,
        0.66, 0.66, 0.77, 1.69, 1.98, 1.98, 2.51;
    MOs = votca::tools::EigenIO_MatrixMarket::ReadMatrix(
        std::string(XTP_TEST_DATA_FOLDER) + "/sigma_exact/MOs.mm");
  }
  Sigma_base::options Options() const {
    Sigma_base::options opt;
    opt.homo = 4;
    opt.qpmin = 1;
    opt.qpmax = 12;
    opt.rpamin = 0;
    opt.rpamax = 16;
    opt.eta = 1e-3;
    return opt;
  }
};

// Sigma_c from the pole form, term by term
double PoleSum(const Sigma_base::DiagPoleForm& pf, double w) {
  const ExactPoles& p = pf.poles;
  const Index B = p.nlevels();
  const double eta2 = p.eta * p.eta;
  double sum = 0.0;
  for (Index s = 0; s < p.nmodes(); ++s) {
    const double c =
        pf.prefactor * (pf.mode_factors ? pf.mode_factors[s] : 1.0);
    for (Index m = 0; m < B; ++m) {
      const double r = pf.residues[m + B * s];
      const double z =
          p.energies(m) + ((m < p.n_occ) ? -p.omegas(s) : p.omegas(s));
      const double t = w - z;
      sum += c * r * r * t / (t * t + eta2);
    }
  }
  return sum;
}

void CheckPoleForm(const Sigma_base& sigma, Index qptotal) {
  BOOST_REQUIRE(sigma.HasPoleForm());
  for (Index n = 0; n < qptotal; ++n) {
    const Sigma_base::DiagPoleForm pf = sigma.PoleForm(n);
    for (double w : {-0.7, -0.25, 0.08, 0.4, 1.3}) {
      const double ref = sigma.CalcCorrelationDiagElement(n, w);
      BOOST_CHECK_CLOSE_FRACTION(PoleSum(pf, w), ref, 1e-10);
    }
  }
}

// Sums over sets of degenerate static weights: inside such a set the frame
// vectors are only fixed up to a rotation, so only sums are comparable.
Eigen::VectorXd GroupSums(const Eigen::VectorXd& w, const Eigen::VectorXd& x) {
  std::vector<double> sums;
  Index i = 0;
  while (i < w.size()) {
    Index j = i + 1;
    while (j < w.size() && std::abs(w(j) - w(i)) < 1e-7 * (1.0 + w(i))) {
      ++j;
    }
    sums.push_back(x.segment(i, j - i).sum());
    i = j;
  }
  return Eigen::Map<Eigen::VectorXd>(sums.data(), Index(sums.size()));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(sigma_dump_test)

BOOST_AUTO_TEST_CASE(pole_form_and_screening) {
  libint2::initialize();
  Setup su;
  Logger log;
  const Sigma_base::options opt = su.Options();
  const Index qptotal = opt.qpmax - opt.qpmin + 1;
  Eigen::VectorXd nu(4);
  nu << 0.0, 0.2, 0.5, 2.0;

  // exact
  TCMatrix_gwbse Mmn;
  Mmn.Initialize(su.aobasis.AOBasisSize(), opt.rpamin, opt.rpamax, opt.rpamin,
                 opt.rpamax);
  Mmn.Fill(su.aobasis, su.aobasis, su.MOs);
  RPA rpa(log, Mmn);
  rpa.configure(opt.homo, opt.rpamin, opt.rpamax);
  rpa.setRPAInputEnergies(su.mo_energy);
  std::unique_ptr<Sigma_base> exact = SigmaFactory().Create("exact", Mmn, rpa);
  exact->configure(opt);
  exact->setKeepScreeningData(true);
  exact->PrepareScreening();
  CheckPoleForm(*exact, qptotal);
  const Sigma_base::ScreeningData sx = exact->ScreeningDump(nu);

  // reference: eps(0) and eps(i nu) from the RPA in the same aux frame
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      rpa.calculate_epsilon_r(0.0));
  const Eigen::MatrixXd phi = es.eigenvectors();
  const Eigen::VectorXd w_ref = 1.0 - es.eigenvalues().array().inverse();
  BOOST_TEST_MESSAGE(
      "max |w - w_ref| = " << (sx.frame_weights - w_ref).cwiseAbs().maxCoeff());
  BOOST_CHECK_SMALL((sx.frame_weights - w_ref).cwiseAbs().maxCoeff(), 1e-7);
  // moment k = 0 is the static weight
  BOOST_CHECK_SMALL((sx.frame_moments.col(0) - w_ref).cwiseAbs().maxCoeff(),
                    1e-7);
  for (Index j = 0; j < nu.size(); ++j) {
    const Eigen::MatrixXd eps_inv = rpa.calculate_epsilon_i(nu(j)).inverse();
    const Eigen::VectorXd ref =
        (phi.transpose() * eps_inv * phi).diagonal().array() - 1.0;
    BOOST_CHECK_SMALL(
        (GroupSums(w_ref, sx.frame_eps_inv_imag.col(j)) - GroupSums(w_ref, ref))
            .cwiseAbs()
            .maxCoeff(),
        1e-8);
  }
  // the mode strengths add up to -Tr[eps^-1(0) - 1]
  BOOST_CHECK_CLOSE_FRACTION(sx.strength.sum(), w_ref.sum(), 1e-7);
  BOOST_CHECK(sx.energies.isApprox(su.mo_energy));

  // plasmon pole: same frame and static weights; the fit point nu = 0.5
  // reproduces the exact eps^-1 diagonal
  TCMatrix_gwbse Mmn2;
  Mmn2.Initialize(su.aobasis.AOBasisSize(), opt.rpamin, opt.rpamax, opt.rpamin,
                  opt.rpamax);
  Mmn2.Fill(su.aobasis, su.aobasis, su.MOs);
  RPA rpa2(log, Mmn2);
  rpa2.configure(opt.homo, opt.rpamin, opt.rpamax);
  rpa2.setRPAInputEnergies(su.mo_energy);
  std::unique_ptr<Sigma_base> ppm = SigmaFactory().Create("ppm", Mmn2, rpa2);
  ppm->configure(opt);
  ppm->PrepareScreening();
  CheckPoleForm(*ppm, qptotal);
  const Sigma_base::ScreeningData sp = ppm->ScreeningDump(nu);
  // modes below the PPM weight cut are zeroed there
  Eigen::VectorXd w_cut = w_ref;
  Eigen::VectorXd x_cut = sx.frame_eps_inv_imag.col(2);
  for (Index i = 0; i < w_ref.size(); ++i) {
    if (w_ref(i) < 1e-5) {
      BOOST_CHECK_EQUAL(sp.frame_weights(i), 0.0);
      w_cut(i) = 0.0;
      x_cut(i) = 0.0;
    }
  }
  BOOST_CHECK_SMALL((sp.frame_weights - w_cut).cwiseAbs().maxCoeff(), 1e-8);
  BOOST_CHECK_SMALL((sp.frame_moments.col(0) - w_cut).cwiseAbs().maxCoeff(),
                    1e-8);
  BOOST_CHECK_SMALL(
      (GroupSums(w_ref, sp.frame_eps_inv_imag.col(2)) - GroupSums(w_ref, x_cut))
          .cwiseAbs()
          .maxCoeff(),
      1e-7);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
