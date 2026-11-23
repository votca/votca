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

// DFTEngine::BSEExcitationGradient against central finite differences of
// the TDA BSE excitation energy (BSE class, static RPA screening) with fixed
// QP corrections, for water with def2-SVP / aux-def2-SVP and PBE0. The
// corrections are arbitrary, level-dependent and different in the BSE
// Hamiltonian and the screening, so that all energy-weight terms are
// exercised.

#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE bsegradient_test

// Standard includes
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/aopotential.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/bse.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/excitonfollow.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/polarsite.h"
#include "votca/xtp/threecenter.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(bsegradient_test)

namespace {

constexpr double kAngstromPerBohr = 0.52917721090380;

std::string DataFile(const std::string& name) {
  return std::string(XTP_TEST_DATA_FOLDER) + "/diabatization/" + name;
}

QMMolecule BuildWater(const Eigen::Matrix3d& pos) {
  QMMolecule mol(" ", 0);
  const char* names[3] = {"O", "H", "H"};
  std::ostringstream xyz;
  xyz.precision(12);
  xyz << "3\n\n";
  for (Index a = 0; a < 3; ++a) {
    xyz << names[a] << " " << pos(a, 0) << " " << pos(a, 1) << " " << pos(a, 2)
        << "\n";
  }
  std::string path = "/tmp/xtp_test_bsegradient_water.xyz";
  std::ofstream out(path);
  out << xyz.str();
  out.close();
  mol.LoadFromFile(path);
  return mol;
}

Eigen::Matrix3d WaterPositions() {
  Eigen::Matrix3d pos;
  pos << 0.0, 0.0, 0.1, 0.96, 0.05, 0.0, -0.24, 0.93, -0.07;
  return pos;
}

votca::xtp::Logger& TestLog() {
  static votca::xtp::Logger log;
  return log;
}

// field: a homogeneous electric field (a.u.), H0 -= F . r (DFTEngine
// externalfield)
void InitializeEngine(DFTEngine& dft, const std::string& functional,
                      const Eigen::Vector3d& field = Eigen::Vector3d::Zero(),
                      bool forces = false) {
  dft.setLogger(&TestLog());
  std::string xml_path = "/tmp/xtp_test_bsegradient_dft.xml";
  std::ofstream xml(xml_path);
  xml << "<dftpackage>\n<spin>1</spin>\n<name>xtp</name>\n<charge>0</charge>\n";
  if (!field.isZero()) {
    xml.precision(17);
    xml << "<externalfield>" << field(0) << " " << field(1) << " " << field(2)
        << "</externalfield>\n";
  }
  xml << "<functional>" << functional << "</functional>\n";
  xml << "<basisset>" << DataFile("def2-svp.xml") << "</basisset>\n";
  xml << "<auxbasisset>" << DataFile("aux-def2-svp.xml") << "</auxbasisset>\n";
  xml << "<initial_guess>independent</initial_guess>\n<xtpdft>\n";
  if (forces) {
    xml << "<compute_forces>true</compute_forces>\n";
  }
  xml << "<screening_eps>1e-12</screening_eps>\n";
  xml << "<fock_matrix_reset>3</fock_matrix_reset>\n<convergence>\n";
  xml << "  <energy>1e-12</energy>\n  <method>DIIS</method>\n";
  xml << "  <DIIS_start>0.002</DIIS_start>\n  <ADIIS_start>0.8</ADIIS_start>\n";
  xml << "  <DIIS_length>20</DIIS_length>\n  <levelshift>0.0</levelshift>\n";
  xml << "  <levelshift_end>0.2</levelshift_end>\n";
  xml << "  <max_iterations>200</max_iterations>\n  <error>1e-10</error>\n";
  xml << "  <DIIS_maxout>false</DIIS_maxout>\n  <mixing>0.7</mixing>\n";
  xml << "  <mixing_max>0.98</mixing_max>\n  <mixing_end>0.8</mixing_end>\n";
  xml << "  <davidson_max_iter>50</davidson_max_iter>\n</convergence>\n";
  xml << "<integration_grid>fine</integration_grid>\n";
  xml << "<max_iterations>200</max_iterations>\n</xtpdft>\n</dftpackage>\n";
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML(xml_path);
  dft.Initialize(prop.get("dftpackage"));
}

struct Windows {
  Index vmin, cmax, rpamin, rpamax;
};

// Fixed corrections, absolute levels: different in the BSE Hamiltonian and
// the screening.
Eigen::VectorXd ShiftBSE(Index nmo, Index homo) {
  Eigen::VectorXd a(nmo);
  for (Index p = 0; p < nmo; ++p) {
    a(p) = (p <= homo) ? -0.04 - 0.01 * double(homo - p)
                       : 0.06 + 0.004 * double(p - homo);
  }
  return a;
}
Eigen::VectorXd ShiftRPA(Index nmo, Index homo) {
  Eigen::VectorXd a(nmo);
  for (Index p = 0; p < nmo; ++p) {
    a(p) = (p <= homo) ? -0.02 + 0.003 * double(homo - p)
                       : 0.03 + 0.002 * double(p - homo);
  }
  return a;
}

// Lowest TDA BSE state for the converged ground state in orb.
double SolveBSE(Orbitals& orb, const Windows& w, bool singlet, bool tda,
                Eigen::VectorXd* X = nullptr, Eigen::VectorXd* Y = nullptr) {
  QMMolecule mol = orb.QMAtoms();
  BasisSet dftset;
  dftset.Load(DataFile("def2-svp.xml"));
  BasisSet auxset;
  auxset.Load(DataFile("aux-def2-svp.xml"));
  AOBasis dftbasis;
  dftbasis.Fill(dftset, mol);
  AOBasis auxbasis;
  auxbasis.Fill(auxset, mol);

  const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
  const Index nmo = eps.size();
  const Index homo = orb.getHomo();
  const Eigen::VectorXd a = ShiftBSE(nmo, homo);
  const Eigen::VectorXd pi = ShiftRPA(nmo, homo);

  TCMatrix_gwbse Mmn;
  Mmn.Initialize(auxbasis.AOBasisSize(), w.rpamin, w.rpamax, w.rpamin,
                 w.rpamax);
  Mmn.Fill(auxbasis, dftbasis, orb.MOs().eigenvectors());

  BSE::options opt;
  opt.useTDA = tda;
  opt.homo = homo;
  opt.rpamin = w.rpamin;
  opt.rpamax = w.rpamax;
  opt.qpmin = w.vmin;
  opt.qpmax = w.cmax;
  opt.vmin = w.vmin;
  opt.cmax = w.cmax;
  opt.nmax = 3;
  opt.davidson_correction = "DPR";
  opt.davidson_tolerance = "lapack";
  opt.davidson_update = "safe";
  opt.davidson_maxiter = 100;
  opt.min_print_weight = 0.5;
  opt.use_Hqp_offdiag = false;
  opt.max_dyn_iter = 0;
  opt.dyn_tolerance = 1e-5;

  const Index nqp = w.cmax - w.vmin + 1;
  Eigen::MatrixXd Hqp =
      (eps.segment(w.vmin, nqp) + a.segment(w.vmin, nqp)).asDiagonal();
  Eigen::VectorXd rpa_energies =
      eps.segment(w.rpamin, w.rpamax - w.rpamin + 1) +
      pi.segment(w.rpamin, w.rpamax - w.rpamin + 1);

  orb.setBSEindices(w.vmin, w.cmax);
  orb.setTDAApprox(tda);
  orb.RPAInputEnergies() = rpa_energies;
  BSE bse(TestLog(), Mmn);
  bse.configure(opt, rpa_energies, Hqp);
  const votca::tools::EigenSystem* es;
  if (singlet) {
    bse.Solve_singlets(orb);
    es = &orb.BSESinglets();
  } else {
    bse.Solve_triplets(orb);
    es = &orb.BSETriplets();
  }
  if (X != nullptr) {
    *X = es->eigenvectors().col(0);
  }
  if (Y != nullptr && !tda) {
    *Y = es->eigenvectors2().col(0);
  }
  return es->eigenvalues()(0);
}

Orbitals RunSCF(const Eigen::Matrix3d& pos, DFTEngine& dft) {
  Orbitals orb;
  orb.QMAtoms() = BuildWater(pos);
  BOOST_REQUIRE(dft.Evaluate(orb));
  return orb;
}

void CheckGradient(const std::string& functional, const Windows& w,
                   bool singlet, bool tda = true) {
  // the gradient needs first derivatives of the one-body, 2- and 3-centre
  // integrals; libint2 builds without them throw in BSEExcitationGradient
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    return;
  }
  DFTEngine dft;
  InitializeEngine(dft, functional);
  Orbitals orb = RunSCF(WaterPositions(), dft);
  Eigen::VectorXd X;
  Eigen::VectorXd Y;
  const double omega_ref = SolveBSE(orb, w, singlet, tda, &X, &Y);

  const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
  const Index nmo = eps.size();
  const Index homo = orb.getHomo();
  DFTEngine::BSEGradientInput in;
  in.singlet = singlet;
  in.vmin = w.vmin;
  in.cmax = w.cmax;
  in.rpamin = w.rpamin;
  in.rpamax = w.rpamax;
  in.X = X;
  in.Y = Y;
  in.energies_bse = eps.segment(w.vmin, w.cmax - w.vmin + 1) +
                    ShiftBSE(nmo, homo).segment(w.vmin, w.cmax - w.vmin + 1);
  in.energies_rpa =
      eps.segment(w.rpamin, w.rpamax - w.rpamin + 1) +
      ShiftRPA(nmo, homo).segment(w.rpamin, w.rpamax - w.rpamin + 1);
  DFTEngine::BSEGradientResult res = dft.BSEExcitationGradient(orb, in);

  BOOST_CHECK_SMALL(res.omega - omega_ref, 1e-9);

  // Directional finite differences: each direction costs two SCF + BSE
  // runs instead of 18 for the full 3x3 gradient, which timed out in CI
  // Debug builds (unit_test_bsegradient > 1500 s on fedora:latest). Fixed,
  // non-symmetric displacements of all atoms, so an error in any gradient
  // component shows up; Debug builds check only the first.
  // Angstrom; smaller steps are limited by SCF and grid noise (~5e-6)
  const double h = 5e-4;
  std::vector<Eigen::Matrix3d> directions(2);
  directions[0] << 0.31, -0.52, 0.17,  //
      -0.44, 0.08, 0.61,               //
      0.23, 0.57, -0.29;
  directions[1] << -0.12, 0.35, 0.48,  //
      0.27, -0.63, 0.09,               //
      0.54, 0.14, -0.38;
#ifndef NDEBUG
  directions.resize(1);
#endif
  const double gnorm = res.gradient.norm();
  double err = 0.0;
  for (const Eigen::Matrix3d& d0 : directions) {
    const Eigen::Matrix3d d = d0 / d0.norm();
    DFTEngine dp;
    InitializeEngine(dp, functional);
    Orbitals op = RunSCF(WaterPositions() + h * d, dp);
    DFTEngine dm;
    InitializeEngine(dm, functional);
    Orbitals om = RunSCF(WaterPositions() - h * d, dm);
    // gradient rows are atoms, columns x,y,z (Hartree/Bohr), as d
    const double fd =
        (SolveBSE(op, w, singlet, tda) - SolveBSE(om, w, singlet, tda)) /
        (2.0 * h) * kAngstromPerBohr;
    const double analytic = (res.gradient.array() * d.array()).sum();
    err = std::max(err, std::abs(analytic - fd) / gnorm);
    std::cout << "  direction derivative: analytic " << analytic
              << ", finite difference " << fd << std::endl;
  }
  std::cout << functional << (singlet ? " singlet" : " triplet")
            << (tda ? " TDA" : " full BSE") << " windows " << w.vmin << "-"
            << w.cmax << " / " << w.rpamin << "-" << w.rpamax << ": Omega "
            << res.omega << ", max rel. error " << err << "\nanalytic\n"
            << res.gradient << std::endl;
  BOOST_CHECK_SMALL(res.gradient.colwise().sum().cwiseAbs().maxCoeff(), 1e-6);
  BOOST_CHECK_SMALL(err, 5e-5);
}

// Gradients of single levels (the QP-level targets) with fixed corrections:
// then they are the derivatives of the KS energies, checked against finite
// differences of the SCF eigenvalues.
void CheckLevelGradients(const std::string& functional, const Windows& w) {
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    return;
  }
  DFTEngine dft;
  InitializeEngine(dft, functional);
  Orbitals orb = RunSCF(WaterPositions(), dft);
  const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
  const Index nmo = eps.size();
  const Index homo = orb.getHomo();
  DFTEngine::BSEGradientInput in;
  in.vmin = w.vmin;
  in.cmax = w.cmax;
  in.rpamin = w.rpamin;
  in.rpamax = w.rpamax;
  in.energies_bse = eps.segment(w.vmin, w.cmax - w.vmin + 1) +
                    ShiftBSE(nmo, homo).segment(w.vmin, w.cmax - w.vmin + 1);
  in.energies_rpa =
      eps.segment(w.rpamin, w.rpamax - w.rpamin + 1) +
      ShiftRPA(nmo, homo).segment(w.rpamin, w.rpamax - w.rpamin + 1);
  in.qp_levels = {homo - 1, homo, homo + 1, homo + 3};
  const DFTEngine::BSEGradientResult res = dft.BSEExcitationGradient(orb, in);
  BOOST_REQUIRE_EQUAL(res.qp_gradients.size(), in.qp_levels.size());
  BOOST_CHECK_EQUAL(res.gradient.size(), 0);

  const double h = 5e-4;
  Eigen::Matrix3d d;
  d << 0.31, -0.52, 0.17,  //
      -0.44, 0.08, 0.61,   //
      0.23, 0.57, -0.29;
  d /= d.norm();
  DFTEngine dp;
  InitializeEngine(dp, functional);
  const Orbitals op = RunSCF(WaterPositions() + h * d, dp);
  DFTEngine dm;
  InitializeEngine(dm, functional);
  const Orbitals om = RunSCF(WaterPositions() - h * d, dm);
  for (std::size_t l = 0; l < in.qp_levels.size(); ++l) {
    const Index p = in.qp_levels[l];
    const double fd = (op.MOs().eigenvalues()(p) - om.MOs().eigenvalues()(p)) /
                      (2.0 * h) * kAngstromPerBohr;
    const Eigen::MatrixXd& grad = res.qp_gradients[l];
    const double analytic = (grad.array() * d.array()).sum();
    std::cout << "  level " << p << ": analytic " << analytic
              << ", finite difference " << fd << std::endl;
    // The finite differences carry an absolute noise of up to ~3e-6
    // Hartree/Bohr (SCF eigenvalues reproducible to ~1e-9 Hartree, depending
    // on the thread count); the HOMO derivative along d is only ~2e-4, so a
    // purely relative check fails on the noise (seen in CI with one thread).
    BOOST_CHECK_SMALL(std::abs(analytic - fd), 5e-5 * grad.norm() + 1e-5);
    BOOST_CHECK_SMALL(grad.colwise().sum().cwiseAbs().maxCoeff(), 1e-6);
  }
}

// The relaxed densities of BSEExcitationGradient: with H0 -= F . r, the
// response of the BSE excitation energy (and of single levels) to the field
// is dOmega/dF_k = -Tr(r_k P), checked against finite fields.
void CheckRelaxedDensities(const std::string& functional, const Windows& w,
                           bool singlet, bool tda) {
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    return;
  }
  DFTEngine dft;
  InitializeEngine(dft, functional);
  Orbitals orb = RunSCF(WaterPositions(), dft);
  Eigen::VectorXd X;
  Eigen::VectorXd Y;
  SolveBSE(orb, w, singlet, tda, &X, &Y);
  const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
  const Index nmo = eps.size();
  const Index homo = orb.getHomo();
  DFTEngine::BSEGradientInput in;
  in.singlet = singlet;
  in.vmin = w.vmin;
  in.cmax = w.cmax;
  in.rpamin = w.rpamin;
  in.rpamax = w.rpamax;
  in.X = X;
  in.Y = Y;
  in.energies_bse = eps.segment(w.vmin, w.cmax - w.vmin + 1) +
                    ShiftBSE(nmo, homo).segment(w.vmin, w.cmax - w.vmin + 1);
  in.energies_rpa =
      eps.segment(w.rpamin, w.rpamax - w.rpamin + 1) +
      ShiftRPA(nmo, homo).segment(w.rpamin, w.rpamax - w.rpamin + 1);
  const DFTEngine::BSEGradientResult res = dft.BSEExcitationGradient(orb, in);
  {
    // density_only: the same density without the integral derivatives
    DFTEngine::BSEGradientInput din = in;
    din.density_only = true;
    const DFTEngine::BSEGradientResult dres =
        dft.BSEExcitationGradient(orb, din);
    BOOST_CHECK_EQUAL(dres.gradient.size(), 0);
    BOOST_CHECK_SMALL(
        (dres.relaxed_density - res.relaxed_density).cwiseAbs().maxCoeff(),
        1e-12);
  }
  DFTEngine::BSEGradientInput lin = in;
  lin.X.resize(0);
  lin.Y.resize(0);
  lin.qp_levels = {homo, homo + 1};
  const DFTEngine::BSEGradientResult lres = dft.BSEExcitationGradient(orb, lin);
  BOOST_REQUIRE_EQUAL(lres.qp_relaxed_densities.size(), 2);

  BasisSet dftset;
  dftset.Load(DataFile("def2-svp.xml"));
  AOBasis dftbasis;
  dftbasis.Fill(dftset, orb.QMAtoms());
  AODipole dipole;
  dipole.setCenter(orb.QMAtoms().getPos());
  dipole.Fill(dftbasis);
  AOOverlap overlap;
  overlap.Fill(dftbasis);
  // difference densities: no change of the electron number
  BOOST_CHECK_SMALL((res.relaxed_density.cwiseProduct(overlap.Matrix())).sum(),
                    1e-8);

  const double f = 2e-4;  // a.u.
#ifdef NDEBUG
  const std::vector<Index> directions = {2, 0};
#else
  const std::vector<Index> directions = {2};  // Debug: one direction
#endif
  for (Index k : directions) {
    const Eigen::MatrixXd& rk = dipole.Matrix()[k];
    Eigen::Vector3d field = Eigen::Vector3d::Zero();
    field(k) = f;
    DFTEngine dp;
    InitializeEngine(dp, functional, field);
    Orbitals op = RunSCF(WaterPositions(), dp);
    DFTEngine dm;
    InitializeEngine(dm, functional, -field);
    Orbitals om = RunSCF(WaterPositions(), dm);
    const double fd =
        (SolveBSE(op, w, singlet, tda) - SolveBSE(om, w, singlet, tda)) /
        (2.0 * f);
    const double analytic = -(rk.cwiseProduct(res.relaxed_density)).sum();
    std::cout << "  field " << k << ": dOmega/dF analytic " << analytic
              << ", finite field " << fd << std::endl;
    BOOST_CHECK_SMALL(std::abs(analytic - fd), 1e-5 + 1e-4 * std::abs(fd));
    for (Index l = 0; l < 2; ++l) {
      const Index p = lin.qp_levels[std::size_t(l)];
      const double fdl =
          (op.MOs().eigenvalues()(p) - om.MOs().eigenvalues()(p)) / (2.0 * f);
      const double al =
          -(rk.cwiseProduct(lres.qp_relaxed_densities[std::size_t(l)])).sum();
      std::cout << "  field " << k << ": level " << p << " analytic " << al
                << ", finite field " << fdl << std::endl;
      BOOST_CHECK_SMALL(std::abs(al - fdl), 1e-5 + 1e-4 * std::abs(fdl));
    }
  }
}

// External sites around water (bohr): a static site with charge, dipole
// and quadrupole, a point polar site with an induced dipole and a
// charge-response site with Gaussian induced charge and dipole; frozen
// moments, positions shifted by h * dir (nsites x 3).
std::vector<std::unique_ptr<StaticSite>> EmbeddingSites(
    const Eigen::MatrixXd& shift = Eigen::MatrixXd::Zero(3, 3)) {
  std::vector<std::unique_ptr<StaticSite>> sites;
  {
    auto s = std::make_unique<StaticSite>(
        0, "O", Eigen::Vector3d(3.6, -1.2, 1.5) + shift.row(0).transpose());
    Vector9d Q;
    Q << 0.30, 0.10, -0.15, 0.20, 0.25, -0.10, 0.15, 0.05, -0.20;
    s->setMultipole(Q, 2);
    sites.push_back(std::move(s));
  }
  {
    auto p = std::make_unique<PolarSite>(
        1, "H", Eigen::Vector3d(-2.8, 2.6, -1.0) + shift.row(1).transpose());
    Vector9d Q = Vector9d::Zero();
    Q(0) = -0.20;
    p->setMultipole(Q, 0);
    p->setInduced_Dipole(Eigen::Vector3d(0.05, -0.08, 0.03));
    sites.push_back(std::move(p));
  }
  {
    auto c = std::make_unique<PolarSite>(
        2, "C", Eigen::Vector3d(0.5, -3.4, -2.2) + shift.row(2).transpose());
    Vector9d Q = Vector9d::Zero();
    Q(0) = 0.10;
    c->setMultipole(Q, 0);
    c->setChargeResponse(1.2, false);
    c->setInduced_Charge(-0.15);
    c->setInduced_Dipole(Eigen::Vector3d(-0.04, 0.06, 0.02));
    sites.push_back(std::move(c));
  }
  return sites;
}

// fixed, non-symmetric directions for the atoms and the sites
Eigen::Matrix3d AtomDirection() {
  Eigen::Matrix3d d;
  d << 0.31, -0.52, 0.17, -0.44, 0.08, 0.61, 0.23, 0.57, -0.29;
  return d / d.norm();
}
Eigen::Matrix3d SiteDirection() {
  Eigen::Matrix3d d;
  d << -0.12, 0.35, 0.48, 0.27, -0.63, 0.09, 0.54, 0.14, -0.38;
  return d / d.norm();
}

}  // namespace

// Debug builds (Eigen without optimisation) are 50 to 60 times slower than
// Release here: with all cases the test ran for more than 1500 s in the CI
// Debug jobs. They keep the level gradients, the two full-BSE gradients,
// which also go through the TDA terms, and the TDA part of the
// root-following test.
#ifdef NDEBUG
const bool kAllCases = true;
#else
const bool kAllCases = false;
#endif

BOOST_AUTO_TEST_CASE(pbe0_level_gradients) {
  libint2::initialize();
  CheckLevelGradients("XC_HYB_GGA_XC_PBEH", Windows{1, 12, 0, 23});
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_relaxed_densities) {
  libint2::initialize();
  CheckRelaxedDensities("XC_HYB_GGA_XC_PBEH", Windows{0, 23, 0, 23}, true,
                        false);
  libint2::finalize();
}

// Ground-state forces on the atoms and on frozen external sites (charge,
// dipole, quadrupole, induced dipole, charge-response Gaussians) against
// directional finite differences of the total energy.
BOOST_AUTO_TEST_CASE(embedded_ground_state_forces) {
  libint2::initialize();
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    libint2::finalize();
    return;
  }
  const std::string xc = "XC_HYB_GGA_XC_PBEH";
  auto energy = [&](const Eigen::Matrix3d& pos, const Eigen::MatrixXd& shift) {
    std::vector<std::unique_ptr<StaticSite>> sites = EmbeddingSites(shift);
    DFTEngine dft;
    InitializeEngine(dft, xc);
    dft.setExternalcharges(&sites);
    return RunSCF(pos, dft).getDFTTotalEnergy();
  };
  std::vector<std::unique_ptr<StaticSite>> sites = EmbeddingSites();
  DFTEngine dft;
  InitializeEngine(dft, xc, Eigen::Vector3d::Zero(), true);
  dft.setExternalcharges(&sites);
  const Orbitals orb = RunSCF(WaterPositions(), dft);
  BOOST_REQUIRE(orb.hasSiteForces());
  const Eigen::MatrixXd& F = orb.getForces();
  const Eigen::MatrixXd& Fs = orb.getSiteForces();
  // forces in Hartree/Bohr, positions here in Angstrom (water) and bohr
  // (sites)
  const double h = 5e-4;
  const Eigen::Matrix3d da = AtomDirection();
  const double fd_atoms =
      (energy(WaterPositions() + h * da, Eigen::Matrix3d::Zero()) -
       energy(WaterPositions() - h * da, Eigen::Matrix3d::Zero())) /
      (2.0 * h) * kAngstromPerBohr;
  const double an_atoms = -(F.array() * da.array()).sum();
  const double hs = 1e-3;
  const Eigen::Matrix3d ds = SiteDirection();
  const double fd_sites =
      (energy(WaterPositions(), hs * ds) - energy(WaterPositions(), -hs * ds)) /
      (2.0 * hs);
  const double an_sites = -(Fs.array() * ds.array()).sum();
  const Eigen::Vector3d net = F.colwise().sum() + Fs.colwise().sum();
  std::cout << "embedded ground state: atoms analytic " << an_atoms << " FD "
            << fd_atoms << "; sites analytic " << an_sites << " FD " << fd_sites
            << "; net force " << net.transpose() << std::endl;
  BOOST_CHECK_SMALL(std::abs(an_atoms - fd_atoms), 5e-5 * F.norm() + 1e-6);
  BOOST_CHECK_SMALL(std::abs(an_sites - fd_sites), 5e-5 * Fs.norm() + 1e-6);
  BOOST_CHECK_SMALL(net.cwiseAbs().maxCoeff(), 1e-5);
  libint2::finalize();
}

// The same for the BSE excitation energy (fixed corrections): the external
// sites enter through the relaxed density.
BOOST_AUTO_TEST_CASE(embedded_bse_gradient) {
  if (!kAllCases) {
    std::cout << "SKIPPING embedded_bse_gradient in a Debug build" << std::endl;
    return;
  }
  libint2::initialize();
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    libint2::finalize();
    return;
  }
  const std::string xc = "XC_HYB_GGA_XC_PBEH";
  const Windows w{0, 23, 0, 23};
  auto omega = [&](const Eigen::Matrix3d& pos, const Eigen::MatrixXd& shift) {
    std::vector<std::unique_ptr<StaticSite>> sites = EmbeddingSites(shift);
    DFTEngine dft;
    InitializeEngine(dft, xc);
    dft.setExternalcharges(&sites);
    Orbitals o = RunSCF(pos, dft);
    return SolveBSE(o, w, true, true);
  };
  std::vector<std::unique_ptr<StaticSite>> sites = EmbeddingSites();
  DFTEngine dft;
  InitializeEngine(dft, xc);
  dft.setExternalcharges(&sites);
  Orbitals orb = RunSCF(WaterPositions(), dft);
  Eigen::VectorXd X;
  Eigen::VectorXd Y;
  const double omega_ref = SolveBSE(orb, w, true, true, &X, &Y);
  const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
  const Index nmo = eps.size();
  const Index homo = orb.getHomo();
  DFTEngine::BSEGradientInput in;
  in.singlet = true;
  in.vmin = w.vmin;
  in.cmax = w.cmax;
  in.rpamin = w.rpamin;
  in.rpamax = w.rpamax;
  in.X = X;
  in.Y = Y;
  in.energies_bse = eps.segment(w.vmin, w.cmax - w.vmin + 1) +
                    ShiftBSE(nmo, homo).segment(w.vmin, w.cmax - w.vmin + 1);
  in.energies_rpa =
      eps.segment(w.rpamin, w.rpamax - w.rpamin + 1) +
      ShiftRPA(nmo, homo).segment(w.rpamin, w.rpamax - w.rpamin + 1);
  const DFTEngine::BSEGradientResult res = dft.BSEExcitationGradient(orb, in);
  BOOST_CHECK_SMALL(res.omega - omega_ref, 1e-9);
  BOOST_REQUIRE_EQUAL(res.site_gradient.rows(), 3);

  const double h = 5e-4;
  const Eigen::Matrix3d da = AtomDirection();
  const double fd_atoms =
      (omega(WaterPositions() + h * da, Eigen::Matrix3d::Zero()) -
       omega(WaterPositions() - h * da, Eigen::Matrix3d::Zero())) /
      (2.0 * h) * kAngstromPerBohr;
  const double an_atoms = (res.gradient.array() * da.array()).sum();
  const double hs = 1e-3;
  const Eigen::Matrix3d ds = SiteDirection();
  const double fd_sites =
      (omega(WaterPositions(), hs * ds) - omega(WaterPositions(), -hs * ds)) /
      (2.0 * hs);
  const double an_sites = (res.site_gradient.array() * ds.array()).sum();
  const Eigen::Vector3d net = res.gradient.colwise().sum().transpose() +
                              res.site_gradient.colwise().sum().transpose();
  std::cout << "embedded BSE: atoms analytic " << an_atoms << " FD " << fd_atoms
            << "; sites analytic " << an_sites << " FD " << fd_sites << "; net "
            << net.transpose() << std::endl;
  BOOST_CHECK_SMALL(std::abs(an_atoms - fd_atoms),
                    5e-5 * res.gradient.norm() + 1e-6);
  BOOST_CHECK_SMALL(std::abs(an_sites - fd_sites),
                    5e-5 * res.site_gradient.norm() + 1e-6);
  BOOST_CHECK_SMALL(net.cwiseAbs().maxCoeff(), 1e-6);
  libint2::finalize();
}

// Polarizable embedding converged self-consistently, as in the iterative
// QM/MM scheme: one polarizable site (charge, isotropic polarizability)
// whose induced dipole mu = alpha F is iterated against the field F of the
// QM density (ground state) or of the relaxed S1 density (with
// qmregion.relaxed_density). The total energy E_QM(R, mu) + mu^2/(2 alpha)
// is then stationary in mu, so the forces at frozen moments are its exact
// derivatives: checked against finite differences of the converged total
// energy (mu re-converged at each displaced geometry).
namespace {
struct PolarizableResult {
  double energy;          // converged total energy (Hartree)
  Eigen::MatrixXd force;  // analytic force at frozen mu (Hartree/Bohr)
  Eigen::Vector3d mu;
};

PolarizableResult ConvergePolarizable(const Eigen::Matrix3d& pos, bool excited,
                                      bool forces) {
  const std::string xc = "XC_HYB_GGA_XC_PBEH";
  const Windows w{0, 23, 0, 23};
  const double alpha = 8.0;
  const Eigen::Vector3d site_pos(3.2, -1.6, 1.4);
  Eigen::Vector3d mu = Eigen::Vector3d::Zero();
  PolarizableResult result;
  for (Index iter = 0; iter < 60; ++iter) {
    std::vector<std::unique_ptr<StaticSite>> sites;
    auto p = std::make_unique<PolarSite>(0, "O", site_pos);
    Vector9d Q = Vector9d::Zero();
    Q(0) = 0.4;
    p->setMultipole(Q, 0);
    p->setpolarization(alpha * Eigen::Matrix3d::Identity());
    p->setInduced_Dipole(mu);
    sites.push_back(std::move(p));
    DFTEngine dft;
    InitializeEngine(dft, xc, Eigen::Vector3d::Zero(), forces);
    dft.setExternalcharges(&sites);
    Orbitals orb = RunSCF(pos, dft);
    double energy = orb.getDFTTotalEnergy();
    Eigen::MatrixXd D = orb.DensityMatrixGroundState();
    Eigen::MatrixXd force = forces ? orb.getForces() : Eigen::MatrixXd();
    if (excited) {
      Eigen::VectorXd X;
      Eigen::VectorXd Y;
      const double omega = SolveBSE(orb, w, true, true, &X, &Y);
      energy += omega;
      const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
      const Index nmo = eps.size();
      const Index homo = orb.getHomo();
      DFTEngine::BSEGradientInput in;
      in.singlet = true;
      in.vmin = w.vmin;
      in.cmax = w.cmax;
      in.rpamin = w.rpamin;
      in.rpamax = w.rpamax;
      in.X = X;
      in.Y = Y;
      in.energies_bse =
          eps.segment(w.vmin, w.cmax - w.vmin + 1) +
          ShiftBSE(nmo, homo).segment(w.vmin, w.cmax - w.vmin + 1);
      in.energies_rpa =
          eps.segment(w.rpamin, w.rpamax - w.rpamin + 1) +
          ShiftRPA(nmo, homo).segment(w.rpamin, w.rpamax - w.rpamin + 1);
      in.density_only = !forces;
      const DFTEngine::BSEGradientResult res =
          dft.BSEExcitationGradient(orb, in);
      D += res.relaxed_density;
      if (forces) {
        force -= res.gradient;
      }
    }
    // F = -dE_QM/dmu: electrons (D) and nuclei, linear in mu
    Eigen::Vector3d field;
    {
      BasisSet bs;
      bs.Load(DataFile("def2-svp.xml"));
      AOBasis basis;
      basis.Fill(bs, orb.QMAtoms());
      eeInteractor interactor;
      for (Index c = 0; c < 3; ++c) {
        StaticSite unit(0, "X", site_pos);
        Vector9d Qu = Vector9d::Zero();
        Qu.segment<3>(1) = Eigen::Vector3d::Unit(c);
        unit.setMultipole(Qu, 1);
        std::vector<std::unique_ptr<StaticSite>> u;
        u.push_back(std::make_unique<StaticSite>(unit));
        AOMultipole m;
        m.FillPotential(basis, u);
        double dE = m.Matrix().cwiseProduct(D).sum();
        for (const QMAtom& atom : orb.QMAtoms()) {
          StaticSite nucleus(atom, double(atom.getNuccharge()));
          dE += interactor.CalcStaticEnergy_site(unit, nucleus);
        }
        field(c) = -dE;
      }
    }
    const Eigen::Vector3d mu_new = alpha * field;
    result.energy = energy + 0.5 * mu.squaredNorm() / alpha;
    result.force = force;
    result.mu = mu;
    if ((mu_new - mu).cwiseAbs().maxCoeff() < 1e-9) {
      break;
    }
    mu = mu_new;
  }
  return result;
}

void CheckPolarizableStationarity(bool excited) {
  const PolarizableResult ref =
      ConvergePolarizable(WaterPositions(), excited, true);
  const double h = 5e-4;
  const Eigen::Matrix3d da = AtomDirection();
  const double ep =
      ConvergePolarizable(WaterPositions() + h * da, excited, false).energy;
  const double em =
      ConvergePolarizable(WaterPositions() - h * da, excited, false).energy;
  const double fd = -(ep - em) / (2.0 * h) * kAngstromPerBohr;
  const double analytic = (ref.force.array() * da.array()).sum();
  std::cout << (excited ? "S1" : "ground state")
            << " in a self-consistent polarizable site: induced dipole "
            << ref.mu.transpose() << ", directional force analytic " << analytic
            << ", FD " << fd << std::endl;
  BOOST_CHECK_SMALL(std::abs(analytic - fd), 5e-5 * ref.force.norm() + 1e-6);
}
}  // namespace

BOOST_AUTO_TEST_CASE(polarizable_embedding_stationarity) {
  if (!kAllCases) {
    std::cout << "SKIPPING polarizable_embedding_stationarity in a Debug build"
              << std::endl;
    return;
  }
  libint2::initialize();
  if (!HasLibint2DerivativeSupport()) {
    std::cout << "SKIPPING: this libint2 build lacks derivative integrals"
              << std::endl;
    libint2::finalize();
    return;
  }
  CheckPolarizableStationarity(false);
  CheckPolarizableStationarity(true);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_singlet_full_windows) {
  if (!kAllCases) {
    std::cout << "SKIPPING pbe0_singlet_full_windows in a Debug build"
              << std::endl;
    return;
  }
  libint2::initialize();
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{0, 23, 0, 23}, true);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_singlet_reduced_windows) {
  if (!kAllCases) {
    std::cout << "SKIPPING pbe0_singlet_reduced_windows in a Debug build"
              << std::endl;
    return;
  }
  libint2::initialize();
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{1, 20, 1, 22}, true);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_triplet_reduced_windows) {
  if (!kAllCases) {
    std::cout << "SKIPPING pbe0_triplet_reduced_windows in a Debug build"
              << std::endl;
    return;
  }
  libint2::initialize();
  // core level outside the BSE window, highest virtuals outside both
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{1, 18, 0, 21}, false);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_singlet_full_bse) {
  libint2::initialize();
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{1, 20, 0, 22}, true, false);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_triplet_full_bse) {
  libint2::initialize();
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{0, 23, 0, 23}, false, false);
  libint2::finalize();
}

// Root following (ExcitonOverlaps): unit overlap with itself and zero with
// the other states at the same geometry; the mixed AO overlap against the
// joint basis of two geometries; the same state stays dominant after a
// small translation and a small distortion.
BOOST_AUTO_TEST_CASE(exciton_overlaps_root_following) {
  libint2::initialize();
  const Windows w{0, 23, 0, 23};
  for (bool tda : {true, false}) {
    if (!tda && !kAllCases) {
      continue;
    }
    DFTEngine d0;
    InitializeEngine(d0, "XC_HYB_GGA_XC_PBEH");
    Orbitals o0 = RunSCF(WaterPositions(), d0);
    SolveBSE(o0, w, true, tda);
    for (Index k = 0; k < 3; ++k) {
      const Eigen::VectorXd ov =
          ExcitonOverlaps(o0, QMState(QMStateType::Singlet, k, false), o0,
                          QMStateType::Singlet);
      for (Index i = 0; i < ov.size(); ++i) {
        BOOST_CHECK_SMALL(ov(i) - (i == k ? 1.0 : 0.0), 1e-8);
      }
    }
    Eigen::Matrix3d shifted = WaterPositions();
    shifted.rowwise() += Eigen::RowVector3d(0.03, -0.02, 0.01);
    DFTEngine d1;
    InitializeEngine(d1, "XC_HYB_GGA_XC_PBEH");
    Orbitals o1 = RunSCF(shifted, d1);
    SolveBSE(o1, w, true, tda);
    {
      // the mixed overlap is the off-diagonal block of the overlap of the
      // basis on both geometries together
      QMMolecule both(" ", 0);
      for (Index i = 0; i < 3; ++i) {
        both.push_back(
            QMAtom(i, o0.QMAtoms()[i].getElement(), o0.QMAtoms()[i].getPos()));
      }
      for (Index i = 0; i < 3; ++i) {
        both.push_back(QMAtom(3 + i, o1.QMAtoms()[i].getElement(),
                              o1.QMAtoms()[i].getPos()));
      }
      BasisSet bs;
      bs.Load(DataFile("def2-svp.xml"));
      AOBasis basis;
      basis.Fill(bs, both);
      AOOverlap joint;
      joint.Fill(basis);
      const Index n = o0.getBasisSetSize();
      BOOST_CHECK_SMALL(
          (MixedAOOverlap(o0, o1) - joint.Matrix().block(0, n, n, n))
              .cwiseAbs()
              .maxCoeff(),
          1e-12);
    }
    Eigen::Matrix3d bent = WaterPositions();
    bent(1, 0) += 0.03;
    bent(2, 1) -= 0.02;
    DFTEngine d2;
    InitializeEngine(d2, "XC_HYB_GGA_XC_PBEH");
    Orbitals o2 = RunSCF(bent, d2);
    SolveBSE(o2, w, true, tda);
    for (Index k = 0; k < 3; ++k) {
      const QMState sk(QMStateType::Singlet, k, false);
      const Eigen::VectorXd t =
          ExcitonOverlaps(o0, sk, o1, QMStateType::Singlet);
      const Eigen::VectorXd b =
          ExcitonOverlaps(o0, sk, o2, QMStateType::Singlet);
      std::cout << (tda ? "TDA" : "full BSE") << " S" << k + 1
                << " overlaps, translated: " << t.transpose()
                << "  distorted: " << b.transpose() << std::endl;
      // the true overlap drops with any displacement (core-like parts of
      // the AOs), but the state stays the dominant one
      Index best = 0;
      t.maxCoeff(&best);
      BOOST_CHECK_EQUAL(best, k);
      BOOST_CHECK_GT(t(k), 0.9);
      b.maxCoeff(&best);
      BOOST_CHECK_EQUAL(best, k);
      BOOST_CHECK_GT(b(k), 0.9);
    }
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
