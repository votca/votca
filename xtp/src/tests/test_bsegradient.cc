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
#include "votca/xtp/basisset.h"
#include "votca/xtp/bse.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/excitonfollow.h"
#include "votca/xtp/orbitals.h"
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

void InitializeEngine(DFTEngine& dft, const std::string& functional) {
  dft.setLogger(&TestLog());
  std::string xml_path = "/tmp/xtp_test_bsegradient_dft.xml";
  std::ofstream xml(xml_path);
  xml << "<dftpackage>\n<spin>1</spin>\n<name>xtp</name>\n<charge>0</charge>\n";
  xml << "<functional>" << functional << "</functional>\n";
  xml << "<basisset>" << DataFile("def2-svp.xml") << "</basisset>\n";
  xml << "<auxbasisset>" << DataFile("aux-def2-svp.xml") << "</auxbasisset>\n";
  xml << "<initial_guess>independent</initial_guess>\n<xtpdft>\n";
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

}  // namespace

BOOST_AUTO_TEST_CASE(pbe0_singlet_full_windows) {
  libint2::initialize();
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{0, 23, 0, 23}, true);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_singlet_reduced_windows) {
  libint2::initialize();
  CheckGradient("XC_HYB_GGA_XC_PBEH", Windows{1, 20, 1, 22}, true);
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(pbe0_triplet_reduced_windows) {
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
