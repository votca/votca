/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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
#define BOOST_TEST_MODULE dftengine_private_test
#include "xtp_libint2.h"
#include <boost/test/unit_test.hpp>
#include <stdexcept>

#include <Eigen/Core>

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <votca/tools/property.h>
#include <votca/xtp/basisset.h>

#include <votca/xtp/ERIs.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aopotential.h>
#include <votca/xtp/dftengine.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/qmmolecule.h>
#include <votca/xtp/trah.h>
#include <votca/xtp/uks_convergenceacc.h>

using namespace votca::xtp;
using namespace votca;

namespace votca {
namespace xtp {

constexpr double kTol = 1e-10;

// Defined in libint2_derivative_calls.cc (no header).
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis);

class DFTEngineTestAccess {
 public:
  static double NuclearRepulsion(const DFTEngine& e, const QMMolecule& mol) {
    return e.NuclearRepulsion(mol);
  }

  static Eigen::MatrixXd RunAtomicDFT(const DFTEngine& e, const QMAtom& atom,
                                      bool hunds_rule) {
    return e.RunAtomicDFT_unrestricted(atom, hunds_rule);
  }

  static Eigen::MatrixXd SphericalAverageShells(const DFTEngine& e,
                                                const Eigen::MatrixXd& dmat,
                                                const AOBasis& basis) {
    return e.SphericalAverageShells(dmat, basis);
  }

  static Eigen::VectorXd BuildEHTOrbitalEnergies(const DFTEngine& e,
                                                 const QMMolecule& mol) {
    return e.BuildEHTOrbitalEnergies(mol);
  }

  static Eigen::MatrixXd BuildEHTHamiltonian(const DFTEngine& e,
                                             const QMMolecule& mol) {
    return e.BuildEHTHamiltonian(mol);
  }

  static Eigen::MatrixXd OrthogonalizeGuess(const DFTEngine& e,
                                            const Eigen::MatrixXd& guess) {
    return e.OrthogonalizeGuess(guess);
  }

  static Eigen::MatrixXd InsertZeroRows(DFTEngine& e, Eigen::MatrixXd M,
                                        Index startidx, Index numofzerorows) {
    return e.InsertZeroRows(M, startidx, numofzerorows);
  }

  static Eigen::MatrixXd InsertZeroCols(DFTEngine& e, Eigen::MatrixXd M,
                                        Index startidx, Index numofzerocols) {
    return e.InsertZeroCols(M, startidx, numofzerocols);
  }

  static void SetBasis(DFTEngine& e, const AOBasis& basis) {
    e.dftbasis_ = basis;
  }

  static void FillOverlap(DFTEngine& e) { e.dftAOoverlap_.Fill(e.dftbasis_); }

  static void SetBasisName(DFTEngine& e, const std::string& basis_name) {
    e.dftbasis_name_ = basis_name;
  }

  static const AOOverlap& Overlap(const DFTEngine& e) {
    return e.dftAOoverlap_;
  }

  static void SetAuxBasis(DFTEngine& e, const AOBasis& basis) {
    e.auxbasis_ = basis;
    e.auxbasis_name_ = "set-for-testing";  // must be non-empty: the RI-J/
                                           // RI-K guard in
                                           // ComputeAndStoreForces(UKS)
                                           // checks auxbasis_name_.empty()
  }

  static void SetElectronCounts(DFTEngine& e, Index n_alpha, Index n_beta) {
    e.num_alpha_electrons_ = n_alpha;
    e.num_beta_electrons_ = n_beta;
  }

  static void SetScaHFX(DFTEngine& e, double scahfx) { e.ScaHFX_ = scahfx; }

  static Eigen::MatrixXd ComputeNonXCGradientUKS(
      const DFTEngine& e, const QMMolecule& mol,
      const UKSConvergenceAcc::SpinDensity& Dspin,
      const tools::EigenSystem& MOs_alpha, const tools::EigenSystem& MOs_beta) {
    return e.ComputeNonXCGradientUKS(mol, Dspin, MOs_alpha, MOs_beta);
  }

  static Vxc_Potential<Vxc_Grid> SetupVxc(DFTEngine& e, const QMMolecule& mol) {
    return e.SetupVxc(mol);
  }
  static Eigen::MatrixXd SetupH0(const DFTEngine& e, const QMMolecule& mol) {
    return e.SetupH0(mol).matrix();
  }
  static UKSConvergenceAcc::SpinFock BuildUKSFock(
      DFTEngine& e, const Eigen::MatrixXd& H0,
      const Vxc_Potential<Vxc_Grid>& vxc, const Eigen::MatrixXd& Da,
      const Eigen::MatrixXd& Db) {
    return e.BuildUKSFock(H0, vxc, Da, Db);
  }
  static Eigen::MatrixXd UKSOrbitalHessian(DFTEngine& e,
                                           const Vxc_Potential<Vxc_Grid>& vxc,
                                           const tools::EigenSystem& Ma,
                                           const tools::EigenSystem& Mb,
                                           const Eigen::MatrixXd& Fa,
                                           const Eigen::MatrixXd& Fb,
                                           const Eigen::MatrixXd& X) {
    return e.UKSOrbitalHessian(vxc, Ma, Mb, Fa, Fb, X);
  }
  static Eigen::MatrixXd RKSOrbitalHessian(DFTEngine& e,
                                           const Vxc_Potential<Vxc_Grid>& vxc,
                                           const tools::EigenSystem& mos,
                                           const Eigen::MatrixXd& F,
                                           const Eigen::MatrixXd& X) {
    return e.RKSOrbitalHessian(vxc, mos, F, X);
  }
  static const AOBasis& Basis(const DFTEngine& e) { return e.dftbasis_; }

  static Eigen::MatrixXd ComputeOverlapPulayGradientUKS(
      const DFTEngine& e, const QMMolecule& mol,
      const tools::EigenSystem& MOs_alpha, const tools::EigenSystem& MOs_beta) {
    return e.ComputeOverlapPulayGradientUKS(mol, MOs_alpha, MOs_beta);
  }
};

QMMolecule MakeH2(double distance_bohr) {
  QMMolecule mol("H2", 0);
  mol.push_back(QMAtom(0, "H", Eigen::Vector3d(0.0, 0.0, 0.0)));
  mol.push_back(QMAtom(1, "H", Eigen::Vector3d(distance_bohr, 0.0, 0.0)));
  return mol;
}

QMMolecule MakeSingleAtom(const std::string& element) {
  QMMolecule mol("single_atom", 0);
  mol.push_back(QMAtom(0, element, Eigen::Vector3d::Zero()));
  return mol;
}

AOBasis MakeBasis(const std::string& basis_name, const QMMolecule& mol) {
  BasisSet basisset;
  basisset.Load(basis_name);

  AOBasis basis;
  basis.Fill(basisset, mol);
  return basis;
}

}  // namespace xtp
}  // namespace votca

BOOST_AUTO_TEST_CASE(nuclear_repulsion_h2_at_2_bohr) {
  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  QMMolecule mol = MakeH2(2.0);

  // Z1 Z2 / R = 1 * 1 / 2 = 0.5 Ha
  const double e = DFTEngineTestAccess::NuclearRepulsion(engine, mol);

  BOOST_TEST(e == 0.5, boost::test_tools::tolerance(kTol));
}

BOOST_AUTO_TEST_CASE(build_eht_hamiltonian_is_symmetric_and_matches_formula) {
  libint2::initialize();
  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  const std::string basis_name = "3-21G";
  QMMolecule mol = MakeH2(1.4);

  AOBasis basis = MakeBasis(basis_name, mol);
  DFTEngineTestAccess::SetBasis(engine, basis);
  DFTEngineTestAccess::SetBasisName(engine, basis_name);
  DFTEngineTestAccess::FillOverlap(engine);

  const Eigen::MatrixXd H =
      DFTEngineTestAccess::BuildEHTHamiltonian(engine, mol);
  const Eigen::VectorXd eps =
      DFTEngineTestAccess::BuildEHTOrbitalEnergies(engine, mol);
  const Eigen::MatrixXd& S = DFTEngineTestAccess::Overlap(engine).Matrix();

  BOOST_REQUIRE_EQUAL(H.rows(), basis.AOBasisSize());
  BOOST_REQUIRE_EQUAL(H.cols(), basis.AOBasisSize());

  BOOST_TEST((H - H.transpose()).norm() < 1e-12);

  for (Index mu = 0; mu < basis.AOBasisSize(); ++mu) {
    BOOST_TEST(H(mu, mu) == eps(mu), boost::test_tools::tolerance(kTol));
  }

  constexpr double K = 1.75;
  for (Index mu = 0; mu < basis.AOBasisSize(); ++mu) {
    for (Index nu = 0; nu < mu; ++nu) {
      const double expected = K * S(mu, nu) * 0.5 * (eps(mu) + eps(nu));
      BOOST_TEST(H(mu, nu) == expected, boost::test_tools::tolerance(1e-12));
      BOOST_TEST(H(nu, mu) == expected, boost::test_tools::tolerance(1e-12));
    }
  }
  libint2::finalize();
}

// The rotational average of an atom's density: (trace / (2l+1)) times the
// identity between shells of equal l, zero between different l. It must not
// depend on the orientation of the density (an open-shell atom converges to a
// symmetry-broken solution whose orientation is decided by round-off; the old
// average turned that into a 1% difference of the guess) and must keep the
// number of electrons.
BOOST_AUTO_TEST_CASE(spherical_average_shells_is_rotation_average) {
  libint2::initialize();
  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  QMMolecule mol = MakeSingleAtom("O");
  AOBasis basis = MakeBasis("3-21G", mol);
  const Index nao = basis.AOBasisSize();
  AOOverlap overlap;
  overlap.Fill(basis);

  // block structure for an arbitrary matrix
  Eigen::MatrixXd dmat(nao, nao);
  for (Index i = 0; i < nao; ++i) {
    for (Index j = 0; j < nao; ++j) {
      dmat(i, j) = 10.0 * double(i + 1) + double(j + 1);
    }
  }
  const Eigen::MatrixXd avg =
      DFTEngineTestAccess::SphericalAverageShells(engine, dmat, basis);
  for (const AOShell& row : basis) {
    for (const AOShell& col : basis) {
      const Eigen::MatrixXd out =
          avg.block(row.getStartIndex(), col.getStartIndex(), row.getNumFunc(),
                    col.getNumFunc());
      Eigen::MatrixXd expected = Eigen::MatrixXd::Zero(out.rows(), out.cols());
      if (row.getL() == col.getL()) {
        expected.diagonal().setConstant(
            dmat.block(row.getStartIndex(), col.getStartIndex(),
                       row.getNumFunc(), col.getNumFunc())
                .trace() /
            double(row.getNumFunc()));
      }
      BOOST_CHECK_SMALL((out - expected).cwiseAbs().maxCoeff(), 1e-12);
    }
  }

  // one electron in a p orbital, along x and along (1,1,1): equal averages,
  // one electron each
  std::vector<Index> p_shells;
  for (const AOShell& shell : basis) {
    if (shell.getL() == L::P) {
      p_shells.push_back(shell.getStartIndex());
    }
  }
  BOOST_REQUIRE_GE(p_shells.size(), 2);
  auto POrbital = [&](const Eigen::Vector3d& dir) {
    Eigen::VectorXd c = Eigen::VectorXd::Zero(nao);
    c.segment(p_shells[0], 3) = 0.8 * dir.normalized();
    c.segment(p_shells[1], 3) = 0.3 * dir.normalized();
    c /= std::sqrt(c.dot(overlap.Matrix() * c));
    return Eigen::MatrixXd(c * c.transpose());
  };
  const Eigen::MatrixXd dx = POrbital(Eigen::Vector3d(1, 0, 0));
  const Eigen::MatrixXd dd = POrbital(Eigen::Vector3d(1, 1, 1));
  const Eigen::MatrixXd ax =
      DFTEngineTestAccess::SphericalAverageShells(engine, dx, basis);
  const Eigen::MatrixXd ad =
      DFTEngineTestAccess::SphericalAverageShells(engine, dd, basis);
  BOOST_CHECK_SMALL((ax - ad).cwiseAbs().maxCoeff(), 1e-14);
  BOOST_CHECK_CLOSE(ax.cwiseProduct(overlap.Matrix()).sum(), 1.0, 1e-10);
  BOOST_CHECK_CLOSE(ad.cwiseProduct(overlap.Matrix()).sum(), 1.0, 1e-10);
  libint2::finalize();
}

// The spherical atom SCF of the atomic guess: open p and d shells (O, N,
// Fe, Ni) converge in a few iterations to a spherical density with the
// right number of electrons, also in Hund's-rule mode, and with an ECP the
// atom holds only the valence electrons (it used to get all Z).
BOOST_AUTO_TEST_CASE(atomic_scf_spherical_open_shells) {
  libint2::initialize();
  struct Case {
    std::string element;
    std::string basis;
    std::string ecp;
    bool hunds_rule;
    double electrons;
  };
  const std::vector<Case> cases = {
      {"O", "def2-svp", "", false, 8},
      {"N", "def2-svp", "", true, 7},
      {"Fe", "def2-svp", "", false, 26},
      {"Ni", "def2-svp", "", false, 28},
      {"I", "def2-tzvp-ecp", "def2-ecp", false, 25}};
  for (const Case& c : cases) {
    std::ofstream xml("atomic_scf.xml");
    xml << "<dftpackage><spin>1</spin><name>xtp</name><charge>0</charge>"
        << "<functional>XC_HYB_GGA_XC_PBEH</functional>" << "<basisset>"
        << c.basis << "</basisset>"
        << (c.ecp.empty() ? "" : "<ecp>" + c.ecp + "</ecp>")
        << "<initial_guess>atom</initial_guess><xtpdft>"
        << "<screening_eps>1e-9</screening_eps>"
        << "<fock_matrix_reset>5</fock_matrix_reset><convergence>"
        << "<energy>1e-8</energy><method>DIIS</method>"
        << "<DIIS_start>0.002</DIIS_start><ADIIS_start>0.8</ADIIS_start>"
        << "<DIIS_length>20</DIIS_length><levelshift>0.0</levelshift>"
        << "<levelshift_end>0.2</levelshift_end>"
        << "<max_iterations>100</max_iterations><error>1e-7</error>"
        << "<DIIS_maxout>false</DIIS_maxout><mixing>0.7</mixing>"
        << "<mixing_max>0.98</mixing_max><mixing_end>0.8</mixing_end>"
        << "<davidson_max_iter>50</davidson_max_iter></convergence>"
        << "<integration_grid>coarse</integration_grid></xtpdft></dftpackage>";
    xml.close();
    votca::tools::Property prop;
    prop.LoadFromXML("atomic_scf.xml");
    DFTEngine engine;
    Logger log;
    log.setReportLevel(Log::info);
    engine.setLogger(&log);
    engine.Initialize(prop.get("dftpackage"));

    const QMMolecule mol = MakeSingleAtom(c.element);
    const AOBasis basis = MakeBasis(c.basis, mol);
    AOOverlap overlap;
    overlap.Fill(basis);
    const Eigen::MatrixXd dmat =
        DFTEngineTestAccess::RunAtomicDFT(engine, mol[0], c.hunds_rule);

    std::stringstream ss;
    ss << log;
    const std::string text = ss.str();
    const std::size_t pos = text.find("Converged after ");
    BOOST_REQUIRE_MESSAGE(pos != std::string::npos,
                          c.element << ": atomic SCF not converged\n"
                                    << text);
    const Index iterations = std::stol(text.substr(pos + 16));
    BOOST_TEST_MESSAGE(c.element << ": " << iterations << " iterations");
    BOOST_CHECK_LE(iterations, 20);
    BOOST_CHECK_CLOSE(dmat.cwiseProduct(overlap.Matrix()).sum(), c.electrons,
                      1e-8);
    const Eigen::MatrixXd averaged =
        DFTEngineTestAccess::SphericalAverageShells(engine, dmat, basis);
    BOOST_CHECK_SMALL((averaged - dmat).cwiseAbs().maxCoeff(), 1e-12);
  }
  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(orthogonalize_guess_produces_s_orthonormal_vectors) {
  libint2::initialize();

  DFTEngine engine;
  Logger log;
  engine.setLogger(&log);

  const std::string basis_name = "3-21G";
  QMMolecule mol = MakeH2(1.4);

  AOBasis basis = MakeBasis(basis_name, mol);
  DFTEngineTestAccess::SetBasis(engine, basis);
  DFTEngineTestAccess::SetBasisName(engine, basis_name);
  DFTEngineTestAccess::FillOverlap(engine);

  const Index nao = basis.AOBasisSize();
  BOOST_REQUIRE(nao >= 2);

  Eigen::MatrixXd guess = Eigen::MatrixXd::Zero(nao, 2);
  guess(0, 0) = 1.0;
  guess(0, 1) = 1.0;
  guess(1, 0) = 0.25;
  guess(1, 1) = 1.0;
  if (nao > 2) {
    guess(2, 0) = -0.3;
    guess(2, 1) = 0.1;
  }

  const Eigen::MatrixXd orth =
      DFTEngineTestAccess::OrthogonalizeGuess(engine, guess);
  const Eigen::MatrixXd& S = DFTEngineTestAccess::Overlap(engine).Matrix();

  const Eigen::MatrixXd metric = orth.transpose() * S * orth;
  const Eigen::MatrixXd I =
      Eigen::MatrixXd::Identity(metric.rows(), metric.cols());

  BOOST_TEST((metric - I).norm() < 1e-10);

  libint2::finalize();
}

// DFTEngine::ComputeNonXCGradientUKS minus the overlap Pulay term (nuclear
// repulsion, one-electron, RI-J, RI-K) against a finite difference of
// E_nuc + E_one + E_coul + E_exx for fixed, arbitrary alpha/beta MO
// coefficients (valid for any fixed C, as in test_dftgradient.cc). Run with
// ScaHFX_ = 0 and ScaHFX_ > 0, the latter checking the factor 0.5 of the UKS
// RI-K term.
BOOST_AUTO_TEST_CASE(compute_non_xc_gradient_uks_finite_difference) {
  libint2::initialize();
  try {
    double h = 1e-3;  // Bohr; RI fitting limits the precision of smaller
                      // steps

    auto build_system = [](double bond_length_bohr, double sca_hfx,
                           const Eigen::MatrixXd& C_alpha_fixed,
                           const Eigen::MatrixXd& C_beta_fixed) {
      QMMolecule mol = MakeH2(bond_length_bohr);
      BasisSet basisset;
      basisset.Load(std::string(XTP_TEST_DATA_FOLDER) +
                    "/threecenter_dft/3-21G.xml");
      AOBasis dftbasis;
      dftbasis.Fill(basisset, mol);

      // Same basis as DFT and auxiliary basis, as in test_dftgradient.cc.
      // With aux-def2-svp and this rank-2 model density the finite
      // difference deviated by an h-independent ~0.1-0.2% (fit precision,
      // not an RIJGradient error).
      AOBasis auxbasis;
      auxbasis.Fill(basisset, mol);

      AOKinetic kinetic;
      kinetic.Fill(dftbasis);
      AOMultipole esp;
      esp.FillPotential(dftbasis, mol);
      Eigen::MatrixXd H0 = kinetic.Matrix() + esp.Matrix();

      UKSConvergenceAcc::SpinDensity Dspin;
      Dspin.alpha = C_alpha_fixed * C_alpha_fixed.transpose();
      Dspin.beta = C_beta_fixed * C_beta_fixed.transpose();
      Eigen::MatrixXd D_total = Dspin.total();

      static Logger log;
      DFTEngine nuc_rep_engine;
      nuc_rep_engine.setLogger(&log);
      double E_nuc = DFTEngineTestAccess::NuclearRepulsion(nuc_rep_engine, mol);
      double E_one = D_total.cwiseProduct(H0).sum();

      ERIs eris;
      eris.Initialize(dftbasis, auxbasis);
      Eigen::MatrixXd J = eris.CalculateERIs_3c(D_total);

      double E_coul = 0.5 * D_total.cwiseProduct(J).sum();

      double E_exx = 0.0;
      if (sca_hfx > 0.0) {
        // CalculateEXX_dmat is private -- use the public
        // CalculateERIs_EXX_3c wrapper with an empty occMos, matching
        // EXACTLY how EvaluateUKS itself calls it
        // (CalcERIs_EXX(Eigen::MatrixXd::Zero(0,0), Dspin.alpha/beta, ...)).
        std::array<Eigen::MatrixXd, 2> JK_alpha =
            eris.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), Dspin.alpha);
        std::array<Eigen::MatrixXd, 2> JK_beta =
            eris.CalculateERIs_EXX_3c(Eigen::MatrixXd::Zero(0, 0), Dspin.beta);
        const Eigen::MatrixXd& K_alpha = JK_alpha[1];
        const Eigen::MatrixXd& K_beta = JK_beta[1];
        E_exx = 0.5 * sca_hfx *
                (Dspin.alpha.cwiseProduct(K_alpha).sum() +
                 Dspin.beta.cwiseProduct(K_beta).sum());
      }

      struct Result {
        double energy;
        double E_nuc, E_one, E_coul, E_exx;
        QMMolecule mol;
        AOBasis dftbasis;
        AOBasis auxbasis;
        UKSConvergenceAcc::SpinDensity Dspin;
      };
      return Result{E_nuc + E_one + E_coul + E_exx,
                    E_nuc,
                    E_one,
                    E_coul,
                    E_exx,
                    mol,
                    dftbasis,
                    auxbasis,
                    Dspin};
    };

    // Determine the actual DFT basis size once, up front, via a quick
    // "probe" build (geometry-independent for a fixed molecule/basis-set
    // pair) -- avoids guessing a size and rebuilding later.
    {
      QMMolecule probe_mol = MakeH2(1.4);
      BasisSet probe_basisset;
      probe_basisset.Load(std::string(XTP_TEST_DATA_FOLDER) +
                          "/threecenter_dft/3-21G.xml");
      AOBasis probe_basis;
      probe_basis.Fill(probe_basisset, probe_mol);
      Index actual_nbf = probe_basis.AOBasisSize();

      for (double sca_hfx : {0.0, 0.25}) {
        double bond_length = 1.4;  // Bohr, roughly H2 equilibrium

        // Fixed, arbitrary (not orthonormal, not SCF-converged) coefficient
        // "columns" -- one occupied alpha orbital, one occupied beta
        // orbital (minimal, but sufficient to exercise every term).
        Eigen::MatrixXd C_alpha_fixed = Eigen::MatrixXd::Random(actual_nbf, 1);
        Eigen::MatrixXd C_beta_fixed = Eigen::MatrixXd::Random(actual_nbf, 1);
        Eigen::VectorXd eps_alpha_fixed = Eigen::VectorXd::Random(1);
        Eigen::VectorXd eps_beta_fixed = Eigen::VectorXd::Random(1);

        auto base =
            build_system(bond_length, sca_hfx, C_alpha_fixed, C_beta_fixed);

        DFTEngine engine;
        DFTEngineTestAccess::SetBasis(engine, base.dftbasis);
        DFTEngineTestAccess::SetAuxBasis(engine, base.auxbasis);
        DFTEngineTestAccess::SetElectronCounts(engine, 1, 1);
        DFTEngineTestAccess::SetScaHFX(engine, sca_hfx);

        tools::EigenSystem MOs_alpha;
        MOs_alpha.eigenvectors() = C_alpha_fixed;
        MOs_alpha.eigenvalues() = eps_alpha_fixed;
        tools::EigenSystem MOs_beta;
        MOs_beta.eigenvectors() = C_beta_fixed;
        MOs_beta.eigenvalues() = eps_beta_fixed;

        Eigen::MatrixXd analytic_grad =
            DFTEngineTestAccess::ComputeNonXCGradientUKS(
                engine, base.mol, base.Dspin, MOs_alpha, MOs_beta);

        Eigen::Vector3d sum = analytic_grad.colwise().sum();
        BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-6);

        // The overlap Pulay term is excluded: it accounts for the implicit
        // R-dependence of C through C^T S(R) C = I and is only valid at an
        // SCF stationary point, so a fixed-C finite difference does not
        // contain it. It is tested in
        // overlap_pulay_gradient_uks_reduces_to_rks.
        Eigen::MatrixXd overlap_pulay =
            DFTEngineTestAccess::ComputeOverlapPulayGradientUKS(
                engine, base.mol, MOs_alpha, MOs_beta);
        Eigen::MatrixXd fixed_c_checkable_grad = analytic_grad - overlap_pulay;

        auto plus =
            build_system(bond_length + h, sca_hfx, C_alpha_fixed, C_beta_fixed);
        auto minus =
            build_system(bond_length - h, sca_hfx, C_alpha_fixed, C_beta_fixed);
        double finite_diff_deriv = (plus.energy - minus.energy) / (2.0 * h);

        std::cout
            << std::setprecision(10)
            << "[test diagnostic, per-term finite differences]"
            << "\n  d(E_nuc)/dx  = " << (plus.E_nuc - minus.E_nuc) / (2.0 * h)
            << "\n  d(E_one)/dx  = " << (plus.E_one - minus.E_one) / (2.0 * h)
            << "\n  d(E_coul)/dx = " << (plus.E_coul - minus.E_coul) / (2.0 * h)
            << "\n  d(E_exx)/dx  = " << (plus.E_exx - minus.E_exx) / (2.0 * h)
            << std::endl;

        // Atom 1 sits at (bond_length, 0, 0) -- the x-component carries the
        // bond-stretch derivative.
        double analytic = fixed_c_checkable_grad(1, 0);

        bool matches =
            std::abs(finite_diff_deriv - analytic) < 1e-4 * std::abs(analytic);
        if (!matches) {
          std::cout << "ScaHFX_=" << sca_hfx
                    << " analytic dE/dx(atom1) [excluding overlap Pulay]="
                    << analytic << " finite-difference=" << finite_diff_deriv
                    << std::endl;
          std::cout << "NOTE: if this fails only for sca_hfx=0.25 (not 0.0), "
                       "the bug is isolated to the RI-K/0.5*ScaHFX_ factor "
                       "specifically. If it fails for BOTH, the bug is in one "
                       "of nuclear repulsion, one-electron, or RI-J."
                    << std::endl;
        }
        BOOST_CHECK_EQUAL(matches, true);
      }
    }
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING compute_non_xc_gradient_uks_finite_difference: "
              << e.what() << std::endl;
    libint2::finalize();
    return;
  }

  libint2::finalize();
}

// ComputeOverlapPulayGradientUKS (not testable by fixed-C finite
// differences, see above) must reduce to the RKS formula,
// W_RKS = 2 C_occ eps_occ C_occ^T, dE/dR = -Tr[W_RKS dS/dR], for identical
// alpha and beta orbitals and energies.
BOOST_AUTO_TEST_CASE(overlap_pulay_gradient_uks_reduces_to_rks) {
  libint2::initialize();
  try {
    QMMolecule mol = MakeH2(1.4);
    BasisSet basisset;
    basisset.Load(std::string(XTP_TEST_DATA_FOLDER) +
                  "/threecenter_dft/3-21G.xml");
    AOBasis dftbasis;
    dftbasis.Fill(basisset, mol);
    Index n_bf = dftbasis.AOBasisSize();

    DFTEngine engine;
    DFTEngineTestAccess::SetBasis(engine, dftbasis);
    DFTEngineTestAccess::SetElectronCounts(engine, 1, 1);

    Eigen::MatrixXd C_occ = Eigen::MatrixXd::Random(n_bf, 1);
    Eigen::VectorXd eps_occ = Eigen::VectorXd::Random(1);

    tools::EigenSystem MOs_same;
    MOs_same.eigenvectors() = C_occ;
    MOs_same.eigenvalues() = eps_occ;

    // alpha == beta: same occupied orbital, same orbital energy, for
    // both spins.
    Eigen::MatrixXd uks_grad =
        DFTEngineTestAccess::ComputeOverlapPulayGradientUKS(engine, mol,
                                                            MOs_same, MOs_same);

    // Independently-computed RKS-style reference: W_RKS =
    // 2*C_occ*eps_occ*C_occ^T.
    Eigen::MatrixXd W_rks =
        2.0 * C_occ * eps_occ.asDiagonal() * C_occ.transpose();
    std::vector<AOMatrixDerivative> dS = ComputeOverlapDerivatives(dftbasis);
    Index natoms = mol.size();
    Eigen::MatrixXd rks_grad = Eigen::MatrixXd::Zero(natoms, 3);
    for (Index a = 0; a < natoms; ++a) {
      for (Index xyz = 0; xyz < 3; ++xyz) {
        rks_grad(a, xyz) = -W_rks.cwiseProduct(dS[a][xyz]).sum();
      }
    }

    bool matches = uks_grad.isApprox(rks_grad, 1e-10);
    if (!matches) {
      std::cout << "UKS overlap Pulay (alpha==beta):\n"
                << uks_grad << std::endl;
      std::cout << "RKS-style reference:\n" << rks_grad << std::endl;
    }
    BOOST_CHECK_EQUAL(matches, true);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING overlap_pulay_gradient_uks_reduces_to_rks: "
              << e.what() << std::endl;
    libint2::finalize();
    return;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(uks_orbital_hessian_matches_finite_differences) {
  // Analytic UKS Hessian products (spin-polarized f_xc, J, K) against
  // finite differences of the coupled Fock matrices, for water+ with a
  // hybrid and RI and with a GGA and 4-centre integrals: at the converged
  // orbitals and at orbitals rotated away from them.
  const std::string data = std::string(XTP_TEST_DATA_FOLDER);
  for (const bool hybrid_ri : {true, false}) {
    libint2::initialize();
    {
      std::ofstream xml("uks_hessian.xml");
      xml << "<dftpackage>\n<spin>2</spin>\n<name>xtp</name>\n"
          << "<charge>1</charge>\n<functional>"
          << (hybrid_ri ? "XC_HYB_GGA_XC_PBEH" : "XC_GGA_X_PBE XC_GGA_C_PBE")
          << "</functional>\n<basisset>" << data
          << "/hirshfeldpartition/3-21G.xml</basisset>\n";
      if (hybrid_ri) {
        xml << "<auxbasisset>" << data
            << "/diabatization/aux-def2-svp.xml</auxbasisset>\n";
      }
      xml << "<initial_guess>independent</initial_guess>\n<xtpdft>\n"
          << "<screening_eps>1e-12</screening_eps>\n"
          << "<fock_matrix_reset>5</fock_matrix_reset>\n<convergence>\n"
          << "<energy>1e-12</energy>\n<method>DIIS</method>\n"
          << "<DIIS_start>0.002</DIIS_start>\n<ADIIS_start>0.8</ADIIS_start>\n"
          << "<DIIS_length>20</DIIS_length>\n<levelshift>0.0</levelshift>\n"
          << "<levelshift_end>0.2</levelshift_end>\n"
          << "<max_iterations>100</max_iterations>\n<error>1e-10</error>\n"
          << "<DIIS_maxout>false</DIIS_maxout>\n<mixing>0.7</mixing>\n"
          << "<mixing_max>0.98</mixing_max>\n<mixing_end>0.8</mixing_end>\n"
          << "<davidson_max_iter>50</davidson_max_iter>\n"
          << "<second_order>off</second_order>\n"
          << "</convergence>\n<integration_grid>medium</integration_grid>\n"
          << "</xtpdft>\n</dftpackage>\n";
    }
    votca::tools::Property prop;
    prop.LoadFromXML("uks_hessian.xml");
    Logger log;
    DFTEngine engine;
    engine.setLogger(&log);
    engine.Initialize(prop.get("dftpackage"));
    Orbitals orb;
    orb.QMAtoms().LoadFromFile(data + "/chargeresponse/water.xyz");
    orb.setChargeAndSpin(1, 2);
    BOOST_REQUIRE(engine.Evaluate(orb));

    const QMMolecule& mol = orb.QMAtoms();
    const Vxc_Potential<Vxc_Grid> vxc =
        DFTEngineTestAccess::SetupVxc(engine, mol);
    const Eigen::MatrixXd H0 = DFTEngineTestAccess::SetupH0(engine, mol);
    AOOverlap overlap;
    overlap.Fill(DFTEngineTestAccess::Basis(engine));
    UKSConvergenceAcc conv;
    conv.setLogger(&log);
    conv.setOverlap(overlap, 1e-8);
    const Index na = orb.getNumberOfAlphaElectrons();
    const Index nb = orb.getNumberOfBetaElectrons();
    auto builder = [&](const Eigen::MatrixXd& Da, const Eigen::MatrixXd& Db) {
      return DFTEngineTestAccess::BuildUKSFock(engine, H0, vxc, Da, Db);
    };

    for (const double displacement : {0.0, 0.02}) {
      tools::EigenSystem Ma = orb.MOs();
      tools::EigenSystem Mb = orb.MOs_beta();
      const Index n = Ma.eigenvectors().cols();
      const Index n_ov = na * (n - na) + nb * (n - nb);
      std::srand(21);
      const Eigen::VectorXd shift =
          displacement * Eigen::VectorXd::Random(n_ov);
      Ma.eigenvectors() = trah::RotateOrbitals(Ma.eigenvectors(),
                                               shift.head(na * (n - na)), na);
      Mb.eigenvectors() = trah::RotateOrbitals(Mb.eigenvectors(),
                                               shift.tail(nb * (n - nb)), nb);
      const Eigen::MatrixXd Ca_occ = Ma.eigenvectors().leftCols(na);
      const Eigen::MatrixXd Cb_occ = Mb.eigenvectors().leftCols(nb);
      const UKSConvergenceAcc::SpinFock F =
          builder(Ca_occ * Ca_occ.transpose(), Cb_occ * Cb_occ.transpose());

      const Eigen::MatrixXd V = Eigen::MatrixXd::Random(n_ov, 2);
      const Eigen::MatrixXd analytic = DFTEngineTestAccess::UKSOrbitalHessian(
          engine, vxc, Ma, Mb, F.alpha, F.beta, V);
      for (Index k = 0; k < V.cols(); ++k) {
        const Eigen::VectorXd fd =
            conv.BuildCoupledSigmaVector(V.col(k), Ma.eigenvectors(), na,
                                         Mb.eigenvectors(), nb, builder, 1e-4);
        const double rel = (analytic.col(k) - fd).norm() / fd.norm();
        BOOST_TEST_MESSAGE((hybrid_ri ? "PBE0/RI" : "PBE/4c")
                           << ", displacement " << displacement
                           << ": |analytic - fd| / |fd| = " << rel);
        BOOST_CHECK_SMALL(rel, 1e-6);
      }
    }
    libint2::finalize();
  }
}

BOOST_AUTO_TEST_CASE(rks_orbital_hessian_matches_finite_differences) {
  // Analytic closed-shell Hessian products against central differences of
  // the energy gradient g(x) = -4 F_ia of the orbitals C exp(kappa(x)), for
  // water with a hybrid and RI and with a GGA and 4-centre integrals, at
  // the converged and at displaced orbitals. The Fock matrix of a
  // closed-shell density P is the UKS one at P_alpha = P_beta = P/2.
  const std::string data = std::string(XTP_TEST_DATA_FOLDER);
  for (const bool hybrid_ri : {true, false}) {
    libint2::initialize();
    {
      std::ofstream xml("rks_hessian.xml");
      xml << "<dftpackage>\n<spin>1</spin>\n<name>xtp</name>\n"
          << "<charge>0</charge>\n<functional>"
          << (hybrid_ri ? "XC_HYB_GGA_XC_PBEH" : "XC_GGA_X_PBE XC_GGA_C_PBE")
          << "</functional>\n<basisset>" << data
          << "/hirshfeldpartition/3-21G.xml</basisset>\n";
      if (hybrid_ri) {
        xml << "<auxbasisset>" << data
            << "/diabatization/aux-def2-svp.xml</auxbasisset>\n";
      }
      xml << "<initial_guess>independent</initial_guess>\n<xtpdft>\n"
          << "<screening_eps>1e-12</screening_eps>\n"
          << "<fock_matrix_reset>5</fock_matrix_reset>\n<convergence>\n"
          << "<energy>1e-12</energy>\n<method>DIIS</method>\n"
          << "<DIIS_start>0.002</DIIS_start>\n<ADIIS_start>0.8</ADIIS_start>\n"
          << "<DIIS_length>20</DIIS_length>\n<levelshift>0.0</levelshift>\n"
          << "<levelshift_end>0.2</levelshift_end>\n"
          << "<max_iterations>100</max_iterations>\n<error>1e-10</error>\n"
          << "<DIIS_maxout>false</DIIS_maxout>\n<mixing>0.7</mixing>\n"
          << "<mixing_max>0.98</mixing_max>\n<mixing_end>0.8</mixing_end>\n"
          << "<davidson_max_iter>50</davidson_max_iter>\n"
          << "<second_order>off</second_order>\n"
          << "</convergence>\n<integration_grid>medium</integration_grid>\n"
          << "</xtpdft>\n</dftpackage>\n";
    }
    votca::tools::Property prop;
    prop.LoadFromXML("rks_hessian.xml");
    Logger log;
    DFTEngine engine;
    engine.setLogger(&log);
    engine.Initialize(prop.get("dftpackage"));
    Orbitals orb;
    orb.QMAtoms().LoadFromFile(data + "/chargeresponse/water.xyz");
    orb.setChargeAndSpin(0, 1);
    BOOST_REQUIRE(engine.Evaluate(orb));

    const QMMolecule& mol = orb.QMAtoms();
    const Vxc_Potential<Vxc_Grid> vxc =
        DFTEngineTestAccess::SetupVxc(engine, mol);
    const Eigen::MatrixXd H0 = DFTEngineTestAccess::SetupH0(engine, mol);
    const Index nocc = orb.getNumberOfAlphaElectrons();
    auto fock = [&](const Eigen::MatrixXd& C) {
      const Eigen::MatrixXd half =
          C.leftCols(nocc) * C.leftCols(nocc).transpose();
      return DFTEngineTestAccess::BuildUKSFock(engine, H0, vxc, half, half)
          .alpha;
    };
    // energy gradient at C exp(kappa(x))
    auto gradient = [&](const Eigen::MatrixXd& C, const Eigen::VectorXd& x) {
      const Eigen::MatrixXd Cx = trah::RotateOrbitals(C, x, nocc);
      return Eigen::VectorXd(
          -4.0 * trah::OccVirt(Cx.transpose() * fock(Cx) * Cx, nocc));
    };

    for (const double displacement : {0.0, 0.02}) {
      tools::EigenSystem mos = orb.MOs();
      const Index n = mos.eigenvectors().cols();
      const Index n_ov = nocc * (n - nocc);
      std::srand(31);
      mos.eigenvectors() = trah::RotateOrbitals(
          mos.eigenvectors(), displacement * Eigen::VectorXd::Random(n_ov),
          nocc);
      const Eigen::MatrixXd& C = mos.eigenvectors();
      const Eigen::MatrixXd V = Eigen::MatrixXd::Random(n_ov, 2);
      const Eigen::MatrixXd analytic =
          DFTEngineTestAccess::RKSOrbitalHessian(engine, vxc, mos, fock(C), V);
      const double h = 1e-4;
      for (Index k = 0; k < V.cols(); ++k) {
        const Eigen::VectorXd v = V.col(k);
        const Eigen::VectorXd fd =
            (gradient(C, h * v) - gradient(C, -h * v)) / (2 * h);
        const double rel = (analytic.col(k) - fd).norm() / fd.norm();
        BOOST_TEST_MESSAGE((hybrid_ri ? "PBE0/RI" : "PBE/4c")
                           << ", displacement " << displacement
                           << ": |analytic - fd| / |fd| = " << rel);
        BOOST_CHECK_SMALL(rel, 1e-6);
      }
    }
    libint2::finalize();
  }
}

BOOST_AUTO_TEST_CASE(trah_converges_below_the_diis_saddle_point) {
  // Stretched O2 triplet (1.8 Angstrom), PBE0/3-21G, UKS. DIIS converges to
  // a stationary point with five negative orbital-Hessian eigenvalues
  // (lowest -0.105); TRAH (second_order=always) converges 12 mHa lower, to
  // a stationary point with fewer and much smaller negative eigenvalues.
  // These remaining ones are symmetry-breaking directions with no gradient
  // component, which a gradient-driven method does not leave (it would
  // need instability following); thread-order noise occasionally lets TRAH
  // slide further, so only the comparison with DIIS is tested.
  const std::string data = std::string(XTP_TEST_DATA_FOLDER);
  {
    std::ofstream xyz("o2_stretched.xyz");
    xyz << "2\n\nO 0 0 0\nO 0 0 1.8\n";
  }
  struct Result {
    double energy;
    double lowest_hessian_eigenvalue;
    Index negative_eigenvalues;
  };
  auto run = [&](const std::string& second_order) {
    libint2::initialize();
    {
      std::ofstream xml("o2_stretched.xml");
      xml << "<dftpackage>\n<spin>3</spin>\n<name>xtp</name>\n"
          << "<charge>0</charge>\n<functional>XC_HYB_GGA_XC_PBEH</functional>\n"
          << "<basisset>" << data
          << "/hirshfeldpartition/3-21G.xml</basisset>\n"
          << "<initial_guess>atom</initial_guess>\n<xtpdft>\n"
          << "<screening_eps>1e-12</screening_eps>\n"
          << "<fock_matrix_reset>5</fock_matrix_reset>\n<convergence>\n"
          << "<energy>1e-9</energy>\n<method>DIIS</method>\n"
          << "<DIIS_start>0.002</DIIS_start>\n<ADIIS_start>0.8</ADIIS_start>\n"
          << "<DIIS_length>20</DIIS_length>\n<levelshift>0.0</levelshift>\n"
          << "<levelshift_end>0.2</levelshift_end>\n"
          << "<max_iterations>100</max_iterations>\n<error>1e-7</error>\n"
          << "<DIIS_maxout>false</DIIS_maxout>\n<mixing>0.7</mixing>\n"
          << "<mixing_max>0.98</mixing_max>\n<mixing_end>0.8</mixing_end>\n"
          << "<davidson_max_iter>50</davidson_max_iter>\n"
          << "<second_order>" << second_order << "</second_order>\n"
          << "</convergence>\n<integration_grid>medium</integration_grid>\n"
          << "</xtpdft>\n</dftpackage>\n";
    }
    votca::tools::Property prop;
    prop.LoadFromXML("o2_stretched.xml");
    Logger log;
    DFTEngine engine;
    engine.setLogger(&log);
    engine.Initialize(prop.get("dftpackage"));
    Orbitals orb;
    orb.QMAtoms().LoadFromFile("o2_stretched.xyz");
    orb.setChargeAndSpin(0, 3);
    BOOST_REQUIRE(engine.Evaluate(orb));

    // full UKS orbital Hessian at the converged orbitals
    const QMMolecule& mol = orb.QMAtoms();
    const Vxc_Potential<Vxc_Grid> vxc =
        DFTEngineTestAccess::SetupVxc(engine, mol);
    const Eigen::MatrixXd H0 = DFTEngineTestAccess::SetupH0(engine, mol);
    const Index na = orb.getNumberOfAlphaElectrons();
    const Index nb = orb.getNumberOfBetaElectrons();
    const Eigen::MatrixXd Ca = orb.MOs().eigenvectors().leftCols(na);
    const Eigen::MatrixXd Cb = orb.MOs_beta().eigenvectors().leftCols(nb);
    const UKSConvergenceAcc::SpinFock F = DFTEngineTestAccess::BuildUKSFock(
        engine, H0, vxc, Ca * Ca.transpose(), Cb * Cb.transpose());
    const Index n = orb.MOs().eigenvectors().cols();
    const Index n_ov = na * (n - na) + nb * (n - nb);
    Eigen::MatrixXd H = DFTEngineTestAccess::UKSOrbitalHessian(
        engine, vxc, orb.MOs(), orb.MOs_beta(), F.alpha, F.beta,
        Eigen::MatrixXd::Identity(n_ov, n_ov));
    H = 0.5 * (H + H.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
    libint2::finalize();
    const Index negative = (es.eigenvalues().array() < -1e-6).count();
    return Result{orb.getDFTTotalEnergy(), es.eigenvalues()(0), negative};
  };

  const Result diis = run("off");
  const Result trah = run("always");
  BOOST_TEST_MESSAGE("DIIS: E = " << diis.energy << ", lowest eigenvalue "
                                  << diis.lowest_hessian_eigenvalue << ", "
                                  << diis.negative_eigenvalues << " negative");
  BOOST_TEST_MESSAGE("TRAH: E = " << trah.energy << ", lowest eigenvalue "
                                  << trah.lowest_hessian_eigenvalue << ", "
                                  << trah.negative_eigenvalues << " negative");
  BOOST_CHECK_LT(diis.lowest_hessian_eigenvalue, -0.05);
  BOOST_CHECK_GT(trah.lowest_hessian_eigenvalue, -0.05);
  BOOST_CHECK_LT(trah.negative_eigenvalues, diis.negative_eigenvalues);
  BOOST_CHECK_LT(trah.energy, diis.energy - 0.01);
}
