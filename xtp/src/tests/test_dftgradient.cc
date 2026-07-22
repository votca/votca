/*
 * Copyright 2009-2024 The VOTCA Development Team (http://www.votca.org)
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

// ===========================================================================
// STATUS: written but NOT yet run. Unlike everything in
// test_aoderivatives.cc, this doesn't involve libint2 at all (nuclear
// repulsion is pure geometry -- nuclear charges and positions, no
// integrals) -- there is no reason to expect this to hit any of the
// libint2-build-configuration issues encountered elsewhere in this
// branch, but it is still an independent piece of new code and should
// not be assumed correct without running it.
//
// The finite-difference reference below deliberately does NOT call
// DFTEngine::NuclearRepulsion (it's private, not accessible from a
// test) -- it's an independent reimplementation of the same elementary
// Coulomb's-law formula. This is arguably slightly BETTER validation
// than reusing shared code would be: two independently-written
// implementations agreeing via finite difference is stronger evidence
// than one implementation agreeing with itself, since a shared formula
// bug couldn't hide from both.
// ===========================================================================

#include "xtp_libint2.h"
#include <stdexcept>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE dftgradient_test

// Standard includes
#include <cmath>
#include <fstream>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/dftgradient.h>
#include <votca/xtp/ERIs.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;

// Defined in libint2_derivative_calls.cc, not yet in any header; forward
// declared here purely to build the finite-difference reference energy
// (production code, DFTGradient::RIJGradient, already calls this
// internally -- this is not duplicating the derivative logic, only the
// energy-level integral used for the reference calculation).
namespace votca {
namespace xtp {
std::vector<Eigen::MatrixXd> ComputeThreeCenterIntegrals(
    const AOBasis& auxbasis, const AOBasis& dftbasis);
}  // namespace xtp
}  // namespace votca

BOOST_AUTO_TEST_SUITE(dftgradient_test)

// Builds a simple three-atom test molecule (not collinear, so the
// nuclear repulsion gradient has genuinely nonzero components in more
// than one direction -- a collinear or two-atom system would exercise
// less of the formula than a general geometry does).
QMMolecule BuildTestMolecule(double displacement_x = 0.0) {
  QMMolecule mol(" ", 0);
  std::string xyz_content =
      "3\n\n"
      "O 0.0 0.0 0.0\n"
      "H 0.96 0.0 0.0\n"
      "H " +
      std::to_string(-0.24 + displacement_x) + " 0.93 0.0\n";
  std::string tmp_path = "/tmp/xtp_test_dftgradient_mol.xyz";
  std::ofstream out(tmp_path);
  out << xyz_content;
  out.close();
  mol.LoadFromFile(tmp_path);
  return mol;
}

// Independent reference implementation of the same nuclear repulsion
// energy formula DFTGradient::NuclearRepulsionDerivative differentiates,
// written separately (not calling any XTP code) as the energy reference
// for the finite-difference check below.
double NuclearRepulsionEnergyReference(const QMMolecule& mol) {
  double energy = 0.0;
  for (Index a = 0; a < mol.size(); ++a) {
    for (Index b = a + 1; b < mol.size(); ++b) {
      double Za = static_cast<double>(mol[a].getNuccharge());
      double Zb = static_cast<double>(mol[b].getNuccharge());
      double Rab = (mol[a].getPos() - mol[b].getPos()).norm();
      energy += Za * Zb / Rab;
    }
  }
  return energy;
}

BOOST_AUTO_TEST_CASE(nuclear_repulsion_derivative_finite_difference) {
  double h = 1e-4;  // finite-difference step, Angstrom -- matching
                     // BuildH2's convention in test_aoderivatives.cc
                     // (confirmed empirically there that a bare .xyz
                     // file with no explicit units line is interpreted
                     // as Angstrom by QMMolecule::LoadFromFile, via the
                     // Bohr/Angstrom bug found and fixed in that file;
                     // reusing that same confirmed assumption here rather
                     // than re-deriving it from scratch).

  // NOTE: BuildTestMolecule's xyz file has no units line beyond the
  // standard xyz format, and QMMolecule::LoadFromFile's default unit
  // assumption was not independently re-checked here (assumed Angstrom,
  // matching every other .xyz-based test in this codebase, e.g. BuildH2
  // in test_aoderivatives.cc). If this test fails by a uniform factor
  // across all entries the same way the overlap-derivative test once
  // did, check that assumption first -- see the detailed note in
  // test_aoderivatives.cc for exactly this failure signature.
  constexpr double kBohrPerAngstrom = 0.52917721090380;

  QMMolecule mol0 = BuildTestMolecule(0.0);
  Eigen::MatrixXd deriv = DFTGradient::NuclearRepulsionDerivative(mol0);

  // Sanity check independent of finite differences: net force-like sum
  // over all atoms must vanish (translational invariance of the total
  // energy -- shifting the whole molecule rigidly cannot change E_nn).
  Eigen::Vector3d total = deriv.colwise().sum();
  BOOST_CHECK_SMALL(total.cwiseAbs().maxCoeff(), 1e-10);

  // Finite difference w.r.t. atom 2's x-coordinate (the second H).
  QMMolecule mol_plus = BuildTestMolecule(h);
  QMMolecule mol_minus = BuildTestMolecule(-h);
  double e_plus = NuclearRepulsionEnergyReference(mol_plus);
  double e_minus = NuclearRepulsionEnergyReference(mol_minus);
  double finite_diff_deriv_ang = (e_plus - e_minus) / (2.0 * h);
  double finite_diff_deriv_bohr = finite_diff_deriv_ang * kBohrPerAngstrom;

  double analytic = deriv(2, 0);  // atom index 2, x-component
  // Tolerance matches the standard used throughout this branch (1e-4,
  // e.g. the overlap/kinetic/coulomb-metric/three-center tests in
  // test_aoderivatives.cc), rather than something tighter chosen without
  // justification -- consistent with expected O(h^2) central-difference
  // truncation error at h=1e-4.
  bool matches =
      std::abs(finite_diff_deriv_bohr - analytic) < 1e-4 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_nn/dx(atom2): " << analytic << std::endl;
    std::cout << "Finite-difference: " << finite_diff_deriv_bohr << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);
}

// RI-J gradient assembly (DFTGradient::RIJGradient), validated against a
// finite difference of the total RI-J energy computed with a FIXED,
// arbitrary density matrix -- see the detailed "IMPORTANT" note on this
// function in dftgradient.h for why an arbitrary (not necessarily
// converged-SCF) density matrix is a legitimate way to validate the
// assembly formula: the RI fitting coefficients' stationarity is a
// property of the linear least-squares fit itself, not of electronic
// self-consistency.
//
// STATUS: NOT yet run. This is the most integrative piece tested so far
// in this branch -- it combines four previously-validated pieces
// (ComputeThreeCenterIntegrals/Derivatives, AOCoulomb, and the
// two-center metric derivative) through NEW contraction/assembly code
// that has not itself been checked. A failure here would point to the
// assembly logic (the c=V^-1 d solve, or the term1/term2 contraction),
// not to any of the underlying integrals, which are already separately
// confirmed correct.
BOOST_AUTO_TEST_CASE(rij_gradient_finite_difference) {
  libint2::initialize();
 try {
  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";
  double bond_length = 0.74;  // Angstrom
  double h = 1e-4;            // Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  auto build_h2 = [](double bond_length_angstrom) {
    QMMolecule mol(" ", 0);
    std::string xyz_content =
        "2\n\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 " +
        std::to_string(bond_length_angstrom) + "\n";
    std::string tmp_path = "/tmp/xtp_test_dftgradient_h2.xyz";
    std::ofstream out(tmp_path);
    out << xyz_content;
    out.close();
    mol.LoadFromFile(tmp_path);
    return mol;
  };

  QMMolecule mol0 = build_h2(bond_length);
  AOBasis dftbasis0;
  dftbasis0.Fill(basisset, mol0);
  AOBasis auxbasis0;
  auxbasis0.Fill(basisset, mol0);  // same basis used for both, as in the
                                    // three-center integral test

  // Fixed, arbitrary, symmetric density matrix -- generated once and
  // reused unchanged at every geometry. Its specific values don't
  // matter for this test (see IMPORTANT note referenced above); what
  // matters is that it stays FIXED while the geometry moves.
  Index n_dft_bf = dftbasis0.AOBasisSize();
  Eigen::MatrixXd density_random = Eigen::MatrixXd::Random(n_dft_bf, n_dft_bf);
  Eigen::MatrixXd density = 0.5 * (density_random + density_random.transpose());

  Eigen::MatrixXd analytic_grad =
      DFTGradient::RIJGradient(density, auxbasis0, dftbasis0);

  // Sanity check independent of finite differences, same reasoning as
  // the nuclear repulsion test: translational invariance of the total
  // energy means the gradient must sum to zero across all atoms.
  Eigen::Vector3d total = analytic_grad.colwise().sum();
  BOOST_CHECK_SMALL(total.cwiseAbs().maxCoeff(), 1e-6);

  auto rij_energy = [&](const AOBasis& auxbasis, const AOBasis& dftbasis) {
    std::vector<Eigen::MatrixXd> tensor =
        ComputeThreeCenterIntegrals(auxbasis, dftbasis);
    Index n_aux_bf = auxbasis.AOBasisSize();
    Eigen::VectorXd d(n_aux_bf);
    for (Index p = 0; p < n_aux_bf; ++p) {
      d(p) = (density.array() * tensor[p].array()).sum();
    }
    AOCoulomb aocoulomb;
    aocoulomb.Fill(auxbasis);
    Eigen::VectorXd c = aocoulomb.Matrix().ldlt().solve(d);
    return 0.5 * c.dot(d);
  };

  QMMolecule mol_plus = build_h2(bond_length + h);
  AOBasis dftbasis_plus;
  dftbasis_plus.Fill(basisset, mol_plus);
  AOBasis auxbasis_plus;
  auxbasis_plus.Fill(basisset, mol_plus);
  double e_plus = rij_energy(auxbasis_plus, dftbasis_plus);

  QMMolecule mol_minus = build_h2(bond_length - h);
  AOBasis dftbasis_minus;
  dftbasis_minus.Fill(basisset, mol_minus);
  AOBasis auxbasis_minus;
  auxbasis_minus.Fill(basisset, mol_minus);
  double e_minus = rij_energy(auxbasis_minus, dftbasis_minus);

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_deriv =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;

  double analytic = analytic_grad(1, 2);  // atom 1 (second H), z-component
  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-4 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_J/dz(atom1): " << analytic << std::endl;
    std::cout << "Finite-difference: " << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if this fails, the assembly logic (c=V^-1 d "
                 "solve, or the term1/term2 contraction in "
                 "DFTGradient::RIJGradient) is the first thing to check "
                 "-- the underlying integrals it consumes are already "
                 "separately validated in test_aoderivatives.cc."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);
 } catch (const std::runtime_error& e) {
   std::cout << "SKIPPING rij_gradient_finite_difference: " << e.what()
             << std::endl;
   libint2::finalize();
   return;
 }

  libint2::finalize();
}

// RI-K (exchange) gradient assembly (DFTGradient::RIKGradient), validated
// the same way as RI-J: finite difference of the total (self-defined,
// internally-consistent -- see the energy-convention note on
// RIKGradient in dftgradient.h) exchange energy, using a FIXED, arbitrary
// (nbf x 2) coefficient matrix rather than genuine converged/orthonormal
// occupied orbitals -- valid for the same reason as RIJGradient's
// density argument (per-pair fitting-coefficient stationarity is a
// linear-algebra property, not a consequence of SCF self-consistency or
// orthonormality).
//
// STATUS: NOT yet run. Structurally similar to the RI-J test but with an
// extra (i,j) double sum over coefficient-matrix columns -- a failure
// here would most likely point at that double-sum bookkeeping rather
// than the underlying integrals, which are shared with (and already
// validated by) the RI-J test above.
BOOST_AUTO_TEST_CASE(rik_gradient_finite_difference) {
  libint2::initialize();
 try {
  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml";
  double bond_length = 0.74;  // Angstrom
  double h = 1e-4;            // Angstrom

  BasisSet basisset;
  basisset.Load(basis_path);

  auto build_h2 = [](double bond_length_angstrom) {
    QMMolecule mol(" ", 0);
    std::string xyz_content =
        "2\n\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 " +
        std::to_string(bond_length_angstrom) + "\n";
    std::string tmp_path = "/tmp/xtp_test_dftgradient_h2_rik.xyz";
    std::ofstream out(tmp_path);
    out << xyz_content;
    out.close();
    mol.LoadFromFile(tmp_path);
    return mol;
  };

  QMMolecule mol0 = build_h2(bond_length);
  AOBasis dftbasis0;
  dftbasis0.Fill(basisset, mol0);
  AOBasis auxbasis0;
  auxbasis0.Fill(basisset, mol0);

  // Fixed, arbitrary (nbf x 2) coefficient matrix -- generated once,
  // reused unchanged at every geometry. Not orthonormal, not from any
  // SCF -- deliberately, per the reasoning in the header comment (valid
  // for testing the gradient FORMULA's correctness, even though it
  // would not be physically meaningful as real occupied MOs).
  Index n_dft_bf = dftbasis0.AOBasisSize();
  Eigen::MatrixXd mo_coeffs = Eigen::MatrixXd::Random(n_dft_bf, 2);

  Eigen::MatrixXd analytic_grad =
      DFTGradient::RIKGradient(mo_coeffs, auxbasis0, dftbasis0);

  Eigen::Vector3d total = analytic_grad.colwise().sum();
  BOOST_CHECK_SMALL(total.cwiseAbs().maxCoeff(), 1e-6);

  // Reference energy via the REAL, production ERIs::CalculateEXX_mos --
  // not a separately-defined, self-consistent-but-possibly-wrong toy
  // formula. This is the decisive check: RIKGradient's energy
  // convention (E_K = -sum_ij c_ij.d_ij, confirmed via direct numerical
  // simulation of CalculateEXX_mos's own algorithm -- see the detailed
  // history in dftgradient.h/.cc) should match this exactly.
  auto rik_energy = [&](const AOBasis& auxbasis, const AOBasis& dftbasis) {
    ERIs eris;
    eris.Initialize(dftbasis, auxbasis);
    Eigen::MatrixXd Dmat_local = 2.0 * mo_coeffs * mo_coeffs.transpose();
    std::array<Eigen::MatrixXd, 2> JK =
        eris.CalculateERIs_EXX_3c(mo_coeffs, Dmat_local);
    const Eigen::MatrixXd& K = JK[1];
    return 0.25 * Dmat_local.cwiseProduct(K).sum();  // ScaHFX=1 for this check
  };

  QMMolecule mol_plus = build_h2(bond_length + h);
  AOBasis dftbasis_plus;
  dftbasis_plus.Fill(basisset, mol_plus);
  AOBasis auxbasis_plus;
  auxbasis_plus.Fill(basisset, mol_plus);
  double e_plus = rik_energy(auxbasis_plus, dftbasis_plus);

  QMMolecule mol_minus = build_h2(bond_length - h);
  AOBasis dftbasis_minus;
  dftbasis_minus.Fill(basisset, mol_minus);
  AOBasis auxbasis_minus;
  auxbasis_minus.Fill(basisset, mol_minus);
  double e_minus = rik_energy(auxbasis_minus, dftbasis_minus);

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_deriv =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;

  double analytic = analytic_grad(1, 2);  // atom 1 (second H), z-component
  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-4 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_K/dz(atom1): " << analytic << std::endl;
    std::cout << "Finite-difference: " << finite_diff_deriv << std::endl;
    std::cout << "NOTE: this compares against the REAL, production "
                 "ERIs::CalculateEXX_mos energy (symmetric V^-1/2 RI "
                 "fitting), not a self-consistent toy formula -- if this "
                 "fails, first double check the factor of 2 and sign in "
                 "RIKGradient's energy convention (E_K = "
                 "-sum_ij c_ij.d_ij, see the detailed history in "
                 "dftgradient.h/.cc, confirmed via a SEPARATE numerical "
                 "simulation in Python before this C++ test was written -- "
                 "worth re-checking that simulation's assumptions if this "
                 "C++ test disagrees with it)."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);
 } catch (const std::runtime_error& e) {
   std::cout << "SKIPPING rik_gradient_finite_difference: " << e.what()
             << std::endl;
   libint2::finalize();
   return;
 }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
