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
// STATUS: written but NOT yet run. This is the combined validation of the
// two XC gradient pieces (Vxc_Potential::PulayGradient and
// ::GridWeightGradient) -- neither piece can be checked against the real
// total XC energy in isolation (see the commit history for both), so
// this test sums them and compares against a genuine finite difference
// of IntegrateVXC's real energy output.
//
// Reuses the EXISTING, already-validated test data from
// test_vxc_potential.cc (methane/3-21G/dmat.mm) rather than a toy
// molecule, specifically because it's a real multi-atom (5-atom) case --
// a good test of the SSW weight derivative's multi-atom bookkeeping,
// which a 2-atom system exercises much less of.
//
// Uses an LDA functional (XC_LDA_X + XC_LDA_C_VWN, standard, should
// exist in any LibXC install) rather than the existing test's GGA
// (XC_GGA_X_PBE + XC_GGA_C_PBE), since PulayGradient/GridWeightGradient
// are explicitly scoped to LDA only -- the additional GGA df_dsigma-driven
// Pulay term is not implemented (see PulayGradient's header comment).
// Using the GGA functional here would compare against a formula that's
// deliberately incomplete for that case and should be expected to fail
// for reasons unrelated to whether the LDA-scoped code itself is correct.
//
// Uses a coarser grid quality ("coarse") than the existing test's
// "medium" (53404 points) to keep runtime reasonable, given
// GridWeightGradient's O(Natoms^3)-per-point cost flagged in its own
// header comment.
// ===========================================================================

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE xcgradient_combined_test

// Standard includes
#include <cmath>
#include <fstream>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>
#include <votca/tools/eigenio_matrixmarket.h>

// Local VOTCA includes
#include "votca/xtp/orbitals.h"
#include "votca/xtp/vxc_grid.h"
#include "votca/xtp/vxc_potential.h"
#include "xtp_libint2.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(xcgradient_combined_test)

// Base methane geometry, read once from the existing test data, so the
// displaced-geometry builder below only needs to shift one coordinate
// rather than hardcode the whole molecule a second time.
std::vector<std::pair<std::string, Eigen::Vector3d>> LoadBaseGeometry() {
  std::ifstream in(std::string(XTP_TEST_DATA_FOLDER) +
                   "/vxc_potential/molecule.xyz");
  std::string line;
  std::getline(in, line);  // atom count
  std::getline(in, line);  // comment line
  std::vector<std::pair<std::string, Eigen::Vector3d>> atoms;
  std::string element;
  double x, y, z;
  while (in >> element >> x >> y >> z) {
    atoms.emplace_back(element, Eigen::Vector3d(x, y, z));
  }
  return atoms;
}

// Builds a QMMolecule with atom `atom_idx`'s z-coordinate shifted by
// `dz` (Angstrom, matching the plain .xyz convention already confirmed
// elsewhere in this branch), all other atoms unchanged.
QMMolecule BuildDisplacedMethane(
    const std::vector<std::pair<std::string, Eigen::Vector3d>>& base_atoms,
    Index atom_idx, double dz) {
  std::string xyz_content =
      std::to_string(base_atoms.size()) + "\n methane displaced\n";
  for (size_t i = 0; i < base_atoms.size(); ++i) {
    Eigen::Vector3d pos = base_atoms[i].second;
    if (static_cast<Index>(i) == atom_idx) {
      pos.z() += dz;
    }
    xyz_content += base_atoms[i].first + " " + std::to_string(pos.x()) + " " +
                   std::to_string(pos.y()) + " " + std::to_string(pos.z()) +
                   "\n";
  }
  std::string tmp_path = "/tmp/xtp_test_xcgradient_methane.xyz";
  std::ofstream out(tmp_path);
  out << xyz_content;
  out.close();

  QMMolecule mol("none", 0);
  mol.LoadFromFile(tmp_path);
  return mol;
}

AOBasis BuildBasis(const QMMolecule& mol) {
  BasisSet basis;
  basis.Load(std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/3-21G.xml");
  AOBasis aobasis;
  aobasis.Fill(basis, mol);
  return aobasis;
}

BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference) {
  libint2::initialize();

  const std::string functional = "XC_LDA_X XC_LDA_C_VWN";
  double h = 1e-4;  // Angstrom
  Index displaced_atom = 1;  // first H (index 1, methane is C,H,H,H,H)

  auto base_atoms = LoadBaseGeometry();

  QMMolecule mol0 = BuildDisplacedMethane(base_atoms, displaced_atom, 0.0);
  AOBasis aobasis0 = BuildBasis(mol0);

  // DIAGNOSTIC SWAP: using an arbitrary random (properly symmetrized, no
  // aliasing) density matrix here instead of the real dmat.mm, keeping
  // everything else (real 5-atom geometry, real grid) identical to the
  // original test -- to test directly whether dmat.mm itself is
  // responsible for the -34.86 blowup, given that H2's test (which
  // behaves normally) also uses an arbitrary random matrix rather than
  // a loaded one. dw (the weight derivative) never touches the density
  // matrix at all -- only prefactor=rho*f_xc does -- so if dmat.mm has
  // some mismatch with a freshly-built AOBasis for this geometry (e.g.
  // a basis-ordering incompatibility), that could produce badly wrong
  // rho values without necessarily breaking translational invariance,
  // if the resulting error pattern happens to be symmetric enough.
  //
  // If this test now passes (small, sane weight term, matching H2's
  // pattern), that confirms dmat.mm specifically. If it still blows up
  // even with an arbitrary matrix, that rules out dmat.mm and points at
  // something about the 5-atom geometry/grid itself.
  Index n_bf_methane = aobasis0.AOBasisSize();
  Eigen::MatrixXd dmat_random_methane =
      Eigen::MatrixXd::Random(n_bf_methane, n_bf_methane);
  Eigen::MatrixXd dmat =
      0.5 * (dmat_random_methane + dmat_random_methane.transpose());

  Vxc_Grid grid0;
  grid0.GridSetup("coarse", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  Eigen::MatrixXd pulay_grad = vxc0.PulayGradient(dmat, aobasis0);
  Eigen::MatrixXd weight_grad = vxc0.GridWeightGradient(dmat, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  // Check translational invariance on EACH term separately, not just
  // the total -- H2's weight term alone satisfied this exactly; if
  // methane's weight term alone does NOT, that's definitive evidence of
  // a real bug specific to this case, rather than something in the
  // shared formula (which H2 didn't expose).
  Eigen::Vector3d pulay_sum_check = pulay_grad.colwise().sum();
  Eigen::Vector3d weight_sum_check = weight_grad.colwise().sum();
  std::cerr << "[xc_gradient_finite_difference diagnostic] "
             << "Pulay alone t.i. sum = " << pulay_sum_check.transpose()
             << "\nWeight alone t.i. sum = " << weight_sum_check.transpose()
             << std::endl;

  // Sanity check independent of finite differences, same reasoning as
  // every other gradient in this branch: translational invariance means
  // the total must sum to zero across all atoms.
  Eigen::Vector3d sum = total_grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    std::cout << "WARNING: translational invariance check failed, sum="
              << sum.transpose() << std::endl;
    std::cout << "Pulay term alone:\n" << pulay_grad << std::endl;
    std::cout << "Weight term alone:\n" << weight_grad << std::endl;
  }
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  // Finite-difference reference: real total XC energy (IntegrateVXC),
  // fixed density matrix, displaced geometry rebuilt from scratch each
  // time (fresh grid, fresh basis) -- exactly mirroring how a genuine
  // geometry change would be evaluated, not an isolated/artificial
  // construction.
  QMMolecule mol_plus = BuildDisplacedMethane(base_atoms, displaced_atom, h);
  AOBasis aobasis_plus = BuildBasis(mol_plus);
  Vxc_Grid grid_plus;
  grid_plus.GridSetup("coarse", mol_plus, aobasis_plus);
  Vxc_Potential<Vxc_Grid> vxc_plus(grid_plus);
  vxc_plus.setXCfunctional(functional);
  double e_plus = vxc_plus.IntegrateVXC(dmat).energy();

  QMMolecule mol_minus = BuildDisplacedMethane(base_atoms, displaced_atom, -h);
  AOBasis aobasis_minus = BuildBasis(mol_minus);
  Vxc_Grid grid_minus;
  grid_minus.GridSetup("coarse", mol_minus, aobasis_minus);
  Vxc_Potential<Vxc_Grid> vxc_minus(grid_minus);
  vxc_minus.setXCfunctional(functional);
  double e_minus = vxc_minus.IntegrateVXC(dmat).energy();

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_deriv =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;

  double analytic = total_grad(displaced_atom, 2);  // z-component
  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_xc/dz(atom" << displaced_atom
              << "): " << analytic << std::endl;
    std::cout << "  Pulay contribution: " << pulay_grad(displaced_atom, 2)
              << std::endl;
    std::cout << "  Weight contribution: " << weight_grad(displaced_atom, 2)
              << std::endl;
    std::cout << "Finite-difference: " << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if analytic and finite-difference disagree "
                 "substantially, try isolating which term is wrong by "
                 "checking the translational-invariance sum above per "
                 "term (not just the total), and by reducing to a "
                 "2-atom toy system where hand-checking is more "
                 "tractable. If the discrepancy is small but exceeds "
                 "tolerance, consider that numerical quadrature grids "
                 "can introduce their own finite-difference noise (the "
                 "grid's discrete structure/pruning doesn't vary "
                 "perfectly smoothly with geometry, unlike an analytic "
                 "functional) -- try a finer grid quality (e.g. "
                 "\"medium\") or a larger h before concluding the "
                 "formula itself is wrong."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// Much smaller reproduction case: H2 instead of methane. Methane's
// per-point formula has now been cross-checked against an independent
// Python implementation for THREE diverse points (different owners,
// including the previously-untested owner==A case) with exact
// agreement every time -- yet the aggregate result is still off by
// ~400x. Continuing to sample individual points at the 5-atom,
// ~28000-point scale isn't converging; H2 has far fewer points and
// only 2 atoms, small enough to reason about (or dump) exhaustively
// rather than by sample, the same way reducing to a minimal case
// resolved the libint2 buffer-ordering question earlier in this branch.
BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference_h2) {
  libint2::initialize();

  const std::string functional = "XC_LDA_X XC_LDA_C_VWN";
  double h = 1e-4;  // Angstrom
  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) +
                "/threecenter_dft/3-21G.xml");

  auto build_h2 = [](double bond_length_angstrom) {
    QMMolecule mol(" ", 0);
    std::string xyz_content =
        "2\n\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 " +
        std::to_string(bond_length_angstrom) + "\n";
    std::string tmp_path = "/tmp/xtp_test_xcgradient_h2.xyz";
    std::ofstream out(tmp_path);
    out << xyz_content;
    out.close();
    mol.LoadFromFile(tmp_path);
    return mol;
  };

  QMMolecule mol0 = build_h2(bond_length);
  BasisSet basis0;
  basis0.Load(std::string(XTP_TEST_DATA_FOLDER) +
              "/threecenter_dft/3-21G.xml");
  AOBasis aobasis0;
  aobasis0.Fill(basis0, mol0);
  Vxc_Grid grid0;
  grid0.GridSetup("coarse", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  // No real converged density for bare H2 at hand here -- reuse the
  // SAME "arbitrary fixed matrix" reasoning already used throughout
  // this branch (RI-J/RI-K): PulayGradient/GridWeightGradient only need
  // A fixed density matrix, not a converged one, to test the gradient
  // ASSEMBLY formula itself. Build a simple fixed, symmetric one.
  Index n_bf = aobasis0.AOBasisSize();
  Eigen::MatrixXd dmat_random = Eigen::MatrixXd::Random(n_bf, n_bf);
  // NOTE: dmat = 0.5*(dmat+dmat.transpose()) is a classic Eigen aliasing
  // trap -- assigning an expression that reads from dmat while
  // simultaneously overwriting it can give a result that's only
  // approximately (not exactly) symmetric, tripping PulayGradient's
  // exact isApprox symmetry assertion. Using a separate source variable
  // (dmat_random) for the RHS avoids the aliasing entirely.
  Eigen::MatrixXd dmat = 0.5 * (dmat_random + dmat_random.transpose());

  Eigen::MatrixXd pulay_grad = vxc0.PulayGradient(dmat, aobasis0);
  Eigen::MatrixXd weight_grad = vxc0.GridWeightGradient(dmat, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  std::cout << "H2 test -- Pulay:\n" << pulay_grad << std::endl;
  std::cout << "H2 test -- Weight:\n" << weight_grad << std::endl;
  std::cout << "H2 test -- Total:\n" << total_grad << std::endl;

  Eigen::Vector3d sum = total_grad.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  QMMolecule mol_plus = build_h2(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basis0, mol_plus);
  Vxc_Grid grid_plus;
  grid_plus.GridSetup("coarse", mol_plus, aobasis_plus);
  Vxc_Potential<Vxc_Grid> vxc_plus(grid_plus);
  vxc_plus.setXCfunctional(functional);
  double e_plus = vxc_plus.IntegrateVXC(dmat).energy();

  QMMolecule mol_minus = build_h2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basis0, mol_minus);
  Vxc_Grid grid_minus;
  grid_minus.GridSetup("coarse", mol_minus, aobasis_minus);
  Vxc_Potential<Vxc_Grid> vxc_minus(grid_minus);
  vxc_minus.setXCfunctional(functional);
  double e_minus = vxc_minus.IntegrateVXC(dmat).energy();

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_deriv =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;
  double analytic = total_grad(1, 2);
  std::cout << "H2 test -- analytic dE_xc/dz(atom1)=" << analytic
             << " finite-difference=" << finite_diff_deriv << std::endl;

  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
