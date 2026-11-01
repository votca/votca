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

// Finite-difference tests of the XC nuclear gradient at fixed density
// matrix: PulayGradient + GridWeightGradient (RKS) and PulayGradientUKS +
// GridWeightGradientUKS, compared with central differences of the XC energy
// of IntegrateVXC/IntegrateVXCSpin, rebuilding basis and grid at each
// displaced geometry. Neither term alone matches the energy derivative, so
// their sum is tested, plus translational invariance. Systems: methane
// (3-21G, dmat.mm from test_vxc_potential.cc, "coarse" grid) and H2 (fixed
// model densities, "medium" grid), with LDA (XC_LDA_X XC_LDA_C_VWN) and GGA
// (XC_GGA_X_PBE XC_GGA_C_PBE).

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

// Methane geometry from the test data, displaced by the builder below.
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

// Copy of the base geometry with atom atom_idx shifted by dz (Angstrom)
// along z.
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
  double h = 1e-4;           // Angstrom
  Index displaced_atom = 1;  // first H (index 1, methane is C,H,H,H,H)

  auto base_atoms = LoadBaseGeometry();

  QMMolecule mol0 = BuildDisplacedMethane(base_atoms, displaced_atom, 0.0);
  AOBasis aobasis0 = BuildBasis(mol0);

  Eigen::MatrixXd dmat = tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/dmat.mm");

  Vxc_Grid grid0;
  grid0.GridSetup("coarse", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  Eigen::MatrixXd pulay_grad = vxc0.PulayGradient(dmat, aobasis0);
  Eigen::MatrixXd weight_grad = vxc0.GridWeightGradient(dmat, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  // translational invariance of each term separately
  Eigen::Vector3d pulay_sum_check = pulay_grad.colwise().sum();
  Eigen::Vector3d weight_sum_check = weight_grad.colwise().sum();
  std::cerr << "[xc_gradient_finite_difference diagnostic] "
            << "Pulay alone t.i. sum = " << pulay_sum_check.transpose()
            << "\nWeight alone t.i. sum = " << weight_sum_check.transpose()
            << std::endl;

  // translational invariance of the total
  Eigen::Vector3d sum = total_grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    std::cout << "WARNING: translational invariance check failed, sum="
              << sum.transpose() << std::endl;
    std::cout << "Pulay term alone:\n" << pulay_grad << std::endl;
    std::cout << "Weight term alone:\n" << weight_grad << std::endl;
  }
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  // Finite-difference reference: XC energy of IntegrateVXC at fixed
  // density matrix, with basis and grid rebuilt at each displaced geometry.
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

  constexpr double kAngstromPerBohr = 0.52917721090380;
  double finite_diff_deriv = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;

  double analytic = total_grad(displaced_atom, 2);  // z-component
  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_xc/dz(atom" << displaced_atom << "): " << analytic
              << std::endl;
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

// H2 (2 atoms, few grid points) with LDA.
BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference_h2) {
  libint2::initialize();

  const std::string functional = "XC_LDA_X XC_LDA_C_VWN";
  double h = 1e-4;            // Angstrom
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
  basis0.Load(std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml");
  AOBasis aobasis0;
  aobasis0.Fill(basis0, mol0);
  Vxc_Grid grid0;
  grid0.GridSetup("medium", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  // Fixed symmetric model density; the gradient formulas hold for any fixed
  // density matrix.
  Index n_bf = aobasis0.AOBasisSize();
  Eigen::MatrixXd dmat_random = Eigen::MatrixXd::Random(n_bf, n_bf);
  // Symmetrize from a separate matrix: dmat = 0.5*(dmat+dmat.transpose())
  // aliases and need not be exactly symmetric (PulayGradient asserts it).
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
  grid_plus.GridSetup("medium", mol_plus, aobasis_plus);
  Vxc_Potential<Vxc_Grid> vxc_plus(grid_plus);
  vxc_plus.setXCfunctional(functional);
  double e_plus = vxc_plus.IntegrateVXC(dmat).energy();

  QMMolecule mol_minus = build_h2(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basis0, mol_minus);
  Vxc_Grid grid_minus;
  grid_minus.GridSetup("medium", mol_minus, aobasis_minus);
  Vxc_Potential<Vxc_Grid> vxc_minus(grid_minus);
  vxc_minus.setXCfunctional(functional);
  double e_minus = vxc_minus.IntegrateVXC(dmat).energy();

  constexpr double kAngstromPerBohr = 0.52917721090380;
  double finite_diff_deriv = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
  double analytic = total_grad(1, 2);
  std::cout << "H2 test -- analytic dE_xc/dz(atom1)=" << analytic
            << " finite-difference=" << finite_diff_deriv << std::endl;

  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// UKS analogue of xc_gradient_finite_difference_h2 (LDA): PulayGradientUKS +
// GridWeightGradientUKS for fixed alpha and beta density matrices against
// finite differences of IntegrateVXCSpin.
BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference_h2_uks) {
  libint2::initialize();

  const std::string functional = "XC_LDA_X XC_LDA_C_VWN";
  double h = 1e-4;            // Angstrom
  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium

  auto build_h2_uks = [](double bond_length_angstrom) {
    QMMolecule mol(" ", 0);
    std::string xyz_content =
        "2\n\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 " +
        std::to_string(bond_length_angstrom) + "\n";
    std::string tmp_path = "/tmp/xtp_test_xcgradient_h2_uks.xyz";
    std::ofstream out(tmp_path);
    out << xyz_content;
    out.close();
    mol.LoadFromFile(tmp_path);
    return mol;
  };

  QMMolecule mol0 = build_h2_uks(bond_length);
  BasisSet basis0;
  basis0.Load(std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml");
  AOBasis aobasis0;
  aobasis0.Fill(basis0, mol0);
  Vxc_Grid grid0;
  grid0.GridSetup("medium", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  Index n_bf = aobasis0.AOBasisSize();
  // Positive semidefinite (C*C^T) spin densities: for arbitrary symmetric
  // matrices rho_a or rho_b can be negative at points that pass the cutoff
  // on rho_a + rho_b, where libxc is not differentiable.
  Eigen::MatrixXd C_alpha = Eigen::MatrixXd::Random(n_bf, 1);
  Eigen::MatrixXd dmat_alpha = C_alpha * C_alpha.transpose();
  Eigen::MatrixXd C_beta = Eigen::MatrixXd::Random(n_bf, 1);
  Eigen::MatrixXd dmat_beta = C_beta * C_beta.transpose();

  Eigen::MatrixXd pulay_grad =
      vxc0.PulayGradientUKS(dmat_alpha, dmat_beta, aobasis0);
  Eigen::MatrixXd weight_grad =
      vxc0.GridWeightGradientUKS(dmat_alpha, dmat_beta, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  std::cout << "H2 UKS test -- Pulay:\n" << pulay_grad << std::endl;
  std::cout << "H2 UKS test -- Weight:\n" << weight_grad << std::endl;
  std::cout << "H2 UKS test -- Total:\n" << total_grad << std::endl;

  Eigen::Vector3d sum = total_grad.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  QMMolecule mol_plus = build_h2_uks(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basis0, mol_plus);
  Vxc_Grid grid_plus;
  grid_plus.GridSetup("medium", mol_plus, aobasis_plus);
  Vxc_Potential<Vxc_Grid> vxc_plus(grid_plus);
  vxc_plus.setXCfunctional(functional);
  double e_plus = vxc_plus.IntegrateVXCSpin(dmat_alpha, dmat_beta).energy;

  QMMolecule mol_minus = build_h2_uks(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basis0, mol_minus);
  Vxc_Grid grid_minus;
  grid_minus.GridSetup("medium", mol_minus, aobasis_minus);
  Vxc_Potential<Vxc_Grid> vxc_minus(grid_minus);
  vxc_minus.setXCfunctional(functional);
  double e_minus = vxc_minus.IntegrateVXCSpin(dmat_alpha, dmat_beta).energy;

  constexpr double kAngstromPerBohr = 0.52917721090380;
  double finite_diff_deriv = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
  double analytic = total_grad(1, 2);
  std::cout << "H2 UKS test -- analytic dE_xc/dz(atom1)=" << analytic
            << " finite-difference=" << finite_diff_deriv << std::endl;

  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  if (!matches) {
    std::cout << "NOTE: if the translational-invariance check above "
                 "failed too, the bug is likely in the shared SSW-"
                 "derivative logic copy-adapted from the restricted "
                 "GridWeightGradient. If T.I. held but this still fails, "
                 "check the basis-type/translation-type prefactor "
                 "construction in PulayGradientUKS specifically."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// As xc_gradient_finite_difference_h2 with GGA (XC_GGA_X_PBE XC_GGA_C_PBE):
// tests the sigma terms of PulayGradient (basis and grid-translation terms,
// using AOShell::EvalAOspaceHessian), which vanish for LDA.
BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference_h2_gga) {
  libint2::initialize();

  const std::string functional = "XC_GGA_X_PBE XC_GGA_C_PBE";
  double h = 1e-4;            // Angstrom
  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium

  auto build_h2_gga = [](double bond_length_angstrom) {
    QMMolecule mol(" ", 0);
    std::string xyz_content =
        "2\n\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 " +
        std::to_string(bond_length_angstrom) + "\n";
    std::string tmp_path = "/tmp/xtp_test_xcgradient_h2_gga.xyz";
    std::ofstream out(tmp_path);
    out << xyz_content;
    out.close();
    mol.LoadFromFile(tmp_path);
    return mol;
  };

  QMMolecule mol0 = build_h2_gga(bond_length);
  BasisSet basis0;
  basis0.Load(std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml");
  AOBasis aobasis0;
  aobasis0.Fill(basis0, mol0);
  Vxc_Grid grid0;
  grid0.GridSetup("medium", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  Index n_bf = aobasis0.AOBasisSize();
  Eigen::MatrixXd dmat_random = Eigen::MatrixXd::Random(n_bf, n_bf);
  Eigen::MatrixXd dmat = 0.5 * (dmat_random + dmat_random.transpose());

  Eigen::MatrixXd pulay_grad = vxc0.PulayGradient(dmat, aobasis0);
  Eigen::MatrixXd weight_grad = vxc0.GridWeightGradient(dmat, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  std::cout << "H2 GGA test -- Pulay:\n" << pulay_grad << std::endl;
  std::cout << "H2 GGA test -- Weight:\n" << weight_grad << std::endl;
  std::cout << "H2 GGA test -- Total:\n" << total_grad << std::endl;

  Eigen::Vector3d sum = total_grad.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  QMMolecule mol_plus = build_h2_gga(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basis0, mol_plus);
  Vxc_Grid grid_plus;
  grid_plus.GridSetup("medium", mol_plus, aobasis_plus);
  Vxc_Potential<Vxc_Grid> vxc_plus(grid_plus);
  vxc_plus.setXCfunctional(functional);
  double e_plus = vxc_plus.IntegrateVXC(dmat).energy();

  QMMolecule mol_minus = build_h2_gga(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basis0, mol_minus);
  Vxc_Grid grid_minus;
  grid_minus.GridSetup("medium", mol_minus, aobasis_minus);
  Vxc_Potential<Vxc_Grid> vxc_minus(grid_minus);
  vxc_minus.setXCfunctional(functional);
  double e_minus = vxc_minus.IntegrateVXC(dmat).energy();

  constexpr double kAngstromPerBohr = 0.52917721090380;
  double finite_diff_deriv = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
  double analytic = total_grad(1, 2);
  std::cout << "H2 GGA test -- analytic dE_xc/dz(atom1)=" << analytic
            << " finite-difference=" << finite_diff_deriv << std::endl;

  bool matches_gga =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  if (!matches_gga) {
    std::cout << "NOTE: if this fails, first check whether the LDA tests "
                 "above still pass -- if they do, the bug is isolated to "
                 "the sigma-dependent terms specifically (Gmat/s_vec/"
                 "Hessian_rho construction in PulayGradient, or the "
                 "underlying EvalAOspaceHessian values themselves, "
                 "already separately validated in test_aoshell_hessian.cc "
                 "but worth re-checking if this fails)."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches_gga, true);

  libint2::finalize();
}

// As xc_gradient_finite_difference with GGA: the sigma terms for several
// atoms whose basis functions and SSW partitions overlap.
BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference_gga) {
  libint2::initialize();

  const std::string functional = "XC_GGA_X_PBE XC_GGA_C_PBE";
  double h = 1e-4;           // Angstrom
  Index displaced_atom = 1;  // first H (index 1, methane is C,H,H,H,H)

  auto base_atoms = LoadBaseGeometry();

  QMMolecule mol0 = BuildDisplacedMethane(base_atoms, displaced_atom, 0.0);
  AOBasis aobasis0 = BuildBasis(mol0);

  // same density matrix as the LDA methane test
  Eigen::MatrixXd dmat = tools::EigenIO_MatrixMarket::ReadMatrix(
      std::string(XTP_TEST_DATA_FOLDER) + "/vxc_potential/dmat.mm");

  Vxc_Grid grid0;
  grid0.GridSetup("coarse", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  Eigen::MatrixXd pulay_grad = vxc0.PulayGradient(dmat, aobasis0);
  Eigen::MatrixXd weight_grad = vxc0.GridWeightGradient(dmat, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  Eigen::Vector3d pulay_sum_check = pulay_grad.colwise().sum();
  Eigen::Vector3d weight_sum_check = weight_grad.colwise().sum();
  std::cerr << "[xc_gradient_finite_difference_gga diagnostic] "
            << "Pulay alone t.i. sum = " << pulay_sum_check.transpose()
            << "\nWeight alone t.i. sum = " << weight_sum_check.transpose()
            << std::endl;

  Eigen::Vector3d sum = total_grad.colwise().sum();
  if (sum.cwiseAbs().maxCoeff() > 1e-4) {
    std::cout << "WARNING: translational invariance check failed, sum="
              << sum.transpose() << std::endl;
    std::cout << "Pulay term alone:\n" << pulay_grad << std::endl;
    std::cout << "Weight term alone:\n" << weight_grad << std::endl;
  }
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

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

  constexpr double kAngstromPerBohr = 0.52917721090380;
  double finite_diff_deriv = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;

  double analytic = total_grad(displaced_atom, 2);  // z-component
  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_xc/dz(atom" << displaced_atom << "): " << analytic
              << std::endl;
    std::cout << "  Pulay contribution: " << pulay_grad(displaced_atom, 2)
              << std::endl;
    std::cout << "  Weight contribution: " << weight_grad(displaced_atom, 2)
              << std::endl;
    std::cout << "Finite-difference: " << finite_diff_deriv << std::endl;
    std::cout << "NOTE: if this fails while the LDA methane test and the "
                 "H2 GGA test both pass, that isolates the bug to "
                 "something specific to MULTIPLE atoms' sigma-dependent "
                 "contributions interacting (e.g. the Gmat/Hessian_rho "
                 "sums spanning several atoms' basis functions at one "
                 "grid point) -- something the 2-atom H2 GGA test cannot "
                 "exercise at all."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// GGA analogue of xc_gradient_finite_difference_h2_uks: tests the
// sigma_aa/sigma_ab/sigma_bb terms of PulayGradientUKS/GridWeightGradientUKS.
// Positive semidefinite spin densities, as in the LDA UKS test.
BOOST_AUTO_TEST_CASE(xc_gradient_finite_difference_h2_uks_gga) {
  libint2::initialize();

  const std::string functional = "XC_GGA_X_PBE XC_GGA_C_PBE";
  double h = 1e-4;            // Angstrom
  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium

  auto build_h2_uks_gga = [](double bond_length_angstrom) {
    QMMolecule mol(" ", 0);
    std::string xyz_content =
        "2\n\n"
        "H 0.0 0.0 0.0\n"
        "H 0.0 0.0 " +
        std::to_string(bond_length_angstrom) + "\n";
    std::string tmp_path = "/tmp/xtp_test_xcgradient_h2_uks_gga.xyz";
    std::ofstream out(tmp_path);
    out << xyz_content;
    out.close();
    mol.LoadFromFile(tmp_path);
    return mol;
  };

  QMMolecule mol0 = build_h2_uks_gga(bond_length);
  BasisSet basis0;
  basis0.Load(std::string(XTP_TEST_DATA_FOLDER) + "/threecenter_dft/3-21G.xml");
  AOBasis aobasis0;
  aobasis0.Fill(basis0, mol0);
  Vxc_Grid grid0;
  grid0.GridSetup("medium", mol0, aobasis0);
  Vxc_Potential<Vxc_Grid> vxc0(grid0);
  vxc0.setXCfunctional(functional);

  Index n_bf = aobasis0.AOBasisSize();
  Eigen::MatrixXd C_alpha = Eigen::MatrixXd::Random(n_bf, 1);
  Eigen::MatrixXd dmat_alpha = C_alpha * C_alpha.transpose();
  Eigen::MatrixXd C_beta = Eigen::MatrixXd::Random(n_bf, 1);
  Eigen::MatrixXd dmat_beta = C_beta * C_beta.transpose();

  Eigen::MatrixXd pulay_grad =
      vxc0.PulayGradientUKS(dmat_alpha, dmat_beta, aobasis0);
  Eigen::MatrixXd weight_grad =
      vxc0.GridWeightGradientUKS(dmat_alpha, dmat_beta, mol0);
  Eigen::MatrixXd total_grad = pulay_grad + weight_grad;

  std::cout << "H2 UKS GGA test -- Pulay:\n" << pulay_grad << std::endl;
  std::cout << "H2 UKS GGA test -- Weight:\n" << weight_grad << std::endl;
  std::cout << "H2 UKS GGA test -- Total:\n" << total_grad << std::endl;

  Eigen::Vector3d sum = total_grad.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  QMMolecule mol_plus = build_h2_uks_gga(bond_length + h);
  AOBasis aobasis_plus;
  aobasis_plus.Fill(basis0, mol_plus);
  Vxc_Grid grid_plus;
  grid_plus.GridSetup("medium", mol_plus, aobasis_plus);
  Vxc_Potential<Vxc_Grid> vxc_plus(grid_plus);
  vxc_plus.setXCfunctional(functional);
  double e_plus = vxc_plus.IntegrateVXCSpin(dmat_alpha, dmat_beta).energy;

  QMMolecule mol_minus = build_h2_uks_gga(bond_length - h);
  AOBasis aobasis_minus;
  aobasis_minus.Fill(basis0, mol_minus);
  Vxc_Grid grid_minus;
  grid_minus.GridSetup("medium", mol_minus, aobasis_minus);
  Vxc_Potential<Vxc_Grid> vxc_minus(grid_minus);
  vxc_minus.setXCfunctional(functional);
  double e_minus = vxc_minus.IntegrateVXCSpin(dmat_alpha, dmat_beta).energy;

  constexpr double kAngstromPerBohr = 0.52917721090380;
  double finite_diff_deriv = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
  double analytic = total_grad(1, 2);
  std::cout << "H2 UKS GGA test -- analytic dE_xc/dz(atom1)=" << analytic
            << " finite-difference=" << finite_diff_deriv << std::endl;

  bool matches =
      std::abs(finite_diff_deriv - analytic) < 1e-3 * std::abs(analytic);
  if (!matches) {
    std::cout << "NOTE: if the LDA UKS test above still passes, the bug "
                 "is isolated to the sigma-dependent terms specifically "
                 "(Gmat_a/Gmat_b/Hessian_rho_a/Hessian_rho_b/V_a/V_b "
                 "construction in PulayGradientUKS, or the real "
                 "sigma_aa/ab/bb computation in GridWeightGradientUKS)."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
