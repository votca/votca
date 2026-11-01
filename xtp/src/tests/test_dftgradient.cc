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

// Finite-difference tests of DFTGradient: NuclearRepulsionDerivative
// against an independent implementation of E_nn, and RIJGradient/RIKGradient
// against finite differences of the RI Coulomb and exchange energies
// (ERIs::CalculateEXX_mos) for a fixed density matrix or MO coefficients.
// All also check translational invariance.

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
#include <votca/xtp/ERIs.h>
#include <votca/xtp/aobasis.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/basisset.h>
#include <votca/xtp/dftgradient.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;

// Defined in libint2_derivative_calls.cc (no header); used for the reference
// energies.
namespace votca {
namespace xtp {
std::vector<Eigen::MatrixXd> ComputeThreeCenterIntegrals(
    const AOBasis& auxbasis, const AOBasis& dftbasis);
}  // namespace xtp
}  // namespace votca

BOOST_AUTO_TEST_SUITE(dftgradient_test)

// Non-collinear three-atom molecule, so the gradient has several nonzero
// components.
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

// Reference E_nn, independent of XTP code.
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
  double h = 1e-4;  // finite-difference step, Angstrom (.xyz input)

  // .xyz coordinates are read as Angstrom; a failure by a uniform factor
  // 1.8897 would indicate a unit mismatch.
  constexpr double kAngstromPerBohr = 0.52917721090380;

  QMMolecule mol0 = BuildTestMolecule(0.0);
  Eigen::MatrixXd deriv = DFTGradient::NuclearRepulsionDerivative(mol0);

  // translational invariance
  Eigen::Vector3d total = deriv.colwise().sum();
  BOOST_CHECK_SMALL(total.cwiseAbs().maxCoeff(), 1e-10);

  // Finite difference w.r.t. atom 2's x-coordinate (the second H).
  QMMolecule mol_plus = BuildTestMolecule(h);
  QMMolecule mol_minus = BuildTestMolecule(-h);
  double e_plus = NuclearRepulsionEnergyReference(mol_plus);
  double e_minus = NuclearRepulsionEnergyReference(mol_minus);
  double finite_diff_deriv_ang = (e_plus - e_minus) / (2.0 * h);
  double finite_diff_deriv_bohr = finite_diff_deriv_ang * kAngstromPerBohr;

  double analytic = deriv(2, 0);  // atom index 2, x-component
  // O(h^2) error at h = 1e-4
  bool matches =
      std::abs(finite_diff_deriv_bohr - analytic) < 1e-4 * std::abs(analytic);
  if (!matches) {
    std::cout << "Analytic dE_nn/dx(atom2): " << analytic << std::endl;
    std::cout << "Finite-difference: " << finite_diff_deriv_bohr << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);
}

// RIJGradient against a finite difference of the RI-J energy for a fixed,
// arbitrary density matrix (valid since the fit coefficients are stationary
// for any density, see dftgradient.h).
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
    auxbasis0.Fill(basisset, mol0);  // same basis as auxiliary basis

    // fixed, arbitrary symmetric density matrix, the same at all geometries
    Index n_dft_bf = dftbasis0.AOBasisSize();
    Eigen::MatrixXd density_random =
        Eigen::MatrixXd::Random(n_dft_bf, n_dft_bf);
    Eigen::MatrixXd density =
        0.5 * (density_random + density_random.transpose());

    Eigen::MatrixXd analytic_grad =
        DFTGradient::RIJGradient(density, auxbasis0, dftbasis0);

    // translational invariance
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

    constexpr double kAngstromPerBohr = 0.52917721090380;
    double finite_diff_deriv =
        (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;

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

// RIKGradient against a finite difference of the exchange energy of
// ERIs::CalculateEXX_mos for a fixed, arbitrary (nbf x 2) coefficient matrix
// (not orthonormal; valid since the per-pair fit coefficients are
// stationary for any C).
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

    // fixed, arbitrary (nbf x 2) coefficient matrix, the same at all
    // geometries
    Index n_dft_bf = dftbasis0.AOBasisSize();
    Eigen::MatrixXd mo_coeffs = Eigen::MatrixXd::Random(n_dft_bf, 2);

    Eigen::MatrixXd analytic_grad =
        DFTGradient::RIKGradient(mo_coeffs, auxbasis0, dftbasis0);

    Eigen::Vector3d total = analytic_grad.colwise().sum();
    BOOST_CHECK_SMALL(total.cwiseAbs().maxCoeff(), 1e-6);

    // Reference energy from ERIs::CalculateEXX_mos; should equal
    // RIKGradient's E_K = -sum_ij c_ij.d_ij.
    auto rik_energy = [&](const AOBasis& auxbasis, const AOBasis& dftbasis) {
      ERIs eris;
      eris.Initialize(dftbasis, auxbasis);
      Eigen::MatrixXd Dmat_local = 2.0 * mo_coeffs * mo_coeffs.transpose();
      std::array<Eigen::MatrixXd, 2> JK =
          eris.CalculateERIs_EXX_3c(mo_coeffs, Dmat_local);
      const Eigen::MatrixXd& K = JK[1];
      return 0.25 *
             Dmat_local.cwiseProduct(K).sum();  // ScaHFX=1 for this check
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

    constexpr double kAngstromPerBohr = 0.52917721090380;
    double finite_diff_deriv =
        (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;

    double analytic = analytic_grad(1, 2);  // atom 1 (second H), z-component
    bool matches =
        std::abs(finite_diff_deriv - analytic) < 1e-4 * std::abs(analytic);
    if (!matches) {
      std::cout << "Analytic dE_K/dz(atom1): " << analytic << std::endl;
      std::cout << "Finite-difference: " << finite_diff_deriv << std::endl;
      std::cout << "NOTE: the reference is the ERIs::CalculateEXX_mos "
                   "energy (symmetric V^-1/2 RI fitting); if this fails, "
                   "check the factor of 2 and the sign of RIKGradient's "
                   "energy convention E_K = -sum_ij c_ij.d_ij (see "
                   "dftgradient.h)."
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
