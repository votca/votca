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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE dftgradient_test

// Standard includes
#include <cmath>
#include <fstream>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include <votca/xtp/dftgradient.h>
#include <votca/xtp/qmmolecule.h>

using namespace votca::xtp;
using namespace votca;

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

BOOST_AUTO_TEST_SUITE_END()
