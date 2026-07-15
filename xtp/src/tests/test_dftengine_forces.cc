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
// REPO NOTE: force computation is opt-in (compute_forces_ defaults to
// false, since it adds real cost to every converged SCF) -- RunSCF
// below explicitly sets <compute_forces>true</compute_forces> in its
// options XML for exactly this reason. If this test ever starts failing
// with hasForces()==false and no other error, check that setting first
// before suspecting a real regression.
//
// STATUS: written but NOT yet run. This is the most decisive validation
// possible for DFTEngine::ComputeAndStoreForces -- a genuine end-to-end
// test: real SCF convergence (not a fixed/arbitrary density matrix, as
// every earlier gradient test in this branch used), real options parsing,
// the actual wiring inside EvaluateClosedShell, the sign convention, and
// the RI-J/non-hybrid scoping guard, all exercised together. Every
// individual gradient term has already been separately validated
// (test_dftgradient.cc, test_xcgradient.cc); this test checks the
// ASSEMBLY and INTEGRATION into DFTEngine specifically, which nothing
// else in this branch has touched.
//
// Runs three FULL SCF calculations (base geometry, +h, -h) -- genuinely
// slow compared to everything else in this branch (which reused a single
// fixed density matrix across displaced geometries rather than
// reconverging SCF each time), but this is the only way to check the
// REAL, self-consistent gradient against a REAL finite difference of the
// REAL total SCF energy, not an XC-only or RI-J-only energy.
//
// Uses a non-hybrid GGA functional (XC_GGA_X_PBE XC_GGA_C_PBE) and a
// real auxiliary basis (diabatization/aux-def2-svp.xml, covers H) to
// stay within ComputeAndStoreForces's explicit scope (RI + non-hybrid) --
// deliberately NOT reusing test_dftengine.cc's existing hybrid-functional
// water setup, which ComputeAndStoreForces would (correctly) skip force
// computation for.
//
// Also includes forces_finite_difference_uks: the same end-to-end
// methodology applied to a genuinely open-shell system (H2+, doublet),
// exercising DFTEngine::ComputeAndStoreForcesUKS and the actual
// EvaluateUKS path for the first time end to end.
// ===========================================================================

#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE dftengine_forces_test

// Standard includes
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/dftengine.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(dftengine_forces_test)

QMMolecule BuildH2(double bond_length_angstrom) {
  QMMolecule mol(" ", 0);
  std::string xyz_content =
      "2\n\n"
      "H 0.0 0.0 0.0\n"
      "H 0.0 0.0 " +
      std::to_string(bond_length_angstrom) + "\n";
  std::string tmp_path = "/tmp/xtp_test_dftengine_forces_h2.xyz";
  std::ofstream out(tmp_path);
  out << xyz_content;
  out.close();
  mol.LoadFromFile(tmp_path);
  return mol;
}

// Writes the options XML and runs one full SCF calculation at the given
// bond length, returning the converged Orbitals object. functional
// defaults to the existing non-hybrid GGA choice; pass a hybrid
// functional (e.g. "XC_HYB_GGA_XC_PBEH", the same name already used and
// confirmed valid in test_dftengine.cc's dft_full test) to exercise the
// RI-K/hybrid path in ComputeAndStoreForces. spin/charge default to the
// existing closed-shell H2 case (singlet, neutral); pass spin=2, charge=1
// for H2+ (doublet, one unpaired electron) to exercise the UKS path.
Orbitals RunSCF(double bond_length_angstrom,
                const std::string& functional = "XC_GGA_X_PBE XC_GGA_C_PBE",
                int spin = 1, int charge = 0) {
  DFTEngine dft;
  Orbitals orb;
  orb.QMAtoms() = BuildH2(bond_length_angstrom);

  // Distinct temp file per (functional, spin, charge), so different test
  // cases' options files can't collide with or shadow each other if
  // tests happen to run concurrently.
  std::string xml_path =
      "/tmp/xtp_test_dftengine_forces_" +
      std::to_string(std::hash<std::string>{}(
          functional + std::to_string(spin) + std::to_string(charge))) +
      ".xml";
  std::ofstream xml(xml_path);
  xml << "<dftpackage>\n";
  xml << "<spin>" << spin << "</spin>\n";
  xml << "<name>xtp</name>\n";
  xml << "<charge>" << charge << "</charge>\n";
  xml << "<functional>" << functional << "</functional>\n";
  xml << "<basisset>" << XTP_TEST_DATA_FOLDER
      << "/threecenter_dft/3-21G.xml</basisset>\n";
  xml << "<auxbasisset>" << XTP_TEST_DATA_FOLDER
      << "/diabatization/aux-def2-svp.xml</auxbasisset>\n";
  xml << "<initial_guess>independent</initial_guess>\n";
  xml << "<xtpdft>\n";
  xml << "<screening_eps>1e-9</screening_eps>\n";
  xml << "<fock_matrix_reset>5</fock_matrix_reset>\n";
  xml << "<convergence>\n";
  xml << "    <energy>1e-9</energy>\n";
  xml << "    <method>DIIS</method>\n";
  xml << "    <DIIS_start>0.002</DIIS_start>\n";
  xml << "    <ADIIS_start>0.8</ADIIS_start>\n";
  xml << "    <DIIS_length>20</DIIS_length>\n";
  xml << "    <levelshift>0.0</levelshift>\n";
  xml << "    <levelshift_end>0.2</levelshift_end>\n";
  xml << "    <max_iterations>100</max_iterations>\n";
  xml << "    <error>1e-8</error>\n";
  xml << "    <DIIS_maxout>false</DIIS_maxout>\n";
  xml << "    <mixing>0.7</mixing>\n";
  xml << "</convergence>\n";
  xml << "<integration_grid>xcoarse</integration_grid>\n";
  xml << "<max_iterations>200</max_iterations>\n";
  xml << "<compute_forces>true</compute_forces>\n";
  xml << "</xtpdft>\n";
  xml << "</dftpackage>\n";
  xml.close();

  votca::tools::Property prop;
  prop.LoadFromXML(xml_path);

  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
  bool converged = dft.Evaluate(orb);
  BOOST_REQUIRE_EQUAL(converged, true);

  return orb;
}

BOOST_AUTO_TEST_CASE(forces_finite_difference) {
  libint2::initialize();

  double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
  double h = 1e-3;  // Angstrom -- larger than the fixed-density-matrix
                     // tests elsewhere in this branch (which could afford
                     // 1e-4 or smaller), since SCF reconvergence at each
                     // displaced geometry introduces its own numerical
                     // noise floor (convergence threshold, grid
                     // discretization) that a too-small h would be
                     // swamped by.

  Orbitals orb0 = RunSCF(bond_length);
  BOOST_REQUIRE_EQUAL(orb0.hasForces(), true);
  Eigen::MatrixXd forces = orb0.getForces();

  std::cout << "H2 SCF forces:\n" << forces << std::endl;
  std::cout << "H2 SCF energy: " << orb0.getDFTTotalEnergy() << std::endl;

  // Translational invariance, independent of the finite-difference
  // comparison below.
  Eigen::Vector3d sum = forces.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  Orbitals orb_plus = RunSCF(bond_length + h);
  Orbitals orb_minus = RunSCF(bond_length - h);

  double e_plus = orb_plus.getDFTTotalEnergy();
  double e_minus = orb_minus.getDFTTotalEnergy();

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_dEdz =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;
  // force = -dE/dR, matching ComputeAndStoreForces's own convention.
  double finite_diff_force = -finite_diff_dEdz;

  double analytic_force = forces(1, 2);  // second H, z-component

  std::cout << "Analytic force on atom 1 (z): " << analytic_force
             << std::endl;
  std::cout << "Finite-difference force on atom 1 (z): "
             << finite_diff_force << std::endl;

  bool matches = std::abs(finite_diff_force - analytic_force) <
                 1e-2 * std::abs(analytic_force);
  if (!matches) {
    std::cout << "NOTE: if this fails, check individually: (1) whether "
                 "orb0.hasForces() was even true (a scoping guard may "
                 "have silently skipped force computation -- check the "
                 "log for a 'Skipping force calculation' message), (2) "
                 "whether the translational-invariance check above "
                 "passed (if not, the bug is likely in the assembly/"
                 "sign convention in ComputeAndStoreForces itself, not "
                 "in any individual already-validated term), (3) "
                 "whether tightening SCF convergence or h changes the "
                 "result substantially (would suggest reconvergence "
                 "noise rather than a real formula bug)."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// Same pattern as forces_finite_difference above, but with a hybrid
// functional (XC_HYB_GGA_XC_PBEH -- same name already used and
// confirmed valid in test_dftengine.cc's dft_full test), to exercise
// the RI-K/exact-exchange path in ComputeAndStoreForces end to end for
// the first time. Everything else in the pipeline (nuclear repulsion,
// one-electron, overlap Pulay force, RI-J, XC) is already validated by
// the non-hybrid test above; this specifically checks the newly-added
// ScaHFX_ * RIKGradient(...) contribution, on top of a real,
// self-consistently-converged hybrid SCF (RIKGradient's underlying
// formula was separately validated against real production code with a
// FIXED, arbitrary MO coefficient matrix in test_dftgradient.cc -- this
// is the first check of it with genuine, converged occupied MOs feeding
// back into its own Fock matrix).
//
// STATUS: written but NOT yet run.
BOOST_AUTO_TEST_CASE(forces_finite_difference_hybrid) {
  libint2::initialize();

  double bond_length = 0.74;  // Angstrom
  double h = 1e-3;            // Angstrom, same reasoning as above
  const std::string functional = "XC_HYB_GGA_XC_PBEH";

  Orbitals orb0 = RunSCF(bond_length, functional);
  BOOST_REQUIRE_EQUAL(orb0.hasForces(), true);
  Eigen::MatrixXd forces = orb0.getForces();

  std::cout << "H2 hybrid SCF forces:\n" << forces << std::endl;
  std::cout << "H2 hybrid SCF energy: " << orb0.getDFTTotalEnergy()
             << std::endl;

  Eigen::Vector3d sum = forces.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  Orbitals orb_plus = RunSCF(bond_length + h, functional);
  Orbitals orb_minus = RunSCF(bond_length - h, functional);

  double e_plus = orb_plus.getDFTTotalEnergy();
  double e_minus = orb_minus.getDFTTotalEnergy();

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_dEdz =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;
  double finite_diff_force = -finite_diff_dEdz;

  double analytic_force = forces(1, 2);

  std::cout << "Analytic force on atom 1 (z): " << analytic_force
             << std::endl;
  std::cout << "Finite-difference force on atom 1 (z): "
             << finite_diff_force << std::endl;

  bool matches = std::abs(finite_diff_force - analytic_force) <
                 1e-2 * std::abs(analytic_force);
  if (!matches) {
    std::cout << "NOTE: since the non-hybrid test above already "
                 "validates everything except the RI-K term, a failure "
                 "here most likely points at the ScaHFX_ * RIKGradient "
                 "contribution specifically -- e.g. an issue with how "
                 "C_occ (built for the overlap Pulay term, reused here) "
                 "is threaded through, or a remaining scale/sign issue "
                 "not caught by the FIXED-MO-matrix test in "
                 "test_dftgradient.cc."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

// Same pattern as forces_finite_difference above, but for a genuinely
// open-shell (UKS) system: H2+ (one electron total, doublet -- spin=2
// meaning multiplicity 2, i.e. one unpaired electron), charge=+1. This
// is the most decisive validation possible for
// DFTEngine::ComputeAndStoreForcesUKS -- real SCF convergence down the
// actual EvaluateUKS path (num_alpha=1, num_beta=0 for H2+, naturally
// unequal, no need for force_uks_path_), the actual wiring inside
// EvaluateUKS, and the newly-added XC gradient (PulayGradientUKS +
// GridWeightGradientUKS) all exercised together for the first time.
// Every individual UKS gradient term has already been separately
// validated (test_dftengine_private.cc, test_xcgradient.cc); this test
// checks the ASSEMBLY and INTEGRATION specifically, the same gap that,
// on the RKS side, turned out to hide two genuinely missing physics
// terms (kinetic+nuclear-attraction, then the overlap Pulay force) that
// no amount of per-term testing could have caught.
//
// STATUS: written but NOT yet run.
BOOST_AUTO_TEST_CASE(forces_finite_difference_uks) {
  libint2::initialize();

  double bond_length = 1.06;  // Angstrom, roughly H2+ equilibrium
                              // (longer than neutral H2's ~0.74 A,
                              // consistent with a weaker one-electron
                              // bond)
  double h = 1e-3;  // Angstrom, same reasoning as the RKS tests above
  const std::string functional = "XC_GGA_X_PBE XC_GGA_C_PBE";
  int spin = 2;    // doublet
  int charge = 1;  // H2+

  Orbitals orb0 = RunSCF(bond_length, functional, spin, charge);
  BOOST_REQUIRE_EQUAL(orb0.hasForces(), true);
  Eigen::MatrixXd forces = orb0.getForces();

  std::cout << "H2+ UKS SCF forces:\n" << forces << std::endl;
  std::cout << "H2+ UKS SCF energy: " << orb0.getDFTTotalEnergy()
             << std::endl;

  Eigen::Vector3d sum = forces.colwise().sum();
  BOOST_CHECK_SMALL(sum.cwiseAbs().maxCoeff(), 1e-4);

  Orbitals orb_plus = RunSCF(bond_length + h, functional, spin, charge);
  Orbitals orb_minus = RunSCF(bond_length - h, functional, spin, charge);

  double e_plus = orb_plus.getDFTTotalEnergy();
  double e_minus = orb_minus.getDFTTotalEnergy();

  constexpr double kBohrPerAngstrom = 0.52917721090380;
  double finite_diff_dEdz =
      (e_plus - e_minus) / (2.0 * h) * kBohrPerAngstrom;
  double finite_diff_force = -finite_diff_dEdz;

  double analytic_force = forces(1, 2);

  std::cout << "Analytic force on atom 1 (z): " << analytic_force
             << std::endl;
  std::cout << "Finite-difference force on atom 1 (z): "
             << finite_diff_force << std::endl;

  bool matches = std::abs(finite_diff_force - analytic_force) <
                 1e-2 * std::abs(analytic_force);
  if (!matches) {
    std::cout << "NOTE: if this fails, check individually: (1) whether "
                 "orb0.hasForces() was even true (a scoping guard may "
                 "have silently skipped force computation), (2) whether "
                 "the translational-invariance check above passed (if "
                 "not, the bug is likely in the assembly/sign convention "
                 "in ComputeAndStoreForcesUKS itself, not in any "
                 "individual already-validated term), (3) whether "
                 "tightening SCF convergence or h changes the result "
                 "substantially (would suggest reconvergence noise "
                 "rather than a real formula bug)."
              << std::endl;
  }
  BOOST_CHECK_EQUAL(matches, true);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
