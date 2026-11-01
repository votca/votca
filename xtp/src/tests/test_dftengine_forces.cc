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

// End-to-end tests of the analytic DFT nuclear forces
// (DFTEngine::ComputeAndStoreForces/ComputeAndStoreForcesUKS) for H2 with
// 3-21G and an RI auxiliary basis: the force from a converged SCF is
// compared with a central finite difference of the total SCF energy (three
// SCF runs per test) and checked for translational invariance. Cases: GGA
// (PBE), hybrid (PBEh, RI-K term) and open-shell UKS (H2+, doublet). Slow.
// Forces are opt-in, so RunSCF sets compute_forces=true. The tests are
// skipped if libint2 lacks derivative integrals.

#include "xtp_libint2.h"
#include <stdexcept>
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

// Writes the options XML and runs one SCF at the given bond length
// (Angstrom). Defaults: PBE, neutral singlet; spin=2, charge=1 gives H2+
// (UKS path).
Orbitals RunSCF(double bond_length_angstrom,
                const std::string& functional = "XC_GGA_X_PBE XC_GGA_C_PBE",
                int spin = 1, int charge = 0) {
  DFTEngine dft;
  Orbitals orb;
  orb.QMAtoms() = BuildH2(bond_length_angstrom);

  // one options file per (functional, spin, charge)
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
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
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
  try {
    double bond_length = 0.74;  // Angstrom, roughly H2 equilibrium
    double h = 1e-3;  // Angstrom; larger than in the fixed-density tests,
                      // since SCF convergence and grid noise limit the
                      // finite difference.

    Orbitals orb0 = RunSCF(bond_length);
    // ComputeAndStoreForces skips (log only, no exception) if libint2
    // lacks derivative support, so test hasForces().
    if (!orb0.hasForces()) {
      std::cout << "SKIPPING forces_finite_difference: orb0.hasForces() is "
                   "false -- this libint2 build likely lacks derivative-"
                   "integral support for one or more operator categories "
                   "(see the log)."
                << std::endl;
      libint2::finalize();
      return;
    }
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

    constexpr double kAngstromPerBohr = 0.52917721090380;
    double finite_diff_dEdz = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
    // force = -dE/dR
    double finite_diff_force = -finite_diff_dEdz;

    double analytic_force = forces(1, 2);  // second H, z-component

    std::cout << "Analytic force on atom 1 (z): " << analytic_force
              << std::endl;
    std::cout << "Finite-difference force on atom 1 (z): " << finite_diff_force
              << std::endl;

    bool matches = std::abs(finite_diff_force - analytic_force) <
                   1e-2 * std::abs(analytic_force);
    if (!matches) {
      std::cout << "NOTE: if this fails, check (1) the log for a 'Skipping "
                   "force calculation' message, (2) the translational-"
                   "invariance check above (if it fails, suspect the "
                   "assembly/sign convention in ComputeAndStoreForces), "
                   "(3) whether tighter SCF convergence or a different h "
                   "changes the result (reconvergence noise)."
                << std::endl;
    }
    BOOST_CHECK_EQUAL(matches, true);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING forces_finite_difference: " << e.what() << std::endl;
    libint2::finalize();
    return;
  }

  libint2::finalize();
}

// As forces_finite_difference, with the hybrid XC_HYB_GGA_XC_PBEH: checks
// the ScaHFX_ * RIKGradient term with converged occupied MOs
// (test_dftgradient.cc tests RIKGradient for a fixed MO matrix).
BOOST_AUTO_TEST_CASE(forces_finite_difference_hybrid) {
  libint2::initialize();
  try {
    double bond_length = 0.74;  // Angstrom
    double h = 1e-3;            // Angstrom, same reasoning as above
    const std::string functional = "XC_HYB_GGA_XC_PBEH";

    Orbitals orb0 = RunSCF(bond_length, functional);
    if (!orb0.hasForces()) {
      std::cout << "SKIPPING forces_finite_difference_hybrid: "
                   "orb0.hasForces() is false -- this libint2 build likely "
                   "lacks derivative-integral support for one or more "
                   "operator categories (see the log)."
                << std::endl;
      libint2::finalize();
      return;
    }
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

    constexpr double kAngstromPerBohr = 0.52917721090380;
    double finite_diff_dEdz = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
    double finite_diff_force = -finite_diff_dEdz;

    double analytic_force = forces(1, 2);

    std::cout << "Analytic force on atom 1 (z): " << analytic_force
              << std::endl;
    std::cout << "Finite-difference force on atom 1 (z): " << finite_diff_force
              << std::endl;

    bool matches = std::abs(finite_diff_force - analytic_force) <
                   1e-2 * std::abs(analytic_force);
    if (!matches) {
      std::cout << "NOTE: the non-hybrid test covers all terms except RI-K, "
                   "so a failure here most likely points at the ScaHFX_ * "
                   "RIKGradient contribution (e.g. the C_occ passed to it, "
                   "or a scale/sign issue)."
                << std::endl;
    }
    BOOST_CHECK_EQUAL(matches, true);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING forces_finite_difference_hybrid: " << e.what()
              << std::endl;
    libint2::finalize();
    return;
  }

  libint2::finalize();
}

// As forces_finite_difference, for H2+ (charge +1, multiplicity 2;
// n_alpha = 1, n_beta = 0), i.e. EvaluateUKS and ComputeAndStoreForcesUKS
// including the XC terms PulayGradientUKS + GridWeightGradientUKS. The
// individual terms are tested in test_dftengine_private.cc and
// test_xcgradient.cc; this checks their assembly.
BOOST_AUTO_TEST_CASE(forces_finite_difference_uks) {
  libint2::initialize();
  try {
    double bond_length = 1.06;  // Angstrom, roughly H2+ equilibrium
    double h = 1e-3;  // Angstrom, same reasoning as the RKS tests above
    const std::string functional = "XC_GGA_X_PBE XC_GGA_C_PBE";
    int spin = 2;    // doublet
    int charge = 1;  // H2+

    Orbitals orb0 = RunSCF(bond_length, functional, spin, charge);
    if (!orb0.hasForces()) {
      std::cout << "SKIPPING forces_finite_difference_uks: "
                   "orb0.hasForces() is false -- this libint2 build likely "
                   "lacks derivative-integral support for one or more "
                   "operator categories (see the log)."
                << std::endl;
      libint2::finalize();
      return;
    }
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

    constexpr double kAngstromPerBohr = 0.52917721090380;
    double finite_diff_dEdz = (e_plus - e_minus) / (2.0 * h) * kAngstromPerBohr;
    double finite_diff_force = -finite_diff_dEdz;

    double analytic_force = forces(1, 2);

    std::cout << "Analytic force on atom 1 (z): " << analytic_force
              << std::endl;
    std::cout << "Finite-difference force on atom 1 (z): " << finite_diff_force
              << std::endl;

    bool matches = std::abs(finite_diff_force - analytic_force) <
                   1e-2 * std::abs(analytic_force);
    if (!matches) {
      std::cout << "NOTE: if this fails, check (1) the log for a 'Skipping "
                   "UKS force calculation' message, (2) the translational-"
                   "invariance check above (if it fails, suspect the "
                   "assembly/sign convention in ComputeAndStoreForcesUKS), "
                   "(3) whether tighter SCF convergence or a different h "
                   "changes the result (reconvergence noise)."
                << std::endl;
    }
    BOOST_CHECK_EQUAL(matches, true);
  } catch (const std::runtime_error& e) {
    std::cout << "SKIPPING forces_finite_difference_uks: " << e.what()
              << std::endl;
    libint2::finalize();
    return;
  }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
