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
// ===========================================================================

#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE dftengine_forces_test

// Standard includes
#include <cmath>
#include <fstream>
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
// bond length, returning the converged Orbitals object.
Orbitals RunSCF(double bond_length_angstrom) {
  DFTEngine dft;
  Orbitals orb;
  orb.QMAtoms() = BuildH2(bond_length_angstrom);

  std::string xml_path = "/tmp/xtp_test_dftengine_forces.xml";
  std::ofstream xml(xml_path);
  xml << "<dftpackage>\n";
  xml << "<spin>1</spin>\n";
  xml << "<name>xtp</name>\n";
  xml << "<charge>0</charge>\n";
  xml << "<functional>XC_GGA_X_PBE XC_GGA_C_PBE</functional>\n";
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

BOOST_AUTO_TEST_SUITE_END()
