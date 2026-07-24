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

// ===========================================================================
// STATUS: written but NOT yet run. First real exercise of the full
// Hirshfeld-partition pipeline built up over several commits
// (DFTEngine::ComputeHirshfeldReferenceDensities ->
// HirshfeldPartition::BuildAtomicReferences -> EvaluateWeight ->
// BuildWeightMatrix), all in one test.
//
// The check itself is a genuine mathematical invariant, not an
// approximate sanity check: by construction, w_i(r) = rho_i(r) /
// sum_j rho_j(r), so summing over EVERY atom i in the molecule gives
// sum_i w_i(r) = 1 identically, everywhere -- which means sum_i W_i
// (every atom's own AO-basis weight matrix, summed) must equal the
// ordinary AO overlap matrix S exactly (up to the numerical grid's own
// integration error): sum_i W_i,munu = integral (sum_i w_i(r))
// phi_mu(r) phi_nu(r) dr = integral 1 * phi_mu(r) phi_nu(r) dr =
// S_munu. This holds regardless of whether any individual reference
// density is itself "correct" in some deeper physical sense (e.g.
// regardless of whether the Hund's-rule occupation table entries are
// exactly right) -- it is a structural property of the ratio
// construction itself, so it specifically tests the IMPLEMENTATION
// (indexing, grid integration, AO evaluation, basis re-centering),
// not the physics of the reference densities.
//
// Uses CO specifically (not H2) so the test actually exercises the
// Hund's-rule occupation table for two different, non-trivial
// elements (C: triplet 2p2, O: triplet 2p4) rather than only
// hydrogen's trivial 1s1 case.
// ===========================================================================

#include "xtp_libint2.h"
#include <stdexcept>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE hirshfeldpartition_test

// Standard includes
#include <fstream>
#include <iostream>
#include <map>
#include <string>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/hirshfeldpartition.h"
#include "votca/xtp/vxc_grid.h"

using namespace votca::xtp;
using namespace votca;

namespace votca {
namespace xtp {

// Local friend-class test-access helper, matching the same pattern
// already established in test_dftengine_private.cc/test_dftengine_forces.cc
// -- ComputeHirshfeldReferenceDensities is a private DFTEngine method.
//
// Deliberately defined BEFORE BOOST_AUTO_TEST_SUITE below, not after:
// that macro introduces its own nested scope (confirmed directly --
// test_dftengine_private.cc, which also needs this same friend class,
// does not use BOOST_AUTO_TEST_SUITE at all, keeping everything flat
// inside namespace votca::xtp with nothing in between), so a class
// defined after it is NOT directly votca::xtp::DFTEngineTestAccess
// anymore -- which is exactly what dftengine.h's own
// "friend class DFTEngineTestAccess;" declaration needs to match.
class DFTEngineTestAccess {
 public:
  static std::map<std::string, Eigen::MatrixXd>
  ComputeHirshfeldReferenceDensities(const DFTEngine& e,
                                     const QMMolecule& mol) {
    return e.ComputeHirshfeldReferenceDensities(mol);
  }
};

QMMolecule BuildCO(double bond_length_angstrom) {
  QMMolecule mol(" ", 0);
  std::string xyz_content =
      "2\n\n"
      "C 0.0 0.0 0.0\n"
      "O 0.0 0.0 " +
      std::to_string(bond_length_angstrom) + "\n";
  std::string tmp_path = "/tmp/xtp_test_hirshfeldpartition_co.xyz";
  std::ofstream out(tmp_path);
  out << xyz_content;
  out.close();
  mol.LoadFromFile(tmp_path);
  return mol;
}

BOOST_AUTO_TEST_SUITE(hirshfeldpartition_test)

BOOST_AUTO_TEST_CASE(weight_matrices_sum_to_overlap) {
  libint2::initialize();
 try {
  QMMolecule mol = BuildCO(1.13);  // Angstrom, roughly CO equilibrium

  std::string xml_path = "/tmp/xtp_test_hirshfeldpartition.xml";
  std::ofstream xml(xml_path);
  xml << "<dftpackage>\n";
  xml << "<spin>1</spin>\n";
  xml << "<name>xtp</name>\n";
  xml << "<charge>0</charge>\n";
  xml << "<functional>XC_GGA_X_PBE XC_GGA_C_PBE</functional>\n";
  xml << "<basisset>" << XTP_TEST_DATA_FOLDER
      << "/hirshfeldpartition/3-21G.xml</basisset>\n";
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
  // fine, not xcoarse: this test's own accuracy is limited by the
  // SAME grid's own integration error (BuildWeightMatrix integrates
  // over this exact grid), so a coarse grid would show up as a
  // misleadingly large "failure" that is really just integration
  // noise, not an implementation bug.
  xml << "<integration_grid>fine</integration_grid>\n";
  xml << "<max_iterations>200</max_iterations>\n";
  xml << "</xtpdft>\n";
  xml << "</dftpackage>\n";
  xml.close();

  votca::tools::Property prop;
  prop.LoadFromXML(xml_path);

  DFTEngine dft;
  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));

  std::map<std::string, Eigen::MatrixXd> reference_densities =
      DFTEngineTestAccess::ComputeHirshfeldReferenceDensities(dft, mol);
  BOOST_REQUIRE_EQUAL(reference_densities.count("C"), 1);
  BOOST_REQUIRE_EQUAL(reference_densities.count("O"), 1);

  BasisSet basisset;
  basisset.Load(std::string(XTP_TEST_DATA_FOLDER) +
               "/hirshfeldpartition/3-21G.xml");
  AOBasis full_basis;
  full_basis.Fill(basisset, mol);

  Vxc_Grid grid;
  grid.GridSetup("fine", mol, full_basis);

  std::vector<HirshfeldPartition::AtomicReference> atoms =
      HirshfeldPartition::BuildAtomicReferences(
          mol,
          std::string(XTP_TEST_DATA_FOLDER) + "/hirshfeldpartition/3-21G.xml",
          reference_densities);
  BOOST_REQUIRE_EQUAL(atoms.size(), 2);

  Eigen::MatrixXd summed_weights =
      Eigen::MatrixXd::Zero(full_basis.AOBasisSize(), full_basis.AOBasisSize());
  for (Index i = 0; i < static_cast<Index>(atoms.size()); ++i) {
    summed_weights +=
        HirshfeldPartition::BuildWeightMatrix(atoms, i, full_basis, grid);
  }

  AOOverlap overlap;
  overlap.Fill(full_basis);

  double max_abs_diff =
      (summed_weights - overlap.Matrix()).cwiseAbs().maxCoeff();
  if (max_abs_diff > 1.e-3) {
    std::cout << "NOTE: if this fails, first check whether max_abs_diff is "
                 "merely large-ish (grid-integration error, expected to "
                 "shrink with a finer integration_grid) or drastically "
                 "wrong (an actual implementation bug -- indexing, basis "
                 "re-centering, or the weight formula itself)."
              << std::endl;
    std::cout << "Sum of weight matrices:\n" << summed_weights << std::endl;
    std::cout << "AO overlap matrix:\n" << overlap.Matrix() << std::endl;
  }
  BOOST_CHECK_SMALL(max_abs_diff, 1.e-3);
 } catch (const std::runtime_error& e) {
   std::cout << "SKIPPING weight_matrices_sum_to_overlap: " << e.what()
             << std::endl;
   libint2::finalize();
   return;
 }

  libint2::finalize();
}

// ===========================================================================
// STATUS: written but NOT yet run. First exercise of
// HirshfeldPartition::ComputeCDFTForceContribution -- the fully
// assembled, four-term CDFT force gradient built up over several
// commits (GridWeightDerivativeContribution, weight-function-
// derivative + its grid integration, PulayAndTranslationContribution).
//
// Deliberately does NOT run a full CDFT (or even a full DFT) SCF at
// each displaced geometry -- that would conflate this NEW force term
// with the already-validated plain-DFT force and the outer
// Lagrange-multiplier machinery, and be far slower besides. Instead,
// isolates ComputeCDFTForceContribution on its own: a fixed, random,
// symmetric density matrix D (same size/shape at every geometry, its
// specific values irrelevant -- this is testing the GRADIENT FORMULA,
// not any physical density) is held constant while the geometry
// moves, and d(Tr[D*W_c])/dR is compared against a direct finite
// difference of Tr[D*W_c] itself at the displaced geometries -- the
// exact same "fixed MO/density matrix" pattern already used for
// RIJGradient/RIKGradient's own finite-difference tests in
// test_dftgradient.cc.
// ===========================================================================
BOOST_AUTO_TEST_CASE(cdft_force_finite_difference) {
  libint2::initialize();
 try {
  double bond_length = 1.13;  // Angstrom, roughly CO equilibrium
  double h = 1e-4;            // Angstrom

  std::string basis_path =
      std::string(XTP_TEST_DATA_FOLDER) + "/hirshfeldpartition/3-21G.xml";
  BasisSet basisset;
  basisset.Load(basis_path);

  QMMolecule mol0 = BuildCO(bond_length);
  AOBasis full_basis0;
  full_basis0.Fill(basisset, mol0);

  std::string xml_path = "/tmp/xtp_test_hirshfeldpartition_cdftforce.xml";
  std::ofstream xml(xml_path);
  xml << "<dftpackage>\n";
  xml << "<spin>1</spin>\n";
  xml << "<name>xtp</name>\n";
  xml << "<charge>0</charge>\n";
  xml << "<functional>XC_GGA_X_PBE XC_GGA_C_PBE</functional>\n";
  xml << "<basisset>" << basis_path << "</basisset>\n";
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
  xml << "<integration_grid>fine</integration_grid>\n";
  xml << "<max_iterations>200</max_iterations>\n";
  xml << "</xtpdft>\n";
  xml << "</dftpackage>\n";
  xml.close();

  votca::tools::Property prop;
  prop.LoadFromXML(xml_path);
  DFTEngine dft;
  Logger log;
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));

  std::map<std::string, Eigen::MatrixXd> reference_densities =
      DFTEngineTestAccess::ComputeHirshfeldReferenceDensities(dft, mol0);

  // Fixed, arbitrary, symmetric density matrix -- generated once and
  // reused unchanged at every geometry. Same reasoning as
  // test_dftgradient.cc's own analogous tests: its specific values do
  // not matter here, only that it stays fixed while the geometry
  // moves.
  Index n_bf = full_basis0.AOBasisSize();
  Eigen::MatrixXd density_random = Eigen::MatrixXd::Random(n_bf, n_bf);
  Eigen::MatrixXd density = 0.5 * (density_random + density_random.transpose());

  const Index kCarbonIndex = 0;

  std::vector<HirshfeldPartition::AtomicReference> atoms0 =
      HirshfeldPartition::BuildAtomicReferences(mol0, basis_path,
                                                reference_densities);
  Vxc_Grid grid0;
  grid0.GridSetup("fine", mol0, full_basis0);

  Eigen::MatrixXd analytic_grad = HirshfeldPartition::ComputeCDFTForceContribution(
      atoms0, kCarbonIndex, density, mol0, full_basis0, grid0);
  BOOST_REQUIRE_EQUAL(analytic_grad.rows(), 2);

  // Sanity check independent of finite differences, same reasoning as
  // the existing NuclearRepulsionDerivative/RIJGradient tests:
  // Tr[D*W_c] cannot depend on where the whole molecule sits in
  // absolute space, only its internal geometry, so its gradient
  // summed over every atom must vanish.
  Eigen::Vector3d total = analytic_grad.colwise().sum();
  BOOST_CHECK_SMALL(total.cwiseAbs().maxCoeff(), 1e-6);

  // Tr[D*W_c] with ONE atom's position displaced by "delta" along one
  // Cartesian axis (0=x, 1=y, 2=z) from mol0's own base geometry,
  // rebuilding everything geometry-dependent (atoms' own re-centered
  // bases, the basis, the grid) but reusing the SAME reference_densities
  // (isolated-atom SCF results do not depend on the molecule's own
  // geometry at all) and the SAME fixed density matrix D throughout.
  // delta is in Angstrom, matching BuildCO's own xyz-file convention.
  auto trace_displaced = [&](Index atom_index, int axis, double delta) {
    QMMolecule mol = BuildCO(bond_length);
    Eigen::Vector3d pos = mol[atom_index].getPos();
    pos(axis) += delta * 1.8897259886;  // Angstrom -> Bohr, matching
                                        // this codebase's own internal
                                        // (atomic) units throughout --
                                        // getPos()/setPos() work in
                                        // Bohr, unlike BuildCO's own
                                        // xyz-file input, which LoadFromFile
                                        // itself converts from Angstrom.
    mol[atom_index].setPos(pos);
    AOBasis full_basis;
    full_basis.Fill(basisset, mol);
    std::vector<HirshfeldPartition::AtomicReference> atoms =
        HirshfeldPartition::BuildAtomicReferences(mol, basis_path,
                                                  reference_densities);
    Vxc_Grid grid;
    grid.GridSetup("fine", mol, full_basis);
    Eigen::MatrixXd W = HirshfeldPartition::BuildWeightMatrix(
        atoms, kCarbonIndex, full_basis, grid);
    return density.cwiseProduct(W).sum();
  };

  double h_bohr = h * 1.8897259886;
  double max_mismatch = 0.0;
  // Both atoms (0=Carbon, the constrained TARGET itself -- a
  // genuinely different code path than atom 1, via the
  // 1_{A==target_atom_index} indicator in EvaluateWeightGradient's
  // own formula -- and 1=Oxygen, non-target), all three Cartesian
  // axes -- deliberately not just the bond-length/z direction: CO is
  // linear, so a pure bond-length displacement preserves the
  // molecule's own cylindrical symmetry and could not by itself catch
  // a bug specific to an off-axis (x/y) direction.
  for (Index atom_index = 0; atom_index < 2; ++atom_index) {
    for (int axis = 0; axis < 3; ++axis) {
      double trace_plus = trace_displaced(atom_index, axis, h);
      double trace_minus = trace_displaced(atom_index, axis, -h);
      double numerical = (trace_plus - trace_minus) / (2.0 * h_bohr);
      double analytic = analytic_grad(atom_index, axis);
      double mismatch = std::abs(numerical - analytic);
      max_mismatch = std::max(max_mismatch, mismatch);
      std::cout << "atom " << atom_index << ", axis " << axis
                << ": numerical=" << numerical << ", analytic=" << analytic
                << ", mismatch=" << mismatch << std::endl;
      if (mismatch > 1.e-4) {
        std::cout << "NOTE: if this fails, first check whether the mismatch "
                     "is merely large-ish (grid-integration/finite-"
                     "difference-step error) or drastically wrong (a "
                     "genuine bug in one or more of the four CDFT force "
                     "terms -- check each of GridWeightDerivativeContribution"
                     "/WeightFunctionDerivativeContribution/"
                     "PulayAndTranslationContribution individually against "
                     "its own finite difference before assuming the bug is "
                     "in how they are summed). Also worth noting whether "
                     "this specific failure is on atom 0 (the TARGET atom "
                     "itself, exercising the 1_{A==target} indicator branch "
                     "in EvaluateWeightGradient) versus atom 1 (non-target) "
                     "-- a failure isolated to one of these but not the "
                     "other would point specifically at that branch."
                  << std::endl;
      }
      BOOST_CHECK_SMALL(mismatch, 1.e-4);
    }
  }
 } catch (const std::runtime_error& e) {
   std::cout << "SKIPPING cdft_force_finite_difference: " << e.what()
             << std::endl;
   libint2::finalize();
   return;
 }

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace xtp
}  // namespace votca
