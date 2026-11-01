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
//
// Tests of PODCoupling, with emphasis on MapAtomsToAOIndices: fragment
// definitions are in general non-contiguous sets of atom indices, and
// index-mapping errors would silently corrupt all results.
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE podcoupling_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/podcoupling.h"
#include "votca/xtp/qmmolecule.h"
#include "votca/xtp/toolfactory.h"
#include "xtp_libint2.h"

#include <cmath>
#include <fstream>
#include <sstream>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(podcoupling_test)

namespace {
// Atom lines (Angstrom) of the co-facial ethylene dimer at the given
// separation, in canonical order: molecule A (C1, C2, H1a, H1b, H2a,
// H2b), then molecule B in the same order. Callers that need another
// atom order permute the list and adapt the fragment indices.
std::vector<std::string> EthyleneDimerAtomLines(double separation_angstrom) {
  std::ostringstream sep;
  sep << separation_angstrom;
  return {
      "C -0.6625 0.0000 0.0000",         // A: C1  (canonical index 0)
      "C 0.6625 0.0000 0.0000",          // A: C2  (canonical index 1)
      "H -1.2332 0.9240 0.0000",         // A: H1a (canonical index 2)
      "H -1.2332 -0.9240 0.0000",        // A: H1b (canonical index 3)
      "H 1.2332 0.9240 0.0000",          // A: H2a (canonical index 4)
      "H 1.2332 -0.9240 0.0000",         // A: H2b (canonical index 5)
      "C -0.6625 0.0000 " + sep.str(),   // B: C1  (canonical index 6)
      "C 0.6625 0.0000 " + sep.str(),    // B: C2  (canonical index 7)
      "H -1.2332 0.9240 " + sep.str(),   // B: H1a (canonical index 8)
      "H -1.2332 -0.9240 " + sep.str(),  // B: H1b (canonical index 9)
      "H 1.2332 0.9240 " + sep.str(),    // B: H2a (canonical index 10)
      "H 1.2332 -0.9240 " + sep.str(),   // B: H2b (canonical index 11)
  };
}

// Runs a DFT SCF on the atoms in the given order and returns |J| [eV]
// of the HOMO-HOMO POD coupling for the given fragment atom indices,
// which must refer to that order.
//
// If half_homo_gap_ev_out is given, it receives 0.5*|E_HOMO - E_HOMO-1|
// [eV] of the same SCF. For a symmetric dimer this equals the coupling
// exactly (Baumeier, Kirkpatrick, Andrienko, PCCP 12, 11103 (2010),
// Eq. 12). If write_orb_file is given, the orbitals are written there.
double RunEthyleneDimerCoupling(const std::vector<std::string>& atom_lines,
                                const std::vector<Index>& fragment_A_atoms,
                                const std::vector<Index>& fragment_B_atoms,
                                double* half_homo_gap_ev_out = nullptr,
                                const std::string* write_orb_file = nullptr) {
  std::string tmp_path = "/tmp/xtp_test_podcoupling_ethylene_dimer.xyz";
  std::ofstream mol_out(tmp_path);
  mol_out << atom_lines.size() << "\n\n";
  for (const std::string& line : atom_lines) {
    mol_out << line << "\n";
  }
  mol_out.close();

  Orbitals orb;
  orb.QMAtoms().LoadFromFile(tmp_path);

  std::ofstream xml("dftengine_ethylene_dimer.xml");
  xml << "<dftpackage>\n";
  xml << "<spin>1</spin>\n";
  xml << "<name>xtp</name>\n";
  xml << "<charge>0</charge>\n";
  xml << "<functional>XC_HYB_GGA_XC_B3LYP</functional>\n";
  xml << "<basisset>def2-svp</basisset>\n";
  xml << "<initial_guess>independent</initial_guess>\n";
  xml << "<xtpdft>\n";
  xml << "<convergence>\n";
  xml << "    <energy>1e-7</energy>\n";
  xml << "    <method>DIIS</method>\n";
  xml << "    <DIIS_start>0.002</DIIS_start>\n";
  xml << "    <ADIIS_start>0.8</ADIIS_start>\n";
  xml << "    <DIIS_length>20</DIIS_length>\n";
  xml << "    <levelshift>0.0</levelshift>\n";
  xml << "    <levelshift_end>0.2</levelshift_end>\n";
  xml << "    <max_iterations>100</max_iterations>\n";
  xml << "    <error>1e-7</error>\n";
  xml << "    <DIIS_maxout>false</DIIS_maxout>\n";
  xml << "    <mixing>0.7</mixing>\n";
  xml << "    <mixing_max>0.98</mixing_max>\n";
  xml << "    <mixing_end>0.8</mixing_end>\n";
  xml << "    <davidson_max_iter>50</davidson_max_iter>\n";
  xml << "</convergence>\n";
  xml << "<integration_grid>medium</integration_grid>\n";
  xml << "</xtpdft>\n";
  xml << "</dftpackage>\n";
  xml.close();
  votca::tools::Property prop;
  prop.LoadFromXML("dftengine_ethylene_dimer.xml");

  DFTEngine dft;
  Logger log(votca::Log::error);
  dft.setLogger(&log);
  dft.Initialize(prop.get("dftpackage"));
  bool converged = dft.Evaluate(orb);
  BOOST_REQUIRE_EQUAL(converged, true);

  // Reused by podcouplingtool_run_end_to_end to avoid another SCF.
  if (write_orb_file != nullptr) {
    orb.WriteToCpt(*write_orb_file);
  }

  if (half_homo_gap_ev_out != nullptr) {
    // Number of doubly occupied orbitals from the total nuclear charge;
    // exact for the neutral closed-shell dimer (32 electrons, 16
    // orbitals).
    double nuccharge = 0.0;
    const QMMolecule& mol = orb.QMAtoms();
    for (Index i = 0; i < mol.size(); ++i) {
      nuccharge += static_cast<double>(mol[i].getNuccharge());
    }
    Index nocc_dimer = static_cast<Index>(std::lround(nuccharge / 2.0));
    const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
    BOOST_REQUIRE_GT(nocc_dimer, 0);
    BOOST_REQUIRE_LT(nocc_dimer, eps.size());
    double homo = eps(nocc_dimer - 1);
    double homo_minus_1 = eps(nocc_dimer - 2);
    *half_homo_gap_ev_out =
        0.5 * std::abs(homo - homo_minus_1) * 27.211386245988;
  }

  // HOMO-HOMO (hole) coupling; numberofstates = 1 gives {HOMO, LUMO}
  // per fragment, as in DFTcoupling.
  PODCoupling pod(orb, &log, fragment_A_atoms, fragment_B_atoms);
  pod.CalculateCouplings(1, 1);
  double coupling_hartree = pod.getCouplingElement(pod.getFragmentAHomoIndex(),
                                                   pod.getFragmentBHomoIndex());

  // Logger only buffers messages; print them (e.g. the PODCoupling
  // diagnostics) explicitly.
  std::cout << log;

  return std::abs(coupling_hartree) * 27.211386245988;  // Hartree -> eV
}
}  // namespace

BOOST_AUTO_TEST_CASE(map_atoms_to_ao_indices_disjoint_out_of_order) {
  // Unequal AO counts per atom, so that prefix-sum errors are detected.
  // Starting AO indices: 0, 3, 4, 8, 9 (total 14).
  //   atom 0: AOs [0, 1, 2]
  //   atom 1: AOs [3]
  //   atom 2: AOs [4, 5, 6, 7]
  //   atom 3: AOs [8]
  //   atom 4: AOs [9, 10, 11, 12, 13]
  std::vector<Index> func_per_atom = {3, 1, 4, 1, 5};

  // Non-contiguous and unsorted: the given order must be preserved.
  std::vector<Index> atom_indices = {3, 0, 4};

  std::vector<Index> result = MapAtomsToAOIndices(atom_indices, func_per_atom);

  // AOs of atom 3, then atom 0, then atom 4.
  std::vector<Index> expected = {8, 0, 1, 2, 9, 10, 11, 12, 13};

  BOOST_REQUIRE_EQUAL(result.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    BOOST_CHECK_EQUAL(result[i], expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(map_atoms_to_ao_indices_single_atom) {
  // A single atom other than the first gets exactly its own AO range.
  std::vector<Index> func_per_atom = {2, 3, 4};
  std::vector<Index> atom_indices = {1};
  std::vector<Index> result = MapAtomsToAOIndices(atom_indices, func_per_atom);
  std::vector<Index> expected = {2, 3, 4};
  BOOST_REQUIRE_EQUAL(result.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    BOOST_CHECK_EQUAL(result[i], expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(map_atoms_to_ao_indices_all_atoms_contiguous_case) {
  // Contiguous fragment.
  std::vector<Index> func_per_atom = {2, 2, 2, 2};
  std::vector<Index> atom_indices = {1, 2};
  std::vector<Index> result = MapAtomsToAOIndices(atom_indices, func_per_atom);
  std::vector<Index> expected = {2, 3, 4, 5};
  BOOST_REQUIRE_EQUAL(result.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    BOOST_CHECK_EQUAL(result[i], expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(map_atoms_to_ao_indices_out_of_range_throws) {
  std::vector<Index> func_per_atom = {2, 3};
  std::vector<Index> atom_indices = {5};
  BOOST_CHECK_THROW(MapAtomsToAOIndices(atom_indices, func_per_atom),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(pod_coupling_ethylene_dimer_consistency_checks) {
  // Two checks on the co-facial ethylene dimer at d = 4.0 Angstrom,
  // combined to limit the number of SCF runs (two in total):
  // 1. Scrambled vs. contiguous atom order, testing the fragment
  //    handling end to end (geometry, basis, Fock/overlap indexing).
  // 2. Coupling of the contiguous case vs. half the HOMO/HOMO-1
  //    splitting, which agree exactly for a symmetric dimer (Baumeier,
  //    Kirkpatrick, Andrienko, PCCP 12, 11103 (2010), Eq. 12).
  libint2::initialize();

  // Scrambled order with non-contiguous, unsorted fragments. In terms of
  // the canonical indices (A: 0-5, B: 6-11):
  //   new position:      0   1  2  3  4  5  6  7  8  9  10 11
  //   canonical index:  11   0  6  2  8  1  7  3  9  4  10 5
  // so fragment A is at {1,5,3,7,9,11} and fragment B at {2,6,4,8,10,0}.
  std::vector<Index> scramble_order = {11, 0, 6, 2, 8, 1, 7, 3, 9, 4, 10, 5};
  std::vector<Index> fragment_A_scrambled = {1, 5, 3, 7, 9, 11};
  std::vector<Index> fragment_B_scrambled = {2, 6, 4, 8, 10, 0};

  double separation = 4.0;
  std::vector<std::string> canonical_lines = EthyleneDimerAtomLines(separation);
  std::vector<std::string> scrambled_lines(canonical_lines.size());
  for (size_t new_pos = 0; new_pos < scramble_order.size(); ++new_pos) {
    scrambled_lines[new_pos] = canonical_lines[size_t(scramble_order[new_pos])];
  }

  std::vector<Index> fragment_A_contiguous = {0, 1, 2, 3, 4, 5};
  std::vector<Index> fragment_B_contiguous = {6, 7, 8, 9, 10, 11};
  double half_homo_gap = 0.0;
  std::string orb_file_path =
      "/tmp/xtp_test_podcoupling_ethylene_dimer_contiguous.orb";
  double coupling_contiguous = RunEthyleneDimerCoupling(
      canonical_lines, fragment_A_contiguous, fragment_B_contiguous,
      &half_homo_gap, &orb_file_path);
  double coupling_scrambled = RunEthyleneDimerCoupling(
      scrambled_lines, fragment_A_scrambled, fragment_B_scrambled);

  std::cout << "PODCoupling ethylene dimer, scrambled-atom-order check: "
               "contiguous |J| = "
            << coupling_contiguous
            << " eV, scrambled |J| = " << coupling_scrambled << " eV"
            << std::endl;
  std::cout << "PODCoupling vs. half-HOMO-gap check: PODCoupling |J| = "
            << coupling_contiguous
            << " eV, 0.5*(HOMO-HOMO-1) = " << half_homo_gap << " eV"
            << std::endl;

  // Same physical system; only numerical noise (SCF convergence, grid,
  // summation order) is allowed. Indexing errors give large deviations.
  BOOST_CHECK_CLOSE(coupling_scrambled, coupling_contiguous, 1.0);

  // Looser tolerance: equal in theory, but computed by different routes
  // (fragment-block diagonalization vs. dimer MO energies).
  BOOST_CHECK_CLOSE(coupling_contiguous, half_homo_gap, 5.0);

  libint2::finalize();
}

BOOST_AUTO_TEST_CASE(podcouplingtool_run_end_to_end) {
  // Runs PodCouplingTool (ParseOptions/Run) end to end on the .orb file
  // written by pod_coupling_ethylene_dimer_consistency_checks, which must
  // run first (Boost.Test runs test cases in declaration order); a
  // missing file makes Orbitals::ReadFromCpt throw.
  //
  // The tool is created through QMToolFactory, as in xtp_tools, since
  // its header in libxtp/tools/ is not on the test include path.
  std::ofstream xml("podcoupling_tool_test.xml");
  xml << "<podcoupling>\n";
  xml << "<job_name>podcoupling_tool_test</job_name>\n";
  xml << "<orb_file>/tmp/xtp_test_podcoupling_ethylene_dimer_contiguous.orb"
         "</orb_file>\n";
  xml << "<fragment_A>0:5</fragment_A>\n";
  xml << "<fragment_B>6:11</fragment_B>\n";
  xml << "<levA>1</levA>\n";
  xml << "<levB>1</levB>\n";
  xml << "<write_cube_files>false</write_cube_files>\n";
  xml << "</podcoupling>\n";
  xml.close();
  tools::Property prop;
  prop.LoadFromXML("podcoupling_tool_test.xml");

  std::unique_ptr<QMTool> tool = QMToolFactory().Create("podcoupling");
  BOOST_REQUIRE(tool != nullptr);
  tool->Initialize(prop.get("podcoupling"));
  bool success = tool->Evaluate();
  BOOST_CHECK_EQUAL(success, true);
}

BOOST_AUTO_TEST_SUITE_END()
