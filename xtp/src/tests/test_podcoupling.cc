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
// Focused specifically on MapAtomsToAOIndices -- the one, genuinely new
// piece of index-mapping logic PODCoupling needs, and the exact place a
// subtle bug could silently corrupt every downstream result (directly
// flagged as the primary risk before any of this was written: real
// fragment definitions, both for this and for the existing CDFT tests,
// use a genuinely DISJOINT, scattered set of atom indices, not a
// contiguous, block-wise range).
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE podcoupling_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/podcoupling.h"
#include "votca/xtp/qmmolecule.h"
#include "xtp_libint2.h"

#include <cmath>
#include <fstream>
#include <sstream>

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(podcoupling_test)

namespace {
// The 12, individual atom lines for the co-facial ethylene dimer at a
// given separation, in a FIXED, "canonical" order (molecule A's own 6
// atoms first: C1, C2, H1a, H1b, H2a, H2b -- then molecule B's own 6,
// in the same internal order) -- see
// pod_coupling_ethylene_dimer_decays_with_distance's own header
// comment, below, for the geometry itself and why these specific bond
// lengths/angle were chosen. Callers that want a different atom
// ORDER (e.g. the scrambled test case below) permute this list
// themselves, keeping track of where each atom ends up so the
// corresponding fragment_A_atoms/fragment_B_atoms can be constructed
// to match.
std::vector<std::string> EthyleneDimerAtomLines(double separation_angstrom) {
  std::ostringstream sep;
  sep << separation_angstrom;
  return {
      "C -0.6625 0.0000 0.0000",          // A: C1  (canonical index 0)
      "C 0.6625 0.0000 0.0000",           // A: C2  (canonical index 1)
      "H -1.2332 0.9240 0.0000",          // A: H1a (canonical index 2)
      "H -1.2332 -0.9240 0.0000",         // A: H1b (canonical index 3)
      "H 1.2332 0.9240 0.0000",           // A: H2a (canonical index 4)
      "H 1.2332 -0.9240 0.0000",          // A: H2b (canonical index 5)
      "C -0.6625 0.0000 " + sep.str(),    // B: C1  (canonical index 6)
      "C 0.6625 0.0000 " + sep.str(),     // B: C2  (canonical index 7)
      "H -1.2332 0.9240 " + sep.str(),    // B: H1a (canonical index 8)
      "H -1.2332 -0.9240 " + sep.str(),   // B: H1b (canonical index 9)
      "H 1.2332 0.9240 " + sep.str(),     // B: H2a (canonical index 10)
      "H 1.2332 -0.9240 " + sep.str(),    // B: H2b (canonical index 11)
  };
}

// Shared by both the contiguous and scrambled ethylene-dimer test
// cases below -- runs a real DFT SCF on the given, already-ordered
// list of atom lines (in EXACTLY the order they will appear in the
// resulting geometry file -- the caller is responsible for whatever
// permutation it wants to test), then computes PODCoupling's own
// HOMO-HOMO coupling for the given fragment atom indices (which must
// correctly correspond to WHATEVER order atom_lines is actually in --
// this is exactly what lets the scrambled test case reuse this same
// function while still validly cross-checking against the contiguous
// case's own result: same physical system, same physical answer,
// different bookkeeping only).
//
// half_homo_gap_ev_out (optional, nullptr by default -- existing
// callers that only want the POD coupling itself are unaffected): if
// given, filled with 0.5*(E_HOMO - E_HOMO-1) of the SAME, already-
// converged dimer SCF this function already ran for PODCoupling
// itself -- essentially free, no additional calculation needed. Per
// the paper's own eqn (12) and its own direct confirmation ("both
// approaches yield identical transfer integrals for this symmetric
// dimer configuration"), this simpler estimate is not merely an
// approximation but EXACTLY equal to the full projective coupling
// for a symmetric dimer specifically (identical monomers, symmetric
// mutual orientation) -- exactly the case this co-facial ethylene
// dimer test already is, so this is a genuinely strong, independent,
// theoretically-exact cross-check, not just another qualitative
// sanity check.
double RunEthyleneDimerCoupling(const std::vector<std::string>& atom_lines,
                                const std::vector<Index>& fragment_A_atoms,
                                const std::vector<Index>& fragment_B_atoms,
                                double* half_homo_gap_ev_out = nullptr) {
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

  if (half_homo_gap_ev_out != nullptr) {
    // Number of occupied (doubly-filled) dimer orbitals -- computed
    // dynamically from the loaded geometry's own total nuclear
    // charge (matching PODCoupling's own convention for a fragment's
    // occupied-orbital count, applied here to the whole dimer instead
    // of a fragment) rather than hardcoded, so this does not silently
    // go stale if the geometry ever changes. For this specific,
    // neutral, closed-shell 24-atom-free (12-atom) ethylene dimer
    // this is exact (32 electrons total, 16 doubly-occupied orbitals),
    // not an approximation -- unlike the analogous per-FRAGMENT count
    // PODCoupling itself uses internally, which is genuinely
    // approximate for a covalently-bonded fragment (see podcoupling.cc's
    // own comment on this) -- there is no such ambiguity for the
    // WHOLE, intact, neutral dimer here at all.
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

  // HOMO-HOMO coupling (hole transfer) -- matching the paper's own,
  // explicit focus ("Throughout the paper we will consider the
  // particular case of hole transport"). numberofstatesA/B=1 gives
  // just {HOMO, LUMO} for each fragment -- matching
  // DFTcoupling::CalculateCouplings' own convention exactly, per
  // direct agreement with the user, rather than this class's own,
  // earlier, single-orbital-pair-only interface.
  PODCoupling pod(orb, &log, fragment_A_atoms, fragment_B_atoms);
  pod.CalculateCouplings(1, 1);
  double coupling_hartree =
      pod.getCouplingElement(pod.getFragmentAHomoIndex(),
                             pod.getFragmentBHomoIndex());

  // Logger only ever BUFFERS its own messages internally -- nothing
  // prints them anywhere on its own. Confirmed directly, from a real
  // run, that this was exactly why an earlier, added S_AB diagnostic
  // (XTP_LOG(...) inside PODCoupling::CalculateCouplings itself) never
  // actually appeared in this test's own output at all -- the same,
  // underlying mechanism already found earlier this session for
  // pyxtp's own capture_standard_output. Explicitly flushed here via
  // the Logger class's own, public friend operator<<(std::ostream&,
  // Logger&), matching the same pattern already confirmed working
  // there.
  std::cout << log;

  return std::abs(coupling_hartree) * 27.211386245988;  // Hartree -> eV
}
}  // namespace

BOOST_AUTO_TEST_CASE(map_atoms_to_ao_indices_disjoint_out_of_order) {
  // Deliberately varied per-atom function counts (NOT all equal -- an
  // off-by-one or wrong-cumulative-sum bug could easily still pass a
  // uniform-size test by accident) for 5 atoms: 3, 1, 4, 1, 5 AOs each.
  // Prefix sum (each atom's own starting AO index): 0, 3, 4, 8, 9, [14].
  //   atom 0: AOs [0, 1, 2]
  //   atom 1: AOs [3]
  //   atom 2: AOs [4, 5, 6, 7]
  //   atom 3: AOs [8]
  //   atom 4: AOs [9, 10, 11, 12, 13]
  std::vector<Index> func_per_atom = {3, 1, 4, 1, 5};

  // Deliberately DISJOINT (skips atoms 1 and 2 entirely) and DELIBERATELY
  // OUT OF NUMERICAL ORDER (3, 0, 4 -- not 0, 3, 4) -- confirms the
  // function preserves the caller's own given order rather than
  // silently, implicitly sorting the atom indices first.
  std::vector<Index> atom_indices = {3, 0, 4};

  std::vector<Index> result = MapAtomsToAOIndices(atom_indices, func_per_atom);

  // Hand-computed expected result, in the SAME order as atom_indices:
  // atom 3's own AOs ([8]), then atom 0's own AOs ([0,1,2]), then atom
  // 4's own AOs ([9,10,11,12,13]).
  std::vector<Index> expected = {8, 0, 1, 2, 9, 10, 11, 12, 13};

  BOOST_REQUIRE_EQUAL(result.size(), expected.size());
  for (size_t i = 0; i < expected.size(); ++i) {
    BOOST_CHECK_EQUAL(result[i], expected[i]);
  }
}

BOOST_AUTO_TEST_CASE(map_atoms_to_ao_indices_single_atom) {
  // Simplest possible case, checked directly: a single atom, not the
  // first one in the molecule, gets exactly its own AO range and
  // nothing else.
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
  // A contiguous, block-wise fragment (the "easy" case this whole
  // function exists to NOT be limited to) should still work correctly
  // -- confirms the general, disjoint-capable implementation does not
  // accidentally break the simpler, contiguous case it needs to
  // subsume.
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
  // Combines two, independent consistency checks on the SAME co-facial
  // ethylene dimer (contiguous atom order, d=4.0 Angstrom) -- merged
  // into one test case, per direct user request to cut this file's
  // own CI runtime (each check needs its own real DFT SCF, and this
  // file previously ran 5 across three separate test cases; this
  // reduces that to 2 total, without losing either check's own
  // coverage): the "distance decays" test case this file used to have
  // was dropped entirely (a genuine, real coverage loss, accepted
  // directly by the user); this merge, in contrast, costs nothing, since
  // the two checks below were already running on the exact same,
  // contiguous d=4.0 geometry as two, separate, fully redundant SCF
  // calculations.
  //
  // 1. Direct, end-to-end validation of the disjoint-fragment handling
  //    (scrambled vs. contiguous atom order) -- exercising the full,
  //    real pipeline (geometry loading, basis construction,
  //    Fock/overlap matrix indexing, AND the fragment gather itself),
  //    not just MapAtomsToAOIndices in isolation (already covered
  //    separately, by map_atoms_to_ao_indices_disjoint_out_of_order at
  //    the top of this file).
  // 2. Cross-check of the CONTIGUOUS case's own result against the
  //    dimer's own half-HOMO-HOMO-1 gap -- exact for a symmetric dimer
  //    per the paper this whole test system is drawn from (Baumeier,
  //    Kirkpatrick, Andrienko, PCCP 2010, 12, 11103, eqn (12): "both
  //    approaches yield identical transfer integrals for this
  //    symmetric dimer configuration"), and free to obtain here since
  //    the contiguous case's own dimer SCF is already being run for
  //    check 1 anyway.
  libint2::initialize();

  // Deliberately scrambled: neither molecule's own 6 atoms are
  // contiguous in this new order, and the resulting fragment indices
  // are not sorted either -- computed precisely (not by hand) as the
  // permutation's own inverse, to avoid exactly the kind of
  // off-by-one bug this test exists to guard against. In terms of
  // the CANONICAL indices from EthyleneDimerAtomLines (A: 0-5, B:
  // 6-11), this new file's own atom order is:
  //   new position:      0   1  2  3  4  5  6  7  8  9  10 11
  //   canonical index:  11   0  6  2  8  1  7  3  9  4  10 5
  // i.e. fragment A's own 6 atoms end up at new positions
  // {1,5,3,7,9,11} and fragment B's own 6 end up at {2,6,4,8,10,0}.
  std::vector<Index> scramble_order = {11, 0, 6, 2, 8, 1, 7, 3, 9, 4, 10, 5};
  std::vector<Index> fragment_A_scrambled = {1, 5, 3, 7, 9, 11};
  std::vector<Index> fragment_B_scrambled = {2, 6, 4, 8, 10, 0};

  double separation = 4.0;
  std::vector<std::string> canonical_lines =
      EthyleneDimerAtomLines(separation);
  std::vector<std::string> scrambled_lines(canonical_lines.size());
  for (size_t new_pos = 0; new_pos < scramble_order.size(); ++new_pos) {
    scrambled_lines[new_pos] =
        canonical_lines[size_t(scramble_order[new_pos])];
  }

  std::vector<Index> fragment_A_contiguous = {0, 1, 2, 3, 4, 5};
  std::vector<Index> fragment_B_contiguous = {6, 7, 8, 9, 10, 11};
  double half_homo_gap = 0.0;
  double coupling_contiguous = RunEthyleneDimerCoupling(
      canonical_lines, fragment_A_contiguous, fragment_B_contiguous,
      &half_homo_gap);
  double coupling_scrambled = RunEthyleneDimerCoupling(
      scrambled_lines, fragment_A_scrambled, fragment_B_scrambled);

  std::cout << "PODCoupling ethylene dimer, scrambled-atom-order check: "
               "contiguous |J| = "
            << coupling_contiguous << " eV, scrambled |J| = "
            << coupling_scrambled << " eV" << std::endl;
  std::cout << "PODCoupling vs. half-HOMO-gap check: PODCoupling |J| = "
            << coupling_contiguous << " eV, 0.5*(HOMO-HOMO-1) = "
            << half_homo_gap << " eV" << std::endl;

  // These describe the exact same physical system (same geometry,
  // same separation, same two fragments -- only the order atoms
  // happen to be listed in, and correspondingly the atom-index
  // bookkeeping, differ) -- a correct implementation must give the
  // same coupling either way, up to ordinary numerical noise (SCF
  // convergence threshold, integration grid, floating-point
  // summation order differences from the different loop iteration
  // order this reordering induces). A real indexing bug in the
  // disjoint-fragment handling would instead show up here as a
  // LARGE, qualitative discrepancy, not a small numerical one.
  BOOST_CHECK_CLOSE(coupling_scrambled, coupling_contiguous, 1.0);

  // A somewhat looser tolerance than the scrambled-vs-contiguous check
  // above (which compares two calls that are, numerically, almost
  // identical up to floating-point summation order) -- these two
  // quantities are computed via genuinely DIFFERENT routes (one via
  // fragment-block Fock diagonalization and an off-diagonal matrix
  // element, the other via the dimer's own, already-diagonalized MO
  // energies directly), so full floating-point-level agreement is not
  // expected even though the paper's own eqn (12) shows them to be
  // mathematically identical in the symmetric-dimer limit -- still
  // tight enough that a genuine, substantial error in either PODCoupling
  // itself or this test's own half-gap calculation would be caught.
  BOOST_CHECK_CLOSE(coupling_contiguous, half_homo_gap, 5.0);

  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
