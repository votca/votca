/*
 *            Copyright 2009-2026 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_PODCOUPLING_H
#define VOTCA_XTP_PODCOUPLING_H

#include "logger.h"
#include "votca/xtp/orbitals.h"
#include <votca/tools/types.h>
#include <votca/xtp/aobasis.h>

namespace votca {
namespace xtp {

// Maps a fragment's own atom indices to the full, concatenated list of
// AO basis-function indices belonging to those atoms, in the SAME
// order the atom indices were given (not sorted) -- confirmed
// necessary directly, with the user, before implementing anything
// else here: the actual, real use case (matching this session's own
// CDFT tests) has fragment atom indices that are a genuinely
// DISJOINT, scattered set, not a contiguous, block-wise range (e.g.
// "atoms 1-5" for one fragment, "atoms 6-11" for the other). Each
// INDIVIDUAL atom's own AOs ARE always contiguous (a direct
// consequence of how AO basis sets are constructed, atom by atom, in
// order -- confirmed via AOBasis::getFuncPerAtom(), not assumed), but
// a FRAGMENT's own AOs, across scattered atoms, generally are not --
// this function does the (otherwise easy to get subtly wrong)
// gathering explicitly, once, in one place, rather than leaving every
// caller to reimplement it.
//
// Deliberately takes func_per_atom directly (AOBasis::getFuncPerAtom()'s
// own return value) rather than a full AOBasis -- makes this function
// trivially unit-testable with a small, synthetic vector, with no need
// to construct any real molecule/basis set at all.
std::vector<Index> MapAtomsToAOIndices(const std::vector<Index>& atom_indices,
                                       const std::vector<Index>& func_per_atom);

// Projection Operator Diabatization, POD2 variant (Ghan et al., avoiding
// the original POD's own global Lowdin-orthogonalization step, which
// was found to make results basis-set-unstable -- see
// arXiv:1512.00200 "Critical analysis of fragment-orbital DFT schemes
// for the calculation of electronic coupling values" for the direct
// comparison against the original POD, FODFT, and CDFT that motivated
// choosing this specific variant).
//
// Unlike CDFT (a self-consistent, constrained calculation on the
// intact supermolecule) or FODFT (two separate, ISOLATED fragment
// calculations, which cannot represent a covalent bond crossing the
// fragment boundary at all), POD2 starts from a SINGLE, already-
// converged, ORDINARY (unconstrained, neutral, ground-state) DFT
// calculation on the intact, covalently-bonded supermolecule -- no
// bond-cutting, no excited-state calculation, no CDFT constraint of
// any kind. The converged Fock matrix's own donor/acceptor AO
// sub-blocks are separately diagonalized (in the ORIGINAL AO basis,
// not a globally Lowdin-orthogonalized one -- this is the specific
// "2" in POD2), giving fragment-localized "diabatic" orbitals
// directly; the coupling is then read off from the off-diagonal
// donor-acceptor block of the same Fock matrix, transformed into this
// new, fragment-block-diagonal basis.
class PODCoupling {
 public:
  PODCoupling(Orbitals& orbitals, Logger* log,
              std::vector<Index> fragment_A_atoms,
              std::vector<Index> fragment_B_atoms);

  // Computes the FULL, pairwise coupling matrix between a symmetric
  // range of fragment A's own orbitals and a symmetric range of
  // fragment B's own orbitals, covering BOTH occupied (hole-transport)
  // and virtual (electron-transport) orbitals in a single call --
  // matching DFTcoupling's own, established convention exactly
  // (numberofstatesA/B, DetermineRangeOfStates, getCouplingElement),
  // per direct agreement with the user rather than inventing a
  // different one for this class.
  //
  // numberofstatesA/B = N means fragment A/B's own range covers N
  // orbitals below its own HOMO (inclusive) through N orbitals above
  // its own LUMO (inclusive) -- i.e. {HOMO-N+1, ..., HOMO, LUMO, ...,
  // LUMO+N-1}, exactly matching
  // DFTcoupling::DetermineRangeOfStates's own definition (minimal =
  // HOMO - numberofstates + 1, maximal = LUMO + numberofstates - 1).
  // N=1 gives just {HOMO, LUMO}.
  void CalculateCouplings(Index numberofstatesA, Index numberofstatesB);

  // Coupling [Hartree] between fragment A's orbital at ABSOLUTE (not
  // fragment-local-range-relative) fragment-orbital index levelA and
  // fragment B's orbital at absolute index levelB -- e.g. levelA =
  // getFragmentAHomoIndex() for fragment A's own HOMO. Must be within
  // the range covered by the most recent CalculateCouplings call, or
  // this throws -- matching DFTcoupling::getCouplingElement's own
  // convention of taking absolute MO indices (there, absolute within
  // the monomer's own, separate Orbitals; here, absolute within the
  // fragment's own, separate eigendecomposition -- there is no
  // separate, standalone "fragment A calculation" with its own MO
  // numbering to refer to at all, since fragment orbitals come from
  // block-diagonalizing sub-blocks of the SAME, single supermolecule
  // Fock matrix -- see this class's own, leading comment above).
  double getCouplingElement(Index levelA, Index levelB) const;

  // Each fragment's own HOMO/LUMO index, in the SAME, absolute index
  // space getCouplingElement expects for levelA/levelB -- lets a
  // caller request e.g. "fragment A's HOMO-1, fragment B's LUMO+1"
  // without needing to know the fragment's own total AO-basis size or
  // do this arithmetic itself. Available immediately after
  // construction -- unlike Range_orbA/B (only meaningful after
  // CalculateCouplings has actually run), each fragment's own HOMO/
  // LUMO index depends only on its own nuclear charge (see
  // CountFragmentElectrons in podcoupling.cc), not on which orbital
  // range was actually requested/computed.
  Index getFragmentAHomoIndex() const { return nocc_A_ - 1; }
  Index getFragmentALumoIndex() const { return nocc_A_; }
  Index getFragmentBHomoIndex() const { return nocc_B_ - 1; }
  Index getFragmentBLumoIndex() const { return nocc_B_; }

  // Returns fragment A's (if fragment_A is true, else fragment B's)
  // orbital at absolute index level -- same index convention as
  // getCouplingElement's own levelA/levelB above -- EMBEDDED into the
  // full, whole-molecule AO basis (i.e. a vector of length equal to
  // the full molecule's own total number of AO basis functions, with
  // zero coefficients everywhere outside this fragment's own AO
  // indices). This is what visualizing a fragment orbital (e.g. via a
  // cube file, which needs a coefficient for every AO of the full
  // molecule the grid is built over) actually requires -- the
  // fragment orbital's own, "native" representation (a much shorter
  // vector, over only this fragment's own AOs) is what getCouplingElement's
  // own internal computation uses, but is not, on its own, something a
  // whole-molecule grid/cube-file writer can consume at all. Must be
  // within the range covered by the most recent CalculateCouplings
  // call, or this throws, matching getCouplingElement's own behavior.
  Eigen::VectorXd GetFragmentOrbital(bool fragment_A, Index level) const;

  // Human-readable breakdown of which AOs actually dominate a given
  // fragment orbital -- the top_n largest-|Mulliken population|
  // AOs, each labeled by which atom (index + element) and shell type
  // (s/p/d/...) it belongs to. Added directly to answer a real,
  // concrete question: large electron density visible on the OTHER
  // fragment's own physical space, in a cube-file rendering of a
  // fragment orbital, does NOT mean a nonzero coefficient on the
  // other fragment's own AOs -- there is never one, by construction
  // (see GetFragmentOrbital's own comment above). It means a large
  // contribution from THIS fragment's own AO, whose own basis
  // function has a spatial tail reaching into the other fragment's
  // space -- this shows directly which AO, on which atom, that
  // actually is.
  //
  // Ranks by Mulliken population (c_i*(S*c)_i), NOT raw |coefficient|
  // -- confirmed directly necessary, from a real, misleading run: in
  // a non-orthogonal AO basis, several strongly-overlapping AOs (e.g.
  // multiple contracted S-shells on nearby atoms) can develop huge,
  // largely-cancelling raw coefficients that dominate a plain
  // |coefficient| ranking while contributing little to the orbital's
  // own, actual normalized shape -- the Mulliken weighting corrects
  // for this by explicitly accounting for overlap between AOs, the
  // same way ordinary Mulliken population analysis does.
  std::string DescribeFragmentOrbitalComposition(bool fragment_A,
                                                  Index level,
                                                  Index top_n = 5) const;

 private:
  Orbitals& orbitals_;
  Logger* pLog_;
  std::vector<Index> fragment_A_atoms_;
  std::vector<Index> fragment_B_atoms_;

  // Each fragment's own, ESTIMATED number of occupied orbitals -- see
  // this class's own .cc file (CountFragmentElectrons) for how this is
  // computed and why it is only ever an estimate for a covalently-
  // bonded fragment specifically. Computed once, in the constructor,
  // since it depends only on each fragment's own atoms/nuclear
  // charge, not on anything CalculateCouplings itself computes.
  Index nocc_A_ = 0;
  Index nocc_B_ = 0;

  // Set by CalculateCouplings; std::pair<start, size>, in each
  // fragment's own, absolute orbital-index space (see
  // getFragmentAHomoIndex/etc.'s own comment above) -- matching
  // DFTcoupling::Range_orbA/Range_orbB's own convention exactly.
  std::pair<Index, Index> Range_orbA_;
  std::pair<Index, Index> Range_orbB_;

  // (Range_orbA_.second, Range_orbB_.second)-sized coupling matrix
  // [Hartree], set by CalculateCouplings -- JAB_(i,j) is the coupling
  // between fragment A's own (Range_orbA_.first + i)-th orbital and
  // fragment B's own (Range_orbB_.first + j)-th orbital. Matches
  // DFTcoupling's own JAB member directly, except sized only to the
  // relevant off-diagonal block (DFTcoupling's own JAB additionally
  // spans the A-A and B-B diagonal blocks in one, larger, combined
  // matrix; POD2's own fragment orbitals are already, separately
  // block-diagonal by construction, so there is no equivalent
  // diagonal-block information to store here at all).
  Eigen::MatrixXd JAB_;

  // Set by CalculateCouplings, needed by GetFragmentOrbital to embed a
  // fragment orbital back into the full, whole-molecule AO basis: each
  // fragment's own eigenvectors (from its own, separate generalized
  // eigenvalue solve -- es_A/es_B in podcoupling.cc), each fragment's
  // own AO indices within the full molecule (the same mapping
  // MapAtomsToAOIndices already computes for CalculateCouplings' own,
  // internal use), and the full molecule's own total AO count (to
  // size the embedded, zero-padded result vector correctly). Cheap to
  // store in full -- each is only ever fragment-sized
  // (n_A x n_A/n_A), not full-molecule-sized.
  Eigen::MatrixXd fragment_A_eigenvectors_;
  Eigen::MatrixXd fragment_B_eigenvectors_;
  std::vector<Index> ao_indices_A_;
  std::vector<Index> ao_indices_B_;
  Index nao_full_ = 0;

  // Each fragment's own AO overlap sub-block (S_AA/S_BB, already
  // computed internally by CalculateCouplings for the generalized
  // eigenvalue solve itself) -- stored here specifically for
  // DescribeFragmentOrbitalComposition's own use: a plain, raw AO
  // coefficient is NOT, on its own, a meaningful measure of that AO's
  // real contribution to a normalized orbital in a non-orthogonal
  // basis (confirmed directly, from a real, misleading run: several
  // strongly-overlapping S-shells developed enormous, largely-
  // cancelling raw coefficients that dominated a plain |coefficient|
  // ranking while contributing little to the orbital's actual shape).
  // The standard, correct fix is a Mulliken-style weighting,
  // c_i*(S*c)_i, which needs this fragment's own S explicitly.
  Eigen::MatrixXd S_AA_;
  Eigen::MatrixXd S_BB_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PODCOUPLING_H
