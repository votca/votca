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

// Returns the AO indices of the given atoms, concatenated in the order
// the atoms are given (not sorted). Fragment atoms need not be a
// contiguous range; the AOs of each single atom are contiguous, those
// of a fragment in general are not.
//
// func_per_atom is AOBasis::getFuncPerAtom(); taking it instead of an
// AOBasis keeps the function testable without a basis set.
std::vector<Index> MapAtomsToAOIndices(const std::vector<Index>& atom_indices,
                                       const std::vector<Index>& func_per_atom);

// Projection Operator Diabatization, POD2 variant (Ghan et al.), which
// avoids the global Lowdin orthogonalization of the original POD that
// makes results basis-set dependent (comparison of POD, FODFT and CDFT:
// arXiv:1512.00200, "Critical analysis of fragment-orbital DFT schemes
// for the calculation of electronic coupling values").
//
// Input is a single converged, unconstrained ground-state DFT calculation
// of the intact supermolecule; unlike FODFT, covalent bonds across the
// fragment boundary are therefore allowed. The donor and acceptor AO
// blocks of the Fock matrix are diagonalized separately in the original
// AO basis, giving fragment-localized diabatic orbitals; the coupling is
// the off-diagonal Fock block transformed into these orbitals.
class PODCoupling {
 public:
  PODCoupling(Orbitals& orbitals, Logger* log,
              std::vector<Index> fragment_A_atoms,
              std::vector<Index> fragment_B_atoms);

  // Computes the couplings between all pairs of orbitals in the ranges
  // of fragments A and B, occupied and virtual, with the conventions of
  // DFTcoupling. numberofstates = N selects {HOMO-N+1, ..., HOMO, LUMO,
  // ..., LUMO+N-1} (as DFTcoupling::DetermineRangeOfStates); N=1 gives
  // {HOMO, LUMO}.
  void CalculateCouplings(Index numberofstatesA, Index numberofstatesB);

  // Coupling [Hartree] between orbitals levelA of fragment A and levelB
  // of fragment B. Levels are absolute indices into the fragment's
  // eigendecomposition (e.g. getFragmentAHomoIndex()), not relative to
  // the computed range, as in DFTcoupling::getCouplingElement. Throws if
  // outside the range of the last CalculateCouplings call.
  double getCouplingElement(Index levelA, Index levelB) const;

  // HOMO/LUMO index of each fragment in the index space of
  // getCouplingElement. Available after construction, since they depend
  // only on the fragment's nuclear charge (CountFragmentElectrons in
  // podcoupling.cc).
  Index getFragmentAHomoIndex() const { return nocc_A_ - 1; }
  Index getFragmentALumoIndex() const { return nocc_A_; }
  Index getFragmentBHomoIndex() const { return nocc_B_ - 1; }
  Index getFragmentBLumoIndex() const { return nocc_B_; }

  // Orbital level (absolute index, as in getCouplingElement) of fragment
  // A (fragment_A true) or B, embedded in the full-molecule AO basis with
  // zero coefficients outside the fragment, e.g. for cube files. Throws
  // if outside the range of the last CalculateCouplings call.
  Eigen::VectorXd GetFragmentOrbital(bool fragment_A, Index level) const;

  // Lists the top_n AOs with the largest |Mulliken population| in a
  // fragment orbital, with atom index, element and shell type. Useful to
  // identify which AO of the fragment causes density in the region of the
  // other fragment (diffuse tails; coefficients on the other fragment's
  // AOs are zero by construction).
  //
  // Ranked by Mulliken population c_i*(S*c)_i rather than |c_i|: in a
  // non-orthogonal basis, strongly overlapping AOs (e.g. several
  // contracted s shells on nearby atoms) can have large, mutually
  // cancelling coefficients that contribute little to the orbital.
  std::string DescribeFragmentOrbitalComposition(bool fragment_A, Index level,
                                                 Index top_n = 5) const;

 private:
  Orbitals& orbitals_;
  Logger* pLog_;
  std::vector<Index> fragment_A_atoms_;
  std::vector<Index> fragment_B_atoms_;

  // Estimated number of occupied orbitals per fragment, set in the
  // constructor (see CountFragmentElectrons in podcoupling.cc).
  Index nocc_A_ = 0;
  Index nocc_B_ = 0;

  // (start, size) of the computed orbital range, in absolute fragment
  // orbital indices; set by CalculateCouplings. As in DFTcoupling.
  std::pair<Index, Index> Range_orbA_;
  std::pair<Index, Index> Range_orbB_;

  // Coupling matrix [Hartree], set by CalculateCouplings: JAB_(i,j)
  // couples orbital Range_orbA_.first + i of A with Range_orbB_.first + j
  // of B. Unlike DFTcoupling's JAB, only the A-B block is stored.
  Eigen::MatrixXd JAB_;

  // Set by CalculateCouplings for GetFragmentOrbital: fragment
  // eigenvectors (all orbitals), fragment AO indices in the full basis,
  // and the full AO count.
  Eigen::MatrixXd fragment_A_eigenvectors_;
  Eigen::MatrixXd fragment_B_eigenvectors_;
  std::vector<Index> ao_indices_A_;
  std::vector<Index> ao_indices_B_;
  Index nao_full_ = 0;

  // Fragment AO overlap blocks, kept for the Mulliken populations in
  // DescribeFragmentOrbitalComposition.
  Eigen::MatrixXd S_AA_;
  Eigen::MatrixXd S_BB_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PODCOUPLING_H
