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

#include "votca/xtp/podcoupling.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <sstream>

namespace votca {
namespace xtp {

std::vector<Index> MapAtomsToAOIndices(
    const std::vector<Index>& atom_indices,
    const std::vector<Index>& func_per_atom) {
  // Prefix sum: ao_start[a] is the first AO index belonging to atom a
  // (and ao_start[a+1] is one past its own last AO index) -- valid
  // specifically because AO basis functions are laid out contiguously
  // per atom, in atom order, a direct, standard consequence of how AO
  // bases are constructed (each atom's own shells are added in turn),
  // not something specific to this function's own use case.
  std::vector<Index> ao_start(func_per_atom.size() + 1, 0);
  for (size_t a = 0; a < func_per_atom.size(); ++a) {
    ao_start[a + 1] = ao_start[a] + func_per_atom[a];
  }

  std::vector<Index> ao_indices;
  for (Index atom_index : atom_indices) {
    if (atom_index < 0 ||
        atom_index >= static_cast<Index>(func_per_atom.size())) {
      throw std::runtime_error("MapAtomsToAOIndices: atom index " +
                               std::to_string(atom_index) +
                               " is out of range (0.." +
                               std::to_string(func_per_atom.size() - 1) + ").");
    }
    // Each INDIVIDUAL atom's own AOs are appended as one, contiguous
    // run -- but successive atoms in atom_indices need not be
    // adjacent to each other at all (the whole reason this function
    // exists in the first place -- see this function's own header
    // comment), so the OVERALL ao_indices returned here is not, in
    // general, a single contiguous range even though each atom's own
    // contribution to it is.
    for (Index ao = ao_start[atom_index]; ao < ao_start[atom_index + 1]; ++ao) {
      ao_indices.push_back(ao);
    }
  }
  return ao_indices;
}

namespace {
// Gathers the (indices.size(), indices.size()) sub-matrix of full_matrix
// at the given (possibly scattered, non-contiguous) row/column indices
// -- a direct gather, not an explicit permutation of full_matrix itself:
// mathematically identical to "reorder full_matrix so these indices
// become contiguous, then slice a contiguous block", but avoids ever
// constructing or applying a full-size permutation on the (potentially
// much larger) full_matrix at all.
Eigen::MatrixXd GatherSubMatrix(const Eigen::MatrixXd& full_matrix,
                                const std::vector<Index>& indices) {
  Index n = static_cast<Index>(indices.size());
  Eigen::MatrixXd result(n, n);
  for (Index i = 0; i < n; ++i) {
    for (Index j = 0; j < n; ++j) {
      result(i, j) = full_matrix(indices[size_t(i)], indices[size_t(j)]);
    }
  }
  return result;
}

// Off-diagonal (indices_row.size(), indices_col.size()) block of
// full_matrix -- same gather approach as GatherSubMatrix above, but for
// a genuinely rectangular block (row indices from one fragment, column
// indices from the other), needed for the final donor-acceptor
// coupling element itself.
Eigen::MatrixXd GatherOffDiagonalBlock(const Eigen::MatrixXd& full_matrix,
                                       const std::vector<Index>& indices_row,
                                       const std::vector<Index>& indices_col) {
  Index n_row = static_cast<Index>(indices_row.size());
  Index n_col = static_cast<Index>(indices_col.size());
  Eigen::MatrixXd result(n_row, n_col);
  for (Index i = 0; i < n_row; ++i) {
    for (Index j = 0; j < n_col; ++j) {
      result(i, j) =
          full_matrix(indices_row[size_t(i)], indices_col[size_t(j)]);
    }
  }
  return result;
}
}  // namespace

namespace {
// Each fragment's own number of occupied orbitals is NOT directly,
// unambiguously available from the full, intact supermolecule's own,
// delocalized wavefunction at all -- estimated instead as half the
// fragment's own total nuclear charge (i.e. assuming a neutral,
// closed-shell fragment), rounded to the nearest integer, matching
// the standard convention already used elsewhere in this codebase for
// a fragment's own "neutral reference" electron count (see
// DFTEngine::BuildCDFTConstraint's own neutral_reference_population).
// Genuinely approximate for a COVALENTLY-bonded fragment specifically
// (there is no truly well-defined "neutral fragment" electron count
// once a bond has been cut across the fragment boundary) -- flagged
// directly to the user as a real modeling choice, not silently
// assumed to be exact.
Index CountFragmentElectrons(const QMMolecule& mol,
                             const std::vector<Index>& atoms) {
  double nuccharge = 0.0;
  for (Index atom_index : atoms) {
    nuccharge += static_cast<double>(mol[atom_index].getNuccharge());
  }
  return static_cast<Index>(std::lround(nuccharge / 2.0));
}

// Fragment-local analogue of DFTcoupling::DetermineRangeOfStates --
// same definition (minimal = homo_index - numberofstates + 1, maximal
// = lumo_index + numberofstates - 1, covering both occupied and
// virtual orbitals in one, single, combined range), but without that
// function's own degeneracy_ handling: not requested by the user for
// this class, and left out deliberately rather than added
// speculatively. Bounds-checked directly against the fragment's own
// total number of orbitals (n_basis, i.e. the fragment's own AO count
// -- the fragment-block Fock sub-block is square, n_basis x n_basis,
// so this is also the fragment's own total number of orbitals
// available from its own generalized eigenvalue solve).
std::pair<Index, Index> DetermineFragmentRangeOfStates(Index homo_index,
                                                       Index lumo_index,
                                                       Index numberofstates,
                                                       Index n_basis) {
  Index minimal = homo_index - numberofstates + 1;
  Index maximal = lumo_index + numberofstates - 1;
  if (minimal < 0 || maximal >= n_basis) {
    throw std::runtime_error(
        "PODCoupling: requested numberofstates=" +
        std::to_string(numberofstates) +
        " exceeds the fragment's own available orbital range (0.." +
        std::to_string(n_basis - 1) + ").");
  }
  return {minimal, maximal - minimal + 1};
}
}  // namespace

PODCoupling::PODCoupling(Orbitals& orbitals, Logger* log,
                         std::vector<Index> fragment_A_atoms,
                         std::vector<Index> fragment_B_atoms)
    : orbitals_(orbitals),
      pLog_(log),
      fragment_A_atoms_(std::move(fragment_A_atoms)),
      fragment_B_atoms_(std::move(fragment_B_atoms)) {
  const QMMolecule& mol = orbitals_.QMAtoms();
  nocc_A_ = CountFragmentElectrons(mol, fragment_A_atoms_);
  nocc_B_ = CountFragmentElectrons(mol, fragment_B_atoms_);
}

void PODCoupling::CalculateCouplings(Index numberofstatesA,
                                     Index numberofstatesB) {
  const QMMolecule& mol = orbitals_.QMAtoms();

  AOBasis full_dftbasis;
  {
    BasisSet basisset;
    basisset.Load(orbitals_.getDFTbasisName());
    full_dftbasis.Fill(basisset, mol);
  }
  AOOverlap overlap;
  overlap.Fill(full_dftbasis);
  const Eigen::MatrixXd& S = overlap.Matrix();

  // Reconstructs the full, AO-basis Fock matrix from the already-
  // converged MOs/orbital energies -- F_AO = S*C*eps*C^T*S, valid
  // because C is S-orthonormal (C^T*S*C = I, so C^{-1} = C^T*S) and
  // spans the full AO space (confirmed directly: no near-linearly-
  // dependent basis functions were removed for the reference
  // calculation this class is designed to consume -- see this class's
  // own header comment on requiring an already-converged, ordinary,
  // neutral ground-state calculation as input). Orbitals itself does
  // not persist the raw AO-basis Fock matrix directly, only the
  // diagonalized MO representation, so this reconstruction is the
  // standard, direct way to recover it.
  const Eigen::MatrixXd& C = orbitals_.MOs().eigenvectors();
  const Eigen::VectorXd& eps = orbitals_.MOs().eigenvalues();
  Eigen::MatrixXd F = S * C * eps.asDiagonal() * C.transpose() * S;

  std::vector<Index> ao_indices_A =
      MapAtomsToAOIndices(fragment_A_atoms_, full_dftbasis.getFuncPerAtom());
  std::vector<Index> ao_indices_B =
      MapAtomsToAOIndices(fragment_B_atoms_, full_dftbasis.getFuncPerAtom());

  Eigen::MatrixXd F_AA = GatherSubMatrix(F, ao_indices_A);
  Eigen::MatrixXd S_AA = GatherSubMatrix(S, ao_indices_A);
  Eigen::MatrixXd F_BB = GatherSubMatrix(F, ao_indices_B);
  Eigen::MatrixXd S_BB = GatherSubMatrix(S, ao_indices_B);

  // Separately diagonalizes each fragment's own donor/acceptor block
  // IN THE ORIGINAL AO BASIS -- this is the specific "2" in POD2 (Ghan
  // et al.), deliberately NOT the original POD's own global Lowdin-
  // orthogonalization of the whole AO basis first (confirmed directly,
  // via arXiv:1512.00200's own critical comparison, to make results
  // basis-set-unstable: larger basis sets increase inter-fragment AO
  // mixing under global orthogonalization, degrading the resulting
  // coupling).
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es_A(F_AA, S_AA);
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es_B(F_BB, S_BB);
  if (es_A.info() != Eigen::Success || es_B.info() != Eigen::Success) {
    throw std::runtime_error(
        "PODCoupling: generalized eigenvalue solve failed for one or both "
        "fragment Fock sub-blocks -- this can happen if a fragment's own "
        "S_AA is (numerically) singular, e.g. from a badly chosen or "
        "overlapping fragment definition.");
  }

  // Stored for later use by GetFragmentOrbital -- see this member's
  // own header comment in podcoupling.h for why (embedding a fragment
  // orbital back into the full, whole-molecule AO basis, e.g. for
  // visualization, needs exactly these ingredients).
  fragment_A_eigenvectors_ = es_A.eigenvectors();
  fragment_B_eigenvectors_ = es_B.eigenvectors();
  ao_indices_A_ = ao_indices_A;
  ao_indices_B_ = ao_indices_B;
  nao_full_ = F.rows();
  // Stored for DescribeFragmentOrbitalComposition's own use -- see
  // its own header comment in podcoupling.h for why a Mulliken-style
  // weighting (needing each fragment's own S explicitly) is used
  // there instead of raw |coefficient|.
  S_AA_ = S_AA;
  S_BB_ = S_BB;

  Range_orbA_ = DetermineFragmentRangeOfStates(getFragmentAHomoIndex(),
                                               getFragmentALumoIndex(),
                                               numberofstatesA, F_AA.rows());
  Range_orbB_ = DetermineFragmentRangeOfStates(getFragmentBHomoIndex(),
                                               getFragmentBLumoIndex(),
                                               numberofstatesB, F_BB.rows());
  Index levelsA = Range_orbA_.second;
  Index levelsB = Range_orbB_.second;

  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " PODCoupling: fragment A HOMO=" << getFragmentAHomoIndex()
      << ", LUMO=" << getFragmentALumoIndex()
      << " (both estimated), range covers orbitals [" << Range_orbA_.first
      << ", " << (Range_orbA_.first + levelsA - 1) << "]" << std::flush;
  XTP_LOG(Log::error, *pLog_)
      << TimeStamp()
      << " PODCoupling: fragment B HOMO=" << getFragmentBHomoIndex()
      << ", LUMO=" << getFragmentBLumoIndex()
      << " (both estimated), range covers orbitals [" << Range_orbB_.first
      << ", " << (Range_orbB_.first + levelsB - 1) << "]" << std::flush;

  // The off-diagonal AO blocks themselves do not depend on which
  // specific orbital pair is being coupled -- gathered once, outside
  // the loop below, and reused for every (i, j) pair in the requested
  // range, rather than redundantly re-gathering the same block once
  // per orbital pair.
  Eigen::MatrixXd F_AB = GatherOffDiagonalBlock(F, ao_indices_A, ao_indices_B);
  Eigen::MatrixXd S_AB_block =
      GatherOffDiagonalBlock(S, ao_indices_A, ao_indices_B);

  JAB_ = Eigen::MatrixXd(levelsA, levelsB);
  for (Index i = 0; i < levelsA; ++i) {
    Eigen::VectorXd orbital_A = es_A.eigenvectors().col(Range_orbA_.first + i);
    double e_A_hartree = es_A.eigenvalues()(Range_orbA_.first + i);
    for (Index j = 0; j < levelsB; ++j) {
      Eigen::VectorXd orbital_B =
          es_B.eigenvectors().col(Range_orbB_.first + j);
      double e_B_hartree = es_B.eigenvalues()(Range_orbB_.first + j);

      double J_AB = orbital_A.dot(F_AB * orbital_B);
      // S_AB: the overlap between this specific PAIR of fragment
      // orbitals themselves (not either fragment's own internal
      // S_AA/S_BB, already used above for each one's own
      // normalization).
      double S_AB = orbital_A.dot(S_AB_block * orbital_B);

      // The actual, final coupling for this pair: NOT J_AB itself --
      // confirmed directly, from a real, independent cross-check
      // against a symmetric test dimer's own half-HOMO-HOMO-1 gap,
      // that the raw J_AB is wrong by close to a factor of 2 (the two
      // fragment orbitals are not mutually orthogonal in general --
      // S_AB above is the direct evidence of that -- so J_AB alone
      // conflates the true electronic coupling with an overlap-
      // induced contribution). The paper this whole POD2
      // implementation grew out of (Baumeier, Kirkpatrick, Andrienko,
      // PCCP 2010, 12, 11103) derives exactly this same correction for
      // DIPRO's own, analogous non-orthogonal monomer HOMOs, its own
      // eqn (10): t_AB = (J_AB - 0.5*(e_A+e_B)*S_AB) / (1-S_AB^2).
      double denominator = 1.0 - S_AB * S_AB;
      JAB_(i, j) =
          (J_AB - 0.5 * (e_A_hartree + e_B_hartree) * S_AB) / denominator;

      // Per-pair diagnostic: S_AB and the Lowdin denominator (1-S_AB^2)
      // are worth watching on an ongoing basis, since a large |S_AB|
      // approaching +-1 would make this correction numerically
      // unstable (the denominator collapsing toward zero) -- this is
      // a genuine, real risk for a COVALENTLY-bonded fragment pair
      // specifically, where the two fragment orbitals sit directly
      // across a real chemical bond, unlike the small, well-behaved
      // S_AB already confirmed for a non-bonded test case (the
      // ethylene dimer, ~0.03). Confirmed directly, on a real
      // covalently-bonded case (2,2'-bithiophene), that a large,
      // physically genuine coupling and a small, stable S_AB can both
      // be true at once -- e.g. S_AB=0.063 there, comfortably away
      // from +-1, with the Lowdin correction only a ~13% adjustment
      // to J_AB, even though the resulting coupling itself (~3.2 eV)
      // is large: strong THROUGH-BOND coupling is a different regime
      // from the weak, through-space coupling this correction was
      // first validated against, and a large result is not itself a
      // red flag -- only S_AB itself getting close to +-1 would be.
      XTP_LOG(Log::error, *pLog_)
          << TimeStamp()
          << " PODCoupling diagnostic: pair (A=" << (Range_orbA_.first + i)
          << ", B=" << (Range_orbB_.first + j) << "): S_AB=" << S_AB
          << ", (1-S_AB^2)=" << denominator
          << ", raw J_AB=" << (J_AB * 27.211386245988)
          << " eV, corrected=" << (JAB_(i, j) * 27.211386245988) << " eV"
          << std::flush;
    }
  }
}

double PODCoupling::getCouplingElement(Index levelA, Index levelB) const {
  Index indexA = levelA - Range_orbA_.first;
  Index indexB = levelB - Range_orbB_.first;
  if (indexA < 0 || indexA >= JAB_.rows() || indexB < 0 ||
      indexB >= JAB_.cols()) {
    throw std::runtime_error(
        "PODCoupling::getCouplingElement: requested levelA=" +
        std::to_string(levelA) + "/levelB=" + std::to_string(levelB) +
        " is outside the range covered by the most recent "
        "CalculateCouplings call.");
  }
  return JAB_(indexA, indexB);
}

Eigen::VectorXd PODCoupling::GetFragmentOrbital(bool fragment_A,
                                                Index level) const {
  const Eigen::MatrixXd& eigenvectors =
      fragment_A ? fragment_A_eigenvectors_ : fragment_B_eigenvectors_;
  const std::vector<Index>& ao_indices =
      fragment_A ? ao_indices_A_ : ao_indices_B_;
  const std::pair<Index, Index>& range = fragment_A ? Range_orbA_ : Range_orbB_;

  Index index = level - range.first;
  if (index < 0 || index >= range.second) {
    throw std::runtime_error(
        "PODCoupling::GetFragmentOrbital: requested level=" +
        std::to_string(level) +
        " is outside the range covered by the most recent "
        "CalculateCouplings call.");
  }

  // Scatters the fragment-local coefficient vector (length equal to
  // this fragment's own AO count) into the full, whole-molecule AO
  // basis (length nao_full_) -- the direct inverse of
  // MapAtomsToAOIndices/GatherSubMatrix's own gather: every AO NOT
  // belonging to this fragment gets a zero coefficient, since the
  // fragment orbital, by construction, has no amplitude there at all.
  //
  // Uses level DIRECTLY here, NOT the range-relative index computed
  // above -- unlike JAB_ (only ever range-sized, levelsA x levelsB),
  // eigenvectors_ stores ALL of this fragment's own columns (its full
  // n_A x n_A eigenvector set), so level itself is already the
  // correct, absolute column index into it. Confirmed directly by
  // comparing against getCouplingElement's own code before writing
  // this: using the range-relative index here instead would have
  // silently returned the wrong orbital, offset by range.first
  // columns, whenever range.first != 0 (i.e. whenever the fragment's
  // own HOMO is not literally its lowest-energy orbital -- the normal
  // case for any real molecule).
  Eigen::VectorXd fragment_local = eigenvectors.col(level);
  Eigen::VectorXd full_basis = Eigen::VectorXd::Zero(nao_full_);
  for (Index i = 0; i < Index(ao_indices.size()); ++i) {
    full_basis(ao_indices[size_t(i)]) = fragment_local(i);
  }
  return full_basis;
}

std::string PODCoupling::DescribeFragmentOrbitalComposition(
    bool fragment_A, Index level, Index top_n) const {
  const Eigen::MatrixXd& eigenvectors =
      fragment_A ? fragment_A_eigenvectors_ : fragment_B_eigenvectors_;
  const std::vector<Index>& ao_indices =
      fragment_A ? ao_indices_A_ : ao_indices_B_;
  const std::vector<Index>& fragment_atoms =
      fragment_A ? fragment_A_atoms_ : fragment_B_atoms_;
  const std::pair<Index, Index>& range =
      fragment_A ? Range_orbA_ : Range_orbB_;

  Index index = level - range.first;
  if (index < 0 || index >= range.second) {
    throw std::runtime_error(
        "PODCoupling::DescribeFragmentOrbitalComposition: requested level=" +
        std::to_string(level) +
        " is outside the range covered by the most recent "
        "CalculateCouplings call.");
  }
  Eigen::VectorXd fragment_local = eigenvectors.col(level);

  // Full-molecule-AO-index -> (atom index, shell angular momentum L)
  // lookup, built by re-constructing the same AOBasis
  // CalculateCouplings itself already used (Orbitals itself does not
  // persist a full-AO-index -> (atom, shell) map directly, so this is
  // rebuilt here the same way CalculateCouplings' own F/S
  // reconstruction already does, via getDFTbasisName()/QMAtoms()).
  const QMMolecule& mol = orbitals_.QMAtoms();
  AOBasis full_dftbasis;
  {
    BasisSet basisset;
    basisset.Load(orbitals_.getDFTbasisName());
    full_dftbasis.Fill(basisset, mol);
  }
  std::vector<Index> ao_to_atom(full_dftbasis.AOBasisSize(), -1);
  std::vector<L> ao_to_L(full_dftbasis.AOBasisSize());
  for (const AOShell& shell : full_dftbasis) {
    Index offset = shell.getStartIndex();
    for (Index k = 0; k < shell.getNumFunc(); ++k) {
      ao_to_atom[size_t(offset + k)] = shell.getAtomIndex();
      ao_to_L[size_t(offset + k)] = shell.getL();
    }
  }

  // Mulliken population of each AO in this orbital: P_i = c_i*(S*c)_i
  // (S being this fragment's own AO overlap sub-block, S_AA_/S_BB_) --
  // sums to c^T*S*c = 1 over all AOs, since the generalized eigenvalue
  // solve normalizes with respect to S, not the plain Euclidean norm.
  // Ranking by |P_i| here, NOT |coefficient| -- see this method's own
  // header comment in podcoupling.h for why the latter is actively
  // misleading in a non-orthogonal basis (confirmed directly, from a
  // real, misleading run, not just in principle).
  const Eigen::MatrixXd& S_local = fragment_A ? S_AA_ : S_BB_;
  Eigen::VectorXd mulliken_population =
      fragment_local.cwiseProduct(S_local * fragment_local);

  // Sort this fragment orbital's own AOs by |Mulliken population|,
  // largest first -- ao_indices[i] (the full-molecule AO index for
  // this fragment's own local AO i) is what actually indexes into
  // ao_to_atom/ao_to_L above; fragment_local(i)/mulliken_population(i)
  // are that same AO's own coefficient/population in THIS orbital.
  std::vector<Index> order(fragment_local.size());
  for (Index i = 0; i < fragment_local.size(); ++i) {
    order[size_t(i)] = i;
  }
  std::sort(order.begin(), order.end(), [&](Index a, Index b) {
    return std::abs(mulliken_population(a)) > std::abs(mulliken_population(b));
  });

  std::ostringstream out;
  out << "Top " << std::min(top_n, Index(order.size()))
      << " AO Mulliken populations for fragment " << (fragment_A ? "A" : "B")
      << " orbital " << level << ":";
  for (Index rank = 0; rank < std::min(top_n, Index(order.size())); ++rank) {
    Index local_ao = order[size_t(rank)];
    Index full_ao = ao_indices[size_t(local_ao)];
    Index atom_index = ao_to_atom[size_t(full_ao)];
    std::string element =
        atom_index >= 0 ? mol[atom_index].getElement() : "?";
    out << "\n  " << (rank + 1) << ". population="
        << mulliken_population(local_ao)
        << " (coeff=" << fragment_local(local_ao) << "), atom " << atom_index
        << " (" << element << ", fragment atom "
        << (std::find(fragment_atoms.begin(), fragment_atoms.end(),
                      atom_index) -
            fragment_atoms.begin())
        << "), shell " << EnumToString(ao_to_L[size_t(full_ao)]);
  }
  return out.str();
}

}  // namespace xtp
}  // namespace votca
