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
#include "votca/xtp/citations.h"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <sstream>

namespace votca {
namespace xtp {

std::vector<Index> MapAtomsToAOIndices(
    const std::vector<Index>& atom_indices,
    const std::vector<Index>& func_per_atom) {
  // Prefix sum: AOs of atom a are [ao_start[a], ao_start[a+1]). Valid
  // because AO bases are built atom by atom, in atom order.
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
    // Contiguous per atom, but not overall for scattered atom indices.
    for (Index ao = ao_start[atom_index]; ao < ao_start[atom_index + 1]; ++ao) {
      ao_indices.push_back(ao);
    }
  }
  return ao_indices;
}

namespace {
// Square sub-matrix of full_matrix at the given (possibly non-contiguous)
// row/column indices, gathered directly without permuting full_matrix.
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

// Rectangular block of full_matrix with rows from one fragment and
// columns from the other (donor-acceptor block).
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
// Number of occupied orbitals of a fragment, estimated as half its
// nuclear charge rounded to the nearest integer (neutral closed-shell
// fragment), as for the neutral reference population in
// DFTEngine::BuildCDFTConstraint. Only an estimate for covalently bonded
// fragments, whose electron count is not well defined.
Index CountFragmentElectrons(const QMMolecule& mol,
                             const std::vector<Index>& atoms) {
  double nuccharge = 0.0;
  for (Index atom_index : atoms) {
    nuccharge += static_cast<double>(mol[atom_index].getNuccharge());
  }
  return static_cast<Index>(std::lround(nuccharge / 2.0));
}

// Fragment-local analogue of DFTcoupling::DetermineRangeOfStates
// (minimal = homo_index - numberofstates + 1, maximal = lumo_index +
// numberofstates - 1), without degeneracy handling. n_basis is the
// fragment's AO count, i.e. its number of orbitals.
std::pair<Index, Index> DetermineFragmentRangeOfStates(Index homo_index,
                                                       Index lumo_index,
                                                       Index numberofstates,
                                                       Index n_basis) {
  Index minimal = homo_index - numberofstates + 1;
  Index maximal = lumo_index + numberofstates - 1;
  if (minimal < 0 || maximal >= n_basis) {
    throw std::runtime_error("PODCoupling: requested numberofstates=" +
                             std::to_string(numberofstates) +
                             " exceeds the fragment's orbital range (0.." +
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
  citations::Cite("Ghan:2020");
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

  // AO Fock matrix from the converged MOs, F_AO = S*C*eps*C^T*S (Orbitals
  // does not store F_AO). Exact if C is S-orthonormal (C^T*S*C = I) and
  // square, i.e. no linearly dependent basis functions were removed.
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

  // Diagonalize each fragment block in the original, non-orthogonal AO
  // basis (generalized problem F c = e S c). This distinguishes POD2 from
  // POD, which Lowdin-orthogonalizes the whole AO basis first; that mixes
  // AOs across fragments increasingly with basis size (arXiv:1512.00200).
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es_A(F_AA, S_AA);
  Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> es_B(F_BB, S_BB);
  if (es_A.info() != Eigen::Success || es_B.info() != Eigen::Success) {
    throw std::runtime_error(
        "PODCoupling: generalized eigenvalue solve failed for one or both "
        "fragment Fock sub-blocks; the fragment overlap matrix may be "
        "numerically singular (check the fragment definition).");
  }

  // Needed by GetFragmentOrbital and DescribeFragmentOrbitalComposition.
  fragment_A_eigenvectors_ = es_A.eigenvectors();
  fragment_B_eigenvectors_ = es_B.eigenvectors();
  ao_indices_A_ = ao_indices_A;
  ao_indices_B_ = ao_indices_B;
  nao_full_ = F.rows();
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

  // Off-diagonal AO blocks, independent of the orbital pair.
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
      // Overlap of the two fragment orbitals.
      double S_AB = orbital_A.dot(S_AB_block * orbital_B);

      // The fragment orbitals are not orthogonal, so J_AB contains an
      // overlap contribution (off by nearly a factor of 2 compared with
      // half the splitting of a symmetric dimer). Lowdin-type correction
      // as in DIPRO (Baumeier, Kirkpatrick, Andrienko, PCCP 12, 11103
      // (2010), Eq. 10): t_AB = (J_AB - 0.5*(e_A+e_B)*S_AB) / (1-S_AB^2).
      double denominator = 1.0 - S_AB * S_AB;
      JAB_(i, j) =
          (J_AB - 0.5 * (e_A_hartree + e_B_hartree) * S_AB) / denominator;

      // Diagnostic: the correction becomes unstable for |S_AB| -> 1,
      // which is more likely for covalently bonded fragments. A large
      // coupling alone is not a problem: for 2,2'-bithiophene S_AB = 0.063
      // and t_AB is about 3.2 eV (ethylene dimer: S_AB about 0.03).
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

  // Scatter the fragment-local coefficients into the full AO basis
  // (zero on AOs outside the fragment). The eigenvector matrix holds all
  // fragment orbitals, so it is indexed with the absolute level, not
  // with the range-relative index used for JAB_.
  Eigen::VectorXd fragment_local = eigenvectors.col(level);
  Eigen::VectorXd full_basis = Eigen::VectorXd::Zero(nao_full_);
  for (Index i = 0; i < Index(ao_indices.size()); ++i) {
    full_basis(ao_indices[size_t(i)]) = fragment_local(i);
  }
  return full_basis;
}

std::string PODCoupling::DescribeFragmentOrbitalComposition(bool fragment_A,
                                                            Index level,
                                                            Index top_n) const {
  const Eigen::MatrixXd& eigenvectors =
      fragment_A ? fragment_A_eigenvectors_ : fragment_B_eigenvectors_;
  const std::vector<Index>& ao_indices =
      fragment_A ? ao_indices_A_ : ao_indices_B_;
  const std::vector<Index>& fragment_atoms =
      fragment_A ? fragment_A_atoms_ : fragment_B_atoms_;
  const std::pair<Index, Index>& range = fragment_A ? Range_orbA_ : Range_orbB_;

  Index index = level - range.first;
  if (index < 0 || index >= range.second) {
    throw std::runtime_error(
        "PODCoupling::DescribeFragmentOrbitalComposition: requested level=" +
        std::to_string(level) +
        " is outside the range covered by the most recent "
        "CalculateCouplings call.");
  }
  Eigen::VectorXd fragment_local = eigenvectors.col(level);

  // AO index -> (atom index, shell L) lookup; the AOBasis is rebuilt as
  // in CalculateCouplings, since Orbitals does not store it.
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

  // Mulliken population P_i = c_i*(S*c)_i with the fragment overlap
  // block S; sums to c^T*S*c = 1 (S-normalized eigenvectors). Ranked by
  // |P_i| rather than |c_i| (see podcoupling.h).
  const Eigen::MatrixXd& S_local = fragment_A ? S_AA_ : S_BB_;
  Eigen::VectorXd mulliken_population =
      fragment_local.cwiseProduct(S_local * fragment_local);

  // Sort local AOs by |P_i|, largest first; ao_indices maps a local AO to
  // its full-molecule index for ao_to_atom/ao_to_L.
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
    std::string element = atom_index >= 0 ? mol[atom_index].getElement() : "?";
    out << "\n  " << (rank + 1)
        << ". population=" << mulliken_population(local_ao)
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
