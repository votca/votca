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

// Local VOTCA includes
#include "votca/xtp/hirshfeldpartition.h"
#include "votca/xtp/aoshell.h"
#include "votca/xtp/basisset.h"

// Standard includes
#include <stdexcept>

namespace votca {
namespace xtp {

std::vector<HirshfeldPartition::AtomicReference>
HirshfeldPartition::BuildAtomicReferences(
    const QMMolecule& mol, const std::string& basisset_name,
    const std::map<std::string, Eigen::MatrixXd>& reference_densities) {
  BasisSet basisset;
  basisset.Load(basisset_name);

  std::vector<AtomicReference> atoms;
  atoms.reserve(mol.size());
  for (const QMAtom& real_atom : mol) {
    auto it = reference_densities.find(real_atom.getElement());
    if (it == reference_densities.end()) {
      // Deliberately fail loudly rather than silently skip: a missing
      // entry here means DFTEngine::ComputeHirshfeldReferenceDensities
      // was not actually run for every element this molecule contains,
      // which would otherwise show up only much later, as a silently
      // wrong (missing-atom) weight denominator.
      throw std::runtime_error(
          "HirshfeldPartition::BuildAtomicReferences: no reference "
          "density found for element '" +
          real_atom.getElement() +
          "' -- was DFTEngine::ComputeHirshfeldReferenceDensities "
          "actually run for every element in this molecule?");
    }

    // Same construction RunAtomicDFT_unrestricted itself uses (a
    // single-atom QMMolecule, then AOBasis::Fill on it) -- except at
    // this REAL atom's own position, not the origin, so the resulting
    // basis functions are already centered where they need to be for
    // EvaluateAtomicDensity/EvaluateWeight to evaluate correctly.
    QMMolecule single_atom("hirshfeld_reference_atom", 0);
    single_atom.push_back(
        QMAtom(0, real_atom.getElement(), real_atom.getPos()));

    AtomicReference ref;
    ref.basis.Fill(basisset, single_atom);
    ref.density = it->second;
    atoms.push_back(std::move(ref));
  }
  return atoms;
}

double HirshfeldPartition::EvaluateAtomicDensity(
    const AOBasis& atom_basis, const Eigen::MatrixXd& reference_density,
    const Eigen::Vector3d& point) {
  Eigen::VectorXd ao_values = Eigen::VectorXd::Zero(atom_basis.AOBasisSize());
  for (const AOShell& shell : atom_basis) {
    AOShell::AOValues vals = shell.EvalAOspace(point);
    ao_values.segment(shell.getStartIndex(), shell.getNumFunc()) =
        vals.values;
  }
  return ao_values.dot(reference_density * ao_values);
}

double HirshfeldPartition::EvaluateWeight(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const Eigen::Vector3d& point) {
  double denominator = 0.0;
  double numerator = 0.0;
  for (Index j = 0; j < static_cast<Index>(atoms.size()); ++j) {
    double rho_j =
        EvaluateAtomicDensity(atoms[j].basis, atoms[j].density, point);
    denominator += rho_j;
    if (j == target_atom_index) {
      numerator = rho_j;
    }
  }
  // Same negligible-denominator guard pattern already used for the SSW
  // grid weights in GridWeightGradient (kNegligibleWOwner there) --
  // points far from every atom, where every rho_j(point) is
  // negligible, would otherwise divide a near-zero numerator by a
  // near-zero denominator.
  constexpr double kNegligibleDenominator = 1.e-12;
  if (denominator < kNegligibleDenominator) {
    return 0.0;
  }
  return numerator / denominator;
}

}  // namespace xtp
}  // namespace votca
