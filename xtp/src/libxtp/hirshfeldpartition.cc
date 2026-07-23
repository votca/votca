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
#include "votca/xtp/gridbox.h"

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

Eigen::Vector3d HirshfeldPartition::EvaluateAtomicDensityGradient(
    const AOBasis& atom_basis, const Eigen::MatrixXd& reference_density,
    const Eigen::Vector3d& point) {
  Index n = atom_basis.AOBasisSize();
  Eigen::VectorXd ao_values = Eigen::VectorXd::Zero(n);
  Eigen::MatrixX3d ao_derivatives = Eigen::MatrixX3d::Zero(n, 3);
  for (const AOShell& shell : atom_basis) {
    AOShell::AOValues vals = shell.EvalAOspace(point);
    ao_values.segment(shell.getStartIndex(), shell.getNumFunc()) =
        vals.values;
    ao_derivatives.block(shell.getStartIndex(), 0, shell.getNumFunc(), 3) =
        vals.derivatives;
  }
  // rho(r) = phi^T P phi, so grad_rho(r) = 2 * derivatives^T * (P * phi) --
  // same quadratic-form pattern as EvaluateAtomicDensity itself, just
  // carrying the extra factor of 2 and the derivative (rather than
  // value) of the second phi in the product rule.
  return 2.0 * ao_derivatives.transpose() * (reference_density * ao_values);
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

Eigen::Vector3d HirshfeldPartition::EvaluateWeightGradient(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    Index differentiate_atom_index, const Eigen::Vector3d& point) {
  double denominator = 0.0;
  for (Index j = 0; j < static_cast<Index>(atoms.size()); ++j) {
    denominator +=
        EvaluateAtomicDensity(atoms[j].basis, atoms[j].density, point);
  }
  // Same negligible-denominator guard as EvaluateWeight itself, for
  // the identical reason -- points far from every atom, where the
  // whole formula is 0/0 in the limit.
  constexpr double kNegligibleDenominator = 1.e-12;
  if (denominator < kNegligibleDenominator) {
    return Eigen::Vector3d::Zero();
  }

  double w_target = EvaluateWeight(atoms, target_atom_index, point);
  Eigen::Vector3d grad_rho_A = EvaluateAtomicDensityGradient(
      atoms[differentiate_atom_index].basis,
      atoms[differentiate_atom_index].density, point);
  double indicator =
      (differentiate_atom_index == target_atom_index) ? 1.0 : 0.0;
  // d w_target/d R_A = [w_target(point) - 1_{A==target}] *
  //                    grad_rho_A(point) / rho_tot(point)
  // -- see this function's own header comment for the derivation.
  return (w_target - indicator) * grad_rho_A / denominator;
}

Eigen::MatrixXd HirshfeldPartition::BuildWeightMatrix(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const AOBasis& full_dftbasis, const Vxc_Grid& grid) {
  Index full_size = full_dftbasis.AOBasisSize();
  Eigen::MatrixXd W = Eigen::MatrixXd::Zero(full_size, full_size);

  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    Eigen::MatrixXd box_contribution =
        Eigen::MatrixXd::Zero(box.Matrixsize(), box.Matrixsize());

    for (Index p = 0; p < static_cast<Index>(points.size()); ++p) {
      double w_i = EvaluateWeight(atoms, target_atom_index, points[p]);
      if (w_i == 0.0) {
        // Purely a performance guard (skip the AO evaluation entirely
        // when it cannot contribute anything), not a correctness one --
        // see this function's own header comment.
        continue;
      }
      AOShell::AOValues ao = box.CalcAOValues(points[p]);
      box_contribution +=
          (weights[p] * w_i) * (ao.values * ao.values.transpose());
    }

    box.AddtoBigMatrix(W, box_contribution);
  }

  // The outer product ao.values * ao.values.transpose() is already
  // exactly symmetric by construction, so this is a numerical no-op --
  // guards only against tiny floating-point asymmetry accumulating
  // across many grid points, matching the same defensive pattern used
  // elsewhere in this branch for operator matrices built this way.
  W = 0.5 * (W + W.transpose());
  return W;
}

}  // namespace xtp
}  // namespace votca
