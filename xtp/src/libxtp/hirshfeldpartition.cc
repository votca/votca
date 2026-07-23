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
#include <cmath>
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

Eigen::Vector3d HirshfeldPartition::EvaluatePointWeightGradient(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const Eigen::Vector3d& point) {
  double denominator = 0.0;
  Eigen::Vector3d grad_denominator = Eigen::Vector3d::Zero();
  for (Index j = 0; j < static_cast<Index>(atoms.size()); ++j) {
    denominator +=
        EvaluateAtomicDensity(atoms[j].basis, atoms[j].density, point);
    grad_denominator +=
        EvaluateAtomicDensityGradient(atoms[j].basis, atoms[j].density, point);
  }
  constexpr double kNegligibleDenominator = 1.e-12;
  if (denominator < kNegligibleDenominator) {
    return Eigen::Vector3d::Zero();
  }

  double w_target = EvaluateWeight(atoms, target_atom_index, point);
  Eigen::Vector3d grad_rho_target = EvaluateAtomicDensityGradient(
      atoms[target_atom_index].basis, atoms[target_atom_index].density,
      point);
  // grad_r w(r) = [grad_r rho_target(r) - w(r)*grad_r rho_tot(r)] /
  //               rho_tot(r) -- see this function's own header comment.
  return (grad_rho_target - w_target * grad_denominator) / denominator;
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

namespace {
// Faithfully duplicated from vxc_potential.cc's own anonymous namespace
// (SSWValue/SSWDerivative/kSSWCutoff there) -- those are genuinely
// inaccessible from this file (confirmed: they are inside an
// anonymous namespace in that translation unit, not exposed via any
// header), and per the explicit decision this grew out of, this is a
// deliberate duplication rather than a refactor of that existing,
// already-validated file. Must stay byte-for-byte consistent with the
// original if it is ever changed there.
constexpr double kSSWAlpha = 1.0 / 0.30;
constexpr double kSSWCutoff = 0.725;
constexpr double kSqrtPi = 1.7724538509055160273;

double SSWValue(double mu) {
  double val = 0.5 * std::erfc(std::abs(mu / (1.0 - mu * mu)) * kSSWAlpha);
  if (mu > 0.0) {
    val = 1.0 - val;
  }
  return val;
}

double SSWDerivative(double mu) {
  double h = std::abs(mu) / (1.0 - mu * mu);
  double sign_mu = (mu > 0.0) ? 1.0 : ((mu < 0.0) ? -1.0 : 0.0);
  double one_minus_mu2 = 1.0 - mu * mu;
  double hprime = sign_mu * (1.0 + mu * mu) / (one_minus_mu2 * one_minus_mu2);
  double d_erf1c = -(kSSWAlpha / kSqrtPi) *
                   std::exp(-(kSSWAlpha * h) * (kSSWAlpha * h)) * hprime;
  return (mu > 0.0) ? -d_erf1c : d_erf1c;
}
}  // namespace

Eigen::MatrixXd HirshfeldPartition::GridWeightDerivativeContribution(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const Eigen::MatrixXd& density_matrix, const QMMolecule& mol,
    const Vxc_Grid& grid) {
  Index natoms = mol.size();
  Eigen::MatrixXd Rij = grid.CalcInverseAtomDist(mol);
  Eigen::MatrixXd force_contribution = Eigen::MatrixXd::Zero(natoms, 3);

  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }

    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

    for (Index pidx = 0; pidx < box.size(); ++pidx) {
      AOShell::AOValues ao = box.CalcAOValues(points[pidx]);
      double rho_molecule = ao.values.dot(DMAT_here * ao.values);
      double weight = weights[pidx];
      if (std::abs(rho_molecule * weight) < 1.e-20) {
        continue;
      }

      Index owner = owner_atoms[pidx];
      if (owner < 0) {
        throw std::runtime_error(
            "GridWeightDerivativeContribution: grid point has no "
            "owner_atom set -- was this grid built via GridSetup after "
            "the owner-atom tracking change?");
      }

      const Eigen::Vector3d& point = points[pidx];
      double w_c = EvaluateWeight(atoms, target_atom_index, point);

      Eigen::VectorXd rq(natoms);
      for (Index k = 0; k < natoms; ++k) {
        rq(k) = (point - mol[k].getPos()).norm();
      }

      Eigen::VectorXd p = Eigen::VectorXd::Ones(natoms);
      Eigen::MatrixXd mu_table = Eigen::MatrixXd::Zero(natoms, natoms);
      Eigen::MatrixXd sk_table = Eigen::MatrixXd::Zero(natoms, natoms);
      Eigen::MatrixXi hard = Eigen::MatrixXi::Zero(natoms, natoms);
      for (Index ii = 1; ii < natoms; ++ii) {
        for (Index jj = 0; jj < ii; ++jj) {
          double mu = (rq(ii) - rq(jj)) * Rij(jj, ii);
          mu_table(jj, ii) = mu;
          if (mu > kSSWCutoff) {
            p(ii) = 0.0;
            hard(jj, ii) = 1;
          } else if (mu < -kSSWCutoff) {
            p(jj) = 0.0;
            hard(jj, ii) = -1;
          } else {
            double sk = SSWValue(mu);
            sk_table(jj, ii) = sk;
            p(jj) *= sk;
            p(ii) *= (1.0 - sk);
          }
        }
      }
      double wsum = p.sum();
      double w_owner = p(owner) / wsum;
      constexpr double kNegligibleWOwner = 1.e-8;
      if (w_owner < kNegligibleWOwner) {
        continue;
      }
      double C_p = weight / w_owner;

      auto d_rq_dR = [&](Index k, Index A) -> Eigen::Vector3d {
        if (A == owner && A == k) {
          return Eigen::Vector3d::Zero();
        } else if (A == owner) {
          return (point - mol[k].getPos()) / rq(k);
        } else if (A == k) {
          return -(point - mol[k].getPos()) / rq(k);
        }
        return Eigen::Vector3d::Zero();
      };
      auto d_Rab_dR = [&](Index a, Index b, Index A) -> Eigen::Vector3d {
        Eigen::Vector3d rvec = mol[a].getPos() - mol[b].getPos();
        double Rab = rvec.norm();
        if (A == a) {
          return rvec / Rab;
        } else if (A == b) {
          return -rvec / Rab;
        }
        return Eigen::Vector3d::Zero();
      };
      auto dmu_dR = [&](Index a, Index b, Index A) -> Eigen::Vector3d {
        Eigen::Vector3d d_rq_b = d_rq_dR(b, A);
        Eigen::Vector3d d_rq_a = d_rq_dR(a, A);
        double Rab = 1.0 / Rij(a, b);
        Eigen::Vector3d dRab = d_Rab_dR(a, b, A);
        double mu = mu_table(a, b);
        return (d_rq_b - d_rq_a) / Rab - (mu / Rab) * dRab;
      };
      auto dp_dR = [&](Index k, Index A) -> Eigen::Vector3d {
        constexpr double kNegligibleP = 1.e-8;
        if (p(k) < kNegligibleP) {
          return Eigen::Vector3d::Zero();
        }
        Eigen::Vector3d total = Eigen::Vector3d::Zero();
        for (Index b = k + 1; b < natoms; ++b) {
          if (hard(k, b) != 0) {
            continue;
          }
          double skv = sk_table(k, b);
          if (skv < kNegligibleP) {
            continue;
          }
          total += (SSWDerivative(mu_table(k, b)) / skv) * dmu_dR(k, b, A);
        }
        for (Index a = 0; a < k; ++a) {
          if (hard(a, k) != 0) {
            continue;
          }
          double skv = sk_table(a, k);
          double one_minus_skv = 1.0 - skv;
          if (one_minus_skv < kNegligibleP) {
            continue;
          }
          total += (-SSWDerivative(mu_table(a, k)) / one_minus_skv) *
                   dmu_dR(a, k, A);
        }
        return p(k) * total;
      };

      // Same prefactor role as GridWeightGradient's own C_p*rho*xc.f_xc,
      // just with w_c(point) in place of xc.f_xc -- this point's
      // contribution to Tr[D*W_c] is weight*w_c*rho_molecule =
      // C_p*w_owner*w_c*rho_molecule, and differentiating w_owner alone
      // (dw below) needs this same missing C_p factor multiplied back
      // in, for the identical reason documented in GridWeightGradient's
      // own comment on this.
      double prefactor = C_p * w_c * rho_molecule;

      for (Index A = 0; A < natoms; ++A) {
        Eigen::Vector3d dp_owner = dp_dR(owner, A);
        Eigen::Vector3d dwsum = Eigen::Vector3d::Zero();
        for (Index k = 0; k < natoms; ++k) {
          dwsum += dp_dR(k, A);
        }
        Eigen::Vector3d dw = dp_owner / wsum - w_owner * dwsum / wsum;
        force_contribution.row(A) += (prefactor * dw).transpose();
      }
    }
  }
  return force_contribution;
}

Eigen::MatrixXd HirshfeldPartition::PulayAndTranslationContribution(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const Eigen::MatrixXd& density_matrix, const AOBasis& full_dftbasis,
    const Vxc_Grid& grid) {
  Index natoms = static_cast<Index>(full_dftbasis.getFuncPerAtom().size());
  Eigen::MatrixXd force_contribution = Eigen::MatrixXd::Zero(natoms, 3);

  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }

    // DMAT_here carries the same factor-of-2 convention as
    // PulayGradient's own DMAT_here -- temp(mu) below is therefore
    // already 2*(D*phi)_mu, matching that function's own contribution
    // formula exactly (no separate factor of 2 needed at the point of
    // use).
    const Eigen::MatrixXd DMAT_here = 2 * box.ReadFromBigMatrix(density_matrix);

    // Same box-local-AO-index -> atom-index bookkeeping as
    // PulayGradient's own, identical reasoning: a box's significant
    // shells are not guaranteed to be grouped contiguously by atom.
    std::vector<Index> local_idx_to_atom(box.Matrixsize());
    const std::vector<const AOShell*>& shells = box.getShells();
    const std::vector<GridboxRange>& ao_ranges = box.getAOranges();
    for (size_t s = 0; s < shells.size(); ++s) {
      Index atom = shells[s]->getAtomIndex();
      for (Index k = 0; k < ao_ranges[s].size; ++k) {
        local_idx_to_atom[ao_ranges[s].start + k] = atom;
      }
    }

    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();
    const std::vector<Index>& owner_atoms = box.getOwnerAtoms();

    for (Index p = 0; p < box.size(); ++p) {
      const Eigen::Vector3d& point = points[p];
      AOShell::AOValues ao = box.CalcAOValues(point);
      Eigen::VectorXd temp = ao.values.transpose() * DMAT_here;
      double rho_molecule = 0.5 * temp.dot(ao.values);
      double weight = weights[p];
      if (std::abs(rho_molecule * weight) < 1.e-20) {
        continue;
      }

      double w_c = EvaluateWeight(atoms, target_atom_index, point);

      // --- Pulay (basis-function) term ---
      // Same sign/accumulation convention as PulayGradient's own LDA
      // term, confirmed directly: contribution = -weight * <potential>
      // * temp(mu) * ao.derivatives.row(mu), with w_c(point) here in
      // place of xc.df_drho there.
      for (Index mu = 0; mu < box.Matrixsize(); ++mu) {
        Index atom = local_idx_to_atom[mu];
        Eigen::Vector3d contribution =
            -weight * w_c * temp(mu) * ao.derivatives.row(mu).transpose();
        force_contribution.row(atom) += contribution.transpose();
      }

      // --- Grid-point-translation term ---
      // Only the point's OWNER atom is affected (r_p moves rigidly
      // with its owner only) -- see EvaluatePointWeightGradient's own
      // header comment for the derivation. Full product rule, since
      // both w_c(r) and rho_molecule(r) vary here (unlike the analogous
      // XC term, which has only one point-varying factor -- xc.df_drho
      // is evaluated AT rho, not itself differentiated in that specific
      // term).
      Index owner = owner_atoms[p];
      if (owner < 0) {
        throw std::runtime_error(
            "PulayAndTranslationContribution: grid point has no "
            "owner_atom set -- was this grid built via GridSetup after "
            "the owner-atom tracking change?");
      }
      Eigen::Vector3d rho_molecule_grad = temp.transpose() * ao.derivatives;
      Eigen::Vector3d grad_w_c =
          EvaluatePointWeightGradient(atoms, target_atom_index, point);
      Eigen::Vector3d translation_term =
          weight * (grad_w_c * rho_molecule + w_c * rho_molecule_grad);
      force_contribution.row(owner) += translation_term.transpose();
    }
  }
  return force_contribution;
}

Eigen::MatrixXd HirshfeldPartition::WeightFunctionDerivativeContribution(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const Eigen::MatrixXd& density_matrix, const AOBasis& full_dftbasis,
    const Vxc_Grid& grid) {
  Index natoms = static_cast<Index>(full_dftbasis.getFuncPerAtom().size());
  Eigen::MatrixXd force_contribution = Eigen::MatrixXd::Zero(natoms, 3);

  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    const GridBox& box = grid[i];
    if (!box.Matrixsize()) {
      continue;
    }
    const Eigen::MatrixXd DMAT_here = box.ReadFromBigMatrix(density_matrix);
    const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
    const std::vector<double>& weights = box.getGridWeights();

    for (Index p = 0; p < box.size(); ++p) {
      const Eigen::Vector3d& point = points[p];
      AOShell::AOValues ao = box.CalcAOValues(point);
      double rho_molecule = ao.values.dot(DMAT_here * ao.values);
      double weight = weights[p];
      if (std::abs(rho_molecule * weight) < 1.e-20) {
        continue;
      }
      double prefactor = weight * rho_molecule;
      for (Index A = 0; A < natoms; ++A) {
        Eigen::Vector3d dw_dR_A =
            EvaluateWeightGradient(atoms, target_atom_index, A, point);
        force_contribution.row(A) += (prefactor * dw_dR_A).transpose();
      }
    }
  }
  return force_contribution;
}

Eigen::MatrixXd HirshfeldPartition::ComputeCDFTForceContribution(
    const std::vector<AtomicReference>& atoms, Index target_atom_index,
    const Eigen::MatrixXd& density_matrix, const QMMolecule& mol,
    const AOBasis& full_dftbasis, const Vxc_Grid& grid) {
  // Linear sum of all four terms -- see this function's own header
  // comment for why the product-rule decomposition means their sum is
  // exactly the full derivative, with no additional cross-terms.
  return GridWeightDerivativeContribution(atoms, target_atom_index,
                                          density_matrix, mol, grid) +
        WeightFunctionDerivativeContribution(atoms, target_atom_index,
                                              density_matrix, full_dftbasis,
                                              grid) +
        PulayAndTranslationContribution(atoms, target_atom_index,
                                        density_matrix, full_dftbasis, grid);
}

}  // namespace xtp
}  // namespace votca
