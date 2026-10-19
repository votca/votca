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

// Standard includes
#include <exception>
#include <sstream>
#include <stdexcept>

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/basisset.h"
#include "votca/xtp/chargeresponsefit.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"

namespace votca {
namespace xtp {

namespace {
// Hirshfeld weights of all atoms at a point (rows: atoms), zero where the
// promolecule density is negligible (as HirshfeldPartition::EvaluateWeight)
Eigen::VectorXd Weights(
    const std::vector<HirshfeldPartition::AtomicReference>& atoms,
    const Eigen::Vector3d& point) {
  Eigen::VectorXd rho(Index(atoms.size()));
  for (Index j = 0; j < rho.size(); ++j) {
    rho(j) = HirshfeldPartition::EvaluateAtomicDensity(
        atoms[std::size_t(j)].basis, atoms[std::size_t(j)].density, point);
  }
  const double denominator = rho.sum();
  constexpr double kNegligibleDenominator = 1.e-12;
  if (denominator < kNegligibleDenominator) {
    return Eigen::VectorXd::Zero(rho.size());
  }
  return rho / denominator;
}
}  // namespace

std::vector<Eigen::MatrixXd> ChargeResponseFit::WeightMatrices(
    const std::vector<HirshfeldPartition::AtomicReference>& atoms,
    const AOBasis& basis, const Vxc_Grid& grid) {
  const Index n_atoms = Index(atoms.size());
  const Index n_ao = basis.AOBasisSize();
  std::vector<Eigen::MatrixXd> W(std::size_t(n_atoms),
                                 Eigen::MatrixXd::Zero(n_ao, n_ao));
  std::exception_ptr eptr = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    try {
      const GridBox& box = grid[i];
      if (!box.Matrixsize()) {
        continue;
      }
      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();
      const Index n_points = Index(points.size());
      Eigen::MatrixXd ao(box.Matrixsize(), n_points);
      Eigen::MatrixXd w(n_points, n_atoms);
      for (Index p = 0; p < n_points; ++p) {
        ao.col(p) = box.CalcAOValues(points[std::size_t(p)]).values;
        w.row(p) = weights[std::size_t(p)] *
                   Weights(atoms, points[std::size_t(p)]).transpose();
      }
      std::vector<Eigen::MatrixXd> local(static_cast<std::size_t>(n_atoms));
      for (Index k = 0; k < n_atoms; ++k) {
        if (w.col(k).cwiseAbs().maxCoeff() == 0.0) {
          continue;
        }
        local[std::size_t(k)] = ao * w.col(k).asDiagonal() * ao.transpose();
      }
#pragma omp critical
      {
        for (Index k = 0; k < n_atoms; ++k) {
          if (local[std::size_t(k)].size() > 0) {
            box.AddtoBigMatrix(W[std::size_t(k)], local[std::size_t(k)]);
          }
        }
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr) {
          eptr = std::current_exception();
        }
      }
    }
  }
  if (eptr) {
    std::rethrow_exception(eptr);
  }
  for (Eigen::MatrixXd& m : W) {
    m = 0.5 * (m + m.transpose()).eval();
  }
  return W;
}

std::vector<Eigen::Matrix3d> ChargeResponseFit::AtomicPolarizabilities(
    const std::vector<HirshfeldPartition::AtomicReference>& atoms,
    const std::vector<Eigen::Vector3d>& positions, const AOBasis& basis,
    const Vxc_Grid& grid, const std::array<Eigen::MatrixXd, 3>& dP_field) {
  const Index n_atoms = Index(atoms.size());
  if (Index(positions.size()) != n_atoms) {
    throw std::runtime_error(
        "ChargeResponseFit: positions and Hirshfeld references differ in "
        "size");
  }
  for (const Eigen::MatrixXd& m : dP_field) {
    if (m.rows() != basis.AOBasisSize() || m.cols() != basis.AOBasisSize()) {
      throw std::runtime_error(
          "ChargeResponseFit: density response does not match the basis");
    }
  }
  std::vector<Eigen::Matrix3d> alpha(std::size_t(n_atoms),
                                     Eigen::Matrix3d::Zero());
  std::exception_ptr eptr = nullptr;
#pragma omp parallel for schedule(guided)
  for (Index i = 0; i < grid.getBoxesSize(); ++i) {
    try {
      const GridBox& box = grid[i];
      if (!box.Matrixsize()) {
        continue;
      }
      std::array<Eigen::MatrixXd, 3> dP_box;
      for (Index c = 0; c < 3; ++c) {
        dP_box[std::size_t(c)] =
            box.ReadFromBigMatrix(dP_field[std::size_t(c)]);
      }
      const std::vector<Eigen::Vector3d>& points = box.getGridPoints();
      const std::vector<double>& weights = box.getGridWeights();
      std::vector<Eigen::Matrix3d> local(std::size_t(n_atoms),
                                         Eigen::Matrix3d::Zero());
      for (std::size_t p = 0; p < points.size(); ++p) {
        const Eigen::VectorXd w = Weights(atoms, points[p]);
        if (w.cwiseAbs().maxCoeff() == 0.0) {
          continue;
        }
        const Eigen::VectorXd ao = box.CalcAOValues(points[p]).values;
        Eigen::Vector3d rho;
        for (Index c = 0; c < 3; ++c) {
          rho(c) = ao.dot(dP_box[std::size_t(c)] * ao);
        }
        for (Index k = 0; k < n_atoms; ++k) {
          const Eigen::Vector3d d = points[p] - positions[std::size_t(k)];
          local[std::size_t(k)] -= (weights[p] * w(k)) * d * rho.transpose();
        }
      }
#pragma omp critical
      {
        for (Index k = 0; k < n_atoms; ++k) {
          alpha[std::size_t(k)] += local[std::size_t(k)];
        }
      }
    } catch (...) {
#pragma omp critical
      {
        if (!eptr) {
          eptr = std::current_exception();
        }
      }
    }
  }
  if (eptr) {
    std::rethrow_exception(eptr);
  }
  return alpha;
}

ChargeResponseFit::Result ChargeResponseFit::Analyse(
    const std::vector<Eigen::Vector3d>& positions,
    const std::vector<Eigen::MatrixXd>& W,
    const std::vector<Eigen::MatrixXd>& dP_potential,
    const std::array<Eigen::MatrixXd, 3>& dipole,
    const std::array<Eigen::MatrixXd, 3>& dP_field,
    const std::vector<Eigen::Matrix3d>& alpha_atoms) {
  const Index n = Index(positions.size());
  if (Index(W.size()) != n || Index(dP_potential.size()) != n ||
      Index(alpha_atoms.size()) != n || n < 2) {
    throw std::runtime_error(
        "ChargeResponseFit::Analyse: inconsistent number of atoms");
  }
  Result r;
  Eigen::MatrixXd pi(n, n);
  for (Index j = 0; j < n; ++j) {
    for (Index k = 0; k < n; ++k) {
      pi(j, k) =
          -W[std::size_t(j)].cwiseProduct(dP_potential[std::size_t(k)]).sum();
    }
  }
  r.pi_asymmetry = (pi - pi.transpose()).norm() / pi.norm();
  const Eigen::MatrixXd projector =
      Eigen::MatrixXd::Identity(n, n) -
      Eigen::MatrixXd::Constant(n, n, 1.0 / double(n));
  r.pi = projector * (0.5 * (pi + pi.transpose())) * projector;
  r.pi = 0.5 * (r.pi + r.pi.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(r.pi,
                                                    Eigen::EigenvaluesOnly);
  r.pi_eigenvalues = es.eigenvalues();
  if (r.pi_eigenvalues(0) < -1e-6 * r.pi_eigenvalues(n - 1)) {
    std::stringstream message;
    message << "ChargeResponseFit: Pi is not positive semidefinite (lowest "
               "eigenvalue "
            << r.pi_eigenvalues(0)
            << "); tighten the SCF convergence or change the step";
    throw std::runtime_error(message.str());
  }

  for (Index a = 0; a < 3; ++a) {
    for (Index c = 0; c < 3; ++c) {
      r.alpha(a, c) =
          -dipole[std::size_t(a)].cwiseProduct(dP_field[std::size_t(c)]).sum();
    }
  }
  r.alpha = 0.5 * (r.alpha + r.alpha.transpose()).eval();

  Eigen::MatrixX3d X(n, 3);
  for (Index k = 0; k < n; ++k) {
    X.row(k) = positions[std::size_t(k)].transpose();
  }
  X.rowwise() -= X.colwise().mean();
  r.alpha_ct = X.transpose() * r.pi * X;
  r.alpha_res = r.alpha - r.alpha_ct;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eres(r.alpha_res,
                                                      Eigen::EigenvaluesOnly);
  if (eres.eigenvalues()(0) <= 0.0) {
    std::stringstream message;
    message << "ChargeResponseFit: alpha - X^T Pi X is not positive definite "
               "(eigenvalues "
            << eres.eigenvalues().transpose()
            << " bohr^3); the charge response over-polarizes";
    throw std::runtime_error(message.str());
  }

  r.alpha_atoms = alpha_atoms;
  r.res_weights.resize(n);
  for (Index k = 0; k < n; ++k) {
    r.res_weights(k) = std::max(alpha_atoms[std::size_t(k)].trace(), 0.0);
  }
  if (!(r.res_weights.sum() > 0.0)) {
    throw std::runtime_error(
        "ChargeResponseFit: no positive Hirshfeld atomic polarizability");
  }
  r.res_weights /= r.res_weights.sum();
  for (Index k = 0; k < n; ++k) {
    r.alpha_res_atoms.push_back(r.res_weights(k) * r.alpha_res);
  }
  return r;
}

ChargeResponseFit::Inversion ChargeResponseFit::FindInversion(
    const QMMolecule& mol, const AOBasis& basis, double tol) {
  Inversion inv;
  const Index n = mol.size();
  Eigen::Vector3d centre = Eigen::Vector3d::Zero();
  for (const QMAtom& atom : mol) {
    centre += atom.getPos();
  }
  centre /= double(n);
  std::vector<Index> map(static_cast<std::size_t>(n), -1);
  for (Index k = 0; k < n; ++k) {
    const Eigen::Vector3d image = 2.0 * centre - mol[k].getPos();
    for (Index j = 0; j < n; ++j) {
      if ((mol[j].getPos() - image).norm() < tol &&
          mol[j].getElement() == mol[k].getElement()) {
        map[std::size_t(k)] = j;
        break;
      }
    }
    if (map[std::size_t(k)] < 0) {
      return inv;
    }
  }
  const Index n_ao = basis.AOBasisSize();
  Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_ao, n_ao);
  for (Index k = 0; k < n; ++k) {
    const std::vector<const AOShell*> from = basis.getShellsofAtom(k);
    const std::vector<const AOShell*> to =
        basis.getShellsofAtom(map[std::size_t(k)]);
    if (from.size() != to.size()) {
      return inv;
    }
    for (std::size_t s = 0; s < from.size(); ++s) {
      if (from[s]->getL() != to[s]->getL()) {
        return inv;
      }
      const double parity = (Index(from[s]->getL()) % 2 == 0) ? 1.0 : -1.0;
      for (Index f = 0; f < from[s]->getNumFunc(); ++f) {
        T(to[s]->getStartIndex() + f, from[s]->getStartIndex() + f) = parity;
      }
    }
  }
  inv.atom_map = map;
  inv.T = T;
  return inv;
}

ChargeResponseFit::Result ChargeResponseFit::Fit(
    const tools::Property& dftoptions, Orbitals& orb, double step, Logger& log,
    bool use_inversion) {
  if (orb.getSpin() != 1) {
    throw std::runtime_error(
        "ChargeResponseFit: closed-shell molecules only (spin 1)");
  }
  if (!(step > 0.0)) {
    throw std::runtime_error("ChargeResponseFit: step must be positive");
  }
  const QMMolecule& mol = orb.QMAtoms();
  const Index n_atoms = mol.size();
  // DFTEngine::Initialize takes a non-const tree
  tools::Property options = dftoptions;

  DFTSetupCache cache;
  // one SCF of mol with the electron potential energy operator v (none if
  // empty), warm-started from the ground state for the perturbed ones
  auto scf = [&](const Eigen::MatrixXd& v, Orbitals& target) {
    DFTEngine engine;
    engine.setLogger(&log);
    engine.Initialize(options);
    engine.setSetupCache(&cache);
    engine.setWarmStart(v.size() > 0);
    if (v.size() > 0) {
      engine.setExternalPotentialMatrix(v);
    }
    const bool ok = engine.Evaluate(target);
    engine.ReturnSetupCache();
    if (!ok) {
      throw std::runtime_error("ChargeResponseFit: SCF did not converge");
    }
  };

  XTP_LOG(Log::error, log) << TimeStamp()
                           << " Charge-response fit: ground state"
                           << std::flush;
  scf(Eigen::MatrixXd(), orb);

  DFTEngine setup;
  setup.setLogger(&log);
  setup.Initialize(options);
  AOBasis basis;
  {
    BasisSet bs;
    bs.Load(options.get("basisset").as<std::string>());
    basis.Fill(bs, mol);
  }
  Vxc_Grid grid;
  grid.GridSetup(setup.getGridName(), mol, basis);
  const std::vector<HirshfeldPartition::AtomicReference> atoms =
      setup.HirshfeldReferences(mol);
  const std::vector<Eigen::MatrixXd> W = WeightMatrices(atoms, basis, grid);
  {
    AOOverlap overlap;
    overlap.Fill(basis);
    Eigen::MatrixXd sum = -overlap.Matrix();
    for (const Eigen::MatrixXd& w : W) {
      sum += w;
    }
    XTP_LOG(Log::error, log)
        << TimeStamp() << " Hirshfeld weight matrices: max |sum_k W_k - S| = "
        << sum.cwiseAbs().maxCoeff() << std::flush;
  }
  AODipole dipole;
  dipole.setCenter(mol.getPos());
  dipole.Fill(basis);

  // central differences: dP = (P(+step) - P(-step)) / (2 step)
  auto response = [&](const Eigen::MatrixXd& h, const std::string& label) {
    Eigen::MatrixXd dP = Eigen::MatrixXd::Zero(h.rows(), h.cols());
    for (double sign : {1.0, -1.0}) {
      XTP_LOG(Log::error, log)
          << TimeStamp() << " Perturbation " << label
          << (sign > 0 ? " +" : " -") << step << std::flush;
      Orbitals perturbed = orb;
      scf(sign * step * h, perturbed);
      dP += sign * perturbed.DensityMatrixGroundState();
    }
    return Eigen::MatrixXd(dP / (2.0 * step));
  };
  Inversion inv;
  if (use_inversion) {
    inv = FindInversion(mol, basis);
    if (!inv.atom_map.empty()) {
      XTP_LOG(Log::error, log)
          << TimeStamp()
          << " Inversion centre: one atom of each pair is perturbed"
          << std::flush;
    }
  }
  std::vector<Eigen::MatrixXd> dP_potential(static_cast<std::size_t>(n_atoms));
  for (Index k = 0; k < n_atoms; ++k) {
    const Index image = inv.atom_map.empty() ? k : inv.atom_map[std::size_t(k)];
    if (image < k) {
      dP_potential[std::size_t(k)] =
          inv.T * dP_potential[std::size_t(image)] * inv.T.transpose();
    } else {
      dP_potential[std::size_t(k)] =
          response(W[std::size_t(k)],
                   "w_" + std::to_string(k) + " (" + mol[k].getElement() + ")");
    }
  }
  std::array<Eigen::MatrixXd, 3> dP_field;
  const std::array<std::string, 3> axes = {"x", "y", "z"};
  for (Index c = 0; c < 3; ++c) {
    dP_field[std::size_t(c)] = response(dipole.Matrix()[std::size_t(c)],
                                        "field " + axes[std::size_t(c)]);
  }

  std::vector<Eigen::Vector3d> positions;
  for (const QMAtom& atom : mol) {
    positions.push_back(atom.getPos());
  }
  const std::vector<Eigen::Matrix3d> alpha_atoms =
      AtomicPolarizabilities(atoms, positions, basis, grid, dP_field);
  return Analyse(positions, W, dP_potential, dipole.Matrix(), dP_field,
                 alpha_atoms);
}

}  // namespace xtp
}  // namespace votca
