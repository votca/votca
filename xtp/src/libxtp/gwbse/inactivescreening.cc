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
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/activeregion.h"
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/inactivescreening.h"
#include "votca/xtp/screening_kernels.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

InactiveSpace InactiveScreening::BuildSpace(const Orbitals& full, Index ncore) {
  const Eigen::MatrixXd& DB = full.getInactiveDensity();
  if (DB.size() == 0) {
    throw std::runtime_error(
        "InactiveScreening: no inactive density; run the projection-based "
        "embedding first");
  }
  const Index nocc = full.getNumberOfAlphaElectrons();
  const Eigen::MatrixXd& C = full.MOs().eigenvectors();
  const Eigen::VectorXd& e = full.MOs().eigenvalues();
  if (DB.rows() != C.rows()) {
    throw std::runtime_error(
        "InactiveScreening: inactive density and MOs have different AO "
        "dimensions (truncated basis is not supported)");
  }

  AOOverlap overlap;
  overlap.Fill(full.getDftBasis());
  const Eigen::MatrixXd SCocc = overlap.Matrix() * C.leftCols(nocc);

  // Projector of the inactive density onto the occupied space; its
  // eigenvalues are 1 for the inactive and 0 for the active orbitals.
  const Eigen::MatrixXd Q = 0.5 * SCocc.transpose() * DB * SCocc;
  const SymmetricEigenSystem q = SymmetricEigen(Q);
  const Index nB = (q.values.array() > 0.5).count();
  const Eigen::MatrixXd UB = q.vectors.rightCols(nB);

  InactiveSpace space;
  space.idempotency_error = (Q * Q - Q).cwiseAbs().maxCoeff();

  // semicanonical: diagonalize the Fock operator within the inactive space
  const Eigen::MatrixXd FB = UB.transpose() * e.head(nocc).asDiagonal() * UB;
  const SymmetricEigenSystem f = SymmetricEigen(FB);
  if (ncore < 0 || ncore > nB) {
    throw std::runtime_error("InactiveScreening: " + std::to_string(ncore) +
                             " core levels requested, but the inactive "
                             "region has only " +
                             std::to_string(nB) + " occupied orbitals");
  }
  const Index nkeep = nB - ncore;
  space.occupied = C.leftCols(nocc) * UB * f.vectors.rightCols(nkeep);
  space.occupied_energies = f.values.tail(nkeep);
  const Index nvirt = C.cols() - nocc;
  space.virtuals = C.rightCols(nvirt);
  space.virtual_energies = e.tail(nvirt);
  return space;
}

std::vector<Eigen::MatrixXd> InactiveScreening::Responses(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& occupied, const Eigen::VectorXd& occupied_energies,
    const Eigen::MatrixXd& virtuals, const Eigen::VectorXd& virtual_energies,
    double chunk_bytes, const std::vector<double>& imaginary_frequencies,
    const Eigen::MatrixXd& inv_sqrt) {
  const Index naux = auxbasis.AOBasisSize();
  const Index nocc = occupied.cols();
  const Index nvirt = virtuals.cols();
  std::vector<Eigen::MatrixXd> X(imaginary_frequencies.size(),
                                 Eigen::MatrixXd::Zero(naux, naux));
  if (nocc == 0 || nvirt == 0) {
    return X;
  }
  const double bytes_per_level = 8.0 * double(nvirt) * double(naux);
  const Index chunk = std::clamp<Index>(
      static_cast<Index>(chunk_bytes / bytes_per_level), 1, nocc);

  // the metric of the auxiliary basis is the same for every chunk
  Eigen::MatrixXd metric = inv_sqrt;
  for (Index first = 0; first < nocc; first += chunk) {
    const Index nc = std::min(chunk, nocc - first);
    Eigen::MatrixXd orbitals(occupied.rows(), nc + nvirt);
    orbitals << occupied.middleCols(first, nc), virtuals;
    // one integral pass for all frequencies
    TCMatrix_gwbse Mmn;
    Mmn.Initialize(naux, 0, nc - 1, nc, nc + nvirt - 1);
    if (metric.size() == 0) {
      Mmn.Fill(auxbasis, dftbasis, orbitals);
      metric = Mmn.InvSqrt();
    } else {
      Mmn.Fill(auxbasis, dftbasis, orbitals, metric);
    }
    for (std::size_t f = 0; f < imaginary_frequencies.size(); ++f) {
      const double w2 = imaginary_frequencies[f] * imaginary_frequencies[f];
      X[f] += WeightedGram::Compute(
          nc, naux,
          [&](Index m, Eigen::VectorXd& w) {
            const Eigen::ArrayXd d =
                virtual_energies.array() - occupied_energies(first + m);
            w = 4.0 * d / (d.square() + w2);
          },
          [&](Index m, const WeightedGram::RowsVisitor& use) { use(Mmn[m]); });
    }
  }
  return X;
}

Eigen::MatrixXd InactiveScreening::Response(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const Eigen::MatrixXd& occupied, const Eigen::VectorXd& occupied_energies,
    const Eigen::MatrixXd& virtuals, const Eigen::VectorXd& virtual_energies,
    double chunk_bytes, double imaginary_frequency) {
  return Responses(auxbasis, dftbasis, occupied, occupied_energies, virtuals,
                   virtual_energies, chunk_bytes, {imaginary_frequency})
      .front();
}

InactiveScreeningInput::Mode InactiveScreeningInput::Parse(
    const std::string& name) {
  if (name == "none") {
    return Mode::None;
  } else if (name == "block") {
    return Mode::Block;
  } else if (name == "exact") {
    return Mode::Exact;
  }
  throw std::runtime_error("inactive_screening: unknown mode '" + name +
                           "' (none, block, exact)");
}

std::string InactiveScreeningInput::Name(Mode mode) {
  switch (mode) {
    case Mode::Block:
      return "block";
    case Mode::Exact:
      return "exact";
    default:
      return "none";
  }
}

InactiveSpace InactiveScreening::NearSpace(const Orbitals& full, double radius,
                                           Index& n_far) {
  const Eigen::MatrixXd& L = full.getLMOs();
  const Index nocc = full.getNumberOfAlphaElectrons();
  if (L.cols() != nocc || L.rows() != full.MOs().eigenvectors().rows()) {
    throw std::runtime_error(
        "inactive_screening_radius needs the localized occupied orbitals of "
        "the full calculation (localize step)");
  }
  AOOverlap overlap;
  overlap.Fill(full.getDftBasis());
  const Eigen::MatrixXd& S = overlap.Matrix();
  const Eigen::MatrixXd SL = S * L;
  // inactive LMOs: weight ~1 in the inactive density
  const Eigen::VectorXd q =
      0.5 * (SL.transpose() * full.getInactiveDensity() * SL).diagonal();
  // charge centroids from Mulliken atom populations
  const QMMolecule& atoms = full.QMAtoms();
  const Eigen::MatrixXd pop = ActiveRegion::AtomPopulations(full, L);
  Eigen::MatrixXd R(3, atoms.size());
  for (Index a = 0; a < atoms.size(); ++a) {
    R.col(a) = atoms[a].getPos();
  }
  const Eigen::MatrixXd centroid =
      R * pop * pop.colwise().sum().cwiseInverse().asDiagonal();
  std::vector<Index> active;
  std::vector<Index> inactive;
  for (Index k = 0; k < nocc; ++k) {
    (q(k) > 0.5 ? inactive : active).push_back(k);
  }
  std::vector<Index> near;
  for (Index k : inactive) {
    double dmin = std::numeric_limits<double>::max();
    for (Index j : active) {
      dmin = std::min(dmin, (centroid.col(k) - centroid.col(j)).norm());
    }
    if (dmin <= radius) {
      near.push_back(k);
    }
  }
  n_far = Index(inactive.size() - near.size());
  Eigen::MatrixXd Ln(L.rows(), Index(near.size()));
  for (Index i = 0; i < Index(near.size()); ++i) {
    Ln.col(i) = L.col(near[std::size_t(i)]);
  }
  // semicanonical among the near ones: F = S C diag(e) C^T S
  const Index n = full.MOs().eigenvectors().cols();
  const Eigen::MatrixXd SC = S * full.MOs().eigenvectors();
  const Eigen::MatrixXd CL = SC.leftCols(nocc).transpose() * Ln;  // MO basis
  const Eigen::MatrixXd FL =
      CL.transpose() * full.MOs().eigenvalues().head(nocc).asDiagonal() * CL;
  const SymmetricEigenSystem f = SymmetricEigen(0.5 * (FL + FL.transpose()));
  InactiveSpace space;
  space.occupied = Ln * f.vectors;
  space.occupied_energies = f.values;
  space.virtuals = full.MOs().eigenvectors().rightCols(n - nocc);
  space.virtual_energies = full.MOs().eigenvalues().tail(n - nocc);
  return space;
}

InactiveScreeningInput InactiveScreening::Prepare(
    const Orbitals& full, InactiveScreeningInput::Mode mode, bool complete_g,
    double radius) {
  InactiveScreeningInput input;
  input.mode = mode;
  input.complete_g = complete_g && mode != InactiveScreeningInput::Mode::None;
  if (mode == InactiveScreeningInput::Mode::None) {
    return input;
  }
  if (full.getCalculationType() == "Truncated") {
    throw std::runtime_error(
        "inactive_screening is not implemented for truncated active "
        "regions (dft_in_dft.truncate_basis).");
  }
  // the inactive space is built in both modes: its size and idempotency
  // are the checks that the full MOs and the inactive density belong
  // together
  InactiveSpace space = BuildSpace(full);
  input.n_inactive = space.size();
  input.idempotency_error = space.idempotency_error;
  if (input.complete_g) {
    input.inactive_occupied = space.occupied;
    input.inactive_energies = space.occupied_energies;
  }
  if (mode == InactiveScreeningInput::Mode::Block && radius > 0.0) {
    InactiveSpace near = NearSpace(full, radius, input.n_far);
    input.radius = radius;
    input.occupied = std::move(near.occupied);
    input.occupied_energies = std::move(near.occupied_energies);
  } else if (mode == InactiveScreeningInput::Mode::Block) {
    input.occupied = std::move(space.occupied);
    input.occupied_energies = std::move(space.occupied_energies);
  } else {
    if (radius > 0.0) {
      throw std::runtime_error(
          "inactive_screening_radius is only available for "
          "inactive_screening block");
    }
    const Index nocc = full.getNumberOfAlphaElectrons();
    input.occupied = full.MOs().eigenvectors().leftCols(nocc);
    input.occupied_energies = full.MOs().eigenvalues().head(nocc);
  }
  input.virtuals = std::move(space.virtuals);
  input.virtual_energies = std::move(space.virtual_energies);
  return input;
}

tools::EigenSystem InactiveScreening::CompleteGOrbitals(
    const tools::EigenSystem& embedded, Index n_active_occupied,
    const InactiveScreeningInput& input) {
  const Index nB = input.inactive_occupied.cols();
  const Index nvirt = input.virtuals.cols();
  const Index nA = n_active_occupied;
  if (nA + nvirt > embedded.eigenvectors().cols()) {
    throw std::runtime_error(
        "CompleteGOrbitals: the embedded MOs have fewer levels than the "
        "active occupied plus virtual ones");
  }
  tools::EigenSystem result;
  result.eigenvectors().resize(embedded.eigenvectors().rows(), nB + nA + nvirt);
  result.eigenvectors() << input.inactive_occupied,
      embedded.eigenvectors().leftCols(nA + nvirt);
  result.eigenvalues().resize(nB + nA + nvirt);
  result.eigenvalues() << input.inactive_energies,
      embedded.eigenvalues().head(nA + nvirt);
  return result;
}

Eigen::MatrixXd InactiveScreening::Kernel(const AOBasis& auxbasis,
                                          const AOBasis& dftbasis,
                                          const InactiveScreeningInput& input,
                                          const tools::EigenSystem& embedded,
                                          Index homo, double chunk_bytes,
                                          double imaginary_frequency,
                                          Index first_active) {
  Eigen::MatrixXd Y = Response(
      auxbasis, dftbasis, input.occupied, input.occupied_energies,
      input.virtuals, input.virtual_energies, chunk_bytes, imaginary_frequency);
  if (input.mode == InactiveScreeningInput::Mode::Exact) {
    const Index nA = homo + 1 - first_active;
    const Index nvirt = input.virtuals.cols();
    if (homo + 1 + nvirt > embedded.eigenvectors().cols()) {
      throw std::runtime_error(
          "inactive_screening: the embedded MOs have fewer levels than the "
          "full calculation's occupied and virtual ones");
    }
    Y -= Response(auxbasis, dftbasis,
                  embedded.eigenvectors().middleCols(first_active, nA),
                  embedded.eigenvalues().segment(first_active, nA),
                  embedded.eigenvectors().middleCols(homo + 1, nvirt),
                  embedded.eigenvalues().segment(homo + 1, nvirt), chunk_bytes,
                  imaginary_frequency);
  }
  return 0.5 * (Y + Y.transpose());
}

std::vector<Eigen::MatrixXd> InactiveScreening::Kernels(
    const AOBasis& auxbasis, const AOBasis& dftbasis,
    const InactiveScreeningInput& input, const TCMatrix_gwbse& active,
    const Eigen::VectorXd& rpa_energies, Index homo, Index rpamin, Index rpamax,
    double chunk_bytes, const std::vector<double>& imaginary_frequencies,
    double virtual_cutoff) {
  // the virtuals of the full calculation inside the RPA range; above them
  // (no complete G) sit the inactive levels the projector shifted up,
  // which are left out of both terms
  Index nvirt = std::min(rpamax - homo, input.virtuals.cols());
  if (virtual_cutoff > 0.0 && nvirt > 0) {
    // only virtuals up to the cutoff above the LUMO, in both terms
    const double top = input.virtual_energies(0) + virtual_cutoff;
    nvirt = std::min(
        nvirt,
        Index((input.virtual_energies.head(nvirt).array() <= top).count()));
  }
  std::vector<Eigen::MatrixXd> Y = Responses(
      auxbasis, dftbasis, input.occupied, input.occupied_energies,
      input.virtuals.leftCols(nvirt), input.virtual_energies.head(nvirt),
      chunk_bytes, imaginary_frequencies, active.InvSqrt());
  if (input.mode == InactiveScreeningInput::Mode::Exact) {
    // X_A from the integrals GW already holds (excluded holes, e.g. the
    // inactive levels of complete G, are left out as in the RPA)
    if (active.HasAdditionalResponse() || active.AdditionalPoleCount() > 0) {
      throw std::runtime_error(
          "inactive_screening: X_A must be built before the inactive "
          "response is attached to the integrals");
    }
    const Index first = rpamin + active.ExcludedHoles();
    const Index mmin = active.get_mmin();
    const Index nmin = active.get_nmin();
    for (std::size_t f = 0; f < Y.size(); ++f) {
      const double w2 = imaginary_frequencies[f] * imaginary_frequencies[f];
      Y[f] -= WeightedGram::Compute(
          homo + 1 - first, active.auxsize(),
          [&](Index m, Eigen::VectorXd& w) {
            const Eigen::ArrayXd d =
                rpa_energies.segment(homo + 1 - rpamin, nvirt).array() -
                rpa_energies(first + m - rpamin);
            w = 4.0 * d / (d.square() + w2);
          },
          [&](Index m, const WeightedGram::RowsVisitor& use) {
            use(active[first + m - mmin].middleRows(homo + 1 - nmin, nvirt));
          });
    }
  }
  for (Eigen::MatrixXd& y : Y) {
    y = 0.5 * (y + y.transpose());
  }
  return Y;
}

InactiveScreening::PoleFit InactiveScreening::FitPoles(
    const Eigen::MatrixXd& Y0, const Eigen::MatrixXd& Yi, double w_i,
    double tolerance) {
  PoleFit fit;
  const SymmetricEigenSystem y = SymmetricEigen(0.5 * (Y0 + Y0.transpose()));
  const double cut = tolerance * std::max(y.values.maxCoeff(), 0.0);
  const Index keep = (y.values.array() > cut).count();
  fit.dropped = y.values.size() - keep;
  if (keep == 0) {
    return fit;
  }
  // Y(0) = G G^T on its range, G = U s^1/2; P = U s^-1/2 = G^+T
  const Eigen::VectorXd s = y.values.tail(keep);
  const Eigen::MatrixXd U = y.vectors.rightCols(keep);
  const Eigen::MatrixXd P = U * s.cwiseSqrt().cwiseInverse().asDiagonal();
  const Eigen::MatrixXd Yp = P.transpose() * Yi * P;
  const SymmetricEigenSystem q = SymmetricEigen(0.5 * (Yp + Yp.transpose()));
  // l = W^2/(W^2 + w_i^2): W from 1e-3 Ha (l_min) to 1e3 Ha (l_max)
  const double l_min = 1e-6 / (1e-6 + w_i * w_i);
  const double l_max = 1e6 / (1e6 + w_i * w_i);
  Eigen::VectorXd l = q.values;
  for (Index k = 0; k < l.size(); ++k) {
    if (l(k) <= 0.0 || l(k) >= 1.0) {
      // no real pole: keep the mode's static response, as a pole far above
      // every excitation (a spurious low-lying pole would be worse for
      // the exact self-energy than a static contribution)
      ++fit.clipped;
      l(k) = l_max;
    } else {
      l(k) = std::clamp(l(k), l_min, l_max);
    }
  }
  fit.energies = w_i * (l.array() / (1.0 - l.array())).sqrt();
  const Eigen::MatrixXd G = U * s.cwiseSqrt().asDiagonal() * q.vectors;
  fit.vectors = G * (0.25 * fit.energies.array()).sqrt().matrix().asDiagonal();
  return fit;
}

Eigen::MatrixXd InactiveScreening::ReactionField(const Eigen::MatrixXd& Y,
                                                 const Eigen::MatrixXd& R_env) {
  const Index n = Y.rows();
  const Eigen::MatrixXd one = Eigen::MatrixXd::Identity(n, n);
  Eigen::MatrixXd u_inv = one + Y;
  if (R_env.size() > 0) {
    if (R_env.rows() != n) {
      throw std::runtime_error(
          "InactiveScreening: environment reaction field and inactive "
          "response have different auxiliary dimensions");
    }
    // (1 + R_env)^-1 + Y; 1 + R_env is positive definite (Check)
    u_inv = InverseSPD(one + R_env) + Y;
  }
  Eigen::MatrixXd R = InverseSPD(u_inv) - one;
  return 0.5 * (R + R.transpose());
}

}  // namespace xtp
}  // namespace votca
