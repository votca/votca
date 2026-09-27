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
#include <cmath>

// Local VOTCA includes
#include "votca/xtp/ewaldreciprocalspacesum.h"

namespace votca {
namespace xtp {

namespace {
constexpr double kPi = 3.14159265358979323846;
}  // namespace

EwaldReciprocalSpaceSum::EwaldReciprocalSpaceSum(const Eigen::Matrix3d& box,
                                                 const EwaldRegistry& registry,
                                                 double alpha, double k_max)
    : box_(box),
      volume_(box.col(0).dot(box.col(1).cross(box.col(2)))),
      registry_(registry),
      alpha_(alpha),
      k_max_(k_max) {
  kvectors_ = GenerateKVectors();
}

Eigen::Matrix3d EwaldReciprocalSpaceSum::SelfFieldMatrix() const {
  // BUG FIX (this session): this used to evaluate the discrete k-lattice
  // sum (4*pi/V) * sum_{k!=0} exp(-k^2/4*alpha^2)/k^2 * k k^T over this
  // class's own k-vector set, on the reasoning that doing so reproduces
  // exactly the self-term this class's own sum numerically produces.
  // That reasoning is wrong, for a reason that has nothing to do with
  // k_max: evaluated at r = 0, that sum is the erf-screened field a site
  // feels from itself AND from every one of its own periodic images. Only
  // the n = 0 piece is the spurious self-interaction that wants removing;
  // the n != 0 image terms are real physical interactions that must be
  // kept. Subtracting the whole lattice sum therefore removes genuine
  // physics along with the artifact.
  //
  // The quantity actually wanted is the r -> 0 limit of the erf-screened
  // dipole field, i.e. the standard analytic Ewald self-term below, which
  // is what legacy applies too (EwdInteractor::FU12_ERF_At_By's own
  // R1 < 1e-2 branch: 4/3 * a3 * rSqrtPi * U1). An earlier version of
  // this comment claimed, from a trace of legacy's reciprocal-space code,
  // that legacy applies no self-correction at all -- that was simply
  // wrong: legacy applies it in real space, as a separate "atomic ERF
  // self-interaction correction" pass (PolarBackground's own step (5)),
  // which a reciprocal-space-only trace never reaches.
  //
  // The two differ badly at realistic box sizes, and NOT by a small
  // amount that tightening k_max would fix. The lattice sum approaches
  // the analytic value only when the k-lattice resolves the Gaussian
  // exp(-k^2/4*alpha^2); with spacing 2*pi/L against width 2*alpha that
  // needs L*alpha/pi >> 1. For a measured case (L = 66.1 bohr, alpha =
  // 0.1058 bohr^-1) there are only ~2.2 k-points per Gaussian width and
  // the lattice sum comes out 1.63% low -- which matched, to five
  // digits, a directly measured 1.63% deficit in this term against
  // legacy across 5000 sites, with the small residual anisotropy of the
  // discrete sum showing up as a cos(angle) of 0.999993 rather than 1.
  //
  // Isotropic by construction, so no anisotropy artifact remains.
  const double self_term =
      (4.0 / 3.0) * alpha_ * alpha_ * alpha_ / std::sqrt(kPi);
  return self_term * Eigen::Matrix3d::Identity();
}

std::vector<EwaldReciprocalSpaceSum::KVector>
    EwaldReciprocalSpaceSum::GenerateKVectors() const {
  // Reciprocal lattice vectors: columns of 2*pi*(box^-1)^T, the standard
  // dual basis (b_i . a_j = 2*pi*delta_ij).
  const Eigen::Matrix3d recip = 2.0 * kPi * box_.inverse().transpose();
  const Eigen::Vector3d b1 = recip.col(0);
  const Eigen::Vector3d b2 = recip.col(1);
  const Eigen::Vector3d b3 = recip.col(2);

  // Generous per-axis bound: enough reciprocal-lattice steps along each
  // primitive direction to reach k_max even along the shortest axis.
  const Index n1 = Index(std::ceil(k_max_ / b1.norm())) + 1;
  const Index n2 = Index(std::ceil(k_max_ / b2.norm())) + 1;
  const Index n3 = Index(std::ceil(k_max_ / b3.norm())) + 1;

  std::vector<KVector> kvecs;
  for (Index i1 = -n1; i1 <= n1; ++i1) {
    for (Index i2 = -n2; i2 <= n2; ++i2) {
      for (Index i3 = -n3; i3 <= n3; ++i3) {
        if (i1 == 0 && i2 == 0 && i3 == 0) {
          continue;  // k=0 term is handled separately (uniform background
                     // / net-charge correction), not part of this sum.
        }
        Eigen::Vector3d k = double(i1) * b1 + double(i2) * b2 + double(i3) * b3;
        double k2 = k.squaredNorm();
        if (k2 <= k_max_ * k_max_) {
          kvecs.push_back({k, k2});
        }
      }
    }
  }
  return kvecs;
}

std::vector<std::complex<double>>
    EwaldReciprocalSpaceSum::TotalStructureFactors(
        EwaldChargeState source_state, const ProgressCallback& progress) const {
  std::vector<std::complex<double>> S(kvectors_.size(),
                                      std::complex<double>(0.0, 0.0));

  // Flatten every source site into contiguous arrays first. Two reasons:
  // the k-vector loop below is parallelized over k rather than over
  // sites, so it needs random access to the site data rather than a
  // nested registry walk; and streaming three packed arrays is far
  // friendlier to cache than chasing segment objects through a std::map
  // for each of ~1e4 k-vectors.
  std::vector<double> q_flat;
  std::vector<Eigen::Vector3d> mu_flat;
  std::vector<Eigen::Vector3d> pos_flat;
  for (Index source_id : registry_.AllIds()) {
    if (!registry_.Has(source_id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(source_id, source_state);
    for (const PolarSite& site : segment) {
      q_flat.push_back(site.getCharge());
      mu_flat.push_back(site.getStaticDipole() + site.getInducedDipole());
      pos_flat.push_back(site.getPos());
    }
  }
  const Index n_sites = Index(q_flat.size());
  const Index n_k = Index(kvectors_.size());

  // Parallelized over K-VECTORS, not over sites. The natural reading --
  // sites outer, k inner -- makes S[idx] a reduction target shared by
  // every thread, needing either atomics or per-thread copies of the
  // whole S array. Inverting the loops gives each thread sole ownership
  // of the S entries it writes, so no reduction, no atomics, and the
  // sum over sites for a given k happens in a fixed order regardless of
  // thread count -- the result is bitwise identical however many
  // threads run it.
  //
  // The k-loop is chunked so the progress callback can still be invoked
  // between chunks, from the serial region: calling it from inside the
  // parallel loop would need the callback itself to be thread-safe,
  // which is not part of its contract.
  const Index chunk = std::max<Index>(1, n_k / 20);
  for (Index k_begin = 0; k_begin < n_k; k_begin += chunk) {
    const Index k_end = std::min(n_k, k_begin + chunk);
#pragma omp parallel for schedule(static)
    for (Index idx = k_begin; idx < k_end; ++idx) {
      const Eigen::Vector3d& k = kvectors_[std::size_t(idx)].k;
      std::complex<double> acc(0.0, 0.0);
      for (Index n = 0; n < n_sites; ++n) {
        const double kr = k.dot(pos_flat[std::size_t(n)]);
        // exp(-i*kr) with |exp| == 1: std::polar avoids the redundant
        // std::exp(0) that std::exp(std::complex) would evaluate.
        const std::complex<double> phase = std::polar(1.0, -kr);
        const double k_dot_mu = k.dot(mu_flat[std::size_t(n)]);
        acc += std::complex<double>(q_flat[std::size_t(n)], -k_dot_mu) * phase;
      }
      S[std::size_t(idx)] = acc;
    }
    if (progress) {
      progress(std::size_t(k_end), std::size_t(n_k));
    }
  }
  return S;
}

double EwaldReciprocalSpaceSum::CalcStaticEnergyBetween(
    const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>>& foreground,
    const std::vector<const PolarSite*>& background_exclusions,
    EwaldChargeState source_state) const {
  // See this method's own declaration for the formula, for why the two
  // structure factors are accumulated separately, and for why the
  // exclusion list is a separate argument rather than being read off the
  // foreground.

  // Foreground moments, taken from the sites themselves and so in the
  // job's own charge state, at the positions they actually occupy. The
  // positions are supplied rather than read off the sites because a
  // foreground segment sits at one particular periodic image and the
  // phase factor must use that image's position.
  std::vector<double> q_fg;
  std::vector<Eigen::Vector3d> mu_fg;
  std::vector<Eigen::Vector3d> pos_fg;
  q_fg.reserve(foreground.size());
  mu_fg.reserve(foreground.size());
  pos_fg.reserve(foreground.size());
  for (const auto& entry : foreground) {
    const PolarSite& site = *entry.first;
    q_fg.push_back(site.getCharge());
    mu_fg.push_back(site.getStaticDipole());
    pos_fg.push_back(entry.second);
  }

  // Background sites: every registered site EXCEPT the ones the caller
  // named. Identity is by address, so a listed site is held out exactly
  // once.
  const std::vector<const PolarSite*>& fg_sites = background_exclusions;
  std::vector<double> q_bg;
  std::vector<Eigen::Vector3d> mu_bg;
  std::vector<Eigen::Vector3d> pos_bg;
  for (Index source_id : registry_.AllIds()) {
    if (!registry_.Has(source_id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(source_id, source_state);
    for (const PolarSite& site : segment) {
      bool is_foreground = false;
      for (const PolarSite* fg : fg_sites) {
        if (fg == &site) {
          is_foreground = true;
          break;
        }
      }
      if (is_foreground) {
        continue;
      }
      q_bg.push_back(site.getCharge());
      mu_bg.push_back(site.getStaticDipole());
      pos_bg.push_back(site.getPos());
    }
  }

  const Index n_fg = Index(q_fg.size());
  const Index n_bg = Index(q_bg.size());
  const Index n_k = Index(kvectors_.size());
  const double prefactor = 4.0 * kPi / volume_;

  // Parallel over k-vectors, as TotalStructureFactors is and for the
  // same reason: each thread owns its own k and the sum over sites for a
  // given k happens in a fixed order, so the result does not depend on
  // thread count.
  double energy = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : energy)
  for (Index idx = 0; idx < n_k; ++idx) {
    const Eigen::Vector3d& k = kvectors_[std::size_t(idx)].k;
    const double k2 = kvectors_[std::size_t(idx)].k2;

    std::complex<double> s_fg(0.0, 0.0);
    for (Index n = 0; n < n_fg; ++n) {
      const double kr = k.dot(pos_fg[std::size_t(n)]);
      const std::complex<double> phase = std::polar(1.0, -kr);
      const double k_dot_mu = k.dot(mu_fg[std::size_t(n)]);
      s_fg += std::complex<double>(q_fg[std::size_t(n)], -k_dot_mu) * phase;
    }

    std::complex<double> s_bg(0.0, 0.0);
    for (Index n = 0; n < n_bg; ++n) {
      const double kr = k.dot(pos_bg[std::size_t(n)]);
      const std::complex<double> phase = std::polar(1.0, -kr);
      const double k_dot_mu = k.dot(mu_bg[std::size_t(n)]);
      s_bg += std::complex<double>(q_bg[std::size_t(n)], -k_dot_mu) * phase;
    }

    const double weight = std::exp(-k2 / (4.0 * alpha_ * alpha_)) / k2;
    energy += prefactor * weight * (std::conj(s_fg) * s_bg).real();
  }
  return energy;
}

double EwaldReciprocalSpaceSum::CalcInducedSourceEnergyBetween(
    const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>>& foreground,
    const std::vector<const PolarSite*>& background_exclusions,
    EwaldChargeState source_state) const {
  // See this method's own declaration. Structurally identical to
  // CalcStaticEnergyBetween above, with one difference: the background
  // contributes its INDUCED dipoles and no charge, rather than its
  // permanent moments.

  // Foreground: permanent moments, at the positions actually occupied.
  std::vector<double> q_fg;
  std::vector<Eigen::Vector3d> mu_fg;
  std::vector<Eigen::Vector3d> pos_fg;
  q_fg.reserve(foreground.size());
  mu_fg.reserve(foreground.size());
  pos_fg.reserve(foreground.size());
  for (const auto& entry : foreground) {
    const PolarSite& site = *entry.first;
    q_fg.push_back(site.getCharge());
    mu_fg.push_back(site.getStaticDipole());
    pos_fg.push_back(entry.second);
  }

  // Background: induced dipoles only. No charge term -- an induced
  // dipole carries none, and the background's permanent charges are
  // already accounted for by CalcStaticEnergyBetween.
  std::vector<Eigen::Vector3d> mu_bg;
  std::vector<Eigen::Vector3d> pos_bg;
  for (Index source_id : registry_.AllIds()) {
    if (!registry_.Has(source_id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(source_id, source_state);
    for (const PolarSite& site : segment) {
      bool is_excluded = false;
      for (const PolarSite* skip : background_exclusions) {
        if (skip == &site) {
          is_excluded = true;
          break;
        }
      }
      if (is_excluded) {
        continue;
      }
      mu_bg.push_back(site.getInducedDipole());
      pos_bg.push_back(site.getPos());
    }
  }

  const Index n_fg = Index(q_fg.size());
  const Index n_bg = Index(mu_bg.size());
  const Index n_k = Index(kvectors_.size());
  const double prefactor = 4.0 * kPi / volume_;

  double energy = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : energy)
  for (Index idx = 0; idx < n_k; ++idx) {
    const Eigen::Vector3d& k = kvectors_[std::size_t(idx)].k;
    const double k2 = kvectors_[std::size_t(idx)].k2;

    std::complex<double> s_fg(0.0, 0.0);
    for (Index n = 0; n < n_fg; ++n) {
      const double kr = k.dot(pos_fg[std::size_t(n)]);
      const std::complex<double> phase = std::polar(1.0, -kr);
      const double k_dot_mu = k.dot(mu_fg[std::size_t(n)]);
      s_fg += std::complex<double>(q_fg[std::size_t(n)], -k_dot_mu) * phase;
    }

    std::complex<double> s_bg(0.0, 0.0);
    for (Index n = 0; n < n_bg; ++n) {
      const double kr = k.dot(pos_bg[std::size_t(n)]);
      const std::complex<double> phase = std::polar(1.0, -kr);
      const double k_dot_mu = k.dot(mu_bg[std::size_t(n)]);
      s_bg += std::complex<double>(0.0, -k_dot_mu) * phase;
    }

    const double weight = std::exp(-k2 / (4.0 * alpha_ * alpha_)) / k2;
    energy += prefactor * weight * (std::conj(s_fg) * s_bg).real();
  }
  return energy;
}

template <enum Estatic CE>
void EwaldReciprocalSpaceSum::AddFieldAt(PolarSite& target,
                                         EwaldChargeState source_state) const {
  std::vector<PolarSite*> targets = {&target};
  AddFieldAtMany<CE>(targets, source_state);
}

template <enum Estatic CE>
void EwaldReciprocalSpaceSum::AddFieldAtMany(
    const std::vector<PolarSite*>& targets, EwaldChargeState source_state,
    const ProgressCallback& progress) const {
  const std::vector<std::complex<double>> S =
      TotalStructureFactors(source_state, progress);
  const double prefactor = 4.0 * kPi / volume_;

  // Parallel over targets. The structure factors S were reduced over
  // every site above and are read-only here, and each iteration writes
  // only to its own target's accumulator, so distinct targets cannot
  // collide. The progress callback is deliberately not invoked from
  // inside this loop -- it is called during TotalStructureFactors, which
  // stays serial.
  const Index n_targets = Index(targets.size());
#pragma omp parallel for schedule(static)
  for (Index t_i = 0; t_i < n_targets; ++t_i) {
    PolarSite& target = *targets[std::size_t(t_i)];

    const Eigen::Vector3d r = target.getPos();
    Eigen::Vector3d field = Eigen::Vector3d::Zero();

    for (std::size_t idx = 0; idx < kvectors_.size(); ++idx) {
      const KVector& kv = kvectors_[idx];
      double weight = std::exp(-kv.k2 / (4.0 * alpha_ * alpha_)) / kv.k2;

      // E(r) = (4*pi/V) * sum_{k!=0} (k/k^2) * exp(-k^2/4a^2) *
      //        Im[ S(k) * exp(i*k.r) ]
      // See class documentation for the derivation.
      const std::complex<double> phase = std::polar(1.0, kv.k.dot(r));
      double im_part = (S[idx] * phase).imag();
      field += prefactor * weight * im_part * kv.k;
    }

    if (CE == Estatic::noE_V) {
      target.V_noE() += field;
    } else {
      target.V() += field;
    }
  }
}

Eigen::VectorXd EwaldReciprocalSpaceSum::PotentialAtMany(
    const std::vector<Eigen::Vector3d>& points, EwaldChargeState source_state,
    const ProgressCallback& progress) const {
  const std::vector<std::complex<double>> S =
      TotalStructureFactors(source_state, progress);
  const double prefactor = 4.0 * kPi / volume_;
  const Index n_points = Index(points.size());
  const std::size_t n_k = kvectors_.size();
  Eigen::VectorXd phi = Eigen::VectorXd::Zero(n_points);

  // Everything that depends on k alone, folded once. The Gaussian weight
  // in particular is an exp() that does not vary over the points, so
  // leaving it in the inner loop costs one transcendental per (point, k)
  // pair -- on a DFT integration grid that is by far the most expensive
  // thing in this routine, and none of it does any work.
  //
  // Kept as separate real and imaginary parts rather than std::complex so
  // the inner loop is two multiplies against a cos/sin pair, with no
  // complex multiply and no temporary.
  std::vector<double> c_re(n_k);
  std::vector<double> c_im(n_k);
  const double inv_four_alpha2 = 1.0 / (4.0 * alpha_ * alpha_);
  for (std::size_t idx = 0; idx < n_k; ++idx) {
    const double weight = prefactor *
                          std::exp(-kvectors_[idx].k2 * inv_four_alpha2) /
                          kvectors_[idx].k2;
    c_re[idx] = weight * S[idx].real();
    c_im[idx] = weight * S[idx].imag();
  }

  // Parallel over points rather than over k, the opposite of
  // TotalStructureFactors above: S is read-only here and each iteration
  // owns its own accumulator, and a DFT grid has far more points than
  // this class has k-vectors.
#pragma omp parallel for schedule(static)
  for (Index p = 0; p < n_points; ++p) {
    const Eigen::Vector3d& r = points[std::size_t(p)];
    double acc = 0.0;
    for (std::size_t idx = 0; idx < n_k; ++idx) {
      // Re[(c_re + i c_im) e^{i theta}] = c_re cos(theta) - c_im sin(theta).
      // Written out rather than through std::polar and a complex multiply,
      // which compute the same two trig calls plus four multiplies.
      const double theta = kvectors_[idx].k.dot(r);
      acc += c_re[idx] * std::cos(theta) - c_im[idx] * std::sin(theta);
    }
    phi[p] = acc;
  }
  return phi;
}

template void EwaldReciprocalSpaceSum::AddFieldAt<Estatic::V>(
    PolarSite&, EwaldChargeState) const;
template void EwaldReciprocalSpaceSum::AddFieldAt<Estatic::noE_V>(
    PolarSite&, EwaldChargeState) const;

template void EwaldReciprocalSpaceSum::AddFieldAtMany<Estatic::V>(
    const std::vector<PolarSite*>&, EwaldChargeState,
    const EwaldReciprocalSpaceSum::ProgressCallback&) const;
template void EwaldReciprocalSpaceSum::AddFieldAtMany<Estatic::noE_V>(
    const std::vector<PolarSite*>&, EwaldChargeState,
    const EwaldReciprocalSpaceSum::ProgressCallback&) const;

}  // namespace xtp
}  // namespace votca
