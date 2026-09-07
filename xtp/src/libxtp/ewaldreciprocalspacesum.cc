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
  Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
  const double prefactor = 4.0 * kPi / volume_;
  for (const KVector& kv : kvectors_) {
    double weight = std::exp(-kv.k2 / (4.0 * alpha_ * alpha_)) / kv.k2;
    M += prefactor * weight * (kv.k * kv.k.transpose());
  }
  return M;
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
  const std::complex<double> i(0.0, 1.0);

  // Count total sites upfront, purely for progress reporting (this is
  // the O(N_k * N_sites) loop -- the dominant cost of a reciprocal-space
  // evaluation -- so it is the one place worth reporting progress from).
  std::size_t total_sites = 0;
  if (progress) {
    for (Index source_id : registry_.AllIds()) {
      if (registry_.Has(source_id, source_state)) {
        total_sites +=
            std::size_t(registry_.Get(source_id, source_state).size());
      }
    }
  }
  std::size_t sites_done = 0;
  // Report at most ~20 times over the whole loop, regardless of how many
  // sites there are, so the callback itself never becomes the bottleneck.
  const std::size_t report_every =
      std::max<std::size_t>(1, total_sites / 20);

  for (Index source_id : registry_.AllIds()) {
    if (!registry_.Has(source_id, source_state)) {
      continue;
    }
    const PolarSegment& segment = registry_.Get(source_id, source_state);
    for (const PolarSite& site : segment) {
      const double q = site.getCharge();
      const Eigen::Vector3d mu =
          site.getStaticDipole() + site.getInducedDipole();
      const Eigen::Vector3d& pos = site.getPos();
      for (std::size_t idx = 0; idx < kvectors_.size(); ++idx) {
        const Eigen::Vector3d& k = kvectors_[idx].k;
        const double kr = k.dot(pos);
        const std::complex<double> phase = std::exp(-i * kr);
        const double k_dot_mu = k.dot(mu);
        S[idx] += (q - i * k_dot_mu) * phase;
      }
      if (progress) {
        ++sites_done;
        if (sites_done % report_every == 0 || sites_done == total_sites) {
          progress(sites_done, total_sites);
        }
      }
    }
  }
  return S;
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
  const std::complex<double> i(0.0, 1.0);
  const double prefactor = 4.0 * kPi / volume_;

  for (PolarSite* target_ptr : targets) {
    PolarSite& target = *target_ptr;

    const Eigen::Vector3d r = target.getPos();
    Eigen::Vector3d field = Eigen::Vector3d::Zero();

    for (std::size_t idx = 0; idx < kvectors_.size(); ++idx) {
      const KVector& kv = kvectors_[idx];
      double weight = std::exp(-kv.k2 / (4.0 * alpha_ * alpha_)) / kv.k2;

      // E(r) = (4*pi/V) * sum_{k!=0} (k/k^2) * exp(-k^2/4a^2) *
      //        Im[ S(k) * exp(i*k.r) ]
      // See class documentation for the derivation.
      std::complex<double> phase = std::exp(i * kv.k.dot(r));
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
