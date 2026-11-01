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
#include <array>
#include <chrono>
#include <fstream>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/ewaldperiodicdipoleoperator.h"

namespace votca {
namespace xtp {

EwaldPeriodicDipoleOperator::EwaldPeriodicDipoleOperator(
    EwaldRegistry& registry, const EwaldRealSpaceSum& real_sum,
    const EwaldReciprocalSpaceSum& recip_sum, const EwaldShapeCorrection& shape,
    std::vector<Index> ids, double alpha_ewald, double thole_a)
    : registry_(registry),
      real_sum_(real_sum),
      recip_sum_(recip_sum),
      shape_(shape),
      ids_(std::move(ids)),
      intra_interactor_(alpha_ewald, thole_a),
      self_field_matrix_(recip_sum.SelfFieldMatrix()) {
  offsets_.reserve(ids_.size() + 1);
  offsets_.push_back(0);
  for (Index id : ids_) {
    if (!registry_.Has(id, EwaldChargeState::Neutral)) {
      throw std::runtime_error(
          "EwaldPeriodicDipoleOperator: segment id not registered at "
          "EwaldChargeState::Neutral");
    }
    Index n_sites = registry_.Get(id, EwaldChargeState::Neutral).size();
    offsets_.push_back(offsets_.back() + 3 * n_sites);
  }
  size_ = offsets_.back();
  // Build every target's real-space neighbour list up front. RawMultiply
  // parallelizes over targets and the cache is written only while a list
  // is being built, so this pass makes that loop race-free. See
  // EwaldRealSpaceSum::PrepareNeighborCache.
  {
    std::vector<std::pair<Index, PolarSite*>> cache_targets;
    cache_targets.reserve(std::size_t(size_ / 3));
    for (std::size_t n = 0; n < ids_.size(); ++n) {
      PolarSegment& segment = registry_.Get(ids_[n], EwaldChargeState::Neutral);
      for (Index s = 0; s < segment.size(); ++s) {
        cache_targets.push_back({ids_[n], &segment[s]});
      }
    }
    real_sum_.PrepareNeighborCache(cache_targets, EwaldChargeState::Neutral);
  }

  // RawMultiply(v) contains a v-independent field from registered segments
  // not in ids_; it is captured once here and subtracted in multiply() so
  // that multiply(0) = 0.
  baseline_ = RawMultiply(Eigen::VectorXd::Zero(size_));
}

std::pair<Index, Index> EwaldPeriodicDipoleOperator::LocateSite(Index i) const {
  // upper_bound finds the first offset strictly greater than i; the
  // segment i belongs to is the one just before that.
  auto it = std::upper_bound(offsets_.begin(), offsets_.end(), i);
  Index seg_idx = Index(it - offsets_.begin()) - 1;
  Index site_idx = (i - offsets_[std::size_t(seg_idx)]) / 3;
  return {seg_idx, site_idx};
}

double EwaldPeriodicDipoleOperator::operator()(Index i, Index j) const {
  auto [seg1, site1] = LocateSite(i);
  auto [seg2, site2] = LocateSite(j);
  Index xyz1 = (i - offsets_[std::size_t(seg1)]) % 3;
  Index xyz2 = (j - offsets_[std::size_t(seg2)]) % 3;

  if (seg1 == seg2 && site1 == site2) {
    const PolarSegment& segment =
        registry_.Get(ids_[std::size_t(seg1)], EwaldChargeState::Neutral);
    return segment[site1].getPInv()(xyz1, xyz2);
  }
  // Cross-site entries (including distinct sites of one segment) are not
  // used by the diagonal preconditioner and are not computed.
  return 0.0;
}

Eigen::VectorXd EwaldPeriodicDipoleOperator::RawMultiply(
    const Eigen::VectorXd& v) const {
  // Phase timers (see RawMultiplyTimings); negligible overhead compared
  // with the pair evaluations of one matvec.
  using clock = std::chrono::steady_clock;
  auto clk = []() { return clock::now(); };
  auto secs = [](clock::time_point a, clock::time_point b) {
    return std::chrono::duration<double>(b - a).count();
  };
  ++timings_.n_calls;
  auto t_phase = clk();

  std::vector<std::pair<Index, PolarSite*>> targets;
  targets.reserve(std::size_t(size_ / 3));

  for (std::size_t n = 0; n < ids_.size(); ++n) {
    PolarSegment& segment = registry_.Get(ids_[n], EwaldChargeState::Neutral);
    Index base = offsets_[n];
    for (Index s = 0; s < segment.size(); ++s) {
      PolarSite& site = segment[s];
      site.setInduced_Dipole(v.segment<3>(base + 3 * s));
      site.Reset();
      targets.push_back({ids_[n], &site});
    }
  }

  timings_.setup += secs(t_phase, clk());
  t_phase = clk();

  // Parallel over targets: each iteration writes only its own target, and
  // the neighbour cache was fully built in the constructor, so AddFieldAt
  // does not modify real_sum_.
  //
  // include_static = false: the permanent-multipole field is independent
  // of v, so it is identical here and in baseline_ = RawMultiply(0), and
  // cancels exactly in multiply().
  {
    const Index n_targets = Index(targets.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (Index i = 0; i < n_targets; ++i) {
      real_sum_.AddFieldAt<Estatic::V>(targets[std::size_t(i)].first,
                                       *targets[std::size_t(i)].second,
                                       EwaldChargeState::Neutral, false);
    }
  }
  timings_.real_space += secs(t_phase, clk());
  t_phase = clk();

  // The reciprocal-space sum excludes nothing, so it needs no segment ids.
  std::vector<PolarSite*> recip_targets;
  recip_targets.reserve(targets.size());
  for (const auto& entry : targets) {
    recip_targets.push_back(entry.second);
  }
  recip_sum_.AddFieldAtMany<Estatic::V>(recip_targets,
                                        EwaldChargeState::Neutral);
  timings_.reciprocal += secs(t_phase, clk());
  t_phase = clk();

  // Shape/surface correction, applied every induction iteration as in
  // legacy (FU12_ShapeField_At_By). EwaldShapeCorrection::TotalDipoleMoment
  // includes charges, static and induced dipoles; the v-independent part
  // cancels against baseline_.
  //
  // The field -(4*pi/3V)*M is the same at every site, so it is computed
  // once with a probe site rather than per target (which would be
  // O(N^2); 0.45 s per matvec for 5000 sites).
  {
    PolarSite shape_probe(-1, "X", Eigen::Vector3d::Zero());
    shape_.AddFieldAt<Estatic::V>(shape_probe, EwaldChargeState::Neutral);
    const Eigen::Vector3d shape_field = shape_probe.V();
    for (const auto& entry : targets) {
      entry.second->V() += shape_field;
    }
  }

  timings_.shape += secs(t_phase, clk());
  t_phase = clk();

  Eigen::VectorXd result(size_);
  for (std::size_t n = 0; n < ids_.size(); ++n) {
    const PolarSegment& segment =
        registry_.Get(ids_[n], EwaldChargeState::Neutral);
    Index base = offsets_[n];
    for (Index s = 0; s < segment.size(); ++s) {
      const PolarSite& site = segment[s];
      // The system is mu = P*(F_perm + FU(mu)), i.e. (P^-1 - C)*mu = F_perm
      // with C*mu = FU the induced-dipole field, so the coupling enters
      // with a minus sign. site.V() contains the spurious reciprocal-space
      // self-field -M*v_i; adding M*v_i to the field (hence -M*v_i here)
      // removes it, as legacy's self-interaction correction does.
      result.segment<3>(base + 3 * s) =
          site.getPInv() * v.segment<3>(base + 3 * s) - site.V() -
          self_field_matrix_ * v.segment<3>(base + 3 * s);
    }
  }
  // Intramolecular coupling is part of C and is subtracted as well;
  // AddIntraSegmentCoupling itself returns the field, not the operator
  // block.
  timings_.assemble += secs(t_phase, clk());
  t_phase = clk();

  Eigen::VectorXd intra = Eigen::VectorXd::Zero(size_);
  AddIntraSegmentCoupling(v, intra);
  result -= intra;

  timings_.intra += secs(t_phase, clk());
  return result;
}

void EwaldPeriodicDipoleOperator::AddIntraSegmentCoupling(
    const Eigen::VectorXd& v, Eigen::VectorXd& result) const {
  // Induced-induced intramolecular pairs are Thole damped, as in legacy
  // FU12_ERFC_At_By (which damps intra- and intermolecular pairs alike).
  // This differs from the permanent field, where legacy has no real-space
  // intramolecular term (see ApplyErfStaticFieldCorrection).
  // Each segment writes only its own slice of result.
#pragma omp parallel for schedule(dynamic, 4)
  for (Index nn = 0; nn < Index(ids_.size()); ++nn) {
    const std::size_t n = std::size_t(nn);
    const PolarSegment& segment =
        registry_.Get(ids_[n], EwaldChargeState::Neutral);
    Index n_sites = segment.size();
    if (n_sites < 2) {
      continue;
    }
    Index base = offsets_[n];
    for (Index i = 0; i < n_sites; ++i) {
      const PolarSite& site_i = segment[i];
      for (Index j = i + 1; j < n_sites; ++j) {
        const PolarSite& site_j = segment[j];
        // r_vec points from site_j (source) to site_i (target), as in
        // ApplyInducedField (r_vec = site2 - site1, site1 = source).
        const Eigen::Vector3d r_vec = site_i.getPos() - site_j.getPos();
        const double r = r_vec.norm();
        const EwaldRealSpaceInteractor::BFunctions b =
            intra_interactor_.ComputeB(r);
        const EwaldRealSpaceInteractor::BFunctions berf =
            intra_interactor_.ComputeErfB(r);
        const EwaldRealSpaceInteractor::TholeFactors t =
            intra_interactor_.ComputeThole(r, site_j, site_i);
        // Same l3*B1 / l5*B2 combination as ApplyInducedField (B2
        // already carries the factor 3). l3/l5 are symmetric in the two
        // sites, so block is symmetric.
        //
        // The pair is excluded from real_sum_ but contained, undamped, in
        // recip_sum_ (B_erf). The total must be the alpha-independent
        // l*B_bare, so this term supplies
        //     l*B_bare - B_erf  ==  l*B_erfc + (l-1)*B_erf.
        // The second form avoids cancellation of large bare terms at
        // small r. Using only l*B_erfc caused a 0.43% alpha dependence
        // of the converged dipoles (alpha = 1.5 ... 3.0 1/nm).
        const double c3 = t.l3 * b.B1 + (t.l3 - 1.0) * berf.B1;
        const double c5 = t.l5 * b.B2 + (t.l5 - 1.0) * berf.B2;
        Eigen::Matrix3d block =
            c5 * (r_vec * r_vec.transpose()) - c3 * Eigen::Matrix3d::Identity();
        result.segment<3>(base + 3 * i) += block * v.segment<3>(base + 3 * j);
        result.segment<3>(base + 3 * j) +=
            block.transpose() * v.segment<3>(base + 3 * i);
      }
    }
  }
}

Eigen::VectorXd EwaldPeriodicDipoleOperator::multiply(
    const Eigen::VectorXd& v) const {
  return RawMultiply(v) - baseline_;
}

}  // namespace xtp
}  // namespace votca
