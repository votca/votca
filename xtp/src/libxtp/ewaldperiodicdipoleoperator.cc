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
#include <fstream>
#include <stdexcept>

// Local VOTCA includes
#include "votca/xtp/ewaldperiodicdipoleoperator.h"

namespace votca {
namespace xtp {

EwaldPeriodicDipoleOperator::EwaldPeriodicDipoleOperator(
    EwaldRegistry& registry, const EwaldRealSpaceSum& real_sum,
    const EwaldReciprocalSpaceSum& recip_sum,
    const EwaldShapeCorrection& shape, std::vector<Index> ids,
    double alpha_ewald, double thole_a, bool apply_shape_correction,
    bool apply_self_field_correction,
    bool apply_thole_damping_intramolecular,
    bool apply_realspace_intermolecular_coupling,
    bool apply_reciprocal_coupling)
    : registry_(registry),
      real_sum_(real_sum),
      recip_sum_(recip_sum),
      shape_(shape),
      apply_shape_correction_(apply_shape_correction),
      ids_(std::move(ids)),
      // Real thole_a now, not a placeholder -- AddIntraSegmentCoupling
      // genuinely calls ComputeThole. See this class's own constructor
      // documentation and AddIntraSegmentCoupling's own documentation
      // for why.
      intra_interactor_(alpha_ewald, thole_a),
      self_field_matrix_(recip_sum.SelfFieldMatrix()),
      apply_self_field_correction_(apply_self_field_correction),
      apply_thole_damping_intramolecular_(
          apply_thole_damping_intramolecular),
      apply_realspace_intermolecular_coupling_(
          apply_realspace_intermolecular_coupling),
      apply_reciprocal_coupling_(apply_reciprocal_coupling) {
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
  // RawMultiply(v) includes a constant "leak" -- the permanent field from
  // every registered segment *not* in ids_ (e.g. a fixed external
  // background) -- since EwaldRealSpaceSum::AddFieldAt /
  // EwaldReciprocalSpaceSum::AddFieldAtMany sum over every other
  // registered segment, not just other ids_ members. That leak is
  // independent of v (RawMultiply(0) already contains all of it, since at
  // v=0 every ids_ site's own induced dipole is zero and so contributes
  // nothing to it), so it is captured once here and subtracted in every
  // multiply() call. Without this, multiply(0) != 0, which silently
  // breaks the linearity ConjugateGradient requires -- confirmed the hard
  // way, by a failing test (see test_ewaldperiodicdipoleoperator.cc).
  baseline_ = RawMultiply(Eigen::VectorXd::Zero(size_));
}

std::pair<Index, Index> EwaldPeriodicDipoleOperator::LocateSite(
    Index i) const {
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
  // Cross-site (including two different sites of the SAME segment -- see
  // class documentation): never actually read by DiagonalPreconditioner
  // (which only keeps entries where row()==col()); 0.0 here rather than
  // the real (expensive, for a cross-segment pair) periodic tensor entry.
  return 0.0;
}

Eigen::VectorXd EwaldPeriodicDipoleOperator::RawMultiply(
    const Eigen::VectorXd& v) const {
  std::vector<std::pair<Index, PolarSite*>> targets;
  targets.reserve(std::size_t(size_ / 3));

  for (std::size_t n = 0; n < ids_.size(); ++n) {
    PolarSegment& segment =
        registry_.Get(ids_[n], EwaldChargeState::Neutral);
    Index base = offsets_[n];
    for (Index s = 0; s < segment.size(); ++s) {
      PolarSite& site = segment[s];
      site.setInduced_Dipole(v.segment<3>(base + 3 * s));
      site.Reset();
      targets.push_back({ids_[n], &site});
    }
  }

  // Gated by apply_realspace_intermolecular_coupling_ -- debug-only, see
  // this class's own constructor documentation.
  if (apply_realspace_intermolecular_coupling_) {
    for (const auto& entry : targets) {
      real_sum_.AddFieldAt<Estatic::V>(entry.first, *entry.second,
                                       EwaldChargeState::Neutral);
    }
  }
  // EwaldReciprocalSpaceSum no longer takes a segment id at all -- it
  // never excludes anything (see its own class documentation, and this
  // class's own documentation for why that matters here specifically) --
  // so only the bare site pointers are needed for this call.
  std::vector<PolarSite*> recip_targets;
  recip_targets.reserve(targets.size());
  for (const auto& entry : targets) {
    recip_targets.push_back(entry.second);
  }
  // Gated by apply_reciprocal_coupling_ -- debug-only, see this class's
  // own constructor documentation.
  if (apply_reciprocal_coupling_) {
    recip_sum_.AddFieldAtMany<Estatic::V>(recip_targets,
                                          EwaldChargeState::Neutral);
  }
  // Shape/surface correction: legacy applies this unconditionally every
  // induction iteration (FU12_ShapeField_At_By, using each site's own
  // CURRENT induced dipole -- see this class's own constructor
  // documentation for the fuller account of why this was missing
  // entirely before). EwaldShapeCorrection::TotalDipoleMoment already
  // sums every registered segment's charge + static dipole + induced
  // dipole together, so the v-dependent part here comes entirely from
  // the induced dipoles just set above (every ids_ site's own, via
  // setInduced_Dipole, plus any other registered segment's induced
  // dipole if set from elsewhere) -- the v-independent part (charges,
  // static dipoles) contributes identically at every call, including
  // v=0, so it is captured once in baseline_ exactly like the real- and
  // reciprocal-space terms above, not something this method needs to
  // handle specially. See this method's own apply_shape_correction_
  // documentation (on this class's own constructor) for why this call
  // is conditional -- debug-only, not a genuine design choice.
  if (apply_shape_correction_) {
    for (const auto& entry : targets) {
      shape_.AddFieldAt<Estatic::V>(*entry.second, EwaldChargeState::Neutral);
    }
  }

  Eigen::VectorXd result(size_);
  for (std::size_t n = 0; n < ids_.size(); ++n) {
    const PolarSegment& segment =
        registry_.Get(ids_[n], EwaldChargeState::Neutral);
    Index base = offsets_[n];
    for (Index s = 0; s < segment.size(); ++s) {
      const PolarSite& site = segment[s];
      // site.V() already includes the true (negative) self-field
      // -self_field_matrix_*v_i, via recip_sum_'s own unconditional
      // (nothing-excluded) sum -- adding +self_field_matrix_*v_i here
      // cancels that leak, matching legacy's own atomic self-interaction
      // correction (see self_field_matrix_'s own class documentation).
      // Gated by apply_self_field_correction_ -- debug-only, see this
      // class's own constructor documentation. Written as a plain if,
      // not a ?: -- Eigen's own Product<> (self_field_matrix_ * v...)
      // and ZeroReturnType (Vector3d::Zero()) are different lazy
      // expression-template types with no implicit common type, so a
      // ternary combining them fails to compile; a real, not
      // hypothetical, error caught by an actual build.
      result.segment<3>(base + 3 * s) =
          site.getPInv() * v.segment<3>(base + 3 * s) + site.V();
      if (apply_self_field_correction_) {
        result.segment<3>(base + 3 * s) +=
            self_field_matrix_ * v.segment<3>(base + 3 * s);
      }
    }
  }
  AddIntraSegmentCoupling(v, result);
  return result;
}

void EwaldPeriodicDipoleOperator::AddIntraSegmentCoupling(
    const Eigen::VectorXd& v, Eigen::VectorXd& result) const {
  // See class documentation for the fuller history: an earlier version
  // of this method used eeInteractor::FillTholeInteraction (a real
  // mistake, since PolarRegion is aperiodic and its "Thole-only"
  // convention was never a considered choice about intramolecular
  // coupling specifically). A LATER version corrected that to erfc-only,
  // UNDAMPED by Thole -- matching, it seemed, legacy's own permanent-
  // field intramolecular treatment (which genuinely has no real-space
  // loop at all, only a reciprocal-space erf compensation -- see
  // EwaldRealSpaceInteractor::ApplyIntramolecularStaticCorrection). That
  // second version was ALSO wrong, for the induced case specifically:
  // legacy's own FU12_ERFC_At_By -- the SAME method used for both
  // intermolecular (confirmed at PolarBackground's own line ~954) and
  // intramolecular (confirmed at its own line ~453-454) induced-induced
  // real-space pairs -- applies Thole damping unconditionally via its
  // own l3/l5 mechanism at both call sites, with no special-casing for
  // same-segment pairs at all. The static and induced cases are
  // genuinely different in legacy's own code (no real-space
  // intramolecular loop at all for statics; a real one, Thole-damped,
  // for induction) -- generalizing from one to the other, twice, in two
  // different directions, was the actual mistake both previous versions
  // made. This version applies ComputeThole exactly as ApplyInducedField
  // itself does, rather than assuming l3=l5=1.0 as prior versions did.
  for (std::size_t n = 0; n < ids_.size(); ++n) {
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
        // r_vec points from site_j (source) to site_i (target), matching
        // ApplyInducedField's own r_vec = site2.getPos() - site1.getPos()
        // convention (site1=source, site2=target) with site_j as source,
        // site_i as target.
        const Eigen::Vector3d r_vec = site_i.getPos() - site_j.getPos();
        const double r = r_vec.norm();
        const EwaldRealSpaceInteractor::BFunctions b =
            intra_interactor_.ComputeB(r);
        // Gated by apply_thole_damping_intramolecular_ -- debug-only,
        // see this class's own constructor documentation. When false,
        // l3=l5=1.0 (fully undamped) is used instead of ComputeThole's
        // own result, reverting to this method's own SECOND (also
        // wrong, see this method's own documentation above) prior
        // version's behavior -- kept only to isolate whether the Thole-
        // damping fix itself is a cause of a real, large PCG
        // divergence found once it, the shape-correction fix, and the
        // self-field-correction fix were all added in the same
        // session; the other two have since been ruled out via their
        // own matching debug flags (see EwaldBackground's own
        // apply_shape_correction_to_induced_/
        // apply_self_field_correction_to_induced_ members for that
        // account), leaving this the last untested mechanism added
        // since the point in that same session where the system still
        // converged, if slowly.
        EwaldRealSpaceInteractor::TholeFactors t;
        if (apply_thole_damping_intramolecular_) {
          t = intra_interactor_.ComputeThole(r, site_j, site_i);
        } else {
          t.l3 = 1.0;
          t.l5 = 1.0;
        }
        // Same l3*B1 / l5*B2 combination ApplyInducedField itself uses
        // (see that method's own comment on why there is no separate
        // 3.0* factor here -- B2 already carries it). Not symmetric in
        // general now that Thole damping is site-pair-dependent (unlike
        // the earlier undamped version, where r_vec (x) r_vec and
        // Identity being individually symmetric made block.transpose()
        // == block trivially) -- ComputeThole(r, site_j, site_i) and
        // ComputeThole(r, site_i, site_j) would give the same l3/l5
        // here specifically (same r, same two sites, order-independent
        // eigendamp product), so the block itself is still symmetric in
        // this case, but that's a property of ComputeThole's own
        // symmetric input, not assumed structurally the way it was
        // before.
        Eigen::Matrix3d block = t.l5 * b.B2 * (r_vec * r_vec.transpose()) -
                               t.l3 * b.B1 * Eigen::Matrix3d::Identity();
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

void EwaldPeriodicDipoleOperator::DumpIntraPairThole(
    const std::string& filename) const {
  // See this method's own declaration for what this is for. Mirrors
  // AddIntraSegmentCoupling's own loop exactly (same pair set, same
  // r_vec convention, same ComputeB/ComputeThole calls) -- deliberately
  // NOT reading anything from that method's own tensor-assembly step,
  // so this dump is independent of whatever that step does with these
  // scalars.
  std::ofstream dump(filename);
  dump.precision(15);
  dump << "segment_id,site_i,site_j,r_bohr,B0,B1,B2,l3,l5\n";
  for (std::size_t n = 0; n < ids_.size(); ++n) {
    const PolarSegment& segment =
        registry_.Get(ids_[n], EwaldChargeState::Neutral);
    Index n_sites = segment.size();
    if (n_sites < 2) {
      continue;
    }
    for (Index i = 0; i < n_sites; ++i) {
      const PolarSite& site_i = segment[i];
      for (Index j = i + 1; j < n_sites; ++j) {
        const PolarSite& site_j = segment[j];
        const Eigen::Vector3d r_vec = site_i.getPos() - site_j.getPos();
        const double r = r_vec.norm();
        const EwaldRealSpaceInteractor::BFunctions b =
            intra_interactor_.ComputeB(r);
        const EwaldRealSpaceInteractor::TholeFactors t =
            intra_interactor_.ComputeThole(r, site_j, site_i);
        dump << ids_[n] << "," << i << "," << j << "," << r << "," << b.B0
             << "," << b.B1 << "," << b.B2 << "," << t.l3 << "," << t.l5
             << "\n";
      }
    }
  }
  dump.close();
}

void EwaldPeriodicDipoleOperator::DumpStagedCoupling(
    const Eigen::VectorXd& v, const std::string& filename) const {
  // See this method's own declaration for what this is for and why.
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

  // Stage A: intramolecular only. AddIntraSegmentCoupling writes into a
  // plain Eigen::VectorXd (result), NOT into site.V() -- unlike every
  // other stage below -- so it's captured separately here, into fu_a,
  // rather than read off site.V() like the rest.
  Eigen::VectorXd fu_a = Eigen::VectorXd::Zero(size_);
  AddIntraSegmentCoupling(v, fu_a);

  auto site_v = [&]() {
    Eigen::VectorXd out(size_);
    Index idx = 0;
    for (std::size_t n = 0; n < ids_.size(); ++n) {
      PolarSegment& segment =
          registry_.Get(ids_[n], EwaldChargeState::Neutral);
      for (Index s = 0; s < segment.size(); ++s) {
        out.segment<3>(idx) = segment[s].V();
        idx += 3;
      }
    }
    return out;
  };

  // Stage B: + real-space intermolecular.
  for (const auto& entry : targets) {
    real_sum_.AddFieldAt<Estatic::V>(entry.first, *entry.second,
                                     EwaldChargeState::Neutral);
  }
  const Eigen::VectorXd fu_b = fu_a + site_v();

  // Stage C: + reciprocal.
  std::vector<PolarSite*> recip_targets;
  recip_targets.reserve(targets.size());
  for (const auto& entry : targets) {
    recip_targets.push_back(entry.second);
  }
  recip_sum_.AddFieldAtMany<Estatic::V>(recip_targets,
                                        EwaldChargeState::Neutral);
  const Eigen::VectorXd fu_c = fu_a + site_v();

  // Stage D: + shape.
  for (const auto& entry : targets) {
    shape_.AddFieldAt<Estatic::V>(*entry.second, EwaldChargeState::Neutral);
  }
  const Eigen::VectorXd fu_d = fu_a + site_v();

  // Stage E: + self-field. Added directly (matches RawMultiply's own
  // self_field_matrix_ * v term), not via site.V().
  Eigen::VectorXd fu_e = fu_d;
  Index idx = 0;
  for (std::size_t n = 0; n < ids_.size(); ++n) {
    PolarSegment& segment = registry_.Get(ids_[n], EwaldChargeState::Neutral);
    for (Index s = 0; s < segment.size(); ++s) {
      fu_e.segment<3>(idx) += self_field_matrix_ * v.segment<3>(idx);
      idx += 3;
    }
  }

  std::ofstream dump(filename);
  dump.precision(15);
  dump << "segment_id,site_index,FUa_x_bohr,FUa_y_bohr,FUa_z_bohr,"
          "FUb_x_bohr,FUb_y_bohr,FUb_z_bohr,FUc_x_bohr,FUc_y_bohr,"
          "FUc_z_bohr,FUd_x_bohr,FUd_y_bohr,FUd_z_bohr,FUe_x_bohr,"
          "FUe_y_bohr,FUe_z_bohr\n";
  idx = 0;
  for (std::size_t n = 0; n < ids_.size(); ++n) {
    const PolarSegment& segment =
        registry_.Get(ids_[n], EwaldChargeState::Neutral);
    for (Index s = 0; s < segment.size(); ++s) {
      dump << ids_[n] << "," << s;
      for (const Eigen::VectorXd* stage :
           std::initializer_list<const Eigen::VectorXd*>{&fu_a, &fu_b, &fu_c,
                                                          &fu_d, &fu_e}) {
        dump << "," << (*stage)(idx) << "," << (*stage)(idx + 1) << ","
             << (*stage)(idx + 2);
      }
      dump << "\n";
      idx += 3;
    }
  }
  dump.close();
}

void EwaldPeriodicDipoleOperator::DumpPerPairIntermolecularField(
    Index target_segment_id, const std::string& filename) const {
  PolarSegment& target_segment =
      registry_.Get(target_segment_id, EwaldChargeState::Neutral);
  real_sum_.DumpPerPairFieldAppend(target_segment_id, target_segment[0],
                                   EwaldChargeState::Neutral, filename);
}

void EwaldPeriodicDipoleOperator::DumpPerPairIntermolecularStaticField(
    Index target_segment_id, const std::string& filename) const {
  PolarSegment& target_segment =
      registry_.Get(target_segment_id, EwaldChargeState::Neutral);
  real_sum_.DumpPerPairStaticFieldAppend(target_segment_id, target_segment[0],
                                        EwaldChargeState::Neutral, filename);
}

void EwaldPeriodicDipoleOperator::DumpCachedNeighborList(
    Index target_segment_id, const std::string& filename) const {
  const PolarSegment& target_segment =
      registry_.Get(target_segment_id, EwaldChargeState::Neutral);
  real_sum_.DumpCachedNeighborListAppend(
      target_segment_id, target_segment[0], EwaldChargeState::Neutral,
      filename);
}

}  // namespace xtp
}  // namespace votca
