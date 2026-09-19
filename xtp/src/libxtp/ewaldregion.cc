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
#include <stdexcept>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/ewaldregion.h"
#include "votca/xtp/ewaldrealspaceinteractor.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/qmregion.h"
#include "votca/xtp/staticregion.h"

namespace votca {
namespace xtp {

namespace {
// Use before Initialize() is a caller error, not a missing feature: the
// background has to be loaded before anything can be asked of it.
[[noreturn]] void NotInitialized(const std::string& what) {
  throw std::runtime_error(
      "EwaldRegion::" + what +
      " was called before Initialize(). The background must be loaded "
      "from its checkpoint first.");
}
}  // namespace

void EwaldRegion::Initialize(const tools::Property& prop) {
  checkpoint_file_ = prop.ifExistsReturnElseReturnDefault<std::string>(
      "checkpoint", "ewaldbackground.hdf5");

  CheckpointFile cpf(checkpoint_file_, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader();
  registry_.ReadFromCpt(r);
  CheckpointReader rp = r.openChild("ewald_parameters");
  params_.ReadFromCpt(rp);
  loaded_ = true;

  XTP_LOG(Log::info, log_)
      << TimeStamp() << " Ewald background read from " << checkpoint_file_
      << ": " << size() << " segments, alpha=" << params_.alpha
      << " bohr^-1, k_max=" << params_.k_max
      << " bohr^-1, box volume=" << params_.box.determinant() << " bohr^3"
      << std::flush;
}

void EwaldRegion::Evaluate(std::vector<std::unique_ptr<Region> >& regions) {
  // Nothing polarizes this region, so all this does is record that the
  // other regions had no effect on it -- the energies are all zero by
  // construction (see the Interactwith* overrides). It is kept rather
  // than skipped so the log still shows the region participating.
  ApplyInfluenceOfOtherRegions(regions);
}

Index EwaldRegion::size() const {
  // The whole periodic cell, not a share of the job's segments -- see
  // this class's own documentation, point 3.
  return Index(registry_.AllIds().size());
}

double EwaldRegion::charge() const {
  if (!loaded_) {
    NotInitialized("charge");
  }
  double q = 0.0;
  for (Index id : registry_.AllIds()) {
    if (!registry_.Has(id, EwaldChargeState::Neutral)) {
      continue;
    }
    for (const PolarSite& site :
         registry_.Get(id, EwaldChargeState::Neutral)) {
      q += site.getCharge();
    }
  }
  return q;
}

double EwaldRegion::Etotal() const {
  // Zero, not a throw, and deliberately so. The background's
  // contribution to the ENERGY is genuinely not implemented -- only its
  // field is -- but throwing here aborts the job after the induction has
  // already converged, destroying the very dipoles the run was for.
  //
  // The omission is announced instead: ApplyFieldTo logs a prominent
  // warning on first use, and AppendResult records it in the job's own
  // output so it survives into the result file rather than living only
  // in a log nobody re-reads.
  return 0.0;
}

void EwaldRegion::WriteToCpt(CheckpointWriter& w) const {
  // Only the path is written, not the background itself. The background
  // is large (thousands of segments) and already lives in its own
  // checkpoint; copying it into every per-iteration job checkpoint would
  // multiply that cost for data that cannot change -- this region is
  // frozen by construction.
  w(checkpoint_file_, "checkpoint_file");
}

void EwaldRegion::ReadFromCpt(CheckpointReader& r) {
  // Reloads from the background checkpoint named at write time -- see
  // WriteToCpt for why the background is not stored inline. If that file
  // has moved since, this fails loudly here rather than silently
  // continuing with an empty background.
  r(checkpoint_file_, "checkpoint_file");

  CheckpointFile cpf(checkpoint_file_, CheckpointAccessLevel::READ);
  CheckpointReader br = cpf.getReader();
  registry_.ReadFromCpt(br);
  CheckpointReader bp = br.openChild("ewald_parameters");
  params_.ReadFromCpt(bp);
  loaded_ = true;

  // The lazily-built sums belong to a particular foreground geometry;
  // after a reload there is no guarantee it is the same one, so they are
  // dropped and rebuilt on next use.
  real_sum_.reset();
  recip_sum_.reset();
  shape_.reset();
  interactor_.reset();
  foreground_copies_.clear();
}

void EwaldRegion::WritePDB(csg::PDBWriter&) const {
  // Deliberately a no-op rather than a throw: the background has no
  // job-local geometry worth writing, and a PDB dump should not be able
  // to abort a run.
}

void EwaldRegion::AppendResult(tools::Property& prop) const {
  prop.add("E_ewald_background", "not_implemented");
  prop.add("segments", std::to_string(size()));
  prop.add("checkpoint", checkpoint_file_);
}


namespace {
Eigen::Vector3d Centroid(const PolarSegment& seg) {
  // Unweighted, matching EwaldRealSpaceSum's own convention (and
  // legacy's PolarSeg::CalcPos). A mass-weighted centre would place the
  // foreground copies fractionally off the positions the real-space sum
  // suppresses, and the erf correction would then remove a copy that was
  // never dropped.
  Eigen::Vector3d pos = Eigen::Vector3d::Zero();
  Index n = 0;
  for (const PolarSite& site : seg) {
    pos += site.getPos();
    ++n;
  }
  return (n > 0) ? Eigen::Vector3d(pos / double(n)) : pos;
}
}  // namespace

void EwaldRegion::BuildSums(const std::vector<PolarSegment>& foreground) const {
  foreground_copies_.clear();
  for (const PolarSegment& seg : foreground) {
    if (!registry_.Has(seg.getId(), EwaldChargeState::Neutral)) {
      throw std::runtime_error(
          "EwaldRegion: the polar region contains segment " +
          std::to_string(seg.getId()) +
          ", which is absent from the periodic background. The foreground "
          "must be carved out of the same system the background was "
          "converged on.");
    }
    // Record the position of the NEUTRAL BACKGROUND COPY this foreground
    // segment displaces -- not the foreground's own centroid.
    //
    // They are not the same point. A charged job maps its central
    // segment with that charge state's own geometry, so the foreground
    // is a relaxed cation (or anion) sitting where a neutral molecule
    // used to be. Its centroid is therefore displaced from the
    // background copy's by far more than the 1e-4 bohr tolerance
    // EwaldRealSpaceSum uses to recognise a copy -- that tolerance is
    // there to absorb round-off, not geometry relaxation.
    //
    // Recording the foreground's own centroid therefore meant the
    // real-space sum found NO match and silently left the neutral copy
    // in the background, while the reciprocal and shape exclusions
    // (which match by address) removed it correctly. The copy was
    // counted in the erfc half and not the erf half: a mismatch that
    // only appears for a charged foreground, and that breaks
    // alpha-independence because erfc and erf shift weight with alpha.
    //
    // Snapping to the nearest lattice image rather than requiring an
    // exact hit also makes this robust to any future geometry that
    // differs between charge states.
    const PolarSegment& bg_copy =
        registry_.Get(seg.getId(), EwaldChargeState::Neutral);
    const Eigen::Vector3d bg_centroid = Centroid(bg_copy);
    const Eigen::Vector3d delta = Centroid(seg) - bg_centroid;
    const Eigen::Vector3d fractional = params_.box.inverse() * delta;
    const Eigen::Vector3d image = params_.box * fractional.array().round().matrix();
    const Eigen::Vector3d residual = delta - image;

    // A residual much larger than a molecule means the foreground
    // segment is not where the background thinks that id lives, which
    // is a mapping error rather than a relaxation. Loud, because the old
    // behaviour for exactly this case was to carry on silently.
    constexpr double kResidualWarn = 5.0;  // bohr
    if (residual.norm() > kResidualWarn) {
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " WARNING: foreground segment " << seg.getId()
          << " sits " << residual.norm()
          << " bohr from the nearest periodic image of its background "
             "copy. That is too far to be geometry relaxation, so the "
             "copy this code is about to suppress may not be the one the "
             "foreground actually displaces."
          << std::flush;
    }

    foreground_copies_.push_back({seg.getId(), bg_centroid + image});
  }

  real_sum_ = std::make_unique<EwaldRealSpaceSum>(
      params_.box, registry_, params_.alpha, params_.thole_a, params_.r_min,
      params_.field_tol, 0.945, 15, params_.screening_factor,
      foreground_copies_);
  recip_sum_ = std::make_unique<EwaldReciprocalSpaceSum>(
      params_.box, registry_, params_.alpha, params_.k_max);
  shape_ = std::make_unique<EwaldShapeCorrection>(
      params_.box.determinant(), registry_, params_.shape);
  interactor_ = std::make_unique<EwaldRealSpaceInteractor>(params_.alpha,
                                                           params_.thole_a);
}

double EwaldRegion::ApplyFieldTo(std::vector<PolarSegment>& foreground) const {
  if (!loaded_) {
    NotInitialized("ApplyFieldTo");
  }
  if (!real_sum_) {
    BuildSums(foreground);
  }

  // Flat target list: the reciprocal sum takes them in one batch, which
  // is what makes its structure factors worth computing once.
  std::vector<PolarSite*> targets;
  std::vector<Index> target_segment_ids;
  for (PolarSegment& seg : foreground) {
    for (PolarSite& site : seg) {
      targets.push_back(&site);
      target_segment_ids.push_back(seg.getId());
    }
  }

  // SIGN CONVENTION. The Ewald code and the region framework store V
  // with OPPOSITE signs, and this is the boundary between them:
  //
  //   EwaldBackground      builds its solver's rhs as b = +V
  //   PolarRegion          builds its solver's rhs as b = -(V + V_noE)
  //
  // Both are internally consistent; neither is wrong on its own. But the
  // sums below write in the Ewald convention, and the polar region is
  // about to negate whatever it finds -- so handing it the field as-is
  // drives the induction BACKWARDS. Measured on a neutral MM/MM job,
  // where the foreground should reproduce the background it was carved
  // from: deviations of 100-700%, flat with distance rather than growing
  // outward, which is what first exposed this.
  //
  // So this region's own contribution is negated before it is handed
  // over. Only the contribution is touched: whatever other regions have
  // already accumulated is preserved, which is why the prior values are
  // recorded rather than assuming V starts at zero.
  std::vector<Eigen::Vector3d> v_before;
  v_before.reserve(targets.size());
  for (const PolarSite* target : targets) {
    v_before.push_back(target->V());
  }

  // (1) real space, foreground copies suppressed
  const Index n_targets = Index(targets.size());
#pragma omp parallel for schedule(dynamic, 16)
  for (Index i = 0; i < n_targets; ++i) {
    real_sum_->AddFieldAt<Estatic::V>(target_segment_ids[std::size_t(i)],
                                      *targets[std::size_t(i)],
                                      EwaldChargeState::Neutral);
  }

  // Permanent-permanent (Q-Q) energy with the background, over exactly
  // the neighbour set the field above used -- same suppression, same
  // cull. Only the PERMANENT part: PolarRegion accounts for the induced
  // contribution itself, as E_polar_ext = sum of mu_ind . V, computed
  // from the very field being delivered here. Adding it again would
  // double-count.
  //
  // Runs after the field loop, not inside it, because the neighbour
  // cache must exist first and CalcStaticEnergyAt deliberately refuses
  // to build one.
  double e_real = 0.0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : e_real)
  for (Index i = 0; i < n_targets; ++i) {
    e_real += real_sum_->CalcStaticEnergyAt(*targets[std::size_t(i)],
                                            EwaldChargeState::Neutral);
  }
  double energy = e_real;

  // (1b) [foreground permanent] x [background INDUCED], real space.
  //      The corner of the permanent/induced product that nothing else
  //      covers: PolarRegion's E_polar_ext contracts the FOREGROUND's
  //      induced dipoles against the delivered field, which gives
  //      [fg induced] x [bg anything]; e_real above gives
  //      [fg permanent] x [bg permanent]. Without this term a job's
  //      induced energy comes out short -- measurably so, by a factor
  //      near two on a neutral foreground, against legacy's _pu channel.
  //
  //      Same neighbour set, same cache, same suppression as e_real.
  double e_real_pu = 0.0;
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : e_real_pu)
  for (Index i = 0; i < n_targets; ++i) {
    e_real_pu += real_sum_->CalcInducedSourceEnergyAt(
        *targets[std::size_t(i)], EwaldChargeState::Neutral);
  }
  energy += e_real_pu;

  // (2) reciprocal space, over the full periodic density
  recip_sum_->AddFieldAtMany<Estatic::V>(targets, EwaldChargeState::Neutral);

  // (3) shape/surface. Target-independent, so evaluated once.
  {
    PolarSite probe(-1, "X", Eigen::Vector3d::Zero());
    shape_->AddFieldAt<Estatic::V>(probe, EwaldChargeState::Neutral);
    const Eigen::Vector3d shape_field = probe.V();
    for (PolarSite* target : targets) {
      target->V() += shape_field;
    }
  }

  // (4) remove the erf-screened field of the neutral foreground copies
  //     that step (2) necessarily put back. See ApplyFieldTo's own
  //     declaration for why these use the background's own multipoles
  //     and dipoles rather than the job's charge state.
#pragma omp parallel for schedule(dynamic, 16)
  for (Index i = 0; i < n_targets; ++i) {
    PolarSite& target = *targets[std::size_t(i)];
    for (const auto& copy : foreground_copies_) {
      const PolarSegment& bg =
          registry_.Get(copy.first, EwaldChargeState::Neutral);
      const Eigen::Vector3d shift = copy.second - Centroid(bg);
      for (const PolarSite& source : bg) {
        interactor_->ApplyErfStaticFieldCorrection<PolarSite, Estatic::V>(
            source, target, shift);
        interactor_->ApplyErfInducedFieldCorrection<Estatic::V>(source, target,
                                                                 shift);
      }
    }
  }
  // (5) reciprocal-space and shape Q-Q energy between the foreground and
  //     the rest of the cell, computed PER FOREGROUND SEGMENT.
  //
  //     Per-segment, not once for the whole foreground, because the
  //     real-space sum's exclusions are per-target and cannot be
  //     expressed in a single structure factor. For a target in
  //     foreground segment i, real space omits
  //
  //       - segment i at EVERY translation (intermolecular scope, see
  //         EwaldRealSpaceSum::AddFieldAt), and
  //       - the COINCIDENT copy of every other foreground segment j,
  //         while keeping segment j's other periodic images, which are
  //         ordinary background molecules.
  //
  //     A structure factor cannot drop a single image: removing a
  //     segment from S_bg removes it from all of them. An earlier
  //     version passed one exclusion list for the whole foreground and
  //     so removed every foreground segment's images too. That is a
  //     different physical system -- a periodic lattice of vacancies
  //     rather than one carved-out cavity -- and it showed up as a
  //     failure of alpha-independence: on a charged 18-segment
  //     foreground the permanent energy moved from -3.1e-4 through
  //     -8.5e-4 to +7.1e-4 eV across alpha = 1.5, 2.0, 3.0 nm^-1, while
  //     the neutral job held to 6e-8. A NEUTRAL foreground cannot see
  //     this, which is why every alpha-independence test in the suite
  //     missed it.
  //
  //     So: exclude only segment i from S_bg (leaving every other
  //     segment, images included), then subtract the erf-screened
  //     energy of the coincident copies of the other foreground
  //     segments. recip + shape together are the erf interaction over
  //     all images -- the shape term IS the k=0 limit the reciprocal
  //     sum omits, which is why its Q0*TrQ2 pieces are not optional --
  //     so subtracting a pair's erf energy removes that pair exactly.
  {
    double e_recip = 0.0;
    double e_shape = 0.0;
    double e_erf = 0.0;

    for (const auto& copy : foreground_copies_) {
      const PolarSegment& bg_self =
          registry_.Get(copy.first, EwaldChargeState::Neutral);
      const Eigen::Vector3d shift_self = copy.second - Centroid(bg_self);

      // This segment's foreground sites, in the job's charge state, at
      // the positions they occupy.
      std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites;
      for (std::size_t n = 0; n < targets.size(); ++n) {
        if (target_segment_ids[n] == copy.first) {
          fg_sites.push_back({targets[n], targets[n]->getPos()});
        }
      }
      if (fg_sites.empty()) {
        continue;
      }

      // Held out of S_bg: this segment only, at every image.
      std::vector<const PolarSite*> excl_self;
      for (const PolarSite& site : bg_self) {
        excl_self.push_back(&site);
      }

      e_recip += recip_sum_->CalcStaticEnergyBetween(fg_sites, excl_self,
                                                     EwaldChargeState::Neutral);
      e_shape += shape_->CalcStaticEnergyBetween(fg_sites, excl_self,
                                                 EwaldChargeState::Neutral);

      // Remove the coincident copies of the OTHER foreground segments,
      // which real space already suppressed. Their images stay.
      for (const auto& other : foreground_copies_) {
        if (other.first == copy.first) {
          continue;
        }
        const PolarSegment& bg_other =
            registry_.Get(other.first, EwaldChargeState::Neutral);
        const Eigen::Vector3d shift_other =
            other.second - Centroid(bg_other);
        for (const PolarSite& source : bg_other) {
          for (const auto& entry : fg_sites) {
            e_erf += interactor_->CalcErfStaticEnergy<PolarSite, PolarSite>(
                source, *entry.first, shift_other);
          }
        }
      }
    }
    energy += e_recip + e_shape - e_erf;

    // (5b), (6b) reciprocal and shape partners of the induced-source
    //      term added at (1b), with exactly the same per-segment
    //      exclusion structure and for exactly the same reason.
    double e_recip_pu = 0.0;
    double e_shape_pu = 0.0;
    double e_erf_pu = 0.0;
    for (const auto& copy : foreground_copies_) {
      const PolarSegment& bg_self =
          registry_.Get(copy.first, EwaldChargeState::Neutral);

      std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites_pu;
      for (std::size_t n = 0; n < targets.size(); ++n) {
        if (target_segment_ids[n] == copy.first) {
          fg_sites_pu.push_back({targets[n], targets[n]->getPos()});
        }
      }
      if (fg_sites_pu.empty()) {
        continue;
      }

      std::vector<const PolarSite*> excl_self;
      for (const PolarSite& site : bg_self) {
        excl_self.push_back(&site);
      }

      e_recip_pu += recip_sum_->CalcInducedSourceEnergyBetween(
          fg_sites_pu, excl_self, EwaldChargeState::Neutral);
      e_shape_pu += shape_->CalcInducedSourceEnergyBetween(
          fg_sites_pu, excl_self, EwaldChargeState::Neutral);

      for (const auto& other : foreground_copies_) {
        if (other.first == copy.first) {
          continue;
        }
        const PolarSegment& bg_other =
            registry_.Get(other.first, EwaldChargeState::Neutral);
        const Eigen::Vector3d shift_other =
            other.second - Centroid(bg_other);
        for (const PolarSite& source : bg_other) {
          for (const auto& entry : fg_sites_pu) {
            e_erf_pu += interactor_->CalcErfInducedSourceEnergy(
                source, *entry.first, shift_other);
          }
        }
      }
    }
    energy += e_recip_pu + e_shape_pu - e_erf_pu;

    // Term-by-term report, for comparison against the legacy `ewald`
    // job calculator's own terms_o block. The correspondence is NOT
    // one-to-one: legacy adds the full background in its K term and
    // then subtracts its C term (the erf-screened fgC<->fgN
    // interaction), where this code excludes the foreground copies up
    // front. So
    //
    //   e_real   <->  legacy R_pp
    //   e_recip  <->  legacy (K_pp - C_pp)
    //   e_shape  <->  legacy J_pp
    //
    // and the TOTAL is what must agree. Printed in eV, because that is
    // the unit legacy reports and the unit the job XML carries.
    const double h2ev = tools::conv::hrt2ev;
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " Ewald energy [eV], permanent x permanent: real = "
        << e_real * h2ev << "  recip = " << e_recip * h2ev
        << "  shape = " << e_shape * h2ev << "  erf = " << e_erf * h2ev
        << std::flush;
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " Ewald energy [eV], fg permanent x bg induced: real = "
        << e_real_pu * h2ev << "  recip = " << e_recip_pu * h2ev
        << "  shape = " << e_shape_pu * h2ev << "  erf = " << e_erf_pu * h2ev
        << std::flush;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Ewald energy [eV], total = " << energy * h2ev
        << std::flush;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Ewald splitting: alpha = " << params_.alpha
        << " 1/bohr (" << params_.alpha * 18.8972612 << " 1/nm)"
        << ", k_max = " << params_.k_max << " 1/bohr ("
        << params_.k_max * 18.8972612 << " 1/nm)"
        << ", r_min = " << params_.r_min << " bohr"
        << ", V = " << params_.box.determinant() << " bohr^3"
        << ", thole_a = " << params_.thole_a << std::flush;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Ewald counts: foreground sites = "
        << targets.size() << ", foreground segments = "
        << foreground_copies_.size() << ", registered segments = "
        << registry_.AllIds().size() << std::flush;

    // Suppression audit. Every target should have had one copy of every
    // OTHER foreground segment suppressed (its own segment is excluded
    // earlier, by the intermolecular-scope rule). A shortfall means some
    // foreground copy was left in the background, which is silent
    // otherwise and was exactly the failure this run is checking for.
    const EwaldRealSpaceSum::NeighborStats stats =
        real_sum_->GetNeighborStats();
    const Index expected =
        Index(targets.size()) * (Index(foreground_copies_.size()) - 1);
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Ewald foreground suppression: " << stats.foreground
        << " of " << expected << " expected"
        << ((stats.foreground == expected) ? " (ok)" : "  <-- MISMATCH")
        << std::flush;
  }

  // Convert this region's contribution into the convention the polar
  // region expects -- see the note where v_before is captured.
  for (std::size_t i = 0; i < targets.size(); ++i) {
    const Eigen::Vector3d ewald_contribution = targets[i]->V() - v_before[i];
    targets[i]->V() = v_before[i] - ewald_contribution;
  }

  if (!warned_no_energy_) {
    warned_no_energy_ = true;
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " NOTE: the Ewald background's permanent (Q-Q) interaction "
           "energy with this foreground is complete -- real space, "
           "reciprocal space and the shape/surface term. Its INDUCED "
           "contribution is not missing either; the polar region "
           "accounts for it itself, from the field delivered here. What "
           "is deliberately absent is the background's own internal "
           "energy, which is a constant of the background and cancels in "
           "any difference taken between charge states of the same "
           "segment."
        << std::flush;
  }
  return energy;
}

}  // namespace xtp
}  // namespace votca
