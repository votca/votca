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
#include <sstream>
#include <stdexcept>

// VOTCA includes
#include <votca/tools/constants.h>

// Local VOTCA includes
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/ewaldrealspaceinteractor.h"
#include "votca/xtp/ewaldregion.h"
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

  XTP_LOG(Log::info, log_) << TimeStamp() << " Ewald background read from "
                           << checkpoint_file_ << ": " << size()
                           << " segments, alpha=" << params_.alpha
                           << " bohr^-1, k_max=" << params_.k_max
                           << " bohr^-1, box volume="
                           << params_.box.determinant() << " bohr^3"
                           << std::flush;
}

void EwaldRegion::Evaluate(std::vector<std::unique_ptr<Region>>& regions) {
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
    for (const PolarSite& site : registry_.Get(id, EwaldChargeState::Neutral)) {
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
  built_foreground_.clear();
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

void EwaldRegion::RegisterForeground(
    const std::vector<std::pair<Index, Eigen::Vector3d>>& foreground) {
  // The union is disjoint by construction -- PartitionRegions marks each
  // segment as it assigns it -- but that guarantee lives far from the
  // code relying on it, and a repeated id would suppress the same
  // background copy twice.
  for (std::size_t i = 0; i < foreground.size(); ++i) {
    for (std::size_t j = i + 1; j < foreground.size(); ++j) {
      if (foreground[i].first == foreground[j].first) {
        std::stringstream message;
        message << "EwaldRegion::RegisterForeground: segment "
                << foreground[i].first
                << " was declared twice. The foreground is the disjoint "
                   "union of the regions that own segments.";
        throw std::runtime_error(message.str());
      }
    }
  }
  registered_foreground_ = foreground;
  // Anything built before the declaration was built from the wrong
  // foreground.
  real_sum_.reset();
  recip_sum_.reset();
  shape_.reset();
  interactor_.reset();
  foreground_copies_.clear();
  built_foreground_.clear();
}

void EwaldRegion::CheckForegroundIsSubset(
    const std::vector<PolarSegment>& foreground) const {
  // The same 1e-4 bohr EwaldRealSpaceSum uses to recognise a copy
  // (kForegroundMatchTol, private there, so restated rather than
  // shared): it absorbs round-off, nothing larger. A foreground's
  // positions do not move between calls within a job, so any difference
  // above this is a different segment, not drift.
  constexpr double kTol = 1e-4;
  for (const PolarSegment& seg : foreground) {
    const Eigen::Vector3d centroid = Centroid(seg);
    bool found = false;
    for (const auto& entry : built_foreground_) {
      if (entry.first == seg.getId() &&
          (centroid - entry.second).norm() <= kTol) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::stringstream message;
      message << "EwaldRegion: asked about segment " << seg.getId()
              << ", which is not part of the foreground these sums were "
                 "built for. The suppression list and the real-space "
                 "neighbour cache belong to that foreground, so answering "
                 "would drop the wrong background copies and subtract erf "
                 "corrections for copies that were never dropped -- a wrong "
                 "energy with no symptom. JobTopology declares the whole "
                 "foreground with RegisterForeground before any region is "
                 "evaluated; if this fires, that declaration is missing or "
                 "incomplete.";
      throw std::runtime_error(message.str());
    }
  }
}

void EwaldRegion::BuildSums(const std::vector<PolarSegment>& fallback) const {
  std::vector<std::pair<Index, Eigen::Vector3d>> source =
      registered_foreground_;
  if (source.empty()) {
    for (const PolarSegment& seg : fallback) {
      source.push_back({seg.getId(), Centroid(seg)});
    }
  }

  foreground_copies_.clear();
  built_foreground_ = source;
  for (const auto& entry : source) {
    const Index seg_id = entry.first;
    if (!registry_.Has(seg_id, EwaldChargeState::Neutral)) {
      throw std::runtime_error(
          "EwaldRegion: the foreground contains segment " +
          std::to_string(seg_id) +
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
        registry_.Get(seg_id, EwaldChargeState::Neutral);
    const Eigen::Vector3d bg_centroid = Centroid(bg_copy);
    const Eigen::Vector3d delta = entry.second - bg_centroid;
    const Eigen::Vector3d fractional = params_.box.inverse() * delta;
    const Eigen::Vector3d image =
        params_.box * fractional.array().round().matrix();
    const Eigen::Vector3d residual = delta - image;

    // A residual much larger than a molecule means the foreground
    // segment is not where the background thinks that id lives, which
    // is a mapping error rather than a relaxation. Loud, because the old
    // behaviour for exactly this case was to carry on silently.
    constexpr double kResidualWarn = 5.0;  // bohr
    if (residual.norm() > kResidualWarn) {
      XTP_LOG(Log::error, log_)
          << TimeStamp() << " WARNING: foreground segment " << seg_id
          << " sits " << residual.norm()
          << " bohr from the nearest periodic image of its background "
             "copy. That is too far to be geometry relaxation, so the "
             "copy this code is about to suppress may not be the one the "
             "foreground actually displaces."
          << std::flush;
    }

    foreground_copies_.push_back({seg_id, bg_centroid + image});
  }

  real_sum_ = std::make_unique<EwaldRealSpaceSum>(
      params_.box, registry_, params_.alpha, params_.thole_a, params_.r_min,
      params_.field_tol, 0.945, 15, params_.screening_factor,
      foreground_copies_);
  recip_sum_ = std::make_unique<EwaldReciprocalSpaceSum>(
      params_.box, registry_, params_.alpha, params_.k_max);
  shape_ = std::make_unique<EwaldShapeCorrection>(params_.box.determinant(),
                                                  registry_, params_.shape);
  interactor_ = std::make_unique<EwaldRealSpaceInteractor>(params_.alpha,
                                                           params_.thole_a);
}

Eigen::VectorXd EwaldRegion::PotentialAt(
    const std::vector<Eigen::Vector3d>& points) const {
  if (!loaded_) {
    NotInitialized("PotentialAt");
  }
  if (registered_foreground_.empty()) {
    throw std::runtime_error(
        "EwaldRegion::PotentialAt: no foreground has been declared. "
        "ApplyFieldTo can fall back on the segments it is handed, but a "
        "list of points carries no segment identity, so the copies to "
        "suppress cannot be inferred. JobTopology declares the foreground "
        "with RegisterForeground before any region is evaluated.");
  }
  if (!real_sum_) {
    BuildSums(std::vector<PolarSegment>());
  }

  // A foreground segment's id. It decides only whether the
  // zero-translation self-pair is skipped, and that skip is disabled for
  // sources that have a foreground copy -- which is what a point that is
  // not a site of its own wants.
  const Index probe_segment_id = built_foreground_.front().first;

  Eigen::VectorXd phi = real_sum_->PotentialAtMany(probe_segment_id, points,
                                                   EwaldChargeState::Neutral);
  phi += recip_sum_->PotentialAtMany(points, EwaldChargeState::Neutral);

  // Shape and the erf removal share a unit probe per point. Neither
  // walks a neighbour list, so both are cheap enough to evaluate through
  // the existing energy routines rather than re-deriving them here --
  // which also keeps them in the same gauge by construction.
  const Index n_points = Index(points.size());
  const std::vector<const PolarSite*> no_exclusions;
#pragma omp parallel for schedule(static)
  for (Index p = 0; p < n_points; ++p) {
    PolarSite probe(0, "H", points[std::size_t(p)]);
    probe.Reset();
    probe.setCharge(1.0);
    probe.setStaticDipole(Eigen::Vector3d::Zero());
    probe.setInduced_Dipole(Eigen::Vector3d::Zero());
    const std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> one{
        {&probe, points[std::size_t(p)]}};

    double extra = shape_->CalcStaticEnergyBetween(one, no_exclusions,
                                                   EwaldChargeState::Neutral) +
                   shape_->CalcInducedSourceEnergyBetween(
                       one, no_exclusions, EwaldChargeState::Neutral);

    // Remove the erf-screened half of the neutral foreground copies that
    // the reciprocal sum necessarily put back -- the same copies, with
    // the same shift, as step (4) of ApplyFieldTo.
    for (const auto& copy : foreground_copies_) {
      const PolarSegment& bg =
          registry_.Get(copy.first, EwaldChargeState::Neutral);
      const Eigen::Vector3d shift = copy.second - Centroid(bg);
      for (const PolarSite& source : bg) {
        extra -= interactor_->CalcErfStaticEnergy<PolarSite, PolarSite>(
            source, probe, shift);
        extra -= interactor_->CalcErfInducedSourceEnergy(source, probe, shift);
      }
    }
    phi[p] += extra;
  }
  return phi;
}

double EwaldRegion::ApplyFieldTo(std::vector<PolarSegment>& foreground) const {
  if (!loaded_) {
    NotInitialized("ApplyFieldTo");
  }
  if (!real_sum_) {
    BuildSums(foreground);
  }
  // Also on the first call: with a registered foreground BuildSums ignores
  // its argument, so this is what catches a client asking about a segment
  // nobody declared.
  CheckForegroundIsSubset(foreground);

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
  //     the rest of the cell.
  //
  //     NOTHING is held out of S_bg. Every segment contributes at every
  //     image; what is removed instead is the erf-screened energy of the
  //     COINCIDENT copy of each foreground segment, because that is
  //     exactly what real space suppresses.
  //
  //     Two earlier versions were wrong in opposite directions and the
  //     answer sits between them. Holding the whole foreground out of
  //     S_bg also deletes its periodic IMAGES -- a lattice of vacancies
  //     rather than one carved-out cavity, worth 1e-3 eV of
  //     alpha-dependence on a charged 18-segment job. Holding out only
  //     the target's own segment fixed that but still deleted each
  //     segment's interaction with its own images, and was needed only
  //     because EwaldRealSpaceSum was skipping that segment at every
  //     translation. With that skip gone, real space suppresses the same
  //     set for every target, so one structure factor is correct.
  //
  //     recip + shape together are the erf interaction over all images
  //     (the shape term IS the k=0 limit the reciprocal sum omits), so
  //     subtracting a pair's erf energy removes that pair exactly. For a
  //     segment's own coincident copy that subtraction is the r -> 0
  //     branch of CalcErfStaticEnergy -- whose sign had to be fixed
  //     before this change was possible.
  {
    double e_recip = 0.0;
    double e_shape = 0.0;
    double e_erf = 0.0;

    std::vector<std::pair<const PolarSite*, Eigen::Vector3d>> fg_sites;
    fg_sites.reserve(targets.size());
    for (std::size_t n = 0; n < targets.size(); ++n) {
      fg_sites.push_back({targets[n], targets[n]->getPos()});
    }
    const std::vector<const PolarSite*> no_exclusions;

    e_recip = recip_sum_->CalcStaticEnergyBetween(fg_sites, no_exclusions,
                                                  EwaldChargeState::Neutral);
    e_shape = shape_->CalcStaticEnergyBetween(fg_sites, no_exclusions,
                                              EwaldChargeState::Neutral);

    // Every foreground copy, the target's own included. Real space
    // suppressed exactly these, so exactly these come back out.
    for (const auto& copy : foreground_copies_) {
      const PolarSegment& bg_copy =
          registry_.Get(copy.first, EwaldChargeState::Neutral);
      const Eigen::Vector3d shift = copy.second - Centroid(bg_copy);
      for (const PolarSite& source : bg_copy) {
        for (const auto& entry : fg_sites) {
          e_erf += interactor_->CalcErfStaticEnergy<PolarSite, PolarSite>(
              source, *entry.first, shift);
        }
      }
    }
    energy += e_recip + e_shape - e_erf;

    // (5b), (6b) reciprocal and shape partners of the induced-source
    //      term added at (1b), with exactly the same exclusion structure
    //      as (5) above and for exactly the same reason.
    double e_recip_pu = recip_sum_->CalcInducedSourceEnergyBetween(
        fg_sites, no_exclusions, EwaldChargeState::Neutral);
    double e_shape_pu = shape_->CalcInducedSourceEnergyBetween(
        fg_sites, no_exclusions, EwaldChargeState::Neutral);
    double e_erf_pu = 0.0;
    for (const auto& copy : foreground_copies_) {
      const PolarSegment& bg_copy =
          registry_.Get(copy.first, EwaldChargeState::Neutral);
      const Eigen::Vector3d shift = copy.second - Centroid(bg_copy);
      for (const PolarSite& source : bg_copy) {
        for (const auto& entry : fg_sites) {
          e_erf_pu += interactor_->CalcErfInducedSourceEnergy(
              source, *entry.first, shift);
        }
      }
    }
    energy += e_recip_pu + e_shape_pu - e_erf_pu;

    // Term-by-term report, for comparison against the legacy `ewald`
    // job calculator's terms_o block. The correspondence is now
    // one-to-one -- real <-> R_pp, recip <-> K_pp, shape <-> J_pp,
    // erf <-> C_pp -- since this code stopped excluding foreground
    // copies from S_bg and started subtracting their erf energy the way
    // legacy does. Measured agreement on an 18-segment job: every term
    // to legacy's six printed figures, on both a rank-0 and an
    // artificially dipolar methane. Printed in eV, the unit legacy
    // reports and the job XML carries.
    const double h2ev = tools::conv::hrt2ev;
    XTP_LOG(Log::error, log_)
        << TimeStamp()
        << " Ewald energy [eV], permanent x permanent: real = " << e_real * h2ev
        << "  recip = " << e_recip * h2ev << "  shape = " << e_shape * h2ev
        << "  erf = " << e_erf * h2ev << std::flush;
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
        << TimeStamp() << " Ewald split: alpha = " << params_.alpha
        << " 1/bohr (" << params_.alpha * 18.8972612 << " 1/nm)"
        << ", k_max = " << params_.k_max << " 1/bohr ("
        << params_.k_max * 18.8972612 << " 1/nm)"
        << ", r_min = " << params_.r_min
        << " bohr, V = " << params_.box.determinant()
        << " bohr^3, thole_a = " << params_.thole_a << std::flush;
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Foreground: " << foreground_copies_.size()
        << " segments, " << targets.size() << " sites, carved from "
        << registry_.AllIds().size() << " registered" << std::flush;

    // Suppression audit. A shortfall means a foreground copy was left in
    // the background -- silent otherwise, and exactly the failure this
    // check exists for. Every copy is suppressed for every target, the
    // target's own segment included, hence size() and not size() - 1.
    const EwaldRealSpaceSum::NeighborStats stats =
        real_sum_->GetNeighborStats();
    const Index expected =
        Index(targets.size()) * Index(foreground_copies_.size());
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Foreground copies suppressed: " << stats.foreground
        << " of " << expected
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
        << " NOTE: the background's own internal energy is deliberately "
           "absent -- it is a constant that cancels in any charge-state "
           "difference. Everything else is here; see ApplyFieldTo's "
           "declaration for the full accounting."
        << std::flush;
  }
  return energy;
}

}  // namespace xtp
}  // namespace votca
