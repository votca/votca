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
    foreground_copies_.push_back({seg.getId(), Centroid(seg)});
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
        << " WARNING: the Ewald background contributes its FIELD to the "
           "polar region, but its contribution to the ENERGY is not "
           "implemented. Induced dipoles are meaningful; reported "
           "energies are not."
        << std::flush;
  }
  return 0.0;
}

}  // namespace xtp
}  // namespace votca
