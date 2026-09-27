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

#pragma once
#ifndef VOTCA_XTP_EWALDREGION_H
#define VOTCA_XTP_EWALDREGION_H

// Local VOTCA includes
#include <memory>

#include "votca/xtp/ewaldparameters.h"
#include "votca/xtp/ewaldrealspacesum.h"
#include "votca/xtp/ewaldreciprocalspacesum.h"
#include "votca/xtp/ewaldregistry.h"
#include "votca/xtp/ewaldshapecorrection.h"
#include "votca/xtp/region.h"

namespace votca {
namespace xtp {
class QMRegion;
class PolarRegion;
class StaticRegion;

/**
 * \brief The periodic Ewald background, as a Region.
 *
 * SKELETON. The region plumbing is complete and correct; the periodic
 * field evaluation is deliberately absent, and every entry point that
 * would need it throws rather than returning a plausible zero. Filling
 * it in is the next step.
 *
 * This region is unlike the others in three ways, all consequences of
 * the same fact -- it represents an environment that was converged
 * once, beforehand, by the ewaldbackground calculator:
 *
 *  1. Its induced dipoles are FROZEN. Other regions do not polarize it,
 *     so InteractwithQMRegion/PolarRegion/StaticRegion genuinely return
 *     0.0. Those are not stubs: "nothing happens" is the correct
 *     physics, and the same one-way relationship legacy has, where the
 *     background is computed once and the foreground alone re-polarizes.
 *
 *  2. It is therefore always Converged(), and Reset() does nothing.
 *     Including it in a job cannot prevent the inter-region SCF loop
 *     from terminating.
 *
 *  3. It does not own segments carved out of the job's topology. It
 *     represents the whole periodic cell, read from the background
 *     calculator's own checkpoint, and the segments a polar region
 *     holds explicitly are SUPPRESSED from it (see
 *     EwaldRealSpaceSum's own foreground documentation) rather than
 *     removed from its ownership. size() therefore reports the
 *     background's segment count, not a share of the job's segments.
 *
 * The influence in the other direction -- this region acting ON a polar
 * or QM region -- is implemented by THAT region's own
 * InteractwithEwaldRegion, following the existing convention in
 * Region::ApplyInfluenceOfOtherRegions, where each region computes the
 * influence of others upon itself.
 */
class EwaldRegion : public Region {
 public:
  EwaldRegion(Index id, Logger& log) : Region(id, log) {}
  ~EwaldRegion() override = default;

  std::string identify() const override { return "ewaldregion"; }

  void Initialize(const tools::Property& prop) override;

  // Always true: this region is frozen, so it can never be the reason a
  // job fails to converge. See this class's own documentation.
  bool Converged() const override { return true; }

  void Evaluate(std::vector<std::unique_ptr<Region>>& regions) override;

  // No-op: there is no per-iteration state to clear.
  void Reset() override {}

  Index size() const override;

  double charge() const override;

  double Etotal() const override;

  void WriteToCpt(CheckpointWriter& w) const override;

  void ReadFromCpt(CheckpointReader& r) override;

  void WritePDB(csg::PDBWriter& writer) const override;

 protected:
  void AppendResult(tools::Property& prop) const override;

  // All three are genuinely 0.0, not placeholders: the background is
  // frozen and no other region polarizes it. See this class's own
  // documentation, point 1.
  double InteractwithQMRegion(const QMRegion&) override { return 0.0; }
  double InteractwithPolarRegion(const PolarRegion&) override { return 0.0; }
  double InteractwithStaticRegion(const StaticRegion&) override { return 0.0; }
  double InteractwithEwaldRegion(const EwaldRegion&) override { return 0.0; }

 public:
  // The converged periodic background, as written by the ewaldbackground
  // calculator. Read-only to consumers: nothing in a job re-polarizes it.
  const EwaldRegistry& Registry() const { return registry_; }

  // The convergence parameters the background was converged WITH, read
  // from the same checkpoint rather than re-specified by the job. See
  // EwaldParameters for why that distinction matters.
  const EwaldParameters& Parameters() const { return params_; }

  // Accumulates this periodic background's field into every site of
  // `foreground`, which is the polar region's own segment list.
  //
  // The four contributions, matching legacy's own ER - EC + EK + E0
  // decomposition:
  //
  //   + real space (erfc), with the foreground's own copies SUPPRESSED
  //     so they are not counted twice -- they are about to be treated
  //     explicitly by the polar region instead
  //   + reciprocal space, over the FULL periodic density. This still
  //     contains the foreground segments in their neutral state: the
  //     k-sum is over the whole lattice and cannot have a hole cut in it
  //   + shape/surface term
  //   - the erf-screened field of exactly those neutral foreground
  //     copies, which removes what the reciprocal sum just put back
  //
  // The last term uses the NEUTRAL multipoles together with the
  // background's own converged induced dipoles, not the job's charge
  // state -- because that is what the reciprocal sum actually placed
  // there. Using the job's state instead would remove something that was
  // never added, and would make a charged job differ from the background
  // for the wrong reason.
  //
  // Returns 0.0: the energy is not implemented yet (see the warning this
  // emits on first use).
  double ApplyFieldTo(std::vector<PolarSegment>& foreground) const;

  // Declares the COMPLETE foreground: the union of every region that owns
  // segments of the job's topology. JobTopology calls this once, after it
  // has built the regions and before any of them is evaluated.
  //
  // Why it cannot be left to ApplyFieldTo's argument. Each client region
  // asks for its own segments, but the background copies that must be
  // suppressed are those of the WHOLE carved-out cluster. In a QM + polar
  // job the QM region holds the central segment, so the polar region
  // hands over everything except it -- and its neutral background copy
  // would be left sitting underneath the QM density, a ghost molecule in
  // every sum, with no symptom to notice it by.
  //
  // Positions are the job-local ones, after JobTopology::ShiftPBC: the
  // background registry is in cell coordinates and the job is recentred,
  // so the two frames differ by a lattice vector. BuildSums snaps that
  // out. Only the lattice image is taken from these centroids, so any
  // reasonable centre of a segment will do.
  void RegisterForeground(
      const std::vector<std::pair<Index, Eigen::Vector3d>>& foreground);

  // The background's electrostatic POTENTIAL at arbitrary points, for a
  // QM region that takes its environment as a potential on a grid rather
  // than as a field on sites.
  //
  // The same four terms ApplyFieldTo assembles -- real space with the
  // foreground copies suppressed, reciprocal space, shape, minus the erf
  // half of the suppressed copies -- and both channels of each, so what
  // comes back is the potential of the permanent background TOGETHER
  // WITH its converged induced dipoles.
  //
  // GAUGE. phi carries an arbitrary additive constant: the k = 0 term is
  // omitted, i.e. the uniform neutralising background. Every term here
  // is evaluated through the same energy routines the classical channels
  // use, with a unit test charge, so the constant is the one every
  // validated number in this code was computed with. A charged region
  // embedded in phi shifts by q * phi_0, so this is not a free choice --
  // see unit_probe_potential_reproduces_the_static_energy, which pins it.
  //
  // Requires a declared foreground (RegisterForeground): a point is not
  // a segment, so there is nothing to fall back on.
  Eigen::VectorXd PotentialAt(const std::vector<Eigen::Vector3d>& points) const;

 private:
  // Built once, on first use, from the calling region's geometry. The
  // foreground's positions are fixed for the whole job even though its
  // dipoles change every SCF iteration, so the sums -- and the neighbour
  // cache inside the real-space one -- stay valid throughout. Mutable
  // for the same reason EwaldRealSpaceSum's own cache is: this is lazily
  // built state behind a const interface, not mutable physics.
  //
  // What is cached is specific to the foreground it was built for:
  // foreground_copies_ decides which background copies the real-space sum
  // suppresses, and that list is baked into real_sum_'s constructor along
  // with a neighbour cache keyed to those positions.
  //
  // Built from registered_foreground_ when RegisterForeground has been
  // called, which is the case inside a JobTopology. The argument is a
  // FALLBACK for direct use of this class without one -- the unit tests,
  // and any single-client setup -- and is only correct when there is
  // exactly one client region. Either way every later call is checked
  // against what was built.
  void BuildSums(const std::vector<PolarSegment>& fallback) const;
  // Every segment a client asks about must be part of the foreground the
  // sums were built for. A subset, not an equality: each client region
  // asks only about its own share of the union.
  void CheckForegroundIsSubset(
      const std::vector<PolarSegment>& foreground) const;

  mutable std::unique_ptr<EwaldRealSpaceSum> real_sum_;
  mutable std::unique_ptr<EwaldReciprocalSpaceSum> recip_sum_;
  mutable std::unique_ptr<EwaldShapeCorrection> shape_;
  mutable std::unique_ptr<EwaldRealSpaceInteractor> interactor_;
  // (segment id, position) of each suppressed foreground copy, so the
  // erf correction removes exactly the copies the real-space sum
  // dropped.
  mutable std::vector<std::pair<Index, Eigen::Vector3d>> foreground_copies_;
  // (segment id, centroid) of the foreground the sums were built for.
  // Keyed on the FOREGROUND's own centroid, not the background copy's:
  // foreground_copies_ holds the latter, which is derived from the id
  // alone, so two different foregrounds sharing a segment id would look
  // identical there.
  mutable std::vector<std::pair<Index, Eigen::Vector3d>> built_foreground_;
  // What RegisterForeground was told. Empty means nobody declared a
  // foreground and BuildSums falls back to its argument.
  std::vector<std::pair<Index, Eigen::Vector3d>> registered_foreground_;
  mutable bool warned_no_energy_ = false;

  std::string checkpoint_file_;
  EwaldRegistry registry_;
  EwaldParameters params_;
  bool loaded_ = false;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDREGION_H
