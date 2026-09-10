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
#ifndef VOTCA_XTP_EWALDPERIODICDIPOLEOPERATOR_H
#define VOTCA_XTP_EWALDPERIODICDIPOLEOPERATOR_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "eeinteractor.h"
#include "eigen.h"
#include "ewaldreciprocalspacesum.h"
#include "ewaldrealspacesum.h"
#include "ewaldregistry.h"
#include "ewaldshapecorrection.h"

/**
 * \brief Matrix-free periodic induced-dipole self-consistent-field
 *        operator, for use with Eigen::ConjugateGradient.
 *
 * This mirrors DipoleDipoleInteraction's own existing interface (see
 * dipoledipoleinteraction.h) exactly -- same required Eigen::EigenBase
 * typedefs/constants, same InnerIterator/operator()/multiply() shape --
 * so it plugs into Eigen::ConjugateGradient the same way
 * PolarRegion::CalcInducedDipolesViaPCG already uses
 * DipoleDipoleInteraction, with the periodic real+reciprocal field
 * (EwaldRealSpaceSum + EwaldReciprocalSpaceSum) standing in for
 * DipoleDipoleInteraction's own aperiodic eeInteractor::
 * FillTholeInteraction wherever the "interaction between two sites"
 * building block is needed.
 *
 * Segments in `ids` may have any number of sites each (this was not
 * always true of this class -- an earlier version required exactly one
 * site per segment; generalizing it was deliberately sequenced after
 * that simpler version was already validated, so that any new bug
 * introduced by the generalization itself would be easy to isolate from
 * bugs in the underlying physics, which by this point had already been
 * separately confirmed). The vector layout any v/x this operator is
 * applied to must follow is: for each id in ids (in order), 3 entries
 * per site of that segment (xyz), in the order EwaldRegistry::Get(id,
 * Neutral)'s own site iteration gives them -- i.e. a variable-width,
 * offset-table ("CSR-style") layout, not a fixed 3-per-segment one.
 *
 * multiply(v) computes A*v, where A is the induced-dipole self-consistent
 * operator restricted to ids: A*v = D*v + (periodic field that induced
 * dipoles v produce at every ids site, from every OTHER ids site's
 * periodic images, at every translation including zero -- see below for
 * why zero-translation same-segment pairs are no longer excluded here)
 * + (intramolecular erfc-screened, Thole-damped coupling between
 * distinct sites of the SAME segment, at zero translation). D is the
 * same block-diagonal inverse-polarizability term (site.getPInv())
 * already used by DipoleDipoleInteraction.
 *
 * The intramolecular term above was NOT present in an earlier version of
 * this class, based on a mistaken belief that DipoleDipoleInteraction
 * excludes intra-segment coupling the same way EwaldRealSpaceSum does for
 * the periodic sum (see that class's own documentation) -- checking
 * DipoleDipoleInteraction::multiply() directly shows this is wrong: its
 * inner loop couples every pair of sites in its own flat sites_ list via
 * FillTholeInteraction, with no segment-awareness at all beyond
 * operator()'s own same-site (not same-segment) diagonal check. Both
 * legacy PolarBackground (see its own explicitly-labeled "Real-space,
 * intramolecular contribution" loop, inside its own induction-iteration
 * section specifically -- NOT its permanent-field section, which has no
 * intramolecular loop at all) and modern PolarRegion/
 * DipoleDipoleInteraction therefore include intramolecular induced-
 * induced coupling as standard practice; omitting it here was an actual
 * correctness gap in this class, not a deliberate scoping choice, caught
 * only once a real full-polarizability comparison against legacy
 * surfaced a large, systematic, atom-specific discrepancy that a static-
 * field-only (Tier 1) comparison could never have exposed, since this
 * term contributes nothing until induction is actually active.
 *
 * This class's own damping choice for that term went through TWO
 * mistaken revisions before landing on Thole-damped, worth recording in
 * full since both were confident and both were wrong in different
 * directions:
 *
 *   1. The first version used eeInteractor::FillTholeInteraction
 *      (Thole-damped) for the intramolecular tensor, copying
 *      DipoleDipoleInteraction's own convention directly -- reasonable
 *      on its face, but DipoleDipoleInteraction is aperiodic (no erfc/
 *      alpha concept at all), so "Thole-only" there isn't a considered
 *      choice about periodic intramolecular coupling specifically, it's
 *      the only mechanism that exists in an aperiodic context at all --
 *      not evidence about what a periodic class should do.
 *
 *   2. The second version, correcting the first, used
 *      EwaldRealSpaceInteractor::ComputeB directly, deliberately
 *      UNDAMPED by Thole -- reasoning by analogy from the PERMANENT-
 *      field case (EwaldRealSpaceInteractor::
 *      ApplyIntramolecularStaticCorrection, which genuinely has no
 *      Thole factor, because legacy's own static intramolecular
 *      treatment has no real-space loop at all, only a reciprocal-space
 *      erf compensation). That analogy does not transfer: legacy's
 *      INDUCED case is structurally different. Its own FU12_ERFC_At_By
 *      -- the SAME method used for both intermolecular (PolarBackground,
 *      confirmed at its own real-space induction loop, roughly line 954)
 *      and intramolecular (confirmed at its own explicitly-labeled
 *      "Real-space, intramolecular contribution" loop, roughly line
 *      453-454) induced-induced real-space pairs -- applies Thole
 *      damping via its own l3/l5 mechanism unconditionally, with no
 *      special-casing for same-segment pairs at either call site. Static
 *      and induced genuinely differ in legacy's own code (no real-space
 *      intramolecular loop at all for statics; a real, Thole-damped one
 *      for induction) -- generalizing from one to the other, in two
 *      different directions across these two mistakes, was the actual
 *      error both times, not a single reversible sign or magnitude slip.
 *
 * A separate, related check was also worth doing directly rather than by
 * analogy: is legacy's OWN Thole damping model here (its L3()/L5(),
 * called via UpdateAllBls from FU12_ERFC_At_By) the same exponential
 * (Thole/van Duijnen-Swart) form ComputeThole already uses, or a
 * different (e.g. linear) one? An earlier investigation concluded
 * "linear" -- but that conclusion was traced from XInteractor, which
 * PolarBackground constructs and uses ONLY inside its own explicitly-
 * labeled "CUTOFF TREATMENT" branches (coulombmethod.method="cutoff",
 * not the "ewald" default this class is actually compared against).
 * Tracing EwdInteractor's own L3()/L5() directly instead -- confirmed as
 * the class PolarBackground actually uses by default, via its own
 * _ewdactor member and the "ewald" default in ewdbgpol.xml's own schema
 * -- shows L3()=1-exp(-ta1*tu3), L5()=1-(1+ta1*tu3)*exp(-ta1*tu3): the
 * SAME exponential form ComputeThole already implements, with ta1
 * traced back to the same polarmethod.aDamp/thole_a parameter. So this
 * class's own Thole model was never the mismatch; AddIntraSegmentCoupling
 * calling ComputeThole now (with the real thole_a -- see this class's
 * own constructor documentation) is a genuine fix, not a model-form
 * change.
 *
 * A FOURTH, related mismatch was found in EwaldReciprocalSpaceSum: an
 * earlier version of that class excluded the target's own segment from
 * its own reciprocal-space sum entirely, modeled on EwaldRealSpaceSum's
 * own real-space exclusion. Tracing legacy's own reciprocal-space code
 * (PolarBackground::KThread::SP_SFactorCalc/FP_KFieldCalc for the
 * permanent case, SU_SFactorCalc/FU_KFieldCalc for the induced case)
 * shows it never excludes anything -- every registered site, including
 * the target's own segment, contributes to the one global structure
 * factor applied back to every site. This matters specifically because
 * of the Ewald-split identity erfc(a*r)/r + [reciprocal-space
 * contribution] = 1/r EXACTLY: if EwaldRealSpaceSum's own real-space sum
 * legitimately excludes a same-segment, zero-translation pair (which it
 * still correctly does, unchanged), the reciprocal-space contribution to
 * that SAME pair must still be included for the two Ewald halves to
 * combine into the correct total wherever they're meant to -- e.g. here,
 * where the erfc-screened intramolecular term above is only the real-
 * space HALF of the intended full-strength (Thole-damped) intramolecular
 * interaction; the other half arrives automatically once
 * EwaldReciprocalSpaceSum stops excluding anything at all (see that
 * class's own documentation for the fuller account). This is why the
 * periodic-field term's own documentation above no longer says
 * "excluding that site's own segment's zero-translation copy" the way
 * an earlier revision of this comment did -- only EwaldRealSpaceSum's
 * own real-space exclusion remains; EwaldReciprocalSpaceSum's own
 * contribution is unconditional now.
 *
 * The periodic-field term is computed by temporarily writing v into
 * registry_'s induced-dipole slots and evaluating EwaldRealSpaceSum::
 * AddFieldAt / EwaldReciprocalSpaceSum::AddFieldAtMany at every ids site;
 * confirmed symmetric for this combination when restricted to ids-only
 * coupling (see test_ewaldperiodicdipolesymmetry.cc -- written against
 * the earlier, segment-excluding version of EwaldReciprocalSpaceSum;
 * worth re-confirming this still holds now that exclusion is gone
 * entirely, though the underlying argument -- a symmetric physical
 * kernel stays symmetric whether or not a term is added to both sides of
 * a pair uniformly -- does not depend on that exclusion having existed),
 * a prerequisite for ConjugateGradient. The intramolecular term is
 * added/transposed into both sites' own result the same block/
 * block.transpose() pairing DipoleDipoleInteraction::multiply() itself
 * uses, which keeps the combined operator symmetric overall. This
 * intramolecular contribution is linear in v and vanishes at v=0, so it
 * does not affect baseline_ (see below) at all.
 *
 * The sign/scale convention relative to
 * DipoleDipoleInteraction's own established one has since been checked
 * against a physically unambiguous case -- a single polarizable site in
 * one fixed external charge's field must develop an induced dipole
 * aligned with that field -- and matches (see
 * test_ewaldperiodicdipoleoperator.cc); this does not by itself prove
 * correctness for larger, genuinely self-consistent (multi-site,
 * multi-iteration) cases, which is worth validating separately before
 * trusting this for anything real -- and is worth re-validating again
 * specifically for the genuinely-multi-site-per-segment case this
 * generalization introduces, since neither existing test exercises that
 * case at all. The intramolecular-coupling addition above is itself
 * NOT YET covered by either existing test (both use single-site
 * segments where the term is identically zero) -- validating it
 * directly, e.g. against a small system where PolarRegion's own
 * DipoleDipoleInteraction can be run for comparison, is a genuine open
 * task, not something this comment should be read as already having
 * confirmed.
 *
 * A genuine bug was caught and fixed here during that same validation:
 * AddFieldAt/AddFieldAtMany sum the field from every OTHER registered
 * segment, not just other members of ids -- so a naive multiply(v) would
 * silently include a constant contribution from any registered segment
 * outside ids (e.g. a fixed external background), independent of v. That
 * breaks the linearity ConjugateGradient requires (multiply(0) must be
 * 0). The fix: baseline_ (computed once at construction, as the same raw
 * computation evaluated at v=0) captures exactly that external leak, and
 * multiply(v) subtracts it, leaving only the genuine ids-internal
 * coupling. Segments outside ids still correctly influence the physics
 * overall -- through the permanent field baked into b before this
 * operator is ever invoked -- just not redundantly inside A itself.
 *
 * IMPORTANT, and easy to get backwards: the right-hand side b to pass to
 * ConjugateGradient::solveWithGuess is b = +V_permanent, not b = -V, even
 * though PolarRegion::CalcInducedDipolesViaPCG itself uses b = -V for the
 * (legacy-derived) DipoleDipoleInteraction operator. The two are not
 * interchangeable conventions to copy blindly: eeInteractor's own V()
 * stores the *negative* of the physical field (confirmed via its own
 * energy expression, e = q*phi + mu.V, which only matches the standard
 * U = q*phi - mu.E if V = -E), whereas EwaldRealSpaceInteractor/
 * EwaldRealSpaceSum/EwaldReciprocalSpaceSum's own V() stores the genuine,
 * un-negated physical field (confirmed via test_ewaldrealspaceinteractor
 * .cc's own explicit CalcStaticEnergy formula, q*phi - mu.E). This sign
 * mismatch was caught the hard way -- as a factor-of-exactly-(-1)
 * discrepancy in test_ewaldperiodicdipoleoperator.cc, after the
 * double-counting bug above was already fixed and could no longer
 * explain it -- rather than reasoned out from first principles, which is
 * exactly why it is worth stating this explicitly here rather than
 * trusting a future caller to rediscover it.
 *
 * operator()(i,j) is used by Eigen's preconditioner setup (via
 * InnerIterator) to extract the literal diagonal only; DiagonalPrecon-
 * ditioner never reads an off-diagonal value even though it is visited
 * during setup. DipoleDipoleInteraction's own operator()(i,j) computes
 * the real interaction tensor for every cross-site pair regardless (its
 * own comment notes this is "not a fast method"), which is tolerable
 * there because that tensor is O(1) per pair (aperiodic). Here, the
 * periodic per-pair tensor would cost a full real+reciprocal evaluation,
 * so cross-site entries return 0.0 unconditionally rather than computing
 * (and discarding) the real value -- same-site entries still return the
 * genuine getPInv() block, exactly as DipoleDipoleInteraction does, since
 * that is both cheap and the only part of operator()'s output the
 * preconditioner actually keeps. "Same-site" here means the literal same
 * PolarSite, not merely the same segment -- two different sites of the
 * same multi-site segment are cross-site entries here too (0.0), even
 * though multiply() now does couple them via the intramolecular term
 * described above: that coupling is never read by
 * DiagonalPreconditioner either (same reasoning as the periodic term),
 * so operator()'s own 0.0 for it is equally safe to leave uncomputed,
 * not an oversight relative to multiply()'s own, more complete behavior.
 *
 * getPInv() alone was NOT actually the true diagonal for a real interval
 * spanning this class's own reciprocal-space-no-exclusion change and the
 * self_field_matrix_ fix below: with the former in place but not yet the
 * latter, a site's own true diagonal was getPInv() + (that site's own
 * spurious reciprocal-space self-field contribution, see
 * self_field_matrix_'s own documentation) -- a real, if modest-looking
 * (~4% relative to getPInv() for the specific system this was first
 * noticed on), mismatch between what operator() reported and what
 * multiply() actually used, silently degrading DiagonalPreconditioner's
 * own effectiveness (observed as PCG iteration counts roughly
 * quadrupling for one real system, alongside the separate, larger
 * reciprocal-space-related slowdown from that same change). Adding
 * self_field_matrix_ back into multiply()'s own diagonal (see below)
 * cancels that same self-field term there too, restoring getPInv() as
 * the genuine true diagonal and operator()'s own reported value as
 * correct again -- not a coincidence: both fixes address the same
 * underlying leak, just in two different places it needed correcting.
 *
 * Registered segments not included in `ids` (the set this operator
 * solves for) are treated as a fixed, non-reacting external background:
 * their own sites still contribute to the periodic field seen by the
 * segments in `ids` (via the same intermolecular-only exclusion already
 * implemented in EwaldRealSpaceSum/EwaldReciprocalSpaceSum), but their
 * own induced dipoles are never themselves solved for or perturbed by
 * multiply(). This is what lets this operator solve for the induced
 * dipoles of, e.g., a growing/adaptive subset of sites without needing
 * to re-solve the whole periodic system's own dipoles at once.
 *
 * registry, real_sum, and recip_sum are all stored by reference and must
 * outlive this object; real_sum and recip_sum must have been constructed
 * against the same registry passed here.
 *
 * Units: bohr, Hartree atomic units throughout.
 */

namespace votca {
namespace xtp {
class EwaldPeriodicDipoleOperator;
}
}  // namespace votca

namespace Eigen {
namespace internal {
// EwaldPeriodicDipoleOperator's traits must be specialized before the
// class itself inherits from Eigen::EigenBase<EwaldPeriodicDipoleOperator>
// below (Eigen::EigenBase's own definition instantiates
// internal::traits<Derived> immediately) -- matching the exact ordering
// DipoleDipoleInteraction's own header (dipoledipoleinteraction.h) uses
// for the same reason.
template <>
struct traits<votca::xtp::EwaldPeriodicDipoleOperator>
    : public Eigen::internal::traits<Eigen::MatrixXd> {};
}  // namespace internal
}  // namespace Eigen

namespace votca {
namespace xtp {

class EwaldPeriodicDipoleOperator
    : public Eigen::EigenBase<EwaldPeriodicDipoleOperator> {
 public:
  // Required typedefs, constants, and method:
  using Scalar = double;
  using RealScalar = double;
  using StorageIndex = votca::Index;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  // ids: the segments this operator solves for, in the order that
  //   defines the layout of any vector this operator is applied to (see
  //   class documentation for the exact, variable-width layout). Every
  //   id must be registered in registry at EwaldChargeState::Neutral.
  // alpha_ewald: Ewald splitting parameter, shared with real_sum/
  //   recip_sum's own construction, needed here to build a matching
  //   erfc-screened intramolecular interactor.
  // thole_a: the SAME Thole damping parameter real_sum's own
  //   construction uses (polarmethod.aDamp in ewdbgpol.xml terms). This
  //   used to be an unused placeholder here (0.39, hardcoded, since
  //   intra_interactor_ was originally only ever used for its own public
  //   ComputeB): AddIntraSegmentCoupling did not apply Thole damping to
  //   the intramolecular term at all, on the -- since corrected --
  //   assumption that legacy's own treatment was undamped there too (see
  //   AddIntraSegmentCoupling's own documentation for the fuller
  //   history). It now genuinely calls ComputeThole, so this parameter
  //   is genuinely used and must be the real value, not a placeholder.
  // shape: the SAME EwaldShapeCorrection instance EwaldBackground's own
  //   permanent-field generation uses. This class previously had NO
  //   shape-correction mechanism at all for the induced case -- a real,
  //   previously-unnoticed gap, found only once EwdInteractor's own
  //   FU12_ShapeField_At_By was traced directly (called unconditionally
  //   inside PolarBackground's own induction-iteration loop, mirroring
  //   FP12_ShapeField_At_By's own permanent-field treatment exactly, but
  //   using each site's own CURRENT induced dipole rather than its
  //   static one). EwaldShapeCorrection::TotalDipoleMoment already sums
  //   charge + static dipole + induced dipole together (see its own
  //   documentation), so no new shape-correction class or formula is
  //   needed here -- only a new call site, added inside RawMultiply
  //   itself (see there for why this is genuinely linear in v and does
  //   not disturb baseline_'s own established v=0 subtraction pattern).
  // apply_shape_correction: debug-only, defaults to true (matches
  //   legacy). Exists ONLY to isolate whether the shape term above is
  //   the cause of a real, large PCG iteration-count regression found
  //   once it and the Thole-damping fix were both added in the same
  //   session -- an open question when this parameter was added, not a
  //   genuine design choice; see EwaldBackground's own
  //   apply_shape_correction_to_induced_ member for the fuller account
  //   and why this has no legacy counterpart or XML schema entry.
  //   Setting this false is a deliberate DEPARTURE from matching legacy
  //   (which always applies its own FU12_ShapeField_At_By
  //   unconditionally) -- never meant to ship as a real option.
  // apply_self_field_correction: debug-only, same purpose and pattern as
  //   apply_shape_correction immediately above -- added alongside it so
  //   the same regression can be isolated against self_field_matrix_
  //   (see that member's own documentation) independently, and against
  //   both together, rather than leaving this one variable untested.
  //   Setting this false is a deliberate DEPARTURE from matching legacy
  //   (which always applies its own atomic self-interaction correction,
  //   FU12_ERF_At_By called with the same site as both arguments) --
  //   never meant to ship as a real option, same as
  //   apply_shape_correction.
  // apply_thole_damping_intramolecular: debug-only, same purpose and
  //   pattern as the two above -- added once both of those were ruled
  //   out (via their own matching debug flags, confirmed to leave a real
  //   PCG divergence entirely unaffected when disabled) as the cause of
  //   a real, large divergence found once this session's own changes
  //   were all in place together, making this the last untested
  //   mechanism from that same session. Gates
  //   AddIntraSegmentCoupling's own use of ComputeThole -- see that
  //   method's own documentation for the fuller history of two earlier,
  //   different mistakes in that same term. Setting this false is a
  //   deliberate DEPARTURE from matching legacy (which applies Thole
  //   damping to intramolecular induced-induced pairs unconditionally,
  //   same as intermolecular ones) -- never meant to ship as a real
  //   option, same as the other two.
  // apply_realspace_intermolecular_coupling / apply_reciprocal_coupling:
  //   debug-only, same purpose and pattern as the three above -- added
  //   after those three turned out to be small contributions on both
  //   this codebase's own side and legacy's own, leaving a real,
  //   substantial coupling-term mismatch found this session unlocalized
  //   any further than "somewhere in {thole, real-space intermolecular,
  //   reciprocal-space} combined". Gate real_sum_.AddFieldAt and
  //   recip_sum_.AddFieldAtMany respectively inside RawMultiply. Setting
  //   either false is a deliberate DEPARTURE from matching legacy
  //   (which always applies its own FX_RealSpace/FX_ReciprocalSpace
  //   unconditionally) -- never meant to ship as a real option, same as
  //   the other three.
  EwaldPeriodicDipoleOperator(EwaldRegistry& registry,
                             const EwaldRealSpaceSum& real_sum,
                             const EwaldReciprocalSpaceSum& recip_sum,
                             const EwaldShapeCorrection& shape,
                             std::vector<Index> ids, double alpha_ewald,
                             double thole_a,
                             bool apply_shape_correction = true,
                             bool apply_self_field_correction = true,
                             bool apply_thole_damping_intramolecular = true,
                             bool apply_realspace_intermolecular_coupling =
                                 true,
                             bool apply_reciprocal_coupling = true);

  // Debug/experimental. Dumps, for every intramolecular pair this
  // class's own AddIntraSegmentCoupling would visit, the raw damping
  // scalars (r, B0, B1, B2, l3, l5) -- not the assembled field tensor
  // (already checked, this session, and found consistent with legacy's
  // own via ApplyInducedField) -- so a per-pair comparison against
  // legacy's own matching EwdInteractor::GetDebugPairState dump can
  // isolate the raw damping computation itself, independent of tensor
  // assembly or vector algebra around it. Depends only on geometry and
  // each site's own polarizability, not on any induced-dipole state, so
  // this can be called at any time (before a solve even starts) and
  // will give the same result throughout.
  void DumpIntraPairThole(const std::string& filename) const;

  // Debug/experimental. Computes op*v's own field contribution to each
  // ids_ target site, staged cumulatively in LEGACY's own order (FUa =
  // intramolecular only, FUb = +real-space intermolecular, FUc =
  // +reciprocal, FUd = +shape, FUe = +self-field), NOT this class's own
  // internal RawMultiply order (real-inter, recip, shape, self-field,
  // THEN intramolecular last) -- added specifically to let a per-stage
  // comparison against legacy's own matching stagedFU dump isolate
  // whether a discrepancy comes from field computation (would show up
  // at the specific stage it enters) or from the induction response
  // built on top of it (would show up equally at every stage, since
  // the same v is used for all five here). All five stages act on the
  // SAME input v (deliberately, to match legacy's own FUa-FUe, which
  // are all computed from the same mu_1) -- this is not RawMultiply's
  // own op*v (which mixes in P^-1*v and is gated by the constructor's
  // own toggles); every stage here is always computed and dumped,
  // regardless of this object's own apply_*_ flags, since the entire
  // point is to see every stage side by side.
  void DumpStagedCoupling(const Eigen::VectorXd& v,
                          const std::string& filename) const;

  // Debug/experimental. Wraps EwaldRealSpaceSum::DumpPerPairFieldAppend
  // for a single target segment's own first site, so the calculator
  // (which has no direct access to real_sum_, a private member here)
  // can trigger it without needing its own accessor for the whole
  // real_sum_ object. Caller is responsible for having already set
  // every site's own induced dipole to whatever state (e.g. x1) the
  // comparison is meant to reflect -- this wrapper does not set any
  // dipole state itself, matching DumpPerPairFieldAppend's own
  // contract.
  void DumpPerPairIntermolecularField(Index target_segment_id,
                                      const std::string& filename) const;

  // Debug/experimental. Wraps EwaldRealSpaceSum::
  // DumpCachedNeighborListAppend the same way
  // DumpPerPairIntermolecularField wraps DumpPerPairFieldAppend -- see
  // that method's own comment for why this wrapper exists. Call this
  // AFTER whatever call is suspected of populating (or not populating,
  // or populating differently than expected) real_sum_'s own neighbor
  // cache for this target -- e.g. right after DumpStagedCoupling.
  void DumpCachedNeighborList(Index target_segment_id,
                              const std::string& filename) const;

  class InnerIterator {
   public:
    InnerIterator(const EwaldPeriodicDipoleOperator& xpr, const Index& id)
        : xpr_(xpr), id_(id) {};

    InnerIterator& operator++() {
      row_++;
      return *this;
    }
    operator bool() const { return row_ < xpr_.size_; }
    double value() const { return xpr_(row_, id_); }
    Index row() const { return row_; }
    Index col() const { return id_; }
    Index index() const { return row(); }

   private:
    const EwaldPeriodicDipoleOperator& xpr_;
    const Index id_;
    Index row_ = 0;
  };

  Index rows() const { return size_; }
  Index cols() const { return size_; }
  Index outerSize() const { return size_; }

  template <typename Vtype>
  Eigen::Product<EwaldPeriodicDipoleOperator, Vtype, Eigen::AliasFreeProduct>
  operator*(const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<EwaldPeriodicDipoleOperator, Vtype,
                          Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // See class documentation: same-site entries return the real getPInv()
  // block; cross-site entries return 0.0 unconditionally (never actually
  // read by DiagonalPreconditioner, and expensive to compute honestly
  // for a periodic pair).
  double operator()(Index i, Index j) const;

  // The actual matrix-vector product; see class documentation for what
  // this computes.
  Eigen::VectorXd multiply(const Eigen::VectorXd& v) const;

 private:
  // Computes getPInv()*v + (periodic field from every OTHER registered
  // segment, including those not in ids_) at every ids_ site. Not itself
  // linear in v (it carries the constant leak from non-ids_ sources) --
  // multiply() below fixes that by subtracting baseline_ = RawMultiply(0).
  Eigen::VectorXd RawMultiply(const Eigen::VectorXd& v) const;

  // Adds the intramolecular erfc-screened, Thole-damped contribution
  // (see class documentation for the fuller history) into result in
  // place, for every segment in ids_ with more than one site. This term
  // is linear in v and contributes 0 at v=0, so it is folded directly
  // into RawMultiply's own result rather than needing any separate
  // baseline_ handling.
  void AddIntraSegmentCoupling(const Eigen::VectorXd& v,
                              Eigen::VectorXd& result) const;

  // Global vector index i's (segment index within ids_, site index
  // within that segment) -- the reverse of the offset table below,
  // found by binary search (std::upper_bound) on offsets_. Used only by
  // operator()(i,j), which is not performance-critical (see class
  // documentation).
  std::pair<Index, Index> LocateSite(Index i) const;

  EwaldRegistry& registry_;
  const EwaldRealSpaceSum& real_sum_;
  const EwaldReciprocalSpaceSum& recip_sum_;
  const EwaldShapeCorrection& shape_;
  // Debug-only, see this class's own constructor documentation.
  bool apply_shape_correction_;
  std::vector<Index> ids_;
  // Used for its own public ComputeB (erfc-screened B-functions) AND
  // ComputeThole (Thole damping factors), constructed with the REAL
  // thole_a passed to this class's own constructor -- unlike an earlier
  // version of this class (see AddIntraSegmentCoupling's own
  // documentation), ComputeThole genuinely is called from this class
  // now, so thole_a here must be the real value, not a placeholder.
  EwaldRealSpaceInteractor intra_interactor_;
  // offsets_[n] is the global vector index (in units of a single scalar,
  // not a 3-vector) where segment ids_[n]'s own block starts;
  // offsets_[n+1]-offsets_[n] = 3 * (that segment's own site count).
  // offsets_.size() == ids_.size()+1, with offsets_.back() == size_.
  std::vector<Index> offsets_;
  Index size_;
  Eigen::VectorXd baseline_;
  // The position-independent 3x3 matrix M such that a site's own trial
  // dipole v_i produces a spurious reciprocal-space self-field -M*v_i at
  // its own position (see EwaldReciprocalSpaceSum::SelfFieldMatrix's own
  // documentation for the derivation). Legacy explicitly removes this
  // same leak for the induced-dipole case (EwdInteractor::
  // FU12_ERF_At_By, called with the SAME site passed as both arguments
  // -- confirmed by tracing PolarBackground's own SOR iteration
  // directly -- with its own explicit "Note the (-): This is a
  // compensation term" comment); this class did not, until this was
  // traced down as the explanation for a large, carbon-specific
  // discrepancy that survived fixing the (separate) intramolecular
  // static-static leak. Computed once here (position-independent, same
  // for every site, depends only on the k-vector set/alpha/volume) and
  // added back into RawMultiply's own result for every site, cancelling
  // that site's own -M*v_i leak already present via site.V() (recip_sum_
  // never excludes anything, self included -- see its own class
  // documentation).
  Eigen::Matrix3d self_field_matrix_;
  // Debug-only, no legacy counterpart, same purpose and pattern as
  // apply_shape_correction_ above -- see this class's own constructor
  // documentation. Gates whether self_field_matrix_ above is added back
  // in RawMultiply; false is a deliberate DEPARTURE from matching legacy
  // (which always applies its own atomic self-interaction correction),
  // never a default.
  bool apply_self_field_correction_;
  // Debug-only, same purpose and pattern as apply_shape_correction_ and
  // apply_self_field_correction_ above, but for AddIntraSegmentCoupling's
  // own Thole-damping term instead -- see that method's own
  // documentation for why this was added, and its own constructor
  // documentation entry immediately below for the fuller account.
  bool apply_thole_damping_intramolecular_;
  // Debug-only, same purpose and pattern as the three above -- added
  // this session specifically because those three (thole, shape, self-
  // field) turned out to be small contributions on both this codebase's
  // own side and legacy's own (confirmed directly: legacy's own FUc,
  // FUd, FUe stages -- +reciprocal, +shape, +self-interaction
  // respectively -- differ only slightly from each other, and toggling
  // shape/self-field off here changed residual2 by only ~0.3% in a real
  // run), meaning a real, substantial mismatch found in a coupling-term
  // comparison against legacy this session could not be localized any
  // further than "somewhere in {thole, real-space intermolecular,
  // reciprocal-space} combined" with only those three toggles. Gates
  // real_sum_.AddFieldAt inside RawMultiply -- see that method's own
  // documentation. Setting this false is a deliberate DEPARTURE from
  // matching legacy (which always applies its own FX_RealSpace
  // unconditionally) -- never meant to ship as a real option, same as
  // the three above.
  bool apply_realspace_intermolecular_coupling_;
  // Debug-only, same purpose and pattern as apply_realspace_
  // intermolecular_coupling_ immediately above -- added alongside it so
  // the two can be isolated independently (real-space intermolecular
  // vs. reciprocal-space) as well as together. Gates
  // recip_sum_.AddFieldAtMany inside RawMultiply. Setting this false is
  // a deliberate DEPARTURE from matching legacy (which always applies
  // its own FX_ReciprocalSpace unconditionally, when !_do_use_cutoff)
  // -- never meant to ship as a real option, same as the others.
  bool apply_reciprocal_coupling_;
};

}  // namespace xtp
}  // namespace votca

namespace Eigen {
namespace internal {
template <typename Vtype>
struct generic_product_impl<votca::xtp::EwaldPeriodicDipoleOperator, Vtype,
                            DenseShape, DenseShape, GemvProduct>
    : generic_product_impl_base<
          votca::xtp::EwaldPeriodicDipoleOperator, Vtype,
          generic_product_impl<votca::xtp::EwaldPeriodicDipoleOperator,
                               Vtype>> {
  typedef
      typename Product<votca::xtp::EwaldPeriodicDipoleOperator, Vtype>::Scalar
          Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst,
                            const votca::xtp::EwaldPeriodicDipoleOperator& op,
                            const Vtype& v, const Scalar& alpha) {
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);
    Eigen::VectorXd temp = op.multiply(v);
    dst = temp.cast<Scalar>();
  }
};
}  // namespace internal
}  // namespace Eigen

#endif  // VOTCA_XTP_EWALDPERIODICDIPOLEOPERATOR_H
