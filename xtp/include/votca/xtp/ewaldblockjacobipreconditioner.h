#ifndef VOTCA_XTP_EWALDBLOCKJACOBIPRECONDITIONER_H
#define VOTCA_XTP_EWALDBLOCKJACOBIPRECONDITIONER_H

// Standard includes
#include <optional>
#include <stdexcept>

// Third party includes
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>

// Local VOTCA includes
#include "ewaldregistry.h"
#include "ewaldrealspaceinteractor.h"

namespace votca {
namespace xtp {

/**
 * \brief Block-Jacobi preconditioner for EwaldPeriodicDipoleOperator's own
 *        PCG solve.
 *
 * USAGE: Eigen::ConjugateGradient owns its own preconditioner as a plain
 * (default-constructed) member (Eigen::IterativeSolverBase<>::
 * m_preconditioner), and cg.compute(op) calls that member's own
 * compute(op.matrix()) -- it never passes a caller-supplied preconditioner
 * instance through. Confirmed directly (a real, not hypothetical, failure
 * caught by an actual compiled-and-run test): calling
 * EwaldBlockJacobiPreconditioner(registry, ids, alpha, thole_a) and then
 * handing that instance to a template parameter alone is NOT enough --
 * cg.compute(op) still default-constructs its own (uninitialized) instance
 * over it. The construct that actually works, verified end to end via
 * Eigen::ConjugateGradient with a genuine block-diagonal test operator
 * (0 iterations to machine precision, since the preconditioner exactly
 * matched the operator there):
 *
 *   Eigen::ConjugateGradient<Op, Lower|Upper, EwaldBlockJacobiPreconditioner> cg;
 *   cg.compute(op);
 *   cg.preconditioner() = EwaldBlockJacobiPreconditioner(registry, ids,
 *                                                        alpha_ewald,
 *                                                        thole_a);
 *   x = cg.solveWithGuess(b, x0);
 *
 * i.e. assign the real, already-built preconditioner via
 * IterativeSolverBase<>::preconditioner()'s own public, non-const
 * accessor, AFTER compute(), not by relying on compute() to build it.
 *
 * This exists because Eigen::DiagonalPreconditioner<double> -- the default
 * used with EwaldPeriodicDipoleOperator -- reads only operator()(i,j)'s own
 * same-site diagonal (getPInv()), with zero knowledge of the intramolecular
 * (Thole-damped) coupling EwaldPeriodicDipoleOperator::
 * AddIntraSegmentCoupling actually adds for multi-site segments. That gap
 * was found, this session, to matter a great deal in practice: enabling
 * Thole damping on the intramolecular term (the physically correct,
 * legacy-matching behavior -- confirmed directly against legacy's own
 * FU12_ERFC_At_By/UpdateAllBls) made a real, large ~5000-site PCG solve
 * diverge outright (residual growing past its own starting value, not
 * merely converging slowly) with the bare-diagonal preconditioner, while
 * disabling that same damping (reverting to the undamped erfc-only tensor)
 * let the identical system converge cleanly in 51 iterations. Since the
 * physics itself is confirmed correct (matching legacy), the fix pursued
 * here is a better preconditioner, not reverting the physics.
 *
 * A direct numerical check (methane geometry, realistic atomic
 * polarizabilities, thole_a=0.6) found the intramolecular block ALONE,
 * isolated from the rest of the system, IS positive-definite when
 * Thole-damped (min eigenvalue ~0.05, comfortably positive) and is NOT
 * when undamped (min eigenvalue ~-0.12) -- the opposite of what a first,
 * trace-based argument suggested (a real correction made mid-session: a
 * negative trace does not by itself imply a negative eigenvalue exists,
 * and here it doesn't). This means the observed divergence is NOT because
 * the local, per-segment intramolecular coupling is itself indefinite --
 * whatever is happening lives either in how that locally-SPD block
 * combines with the periodic (intermolecular, shape, self-field) terms
 * across the full ~5000-site system, or is purely a preconditioner-
 * tracking problem (the true global operator staying SPD throughout, but
 * DiagonalPreconditioner's blindness to the intramolecular term causing
 * CG's implicit preconditioned-residual bookkeeping to behave as though
 * it doesn't). This class was not built having resolved which of those is
 * true -- only having confirmed the local block itself is a valid (SPD)
 * building block either way, and that a preconditioner genuinely
 * reflecting the intramolecular coupling should track the true operator
 * far better than the bare diagonal does regardless of which explanation
 * is correct. Whether this actually resolves the real ~5000-site
 * divergence, as opposed to just being well-motivated, has not yet been
 * tested against that real system -- only the mechanics above have been
 * verified, with a synthetic operator, not the real EwaldPeriodicDipole
 * Operator/real molecular system.
 *
 * What this computes: for every segment in ids with 2+ sites, the full
 * dense 3*N_site x 3*N_site local block -- getPInv() on each site's own
 * diagonal, PLUS the same Thole-damped intramolecular coupling
 * AddIntraSegmentCoupling itself adds (identical B-function/ComputeThole
 * calls, same sign convention) -- factorized ONCE via LDLT (chosen over
 * LLT/Cholesky specifically because the local block was found to be only
 * marginally SPD at higher thole_a, e.g. min eigenvalue ~0.0016 at
 * thole_a=1.0 in the same methane test above -- LDLT degrades gracefully
 * on a near-singular or, for some future geometry/polarizability
 * combination, genuinely indefinite block, where LLT would simply fail).
 * Single-site segments fall back to the plain diagonal (getPInv()),
 * identical to what DiagonalPreconditioner itself would do for them.
 *
 * solve(b) applies each segment's own precomputed factorization to that
 * segment's own slice of b independently -- a genuine block-Jacobi
 * scheme: exact within each segment's own local block, still ignoring
 * inter-segment (periodic, real+reciprocal-space) coupling entirely, the
 * same simplification DiagonalPreconditioner itself already makes at the
 * single-site level. The intermolecular coupling is left to PCG's own
 * outer iteration to resolve, same as before -- only the intramolecular
 * piece moves from "completely ignored" to "handled exactly."
 */
class EwaldBlockJacobiPreconditioner {
  using Scalar = double;
  using Vector = Eigen::VectorXd;

 public:
  using StorageIndex = votca::Index;
  enum { ColsAtCompileTime = Eigen::Dynamic, MaxColsAtCompileTime = Eigen::Dynamic };

  EwaldBlockJacobiPreconditioner() : is_initialized_(false) {}

  // registry/ids/alpha_ewald/thole_a mirror EwaldPeriodicDipoleOperator's
  // own constructor exactly -- same segments, same layout, same
  // intramolecular interactor parameters -- so a caller already holding
  // those values for the operator can pass the identical ones here.
  EwaldBlockJacobiPreconditioner(EwaldRegistry& registry,
                                 std::vector<Index> ids, double alpha_ewald,
                                 double thole_a)
      : registry_(&registry), ids_(std::move(ids)) {
    intra_interactor_.emplace(alpha_ewald, thole_a);
    BuildOffsets();
    FactorizeBlocks();
    is_initialized_ = true;
  }

  Index rows() const { return size_; }
  Index cols() const { return size_; }

  // Matches DiagonalPreconditioner's own analyzePattern/factorize/compute
  // split: this class's own local blocks depend only on registry_/ids_/
  // intra_interactor_ (fixed at construction, never on the operator
  // matrix passed in here), so all of these are no-ops that just
  // return *this -- the real work already happened in the constructor.
  // Present only because Eigen::ConjugateGradient's own compute() path
  // requires them to exist on whatever preconditioner type it's given.
  template <typename MatType>
  EwaldBlockJacobiPreconditioner& analyzePattern(const MatType&) {
    return *this;
  }
  template <typename MatType>
  EwaldBlockJacobiPreconditioner& factorize(const MatType&) {
    return *this;
  }
  template <typename MatType>
  EwaldBlockJacobiPreconditioner& compute(const MatType&) {
    return *this;
  }

  template <typename Rhs, typename Dest>
  void _solve_impl(const Rhs& b, Dest& x) const {
    for (std::size_t n = 0; n < ids_.size(); ++n) {
      const Index base = offsets_[n];
      const Index width = offsets_[n + 1] - base;
      // LDLT::solve works identically regardless of block size -- no
      // special case needed for single-site (width==3) segments; their
      // own factorizations_[n] already holds a plain 3x3 LDLT of
      // getPInv() alone (see FactorizeBlocks), the same result
      // DiagonalPreconditioner itself would give them.
      x.segment(base, width) = factorizations_[n].solve(b.segment(base, width));
    }
  }

  template <typename Rhs>
  inline const Eigen::Solve<EwaldBlockJacobiPreconditioner, Rhs> solve(
      const Eigen::MatrixBase<Rhs>& b) const {
    eigen_assert(is_initialized_ &&
                "EwaldBlockJacobiPreconditioner is not initialized.");
    eigen_assert(size_ == b.rows() &&
                "EwaldBlockJacobiPreconditioner::solve(): invalid number "
                "of rows of the right hand side matrix b");
    return Eigen::Solve<EwaldBlockJacobiPreconditioner, Rhs>(*this,
                                                             b.derived());
  }

  Eigen::ComputationInfo info() const { return Eigen::Success; }

 private:
  void BuildOffsets() {
    offsets_.reserve(ids_.size() + 1);
    offsets_.push_back(0);
    for (Index id : ids_) {
      if (!registry_->Has(id, EwaldChargeState::Neutral)) {
        throw std::runtime_error(
            "EwaldBlockJacobiPreconditioner: segment id not registered at "
            "EwaldChargeState::Neutral");
      }
      Index n_sites = registry_->Get(id, EwaldChargeState::Neutral).size();
      offsets_.push_back(offsets_.back() + 3 * n_sites);
    }
    size_ = offsets_.back();
  }

  // Builds and factorizes every segment's own local block, matching
  // EwaldPeriodicDipoleOperator's own intramolecular term: the same
  // B-function/ComputeThole calls, and the same sign it enters the
  // OPERATOR with -- i.e. subtracted, since RawMultiply forms
  // A = P^-1 - C. Note this is the opposite sign to
  // AddIntraSegmentCoupling's own internal accumulation, which builds
  // the coupling field itself (a positive quantity) that RawMultiply
  // then subtracts. An earlier version of this comment claimed to
  // mirror that method's sign convention "exactly", and the code did --
  // which is precisely why this was wrong, and why it survived the
  // operator's own sign fix without being noticed. Compare against the
  // operator's assembled block, not against AddIntraSegmentCoupling in
  // isolation.
  void FactorizeBlocks() {
    factorizations_.reserve(ids_.size());
    for (std::size_t n = 0; n < ids_.size(); ++n) {
      const PolarSegment& segment =
          registry_->Get(ids_[n], EwaldChargeState::Neutral);
      const Index n_sites = segment.size();
      const Index width = 3 * n_sites;
      Eigen::MatrixXd block = Eigen::MatrixXd::Zero(width, width);

      for (Index s = 0; s < n_sites; ++s) {
        block.block<3, 3>(3 * s, 3 * s) = segment[s].getPInv();
      }

      if (n_sites >= 2) {
        for (Index i = 0; i < n_sites; ++i) {
          const PolarSite& site_i = segment[i];
          for (Index j = i + 1; j < n_sites; ++j) {
            const PolarSite& site_j = segment[j];
            const Eigen::Vector3d r_vec = site_i.getPos() - site_j.getPos();
            const double r = r_vec.norm();
            const EwaldRealSpaceInteractor::BFunctions b =
                intra_interactor_->ComputeB(r);
            const EwaldRealSpaceInteractor::TholeFactors t =
                intra_interactor_->ComputeThole(r, site_j, site_i);
            const Eigen::Matrix3d coupling =
                t.l5 * b.B2 * (r_vec * r_vec.transpose()) -
                t.l3 * b.B1 * Eigen::Matrix3d::Identity();
            // BUG FIX (this session): MINUS, not plus. The operator
            // being preconditioned is A = P^-1 - C: RawMultiply builds
            // the intramolecular term with exactly the `coupling`
            // expression above and then subtracts it (result -= intra).
            // This class assembled it with a plus, so it was
            // factorizing P^-1 + C_intra -- the operator as it stood
            // BEFORE the coupling-sign fix, which this file was written
            // against and which never propagated here.
            //
            // The consequence was not a wrong answer (a preconditioner
            // cannot change the fixed point, only the path to it) but a
            // markedly worse one: on a real 5000-site solve this took
            // PCG from 16 iterations to 29, with a visibly
            // non-monotonic p.A.p curvature trace, versus a smoothly
            // decreasing one unpreconditioned.
            block.block<3, 3>(3 * i, 3 * j) -= coupling;
            block.block<3, 3>(3 * j, 3 * i) -= coupling.transpose();
          }
        }
      }

      Eigen::LDLT<Eigen::MatrixXd> ldlt(block);
      factorizations_.push_back(std::move(ldlt));
    }
  }

  EwaldRegistry* registry_ = nullptr;
  std::vector<Index> ids_;
  std::optional<EwaldRealSpaceInteractor> intra_interactor_;
  std::vector<Index> offsets_;
  Index size_ = 0;
  std::vector<Eigen::LDLT<Eigen::MatrixXd>> factorizations_;
  bool is_initialized_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EWALDBLOCKJACOBIPRECONDITIONER_H
