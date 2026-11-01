#ifndef VOTCA_XTP_EWALDBLOCKJACOBIPRECONDITIONER_H
#define VOTCA_XTP_EWALDBLOCKJACOBIPRECONDITIONER_H

// Standard includes
#include <optional>
#include <stdexcept>

// Third party includes
#include <Eigen/Cholesky>
#include <Eigen/IterativeLinearSolvers>

// Local VOTCA includes
#include "ewaldrealspaceinteractor.h"
#include "ewaldregistry.h"

namespace votca {
namespace xtp {

/**
 * \brief Block-Jacobi preconditioner for the PCG solve of
 *        EwaldPeriodicDipoleOperator.
 *
 * Usage: Eigen::ConjugateGradient default-constructs its preconditioner
 * member and cg.compute(op) only calls that member's compute(), so a
 * constructed instance must be assigned after compute():
 *
 *   Eigen::ConjugateGradient<Op, Lower|Upper, EwaldBlockJacobiPreconditioner>
 * cg; cg.compute(op); cg.preconditioner() =
 * EwaldBlockJacobiPreconditioner(registry, ids, alpha_ewald, thole_a); x =
 * cg.solveWithGuess(b, x0);
 *
 * Eigen::DiagonalPreconditioner only sees the same-site getPInv() blocks
 * and ignores the Thole-damped intramolecular coupling of
 * EwaldPeriodicDipoleOperator::AddIntraSegmentCoupling. With that
 * preconditioner a ~5000-site PCG solve diverged once intramolecular
 * Thole damping was enabled. The isolated intramolecular block itself is
 * SPD when Thole damped (methane, thole_a = 0.6: min eigenvalue ~0.05;
 * undamped ~-0.12), so it is a valid preconditioner block.
 *
 * For every segment in ids the dense 3*N_site x 3*N_site block of
 * getPInv() on the diagonal minus the intramolecular coupling is built
 * and factorized once with LDLT. LDLT rather than LLT because the block
 * can be nearly singular (min eigenvalue ~0.0016 at thole_a = 1.0 for
 * methane). Single-site segments reduce to the plain diagonal block.
 *
 * solve(b) applies each segment's factorization to its slice of b
 * (block Jacobi): intramolecular coupling is treated exactly,
 * intermolecular (periodic) coupling is left to the outer PCG iteration.
 */
class EwaldBlockJacobiPreconditioner {
  using Scalar = double;
  using Vector = Eigen::VectorXd;

 public:
  using StorageIndex = votca::Index;
  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic
  };

  EwaldBlockJacobiPreconditioner() : is_initialized_(false) {}

  // registry/ids/alpha_ewald/thole_a must be the same as those passed to
  // EwaldPeriodicDipoleOperator (same segments, layout and interactor).
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

  // Required by the Eigen preconditioner interface. No-ops: the blocks
  // depend only on registry_/ids_/intra_interactor_ and are factorized
  // in the constructor.
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
    // Segments own disjoint slices of x, so they are independent.
#pragma omp parallel for schedule(dynamic, 4)
    for (Index nn = 0; nn < Index(ids_.size()); ++nn) {
      const std::size_t n = std::size_t(nn);
      const Index base = offsets_[n];
      const Index width = offsets_[n + 1] - base;
      // Single-site segments hold a 3x3 LDLT of getPInv(), so no special
      // case is needed.
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

  // Builds and factorizes the local block of every segment. The
  // intramolecular coupling enters with the sign it has in the operator
  // A = P^-1 - C (subtracted), i.e. opposite to the field accumulated by
  // AddIntraSegmentCoupling.
  void FactorizeBlocks() {
    // Blocks are independent and factorized in parallel.
    factorizations_.resize(ids_.size());
#pragma omp parallel for schedule(dynamic, 4)
    for (Index nn = 0; nn < Index(ids_.size()); ++nn) {
      const std::size_t n = std::size_t(nn);
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
            // The operator couples an intramolecular pair through its
            // explicit term l*B_erfc + (l-1)*B_erf plus the B_erf part of
            // the reciprocal sum, in total the Thole-damped bare coupling
            // l*(B_erfc + B_erf), which is used here.
            const EwaldRealSpaceInteractor::BFunctions b =
                intra_interactor_->ComputeB(r);
            const EwaldRealSpaceInteractor::BFunctions berf =
                intra_interactor_->ComputeErfB(r);
            const EwaldRealSpaceInteractor::TholeFactors t =
                intra_interactor_->ComputeThole(r, site_j, site_i);
            const Eigen::Matrix3d coupling =
                t.l5 * (b.B2 + berf.B2) * (r_vec * r_vec.transpose()) -
                t.l3 * (b.B1 + berf.B1) * Eigen::Matrix3d::Identity();
            // Subtracted, as in A = P^-1 - C. With the wrong sign
            // (P^-1 + C_intra) a 5000-site solve needed 29 instead of 16
            // PCG iterations.
            block.block<3, 3>(3 * i, 3 * j) -= coupling;
            block.block<3, 3>(3 * j, 3 * i) -= coupling.transpose();
          }
        }
      }

      factorizations_[n].compute(block);
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
