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
#include "ewaldrealspacesum.h"
#include "ewaldreciprocalspacesum.h"
#include "ewaldregistry.h"
#include "ewaldshapecorrection.h"

/**
 * \brief Matrix-free periodic induced-dipole self-consistent-field
 *        operator, for use with Eigen iterative solvers.
 *
 * Follows the Eigen::EigenBase interface of DipoleDipoleInteraction
 * (dipoledipoleinteraction.h), with the periodic real- plus
 * reciprocal-space field (EwaldRealSpaceSum, EwaldReciprocalSpaceSum,
 * EwaldShapeCorrection) in place of the aperiodic
 * eeInteractor::FillTholeInteraction.
 *
 * Vector layout: for each id in `ids` (in order), 3 entries (xyz) per site
 * of that segment, in the site order of EwaldRegistry::Get(id, Neutral).
 * Segments may have any number of sites, so the layout is variable-width
 * (see offsets_).
 *
 * The induced dipoles solve mu = P*(F_perm + C*mu), i.e.
 *   A*mu = (P^-1 - C)*mu = F_perm,
 * where C*mu is the field produced by the induced dipoles of the `ids`
 * sites: the periodic real-space sum (excluding zero-translation
 * same-segment pairs), the reciprocal-space sum (no exclusions), the
 * shape/surface term, and the intramolecular Thole-damped coupling
 * between distinct sites of the same segment (AddIntraSegmentCoupling).
 * Like legacy PolarBackground and DipoleDipoleInteraction, intramolecular
 * induced-induced coupling is included; it is Thole damped with the same
 * exponential (van Duijnen-Swart) form as legacy EwdInteractor::L3()/L5()
 * (FU12_ERFC_At_By damps intra- and intermolecular pairs alike).
 *
 * Because the reciprocal-space sum excludes nothing, it also contains
 * each site's spurious self-field -M*v_i, which is cancelled by
 * self_field_matrix_ (as in legacy FU12_ERF_At_By with identical sites).
 * With that correction getPInv() is the exact diagonal block of A.
 *
 * Sign convention: b = +F_perm, not b = -V as in
 * PolarRegion::CalcInducedDipolesViaPCG. eeInteractor::V() stores the
 * negative of the field (e = q*phi + mu.V), whereas the Ewald classes'
 * V() stores the physical field E (energy q*phi - mu.E).
 *
 * The field sums include every registered segment, not only `ids`, so
 * RawMultiply(v) contains a v-independent contribution from segments
 * outside `ids`. multiply(v) = RawMultiply(v) - baseline_, with
 * baseline_ = RawMultiply(0), which restores linearity (multiply(0) = 0).
 * Segments outside `ids` thus act as a fixed, non-responding background;
 * their influence enters only through the permanent field in b.
 *
 * The operator is symmetric when restricted to `ids` (see
 * test_ewaldperiodicdipolesymmetry.cc); the intramolecular term is added
 * as a block/block.transpose() pair to preserve this. The sign/scale
 * convention is tested on a single polarizable site in a point-charge
 * field (test_ewaldperiodicdipoleoperator.cc), multi-site segments in
 * test_ewaldperiodicdipoleoperatormultisite.cc.
 *
 * operator()(i,j) is only used by Eigen's preconditioner setup, which
 * keeps the diagonal. Same-site entries return the getPInv() block;
 * all other entries (including distinct sites of one segment) return
 * 0.0, since a periodic pair tensor would require a full Ewald
 * evaluation.
 *
 * registry, real_sum, recip_sum and shape are stored by reference and
 * must outlive this object; the sums must be built on the same registry.
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
// The traits must be specialized before the class derives from
// Eigen::EigenBase, which instantiates internal::traits of the derived
// type immediately (same ordering as dipoledipoleinteraction.h).
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

  // ids: the segments solved for; their order defines the vector layout
  //   (see class documentation). Every id must be registered at
  //   EwaldChargeState::Neutral.
  // alpha_ewald: Ewald splitting parameter, the same as used for
  //   real_sum/recip_sum; needed for the intramolecular interactor.
  // thole_a: Thole damping parameter, the same as used for real_sum
  //   (polarmethod.aDamp in ewdbgpol.xml); used by
  //   AddIntraSegmentCoupling.
  // shape: the EwaldShapeCorrection also used for the permanent field.
  //   Applied to the induced dipoles in every product, as legacy does
  //   with FU12_ShapeField_At_By in each induction iteration.
  EwaldPeriodicDipoleOperator(EwaldRegistry& registry,
                              const EwaldRealSpaceSum& real_sum,
                              const EwaldReciprocalSpaceSum& recip_sum,
                              const EwaldShapeCorrection& shape,
                              std::vector<Index> ids, double alpha_ewald,
                              double thole_a);

  // Accumulated wall-clock time (seconds) spent in RawMultiply, split by
  // phase, over the lifetime of this object. n_calls is the number of
  // RawMultiply calls (one per solver iteration plus one for baseline_).
  struct RawMultiplyTimings {
    double setup = 0.0;       // dipole set + Reset over every site
    double real_space = 0.0;  // EwaldRealSpaceSum::AddFieldAt
    double reciprocal = 0.0;  // EwaldReciprocalSpaceSum::AddFieldAtMany
    double shape = 0.0;       // EwaldShapeCorrection::AddFieldAt
    double assemble = 0.0;    // P^-1*v - V - M*v result assembly
    double intra = 0.0;       // AddIntraSegmentCoupling
    Index n_calls = 0;
    double total() const {
      return setup + real_space + reciprocal + shape + assemble + intra;
    }
  };
  const RawMultiplyTimings& Timings() const { return timings_; }

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

  // Same-site entries return the getPInv() block; cross-site entries
  // return 0.0 (only the diagonal is used by the preconditioner).
  double operator()(Index i, Index j) const;

  // Matrix-vector product A*v (see class documentation).
  Eigen::VectorXd multiply(const Eigen::VectorXd& v) const;

 private:
  // Computes getPInv()*v minus the field at every ids_ site from all
  // registered segments (including those not in ids_). Affine rather
  // than linear in v; multiply() subtracts baseline_ = RawMultiply(0).
  Eigen::VectorXd RawMultiply(const Eigen::VectorXd& v) const;

  // Adds the intramolecular Thole-damped induced-dipole field (the part
  // not already contained in the reciprocal-space sum) to result, for
  // every segment in ids_ with more than one site. Linear in v.
  void AddIntraSegmentCoupling(const Eigen::VectorXd& v,
                               Eigen::VectorXd& result) const;

  // Maps a global vector index to (segment index in ids_, site index in
  // that segment) by binary search on offsets_. Used by operator()(i,j).
  std::pair<Index, Index> LocateSite(Index i) const;

  EwaldRegistry& registry_;
  const EwaldRealSpaceSum& real_sum_;
  const EwaldReciprocalSpaceSum& recip_sum_;
  const EwaldShapeCorrection& shape_;
  std::vector<Index> ids_;
  // Provides ComputeB/ComputeErfB and ComputeThole for the
  // intramolecular term.
  EwaldRealSpaceInteractor intra_interactor_;
  // offsets_[n] is the scalar index where segment ids_[n] starts;
  // offsets_[n+1]-offsets_[n] = 3 * (number of sites in that segment).
  // offsets_.size() == ids_.size()+1, with offsets_.back() == size_.
  std::vector<Index> offsets_;
  Index size_;
  Eigen::VectorXd baseline_;
  mutable RawMultiplyTimings timings_;
  // 3x3 matrix M such that a dipole v_i produces a spurious
  // reciprocal-space self-field -M*v_i at its own position (see
  // EwaldReciprocalSpaceSum::SelfFieldMatrix). Depends only on the
  // k-vectors, alpha and volume. RawMultiply cancels this self-field for
  // every site, as legacy does via EwdInteractor::FU12_ERF_At_By with
  // the same site as both arguments.
  Eigen::Matrix3d self_field_matrix_;
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
