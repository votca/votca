/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_DIPOLEDIPOLEINTERACTION_H
#define VOTCA_XTP_DIPOLEDIPOLEINTERACTION_H

// Local VOTCA includes
#include "eeinteractor.h"
#include "eigen.h"

namespace votca {
namespace xtp {
class DipoleDipoleInteraction;
}
}  // namespace votca
namespace Eigen {
namespace internal {
// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
template <>
struct traits<votca::xtp::DipoleDipoleInteraction>
    : public Eigen::internal::traits<Eigen::MatrixXd> {};
}  // namespace internal
}  // namespace Eigen

namespace votca {
namespace xtp {

class DipoleDipoleInteraction
    : public Eigen::EigenBase<DipoleDipoleInteraction> {
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

  DipoleDipoleInteraction(const eeInteractor& interactor,
                          const std::vector<PolarSegment>& segs)
      : interactor_(interactor) {
    size_ = 0;
    for (const PolarSegment& seg : segs) {
      size_ += 3 * seg.size();
    }
    sites_.reserve(size_ / 3);
    for (const PolarSegment& seg : segs) {
      for (const PolarSite& site : seg) {
        sites_.push_back(&site);
      }
    }
  }

  class InnerIterator {
   public:
    InnerIterator(const DipoleDipoleInteraction& xpr, const Index& id)
        : xpr_(xpr), id_(id){};

    InnerIterator& operator++() {
      row_++;
      return *this;
    }
    operator bool() const {
      return row_ < xpr_.size_;
    }  // DO not use the size method, it returns linear dimension*linear
       // dimension i.e  size_^2
    double value() const { return xpr_(row_, id_); }
    Index row() const { return row_; }
    Index col() const { return id_; }
    Index index() const { return row(); }

   private:
    const DipoleDipoleInteraction& xpr_;
    const Index id_;
    Index row_ = 0;
  };

  Index rows() const { return this->size_; }
  Index cols() const { return this->size_; }
  Index outerSize() const { return this->size_; }

  template <typename Vtype>
  Eigen::Product<DipoleDipoleInteraction, Vtype, Eigen::AliasFreeProduct>
      operator*(const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<DipoleDipoleInteraction, Vtype,
                          Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // this is not a fast method
  double operator()(const Index i, const Index j) const {
    Index seg1id = Index(i / 3);
    Index xyz1 = Index(i % 3);
    Index seg2id = Index(j / 3);
    Index xyz2 = Index(j % 3);

    if (seg1id == seg2id) {
      return sites_[seg1id]->getPInv()(xyz1, xyz2);
    } else {
      const PolarSite& site1 = *sites_[seg1id];
      const PolarSite& site2 = *sites_[seg2id];
      return interactor_.FillTholeInteraction(site1, site2)(xyz1, xyz2);
    }
  };

  Eigen::VectorXd multiply(const Eigen::VectorXd& v) const {
    assert(v.size() == size_ &&
           "input vector has the wrong size for multiply with operator");
    const Index segment_size = Index(sites_.size());
    Eigen::VectorXd result = Eigen::VectorXd::Zero(size_);
#pragma omp parallel for schedule(dynamic) reduction(+ : result)
    for (Index i = 0; i < segment_size; i++) {
      const PolarSite& site1 = *sites_[i];
      result.segment<3>(3 * i) += site1.getPInv() * v.segment<3>(3 * i);
      for (Index j = i + 1; j < segment_size; j++) {
        const PolarSite& site2 = *sites_[j];
        Eigen::Matrix3d block = interactor_.FillTholeInteraction(site1, site2);
        result.segment<3>(3 * i) += block * v.segment<3>(3 * j);
        result.segment<3>(3 * j) += block.transpose() * v.segment<3>(3 * i);
      }
    }

    return result;
  }

 private:
  const eeInteractor& interactor_;
  std::vector<const PolarSite*> sites_;
  Index size_;
};
}  // namespace xtp
}  // namespace votca

namespace Eigen {

namespace internal {

// replacement of the mat*vect operation
template <typename Vtype>
struct generic_product_impl<votca::xtp::DipoleDipoleInteraction, Vtype,
                            DenseShape, DenseShape, GemvProduct>
    : generic_product_impl_base<
          votca::xtp::DipoleDipoleInteraction, Vtype,
          generic_product_impl<votca::xtp::DipoleDipoleInteraction, Vtype>> {

  typedef typename Product<votca::xtp::DipoleDipoleInteraction, Vtype>::Scalar
      Scalar;

  template <typename Dest>
  static void scaleAndAddTo(Dest& dst,
                            const votca::xtp::DipoleDipoleInteraction& op,
                            const Vtype& v, const Scalar& alpha) {
    // returns dst = alpha * op * v
    // alpha must be 1 here
    assert(alpha == Scalar(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);
    Eigen::VectorXd temp = op.multiply(v);
    dst = temp.cast<Scalar>();  // tumbleweed fix do not delete
  }
};

}  // namespace internal
}  // namespace Eigen

#endif  // VOTCA_XTP_DIPOLEDIPOLEINTERACTION_H
