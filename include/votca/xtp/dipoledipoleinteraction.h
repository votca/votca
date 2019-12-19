/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#ifndef __VOTCA_DIPOLEDIPOLEINTERACTION_H
#define __VOTCA_DIPOLEDIPOLEINTERACTION_H
#include <votca/xtp/eeinteractor.h>
#include <votca/xtp/eigen.h>
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
      : _interactor(interactor) {
    _size = 0;
    for (const PolarSegment seg : segs) {
      _size += 3 * seg.size();
    }
    _sites.reserve(_size / 3);
    for (const PolarSegment& seg : segs) {
      for (const PolarSite& site : seg) {
        _sites.push_back(&site);
      }
    }
  }

  class InnerIterator {
   public:
    InnerIterator(const DipoleDipoleInteraction& xpr, const Index& id)
        : _xpr(xpr), _id(id){};

    InnerIterator& operator++() {
      _row++;
      return *this;
    }
    operator bool() const {
      return _row < _xpr._size;
    }  // DO not use the size method, it returns linear dimension*linear
       // dimension i.e _size^2
    double value() const { return _xpr(_row, _id); }
    Index row() const { return _row; }
    Index col() const { return _id; }
    Index index() const { return row(); }

   private:
    const DipoleDipoleInteraction& _xpr;
    const Index _id;
    Index _row = 0;
  };

  Index rows() const { return this->_size; }
  Index cols() const { return this->_size; }
  Index outerSize() const { return this->_size; }

  template <typename Vtype>
  Eigen::Product<DipoleDipoleInteraction, Vtype, Eigen::AliasFreeProduct>
      operator*(const Eigen::MatrixBase<Vtype>& x) const {
    return Eigen::Product<DipoleDipoleInteraction, Vtype,
                          Eigen::AliasFreeProduct>(*this, x.derived());
  }

  // this is not a fast method
  const double& operator()(const Index i, const Index j) const {
    Index seg1id = Index(i / 3);
    Index xyz1 = Index(i % 3);
    Index seg2id = Index(j / 3);
    Index xyz2 = Index(j % 3);

    if (seg1id == seg2id) {
      return _sites[seg1id]->getPInv()(xyz1, xyz2);
    } else {
      const PolarSite& site1 = *_sites[seg1id];
      const PolarSite& site2 = *_sites[seg2id];
      return _interactor.FillTholeInteraction(site1, site2)(xyz1, xyz2);
    }
  };

  Eigen::VectorXd multiply(const Eigen::VectorXd& v) const {
    assert(v.size() == _size &&
           "input vector has the wrong size for multiply with operator");
    const Index segment_size = Index(_sites.size());
    std::vector<Eigen::VectorXd> thread_result(OPENMP::getMaxThreads(),
                                               Eigen::VectorXd::Zero(_size));
#pragma omp parallel for schedule(dynamic)
    for (Index i = 0; i < segment_size; i++) {
      const PolarSite& site1 = *_sites[i];
      for (Index j = i + 1; j < segment_size; j++) {
        const PolarSite& site2 = *_sites[j];
        Eigen::Matrix3d block = _interactor.FillTholeInteraction(site1, site2);
        thread_result[OPENMP::getThreadId()].segment<3>(3 * i) +=
            block * v.segment<3>(3 * j);
        thread_result[OPENMP::getThreadId()].segment<3>(3 * j) +=
            block.transpose() * v.segment<3>(3 * i);
      }
    }
#pragma omp parallel for schedule(dynamic)
    for (Index i = 0; i < segment_size; i++) {
      const PolarSite& site = *_sites[i];
      thread_result[OPENMP::getThreadId()].segment<3>(3 * i) +=
          site.getPInv() * v.segment<3>(3 * i);
    }

    return std::accumulate(thread_result.begin(), thread_result.end(),
                           Eigen::VectorXd::Zero(_size).eval());
  }

 private:
  const eeInteractor& _interactor;
  std::vector<const PolarSite*> _sites;
  Index _size;
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
    dst = op.multiply(v);
  }
};

}  // namespace internal
}  // namespace Eigen

#endif  //__VOTCA_DIPOLEDIPOLEINTERACTION_H