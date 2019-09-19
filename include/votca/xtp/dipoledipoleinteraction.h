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
  typedef double Scalar;
  typedef double RealScalar;
  typedef int StorageIndex;
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
  const double& operator()(const size_t i, const size_t j) const {
    int seg1id = i / 3;
    int xyz1 = i % 3;
    int seg2id = j / 3;
    int xyz2 = j % 3;

    if (seg1id == seg2id) {
      return _sites[seg1id]->getPInv()(xyz1, xyz2);
    } else {
      const PolarSite& site1 = *_sites[seg1id];
      const PolarSite& site2 = *_sites[seg2id];
      return _interactor.FillTholeInteraction(site1, site2)(xyz1, xyz2);
    }
  };

  Eigen::Vector3d Block(int i, const Eigen::VectorXd& v) const {
    Eigen::Vector3d result = Eigen::Vector3d::Zero();
    const PolarSite& site1 = *_sites[i];
    const int segment_size = _sites.size();
    for (int j = 0; j < segment_size; j++) {
      const Eigen::Vector3d v_small = v.segment<3>(3 * j);

      if (i == j) {
        result += site1.getPInv() * v_small;
      } else {
        const PolarSite& site2 = *_sites[j];
        result += _interactor.VThole(site1, site2, v_small);
      }
    }
    return result;
  }

 private:
  const eeInteractor& _interactor;
  std::vector<const PolarSite*> _sites;
  int _size;
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
    int sites = op.rows() / 3;
// make the mat vect product
#pragma omp parallel for
    for (int i = 0; i < sites; i++) {
      const Eigen::Vector3d result = op.Block(i, v);
      dst.segment(3 * i, 3) = result;
    }
  }
};

}  // namespace internal
}  // namespace Eigen

#endif  //__VOTCA_DIPOLEDIPOLEINTERACTION_H