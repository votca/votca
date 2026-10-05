/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Local VOTCA includes
#include "votca/xtp/bse_operator.h"
#include "votca/xtp/bse_operator_kernels.h"
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR<cqp, cx, cd, cd2>::configure(BSEOperator_Options opt) {
  opt_ = opt;
  Index bse_vmax = opt_.homo;
  bse_cmin_ = opt_.homo + 1;
  bse_vtotal_ = bse_vmax - opt_.vmin + 1;
  bse_ctotal_ = opt_.cmax - bse_cmin_ + 1;
  bse_size_ = bse_vtotal_ * bse_ctotal_;
  this->set_size(bse_size_);
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::MatrixXd BSE_OPERATOR<cqp, cx, cd, cd2>::matmul(
    const Eigen::MatrixXd& input) const {

  static_assert(!(cd2 != 0 && cd != 0),
                "Hamiltonian cannot contain Hd and Hd2 at the same time");

  Index auxsize = Mmn_.auxsize();
  vc2index vc = vc2index(0, 0, bse_ctotal_);

  Index vmin = opt_.vmin - opt_.rpamin;
  Index cmin = bse_cmin_ - opt_.rpamin;

  Eigen::MatrixXd result;
  if ((cd != 0 || cd2 != 0) && !direct_built_ && direct_cache_limit_ > 0.0 &&
      8.0 * double(bse_size_) * double(bse_size_) <= direct_cache_limit_) {
    BuildDirectMatrix();
  }
  if (direct_built_) {
    result = direct_ * input;
  } else if (cd != 0 || cd2 != 0) {
    // screened direct terms, one Hamiltonian row per (v1,c1)
    OpenMP_CUDA transform;
    if (cd != 0) {
      transform.createTemporaries(epsilon_0_inv_, input, bse_ctotal_,
                                  bse_vtotal_, auxsize);
    } else {
      transform.createTemporaries(epsilon_0_inv_, input, bse_vtotal_,
                                  bse_ctotal_, auxsize);
    }

#pragma omp parallel
    {
      Index threadid = OPENMP::getThreadId();
#pragma omp for schedule(dynamic)
      for (Index c1 = 0; c1 < bse_ctotal_; c1++) {

        // Temp matrix has to stay in this scope, because transform only
        // holds a reference to it
        Eigen::MatrixXd Temp;
        if (cd != 0) {
          Temp = -cd * (Mmn_[c1 + cmin].middleRows(cmin, bse_ctotal_));
          transform.PrepareMatrix1(Temp, threadid);
        } else if (cd2 != 0) {
          Temp = -cd2 * (Mmn_[c1 + cmin].middleRows(vmin, bse_vtotal_));
          transform.PrepareMatrix1(Temp, threadid);
        }

        for (Index v1 = 0; v1 < bse_vtotal_; v1++) {
          transform.SetTempZero(threadid);
          if (cd != 0) {
            transform.PrepareMatrix2(
                Mmn_[v1 + vmin].middleRows(vmin, bse_vtotal_), cd2 != 0,
                threadid);
          }
          if (cd2 != 0) {
            transform.PrepareMatrix2(
                Mmn_[v1 + vmin].middleRows(cmin, bse_ctotal_), cd2 != 0,
                threadid);
          }
          transform.MultiplyRow(vc.I(v1, c1), threadid);
        }
      }
    }
    result = transform.getReductionVar();
  } else {
    result = Eigen::MatrixXd::Zero(bse_size_, input.cols());
  }

  if (cqp != 0) {
    result += ApplyBSEQuasiparticle(Hqp_, bse_vtotal_, bse_ctotal_, input,
                                    static_cast<double>(cqp));
  }
  if (cx != 0) {
    const BSEWindow window{vmin, cmin, bse_vtotal_, bse_ctotal_};
    result += ApplyBSEExchange(Mmn_, window, Mmn_, window, input,
                               static_cast<double>(cx));
  }
  return result;
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR<cqp, cx, cd, cd2>::BuildDirectMatrix() const {
  const Index vmin = opt_.vmin - opt_.rpamin;
  const Index cmin = bse_cmin_ - opt_.rpamin;
  vc2index vc = vc2index(0, 0, bse_ctotal_);
  direct_.resize(bse_size_, bse_size_);
  // Same contractions as the row loop in matmul; rows are contiguous in the
  // row-major matrix.
#pragma omp parallel for schedule(dynamic)
  for (Index c1 = 0; c1 < bse_ctotal_; c1++) {
    Eigen::MatrixXd Temp;
    if (cd != 0) {
      Temp = -double(cd) * (Mmn_[c1 + cmin].middleRows(cmin, bse_ctotal_)) *
             epsilon_0_inv_.asDiagonal();
    } else {
      Temp = -double(cd2) * (Mmn_[c1 + cmin].middleRows(vmin, bse_vtotal_)) *
             epsilon_0_inv_.asDiagonal();
    }
    for (Index v1 = 0; v1 < bse_vtotal_; v1++) {
      Eigen::Map<Eigen::MatrixXd> row(direct_.row(vc.I(v1, c1)).data(),
                                      bse_ctotal_, bse_vtotal_);
      if (cd != 0) {
        // (c2, v2) block of row (v1, c1)
        row.noalias() =
            Temp * Mmn_[v1 + vmin].middleRows(vmin, bse_vtotal_).transpose();
      } else {
        row.noalias() =
            Mmn_[v1 + vmin].middleRows(cmin, bse_ctotal_) * Temp.transpose();
      }
    }
  }
  direct_built_ = true;
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::diagonal() const {

  static_assert(!(cd2 != 0 && cd != 0),
                "Hamiltonian cannot contain Hd and Hd2 at the same time");

  vc2index vc = vc2index(0, 0, bse_ctotal_);
  Index vmin = opt_.vmin - opt_.rpamin;
  Index cmin = bse_cmin_ - opt_.rpamin;

  Eigen::VectorXd result = Eigen::VectorXd::Zero(bse_size_);

#pragma omp parallel for schedule(dynamic) reduction(+ : result)
  for (Index v = 0; v < bse_vtotal_; v++) {
    for (Index c = 0; c < bse_ctotal_; c++) {

      double entry = 0.0;
      if (cx != 0) {
        entry += cx * Mmn_[v + vmin].row(cmin + c).squaredNorm();
      }

      if (cqp != 0) {
        Index cmin_qp = bse_vtotal_;
        entry += cqp * (Hqp_(c + cmin_qp, c + cmin_qp) - Hqp_(v, v));
      }
      if (cd != 0) {
        entry -=
            cd * (Mmn_[c + cmin].row(c + cmin) * epsilon_0_inv_.asDiagonal() *
                  Mmn_[v + vmin].row(v + vmin).transpose())
                     .value();
      }
      if (cd2 != 0) {
        entry -=
            cd2 * (Mmn_[c + cmin].row(v + vmin) * epsilon_0_inv_.asDiagonal() *
                   Mmn_[v + vmin].row(c + cmin).transpose())
                      .value();
      }

      result(vc.I(v, c)) = entry;
    }
  }
  return result;
}

template class BSE_OPERATOR<1, 2, 1, 0>;
template class BSE_OPERATOR<1, 0, 1, 0>;

template class BSE_OPERATOR<1, 0, 0, 0>;
template class BSE_OPERATOR<0, 1, 0, 0>;
template class BSE_OPERATOR<0, 0, 1, 0>;
template class BSE_OPERATOR<0, 0, 0, 1>;

template class BSE_OPERATOR<0, 2, 0, 1>;

template class BSE_OPERATOR<1, 1, 1, 0>;
template class BSE_OPERATOR<0, 1, 0, 1>;

}  // namespace xtp

}  // namespace votca
