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
    result = row_kernel_ ? ApplyDirectRows(input) : ApplyDirectBlocks(input);
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
Eigen::MatrixXd BSE_OPERATOR<cqp, cx, cd, cd2>::ApplyDirectRows(
    const Eigen::MatrixXd& input) const {
  vc2index vc = vc2index(0, 0, bse_ctotal_);
  Index vmin = opt_.vmin - opt_.rpamin;
  Index cmin = bse_cmin_ - opt_.rpamin;
  // screened direct terms, one Hamiltonian row per (v1,c1):
  //   cd:  K(c2,v2) = sum_P S(c2,P) M_{v1}(v2,P), S = -cd  M_{c1}(c,:) diag(w)
  //   cd2: K(c2,v2) = sum_P M_{v1}(c2,P) S(v2,P), S = -cd2 M_{c1}(v,:) diag(w)
  // with K flattened column-major as the row's (v2,c2) index v2*n_c + c2.
  Eigen::MatrixXd result(bse_size_, input.cols());
#pragma omp parallel for schedule(dynamic)
  for (Index c1 = 0; c1 < bse_ctotal_; c1++) {
    const Eigen::MatrixXd S =
        (cd != 0)
            ? Eigen::MatrixXd(-double(cd) *
                              Mmn_[c1 + cmin].middleRows(cmin, bse_ctotal_) *
                              epsilon_0_inv_.asDiagonal())
            : Eigen::MatrixXd(-double(cd2) *
                              Mmn_[c1 + cmin].middleRows(vmin, bse_vtotal_) *
                              epsilon_0_inv_.asDiagonal());
    Eigen::MatrixXd K(bse_ctotal_, bse_vtotal_);
    for (Index v1 = 0; v1 < bse_vtotal_; v1++) {
      if (cd != 0) {
        K.noalias() =
            S * Mmn_[v1 + vmin].middleRows(vmin, bse_vtotal_).transpose();
      } else {
        K.noalias() =
            Mmn_[v1 + vmin].middleRows(cmin, bse_ctotal_) * S.transpose();
      }
      const Eigen::Map<const Eigen::VectorXd> row(K.data(), K.size());
      result.row(vc.I(v1, c1)).noalias() = row.transpose() * input;
    }
  }
  return result;
}

// The screened direct term in blocks of Hamiltonian rows,
//   cd:  K(v1c1, v2c2) = -cd  sum_P M_{c1}(c2,P) w_P M_{v1}(v2,P)
//   cd2: K(v1c1, v2c2) = -cd2 sum_P M_{c1}(v2,P) w_P M_{v1}(c2,P)
// cd:  for each c1, the rows (v1, c1) of all v1 as one (n_v x N) matrix
//      R = stacked_vv * S^T, S = -cd M_{c1}(c,:) diag(w);
// cd2: for each v1, the rows (v1, c1) of all c1 as one (n_c x N) matrix
//      R = stacked_cv * S^T, S = -cd2 M_{v1}(c,:) diag(w).
// R comes out row-major in the order of the BSE vectors and multiplies all
// input vectors at once: the same flops as the row kernel, but in large
// products, and the input is read once per block.
template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::MatrixXd BSE_OPERATOR<cqp, cx, cd, cd2>::ApplyDirectBlocks(
    const Eigen::MatrixXd& input) const {
  using RowMajorMatrix =
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  const Index vmin = opt_.vmin - opt_.rpamin;
  const Index cmin = bse_cmin_ - opt_.rpamin;
  const Index nv = bse_vtotal_;
  const Index nc = bse_ctotal_;
  const Index aux = Mmn_.auxsize();
  const double prefactor = -double((cd != 0) ? cd : cd2);

  // cd: the vv blocks of the occupied slices, row v1*n_v + v2;
  // cd2: the cv blocks of the virtual slices, row c1*n_v + v2
  const Index nblocks = (cd != 0) ? nc : nv;
  const Index blockrows = (cd != 0) ? nv : nc;
  if (stacked_.size() == 0) {
    const Index nslice = (cd != 0) ? nv : nc;
    const Index first_slice = (cd != 0) ? vmin : cmin;
    stacked_.resize(nslice * nv, aux);
#pragma omp parallel for
    for (Index m = 0; m < nslice; m++) {
      stacked_.middleRows(m * nv, nv) =
          Mmn_[m + first_slice].middleRows(vmin, nv);
    }
  }

  Eigen::MatrixXd result(bse_size_, input.cols());
  Eigen::MatrixXd S(stacked_.rows() / blockrows, aux);
  RowMajorMatrix T(stacked_.rows(), S.rows());
  Eigen::MatrixXd out(blockrows, input.cols());
  for (Index b = 0; b < nblocks; b++) {
    // cd: b = c1, S = M_{c1}(c,:); cd2: b = v1, S = M_{v1}(c,:)
    const Index slice = (cd != 0) ? b + cmin : b + vmin;
    S.noalias() = prefactor * Mmn_[slice].middleRows(cmin, nc) *
                  epsilon_0_inv_.asDiagonal();
    T.noalias() = stacked_ * S.transpose();
    const Eigen::Map<const RowMajorMatrix> R(T.data(), blockrows, bse_size_);
    out.noalias() = R * input;
    if (cd != 0) {
      for (Index v1 = 0; v1 < nv; v1++) {
        result.row(v1 * nc + b) = out.row(v1);
      }
    } else {
      result.middleRows(b * nc, nc) = out;
    }
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
