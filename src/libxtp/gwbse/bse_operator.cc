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
#include "votca/xtp/openmp_cuda.h"
#include "votca/xtp/vc2index.h"

namespace votca {
namespace xtp {

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR<cqp, cx, cd, cd2>::configure(BSEOperator_Options opt) {
  _opt = opt;
  Index bse_vmax = _opt.homo;
  _bse_cmin = _opt.homo + 1;
  _bse_vtotal = bse_vmax - _opt.vmin + 1;
  _bse_ctotal = _opt.cmax - _bse_cmin + 1;
  _bse_size = _bse_vtotal * _bse_ctotal;
  this->set_size(_bse_size);
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::MatrixXd BSE_OPERATOR<cqp, cx, cd, cd2>::matmul(
    const Eigen::MatrixXd& input) const {

  static_assert(!(cd2 != 0 && cd != 0),
                "Hamiltonian cannot contain Hd and Hd2 at the same time");

  Index auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);

  Index vmin = _opt.vmin - _opt.rpamin;
  Index cmin = _bse_cmin - _opt.rpamin;

  OpenMP_CUDA transform;
  if (cd != 0) {
    transform.createTemporaries(_epsilon_0_inv, input, _bse_ctotal, _bse_vtotal,
                                auxsize);
  } else {
    transform.createTemporaries(_epsilon_0_inv, input, _bse_vtotal, _bse_ctotal,
                                auxsize);
  }

#pragma omp parallel for schedule(dynamic)
  for (Index c1 = 0; c1 < _bse_ctotal; c1++) {

    Eigen::MatrixXd Temp;
    if (cd != 0) {
      Temp = -cd * (_Mmn[c1 + cmin].block(cmin, 0, _bse_ctotal, auxsize));
      transform.PrepareMatrix1(Temp);
    }
    if (cd2 != 0) {
      Temp = -cd2 * (_Mmn[c1 + cmin].block(vmin, 0, _bse_vtotal, auxsize));
      transform.PrepareMatrix1(Temp);
    }

    for (Index v1 = 0; v1 < _bse_vtotal; v1++) {
      transform.SetTempZero();
      Eigen::MatrixXd Temp2;
      if (cd != 0) {
        Temp2 = _Mmn[v1 + vmin].block(vmin, 0, _bse_vtotal, auxsize);
        transform.PrepareMatrix2(Temp2, cd2 != 0);
      }
      if (cd2 != 0) {
        Temp2 = _Mmn[v1 + vmin].block(cmin, 0, _bse_ctotal, auxsize);
        transform.PrepareMatrix2(Temp2, cd2 != 0);
      }
      if (cqp != 0) {
        Eigen::VectorXd vec = Hqp_row(v1, c1);
        transform.Addvec(vec);
      }
      transform.MultiplyRow(vc.I(v1, c1));
    }
  }
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(input.rows(), input.cols());
  if (cx > 0) {
#pragma omp parallel for schedule(dynamic) reduction(+ : result)
    for (Index v1 = 0; v1 < _bse_vtotal; v1++) {
      Index va = v1 + vmin;
      Eigen::MatrixXd Mmn1 = cx * _Mmn[va].block(cmin, 0, _bse_ctotal, auxsize);
      for (Index v2 = v1; v2 < _bse_vtotal; v2++) {
        Index vb = v2 + vmin;
        const Eigen::MatrixXd blockmat =
            Mmn1 * _Mmn[vb].block(cmin, 0, _bse_ctotal, auxsize).transpose();
        result.block(v1 * _bse_ctotal, 0, _bse_ctotal, result.cols()) +=
            blockmat *
            input.block(v2 * _bse_ctotal, 0, _bse_ctotal, input.cols());
        if (v1 != v2) {
          result.block(v2 * _bse_ctotal, 0, _bse_ctotal, result.cols()) +=
              blockmat.transpose() *
              input.block(v1 * _bse_ctotal, 0, _bse_ctotal, input.cols());
        }
      }
    }
  }

  return result + transform.getReductionVar();
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hqp_row(Index v1,
                                                           Index c1) const {
  Eigen::MatrixXd Result = Eigen::MatrixXd::Zero(_bse_ctotal, _bse_vtotal);
  Index cmin = _bse_vtotal;
  // v->c
  Result.col(v1) += cqp * _Hqp.col(c1 + cmin).segment(cmin, _bse_ctotal);
  // c-> v
  Result.row(c1) -= cqp * _Hqp.col(v1).head(_bse_vtotal);
  return Eigen::Map<Eigen::VectorXd>(Result.data(), Result.size());
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::diagonal() const {

  static_assert(!(cd2 != 0 && cd != 0),
                "Hamiltonian cannot contain Hd and Hd2 at the same time");

  vc2index vc = vc2index(0, 0, _bse_ctotal);
  Index vmin = _opt.vmin - _opt.rpamin;
  Index cmin = _bse_cmin - _opt.rpamin;

  Eigen::VectorXd result = Eigen::VectorXd::Zero(_bse_size);

#pragma omp parallel for schedule(dynamic) reduction(+ : result)
  for (Index v = 0; v < _bse_vtotal; v++) {
    for (Index c = 0; c < _bse_ctotal; c++) {

      double entry = 0.0;
      if (cx != 0) {
        entry += cx * _Mmn[v + vmin].row(cmin + c).squaredNorm();
      }

      if (cqp != 0) {
        Index cmin_qp = _bse_vtotal;
        entry += cqp * (_Hqp(c + cmin_qp, c + cmin_qp) - _Hqp(v, v));
      }
      if (cd != 0) {
        entry -=
            cd * (_Mmn[c + cmin].row(c + cmin) * _epsilon_0_inv.asDiagonal() *
                  _Mmn[v + vmin].row(v + vmin).transpose())
                     .value();
      }
      if (cd2 != 0) {
        entry -=
            cd2 * (_Mmn[c + cmin].row(v + vmin) * _epsilon_0_inv.asDiagonal() *
                   _Mmn[v + vmin].row(c + cmin).transpose())
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

}  // namespace xtp

}  // namespace votca
