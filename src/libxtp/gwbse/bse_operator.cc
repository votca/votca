/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/xtp/bse_operator.h>
#include <votca/xtp/vc2index.h>

namespace votca {
namespace xtp {

template <int cqp, int cx, int cd, int cd2>
void BSE_OPERATOR<cqp, cx, cd, cd2>::configure(BSEOperator_Options opt) {
  _opt = opt;
  int bse_vmax = _opt.homo;
  _bse_cmin = _opt.homo + 1;
  _bse_vtotal = bse_vmax - _opt.vmin + 1;
  _bse_ctotal = _opt.cmax - _bse_cmin + 1;
  _bse_size = _bse_vtotal * _bse_ctotal;
  this->set_size(_bse_size);
}

template <int cqp, int cx, int cd, int cd2>
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::OperatorRow(
    long index) const {
  Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(_bse_size);
  if (cd != 0) {
    row += cd * Hd_row(index);
  }
  if (cd2 != 0) {
    row += cd2 * Hd2_row(index);
  }
  if (cqp != 0) {
    row += cqp * Hqp_row(index);
  }
  return row;
}

template <int cqp, int cx, int cd, int cd2>
Eigen::MatrixXd BSE_OPERATOR<cqp, cx, cd, cd2>::OperatorBlock(long row,
                                                              long col) const {
  return cx * HxBlock(row, col);
}

template <int cqp, int cx, int cd, int cd2>
Eigen::MatrixXd BSE_OPERATOR<cqp, cx, cd, cd2>::HxBlock(long row,
                                                        long col) const {
  int auxsize = _Mmn.auxsize();
  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = int(row) + vmin;
  int v2 = int(col) + vmin;
  return _Mmn[v1].block(cmin, 0, _bse_ctotal, auxsize) *
         _Mmn[v2].block(cmin, 0, _bse_ctotal, auxsize).transpose();
}

template <int cqp, int cx, int cd, int cd2>
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hd_row(long index) const {
  int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);

  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = vc.v(index);
  int c1 = vc.c(index);

  const Eigen::MatrixXd Mmn1T =
      -(_Mmn[v1 + vmin].block(vmin, 0, _bse_vtotal, auxsize) *
        _epsilon_0_inv.asDiagonal())
           .transpose();
  const Eigen::MatrixXd& Mmn2 = _Mmn[c1 + cmin];
  Eigen::MatrixXd Mmn2xMmn1T =
      Mmn2.block(cmin, 0, _bse_ctotal, auxsize) * Mmn1T;
  return Eigen::Map<Eigen::RowVectorXd>(Mmn2xMmn1T.data(), Mmn2xMmn1T.size());
}

template <int cqp, int cx, int cd, int cd2>
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hqp_row(long index) const {
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  int v1 = vc.v(index);
  int c1 = vc.c(index);
  Eigen::MatrixXd Result = Eigen::MatrixXd::Zero(_bse_ctotal, _bse_vtotal);
  int cmin = _bse_cmin - _opt.qpmin;
  // v->c
  Result.col(v1) = _Hqp.col(c1 + cmin).segment(cmin, _bse_ctotal);
  // c-> v
  int vmin = _opt.vmin - _opt.qpmin;
  Result.row(c1) -= _Hqp.col(v1 + vmin).segment(vmin, _bse_vtotal);
  return Eigen::Map<Eigen::RowVectorXd>(Result.data(), Result.size());
}

template <int cqp, int cx, int cd, int cd2>
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hd2_row(long index) const {

  int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = vc.v(index);
  int c1 = vc.c(index);

  const Eigen::MatrixXd Mmn2T =
      -(_Mmn[c1 + cmin].block(vmin, 0, _bse_vtotal, auxsize) *
        _epsilon_0_inv.asDiagonal())
           .transpose();
  const Eigen::MatrixXd& Mmn1 = _Mmn[v1 + vmin];
  Eigen::MatrixXd Mmn1xMmn2T =
      Mmn1.block(cmin, 0, _bse_ctotal, auxsize) * Mmn2T;
  return Eigen::Map<Eigen::RowVectorXd>(Mmn1xMmn2T.data(), Mmn1xMmn2T.size());
}

template class BSE_OPERATOR<1, 2, 1, 0>;
template class BSE_OPERATOR<1, 0, 1, 0>;

template class BSE_OPERATOR<1, 4, 1, 1>;
template class BSE_OPERATOR<1, 0, 1, 1>;
template class BSE_OPERATOR<1, 0, 1, -1>;

template class BSE_OPERATOR<1, 0, 0, 0>;
template class BSE_OPERATOR<0, 1, 0, 0>;
template class BSE_OPERATOR<0, 0, 1, 0>;
template class BSE_OPERATOR<0, 0, 0, 1>;

template class BSE_OPERATOR<0, 2, 0, 1>;

}  // namespace xtp
}  // namespace votca
