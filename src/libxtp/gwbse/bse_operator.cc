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

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

template <int cqp, int cx, int cd, int cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::col(int index) const {
  Eigen::VectorXd col = Eigen::VectorXd::Zero(_bse_size);
  if (cx != 0) {
    col += cx * Hx_col(index);
  }
  if (cd != 0) {
    col += cd * Hd_col(index);
  }
  if (cd2 != 0) {
    col += cd2 * Hd2_col(index);
  }
  if (cqp != 0) {
    col += cqp * Hqp_col(index);
  }
  return col;
}

template <int cqp, int cx, int cd, int cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hx_col(int index) const {
  int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  Eigen::VectorXd Hcol = Eigen::VectorXd::Zero(_bse_size);
  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = vc.v(index);
  int c1 = vc.c(index);
  const Eigen::VectorXd Mmn1Tc1 = _Mmn[v1 + vmin].row(c1 + cmin);
  for (int v2 = 0; v2 < _bse_vtotal; v2++) {
    const Eigen::MatrixXd& Mmn2 = _Mmn[v2 + vmin];
    auto Mmnx2 = Mmn2.block(cmin, 0, _bse_ctotal, auxsize) * Mmn1Tc1;
    int i2 = vc.I(v2, 0);
    Hcol.segment(i2, _bse_ctotal).noalias() = Mmnx2;
  }
  return Hcol;
}

template <int cqp, int cx, int cd, int cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hd_col(int index) const {
  int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  Eigen::VectorXd Hcol = Eigen::VectorXd::Zero(_bse_size);
  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = vc.v(index);
  int c1 = vc.c(index);

  const Eigen::MatrixXd Mmn1T =
      (_Mmn[v1 + vmin].block(vmin, 0, _bse_vtotal, auxsize) *
       _epsilon_0_inv.asDiagonal())
          .transpose();
  const Eigen::MatrixXd& Mmn2 = _Mmn[c1 + cmin];
  const Eigen::MatrixXd Mmn2xMmn1T =
      Mmn2.block(cmin, 0, _bse_ctotal, auxsize) * Mmn1T;

  for (int v2 = 0; v2 < _bse_vtotal; v2++) {
    int i2 = vc.I(v2, 0);
    Hcol.segment(i2, _bse_ctotal) = -Mmn2xMmn1T.col(v2);
  }
  return Hcol;
}

template <int cqp, int cx, int cd, int cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hqp_col(int index) const {
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  int v1 = vc.v(index);
  int c1 = vc.c(index);
  int index_vc = vc.I(v1, c1);
  Eigen::VectorXd Hcol = Eigen::VectorXd::Zero(_bse_size);

  // v->c
  for (int c2 = 0; c2 < _bse_ctotal; c2++) {
    int index_vc2 = vc.I(v1, c2);
    Hcol(index_vc2) +=
        _Hqp(c2 + _bse_vtotal - _opt.qpmin, c1 + _bse_vtotal - _opt.qpmin);
  }
  // c-> v
  for (int v2 = 0; v2 < _bse_vtotal; v2++) {
    int index_vc2 = vc.I(v2, c1);
    Hcol(index_vc2) -= _Hqp(v2 - _opt.qpmin, v1 - _opt.qpmin);
  }

  return Hcol;
}

template <int cqp, int cx, int cd, int cd2>
Eigen::VectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hd2_col(int index) const {

  int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = vc.v(index);
  int c1 = vc.c(index);

  Eigen::VectorXd Hcol = Eigen::VectorXd::Zero(_bse_size);

  const Eigen::MatrixXd Mmn2T =
      (_Mmn[c1 + cmin].block(vmin, 0, _bse_vtotal, auxsize) *
       _epsilon_0_inv.asDiagonal())
          .transpose();
  const Eigen::MatrixXd& Mmn1 = _Mmn[v1 + vmin];
  Eigen::MatrixXd Mmn1xMmn2T =
      Mmn1.block(cmin, 0, _bse_ctotal, auxsize) * Mmn2T;

  for (int v2 = 0; v2 < _bse_vtotal; v2++) {
    int i2 = vc.I(v2, 0);
    Hcol.segment(i2, _bse_ctotal) = -Mmn1xMmn2T.col(v2);
  }

  return Hcol;
}

template class BSE_OPERATOR<1, 2, 1, 0>;
template class BSE_OPERATOR<1, 0, 1, 0>;

template class BSE_OPERATOR<1, 4, 1, 1>;
template class BSE_OPERATOR<1, 0, 1, -1>;

template class BSE_OPERATOR<1, 0, 0, 0>;
template class BSE_OPERATOR<0, 1, 0, 0>;
template class BSE_OPERATOR<0, 0, 1, 0>;
template class BSE_OPERATOR<0, 0, 0, 1>;

}  // namespace xtp
}  // namespace votca
