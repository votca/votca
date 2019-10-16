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

  if (cx != 0) {
    _Hx_cache = std::vector<cache_block>(OPENMP::getMaxThreads());
  }
}

template <int cqp, int cx, int cd, int cd2>
void BSE_OPERATOR<cqp, cx, cd, cd2>::cache_block::FillCache(
    const Eigen::MatrixXd& matrix, int index) {
  _index = index + 1;
  _size = matrix.cols() - 1;
  _values.resize(_size);
  for (int i = 0; i < _size; i++) {
    _values[i] = matrix.col(i + 1).transpose();
  }
}

template <int cqp, int cx, int cd, int cd2>
bool BSE_OPERATOR<cqp, cx, cd, cd2>::cache_block::hasValue(int index) const {
  if (index >= _index && index < (_index + _size)) {
    return true;
  } else {
    return false;
  }
}

template <int cqp, int cx, int cd, int cd2>
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::row(int index) const {
  Eigen::RowVectorXd row = Eigen::RowVectorXd::Zero(_bse_size);
  if (cx != 0) {
    row += cx * Hx_row(index);
  }
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
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hx_row(int index) const {
  int thread_id = OPENMP::getThreadId();

  if (_Hx_cache[thread_id].hasValue(index)) {
    return _Hx_cache[thread_id].getValue(index);
  }

  int auxsize = _Mmn.auxsize();
  vc2index vc = vc2index(0, 0, _bse_ctotal);

  const int vmin = _opt.vmin - _opt.rpamin;
  const int cmin = _bse_cmin - _opt.rpamin;
  int v1 = vc.v(index);
  int c1 = vc.c(index);
  int cache_size = 50;

  if ((_bse_ctotal - c1) < cache_size) {
    cache_size = _bse_ctotal - c1;
  }

  Eigen::MatrixXd H_cache = Eigen::MatrixXd::Zero(_bse_size, cache_size);
  const Eigen::MatrixXd Mmn1T =
      _Mmn[v1 + vmin].block(c1 + cmin, 0, cache_size, auxsize).transpose();
  for (int v2 = 0; v2 < _bse_vtotal; v2++) {
    const Eigen::MatrixXd& Mmn2 = _Mmn[v2 + vmin];
    int i2 = vc.I(v2, 0);
    H_cache.block(i2, 0, _bse_ctotal, cache_size) =
        Mmn2.block(cmin, 0, _bse_ctotal, auxsize) * Mmn1T;
  }
  _Hx_cache[thread_id].FillCache(H_cache, index);
  return H_cache.col(0).transpose();
}

template <int cqp, int cx, int cd, int cd2>
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hd_row(int index) const {
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
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hqp_row(int index) const {
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
Eigen::RowVectorXd BSE_OPERATOR<cqp, cx, cd, cd2>::Hd2_row(int index) const {

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
