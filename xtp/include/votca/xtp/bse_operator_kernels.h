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
#ifndef VOTCA_XTP_BSE_OPERATOR_KERNELS_H
#define VOTCA_XTP_BSE_OPERATOR_KERNELS_H

// Local VOTCA includes
#include "eigen.h"
#include "threecenter.h"

namespace votca {
namespace xtp {

/// A window of particle-hole pairs (v,c) in a TCMatrix_gwbse: vmin and cmin
/// are row/slice indices of the tensor (counted from its first level), the
/// compound index of a pair is v * ctotal + c, as in vc2index.
struct BSEWindow {
  Index vmin;
  Index cmin;
  Index vtotal;
  Index ctotal;
  Index size() const { return vtotal * ctotal; }
};

/// y = prefactor * K_x x with the bare exchange
///   K_x(v1 c1, v2 c2) = sum_P M_out[v1](c1,P) M_in[v2](c2,P).
/// K_x has rank auxsize, so it is applied as M_out (M_in^T x): about
/// 4 N auxsize k operations instead of building the N x N blocks.
inline Eigen::MatrixXd ApplyBSEExchange(
    const TCMatrix_gwbse& Mout, const BSEWindow& out, const TCMatrix_gwbse& Min,
    const BSEWindow& in, const Eigen::MatrixXd& x, double prefactor) {
  const Index auxsize = Min.auxsize();
  const Index k = x.cols();
  // B = sum_v2 M_in[v2]_c^T x_v2, accumulated per thread
  Eigen::MatrixXd B = Eigen::MatrixXd::Zero(auxsize, k);
#pragma omp parallel
  {
    Eigen::MatrixXd B_thread = Eigen::MatrixXd::Zero(auxsize, k);
#pragma omp for schedule(dynamic)
    for (Index v2 = 0; v2 < in.vtotal; ++v2) {
      B_thread.noalias() +=
          Min[v2 + in.vmin].middleRows(in.cmin, in.ctotal).transpose() *
          x.middleRows(v2 * in.ctotal, in.ctotal);
    }
#pragma omp critical
    B += B_thread;
  }
  Eigen::MatrixXd y(out.size(), k);
#pragma omp parallel for schedule(dynamic)
  for (Index v1 = 0; v1 < out.vtotal; ++v1) {
    y.middleRows(v1 * out.ctotal, out.ctotal).noalias() =
        prefactor * Mout[v1 + out.vmin].middleRows(out.cmin, out.ctotal) * B;
  }
  return y;
}

/// y = prefactor * H_qp x with
///   H_qp(v1 c1, v2 c2) = delta_v1v2 Hqp(c1,c2) - delta_c1c2 Hqp(v1,v2),
/// Hqp indexed from the first valence level of the window (the conduction
/// levels start at vtotal). For each column, x is the ctotal x vtotal matrix
/// X(c,v), and H_qp x = Hqp_cc X - X Hqp_vv.
inline Eigen::MatrixXd ApplyBSEQuasiparticle(const Eigen::MatrixXd& Hqp,
                                             Index vtotal, Index ctotal,
                                             const Eigen::MatrixXd& x,
                                             double prefactor) {
  const Eigen::MatrixXd Hvv = Hqp.topLeftCorner(vtotal, vtotal);
  const Eigen::MatrixXd Hcc = Hqp.block(vtotal, vtotal, ctotal, ctotal);
  Eigen::MatrixXd y(x.rows(), x.cols());
  // all columns at once for Hcc: x as ctotal x (vtotal * k)
  Eigen::Map<const Eigen::MatrixXd> X_all(x.data(), ctotal, vtotal * x.cols());
  Eigen::Map<Eigen::MatrixXd> Y_all(y.data(), ctotal, vtotal * x.cols());
  Y_all.noalias() = prefactor * Hcc * X_all;
  for (Index col = 0; col < x.cols(); ++col) {
    Eigen::Map<const Eigen::MatrixXd> X(x.col(col).data(), ctotal, vtotal);
    Eigen::Map<Eigen::MatrixXd> Y(y.col(col).data(), ctotal, vtotal);
    Y.noalias() -= prefactor * X * Hvv;
  }
  return y;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_BSE_OPERATOR_KERNELS_H
