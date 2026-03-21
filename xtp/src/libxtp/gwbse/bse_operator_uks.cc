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

 #include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "votca/xtp/bse_operator_uks.h"

namespace votca {
namespace xtp {

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::setup_block(SpinBlockInfo& blk,
                                                     Index homo, Index offset) {
  blk.homo = homo;
  blk.vmin_rpa = opt_.vmin - opt_.rpamin;
  blk.cmin_rpa = homo + 1 - opt_.rpamin;
  blk.vtotal = homo - opt_.vmin + 1;
  blk.ctotal = opt_.cmax - (homo + 1) + 1;
  blk.size = blk.vtotal * blk.ctotal;
  blk.offset = offset;
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::configure(BSEOperatorUKS_Options opt) {
  opt_ = opt;
  setup_block(alpha_, opt_.homo_alpha, 0);
  setup_block(beta_, opt_.homo_beta, alpha_.size);
  size_total_ = alpha_.size + beta_.size;
  this->set_size(size_total_);
}


template <Index cqp, Index cx, Index cd, Index cd2>
std::string BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::spin_block_info_string(
    const SpinBlockInfo& blk) const {
  std::ostringstream oss;
  oss << "{homo=" << blk.homo
      << ", vmin_rpa=" << blk.vmin_rpa
      << ", cmin_rpa=" << blk.cmin_rpa
      << ", vtotal=" << blk.vtotal
      << ", ctotal=" << blk.ctotal
      << ", size=" << blk.size
      << ", offset=" << blk.offset
      << "}";
  return oss.str();
}

template <Index cqp, Index cx, Index cd, Index cd2>
std::string BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::matrix_label(
    const TCMatrix_gwbse& M) const {
  if (&M == &Mmn_.alpha) {
    return "Mmn_.alpha";
  }
  if (&M == &Mmn_.beta) {
    return "Mmn_.beta";
  }
  return "unknown";
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::log_add_direct2_call(
    const std::string& block_label, const SpinBlockInfo& out_blk,
    const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
    const TCMatrix_gwbse& Min, double prefactor) const {
  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  oss << std::setprecision(16);
  oss << "[UKS Hd2 call] block=" << block_label
      << " prefactor=" << prefactor
      << " Mout=" << matrix_label(Mout)
      << " Min=" << matrix_label(Min)
      << " out_blk=" << spin_block_info_string(out_blk)
      << " in_blk=" << spin_block_info_string(in_blk);
  std::cout << oss.str() << std::endl;
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::VectorXd BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::Hqp_row(
    const Eigen::MatrixXd& Hqp, const SpinBlockInfo& blk, Index v1,
    Index c1) const {
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(blk.ctotal, blk.vtotal);
  Index cmin_qp = blk.vtotal;
  result.col(v1) += Hqp.col(c1 + cmin_qp).segment(cmin_qp, blk.ctotal);
  result.row(c1) -= Hqp.col(v1).head(blk.vtotal).transpose();
  return Eigen::Map<Eigen::VectorXd>(result.data(), result.size());
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::add_qp_block(
    Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const SpinBlockInfo& blk,
    const Eigen::MatrixXd& Hqp) const {
  if (cqp == 0) {
    return;
  }

  for (Index c1 = 0; c1 < blk.ctotal; ++c1) {
    for (Index v1 = 0; v1 < blk.vtotal; ++v1) {
      const Index out_idx = v1 * blk.ctotal + c1;
      Eigen::VectorXd row = Hqp_row(Hqp, blk, v1, c1);
      y.row(out_idx) += row.transpose() * x;
    }
  }
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::add_exchange_block(
    Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const SpinBlockInfo& out_blk,
    const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
    const TCMatrix_gwbse& Min, double prefactor) const {
  if (cx == 0 || prefactor == 0.0) {
    return;
  }

  for (Index v1 = 0; v1 < out_blk.vtotal; ++v1) {
    const Eigen::MatrixXd left =
        prefactor * Mout[v1 + out_blk.vmin_rpa].middleRows(out_blk.cmin_rpa,
                                                           out_blk.ctotal);

    for (Index v2 = 0; v2 < in_blk.vtotal; ++v2) {
      const Eigen::MatrixXd right =
          Min[v2 + in_blk.vmin_rpa].middleRows(in_blk.cmin_rpa, in_blk.ctotal);

      const Eigen::MatrixXd block = left * right.transpose();
      const Index out_row0 = v1 * out_blk.ctotal;
      const Index in_row0 = v2 * in_blk.ctotal;
      y.middleRows(out_row0, out_blk.ctotal) +=
          block * x.middleRows(in_row0, in_blk.ctotal);
    }
  }
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::add_direct_block(
    Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const SpinBlockInfo& out_blk,
    const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
    const TCMatrix_gwbse& Min, double prefactor) const {
  if (cd == 0 || prefactor == 0.0) {
    return;
  }

  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> eps =
      epsilon_0_inv_.asDiagonal();

  for (Index v1 = 0; v1 < out_blk.vtotal; ++v1) {
    for (Index c1 = 0; c1 < out_blk.ctotal; ++c1) {

      // rows correspond to c2, columns correspond to v2
      const Eigen::MatrixXd left =
          prefactor * Mout[c1 + out_blk.cmin_rpa].middleRows(in_blk.cmin_rpa,
                                                             in_blk.ctotal);

      const Eigen::MatrixXd right =
          Min[v1 + out_blk.vmin_rpa].middleRows(in_blk.vmin_rpa, in_blk.vtotal);

      // (ctotal_in x naux) * (naux x vtotal_in) = (ctotal_in x vtotal_in)
      const Eigen::MatrixXd block = left * eps * right.transpose();

      // Flattening column-major yields index order v2*ctotal + c2
      const Eigen::VectorXd row =
          Eigen::Map<const Eigen::VectorXd>(block.data(), block.size());

      const Index out_idx = v1 * out_blk.ctotal + c1;
      y.row(out_idx) += row.transpose() * x;
    }
  }
}

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::add_direct2_block(
    Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const SpinBlockInfo& out_blk,
    const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
    const TCMatrix_gwbse& Min, double prefactor) const {
  if (cd2 == 0 || prefactor == 0.0) {
    return;
  }

  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> eps =
      epsilon_0_inv_.asDiagonal();

  for (Index c1 = 0; c1 < out_blk.ctotal; ++c1) {
    const Eigen::MatrixXd left =
        prefactor *
        Mout[c1 + out_blk.cmin_rpa].middleRows(in_blk.vmin_rpa, in_blk.vtotal);
    // left: (vtotal_in x naux)

    for (Index v1 = 0; v1 < out_blk.vtotal; ++v1) {
      const Eigen::MatrixXd right =
          Min[v1 + out_blk.vmin_rpa].middleRows(in_blk.cmin_rpa, in_blk.ctotal);
      // right: (ctotal_in x naux)

      // swapped contraction trial
      const Eigen::MatrixXd block = left * eps * right.transpose();
      // block: (vtotal_in x ctotal_in)

      Eigen::VectorXd row(in_blk.vtotal * in_blk.ctotal);
      for (Index v2 = 0; v2 < in_blk.vtotal; ++v2) {
        for (Index c2 = 0; c2 < in_blk.ctotal; ++c2) {
          const Index idx = v2 * in_blk.ctotal + c2;
          row(idx) = block(v2, c2);
        }
      }

      if (v1 < 2 && c1 < 2) {
        std::ostringstream oss;
        oss.setf(std::ios::scientific);
        oss.precision(16);
        oss << "Hd2 row v1=" << v1 << " c1=" << c1 << ":";
        for (Index k = 0; k < std::min<Index>(8, row.size()); ++k) {
          oss << " " << row(k);
        }
        std::cout << " " << oss.str() << std::endl;
      }

      const Index out_idx = v1 * out_blk.ctotal + c1;
      y.row(out_idx) += row.transpose() * x;
    }
  }
}

/*template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::add_direct2_block(
    Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const SpinBlockInfo& out_blk,
    const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
    const TCMatrix_gwbse& Min, double prefactor) const {
  if (cd2 == 0 || prefactor == 0.0) {
    return;
  }

  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> eps =
      epsilon_0_inv_.asDiagonal();

  for (Index c1 = 0; c1 < out_blk.ctotal; ++c1) {
    const Eigen::MatrixXd left =
        prefactor *
        Mout[c1 + out_blk.cmin_rpa].middleRows(in_blk.vmin_rpa, in_blk.vtotal);

    for (Index v1 = 0; v1 < out_blk.vtotal; ++v1) {
      const Eigen::MatrixXd right =
          Min[v1 + out_blk.vmin_rpa].middleRows(in_blk.cmin_rpa, in_blk.ctotal);

      // We need matrix dimensions (ctotal_in x vtotal_in) so that the
      // flattened memory order matches the existing vc = v*ctotal + c layout.
      const Eigen::MatrixXd block = right * eps * left.transpose();
      const Eigen::VectorXd row =
          Eigen::Map<const Eigen::VectorXd>(block.data(), block.size());

      const Index out_idx = v1 * out_blk.ctotal + c1;
      y.row(out_idx) += row.transpose() * x;
    }
  }
}*/

template <Index cqp, Index cx, Index cd, Index cd2>
void BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::add_direct_cross_tda_block(
    Eigen::MatrixXd& y, const Eigen::MatrixXd& x, const SpinBlockInfo& out_blk,
    const SpinBlockInfo& in_blk, const TCMatrix_gwbse& Mout,
    const TCMatrix_gwbse& Min, double prefactor) const {
  if (cd == 0 || prefactor == 0.0) {
    return;
  }

  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> eps =
      epsilon_0_inv_.asDiagonal();

  for (Index v1 = 0; v1 < out_blk.vtotal; ++v1) {
    for (Index c1 = 0; c1 < out_blk.ctotal; ++c1) {

      // Transition density for the output excitation (v1 -> c1)
      const Eigen::RowVectorXd tout =
          prefactor * Mout[c1 + out_blk.cmin_rpa].row(v1 + out_blk.vmin_rpa);

      const Index out_idx = v1 * out_blk.ctotal + c1;

      // Build the row blockwise in the input excitation space:
      // input ordering is (v2 * ctotal + c2), i.e. c2 runs fastest.
      for (Index v2 = 0; v2 < in_blk.vtotal; ++v2) {
        const Eigen::MatrixXd Tin = Min[v2 + in_blk.vmin_rpa].middleRows(
            in_blk.cmin_rpa, in_blk.ctotal);

        // Tin rows correspond to c2, and column-major flattening over
        // successive v2 blocks is therefore consistent with vc = v*ctotal + c.
        const Eigen::VectorXd row_block = Tin * eps * tout.transpose();

        const Index in_row0 = v2 * in_blk.ctotal;
        y.row(out_idx) +=
            row_block.transpose() * x.middleRows(in_row0, in_blk.ctotal);
      }
    }
  }
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::MatrixXd BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::matmul(
    const Eigen::MatrixXd& input) const {

  static_assert(!(cd != 0 && cd2 != 0),
                "Hamiltonian cannot contain Hd and Hd2 at the same time");

  Eigen::MatrixXd y = Eigen::MatrixXd::Zero(size_total_, input.cols());

  const Eigen::MatrixXd x_alpha = input.topRows(alpha_.size);
  const Eigen::MatrixXd x_beta = input.bottomRows(beta_.size);

  Eigen::MatrixXd y_alpha = Eigen::MatrixXd::Zero(alpha_.size, input.cols());
  Eigen::MatrixXd y_beta = Eigen::MatrixXd::Zero(beta_.size, input.cols());

  // Same-spin alpha-alpha
  add_qp_block(y_alpha, x_alpha, alpha_, Hqp_alpha_);
  add_exchange_block(y_alpha, x_alpha, alpha_, alpha_, Mmn_.alpha, Mmn_.alpha,
                     static_cast<double>(cx));
  add_direct_block(y_alpha, x_alpha, alpha_, alpha_, Mmn_.alpha, Mmn_.alpha,
                   -static_cast<double>(cd));
  log_add_direct2_call("aa", alpha_, alpha_, Mmn_.alpha, Mmn_.alpha,
                       -static_cast<double>(cd2));
  add_direct2_block(y_alpha, x_alpha, alpha_, alpha_, Mmn_.alpha, Mmn_.alpha,
                    -static_cast<double>(cd2));

  // Same-spin beta-beta
  add_qp_block(y_beta, x_beta, beta_, Hqp_beta_);
  add_exchange_block(y_beta, x_beta, beta_, beta_, Mmn_.beta, Mmn_.beta,
                     static_cast<double>(cx));
  add_direct_block(y_beta, x_beta, beta_, beta_, Mmn_.beta, Mmn_.beta,
                   -static_cast<double>(cd));
  log_add_direct2_call("bb", beta_, beta_, Mmn_.beta, Mmn_.beta,
                       -static_cast<double>(cd2));
  add_direct2_block(y_beta, x_beta, beta_, beta_, Mmn_.beta, Mmn_.beta,
                    -static_cast<double>(cd2));

  // Cross-spin TDA coupling: use transition-density form.
  add_direct_cross_tda_block(y_alpha, x_beta, alpha_, beta_, Mmn_.alpha,
                             Mmn_.beta, -static_cast<double>(cd));
  add_direct_cross_tda_block(y_beta, x_alpha, beta_, alpha_, Mmn_.beta,
                             Mmn_.alpha, -static_cast<double>(cd));

  // Cross-spin full-BSE B-block coupling is not the same object as the TDA
  // cross block above; keep the existing Hd2-style contraction for now.
   log_add_direct2_call("ab", alpha_, beta_, Mmn_.alpha, Mmn_.beta,
                       -static_cast<double>(cd2));
  add_direct2_block(y_alpha, x_beta, alpha_, beta_, Mmn_.alpha, Mmn_.beta,
                    -static_cast<double>(cd2));

  log_add_direct2_call("ba", beta_, alpha_, Mmn_.beta, Mmn_.alpha,
                       -static_cast<double>(cd2));
  add_direct2_block(y_beta, x_alpha, beta_, alpha_, Mmn_.beta, Mmn_.alpha,
                    -static_cast<double>(cd2));

  y.topRows(alpha_.size) = y_alpha;
  y.bottomRows(beta_.size) = y_beta;
  return y;
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::MatrixXd BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::dense_matrix() const {
  Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(rows(), cols());
  for (Index i = 0; i < cols(); ++i) {
    Eigen::MatrixXd e = Eigen::MatrixXd::Zero(rows(), 1);
    e(i, 0) = 1.0;
    dense.col(i) = matmul(e);
  }
  return dense;
}

template <Index cqp, Index cx, Index cd, Index cd2>
Eigen::VectorXd BSE_OPERATOR_UKS<cqp, cx, cd, cd2>::diagonal() const {
  static_assert(!(cd != 0 && cd2 != 0),
                "Hamiltonian cannot contain Hd and Hd2 at the same time");

  Eigen::VectorXd result = Eigen::VectorXd::Zero(size_total_);

  const Eigen::DiagonalMatrix<double, Eigen::Dynamic> eps =
      epsilon_0_inv_.asDiagonal();

  // alpha block
  for (Index v = 0; v < alpha_.vtotal; ++v) {
    for (Index c = 0; c < alpha_.ctotal; ++c) {
      double entry = 0.0;

      if (cx != 0) {
        entry += cx * Mmn_.alpha[v + alpha_.vmin_rpa]
                          .row(alpha_.cmin_rpa + c)
                          .squaredNorm();
      }
      if (cqp != 0) {
        Index cmin_qp = alpha_.vtotal;
        entry += Hqp_alpha_(c + cmin_qp, c + cmin_qp) - Hqp_alpha_(v, v);
      }
      if (cd != 0) {
        entry -=
            (Mmn_.alpha[c + alpha_.cmin_rpa].row(alpha_.cmin_rpa + c) * eps *
             Mmn_.alpha[v + alpha_.vmin_rpa]
                 .row(alpha_.vmin_rpa + v)
                 .transpose())
                .value();
      }
      if (cd2 != 0) {
        entry -=
            (Mmn_.alpha[c + alpha_.cmin_rpa].row(alpha_.vmin_rpa + v) * eps *
             Mmn_.alpha[v + alpha_.vmin_rpa]
                 .row(alpha_.cmin_rpa + c)
                 .transpose())
                .value();
      }

      result(alpha_.offset + v * alpha_.ctotal + c) = entry;
    }
  }

  // beta block
  for (Index v = 0; v < beta_.vtotal; ++v) {
    for (Index c = 0; c < beta_.ctotal; ++c) {
      double entry = 0.0;

      if (cx != 0) {
        entry +=
            cx *
            Mmn_.beta[v + beta_.vmin_rpa].row(beta_.cmin_rpa + c).squaredNorm();
      }
      if (cqp != 0) {
        Index cmin_qp = beta_.vtotal;
        entry += Hqp_beta_(c + cmin_qp, c + cmin_qp) - Hqp_beta_(v, v);
      }
      if (cd != 0) {
        entry -=
            (Mmn_.beta[c + beta_.cmin_rpa].row(beta_.cmin_rpa + c) * eps *
             Mmn_.beta[v + beta_.vmin_rpa].row(beta_.vmin_rpa + v).transpose())
                .value();
      }
      if (cd2 != 0) {
        entry -=
            (Mmn_.beta[c + beta_.cmin_rpa].row(beta_.vmin_rpa + v) * eps *
             Mmn_.beta[v + beta_.vmin_rpa].row(beta_.cmin_rpa + c).transpose())
                .value();
      }

      result(beta_.offset + v * beta_.ctotal + c) = entry;
    }
  }

  return result;
}

template class BSE_OPERATOR_UKS<1, 1, 1, 0>;
template class BSE_OPERATOR_UKS<0, 1, 0, 1>;
template class BSE_OPERATOR_UKS<0, 0, 1, 0>;
template class BSE_OPERATOR_UKS<0, 0, 0, 1>;

template class BSE_OPERATOR_UKS<1, 0, 0, 0>;
template class BSE_OPERATOR_UKS<0, 1, 0, 0>;

}  // namespace xtp
}  // namespace votca