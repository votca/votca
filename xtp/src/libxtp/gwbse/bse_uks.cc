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

#include <chrono>
#include <iostream>

#include <votca/tools/linalg.h>

#include "votca/xtp/bse_uks.h"
#include "votca/xtp/bse_operator_uks.h"
#include "votca/xtp/bseoperator_btda.h"
#include "votca/xtp/davidsonsolver.h"

using std::flush;

namespace votca {
namespace xtp {

void BSE_UKS::configure_with_precomputed_screening(
    const options& opt,
    Index homo_alpha, Index homo_beta,
    const Eigen::VectorXd& RPAInputEnergiesAlpha,
    const Eigen::VectorXd& RPAInputEnergiesBeta,
    const Eigen::MatrixXd& Hqp_alpha_in,
    const Eigen::MatrixXd& Hqp_beta_in,
    const Eigen::VectorXd& epsilon_0_inv) {

  opt_ = opt;
  homo_alpha_ = homo_alpha;
  homo_beta_ = homo_beta;

  alpha_vtotal_ = homo_alpha_ - opt_.vmin + 1;
  alpha_ctotal_ = opt_.cmax - (homo_alpha_ + 1) + 1;
  alpha_size_ = alpha_vtotal_ * alpha_ctotal_;

  beta_vtotal_ = homo_beta_ - opt_.vmin + 1;
  beta_ctotal_ = opt_.cmax - (homo_beta_ + 1) + 1;
  beta_size_ = beta_vtotal_ * beta_ctotal_;

  if (opt_.use_Hqp_offdiag) {
    Hqp_alpha_ = AdjustHqpSize(Hqp_alpha_in, RPAInputEnergiesAlpha, homo_alpha_);
    Hqp_beta_ = AdjustHqpSize(Hqp_beta_in, RPAInputEnergiesBeta, homo_beta_);
  } else {
    Hqp_alpha_ =
        AdjustHqpSize(Hqp_alpha_in, RPAInputEnergiesAlpha, homo_alpha_)
            .diagonal()
            .asDiagonal();
    Hqp_beta_ =
        AdjustHqpSize(Hqp_beta_in, RPAInputEnergiesBeta, homo_beta_)
            .diagonal()
            .asDiagonal();
  }

  epsilon_0_inv_ = epsilon_0_inv;
}

Eigen::MatrixXd BSE_UKS::AdjustHqpSize(const Eigen::MatrixXd& Hqp,
                                       const Eigen::VectorXd& RPAInputEnergies,
                                       Index homo) const {
  Index bse_vtotal = homo - opt_.vmin + 1;
  Index bse_ctotal = opt_.cmax - (homo + 1) + 1;
  Index hqp_size = bse_vtotal + bse_ctotal;
  Index gwsize = opt_.qpmax - opt_.qpmin + 1;
  Index RPAoffset = opt_.vmin - opt_.rpamin;
  Eigen::MatrixXd Hqp_BSE = Eigen::MatrixXd::Zero(hqp_size, hqp_size);

  if (opt_.vmin >= opt_.qpmin) {
    Index start = opt_.vmin - opt_.qpmin;
    if (opt_.cmax <= opt_.qpmax) {
      Hqp_BSE = Hqp.block(start, start, hqp_size, hqp_size);
    } else {
      Index virtoffset = gwsize - start;
      Hqp_BSE.topLeftCorner(virtoffset, virtoffset) =
          Hqp.block(start, start, virtoffset, virtoffset);

      Index virt_extra = opt_.cmax - opt_.qpmax;
      Hqp_BSE.diagonal().tail(virt_extra) =
          RPAInputEnergies.segment(RPAoffset + virtoffset, virt_extra);
    }
  }

  if (opt_.vmin < opt_.qpmin) {
    Index occ_extra = opt_.qpmin - opt_.vmin;
    Hqp_BSE.diagonal().head(occ_extra) =
        RPAInputEnergies.segment(RPAoffset, occ_extra);

    Hqp_BSE.block(occ_extra, occ_extra, gwsize, gwsize) = Hqp;

    if (opt_.cmax > opt_.qpmax) {
      Index virtoffset = occ_extra + gwsize;
      Index virt_extra = opt_.cmax - opt_.qpmax;
      Hqp_BSE.diagonal().tail(virt_extra) =
          RPAInputEnergies.segment(RPAoffset + virtoffset, virt_extra);
    }
  }

  return Hqp_BSE;
}

tools::EigenSystem BSE_UKS::Solve_excitons_uks_TDA() const {
  ExcitonUKSOperator_TDA H(epsilon_0_inv_, Mmn_, Hqp_alpha_, Hqp_beta_);
  BSEOperatorUKS_Options opt;
  opt.homo_alpha = homo_alpha_;
  opt.homo_beta = homo_beta_;
  opt.rpamin = opt_.rpamin;
  opt.qpmin = opt_.qpmin;
  opt.vmin = opt_.vmin;
  opt.cmax = opt_.cmax;
  H.configure(opt);
  Eigen::MatrixXd dense = Eigen::MatrixXd::Zero(H.rows(), H.cols());
  for (Index i = 0; i < H.cols(); ++i) {
    Eigen::MatrixXd e = Eigen::MatrixXd::Zero(H.rows(), 1);
    e(i, 0) = 1.0;
    dense.col(i) = H.matmul(e);
  }

  double asym = (dense - dense.transpose()).norm();
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " TDA asymmetry ||A-A^T|| = " << asym << flush;
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Setup combined UKS TDA exciton hamiltonian " << flush;
  return solve_hermitian(H);
}

tools::EigenSystem BSE_UKS::Solve_excitons_uks_BTDA() const {
  ExcitonUKSOperator_TDA A(epsilon_0_inv_, Mmn_, Hqp_alpha_, Hqp_beta_);
  ExcitonUKSOperator_BTDA_B B(epsilon_0_inv_, Mmn_, Hqp_alpha_, Hqp_beta_);

  BSEOperatorUKS_Options opt;
  opt.homo_alpha = homo_alpha_;
  opt.homo_beta = homo_beta_;
  opt.rpamin = opt_.rpamin;
  opt.qpmin = opt_.qpmin;
  opt.vmin = opt_.vmin;
  opt.cmax = opt_.cmax;

  A.configure(opt);
  B.configure(opt);

  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Setup combined UKS full BSE exciton hamiltonian "
      << flush;
  return Solve_nonhermitian_Davidson(A, B);
}

void BSE_UKS::Solve_excitons_uks(Orbitals& orb) const {
  orb.setTDAApprox(opt_.useTDA);
  if (opt_.useTDA) {
    orb.BSEUKS() = Solve_excitons_uks_TDA();
  } else {
    orb.BSEUKS() = Solve_excitons_uks_BTDA();
  }
}

template <typename BSE_OPERATOR>
tools::EigenSystem BSE_UKS::solve_hermitian(BSE_OPERATOR& h) const {
  std::chrono::time_point<std::chrono::system_clock> start =
      std::chrono::system_clock::now();

  tools::EigenSystem result;

  DavidsonSolver DS(log_);
  DS.set_correction(opt_.davidson_correction);
  DS.set_tolerance(opt_.davidson_tolerance);
  DS.set_size_update(opt_.davidson_update);
  DS.set_iter_max(opt_.davidson_maxiter);
  DS.set_max_search_space(10 * opt_.nmax);
  DS.solve(h, opt_.nmax);

  result.eigenvalues() = DS.eigenvalues();
  result.eigenvectors() = DS.eigenvectors();

  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;

  XTP_LOG(Log::info, log_) << TimeStamp() << " Diagonalization done in "
                           << elapsed_time.count() << " secs" << flush;

  return result;
}

template <typename BSE_OPERATOR_A, typename BSE_OPERATOR_B>
tools::EigenSystem BSE_UKS::Solve_nonhermitian_Davidson(
    BSE_OPERATOR_A& Aop, BSE_OPERATOR_B& Bop) const {
  std::chrono::time_point<std::chrono::system_clock> start =
      std::chrono::system_clock::now();

  HamiltonianOperator<BSE_OPERATOR_A, BSE_OPERATOR_B> Hop(Aop, Bop);

  DavidsonSolver DS(log_);
  DS.set_correction(opt_.davidson_correction);
  DS.set_tolerance(opt_.davidson_tolerance);
  DS.set_size_update(opt_.davidson_update);
  DS.set_iter_max(opt_.davidson_maxiter);
  DS.set_max_search_space(10 * opt_.nmax);
  DS.set_matrix_type("HAM");
  DS.solve(Hop, opt_.nmax);

  tools::EigenSystem result;
  result.eigenvalues() = DS.eigenvalues();

  Eigen::MatrixXd tmpX = DS.eigenvectors().topRows(Aop.rows());
  Eigen::MatrixXd tmpY = DS.eigenvectors().bottomRows(Bop.rows());

  Eigen::VectorXd normX = tmpX.colwise().squaredNorm();
  Eigen::VectorXd normY = tmpY.colwise().squaredNorm();
  Eigen::ArrayXd sqinvnorm = (normX - normY).array().inverse().cwiseSqrt();

  result.eigenvectors() = tmpX * sqinvnorm.matrix().asDiagonal();
  result.eigenvectors2() = tmpY * sqinvnorm.matrix().asDiagonal();

  std::chrono::time_point<std::chrono::system_clock> end =
      std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_time = end - start;

  XTP_LOG(Log::info, log_) << TimeStamp() << " Diagonalization done in "
                           << elapsed_time.count() << " secs" << flush;

  return result;
}

void BSE_UKS::PrintWeightsUKS(const Eigen::VectorXd& coeffs) const {
  const Eigen::VectorXd weights = coeffs.array().square();

  double alpha_weight = weights.head(alpha_size_).sum();
  double beta_weight = weights.tail(beta_size_).sum();

  XTP_LOG(Log::error, log_)
      << boost::format("           alpha-sector: %1$5.1f%%   beta-sector: %2$5.1f%%")
             % (100.0 * alpha_weight) % (100.0 * beta_weight)
      << flush;

  for (Index i = 0; i < alpha_size_; ++i) {
    if (weights(i) > opt_.min_print_weight) {
      Index v = i / alpha_ctotal_;
      Index c = i % alpha_ctotal_;
      XTP_LOG(Log::error, log_)
          << boost::format("           [alpha] HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                 % (homo_alpha_ - (opt_.vmin + v))
                 % ((homo_alpha_ + 1 + c) - (homo_alpha_ + 1))
                 % (100.0 * weights(i))
          << flush;
    }
  }

  for (Index i = 0; i < beta_size_; ++i) {
    if (weights(alpha_size_ + i) > opt_.min_print_weight) {
      Index v = i / beta_ctotal_;
      Index c = i % beta_ctotal_;
      XTP_LOG(Log::error, log_)
          << boost::format("           [beta ] HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%")
                 % (homo_beta_ - (opt_.vmin + v))
                 % ((homo_beta_ + 1 + c) - (homo_beta_ + 1))
                 % (100.0 * weights(alpha_size_ + i))
          << flush;
    }
  }
}

void BSE_UKS::Analyze_excitons_uks(
    std::vector<QMFragment<BSE_Population>> fragments,
    const Orbitals& orb) const {
  (void)fragments;

  XTP_LOG(Log::error, log_)
      << "  ====== combined UKS exciton energies (eV) ====== " << flush;

  const tools::EigenSystem& es = orb.BSEUKS();
  for (Index i = 0; i < std::min<Index>(opt_.nmax, es.eigenvalues().size());
       ++i) {
    XTP_LOG(Log::error, log_)
        << boost::format("  XU%-4d %+1.6f") % (i + 1)
               % (es.eigenvalues()(i) * tools::conv::hrt2ev)
        << flush;
    PrintWeightsUKS(es.eigenvectors().col(i));
  }
}

template tools::EigenSystem BSE_UKS::solve_hermitian<ExcitonUKSOperator_TDA>(
    ExcitonUKSOperator_TDA&) const;

template tools::EigenSystem
BSE_UKS::Solve_nonhermitian_Davidson<ExcitonUKSOperator_TDA,
                                     ExcitonUKSOperator_BTDA_B>(
    ExcitonUKSOperator_TDA&, ExcitonUKSOperator_BTDA_B&) const;

}  // namespace xtp
}  // namespace votca