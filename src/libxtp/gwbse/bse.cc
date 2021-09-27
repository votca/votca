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

// Standard includes
#include <chrono>
#include <iostream>

// VOTCA includes
#include <votca/tools/linalg.h>

// Local VOTCA includes
#include "votca/xtp/bse.h"
#include "votca/xtp/bse_operator.h"
#include "votca/xtp/bseoperator_btda.h"
#include "votca/xtp/davidsonsolver.h"
#include "votca/xtp/populationanalysis.h"
#include "votca/xtp/qmfragment.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/vc2index.h"

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void BSE::configure(const options& opt, const Eigen::VectorXd& RPAInputEnergies,
                    const Eigen::MatrixXd& Hqp_in) {
  opt_ = opt;
  bse_vmax_ = opt_.homo;
  bse_cmin_ = opt_.homo + 1;
  bse_vtotal_ = bse_vmax_ - opt_.vmin + 1;
  bse_ctotal_ = opt_.cmax - bse_cmin_ + 1;
  bse_size_ = bse_vtotal_ * bse_ctotal_;
  max_dyn_iter_ = opt_.max_dyn_iter;
  dyn_tolerance_ = opt_.dyn_tolerance;
  if (opt_.use_Hqp_offdiag) {
    Hqp_ = AdjustHqpSize(Hqp_in, RPAInputEnergies);
  } else {
    Hqp_ = AdjustHqpSize(Hqp_in, RPAInputEnergies).diagonal().asDiagonal();
  }
  SetupDirectInteractionOperator(RPAInputEnergies, 0.0);
}

Eigen::MatrixXd BSE::AdjustHqpSize(const Eigen::MatrixXd& Hqp,
                                   const Eigen::VectorXd& RPAInputEnergies) {

  Index hqp_size = bse_vtotal_ + bse_ctotal_;
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

void BSE::SetupDirectInteractionOperator(
    const Eigen::VectorXd& RPAInputEnergies, double energy) {
  RPA rpa = RPA(log_, Mmn_);
  rpa.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  rpa.setRPAInputEnergies(RPAInputEnergies);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      rpa.calculate_epsilon_r(energy));
  Mmn_.MultiplyRightWithAuxMatrix(es.eigenvectors());

  epsilon_0_inv_ = Eigen::VectorXd::Zero(es.eigenvalues().size());
  for (Index i = 0; i < es.eigenvalues().size(); ++i) {
    if (es.eigenvalues()(i) > 1e-8) {
      epsilon_0_inv_(i) = 1 / es.eigenvalues()(i);
    }
  }
}

template <typename BSE_OPERATOR>
void BSE::configureBSEOperator(BSE_OPERATOR& H) const {
  BSEOperator_Options opt;
  opt.cmax = opt_.cmax;
  opt.homo = opt_.homo;
  opt.qpmin = opt_.qpmin;
  opt.rpamin = opt_.rpamin;
  opt.vmin = opt_.vmin;
  H.configure(opt);
}

tools::EigenSystem BSE::Solve_triplets_TDA() const {

  TripletOperator_TDA Ht(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(Ht);
  return solve_hermitian(Ht);
}

void BSE::Solve_singlets(Orbitals& orb) const {
  orb.setTDAApprox(opt_.useTDA);
  if (opt_.useTDA) {
    orb.BSESinglets() = Solve_singlets_TDA();
  } else {
    orb.BSESinglets() = Solve_singlets_BTDA();
  }
  orb.CalcCoupledTransition_Dipoles();
}

void BSE::Solve_triplets(Orbitals& orb) const {
  orb.setTDAApprox(opt_.useTDA);
  if (opt_.useTDA) {
    orb.BSETriplets() = Solve_triplets_TDA();
  } else {
    orb.BSETriplets() = Solve_triplets_BTDA();
  }
}

tools::EigenSystem BSE::Solve_singlets_TDA() const {

  SingletOperator_TDA Hs(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(Hs);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Setup TDA singlet hamiltonian " << flush;
  return solve_hermitian(Hs);
}

SingletOperator_TDA BSE::getSingletOperator_TDA() const {

  SingletOperator_TDA Hs(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(Hs);
  return Hs;
}

TripletOperator_TDA BSE::getTripletOperator_TDA() const {

  TripletOperator_TDA Ht(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(Ht);
  return Ht;
}

template <typename BSE_OPERATOR>
tools::EigenSystem BSE::solve_hermitian(BSE_OPERATOR& h) const {

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

tools::EigenSystem BSE::Solve_singlets_BTDA() const {
  SingletOperator_TDA A(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(A);
  SingletOperator_BTDA_B B(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(B);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Setup Full singlet hamiltonian " << flush;
  return Solve_nonhermitian_Davidson(A, B);
}

tools::EigenSystem BSE::Solve_triplets_BTDA() const {
  TripletOperator_TDA A(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(A);
  Hd2Operator B(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(B);
  XTP_LOG(Log::error, log_)
      << TimeStamp() << " Setup Full triplet hamiltonian " << flush;
  return Solve_nonhermitian_Davidson(A, B);
}

template <typename BSE_OPERATOR_A, typename BSE_OPERATOR_B>
tools::EigenSystem BSE::Solve_nonhermitian_Davidson(BSE_OPERATOR_A& Aop,
                                                    BSE_OPERATOR_B& Bop) const {
  std::chrono::time_point<std::chrono::system_clock> start =
      std::chrono::system_clock::now();

  // operator
  HamiltonianOperator<BSE_OPERATOR_A, BSE_OPERATOR_B> Hop(Aop, Bop);

  // Davidson solver
  DavidsonSolver DS(log_);
  DS.set_correction(opt_.davidson_correction);
  DS.set_tolerance(opt_.davidson_tolerance);
  DS.set_size_update(opt_.davidson_update);
  DS.set_iter_max(opt_.davidson_maxiter);
  DS.set_max_search_space(10 * opt_.nmax);
  DS.set_matrix_type("HAM");
  DS.solve(Hop, opt_.nmax);

  // results
  tools::EigenSystem result;
  result.eigenvalues() = DS.eigenvalues();
  Eigen::MatrixXd tmpX = DS.eigenvectors().topRows(Aop.rows());
  Eigen::MatrixXd tmpY = DS.eigenvectors().bottomRows(Bop.rows());

  // // normalization so that eigenvector^2 - eigenvector2^2 = 1
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

void BSE::printFragInfo(const std::vector<QMFragment<BSE_Population> >& frags,
                        Index state) const {
  for (const QMFragment<BSE_Population>& frag : frags) {
    double dq = frag.value().H[state] + frag.value().E[state];
    double qeff = dq + frag.value().Gs;
    XTP_LOG(Log::error, log_)
        << format(
               "           Fragment %1$4d -- hole: %2$5.1f%%  electron: "
               "%3$5.1f%%  dQ: %4$+5.2f  Qeff: %5$+5.2f") %
               int(frag.getId()) % (100.0 * frag.value().H[state]) %
               (-100.0 * frag.value().E[state]) % dq % qeff
        << flush;
  }
  return;
}

void BSE::PrintWeights(const Eigen::VectorXd& weights) const {
  vc2index vc = vc2index(opt_.vmin, bse_cmin_, bse_ctotal_);
  for (Index i_bse = 0; i_bse < bse_size_; ++i_bse) {
    double weight = weights(i_bse);
    if (weight > opt_.min_print_weight) {
      XTP_LOG(Log::error, log_)
          << format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") %
                 (opt_.homo - vc.v(i_bse)) % (vc.c(i_bse) - opt_.homo - 1) %
                 (100.0 * weight)
          << flush;
    }
  }
  return;
}

void BSE::Analyze_singlets(std::vector<QMFragment<BSE_Population> > fragments,
                           const Orbitals& orb) const {

  QMStateType singlet = QMStateType(QMStateType::Singlet);

  Eigen::VectorXd oscs = orb.Oscillatorstrengths();

  Interaction act = Analyze_eh_interaction(singlet, orb);

  if (fragments.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(fragments, orb, singlet);
  }

  const Eigen::VectorXd& energies = orb.BSESinglets().eigenvalues();

  double hrt2ev = tools::conv::hrt2ev;
  XTP_LOG(Log::error, log_)
      << "  ====== singlet energies (eV) ====== " << flush;
  for (Index i = 0; i < opt_.nmax; ++i) {
    Eigen::VectorXd weights =
        orb.BSESinglets().eigenvectors().col(i).cwiseAbs2();
    if (!orb.getTDAApprox()) {
      weights -= orb.BSESinglets().eigenvectors2().col(i).cwiseAbs2();
    }

    double osc = oscs[i];
    XTP_LOG(Log::error, log_)
        << format(
               "  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> "
               "= %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f") %
               (i + 1) % (hrt2ev * energies(i)) %
               (1240.0 / (hrt2ev * energies(i))) %
               (hrt2ev * act.qp_contrib(i)) %
               (hrt2ev * act.exchange_contrib(i)) %
               (hrt2ev * act.direct_contrib(i))
        << flush;

    const Eigen::Vector3d& trdip = orb.TransitionDipoles()[i];
    XTP_LOG(Log::error, log_)
        << format(
               "           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = "
               "%2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f") %
               trdip[0] % trdip[1] % trdip[2] % (trdip.squaredNorm()) % osc
        << flush;

    PrintWeights(weights);
    if (fragments.size() > 0) {
      printFragInfo(fragments, i);
    }

    XTP_LOG(Log::error, log_) << flush;
  }
  return;
}

void BSE::Analyze_triplets(std::vector<QMFragment<BSE_Population> > fragments,
                           const Orbitals& orb) const {

  QMStateType triplet = QMStateType(QMStateType::Triplet);
  Interaction act = Analyze_eh_interaction(triplet, orb);

  if (fragments.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(fragments, orb, triplet);
  }

  const Eigen::VectorXd& energies = orb.BSETriplets().eigenvalues();
  XTP_LOG(Log::error, log_)
      << "  ====== triplet energies (eV) ====== " << flush;
  for (Index i = 0; i < opt_.nmax; ++i) {
    Eigen::VectorXd weights =
        orb.BSETriplets().eigenvectors().col(i).cwiseAbs2();
    if (!orb.getTDAApprox()) {
      weights -= orb.BSETriplets().eigenvectors2().col(i).cwiseAbs2();
    }

    XTP_LOG(Log::error, log_)
        << format(
               "  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> "
               "= %4$+1.4f <K_d> = %5$+1.4f") %
               (i + 1) % (tools::conv::hrt2ev * energies(i)) %
               (1240.0 / (tools::conv::hrt2ev * energies(i))) %
               (tools::conv::hrt2ev * act.qp_contrib(i)) %
               (tools::conv::hrt2ev * act.direct_contrib(i))
        << flush;

    PrintWeights(weights);
    if (fragments.size() > 0) {
      printFragInfo(fragments, i);
    }
    XTP_LOG(Log::error, log_) << format("   ") << flush;
  }

  return;
}
template <class OP>
Eigen::VectorXd ExpValue(const Eigen::MatrixXd& state1, OP OPxstate2) {
  return state1.cwiseProduct(OPxstate2.eval()).colwise().sum().transpose();
}

Eigen::VectorXd ExpValue(const Eigen::MatrixXd& state1,
                         const Eigen::MatrixXd& OPxstate2) {
  return state1.cwiseProduct(OPxstate2).colwise().sum().transpose();
}

template <typename BSE_OPERATOR>
BSE::ExpectationValues BSE::ExpectationValue_Operator(
    const QMStateType& type, const Orbitals& orb, const BSE_OPERATOR& H) const {

  const tools::EigenSystem& BSECoefs =
      (type == QMStateType::Singlet) ? orb.BSESinglets() : orb.BSETriplets();

  ExpectationValues expectation_values;

  const Eigen::MatrixXd temp = H * BSECoefs.eigenvectors();

  expectation_values.direct_term = ExpValue(BSECoefs.eigenvectors(), temp);
  if (!orb.getTDAApprox()) {
    expectation_values.direct_term +=
        ExpValue(BSECoefs.eigenvectors2(), H * BSECoefs.eigenvectors2());
    expectation_values.cross_term =
        2 * ExpValue(BSECoefs.eigenvectors2(), temp);
  } else {
    expectation_values.cross_term = Eigen::VectorXd::Zero(0);
  }
  return expectation_values;
}

template <typename BSE_OPERATOR>
BSE::ExpectationValues BSE::ExpectationValue_Operator_State(
    const QMState& state, const Orbitals& orb, const BSE_OPERATOR& H) const {

  const tools::EigenSystem& BSECoefs = (state.Type() == QMStateType::Singlet)
                                           ? orb.BSESinglets()
                                           : orb.BSETriplets();

  ExpectationValues expectation_values;

  const Eigen::MatrixXd BSECoefs_state =
      BSECoefs.eigenvectors().col(state.StateIdx());

  const Eigen::MatrixXd temp = H * BSECoefs_state;

  expectation_values.direct_term = ExpValue(BSECoefs_state, temp);

  if (!orb.getTDAApprox()) {
    const Eigen::MatrixXd BSECoefs2_state =
        BSECoefs.eigenvectors2().col(state.StateIdx());

    expectation_values.direct_term +=
        ExpValue(BSECoefs2_state, H * BSECoefs2_state);
    expectation_values.cross_term = 2 * ExpValue(BSECoefs2_state, temp);
  } else {
    expectation_values.cross_term = Eigen::VectorXd::Zero(0);
  }

  return expectation_values;
}

// Composition of the excitation energy in terms of QP, direct (screened),
// and exchance contributions in the BSE
// Full BSE:
//
// |  A* | |  H  K | | A |
// | -B* | | -K -H | | B | = A*.H.A + B*.H.B + 2A*.K.B
//
// with: H = H_qp + H_d  + eta.H_x
//       K =        H_d2 + eta.H_x
//
// reports composition for FULL BSE as
//  <FT> = A*.H_qp.A + B*.H_qp.B
//  <Kx> = eta.(A*.H_x.A + B*.H_x.B + 2A*.H_x.B)
//  <Kd> = A*.H_d.A + B*.H_d.B + 2A*.H_d2.B
BSE::Interaction BSE::Analyze_eh_interaction(const QMStateType& type,
                                             const Orbitals& orb) const {
  Interaction analysis;
  {
    HqpOperator hqp(epsilon_0_inv_, Mmn_, Hqp_);
    configureBSEOperator(hqp);
    ExpectationValues expectation_values =
        ExpectationValue_Operator(type, orb, hqp);
    analysis.qp_contrib = expectation_values.direct_term;
  }
  {
    HdOperator hd(epsilon_0_inv_, Mmn_, Hqp_);
    configureBSEOperator(hd);
    ExpectationValues expectation_values =
        ExpectationValue_Operator(type, orb, hd);
    analysis.direct_contrib = expectation_values.direct_term;
  }
  if (!orb.getTDAApprox()) {
    Hd2Operator hd2(epsilon_0_inv_, Mmn_, Hqp_);
    configureBSEOperator(hd2);
    ExpectationValues expectation_values =
        ExpectationValue_Operator(type, orb, hd2);
    analysis.direct_contrib += expectation_values.cross_term;
  }

  if (type == QMStateType::Singlet) {
    HxOperator hx(epsilon_0_inv_, Mmn_, Hqp_);
    configureBSEOperator(hx);
    ExpectationValues expectation_values =
        ExpectationValue_Operator(type, orb, hx);
    analysis.exchange_contrib = 2.0 * expectation_values.direct_term;
    if (!orb.getTDAApprox()) {
      analysis.exchange_contrib += 2.0 * expectation_values.cross_term;
    }
  } else {
    analysis.exchange_contrib = Eigen::VectorXd::Zero(0);
  }

  return analysis;
}

// Dynamical Screening in BSE as perturbation to static excitation energies
// as in Phys. Rev. B 80, 241405 (2009) for the TDA case
void BSE::Perturbative_DynamicalScreening(const QMStateType& type,
                                          Orbitals& orb) {

  const tools::EigenSystem& BSECoefs =
      (type == QMStateType::Singlet) ? orb.BSESinglets() : orb.BSETriplets();

  const Eigen::VectorXd& RPAInputEnergies = orb.RPAInputEnergies();

  // static case as reference
  SetupDirectInteractionOperator(RPAInputEnergies, 0.0);
  HdOperator Hd_static(epsilon_0_inv_, Mmn_, Hqp_);
  configureBSEOperator(Hd_static);
  ExpectationValues expectation_values =
      ExpectationValue_Operator(type, orb, Hd_static);
  Eigen::VectorXd Hd_static_contribution = expectation_values.direct_term;
  if (!orb.getTDAApprox()) {
    Hd2Operator Hd2_static(epsilon_0_inv_, Mmn_, Hqp_);
    configureBSEOperator(Hd2_static);
    expectation_values = ExpectationValue_Operator(type, orb, Hd2_static);
    Hd_static_contribution += expectation_values.cross_term;
  }

  const Eigen::VectorXd& BSEenergies = BSECoefs.eigenvalues();

  // initial copy of static BSE energies to dynamic
  Eigen::VectorXd BSEenergies_dynamic = BSEenergies;

  // recalculate Hd at the various energies
  for (Index i_exc = 0; i_exc < BSEenergies.size(); i_exc++) {
    XTP_LOG(Log::info, log_) << "Dynamical Screening BSE, Excitation " << i_exc
                             << " static " << BSEenergies(i_exc) << flush;

    for (Index iter = 0; iter < max_dyn_iter_; iter++) {

      // setup the direct operator with the last energy as screening frequency
      double old_energy = BSEenergies_dynamic(i_exc);
      SetupDirectInteractionOperator(RPAInputEnergies, old_energy);
      HdOperator Hd_dyn(epsilon_0_inv_, Mmn_, Hqp_);
      configureBSEOperator(Hd_dyn);

      // get the contribution of Hd for the dynamic case
      QMState state(type, i_exc, false);
      expectation_values = ExpectationValue_Operator_State(state, orb, Hd_dyn);
      Eigen::VectorXd Hd_dynamic_contribution = expectation_values.direct_term;
      if (!orb.getTDAApprox()) {
        Hd2Operator Hd2_dyn(epsilon_0_inv_, Mmn_, Hqp_);
        configureBSEOperator(Hd2_dyn);
        expectation_values =
            ExpectationValue_Operator_State(state, orb, Hd2_dyn);
        Hd_dynamic_contribution += expectation_values.cross_term;
      }

      // new energy perturbatively
      BSEenergies_dynamic(i_exc) = BSEenergies(i_exc) +
                                   Hd_static_contribution(i_exc) -
                                   Hd_dynamic_contribution(0);

      XTP_LOG(Log::info, log_)
          << "Dynamical Screening BSE, excitation " << i_exc << " iteration "
          << iter << " dynamic " << BSEenergies_dynamic(i_exc) << flush;

      // check tolerance
      if (std::abs(BSEenergies_dynamic(i_exc) - old_energy) < dyn_tolerance_) {
        break;
      }
    }
  }

  double hrt2ev = tools::conv::hrt2ev;

  if (type == QMStateType::Singlet) {
    orb.BSESinglets_dynamic() = BSEenergies_dynamic;
    XTP_LOG(Log::error, log_) << "  ====== singlet energies with perturbative "
                                 "dynamical screening (eV) ====== "
                              << flush;
    Eigen::VectorXd oscs = orb.Oscillatorstrengths();
    for (Index i = 0; i < opt_.nmax; ++i) {
      double osc = oscs[i];
      XTP_LOG(Log::error, log_)
          << format(
                 "  S(dynamic) = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f "
                 "nm f "
                 "= %4$+1.4f") %
                 (i + 1) % (hrt2ev * BSEenergies_dynamic(i)) %
                 (1240.0 / (hrt2ev * BSEenergies_dynamic(i))) %
                 (osc * BSEenergies_dynamic(i) / BSEenergies(i))
          << flush;
    }

  } else {
    orb.BSETriplets_dynamic() = BSEenergies_dynamic;
    XTP_LOG(Log::error, log_) << "  ====== triplet energies with perturbative "
                                 "dynamical screening (eV) ====== "
                              << flush;
    for (Index i = 0; i < opt_.nmax; ++i) {
      XTP_LOG(Log::error, log_)
          << format(
                 "  T(dynamic) = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f "
                 "nm ") %
                 (i + 1) % (hrt2ev * BSEenergies_dynamic(i)) %
                 (1240.0 / (hrt2ev * BSEenergies_dynamic(i)))
          << flush;
    }
  }
}

}  // namespace xtp
}  // namespace votca
