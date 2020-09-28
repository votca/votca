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
  _opt = opt;
  _bse_vmax = _opt.homo;
  _bse_cmin = _opt.homo + 1;
  _bse_vtotal = _bse_vmax - _opt.vmin + 1;
  _bse_ctotal = _opt.cmax - _bse_cmin + 1;
  _bse_size = _bse_vtotal * _bse_ctotal;
  _max_dyn_iter = _opt.max_dyn_iter;
  _dyn_tolerance = _opt.dyn_tolerance;
  if (_opt.use_Hqp_offdiag) {
    _Hqp = AdjustHqpSize(Hqp_in, RPAInputEnergies);
  } else {
    _Hqp = AdjustHqpSize(Hqp_in, RPAInputEnergies).diagonal().asDiagonal();
  }
  SetupDirectInteractionOperator(RPAInputEnergies, 0.0);
}

Eigen::MatrixXd BSE::AdjustHqpSize(const Eigen::MatrixXd& Hqp,
                                   const Eigen::VectorXd& RPAInputEnergies) {

  Index hqp_size = _bse_vtotal + _bse_ctotal;
  Index gwsize = _opt.qpmax - _opt.qpmin + 1;
  Index RPAoffset = _opt.vmin - _opt.rpamin;
  Eigen::MatrixXd Hqp_BSE = Eigen::MatrixXd::Zero(hqp_size, hqp_size);

  if (_opt.vmin >= _opt.qpmin) {
    Index start = _opt.vmin - _opt.qpmin;
    if (_opt.cmax <= _opt.qpmax) {
      Hqp_BSE = Hqp.block(start, start, hqp_size, hqp_size);
    } else {
      Index virtoffset = gwsize - start;
      Hqp_BSE.topLeftCorner(virtoffset, virtoffset) =
          Hqp.block(start, start, virtoffset, virtoffset);

      Index virt_extra = _opt.cmax - _opt.qpmax;
      Hqp_BSE.diagonal().tail(virt_extra) =
          RPAInputEnergies.segment(RPAoffset + virtoffset, virt_extra);
    }
  }

  if (_opt.vmin < _opt.qpmin) {
    Index occ_extra = _opt.qpmin - _opt.vmin;
    Hqp_BSE.diagonal().head(occ_extra) =
        RPAInputEnergies.segment(RPAoffset, occ_extra);

    Hqp_BSE.block(occ_extra, occ_extra, gwsize, gwsize) = Hqp;

    if (_opt.cmax > _opt.qpmax) {
      Index virtoffset = occ_extra + gwsize;
      Index virt_extra = _opt.cmax - _opt.qpmax;
      Hqp_BSE.diagonal().tail(virt_extra) =
          RPAInputEnergies.segment(RPAoffset + virtoffset, virt_extra);
    }
  }

  return Hqp_BSE;
}

void BSE::SetupDirectInteractionOperator(
    const Eigen::VectorXd& RPAInputEnergies, double energy) {
  RPA rpa = RPA(_log, _Mmn);
  rpa.configure(_opt.homo, _opt.rpamin, _opt.rpamax);
  rpa.setRPAInputEnergies(RPAInputEnergies);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(
      rpa.calculate_epsilon_r(energy));
  _Mmn.MultiplyRightWithAuxMatrix(es.eigenvectors());

  _epsilon_0_inv = Eigen::VectorXd::Zero(es.eigenvalues().size());
  for (Index i = 0; i < es.eigenvalues().size(); ++i) {
    if (es.eigenvalues()(i) > 1e-8) {
      _epsilon_0_inv(i) = 1 / es.eigenvalues()(i);
    }
  }
}

template <typename BSE_OPERATOR>
void BSE::configureBSEOperator(BSE_OPERATOR& H) const {
  BSEOperator_Options opt;
  opt.cmax = _opt.cmax;
  opt.homo = _opt.homo;
  opt.qpmin = _opt.qpmin;
  opt.rpamin = _opt.rpamin;
  opt.vmin = _opt.vmin;
  H.configure(opt);
}

tools::EigenSystem BSE::Solve_triplets_TDA() const {

  TripletOperator_TDA Ht(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Ht);
  return solve_hermitian(Ht);
}

void BSE::Solve_singlets(Orbitals& orb) const {
  orb.setTDAApprox(_opt.useTDA);
  if (_opt.useTDA) {
    orb.BSESinglets() = Solve_singlets_TDA();
  } else {
    orb.BSESinglets() = Solve_singlets_BTDA();
  }
  orb.CalcCoupledTransition_Dipoles();
}

void BSE::Solve_triplets(Orbitals& orb) const {
  orb.setTDAApprox(_opt.useTDA);
  if (_opt.useTDA) {
    orb.BSETriplets() = Solve_triplets_TDA();
  } else {
    orb.BSETriplets() = Solve_triplets_BTDA();
  }
}

tools::EigenSystem BSE::Solve_singlets_TDA() const {

  SingletOperator_TDA Hs(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hs);
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Setup TDA singlet hamiltonian " << flush;
  return solve_hermitian(Hs);
}

SingletOperator_TDA BSE::getSingletOperator_TDA() const {

  SingletOperator_TDA Hs(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hs);
  return Hs;
}

TripletOperator_TDA BSE::getTripletOperator_TDA() const {

  TripletOperator_TDA Ht(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Ht);
  return Ht;
}

template <typename BSE_OPERATOR>
tools::EigenSystem BSE::solve_hermitian(BSE_OPERATOR& h) const {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::time_point<std::chrono::system_clock> hstart, hend;
  std::chrono::duration<double> elapsed_time;
  start = std::chrono::system_clock::now();

  tools::EigenSystem result;
  if (_opt.davidson) {

    DavidsonSolver DS(_log);

    DS.set_correction(_opt.davidson_correction);
    DS.set_tolerance(_opt.davidson_tolerance);
    DS.set_ortho(_opt.davidson_ortho);
    DS.set_size_update(_opt.davidson_update);
    DS.set_iter_max(_opt.davidson_maxiter);
    DS.set_max_search_space(10 * _opt.nmax);

    if (_opt.matrixfree) {
      XTP_LOG(Log::error, _log)
          << TimeStamp() << " Using matrix free method" << flush;
      DS.solve(h, _opt.nmax);
    } else {
      XTP_LOG(Log::error, _log)
          << TimeStamp() << " Using full matrix method" << flush;

      // get the full matrix
      hstart = std::chrono::system_clock::now();
      Eigen::MatrixXd hfull = h.get_full_matrix();
      hend = std::chrono::system_clock::now();

      elapsed_time = hend - hstart;

      XTP_LOG(Log::info, _log) << TimeStamp() << " Full matrix assembled in "
                               << elapsed_time.count() << " secs" << flush;

      // solve theeigenalue problem
      hstart = std::chrono::system_clock::now();
      DS.solve(hfull, _opt.nmax);
      hend = std::chrono::system_clock::now();

      elapsed_time = hend - hstart;
      XTP_LOG(Log::info, _log) << TimeStamp() << " Davidson solve done in "
                               << elapsed_time.count() << " secs" << flush;
    }
    result.eigenvalues() = DS.eigenvalues();
    result.eigenvectors() = DS.eigenvectors();

  } else {

    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Lapack Diagonalization" << flush;

    hstart = std::chrono::system_clock::now();
    Eigen::MatrixXd hfull = h.get_full_matrix();
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;
    XTP_LOG(Log::info, _log) << TimeStamp() << " Full matrix assembled in "
                             << elapsed_time.count() << " secs" << flush;
    hstart = std::chrono::system_clock::now();
    result = tools::linalg_eigenvalues(hfull, _opt.nmax);
    hend = std::chrono::system_clock::now();
    elapsed_time = hend - hstart;
    XTP_LOG(Log::info, _log) << TimeStamp() << " Lapack solve done in "
                             << elapsed_time.count() << " secs" << flush;
  }
  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  XTP_LOG(Log::info, _log) << TimeStamp() << " Diagonalization done in "
                           << elapsed_time.count() << " secs" << flush;

  return result;
}

tools::EigenSystem BSE::Solve_singlets_BTDA() const {
  if (_opt.davidson) {
    SingletOperator_TDA A(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(A);

    SingletOperator_BTDA_B B(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(B);

    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Setup Full singlet hamiltonian " << flush;
    return Solve_nonhermitian_Davidson(A, B);
  } else {
    SingletOperator_BTDA_ApB Hs_ApB(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(Hs_ApB);
    Operator_BTDA_AmB Hs_AmB(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(Hs_AmB);
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Setup Full singlet hamiltonian " << flush;
    return Solve_nonhermitian(Hs_ApB, Hs_AmB);
  }
}

tools::EigenSystem BSE::Solve_triplets_BTDA() const {
  if (_opt.davidson) {
    TripletOperator_TDA A(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(A);
    Hd2Operator B(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(B);
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Setup Full triplet hamiltonian " << flush;
    return Solve_nonhermitian_Davidson(A, B);
  } else {
    TripletOperator_BTDA_ApB Ht_ApB(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(Ht_ApB);
    Operator_BTDA_AmB Ht_AmB(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(Ht_AmB);
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Setup Full triplet hamiltonian " << flush;
    return Solve_nonhermitian(Ht_ApB, Ht_AmB);
  }
}

template <typename BSE_OPERATOR_ApB, typename BSE_OPERATOR_AmB>
tools::EigenSystem BSE::Solve_nonhermitian(BSE_OPERATOR_ApB& apb,
                                           BSE_OPERATOR_AmB& amb) const {

  // For details of the method, see EPL,78(2007)12001,
  // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.
  // setup resonant (A) and RARC blocks (B)
  // corresponds to
  // _ApB = (_eh_d +_eh_qp + _eh_d2 + 4.0 * _eh_x);
  // _AmB = (_eh_d +_eh_qp - _eh_d2);

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::time_point<std::chrono::system_clock> hstart, hend;
  std::chrono::duration<double> elapsed_time;
  start = std::chrono::system_clock::now();

  Eigen::MatrixXd ApB = apb.get_full_matrix();
  Eigen::MatrixXd AmB = amb.get_full_matrix();
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Setup singlet hamiltonian " << flush;
  XTP_LOG(Log::error, _log)
      << TimeStamp() << " Lapack Diagonalization" << flush;

  // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not
  // positive definite
  //(A-B) is not needed any longer and can be overwritten
  XTP_LOG(Log::info, _log) << TimeStamp()
                           << " Trying Cholesky decomposition of KAA-KAB"
                           << flush;
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd> > L(AmB);

  for (Index i = 0; i < AmB.rows(); ++i) {
    for (Index j = i + 1; j < AmB.cols(); ++j) {
      AmB(i, j) = 0;
    }
  }
  if (L.info() != Eigen::ComputationInfo::Success) {
    XTP_LOG(Log::error, _log)
        << TimeStamp()
        << " Cholesky decomposition of KAA-KAB was unsucessful. Try a smaller "
           "basisset. This can indicate a triplet instability."
        << flush;
    throw std::runtime_error("Cholesky decompostion failed");
  } else {
    XTP_LOG(Log::info, _log)
        << TimeStamp() << " Cholesky decomposition of KAA-KAB was successful"
        << flush;
  }

  ApB = AmB.transpose() * ApB * AmB;

  XTP_LOG(Log::info, _log) << TimeStamp() << " Calculated H = L^T(A+B)L "
                           << flush;

  XTP_LOG(Log::error, _log) << TimeStamp() << " Solving for first " << _opt.nmax
                            << " eigenvectors" << flush;

  hstart = std::chrono::system_clock::now();
  tools::EigenSystem eigensys = tools::linalg_eigenvalues(ApB, _opt.nmax);
  hend = std::chrono::system_clock::now();
  elapsed_time = hend - hstart;
  XTP_LOG(Log::info, _log) << TimeStamp() << " Lapack solve done in "
                           << elapsed_time.count() << " secs" << flush;

  if (eigensys.info() != Eigen::ComputationInfo::Success) {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Could not solve problem" << flush;
  } else {
    XTP_LOG(Log::info, _log)
        << TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
  }
  ApB.resize(0, 0);
  if ((eigensys.eigenvalues().array() < 0).any()) {
    throw std::runtime_error("Negative eigenvalues in BTDA");
  }
  Eigen::VectorXd energies =
      eigensys.eigenvalues().cwiseSqrt();  // has to stay otherwise mkl
                                           // complains

  tools::EigenSystem result;
  result.eigenvalues() = energies;
  // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 +
  // 1/sqrt(eps_l)L ] R_l
  //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 -
  //                               1/sqrt(eps_l)L ] R_l
  //                               1/sqrt(eps_l)L ] R_l
  // determine inverse of L^T
  Eigen::MatrixXd LmT = AmB.inverse().transpose();
  Index dim = LmT.rows();
  result.eigenvectors().resize(dim, _opt.nmax);  // resonant part (_X_evec)
  result.eigenvectors2().resize(dim,
                                _opt.nmax);  // anti-resonant part (_Y_evec)
  for (Index level = 0; level < _opt.nmax; level++) {
    double sqrt_eval = std::sqrt(energies(level));
    // get l-th reduced EV
    result.eigenvectors().col(level) =
        (0.5 / sqrt_eval * (energies(level) * LmT + AmB) *
         eigensys.eigenvectors().col(level));
    result.eigenvectors2().col(level) =
        (0.5 / sqrt_eval * (energies(level) * LmT - AmB) *
         eigensys.eigenvectors().col(level));
  }

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  XTP_LOG(Log::info, _log) << TimeStamp() << " Diagonalization done in "
                           << elapsed_time.count() << " secs" << flush;

  return result;
}

template <typename BSE_OPERATOR_A, typename BSE_OPERATOR_B>
tools::EigenSystem BSE::Solve_nonhermitian_Davidson(BSE_OPERATOR_A& Aop,
                                                    BSE_OPERATOR_B& Bop) const {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::time_point<std::chrono::system_clock> hstart, hend;
  std::chrono::duration<double> elapsed_time;
  start = std::chrono::system_clock::now();

  // operator
  HamiltonianOperator<BSE_OPERATOR_A, BSE_OPERATOR_B> Hop(Aop, Bop);

  // Davidson solver
  DavidsonSolver DS(_log);
  DS.set_correction(_opt.davidson_correction);
  DS.set_tolerance(_opt.davidson_tolerance);
  DS.set_ortho(_opt.davidson_ortho);
  DS.set_size_update(_opt.davidson_update);
  DS.set_iter_max(_opt.davidson_maxiter);
  DS.set_max_search_space(10 * _opt.nmax);
  DS.set_matrix_type("HAM");

  if (_opt.matrixfree) {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Using matrix free method" << flush;
    DS.solve(Hop, _opt.nmax);

  } else {
    XTP_LOG(Log::error, _log)
        << TimeStamp() << " Using full matrix method" << flush;
    // get the full matrix
    hstart = std::chrono::system_clock::now();
    Eigen::MatrixXd hfull = Hop.get_full_matrix();
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;

    XTP_LOG(Log::info, _log) << TimeStamp() << " Full matrix assembled in "
                             << elapsed_time.count() << " secs" << flush;

    // solve theeigenalue problem
    hstart = std::chrono::system_clock::now();
    DS.solve(hfull, _opt.nmax);
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;
    XTP_LOG(Log::info, _log) << TimeStamp() << " Davidson solve done in "
                             << elapsed_time.count() << " secs" << flush;
  }

  // results
  tools::EigenSystem result;
  result.eigenvalues() = DS.eigenvalues();
  Eigen::MatrixXd _tmpX = DS.eigenvectors().topRows(Aop.size());
  Eigen::MatrixXd _tmpY = DS.eigenvectors().bottomRows(Aop.size());

  // // normalization so that eigenvector^2 - eigenvector2^2 = 1
  Eigen::VectorXd normX = _tmpX.colwise().squaredNorm();
  Eigen::VectorXd normY = _tmpY.colwise().squaredNorm();

  Eigen::ArrayXd sqinvnorm = (normX - normY).array().inverse().cwiseSqrt();

  result.eigenvectors() = _tmpX * sqinvnorm.matrix().asDiagonal();
  result.eigenvectors2() = _tmpY * sqinvnorm.matrix().asDiagonal();

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  XTP_LOG(Log::info, _log) << TimeStamp() << " Diagonalization done in "
                           << elapsed_time.count() << " secs" << flush;

  return result;
}

void BSE::printFragInfo(const std::vector<QMFragment<BSE_Population> >& frags,
                        Index state) const {
  for (const QMFragment<BSE_Population>& frag : frags) {
    double dq = frag.value().H[state] + frag.value().E[state];
    double qeff = dq + frag.value().Gs;
    XTP_LOG(Log::error, _log)
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
  vc2index vc = vc2index(_opt.vmin, _bse_cmin, _bse_ctotal);
  for (Index i_bse = 0; i_bse < _bse_size; ++i_bse) {
    double weight = weights(i_bse);
    if (weight > _opt.min_print_weight) {
      XTP_LOG(Log::error, _log)
          << format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") %
                 (_opt.homo - vc.v(i_bse)) % (vc.c(i_bse) - _opt.homo - 1) %
                 (100.0 * weight)
          << flush;
    }
  }
  return;
}

void BSE::Analyze_singlets(std::vector<QMFragment<BSE_Population> > fragments,
                           const Orbitals& orb) const {

  Interaction act;
  QMStateType singlet = QMStateType(QMStateType::Singlet);

  Eigen::VectorXd oscs = orb.Oscillatorstrengths();

  act = Analyze_eh_interaction(singlet, orb);

  if (fragments.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(fragments, orb, singlet);
  }

  const Eigen::VectorXd& energies = orb.BSESinglets().eigenvalues();

  double hrt2ev = tools::conv::hrt2ev;
  XTP_LOG(Log::error, _log)
      << "  ====== singlet energies (eV) ====== " << flush;
  for (Index i = 0; i < _opt.nmax; ++i) {
    Eigen::VectorXd weights =
        orb.BSESinglets().eigenvectors().col(i).cwiseAbs2();
    if (!orb.getTDAApprox()) {
      weights -= orb.BSESinglets().eigenvectors2().col(i).cwiseAbs2();
    }

    double osc = oscs[i];
    XTP_LOG(Log::error, _log)
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
    XTP_LOG(Log::error, _log)
        << format(
               "           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = "
               "%2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f") %
               trdip[0] % trdip[1] % trdip[2] % (trdip.squaredNorm()) % osc
        << flush;

    PrintWeights(weights);
    if (fragments.size() > 0) {
      printFragInfo(fragments, i);
    }

    XTP_LOG(Log::error, _log) << flush;
  }
  return;
}

void BSE::Analyze_triplets(std::vector<QMFragment<BSE_Population> > fragments,
                           const Orbitals& orb) const {

  Interaction act;
  QMStateType triplet = QMStateType(QMStateType::Triplet);
  act = Analyze_eh_interaction(triplet, orb);

  if (fragments.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(fragments, orb, triplet);
  }

  const Eigen::VectorXd& energies = orb.BSETriplets().eigenvalues();
  XTP_LOG(Log::error, _log)
      << "  ====== triplet energies (eV) ====== " << flush;
  for (Index i = 0; i < _opt.nmax; ++i) {
    Eigen::VectorXd weights =
        orb.BSETriplets().eigenvectors().col(i).cwiseAbs2();
    if (!orb.getTDAApprox()) {
      weights -= orb.BSETriplets().eigenvectors2().col(i).cwiseAbs2();
    }

    XTP_LOG(Log::error, _log)
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
    XTP_LOG(Log::error, _log) << format("   ") << flush;
  }

  return;
}

template <typename BSE_OPERATOR>
BSE::ExpectationValues BSE::ExpectationValue_Operator(
    const QMStateType& type, const Orbitals& orb, const BSE_OPERATOR& H) const {

  const tools::EigenSystem& BSECoefs =
      (type == QMStateType::Singlet) ? orb.BSESinglets() : orb.BSETriplets();

  ExpectationValues expectation_values;

  expectation_values.direct_term =
      BSECoefs.eigenvectors()
          .cwiseProduct((H * BSECoefs.eigenvectors()))
          .colwise()
          .sum()
          .transpose();

  if (!orb.getTDAApprox()) {
    expectation_values.direct_term +=
        BSECoefs.eigenvectors2()
            .cwiseProduct((H * BSECoefs.eigenvectors2()))
            .colwise()
            .sum()
            .transpose();

    expectation_values.cross_term =
        2 * BSECoefs.eigenvectors2()
                .cwiseProduct((H * BSECoefs.eigenvectors()))
                .colwise()
                .sum()
                .transpose();
  } else {
    expectation_values.cross_term = Eigen::VectorXd::Zero(0);
  }

  return expectation_values;
}

template <typename BSE_OPERATOR>
BSE::ExpectationValues BSE::ExpectationValue_Operator_State(
    const QMStateType& type, const Orbitals& orb, const BSE_OPERATOR& H,
    const Index state) const {

  const tools::EigenSystem& BSECoefs =
      (type == QMStateType::Singlet) ? orb.BSESinglets() : orb.BSETriplets();

  ExpectationValues expectation_values;

  const Eigen::MatrixXd& BSECoefs_state = BSECoefs.eigenvectors().col(state);

  expectation_values.direct_term =
      BSECoefs_state.cwiseProduct((H * BSECoefs_state))
          .colwise()
          .sum()
          .transpose();

  if (!orb.getTDAApprox()) {
    const Eigen::MatrixXd& BSECoefs2_state =
        BSECoefs.eigenvectors2().col(state);

    expectation_values.direct_term +=
        BSECoefs2_state.cwiseProduct((H * BSECoefs2_state))
            .colwise()
            .sum()
            .transpose();

    expectation_values.cross_term =
        2 * BSECoefs2_state.cwiseProduct((H * BSECoefs_state))
                .colwise()
                .sum()
                .transpose();
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
  ExpectationValues expectation_values;

  HqpOperator hqp(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(hqp);
  expectation_values = ExpectationValue_Operator(type, orb, hqp);
  analysis.qp_contrib = expectation_values.direct_term;

  HdOperator hd(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(hd);
  expectation_values = ExpectationValue_Operator(type, orb, hd);
  analysis.direct_contrib = expectation_values.direct_term;

  if (!orb.getTDAApprox()) {
    Hd2Operator hd2(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(hd2);
    expectation_values = ExpectationValue_Operator(type, orb, hd2);
    analysis.direct_contrib += expectation_values.cross_term;
  }

  if (type == QMStateType::Singlet) {
    HxOperator hx(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(hx);
    expectation_values = ExpectationValue_Operator(type, orb, hx);
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

  ExpectationValues expectation_values;

  // static case as reference
  SetupDirectInteractionOperator(RPAInputEnergies, 0.0);
  HdOperator Hd_static(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hd_static);
  expectation_values = ExpectationValue_Operator(type, orb, Hd_static);
  Eigen::VectorXd Hd_static_contribution = expectation_values.direct_term;
  if (!orb.getTDAApprox()) {
    Hd2Operator Hd2_static(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(Hd2_static);
    expectation_values = ExpectationValue_Operator(type, orb, Hd2_static);
    Hd_static_contribution += expectation_values.cross_term;
  }

  const Eigen::VectorXd& BSEenergies = BSECoefs.eigenvalues();

  // initial copy of static BSE energies to dynamic
  Eigen::VectorXd BSEenergies_dynamic = BSEenergies;

  // recalculate Hd at the various energies
  for (Index i_exc = 0; i_exc < BSEenergies.size(); i_exc++) {
    XTP_LOG(Log::info, _log) << "Dynamical Screening BSE, Excitation " << i_exc
                             << " static " << BSEenergies(i_exc) << flush;

    for (Index iter = 0; iter < _max_dyn_iter; iter++) {

      // setup the direct operator with the last energy as screening frequency
      double old_energy = BSEenergies_dynamic(i_exc);
      SetupDirectInteractionOperator(RPAInputEnergies, old_energy);
      HdOperator Hd_dyn(_epsilon_0_inv, _Mmn, _Hqp);
      configureBSEOperator(Hd_dyn);

      // get the contribution of Hd for the dynamic case
      expectation_values =
          ExpectationValue_Operator_State(type, orb, Hd_dyn, i_exc);
      Eigen::VectorXd Hd_dynamic_contribution = expectation_values.direct_term;
      if (!orb.getTDAApprox()) {
        Hd2Operator Hd2_dyn(_epsilon_0_inv, _Mmn, _Hqp);
        configureBSEOperator(Hd2_dyn);
        expectation_values =
            ExpectationValue_Operator_State(type, orb, Hd2_dyn, i_exc);
        Hd_dynamic_contribution += expectation_values.cross_term;
      }

      // new energy perturbatively
      BSEenergies_dynamic(i_exc) = BSEenergies(i_exc) +
                                   Hd_static_contribution(i_exc) -
                                   Hd_dynamic_contribution(0);

      XTP_LOG(Log::info, _log)
          << "Dynamical Screening BSE, excitation " << i_exc << " iteration "
          << iter << " dynamic " << BSEenergies_dynamic(i_exc) << flush;

      // check tolerance
      if (std::abs(BSEenergies_dynamic(i_exc) - old_energy) < _dyn_tolerance) {
        break;
      }
    }
  }

  double hrt2ev = tools::conv::hrt2ev;

  if (type == QMStateType::Singlet) {
    orb.BSESinglets_dynamic() = BSEenergies_dynamic;
    XTP_LOG(Log::error, _log) << "  ====== singlet energies with perturbative "
                                 "dynamical screening (eV) ====== "
                              << flush;
    Eigen::VectorXd oscs = orb.Oscillatorstrengths();
    for (Index i = 0; i < _opt.nmax; ++i) {
      double osc = oscs[i];
      XTP_LOG(Log::error, _log)
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
    XTP_LOG(Log::error, _log) << "  ====== triplet energies with perturbative "
                                 "dynamical screening (eV) ====== "
                              << flush;
    for (Index i = 0; i < _opt.nmax; ++i) {
      XTP_LOG(Log::error, _log)
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
