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

#include "votca/xtp/vc2index.h"
#include <iostream>
#include <votca/tools/linalg.h>
#include <votca/xtp/bse.h>
#include <votca/xtp/bse_operator.h>
#include <votca/xtp/davidsonsolver.h>
#include <votca/xtp/populationanalysis.h>
#include <votca/xtp/qmfragment.h>
#include <votca/xtp/rpa.h>

#include <chrono>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

void BSE::configure(const options& opt) {
  _opt = opt;
  _bse_vmax = _opt.homo;
  _bse_cmin = _opt.homo + 1;
  _bse_vtotal = _bse_vmax - _opt.vmin + 1;
  _bse_ctotal = _opt.cmax - _bse_cmin + 1;
  _bse_size = _bse_vtotal * _bse_ctotal;
  SetupDirectInteractionOperator();
}

void BSE::SetupDirectInteractionOperator() {
  RPA rpa = RPA(_Mmn);
  rpa.configure(_opt.homo, _opt.rpamin, _opt.rpamax);
  rpa.UpdateRPAInputEnergies(_orbitals.MOEnergies(), _Hqp.diagonal(),
                             _opt.qpmin);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(rpa.calculate_epsilon_r(0));
  _Mmn.MultiplyRightWithAuxMatrix(es.eigenvectors());

  _epsilon_0_inv = Eigen::VectorXd::Zero(es.eigenvalues().size());
  for (int i = 0; i < es.eigenvalues().size(); ++i) {
    if (es.eigenvalues()(i) > 1e-8) {
      _epsilon_0_inv(i) = 1 / es.eigenvalues()(i);
    }
  }
}

template <typename BSE_OPERATOR>
void BSE::configureBSEOperator(BSE_OPERATOR& H) {
  BSEOperator_Options opt;
  opt.cmax = _opt.cmax;
  opt.homo = _opt.homo;
  opt.qpmin = _opt.qpmin;
  opt.rpamin = _opt.rpamin;
  opt.vmin = _opt.vmin;
  H.configure(opt);
}

void BSE::Solve_triplets_TDA() {

  TripletOperator_TDA Ht(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Ht);
  solve_hermitian(Ht, _bse_triplet_energies, _bse_triplet_coefficients);

  return;
}

void BSE::Solve_singlets() {
  if (_opt.useTDA) {
    Solve_singlets_TDA();
  } else {
    if (_opt.davidson) {
      XTP_LOG(logDEBUG, _log)
          << TimeStamp()
          << " Davidson solver not implemented for BTDA. Full diagonalization."
          << flush;
      _opt.davidson = 0;
    }
    Solve_singlets_BTDA();
  }
}

void BSE::Solve_triplets() {
  if (_opt.useTDA) {
    Solve_triplets_TDA();
  } else {
    if (_opt.davidson) {
      XTP_LOG(logDEBUG, _log)
          << TimeStamp()
          << " Davidson solver not implemented for BTDA. Full diagonalization."
          << flush;
      _opt.davidson = 0;
    }
    Solve_triplets_BTDA();
  }
}

void BSE::Solve_singlets_TDA() {

  SingletOperator_TDA Hs(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hs);
  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Setup TDA singlet hamiltonian "
                          << flush;

  solve_hermitian(Hs, _bse_singlet_energies, _bse_singlet_coefficients);
}

SingletOperator_TDA BSE::getSingletOperator_TDA() {

  SingletOperator_TDA Hs(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hs);
  return Hs;
}

TripletOperator_TDA BSE::getTripletOperator_TDA() {

  TripletOperator_TDA Ht(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Ht);
  return Ht;
}

template <typename BSE_OPERATOR>
void BSE::solve_hermitian(BSE_OPERATOR& h, Eigen::VectorXd& energies,
                          Eigen::MatrixXd& coefficients) {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::time_point<std::chrono::system_clock> hstart, hend;
  std::chrono::duration<double> elapsed_time;
  start = std::chrono::system_clock::now();

  if (_opt.davidson) {

    DavidsonSolver DS(_log);

    DS.set_correction(_opt.davidson_correction);
    DS.set_tolerance(_opt.davidson_tolerance);
    DS.set_ortho(_opt.davidson_ortho);
    DS.set_size_update(_opt.davidson_update);
    DS.set_iter_max(_opt.davidson_maxiter);
    DS.set_max_search_space(10 * _opt.nmax);

    if (_opt.matrixfree) {
      XTP_LOG(logDEBUG, _log)
          << TimeStamp() << " Using matrix free method" << flush;
      DS.solve(h, _opt.nmax);
    }

    else {
      XTP_LOG(logDEBUG, _log)
          << TimeStamp() << " Using full matrix method" << flush;

      // get the full matrix
      hstart = std::chrono::system_clock::now();
      Eigen::MatrixXd hfull = h.get_full_matrix();
      hend = std::chrono::system_clock::now();

      elapsed_time = hend - hstart;

      XTP_LOG(logDEBUG, _log) << TimeStamp() << " Full matrix assembled in "
                              << elapsed_time.count() << " secs" << flush;

      // solve theeigenalue problem
      hstart = std::chrono::system_clock::now();
      DS.solve(hfull, _opt.nmax);
      hend = std::chrono::system_clock::now();

      elapsed_time = hend - hstart;
      XTP_LOG(logDEBUG, _log) << TimeStamp() << " Davidson solve done in "
                              << elapsed_time.count() << " secs" << flush;
    }

    energies = DS.eigenvalues();
    coefficients = DS.eigenvectors();

  }

  else {

    XTP_LOG(logDEBUG, _log)
        << TimeStamp() << " Lapack Diagonalization" << flush;

    hstart = std::chrono::system_clock::now();
    Eigen::MatrixXd hfull = h.get_full_matrix();
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;
    XTP_LOG(logDEBUG, _log) << TimeStamp() << " Full matrix assembled in "
                            << elapsed_time.count() << " secs" << flush;

    hstart = std::chrono::system_clock::now();
    tools::linalg_eigenvalues(hfull, energies, coefficients, _opt.nmax);
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;
    XTP_LOG(logDEBUG, _log) << TimeStamp() << " Lapack solve done in "
                            << elapsed_time.count() << " secs" << flush;
  }

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Diagonalization done in "
                          << elapsed_time.count() << " secs" << flush;

  return;
}

void BSE::Solve_singlets_BTDA() {
  SingletOperator_BTDA_ApB Hs_ApB(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hs_ApB);
  Operator_BTDA_AmB Hs_AmB(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Hs_AmB);
  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Setup Full singlet hamiltonian "
                          << flush;
  Solve_antihermitian(Hs_ApB, Hs_AmB, _bse_singlet_energies,
                      _bse_singlet_coefficients, _bse_singlet_coefficients_AR);
}

void BSE::Solve_triplets_BTDA() {
  TripletOperator_BTDA_ApB Ht_ApB(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Ht_ApB);
  Operator_BTDA_AmB Ht_AmB(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(Ht_AmB);
  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Setup Full triplet hamiltonian "
                          << flush;

  Solve_antihermitian(Ht_ApB, Ht_AmB, _bse_triplet_energies,
                      _bse_triplet_coefficients, _bse_triplet_coefficients_AR);
}

template <typename BSE_OPERATOR_ApB, typename BSE_OPERATOR_AmB>
void BSE::Solve_antihermitian(BSE_OPERATOR_ApB& apb, BSE_OPERATOR_AmB& amb,
                              Eigen::VectorXd& energies,
                              Eigen::MatrixXd& coefficients,
                              Eigen::MatrixXd& coefficients_AR) {

  // For details of the method, see EPL,78(2007)12001,
  // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.
  // setup resonant (A) and RARC blocks (B)
  // corresponds to
  // _ApB = (_eh_d +_eh_qp + _eh_d2 + 4.0 * _eh_x);
  // _AmB = (_eh_d +_eh_qp - _eh_d2);

  Eigen::MatrixXd ApB = apb.get_full_matrix();
  Eigen::MatrixXd AmB = amb.get_full_matrix();
  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Setup singlet hamiltonian "
                          << flush;

  // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not
  // positive definite
  //(A-B) is not needed any longer and can be overwritten
  XTP_LOG(logDEBUG, _log) << TimeStamp()
                          << " Trying Cholesky decomposition of KAA-KAB"
                          << flush;
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd> > L(AmB);

  for (int i = 0; i < AmB.rows(); ++i) {
    for (int j = i + 1; j < AmB.cols(); ++j) {
      AmB(i, j) = 0;
    }
  }
  if (L.info() != 0) {
    XTP_LOG(logDEBUG, _log)
        << TimeStamp()
        << " Cholesky decomposition of KAA-KAB was unsucessful. Try a smaller "
           "basisset. This can indicate a triplet instability."
        << flush;
    throw std::runtime_error("Cholesky decompostion failed");
  } else {
    XTP_LOG(logDEBUG, _log)
        << TimeStamp() << " Cholesky decomposition of KAA-KAB was successful"
        << flush;
  }
  Eigen::MatrixXd temp = ApB * AmB;
  ApB.noalias() = AmB.transpose() * temp;
  temp.resize(0, 0);
  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Calculated H = L^T(A+B)L "
                          << flush;
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors;
  XTP_LOG(logDEBUG, _log) << TimeStamp() << " Solving for first " << _opt.nmax
                          << " eigenvectors" << flush;
  bool success_diag =
      tools::linalg_eigenvalues(ApB, eigenvalues, eigenvectors, _opt.nmax);
  if (!success_diag) {
    XTP_LOG(logDEBUG, _log)
        << TimeStamp() << " Could not solve problem" << flush;
  } else {
    XTP_LOG(logDEBUG, _log)
        << TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
  }
  ApB.resize(0, 0);
  if ((eigenvalues.array() < 0).any()) {
    throw std::runtime_error("Negative eigenvalues in BTDA");
  }
  Eigen::VectorXd tempvec =
      eigenvalues.cwiseSqrt();  // has to stay otherwise mkl complains
  energies = tempvec;
  // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 +
  // 1/sqrt(eps_l)L ] R_l
  //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 -
  //                               1/sqrt(eps_l)L ] R_l
  //                               1/sqrt(eps_l)L ] R_l
  // determine inverse of L^T
  Eigen::MatrixXd LmT = AmB.inverse().transpose();
  int dim = LmT.rows();
  coefficients.resize(dim, _opt.nmax);     // resonant part (_X_evec)
  coefficients_AR.resize(dim, _opt.nmax);  // anti-resonant part (_Y_evec)
  for (int level = 0; level < _opt.nmax; level++) {
    double sqrt_eval = std::sqrt(energies(level));
    // get l-th reduced EV
    coefficients.col(level) = (0.5 / sqrt_eval * (energies(level) * LmT + AmB) *
                               eigenvectors.col(level));
    coefficients_AR.col(level) =
        (0.5 / sqrt_eval * (energies(level) * LmT - AmB) *
         eigenvectors.col(level));
  }
  return;
}

void BSE::printFragInfo(const std::vector<QMFragment<BSE_Population> >& frags,
                        int state) const {
  for (const QMFragment<BSE_Population>& frag : frags) {
    double dq = frag.value().H[state] - frag.value().E[state];
    double qeff = dq + frag.value().Gs;
    XTP_LOG(logINFO, _log)
        << format(
               "           Fragment %1$8s%% -- hole: %2$5.1f%%  electron: "
               "%3$5.1f%%  dQ: %4$+5.2f  Qeff: %5$+5.2f") %
               frag.name() % (100.0 * frag.value().H[state]) %
               (100.0 * frag.value().E[state]) % dq % qeff
        << flush;
  }
  return;
}

void BSE::PrintWeight(int i, int i_bse, QMStateType state) {
  const Eigen::MatrixXd& BSECoefs_AR = (state == QMStateType::Singlet)
                                           ? _bse_singlet_coefficients_AR
                                           : _bse_triplet_coefficients_AR;
  const Eigen::MatrixXd& BSECoefs = (state == QMStateType::Singlet)
                                        ? _bse_singlet_coefficients
                                        : _bse_triplet_coefficients;
  double weight = std::pow(BSECoefs(i_bse, i), 2);
  if (!_opt.useTDA) {
    weight -= std::pow(BSECoefs_AR(i_bse, i), 2);
  }
  vc2index vc = vc2index(_opt.vmin, _bse_cmin, _bse_ctotal);
  if (weight > _opt.min_print_weight) {
    XTP_LOG(logINFO, _log)
        << format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") %
               (_opt.homo - vc.v(i_bse)) % (vc.c(i_bse) - _opt.homo - 1) %
               (100.0 * weight)
        << flush;
  }
  return;
}

void BSE::Analyze_singlets(std::vector<QMFragment<BSE_Population> >& singlets) {

  Interaction act;
  QMStateType singlet = QMStateType(QMStateType::Singlet);
  _orbitals.CalcCoupledTransition_Dipoles();
  std::vector<double> oscs = _orbitals.Oscillatorstrengths();
  if (tools::globals::verbose) {
    act = Analyze_eh_interaction(singlet);
  }
  if (singlets.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(singlets, _orbitals, singlet);
  }

  double hrt2ev = tools::conv::hrt2ev;
  XTP_LOG(logINFO, _log) << "  ====== singlet energies (eV) ====== " << flush;
  for (int i = 0; i < _opt.nmax; ++i) {
    double osc = oscs[i];
    if (tools::globals::verbose) {
      XTP_LOG(logINFO, _log)
          << format(
                 "  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> "
                 "= %4$+1.4f <K_x> = %5$+1.4f <K_d> = %6$+1.4f") %
                 (i + 1) % (hrt2ev * _bse_singlet_energies(i)) %
                 (1240.0 / (hrt2ev * _bse_singlet_energies(i))) %
                 (hrt2ev * act.qp_contrib(i)) %
                 (hrt2ev * act.exchange_contrib(i)) %
                 (hrt2ev * act.direct_contrib(i))
          << flush;
    } else {
      XTP_LOG(logINFO, _log)
          << format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm") %
                 (i + 1) % (hrt2ev * _bse_singlet_energies(i)) %
                 (1240.0 / (hrt2ev * _bse_singlet_energies(i)))
          << flush;
    }
    const Eigen::Vector3d& trdip = _orbitals.TransitionDipoles()[i];
    XTP_LOG(logINFO, _log)
        << format(
               "           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = "
               "%2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f") %
               trdip[0] % trdip[1] % trdip[2] % (trdip.squaredNorm()) % osc
        << flush;
    for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
      // if contribution is larger than 0.2, print
      PrintWeight(i, i_bse, singlet);
    }
    // results of fragment population analysis
    if (singlets.size() > 0) {
      printFragInfo(singlets, i);
    }

    XTP_LOG(logINFO, _log) << flush;
  }
  return;
}

void BSE::Analyze_triplets(std::vector<QMFragment<BSE_Population> >& triplets) {

  Interaction act;
  QMStateType triplet = QMStateType(QMStateType::Triplet);
  if (tools::globals::verbose) {
    act = Analyze_eh_interaction(triplet);
  }
  if (triplets.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(triplets, _orbitals, triplet);
  }

  XTP_LOG(logINFO, _log) << "  ====== triplet energies (eV) ====== " << flush;
  for (int i = 0; i < _opt.nmax; ++i) {
    if (tools::globals::verbose) {
      XTP_LOG(logINFO, _log)
          << format(
                 "  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> "
                 "= %4$+1.4f <K_d> = %5$+1.4f") %
                 (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) %
                 (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i))) %
                 (tools::conv::hrt2ev * act.qp_contrib(i)) %
                 (tools::conv::hrt2ev * act.direct_contrib(i))
          << flush;
    } else {
      XTP_LOG(logINFO, _log)
          << format("  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm") %
                 (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) %
                 (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i)))
          << flush;
    }
    for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
      // if contribution is larger than 0.2, print
      PrintWeight(i, i_bse, triplet);
    }
    // results of fragment population analysis
    if (triplets.size() > 0) {
      printFragInfo(triplets, i);
    }
    XTP_LOG(logINFO, _log) << format("   ") << flush;
  }
  // storage to orbitals object

  return;
}

template <typename BSE_OPERATOR>
Eigen::VectorXd BSE::Analyze_IndividualContribution(const QMStateType& type,
                                                    const BSE_OPERATOR& H) {
  Eigen::VectorXd contrib = Eigen::VectorXd::Zero(_opt.nmax);
  const Eigen::MatrixXd& BSECoefs_AR = (type == QMStateType::Singlet)
                                           ? _bse_singlet_coefficients_AR
                                           : _bse_triplet_coefficients_AR;
  const Eigen::MatrixXd& BSECoefs = (type == QMStateType::Singlet)
                                        ? _bse_singlet_coefficients
                                        : _bse_triplet_coefficients;

  for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
    Eigen::MatrixXd slice_R = BSECoefs.block(0, i_exc, _bse_size, 1);
    contrib(i_exc) = (slice_R.transpose() * (H * slice_R)).value();
    if (!_opt.useTDA) {
      Eigen::MatrixXd slice_AR = BSECoefs_AR.block(0, i_exc, _bse_size, 1);
      // get anti-resonant contribution from direct Keh
      contrib(i_exc) -= (slice_AR.transpose() * (H * slice_AR)).value();
    }
  }
  return contrib;
}

BSE::Interaction BSE::Analyze_eh_interaction(const QMStateType& type) {
  Interaction analysis;

  HqpOperator hqp(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(hqp);
  analysis.qp_contrib = Analyze_IndividualContribution(type, hqp);

  HdOperator hd(_epsilon_0_inv, _Mmn, _Hqp);
  configureBSEOperator(hd);
  analysis.direct_contrib = Analyze_IndividualContribution(type, hd);

  if (type == QMStateType::Singlet) {
    HxOperator hx(_epsilon_0_inv, _Mmn, _Hqp);
    configureBSEOperator(hx);
    analysis.exchange_contrib = Analyze_IndividualContribution(type, hx);
  } else {
    analysis.exchange_contrib = Eigen::VectorXd::Zero(0);
  }

  return analysis;
}

}  // namespace xtp
};  // namespace votca
