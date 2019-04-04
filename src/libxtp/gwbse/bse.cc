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

#include <votca/tools/linalg.h>
#include <votca/xtp/bse.h>

#include <votca/xtp/bse_operator.h>
#include <votca/xtp/davidsonsolver.h>

#include "votca/xtp/qmstate.h"
#include "votca/xtp/vc2index.h"

#include <chrono>

using boost::format;
using std::flush;

namespace votca {
namespace xtp {

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

  TripletOperator_TDA Ht(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Ht);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Setup TDA triplet hamiltonian " << flush;
  solve_hermitian(Ht, _bse_triplet_energies, _bse_triplet_coefficients);

  return;
}

void BSE::Solve_singlets() {
  if (_opt.useTDA) {
    Solve_singlets_TDA();
  } else {
    if (_opt.davidson) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp()
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
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp()
          << " Davidson solver not implemented for BTDA. Full diagonalization."
          << flush;
      _opt.davidson = 0;
    }
    Solve_triplets_BTDA();
  }
}

void BSE::Solve_singlets_TDA() {

  SingletOperator_TDA Hs(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Hs);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Setup TDA singlet hamiltonian " << flush;

  solve_hermitian(Hs, _bse_singlet_energies, _bse_singlet_coefficients);
}

void BSE::SetupHs() {

  SingletOperator_TDA Hs(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Hs);
  _eh_s = Hs.get_full_matrix();
}

void BSE::SetupHt() {

  TripletOperator_TDA Ht(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Ht);
  _eh_t = Ht.get_full_matrix();
}

template <typename BSE_OPERATOR>
void BSE::solve_hermitian(BSE_OPERATOR& h, Eigen::VectorXd& energies,
                          Eigen::MatrixXd& coefficients) {

  std::chrono::time_point<std::chrono::system_clock> start, end;
  std::chrono::time_point<std::chrono::system_clock> hstart, hend;
  std::chrono::duration<double> elapsed_time;
  start = std::chrono::system_clock::now();

  CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Solving for first "
                               << _opt.nmax << " eigenvectors" << flush;

  if (_opt.davidson) {

    DavidsonSolver DS(_log);
    DS.set_correction(_opt.davidson_correction);
    DS.set_tolerance(_opt.davidson_tolerance);
    DS.set_max_search_space(10 * _opt.nmax);

    if (_opt.matrixfree) {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Using matrix free method" << flush;
      DS.solve(h, _opt.nmax);
    }

    else {
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Using full matrix method" << flush;

      // get the full matrix
      hstart = std::chrono::system_clock::now();
      Eigen::MatrixXd hfull = h.get_full_matrix();
      hend = std::chrono::system_clock::now();

      elapsed_time = hend - hstart;

      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Full matrix assembled in "
          << elapsed_time.count() << " secs" << flush;

      // solve theeigenalue problem
      hstart = std::chrono::system_clock::now();
      DS.solve(hfull, _opt.nmax);
      hend = std::chrono::system_clock::now();

      elapsed_time = hend - hstart;
      CTP_LOG(ctp::logDEBUG, _log)
          << ctp::TimeStamp() << " Davidson solve done in "
          << elapsed_time.count() << " secs" << flush;
    }

    energies = DS.eigenvalues();
    coefficients = DS.eigenvectors();

  }

  else {

    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Lapack Diagonalization" << flush;

    hstart = std::chrono::system_clock::now();
    Eigen::MatrixXd hfull = h.get_full_matrix();
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Full matrix assembled in "
        << elapsed_time.count() << " secs" << flush;

    hstart = std::chrono::system_clock::now();
    tools::linalg_eigenvalues(hfull, energies, coefficients, _opt.nmax);
    hend = std::chrono::system_clock::now();

    elapsed_time = hend - hstart;
    CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Lapack solve done in "
                                 << elapsed_time.count() << " secs" << flush;
  }

  end = std::chrono::system_clock::now();
  elapsed_time = end - start;

  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Diagonalization done in " << elapsed_time.count()
      << " secs" << flush;

  return;
}

void BSE::Solve_singlets_BTDA() {
  SingletOperator_BTDA_ApB Hs_ApB(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Hs_ApB);
  Operator_BTDA_AmB Hs_AmB(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Hs_AmB);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Setup Full singlet hamiltonian " << flush;
  Solve_antihermitian(Hs_ApB, Hs_AmB, _bse_singlet_energies,
                      _bse_singlet_coefficients, _bse_singlet_coefficients_AR);
}

void BSE::Solve_triplets_BTDA() {
  TripletOperator_BTDA_ApB Ht_ApB(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Ht_ApB);
  Operator_BTDA_AmB Ht_AmB(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(Ht_AmB);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Setup Full triplet hamiltonian " << flush;

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

  // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not
  // positive definite
  //(A-B) is not needed any longer and can be overwritten
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Trying Cholesky decomposition of KAA-KAB"
      << flush;
  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd> > L(AmB);

  for (int i = 0; i < AmB.rows(); ++i) {
    for (int j = i + 1; j < AmB.cols(); ++j) {
      AmB(i, j) = 0;
    }
  }
  if (L.info() != 0) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp()
        << " Cholesky decomposition of KAA-KAB was unsucessful. Try a smaller "
           "basisset. This can indicate a triplet instability."
        << flush;
    throw std::runtime_error("Cholesky decompostion failed");
  } else {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp()
        << " Cholesky decomposition of KAA-KAB was successful" << flush;
  }

  Eigen::MatrixXd temp = ApB * AmB;
  ApB.noalias() = AmB.transpose() * temp;
  temp.resize(0, 0);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Calculated H = L^T(A+B)L " << flush;
  Eigen::VectorXd eigenvalues;
  Eigen::MatrixXd eigenvectors;
  CTP_LOG(ctp::logDEBUG, _log) << ctp::TimeStamp() << " Solving for first "
                               << _opt.nmax << " eigenvectors" << flush;
  bool success_diag =
      tools::linalg_eigenvalues(ApB, eigenvalues, eigenvectors, _opt.nmax);
  if (!success_diag) {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Could not solve problem" << flush;
  } else {
    CTP_LOG(ctp::logDEBUG, _log)
        << ctp::TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
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

void BSE::printFragInfo(const Population& pop, int i) {
  CTP_LOG(ctp::logINFO, _log)
      << format(
             "           Fragment A -- hole: %1$5.1f%%  electron: %2$5.1f%%  "
             "dQ: %3$+5.2f  Qeff: %4$+5.2f") %
             (100.0 * pop.popH[i](0)) % (100.0 * pop.popE[i](0)) %
             (pop.Crgs[i](0)) % (pop.Crgs[i](0) + pop.popGs(0))
      << flush;
  CTP_LOG(ctp::logINFO, _log)
      << format(
             "           Fragment B -- hole: %1$5.1f%%  electron: %2$5.1f%%  "
             "dQ: %3$+5.2f  Qeff: %4$+5.2f") %
             (100.0 * pop.popH[i](1)) % (100.0 * pop.popE[i](1)) %
             (pop.Crgs[i](1)) % (pop.Crgs[i](1) + pop.popGs(1))
      << flush;
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
    CTP_LOG(ctp::logINFO, _log)
        << format("           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") %
               (_opt.homo - vc.v(i_bse)) % (vc.c(i_bse) - _opt.homo - 1) %
               (100.0 * weight)
        << flush;
  }
  return;
}

void BSE::Analyze_singlets(const AOBasis& dftbasis) {

  Interaction act;
  Population pop;
  QMStateType singlet = QMStateType(QMStateType::Singlet);
  std::vector<tools::vec> transition_dipoles =
      CalcCoupledTransition_Dipoles(dftbasis);
  _orbitals.TransitionDipoles() = transition_dipoles;
  std::vector<double> oscs = _orbitals.Oscillatorstrengths();

  if (tools::globals::verbose) {
    act = Analyze_eh_interaction(singlet);
  }
  if (dftbasis.getAOBasisFragA() > 0 && dftbasis.getAOBasisFragB() > 0) {
    pop = FragmentPopulations(singlet, dftbasis);
    _orbitals.setFragmentChargesSingEXC(pop.Crgs);
    _orbitals.setFragment_E_localisation_singlet(pop.popE);
    _orbitals.setFragment_H_localisation_singlet(pop.popH);
    _orbitals.setFragmentChargesGS(pop.popGs);
  }

  double hrt2ev = tools::conv::hrt2ev;
  CTP_LOG(ctp::logINFO, _log)
      << "  ====== singlet energies (eV) ====== " << flush;
  for (int i = 0; i < _opt.nmax; ++i) {
    const tools::vec& trdip = transition_dipoles[i];
    double osc = oscs[i];
    if (tools::globals::verbose) {
      CTP_LOG(ctp::logINFO, _log)
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
      CTP_LOG(ctp::logINFO, _log)
          << format("  S = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm") %
                 (i + 1) % (hrt2ev * _bse_singlet_energies(i)) %
                 (1240.0 / (hrt2ev * _bse_singlet_energies(i)))
          << flush;
    }
    CTP_LOG(ctp::logINFO, _log)
        << format(
               "           TrDipole length gauge[e*bohr]  dx = %1$+1.4f dy = "
               "%2$+1.4f dz = %3$+1.4f |d|^2 = %4$+1.4f f = %5$+1.4f") %
               trdip.getX() % trdip.getY() % trdip.getZ() % (trdip * trdip) %
               osc
        << flush;
    for (int i_bse = 0; i_bse < _bse_size; ++i_bse) {
      // if contribution is larger than 0.2, print
      PrintWeight(i, i_bse, singlet);
    }
    // results of fragment population analysis
    if (dftbasis.getAOBasisFragA() > 0 && dftbasis.getAOBasisFragB() > 0) {
      printFragInfo(pop, i);
    }

    CTP_LOG(ctp::logINFO, _log) << flush;
  }
  return;
}

void BSE::Analyze_triplets(const AOBasis& dftbasis) {

  Interaction act;
  Population pop;
  QMStateType triplet = QMStateType(QMStateType::Triplet);
  if (tools::globals::verbose) {
    act = Analyze_eh_interaction(triplet);
  }
  if (dftbasis.getAOBasisFragA() > 0) {
    pop = FragmentPopulations(triplet, dftbasis);
    _orbitals.setFragmentChargesTripEXC(pop.Crgs);
    _orbitals.setFragment_E_localisation_triplet(pop.popE);
    _orbitals.setFragment_H_localisation_triplet(pop.popH);
    _orbitals.setFragmentChargesGS(pop.popGs);
  }
  CTP_LOG(ctp::logINFO, _log)
      << "  ====== triplet energies (eV) ====== " << flush;
  for (int i = 0; i < _opt.nmax; ++i) {
    if (tools::globals::verbose) {
      CTP_LOG(ctp::logINFO, _log)
          << format(
                 "  T = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f nm <FT> "
                 "= %4$+1.4f <K_d> = %5$+1.4f") %
                 (i + 1) % (tools::conv::hrt2ev * _bse_triplet_energies(i)) %
                 (1240.0 / (tools::conv::hrt2ev * _bse_triplet_energies(i))) %
                 (tools::conv::hrt2ev * act.qp_contrib(i)) %
                 (tools::conv::hrt2ev * act.direct_contrib(i))
          << flush;
    } else {
      CTP_LOG(ctp::logINFO, _log)
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
    if (dftbasis.getAOBasisFragA() > 0) {
      printFragInfo(pop, i);
    }
    CTP_LOG(ctp::logINFO, _log) << format("   ") << flush;
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

  HqpOperator hqp(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(hqp);
  analysis.qp_contrib = Analyze_IndividualContribution(type, hqp);

  HdOperator hd(_epsilon_0_inv, _log, _Mmn, _Hqp);
  configureBSEOperator(hd);
  analysis.direct_contrib = Analyze_IndividualContribution(type, hd);

  if (type == QMStateType::Singlet) {
    HxOperator hx(_epsilon_0_inv, _log, _Mmn, _Hqp);
    configureBSEOperator(hx);
    analysis.exchange_contrib = Analyze_IndividualContribution(type, hx);
  } else {
    analysis.exchange_contrib = Eigen::VectorXd::Zero(0);
  }

  return analysis;
}

BSE::Population BSE::FragmentPopulations(const QMStateType& type,
                                         const AOBasis& dftbasis) {
  Population pop;
  // Mulliken fragment population analysis
  AOOverlap dftoverlap;
  dftoverlap.Fill(dftbasis);
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Filled DFT Overlap matrix of dimension: "
      << dftoverlap.Matrix().rows() << flush;
  // ground state populations
  Eigen::MatrixXd DMAT = _orbitals.DensityMatrixGroundState();
  Eigen::VectorXd nuccharges =
      _orbitals.FragmentNuclearCharges(dftbasis.getAOBasisFragA());
  Eigen::VectorXd pops = _orbitals.LoewdinPopulation(
      DMAT, dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
  pop.popGs = nuccharges - pops;
  // population to electron charges and add nuclear charges
  for (int i_state = 0; i_state < _opt.nmax; i_state++) {
    QMState state = QMState(type, i_state, false);
    // checking Density Matrices
    std::vector<Eigen::MatrixXd> DMAT =
        _orbitals.DensityMatrixExcitedState(state);
    // hole part
    Eigen::VectorXd popsH = _orbitals.LoewdinPopulation(
        DMAT[0], dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
    pop.popH.push_back(popsH);
    // electron part
    Eigen::VectorXd popsE = _orbitals.LoewdinPopulation(
        DMAT[1], dftoverlap.Matrix(), dftbasis.getAOBasisFragA());
    pop.popE.push_back(popsE);
    // update effective charges
    Eigen::VectorXd diff = popsH - popsE;
    pop.Crgs.push_back(diff);
  }
  CTP_LOG(ctp::logDEBUG, _log)
      << ctp::TimeStamp() << " Ran Excitation fragment population analysis "
      << flush;

  return pop;
}

std::vector<Eigen::MatrixXd> BSE::CalcFreeTransition_Dipoles(
    const AOBasis& dftbasis) {
  const Eigen::MatrixXd& dft_orbitals = _orbitals.MOCoefficients();
  // Testing electric dipole AOMatrix
  AODipole dft_dipole;
  dft_dipole.Fill(dftbasis);

  // now transition dipole elements for free interlevel transitions
  std::vector<Eigen::MatrixXd> interlevel_dipoles;

  Eigen::MatrixXd empty =
      dft_orbitals.block(0, _bse_cmin, dftbasis.AOBasisSize(), _bse_ctotal);
  Eigen::MatrixXd occ =
      dft_orbitals.block(0, _opt.vmin, dftbasis.AOBasisSize(), _bse_vtotal);
  for (int i_comp = 0; i_comp < 3; i_comp++) {
    interlevel_dipoles.push_back(occ.transpose() * dft_dipole.Matrix()[i_comp] *
                                 empty);
  }
  return interlevel_dipoles;
}

std::vector<tools::vec> BSE::CalcCoupledTransition_Dipoles(
    const AOBasis& dftbasis) {
  std::vector<Eigen::MatrixXd> interlevel_dipoles =
      CalcFreeTransition_Dipoles(dftbasis);
  vc2index vc = vc2index(0, 0, _bse_ctotal);
  std::vector<tools::vec> dipols;
  const double sqrt2 = sqrt(2.0);
  for (int i_exc = 0; i_exc < _opt.nmax; i_exc++) {
    tools::vec tdipole = tools::vec(0, 0, 0);
    for (int c = 0; c < _bse_ctotal; c++) {
      for (int v = 0; v < _bse_vtotal; v++) {
        int index_vc = vc.I(v, c);
        double factor = _bse_singlet_coefficients(index_vc, i_exc);
        if (_bse_singlet_coefficients_AR.rows() > 0) {
          factor += _bse_singlet_coefficients_AR(index_vc, i_exc);
        }
        // The Transition dipole is sqrt2 bigger because of the spin, the
        // excited state is a linear combination of 2 slater determinants, where
        // either alpha or beta spin electron is excited
        tdipole.x() += factor * interlevel_dipoles[0](v, c);
        tdipole.y() += factor * interlevel_dipoles[1](v, c);
        tdipole.z() += factor * interlevel_dipoles[2](v, c);
      }
    }

    dipols.push_back(-sqrt2 * tdipole);  //- because electrons are negative
  }
  return dipols;
}

}  // namespace xtp
};  // namespace votca
