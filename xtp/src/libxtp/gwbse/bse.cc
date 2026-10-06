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
#include "votca/xtp/bse_fullsolver.h"
#include "votca/xtp/bse_initialization.h"
#include "votca/xtp/bse_operator.h"
#include "votca/xtp/bseoperator_btda.h"
#include "votca/xtp/davidsonsolver.h"
#include "votca/xtp/environmentscreening.h"
#include "votca/xtp/populationanalysis.h"
#include "votca/xtp/qmfragment.h"
#include "votca/xtp/rpa.h"
#include "votca/xtp/vc2index.h"

using std::flush;

namespace votca {
namespace xtp {

namespace {
double DavidsonToleranceValue(const std::string& tol) {
  if (tol == "loose") {
    return 1e-3;
  } else if (tol == "normal") {
    return 1e-4;
  } else if (tol == "strict") {
    return 1e-5;
  } else if (tol == "lapack") {
    return 1e-9;
  }
  throw std::runtime_error(tol + " is not a valid Davidson tolerance");
}
}  // namespace

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

void BSE::configure_with_precomputed_screening(
    const options& opt, const Eigen::VectorXd& RPAInputEnergies,
    const Eigen::MatrixXd& Hqp_in, const Eigen::VectorXd& epsilon_0_inv) {
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
  epsilon_0_inv_ = epsilon_0_inv;
}

tools::EigenSystem& BSE::GetBSEEigenSystem(const QMStateType& type,
                                           Orbitals& orb) const {
  if (type == QMStateType::Singlet) {
    return orb.BSESinglets();
  } else if (type == QMStateType::Triplet) {
    return orb.BSETriplets();
  } else {
    throw std::runtime_error(
        "Unsupported QMStateType in BSE::GetBSEEigenSystem");
  }
}

const tools::EigenSystem& BSE::GetBSEEigenSystem(const QMStateType& type,
                                                 const Orbitals& orb) const {
  if (type == QMStateType::Singlet) {
    return orb.BSESinglets();
  } else if (type == QMStateType::Triplet) {
    return orb.BSETriplets();
  } else {
    throw std::runtime_error(
        "Unsupported QMStateType in BSE::GetBSEEigenSystem");
  }
}

std::string BSE::StateEnergiesHeader(const QMStateType& type) const {
  if (type == QMStateType::Singlet) {
    return "  ====== singlet energies (eV) ====== ";
  } else if (type == QMStateType::Triplet) {
    return "  ====== triplet energies (eV) ====== ";
  } else {
    throw std::runtime_error(
        "Unsupported QMStateType in BSE::StateEnergiesHeader");
  }
}

std::string BSE::StateShortLabel(const QMStateType& type) const {
  if (type == QMStateType::Singlet) {
    return "S";
  } else if (type == QMStateType::Triplet) {
    return "T";
  } else {
    throw std::runtime_error("Unsupported QMStateType in BSE::StateShortLabel");
  }
}

std::string BSE::StateDynamicLabel(const QMStateType& type) const {
  if (type == QMStateType::Singlet) {
    return "singlet";
  } else if (type == QMStateType::Triplet) {
    return "triplet";
  } else {
    throw std::runtime_error(
        "Unsupported QMStateType in BSE::StateDynamicLabel");
  }
}

double BSE::ExchangePrefactor(const QMStateType& type) const {
  if (type == QMStateType::Singlet) {
    return 2.0;
  } else if (type == QMStateType::Triplet) {
    return 0.0;
  } else {
    throw std::runtime_error(
        "Unsupported QMStateType in BSE::ExchangePrefactor");
  }
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

void BSE::setReactionField(const Eigen::MatrixXd& R, bool include_kreac) {
  reaction_field_ = R;
  include_kreac_ = include_kreac;
  dressing_ = (R.size() > 0) ? EnvironmentScreening::DressingMatrix(R)
                             : Eigen::MatrixXd();
}

// The BSE operators read slices and rows up to the highest conduction level;
// the dielectric matrix reads all rows of the occupied slices. Only those
// are rotated (see TCMatrix_gwbse::MultiplyRightWithAuxMatrixLeading).
void BSE::RotateIntoScreeningFrame(const Eigen::MatrixXd& U) {
  const Index occupied_slices = opt_.homo + 1 - opt_.rpamin;
  const Index lead = opt_.cmax + 1 - opt_.rpamin;
  Mmn_.MultiplyRightWithAuxMatrixLeading(U, occupied_slices, lead);
}

void BSE::MatchDressingToRoute() {
  const bool environment = reaction_field_.size() > 0;
  const bool dressed_route = environment && include_kreac_;
  // Dressed: every pairing of two M's is through u = v + v_reac, so eps
  // below is the environment-screened one and Kx carries K_reac. Bare:
  // Kx stays v, and W_tot is assembled explicitly below. GW may have left
  // the integrals either way.
  if (dressed_route && !Mmn_.Dressed()) {
    Mmn_.DressAuxIndex(dressing_);
  } else if (!dressed_route && Mmn_.Dressed()) {
    Mmn_.UndressAuxIndex();
  }
}

// W = vectors diag(values) vectors^T in the current frame of the
// integrals, as the operators use it: the inverse of epsilon without its
// non-positive modes, or W_tot of the bare environment route.
SymmetricEigenSystem BSE::ScreenedInteraction(
    const Eigen::VectorXd& RPAInputEnergies, double energy,
    DFTTimings& timings) const {
  RPA rpa = RPA(log_, Mmn_);
  rpa.configure(opt_.homo, opt_.rpamin, opt_.rpamax);
  rpa.setRPAInputEnergies(RPAInputEnergies);
  rpa.setTimings(&timings);
  const Eigen::MatrixXd screening = rpa.calculate_epsilon_r(energy);

  auto t = timings.Measure("eigensolver");
  if (reaction_field_.size() > 0 && !include_kreac_) {
    // Operators expect M U with U orthogonal and diag(epsilon_0_inv_) the
    // screened interaction in that frame: Kd = M U d U^T M^T and
    // Kx = M U U^T M^T = M M^T = v. So diagonalize W_tot, not eps, and
    // hand its eigenvalues over as epsilon_0_inv_.
    const Index n = screening.rows();
    const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);
    // 1 + R and u^-1 + eps - 1 are positive definite (Cholesky, LU if not)
    const Eigen::MatrixXd u_inv =
        InverseSPD(I + Mmn_.ToCurrentAuxFrame(reaction_field_));
    Eigen::MatrixXd W = InverseSPD(u_inv + screening - I);
    W = 0.5 * (W + W.transpose()).eval();
    return SymmetricEigen(W);
  }

  SymmetricEigenSystem es = SymmetricEigen(screening);
  for (Index i = 0; i < es.values.size(); ++i) {
    es.values(i) = (es.values(i) > 1e-8) ? 1 / es.values(i) : 0.0;
  }
  return es;
}

void BSE::SetupDirectInteractionOperator(
    const Eigen::VectorXd& RPAInputEnergies, double energy) {
  MatchDressingToRoute();
  DFTTimings timings;
  SymmetricEigenSystem W =
      ScreenedInteraction(RPAInputEnergies, energy, timings);
  {
    auto t = timings.Measure("rotation");
    RotateIntoScreeningFrame(W.vectors);
  }
  epsilon_0_inv_ = std::move(W.values);
  XTP_LOG(Log::error, log_) << TimeStamp() << " BSE screening at " << energy
                            << " Ha: " << timings.Format() << flush;
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
  H.set_direct_cache_limit(opt_.direct_cache_gb * 1e9);
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

  if (reaction_field_.size() > 0) {
    const Eigen::VectorXd kreac =
        ReactionFieldExchange(orb.BSESinglets(), opt_.useTDA);
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Environment K_reac (linear-response part), "
        << (include_kreac_ ? "included" : "NOT included, estimate only")
        << ", first order [eV]:" << flush;
    for (Index i = 0; i < kreac.size(); ++i) {
      XTP_LOG(Log::error, log_)
          << (boost::format("  S = %1$4d Omega = %2$+1.4f eV  <K_reac> = "
                            "%3$+1.4f eV") %
              (i + 1) %
              (orb.BSESinglets().eigenvalues()(i) * tools::conv::hrt2ev) %
              (kreac(i) * tools::conv::hrt2ev))
                 .str()
          << flush;
    }
  }
}

Eigen::VectorXd BSE::ReactionFieldExchange(const tools::EigenSystem& es,
                                           bool tda) const {
  // The transition density of each state in the auxiliary basis,
  // d = sum_vc (X+Y)_vc M_v,c, and <K_reac> = 2 d^T R d: the singlet
  // spin factor of Kx, which K_reac shares, times the first-order change
  // (X+Y)^T dK (X+Y) that a perturbation entering A and B alike produces.
  const Eigen::MatrixXd Rc = Mmn_.ToCurrentAuxFrame(reaction_field_);
  const Index vmin = opt_.vmin - opt_.rpamin;
  const Index cmin = bse_cmin_ - opt_.rpamin;
  const Index nstates = es.eigenvalues().size();
  Eigen::VectorXd result = Eigen::VectorXd::Zero(nstates);
#pragma omp parallel for schedule(dynamic)
  for (Index s = 0; s < nstates; ++s) {
    Eigen::VectorXd coeffs = es.eigenvectors().col(s);
    if (!tda) {
      coeffs += es.eigenvectors2().col(s);
    }
    const Eigen::Map<const Eigen::MatrixXd> mat(coeffs.data(), bse_ctotal_,
                                                bse_vtotal_);
    Eigen::VectorXd d = Eigen::VectorXd::Zero(Mmn_.auxsize());
    for (Index v = 0; v < bse_vtotal_; ++v) {
      d +=
          Mmn_[v + vmin].middleRows(cmin, bse_ctotal_).transpose() * mat.col(v);
    }
    result(s) = 2.0 * d.dot(Rc * d);
  }
  return result;
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
  // Hermitian (A-B)(A+B) form: one application of A and of B per new
  // vector. Falls back to the general non-Hermitian Davidson if A-B or A+B
  // is not positive definite.
  try {
    FullBSEDavidson::Options fopt;
    fopt.tolerance = DavidsonToleranceValue(opt_.davidson_tolerance);
    fopt.max_iterations = opt_.davidson_maxiter;
    fopt.max_subspace = 10 * opt_.nmax;
    FullBSEDavidson solver(log_, fopt);
    return solver.Solve(Aop, Bop, opt_.nmax);
  } catch (const std::runtime_error& error) {
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " " << error.what()
        << "; using the non-Hermitian Davidson solver" << flush;
  }

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
  Eigen::MatrixXd initial_guess = BuildFullBSEXRankedInitialGuess(
      Aop.diagonal(), Bop.diagonal(), opt_.nmax);

  DS.solve(Hop, opt_.nmax, initial_guess);

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

void BSE::printFragInfo(const std::vector<QMFragment<BSE_Population>>& frags,
                        Index state) const {
  for (const QMFragment<BSE_Population>& frag : frags) {
    double dq = frag.value().H[state] + frag.value().E[state];
    double qeff = dq + frag.value().Gs;
    XTP_LOG(Log::error, log_)
        << boost::format(
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
          << boost::format(
                 "           HOMO-%1$-3d -> LUMO+%2$-3d  : %3$3.1f%%") %
                 (opt_.homo - vc.v(i_bse)) % (vc.c(i_bse) - opt_.homo - 1) %
                 (100.0 * weight)
          << flush;
    }
  }
  return;
}

void BSE::Analyze_singlets(std::vector<QMFragment<BSE_Population>> fragments,
                           const Orbitals& orb) const {

  QMStateType type = QMStateType(QMStateType::Singlet);

  Eigen::VectorXd oscs = orb.Oscillatorstrengths();
  Interaction act = Analyze_eh_interaction(type, orb);

  if (fragments.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(fragments, orb, type);
  }

  const tools::EigenSystem& bse = GetBSEEigenSystem(type, orb);
  const Eigen::VectorXd& energies = bse.eigenvalues();

  double hrt2ev = tools::conv::hrt2ev;
  XTP_LOG(Log::error, log_) << StateEnergiesHeader(type) << flush;
  for (Index i = 0; i < opt_.nmax; ++i) {
    Eigen::VectorXd weights = bse.eigenvectors().col(i).cwiseAbs2();
    if (!orb.getTDAApprox()) {
      weights -= bse.eigenvectors2().col(i).cwiseAbs2();
    }

    double osc = oscs[i];
    // On the dressed route (environment with K_reac) the exchange operator
    // acts on integrals dressed with (1+R)^1/2, so its expectation value is
    // that of K_x + K_reac. Label it as such.
    const bool dressed = reaction_field_.size() > 0 && include_kreac_;
    const std::string kx_label = dressed ? "<K_x+K_reac>" : "<K_x>";
    XTP_LOG(Log::error, log_)
        << boost::format(
               "  %1$2s = %2$4d Omega = %3$+1.12f eV  lamdba = %4$+3.2f nm "
               "<FT> = %5$+1.4f %8$s = %6$+1.4f <K_d> = %7$+1.4f") %
               StateShortLabel(type) % (i + 1) % (hrt2ev * energies(i)) %
               (1240.0 / (hrt2ev * energies(i))) %
               (hrt2ev * act.qp_contrib(i)) %
               (hrt2ev * act.exchange_contrib(i)) %
               (hrt2ev * act.direct_contrib(i)) % kx_label
        << flush;

    const Eigen::Vector3d& trdip = orb.TransitionDipoles()[i];
    XTP_LOG(Log::error, log_)
        << boost::format(
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

void BSE::Analyze_triplets(std::vector<QMFragment<BSE_Population>> fragments,
                           const Orbitals& orb) const {

  QMStateType type = QMStateType(QMStateType::Triplet);
  Interaction act = Analyze_eh_interaction(type, orb);

  if (fragments.size() > 0) {
    Lowdin low;
    low.CalcChargeperFragment(fragments, orb, type);
  }

  const tools::EigenSystem& bse = GetBSEEigenSystem(type, orb);
  const Eigen::VectorXd& energies = bse.eigenvalues();

  XTP_LOG(Log::error, log_) << StateEnergiesHeader(type) << flush;
  for (Index i = 0; i < opt_.nmax; ++i) {
    Eigen::VectorXd weights = bse.eigenvectors().col(i).cwiseAbs2();
    if (!orb.getTDAApprox()) {
      weights -= bse.eigenvectors2().col(i).cwiseAbs2();
    }

    XTP_LOG(Log::error, log_)
        << boost::format(
               "  %1$2s = %2$4d Omega = %3$+1.12f eV  lamdba = %4$+3.2f nm "
               "<FT> = %5$+1.4f <K_d> = %6$+1.4f") %
               StateShortLabel(type) % (i + 1) %
               (tools::conv::hrt2ev * energies(i)) %
               (1240.0 / (tools::conv::hrt2ev * energies(i))) %
               (tools::conv::hrt2ev * act.qp_contrib(i)) %
               (tools::conv::hrt2ev * act.direct_contrib(i))
        << flush;

    PrintWeights(weights);
    if (fragments.size() > 0) {
      printFragInfo(fragments, i);
    }
    XTP_LOG(Log::error, log_) << boost::format("   ") << flush;
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

  const tools::EigenSystem& BSECoefs = GetBSEEigenSystem(type, orb);

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

  double xpref = ExchangePrefactor(type);
  if (xpref != 0.0) {
    HxOperator hx(epsilon_0_inv_, Mmn_, Hqp_);
    configureBSEOperator(hx);
    ExpectationValues expectation_values =
        ExpectationValue_Operator(type, orb, hx);
    analysis.exchange_contrib = xpref * expectation_values.direct_term;
    if (!orb.getTDAApprox()) {
      analysis.exchange_contrib += xpref * expectation_values.cross_term;
    }
  } else {
    analysis.exchange_contrib =
        Eigen::VectorXd::Zero(analysis.direct_contrib.size());
  }

  return analysis;
}

double BSE::state_matrix_bytes_ = 1e9;

namespace {
using RowMajorMatrix =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

// Number of hole levels per pass, so that the stacked intermediates of one
// pass stay within the given bytes.
Index HoleChunk(Index vtotal, Index ctotal, Index auxsize, double bytes) {
  const double per_level =
      8.0 * double(auxsize) * double(ctotal + 2 * std::max(vtotal, ctotal));
  return std::clamp<Index>(Index(bytes / per_level), 1, vtotal);
}

// sum_PQ W_PQ G_PQ for W = vectors diag(values) vectors^T
double ContractScreened(const SymmetricEigenSystem& W,
                        const Eigen::MatrixXd& G) {
  const Eigen::MatrixXd GPhi = G * W.vectors;
  return (W.vectors.cwiseProduct(GPhi).colwise().sum().transpose().array() *
          W.values.array())
      .sum();
}
}  // namespace

// G with <A|Hd|B> = -sum_PQ W_PQ G_PQ, summed over the pairs (A, B), for the
// screened direct term Hd(v1c1, v2c2) = -M[c1](c2,:) W M[v1](v2,:)^T:
//   G = sum_{v1 v2} b_{v1v2} a_{v1v2}^T,  a_{v1v2} = M[v1](v2,:),
//   b_{v1v2} = sum_{c1 c2} A(v1,c1) B(v2,c2) M[c1](c2,:).
Eigen::MatrixXd BSE::DirectTermDensity(
    const std::vector<std::pair<const Eigen::VectorXd*,
                                const Eigen::VectorXd*>>& pairs) const {
  const Index vmin = opt_.vmin - opt_.rpamin;
  const Index cmin = bse_cmin_ - opt_.rpamin;
  const Index vt = bse_vtotal_;
  const Index ct = bse_ctotal_;
  const Index N = Mmn_.auxsize();
  const Index chunk = HoleChunk(vt, ct, N, state_matrix_bytes_);
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N, N);
  for (Index v0 = 0; v0 < vt; v0 += chunk) {
    const Index nc = std::min(chunk, vt - v0);
    // rows j * vt + v1: b_{v1, v0 + j}
    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(nc * vt, N);
    for (const auto& pair : pairs) {
      // coefficient vectors as (c, v) matrices, index v * ct + c
      Eigen::Map<const Eigen::MatrixXd> A(pair.first->data(), ct, vt);
      Eigen::Map<const Eigen::MatrixXd> B(pair.second->data(), ct, vt);
      const Eigen::MatrixXd Bchunk = B.middleCols(v0, nc).transpose();
      // T[j](c1, :) = sum_c2 B(v0 + j, c2) M[c1](c2, :)
      std::vector<RowMajorMatrix> T(std::size_t(nc), RowMajorMatrix(ct, N));
#pragma omp parallel for schedule(dynamic)
      for (Index c1 = 0; c1 < ct; ++c1) {
        const Eigen::MatrixXd t = Bchunk * Mmn_[c1 + cmin].middleRows(cmin, ct);
        for (Index j = 0; j < nc; ++j) {
          T[std::size_t(j)].row(c1) = t.row(j);
        }
      }
#pragma omp parallel for schedule(dynamic)
      for (Index j = 0; j < nc; ++j) {
        b.middleRows(j * vt, vt).noalias() += A.transpose() * T[std::size_t(j)];
      }
    }
    Eigen::MatrixXd a(nc * vt, N);
#pragma omp parallel for schedule(dynamic)
    for (Index v1 = 0; v1 < vt; ++v1) {
      for (Index j = 0; j < nc; ++j) {
        a.row(j * vt + v1) = Mmn_[v1 + vmin].row(vmin + v0 + j);
      }
    }
    G.noalias() += b.transpose() * a;
  }
  return G;
}

// G with <A|Hd2|B> = -sum_PQ W_PQ G_PQ for the screened coupling term
// Hd2(v1c1, v2c2) = -M[v1](c2,:) W M[c1](v2,:)^T:
//   G = sum_{c1 v2} b_{c1v2} a_{c1v2}^T,  a_{c1v2} = M[c1](v2,:),
//   b_{c1v2} = sum_{v1 c2} A(v1,c1) B(v2,c2) M[v1](c2,:).
Eigen::MatrixXd BSE::CouplingTermDensity(const Eigen::VectorXd& Avec,
                                         const Eigen::VectorXd& Bvec) const {
  const Index vmin = opt_.vmin - opt_.rpamin;
  const Index cmin = bse_cmin_ - opt_.rpamin;
  const Index vt = bse_vtotal_;
  const Index ct = bse_ctotal_;
  const Index N = Mmn_.auxsize();
  const Index chunk = HoleChunk(vt, ct, N, state_matrix_bytes_);
  Eigen::Map<const Eigen::MatrixXd> A(Avec.data(), ct, vt);
  Eigen::Map<const Eigen::MatrixXd> B(Bvec.data(), ct, vt);
  Eigen::MatrixXd G = Eigen::MatrixXd::Zero(N, N);
  for (Index v0 = 0; v0 < vt; v0 += chunk) {
    const Index nc = std::min(chunk, vt - v0);
    const Eigen::MatrixXd Bchunk = B.middleCols(v0, nc).transpose();
    // U[j](v1, :) = sum_c2 B(v0 + j, c2) M[v1](c2, :)
    std::vector<RowMajorMatrix> U(std::size_t(nc), RowMajorMatrix(vt, N));
#pragma omp parallel for schedule(dynamic)
    for (Index v1 = 0; v1 < vt; ++v1) {
      const Eigen::MatrixXd u = Bchunk * Mmn_[v1 + vmin].middleRows(cmin, ct);
      for (Index j = 0; j < nc; ++j) {
        U[std::size_t(j)].row(v1) = u.row(j);
      }
    }
    // rows j * ct + c1: b_{c1, v0 + j} and a_{c1, v0 + j}
    Eigen::MatrixXd b(nc * ct, N);
#pragma omp parallel for schedule(dynamic)
    for (Index j = 0; j < nc; ++j) {
      b.middleRows(j * ct, ct).noalias() = A * U[std::size_t(j)];
    }
    Eigen::MatrixXd a(nc * ct, N);
#pragma omp parallel for schedule(dynamic)
    for (Index c1 = 0; c1 < ct; ++c1) {
      for (Index j = 0; j < nc; ++j) {
        a.row(j * ct + c1) = Mmn_[c1 + cmin].row(vmin + v0 + j);
      }
    }
    G.noalias() += b.transpose() * a;
  }
  return G;
}

// Dynamical Screening in BSE as perturbation to static excitation energies
// as in Phys. Rev. B 80, 241405 (2009) for the TDA case
//
// The integrals stay in the frame of the static screening set up in
// configure, where W(0) = diag(epsilon_0_inv_). For each state, the
// screened terms <X|Hd|X> + <Y|Hd|Y> + 2<Y|Hd2|X> are -sum_PQ W_PQ G_PQ with
// state matrices G built once, so W(w) at the next frequency costs one
// dielectric matrix and its eigenvectors: no pass over the BSE Hamiltonian
// and no rotation of the integrals per state and step.
void BSE::Perturbative_DynamicalScreening(const QMStateType& type,
                                          Orbitals& orb) {

  const tools::EigenSystem& BSECoefs = GetBSEEigenSystem(type, orb);
  const Eigen::VectorXd& RPAInputEnergies = orb.RPAInputEnergies();
  const bool full = !orb.getTDAApprox();

  const Eigen::VectorXd& BSEenergies = BSECoefs.eigenvalues();

  // initial copy of static BSE energies to dynamic
  Eigen::VectorXd BSEenergies_dynamic = BSEenergies;

  for (Index i_exc = 0; i_exc < BSEenergies.size(); i_exc++) {
    XTP_LOG(Log::info, log_) << "Dynamical Screening BSE, Excitation " << i_exc
                             << " static " << BSEenergies(i_exc) << flush;

    const Eigen::VectorXd X = BSECoefs.eigenvectors().col(i_exc);
    Eigen::VectorXd Y;
    std::vector<std::pair<const Eigen::VectorXd*, const Eigen::VectorXd*>>
        pairs{{&X, &X}};
    if (full) {
      Y = BSECoefs.eigenvectors2().col(i_exc);
      pairs.emplace_back(&Y, &Y);
    }
    const Eigen::MatrixXd Gd = DirectTermDensity(pairs);
    Eigen::MatrixXd G2;
    if (full) {
      G2 = CouplingTermDensity(Y, X);
    }
    // screened direct contribution for W(w), or for the static W(0)
    auto Hd = [&](const SymmetricEigenSystem* W) {
      double value = (W != nullptr) ? ContractScreened(*W, Gd)
                                    : epsilon_0_inv_.dot(Gd.diagonal());
      if (full) {
        value += 2 * ((W != nullptr) ? ContractScreened(*W, G2)
                                     : epsilon_0_inv_.dot(G2.diagonal()));
      }
      return -value;
    };
    const double Hd_static_contribution = Hd(nullptr);

    for (Index iter = 0; iter < max_dyn_iter_; iter++) {

      // screened interaction at the last energy as screening frequency
      double old_energy = BSEenergies_dynamic(i_exc);
      DFTTimings timings;
      const SymmetricEigenSystem W =
          ScreenedInteraction(RPAInputEnergies, old_energy, timings);

      // new energy perturbatively
      BSEenergies_dynamic(i_exc) =
          BSEenergies(i_exc) + Hd_static_contribution - Hd(&W);

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
          << boost::format(
                 "  S(dynamic) = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f "
                 "nm f = %4$+1.4f") %
                 (i + 1) % (hrt2ev * BSEenergies_dynamic(i)) %
                 (1240.0 / (hrt2ev * BSEenergies_dynamic(i))) %
                 (osc * BSEenergies_dynamic(i) / BSEenergies(i))
          << flush;
    }

  } else if (type == QMStateType::Triplet) {
    orb.BSETriplets_dynamic() = BSEenergies_dynamic;
    XTP_LOG(Log::error, log_) << "  ====== triplet energies with perturbative "
                                 "dynamical screening (eV) ====== "
                              << flush;
    for (Index i = 0; i < opt_.nmax; ++i) {
      XTP_LOG(Log::error, log_)
          << boost::format(
                 "  T(dynamic) = %1$4d Omega = %2$+1.12f eV  lamdba = %3$+3.2f "
                 "nm ") %
                 (i + 1) % (hrt2ev * BSEenergies_dynamic(i)) %
                 (1240.0 / (hrt2ev * BSEenergies_dynamic(i)))
          << flush;
    }

  } else {
    throw std::runtime_error(
        "Unsupported QMStateType in BSE::Perturbative_DynamicalScreening");
  }
}

}  // namespace xtp
}  // namespace votca
