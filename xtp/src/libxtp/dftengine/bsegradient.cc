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

// Analytic nuclear gradient of a BSE excitation energy (TDA or full) with
// static RPA screening and fixed (geometry-independent) QP corrections.
//
// Notation. p, q, r: MOs; i, j, k: occupied, a, b: virtual (all MOs of the
// closed-shell KS reference); v, c: occupied/virtual levels of the BSE window
// [vmin, cmax]; k, d: occupied/virtual levels of the RPA window
// [rpamin, rpamax]. B^P_pq = sum_R L_RP (R|pq) with L L^T = V^-1 (here
// L = V^-1/2; the energy is invariant under orthogonal rotations of the aux
// index, so VOTCA's Pseudo_InvSqrt_GWBSE frame gives the same Omega).
//
// Energy, written here for the TDA (X normalized, kappa = 1 for singlets, 0
// for triplets; the full BSE adds Y, see the amplitude section below):
//   Omega = sum_vc X_vc^2 (e_c - e_v) + 2 kappa |x|^2
//           - sum_{vc,v'c'} X_vc X_v'c' B_vv'^T M B_cc',
//   x^P = sum_vc X_vc B^P_vc,  M = eps^-1,
//   eps = 1 + sum_kd f(Delta_kd) b_kd b_kd^T,  f(D) = 4 D / (D^2 + eta^2),
//   Delta_kd = e'_d - e'_k,
// with e_p = eps_p + a_p (BSE) and e'_p = eps_p + pi_p (RPA), a and pi fixed.
// By Hellmann-Feynman dOmega = X^T dH X, so Omega depends on the geometry
// through B (integrals and MO coefficients) and the KS energies eps_p:
//   dOmega = sum_P g^P . dB^P + sum_p w_p deps_p,
// g^P = dOmega/dB^P (exchange, direct, and screening response through M),
// w_p = dOmega/deps_p (diagonal of the BSE difference density and of the
// screening response density, cf. Holzer, arXiv:2611.12617).
//
// dB^P has an integral part (three-centre and metric derivatives at fixed
// MO coefficients, ContractThreeCenterDerivatives and
// ContractCoulombMetricDerivatives with
// Gamma^R = C (sum_P L_RP g^P) C^T and G = -1/2 L (g . B)_sym L^T) and an
// orbital part sum_rp U_rp Q_rp with Q_rp = sum_{P,q} B^P_rq (g^P_pq +
// g^P_qp) for the MO response C -> C (1 + U). The canonical KS orbitals
// respond as
//   same space, r != p:  U_rp = (F'_rp - eps_p S_rp) / (eps_p - eps_r)
//   diagonal:            U_pp = -S_pp / 2,  deps_p = F'_pp - eps_p S_pp
//   occupied-virtual:    U_ia = -S_ia - U_ai, U_ai from coupled-perturbed KS
// with S_pq = C_p^T dS C_q, F' = F^(R) + G[dD], F^(R) the derivative of the
// Fock matrix at fixed AO density and G the closed-shell Fock response
// (J + 1/2 a_x K + f_xc) to the density change dD (dD_ai = dD_ia = 2 U_ai,
// dD_jk = -2 S_jk). Collecting terms:
//   D_rp = (Q_rp - Q_pr) / (2 (eps_p - eps_r))    same space, r != p
//   D_pp = w_p
//   L_ai = Q_ai - Q_ia + 4 G_ai[D],   (A + B) Z = L  (Z-vector)
//   Ptilde = D (occ-occ, virt-virt),  Ptilde_ai = Ptilde_ia = -Z_ai / 2
//   W_rp = (Q_rp eps_p - Q_pr eps_r) / (2 (eps_p - eps_r))   same space
//   W_pp = Q_pp / 2 + w_p eps_p,  W_jk += 2 G_jk[Ptilde],
//   W_ai = W_ia = (Q_ia - Z_ai eps_i) / 2
// and
//   dOmega/dR = (integral part) + sum Ptilde_AO dF/dR|_D - sum W_AO dS/dR.
// The Fock derivative at fixed density contracted with Ptilde is split into
// the one-electron part, RI-J and RI-K (exact, both are bilinear) and XC
// (d/dt of the ground-state XC gradient at D + t Ptilde, analytic in
// Vxc_Potential::GradientResponse), which contains the f_xc term of the moving
// basis functions and grid. No g_xc is needed: the BSE kernel has no XC.
//
// Degenerate KS levels within a space get D_rp = 0 (the antisymmetric part
// of Q vanishes for symmetry-related partners, not for accidental
// degeneracies).

// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/dftengine.h"
#include "votca/xtp/dftgradient.h"
#include "votca/xtp/eigen.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/out_of_window_shift.h"
#include "votca/xtp/riderivatives.h"
#include "votca/xtp/screening_kernels.h"
#include <chrono>
#include <functional>

namespace votca {
namespace xtp {

// Defined in libint2_derivative_calls.cc.
using AOMatrixDerivative = std::array<Eigen::MatrixXd, 3>;
std::vector<AOMatrixDerivative> ComputeOverlapDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeKineticDerivatives(
    const AOBasis& aobasis);
std::vector<AOMatrixDerivative> ComputeNuclearAttractionDerivatives(
    const AOBasis& aobasis, const QMMolecule& mol);
bool HasLibint2DerivativeSupport();

namespace {

Eigen::MatrixXd ContractAODerivatives(const std::vector<AOMatrixDerivative>& d,
                                      const Eigen::MatrixXd& P) {
  Eigen::MatrixXd grad = Eigen::MatrixXd::Zero(Index(d.size()), 3);
  for (Index a = 0; a < Index(d.size()); ++a) {
    for (Index xyz = 0; xyz < 3; ++xyz) {
      grad(a, xyz) = P.cwiseProduct(d[a][xyz]).sum();
    }
  }
  return grad;
}

// A <- A B for a square B, by blocks of rows (scratch of one block)
void RightMultiplyInPlace(Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  const Index rows = A.rows();
  const Index block = std::max<Index>(1, std::min<Index>(rows, 4096));
  for (Index r0 = 0; r0 < rows; r0 += block) {
    const Index n = std::min(block, rows - r0);
    Eigen::MatrixXd tmp = A.middleRows(r0, n) * B;
    A.middleRows(r0, n) = tmp;
  }
}

}  // namespace

DFTEngine::BSEGradientResult DFTEngine::BSEExcitationGradient(
    const Orbitals& orb, const BSEGradientInput& in) {
  if (auxbasis_name_.empty()) {
    throw std::runtime_error(
        "BSEExcitationGradient: needs the RI ground state (auxbasisset)");
  }
  if (!ecp_name_.empty()) {
    throw std::runtime_error("BSEExcitationGradient: ECPs are not supported");
  }
  if (!HasLibint2DerivativeSupport()) {
    throw std::runtime_error(
        "BSEExcitationGradient: libint2 lacks derivative integrals");
  }
  Orbitals work = orb;
  Prepare(work);
  if (num_alpha_electrons_ != num_beta_electrons_) {
    throw std::runtime_error(
        "BSEExcitationGradient: closed-shell ground states only");
  }
  const QMMolecule& mol = work.QMAtoms();
  const Vxc_Potential<Vxc_Grid> vxc = SetupVxc(mol);
  // wall time of the phases, reported at the end
  std::vector<std::pair<std::string, double>> phase_times;
  auto phase_start = std::chrono::steady_clock::now();
  auto phase = [&](const std::string& name) {
    const auto now = std::chrono::steady_clock::now();
    phase_times.emplace_back(
        name, std::chrono::duration<double>(now - phase_start).count());
    phase_start = now;
  };

  const Eigen::MatrixXd& C = orb.MOs().eigenvectors();
  const Eigen::VectorXd& eps = orb.MOs().eigenvalues();
  const Index nmo = C.cols();
  const Index nocc = num_docc_;
  const Index nvirt = nmo - nocc;
  const Index homo = nocc - 1;
  const Eigen::MatrixXd C_occ = C.leftCols(nocc);
  const Eigen::MatrixXd C_virt = C.rightCols(nvirt);
  const Eigen::MatrixXd D0 = 2.0 * C_occ * C_occ.transpose();

  // windows
  if (in.vmin < 0 || in.vmin > homo || in.cmax <= homo || in.cmax >= nmo ||
      in.rpamin < 0 || in.rpamin > homo || in.rpamax <= homo ||
      in.rpamax >= nmo) {
    throw std::runtime_error("BSEExcitationGradient: invalid windows");
  }
  const Index vt = homo - in.vmin + 1;
  const Index ct = in.cmax - homo;
  if (in.X.size() != vt * ct || (in.Y.size() != 0 && in.Y.size() != vt * ct) ||
      in.energies_bse.size() != in.cmax - in.vmin + 1 ||
      in.energies_rpa.size() != in.rpamax - in.rpamin + 1) {
    throw std::runtime_error("BSEExcitationGradient: input sizes");
  }
  if (in.qp_response &&
      (in.qpmin < in.rpamin || in.qpmax > in.rpamax || in.qpmin > in.qpmax ||
       in.qp_energies.size() != in.qpmax - in.qpmin + 1)) {
    throw std::runtime_error(
        "BSEExcitationGradient: QP window must lie in the RPA window");
  }
  const Index lo = std::min(in.vmin, in.rpamin);
  const Index hi = std::max(in.cmax, in.rpamax);
  const Index nw = hi - lo + 1;
  const Eigen::MatrixXd Cw = C.middleCols(lo, nw);
  const Index v0 = in.vmin - lo;   // window offsets
  const Index c0 = homo + 1 - lo;  // first virtual
  const Index k0 = in.rpamin - lo;
  const Index nk = homo - in.rpamin + 1;
  const Index nd = in.rpamax - homo;
  const double kappa = in.singlet ? 1.0 : 0.0;

  // energies used in the Hamiltonian, by window position
  Eigen::VectorXd eA = Eigen::VectorXd::Zero(nw);
  eA.segment(v0, in.cmax - in.vmin + 1) = in.energies_bse;
  Eigen::VectorXd eP = Eigen::VectorXd::Zero(nw);
  eP.segment(k0, in.rpamax - in.rpamin + 1) = in.energies_rpa;

  // ---- B^P_pq over the window, L = V^-1/2
  AOCoulomb aocoulomb;
  aocoulomb.Fill(auxbasis_);
  const Index naux = auxbasis_.AOBasisSize();
  Eigen::MatrixXd L;
  {
    const SymmetricEigenSystem es = SymmetricEigen(aocoulomb.Matrix());
    if (es.values.minCoeff() <= 0.0) {
      throw std::runtime_error(
          "BSEExcitationGradient: aux Coulomb metric not positive definite");
    }
    L = es.vectors * es.values.cwiseInverse().cwiseSqrt().asDiagonal() *
        es.vectors.transpose();
  }
  // (R|pq) first, then B = J L in place by blocks of rows: one array of
  // nw^2 x Naux (the size of the GW-BSE integrals) instead of two
  Eigen::MatrixXd Bmat(nw * nw, naux);
  ForEachThreeCenterIntegral(
      auxbasis_, dftbasis_,
      [&Cw, &Bmat, nw](Index R, const Eigen::MatrixXd& I) {
        Eigen::MatrixXd m = Cw.transpose() * (I * Cw);
        Bmat.col(R) = Eigen::Map<const Eigen::VectorXd>(m.data(), nw * nw);
      });
  RightMultiplyInPlace(Bmat, L);
  auto B = [&Bmat, nw](Index P) {
    return Eigen::Map<const Eigen::MatrixXd>(Bmat.col(P).data(), nw, nw);
  };

  // X (and Y for the full BSE) as (v, c); the energy functional
  //   Omega = (X A X + Y A Y + 2 X B Y) / (X X - Y Y)
  // is stationary in X, Y, so dOmega = (X dA X + Y dA Y + 2 X dB Y) / N
  const bool full = in.Y.size() > 0;
  const Eigen::MatrixXd Xvc =
      Eigen::Map<const Eigen::MatrixXd>(in.X.data(), ct, vt).transpose();
  const Eigen::MatrixXd Yvc =
      full ? Eigen::MatrixXd(
                 Eigen::Map<const Eigen::MatrixXd>(in.Y.data(), ct, vt)
                     .transpose())
           : Eigen::MatrixXd::Zero(vt, ct);
  const double norm = Xvc.squaredNorm() - Yvc.squaredNorm();
  if (norm <= 0.0) {
    throw std::runtime_error("BSEExcitationGradient: X X - Y Y <= 0");
  }

  phase("three-index integrals B");
  // ---- energy and its partial derivatives g^P (nw x nw per P), w
  Eigen::MatrixXd gmat = Eigen::MatrixXd::Zero(nw * nw, naux);
  auto g = [&gmat, nw](Index P) {
    return Eigen::Map<Eigen::MatrixXd>(gmat.col(P).data(), nw, nw);
  };
  Eigen::VectorXd w = Eigen::VectorXd::Zero(nw);  // dOmega/deps, window

  // dOmega/de of the BSE energies (wA) and of the RPA energies (w)
  Eigen::VectorXd wA = Eigen::VectorXd::Zero(nw);
  double omega_qp = 0.0;
  for (Index v = 0; v < vt; ++v) {
    for (Index c = 0; c < ct; ++c) {
      const double x2 = Xvc(v, c) * Xvc(v, c) + Yvc(v, c) * Yvc(v, c);
      omega_qp += x2 * (eA(c0 + c) - eA(v0 + v));
      wA(c0 + c) += x2;
      wA(v0 + v) -= x2;
    }
  }

  // exchange, A and B: 2 kappa |x + y|^2
  const Eigen::MatrixXd XpY = Xvc + Yvc;
  Eigen::VectorXd x(naux);
  for (Index P = 0; P < naux; ++P) {
    x(P) = XpY.cwiseProduct(B(P).block(v0, c0, vt, ct)).sum();
  }
  const double omega_x = 2.0 * kappa * x.squaredNorm();
  for (Index P = 0; P < naux; ++P) {
    g(P).block(v0, c0, vt, ct) += 4.0 * kappa * x(P) * XpY;
  }

  // static RPA screening
  Eigen::MatrixXd Bkd(nk * nd, naux);
  Eigen::VectorXd f(nk * nd);
  Eigen::VectorXd fprime(nk * nd);
  const double eta2 = in.eta * in.eta;
  for (Index d = 0; d < nd; ++d) {
    for (Index k = 0; k < nk; ++k) {
      const double Delta = eP(c0 + d) - eP(k0 + k);
      const double den = Delta * Delta + eta2;
      f(k + nk * d) = 4.0 * Delta / den;
      fprime(k + nk * d) = 4.0 * (eta2 - Delta * Delta) / (den * den);
    }
  }
  for (Index P = 0; P < naux; ++P) {
    const auto b = B(P);
    for (Index d = 0; d < nd; ++d) {
      Bkd.col(P).segment(nk * d, nk) = b.block(k0, c0 + d, nk, 1);
    }
  }
  Eigen::MatrixXd epsilon = Bkd.transpose() * f.asDiagonal() * Bkd;
  epsilon.diagonal().array() += 1.0;
  Eigen::MatrixXd M;
  {
    SymmetricEigenSystem es = SymmetricEigen(epsilon);
    for (Index i = 0; i < es.values.size(); ++i) {
      es.values(i) = (es.values(i) > 1e-8) ? 1.0 / es.values(i) : 0.0;
    }
    M = es.vectors * es.values.asDiagonal() * es.vectors.transpose();
  }

  // direct term of A: Y^Q = X B^Q_cc X^T + Y B^Q_cc Y^T (v x v),
  // Z^P = X^T B^P_vv X + Y^T B^P_vv Y (c x c)
  Eigen::MatrixXd Bvv(vt * vt, naux);
  Eigen::MatrixXd Ymat(vt * vt, naux);
  Eigen::MatrixXd Zmat(ct * ct, naux);
  for (Index P = 0; P < naux; ++P) {
    const auto b = B(P);
    const Eigen::MatrixXd bvv = b.block(v0, v0, vt, vt);
    const Eigen::MatrixXd bcc = b.block(c0, c0, ct, ct);
    Eigen::MatrixXd Yq = Xvc * bcc * Xvc.transpose();
    Eigen::MatrixXd Zc = Xvc.transpose() * bvv * Xvc;
    if (full) {
      Yq += Yvc * bcc * Yvc.transpose();
      Zc += Yvc.transpose() * bvv * Yvc;
    }
    Bvv.col(P) = Eigen::Map<const Eigen::VectorXd>(bvv.data(), vt * vt);
    Ymat.col(P) = Eigen::Map<const Eigen::VectorXd>(Yq.data(), vt * vt);
    Zmat.col(P) = Eigen::Map<const Eigen::VectorXd>(Zc.data(), ct * ct);
  }
  Eigen::MatrixXd GammaW = Bvv.transpose() * Ymat;  // G_PQ = B^P . Y^Q
  {
    const Eigen::MatrixXd MY = Ymat * M;
    const Eigen::MatrixXd MZ = Zmat * M;
    for (Index P = 0; P < naux; ++P) {
      g(P).block(v0, v0, vt, vt) -=
          Eigen::Map<const Eigen::MatrixXd>(MY.col(P).data(), vt, vt);
      g(P).block(c0, c0, ct, ct) -=
          Eigen::Map<const Eigen::MatrixXd>(MZ.col(P).data(), ct, ct);
    }
  }
  // direct term of B (full BSE): 2 X Hd2 Y with
  // Hd2(v1c1, v2c2) = -B_{v1c2}^T M B_{v2c1}, i.e.
  //   -2 sum_PQ M_PQ B^P_vc . T^Q,  T^Q = X (B^Q_vc)^T Y  (v x c)
  // and dT^Q/dB^Q gives T'^P = Y (B^P_vc)^T X.
  if (full) {
    Eigen::MatrixXd Bvc(vt * ct, naux);
    Eigen::MatrixXd Tmat(vt * ct, naux);
    Eigen::MatrixXd Tpmat(vt * ct, naux);
    for (Index P = 0; P < naux; ++P) {
      const Eigen::MatrixXd bvc = B(P).block(v0, c0, vt, ct);
      const Eigen::MatrixXd T = Xvc * bvc.transpose() * Yvc;
      const Eigen::MatrixXd Tp = Yvc * bvc.transpose() * Xvc;
      Bvc.col(P) = Eigen::Map<const Eigen::VectorXd>(bvc.data(), vt * ct);
      Tmat.col(P) = Eigen::Map<const Eigen::VectorXd>(T.data(), vt * ct);
      Tpmat.col(P) = Eigen::Map<const Eigen::VectorXd>(Tp.data(), vt * ct);
    }
    GammaW += 2.0 * Bvc.transpose() * Tmat;
    const Eigen::MatrixXd MT = (Tmat + Tpmat) * M;
    for (Index P = 0; P < naux; ++P) {
      g(P).block(v0, c0, vt, ct) -=
          2.0 * Eigen::Map<const Eigen::MatrixXd>(MT.col(P).data(), vt, ct);
    }
  }
  const double omega_d = -(M.cwiseProduct(GammaW)).sum();
  GammaW = 0.5 * (GammaW + GammaW.transpose()).eval();
  // screening response: dOmega = Tr(Gbar d eps), Gbar = M GammaW M
  {
    const Eigen::MatrixXd Gbar = M * GammaW * M;
    const Eigen::MatrixXd GB = Bkd * Gbar;  // (kd x P)
    for (Index P = 0; P < naux; ++P) {
      for (Index d = 0; d < nd; ++d) {
        g(P).block(k0, c0 + d, nk, 1) +=
            2.0 *
            f.segment(nk * d, nk).cwiseProduct(GB.col(P).segment(nk * d, nk));
      }
    }
    const Eigen::VectorXd bGb = Bkd.cwiseProduct(GB).rowwise().sum();
    for (Index d = 0; d < nd; ++d) {
      for (Index k = 0; k < nk; ++k) {
        const double t = fprime(k + nk * d) * bGb(k + nk * d);
        w(c0 + d) += t;
        w(k0 + k) -= t;
      }
    }
  }

  phase("BSE energy and partial derivatives");
  // normalization of the full BSE (1 for TDA and normalized vectors)
  gmat /= norm;
  w /= norm;
  wA /= norm;
  BSEGradientResult result;
  result.omega = (omega_qp + omega_x + omega_d) / norm;

  // u_p: coefficients of sum_p u_p <p|V_xc|p> in the energy (QP response)
  Eigen::VectorXd u = Eigen::VectorXd::Zero(nmo);
  if (!in.qp_response) {
    // all energies are KS energies plus fixed shifts
    w += wA;
  } else {
    // ---- response of the QP corrections, plasmon-pole model of VOTCA
    // (Sigma_PPM, PPM::PPM_construct_parameters). The QP energies E of
    // qpmin..qpmax solve E_p = eps_p + sigma_p(E_p; R), sigma = (1 - a_x)
    // Sigma_x - v_xc + Sigma_c, with the screening and the poles of G at
    // the RPA energies R: G0W0: R = KS energies plus fixed shifts; evGW:
    // R = E in the QP window and, outside, KS energies plus the
    // out-of-window shift of the corrections E - eps (OutOfWindowShift).
    // With the weights W of dOmega on E, the adjoint (1 - J)^T v = W,
    // J = dsigma/dE (through omega_p = E_p and, for evGW, R), turns
    // W . dE into sum_p v_p (deps_p + dsigma_p|_{omega, R}) plus the
    // explicit dependence of R on eps; G0W0: v_p = Z_p W_p.
    // Sigma_c(omega) = sum_i a_i sum_n (Bt^i_pn)^2 h(omega - e_n +- Omega_i),
    // h(t) = t / (t^2 + eta^2), a_i = w_i Omega_i / 2, Bt = B Phi with
    // eps(0) = Phi diag(lambda) Phi^T, w_i = 1 - 1/lambda_i,
    // Omega_i^2 = -nu^2 (y_i - 1) / (y_i - 1 + w_i), y_i = phi_i^T eps(i nu)^-1
    // phi_i; + for occupied, - for virtual n. The reverse pass goes through
    // all of these down to B and R.
    const Index q0 = in.qpmin - lo;
    const Index nq = in.qpmax - in.qpmin + 1;
    const Index r0 = k0;
    const Index nr = in.rpamax - in.rpamin + 1;
    const Index qr = in.qpmin - in.rpamin;  // QP level p is RPA level qr + p
    const double nu2 = in.ppm_imag * in.ppm_imag;
    const double seta2 = in.sigma_eta * in.sigma_eta;
    const double ax = ScaHFX_;

    using ReversePass =
        std::function<Eigen::VectorXd(const Eigen::VectorXd&, bool)>;
    // the adjoint and the bookkeeping of the energy weights, the same for
    // all self-energy models: dsigma = Sigma_c'(E_p), reverse(v, accumulate)
    // the reverse pass of sum_p v_p sigma_p(omega_p; R) at fixed omega
    // (returns d/dR; with accumulate also adds d/dB, d/deps_p and v_xc),
    // poles (nr x nq): column p the part of d sigma_p/dR through the pole
    // energies of G alone (cheap, preconditioner of the evGW adjoint)
    auto finish = [&](const Eigen::VectorXd& dsigma, const ReversePass& reverse,
                      const Eigen::MatrixXd& poles) {
      // weights of dOmega: wE on the QP energies, wR on the RPA energies
      // (screening of the BSE; BSE levels outside the QP window use R)
      Eigen::VectorXd wR = w.segment(r0, nr);
      w.setZero();
      Eigen::VectorXd wE = wA.segment(q0, nq);
      for (Index j = 0; j < nw; ++j) {
        if (j >= q0 && j < q0 + nq) {
          continue;
        }
        if (j < r0 || j >= r0 + nr) {
          if (wA(j) != 0.0) {
            throw std::runtime_error(
                "BSEExcitationGradient: BSE level outside the RPA window");
          }
          continue;
        }
        wR(j - r0) += wA(j);
      }

      if (!in.qp_selfconsistent) {
        // G0W0: R = eps + fixed shifts
        const Eigen::VectorXd v =
            wE.cwiseQuotient((1.0 - dsigma.array()).matrix());
        const Eigen::VectorXd rho = reverse(v, true);
        w.segment(r0, nr) += wR + rho;
      } else {
        // evGW: Jacobian of the out-of-window RPA energies with respect to
        // the corrections c = E - eps of the QP window (piecewise linear in c)
        OutOfWindowShift shift;
        shift.mode = OutOfWindowShift::Parse(in.out_of_window_shift);
        shift.fit_fraction = in.out_of_window_fit_fraction;
        // R(c, eps): the linear shift also depends on the KS energies (the
        // abscissae of its fit), so both Jacobians are taken
        auto rpa_energies_eps = [&](const Eigen::VectorXd& c,
                                    const Eigen::VectorXd& e) {
          Eigen::VectorXd R = e.segment(in.rpamin, nr);
          R.segment(qr, nq) += c;
          shift.Apply(R, e, in.rpamin, in.rpamax, homo, in.qpmin, in.qpmax);
          return R;
        };
        auto rpa_energies = [&](const Eigen::VectorXd& c) {
          return rpa_energies_eps(c, eps);
        };
        const Eigen::VectorXd c_ref =
            in.qp_energies - eps.segment(in.qpmin, nq);
        const Eigen::VectorXd R_ref = rpa_energies(c_ref);
        result.rpa_rebuild_error =
            (R_ref - in.energies_rpa).cwiseAbs().maxCoeff();
        Eigen::MatrixXd Js = Eigen::MatrixXd::Zero(nr, nq);  // dR/dc
        const double dc = 1e-7;
        for (Index q = 0; q < nq; ++q) {
          Eigen::VectorXd cp = c_ref;
          cp(q) += dc;
          Eigen::VectorXd cm = c_ref;
          cm(q) -= dc;
          Js.col(q) = (rpa_energies(cp) - rpa_energies(cm)) / (2.0 * dc);
        }
        Js.middleRows(qr, nq).setZero();  // inside the window R = E itself
        // W . dE with W_q = wE_q + wR_{window q} + sum_o wR_o Js_oq
        const Eigen::VectorXd Wbar =
            wE + wR.segment(qr, nq) + Js.transpose() * wR;
        // adjoint: v_q (1 - Sigma_c'_q) - rho_q(v) - sum_o rho_o(v) Js_oq =
        // Wbar_q, linear in v. A fixed-point iteration need not converge (it
        // mirrors the evGW iteration, which VOTCA stabilizes by mixing), so
        // GMRES, right-preconditioned with the same operator restricted to
        // the pole energies of G (screening frozen); one reverse pass per
        // iteration. Fallback: the matrix column by column (nq passes).
        const Eigen::VectorXd diag1 = (1.0 - dsigma.array()).matrix();
        auto apply = [&](const Eigen::VectorXd& x) {
          const Eigen::VectorXd rho = reverse(x, false);
          return Eigen::VectorXd(x.cwiseProduct(diag1) - rho.segment(qr, nq) -
                                 Js.transpose() * rho);
        };
        Eigen::MatrixXd M = -poles.middleRows(qr, nq) - Js.transpose() * poles;
        M.diagonal() += diag1;
        const Eigen::PartialPivLU<Eigen::MatrixXd> luM(M);
        const double bnorm = Wbar.norm();
        const double tol = 1e-11 * std::max(1.0, bnorm);
        Eigen::VectorXd v = Eigen::VectorXd::Zero(nq);
        Index iterations = 0;
        bool converged = (bnorm == 0.0);
        if (!converged) {
          const Index max_krylov = std::min<Index>(nq, 60);
          Eigen::MatrixXd Vk(nq, max_krylov + 1);
          Eigen::MatrixXd Hk =
              Eigen::MatrixXd::Zero(max_krylov + 1, max_krylov);
          Vk.col(0) = Wbar / bnorm;
          for (Index j = 0; j < max_krylov; ++j) {
            Eigen::VectorXd z = apply(luM.solve(Vk.col(j)));
            for (Index pass = 0; pass < 2; ++pass) {  // Gram-Schmidt, twice
              for (Index i = 0; i <= j; ++i) {
                const double hij = Vk.col(i).dot(z);
                Hk(i, j) += hij;
                z -= hij * Vk.col(i);
              }
            }
            Hk(j + 1, j) = z.norm();
            iterations = j + 1;
            // least squares min |bnorm e1 - H y|
            Eigen::VectorXd rhs = Eigen::VectorXd::Zero(j + 2);
            rhs(0) = bnorm;
            const Eigen::MatrixXd Hj = Hk.topLeftCorner(j + 2, j + 1);
            const Eigen::VectorXd y = Hj.colPivHouseholderQr().solve(rhs);
            const double res = (rhs - Hj * y).norm();
            if (res < tol || Hk(j + 1, j) < 1e-14 * bnorm) {
              v = luM.solve(Vk.leftCols(j + 1) * y);
              converged = true;
              break;
            }
            Vk.col(j + 1) = z / Hk(j + 1, j);
          }
        }
        if (converged &&
            (apply(v) - Wbar).norm() > 1e-8 * std::max(1.0, bnorm)) {
          converged = false;
        }
        if (converged) {
          XTP_LOG(Log::error, *pLog_)
              << TimeStamp() << " evGW adjoint: GMRES converged in "
              << iterations << " iterations" << std::flush;
        } else {
          XTP_LOG(Log::error, *pLog_)
              << TimeStamp() << " evGW adjoint: GMRES did not converge in "
              << iterations << " iterations, building the matrix" << std::flush;
          Eigen::MatrixXd A(nq, nq);
          for (Index q = 0; q < nq; ++q) {
            Eigen::VectorXd e = Eigen::VectorXd::Zero(nq);
            e(q) = 1.0;
            A.col(q) = apply(e);
          }
          const Eigen::PartialPivLU<Eigen::MatrixXd> lu(A);
          v = lu.solve(Wbar);
          if ((A * v - Wbar).cwiseAbs().maxCoeff() >
              1e-8 * std::max(1.0, Wbar.cwiseAbs().maxCoeff())) {
            throw std::runtime_error(
                "BSEExcitationGradient: the evGW adjoint equation is "
                "singular");
          }
        }
        const Eigen::VectorXd rho = reverse(v, true);
        // explicit eps dependence of the out-of-window R at fixed E (the dE
        // part is in the adjoint): K_ok = dR_o/deps_k with c = E - eps
        Eigen::VectorXd wRo = wR + rho;
        wRo.segment(qr, nq).setZero();
        Eigen::MatrixXd K = Eigen::MatrixXd::Zero(nr, nr);
        const double de = 1e-7;
        for (Index k = 0; k < nr; ++k) {
          Eigen::VectorXd ep = eps;
          ep(in.rpamin + k) += de;
          Eigen::VectorXd em = eps;
          em(in.rpamin + k) -= de;
          const Eigen::VectorXd cp = in.qp_energies - ep.segment(in.qpmin, nq);
          const Eigen::VectorXd cm = in.qp_energies - em.segment(in.qpmin, nq);
          K.col(k) = (rpa_energies_eps(cp, ep) - rpa_energies_eps(cm, em)) /
                     (2.0 * de);
        }
        w.segment(r0, nr) += K.transpose() * wRo;
      }
    };

    if (in.sigma_integrator == "ppm") {

      const SymmetricEigenSystem es0 = SymmetricEigen(epsilon);
      const Eigen::MatrixXd& Phi = es0.vectors;
      const Eigen::VectorXd& lam = es0.values;
      Eigen::VectorXd fi(nk * nd);
      Eigen::VectorXd fiprime(nk * nd);
      for (Index d = 0; d < nd; ++d) {
        for (Index k = 0; k < nk; ++k) {
          const double Delta = eP(c0 + d) - eP(k0 + k);
          const double den = Delta * Delta + nu2;
          fi(k + nk * d) = 4.0 * Delta / den;
          fiprime(k + nk * d) = 4.0 * (nu2 - Delta * Delta) / (den * den);
        }
      }
      Eigen::MatrixXd epsi = Bkd.transpose() * fi.asDiagonal() * Bkd;
      epsi.diagonal().array() += 1.0;
      const Eigen::MatrixXd EiPhi = InverseSPD(epsi) * Phi;

      Eigen::VectorXd wppm = Eigen::VectorXd::Zero(naux);
      Eigen::VectorXd Om = Eigen::VectorXd::Zero(naux);
      Eigen::VectorXd a = Eigen::VectorXd::Zero(naux);
      Eigen::VectorXd nom = Eigen::VectorXd::Zero(naux);
      Eigen::VectorXd sgn = Eigen::VectorXd::Ones(naux);
      for (Index i = 0; i < naux; ++i) {
        const double wi = 1.0 - 1.0 / lam(i);
        if (wi < 1e-5) {
          continue;  // mode dropped, as in PPM_construct_parameters
        }
        wppm(i) = wi;
        nom(i) = Phi.col(i).dot(EiPhi.col(i)) - 1.0;
        const double frac = -nom(i) / (nom(i) + wi) * nu2;
        sgn(i) = (frac < 0.0) ? -1.0 : 1.0;
        Om(i) = std::sqrt(std::abs(frac));
        a(i) = 0.5 * wi * Om(i);
      }
      const Eigen::VectorXd epole = eP.segment(r0, nr);
      Eigen::VectorXd s_n = -Eigen::VectorXd::Ones(nr);
      s_n.head(nk).setOnes();

      const Eigen::MatrixXd Vxc_mo =
          C.transpose() * vxc.IntegrateVXC(D0).matrix() * C;

      // per QP level: B_pn (rows of B), Bt = B_p Phi (stored, nr x Naux per
      // level), pole factors h, h' (recomputed when needed), Sigma_c'(E_p)
      auto gatherB = [&](Index p) {
        Eigen::MatrixXd Bp(nr, naux);
        for (Index n = 0; n < nr; ++n) {
          Bp.row(n) = Bmat.row((q0 + p) + nw * (r0 + n));
        }
        return Bp;
      };
      auto poleFactors = [&](Index p, Eigen::ArrayXXd& H, Eigen::ArrayXXd& Hp) {
        const double omega_p = in.qp_energies(p);
        H.resize(nr, naux);
        Hp.resize(nr, naux);
        for (Index i = 0; i < naux; ++i) {
          for (Index n = 0; n < nr; ++n) {
            const double t = omega_p - epole(n) + s_n(n) * Om(i);
            const double den = t * t + seta2;
            H(n, i) = t / den;
            Hp(n, i) = (seta2 - t * t) / (den * den);
          }
        }
      };
      std::vector<Eigen::MatrixXd> Bts(static_cast<std::size_t>(nq));
      Eigen::VectorXd dsigma(nq);
      Eigen::MatrixXd poles(nr, nq);
      for (Index p = 0; p < nq; ++p) {
        const Index pabs = in.qpmin + p;
        const double omega_p = in.qp_energies(p);
        const Eigen::MatrixXd Bp = gatherB(p);
        Bts[std::size_t(p)] = Bp * Phi;
        Eigen::ArrayXXd H, Hp;
        poleFactors(p, H, Hp);
        const Eigen::ArrayXXd aBt2 =
            Bts[std::size_t(p)].array().square().rowwise() *
            a.transpose().array();
        const double sigma_c = (aBt2 * H).sum();
        dsigma(p) = (aBt2 * Hp).sum();
        poles.col(p) = -(aBt2 * Hp).rowwise().sum().matrix();
        const double sigma_x = -(1.0 - ax) * Bp.topRows(nk).squaredNorm();
        const double e_rebuilt =
            eps(pabs) + sigma_x - Vxc_mo(pabs, pabs) + sigma_c;
        result.qp_rebuild_error =
            std::max(result.qp_rebuild_error, std::abs(e_rebuilt - omega_p));
      }

      // reverse pass of sum_p v_p sigma_p(omega_p; R) at fixed omega: returns
      // its derivative with respect to R (rho, RPA window); with accumulate,
      // also adds the derivatives with respect to B (g), eps_p and v_xc
      auto reverse = [&](const Eigen::VectorXd& vq, bool accumulate) {
        Eigen::VectorXd rho = Eigen::VectorXd::Zero(nr);
        Eigen::VectorXd gOm = Eigen::VectorXd::Zero(naux);
        Eigen::VectorXd gw = Eigen::VectorXd::Zero(naux);
        Eigen::MatrixXd GPhi = Eigen::MatrixXd::Zero(naux, naux);
        for (Index p = 0; p < nq; ++p) {
          const double v = vq(p);
          if (v == 0.0) {
            continue;
          }
          const Eigen::MatrixXd Bp = gatherB(p);
          const Eigen::MatrixXd& Bt = Bts[std::size_t(p)];
          Eigen::ArrayXXd H, Hp;
          poleFactors(p, H, Hp);
          const Eigen::ArrayXXd Bt2 = Bt.array().square();
          const Eigen::ArrayXXd aBt2 = Bt2.rowwise() * a.transpose().array();
          const Eigen::ArrayXXd BtH = Bt.array() * H;
          const Eigen::MatrixXd GBt =
              (2.0 * v) * (BtH.rowwise() * a.transpose().array()).matrix();
          if (accumulate) {
            w(q0 + p) += v;        // eps_p
            u(in.qpmin + p) -= v;  // -v_xc
            const Eigen::MatrixXd gBp = GBt * Phi.transpose();
            for (Index n = 0; n < nr; ++n) {
              for (Index P = 0; P < naux; ++P) {
                double gpn = gBp(n, P);
                if (n < nk) {  // Sigma_x: -(1 - a_x) sum_{i occ} (B_pi)^2
                  gpn += -2.0 * (1.0 - ax) * v * Bp(n, P);
                }
                g(P)(q0 + p, r0 + n) += gpn;
              }
            }
          }
          GPhi.noalias() += Bp.transpose() * GBt;
          rho += v * poles.col(p);
          for (Index i = 0; i < naux; ++i) {
            if (wppm(i) == 0.0) {
              continue;
            }
            double sH = 0.0;
            double sHp = 0.0;
            for (Index n = 0; n < nr; ++n) {
              sH += Bt2(n, i) * H(n, i);
              sHp += Bt2(n, i) * s_n(n) * Hp(n, i);
            }
            gOm(i) += v * (0.5 * wppm(i) * sH + a(i) * sHp);
            gw(i) += v * 0.5 * Om(i) * sH;
          }
        }
        // PPM parameters -> eps(0), eps(i nu)
        Eigen::VectorXd gy = Eigen::VectorXd::Zero(naux);
        Eigen::VectorXd glam = Eigen::VectorXd::Zero(naux);
        for (Index i = 0; i < naux; ++i) {
          if (wppm(i) == 0.0) {
            continue;
          }
          const double d2 = (nom(i) + wppm(i)) * (nom(i) + wppm(i));
          const double dOm_dfrac = sgn(i) / (2.0 * Om(i));
          gy(i) = gOm(i) * dOm_dfrac * (-nu2 * wppm(i) / d2);
          const double gwi = gw(i) + gOm(i) * dOm_dfrac * (nu2 * nom(i) / d2);
          glam(i) = gwi / (lam(i) * lam(i));
          GPhi.col(i) += 2.0 * gy(i) * EiPhi.col(i);
        }
        const Eigen::MatrixXd Gepsi =
            -(EiPhi * gy.asDiagonal() * EiPhi.transpose());
        Eigen::MatrixXd Mx = Phi.transpose() * GPhi;
        for (Index j = 0; j < naux; ++j) {
          for (Index i = 0; i < naux; ++i) {
            if (i == j) {
              Mx(i, i) = glam(i);
              continue;
            }
            const double dl = lam(j) - lam(i);
            Mx(i, j) = (std::abs(dl) > 1e-12) ? Mx(i, j) / dl : 0.0;
          }
        }
        Eigen::MatrixXd Geps0 = Phi * Mx * Phi.transpose();
        Geps0 = (0.5 * (Geps0 + Geps0.transpose())).eval();
        // eps = 1 + sum_kd f(Delta_kd) b b^T for both frequencies
        const Eigen::MatrixXd GB0 = Bkd * Geps0;
        const Eigen::MatrixXd GBi = Bkd * Gepsi;
        if (accumulate) {
          for (Index P = 0; P < naux; ++P) {
            for (Index d = 0; d < nd; ++d) {
              g(P).block(k0, c0 + d, nk, 1) +=
                  2.0 * (f.segment(nk * d, nk)
                             .cwiseProduct(GB0.col(P).segment(nk * d, nk)) +
                         fi.segment(nk * d, nk)
                             .cwiseProduct(GBi.col(P).segment(nk * d, nk)));
            }
          }
        }
        const Eigen::VectorXd bGb0 = Bkd.cwiseProduct(GB0).rowwise().sum();
        const Eigen::VectorXd bGbi = Bkd.cwiseProduct(GBi).rowwise().sum();
        for (Index d = 0; d < nd; ++d) {
          for (Index k = 0; k < nk; ++k) {
            const Index kd = k + nk * d;
            const double t = fprime(kd) * bGb0(kd) + fiprime(kd) * bGbi(kd);
            rho(nk + d) += t;
            rho(k) -= t;
          }
        }
        return rho;
      };

      finish(dsigma, reverse, poles);
    } else if (in.sigma_integrator == "exact") {
      // Sigma_c of Sigma_Exact: the RPA pair eigenproblem
      // C = Delta^1/2 (Delta + 4 b b^T) Delta^1/2 = F diag(Omega^2) F^T over
      // the RPA transitions (k, d), modes Z_s = sum_kd b_kd (X+Y)_kd,s with
      // (X+Y)_s = Delta^1/2 F_s Omega_s^-1/2, and
      //   Sigma_c(omega) = 2 sum_s sum_n (B_pn . Z_s)^2 h(omega - e_n +-
      //   Omega_s)
      // = 2 sum_n k_n^T phi_n(C) k_n, k_n = (B_pn . b_kd) Delta_kd^1/2,
      // phi_n(lambda) = lambda^-1/2 h(omega - e_n +- lambda^1/2): a matrix
      // function of C, differentiated with the divided differences of
      // phi_n in the eigenbasis of C (Daleckii-Krein).
      const Index S = nk * nd;
      Eigen::VectorXd Delta(S);
      for (Index d = 0; d < nd; ++d) {
        for (Index k = 0; k < nk; ++k) {
          Delta(k + nk * d) = eP(c0 + d) - eP(k0 + k);
        }
      }
      const Eigen::VectorXd sD = Delta.cwiseSqrt();
      Eigen::MatrixXd Cm = 4.0 * Bkd * Bkd.transpose();
      Cm.diagonal() += Delta;
      Cm = (sD.asDiagonal() * Cm * sD.asDiagonal()).eval();
      const SymmetricEigenSystem esC = SymmetricEigen(Cm);
      const Eigen::MatrixXd& F = esC.vectors;
      const Eigen::VectorXd& lamC = esC.values;
      if (lamC.minCoeff() <= 0.0) {
        throw std::runtime_error(
            "BSEExcitationGradient: non-positive RPA eigenvalue");
      }
      const Eigen::ArrayXd OmC = lamC.array().sqrt();
      // 1/(lambda_s - lambda_t), zero for (near) degenerate pairs, which
      // are handled with the derivative
      Eigen::MatrixXd invdl(S, S);
      std::vector<std::pair<Index, Index>> degenerate;
      for (Index t = 0; t < S; ++t) {
        for (Index s2 = 0; s2 < S; ++s2) {
          const double dl = lamC(s2) - lamC(t);
          if (s2 != t && std::abs(dl) > 1e-10 * std::max(1.0, lamC(s2))) {
            invdl(s2, t) = 1.0 / dl;
          } else {
            invdl(s2, t) = 0.0;
            if (s2 != t) {
              degenerate.emplace_back(s2, t);
            }
          }
        }
      }
      const Eigen::VectorXd epole = eP.segment(r0, nr);
      Eigen::ArrayXd s_n = -Eigen::ArrayXd::Ones(nr);
      s_n.head(nk).setOnes();
      const Eigen::MatrixXd Vxc_mo =
          C.transpose() * vxc.IntegrateVXC(D0).matrix() * C;

      // per QP level: a_n = F^T k_n stored (nr x S); B_pn, k_n and the
      // pole terms h, h', phi_n(lambda_s), phi_n'(lambda_s) at omega_p are
      // recomputed when needed
      auto gatherB = [&](Index p) {
        Eigen::MatrixXd Bp(nr, naux);
        for (Index n = 0; n < nr; ++n) {
          Bp.row(n) = Bmat.row((q0 + p) + nw * (r0 + n));
        }
        return Bp;
      };
      struct PoleTerms {
        Eigen::ArrayXXd phi, dphi, hp;
      };
      auto poleTerms = [&](Index p) {
        const double omega_p = in.qp_energies(p);
        PoleTerms pt;
        pt.phi.resize(nr, S);
        pt.dphi.resize(nr, S);
        pt.hp.resize(nr, S);
        for (Index s2 = 0; s2 < S; ++s2) {
          const double om = OmC(s2);
          for (Index n = 0; n < nr; ++n) {
            const double t = omega_p - epole(n) + s_n(n) * om;
            const double den = t * t + seta2;
            const double h = t / den;
            const double hp = (seta2 - t * t) / (den * den);
            pt.phi(n, s2) = h / om;
            pt.hp(n, s2) = hp / om;  // h' / Omega
            pt.dphi(n, s2) =
                -h / (2.0 * om * om * om) + s_n(n) * hp / (2.0 * om * om);
          }
        }
        return pt;
      };
      std::vector<Eigen::MatrixXd> Aps(static_cast<std::size_t>(nq));
      Eigen::VectorXd dsigma(nq);
      Eigen::MatrixXd poles(nr, nq);
      for (Index p = 0; p < nq; ++p) {
        const Index pabs = in.qpmin + p;
        const double omega_p = in.qp_energies(p);
        const Eigen::MatrixXd Bp = gatherB(p);
        Aps[std::size_t(p)] =
            ((Bp * Bkd.transpose()) * sD.asDiagonal()) * F;  // a_n = F^T k_n
        const PoleTerms pt = poleTerms(p);
        const Eigen::ArrayXXd A2 = Aps[std::size_t(p)].array().square();
        const double sigma_c = 2.0 * (A2 * pt.phi).sum();
        dsigma(p) = 2.0 * (A2 * pt.hp).sum();
        poles.col(p) = -2.0 * (A2 * pt.hp).rowwise().sum().matrix();
        const double sigma_x = -(1.0 - ax) * Bp.topRows(nk).squaredNorm();
        const double e_rebuilt =
            eps(pabs) + sigma_x - Vxc_mo(pabs, pabs) + sigma_c;
        result.qp_rebuild_error =
            std::max(result.qp_rebuild_error, std::abs(e_rebuilt - omega_p));
      }

      auto reverse = [&](const Eigen::VectorXd& vq, bool accumulate) {
        Eigen::VectorXd rho = Eigen::VectorXd::Zero(nr);
        Eigen::MatrixXd Gt = Eigen::MatrixXd::Zero(S, S);     // eigenbasis of C
        Eigen::MatrixXd gB = Eigen::MatrixXd::Zero(S, naux);  // d/db_kd
        Eigen::VectorXd gDelta = Eigen::VectorXd::Zero(S);
        for (Index p = 0; p < nq; ++p) {
          const double v = vq(p);
          if (v == 0.0) {
            continue;
          }
          const Eigen::MatrixXd Bp = gatherB(p);
          const Eigen::MatrixXd Kraw = Bp * Bkd.transpose();  // (n, kd)
          const Eigen::MatrixXd& A = Aps[std::size_t(p)];
          const PoleTerms pt = poleTerms(p);
          const Eigen::ArrayXXd A2 = A.array().square();
          // d/dC through phi_n(C): divided differences in the eigenbasis
          const Eigen::MatrixXd Aphi = (A.array() * pt.phi).matrix();
          Gt.noalias() +=
              (2.0 * v) * ((Aphi.transpose() * A - A.transpose() * Aphi)
                               .cwiseProduct(invdl));
          Gt.diagonal() += (2.0 * v) * (A2 * pt.dphi).colwise().sum().matrix();
          for (const auto& st : degenerate) {
            Gt(st.first, st.second) +=
                2.0 * v *
                (A.col(st.first).array() * A.col(st.second).array() *
                 pt.dphi.col(st.first))
                    .sum();
          }
          // pole energies e_n
          rho += v * poles.col(p);
          // k_n = Kraw_n * Delta^1/2, a_n = F^T k_n
          const Eigen::MatrixXd gK =
              ((4.0 * v) * (A.array() * pt.phi)).matrix() * F.transpose();
          gDelta += (gK.array() * Kraw.array())
                        .colwise()
                        .sum()
                        .matrix()
                        .transpose()
                        .cwiseQuotient(2.0 * sD);
          const Eigen::MatrixXd gKraw = gK * sD.asDiagonal();  // d/d(B_pn.b)
          gB.noalias() += gKraw.transpose() * Bp;
          if (accumulate) {
            w(q0 + p) += v;                           // eps_p
            u(in.qpmin + p) -= v;                     // -v_xc
            const Eigen::MatrixXd gBp = gKraw * Bkd;  // d/dB_pn
            for (Index n = 0; n < nr; ++n) {
              for (Index P = 0; P < naux; ++P) {
                double gpn = gBp(n, P);
                if (n < nk) {  // Sigma_x
                  gpn += -2.0 * (1.0 - ax) * v * Bp(n, P);
                }
                g(P)(q0 + p, r0 + n) += gpn;
              }
            }
          }
        }
        // C = Delta^1/2 (Delta + 4 b b^T) Delta^1/2
        Eigen::MatrixXd GC = F * Gt * F.transpose();
        GC = (0.5 * (GC + GC.transpose())).eval();
        gDelta += (GC.cwiseProduct(Cm)).rowwise().sum().cwiseQuotient(Delta) +
                  GC.diagonal().cwiseProduct(Delta);
        gB.noalias() += 8.0 * (sD.asDiagonal() * GC * sD.asDiagonal()) * Bkd;
        if (accumulate) {
          for (Index P = 0; P < naux; ++P) {
            for (Index d = 0; d < nd; ++d) {
              g(P).block(k0, c0 + d, nk, 1) += gB.col(P).segment(nk * d, nk);
            }
          }
        }
        for (Index d = 0; d < nd; ++d) {
          for (Index k = 0; k < nk; ++k) {
            const double t = gDelta(k + nk * d);
            rho(nk + d) += t;
            rho(k) -= t;
          }
        }
        return rho;
      };

      finish(dsigma, reverse, poles);
    } else {
      throw std::runtime_error(
          "BSEExcitationGradient: QP response only for the ppm and exact "
          "self-energies");
    }
  }

  phase("QP-correction response");
  // ---- integral part of sum g . dB
  // metric term first (needs g and B), then B is released and g becomes
  // gL (column R: sum_P L_RP g^P) in place; the integral derivatives and
  // the orbital part Q (with s^R = gL^R + gL^R^T) both use gL
  Eigen::MatrixXd grad;
  {
    Eigen::MatrixXd G = gmat.transpose() * Bmat;
    Bmat.resize(0, 0);
    G = (0.5 * (G + G.transpose())).eval();
    grad = -0.5 * ContractCoulombMetricDerivatives(auxbasis_, L * G * L);
  }
  RightMultiplyInPlace(gmat, L);
  const Eigen::MatrixXd& gL = gmat;
  grad += ContractThreeCenterDerivatives(
      auxbasis_, dftbasis_, [&gL, &Cw, nw](Index R, Eigen::MatrixXd& Gamma) {
        Gamma =
            Cw * (Eigen::Map<const Eigen::MatrixXd>(gL.col(R).data(), nw, nw) *
                  Cw.transpose());
      });

  phase("integral derivatives of B");
  // ---- orbital part: Q_rp = sum_{P,q} B^P_rq s^P_qp, s^P = g^P + g^P^T,
  // for all r (all MOs), p in the window
  Eigen::MatrixXd Q = Eigen::MatrixXd::Zero(nmo, nmo);
  {
    Eigen::MatrixXd Qw = Eigen::MatrixXd::Zero(nmo, nw);
    ForEachThreeCenterIntegral(
        auxbasis_, dftbasis_, [&](Index R, const Eigen::MatrixXd& I) {
          const Eigen::Map<const Eigen::MatrixXd> gLR(gL.col(R).data(), nw, nw);
          const Eigen::MatrixXd sLR = gLR + gLR.transpose();
          const Eigen::MatrixXd t = C.transpose() * (I * (Cw * sLR));
#pragma omp critical(bsegradient_Q)
          { Qw += t; }
        });
    Q.middleCols(lo, nw) = Qw;
  }
  gmat.resize(0, 0);
  // v_xc expectation values (QP response): orbital part 2 u_p Vxc_rp
  const bool has_u = u.cwiseAbs().maxCoeff() > 0.0;
  if (has_u) {
    const Eigen::MatrixXd Vxc_mo =
        C.transpose() * vxc.IntegrateVXC(D0).matrix() * C;
    Q += 2.0 * Vxc_mo * u.asDiagonal();
  }
  Eigen::VectorXd wfull = Eigen::VectorXd::Zero(nmo);
  wfull.segment(lo, nw) = w;

  phase("orbital part Q");
  // ---- D, Theta (same-space blocks and diagonal)
  const double kDegenerate = 1e-8;
  Eigen::MatrixXd Dmo = Eigen::MatrixXd::Zero(nmo, nmo);
  Eigen::MatrixXd Wmo = Eigen::MatrixXd::Zero(nmo, nmo);
  auto same_space = [&](Index first, Index n) {
    for (Index r = first; r < first + n; ++r) {
      for (Index p = first; p < first + n; ++p) {
        if (r == p) {
          Dmo(p, p) = wfull(p);
          Wmo(p, p) = 0.5 * Q(p, p) + wfull(p) * eps(p);
          continue;
        }
        const double de = eps(p) - eps(r);
        if (std::abs(de) < kDegenerate) {
          Wmo(r, p) = 0.25 * (Q(r, p) + Q(p, r));
          continue;
        }
        Dmo(r, p) = (Q(r, p) - Q(p, r)) / (2.0 * de);
        Wmo(r, p) = (Q(r, p) * eps(p) - Q(p, r) * eps(r)) / (2.0 * de);
      }
    }
  };
  same_space(0, nocc);
  same_space(nocc, nvirt);

  // closed-shell Fock response G[X] (MO in, MO out) to a symmetric MO
  // density change X
  auto fock_response = [&](const Eigen::MatrixXd& Xmo) {
    const Eigen::MatrixXd Xao = C * Xmo * C.transpose();
    Eigen::MatrixXd G = vxc.IntegrateFXC(D0, {Xao})[0];
    G += ERIs_.CalculateERIs_3c(Xao);
    if (ScaHFX_ > 0) {
      G += 0.5 * ScaHFX_ * ERIs_.CalculateEXX_3c(Eigen::MatrixXd(0, 0), Xao);
    }
    return Eigen::MatrixXd(C.transpose() * G * C);
  };

  // XC-kernel response alone (the v_xc expectation values), MO in and out
  Eigen::MatrixXd Gu;
  if (has_u) {
    const Eigen::MatrixXd Uao = C * u.asDiagonal() * C.transpose();
    Gu = C.transpose() * vxc.IntegrateFXC(D0, {Uao})[0] * C;
  }

  phase("Fock responses");
  // ---- Z-vector: (A + B) Z = L
  Eigen::MatrixXd Zvo;
  {
    const Eigen::MatrixXd GD = fock_response(Dmo);
    Eigen::MatrixXd Lvo = Q.block(nocc, 0, nvirt, nocc) -
                          Q.block(0, nocc, nocc, nvirt).transpose() +
                          4.0 * GD.block(nocc, 0, nvirt, nocc);
    if (has_u) {
      Lvo += 4.0 * Gu.block(nocc, 0, nvirt, nocc);
    }
    // StaticResponse solves (A + B) U = -C_v^T h C_o
    AOOverlap overlap;
    overlap.Fill(dftbasis_);
    const Eigen::MatrixXd& S = overlap.Matrix();
    Eigen::MatrixXd h = -(S * C_virt * Lvo * C_occ.transpose() * S);
    h = (h + h.transpose()).eval();
    // the gradient error is linear in the Z residual; 1e-8 relative is
    // well below the other errors (1e-10 cost four more iterations)
    std::vector<Eigen::MatrixXd> dP = StaticResponse(orb, {h}, 1e-8, 200);
    Zvo = 0.5 * C_virt.transpose() * S * dP[0] * S * C_occ;
  }

  phase("Z-vector");
  // ---- relaxed density Ptilde and energy-weighted density W
  Eigen::MatrixXd Pmo = Dmo;
  Pmo.block(nocc, 0, nvirt, nocc) = -0.5 * Zvo;
  Pmo.block(0, nocc, nocc, nvirt) = -0.5 * Zvo.transpose();
  {
    const Eigen::MatrixXd GP = fock_response(Pmo);
    Wmo.topLeftCorner(nocc, nocc) += 2.0 * GP.topLeftCorner(nocc, nocc);
    if (has_u) {
      Wmo.topLeftCorner(nocc, nocc) += 2.0 * Gu.topLeftCorner(nocc, nocc);
    }
  }
  for (Index i = 0; i < nocc; ++i) {
    for (Index a = 0; a < nvirt; ++a) {
      const double wai = 0.5 * (Q(i, nocc + a) - Zvo(a, i) * eps(i));
      Wmo(nocc + a, i) = wai;
      Wmo(i, nocc + a) = wai;
    }
  }
  const Eigen::MatrixXd Pao = C * Pmo * C.transpose();
  const Eigen::MatrixXd Wao = C * Wmo * C.transpose();

  phase("energy-weighted density");
  // ---- sum Ptilde dF/dR|_D - sum W dS/dR
  grad += ContractAODerivatives(ComputeKineticDerivatives(dftbasis_), Pao);
  grad += ContractAODerivatives(
      ComputeNuclearAttractionDerivatives(dftbasis_, mol), Pao);
  grad -= ContractAODerivatives(ComputeOverlapDerivatives(dftbasis_), Wao);
  phase("one-electron and overlap derivatives");
  // RI-J, exact: the Hartree energy is quadratic in the density
  grad += 0.5 * (DFTGradient::RIJGradient(D0 + Pao, auxbasis_, dftbasis_) -
                 DFTGradient::RIJGradient(D0 - Pao, auxbasis_, dftbasis_));
  if (ScaHFX_ > 0) {
    const SymmetricEigenSystem es = SymmetricEigen(Pao);
    const double cutoff = 1e-12 * es.values.cwiseAbs().maxCoeff();
    std::vector<Index> keep;
    for (Index k = 0; k < es.values.size(); ++k) {
      if (std::abs(es.values(k)) > cutoff) {
        keep.push_back(k);
      }
    }
    Eigen::MatrixXd U(es.vectors.rows(), Index(keep.size()));
    Eigen::VectorXd lambda(Index(keep.size()));
    for (Index k = 0; k < Index(keep.size()); ++k) {
      U.col(k) = es.vectors.col(keep[std::size_t(k)]);
      lambda(k) = es.values(keep[std::size_t(k)]);
    }
    grad -= ScaHFX_ *
            DFTGradient::RIPairGradient(U, lambda, C_occ, auxbasis_, dftbasis_);
  }
  phase("RI-J and RI-K");
  // XC: d/dt of the ground-state XC gradient at D0 + t (Ptilde + U), with
  // U = sum_p u_p |p><p| from the v_xc expectation values (analytic, one
  // grid pass)
  {
    const Eigen::MatrixXd Pxc = Pao + C * u.asDiagonal() * C.transpose();
    grad += vxc.GradientResponse(D0, Pxc, mol);
  }
  phase("XC");
  {
    double total = 0.0;
    for (const auto& pt : phase_times) {
      total += pt.second;
    }
    XTP_LOG(Log::error, *pLog_)
        << TimeStamp() << " BSE gradient timing (s), total " << total
        << std::flush;
    for (const auto& pt : phase_times) {
      XTP_LOG(Log::error, *pLog_)
          << "   " << pt.first << ": " << pt.second << std::flush;
    }
  }
  result.gradient = grad;
  return result;
}

DFTEngine::BSEGradientInput DFTEngine::BSEGradientInputFromOrbitals(
    const Orbitals& orb, bool singlet, Index index) {
  if (orb.hasBetaMOs()) {
    throw std::runtime_error(
        "BSE gradient: closed-shell (RKS) GW-BSE runs only");
  }
  if (orb.GetFlagUseHqpOffdiag()) {
    throw std::runtime_error(
        "BSE gradient: needs the diagonal QP Hamiltonian "
        "(bse.use_Hqp_offdiag false)");
  }
  if (orb.isQSGW()) {
    throw std::runtime_error("BSE gradient: QSGW is not supported");
  }
  const tools::EigenSystem& es =
      singlet ? orb.BSESinglets() : orb.BSETriplets();
  if (index < 0 || index >= es.eigenvalues().size()) {
    throw std::runtime_error("BSE gradient: state " + std::to_string(index) +
                             " was not computed");
  }
  BSEGradientInput in;
  in.singlet = singlet;
  in.vmin = orb.getBSEvmin();
  in.cmax = orb.getBSEcmax();
  in.rpamin = orb.getRPAmin();
  in.rpamax = orb.getRPAmax();
  in.X = es.eigenvectors().col(index);
  if (!orb.getTDAApprox()) {
    in.Y = es.eigenvectors2().col(index);
  }

  const Eigen::VectorXd& rpa = orb.RPAInputEnergies();
  if (rpa.size() != in.rpamax - in.rpamin + 1) {
    throw std::runtime_error("BSE gradient: no RPA input energies in orb");
  }
  in.energies_rpa = rpa;
  // diagonal of Hqp = V diag(e) V^T (QPdiag), as GWBSE builds it
  const Eigen::MatrixXd& V = orb.QPdiag().eigenvectors();
  const Eigen::VectorXd hqp_diag =
      V.cwiseProduct(V) * orb.QPdiag().eigenvalues();
  const Index qpmin = orb.getGWAmin();
  const Index qpmax = orb.getGWAmax();
  if (hqp_diag.size() != qpmax - qpmin + 1) {
    throw std::runtime_error("BSE gradient: no QP Hamiltonian in orb");
  }
  in.energies_bse.resize(in.cmax - in.vmin + 1);
  for (Index p = in.vmin; p <= in.cmax; ++p) {
    in.energies_bse(p - in.vmin) =
        (p >= qpmin && p <= qpmax) ? hqp_diag(p - qpmin) : rpa(p - in.rpamin);
  }
  return in;
}

}  // namespace xtp
}  // namespace votca
