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

#pragma once
#ifndef VOTCA_XTP_VXC_POTENTIAL_H
#define VOTCA_XTP_VXC_POTENTIAL_H

#include <xc.h>

#include <array>
#include <vector>

#include "grid_containers.h"
#include "gridbox.h"

#undef LOG

namespace votca {
namespace xtp {

class AOBasis;
class QMMolecule;

template <class Grid>
class Vxc_Potential {
 public:
  struct SpinResult {
    double energy = 0.0;
    Eigen::MatrixXd vxc_alpha;
    Eigen::MatrixXd vxc_beta;
  };

  explicit Vxc_Potential(const Grid& grid) : grid_(grid) {};
  ~Vxc_Potential();

  static double getExactExchange(const std::string& functional);
  void setXCfunctional(const std::string& functional);

  Mat_p_Energy IntegrateVXC(const Eigen::MatrixXd& density_matrix) const;

  /// The xc kernel of the closed-shell density P (density_matrix, total)
  /// applied to symmetric first-order density matrices dP_k (total): the
  /// first-order change of IntegrateVXC, int f_xc[rho] drho_k chi chi, with
  /// the gradient terms for GGAs. The AO values are evaluated once for all
  /// k. Throws if the functional has no second derivatives in libxc.
  std::vector<Eigen::MatrixXd> IntegrateFXC(
      const Eigen::MatrixXd& density_matrix,
      const std::vector<Eigen::MatrixXd>& response_matrices) const;
  SpinResult IntegrateVXCSpin(const Eigen::MatrixXd& dmat_alpha,
                              const Eigen::MatrixXd& dmat_beta) const;

  /// Spin-polarized analogue of IntegrateFXC: the first-order changes
  /// {dV_alpha, dV_beta} of IntegrateVXCSpin at the spin densities
  /// (dmat_alpha, dmat_beta) for symmetric first-order spin densities
  /// response[k] = {dP_alpha, dP_beta}, from the polarized second
  /// derivatives of libxc (LDA/GGA). The AO values are evaluated once for
  /// all k. Throws if the functional has no second derivatives in libxc.
  std::vector<std::array<Eigen::MatrixXd, 2>> IntegrateFXCSpin(
      const Eigen::MatrixXd& dmat_alpha, const Eigen::MatrixXd& dmat_beta,
      const std::vector<std::array<Eigen::MatrixXd, 2>>& response) const;

  // ===========================================================================
  // XC nuclear gradient at fixed density matrix, part 1: basis-function
  // (Pulay) and grid-point translation terms, for LDA and GGA (the GGA sigma
  // terms use AOShell::EvalAOspaceHessian and vanish for LDA, where
  // df_dsigma stays 0). Despite the name it computes both terms, which
  // share one AO evaluation per grid point. Derivation in vxc_potential.cc.
  // ===========================================================================
  Eigen::MatrixXd PulayGradient(const Eigen::MatrixXd& density_matrix,
                                const AOBasis& dftbasis) const;

  // ===========================================================================
  // XC nuclear gradient, part 2: derivative of the SSW grid weights. Only the
  // sum PulayGradient + GridWeightGradient equals the derivative of the XC
  // energy of IntegrateVXC (checked by finite differences in
  // test_xcgradient.cc). Cost is O(Natoms^3) per grid point, one power of
  // Natoms more than the SSW partition itself.
  // ===========================================================================
  Eigen::MatrixXd GridWeightGradient(const Eigen::MatrixXd& density_matrix,
                                     const QMMolecule& atoms) const;

  // ===========================================================================
  // UKS (spin-polarized) analogs of PulayGradient/GridWeightGradient
  // above, for LDA and GGA.
  //
  // Basis-type and translation-type rho terms, for each spin channel
  // s in {a,b}:
  //   d(rho_s,p)/dR_A|_basis = -2 * sum_{mu in A} temp_s,mu * grad(chi_mu)
  //   d(rho_s,p)/dR_A|_translation = grad_s,p   if A == owner(p), else 0
  // with temp_s = dmat_s * ao.values and grad_s = 2*(ao.derivatives^T *
  // temp_s), as in IntegrateVXCSpin (UKS spin density matrices are not
  // doubled, unlike DMAT_here = 2*P for RKS), contracted with vrho_a/vrho_b
  // and summed over both spins.
  //
  // GGA: the three sigma_aa/sigma_ab/sigma_bb terms reduce to the per-spin
  // machinery of the restricted GGA case (Gmat_s, Hessian_rho_s) applied
  // once per spin, with the effective gradient vectors
  //   V_a = 2*vsigma_aa*rho_a_grad + vsigma_ab*rho_b_grad
  //   V_b = 2*vsigma_bb*rho_b_grad + vsigma_ab*rho_a_grad
  // Basis-type sigma contribution per mu in A:
  //   -[grad(chi_mu)*(Gmat_a.V_a)_mu + temp_a,mu*(Hessian_mu.V_a)
  //     + grad(chi_mu)*(Gmat_b.V_b)_mu + temp_b,mu*(Hessian_mu.V_b)]
  // Translation-type sigma contribution (A==owner(p) only):
  //   Hessian_rho_a . V_a + Hessian_rho_b . V_b
  // where Gmat_s = 2*dmat_s*ao.derivatives and Hessian_rho_s is the
  // restricted Hessian_rho per spin channel.
  //
  // The weight-derivative term uses the same dw/dR as GridWeightGradient;
  // only the energy density changes to (rho_a+rho_b)*f_xc, evaluated with
  // sigma_aa/sigma_ab/sigma_bb from the density gradients, as in the energy
  // of IntegrateVXCSpin.
  Eigen::MatrixXd PulayGradientUKS(const Eigen::MatrixXd& dmat_alpha,
                                   const Eigen::MatrixXd& dmat_beta,
                                   const AOBasis& dftbasis) const;
  Eigen::MatrixXd GridWeightGradientUKS(const Eigen::MatrixXd& dmat_alpha,
                                        const Eigen::MatrixXd& dmat_beta,
                                        const QMMolecule& atoms) const;
  // ===========================================================================

 private:
  struct XC_entry {
    double f_xc = 0;
    double df_drho = 0;
    double df_dsigma = 0;
  };

  struct XC_entry_spin {
    double f_xc = 0;
    double vrho_a = 0;
    double vrho_b = 0;
    double vsigma_aa = 0;
    double vsigma_ab = 0;
    double vsigma_bb = 0;
  };

  XC_entry EvaluateXC(double rho, double sigma) const;
  XC_entry_spin EvaluateXCSpin(double rho_a, double rho_b, double sigma_aa,
                               double sigma_ab, double sigma_bb) const;

  // libxc for n points at once, exchange and correlation summed. Unpolarized:
  // rho, sigma, exc, vrho, vsigma have n entries. Polarized (libxc layout):
  // rho and vrho 2n, sigma and vsigma 3n, exc n. vsigma is zero for LDA.
  void EvaluateXCBatch(Index n, const double* rho, const double* sigma,
                       double* exc, double* vrho, double* vsigma) const;
  void EvaluateXCSpinBatch(Index n, const double* rho, const double* sigma,
                           double* exc, double* vrho, double* vsigma) const;
  // first and second derivatives for n points, unpolarized, exchange and
  // correlation summed; vsigma, v2rhosigma, v2sigma2 are zero for LDA
  void EvaluateFXCBatch(Index n, const double* rho, const double* sigma,
                        double* vsigma, double* v2rho2, double* v2rhosigma,
                        double* v2sigma2) const;
  // the same, polarized (libxc layout): rho 2n, sigma and vsigma 3n,
  // v2rho2 3n (aa, ab, bb), v2rhosigma 6n (a_aa, a_ab, a_bb, b_aa, b_ab,
  // b_bb), v2sigma2 6n (aa_aa, aa_ab, aa_bb, ab_ab, ab_bb, bb_bb)
  void EvaluateFXCSpinBatch(Index n, const double* rho, const double* sigma,
                            double* vsigma, double* v2rho2, double* v2rhosigma,
                            double* v2sigma2) const;

  // Grid points processed together inside a box.
  static constexpr Index kPointBlock = 128;

  const Grid grid_;
  int xfunc_id;
  bool setXC_ = false;
  bool use_separate_;
  int cfunc_id;
  xc_func_type xfunc;
  xc_func_type cfunc;
  // Polarized handles for the UKS path, set up once with the unpolarized
  // ones instead of per grid point.
  xc_func_type xfunc_pol_;
  xc_func_type cfunc_pol_;
  bool is_gga_ = false;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_VXC_POTENTIAL_H
