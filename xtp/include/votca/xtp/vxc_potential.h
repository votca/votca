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
  SpinResult IntegrateVXCSpin(const Eigen::MatrixXd& dmat_alpha,
                              const Eigen::MatrixXd& dmat_beta) const;

  // ===========================================================================
  // STATUS: originally basis-function/Pulay term only; a second, distinct
  // term (grid-point translation) was added later after the
  // GridWeightGradient C_p fix substantially improved but did not fully
  // resolve a residual discrepancy against finite differences -- see
  // vxc_potential.cc for the full derivation of both terms. The name
  // "PulayGradient" is now a slight misnomer (it computes basis-function
  // + grid-point-translation contributions together, since both loop
  // over the same grid points and share the same expensive AO
  // evaluation) -- kept for now rather than renaming mid-branch, but
  // worth reconsidering once the whole XC gradient is confirmed correct.
  //
  // Valid in FULL for LDA functionals for the basis-function part; for
  // GGA functionals it captures only the df_drho-driven part of that
  // term, NOT the additional df_dsigma-driven Pulay contribution, which
  // needs second derivatives of basis functions w.r.t. electron position
  // (not yet available anywhere in this codebase). The grid-point-
  // translation term added later has the same LDA-only caveat (uses
  // df_drho only, not df_dsigma). The grid-weight (SSW partition)
  // derivative term is a separate, third piece, computed in
  // GridWeightGradient below.
  // ===========================================================================
  Eigen::MatrixXd PulayGradient(const Eigen::MatrixXd& density_matrix,
                                const AOBasis& dftbasis) const;

  // ===========================================================================
  // STATUS: written but NOT yet run/tested. This is the second (and
  // harder) of the two pieces needed for the XC gradient, alongside
  // PulayGradient above -- the SSW grid-weight nuclear derivative.
  //
  // The underlying formula was derived and independently verified
  // NUMERICALLY in Python (two separate multi-atom test configurations,
  // matching finite differences to ~1e-9/1e-10) before writing any of
  // this C++ -- see conversation history for the full derivation. This
  // C++ is a direct translation of that verified Python, not a fresh
  // derivation done blind.
  //
  // Combined with PulayGradient, this should give the complete LDA-level
  // XC gradient (GGA's additional df_dsigma-driven Pulay term is still
  // out of scope, per the note on PulayGradient). Validate the SUM of
  // the two against a finite difference of the real total XC energy
  // (IntegrateVXC's energy output) -- neither piece alone can be checked
  // against the total energy in isolation.
  //
  // SCALING NOTE: this is O(Natoms^3) per grid point (derivative w.r.t.
  // every atom, each needing a sum over every other atom's pairwise
  // term) -- consistent with the SSW partition's inherent O(Natoms^2)
  // energy-level cost, but one power of Natoms worse for the gradient.
  // Fine for validation-scale systems; would need optimizing before use
  // on anything large.
  // ===========================================================================
  Eigen::MatrixXd GridWeightGradient(const Eigen::MatrixXd& density_matrix,
                                     const QMMolecule& atoms) const;

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

  const Grid grid_;
  int xfunc_id;
  bool setXC_ = false;
  bool use_separate_;
  int cfunc_id;
  xc_func_type xfunc;
  xc_func_type cfunc;
};

}  // namespace xtp
}  // namespace votca
#endif  // VOTCA_XTP_VXC_POTENTIAL_H
