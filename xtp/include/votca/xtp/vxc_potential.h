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
