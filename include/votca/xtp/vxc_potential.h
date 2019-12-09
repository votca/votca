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

#pragma once
#ifndef XTP_VXC_POTENTIAL_H
#define XTP_VXC_POTENTIAL_H

#include <votca/xtp/grid_containers.h>
#include <votca/xtp/gridbox.h>

#include <xc.h>
#undef LOG

namespace votca {
namespace xtp {

template <class Grid>
class Vxc_Potential {
 public:
  explicit Vxc_Potential(const Grid& grid) : _grid(grid){};
  ~Vxc_Potential();

  static double getExactExchange(const std::string& functional);
  void setXCfunctional(const std::string& functional);
  Mat_p_Energy IntegrateVXC(const Eigen::MatrixXd& density_matrix) const;

 private:
  struct XC_entry {
    double f_xc = 0;  // E_xc[n] = int{n(r)*eps_xc[n(r)] d3r} = int{ f_xc(r) d3r
    double df_drho = 0;    // v_xc_rho(r) = df/drho
    double df_dsigma = 0;  // df/dsigma ( df/dgrad(rho) = df/dsigma *
                           // dsigma/dgrad(rho) = df/dsigma * 2*grad(rho))
  };

  XC_entry EvaluateXC(double rho, double sigma) const;

  const Grid _grid;
  int xfunc_id;
  bool _setXC = false;
  bool _use_separate;
  int cfunc_id;
  xc_func_type xfunc;  // handle for exchange functional
  xc_func_type cfunc;  // handle for correlation functional
};

}  // namespace xtp
}  // namespace votca
#endif  // XTP_NUMERICAL_INTEGRATION_H
