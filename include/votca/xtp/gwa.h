/*
 *            Copyright 2009-2017 The VOTCA Development Team
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

#ifndef _VOTCA_XTP_GWA_H
#define _VOTCA_XTP_GWA_H
#include <votca/xtp/votca_config.h>

namespace votca {
namespace xtp {


class GWA {
 public:
  GWA(Orbitals* orbitals)
      : _qp_diag_energies(orbitals->QPdiagEnergies()),
        _qp_diag_coefficients(orbitals->QPdiagCoefficients()){};
        


 private:
 
  // program tasks
  bool _do_qp_diag;
 

  // storage tasks
  bool _store_qp_pert;
  bool _store_qp_diag;
  

  // iterate G and W and not only G
  bool _iterate_gw;

 
  int _openmp_threads;

 


 
  unsigned int _homo;   // HOMO index
  unsigned int _rpamin;
  unsigned int _rpamax;
  double _rpamaxfactor;  // RPA level range
  unsigned int _qpmin;
  unsigned int _qpmax;
  unsigned int _qptotal;
  double _qpminfactor;
  double _qpmaxfactor;  // QP level range
  

  double _g_sc_limit;  // convergence criteria for g iteration [Hartree]]
  unsigned int _g_sc_max_iterations;
  unsigned int _gw_sc_max_iterations;
  double _gw_sc_limit;  // convergence criteria for gw iteration [Hartree]]
  
  double _shift;  // pre-shift of DFT energies
 

  

 
  // Sigma related variables and functions
  Eigen::MatrixXd _sigma_x;  // exchange term
  Eigen::MatrixXd _sigma_c;  // correlation term

  void sigma_prepare_threecenters(TCMatrix_gwbse& _Mmn);

  void sigma_diag(const TCMatrix_gwbse& _Mmn);
  void sigma_offdiag(const TCMatrix_gwbse& _Mmn);

  // QP variables and functions
  Eigen::VectorXd _qp_energies;
  Eigen::VectorXd& _qp_diag_energies;      // stored in orbitals object
  Eigen::MatrixXd& _qp_diag_coefficients;  // dito
  void FullQPHamiltonian();

};
}
}

#endif /* _VOTCA_XTP_GWA_H */
