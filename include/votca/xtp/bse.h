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

#ifndef _VOTCA_XTP_BSE_H
#define _VOTCA_XTP_BSE_H
#include <votca/xtp/votca_config.h>

#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/threecenter.h>



namespace votca {
namespace xtp {




class BSE {
 public:
  BSE(Orbitals* orbitals):
        _eh_x(orbitals->eh_x()),
        _eh_d(orbitals->eh_d()),
        _bse_singlet_energies(orbitals->BSESingletEnergies()),
        _bse_singlet_coefficients(orbitals->BSESingletCoefficients()),
        _bse_singlet_coefficients_AR(orbitals->BSESingletCoefficientsAR()),
        _bse_triplet_energies(orbitals->BSETripletEnergies()),
        _bse_triplet_coefficients(orbitals->BSETripletCoefficients()){};

  ~BSE(){};

  

 

 private:
 

  // BSE variant
  bool _do_full_BSE;


  
  unsigned int _bse_vmin;
  unsigned int _bse_vmax;
  unsigned int _bse_cmin;
  unsigned int _bse_cmax;
  unsigned int _bse_size;
  unsigned int _bse_vtotal;
  unsigned int _bse_ctotal;
  int _bse_nmax;
  int _bse_nprint;
  double _min_print_weight;

  

  // BSE variables and functions
  MatrixXfd& _eh_x;  // stored in orbitals object
  MatrixXfd& _eh_d;  // stored in orbitals object
  MatrixXfd _eh_d2;  // because it is not stored in orbitals object
  MatrixXfd _eh_qp;

  VectorXfd& _bse_singlet_energies;  // stored in orbitals object
  MatrixXfd& _bse_singlet_coefficients;  // stored in orbitals
                                                      // object
  MatrixXfd& _bse_singlet_coefficients_AR;  // stored in orbitals
                                                         // object
  VectorXfd& _bse_triplet_energies;  // stored in orbitals object
  MatrixXfd& _bse_triplet_coefficients;  // stored in orbitals
                                                      // object

  std::vector<Eigen::MatrixXd > _interlevel_dipoles;
  std::vector<Eigen::MatrixXd > _interlevel_dipoles_electrical;
  void BSE_x_setup(TCMatrix_gwbse& _Mmn);
  void BSE_d_setup(TCMatrix_gwbse& _Mmn);
  void BSE_d2_setup(TCMatrix_gwbse& _Mmn);
  void BSE_qp_setup();
  void BSE_Add_qp2H(MatrixXfd& qp);
  void BSE_solve_triplets();
  void BSE_solve_singlets();
  void BSE_solve_singlets_BTDA();

  void Solve_nonhermitian(Eigen::MatrixXd& H, Eigen::MatrixXd& L);
  std::vector<int> _index2v;
  std::vector<int> _index2c;

  // some cleaner analysis
  void BSE_analyze_triplets();
  void BSE_analyze_singlets();
 

  void BSE_analyze_eh_interaction_Triplet(std::vector<real_gwbse>& _c_d,
                                          std::vector<real_gwbse>& _c_qp);
  
  void BSE_analyze_eh_interaction_Singlet(std::vector<real_gwbse>& _c_x,
                                               std::vector<real_gwbse>& _c_d,
                                               std::vector<real_gwbse>& _c_qp);

  void BSE_FragmentPopulations(const string& spin,
                               std::vector<Eigen::VectorXd >& popH,
                               std::vector<Eigen::VectorXd >& popE,
                               std::vector<Eigen::VectorXd >& Crgs);

  void BSE_FreeTransition_Dipoles();

  void BSE_CoupledTransition_Dipoles();
};
}
}

#endif /* _VOTCA_XTP_BSE_H */
