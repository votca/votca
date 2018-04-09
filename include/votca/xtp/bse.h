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
#include <votca/xtp/orbitals.h>
#include <votca/xtp/threecenter.h>



namespace votca {
namespace xtp {




class BSE {
 public:
  BSE(Orbitals* orbitals):
        _orbitals(orbitals),
        _eh_x(orbitals->eh_x()),
        _eh_d(orbitals->eh_d()),
        _bse_singlet_energies(orbitals->BSESingletEnergies()),
        _bse_singlet_coefficients(orbitals->BSESingletCoefficients()),
        _bse_singlet_coefficients_AR(orbitals->BSESingletCoefficientsAR()),
        _bse_triplet_energies(orbitals->BSETripletEnergies()),
        _bse_triplet_coefficients(orbitals->BSETripletCoefficients()){};

  ~BSE(){};
  
  void setBSEindices(int vmin, int vmax, int cmin, int cmax, int nmax) {
                _bse_vmin = vmin;
                _bse_vmax = vmax;
                _bse_cmin = cmin;
                _bse_cmax = cmax;
                _bse_nmax = nmax;
                _bse_vtotal = _bse_vmax - _bse_vmin + 1;
                _bse_ctotal = _bse_cmax - _bse_cmin + 1;
                _bse_size = _bse_vtotal * _bse_ctotal;
                return;
            }

  void Setup_Hx(TCMatrix_gwbse& _Mmn);
  void Setup_Hd(TCMatrix_gwbse& _Mmn);
  void Setup_Hd_BTDA(TCMatrix_gwbse& _Mmn);
  void Add_HqpToHd(const Eigen::MatrixXd& Hqp );
  
  void Solve_triplets();
  void Solve_singlets();
  void Solve_singlets_BTDA();
  void Analyze_triplets();
  void Analyze_singlets();
  
  void FreeMatrices(){
      _eh_d.resize(0, 0);
      _eh_x.resize(0, 0);
  }
 

 private:
 
  
  unsigned  _bse_vmin;
  unsigned  _bse_vmax;
  unsigned  _bse_cmin;
  unsigned  _bse_cmax;
  unsigned  _bse_size;
  unsigned  _bse_vtotal;
  unsigned  _bse_ctotal;
  int _bse_nmax;
  int _bse_nprint;
  double _min_print_weight;

  Orbitals* _orbitals;
  

  // BSE variables and functions
  MatrixXfd& _eh_x;  // stored in orbitals object
  MatrixXfd& _eh_d;  // stored in orbitals object
  MatrixXfd _eh_d2;  // because it is not stored in orbitals object

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
  

  void Solve_nonhermitian(Eigen::MatrixXd& H, Eigen::MatrixXd& L);
  std::vector<int> _index2v;
  std::vector<int> _index2c;

  // some cleaner analysis
  
 

  void Analyze_eh_interaction_Triplet(std::vector<real_gwbse>& _c_d,
                                          std::vector<real_gwbse>& _c_qp);
  
  void Analyze_eh_interaction_Singlet(std::vector<real_gwbse>& _c_x,
                                               std::vector<real_gwbse>& _c_d,
                                               std::vector<real_gwbse>& _c_qp);

  void FragmentPopulations(const string& spin,
                               std::vector<Eigen::VectorXd >& popH,
                               std::vector<Eigen::VectorXd >& popE,
                               std::vector<Eigen::VectorXd >& Crgs);

  void CalcFreeTransition_Dipoles();

  void CalcCoupledTransition_Dipoles();
};
}
}

#endif /* _VOTCA_XTP_BSE_H */
