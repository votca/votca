/*
 *            Copyright 2009-2018 The VOTCA Development Team
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



#include <votca/xtp/rpa.h>
#include <votca/xtp/aomatrix.h>
#include "votca/xtp/threecenter.h"



namespace votca {
  namespace xtp {

 void RPA::calculate_epsilon(const Eigen::VectorXd& qp_energies,const TCMatrix_gwbse& _Mmn_full) {
const int _size = _Mmn_full.getAuxDimension(); // size of gwbasis
            for (auto& matrix : _epsilon_r) {
                matrix = Eigen::MatrixXd::Identity(_size,_size);
            }
            for (auto& matrix : _epsilon_i) {
                matrix = Eigen::MatrixXd::Identity(_size,_size);
            }

            int lumo=_homo+1;
            int n_occ=lumo-_rpamin;
            int n_unocc=_rpamax-_homo;
            
#pragma omp parallel for 
            for (int _m_level = 0; _m_level < n_occ; _m_level++) {
                const double _qp_energy_m = qp_energies(_m_level + _rpamin);
#if (GWBSE_DOUBLE)
                const Eigen::MatrixXd Mmn_RPA = _Mmn_full[ _m_level ].block(n_occ, 0,n_unocc, _size );
#else
                const Eigen::MatrixXd Mmn_RPA = _Mmn_full[ _m_level ].block(n_occ,0,  n_unocc, _size).cast<double>();       
#endif
                Eigen::MatrixXd tempresult=Eigen::MatrixXd::Zero(_size,_size);
                Eigen::MatrixXd denom_x_Mmn_RPA=Eigen::MatrixXd::Zero(n_unocc,_size);
                for (int i = 0; i < screen_freq_i.size(); ++i) {   
                    // a temporary matrix, that will get filled in empty levels loop
                    const double screen_freq2 = screen_freq_i(i) * screen_freq_i(i);
                    for (int _n_level = 0; _n_level < n_unocc; _n_level++) {
                        const double _deltaE = qp_energies(_n_level + lumo) - _qp_energy_m;
                        const double denom=4.0 * _deltaE / (_deltaE * _deltaE + screen_freq2);  
                        denom_x_Mmn_RPA.row(_n_level)=Mmn_RPA.row(_n_level)*denom; //hartree    
                    }
                    tempresult.noalias() = Mmn_RPA.transpose() * denom_x_Mmn_RPA;

#pragma omp critical
                    {
                        _epsilon_i[i] += tempresult;
                    }
                }

                //real parts
                for (int i = 0; i < screen_freq_r.size(); ++i) {
                    for (int _n_level = 0;  _n_level < n_unocc; _n_level++) {
                        const double _deltaE = qp_energies(_n_level + lumo) - _qp_energy_m;
                        const double denom=2.0 * (1.0 / (_deltaE - screen_freq_r(i)) + 1.0 / (_deltaE + screen_freq_r(i)));
                        denom_x_Mmn_RPA.row(_n_level)=Mmn_RPA.row(_n_level)*denom; //hartree    
                    }
                    tempresult.noalias() = Mmn_RPA.transpose() * denom_x_Mmn_RPA;

#pragma omp critical
                    {
                        _epsilon_r[i] += tempresult;
                    }
                }

            } // occupied levels


            return;
        }

   
  }
};
