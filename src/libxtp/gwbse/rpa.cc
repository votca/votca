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



#include <votca/xtp/rpa.h>
#include <votca/xtp/aomatrix.h>

#include "votca/xtp/threecenter.h"



namespace votca {
    namespace xtp {


                
        
        //imaginary
        Eigen::MatrixXd RPA::RPA_imaginary(const Eigen::VectorXd& qp_energies,const double screening_freq) {


            const int _size = _Mmn_RPA.getAuxDimension(); // size of gwbasis
            const int index_n = _Mmn_RPA.get_nmin();
            const int index_m = _Mmn_RPA.get_mmin();
            const double screenf2=screening_freq * screening_freq;
            Eigen::MatrixXd result=Eigen::MatrixXd::Zero(_size,_size);
            
            #pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {
                const double _qp_energy_m=qp_energies(_m_level + index_m);
#if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const Eigen::MatrixXd Mmn_RPA = _Mmn_RPA[ _m_level ].cast<double>();
#endif
                // a temporary matrix, that will get filled in empty levels loop
                Eigen::MatrixXd _temp = Eigen::MatrixXd(_Mmn_RPA.get_ntot(), _size);

                // loop over empty levels
                for (int _n_level = 0; _n_level < _Mmn_RPA.get_ntot(); _n_level++) {
                    

                    const double _deltaE = qp_energies(_n_level + index_n) -_qp_energy_m ; // get indices and units right!!!

                    // this only works, if we have either purely real or purely imaginary frequencies

                    // purely imaginary
                    const double _energy_factor = 4.0 * _deltaE / (_deltaE * _deltaE + screenf2);//hartree
                    for (int _i_gw = 0; _i_gw < _size; _i_gw++) {
                        _temp(_n_level, _i_gw) = _energy_factor * Mmn_RPA(_i_gw, _n_level);
                    } // matrix size

                } // empty levels
                _temp=Mmn_RPA*_temp;
                // now multiply and add to epsilon
                #pragma omp critical
                {
                result+=_temp;
                }
            } // occupied levels
              
            return result;
        }
        //real

        Eigen::MatrixXd RPA::RPA_real(const Eigen::VectorXd& qp_energies,const double screening_freq) {
            const int _size = _Mmn_RPA.getAuxDimension(); // size of gwbasis
            const int index_n = _Mmn_RPA.get_nmin();
            const int index_m = _Mmn_RPA.get_mmin();
            Eigen::MatrixXd result=Eigen::MatrixXd::Zero(_size,_size);
           
            #pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {
                const double _qp_energy_m=qp_energies(_m_level + index_m);
                
                
#if (GWBSE_DOUBLE)
                const Eigen::MatrixXd& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const Eigen::MatrixXd Mmn_RPA = _Mmn_RPA[ _m_level ].cast<double>();
#endif

                // a temporary matrix, that will get filled in empty levels loop
                Eigen::MatrixXd _temp =Eigen::MatrixXd(_Mmn_RPA.get_ntot(), _size);
                

                // loop over empty levels
                for (int _n_level = 0; _n_level < _Mmn_RPA.get_ntot(); _n_level++) {
                   
                    const double _deltaE = qp_energies(_n_level + index_n) - _qp_energy_m; // get indices and units right!!!

                    // this only works, if we have either purely real or purely imaginary frequencies

                    // purely real
                    const double _energy_factor =2.0*  (1.0 / (_deltaE - screening_freq) + 1.0 / (_deltaE + screening_freq));//hartree

                    for (int _i_gw = 0; _i_gw < _size; _i_gw++) {
                        _temp(_n_level, _i_gw) = _energy_factor * Mmn_RPA(_i_gw, _n_level);
                    } // matrix size

                } // empty levels

                // now multiply and add to epsilon
               _temp=Mmn_RPA*_temp;
                #pragma omp critical
                {
                result+=_temp;
                }
            } // occupied levels
   
            return result;
        }
        
void RPA::RPA_calculate_epsilon(const Eigen::VectorXd& qp_energies) {
  
 
      for (auto& matrix : _epsilon) {
        matrix = Eigen::MatrixXd::Identity(_Mmn_RPA.getAuxDimension(),_Mmn_RPA.getAuxDimension());
      }
      
      // loop over frequencies
      for (unsigned _i_freq = 0; _i_freq < _screening_freq.rows(); _i_freq++) {

        if (_screening_freq(_i_freq, 0) == 0.0) {
          _epsilon[ _i_freq ] += RPA_imaginary(qp_energies, _screening_freq(_i_freq, 1));
        }
        else if (_screening_freq(_i_freq, 1) == 0.0) {
          // purely real
          _epsilon[ _i_freq ] += RPA_real(qp_energies, _screening_freq(_i_freq, 0));
        }
        else {
          // mixed -> FAIL
          cerr << " mixed frequency! real part: " << _screening_freq(_i_freq, 0) << " imaginary part: " << _screening_freq(_i_freq, 1) << flush;
          exit(1);
        }

      } // loop over frequencies

      return;
    }
        
        
   
    void RPA::RPA_prepare_threecenters(const TCMatrix_gwbse& _Mmn_full){
      //TODO maybe remove and instead not make a copy but use full Mmn with restricted indices
        _Mmn_RPA.Initialize(_Mmn_full.getAuxDimension(), _rpamin, _homo, _homo + 1,
                      _rpamax);      
        unsigned start=_Mmn_RPA.get_nmin() - _Mmn_full.get_nmin();
              
            // loop over m-levels in _Mmn_RPA
            #pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {  
                // copy to _Mmn_RPA
                _Mmn_RPA[ _m_level ] =_Mmn_full[ _m_level ].block(0,start,_Mmn_full.getAuxDimension(),_Mmn_RPA.get_ntot());           
            }// loop m-levels
     return;   
    } 

    
    
    
 
}};
