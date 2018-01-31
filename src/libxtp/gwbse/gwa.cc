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



#include <votca/xtp/gwbse.h>

#include <boost/format.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/constants.h>

//#include "mathimf.h"

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // MBPT MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        
        

        void GWBSE::FullQPHamiltonian(){
            
            // constructing full QP Hamiltonian, storage in vxc
            _vxc = -_vxc + _sigma_x + _sigma_c;
            // diagonal elements are given by _qp_energies
            for (unsigned _m = 0; _m < _vxc.size1(); _m++ ){
              _vxc( _m,_m ) = _qp_energies( _m + _qpmin );
            }

             // sigma matrices can be freed
            _sigma_x.resize(0);
            _sigma_c.resize(0);
            
            
            
            if ( _do_qp_diag ){
                _qp_diag_energies.resize(_vxc.size1());
                _qp_diag_coefficients.resize(_vxc.size1(), _vxc.size1());
                linalg_eigenvalues(_vxc, _qp_diag_energies, _qp_diag_coefficients);
            }
           return; 
        }
        
        
        void GWBSE::sigma_diag(const TCMatrix& _Mmn){
            
            unsigned _levelsum = _Mmn[0].size2(); // total number of bands
            unsigned _gwsize = _Mmn[0].size1(); // size of the GW basis
            const double pi = boost::math::constants::pi<double>();

            
             #pragma omp parallel for
                for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
                    const ub::matrix<real_gwbse>& Mmn = _Mmn[ _gw_level + _qpmin ];
                    double sigma_x=0;
                        for ( unsigned _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){
                            // loop over all occupied bands used in screening
                            for ( unsigned _i_occ = 0 ; _i_occ <= _homo ; _i_occ++ ){
                                sigma_x-= Mmn( _i_gw , _i_occ ) * Mmn( _i_gw , _i_occ );
                            } // occupied bands
                        } // gwbasis functions             
                    _sigma_x(_gw_level,_gw_level)=( 1.0 - _ScaHFX ) * sigma_x; 
                }
            
                if(_g_sc_max_iterations==0) {_g_sc_max_iterations=1;}
                ub::vector<double>& dftenergies=_orbitals->MOEnergies();
                // initial _qp_energies are dft energies
                ub::vector<double>_qp_old=_qp_energies;

                bool energies_converged=false;

            
	    // only diagonal elements except for in final iteration
            for (unsigned _g_iter = 0; _g_iter < _g_sc_max_iterations; _g_iter++) {
                // loop over all GW levels

                #pragma omp parallel for
                for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
                    const ub::matrix<real_gwbse>& Mmn = _Mmn[ _gw_level + _qpmin ];
                    const double qpmin = _qp_old(_gw_level + _qpmin);
                    
                    double sigma_c=0.0;
                    
                    // loop over all functions in GW basis
                    for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
                        // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
                        if (_ppm_weight(_i_gw) < 1.e-9) { continue;}
                        const double ppm_freq = _ppm_freq(_i_gw);
                        const double fac = _ppm_weight(_i_gw) * ppm_freq;
                        // loop over all bands
                        for (unsigned _i = 0; _i < _levelsum; _i++) {

                            double occ = 1.0;
                            if (_i > _homo) occ = -1.0; // sign for empty levels

                            // energy denominator
                            const double _denom = qpmin - _qp_old(_i) + occ*ppm_freq;

                            double _stab = 1.0;
                            if (std::abs(_denom) < 0.25) {
                                 
                                _stab = 0.5 * (1.0 - std::cos(4.0 * pi * std::abs(_denom)));
                               
                            }
                            
                            const double _factor =0.5* fac * _stab / _denom; //Hartree

                            // sigma_c diagonal elements
                            sigma_c += _factor * Mmn(_i_gw, _i) * Mmn(_i_gw, _i);
                           
                        }// bands

                    }// GW functions
                    _sigma_c(_gw_level, _gw_level)=sigma_c;
                    // update _qp_energies
                   
                    _qp_energies(_gw_level + _qpmin) = dftenergies(_gw_level + _qpmin) + sigma_c + _sigma_x(_gw_level, _gw_level) - _vxc(_gw_level, _gw_level);

                }// all bands
                ub::vector<double> diff= _qp_old - _qp_energies;
                energies_converged = true;
                double diff_max=0;
                unsigned state_max=0;
                for (unsigned l = 0; l < diff.size(); l++) {
                    if(std::abs(diff(l))>std::abs(diff_max)){
                            diff_max=diff(l);
                            state_max=l;
                    }
                    if (std::abs(diff(l))>_g_sc_limit) {
                        energies_converged = false;   
                    }
                }
                if(tools::globals::verbose){
                    double _DFTgap =dftenergies(_homo + 1) - dftenergies(_homo);
                    double _QPgap = _qp_energies( _homo +1 ) - _qp_energies( _homo  );
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " QP_Iteration: " << _g_iter+1 <<" shift="<<_QPgap - _DFTgap <<" E_diff max="<<diff_max<<" StateNo:"<<state_max << flush;
                }
                double alpha=0.0;
                _qp_energies=(1-alpha)*_qp_energies+alpha*_qp_old;
                
                if (energies_converged) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converged after " << _g_iter+1 << " G iterations." << flush;
                    break;
                }else if(_g_iter==_g_sc_max_iterations-1){
                        CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " G-self-consistency cycle not converged after " << _g_sc_max_iterations << " iterations." << flush;  
                        break;
                    
                } else {
                    _qp_old = _qp_energies;
                }

            } // iterations

            return;
        }

      
       
         void GWBSE::sigma_offdiag(const TCMatrix& _Mmn) {
            unsigned _levelsum = _Mmn[0].size2(); // total number of bands
            unsigned _gwsize = _Mmn[0].size1(); // size of the GW basis
            const double pi = boost::math::constants::pi<double>();
           
            #pragma omp parallel for
            for (unsigned _gw_level1 = 0; _gw_level1 < _qptotal; _gw_level1++) {
                const ub::matrix<real_gwbse>& Mmn1 =  _Mmn[ _gw_level1 + _qpmin ];
                for (unsigned _gw_level2 = 0; _gw_level2 < _gw_level1; _gw_level2++) {
                    const ub::matrix<real_gwbse>& Mmn2 =  _Mmn[ _gw_level2 + _qpmin ];
                    double sigma_x=0;
                    for ( unsigned _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){
                        // loop over all occupied bands used in screening
                        for ( unsigned _i_occ = 0 ; _i_occ <= _homo ; _i_occ++ ){
                            sigma_x -= Mmn1( _i_gw , _i_occ ) * Mmn2( _i_gw , _i_occ );
                        } // occupied bands
                    } // gwbasis functions
                    _sigma_x(_gw_level1, _gw_level2)=( 1.0 - _ScaHFX ) * sigma_x;
                }
            }
            
            #pragma omp parallel for
            for (unsigned _gw_level1 = 0; _gw_level1 < _qptotal; _gw_level1++) {
                const double qpmin1 = _qp_energies(_gw_level1 + _qpmin);
                const ub::matrix<real_gwbse>& Mmn1 = _Mmn[ _gw_level1 + _qpmin ];
                for (unsigned _gw_level2 = 0; _gw_level2 < _gw_level1; _gw_level2++) {
                    const double qpmin2 = _qp_energies(_gw_level1 + _qpmin);
                    const ub::matrix<real_gwbse>& Mmn2 = _Mmn[ _gw_level2 + _qpmin ];
                    double sigma_c = 0;
                    for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
                        // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
                        if (_ppm_weight(_i_gw) < 1.e-9) {
                            continue;
                        }
                        const double ppm_freq = _ppm_freq(_i_gw);
                        const double fac = _ppm_weight(_i_gw) * ppm_freq;
                        // loop over all screening levels
                        for (unsigned _i = 0; _i < _levelsum; _i++) {

                            double occ = 1.0;
                            if (_i > _homo) occ = -1.0; // sign for empty levels

                            // energy denominator
                            const double _denom1 = qpmin1 - _qp_energies(_i) + occ * ppm_freq;
                            const double _denom2 = qpmin2 - _qp_energies(_i) + occ * ppm_freq;

                            double _stab1 = 1.0;
                            if (std::abs(_denom1) < 0.25) {
                                _stab1 = 0.5 * (1.0 - std::cos(4.0 * pi * std::abs(_denom1)));
                            }
                            const double factor1 = 0.5 * fac * _stab1 / _denom1; //Hartree
                            double _stab2 = 1.0;
                            if (std::abs(_denom2) < 0.25) {
                                _stab2 = 0.5 * (1.0 - std::cos(4.0 * pi * std::abs(_denom2)));
                            }
                            const double factor2 = 0.5 * fac * _stab2 / _denom2; //Hartree
                            sigma_c +=Mmn1(_i_gw, _i) * Mmn2(_i_gw, _i)*0.5*(factor1+factor2);
                        }// screening levels 
                    }// GW functions 

                    _sigma_c(_gw_level1, _gw_level2) = sigma_c;


                }// GW row             
            }//GW col
         
        return;
        } 


        void GWBSE::sigma_prepare_threecenters(TCMatrix& _Mmn){
            #if (GWBSE_DOUBLE)
                const ub::matrix<double>& ppm_phi=_ppm_phi;
            #else
                const ub::matrix<float> ppm_phi=_ppm_phi;        
            #endif
            
            
            #pragma omp parallel for
            for ( int _m_level = 0 ; _m_level < _Mmn.get_mtot(); _m_level++ ){
                // get Mmn for this _m_level
                // and multiply with _ppm_phi = eigenvectors of epsilon
              _Mmn[ _m_level ] = ub::prod(  ppm_phi , _Mmn[_m_level] );
            }
            return;
        }        
        


    }
    
 
};
