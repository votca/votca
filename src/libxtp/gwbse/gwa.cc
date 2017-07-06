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

        void GWBSE::sigma_c_setup(const TCMatrix& _Mmn) {

            // iterative refinement of qp energies
            int _max_iter = 15;
            unsigned _levelsum = _Mmn[0].size2(); // total number of bands
            unsigned _gwsize = _Mmn[0].size1(); // size of the GW basis
            const double pi = boost::math::constants::pi<double>();

            ub::vector<double>& dftenergies=_orbitals->MOEnergies();
            // initial _qp_energies are dft energies
            ub::vector<double>_qp_old_rpa=_qp_energies;
            ub::vector<double>_qp_old=_qp_energies;
            
            bool energies_converged=false;
            _sigma_c.resize(_qptotal);
	    // only diagonal elements except for in final iteration
            for (int _i_iter = 0; _i_iter < _max_iter - 1; _i_iter++) {
                // loop over all GW levels

                #pragma omp parallel for
                for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
                    
                    const double qpmin = _qp_old(_gw_level + _qpmin);
                    const ub::matrix<real_gwbse>& Mmn = _Mmn[ _gw_level + _qpmin ];
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
                                _stab = 0.25 * (1.0 - std::cos(4.0 * pi * std::abs(_denom)));
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
                
                _qp_old = _qp_old - _qp_energies;
                energies_converged = true;
                for (unsigned l = 0; l < _qp_old.size(); l++) {
                    if (std::abs(_qp_old(l)) > _qp_limit) {
                        energies_converged = false;
                        break;
                    }
                }
                if (energies_converged) {
                    CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Converged after " << _i_iter+1 << " qp_energy iterations." << flush;
                    break;
                } else {
                    _qp_old = _qp_energies;
                }

            } // iterations

            _qp_converged = true;
            _qp_old_rpa= _qp_old_rpa - _qp_energies;
            for (unsigned l = 0; l < _qp_old_rpa.size(); l++) {
                    if (std::abs(_qp_old_rpa(l)) > _shift_limit) {
                        _qp_converged = false;
                        break;
                    }
                } 
            double _DFTgap =dftenergies(_homo + 1) - dftenergies(_homo);
            double _QPgap = _qp_energies( _homo +1 ) - _qp_energies( _homo  );
            _shift = _QPgap - _DFTgap;
            
            
            if(_iterate_qp){
            // qp energies outside the update range are simply shifted. 
            for(unsigned i=_qpmax+1;i<dftenergies.size();++i){
                _qp_energies(i)=dftenergies(i)+_shift;
            }
            }else{
                _qp_converged = true;
            }


            // only if _shift is converged
            if (_qp_converged) {
                // in final step, also calc offdiagonal elements
                // initialize sigma_c to zero at the beginning
                
                
            
                
      
            //this is not the fastest algorithm but faster ones throw igwbse off, so this is good enough.    
            #pragma omp parallel for
            for (unsigned _gw_level = 0; _gw_level < _qptotal; _gw_level++) {
                const double qpmin=_qp_energies(_gw_level + _qpmin);

                const ub::matrix<real_gwbse>& Mmn = _Mmn[ _gw_level + _qpmin ];
                for (unsigned _m = 0; _m < _gw_level; _m++) {
                     double sigma_c = 0;
                    const ub::matrix<real_gwbse>& Mmn2 = _Mmn[_m + _qpmin];

                    // loop over all functions in GW basis
                    for (unsigned _i_gw = 0; _i_gw < _gwsize; _i_gw++) {
                        // the ppm_weights smaller 1.e-5 are set to zero in rpa.cc PPM_construct_parameters
                        if (_ppm_weight(_i_gw) < 1.e-9) { continue;}
                        const double ppm_freq = _ppm_freq(_i_gw);
                        const double fac = _ppm_weight(_i_gw) * ppm_freq;
                        // loop over all screening levels
                        for (unsigned _i = 0; _i < _levelsum; _i++) {

                            double occ = 1.0;
                            if (_i > _homo) occ = -1.0; // sign for empty levels

                            // energy denominator
                            const double _denom = qpmin - _qp_energies(_i) + occ * ppm_freq;

                            double _stab = 1.0;
                            if (std::abs(_denom) < 0.25) {
                                _stab = 0.25 * (1.0 - std::cos(4.0 * pi * std::abs(_denom)));
                            }

                            const double _factor = 0.5*fac * Mmn(_i_gw, _i) * _stab / _denom; //Hartree

                            sigma_c+= _factor * Mmn2(_i_gw, _i);


                        }// screening levels 
                    }// GW functions 
                    _sigma_c(_gw_level, _m)=sigma_c;
                }// GW row 
            } // GW col 
        } 
            
        return;
        } // sigma_c_setup


        void GWBSE::sigma_x_setup(const TCMatrix& _Mmn){
        
            // initialize sigma_x
            _sigma_x.resize(_qptotal);
            int _size  = _Mmn[0].size1();

            // band 1 loop over all GW levels
            #pragma omp parallel for
            for ( unsigned _m1 = 0 ; _m1 < _qptotal ; _m1++ ){
                
                const ub::matrix<real_gwbse>& M1mn =  _Mmn[ _m1 + _qpmin ];
                
                // band 2 loop over all GW levels
                //for ( int _m2 = _qpmin ; _m2 <= _qpmax ; _m2++ ){
                for ( unsigned _m2 = 0 ; _m2 <= _m1 ; _m2++ ){
                    _sigma_x( _m1, _m2 )=0;
                    const ub::matrix<real_gwbse>& M2mn =  _Mmn[ _m2 + _qpmin ];
                    
                    // loop over all basis functions
                    for ( int _i_gw = 0 ; _i_gw < _size ; _i_gw++ ){
                        // loop over all occupied bands used in screening
                        for ( unsigned _i_occ = 0 ; _i_occ <= _homo ; _i_occ++ ){
                            _sigma_x( _m1, _m2 ) -= M1mn( _i_gw , _i_occ ) * M2mn( _i_gw , _i_occ );
                        } // occupied bands
                    } // gwbasis functions
                } // level 2
            } // level 1

	    // factor for hybrid DFT
	    _sigma_x = ( 1.0 - _ScaHFX ) * _sigma_x;
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
