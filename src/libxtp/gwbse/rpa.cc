/*
 *            Copyright 2009-2016 The VOTCA Development Team
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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/xtp/votca_xtp_config.h>

#include <votca/xtp/gwbse.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/threecenters.h>
// #include <votca/xtp/logger.h>
#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>
#include <votca/xtp/aoshell.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        
        void GWBSE::PPM_construct_parameters(  ub::matrix<double>& _overlap_cholesky_inverse ){
            
            // multiply with L-1^t from the right
            ub::matrix<double> _overlap_cholesky_inverse_transposed = ub::trans( _overlap_cholesky_inverse );
            ub::matrix<double> _temp = ub::prod( _epsilon[0] , _overlap_cholesky_inverse_transposed );
            // multiply with L-1 from the left
            _temp = ub::prod( _overlap_cholesky_inverse, _temp );
            
            // get eigenvalues and eigenvectors of this matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
            _eigenvalues.resize(_temp.size1());
            _eigenvectors.resize(_temp.size1(), _temp.size1());
            linalg_eigenvalues(_temp, _eigenvalues, _eigenvectors);
            
            // multiply eigenvectors with overlap_cholesky_inverse_transpose and store as eigenvalues of epsilon
            _ppm_phi = ub::prod( _overlap_cholesky_inverse_transposed , _eigenvectors ); 

 
            
            // store PPM weights from eigenvalues
            _ppm_weight.resize( _eigenvalues.size() );
            for ( unsigned _i = 0 ; _i <  _eigenvalues.size(); _i++   ){
                _ppm_weight(_i) = 1.0 - 1.0/_eigenvalues(_i);
            }

            // determine PPM frequencies
            _ppm_freq.resize( _eigenvalues.size() );
            // a) phi^t * epsilon(1) * phi 
            _temp = ub::prod( ub::trans( _ppm_phi ) , _epsilon[1] );
            _eigenvectors  = ub::prod( _temp ,  _ppm_phi  );
            // b) invert
            _temp = ub::zero_matrix<double>( _eigenvalues.size(),_eigenvalues.size() )  ;
            linalg_invert( _eigenvectors , _temp ); //eigenvectors is destroyed after!
            // c) PPM parameters -> diagonal elements
            for ( unsigned _i = 0 ; _i <  _eigenvalues.size(); _i++   ){
                
                double _nom  = _temp( _i, _i ) - 1.0;
                
                // only purely imaginary frequency assumed
                if ( _screening_freq(1,0) != 0.0 ) {
                    cerr << " mixed frequency! real part: " << _screening_freq( 1, 0 ) << " imaginary part: "  << _screening_freq( 1 , 1 ) << flush;
                    exit(1);
                } else {
                    
                    double _frac = -1.0 * _nom/(_nom + _ppm_weight( _i )) * _screening_freq(1,1) * _screening_freq(1,1) ;
                    _ppm_freq ( _i ) =  sqrt( std::abs(_frac )) ;

		    if ( _ppm_weight(_i) < 1.e-5 ){
		      _ppm_weight(_i) = 0.0;
		      _ppm_freq(_i)   = 1.0;

		    }

                }

            }
            
            // will be needed transposed later
            _ppm_phi = ub::trans( _ppm_phi );
            
            // epsilon can be deleted
            _epsilon[0].resize(0,0);
            _epsilon[1].resize(0,0);
                   
            
        }
        
        

        void GWBSE::RPA_calculate_epsilon(TCMatrix& _Mmn_RPA, ub::matrix<double> _screening_freq, double _shift, ub::vector<double>& _dft_energies){
            
            int _size = _Mmn_RPA[0].size1(); // size of gwbasis
            
            // loop over frequencies
            for ( unsigned _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
                
                // loop over occupied levels -> vector index of _Mmn_RPA
                // by default all variable shared, except for one defined in parallel region
                #pragma omp parallel for 
                for ( int _m_level = 0; _m_level < _Mmn_RPA.get_mtot() ; _m_level++ ){
                    //cout << " act threads: " << omp_get_thread_num( ) << " total threads " << omp_get_num_threads( ) << " max threads " << omp_get_max_threads( ) <<endl;
                    int index_m = _Mmn_RPA.get_mmin();
                    const ub::matrix<double>& Mmn_RPA =  _Mmn_RPA[ _m_level ];

                    // a temporary matrix, that will get filled in empty levels loop
                    ub::matrix<double> _temp = ub::zero_matrix<double>( _Mmn_RPA.get_ntot(), _size );
                    
                        
                    // loop over empty levels
                    for ( int _n_level = 0 ; _n_level < _Mmn_RPA.get_ntot() ; _n_level++ ){
                        int index_n = _Mmn_RPA.get_nmin();
                        
                        
                        double _deltaE = _shift + _dft_energies( _n_level + index_n ) - _dft_energies( _m_level + index_m ); // get indices and units right!!!
                        double _energy_factor;
                        // this only works, if we have either purely real or purely imaginary frequencies
                        if ( _screening_freq( _i_freq, 0) == 0.0 ) {
                            // purely imaginary
                            _energy_factor = 8.0 * _deltaE / (_deltaE*_deltaE + _screening_freq( _i_freq, 1) *  _screening_freq( _i_freq, 1 ));
                        } else if ( _screening_freq( _i_freq, 1) == 0.0  ) {
                            // purely real
                            _energy_factor = 4.0 * (1.0 / (_deltaE - _screening_freq( _i_freq, 0 ) ) +  1.0 / (_deltaE + _screening_freq( _i_freq, 0 ) ) );
                        } else {
                            // mixed -> FAIL
                            cerr << " mixed frequency! real part: " << _screening_freq( _i_freq, 0 ) << " imaginary part: "  << _screening_freq( _i_freq, 1 ) << flush;
                            exit(1);
                        }

                        
                        
                        // _temp = _energy_factor * ub::trans( Mmn_RPA );
                        for ( int _i_gw = 0 ; _i_gw < _size ; _i_gw++ ){
                            _temp( _n_level , _i_gw ) = _energy_factor * Mmn_RPA( _i_gw , _n_level );
                         } // matrix size
                        
                    } // empty levels

                   // now multiply and add to epsilon
                    ub::matrix<double> _add = ub::prod( Mmn_RPA , _temp  );
                   #pragma omp critical
                    {
                   _epsilon[ _i_freq ] += _add;// ub::prod( Mmn_RPA , _temp  );
                    }
                } // occupied levels
                
            } // loop over frequencies
            
            
        }
        
        
   
    
    void GWBSE::RPA_prepare_threecenters( TCMatrix& _Mmn_RPA, TCMatrix& _Mmn_full, AOBasis& gwbasis, AOMatrix& gwoverlap, AOMatrix& gwoverlap_inverse     ){
        

      //ub::matrix<double> _temp;
      //ub::matrix<double> _temp2;


        // loop over m-levels in _Mmn_RPA
        #pragma omp parallel for 
        for ( int _m_level = 0; _m_level < _Mmn_RPA.size() ; _m_level++ ){
        
	  //ub::matrix<double> _temp = ub::prod( gwoverlap_inverse._aomatrix , _Mmn_full[ _m_level ] );
	  // try casting for efficient prod() overloading
	  // cast _Mmn_full to double
	  ub::matrix<double> _Mmn_double = _Mmn_full[ _m_level ];
	  ub::matrix<double> _temp = ub::prod( gwoverlap_inverse._aomatrix , _Mmn_double );
	 

            // loop over n-levels in _Mmn_full 
            for ( int _n_level = 0; _n_level < _Mmn_full.get_ntot() ; _n_level++ ){

                double sc_plus  = 0.0;
                double sc_minus = 0.0;
                
                // loop over gwbasis shells
                for (vector< AOShell* >::iterator _is = gwbasis.firstShell(); _is != gwbasis.lastShell(); _is++) {
                    AOShell* _shell = gwbasis.getShell(_is);
                    double decay = (*_shell->firstGaussian())->decay;
                    //int _lmax    = _shell->getLmax();
                    int _size    = _shell->getNumFunc();
                    int _start  = _shell->getStartIndex();

                    const double pi = boost::math::constants::pi<double>();
                    double _factor = pow((2.0 *pi/decay),0.75);
                    vector<double> chi( _size, 0.0 );
                    chi[0] = _factor;

                    // some block from the fortran code that I'm not sure we need 
                    /*
                                  if ( lmax .ge. 0 ) then    
                      if(lmax .ge. 2 ) then
                       chi(10)= 6.d0*factor/sqrt(15.d0)   
                       if( lmax .ge. 4) then
                          fak = 0.25d0*factor*sqrt(beta_gwa) 
                          ! xxxx, yyyy, zzzz
                          chi(21) = fak*3.d0
                          chi(22) = chi(21)
                          chi(23) = chi(21)
                           ! xxzz, yyzz, xxyy
                           chi(30) = fak
                          chi(31) = fak
                          chi(32) = fak
                       end if
                    endif
                   end if
                     
                     */

                    // loop over all functions in shell
                    for ( int _i_gw = 0; _i_gw < _size ; _i_gw++ ){
                        double _test = _temp( _i_gw + _start, _n_level   );
                        if ( _test > 0.0  ){
                            sc_plus += chi[ _i_gw ]* _test;
                        } else if ( _test < 0.0 ){
                            sc_minus -= chi[ _i_gw ]* _test;
                        }
                    } // end loop over functions in shell

                } // end loop over all shells

                if ( _m_level <= _Mmn_RPA.get_mmax() && _n_level >= _Mmn_RPA.get_nmin()  ){
                    
                    double target = sqrt( sc_plus * sc_minus );
                    sc_plus  = target / sc_plus;
                    sc_minus = target / sc_minus;

                    // loop over gwbasis shells
                    for (vector< AOShell* >::iterator _is = gwbasis.firstShell(); _is != gwbasis.lastShell(); _is++) {
                        AOShell* _shell = gwbasis.getShell(_is);
                        double decay = (*_shell->firstGaussian())->decay;
                        int _size    = _shell->getNumFunc();
                        int _start  = _shell->getStartIndex();
                        vector<double> chi( _size, 0.0 );
                        const double pi = boost::math::constants::pi<double>();
                        double _factor = pow((2.0 *pi/decay),0.75);
                        chi[0] = _factor;
                        // loop over all functions in shell
                        for ( int _i_gw = 0; _i_gw < _size ; _i_gw++ ){
                            double _test = _temp( _i_gw + _start, _n_level   );
                            if ( _test > 0.0 && std::abs( chi[_i_gw] ) > 1.e-10 ){
                               _temp( _i_gw + _start, _n_level   ) = _temp( _i_gw + _start, _n_level   ) * sc_plus;
                            } else if ( _test < 0.0 && std::abs( chi[_i_gw] ) > 1.e-10  ){
                               _temp( _i_gw + _start, _n_level   ) = _temp( _i_gw + _start, _n_level   ) * sc_minus;
                            }
                        } // end loop over functions in shell
                    } // end loop over all shells
                    
                }                
                
            }// loop n-levels

            // multiply _temp with overlap
            ub::matrix<float> _temp2 = ub::prod( gwoverlap._aomatrix , _temp );
	    //_temp2 = ub::prod( gwoverlap._aomatrix , _temp );
            // copy to _Mmn_RPA
                      
            //ub::matrix<float> _cut = ub::project( _temp2, ub::range(0, gwbasis._AOBasisSize) , ub::range(_Mmn_RPA.get_nmin() - _Mmn_full.get_nmin()  , _Mmn_RPA.get_nmax() - _Mmn_full.get_nmin() +1 ));
            
            _Mmn_RPA[ _m_level ] = ub::project( _temp2, ub::range(0, gwbasis._AOBasisSize) , ub::range(_Mmn_RPA.get_nmin() - _Mmn_full.get_nmin()  , _Mmn_RPA.get_nmax() - _Mmn_full.get_nmin() +1 ));
            //_Mmn_RPA[ _m_level ] =  ub::zero_matrix<float>(gwbasis._AOBasisSize,_Mmn_RPA[0].size2() );
            //_Mmn_RPA[ _m_level ] = _cut;
            
        }// loop m-levels
        
    } // end RPA_prepare_threecenters


    }
    
 
};
