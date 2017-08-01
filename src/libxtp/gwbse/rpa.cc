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
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/threecenters.h>
// #include <votca/xtp/logger.h>
#include <votca/xtp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/xtp/aoshell.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

        void GWBSE::PPM_construct_parameters(const ub::matrix<double>& _overlap_cholesky_inverse,const ub::matrix<double>& transcholeskyinverse) {
            
            // multiply with L-1^t from the right
             
            ub::matrix<double> _temp = ub::prod(_epsilon[0], ub::trans(_overlap_cholesky_inverse));
            // multiply with L-1 from the left
            
            _temp = ub::prod(_overlap_cholesky_inverse, _temp);

            // get eigenvalues and eigenvectors of this matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
           
            linalg_eigenvalues(_temp, _eigenvalues, _eigenvectors);
            
            // multiply eigenvectors with overlap_cholesky_inverse and store as eigenvalues of epsilon
            _ppm_phi = ub::prod(transcholeskyinverse, _eigenvectors);
            
            // store PPM weights from eigenvalues
            _ppm_weight.resize(_eigenvalues.size());
            for (unsigned _i = 0; _i < _eigenvalues.size(); _i++) {
                _ppm_weight(_i) = 1.0 - 1.0 / _eigenvalues(_i);
            }

            // determine PPM frequencies
            _ppm_freq.resize(_eigenvalues.size());
            // a) phi^t * epsilon(1) * phi 
            _temp=ub::trans(_ppm_phi);
            _temp = ub::prod(_temp, _epsilon[1]);
            
            _eigenvectors = ub::prod(_temp, _ppm_phi);
            // b) invert
            
            _temp = ub::zero_matrix<double>(_eigenvalues.size(), _eigenvalues.size());
            
            linalg_invert(_eigenvectors, _temp); //eigenvectors is destroyed after!
            // c) PPM parameters -> diagonal elements
            
            #pragma omp parallel for 
            for (unsigned _i = 0; _i < _eigenvalues.size(); _i++) {

                if (_screening_freq(1, 0) == 0.0) {
                    if (_ppm_weight(_i) < 1.e-5) {
                        _ppm_weight(_i) = 0.0;
                        _ppm_freq(_i) = 0.5;//Hartree
                        continue;
                    } else {
                        double _nom = _temp(_i, _i) - 1.0;
                        double _frac = -1.0 * _nom / (_nom + _ppm_weight(_i)) * _screening_freq(1, 1) * _screening_freq(1, 1);
                        _ppm_freq(_i) = sqrt(std::abs(_frac));
                    }

                } else {
                    // only purely imaginary frequency assumed
                    cerr << " mixed frequency! real part: " << _screening_freq(1, 0) << " imaginary part: " << _screening_freq(1, 1) << flush;
                    exit(1);
                }

            }

            
            // will be needed transposed later
            _ppm_phi = ub::trans(_ppm_phi);
            
            // epsilon can be deleted
            _epsilon[0].resize(0, 0);
            _epsilon[1].resize(0, 0);
            

            return;
        }
        
        
        //imaginary
        void GWBSE::RPA_imaginary(ub::matrix<double>& result, const TCMatrix& _Mmn_RPA,const double screening_freq) {
            int _size = _Mmn_RPA[0].size1(); // size of gwbasis
            
            #pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {
                int index_m = _Mmn_RPA.get_mmin();
#if (GWBSE_DOUBLE)
                const ub::matrix<double>& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const ub::matrix<double> Mmn_RPA = _Mmn_RPA[ _m_level ];
#endif
                // a temporary matrix, that will get filled in empty levels loop
                ub::matrix<double> _temp = ub::zero_matrix<double>(_Mmn_RPA.get_ntot(), _size);

                // loop over empty levels
                for (int _n_level = 0; _n_level < _Mmn_RPA.get_ntot(); _n_level++) {
                    int index_n = _Mmn_RPA.get_nmin();

                    double _deltaE = _qp_energies(_n_level + index_n) - _qp_energies(_m_level + index_m); // get indices and units right!!!

                    // this only works, if we have either purely real or purely imaginary frequencies

                    // purely imaginary
                    double _energy_factor = 4.0 * _deltaE / (_deltaE * _deltaE + screening_freq * screening_freq);//hartree
                    for (int _i_gw = 0; _i_gw < _size; _i_gw++) {
                        _temp(_n_level, _i_gw) = _energy_factor * Mmn_RPA(_i_gw, _n_level);
                    } // matrix size

                } // empty levels

                // now multiply and add to epsilon
                ub::matrix<double> _add = ub::prod(Mmn_RPA, _temp);
                #pragma omp critical
                {
                    result += _add; 
                }
            } // occupied levels
            return;
        }
        //real

        void GWBSE::RPA_real(ub::matrix<double>& result, const TCMatrix& _Mmn_RPA, const double screening_freq) {
            int _size = _Mmn_RPA[0].size1(); // size of gwbasis
           
            #pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.get_mtot(); _m_level++) {
                
                int index_m = _Mmn_RPA.get_mmin();
                
#if (GWBSE_DOUBLE)
                const ub::matrix<double>& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const ub::matrix<double> Mmn_RPA = _Mmn_RPA[ _m_level ];
#endif

                // a temporary matrix, that will get filled in empty levels loop
                ub::matrix<double> _temp = ub::zero_matrix<double>(_Mmn_RPA.get_ntot(), _size);


                // loop over empty levels
                for (int _n_level = 0; _n_level < _Mmn_RPA.get_ntot(); _n_level++) {
                    int index_n = _Mmn_RPA.get_nmin();


                    double _deltaE = _qp_energies(_n_level + index_n) - _qp_energies(_m_level + index_m); // get indices and units right!!!

                    // this only works, if we have either purely real or purely imaginary frequencies

                    // purely real
                    double _energy_factor =2.0*  (1.0 / (_deltaE - screening_freq) + 1.0 / (_deltaE + screening_freq));//hartree

                    for (int _i_gw = 0; _i_gw < _size; _i_gw++) {
                        _temp(_n_level, _i_gw) = _energy_factor * Mmn_RPA(_i_gw, _n_level);
                    } // matrix size

                } // empty levels

                // now multiply and add to epsilon
                ub::matrix<double> _add = ub::prod(Mmn_RPA, _temp);
                #pragma omp critical
                {
                    result += _add; 
                }
            } // occupied levels
            return;
        }
        

    void GWBSE::RPA_calculate_epsilon(const TCMatrix& _Mmn_RPA){

        // loop over frequencies
        for ( unsigned _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
           
             if ( _screening_freq( _i_freq, 0) == 0.0 ) {         
                 RPA_imaginary(_epsilon[ _i_freq ],_Mmn_RPA, _screening_freq( _i_freq, 1));
             }

             else if ( _screening_freq( _i_freq, 1) == 0.0  ) {
                  // purely real
                  RPA_real(_epsilon[ _i_freq ],_Mmn_RPA, _screening_freq( _i_freq, 0));
                    } 
             else {
                    // mixed -> FAIL
                    cerr << " mixed frequency! real part: " << _screening_freq( _i_freq, 0 ) << " imaginary part: "  << _screening_freq( _i_freq, 1 ) << flush;
                    exit(1);
                    }   

        } // loop over frequencies

        return;
    }
        
        
   
    // this is apparently the fourier transform of the coulomb matrix, I am not sure according to old gwbse code it restores  orthogonality
    void GWBSE::RPA_prepare_threecenters(TCMatrix& _Mmn_RPA,const TCMatrix& _Mmn_full){
        
        const double pi = boost::math::constants::pi<double>();
        ub::range full=ub::range(0, _Mmn_full.get_beta());
        ub::range RPA_cut=ub::range(_Mmn_RPA.get_nmin() - _Mmn_full.get_nmin(), _Mmn_RPA.get_nmax() - _Mmn_full.get_nmin() + 1);
            // loop over m-levels in _Mmn_RPA
            #pragma omp parallel for 
            for (int _m_level = 0; _m_level < _Mmn_RPA.size(); _m_level++) {
          
                // copy to _Mmn_RPA
                _Mmn_RPA[ _m_level ] = ub::project(_Mmn_full[ _m_level ], full, RPA_cut);
              

            }// loop m-levels
     return;   
    } 

    
    
    
 
}};
