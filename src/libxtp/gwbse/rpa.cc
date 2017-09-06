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

        void GWBSE::PPM_construct_parameters(const ub::matrix<double>& _overlap_cholesky_inverse) {
            
            // multiply with L-1^t from the right
             
            ub::matrix<double> _temp = ub::prod(_epsilon[0], ub::trans(_overlap_cholesky_inverse));
            // multiply with L-1 from the left
            
            _temp = ub::prod(_overlap_cholesky_inverse, _temp);

            // get eigenvalues and eigenvectors of this matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
           
            linalg_eigenvalues(_temp, _eigenvalues, _eigenvectors);
            _eigenvectors=ub::trans(_eigenvectors);
            // multiply eigenvectors with overlap_cholesky_inverse and store as eigenvalues of epsilon
            _ppm_phi = ub::prod(_eigenvectors, _overlap_cholesky_inverse);
            
            // store PPM weights from eigenvalues
            _ppm_weight.resize(_eigenvalues.size());
            for (unsigned _i = 0; _i < _eigenvalues.size(); _i++) {
                _ppm_weight(_i) = 1.0 - 1.0 / _eigenvalues(_i);
            }

            // determine PPM frequencies
            _ppm_freq.resize(_eigenvalues.size());
            // a) phi^t * epsilon(1) * phi 
           
            _temp = ub::prod(_ppm_phi, _epsilon[1]);
            
            _eigenvectors = ub::prod(_temp, ub::trans(_ppm_phi));
            // b) invert
            
            linalg_invert(_eigenvectors, _temp);
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
 
            // epsilon can be deleted
            _epsilon[0].resize(0, 0);
            _epsilon[1].resize(0, 0);
            

            return;
        }
        
        
        //imaginary
        ub::matrix<double> GWBSE::RPA_imaginary(const TCMatrix& _Mmn_RPA,const double screening_freq) {
            const int _size = _Mmn_RPA.get_beta(); // size of gwbasis
            const int index_n = _Mmn_RPA.get_nmin();
            const int index_m = _Mmn_RPA.get_mmin();
            const double screenf2=screening_freq * screening_freq;
            
            std::vector<ub::matrix<double> > result_thread;
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
            const ub::vector<double>& qp_energies=   _qp_energies;
              
            for(unsigned i=0;i<nthreads;++i){
                result_thread.push_back(ub::zero_matrix<double>(_size));
            }
            
            #pragma omp parallel for 
            for(unsigned thread=0;thread<nthreads;++thread){
            for (int _m_level = thread; _m_level < _Mmn_RPA.get_mtot(); _m_level+=nthreads) {
                const double _qp_energy_m=qp_energies(_m_level + index_m);
#if (GWBSE_DOUBLE)
                const ub::matrix<double>& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const ub::matrix<double> Mmn_RPA = _Mmn_RPA[ _m_level ];
#endif
                // a temporary matrix, that will get filled in empty levels loop
                ub::matrix<double> _temp = ub::matrix<double>(_Mmn_RPA.get_ntot(), _size);

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

                // now multiply and add to epsilon
                result_thread[thread] += ub::prod(Mmn_RPA, _temp);
            } // occupied levels
            }
            ub::matrix<double> result=ub::zero_matrix<double>(_size);
               for(unsigned thread=0;thread<nthreads;++thread){
                   result+=result_thread[thread];
               }
            return result;
        }
        //real

        ub::matrix<double> GWBSE::RPA_real(const TCMatrix& _Mmn_RPA, const double screening_freq) {
            const int _size = _Mmn_RPA.get_beta(); // size of gwbasis
            const int index_n = _Mmn_RPA.get_nmin();
            const int index_m = _Mmn_RPA.get_mmin();
            const ub::vector<double>& qp_energies=   _qp_energies;
            
            std::vector<ub::matrix<double> > result_thread;
            unsigned nthreads = 1;
            #ifdef _OPENMP
               nthreads = omp_get_max_threads();
            #endif
               
              
            for(unsigned i=0;i<nthreads;++i){
                result_thread.push_back(ub::zero_matrix<double>(_size));
            }
            
            #pragma omp parallel for 
            for(unsigned thread=0;thread<nthreads;++thread){
            for (int _m_level = thread; _m_level < _Mmn_RPA.get_mtot(); _m_level+=nthreads) {
                const double _qp_energy_m=qp_energies(_m_level + index_m);
                
                
#if (GWBSE_DOUBLE)
                const ub::matrix<double>& Mmn_RPA = _Mmn_RPA[ _m_level ];
#else
                const ub::matrix<double> Mmn_RPA = _Mmn_RPA[ _m_level ];
#endif

                // a temporary matrix, that will get filled in empty levels loop
                ub::matrix<double> _temp = ub::matrix<double>(_Mmn_RPA.get_ntot(), _size);
                

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
               
               
                result_thread[thread] += ub::prod(Mmn_RPA, _temp);
                
            } // occupied levels
            }
               ub::matrix<double> result=ub::zero_matrix<double>(_size);
               for(unsigned thread=0;thread<nthreads;++thread){
                   result+=result_thread[thread];
               }
              
            return result;
        }
        

    void GWBSE::RPA_calculate_epsilon(const TCMatrix& _Mmn_RPA){

        // loop over frequencies
        for ( unsigned _i_freq = 0 ; _i_freq < _screening_freq.size1() ; _i_freq++ ){
           
             if ( _screening_freq( _i_freq, 0) == 0.0 ) {         
                 _epsilon[ _i_freq ]+=RPA_imaginary(_Mmn_RPA, _screening_freq( _i_freq, 1));
             }

             else if ( _screening_freq( _i_freq, 1) == 0.0  ) {
                  // purely real
                 _epsilon[ _i_freq ]+= RPA_real(_Mmn_RPA, _screening_freq( _i_freq, 0));
                    } 
             else {
                    // mixed -> FAIL
                    cerr << " mixed frequency! real part: " << _screening_freq( _i_freq, 0 ) << " imaginary part: "  << _screening_freq( _i_freq, 1 ) << flush;
                    exit(1);
                    }   

        } // loop over frequencies

        return;
    }
        
        
   
    void GWBSE::RPA_prepare_threecenters(TCMatrix& _Mmn_RPA,const TCMatrix& _Mmn_full){
        
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
