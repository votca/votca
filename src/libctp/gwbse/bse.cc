/*
 *            Copyright 2009-2012 The VOTCA Development Team
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

#define NDEBUG

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/gwbse.h>

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>
// #include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // MBPT MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //


        
        void GWBSE::BSE_qp_setup(){
            _eh_qp = ub::zero_matrix<float>( _bse_size , _bse_size );
            
            
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    _eh_qp( _index_vc , _index_vc ) += _vxc(_c1 + _bse_vtotal ,_c1 + _bse_vtotal ) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            _eh_qp( _index_vc , _index_vc2 ) += _vxc(_c1+ _bse_vtotal ,_c2 + _bse_vtotal );
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _eh_qp( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
            
            
        }
        
    


        void GWBSE::BSE_solve_triplets(){
            
            ub::matrix<float> _bse =  -_eh_d;
            
            // add full QP Hamiltonian contributions to free transitions
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    //_bse( _index_vc , _index_vc ) += _vxc(_c1 + _bse_cmin ,_c1 + _bse_cmin ) - _vxc(_v1,_v1);
                    _bse( _index_vc , _index_vc ) += _vxc(_c1 + _bse_vtotal ,_c1 + _bse_vtotal ) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            //_bse( _index_vc , _index_vc2 ) += _vxc(_c1+ _bse_cmin ,_c2 + _bse_cmin );
                            _bse( _index_vc , _index_vc2 ) += _vxc(_c1+ _bse_vtotal ,_c2 + _bse_vtotal );
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _bse( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
            
            
            
            linalg_eigenvalues( _bse, _bse_triplet_energies, _bse_triplet_coefficients, _bse_nmax);
        }
        

   
        
      void GWBSE::BSE_solve_singlets(){
            
            ub::matrix<float> _bse = -_eh_d + 2.0 * _eh_x;


            // add full QP Hamiltonian contributions to free transitions
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    _bse( _index_vc , _index_vc ) += _vxc(_c1 + _bse_vtotal ,_c1 + _bse_vtotal) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            _bse( _index_vc , _index_vc2 ) += _vxc(_c1 + _bse_vtotal ,_c2 + _bse_vtotal);
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _bse( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
            
            // _bse_singlet_energies.resize(_bse_singlet_coefficients.size1());
            linalg_eigenvalues(_bse, _bse_singlet_energies, _bse_singlet_coefficients, _bse_nmax);
            
       
            
        } 
        
        
        void GWBSE::BSE_d_setup ( TCMatrix& _Mmn){
            // gwbasis size
            size_t _gwsize = _Mmn[_homo].size1();

            // messy procedure, first get two matrices for occ and empty subbparts
            // store occs directly transposed
            ub::matrix<float> _storage_v = ub::zero_matrix<float>(  _bse_vtotal * _bse_vtotal , _gwsize );
            #pragma omp parallel for
            for ( size_t _v1 = 0; _v1 < _bse_vtotal; _v1++){
                const ub::matrix<float>& Mmn = _Mmn[_v1 + _bse_vmin ];
                for ( size_t _v2 = 0; _v2 < _bse_vtotal; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_v( _index_vv , _i_gw ) = Mmn( _i_gw , _v2 + _bse_vmin );
                    }
                }
            }
            
            
            ub::matrix<float> _storage_c = ub::zero_matrix<float>( _gwsize, _bse_ctotal * _bse_ctotal );
            #pragma omp parallel for
            for ( size_t _c1 = 0; _c1 < _bse_ctotal; _c1++){
                const ub::matrix<float>& Mmn = _Mmn[_c1 + _bse_cmin];
                for ( size_t _c2 = 0; _c2 < _bse_ctotal; _c2++){
                    size_t _index_cc = _bse_ctotal * _c1 + _c2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_c( _i_gw , _index_cc ) = Mmn( _i_gw , _c2 + _bse_cmin );
                    }
                }
            }
            
            if ( ! _do_bse_singlets )  _Mmn.Cleanup();
            
            // store elements in a vtotal^2 x ctotal^2 matrix
            // cout << "BSE_d_setup 1 [" << _storage_v.size1() << "x" << _storage_v.size2() << "]\n" << std::flush;
            ub::matrix<float> _storage_prod = ub::prod( _storage_v , _storage_c );
            
            
            
            // now patch up _storage for screened interaction
            #pragma omp parallel for
            for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){  
                double _ppm_factor = sqrt( _ppm_weight( _i_gw ));
                for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_vtotal) ; _v++){
                    _storage_v( _v , _i_gw ) = _ppm_factor * _storage_v(_v , _i_gw );
                }
                for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_ctotal) ; _c++){
                    _storage_c( _i_gw , _c ) = _ppm_factor * _storage_c( _i_gw , _c  );
                }
            }
            
            // multiply and subtract from _storage_prod
            // cout << "BSE_d_setup 2 [" << _storage_c.size1() << "x" << _storage_c.size2() << "]\n" << std::flush;
            _storage_prod -= ub::prod( _storage_v , _storage_c );
            
            // free storage_v and storage_c
            _storage_c.resize(0,0);
            _storage_v.resize(0,0);
            
            
            // finally resort into _eh_d and multiply by to for Rydbergs
            // can be limited to upper diagonal !
            _eh_d = ub::zero_matrix<float>( _bse_size , _bse_size );
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_vc1 = _bse_ctotal * _v1 + _c1 ;
                              
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_vc2 = _bse_ctotal * _v2 + _c2 ;
                            size_t _index_cc  = _bse_ctotal * _c1 + _c2;

                            _eh_d( _index_vc1 , _index_vc2 ) = 2.0 * _storage_prod( _index_vv , _index_cc ); 

                            
                        }
                    }
                }
            }
            
            
        }
        
        
        
        void GWBSE::BSE_x_setup( TCMatrix& _Mmn){
            
            /* unlike the fortran code, we store eh interaction directly in
             * a suitable matrix form instead of a four-index array
             */
            
            
            // gwbasis size
            size_t _gwsize = _Mmn[_homo].size1();
            
            // get a different storage for 3-center integrals we need
            ub::matrix<float> _storage = ub::zero_matrix<float>( _gwsize , _bse_size);


            
            // occupied levels
            #pragma omp parallel for
            for ( size_t _v = 0; _v < _bse_vtotal ; _v++ ){
                // cout << " act threads: " << omp_get_thread_num( ) << " total threads " << omp_get_num_threads( ) << " max threads " << omp_get_max_threads( ) <<endl;
                ub::matrix<float>& Mmn = _Mmn[_v + _bse_vmin];
                // empty levels
                for (size_t _c =0 ; _c < _bse_ctotal ; _c++ ){
                    size_t _index_vc = _bse_ctotal * _v + _c ;
                    for (size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){
                        _storage( _i_gw, _index_vc ) = Mmn( _i_gw, _c + _bse_cmin);
                    }
                }
            }

            
            _Mmn.Cleanup();
            

            
            // with this storage, _eh_x is obtained by matrix multiplication
            
	    _eh_x = ub::prod( ub::trans( _storage ), _storage ); 
            _eh_x = 2.0 * _eh_x; // Rydberg
  
            
            
            
            
            
            
            
            
            
        }
        
        
        

    }
    
 
};
