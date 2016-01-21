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

#define NDEBUG

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

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace xtp {
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
            
            
            
        /*      cout << endl;
        for ( int i = 0; i<_eh_d.size1(); i++ ){
            for ( int j = 0; j<_eh_d.size2(); j++ ){
                
                
                cout << "eh D [" << i << " : " << j << "]: " << _bse(i,j) << endl; 
                
            }
            
        } */
            
            
            linalg_eigenvalues( _bse, _bse_triplet_energies, _bse_triplet_coefficients, _bse_nmax);
        }
        
        
        
      void GWBSE::BSE_solve_singlets_BTDA(){
        
          
        // For details of the method, see EPL,78(2007)12001,
        // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.
        
          // setup resonant (A) and RARC blocks (B)

            ub::matrix<float> _A = -_eh_d + 2.0 * _eh_x;


            // add full QP Hamiltonian contributions to free transitions
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;

                    // diagonal
                    _A( _index_vc , _index_vc ) += _vxc(_c1 + _bse_vtotal ,_c1 + _bse_vtotal) - _vxc(_v1,_v1);

                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            _A( _index_vc , _index_vc2 ) += _vxc(_c1 + _bse_vtotal ,_c2 + _bse_vtotal);
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            _A( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    }
                    
                    
                }
            }
          
            ub::matrix<float> _B = -_eh_d2 + 2.0 * _eh_x;

            ub::matrix<float> _ApB = _A + _B;
            ub::matrix<float> _AmB = _A - _B;
   
          // check for positive definiteness of A-B
            ub::vector<float> _eigenvalues;
            ub::matrix<float> _eigenvectors;
            int dim = _AmB.size1();
            ub::matrix<float> _test = _AmB;
            linalg_eigenvalues(_test, _eigenvalues, _eigenvectors, dim);
            
            for ( unsigned _i = 0; _i < _eigenvalues.size(); _i++ ) {
                
                if ( _eigenvalues(_i) < 0.0 ) {
                    LOG(logDEBUG, *_pLog) << TimeStamp() << " KAA-KAB has negative eigenvalues! : " << _i << " : " << _eigenvalues[_i] << flush;
                }

            }
          // if, positive definite, calculate Cholesky decomposition of A-B = LL^T (A-B) is not needed any longer and can be overwritten
            
            
            
            linalg_cholesky_decompose( _AmB );
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Cholesky decomposition of KAA-KAB" << flush;

            // remove L^T from Cholesky
            for (unsigned i =0; i < _AmB.size1(); i++ ){
                for (unsigned j = i+1; j < _AmB.size1(); j++ ){
                    _AmB(i,j) = 0.0;
                }
            }
            
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Removed L^T" << flush;
          // determine H = L^T(A-B)L
          ub::matrix<float> _temp = ub::prod( _ApB , _AmB );
          ub::matrix<float> _H = ub::prod( ub::trans(_AmB), _temp );
          _temp.resize(0,0);
          LOG(logDEBUG, *_pLog) << TimeStamp() << " Calculated H = L^T(A-B)L " << flush;
          
          // solve eigenvalue problem: HR_l = eps_l^2 R_l
          linalg_eigenvalues(_H, _eigenvalues, _eigenvectors, dim);
          LOG(logDEBUG, *_pLog) << TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;

          // reconstruct real eigenvalues eps_l = sqrt(eps_l^2)
          for ( unsigned _i = 0; _i < _eigenvalues.size(); _i++) {
            LOG(logDEBUG, *_pLog) << TimeStamp() << " Eigenvalue: " << _eigenvalues(_i) << " : " << sqrt(_eigenvalues(_i))*13.6058 <<  flush;  
          }
          
          
          // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 + 1/sqrt(eps_l)L ] R_l
          //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 - 1/sqrt(eps_l)L ] R_l

          // determine inverse of L^T
          ub::matrix<float> _cholesky_transposed_invert;
          ub::matrix<float> _cholesky_transposed =ub::trans(_AmB);
          linalg_invert( _cholesky_transposed , _cholesky_transposed_invert );
          
          ub::matrix<float> _X_evec(dim,dim);
          ub::matrix<float> _Y_evec(dim,dim);
          
          for ( int _i = 0; _i < dim; _i++) {
              
              float sqrt_eval = sqrt(_eigenvalues(_i));
              // get l-th reduced EV
              ub::matrix<float> _reduced_evec = ub::project(_eigenvectors,  ub::range(0, dim), ub::range(_i, _i + 1)); // potentially col<->row

              ub::matrix<float> _transform = 0.5*( sqrt_eval * _cholesky_transposed_invert + 1.0/sqrt_eval * _AmB  );
              ub::project( _X_evec, ub::range (0, dim ), ub::range ( _i, _i+1 ) ) = ub::prod(_transform,_reduced_evec);
              _transform = 0.5*( sqrt_eval * _cholesky_transposed_invert - 1.0/sqrt_eval * _AmB  );
              ub::project( _Y_evec, ub::range (0, dim ), ub::range ( _i, _i+1 ) ) = ub::prod(_transform,_reduced_evec);
              
             /* // check normalization:
              ub::matrix<float> X = ub::project( _X_evec, ub::range (0, dim ), ub::range ( _i, _i+1 ) );
              ub::matrix<float> Y = ub::project( _Y_evec, ub::range (0, dim ), ub::range ( _i, _i+1 ) );
              ub::matrix<float> norm = ub::prod(ub::trans(X),X) - ub::prod(ub::trans(Y),Y);
              
              LOG(logDEBUG, *_pLog) << TimeStamp() << " Norm of Eigenvector: " << norm(0,0) << "[" << norm.size1() << ":" << norm.size2() <<  flush;  */
              
              
          }
          
          
          
          
          // check orthogonality
          for ( int _i = 0; _i < dim; _i++) {
              
              for ( int _j = 0; _j < dim; _j++) {
                   // check normalization:
              ub::matrix<float> X1 = ub::project( _X_evec, ub::range (0, dim ), ub::range ( _i, _i+1 ) );
              ub::matrix<float> X2 = ub::project( _X_evec, ub::range (0, dim ), ub::range ( _j, _j+1 ) );
              ub::matrix<float> Y1 = ub::project( _Y_evec, ub::range (0, dim ), ub::range ( _i, _i+1 ) );
              ub::matrix<float> Y2 = ub::project( _Y_evec, ub::range (0, dim ), ub::range ( _j, _j+1 ) );
              ub::matrix<float> norm = ub::prod(ub::trans(X1),X2) - ub::prod(ub::trans(Y1),Y2);
              
              LOG(logDEBUG, *_pLog) << TimeStamp() << " Orthogonality: (" << _i << " : " << _j << ") = " <<  norm(0,0) <<  flush;  
                  
              }
          }
         
          
          exit(0);
          
          

          
          
          
          
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
            
            // print check
            
      /*                 for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_vc1 = _bse_ctotal * _v1 + _c1 ;
                              
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_vc2 = _bse_ctotal * _v2 + _c2 ;
                            size_t _index_cc  = _bse_ctotal * _c1 + _c2;

                            
                            
                            cout << " eh_d1 : " << _index_vc1 << " : " << _index_vc2 << " = " <<  _eh_d( _index_vc1 , _index_vc2 ) << endl;
     
                            
                        }
                    }
                }
            } */
            
            
            
        }
        
        
         void GWBSE::BSE_d2_setup ( TCMatrix& _Mmn){
            // gwbasis size
            size_t _gwsize = _Mmn[_homo].size1();

            // messy procedure, first get two matrices for occ and empty subbparts
            // store occs directly transposed
            ub::matrix<float> _storage_cv = ub::zero_matrix<float>(  _bse_vtotal * _bse_ctotal , _gwsize );
            #pragma omp parallel for
            for ( size_t _c1 = 0; _c1 < _bse_ctotal; _c1++){
                const ub::matrix<float>& Mmn = _Mmn[_c1 + _bse_cmin ];
                for ( size_t _v2 = 0; _v2 < _bse_vtotal; _v2++){
                    size_t _index_cv = _bse_vtotal * _c1 + _v2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_cv( _index_cv , _i_gw ) = Mmn( _i_gw , _v2 + _bse_vmin );
                    }
                }
            }
            
            cout << " ======= check 1 ======= " << endl;
                    
            
            
            
            ub::matrix<float> _storage_vc = ub::zero_matrix<float>( _gwsize, _bse_vtotal * _bse_ctotal );
            #pragma omp parallel for
            for ( size_t _v1 = 0; _v1 < _bse_vtotal; _v1++){
                const ub::matrix<float>& Mmn = _Mmn[_v1 + _bse_vmin];
                for ( size_t _c2 = 0; _c2 < _bse_ctotal; _c2++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_vc( _i_gw , _index_vc ) = Mmn( _i_gw , _c2 + _bse_cmin );
                    }
                }
            }
            
            
            cout << " ======= check 2 ========= " << endl;
            
            
            if ( ! _do_bse_singlets )  _Mmn.Cleanup();
            
            // store elements in a vtotal^2 x ctotal^2 matrix
            // cout << "BSE_d_setup 1 [" << _storage_v.size1() << "x" << _storage_v.size2() << "]\n" << std::flush;
            ub::matrix<float> _storage_prod = ub::prod( _storage_cv , _storage_vc );
            
            cout << " ======= check 3 ========== PASSED" << endl;
            
            // now patch up _storage for screened interaction
            #pragma omp parallel for
            for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){  
                double _ppm_factor = sqrt( _ppm_weight( _i_gw ));
                for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_vtotal) ; _v++){
                    _storage_vc(  _i_gw , _v ) = _ppm_factor * _storage_vc( _i_gw , _v );
                }
                for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_ctotal) ; _c++){
                    _storage_cv( _c, _i_gw  ) = _ppm_factor * _storage_cv( _c , _i_gw  );
                }
            }
            
            // multiply and subtract from _storage_prod
            // cout << "BSE_d_setup 2 [" << _storage_c.size1() << "x" << _storage_c.size2() << "]\n" << std::flush;
            _storage_prod -= ub::prod( _storage_cv , _storage_vc );
            
            // free storage_v and storage_c
            _storage_cv.resize(0,0);
            _storage_vc.resize(0,0);
            
            
            cout << " ==== check 4 ====== " << endl;
            
            
            // finally resort into _eh_d and multiply by to for Rydbergs
            // can be limited to upper diagonal !
            _eh_d2 = ub::zero_matrix<float>( _bse_size , _bse_size );
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                    
                    
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_v1c1 = _bse_ctotal * _v1 + _c1 ;

                        //// OK?

                        size_t _index_c1v2 =_bse_vtotal * _c1 + _v2;
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_v2c2 = _bse_ctotal * _v2 + _c2 ;
                            size_t _index_v1c2 = _bse_ctotal * _v1 + _c2;
                            
                            
                            
                            _eh_d2( _index_v1c1 , _index_v2c2 ) = 2.0 * _storage_prod( _index_c1v2 , _index_v1c2 ); 

                            
                           
                            
                        }
                    }
                }
            }
                        cout << " ==== check 5 ====== " << endl;
            
             for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                    
                    
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_v1c1 = _bse_ctotal * _v1 + _c1 ;

                        //// OK?

                        
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_v2c2 = _bse_ctotal * _v2 + _c2 ;
                         
            
            
            
             cout << " eh_d2 : " << _index_v1c1 << " : " << _index_v2c2 << " = " <<  _eh_d2( _index_v1c1 , _index_v2c2 ) << endl;
            
            
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
