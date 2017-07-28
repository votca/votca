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
            _eh_qp = ub::zero_matrix<real_gwbse>( _bse_size , _bse_size );
            BSE_Add_qp2H( _eh_qp );
            return;
        }
        
    


        void GWBSE::BSE_solve_triplets(){
            
            // add full QP Hamiltonian contributions to free transitions
           
            ub::matrix<real_gwbse> _bse=_eh_d;
            
            linalg_eigenvalues(  _bse, _bse_triplet_energies, _bse_triplet_coefficients, _bse_nmax);
            return;
        }
        
        
        void GWBSE::Solve_nonhermitian(ub::matrix<double>& H, ub::matrix<double>& LT) {

            // remove stuff from Cholesky and Calculated L^T,, because more efficient for mat prods 
            #pragma omp parallel for
            for (unsigned i = 0; i < LT.size1(); i++) {
                for (unsigned j = i + 1; j < LT.size1(); j++) {
                    LT(i, j) = LT(j, i);
                    LT(j, i) = 0;
                }
            }

            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Removed non referenced part of Cholesky decompostion" << flush;
            // determine H = L^T(A-B)L
            ub::matrix<double> _temp = ub::prod(H, ub::trans(LT));
            H= ub::prod(LT, _temp);
            _temp.resize(0, 0);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Calculated H = L^T(A+B)L " << flush;

            // solve eigenvalue problem: HR_l = eps_l^2 R_l
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;

            linalg_eigenvalues(H, _eigenvalues, _eigenvectors, _bse_nmax);
            CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Solved HR_l = eps_l^2 R_l " << flush;
            H.resize(0, 0);
            // reconstruct real eigenvalues eps_l = sqrt(eps_l^2)
            _bse_singlet_energies.resize(_bse_nmax);
            for (int _i = 0; _i < _bse_nmax; _i++) {
                _bse_singlet_energies(_i) = sqrt(_eigenvalues(_i)); // store positive energies in orbitals objects
            }


            // reconstruct real eigenvectors X_l = 1/2 [sqrt(eps_l) (L^T)^-1 + 1/sqrt(eps_l)L ] R_l
            //                               Y_l = 1/2 [sqrt(eps_l) (L^T)^-1 - 1/sqrt(eps_l)L ] R_l

            // determine inverse of L^T
            ub::matrix<double> _cholesky_transposed_invert;
            ub::matrix<double> L = ub::trans(LT);
            linalg_invert(LT, _cholesky_transposed_invert);

            int dim = L.size1();
            _bse_singlet_coefficients.resize(dim, _bse_nmax); // resonant part (_X_evec)
            _bse_singlet_coefficients_AR.resize(dim, _bse_nmax); // anti-resonant part (_Y_evec)


            for (int _i = 0; _i < _bse_nmax; _i++) {
                //real_gwbse sqrt_eval = sqrt(_eigenvalues(_i));
                double sqrt_eval = sqrt(_bse_singlet_energies(_i));
                // get l-th reduced EV
                ub::matrix<double> _reduced_evec = ub::project(_eigenvectors, ub::range(0, dim), ub::range(_i, _i + 1)); // potentially col<->row

                ub::matrix<double> _transform = 0.5 * (sqrt_eval * _cholesky_transposed_invert + 1.0 / sqrt_eval * L);
                ub::project(_bse_singlet_coefficients, ub::range(0, dim), ub::range(_i, _i + 1)) = ub::prod(_transform, _reduced_evec);
                _transform = 0.5 * (sqrt_eval * _cholesky_transposed_invert - 1.0 / sqrt_eval * L);
                ub::project(_bse_singlet_coefficients_AR, ub::range(0, dim), ub::range(_i, _i + 1)) = ub::prod(_transform, _reduced_evec);
            }
            return;
        }
        
      void GWBSE::BSE_solve_singlets_BTDA(){
        
          
        // For details of the method, see EPL,78(2007)12001,
        // Nuclear Physics A146(1970)449, Nuclear Physics A163(1971)257.
        
          // setup resonant (A) and RARC blocks (B)
          // TOCHECK: Isn't that memory overkill here? _A and _B are never needed again?
           // ub::matrix<real_gwbse> _A = _eh_d + 2.0 * _eh_x;
           // ub::matrix<real_gwbse> _B = _eh_d2 + 2.0 * _eh_x;
          ub::matrix<double> _ApB = _eh_d + _eh_d2 + 4.0 * _eh_x;
          ub::matrix<double> _AmB = _eh_d - _eh_d2;

            
        
            
          // calculate Cholesky decomposition of A-B = LL^T. It throws an error if not positive definite
            //(A-B) is not needed any longer and can be overwritten
          
          bool positive_definite=true;
          CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << " Trying Cholesky decomposition of KAA-KAB" << flush;
          try
            {
            linalg_cholesky_decompose( _AmB );
            }
            catch (const std::runtime_error& error)
            {
                positive_definite=false;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() <<error.what()<<endl;
                CTP_LOG(ctp::logDEBUG, *_pLog) << ctp::TimeStamp() << "Trying Cholesky decomposition of KAA+KAB" << flush;
                linalg_cholesky_decompose( _ApB );
                _AmB = _eh_d - _eh_d2;
            }
            
          if(positive_definite){
              Solve_nonhermitian(_ApB,_AmB);
          }
          else{
              Solve_nonhermitian(_AmB,_ApB);
          }
          
          return;
      }
        
        
      void GWBSE::BSE_Add_qp2H( ub::matrix<real_gwbse>& qp ){
              
          #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c1;
                    // diagonal
                    qp( _index_vc , _index_vc ) += _vxc(_c1 + _bse_vtotal ,_c1 + _bse_vtotal ) - _vxc(_v1,_v1);
                    // v->c
                    for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                        size_t _index_vc2 = _bse_ctotal * _v1 + _c2;
                        if ( _c1 != _c2 ){
                            qp( _index_vc , _index_vc2 ) += _vxc(_c1+ _bse_vtotal ,_c2 + _bse_vtotal );
                        }
                    }
                    
                    // c-> v
                    for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                        size_t _index_vc2 = _bse_ctotal * _v2 + _c1;
                        if ( _v1 != _v2 ){
                            qp( _index_vc , _index_vc2 ) -= _vxc(_v1,_v2);
                        }
                    } 
                }
            }
            return;
      }
   
        
      void GWBSE::BSE_solve_singlets(){
            
            ub::matrix<real_gwbse> _bse = _eh_d + 2.0 * _eh_x;
            
            // _bse_singlet_energies.resize(_bse_singlet_coefficients.size1());
            linalg_eigenvalues(_bse, _bse_singlet_energies, _bse_singlet_coefficients, _bse_nmax);
            return;
        } 
        
        
        void GWBSE::BSE_d_setup ( TCMatrix& _Mmn){
            // gwbasis size
            size_t _gwsize = _Mmn[_homo].size1();

            // messy procedure, first get two matrices for occ and empty subbparts
            // store occs directly transposed
            ub::matrix<real_gwbse> _storage_v = ub::zero_matrix<real_gwbse>(  _bse_vtotal * _bse_vtotal , _gwsize );
            #pragma omp parallel for
            for ( size_t _v1 = 0; _v1 < _bse_vtotal; _v1++){
                const ub::matrix<real_gwbse>& Mmn = _Mmn[_v1 + _bse_vmin ];
                for ( size_t _v2 = 0; _v2 < _bse_vtotal; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_v( _index_vv , _i_gw ) = Mmn( _i_gw , _v2 + _bse_vmin );
                    }
                }
            }
            
            
            ub::matrix<real_gwbse> _storage_c = ub::zero_matrix<real_gwbse>( _gwsize, _bse_ctotal * _bse_ctotal );
            #pragma omp parallel for
            for ( size_t _c1 = 0; _c1 < _bse_ctotal; _c1++){
                const ub::matrix<real_gwbse>& Mmn = _Mmn[_c1 + _bse_cmin];
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
            ub::matrix<real_gwbse> _storage_prod = ub::prod( _storage_v , _storage_c );
            

            // now patch up _storage for screened interaction
            #pragma omp parallel for
            for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){  
                if (_ppm_weight(_i_gw) < 1.e-9) {                    
                    for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_vtotal) ; _v++){
                        _storage_v( _v , _i_gw ) = 0;
                    }
                    for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_ctotal) ; _c++){
                        _storage_c( _i_gw , _c ) =0;
                    }
                
                }else{
                    double _ppm_factor = sqrt( _ppm_weight( _i_gw ));
                    for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_vtotal) ; _v++){
                        _storage_v( _v , _i_gw ) = _ppm_factor * _storage_v(_v , _i_gw );
                    }
                    for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_ctotal) ; _c++){
                        _storage_c( _i_gw , _c ) = _ppm_factor * _storage_c( _i_gw , _c  );
                    }
                }
            }
            
            // multiply and subtract from _storage_prod
         
            _storage_prod -= ub::prod( _storage_v , _storage_c );
            
            // free storage_v and storage_c
            _storage_c.resize(0,0);
            _storage_v.resize(0,0);
            
            // finally resort into _eh_d
            // can be limited to upper diagonal !
            _eh_d = ub::zero_matrix<real_gwbse>( _bse_size , _bse_size );
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){
                    size_t _index_vv = _bse_vtotal * _v1 + _v2;
                    
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_vc1 = _bse_ctotal * _v1 + _c1 ;
                              
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_vc2 = _bse_ctotal * _v2 + _c2 ;
                            size_t _index_cc  = _bse_ctotal * _c1 + _c2;

                            _eh_d( _index_vc1 , _index_vc2 ) = -_storage_prod( _index_vv , _index_cc ); 
                        }
                    }
                }
            }
            
            return;
        }
        
        
         void GWBSE::BSE_d2_setup ( TCMatrix& _Mmn){
            // gwbasis size
            size_t _gwsize = _Mmn[_homo].size1();

            // messy procedure, first get two matrices for occ and empty subbparts
            // store occs directly transposed
            ub::matrix<real_gwbse> _storage_cv = ub::zero_matrix<real_gwbse>(  _bse_vtotal * _bse_ctotal , _gwsize );
            #pragma omp parallel for
            for ( size_t _c1 = 0; _c1 < _bse_ctotal; _c1++){
                const ub::matrix<real_gwbse>& Mmn = _Mmn[_c1 + _bse_cmin ];
                for ( size_t _v2 = 0; _v2 < _bse_vtotal; _v2++){
                    size_t _index_cv = _bse_vtotal * _c1 + _v2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_cv( _index_cv , _i_gw ) = Mmn( _i_gw , _v2 + _bse_vmin );
                    }
                }
            }
         
            ub::matrix<real_gwbse> _storage_vc = ub::zero_matrix<real_gwbse>( _gwsize, _bse_vtotal * _bse_ctotal );
            #pragma omp parallel for
            for ( size_t _v1 = 0; _v1 < _bse_vtotal; _v1++){
                const ub::matrix<real_gwbse>& Mmn = _Mmn[_v1 + _bse_vmin];
                for ( size_t _c2 = 0; _c2 < _bse_ctotal; _c2++){
                    size_t _index_vc = _bse_ctotal * _v1 + _c2;
                    for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++) {
                        _storage_vc( _i_gw , _index_vc ) = Mmn( _i_gw , _c2 + _bse_cmin );
                    }
                }
            }
            
            if ( ! _do_bse_singlets )  _Mmn.Cleanup();
            
            // store elements in a vtotal^2 x ctotal^2 matrix
            ub::matrix<real_gwbse> _storage_prod = ub::prod( _storage_cv , _storage_vc );
       
            
            // now patch up _storage for screened interaction
            #pragma omp parallel for
            for ( size_t _i_gw = 0 ; _i_gw < _gwsize ; _i_gw++ ){  
                if (_ppm_weight(_i_gw) < 1.e-9) {
                    for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_ctotal) ; _v++){
                    _storage_vc(  _i_gw , _v ) =0;
                }
                for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_vtotal) ; _c++){
                    _storage_cv( _c, _i_gw  ) = 0;
                }
                }else{
                double _ppm_factor = sqrt( _ppm_weight( _i_gw ));
                for ( size_t _v = 0 ; _v < (_bse_vtotal* _bse_ctotal) ; _v++){
                    _storage_vc(  _i_gw , _v ) = _ppm_factor * _storage_vc( _i_gw , _v );
                }
                for ( size_t _c = 0 ; _c < (_bse_ctotal*_bse_vtotal) ; _c++){
                    _storage_cv( _c, _i_gw  ) = _ppm_factor * _storage_cv( _c , _i_gw  );
                }
                    }
            }
         
            // multiply and subtract from _storage_prod
            _storage_prod -= ub::prod( _storage_cv , _storage_vc );
            
            // free storage_v and storage_c
            _storage_cv.resize(0,0);
            _storage_vc.resize(0,0);
            // finally resort into _eh_d
            // can be limited to upper diagonal !
            _eh_d2 = ub::zero_matrix<real_gwbse>( _bse_size , _bse_size );
            #pragma omp parallel for
            for ( size_t _v1 = 0 ; _v1 < _bse_vtotal ; _v1++){
                for ( size_t _v2 = 0 ; _v2 < _bse_vtotal ; _v2++){ 
                    for ( size_t _c1 = 0 ; _c1 < _bse_ctotal ; _c1++){
                        size_t _index_v1c1 = _bse_ctotal * _v1 + _c1 ;

                        size_t _index_c1v2 =_bse_vtotal * _c1 + _v2;
                        
                        for ( size_t _c2 = 0 ; _c2 < _bse_ctotal ; _c2++){
                            size_t _index_v2c2 = _bse_ctotal * _v2 + _c2 ;
                            size_t _index_v1c2 = _bse_ctotal * _v1 + _c2;

                            _eh_d2( _index_v1c1 , _index_v2c2 ) = -_storage_prod( _index_c1v2 , _index_v1c2 ); 

                        }
                    }
                }
            }
         return;   
        }
        
        
        
        void GWBSE::BSE_x_setup( TCMatrix& _Mmn){
            
            /* unlike the fortran code, we store eh interaction directly in
             * a suitable matrix form instead of a four-index array
             */
                        
            // gwbasis size
            size_t _gwsize = _Mmn[_homo].size1();
            
            // get a different storage for 3-center integrals we need
            //cout<< "Starting to set up "<< endl;
            ub::matrix<real_gwbse> _storage = ub::zero_matrix<real_gwbse>( _gwsize , _bse_size);
            //cout<< "Storage set up"<< endl;
         
            // occupied levels
            #pragma omp parallel for
            for ( size_t _v = 0; _v < _bse_vtotal ; _v++ ){
                // cout << " act threads: " << omp_get_thread_num( ) << " total threads " << omp_get_num_threads( ) << " max threads " << omp_get_max_threads( ) <<endl;
                ub::matrix<real_gwbse>& Mmn = _Mmn[_v + _bse_vmin];
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
            return;    
        }
        
        
        

    }
    
 
};
