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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/xtp/votca_xtp_config.h>

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>
#include <votca/xtp/logger.h>
#include <votca/tools/linalg.h>
// #include <omp.h>

using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;

    void AOSuperMatrix::PrintIndexToFunction( AOBasis* aobasis){
        for (vector< AOShell* >::iterator _row = aobasis->firstShell(); _row != aobasis->lastShell() ; _row++ ) {
            AOShell* _shell_row = aobasis->getShell( _row );
            int _row_start = _shell_row->getStartIndex();
            string type = _shell_row->getType();
            cout << "Shell " << type << "starts at " << _row_start+1 << endl;
        }
    }
    
    void AOMatrix::Fill( AOBasis* aobasis, ub::vector<double> r, AOBasis* ecp ) {
        // cout << "I'm supposed to fill out the AO overlap matrix" << endl;

          //      cout << aobasis->_aoshells.size();
      
        _gridpoint = r;
        // loop row
        #pragma omp parallel for
        for (unsigned _row = 0; _row <  aobasis->_aoshells.size() ; _row++ ){
        //for (vector< AOShell* >::iterator _row = aobasis->firstShell(); _row != aobasis->lastShell() ; _row++ ) {
            //cout << " act threads: " << omp_get_thread_num( ) << " total threads " << omp_get_num_threads( ) << " max threads " << omp_get_max_threads( ) <<endl;
            AOShell* _shell_row = aobasis->getShell( _row );
            int _row_start = _shell_row->getStartIndex();
            int _row_end   = _row_start + _shell_row->getNumFunc();
           
            // AOMatrix is symmetric, restrict explicit calculation to triangular matrix
            for ( unsigned _col = 0; _col <= _row ; _col++ ){

                AOShell* _shell_col = aobasis->getShell( _col );
                
                // figure out the submatrix
                int _col_start = _shell_col->getStartIndex();
                int _col_end   = _col_start + _shell_col->getNumFunc();
                //cout << _row << ":" << _row_start << ":" << _row_end << "/" << _col << ":" <<  _col_start << ":" << _col_end << endl;
                ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange(this->_aomatrix, _row_start, _row_end, _col_start, _col_end);

                // Fill block
                FillBlock( _submatrix, _shell_row, _shell_col, ecp );

            }
        }
        
        // Fill whole matrix by copying
        for ( unsigned _i=0; _i < _aomatrix.size1(); _i++){
            for ( unsigned _j=0; _j < _i; _j++){
               _aomatrix(_j,_i) = _aomatrix(_i,_j); 
                       
            }
        }
     
 
        
      
        // check symmetry
         bool _is_symmetric = true;
        
        // Copy stuff to fill lower triangular part
         for ( unsigned _i=0; _i < this->_aomatrix.size1(); _i++){
            for (unsigned _j=0; _j <= _i; _j++){
         
                if ( std::abs(this->_aomatrix(_i,_j) - this->_aomatrix(_j,_i) ) > 1e-4 ) {
                    
                    cerr << _i << ":" << _j << " == " << this->_aomatrix(_i,_j) << " vs " <<  this->_aomatrix(_j,_i) << endl;
                    _is_symmetric = false;
                }
                
            }
        }
        if ( !_is_symmetric) {cerr << " Error: AOMatrix is not symmetric! "; exit(1);}
        
       
        
    }
    
    
    void AOMatrix3D::Fill( AOBasis* aobasis ) {
        // cout << "I'm supposed to fill out the AO overlap matrix" << endl;
        
        // loop row
        #pragma omp parallel for
        for ( unsigned _row = 0; _row <  aobasis->_aoshells.size() ; _row++ ){
        // for (vector< AOShell* >::iterator _row = aobasis->firstShell(); _row != aobasis->lastShell() ; _row++ ) {
            AOShell* _shell_row = aobasis->getShell( _row );
            int _row_start = _shell_row->getStartIndex();
            int _row_end   = _row_start + _shell_row->getNumFunc();

            // loop column
            for (vector< AOShell* >::iterator _col = aobasis->firstShell(); _col != aobasis->lastShell() ; _col++ ) {
                AOShell* _shell_col = aobasis->getShell( _col );
                
                // figure out the submatrix
                int _col_start = _shell_col->getStartIndex();
                int _col_end   = _col_start + _shell_col->getNumFunc();
                std::vector< ub::matrix_range< ub::matrix<double> > > _submatrix;
                for ( int _i = 0; _i < 3; _i++){
                   _submatrix.push_back(   ub::subrange(this->_aomatrix[_i], _row_start, _row_end, _col_start, _col_end) );
                //ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange(this->_aomatrix, _row_start, _row_end, _col_start, _col_end);
                }
                // Fill block
                FillBlock( _submatrix, _shell_row, _shell_col);

            }
        }
    }
    
    void AOMatrix3D::Cleanup(){
        
        for (int i = 0; i < 3; i++){
            
            _aomatrix[i].resize(0,0);
            
        }
        _aomatrix.clear();
        
    }
    
    
    
    
    
    
    
    
    
    
    //    AOMatrix::~AOMatrix() {
    // ;
      //_aomatrix.clear();
    //    _aomatrix.resize(0,0);
    // };
    
    
    void AOMatrix::Print( string _ident){
        cout << "\n" << endl;
        std::cout.precision(12);
        for ( unsigned i =0; i< this->_aomatrix.size1(); i++){
            for ( unsigned j =0; j< this->_aomatrix.size2(); j++){
                cout << _ident << "[" << i+1 << ":" << j+1 << "] " << scientific << this->_aomatrix(i,j) << endl;
            }
        }
    }
    
    
       void AOMatrix3D::Print( string _ident){
        cout << "\n" << endl;
        for ( unsigned i =0; i< this->_aomatrix[0].size1(); i++){
            for ( unsigned j =0; j< this->_aomatrix[0].size2(); j++){
                cout << _ident << "[" << i+1 << ":" << j+1 << "] " <<  this->_aomatrix[0](i,j) << " : " <<  this->_aomatrix[1](i,j) << " : " <<  this->_aomatrix[2](i,j)  << endl;
            }
        }
    }
       
       void AOSuperMatrix::getTrafo(ub::matrix<double>& _trafo, int _lmax, const double& _decay, std::vector<double> contractions){
        // s-functions
        _trafo(0,0) = 1.0*contractions[0]; // s
       
        // p-functions
        if ( _lmax > 0 ){
            //cout << _trafo_row.size1() << ":" << _trafo_row.size2() << endl;
            _trafo(1,1) = 2.0*sqrt(_decay)*contractions[1];
            _trafo(2,2) = 2.0*sqrt(_decay)*contractions[1];
            _trafo(3,3) = 2.0*sqrt(_decay)*contractions[1];
        }
        //votca order is dxz dyz dxy d3z2-r2 dx2-y2
        // d-functions
        if ( _lmax > 1 ){
            _trafo(4,5) = 4.0*_decay*contractions[2];             // dxz
            _trafo(5,6) = _trafo(4,5);            // dyz
            _trafo(6,4) = _trafo(4,5);            // dxy
            _trafo(7,7) = -2.0*_decay/sqrt(3.0)*contractions[2];  // d3z2-r2 (dxx)
            _trafo(7,8) = _trafo(7,7);            // d3z2-r2 (dyy)
            _trafo(7,9) = -2.0*_trafo(7,7);       // d3z2-r2 (dzz)
            _trafo(8,7) = 2.0*_decay*contractions[2];             // dx2-y2 (dxx)
            _trafo(8,8) = -_trafo(8,7);           // dx2-y2 (dzz)
        }
        
        // f-functions
        if ( _lmax > 2 ){
            _trafo(9,12) = 4.0 * 2.0 *pow(_decay,1.5)/sqrt(15.)*contractions[3]; // f1 (f??)
            _trafo(9,15) = -1.5 * _trafo(9,12);        // f1 (f??)
            _trafo(9,17) = _trafo(9,15);               // f1 (f??)
            
            _trafo(10,16) = 4.0 * 2.0 * sqrt(2.0)/sqrt(5.0) * pow(_decay,1.5)*contractions[3]; // f2 (f??)
            _trafo(10,10) = -0.25 * _trafo(10,16);                             // f2 f(??)
            _trafo(10,14) = _trafo(10,10);                                     // f2 f(??)
            
            _trafo(11,18) = _trafo(10,16);                                     // f3 (f??)
            _trafo(11,13) = -0.25 * _trafo(11,18);                             // f3 f(??)
            _trafo(11,11) = _trafo(11,13);                                     // f3 f(??)            
                   
            _trafo(12,13) = 3.0 * 2.0 * sqrt(2.0)/sqrt(3.0) * pow(_decay,1.5)*contractions[3]; // f4 (f??)
            _trafo(12,11) = -_trafo(12,13)/3.0;                                // f4 (f??)
            
            _trafo(13,10) = -_trafo(12,11);                                    // f5 (f??)
            _trafo(13,14) = -_trafo(12,13);                                    // f5 (f??)
            
            _trafo(14,19) = 8.0 * pow(_decay,1.5)*contractions[3];                             // f6 (f??)
            
            _trafo(15,15) = 0.5 * _trafo(14,19);                               // f7 (f??)
            _trafo(15,17) = -_trafo(15,15);                                    // f7 (f??)
        }
        
        // g-functions
        if ( _lmax > 3 ){
            _trafo(16,22) = 8.0 * 2.0/sqrt(105.0) * pow(_decay,2.0)*contractions[4];
            _trafo(16,21) = 3.0 * 2.0/sqrt(105.0) * pow(_decay,2.0)*contractions[4];
            _trafo(16,20) = _trafo(16,21);
            _trafo(16,29) = -3.0 * _trafo(16,22);
            _trafo(16,31) = 2.0 * _trafo(16,21);
            _trafo(16,30) = _trafo(16,29);
            _trafo(16,5)  = _trafo(16,31);
            
             /* vv(17,:) =  (/   23,  22, 21, 30, 32, 31,   6 /) ! g
                cc(17,:) =  (/    8,  3, 3, -24, 6, -24,    6 /)
                normConst(17,:) = (/ 2.d0/sqrt(105.d0) ,2.d0  /)
              */
            _trafo(17,26) = 4.0 * 4.0*sqrt(2.0)/sqrt(21.0) * pow(_decay,2.0)*contractions[4];
            _trafo(17,25) = -0.75 * _trafo(17,26);
            _trafo(17,33) = _trafo(17,25);
             
             /* vv(18,:) =  (/   27,  26, 34,  0,  0,  0,   3 /) ! g
                cc(18,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
                normConst(18,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
              */
            
            _trafo(18,28) = _trafo(17,26);
            _trafo(18,32) = _trafo(17,25);
            _trafo(18,27) = _trafo(17,25);
             
            /* vv(19,:) =  (/   29,  33, 28,  0,  0,  0,   3 /) ! g 
               cc(19,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
               normConst(19,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
             */
     
            _trafo(19,34) = 6.0 * 8.0/sqrt(21.0) * pow(_decay,2.0)*contractions[4];
            _trafo(19,23) = -_trafo(19,34)/6.0;
            _trafo(19,24) = _trafo(19,23);
             
            /* vv(20,:) =  (/   35,  24, 25,  0,  0,  0,   3 /) ! g
               cc(20,:) =  (/    6,  -1, -1,  0,  0,  0,   3 /)
               normConst(20,:) = (/ 8.d0/sqrt(21.d0) ,2.d0  /)
             */
    
            _trafo(20,29) = 6.0 * 4.0/sqrt(21.0) * pow(_decay,2.0)*contractions[4];
            _trafo(20,20) = -_trafo(20,29)/6.0;
            _trafo(20,30) = -_trafo(20,29);
            _trafo(20,21) = -_trafo(20,20);

            /* vv(21,:) =  (/   30,  21, 31, 22,  0,  0,   4 /) ! g
               cc(21,:) =  (/    6,  -1, -6, 1,  0,  0,    4 /)
               normConst(21,:) = (/ 4.d0/sqrt(21.d0) ,2.d0  /)
             */
    
            _trafo(21,25) = 4.0 * sqrt(2.0)/sqrt(3.0) * pow(_decay,2.0)*contractions[4];
            _trafo(21,33) = -3.0 * _trafo(21,25);
             
            /* vv(22,:) =  (/   26,  34,  0,  0,  0,  0,   2 /) ! g
               cc(22,:) =  (/    1,  -3,  0,  0,  0,  0,   2 /)
               normConst(22,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
             */
    
            _trafo(22,32) = -_trafo(21,33);
            _trafo(22,27) = -_trafo(21,25);
            
            /* vv(23,:) =  (/   33,  28,  0,  0,  0,  0,   2 /) ! g
               cc(23,:) =  (/    3,  -1,  0,  0,  0,  0,   2 /)
               normConst(23,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
             */
    
            _trafo(23,23) = 8.0/sqrt(3.0) * pow(_decay,2.0)*contractions[4];
            _trafo(23,24) = -_trafo(23,23);
             
            /* vv(24,:) =  (/   24,  25,  0,  0,  0,  0,   2 /) ! g 
               cc(24,:) =  (/    1,  -1,  0,  0,  0,  0,   2 /)
               normConst(24,:) = (/ 8.d0/sqrt(3.d0) ,2.d0  /)
             */
    
            _trafo(24,20) = 2.0/sqrt(3.0) * pow(_decay,2.0)*contractions[4];
            _trafo(24,21) = _trafo(24,20);
            _trafo(24,31) = -6.0 * _trafo(24,20);
             
            /* vv(25,:) =  (/   21,  22, 32,  0,  0,  0,   3 /) ! g
               cc(25,:) =  (/    1,  1, -6,  0,  0,  0,   3  /)
               normConst(25,:) = (/ 2.d0/sqrt(3.d0) ,2.d0  /)
             */
           
       
       }
       }
       
    void AOSuperMatrix::getTrafo(ub::matrix<double>& _trafo, int _lmax, const double& _decay) {
        // s-functions
        _trafo(0,0) = 1.0; // s
       
        // p-functions
        if ( _lmax > 0 ){
            //cout << _trafo_row.size1() << ":" << _trafo_row.size2() << endl;
            _trafo(1,1) = 2.0*sqrt(_decay);
            _trafo(2,2) = 2.0*sqrt(_decay);
            _trafo(3,3) = 2.0*sqrt(_decay);
        }

        // d-functions
        if ( _lmax > 1 ){
            _trafo(4,5) = 4.0*_decay;             // dxz
            _trafo(5,6) = _trafo(4,5);            // dyz
            _trafo(6,4) = _trafo(4,5);            // dxy
            _trafo(7,7) = -2.0*_decay/sqrt(3.0);  // d3z2-r2 (dxx)
            _trafo(7,8) = _trafo(7,7);            // d3z2-r2 (dyy)
            _trafo(7,9) = -2.0*_trafo(7,7);       // d3z2-r2 (dzz)
            _trafo(8,7) = 2.0*_decay;             // dx2-y2 (dxx)
            _trafo(8,8) = -_trafo(8,7);           // dx2-y2 (dzz)
        }
        
        // f-functions
        if ( _lmax > 2 ){
            _trafo(9,12) = 4.0 * 2.0 *pow(_decay,1.5)/sqrt(15.); // f1 (f??)
            _trafo(9,15) = -1.5 * _trafo(9,12);        // f1 (f??)
            _trafo(9,17) = _trafo(9,15);               // f1 (f??)
            
            _trafo(10,16) = 4.0 * 2.0 * sqrt(2.0)/sqrt(5.0) * pow(_decay,1.5); // f2 (f??)
            _trafo(10,10) = -0.25 * _trafo(10,16);                             // f2 f(??)
            _trafo(10,14) = _trafo(10,10);                                     // f2 f(??)
            
            _trafo(11,18) = _trafo(10,16);                                     // f3 (f??)
            _trafo(11,13) = -0.25 * _trafo(11,18);                             // f3 f(??)
            _trafo(11,11) = _trafo(11,13);                                     // f3 f(??)            
                   
            _trafo(12,13) = 3.0 * 2.0 * sqrt(2.0)/sqrt(3.0) * pow(_decay,1.5); // f4 (f??)
            _trafo(12,11) = -_trafo(12,13)/3.0;                                // f4 (f??)
            
            _trafo(13,10) = -_trafo(12,11);                                    // f5 (f??)
            _trafo(13,14) = -_trafo(12,13);                                    // f5 (f??)
            
            _trafo(14,19) = 8.0 * pow(_decay,1.5);                             // f6 (f??)
            
            _trafo(15,15) = 0.5 * _trafo(14,19);                               // f7 (f??)
            _trafo(15,17) = -_trafo(15,15);                                    // f7 (f??)
        }
        
        // g-functions
        if ( _lmax > 3 ){
            _trafo(16,22) = 8.0 * 2.0/sqrt(105.0) * pow(_decay,2.0);
            _trafo(16,21) = 3.0 * 2.0/sqrt(105.0) * pow(_decay,2.0);
            _trafo(16,20) = _trafo(16,21);
            _trafo(16,29) = -3.0 * _trafo(16,22);
            _trafo(16,31) = 2.0 * _trafo(16,21);
            _trafo(16,30) = _trafo(16,29);
            _trafo(16,5)  = _trafo(16,31);
            
             /* vv(17,:) =  (/   23,  22, 21, 30, 32, 31,   6 /) ! g
                cc(17,:) =  (/    8,  3, 3, -24, 6, -24,    6 /)
                normConst(17,:) = (/ 2.d0/sqrt(105.d0) ,2.d0  /)
              */
            _trafo(17,26) = 4.0 * 4.0*sqrt(2.0)/sqrt(21.0) * pow(_decay,2.0);
            _trafo(17,25) = -0.75 * _trafo(17,26);
            _trafo(17,33) = _trafo(17,25);
             
             /* vv(18,:) =  (/   27,  26, 34,  0,  0,  0,   3 /) ! g
                cc(18,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
                normConst(18,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
              */
            
            _trafo(18,28) = _trafo(17,26);
            _trafo(18,32) = _trafo(17,25);
            _trafo(18,27) = _trafo(17,25);
             
            /* vv(19,:) =  (/   29,  33, 28,  0,  0,  0,   3 /) ! g 
               cc(19,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
               normConst(19,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
             */
     
            _trafo(19,34) = 6.0 * 8.0/sqrt(21.0) * pow(_decay,2.0);
            _trafo(19,23) = -_trafo(19,34)/6.0;
            _trafo(19,24) = _trafo(19,23);
             
            /* vv(20,:) =  (/   35,  24, 25,  0,  0,  0,   3 /) ! g
               cc(20,:) =  (/    6,  -1, -1,  0,  0,  0,   3 /)
               normConst(20,:) = (/ 8.d0/sqrt(21.d0) ,2.d0  /)
             */
    
            _trafo(20,29) = 6.0 * 4.0/sqrt(21.0) * pow(_decay,2.0);
            _trafo(20,20) = -_trafo(20,29)/6.0;
            _trafo(20,30) = -_trafo(20,29);
            _trafo(20,21) = -_trafo(20,20);

            /* vv(21,:) =  (/   30,  21, 31, 22,  0,  0,   4 /) ! g
               cc(21,:) =  (/    6,  -1, -6, 1,  0,  0,    4 /)
               normConst(21,:) = (/ 4.d0/sqrt(21.d0) ,2.d0  /)
             */
    
            _trafo(21,25) = 4.0 * sqrt(2.0)/sqrt(3.0) * pow(_decay,2.0);
            _trafo(21,33) = -3.0 * _trafo(21,25);
             
            /* vv(22,:) =  (/   26,  34,  0,  0,  0,  0,   2 /) ! g
               cc(22,:) =  (/    1,  -3,  0,  0,  0,  0,   2 /)
               normConst(22,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
             */
    
            _trafo(22,32) = -_trafo(21,33);
            _trafo(22,27) = -_trafo(21,25);
            
            /* vv(23,:) =  (/   33,  28,  0,  0,  0,  0,   2 /) ! g
               cc(23,:) =  (/    3,  -1,  0,  0,  0,  0,   2 /)
               normConst(23,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
             */
    
            _trafo(23,23) = 8.0/sqrt(3.0) * pow(_decay,2.0);
            _trafo(23,24) = -_trafo(23,23);
             
            /* vv(24,:) =  (/   24,  25,  0,  0,  0,  0,   2 /) ! g 
               cc(24,:) =  (/    1,  -1,  0,  0,  0,  0,   2 /)
               normConst(24,:) = (/ 8.d0/sqrt(3.d0) ,2.d0  /)
             */
    
            _trafo(24,20) = 2.0/sqrt(3.0) * pow(_decay,2.0);
            _trafo(24,21) = _trafo(24,20);
            _trafo(24,31) = -6.0 * _trafo(24,20);
             
            /* vv(25,:) =  (/   21,  22, 32,  0,  0,  0,   3 /) ! g
               cc(25,:) =  (/    1,  1, -6,  0,  0,  0,   3  /)
               normConst(25,:) = (/ 2.d0/sqrt(3.d0) ,2.d0  /)
             */
               
        }
        
        
    }

    void AOMatrix::XIntegrate(vector<double>& _FmT, const double& _T  ){
        
        const int _mm = _FmT.size() - 1;
        const double pi = boost::math::constants::pi<double>();
        if ( _mm < 0 || _mm > 10){
            cerr << "mm is: " << _mm << " This should not have happened!" << flush;
            exit(1);
        }
        
        if ( _T < 0.0 ) {
            cerr << "T is: " << _T << " This should not have happened!" << flush;
            exit(1);
        }
  
        if ( _T >= 10.0 ) {
            // forward iteration
            _FmT[0]=0.50*sqrt(pi/_T)* erf(sqrt(_T));

            for (unsigned m = 1; m < _FmT.size(); m++ ){
                _FmT[m] = (2*m-1) * _FmT[m-1]/(2.0*_T) - exp(-_T)/(2.0*_T) ;
            }
        }

        if ( _T < 1e-10 ){
           for ( unsigned m=0; m < _FmT.size(); m++){
               _FmT[m] = 1.0/(2.0*m+1.0) - _T/(2.0*m+3.0); 
           }
        }

        
        if ( _T >= 1e-10 && _T < 10.0 ){
            // backward iteration
            double fm = 0.0;
            for ( int m = 60; m >= _mm; m--){
                fm = (2.0*_T)/(2.0*m+1.0) * ( fm + exp(-_T)/(2.0*_T));
            } 
            _FmT[_mm] = fm;
            for (int m = _mm-1 ; m >= 0; m--){
                _FmT[m] = (2.0*_T)/(2.0*m+1.0) * (_FmT[m+1] + exp(-_T)/(2.0*_T));
            }
        }
        

    }
    
    
    int AOSuperMatrix::getBlockSize(int _lmax){
        int _block_size;
        if ( _lmax == 0 ) { _block_size = 1  ;}  // s
        if ( _lmax == 1 ) { _block_size = 4  ;}  // p
        if ( _lmax == 2 ) { _block_size = 10 ;}  // d
        if ( _lmax == 3 ) { _block_size = 20 ;}  // f
        if ( _lmax == 4 ) { _block_size = 35 ;}  // g
        if ( _lmax == 5 ) { _block_size = 56 ;}  // h
        if ( _lmax <= 6 ) { _block_size = 84 ;}  // i
        
        return _block_size;
    }
    
    
    
    
    
}}

