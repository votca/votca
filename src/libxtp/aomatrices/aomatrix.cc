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
#include <votca/tools/linalg.h>

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>

#include <vector>



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
    
    void AOMatrix::Fill( AOBasis* aobasis, vec r, AOBasis* ecp ) {
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
                //cout << _submatrix<<endl;
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
        _trafo(0,0) = 1.0*contractions[0]; //  // s  Y 0,0
       // p-functions
        if ( _lmax > 0 ){ // order of functions changed
          double factor = 2.*sqrt(_decay)*contractions[1];
          _trafo(1,3) = factor;  // Y 1,0
          _trafo(2,2) = factor;  // Y 1,-1
          _trafo(3,1) = factor;  // Y 1,1
        }

        // d-functions
        if ( _lmax > 1 ){ // order of functions changed
          double factor = 2.*_decay*contractions[2];
          double factor_1 =  factor/sqrt(3.);
          _trafo(4,Cart::xx) = -factor_1;    // d3z2-r2 (dxx)
          _trafo(4,Cart::yy) = -factor_1;    // d3z2-r2 (dyy)  Y 2,0
          _trafo(4,Cart::zz) = 2.*factor_1;  // d3z2-r2 (dzz)

          _trafo(5,Cart::yz) = 2.*factor;     // dyz           Y 2,-1

          _trafo(6,Cart::xz) = 2.*factor;     // dxz           Y 2,1

          _trafo(7,Cart::xy) = 2.*factor;     // dxy           Y 2,-2

          _trafo(8,Cart::xx) = factor;       // dx2-y2 (dxx)   Y 2,2
          _trafo(8,Cart::yy) = -factor;      // dx2-y2 (dzz)
        }
       
        // f-functions
        if ( _lmax > 2 ){ // order of functions changed
          double factor = 2.*pow(_decay,1.5)*contractions[3];
          double factor_1 = factor*2./sqrt(15.);
          double factor_2 = factor*sqrt(2.)/sqrt(5.);
          double factor_3 = factor*sqrt(2.)/sqrt(3.);

          _trafo(9,Cart::xxz) = -3.*factor_1;        // f1 (f??) xxz 13
          _trafo(9,Cart::yyz) = -3.*factor_1;        // f1 (f??) yyz 15        Y 3,0
          _trafo(9,Cart::zzz) = 2.*factor_1;         // f1 (f??) zzz 19

          _trafo(10,Cart::xxy) = -factor_2;          // f3 xxy 10
          _trafo(10,Cart::yyy) = -factor_2;          // f3 yyy 18   Y 3,-1
          _trafo(10,Cart::yzz) = 4.*factor_2;        // f3 yzz 16

          _trafo(11,Cart::xxx) = -factor_2;          // f2 xxx 17
          _trafo(11,Cart::xyy) = -factor_2;          // f2 xyy 11   Y 3,1
          _trafo(11,Cart::xzz) = 4.*factor_2;        // f2 xzz 14

          _trafo(12,Cart::xyz) = 4.*factor;          // f6 xyz 12     Y 3,-2

          _trafo(13,Cart::xxz) = 2.*factor;          // f7 (f??)   xxz   13
          _trafo(13,Cart::yyz) = -2.*factor;         // f7 (f??)   yyz   15   Y 3,2

          _trafo(14,Cart::xxy) = 3.*factor_3;        // f4 xxy 10
          _trafo(14,Cart::yyy) = -factor_3;          // f4 yyy 18   Y 3,-3

          _trafo(15,Cart::xxx) = factor_3;           // f5 (f??) xxx 17
          _trafo(15,Cart::xyy) = -3.*factor_3;       // f5 (f??) xyy 11     Y 3,3
        }

        // g-functions
        if ( _lmax > 3 ) {
          double factor = 2./sqrt(3.)*_decay*_decay*contractions[4];
          double factor_1 = factor/sqrt(35.);
          double factor_2 = factor*4./sqrt(14.);
          double factor_3 = factor*2./sqrt(7.);
          double factor_4 = factor*2.*sqrt(2.);

          _trafo(16,Cart::xxxx) = 3.*factor_1;   /// Y 4,0
          _trafo(16,Cart::xxyy) = 6.*factor_1;
          _trafo(16,Cart::xxzz) = -24.*factor_1;
          _trafo(16,Cart::yyyy) = 3.*factor_1;
          _trafo(16,Cart::yyzz) = -24.*factor_1;
          _trafo(16,Cart::zzzz) = 8.*factor_1;

          _trafo(17,Cart::xxyz) = -3.*factor_2;  /// Y 4,-1
          _trafo(17,Cart::yyyz) = -3.*factor_2;
          _trafo(17,Cart::yzzz) = 4.*factor_2;

          _trafo(18,Cart::xxxz) = -3.*factor_2;  /// Y 4,1
          _trafo(18,Cart::xyyz) = -3.*factor_2;
          _trafo(18,Cart::xzzz) = 4.*factor_2;

          _trafo(19,Cart::xxxy) = -2.*factor_3;  /// Y 4,-2
          _trafo(19,Cart::xyyy) = -2.*factor_3;
          _trafo(19,Cart::xyzz) = 12.*factor_3;

          _trafo(20,Cart::xxxx) = -factor_3;     /// Y 4,2
          _trafo(20,Cart::xxzz) = 6.*factor_3;
          _trafo(20,Cart::yyyy) = factor_3;
          _trafo(20,Cart::yyzz) = -6.*factor_3;

          _trafo(21,Cart::xxyz) = 3.*factor_4;   /// Y 4,-3
          _trafo(21,Cart::yyyz) = -factor_4;

          _trafo(22,Cart::xxxz) = factor_4;      /// Y 4,3
          _trafo(22,Cart::xyyz) = -3.*factor_4;

          _trafo(23,Cart::xxxy) = 4.*factor;     /// Y 4,-4
          _trafo(23,Cart::xyyy) = -4.*factor;

          _trafo(24,Cart::xxxx) = factor;        /// Y 4,4
          _trafo(24,Cart::xxyy) = -6.*factor;
          _trafo(24,Cart::yyyy) = factor;
        }
        return;
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
        else if ( _lmax == 1 ) { _block_size = 4  ;}  // p
        else if ( _lmax == 2 ) { _block_size = 10 ;}  // d
        else if ( _lmax == 3 ) { _block_size = 20 ;}  // f
        else if ( _lmax == 4 ) { _block_size = 35 ;}  // g
        else if ( _lmax == 5 ) { _block_size = 56 ;}  // h
        else if ( _lmax == 6 ) { _block_size = 84 ;}  // i
        else if ( _lmax == 7 ) { _block_size = 120 ;} // j //////
        else if ( _lmax == 8 ) { _block_size = 165 ;} // k //////
        else{
            cerr << "GetBlocksize for l greater 8 not implemented!" << flush;
            exit(1);
        }
        return _block_size;
    }
    
    
    
    
    
}}

