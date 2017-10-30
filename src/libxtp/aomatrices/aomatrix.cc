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

    void AOSuperMatrix::PrintIndexToFunction(const AOBasis& aobasis){
        for (AOBasis::AOShellIterator _row = aobasis.firstShell(); _row != aobasis.lastShell() ; _row++ ) {
            const AOShell* _shell_row = aobasis.getShell( _row );
            int _row_start = _shell_row->getStartIndex();
            string type = _shell_row->getType();
            cout << "Shell " << type << "starts at " << _row_start+1 << endl;
        }
        return;
    }
    
    void AOMatrix::Fill(const AOBasis& aobasis,vec r, AOBasis* ecp ) {
        _aomatrix = ub::zero_matrix<double>(aobasis.AOBasisSize());
        _gridpoint = r;
        // loop row
        #pragma omp parallel for
        for (unsigned _row = 0; _row <  aobasis.getNumofShells() ; _row++ ){
       
            const AOShell* _shell_row = aobasis.getShell( _row );
            int _row_start = _shell_row->getStartIndex();
            int _row_end   = _row_start + _shell_row->getNumFunc();
           
            // AOMatrix is symmetric, restrict explicit calculation to triangular matrix
            for ( unsigned _col = 0; _col <= _row ; _col++ ){

                const AOShell* _shell_col = aobasis.getShell( _col );
                
                // figure out the submatrix
                int _col_start = _shell_col->getStartIndex();
                int _col_end   = _col_start + _shell_col->getNumFunc();
                
                ub::matrix_range< ub::matrix<double> > _submatrix = ub::subrange(_aomatrix, _row_start, _row_end, _col_start, _col_end);
                
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
       
        return;
    }
    
    
    void AOMatrix3D::Fill(const AOBasis& aobasis ) {
        // cout << "I'm supposed to fill out the AO overlap matrix" << endl;
        _aomatrix.resize(3);
        for (int i = 0; i < 3 ; i++){
          _aomatrix[ i ] = ub::zero_matrix<double>(aobasis.AOBasisSize());
        }
        
        // loop row
        #pragma omp parallel for
        for ( unsigned _row = 0; _row <  aobasis.getNumofShells() ; _row++ ){
       
            const AOShell* _shell_row = aobasis.getShell( _row );
            int _row_start = _shell_row->getStartIndex();
            int _row_end   = _row_start + _shell_row->getNumFunc();

            // loop column
            for (AOBasis::AOShellIterator _col = aobasis.firstShell(); _col != aobasis.lastShell() ; _col++ ) {
                const AOShell* _shell_col = aobasis.getShell( _col );
                
                // figure out the submatrix
                int _col_start = _shell_col->getStartIndex();
                int _col_end   = _col_start + _shell_col->getNumFunc();
                std::vector< ub::matrix_range< ub::matrix<double> > > _submatrix;
                for ( int _i = 0; _i < 3; _i++){
                   _submatrix.push_back(   ub::subrange(_aomatrix[_i], _row_start, _row_end, _col_start, _col_end) );
               
                }
                // Fill block
                FillBlock( _submatrix, _shell_row, _shell_col);

            }
        }
        return;
    }
    
    void AOMatrix3D::Cleanup(){
        
        for (int i = 0; i < 3; i++){
            
            _aomatrix[i].resize(0,0);
            
        }
        _aomatrix.clear();
       return; 
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
       
       ub::matrix<double> AOSuperMatrix::getTrafo(const AOGaussianPrimitive& gaussian){
         ///         0    1  2  3    4  5  6  7  8  9   10  11  12  13  14  15  16  17  18  19       20    21    22    23    24    25    26    27    28    29    30    31    32    33    34 
        ///         s,   x, y, z,   xy xz yz xx yy zz, xxy xyy xyz xxz xzz yyz yzz xxx yyy zzz,    xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyz, yyzz, yzzz, xxxx, yyyy, zzzz,
         const AOShell* shell=gaussian.getShell();
         const int ntrafo = shell->getNumFunc() + shell->getOffset();
         const double _decay=gaussian.getDecay();
         const int _lmax=shell->getLmax();
         const int n=getBlockSize( _lmax );
         ub::matrix<double> _trafo=ub::zero_matrix<double>(ntrafo,n); 
         const std::vector<double>& contractions=gaussian.getContraction();
         
        // s-functions
        _trafo(0,0) = contractions[0]; //  // s  Y 0,0
       // p-functions
        if ( _lmax > 0 ){ // order of functions changed
          const double factor = 2.*sqrt(_decay)*contractions[1];
          _trafo(1,3) = factor;  // Y 1,0
          _trafo(2,2) = factor;  // Y 1,-1
          _trafo(3,1) = factor;  // Y 1,1
        }

        // d-functions
        if ( _lmax > 1 ){ // order of functions changed
          const double factor = 2.*_decay*contractions[2];
          const double factor_1 =  factor/sqrt(3.);
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
          const double factor = 2.*pow(_decay,1.5)*contractions[3];
          const double factor_1 = factor*2./sqrt(15.);
          const double factor_2 = factor*sqrt(2.)/sqrt(5.);
          const double factor_3 = factor*sqrt(2.)/sqrt(3.);

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
          const double factor = 2./sqrt(3.)*_decay*_decay*contractions[4];
          const double factor_1 = factor/sqrt(35.);
          const double factor_2 = factor*4./sqrt(14.);
          const double factor_3 = factor*2./sqrt(7.);
          const double factor_4 = factor*2.*sqrt(2.);

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
        // h-functions
        if (_lmax > 4) {
          const double factor = (2./3.)*pow(_decay,2.5)*contractions[5];
          const double factor_1 = factor*2./sqrt(105.);
          const double factor_2 = factor*2./sqrt(7.);
          const double factor_3 = factor*sqrt(6.)/3.;
          const double factor_4 = factor*2.*sqrt(3.);
          const double factor_5 = factor*.2*sqrt(30.);

          _trafo(25,Cart::xxxxz) = 15.*factor_1;      /// Y 5,0
          _trafo(25,Cart::xxyyz) = 30.*factor_1;
          _trafo(25,Cart::xxzzz) = -40.*factor_1;
          _trafo(25,Cart::yyyyz) = 15.*factor_1;
          _trafo(25,Cart::yyzzz) = -40.*factor_1;
          _trafo(25,Cart::zzzzz) = 8.*factor_1;

          _trafo(26,Cart::xxxxy) = factor_2;          /// Y 5,-1
          _trafo(26,Cart::xxyyy) = 2.*factor_2;
          _trafo(26,Cart::xxyzz) = -12.*factor_2;
          _trafo(26,Cart::yyyyy) = factor_2;
          _trafo(26,Cart::yyyzz) = -12.*factor_2;
          _trafo(26,Cart::yzzzz) = 8.*factor_2;

          _trafo(27,Cart::xxxxx) = factor_2;          /// Y 5,1
          _trafo(27,Cart::xxxyy) = 2.*factor_2;
          _trafo(27,Cart::xxxzz) = -12.*factor_2;
          _trafo(27,Cart::xyyyy) = factor_2;
          _trafo(27,Cart::xyyzz) = -12.*factor_2;
          _trafo(27,Cart::xzzzz) = 8.*factor_2;

          _trafo(28,Cart::xxxyz) = -8.*factor;        /// Y 5,-2
          _trafo(28,Cart::xyyyz) = -8.*factor;
          _trafo(28,Cart::xyzzz) = 16.*factor;

          _trafo(29,Cart::xxxxz) = -4.*factor;        /// Y 5,2
          _trafo(29,Cart::xxzzz) = 8.*factor;
          _trafo(29,Cart::yyyyz) = 4.*factor;
          _trafo(29,Cart::yyzzz) = -8.*factor;

          _trafo(30,Cart::xxxxy) = -3.*factor_3;      /// Y 5,-3
          _trafo(30,Cart::xxyyy) = -2.*factor_3;
          _trafo(30,Cart::xxyzz) = 24.*factor_3;
          _trafo(30,Cart::yyyyy) = factor_3;
          _trafo(30,Cart::yyyzz) = -8.*factor_3;

          _trafo(31,Cart::xxxxx) = -factor_3;         /// Y 5,3
          _trafo(31,Cart::xxxyy) = 2.*factor_3;
          _trafo(31,Cart::xxxzz) = 8.*factor_3;
          _trafo(31,Cart::xyyyy) = 3.*factor_3;
          _trafo(31,Cart::xyyzz) = -24.*factor_3;

          _trafo(32,Cart::xxxyz) = 4.*factor_4;       /// Y 5,-4
          _trafo(32,Cart::xyyyz) = -4.*factor_4;

          _trafo(33,Cart::xxxxz) = factor_4;          /// Y 5,4
          _trafo(33,Cart::xxyyz) = -6.*factor_4;
          _trafo(33,Cart::yyyyz) = factor_4;

          _trafo(34,Cart::xxxxy) = 5.*factor_5;       /// Y 5,-5
          _trafo(34,Cart::xxyyy) = -10.*factor_5;
          _trafo(34,Cart::yyyyy) = factor_5;

          _trafo(35,Cart::xxxxx) = factor_5;          /// Y 5,5
          _trafo(35,Cart::xxxyy) = -10.*factor_5;
          _trafo(35,Cart::xyyyy) = 5.*factor_5;
        }

        // i-functions
        if (_lmax > 5) {
          const double factor = (2./3.)*_decay*_decay*_decay*contractions[6];
          const double factor_1 = factor*2./sqrt(1155.);
          const double factor_2 = factor*4./sqrt(55.);
          const double factor_3 = factor*sqrt(22.)/11.;
          const double factor_4 = factor*2.*sqrt(165.)/55.;
          const double factor_5 = factor*.4*sqrt(30.);
          const double factor_6 = factor*.2*sqrt(10.);

          _trafo(36,Cart::xxxxxx) = -5.*factor_1;     /// Y 6,0
          _trafo(36,Cart::xxxxyy) = -15.*factor_1;
          _trafo(36,Cart::xxxxzz) = 90.*factor_1;
          _trafo(36,Cart::xxyyyy) = -15.*factor_1;
          _trafo(36,Cart::xxyyzz) = 180.*factor_1;
          _trafo(36,Cart::xxzzzz) = -120.*factor_1;
          _trafo(36,Cart::yyyyyy) = -5.*factor_1;
          _trafo(36,Cart::yyyyzz) = 90.*factor_1;
          _trafo(36,Cart::yyzzzz) = -120.*factor_1;
          _trafo(36,Cart::zzzzzz) = 16.*factor_1;

          _trafo(37,Cart::xxxxyz) = 5.*factor_2;      /// Y 6,-1
          _trafo(37,Cart::xxyyyz) = 10.*factor_2;
          _trafo(37,Cart::xxyzzz) = -20.*factor_2;
          _trafo(37,Cart::yyyyyz) = 5.*factor_2;
          _trafo(37,Cart::yyyzzz) = -20.*factor_2;
          _trafo(37,Cart::yzzzzz) = 8.*factor_2;

          _trafo(38,Cart::xxxxxz) = 5.*factor_2;      /// Y 6,1
          _trafo(38,Cart::xxxyyz) = 10.*factor_2;
          _trafo(38,Cart::xxxzzz) = -20.*factor_2;
          _trafo(38,Cart::xyyyyz) = 5.*factor_2;
          _trafo(38,Cart::xyyzzz) = -20.*factor_2;
          _trafo(38,Cart::xzzzzz) = 8.*factor_2;

          _trafo(39,Cart::xxxxxy) = 2.*factor_3;      /// Y 6,-2
          _trafo(39,Cart::xxxyyy) = 4.*factor_3;
          _trafo(39,Cart::xxxyzz) = -32.*factor_3;
          _trafo(39,Cart::xyyyyy) = 2.*factor_3;
          _trafo(39,Cart::xyyyzz) = -32.*factor_3;
          _trafo(39,Cart::xyzzzz) = 32.*factor_3;

          _trafo(40,Cart::xxxxxy) = factor_3;         /// Y 6,2
          _trafo(40,Cart::xxxxyy) = factor_3;
          _trafo(40,Cart::xxxxzz) = -16.*factor_3;
          _trafo(40,Cart::xxyyyy) = -factor_3;
          _trafo(40,Cart::xxzzzz) = 16.*factor_3;
          _trafo(40,Cart::yyyyyy) = -factor_3;
          _trafo(40,Cart::yyyyzz) = 16.*factor_3;
          _trafo(40,Cart::yyzzzz) = -16.*factor_3;

          _trafo(41,Cart::xxxxyz) = -18.*factor_3;    /// Y 6,-3
          _trafo(41,Cart::xxyyyz) = -12.*factor_3;
          _trafo(41,Cart::xxyzzz) = 48.*factor_3;
          _trafo(41,Cart::yyyyyz) = 6.*factor_3;
          _trafo(41,Cart::yyyzzz) = -16.*factor_3;

          _trafo(42,Cart::xxxxxz) = -6.*factor_3;     /// Y 6,3
          _trafo(42,Cart::xxxyyz) = 12.*factor_3;
          _trafo(42,Cart::xxxzzz) = 16.*factor_3;
          _trafo(42,Cart::xyyyyz) = 18.*factor_3;
          _trafo(42,Cart::xyyzzz) = -48.*factor_3;

          _trafo(43,Cart::xxxxxy) = -4.*factor_4;     /// Y 6,-4
          _trafo(43,Cart::xxxyzz) = 40.*factor_4;
          _trafo(43,Cart::xyyyyy) = 4.*factor_4;
          _trafo(43,Cart::xyyyzz) = -40.*factor_4;

          _trafo(44,Cart::xxxxxx) = -factor_4;        /// Y 6,4
          _trafo(44,Cart::xxxxyy) = 5.*factor_4;
          _trafo(44,Cart::xxxxzz) = 10.*factor_4;
          _trafo(44,Cart::xxyyyy) = 5.*factor_4;
          _trafo(44,Cart::xxyyzz) = -60.*factor_4;
          _trafo(44,Cart::yyyyyy) = -factor_4;
          _trafo(44,Cart::yyyyzz) = 10.*factor_4;

          _trafo(45,Cart::xxxxyz) = 5.*factor_5;      /// Y 6,-5
          _trafo(45,Cart::xxyyyz) = -10.*factor_5;
          _trafo(45,Cart::yyyyyz) = factor_5;

          _trafo(46,Cart::xxxxxz) = factor_5;         /// Y 6,5
          _trafo(46,Cart::xxxyyz) = -10.*factor_5;
          _trafo(46,Cart::xyyyyz) = 5.*factor_5;

          _trafo(47,Cart::xxxxxy) = 6.*factor_6;      /// Y 6,-6
          _trafo(47,Cart::xxxyyy) = -20.*factor_6;
          _trafo(47,Cart::xyyyyy) = 6.*factor_6;

          _trafo(48,Cart::xxxxxx) = factor_6;         /// Y 6,6
          _trafo(48,Cart::xxxxyy) = -15.*factor_6;
          _trafo(48,Cart::xxyyyy) = 15.*factor_6;
          _trafo(48,Cart::yyyyyy) = -factor_6;
        }
        return _trafo;
       }
       
    

    std::vector<double> AOMatrix::XIntegrate(int _n, double _T  ){
        std::vector<double> _FmT=std::vector<double>(_n,0.0);
        const int _mm = _FmT.size() - 1;
        const double pi = boost::math::constants::pi<double>();
        if ( _mm < 0){
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
        
        return _FmT;
    }
    
    
    int AOSuperMatrix::getBlockSize( int _lmax){
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

