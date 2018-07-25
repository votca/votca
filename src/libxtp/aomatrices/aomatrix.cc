/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>

#include <vector>



namespace votca {
    namespace xtp {

        void AOSuperMatrix::PrintIndexToFunction(const AOBasis& aobasis) {
            for ( const AOShell* shell_row:aobasis) {
                int _row_start = shell_row->getStartIndex();
                std::string type = shell_row->getType();
                std::cout << "Shell " << type << "starts at " << _row_start + 1 << std::endl;
            }
            return;
        }

        
        template< class T> 
        void AOMatrix<T>::Fill(const AOBasis& aobasis) {
            _aomatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
            // loop row
#pragma omp parallel for schedule(guided)
            for (unsigned _row = 0; _row < aobasis.getNumofShells(); _row++) {

                const AOShell* _shell_row = aobasis.getShell(_row);
                int _row_start = _shell_row->getStartIndex();
      

                // AOMatrix is symmetric, restrict explicit calculation to triangular matrix
                for (unsigned _col = _row; _col <  aobasis.getNumofShells(); _col++) {

                    const AOShell* _shell_col = aobasis.getShell(_col);

                    // figure out the submatrix
                    int _col_start = _shell_col->getStartIndex();
                    Eigen::Block< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > block=
                            _aomatrix.block(_row_start,_col_start, _shell_row->getNumFunc(),_shell_col->getNumFunc());
                    // Fill block
                    FillBlock(block, _shell_row, _shell_col);

                }
            }

            // Fill whole matrix by copying
        for ( unsigned _i=0; _i < _aomatrix.rows(); _i++){
                for (unsigned _j = 0; _j < _i; _j++) {
                    _aomatrix(_i, _j) = _aomatrix(_j, _i);

                }
            }

            return;
        }


        void AOMatrix3D::Fill(const AOBasis& aobasis) {
            // cout << "I'm supposed to fill out the AO overlap matrix" << endl;
            _aomatrix.resize(3);
            for (int i = 0; i < 3; i++) {
          _aomatrix[ i ] = Eigen::MatrixXd::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
            }
            // loop row
#pragma omp parallel for
            for (unsigned _row = 0; _row < aobasis.getNumofShells(); _row++) {

                const AOShell* shell_row = aobasis.getShell(_row);
                int _row_start = shell_row->getStartIndex();

                // loop column
                for (const AOShell* shell_col:aobasis) {
                    // figure out the submatrix
                    int _col_start = shell_col->getStartIndex();
                std::vector< Eigen::Block<Eigen::MatrixXd> > _submatrix;
                    for (int _i = 0; _i < 3; _i++) {
                   Eigen::Block<Eigen::MatrixXd> block=_aomatrix[_i].block( _row_start,_col_start,shell_row->getNumFunc(),shell_col->getNumFunc());
                   _submatrix.push_back(block );

                    }
                    // Fill block
                    FillBlock(_submatrix, shell_row, shell_col);

                }
            }
            return;
        }

        void AOMatrix3D::FreeMatrix() {

            for (int i = 0; i < 3; i++) {

                _aomatrix[i].resize(0, 0);

            }
            _aomatrix.clear();
            return;
        }

        template<class T> void AOMatrix<T>::Print(std::string _ident) {
            std::cout << "\n" << std::endl;
            std::cout.precision(12);
        for ( unsigned i =0; i< _aomatrix.rows(); i++){
            for ( unsigned j =0; j< _aomatrix.cols(); j++){
                std::cout << _ident << "[" << i+1 << ":" << j+1 << "] " << std::scientific <<_aomatrix(i,j) << std::endl;
                }
            }
    }   

        void AOMatrix3D::Print(std::string _ident) {
            std::cout << "\n" << std::endl;
        for ( unsigned i =0; i< _aomatrix[0].rows(); i++){
            for ( unsigned j =0; j< _aomatrix[0].cols(); j++){
                std::cout << _ident << "[" << i+1 << ":" << j+1 << "] " <<  _aomatrix[0](i,j) << " : " <<  _aomatrix[1](i,j) << " : " <<  _aomatrix[2](i,j)  << std::endl;
                }
            }
        }

       Eigen::MatrixXd AOSuperMatrix::getTrafo(const AOGaussianPrimitive& gaussian){
            ///         0    1  2  3    4  5  6  7  8  9   10  11  12  13  14  15  16  17  18  19       20    21    22    23    24    25    26    27    28    29    30    31    32    33    34 
            ///         s,   x, y, z,   xy xz yz xx yy zz, xxy xyy xyz xxz xzz yyz yzz xxx yyy zzz,    xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyz, yyzz, yzzz, xxxx, yyyy, zzzz,
            const AOShell* shell = gaussian.getShell();
            const int ntrafo = shell->getNumFunc() + shell->getOffset();
            const double _decay = gaussian.getDecay();
            const int _lmax = shell->getLmax();
            const int n = getBlockSize(_lmax);
         Eigen::MatrixXd _trafo=Eigen::MatrixXd::Zero(n,ntrafo); 
            const std::vector<double>& contractions = gaussian.getContraction();

            // s-functions
            _trafo(0, 0) = contractions[0]; //  // s  Y 0,0
            // p-functions
            if (_lmax > 0) { // order of functions changed
                const double factor = 2. * sqrt(_decay) * contractions[1];
                _trafo(3, 1) = factor; // Y 1,0
                _trafo(2, 2) = factor; // Y 1,-1
                _trafo(1, 3) = factor; // Y 1,1
            }

            // d-functions
            if (_lmax > 1) { // order of functions changed
                const double factor = 2. * _decay * contractions[2];
                const double factor_1 = factor / sqrt(3.);
                _trafo(Cart::xx,4) = -factor_1; // d3z2-r2 (dxx)
                _trafo(Cart::yy,4) = -factor_1; // d3z2-r2 (dyy)  Y 2,0
                _trafo(Cart::zz,4) = 2. * factor_1; // d3z2-r2 (dzz)

                _trafo(Cart::yz,5) = 2. * factor; // dyz           Y 2,-1

                _trafo(Cart::xz,6) = 2. * factor; // dxz           Y 2,1

                _trafo(Cart::xy,7) = 2. * factor; // dxy           Y 2,-2

                _trafo(Cart::xx,8) = factor; // dx2-y2 (dxx)   Y 2,2
                _trafo(Cart::yy,8) = -factor; // dx2-y2 (dzz)
            }

            // f-functions
            if (_lmax > 2) { // order of functions changed
                const double factor = 2. * pow(_decay, 1.5) * contractions[3];
                const double factor_1 = factor * 2. / sqrt(15.);
                const double factor_2 = factor * sqrt(2.) / sqrt(5.);
                const double factor_3 = factor * sqrt(2.) / sqrt(3.);

                _trafo(Cart::xxz,9) = -3. * factor_1; // f1 (f??) xxz 13
                _trafo(Cart::yyz,9) = -3. * factor_1; // f1 (f??) yyz 15        Y 3,0
                _trafo(Cart::zzz,9) = 2. * factor_1; // f1 (f??) zzz 19

                _trafo(Cart::xxy) = -factor_2; // f3 xxy 10
                _trafo(Cart::yyy,10) = -factor_2; // f3 yyy 18   Y 3,-1
                _trafo(Cart::yzz,10) = 4. * factor_2; // f3 yzz 16

                _trafo(Cart::xxx,11) = -factor_2; // f2 xxx 17
                _trafo(Cart::xyy,11) = -factor_2; // f2 xyy 11   Y 3,1
                _trafo(Cart::xzz,11) = 4. * factor_2; // f2 xzz 14

                _trafo(Cart::xyz,12) = 4. * factor; // f6 xyz 12     Y 3,-2

                _trafo(Cart::xxz,13) = 2. * factor; // f7 (f??)   xxz   13
                _trafo(Cart::yyz,13) = -2. * factor; // f7 (f??)   yyz   15   Y 3,2

                _trafo(Cart::xxy,14) = 3. * factor_3; // f4 xxy 10
                _trafo(Cart::yyy,14) = -factor_3; // f4 yyy 18   Y 3,-3

                _trafo(Cart::xxx,15) = factor_3; // f5 (f??) xxx 17
                _trafo(Cart::xyy,15) = -3. * factor_3; // f5 (f??) xyy 11     Y 3,3
            }

            // g-functions
            if (_lmax > 3) {
                const double factor = 2. / sqrt(3.) * _decay * _decay * contractions[4];
                const double factor_1 = factor / sqrt(35.);
                const double factor_2 = factor * 4. / sqrt(14.);
                const double factor_3 = factor * 2. / sqrt(7.);
                const double factor_4 = factor * 2. * sqrt(2.);

                _trafo(Cart::xxxx,16) = 3. * factor_1; /// Y 4,0
                _trafo(Cart::xxyy,16) = 6. * factor_1;
                _trafo(Cart::xxzz,16) = -24. * factor_1;
                _trafo(Cart::yyyy,16) = 3. * factor_1;
                _trafo(Cart::yyzz,16) = -24. * factor_1;
                _trafo(Cart::zzzz,16) = 8. * factor_1;

                _trafo(Cart::xxyz,17) = -3. * factor_2; /// Y 4,-1
                _trafo(Cart::yyyz,17) = -3. * factor_2;
                _trafo(Cart::yzzz,17) = 4. * factor_2;

                _trafo(Cart::xxxz,18) = -3. * factor_2; /// Y 4,1
                _trafo(Cart::xyyz,18) = -3. * factor_2;
                _trafo(Cart::xzzz,18) = 4. * factor_2;

                _trafo(Cart::xxxy,19) = -2. * factor_3; /// Y 4,-2
                _trafo(Cart::xyyy,19) = -2. * factor_3;
                _trafo(Cart::xyzz,19) = 12. * factor_3;

                _trafo(Cart::xxxx,20) = -factor_3; /// Y 4,2
                _trafo(Cart::xxzz,20) = 6. * factor_3;
                _trafo(Cart::yyyy,20) = factor_3;
                _trafo(Cart::yyzz,20) = -6. * factor_3;

                _trafo(Cart::xxyz,21) = 3. * factor_4; /// Y 4,-3
                _trafo(Cart::yyyz,21) = -factor_4;

                _trafo(Cart::xxxz,22) = factor_4; /// Y 4,3
                _trafo(Cart::xyyz,22) = -3. * factor_4;

                _trafo(Cart::xxxy,23) = 4. * factor; /// Y 4,-4
                _trafo(Cart::xyyy,23) = -4. * factor;

                _trafo(Cart::xxxx,24) = factor; /// Y 4,4
                _trafo(Cart::xxyy,24) = -6. * factor;
                _trafo(Cart::yyyy,24) = factor;
            }
            // h-functions
            if (_lmax > 4) {
                const double factor = (2. / 3.) * pow(_decay, 2.5) * contractions[5];
                const double factor_1 = factor * 2. / sqrt(105.);
                const double factor_2 = factor * 2. / sqrt(7.);
                const double factor_3 = factor * sqrt(6.) / 3.;
                const double factor_4 = factor * 2. * sqrt(3.);
                const double factor_5 = factor * .2 * sqrt(30.);

                _trafo(Cart::xxxxz,25) = 15. * factor_1; /// Y 5,0
                _trafo(Cart::xxyyz,25) = 30. * factor_1;
                _trafo(Cart::xxzzz,25) = -40. * factor_1;
                _trafo(Cart::yyyyz,25) = 15. * factor_1;
                _trafo(Cart::yyzzz,25) = -40. * factor_1;
                _trafo(Cart::zzzzz,25) = 8. * factor_1;

                _trafo(Cart::xxxxy,26) = factor_2; /// Y 5,-1
                _trafo(Cart::xxyyy,26) = 2. * factor_2;
                _trafo(Cart::xxyzz,26) = -12. * factor_2;
                _trafo(Cart::yyyyy,26) = factor_2;
                _trafo(Cart::yyyzz,26) = -12. * factor_2;
                _trafo(Cart::yzzzz,26) = 8. * factor_2;

                _trafo(Cart::xxxxx,27) = factor_2; /// Y 5,1
                _trafo(Cart::xxxyy,27) = 2. * factor_2;
                _trafo(Cart::xxxzz,27) = -12. * factor_2;
                _trafo(Cart::xyyyy,27) = factor_2;
                _trafo(Cart::xyyzz,27) = -12. * factor_2;
                _trafo(Cart::xzzzz,27) = 8. * factor_2;

                _trafo(Cart::xxxyz,28) = -8. * factor; /// Y 5,-2
                _trafo(Cart::xyyyz,28) = -8. * factor;
                _trafo(Cart::xyzzz,28) = 16. * factor;

                _trafo(Cart::xxxxz,29) = -4. * factor; /// Y 5,2
                _trafo(Cart::xxzzz,29) = 8. * factor;
                _trafo(Cart::yyyyz,29) = 4. * factor;
                _trafo(Cart::yyzzz,29) = -8. * factor;

                _trafo(Cart::xxxxy,30) = -3. * factor_3; /// Y 5,-3
                _trafo(Cart::xxyyy,30) = -2. * factor_3;
                _trafo(Cart::xxyzz,30) = 24. * factor_3;
                _trafo(Cart::yyyyy,30) = factor_3;
                _trafo(Cart::yyyzz,30) = -8. * factor_3;

                _trafo(Cart::xxxxx,31) = -factor_3; /// Y 5,3
                _trafo(Cart::xxxyy,31) = 2. * factor_3;
                _trafo(Cart::xxxzz,31) = 8. * factor_3;
                _trafo(Cart::xyyyy,31) = 3. * factor_3;
                _trafo(Cart::xyyzz,31) = -24. * factor_3;

                _trafo(Cart::xxxyz,32) = 4. * factor_4; /// Y 5,-4
                _trafo(Cart::xyyyz,32) = -4. * factor_4;

                _trafo(Cart::xxxxz,33) = factor_4; /// Y 5,4
                _trafo(Cart::xxyyz,33) = -6. * factor_4;
                _trafo(Cart::yyyyz,33) = factor_4;

                _trafo(Cart::xxxxy,34) = 5. * factor_5; /// Y 5,-5
                _trafo(Cart::xxyyy,34) = -10. * factor_5;
                _trafo(Cart::yyyyy,34) = factor_5;

                _trafo(Cart::xxxxx,35) = factor_5; /// Y 5,5
                _trafo(Cart::xxxyy,35) = -10. * factor_5;
                _trafo(Cart::xyyyy,35) = 5. * factor_5;
            }

            // i-functions
            if (_lmax > 5) {
                const double factor = (2. / 3.) * _decay * _decay * _decay * contractions[6];
                const double factor_1 = factor * 2. / sqrt(1155.);
                const double factor_2 = factor * 4. / sqrt(55.);
                const double factor_3 = factor * sqrt(22.) / 11.;
                const double factor_4 = factor * 2. * sqrt(165.) / 55.;
                const double factor_5 = factor * .4 * sqrt(30.);
                const double factor_6 = factor * .2 * sqrt(10.);

                _trafo(Cart::xxxxxx,36) = -5. * factor_1; /// Y 6,0
                _trafo(Cart::xxxxyy,36) = -15. * factor_1;
                _trafo(Cart::xxxxzz,36) = 90. * factor_1;
                _trafo(Cart::xxyyyy,36) = -15. * factor_1;
                _trafo(Cart::xxyyzz,36) = 180. * factor_1;
                _trafo(Cart::xxzzzz,36) = -120. * factor_1;
                _trafo(Cart::yyyyyy,36) = -5. * factor_1;
                _trafo(Cart::yyyyzz,36) = 90. * factor_1;
                _trafo(Cart::yyzzzz,36) = -120. * factor_1;
                _trafo(Cart::zzzzzz,36) = 16. * factor_1;

                _trafo(Cart::xxxxyz,37) = 5. * factor_2; /// Y 6,-1
                _trafo(Cart::xxyyyz,37) = 10. * factor_2;
                _trafo(Cart::xxyzzz,37) = -20. * factor_2;
                _trafo(Cart::yyyyyz,37) = 5. * factor_2;
                _trafo(Cart::yyyzzz,37) = -20. * factor_2;
                _trafo(Cart::yzzzzz,37) = 8. * factor_2;

                _trafo(Cart::xxxxxz,38) = 5. * factor_2; /// Y 6,1
                _trafo(Cart::xxxyyz,38) = 10. * factor_2;
                _trafo(Cart::xxxzzz,38) = -20. * factor_2;
                _trafo(Cart::xyyyyz,38) = 5. * factor_2;
                _trafo(Cart::xyyzzz,38) = -20. * factor_2;
                _trafo(Cart::xzzzzz,38) = 8. * factor_2;

                _trafo(Cart::xxxxxy,39) = 2. * factor_3; /// Y 6,-2
                _trafo(Cart::xxxyyy,39) = 4. * factor_3;
                _trafo(Cart::xxxyzz,39) = -32. * factor_3;
                _trafo(Cart::xyyyyy,39) = 2. * factor_3;
                _trafo(Cart::xyyyzz,39) = -32. * factor_3;
                _trafo(Cart::xyzzzz,39) = 32. * factor_3;

                _trafo(Cart::xxxxxy,40) = factor_3; /// Y 6,2
                _trafo(Cart::xxxxyy,40) = factor_3;
                _trafo(Cart::xxxxzz,40) = -16. * factor_3;
                _trafo(Cart::xxyyyy,40) = -factor_3;
                _trafo(Cart::xxzzzz,40) = 16. * factor_3;
                _trafo(Cart::yyyyyy,40) = -factor_3;
                _trafo(Cart::yyyyzz,40) = 16. * factor_3;
                _trafo(Cart::yyzzzz,40) = -16. * factor_3;

                _trafo(Cart::xxxxyz,41) = -18. * factor_3; /// Y 6,-3
                _trafo(Cart::xxyyyz,41) = -12. * factor_3;
                _trafo(Cart::xxyzzz,41) = 48. * factor_3;
                _trafo(Cart::yyyyyz,41) = 6. * factor_3;
                _trafo(Cart::yyyzzz,41) = -16. * factor_3;

                _trafo(Cart::xxxxxz,42) = -6. * factor_3; /// Y 6,3
                _trafo(Cart::xxxyyz,42) = 12. * factor_3;
                _trafo(Cart::xxxzzz,42) = 16. * factor_3;
                _trafo(Cart::xyyyyz,42) = 18. * factor_3;
                _trafo(Cart::xyyzzz,42) = -48. * factor_3;

                _trafo(Cart::xxxxxy,43) = -4. * factor_4; /// Y 6,-4
                _trafo(Cart::xxxyzz,43) = 40. * factor_4;
                _trafo(Cart::xyyyyy,43) = 4. * factor_4;
                _trafo(Cart::xyyyzz,43) = -40. * factor_4;

                _trafo(Cart::xxxxxx,44) = -factor_4; /// Y 6,4
                _trafo(Cart::xxxxyy,44) = 5. * factor_4;
                _trafo(Cart::xxxxzz,44) = 10. * factor_4;
                _trafo(Cart::xxyyyy,44) = 5. * factor_4;
                _trafo(Cart::xxyyzz,44) = -60. * factor_4;
                _trafo(Cart::yyyyyy,44) = -factor_4;
                _trafo(Cart::yyyyzz,44) = 10. * factor_4;

                _trafo(Cart::xxxxyz,45) = 5. * factor_5; /// Y 6,-5
                _trafo(Cart::xxyyyz,45) = -10. * factor_5;
                _trafo(Cart::yyyyyz,45) = factor_5;

                _trafo(Cart::xxxxxz,46) = factor_5; /// Y 6,5
                _trafo(Cart::xxxyyz,46) = -10. * factor_5;
                _trafo(Cart::xyyyyz,46) = 5. * factor_5;

                _trafo(Cart::xxxxxy,47) = 6. * factor_6; /// Y 6,-6
                _trafo(Cart::xxxyyy,47) = -20. * factor_6;
                _trafo(Cart::xyyyyy,47) = 6. * factor_6;

                _trafo(Cart::xxxxxx,48) = factor_6; /// Y 6,6
                _trafo(Cart::xxxxyy,48) = -15. * factor_6;
                _trafo(Cart::xxyyyy,48) = 15. * factor_6;
                _trafo(Cart::yyyyyy,48) = -factor_6;
            }
            return _trafo;
        }


        template<class T> 
        std::vector<double> AOMatrix<T>::XIntegrate(int _n, double _T) {
            std::vector<double> _FmT = std::vector<double>(_n, 0.0);
            const int _mm = _FmT.size() - 1;
            const double pi = boost::math::constants::pi<double>();
            if (_mm < 0) {
                std::cerr << "mm is: " << _mm << " This should not have happened!" << std::flush;
                exit(1);
            }

            if (_T < 0.0) {
                std::cerr << "T is: " << _T << " This should not have happened!" << std::flush;
                exit(1);
            }

            if (_T >= 10.0) {
                // forward iteration
                _FmT[0] = 0.50 * sqrt(pi / _T) * erf(sqrt(_T));

                for (unsigned m = 1; m < _FmT.size(); m++) {
                    _FmT[m] = (2.0 * m - 1) * _FmT[m - 1] / (2.0 * _T) - exp(-_T) / (2.0 * _T);
                }
            }

            if (_T < 1e-10) {
                for (unsigned m = 0; m < _FmT.size(); m++) {
                    _FmT[m] = 1.0 / (2.0 * m + 1.0) - _T / (2.0 * m + 3.0);
                }
            }


            if (_T >= 1e-10 && _T < 10.0) {
                // backward iteration
                double fm = 0.0;
                for (int m = 60; m >= _mm; m--) {
                    fm = (2.0 * _T) / (2.0 * m + 1.0) * (fm + exp(-_T) / (2.0 * _T));
                }
                _FmT[_mm] = fm;
                for (int m = _mm - 1; m >= 0; m--) {
                    _FmT[m] = (2.0 * _T) / (2.0 * m + 1.0) * (_FmT[m + 1] + exp(-_T) / (2.0 * _T));
                }
            }

            return _FmT;
        }
        
int AOSuperMatrix::getBlockSize(int _lmax) {
      //Each cartesian shells has (l+1)(l+2)/2 elements
      //Sum of all shells up to _lmax leads to blocksize=1+11/6 l+l^2+1/6 l^3
      int blocksize = 6 + 11 * _lmax + 6 * _lmax * _lmax + _lmax * _lmax*_lmax;
      blocksize /= 6;
      return blocksize;
    }



template class AOMatrix<double>;
template class AOMatrix< std::complex<double> >;

    }
}

