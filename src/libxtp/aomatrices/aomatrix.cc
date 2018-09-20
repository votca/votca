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
                int row_start = shell_row->getStartIndex();
                std::string type = shell_row->getType();
                std::cout << "Shell " << type << "starts at " << row_start + 1 << std::endl;
            }
            return;
        }

        
        template< class T> 
        void AOMatrix<T>::Fill(const AOBasis& aobasis) {
            _aomatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>::Zero(aobasis.AOBasisSize(),aobasis.AOBasisSize());
            // loop row
#pragma omp parallel for schedule(guided)
            for (unsigned row = 0; row < aobasis.getNumofShells(); row++) {
                const AOShell* shell_row = aobasis.getShell(row);
                int row_start = shell_row->getStartIndex();

                // AOMatrix is symmetric, restrict explicit calculation to triangular matrix
                for (unsigned col = row; col <  aobasis.getNumofShells(); col++) {
                    const AOShell* shell_col = aobasis.getShell(col);
                    // figure out the submatrix
                    int col_start = shell_col->getStartIndex();
                    Eigen::Block< Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > block=
                            _aomatrix.block(row_start,col_start, shell_row->getNumFunc(),shell_col->getNumFunc());
                    // Fill block
                    FillBlock(block, shell_row, shell_col);
                }
            }
            // Fill whole matrix by copying
        for ( unsigned i=0; i < _aomatrix.rows(); i++){
                for (unsigned j = 0; j < i; j++) {
                    _aomatrix(i, j) = _aomatrix(j, i);
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
            for (unsigned row = 0; row < aobasis.getNumofShells(); row++) {
                const AOShell* shell_row = aobasis.getShell(row);
                int row_start = shell_row->getStartIndex();
                // loop column
                for (const AOShell* shell_col:aobasis) {
                    // figure out the submatrix
                    int col_start = shell_col->getStartIndex();
                std::vector< Eigen::Block<Eigen::MatrixXd> > submatrix;
                    for (int i = 0; i < 3; i++) {
                   Eigen::Block<Eigen::MatrixXd> block=_aomatrix[i].block( row_start,col_start,shell_row->getNumFunc(),shell_col->getNumFunc());
                   submatrix.push_back(block );
                    }
                    // Fill block
                    FillBlock(submatrix, shell_row, shell_col);
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

        template<class T> void AOMatrix<T>::Print(std::string ident) {
            std::cout << "\n" << std::endl;
            std::cout.precision(12);
        for ( unsigned i =0; i< _aomatrix.rows(); i++){
            for ( unsigned j =0; j< _aomatrix.cols(); j++){
                std::cout << ident << "[" << i+1 << ":" << j+1 << "] " << std::scientific <<_aomatrix(i,j) << std::endl;
                }
            }
    }   

        void AOMatrix3D::Print(std::string ident) {
            std::cout << "\n" << std::endl;
        for ( unsigned i =0; i< _aomatrix[0].rows(); i++){
            for ( unsigned j =0; j< _aomatrix[0].cols(); j++){
                std::cout << ident << "[" << i+1 << ":" << j+1 << "] " <<  _aomatrix[0](i,j) << " : " <<  _aomatrix[1](i,j) << " : " <<  _aomatrix[2](i,j)  << std::endl;
                }
            }
        }

       Eigen::MatrixXd AOSuperMatrix::getTrafo(const AOGaussianPrimitive& gaussian){
            ///         0    1  2  3    4  5  6  7  8  9   10  11  12  13  14  15  16  17  18  19       20    21    22    23    24    25    26    27    28    29    30    31    32    33    34 
            ///         s,   x, y, z,   xy xz yz xx yy zz, xxy xyy xyz xxz xzz yyz yzz xxx yyy zzz,    xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyz, yyzz, yzzz, xxxx, yyyy, zzzz,
            const AOShell* shell = gaussian.getShell();
            const int ntrafo = shell->getNumFunc() + shell->getOffset();
            const double decay = gaussian.getDecay();
            const int lmax = shell->getLmax();
            const int n = getBlockSize(lmax);
         Eigen::MatrixXd trafo=Eigen::MatrixXd::Zero(n,ntrafo); 
            const std::vector<double>& contractions = gaussian.getContraction();

            // s-functions
            trafo(0, 0) = contractions[0]; //  // s  Y 0,0
            // p-functions
            if (lmax > 0) { // order of functions changed
                const double factor = 2. * sqrt(decay) * contractions[1];
                trafo(3, 1) = factor; // Y 1,0
                trafo(2, 2) = factor; // Y 1,-1
                trafo(1, 3) = factor; // Y 1,1
            }

            // d-functions
            if (lmax > 1) { // order of functions changed
                const double factor = 2. * decay * contractions[2];
                const double factor_1 = factor / sqrt(3.);
                trafo(Cart::xx,4) = -factor_1; // d3z2-r2 (dxx)
                trafo(Cart::yy,4) = -factor_1; // d3z2-r2 (dyy)  Y 2,0
                trafo(Cart::zz,4) = 2. * factor_1; // d3z2-r2 (dzz)

                trafo(Cart::yz,5) = 2. * factor; // dyz           Y 2,-1

                trafo(Cart::xz,6) = 2. * factor; // dxz           Y 2,1

                trafo(Cart::xy,7) = 2. * factor; // dxy           Y 2,-2

                trafo(Cart::xx,8) = factor; // dx2-y2 (dxx)   Y 2,2
                trafo(Cart::yy,8) = -factor; // dx2-y2 (dzz)
            }

            // f-functions
            if (lmax > 2) { // order of functions changed
                const double factor = 2. * pow(decay, 1.5) * contractions[3];
                const double factor_1 = factor * 2. / sqrt(15.);
                const double factor_2 = factor * sqrt(2.) / sqrt(5.);
                const double factor_3 = factor * sqrt(2.) / sqrt(3.);

                trafo(Cart::xxz,9) = -3. * factor_1; // f1 (f??) xxz 13
                trafo(Cart::yyz,9) = -3. * factor_1; // f1 (f??) yyz 15        Y 3,0
                trafo(Cart::zzz,9) = 2. * factor_1; // f1 (f??) zzz 19

                trafo(Cart::xxy,10) = -factor_2; // f3 xxy 10
                trafo(Cart::yyy,10) = -factor_2; // f3 yyy 18   Y 3,-1
                trafo(Cart::yzz,10) = 4. * factor_2; // f3 yzz 16

                trafo(Cart::xxx,11) = -factor_2; // f2 xxx 17
                trafo(Cart::xyy,11) = -factor_2; // f2 xyy 11   Y 3,1
                trafo(Cart::xzz,11) = 4. * factor_2; // f2 xzz 14

                trafo(Cart::xyz,12) = 4. * factor; // f6 xyz 12     Y 3,-2

                trafo(Cart::xxz,13) = 2. * factor; // f7 (f??)   xxz   13
                trafo(Cart::yyz,13) = -2. * factor; // f7 (f??)   yyz   15   Y 3,2

                trafo(Cart::xxy,14) = 3. * factor_3; // f4 xxy 10
                trafo(Cart::yyy,14) = -factor_3; // f4 yyy 18   Y 3,-3

                trafo(Cart::xxx,15) = factor_3; // f5 (f??) xxx 17
                trafo(Cart::xyy,15) = -3. * factor_3; // f5 (f??) xyy 11     Y 3,3
            }

            // g-functions
            if (lmax > 3) {
                const double factor = 2. / sqrt(3.) * decay * decay * contractions[4];
                const double factor_1 = factor / sqrt(35.);
                const double factor_2 = factor * 4. / sqrt(14.);
                const double factor_3 = factor * 2. / sqrt(7.);
                const double factor_4 = factor * 2. * sqrt(2.);

                trafo(Cart::xxxx,16) = 3. * factor_1; /// Y 4,0
                trafo(Cart::xxyy,16) = 6. * factor_1;
                trafo(Cart::xxzz,16) = -24. * factor_1;
                trafo(Cart::yyyy,16) = 3. * factor_1;
                trafo(Cart::yyzz,16) = -24. * factor_1;
                trafo(Cart::zzzz,16) = 8. * factor_1;

                trafo(Cart::xxyz,17) = -3. * factor_2; /// Y 4,-1
                trafo(Cart::yyyz,17) = -3. * factor_2;
                trafo(Cart::yzzz,17) = 4. * factor_2;

                trafo(Cart::xxxz,18) = -3. * factor_2; /// Y 4,1
                trafo(Cart::xyyz,18) = -3. * factor_2;
                trafo(Cart::xzzz,18) = 4. * factor_2;

                trafo(Cart::xxxy,19) = -2. * factor_3; /// Y 4,-2
                trafo(Cart::xyyy,19) = -2. * factor_3;
                trafo(Cart::xyzz,19) = 12. * factor_3;

                trafo(Cart::xxxx,20) = -factor_3; /// Y 4,2
                trafo(Cart::xxzz,20) = 6. * factor_3;
                trafo(Cart::yyyy,20) = factor_3;
                trafo(Cart::yyzz,20) = -6. * factor_3;

                trafo(Cart::xxyz,21) = 3. * factor_4; /// Y 4,-3
                trafo(Cart::yyyz,21) = -factor_4;

                trafo(Cart::xxxz,22) = factor_4; /// Y 4,3
                trafo(Cart::xyyz,22) = -3. * factor_4;

                trafo(Cart::xxxy,23) = 4. * factor; /// Y 4,-4
                trafo(Cart::xyyy,23) = -4. * factor;

                trafo(Cart::xxxx,24) = factor; /// Y 4,4
                trafo(Cart::xxyy,24) = -6. * factor;
                trafo(Cart::yyyy,24) = factor;
            }
            // h-functions
            if (lmax > 4) {
                const double factor = (2. / 3.) * pow(decay, 2.5) * contractions[5];
                const double factor_1 = factor * 2. / sqrt(105.);
                const double factor_2 = factor * 2. / sqrt(7.);
                const double factor_3 = factor * sqrt(6.) / 3.;
                const double factor_4 = factor * 2. * sqrt(3.);
                const double factor_5 = factor * .2 * sqrt(30.);

                trafo(Cart::xxxxz,25) = 15. * factor_1; /// Y 5,0
                trafo(Cart::xxyyz,25) = 30. * factor_1;
                trafo(Cart::xxzzz,25) = -40. * factor_1;
                trafo(Cart::yyyyz,25) = 15. * factor_1;
                trafo(Cart::yyzzz,25) = -40. * factor_1;
                trafo(Cart::zzzzz,25) = 8. * factor_1;

                trafo(Cart::xxxxy,26) = factor_2; /// Y 5,-1
                trafo(Cart::xxyyy,26) = 2. * factor_2;
                trafo(Cart::xxyzz,26) = -12. * factor_2;
                trafo(Cart::yyyyy,26) = factor_2;
                trafo(Cart::yyyzz,26) = -12. * factor_2;
                trafo(Cart::yzzzz,26) = 8. * factor_2;

                trafo(Cart::xxxxx,27) = factor_2; /// Y 5,1
                trafo(Cart::xxxyy,27) = 2. * factor_2;
                trafo(Cart::xxxzz,27) = -12. * factor_2;
                trafo(Cart::xyyyy,27) = factor_2;
                trafo(Cart::xyyzz,27) = -12. * factor_2;
                trafo(Cart::xzzzz,27) = 8. * factor_2;

                trafo(Cart::xxxyz,28) = -8. * factor; /// Y 5,-2
                trafo(Cart::xyyyz,28) = -8. * factor;
                trafo(Cart::xyzzz,28) = 16. * factor;

                trafo(Cart::xxxxz,29) = -4. * factor; /// Y 5,2
                trafo(Cart::xxzzz,29) = 8. * factor;
                trafo(Cart::yyyyz,29) = 4. * factor;
                trafo(Cart::yyzzz,29) = -8. * factor;

                trafo(Cart::xxxxy,30) = -3. * factor_3; /// Y 5,-3
                trafo(Cart::xxyyy,30) = -2. * factor_3;
                trafo(Cart::xxyzz,30) = 24. * factor_3;
                trafo(Cart::yyyyy,30) = factor_3;
                trafo(Cart::yyyzz,30) = -8. * factor_3;

                trafo(Cart::xxxxx,31) = -factor_3; /// Y 5,3
                trafo(Cart::xxxyy,31) = 2. * factor_3;
                trafo(Cart::xxxzz,31) = 8. * factor_3;
                trafo(Cart::xyyyy,31) = 3. * factor_3;
                trafo(Cart::xyyzz,31) = -24. * factor_3;

                trafo(Cart::xxxyz,32) = 4. * factor_4; /// Y 5,-4
                trafo(Cart::xyyyz,32) = -4. * factor_4;

                trafo(Cart::xxxxz,33) = factor_4; /// Y 5,4
                trafo(Cart::xxyyz,33) = -6. * factor_4;
                trafo(Cart::yyyyz,33) = factor_4;

                trafo(Cart::xxxxy,34) = 5. * factor_5; /// Y 5,-5
                trafo(Cart::xxyyy,34) = -10. * factor_5;
                trafo(Cart::yyyyy,34) = factor_5;

                trafo(Cart::xxxxx,35) = factor_5; /// Y 5,5
                trafo(Cart::xxxyy,35) = -10. * factor_5;
                trafo(Cart::xyyyy,35) = 5. * factor_5;
            }

            // i-functions
            if (lmax > 5) {
                const double factor = (2. / 3.) * decay * decay * decay * contractions[6];
                const double factor_1 = factor * 2. / sqrt(1155.);
                const double factor_2 = factor * 4. / sqrt(55.);
                const double factor_3 = factor * sqrt(22.) / 11.;
                const double factor_4 = factor * 2. * sqrt(165.) / 55.;
                const double factor_5 = factor * .4 * sqrt(30.);
                const double factor_6 = factor * .2 * sqrt(10.);

                trafo(Cart::xxxxxx,36) = -5. * factor_1; /// Y 6,0
                trafo(Cart::xxxxyy,36) = -15. * factor_1;
                trafo(Cart::xxxxzz,36) = 90. * factor_1;
                trafo(Cart::xxyyyy,36) = -15. * factor_1;
                trafo(Cart::xxyyzz,36) = 180. * factor_1;
                trafo(Cart::xxzzzz,36) = -120. * factor_1;
                trafo(Cart::yyyyyy,36) = -5. * factor_1;
                trafo(Cart::yyyyzz,36) = 90. * factor_1;
                trafo(Cart::yyzzzz,36) = -120. * factor_1;
                trafo(Cart::zzzzzz,36) = 16. * factor_1;

                trafo(Cart::xxxxyz,37) = 5. * factor_2; /// Y 6,-1
                trafo(Cart::xxyyyz,37) = 10. * factor_2;
                trafo(Cart::xxyzzz,37) = -20. * factor_2;
                trafo(Cart::yyyyyz,37) = 5. * factor_2;
                trafo(Cart::yyyzzz,37) = -20. * factor_2;
                trafo(Cart::yzzzzz,37) = 8. * factor_2;

                trafo(Cart::xxxxxz,38) = 5. * factor_2; /// Y 6,1
                trafo(Cart::xxxyyz,38) = 10. * factor_2;
                trafo(Cart::xxxzzz,38) = -20. * factor_2;
                trafo(Cart::xyyyyz,38) = 5. * factor_2;
                trafo(Cart::xyyzzz,38) = -20. * factor_2;
                trafo(Cart::xzzzzz,38) = 8. * factor_2;

                trafo(Cart::xxxxxy,39) = 2. * factor_3; /// Y 6,-2
                trafo(Cart::xxxyyy,39) = 4. * factor_3;
                trafo(Cart::xxxyzz,39) = -32. * factor_3;
                trafo(Cart::xyyyyy,39) = 2. * factor_3;
                trafo(Cart::xyyyzz,39) = -32. * factor_3;
                trafo(Cart::xyzzzz,39) = 32. * factor_3;

                trafo(Cart::xxxxxy,40) = factor_3; /// Y 6,2
                trafo(Cart::xxxxyy,40) = factor_3;
                trafo(Cart::xxxxzz,40) = -16. * factor_3;
                trafo(Cart::xxyyyy,40) = -factor_3;
                trafo(Cart::xxzzzz,40) = 16. * factor_3;
                trafo(Cart::yyyyyy,40) = -factor_3;
                trafo(Cart::yyyyzz,40) = 16. * factor_3;
                trafo(Cart::yyzzzz,40) = -16. * factor_3;

                trafo(Cart::xxxxyz,41) = -18. * factor_3; /// Y 6,-3
                trafo(Cart::xxyyyz,41) = -12. * factor_3;
                trafo(Cart::xxyzzz,41) = 48. * factor_3;
                trafo(Cart::yyyyyz,41) = 6. * factor_3;
                trafo(Cart::yyyzzz,41) = -16. * factor_3;

                trafo(Cart::xxxxxz,42) = -6. * factor_3; /// Y 6,3
                trafo(Cart::xxxyyz,42) = 12. * factor_3;
                trafo(Cart::xxxzzz,42) = 16. * factor_3;
                trafo(Cart::xyyyyz,42) = 18. * factor_3;
                trafo(Cart::xyyzzz,42) = -48. * factor_3;

                trafo(Cart::xxxxxy,43) = -4. * factor_4; /// Y 6,-4
                trafo(Cart::xxxyzz,43) = 40. * factor_4;
                trafo(Cart::xyyyyy,43) = 4. * factor_4;
                trafo(Cart::xyyyzz,43) = -40. * factor_4;

                trafo(Cart::xxxxxx,44) = -factor_4; /// Y 6,4
                trafo(Cart::xxxxyy,44) = 5. * factor_4;
                trafo(Cart::xxxxzz,44) = 10. * factor_4;
                trafo(Cart::xxyyyy,44) = 5. * factor_4;
                trafo(Cart::xxyyzz,44) = -60. * factor_4;
                trafo(Cart::yyyyyy,44) = -factor_4;
                trafo(Cart::yyyyzz,44) = 10. * factor_4;

                trafo(Cart::xxxxyz,45) = 5. * factor_5; /// Y 6,-5
                trafo(Cart::xxyyyz,45) = -10. * factor_5;
                trafo(Cart::yyyyyz,45) = factor_5;

                trafo(Cart::xxxxxz,46) = factor_5; /// Y 6,5
                trafo(Cart::xxxyyz,46) = -10. * factor_5;
                trafo(Cart::xyyyyz,46) = 5. * factor_5;

                trafo(Cart::xxxxxy,47) = 6. * factor_6; /// Y 6,-6
                trafo(Cart::xxxyyy,47) = -20. * factor_6;
                trafo(Cart::xyyyyy,47) = 6. * factor_6;

                trafo(Cart::xxxxxx,48) = factor_6; /// Y 6,6
                trafo(Cart::xxxxyy,48) = -15. * factor_6;
                trafo(Cart::xxyyyy,48) = 15. * factor_6;
                trafo(Cart::yyyyyy,48) = -factor_6;
            }
            return trafo;
        }


        template<class T> 
        std::vector<double> AOMatrix<T>::XIntegrate(int size, double U) {
            std::vector<double> FmU = std::vector<double>(size, 0.0);
            const int mm = FmU.size() - 1;
            const double pi = boost::math::constants::pi<double>();
            if (mm < 0) {
                std::cerr << "mm is: " << mm << " This should not have happened!" << std::flush;
                exit(1);
            }

            if (U < 0.0) {
                std::cerr << "U is: " << U << " This should not have happened!" << std::flush;
                exit(1);
            }

            if (U >= 10.0) {
                // forward iteration
                FmU[0] = 0.50 * sqrt(pi / U) * erf(sqrt(U));

                for (unsigned m = 1; m < FmU.size(); m++) {
                    FmU[m] = (2.0 * m - 1) * FmU[m - 1] / (2.0 * U) - exp(-U) / (2.0 * U);
                }
            }

            if (U < 1e-10) {
                for (unsigned m = 0; m < FmU.size(); m++) {
                    FmU[m] = 1.0 / (2.0 * m + 1.0) - U / (2.0 * m + 3.0);
                }
            }


            if (U >= 1e-10 && U < 10.0) {
                // backward iteration
                double fm = 0.0;
                for (int m = 60; m >= mm; m--) {
                    fm = (2.0 * U) / (2.0 * m + 1.0) * (fm + exp(-U) / (2.0 * U));
                }
                FmU[mm] = fm;
                for (int m = mm - 1; m >= 0; m--) {
                    FmU[m] = (2.0 * U) / (2.0 * m + 1.0) * (FmU[m + 1] + exp(-U) / (2.0 * U));
                }
            }

            return FmU;
        }
        
int AOSuperMatrix::getBlockSize(int lmax) {
      //Each cartesian shells has (l+1)(l+2)/2 elements
      //Sum of all shells up to _lmax leads to blocksize=1+11/6 l+l^2+1/6 l^3
      int blocksize = 6 + 11 * lmax + 6 * lmax * lmax + lmax * lmax*lmax;
      blocksize /= 6;
      return blocksize;
    }



template class AOMatrix<double>;
template class AOMatrix< std::complex<double> >;

    }
}

