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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aomatrix.h>


namespace votca {
    namespace xtp {
       

        void AOPlanewave::FillBlock(Eigen::Block<Eigen::MatrixXcd>& _matrix, const AOShell* _shell_row, const AOShell* _shell_col) {

            // shell info, only lmax tells how far to go
            int _lmax_row = _shell_row->getLmax();
            int _lmax_col = _shell_col->getLmax();
            // set size of internal block for recursion
            int _nrows = this->getBlockSize(_lmax_row);
            int _ncols = this->getBlockSize(_lmax_col);
            if (_lmax_col > 6 || _lmax_row > 6) {
                std::cerr << "Orbitals higher than i are not yet implemented. This should not have happened!" << std::flush;
                exit(1);
            }
            // get shell positions
            const tools::vec& _pos_row = _shell_row->getPos(); //get position R_{i}
            const tools::vec& _pos_col = _shell_col->getPos(); //get position R_{j}
            const tools::vec _diff = _pos_row - _pos_col; //get difference r_{ij}
            double _distsq = (_diff * _diff); //get |R_{ij}|^2
            // get kvector modulus
            double _kmodulus = (_k * _k); //get |k|^2

            // cout << "k is " << _k << endl;

            int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};

            int nx[] = {0,
                1, 0, 0,
                2, 1, 1, 0, 0, 0,
                3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
                4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0,
                5, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0};

            int ny[] = {0,
                0, 1, 0,
                0, 1, 0, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
                0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0};

            int nz[] = {0,
                0, 0, 1,
                0, 0, 1, 0, 1, 2,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
                0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6};


            int i_less_x[] = {0,
                0, 0, 0,
                1, 2, 3, 0, 0, 0,
                4, 5, 6, 7, 8, 9, 0, 0, 0, 0,
                10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 0, 0, 0, 0, 0,
                20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 0, 0, 0, 0, 0, 0,
                35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 0, 0, 0, 0, 0, 0, 0};

            int i_less_y[] = {0,
                0, 0, 0,
                0, 1, 0, 2, 3, 0,
                0, 4, 0, 5, 6, 0, 7, 8, 9, 0,
                0, 10, 0, 11, 12, 0, 13, 14, 15, 0, 16, 17, 18, 19, 0,
                0, 20, 0, 21, 22, 0, 23, 24, 25, 0, 26, 27, 28, 29, 0, 30, 31, 32, 33, 34, 0,
                0, 35, 0, 36, 37, 0, 38, 39, 40, 0, 41, 42, 43, 44, 0, 45, 46, 47, 48, 49, 0, 50, 51, 52, 53, 54, 55, 0};

            int i_less_z[] = {0,
                0, 0, 0,
                0, 0, 1, 0, 2, 3,
                0, 0, 4, 0, 5, 6, 0, 7, 8, 9,
                0, 0, 10, 0, 11, 12, 0, 13, 14, 15, 0, 16, 17, 18, 19,
                0, 0, 20, 0, 21, 22, 0, 23, 24, 25, 0, 26, 27, 28, 29, 0, 30, 31, 32, 33, 34,
                0, 0, 35, 0, 36, 37, 0, 38, 39, 40, 0, 41, 42, 43, 44, 0, 45, 46, 47, 48, 49, 0, 50, 51, 52, 53, 54, 55};
            // iterate over Gaussians in this _shell_row   
            for (AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr) {
                // iterate over Gaussians in this _shell_col
                // get decay constant
                const double _decay_row = itr->getDecay();

                for (AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc) {
                    //get decay constant
                    const double _decay_col = itc->getDecay();

                    // some helpers

                    const double _fak = 0.5 / (_decay_row + _decay_col);
                    const double _fak2 = 2.0 * _fak;
                    double _exparg = _fak2 * _decay_row * _decay_col *_distsq;

                    // check if distance between postions is big, then skip step   

                    if (_exparg > 30.0) {
                        continue;
                    }

                    // initialize local matrix block for unnormalized cartesians
                    Eigen::MatrixXcd _olk = Eigen::MatrixXcd::Zero(_nrows, _ncols);

                    typedef std::complex<double> COMPLEX; // Define an abbreviation for complex numbers 

                    const COMPLEX PmA0(_fak2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_row.getX(), _fak * _k.getX());
                    const COMPLEX PmA1(_fak2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_row.getY(), _fak * _k.getY());
                    const COMPLEX PmA2(_fak2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_row.getZ(), _fak * _k.getZ());

                    const COMPLEX PmB0(_fak2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_col.getX(), _fak * _k.getX());
                    const COMPLEX PmB1(_fak2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_col.getY(), _fak * _k.getY());
                    const COMPLEX PmB2(_fak2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_col.getZ(), _fak * _k.getZ());

                    const COMPLEX _cfak(_fak, 0.0);
                    const COMPLEX _cfak2(_fak2, 0.0);


                    // calculate s-s- overlap matrix element
                    COMPLEX _ssol(pow(4.0 * _decay_row*_decay_col, 0.75) * pow(_fak2, 1.5) * exp(-_exparg), 0.0); // s-s element

                    // calculate s-W-s matrix element
                    double _kdotr_row = (_k * _pos_row);
                    double _kdotr_col = (_k * _pos_col);
                    COMPLEX _kexparg(_fak2 * (-0.25)*(_kmodulus), _fak2 * (_decay_row)*(_kdotr_row) + _fak2 * (_decay_col)*(_kdotr_col));

                    _olk(0, 0) = _ssol * (std::exp(_kexparg)); // s-W-s element


                    // Integral p-W-s
                    if (_lmax_row > 0) {
                        _olk(Cart::x, 0) = PmA0 * _olk(0, 0);
                        _olk(Cart::y, 0) = PmA1 * _olk(0, 0);
                        _olk(Cart::z, 0) = PmA2 * _olk(0, 0);
                    }
                    //------------------------------------------------------

                    //Integrals     d - W - s
                    if (_lmax_row > 1) {
                        COMPLEX term = (_cfak)*(_olk(0, 0));
                        _olk(Cart::xx, 0) = PmA0 * _olk(Cart::x, 0) + term;
                        _olk(Cart::xy, 0) = PmA0 * _olk(Cart::y, 0);
                        _olk(Cart::xz, 0) = PmA0 * _olk(Cart::z, 0);
                        _olk(Cart::yy, 0) = PmA1 * _olk(Cart::y, 0) + term;
                        _olk(Cart::yz, 0) = PmA1 * _olk(Cart::z, 0);
                        _olk(Cart::zz, 0) = PmA2 * _olk(Cart::z, 0) + term;
                    }
                    //------------------------------------------------------
                    //Integrals     f - W - s
                    if (_lmax_row > 2) {
                        _olk(Cart::xxx, 0) = PmA0 * _olk(Cart::xx, 0) + _cfak2 * _olk(Cart::x, 0);
                        _olk(Cart::xxy, 0) = PmA1 * _olk(Cart::xx, 0);
                        _olk(Cart::xxz, 0) = PmA2 * _olk(Cart::xx, 0);
                        _olk(Cart::xyy, 0) = PmA0 * _olk(Cart::yy, 0);
                        _olk(Cart::xyz, 0) = PmA0 * _olk(Cart::yz, 0);
                        _olk(Cart::xzz, 0) = PmA0 * _olk(Cart::zz, 0);
                        _olk(Cart::yyy, 0) = PmA1 * _olk(Cart::yy, 0) + _cfak2 * _olk(Cart::y, 0);
                        _olk(Cart::yyz, 0) = PmA2 * _olk(Cart::yy, 0);
                        _olk(Cart::yzz, 0) = PmA1 * _olk(Cart::zz, 0);
                        _olk(Cart::zzz, 0) = PmA2 * _olk(Cart::zz, 0) + _cfak2 * _olk(Cart::z, 0);
                    }
                    //------------------------------------------------------
                    //Integrals     g - W - s
                    if (_lmax_row > 3) {
                        COMPLEX term_xx = (_cfak)*(_olk(Cart::xx, 0));
                        COMPLEX term_yy = (_cfak)*(_olk(Cart::yy, 0));
                        COMPLEX term_zz = (_cfak)*(_olk(Cart::zz, 0));
                        _olk(Cart::xxxx, 0) = PmA0 * _olk(Cart::xxx, 0) + 3.0 * term_xx;
                        _olk(Cart::xxxy, 0) = PmA1 * _olk(Cart::xxx, 0);
                        _olk(Cart::xxxz, 0) = PmA2 * _olk(Cart::xxx, 0);
                        _olk(Cart::xxyy, 0) = PmA0 * _olk(Cart::xyy, 0) + term_yy;
                        _olk(Cart::xxyz, 0) = PmA1 * _olk(Cart::xxz, 0);
                        _olk(Cart::xxzz, 0) = PmA0 * _olk(Cart::xzz, 0) + term_zz;
                        _olk(Cart::xyyy, 0) = PmA0 * _olk(Cart::yyy, 0);
                        _olk(Cart::xyyz, 0) = PmA0 * _olk(Cart::yyz, 0);
                        _olk(Cart::xyzz, 0) = PmA0 * _olk(Cart::yzz, 0);
                        _olk(Cart::xzzz, 0) = PmA0 * _olk(Cart::zzz, 0);
                        _olk(Cart::yyyy, 0) = PmA1 * _olk(Cart::yyy, 0) + 3.0 * term_yy;
                        _olk(Cart::yyyz, 0) = PmA2 * _olk(Cart::yyy, 0);
                        _olk(Cart::yyzz, 0) = PmA1 * _olk(Cart::yzz, 0) + term_zz;
                        _olk(Cart::yzzz, 0) = PmA1 * _olk(Cart::zzz, 0);
                        _olk(Cart::zzzz, 0) = PmA2 * _olk(Cart::zzz, 0) + 3.0 * term_zz;
                    }
                    //------------------------------------------------------        
                    //Integrals     h - W - s
                    if (_lmax_row > 4) {
                        COMPLEX term_xxx = (_cfak)*(_olk(Cart::xxx, 0));
                        COMPLEX term_yyy = (_cfak)*(_olk(Cart::yyy, 0));
                        COMPLEX term_zzz = (_cfak)*(_olk(Cart::zzz, 0));
                        _olk(Cart::xxxxx, 0) = PmA0 * _olk(Cart::xxxx, 0) + 4.0 * term_xxx;
                        _olk(Cart::xxxxy, 0) = PmA1 * _olk(Cart::xxxx, 0);
                        _olk(Cart::xxxxz, 0) = PmA2 * _olk(Cart::xxxx, 0);
                        _olk(Cart::xxxyy, 0) = PmA1 * _olk(Cart::xxxy, 0) + term_xxx;
                        _olk(Cart::xxxyz, 0) = PmA1 * _olk(Cart::xxxz, 0);
                        _olk(Cart::xxxzz, 0) = PmA2 * _olk(Cart::xxxz, 0) + term_xxx;
                        _olk(Cart::xxyyy, 0) = PmA0 * _olk(Cart::xyyy, 0) + term_yyy;
                        _olk(Cart::xxyyz, 0) = PmA2 * _olk(Cart::xxyy, 0);
                        _olk(Cart::xxyzz, 0) = PmA1 * _olk(Cart::xxzz, 0);
                        _olk(Cart::xxzzz, 0) = PmA0 * _olk(Cart::xzzz, 0) + term_zzz;
                        _olk(Cart::xyyyy, 0) = PmA0 * _olk(Cart::yyyy, 0);
                        _olk(Cart::xyyyz, 0) = PmA0 * _olk(Cart::yyyz, 0);
                        _olk(Cart::xyyzz, 0) = PmA0 * _olk(Cart::yyzz, 0);
                        _olk(Cart::xyzzz, 0) = PmA0 * _olk(Cart::yzzz, 0);
                        _olk(Cart::xzzzz, 0) = PmA0 * _olk(Cart::zzzz, 0);
                        _olk(Cart::yyyyy, 0) = PmA1 * _olk(Cart::yyyy, 0) + 4.0 * term_yyy;
                        _olk(Cart::yyyyz, 0) = PmA2 * _olk(Cart::yyyy, 0);
                        _olk(Cart::yyyzz, 0) = PmA2 * _olk(Cart::yyyz, 0) + term_yyy;
                        _olk(Cart::yyzzz, 0) = PmA1 * _olk(Cart::yzzz, 0) + term_zzz;
                        _olk(Cart::yzzzz, 0) = PmA1 * _olk(Cart::zzzz, 0);
                        _olk(Cart::zzzzz, 0) = PmA2 * _olk(Cart::zzzz, 0) + 4.0 * term_zzz;
                    }
                    //------------------------------------------------------
                    //Integrals     i -W - s
                    if (_lmax_row > 5) {
                        COMPLEX term_xxxx = (_cfak)*(_olk(Cart::xxxx, 0));
                        COMPLEX term_xyyy = (_cfak)*(_olk(Cart::xyyy, 0));
                        COMPLEX term_xzzz = (_cfak)*(_olk(Cart::xzzz, 0));
                        COMPLEX term_yyyy = (_cfak)*(_olk(Cart::yyyy, 0));
                        COMPLEX term_yyzz = (_cfak)*(_olk(Cart::yyzz, 0));
                        COMPLEX term_yzzz = (_cfak)*(_olk(Cart::yzzz, 0));
                        COMPLEX term_zzzz = (_cfak)*(_olk(Cart::zzzz, 0));
                        _olk(Cart::xxxxxx, 0) = PmA0 * _olk(Cart::xxxxx, 0) + 5.0 * term_xxxx;
                        _olk(Cart::xxxxxy, 0) = PmA1 * _olk(Cart::xxxxx, 0);
                        _olk(Cart::xxxxxz, 0) = PmA2 * _olk(Cart::xxxxx, 0);
                        _olk(Cart::xxxxyy, 0) = PmA1 * _olk(Cart::xxxxy, 0) + term_xxxx;
                        _olk(Cart::xxxxyz, 0) = PmA1 * _olk(Cart::xxxxz, 0);
                        _olk(Cart::xxxxzz, 0) = PmA2 * _olk(Cart::xxxxz, 0) + term_xxxx;
                        _olk(Cart::xxxyyy, 0) = PmA0 * _olk(Cart::xxyyy, 0) + 2.0 * term_xyyy;
                        _olk(Cart::xxxyyz, 0) = PmA2 * _olk(Cart::xxxyy, 0);
                        _olk(Cart::xxxyzz, 0) = PmA1 * _olk(Cart::xxxzz, 0);
                        _olk(Cart::xxxzzz, 0) = PmA0 * _olk(Cart::xxzzz, 0) + 2.0 * term_xzzz;
                        _olk(Cart::xxyyyy, 0) = PmA0 * _olk(Cart::xyyyy, 0) + term_yyyy;
                        _olk(Cart::xxyyyz, 0) = PmA2 * _olk(Cart::xxyyy, 0);
                        _olk(Cart::xxyyzz, 0) = PmA0 * _olk(Cart::xyyzz, 0) + term_yyzz;
                        _olk(Cart::xxyzzz, 0) = PmA1 * _olk(Cart::xxzzz, 0);
                        _olk(Cart::xxzzzz, 0) = PmA0 * _olk(Cart::xzzzz, 0) + term_zzzz;
                        _olk(Cart::xyyyyy, 0) = PmA0 * _olk(Cart::yyyyy, 0);
                        _olk(Cart::xyyyyz, 0) = PmA0 * _olk(Cart::yyyyz, 0);
                        _olk(Cart::xyyyzz, 0) = PmA0 * _olk(Cart::yyyzz, 0);
                        _olk(Cart::xyyzzz, 0) = PmA0 * _olk(Cart::yyzzz, 0);
                        _olk(Cart::xyzzzz, 0) = PmA0 * _olk(Cart::yzzzz, 0);
                        _olk(Cart::xzzzzz, 0) = PmA0 * _olk(Cart::zzzzz, 0);
                        _olk(Cart::yyyyyy, 0) = PmA1 * _olk(Cart::yyyyy, 0) + 5.0 * term_yyyy;
                        _olk(Cart::yyyyyz, 0) = PmA2 * _olk(Cart::yyyyy, 0);
                        _olk(Cart::yyyyzz, 0) = PmA2 * _olk(Cart::yyyyz, 0) + term_yyyy;
                        _olk(Cart::yyyzzz, 0) = PmA1 * _olk(Cart::yyzzz, 0) + 2.0 * term_yzzz;
                        _olk(Cart::yyzzzz, 0) = PmA1 * _olk(Cart::yzzzz, 0) + term_zzzz;
                        _olk(Cart::yzzzzz, 0) = PmA1 * _olk(Cart::zzzzz, 0);
                        _olk(Cart::zzzzzz, 0) = PmA2 * _olk(Cart::zzzzz, 0) + 5.0 * term_zzzz;
                    }
                    //------------------------------------------------------
                    if (_lmax_col > 0) {

                        //Integrals     s - W - p
                        _olk(0, Cart::x) = PmB0 * _olk(0, 0);
                        _olk(0, Cart::y) = PmB1 * _olk(0, 0);
                        _olk(0, Cart::z) = PmB2 * _olk(0, 0);
                        //------------------------------------------------------

                        //Integrals     p - W - p     d - W - p     f - W - p     g - W - p     h - W - p     i - W - p
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            //COMPLEX cnx(nx[_i] * _fak, 0.0);
                            _olk(_i, Cart::x) = PmB0 * _olk(_i, 0) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], 0);
                            //COMPLEX cny = (ny[_i] * _fak, 0.0);
                            _olk(_i, Cart::y) = PmB1 * _olk(_i, 0) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], 0);
                            //COMPLEX cnz = (nz[_i] * _fak, 0.0);
                            _olk(_i, Cart::z) = PmB2 * _olk(_i, 0) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], 0);
                        }
                        //------------------------------------------------------


                    } // end if (_lmax_col > 0)
                    if (_lmax_col > 1) {

                        //Integrals     s - W - d
                        COMPLEX term = _cfak * _olk(0, 0);
                        _olk(0, Cart::xx) = PmB0 * _olk(0, Cart::x) + term;
                        _olk(0, Cart::xy) = PmB0 * _olk(0, Cart::y);
                        _olk(0, Cart::xz) = PmB0 * _olk(0, Cart::z);
                        _olk(0, Cart::yy) = PmB1 * _olk(0, Cart::y) + term;
                        _olk(0, Cart::yz) = PmB1 * _olk(0, Cart::z);
                        _olk(0, Cart::zz) = PmB2 * _olk(0, Cart::z) + term;
                        //------------------------------------------------------

                        //Integrals     p - W - d     d - W - d     f - W - d     g - W - d     h - W - d     i - W - d
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            COMPLEX term = _cfak * _olk(_i, 0);
                            _olk(_i, Cart::xx) = PmB0 * _olk(_i, Cart::x) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::x) + term;
                            _olk(_i, Cart::xy) = PmB0 * _olk(_i, Cart::y) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::y);
                            _olk(_i, Cart::xz) = PmB0 * _olk(_i, Cart::z) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::z);
                            _olk(_i, Cart::yy) = PmB1 * _olk(_i, Cart::y) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::y) + term;
                            _olk(_i, Cart::yz) = PmB1 * _olk(_i, Cart::z) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::z);
                            _olk(_i, Cart::zz) = PmB2 * _olk(_i, Cart::z) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::z) + term;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 1)

                    if (_lmax_col > 2) {

                        //Integrals     s - W - f
                        _olk(0, Cart::xxx) = PmB0 * _olk(0, Cart::xx) + 2.0 * _cfak * _olk(0, Cart::x);
                        _olk(0, Cart::xxy) = PmB1 * _olk(0, Cart::xx);
                        _olk(0, Cart::xxz) = PmB2 * _olk(0, Cart::xx);
                        _olk(0, Cart::xyy) = PmB0 * _olk(0, Cart::yy);
                        _olk(0, Cart::xyz) = PmB0 * _olk(0, Cart::yz);
                        _olk(0, Cart::xzz) = PmB0 * _olk(0, Cart::zz);
                        _olk(0, Cart::yyy) = PmB1 * _olk(0, Cart::yy) + 2.0 * _cfak * _olk(0, Cart::y);
                        _olk(0, Cart::yyz) = PmB2 * _olk(0, Cart::yy);
                        _olk(0, Cart::yzz) = PmB1 * _olk(0, Cart::zz);
                        _olk(0, Cart::zzz) = PmB2 * _olk(0, Cart::zz) + 2.0 * _cfak * _olk(0, Cart::z);
                        //------------------------------------------------------

                        //Integrals     p - f     d - f     f - f     g - f     h - f     i - f
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            COMPLEX term_x = 2.0 * _cfak * _olk(_i, Cart::x);
                            COMPLEX term_y = 2.0 * _cfak * _olk(_i, Cart::y);
                            COMPLEX term_z = 2.0 * _cfak * _olk(_i, Cart::z);
                            _olk(_i, Cart::xxx) = PmB0 * _olk(_i, Cart::xx) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xx) + term_x;
                            _olk(_i, Cart::xxy) = PmB1 * _olk(_i, Cart::xx) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xx);
                            _olk(_i, Cart::xxz) = PmB2 * _olk(_i, Cart::xx) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xx);
                            _olk(_i, Cart::xyy) = PmB0 * _olk(_i, Cart::yy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yy);
                            _olk(_i, Cart::xyz) = PmB0 * _olk(_i, Cart::yz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yz);
                            _olk(_i, Cart::xzz) = PmB0 * _olk(_i, Cart::zz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::zz);
                            _olk(_i, Cart::yyy) = PmB1 * _olk(_i, Cart::yy) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yy) + term_y;
                            _olk(_i, Cart::yyz) = PmB2 * _olk(_i, Cart::yy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::yy);
                            _olk(_i, Cart::yzz) = PmB1 * _olk(_i, Cart::zz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::zz);
                            _olk(_i, Cart::zzz) = PmB2 * _olk(_i, Cart::zz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::zz) + term_z;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 2)         

                    if (_lmax_col > 3) {

                        //Integrals     s - W - g
                        COMPLEX term_xx = _cfak * _olk(0, Cart::xx);
                        COMPLEX term_yy = _cfak * _olk(0, Cart::yy);
                        COMPLEX term_zz = _cfak * _olk(0, Cart::zz);
                        _olk(0, Cart::xxxx) = PmB0 * _olk(0, Cart::xxx) + 3.0 * term_xx;
                        _olk(0, Cart::xxxy) = PmB1 * _olk(0, Cart::xxx);
                        _olk(0, Cart::xxxz) = PmB2 * _olk(0, Cart::xxx);
                        _olk(0, Cart::xxyy) = PmB0 * _olk(0, Cart::xyy) + term_yy;
                        _olk(0, Cart::xxyz) = PmB1 * _olk(0, Cart::xxz);
                        _olk(0, Cart::xxzz) = PmB0 * _olk(0, Cart::xzz) + term_zz;
                        _olk(0, Cart::xyyy) = PmB0 * _olk(0, Cart::yyy);
                        _olk(0, Cart::xyyz) = PmB0 * _olk(0, Cart::yyz);
                        _olk(0, Cart::xyzz) = PmB0 * _olk(0, Cart::yzz);
                        _olk(0, Cart::xzzz) = PmB0 * _olk(0, Cart::zzz);
                        _olk(0, Cart::yyyy) = PmB1 * _olk(0, Cart::yyy) + 3.0 * term_yy;
                        _olk(0, Cart::yyyz) = PmB2 * _olk(0, Cart::yyy);
                        _olk(0, Cart::yyzz) = PmB1 * _olk(0, Cart::yzz) + term_zz;
                        _olk(0, Cart::yzzz) = PmB1 * _olk(0, Cart::zzz);
                        _olk(0, Cart::zzzz) = PmB2 * _olk(0, Cart::zzz) + 3.0 * term_zz;
                        //------------------------------------------------------

                        //Integrals     p - g     d - g     f - g     g - g     h - g     i - g
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            COMPLEX term_xx = _cfak * _olk(_i, Cart::xx);
                            COMPLEX term_yy = _cfak * _olk(_i, Cart::yy);
                            COMPLEX term_zz = _cfak * _olk(_i, Cart::zz);
                            _olk(_i, Cart::xxxx) = PmB0 * _olk(_i, Cart::xxx) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xxx) + 3.0 * term_xx;
                            _olk(_i, Cart::xxxy) = PmB1 * _olk(_i, Cart::xxx) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxx);
                            _olk(_i, Cart::xxxz) = PmB2 * _olk(_i, Cart::xxx) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxx);
                            _olk(_i, Cart::xxyy) = PmB0 * _olk(_i, Cart::xyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xyy) + term_yy;
                            _olk(_i, Cart::xxyz) = PmB1 * _olk(_i, Cart::xxz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxz);
                            _olk(_i, Cart::xxzz) = PmB0 * _olk(_i, Cart::xzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xzz) + term_zz;
                            _olk(_i, Cart::xyyy) = PmB0 * _olk(_i, Cart::yyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyy);
                            _olk(_i, Cart::xyyz) = PmB0 * _olk(_i, Cart::yyz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyz);
                            _olk(_i, Cart::xyzz) = PmB0 * _olk(_i, Cart::yzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yzz);
                            _olk(_i, Cart::xzzz) = PmB0 * _olk(_i, Cart::zzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::zzz);
                            _olk(_i, Cart::yyyy) = PmB1 * _olk(_i, Cart::yyy) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yyy) + 3.0 * term_yy;
                            _olk(_i, Cart::yyyz) = PmB2 * _olk(_i, Cart::yyy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::yyy);
                            _olk(_i, Cart::yyzz) = PmB1 * _olk(_i, Cart::yzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yzz) + term_zz;
                            _olk(_i, Cart::yzzz) = PmB1 * _olk(_i, Cart::zzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::zzz);
                            _olk(_i, Cart::zzzz) = PmB2 * _olk(_i, Cart::zzz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::zzz) + 3.0 * term_zz;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 3)        

                    if (_lmax_col > 4) {

                        //Integrals     s - h
                        COMPLEX term_xxx = _cfak * _olk(0, Cart::xxx);
                        COMPLEX term_yyy = _cfak * _olk(0, Cart::yyy);
                        COMPLEX term_zzz = _cfak * _olk(0, Cart::zzz);
                        _olk(0, Cart::xxxxx) = PmB0 * _olk(0, Cart::xxxx) + 4.0 * term_xxx;
                        _olk(0, Cart::xxxxy) = PmB1 * _olk(0, Cart::xxxx);
                        _olk(0, Cart::xxxxz) = PmB2 * _olk(0, Cart::xxxx);
                        _olk(0, Cart::xxxyy) = PmB1 * _olk(0, Cart::xxxy) + term_xxx;
                        _olk(0, Cart::xxxyz) = PmB1 * _olk(0, Cart::xxxz);
                        _olk(0, Cart::xxxzz) = PmB2 * _olk(0, Cart::xxxz) + term_xxx;
                        _olk(0, Cart::xxyyy) = PmB0 * _olk(0, Cart::xyyy) + term_yyy;
                        _olk(0, Cart::xxyyz) = PmB2 * _olk(0, Cart::xxyy);
                        _olk(0, Cart::xxyzz) = PmB1 * _olk(0, Cart::xxzz);
                        _olk(0, Cart::xxzzz) = PmB0 * _olk(0, Cart::xzzz) + term_zzz;
                        _olk(0, Cart::xyyyy) = PmB0 * _olk(0, Cart::yyyy);
                        _olk(0, Cart::xyyyz) = PmB0 * _olk(0, Cart::yyyz);
                        _olk(0, Cart::xyyzz) = PmB0 * _olk(0, Cart::yyzz);
                        _olk(0, Cart::xyzzz) = PmB0 * _olk(0, Cart::yzzz);
                        _olk(0, Cart::xzzzz) = PmB0 * _olk(0, Cart::zzzz);
                        _olk(0, Cart::yyyyy) = PmB1 * _olk(0, Cart::yyyy) + 4.0 * term_yyy;
                        _olk(0, Cart::yyyyz) = PmB2 * _olk(0, Cart::yyyy);
                        _olk(0, Cart::yyyzz) = PmB2 * _olk(0, Cart::yyyz) + term_yyy;
                        _olk(0, Cart::yyzzz) = PmB1 * _olk(0, Cart::yzzz) + term_zzz;
                        _olk(0, Cart::yzzzz) = PmB1 * _olk(0, Cart::zzzz);
                        _olk(0, Cart::zzzzz) = PmB2 * _olk(0, Cart::zzzz) + 4.0 * term_zzz;
                        //------------------------------------------------------

                        //Integrals     p - h     d - h     f - h     g - h     h - h     i - h
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            COMPLEX term_xxx = _cfak * _olk(_i, Cart::xxx);
                            COMPLEX term_yyy = _cfak * _olk(_i, Cart::yyy);
                            COMPLEX term_zzz = _cfak * _olk(_i, Cart::zzz);
                            _olk(_i, Cart::xxxxx) = PmB0 * _olk(_i, Cart::xxxx) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xxxx) + 4.0 * term_xxx;
                            _olk(_i, Cart::xxxxy) = PmB1 * _olk(_i, Cart::xxxx) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxx);
                            _olk(_i, Cart::xxxxz) = PmB2 * _olk(_i, Cart::xxxx) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxxx);
                            _olk(_i, Cart::xxxyy) = PmB1 * _olk(_i, Cart::xxxy) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxy) + term_xxx;
                            _olk(_i, Cart::xxxyz) = PmB1 * _olk(_i, Cart::xxxz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxz);
                            _olk(_i, Cart::xxxzz) = PmB2 * _olk(_i, Cart::xxxz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxxz) + term_xxx;
                            _olk(_i, Cart::xxyyy) = PmB0 * _olk(_i, Cart::xyyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xyyy) + term_yyy;
                            _olk(_i, Cart::xxyyz) = PmB2 * _olk(_i, Cart::xxyy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxyy);
                            _olk(_i, Cart::xxyzz) = PmB1 * _olk(_i, Cart::xxzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxzz);
                            _olk(_i, Cart::xxzzz) = PmB0 * _olk(_i, Cart::xzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xzzz) + term_zzz;
                            _olk(_i, Cart::xyyyy) = PmB0 * _olk(_i, Cart::yyyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyyy);
                            _olk(_i, Cart::xyyyz) = PmB0 * _olk(_i, Cart::yyyz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyyz);
                            _olk(_i, Cart::xyyzz) = PmB0 * _olk(_i, Cart::yyzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyzz);
                            _olk(_i, Cart::xyzzz) = PmB0 * _olk(_i, Cart::yzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yzzz);
                            _olk(_i, Cart::xzzzz) = PmB0 * _olk(_i, Cart::zzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::zzzz);
                            _olk(_i, Cart::yyyyy) = PmB1 * _olk(_i, Cart::yyyy) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yyyy) + 4.0 * term_yyy;
                            _olk(_i, Cart::yyyyz) = PmB2 * _olk(_i, Cart::yyyy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::yyyy);
                            _olk(_i, Cart::yyyzz) = PmB2 * _olk(_i, Cart::yyyz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::yyyz) + term_yyy;
                            _olk(_i, Cart::yyzzz) = PmB1 * _olk(_i, Cart::yzzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yzzz) + term_zzz;
                            _olk(_i, Cart::yzzzz) = PmB1 * _olk(_i, Cart::zzzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::zzzz);
                            _olk(_i, Cart::zzzzz) = PmB2 * _olk(_i, Cart::zzzz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::zzzz) + 4.0 * term_zzz;
                        }
                        //------------------------------------------------------

                    } // end if (_lmax_col > 4)

                    if (_lmax_col > 5) {

                        //Integrals     s - W -i
                        COMPLEX term_xxxx = _cfak * _olk(0, Cart::xxxx);
                        COMPLEX term_xyyy = _cfak * _olk(0, Cart::xyyy);
                        COMPLEX term_xzzz = _cfak * _olk(0, Cart::xzzz);
                        COMPLEX term_yyyy = _cfak * _olk(0, Cart::yyyy);
                        COMPLEX term_yyzz = _cfak * _olk(0, Cart::yyzz);
                        COMPLEX term_yzzz = _cfak * _olk(0, Cart::yzzz);
                        COMPLEX term_zzzz = _cfak * _olk(0, Cart::zzzz);
                        _olk(0, Cart::xxxxxx) = PmB0 * _olk(0, Cart::xxxxx) + 5.0 * term_xxxx;
                        _olk(0, Cart::xxxxxy) = PmB1 * _olk(0, Cart::xxxxx);
                        _olk(0, Cart::xxxxxz) = PmB2 * _olk(0, Cart::xxxxx);
                        _olk(0, Cart::xxxxyy) = PmB1 * _olk(0, Cart::xxxxy) + term_xxxx;
                        _olk(0, Cart::xxxxyz) = PmB1 * _olk(0, Cart::xxxxz);
                        _olk(0, Cart::xxxxzz) = PmB2 * _olk(0, Cart::xxxxz) + term_xxxx;
                        _olk(0, Cart::xxxyyy) = PmB0 * _olk(0, Cart::xxyyy) + 2.0 * term_xyyy;
                        _olk(0, Cart::xxxyyz) = PmB2 * _olk(0, Cart::xxxyy);
                        _olk(0, Cart::xxxyzz) = PmB1 * _olk(0, Cart::xxxzz);
                        _olk(0, Cart::xxxzzz) = PmB0 * _olk(0, Cart::xxzzz) + 2.0 * term_xzzz;
                        _olk(0, Cart::xxyyyy) = PmB0 * _olk(0, Cart::xyyyy) + term_yyyy;
                        _olk(0, Cart::xxyyyz) = PmB2 * _olk(0, Cart::xxyyy);
                        _olk(0, Cart::xxyyzz) = PmB0 * _olk(0, Cart::xyyzz) + term_yyzz;
                        _olk(0, Cart::xxyzzz) = PmB1 * _olk(0, Cart::xxzzz);
                        _olk(0, Cart::xxzzzz) = PmB0 * _olk(0, Cart::xzzzz) + term_zzzz;
                        _olk(0, Cart::xyyyyy) = PmB0 * _olk(0, Cart::yyyyy);
                        _olk(0, Cart::xyyyyz) = PmB0 * _olk(0, Cart::yyyyz);
                        _olk(0, Cart::xyyyzz) = PmB0 * _olk(0, Cart::yyyzz);
                        _olk(0, Cart::xyyzzz) = PmB0 * _olk(0, Cart::yyzzz);
                        _olk(0, Cart::xyzzzz) = PmB0 * _olk(0, Cart::yzzzz);
                        _olk(0, Cart::xzzzzz) = PmB0 * _olk(0, Cart::zzzzz);
                        _olk(0, Cart::yyyyyy) = PmB1 * _olk(0, Cart::yyyyy) + 5.0 * term_yyyy;
                        _olk(0, Cart::yyyyyz) = PmB2 * _olk(0, Cart::yyyyy);
                        _olk(0, Cart::yyyyzz) = PmB2 * _olk(0, Cart::yyyyz) + term_yyyy;
                        _olk(0, Cart::yyyzzz) = PmB1 * _olk(0, Cart::yyzzz) + 2.0 * term_yzzz;
                        _olk(0, Cart::yyzzzz) = PmB1 * _olk(0, Cart::yzzzz) + term_zzzz;
                        _olk(0, Cart::yzzzzz) = PmB1 * _olk(0, Cart::zzzzz);
                        _olk(0, Cart::zzzzzz) = PmB2 * _olk(0, Cart::zzzzz) + 5.0 * term_zzzz;
                        //------------------------------------------------------

                        //Integrals     p - W - i     d - W - i     f - W - i     g - W -i     h - W - i     i - W - i
                        for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
                            COMPLEX term_xxxx = _cfak * _olk(_i, Cart::xxxx);
                            COMPLEX term_xyyy = _cfak * _olk(_i, Cart::xyyy);
                            COMPLEX term_xzzz = _cfak * _olk(_i, Cart::xzzz);
                            COMPLEX term_yyyy = _cfak * _olk(_i, Cart::yyyy);
                            COMPLEX term_yyzz = _cfak * _olk(_i, Cart::yyzz);
                            COMPLEX term_yzzz = _cfak * _olk(_i, Cart::yzzz);
                            COMPLEX term_zzzz = _cfak * _olk(_i, Cart::zzzz);
                            _olk(_i, Cart::xxxxxx) = PmB0 * _olk(_i, Cart::xxxxx) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xxxxx) + 5.0 * term_xxxx;
                            _olk(_i, Cart::xxxxxy) = PmB1 * _olk(_i, Cart::xxxxx) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxxx);
                            _olk(_i, Cart::xxxxxz) = PmB2 * _olk(_i, Cart::xxxxx) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxxxx);
                            _olk(_i, Cart::xxxxyy) = PmB1 * _olk(_i, Cart::xxxxy) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxxy) + term_xxxx;
                            _olk(_i, Cart::xxxxyz) = PmB1 * _olk(_i, Cart::xxxxz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxxz);
                            _olk(_i, Cart::xxxxzz) = PmB2 * _olk(_i, Cart::xxxxz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxxxz) + term_xxxx;
                            _olk(_i, Cart::xxxyyy) = PmB0 * _olk(_i, Cart::xxyyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xxyyy) + 2.0 * term_xyyy;
                            _olk(_i, Cart::xxxyyz) = PmB2 * _olk(_i, Cart::xxxyy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxxyy);
                            _olk(_i, Cart::xxxyzz) = PmB1 * _olk(_i, Cart::xxxzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxxzz);
                            _olk(_i, Cart::xxxzzz) = PmB0 * _olk(_i, Cart::xxzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xxzzz) + 2.0 * term_xzzz;
                            _olk(_i, Cart::xxyyyy) = PmB0 * _olk(_i, Cart::xyyyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xyyyy) + term_yyyy;
                            _olk(_i, Cart::xxyyyz) = PmB2 * _olk(_i, Cart::xxyyy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::xxyyy);
                            _olk(_i, Cart::xxyyzz) = PmB0 * _olk(_i, Cart::xyyzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xyyzz) + term_yyzz;
                            _olk(_i, Cart::xxyzzz) = PmB1 * _olk(_i, Cart::xxzzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::xxzzz);
                            _olk(_i, Cart::xxzzzz) = PmB0 * _olk(_i, Cart::xzzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::xzzzz) + term_zzzz;
                            _olk(_i, Cart::xyyyyy) = PmB0 * _olk(_i, Cart::yyyyy) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyyyy);
                            _olk(_i, Cart::xyyyyz) = PmB0 * _olk(_i, Cart::yyyyz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyyyz);
                            _olk(_i, Cart::xyyyzz) = PmB0 * _olk(_i, Cart::yyyzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyyzz);
                            _olk(_i, Cart::xyyzzz) = PmB0 * _olk(_i, Cart::yyzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yyzzz);
                            _olk(_i, Cart::xyzzzz) = PmB0 * _olk(_i, Cart::yzzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::yzzzz);
                            _olk(_i, Cart::xzzzzz) = PmB0 * _olk(_i, Cart::zzzzz) + double(nx[_i]) * _cfak * _olk(i_less_x[_i], Cart::zzzzz);
                            _olk(_i, Cart::yyyyyy) = PmB1 * _olk(_i, Cart::yyyyy) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yyyyy) + 5.0 * term_yyyy;
                            _olk(_i, Cart::yyyyyz) = PmB2 * _olk(_i, Cart::yyyyy) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::yyyyy);
                            _olk(_i, Cart::yyyyzz) = PmB2 * _olk(_i, Cart::yyyyz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::yyyyz) + term_yyyy;
                            _olk(_i, Cart::yyyzzz) = PmB1 * _olk(_i, Cart::yyzzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yyzzz) + 2.0 * term_yzzz;
                            _olk(_i, Cart::yyzzzz) = PmB1 * _olk(_i, Cart::yzzzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::yzzzz) + term_zzzz;
                            _olk(_i, Cart::yzzzzz) = PmB1 * _olk(_i, Cart::zzzzz) + double(ny[_i]) * _cfak * _olk(i_less_y[_i], Cart::zzzzz);
                            _olk(_i, Cart::zzzzzz) = PmB2 * _olk(_i, Cart::zzzzz) + double(nz[_i]) * _cfak * _olk(i_less_z[_i], Cart::zzzzz) + 5.0 * term_zzzz;
                        }
                        //------------------------------------------------------


                    } // end if (_lmax_col > 5)        

                    Eigen::MatrixXd _trafo_row = getTrafo(*itr);
                    Eigen::MatrixXd _trafo_col = getTrafo(*itc);

                    // cartesian -> spherical
                    Eigen::MatrixXcd _olk_sph=_trafo_row.transpose()*_olk*_trafo_col;
                   

                    // save to _matrix
                    for (int i = 0; i < _matrix.rows(); i++) {
                        for (int j = 0; j < _matrix.cols(); j++) {
                            _matrix(i, j) += _olk_sph(i + _shell_row->getOffset(), j + _shell_col->getOffset());
                        }
                    }


                } // close Gaussian _shell_col     

            } // close Gaussian _shell_row


        } // End AOPlanewave
        
        
        
        void AOPlanewave::Fillextpotential(const AOBasis& aobasis, const std::vector< tools::vec>& _kpoints) {
            
            _externalpotential = Eigen::MatrixXcd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

            for (const auto& kpoint:_kpoints) {
                    _aomatrix = Eigen::MatrixXcd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
                    setkVector(kpoint);
                    Fill(aobasis);
                    _externalpotential+=_aomatrix;     
                }
            
            return;
        }   

    }
}
