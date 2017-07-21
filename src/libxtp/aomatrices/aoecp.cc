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
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <votca/tools/linalg.h>
#include <votca/xtp/elements.h>






namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    using namespace votca::tools;
    
    
    void AOECP::FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, const AOShell* _shell_row, const AOShell* _shell_col, AOBasis* ecp) {

            // get shell positions
            
            int _lmax_row = _shell_row->getLmax();
            std::vector<double> _contractions_row_full((_lmax_row + 1)*(_lmax_row + 1));

            
            int _lmax_col = _shell_col->getLmax();
            std::vector<double> _contractions_col_full((_lmax_col + 1)*(_lmax_col + 1));

            
            
            const vec& _pos_row = _shell_row->getPos();
            const vec& _pos_col = _shell_col->getPos();
            const vec _diff = _pos_row - _pos_col;
            // initialize some helper
            double _distsq = _diff*_diff;


            // iterate over Gaussians in this _shell_row
            for (AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr) {
                // iterate over Gaussians in this _shell_col
                // get decay constant
                const double _decay_row = (*itr)->getDecay();

                //if ( _decay_row > 0.08 ) continue;

                const std::vector<double>& _contractions_row = (*itr)->getContraction();
                // shitty magic
                for (int L = 0; L <= _lmax_row; L++) {
                    for (int M = L * L; M < (L + 1)*(L + 1); M++) {
                        _contractions_row_full[M] = _contractions_row[L];
                    }
                }


                for (AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc) {
                    //get decay constant
                    const double _decay_col = (*itc)->getDecay();
                    const double _fak = 0.5 / (_decay_row + _decay_col);
                    const double _fak2 = 2.0 * _fak;


                    double _exparg = _fak2 * _decay_row * _decay_col *_distsq;

                    // check if distance between postions is big, then skip step   

                    if (_exparg > 30.0) {
                        continue;
                    }

                    const std::vector<double>& _contractions_col = (*itc)->getContraction();
                    for (int L = 0; L <= _lmax_col; L++) {
                        for (int M = L * L; M < (L + 1)*(L + 1); M++) {
                            _contractions_col_full[M] = _contractions_col[L];
                        }
                    }
                    // for each atom and its pseudopotential, get a matrix
                    int _atomidx = 0;

                    ub::matrix<int> _power_matrix = ub::zero_matrix<int>(5, 4); ///// // max 12 fit components, max non-local ECP l=0,1,2,3 ///////////////
                    ub::matrix<double> _decay_matrix = ub::zero_matrix<double>(5, 4); ///// // max 12 fit components, max non-local ECP l=0,1,2,3
                    ub::matrix<double> _coef_matrix = ub::zero_matrix<double>(5, 4); ///// // max 12 fit components, max non-local ECP l=0,1,2,3

                    AOBasis::AOShellIterator final_iter = ecp->lastShell();
                    --final_iter;
                    vec _ecp_eval_pos = vec(0.0);
                    int _lmax_ecp_act = 0;
                    for (AOBasis::AOShellIterator _ecp = ecp->firstShell(); _ecp != ecp->lastShell(); ++_ecp) {

                        const AOShell* _shell_ecp = ecp->getShell(_ecp);
                        const vec& _ecp_pos = _shell_ecp->getPos();

                        int this_atom = _shell_ecp->getIndex();

                        const int _ecp_l = _shell_ecp->getOffset(); //  angular momentum l is stored in offset for ECP

                        // only do the non-local parts
                        if (_ecp_l < _shell_ecp->getNumFunc()) {
                            int _lmax_ecp_old = _lmax_ecp_act; ///////
                            _lmax_ecp_act = _shell_ecp->getNumFunc() - 1; ///////
                            int i_fit = -1;
                            for (AOShell::GaussianIterator itecp = _shell_ecp->firstGaussian(); itecp != _shell_ecp->lastGaussian(); ++itecp) {
                                i_fit++;

                                // get info for this angular momentum shell
                                const int _power_ecp = (*itecp)->getPower(); ///////////////
                                const double _decay_ecp = (*itecp)->getDecay();
                                const double _contraction_ecp = (*itecp)->getContraction()[0];

                                // collect atom ECP
                                if (this_atom == _atomidx) {
                                    _ecp_eval_pos = _ecp_pos;
                                    _power_matrix(i_fit, _ecp_l) = _power_ecp; ////////////////////
                                    _decay_matrix(i_fit, _ecp_l) = _decay_ecp;
                                    _coef_matrix(i_fit, _ecp_l) = _contraction_ecp;

                                }
                                if ((this_atom != _atomidx) || (_ecp == final_iter)) {

                                    // evaluate collected data, returns a (10x10) matrix of already normalized matrix elements
                                    ub::matrix<double> VNL_ECP = calcVNLmatrix(_lmax_ecp_old, _ecp_eval_pos, *itr, *itc,  _power_matrix,_decay_matrix, _coef_matrix); ////////////////

                                    // consider contractions
                                    // cut out block that is needed. sum
                                    //                                   cout << "_matrix.size1,2()   " << _matrix.size1() << "    " << _matrix.size2() << endl;
                                    for (unsigned i = 0; i < _matrix.size1(); i++) {
                                        for (unsigned j = 0; j < _matrix.size2(); j++) {
                                            _matrix(i, j) += VNL_ECP(i + _shell_row->getOffset(), j + _shell_col->getOffset()) * _contractions_row_full[i + _shell_row->getOffset()] * _contractions_col_full[j + _shell_col->getOffset()];
                                        }
                                    }


                                    // reset atom ECP containers
                                    _power_matrix = ub::zero_matrix<int>(5,4); ///// // max 12 fit components, max non-local ECP l=0,1,2 /////////////////////
                                    _decay_matrix = ub::zero_matrix<double>(5, 4); ///// // max 12 fit components, max non-local ECP l=0,1,2
                                    _coef_matrix = ub::zero_matrix<double>(5, 4); ///// // max 12 fit components, max non-local ECP l=0,1,2
                                    _atomidx++;
                                    i_fit = 0;
                                    //cout << "setting new matrix " << i_fit << " l " << _ecp_l << " alpha  " << _decay_ecp <<  " pref " << _contraction_ecp << endl;
                                    _power_matrix(i_fit, _ecp_l) = _power_ecp; //////////////////
                                    _decay_matrix(i_fit, _ecp_l) = _decay_ecp;
                                    _coef_matrix(i_fit, _ecp_l) = _contraction_ecp;
                                } // evaluate if new atom is found

                            } // all Gaussians in ecp_shell
                        } // only for non local parts

                    } // all ecp_shells

                }// _shell_col Gaussians
            }// _shell_row Gaussians

            return;
        }

        ub::matrix<double> AOECP::calcVNLmatrix(int _lmax_ecp, const vec& posC, const AOGaussianPrimitive* _g_row, const AOGaussianPrimitive* _g_col,const ub::matrix<int>& _power_ecp, const ub::matrix<double>& _gamma_ecp,const ub::matrix<double>& _pref_ecp) {


            const double pi = boost::math::constants::pi<double>();
            const double conv = 1.e-8;
            /* calculate the contribution of the nonlocal 
             *     ECP of atom at posC with 
             *       decay constants in _gamma_ecp
             *       coefficients in    _pref_ecp
             *       with angular momentum of max 3
             * 
             * to DFT basis shell pair 
             *    with decay alpha at posA
             *         decay beta  at posB
  
             */

            double SQ2 = sqrt(2.);
            double SQ3 = sqrt(3.);
            double SQ5 = sqrt(5.);
            double SQ7 = sqrt(7.);

            
            double alpha = _g_row->getDecay();
            double beta = _g_col->getDecay();
            const vec& posA = _g_row->getShell()->getPos();
            const vec& posB = _g_col->getShell()->getPos();
            int _lmax_row = _g_row->getShell()->getLmax();
            int _lmax_col = _g_col->getShell()->getLmax();
            
            int _lmax_dft = max(_lmax_row, _lmax_col);
            int _lmax_dft_ecp = max(_lmax_dft, _lmax_ecp);
            int _nsph_row = (_lmax_row + 1) * (_lmax_row + 1);
            int _nsph_col = (_lmax_col + 1) * (_lmax_col + 1);
            //double DIL = (alpha + beta + _gamma_ecp(1, 2));
            //ub::vector<double> Int_r_exp = CalcInt_r_exp(3, DIL);

            ub::matrix<double> matrix = ub::zero_matrix<double>(_nsph_row,_nsph_col);
            const int nnonsep = _gamma_ecp.size1();
            const int nmax = 90;
            ub::matrix<double> XI(_lmax_ecp + 1, nmax);

            double f_even_r0 = .5 * sqrt(2.) * sqrt(.5 * pi);
            double f_even_r1 = .5;
            double f_even_r2 = .5 * f_even_r0;
            double f_odd_r0 = .5;
            double f_odd_r1 = .25 * sqrt(2.) * sqrt(.5 * pi);
            double f_odd_r2 = 1.;
            double DGAMAF_r0;
            double DGAMAF_r1;
            double DGAMAF_r2;

            for (int N = 0; N < nmax ; N++) {

                if ((N % 2) == 0) { // N even

                    if (N > 0) {
                        f_even_r0 = f_even_r2;
                        f_even_r1 = f_even_r1 * double(N / 2);
                        f_even_r2 = .5 * f_even_r0 * double(N + 1);
                    }
                    DGAMAF_r0 = f_even_r0;
                    DGAMAF_r1 = f_even_r1;
                    DGAMAF_r2 = f_even_r2;

                } else { // N odd

                    if (N > 1) {
                        f_odd_r0 = f_odd_r2;
                        f_odd_r1 = .5 * f_odd_r1 * double(N);
                        f_odd_r2 = f_odd_r0 * double((N + 1) / 2);
                    }
                    DGAMAF_r0 = f_odd_r0;
                    DGAMAF_r1 = f_odd_r1;
                    DGAMAF_r2 = f_odd_r2;
                }

                double DFAK_r0 = .5 * double(N + 1);
                double DFAK_r1 = .5 * double(N + 2);
                double DFAK_r2 = .5 * double(N + 3);

                for (int L = 0; L <= _lmax_ecp; L++) {

                    XI(L, N) = 0.0;
                    for (int I = 0; I < nnonsep; I++) {
                        int power = _power_ecp(I, L);
                        double DLI = (alpha + beta + _gamma_ecp(I, L));
                        if (power == 2) {
                            XI(L, N) += DGAMAF_r2 * _pref_ecp(I, L) / pow(DLI, DFAK_r2); // r^2 terms
                        } else if (power == 0) {
                            XI(L, N) += DGAMAF_r0 * _pref_ecp(I, L) / pow(DLI, DFAK_r0); // r^0 terms
                        } else if (power == 1) {
                            XI(L, N) += DGAMAF_r1 * _pref_ecp(I, L) / pow(DLI, DFAK_r1); // r^1 terms
                        }
                    }

                }

            }

            /****** ORIGINAL CKO SUBROUTINE **********/
            // get a multi dimensional array
            typedef boost::multi_array<double, 4> ma_type;
            typedef boost::multi_array_types::extent_range range;
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type COEF;
            COEF.resize(extents[ range(0, 4) ][ range(0, 4) ][ range(0, 7)][range(0, 42)]);

            // init it all to 0
            for (index i1 = 0; i1 < 4; i1++) { //////////
                for (index i2 = 0; i2 < 4; i2++) { ///////////////
                    for (index i3 = 0; i3 < 7; i3++) {
                        for (index i4 = 0; i4 < 42; i4++) {
                            COEF[i1][i2][i3][i4] = 0.0;
                        }
                    }
                }
            }
            for (index i4 = 0; i4 < 42; i4++) {

                int NU = i4 % 2; ///
                int NG = (i4 + 1) % 2;
                double FN1 = double(i4 + 1);
                double FN2 = double(i4 + 2);
                double FN3 = double(i4 + 3);
                double FN4 = double(i4 + 4);
                double FN5 = double(i4 + 5);
                double FN6 = double(i4 + 6);
                double FN7 = double(i4 + 7);

                COEF[0][0][3][i4] = NG / FN1; /////////   M0(x)

                if (_lmax_dft_ecp > 0) {

                    COEF[0][1][3][i4] = NU / FN2*SQ3; ////////  SQ(3) * M1(x)
                    COEF[1][0][3][i4] = COEF[0][1][3][i4];

                    COEF[1][1][3][i4] = NG * 3.0 / FN3; ///////    M0(x) + 2 * M2(x)
                    COEF[1][1][4][i4] = 3.0 / 2.0 * NG * (1.0 / FN1 - 1.0 / FN3); ////////    M0(x) - M2(x) 
                    COEF[1][1][2][i4] = COEF[1][1][4][i4];

                }

                if (_lmax_dft_ecp > 1) {

                    COEF[0][2][3][i4] = NG / 2.0 * SQ5 * (3.0 / FN3 - 1.0 / FN1); //////   SQ(5) * M2(x)
                    COEF[2][0][3][i4] = COEF[0][2][3][i4];

                    COEF[1][2][3][i4] = SQ3 * SQ5 / 2.0 * NU * (3.0 / FN4 - 1.0 / FN2); ///////   (2/5) * SQ(15) * ( M1(x) + (3/2) * M3(x) )
                    COEF[2][1][3][i4] = COEF[1][2][3][i4];
                    COEF[1][2][4][i4] = 3. * SQ5 / 2.0 * NU * (1.0 / FN2 - 1.0 / FN4); ///////   (SQ(45)/5) * ( M1(x) - M3(x) )
                    COEF[1][2][2][i4] = COEF[1][2][4][i4];
                    COEF[2][1][4][i4] = COEF[1][2][4][i4];
                    COEF[2][1][2][i4] = COEF[1][2][4][i4];

                    COEF[2][2][3][i4] = 5.0 / 4.0 * NG * (9.0 / FN5 - 6.0 / FN3 + 1.0 / FN1); ///////  M0(x) + (10/7) * M2(x) + (18/7) * M4(x)
                    COEF[2][2][4][i4] = NG * 15.0 / 2.0 * (1.0 / FN3 - 1.0 / FN5); ///////  M0(x) + (5/7) * M2(x) - (12/7) * M4(x)
                    COEF[2][2][2][i4] = COEF[2][2][4][i4];
                    COEF[2][2][5][i4] = 15.0 / 8.0 * NG * (1.0 / FN1 - 2.0 / FN3 + 1.0 / FN5); ///////  M0(x) - (10/7) * M2(x) + (3/7) * M4(x) 
                    COEF[2][2][1][i4] = COEF[2][2][5][i4];
                }

                if (_lmax_dft_ecp > 2) {

                    ////         cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& COEF[][][]  FFFFFFFFFF" << endl;
                    double COEFF = NU * .5 * SQ7 * (5. / FN4 - 3. / FN2); ///  SQ(7) * M3(x)
                    COEF[0][3][3][i4] = COEFF;
                    COEF[3][0][3][i4] = COEFF;

                    COEFF = NG * .5 * SQ3 * SQ7 * (5. / FN5 - 3. / FN3); ///  SQ(3/7) * ( 3 * M2(x) + 4 * M4(x) )
                    COEF[1][3][3][i4] = COEFF;
                    COEF[3][1][3][i4] = COEFF;
                    COEFF = NG * .375 * SQ2 * SQ7 * (-5. / FN5 + 6. / FN3 - 1. / FN1); ///  3 * SQ(2/7) * ( M2(x) - M4(x) )
                    COEF[1][3][2][i4] = COEFF;
                    COEF[1][3][4][i4] = COEFF;
                    COEF[3][1][2][i4] = COEFF;
                    COEF[3][1][4][i4] = COEFF;

                    COEFF = NU * .25 * SQ5 * SQ7 * (15. / FN6 - 14. / FN4 + 3. / FN2); ///  ( 1/(3*SQ(35)) ) * ( 27 * M1(x) + 28 * M3(x) + 50 * M5(x) )
                    COEF[2][3][3][i4] = COEFF;
                    COEF[3][2][3][i4] = COEFF;
                    COEFF = NU * .375 * SQ2 * SQ5 * SQ7 * (-5. / FN6 + 6. / FN4 - 1. / FN2); ///  ( SQ(2/35)/3 ) * ( 18 * M1(x) + 7 * M3(x) - 25 * M5(x) )
                    COEF[2][3][2][i4] = COEFF;
                    COEF[2][3][4][i4] = COEFF;
                    COEF[3][2][2][i4] = COEFF;
                    COEF[3][2][4][i4] = COEFF;
                    COEFF = NU * 1.875 * SQ7 * (1. / FN6 - 2. / FN4 + 1. / FN2); ///  ( 1/(3*SQ(7)) ) * ( 9 * M1(x) - 14 * M3(x) + 5 * M5(x) )
                    COEF[2][3][1][i4] = COEFF;
                    COEF[2][3][5][i4] = COEFF;
                    COEF[3][2][1][i4] = COEFF;
                    COEF[3][2][5][i4] = COEFF;

                    COEF[3][3][3][i4] = NG * 1.75 * (50. / FN7 - 30. / FN5 + 9. / FN3); ///  (1/33) * ( 33 * M0(x) + 44 * M2(x) + 54 * M4(x) + 100 * M6(x) )
                    COEFF = NG * 1.3125 * (-25. / FN7 + 35. / FN5 - 11. / FN3 + 1. / FN1); ///  (1/11) * ( 11 * M0(x) + 11 * M2(x) + 3 * M4(x) - 25 * M6(x) )
                    COEF[3][3][2][i4] = COEFF;
                    COEF[3][3][4][i4] = COEFF;
                    COEFF = NG * 105 * .125 * (1. / FN7 - 2. / FN5 + 1. / FN3); ///  (1/11) * ( 11 * M0(x) - 21 * M4(x) + 10 * M6(x) )
                    COEF[3][3][1][i4] = COEFF;
                    COEF[3][3][5][i4] = COEFF;
                    COEFF = NG * 35 * .0625 * (-1. / FN7 + 3. / FN5 - 3. / FN3 + 1. / FN1); ///  (1/33) * ( 33 * M0(x) - 55 * M2(x) + 27 * M4(x) - 5 * M6(x) )
                    COEF[3][3][0][i4] = COEFF;
                    COEF[3][3][6][i4] = COEFF;

                }

            } // i4 loop (== CKO )

            /**** PREPARATIONS DONE, NOW START ******/
            vec AVS = posA - posC;
            vec BVS = posB - posC;

            double AVS2 = AVS*AVS;
            double BVS2 = BVS*BVS;

            double AVSSQ = sqrt(AVS2);
            double BVSSQ = sqrt(BVS2);
            double GAUSS = exp(-alpha * AVS2 - beta * BVS2);

            // some limit determinations
            double G1 = exp(-alpha * AVS2);
            double G2 = exp(-beta * BVS2);

            int NMAX1 = 0;
            int NMAX2 = 0;

            if (AVSSQ <= 0.1) {

                NMAX1 = 1;

            } else {

                double AMAX = 0.0;
                double fak = 2.0 * alpha*AVSSQ;
                double Pow = 1;
                double factorialNN = 1;
                for (int NN = 0; NN < 42; NN++) {

                    if (NN != 0) {
                        Pow = Pow*fak;
                        factorialNN = factorialNN*NN;
                    }

                    double AF = Pow / factorialNN*G1;
                    double AF1 = std::abs(AF * XI(0, NN + 4));
                    double AF2 = std::abs(AF * XI(1, NN + 4));
                    double AF3 = std::abs(AF * XI(2, NN + 4));
                    AMAX = std::max(AF1, AF2);
                    AMAX = std::max(AMAX, AF3);

                    if (NMAX1 == 0 && AMAX <= conv) NMAX1 = NN + 1;
                    if (NMAX1 != 0 && AMAX > conv) NMAX1 = 0;

                }
                if (NMAX1 == 0 && AMAX > conv) NMAX1 = 42;
            }

            // same story for B
            if (BVSSQ <= 0.1) {

                NMAX2 = 1;

            } else {

                double BMAX = 0.0;
                double fak = 2.0 * beta*BVSSQ;
                double Pow = 1;
                double factorialNN = 1;
                for (int NN = 0; NN < 42; NN++) {

                    if (NN != 0) {
                        Pow = Pow*fak;
                        factorialNN = factorialNN*NN;
                    }
                    double BF = Pow / factorialNN*G2;
                    double BF1 = std::abs(BF * XI(0, NN + 4));
                    double BF2 = std::abs(BF * XI(1, NN + 4));
                    double BF3 = std::abs(BF * XI(2, NN + 4));
                    BMAX = std::max(BF1, BF2);
                    BMAX = std::max(BMAX, BF3);

                    if (NMAX2 == 0 && BMAX <= conv) NMAX2 = NN + 1;
                    if (NMAX2 != 0 && BMAX > conv) NMAX2 = 0;

                }

                if (NMAX2 == 0 && BMAX > conv) NMAX2 = 42;
            }


            int INULL = 1;
            if (AVSSQ <= 1.e-1) INULL = 3;
            if (BVSSQ <= 1.e-1) INULL = 4;
            if (AVSSQ <= 1.e-1 && BVSSQ <= 1.e-1) INULL = 2;

            type_3D BLMA;
            type_3D CA;
            getBLMCOF(_lmax_ecp, _lmax_row, AVS, BLMA, CA);


            type_3D BLMB;
            type_3D CB;
            getBLMCOF(_lmax_ecp, _lmax_col, BVS, BLMB, CB);

            typedef boost::multi_array_types::extent_range range;
            typedef type_3D::index index;
            type_3D::extent_gen extents3D;

            type_3D CC;
            CC.resize(extents3D[ range(0, 4)][range(0, 7)][range(0, 7)]); ///////
            for (index L = 0;L <= _lmax_ecp; L++) { ////
                for (index M1 = 3 - L; M1 <= 3 + L; M1++) {
                    for (index M2 = 3 - L; M2 <= 3 + L; M2++) {

                        CC[L][M1][M2] = 0.0;
                        for (index M = 3 - L; M <= 3 + L; M++) {

                            CC[L][M1][M2] += CA[L][M][M1] * CB[L][M][M2]; /////

                        }
                    }
                }
            }

            typedef boost::multi_array<double, 5> type_5D;
            type_5D::extent_gen extents5D;
            type_5D SUMCI;
            SUMCI.resize(extents5D[range(0, 4)][range(0, 4)][ range(0, 4)][range(0, 7)][range(0, 7)]);
            type_3D SUMCI3;
            SUMCI3.resize(extents3D[range(0, 4)][range(0, 4)][range(0, 7)]);


            switch (INULL) {

                case 1:
                {

                    for (index L = 0; L <= _lmax_ecp; L++) {
                        for (index L1 = 0; L1 <= _lmax_row; L1++) { /////                
                            for (index L2 = 0; L2 <= _lmax_col; L2++) { /////          
                                for (index M1 = 3 - L; M1 <= 3 + L; M1++) {
                                    for (index M2 = 3 - L; M2 <= 3 + L; M2++) {

                                        SUMCI[L][L1][L2][M1][M2] = 0.0;

                                        double fak1 = 2.0 * alpha*AVSSQ;
                                        double pow1 = 1;
                                        double factorialN = 1;

                                        for (int N = 0; N < NMAX1-1; N++) { ///

                                            if (N != 0) {
                                                pow1 = pow1*fak1;
                                                factorialN = factorialN*N;
                                            }

                                            double VAR1 = COEF[L][L1][M1][N] * pow1 / factorialN; ///////////
                                            double VAR2 = 0.0;
                                            double fak2 = 2.0 * beta*BVSSQ;
                                            double pow2 = 1;
                                            double factorialNN = 1;

                                            for (int NN = 0; NN < NMAX2-1; NN++) {

                                                if (NN != 0) {
                                                    pow2 = pow2*fak2;
                                                    factorialNN = factorialNN*NN;
                                                }
                                                double XDUM = COEF[L][L2][M2][NN] * pow2 / factorialNN; //////
                                                VAR2 += XDUM * XI(L, N + NN + L1 + L2); // L index of XI starts with 0 !! /////////

                                            }

                                            SUMCI[L][L1][L2][M1][M2] += VAR1*VAR2; /////

                                        }

                                    } // end M2
                                } // end M1
                            } // end L2
                        } // end L1 
                    } // end L

                    // now finally calculate matrix

                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {

                            for (index L = 0; L <= _lmax_ecp; L++) {
                                for (index L1 = 0; L1 <= _lmax_row; L1++) {
                                    for (index L2 = 0; L2 <= _lmax_col; L2++) {

                                        for (index M1 = 3 - L; M1 <= 3 + L; M1++) {
                                            for (index M2 = 3 - L; M2 <= 3 + L; M2++) {
                                                matrix(i, j) += BLMA[i][L1][M1] * BLMB[j][L2][M2] * SUMCI[L][L1][L2][M1][M2] * CC[L][M1][M2];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }

                    break;
                }


                case 2: //  AVSSQ<=1.e-1 && BVSSQ <= 1.e-1
                {
                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {
                            for (index L = 0; L <= _lmax_ecp; L++) {
                                double XI_L = XI(L, L + L);
                                for (index M = 3 - L; M <= 3 + L; M++) {
                                    matrix(i, j) += BLMA[i][L][M] * BLMB[j][L][M] * XI_L;
                                }
                            }
                        }
                    }
                    break;
                }



                case 3: //  AVSSQ <= 1.e-1
                {

                    for (index L = 0; L <= _lmax_ecp; L++) { /// for (index L = 0; L <= 2; L++) {
                        for (index L2 = 0; L2 <= _lmax_col; L2++) {
                            for (index M2 = 3 - L; M2 <= 3 + L; M2++) {

                                double VAR2 = 0.0;
                                double fak = 2.0 * beta*BVSSQ;
                                double pow = 1;
                                double factorialNN = 1;

                                for (int NN = 0; NN <= NMAX2 - 1; NN++) {

                                    if (NN != 0) {
                                        pow = pow*fak;
                                        factorialNN = factorialNN*NN;
                                    }

                                    double XDUM = COEF[L][L2][M2][NN] * pow / factorialNN;
                                    VAR2 += XDUM * XI(L, NN + L + L2);

                                }

                                SUMCI3[L][L2][M2] = VAR2;

                            } // end M2
                        } // end L2
                    } // end L


                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {

                            for (index L = 0; L <= _lmax_ecp; L++) { /// for (index L = 0; L <= 2; L++) {
                                for (index L2 = 0; L2 <= _lmax_col; L2++) {

                                    for (index M1 = 3 - L; M1 <= 3 + L; M1++) {
                                        for (index M2 = 3 - L; M2 <= 3 + L; M2++) {

                                            matrix(i, j) += BLMA[i][L][M1] * BLMB[j][L2][M2] * SUMCI3[L][L2][M2] * CC[L][M1][M2];

                                        }
                                    }
                                }
                            }
                        }
                    }

                    break;
                }


                case 4: //  BVSSQ <= 1.e-1
                {

                    for (index L = 0; L <= _lmax_ecp; L++) { /// for (index L = 0; L <= 2; L++) {
                        for (index L1 = 0; L1 <= _lmax_row; L1++) {
                            for (index M1 = 3 - L; M1 <= 3 + L; M1++) {

                                double VAR1 = 0.0;
                                double fak = 2.0 * alpha*AVSSQ;
                                double pow = 1;
                                double factorialN = 1;

                                for (int N = 0; N <= NMAX1 - 1; N++) {

                                    if (N != 0) {
                                        pow = pow*fak;
                                        factorialN = factorialN*N;
                                    }

                                    double XDUM = COEF[L][L1][M1][N] * pow / factorialN;
                                    VAR1 += XDUM * XI(L, N + L1 + L);

                                }

                                SUMCI3[L][L1][M1] = VAR1;

                            } // end M1
                        } // end L1
                    } // end L


                    for (int i = 0; i < _nsph_row; i++) {
                        for (int j = 0; j < _nsph_col; j++) {

                            for (index L = 0; L <= _lmax_ecp; L++) { /// for (index L = 0; L <= 2; L++) {
                                for (index L1 = 0; L1 <= _lmax_row; L1++) {

                                    for (index M1 = 3 - L; M1 <= 3 + L; M1++) {
                                        for (index M2 = 3 - L; M2 <= 3 + L; M2++) {

                                            matrix(i, j) += BLMA[i][L1][M1] * BLMB[j][L][M2] * SUMCI3[L][L1][M1] * CC[L][M1][M2];

                                        }
                                    }
                                }
                            }
                        }
                    }

                    break;
                }


                default:
                    cout << "Wrong ECP summation mode";
                    exit(1);
            } // switch


            // GET TRAFO HERE ALREADY         
            ub::vector<double> NormA = CalcNorms(alpha,_nsph_row);
            ub::vector<double> NormB = CalcNorms(beta,_nsph_col);


            for (int i = 0; i < _nsph_row; i++) {
                for (int j = 0; j < _nsph_col; j++) {

                    matrix(i, j) = matrix(i, j) * GAUSS * NormA[i] * NormB[j];

                }
            }


            return matrix;
        }

        ub::vector<double> AOECP::CalcNorms(double decay, int size) {
            ub::vector<double> Norms = ub::vector<double>(size);
            const double PI = boost::math::constants::pi<double>();
            double SQ3 = sqrt(3.);

            double Norm_S = pow(2.0 * decay / PI, 0.75);
            double Norm_P;
            Norms[0] = Norm_S; ///  Y 00

            if (size > 1) {
                Norm_P = 2.0 * sqrt(decay) * Norm_S;
                Norms[1] = Norm_P; ///  Y 10
                Norms[2] = Norm_P; ///  Y 1-1
                Norms[3] = Norm_P; ///  Y 11
            }

            if (size > 4) {
                double Norm_D = 4.00 * decay * Norm_S;
                Norms[4] = .5 * Norm_D / SQ3; ///  Y 20
                Norms[5] = Norm_D; ///  Y 2-1
                Norms[6] = Norm_D; ///  Y 21
                Norms[7] = Norm_D; ///  Y 2-2
                Norms[8] = .5 * Norm_D; ///  Y 22
            }

            if (size > 9) {

                double SQ2 = sqrt(2.);
                double SQ5 = sqrt(5.);
                double Norm_F = 4.00 * decay * Norm_P;
                double Norm_F_1 = .5 * Norm_F / (SQ2 * SQ5);
                double Norm_F_3 = .5 * Norm_F / (SQ2 * SQ3);
                Norms[9] = .5 * Norm_F / (SQ3 * SQ5); ///  Y 30
                Norms[10] = Norm_F_1; ///  Y 3-1
                Norms[11] = Norm_F_1; ///  Y 31
                Norms[12] = Norm_F; ///  Y 3-2
                Norms[13] = .5 * Norm_F; ///  Y 32
                Norms[14] = Norm_F_3; ///  Y 3-3
                Norms[15] = Norm_F_3; ///  Y 33
            }
            return Norms;
        }

        void AOECP::getBLMCOF(int _lmax_ecp, int _lmax_dft, const vec& pos, type_3D& BLC, type_3D& C) {

            typedef boost::multi_array_types::extent_range range;
            typedef type_3D::index index;
            type_3D::extent_gen extents;

            int _nsph = (_lmax_dft + 1) * (_lmax_dft + 1);
            int _lmax_dft_ecp = max(_lmax_ecp, _lmax_dft);

            BLC.resize(extents[ range(0, _nsph) ][ range(0, 4) ][ range(0, 7)]);
            C.resize(extents[ range(0, 4) ][ range(0, 7) ][ range(0, 7)]);

            type_3D BLM;
            BLM.resize(extents[ range(0, _nsph) ][ range(0, 4) ][ range(0, 7)]);

            const double PI = boost::math::constants::pi<double>();
            double SPI = sqrt(PI);
            double S2PI = sqrt(2. * PI);
            double SQ2 = sqrt(2.);
            double SQ3 = sqrt(3.);
            double SQ5 = sqrt(5.);
            double XS = 2.0 * SPI; ///////       2 * SQ(pi)         
            double XP = XS / SQ3; /////////////   2 * SQ(pi/3)
            double XD = XP / SQ5; //////////     2 * SQ(pi/15)
            double XD_0 = 4.0 * SPI / SQ5; ///   4 * SQ(pi/5)         Y 20
            double XD_p2 = 2. * XD; //////         4 * SQ(pi/15)        Y 22

            for (index I = 0; I < _nsph; I++) {
                for (index L = 0; L <= _lmax_dft; L++) {
                    for (index M = 3 - L; M <= 3 + L; M++) {
                        BLM[I][L][M] = 0.0;
                    }
                }
            }

            double BVS_X = pos.getX();
            double BVS_Y = pos.getY();
            double BVS_Z = pos.getZ();
            double BVS_XX = BVS_X * BVS_X;
            double BVS_YY = BVS_Y * BVS_Y;
            double BVS_ZZ = BVS_Z * BVS_Z;
            double BVS_XY, BVS_XZ, BVS_YZ;

            BLM[0][0][3] = XS; ///  Y 00

            if (_lmax_dft > 0) {

                BLM[1][0][3] = -BVS_Z * XS; ///  Y 10
                BLM[1][1][3] = XP;

                BLM[2][0][3] = -BVS_Y * XS; ///  Y 1 -1
                BLM[2][1][2] = XP;

                BLM[3][0][3] = -BVS_X * XS; ///  Y 11
                BLM[3][1][4] = XP;

            }

            if (_lmax_dft > 1) {

                BVS_XY = BVS_X * BVS_Y;
                BVS_XZ = BVS_X * BVS_Z;
                BVS_YZ = BVS_Y * BVS_Z;

                BLM[4][0][3] = (2.0 * BVS_ZZ - BVS_XX - BVS_YY) * XS; ///  Y 20
                BLM[4][1][4] = 2.0 * BVS_X * XP;
                BLM[4][1][2] = 2.0 * BVS_Y * XP;
                BLM[4][1][3] = -4.0 * BVS_Z * XP;
                BLM[4][2][3] = XD_0;

                BLM[5][0][3] = BVS_YZ * XS; ///  Y 2 -1
                BLM[5][1][2] = -BVS_Z * XP;
                BLM[5][1][3] = -BVS_Y * XP;
                BLM[5][2][2] = XD;

                BLM[6][0][3] = BVS_XZ * XS; ///  Y 21
                BLM[6][1][4] = -BVS_Z * XP;
                BLM[6][1][3] = -BVS_X * XP;
                BLM[6][2][4] = XD;

                BLM[7][0][3] = BVS_XY * XS; ///  Y 2 -2
                BLM[7][1][4] = -BVS_Y * XP;
                BLM[7][1][2] = -BVS_X * XP;
                BLM[7][2][1] = XD;

                BLM[8][0][3] = (BVS_XX - BVS_YY) * XS; ///  Y 22
                BLM[8][1][4] = -2.0 * BVS_X * XP;
                BLM[8][1][2] = 2.0 * BVS_Y * XP;
                BLM[8][2][5] = XD_p2;
            }

            if (_lmax_dft > 2) {

                double SQ7 = sqrt(7.);
                double XF_0 = 4. * SPI / SQ7;
                double XF_3 = 4. * S2PI / (SQ5 * SQ7);
                double XF_m2 = 2. * SPI / (SQ3 * SQ5 * SQ7);
                double XF_p2 = 2. * XF_m2;
                double XF_1 = 4. * S2PI / (SQ3 * SQ7);

                BLM[9][0][3] = (3. * (BVS_XX + BVS_YY) - 2. * BVS_ZZ) * BVS_Z * XS; ///  Y 30
                BLM[9][1][2] = -6. * BVS_YZ * XP;
                BLM[9][1][3] = 3. * (2. * BVS_ZZ - BVS_XX - BVS_YY) * XP;
                BLM[9][1][4] = -6. * BVS_XZ * XP;
                BLM[9][2][2] = 6. * BVS_Y * XD;
                BLM[9][2][3] = -3. * BVS_Z * XD_0;
                BLM[9][2][4] = 6. * BVS_X * XD;
                BLM[9][3][3] = XF_0;

                BLM[10][0][3] = (BVS_XX + BVS_YY - 4. * BVS_ZZ) * BVS_Y * XS; ///  Y 3 -1
                BLM[10][1][2] = (4. * BVS_ZZ - BVS_XX - 3. * BVS_YY) * XP;
                BLM[10][1][3] = 8. * BVS_YZ * XP;
                BLM[10][1][4] = -2. * BVS_XY * XP;
                BLM[10][2][1] = 2. * BVS_X * XD;
                BLM[10][2][2] = -8. * BVS_Z * XD;
                BLM[10][2][3] = -2. * BVS_Y * XD_0;
                BLM[10][2][5] = -BVS_Y * XD_p2;
                BLM[10][3][2] = XF_1;

                BLM[11][0][3] = (BVS_XX + BVS_YY - 4. * BVS_ZZ) * BVS_X * XS; ///  Y 31
                BLM[11][1][2] = -2. * BVS_XY * XP;
                BLM[11][1][3] = 8. * BVS_XZ * XP;
                BLM[11][1][4] = (4. * BVS_ZZ - 3. * BVS_XX - BVS_YY) * XP;
                BLM[11][2][1] = 2. * BVS_Y * XD;
                BLM[11][2][3] = -2. * BVS_X * XD_0;
                BLM[11][2][4] = -8. * BVS_Z * XD;
                BLM[11][2][5] = BVS_X * XD_p2;
                BLM[11][3][4] = XF_1;

                BLM[12][0][3] = -BVS_XY * BVS_Z * XS; ///  Y 3 -2
                BLM[12][1][2] = BVS_XZ * XP;
                BLM[12][1][3] = BVS_XY * XP;
                BLM[12][1][4] = BVS_YZ * XP;
                BLM[12][2][1] = -BVS_Z * XD;
                BLM[12][2][2] = -BVS_X * XD;
                BLM[12][2][4] = -BVS_Y * XD;
                BLM[12][3][1] = XF_m2;

                BLM[13][0][3] = (BVS_YY - BVS_XX) * BVS_Z * XS; ///  Y 32
                BLM[13][1][2] = -2. * BVS_YZ * XP;
                BLM[13][1][3] = (BVS_XX - BVS_YY) * XP;
                BLM[13][1][4] = 2. * BVS_XZ * XP;
                BLM[13][2][2] = 2. * BVS_Y * XD;
                BLM[13][2][4] = -2. * BVS_X * XD;
                BLM[13][2][5] = -BVS_Z * XD_p2;
                BLM[13][3][5] = XF_p2;

                BLM[14][0][3] = (BVS_YY - 3. * BVS_XX) * BVS_Y * XS; ///  Y 3 -3
                BLM[14][1][2] = 3. * (BVS_XX - BVS_YY) * XP;
                BLM[14][1][4] = 6. * BVS_XY * XP;
                BLM[14][2][1] = -6. * BVS_X * XD;
                BLM[14][2][5] = -3. * BVS_Y * XD_p2;
                BLM[14][3][0] = XF_3;

                BLM[15][0][3] = (3. * BVS_YY - BVS_XX) * BVS_X * XS; ///  Y 33
                BLM[15][1][2] = -6. * BVS_XY * XP;
                BLM[15][1][4] = 3. * (BVS_XX - BVS_YY) * XP;
                BLM[15][2][1] = 6. * BVS_Y * XD;
                BLM[15][2][5] = -3. * BVS_X * XD_p2;
                BLM[15][3][6] = XF_3;

            }



            for (index L = 0; L <= _lmax_dft_ecp; L++) {
                for (index M = 3 - L; M <= 3 + L; M++) {
                    for (index MM = 3 - L; MM <= 3 + L; MM++) {
                        C[L][M][MM] = 0.0;
                    }
                }
            }
            double XY = BVS_XX + BVS_YY;
            double XYZ = XY + BVS_ZZ;
            double SXY = sqrt(XY); //// SXY = r * sin(theta)
            double SXYZ = sqrt(XYZ); //// SXYZ = r
            double CP = 1.0;
            double SP = 0.0;

            if (SXY > 1.e-4) {

                CP = BVS_X / SXY; //// CP = cos(phi)
                SP = BVS_Y / SXY; //// SP = sin(phi)

            }


            if (SXYZ > 1.e-4) {

                double CT = BVS_Z / SXYZ; /// CT = cos(theta)
                double ST = SXY / SXYZ; /// ST = sin(theta)

                C[0][3][3] = 1.0; // 2*SQ(pi) * (Y 00)

                if (_lmax_dft_ecp > 0) {

                    C[1][2][2] = CP; // 2*SQ(pi/3) * (Z 1-1)
                    C[1][2][3] = ST * SP; // 2*SQ(pi/3) * (Y 1-1)
                    C[1][2][4] = CT * SP; // 2*SQ(pi/3) * (Y 1-1)'

                    C[1][3][2] = 0.0;
                    C[1][3][3] = CT; // 2*SQ(pi/3) * Y 10
                    C[1][3][4] = -ST; // 2*SQ(pi/3) * (Y 10)'

                    C[1][4][2] = -SP; // 2*SQ(pi/3) * (Z 11)
                    C[1][4][3] = ST * CP; // 2*SQ(pi/3) * (Y 11)
                    C[1][4][4] = CT * CP; // 2*SQ(pi/3) * (Y 11)'
                }

                if (_lmax_dft_ecp > 1) {

                  
                    C[2][1][1] = CT * (2.0 * CP * CP - 1.0); // 2*SQ(pi/15) * (Z 2-2)'
                    C[2][1][2] = ST * (2.0 * CP * CP - 1.0); // 2*SQ(pi/15) * (Z 2-2)
                    C[2][1][3] = SQ3 * ST * ST * CP * SP; // 2*SQ(pi/5) * (Y 2-2)
                    C[2][1][4] = 2.0 * CT * ST * CP * SP; // 2*SQ(pi/15) * (Y 2-2)'
                    C[2][1][5] = (1.0 + CT * CT) * CP * SP; // 2*SQ(pi/15) * ( (Y 2-2)'' + 3*(Y 2-2) )

                    C[2][2][1] = -ST * CP; // 2*SQ(pi/15) * (Z 2-1)'
                    C[2][2][2] = CT * CP; // 2*SQ(pi/15) * (Z 2-1)
                    C[2][2][3] = SQ3 * CT * ST * SP; // 2*SQ(pi/5) * (Y 2-1)
                    C[2][2][4] = (2.0 * CT * CT - 1.0) * SP; // 2*SQ(pi/15) * (Y 2-1)'
                    C[2][2][5] = -CT * ST * SP; // 2*SQ(pi/15) * ( (Y 2-1)'' + 3*(Y 2-1) )

                    C[2][3][1] = 0.0;
                    C[2][3][2] = 0.0;
                    C[2][3][3] = 1.5 * CT * CT - 0.5; // 2*SQ(pi/5) * (Y 20)
                    C[2][3][4] = -SQ3 * CT * ST; // 2*SQ(pi/15) * (Y 20)'
                    C[2][3][5] = .5 * SQ3 * ST * ST; // 2*SQ(pi/15) * ( (Y 20)'' + 3*(Y 20) )         .5 * SQ3 * (1.0 - CT * CT)

                    C[2][4][1] = ST * SP; // 2*SQ(pi/15) * (Z 21)'
                    C[2][4][2] = -CT * SP; // 2*SQ(pi/15) * (Z 21)
                    C[2][4][3] = SQ3 * CT * ST * CP; // 2*SQ(pi/5) * (Y 21)
                    C[2][4][4] = (2.0 * CT * CT - 1.0) * CP; // 2*SQ(pi/15) * (Y 21)'
                    C[2][4][5] = -CT * ST * CP; // 2*SQ(pi/15) * ( (Y 21)'' + 3*(Y 21) )

                    C[2][5][1] = -2.0 * CT * CP * SP; // 2*SQ(pi/15) * (Z 22)'
                    C[2][5][2] = -2.0 * ST * CP * SP; // 2*SQ(pi/15) * (Z 22)
                    C[2][5][3] = SQ3 * ST * ST * (CP * CP - .5); // 2*SQ(pi/5) * (Y 22)         0.5 * SQ3 * (CT * CT * (1.0 - 2.0 * CP * CP) + 2.0 * CP * CP - 1.0)
                    C[2][5][4] = CT * ST * (2.0 * CP * CP - 1.); // 2*SQ(pi/15) * (Y 22)'
                    C[2][5][5] = (1. + CT * CT) * (CP * CP - .5); // 2*SQ(pi/15) * ( (Y 22)'' + 3*(Y 22) )       0.5 * (CP * CP * (2.0 + 2.0 * CT * CT) - CT * CT - 1.0)
                }

                if (_lmax_dft_ecp > 2) {


                    double f_phi = (4. * CP * CP - 1.) * SP;
                    double df_dphi = 3. * (1. - 4. * SP * SP) * CP;
                    C[3][0][0] = (SQ2 / (3. * SQ5))*(.5 * SQ5 / SQ2) * (.5 + 1.5 * CT * CT) * df_dphi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Z 3-3) + (Z 3-3)'' )
                    C[3][0][1] = (1. / (SQ3 * SQ5))*(.5 * SQ5 / SQ2) * 2. * ST * CT * df_dphi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * (Z 3-3)'
                    C[3][0][2] = (1. / (SQ2 * SQ3))*(.5 * SQ5 / SQ2) * ST * ST * df_dphi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Z 3-3)
                    C[3][0][3] = (.5 * SQ5 / SQ2) * ST * ST * ST * f_phi; // 2*SQ(pi/7) * (Y 3-3)
                    C[3][0][4] = (1. / (SQ2 * SQ3))*(.5 * SQ5 / SQ2) * 3. * ST * ST * CT * f_phi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 3-3)'
                    C[3][0][5] = (1. / (SQ3 * SQ5))*(.5 * SQ5 / SQ2) * 3. * (1. + CT * CT) * ST * f_phi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 3-3)'' + 6*(Y 3-3) )
                    C[3][0][6] = (SQ2 / (3. * SQ5))*(.5 * SQ5 / SQ2) * 1.5 * (3. + CT * CT) * CT * f_phi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 3-3)' + ( (Y 3-3)'' + 6*(Y 3-3) )' 

                    f_phi = CP * SP;
                    df_dphi = 1. - 2. * SP * SP;
                    C[3][1][0] = (SQ2 / (3. * SQ5))*(SQ3 * SQ5) * (-1.5) * ST * CT * df_dphi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Z 3-2) + (Z 3-2)'' )
                    C[3][1][1] = (1. / (SQ3 * SQ5))*(SQ3 * SQ5) * (1. - 2. * ST * ST) * df_dphi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * (Z 3-2)'
                    C[3][1][2] = (1. / (SQ2 * SQ3))*(SQ3 * SQ5) * ST * CT * df_dphi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Z 3-2)
                    C[3][1][3] = (SQ3 * SQ5) * ST * ST * CT * f_phi; // 2*SQ(pi/7) * (Y 3-2)
                    C[3][1][4] = (1. / (SQ2 * SQ3))*(SQ3 * SQ5) * (3. * CT * CT - 1.) * ST * f_phi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 3-2)'
                    C[3][1][5] = (1. / (SQ3 * SQ5))*(SQ3 * SQ5) * (3. * CT * CT - 1.) * CT * f_phi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 3-2)'' + 6*(Y 3-2) )
                    C[3][1][6] = (SQ2 / (3. * SQ5))*(SQ3 * SQ5) * (-1.5) * (1. + CT * CT) * ST * f_phi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 3-2)' + ( (Y 3-2)'' + 6*(Y 3-2) )'

                    f_phi = SP;
                    df_dphi = CP;
                    C[3][2][0] = (SQ2 / (3. * SQ5))*(.5 * SQ3 / SQ2) * 7.5 * ST * ST * df_dphi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Z 3-1) + (Z 3-1)'' )
                    C[3][2][1] = (1. / (SQ3 * SQ5))*(.5 * SQ3 / SQ2) * (-10.) * ST * CT * df_dphi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * (Z 3-1)'
                    C[3][2][2] = (1. / (SQ2 * SQ3))*(.5 * SQ3 / SQ2) * (5. * CT * CT - 1.) * df_dphi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Z 3-1)
                    C[3][2][3] = (.5 * SQ3 / SQ2) * (5. * CT * CT - 1.) * ST * f_phi; // 2*SQ(pi/7) * (Y 3-1)
                    C[3][2][4] = (1. / (SQ2 * SQ3))*(.5 * SQ3 / SQ2) * (4. - 15. * ST * ST) * CT * f_phi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 3-1)'
                    C[3][2][5] = (1. / (SQ3 * SQ5))*(.5 * SQ3 / SQ2) * (5. - 15. * CT * CT) * ST * f_phi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 3-1)'' + 6*(Y 3-1) )
                    C[3][2][6] = (SQ2 / (3. * SQ5))*(.5 * SQ3 / SQ2) * 7.5 * ST * ST * CT * f_phi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 3-1)' + ( (Y 3-1)'' + 6*(Y 3-1) )'

                    C[3][3][0] = 0.;
                    C[3][3][1] = 0.;
                    C[3][3][2] = 0.;
                    C[3][3][3] = (2.5 * CT * CT - 1.5) * CT; // 2*SQ(pi/7) * (Y 30)
                    C[3][3][4] = (1. / (SQ2 * SQ3))* (1.5 - 7.5 * CT * CT) * ST; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 30)'
                    C[3][3][5] = (1. / (SQ3 * SQ5))* 7.5 * ST * ST * CT; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 30)'' + 6*(Y 30) )
                    C[3][3][6] = (SQ2 / (3. * SQ5))* (-3.75 * ST * ST * ST); // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 30)' + ( (Y 30)'' + 6*(Y 30) )' )

                    f_phi = CP;
                    df_dphi = -SP;
                    C[3][4][0] = (SQ2 / (3. * SQ5))*(.5 * SQ3 / SQ2) * 7.5 * ST * ST * df_dphi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Z 31) + (Z 31)'' )
                    C[3][4][1] = (1. / (SQ3 * SQ5))*(.5 * SQ3 / SQ2) * (-10.) * ST * CT * df_dphi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * (Z 31)'
                    C[3][4][2] = (1. / (SQ2 * SQ3))*(.5 * SQ3 / SQ2) * (5. * CT * CT - 1.) * df_dphi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Z 31)
                    C[3][4][3] = (.5 * SQ3 / SQ2) * (5. * CT * CT - 1.) * ST * f_phi; // 2*SQ(pi/7) * (Y 31)
                    C[3][4][4] = (1. / (SQ2 * SQ3))*(.5 * SQ3 / SQ2) * (4. - 15. * ST * ST) * CT * f_phi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 31)'
                    C[3][4][5] = (1. / (SQ3 * SQ5))*(.5 * SQ3 / SQ2) * (5. - 15. * CT * CT) * ST * f_phi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 31)'' + 6*(Y 31) )
                    C[3][4][6] = (SQ2 / (3. * SQ5))*(.5 * SQ3 / SQ2) * 7.5 * ST * ST * CT * f_phi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 31)' + ( (Y 31)'' + 6*(Y 31) )'

                    f_phi = 1. - 2. * SP * SP;
                    df_dphi = -4. * CP * SP;
                    C[3][5][0] = (SQ2 / (3. * SQ5))*(.5 * SQ3 * SQ5) * (-1.5) * ST * CT * df_dphi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Z 32) + (Z 32)'' )
                    C[3][5][1] = (1. / (SQ3 * SQ5))*(.5 * SQ3 * SQ5) * (1. - 2. * ST * ST) * df_dphi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * (Z 32)'
                    C[3][5][2] = (1. / (SQ2 * SQ3))*(.5 * SQ3 * SQ5) * ST * CT * df_dphi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Z 32)
                    C[3][5][3] = (.5 * SQ3 * SQ5) * ST * ST * CT * f_phi; // 2*SQ(pi/7) * (Y 32)
                    C[3][5][4] = (1. / (SQ2 * SQ3))*(.5 * SQ3 * SQ5) * (3. * CT * CT - 1.) * ST * f_phi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 32)'
                    C[3][5][5] = (1. / (SQ3 * SQ5))*(.5 * SQ3 * SQ5) * (3. * CT * CT - 1.) * CT * f_phi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 32)'' + 6*(Y 32) )
                    C[3][5][6] = (SQ2 / (3. * SQ5))*(.5 * SQ3 * SQ5) * (-1.5) * (1. + CT * CT) * ST * f_phi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 32)' + ( (Y 32)'' + 6*(Y 32) )'

                    f_phi = (1. - 4. * SP * SP) * CP;
                    df_dphi = 3. * (1. - 4. * CP * CP) * SP;
                    C[3][6][0] = (SQ2 / (3. * SQ5))*(.5 * SQ5 / SQ2) * (.5 + 1.5 * CT * CT) * df_dphi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Z 33) + (Z 33)'' )
                    C[3][6][1] = (1. / (SQ3 * SQ5))*(.5 * SQ5 / SQ2) * 2. * ST * CT * df_dphi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * (Z 33)'
                    C[3][6][2] = (1. / (SQ2 * SQ3))*(.5 * SQ5 / SQ2) * ST * ST * df_dphi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Z 33)
                    C[3][6][3] = (.5 * SQ5 / SQ2) * ST * ST * ST * f_phi; // 2*SQ(pi/7) * (Y 33)
                    C[3][6][4] = (1. / (SQ2 * SQ3))*(.5 * SQ5 / SQ2) * 3. * ST * ST * CT * f_phi; // (1./(SQ2*SQ3)) * 2*SQ(pi/7) * (Y 33)'
                    C[3][6][5] = (1. / (SQ3 * SQ5))*(.5 * SQ5 / SQ2) * 3. * (1. + CT * CT) * ST * f_phi; // (1./(SQ3*SQ5)) * 2*SQ(pi/7) * ( (Y 33)'' + 6*(Y 33) )
                    C[3][6][6] = (SQ2 / (3. * SQ5))*(.5 * SQ5 / SQ2) * 1.5 * (3. + CT * CT) * CT * f_phi; // (SQ2/(3.*SQ5)) * 2*SQ(pi/7) * ( (5/2)*(Y 33)' + ( (Y 33)'' + 6*(Y 33) )'
                }


            } else {


                for (index L = 0; L <= _lmax_dft_ecp; L++) {
                    for (index M = 3 - L; M <= 3 + L; M++) {

                        C[L][M][M] = 1.;

                    }
                }
            }

            for (index I = 0; I < _nsph; I++) {
                for (index L = 0; L <= _lmax_dft; L++) {
                    for (index M = 3 - L; M <= 3 + L; M++) {
                        BLC[I][L][M] = 0.0;
                        for (index M1 = 3 - L; M1 <= 3 + L; M1++) {
                            BLC[I][L][M] += BLM[I][L][M1] * C[L][M1][M];
                        }
                    }
                }

            }

            return;
        } // getBLMCOF

        ub::vector<double> AOECP::CalcInt_r_exp(int nmax, double decay) {

            ub::vector<double> result = ub::zero_vector<double>(nmax + 1);
            const double PI = boost::math::constants::pi<double>();
            double factor = .5 / decay;

            result[0] = .5 * sqrt(PI / decay);
            result[1] = factor;
            for (int i = 2; i <= nmax; i++) {
                result[i] = (i - 1) * factor * result[i - 2];
            }

            return result;
        }    
    
    
    
}}

