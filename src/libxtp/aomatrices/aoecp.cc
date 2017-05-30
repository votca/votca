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
//#include <boost/timer/timer.hpp>


using namespace votca::tools;



namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
    

    void AOECP::FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col, AOBasis* ecp) {

            // get shell positions
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
                std::vector<double> _contractions_row_full(9);
                _contractions_row_full[0] = _contractions_row[0];
                _contractions_row_full[1] = _contractions_row[1];
                _contractions_row_full[2] = _contractions_row[1];
                _contractions_row_full[3] = _contractions_row[1];

                _contractions_row_full[4] = _contractions_row[2];
                _contractions_row_full[5] = _contractions_row[2];
                _contractions_row_full[6] = _contractions_row[2];
                _contractions_row_full[7] = _contractions_row[2];
                _contractions_row_full[8] = _contractions_row[2];


                for (AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc) {
                    //get decay constant
                    const double _decay_col = (*itc)->getDecay();
                    // if (_decay_col > 0.16) continue;
                    const double _fak = 0.5 / (_decay_row + _decay_col);
                    const double _fak2 = 2.0 * _fak;


                    double _exparg = _fak2 * _decay_row * _decay_col *_distsq;

                    // check if distance between postions is big, then skip step   

                    if (_exparg > 30.0) {
                        continue;
                    }

                    const std::vector<double>& _contractions_col = (*itc)->getContraction();
                    // shitty magic
                    std::vector<double> _contractions_col_full(9);
                    _contractions_col_full[0] = _contractions_col[0];
                    _contractions_col_full[1] = _contractions_col[1];
                    _contractions_col_full[2] = _contractions_col[1];
                    _contractions_col_full[3] = _contractions_col[1];

                    _contractions_col_full[4] = _contractions_col[2];
                    _contractions_col_full[5] = _contractions_col[2];
                    _contractions_col_full[6] = _contractions_col[2];
                    _contractions_col_full[7] = _contractions_col[2];
                    _contractions_col_full[8] = _contractions_col[2];
                    // for each atom and its pseudopotential, get a matrix
                    int _atomidx = 0;
                    
                    ub::matrix<int> _power_matrix = ub::zero_matrix<int>(5, 3); ///// // max 12 fit components, max non-local ECP l=0,1,2 ///////////////
                    ub::matrix<double> _decay_matrix = ub::zero_matrix<double>(5, 3); ///// // max 12 fit components, max non-local ECP l=0,1,2
                    ub::matrix<double> _coef_matrix = ub::zero_matrix<double>(5, 3); ///// // max 12 fit components, max non-local ECP l=0,1,2

                    AOBasis::AOShellIterator final_iter = ecp->lastShell();
                    --final_iter;
                    vec _ecp_eval_pos;
                    for (AOBasis::AOShellIterator _ecp = ecp->firstShell(); _ecp != ecp->lastShell(); _ecp++) {

                        const AOShell* _shell_ecp = ecp->getShell(_ecp);
                        const vec& _ecp_pos = _shell_ecp->getPos();

                        int this_atom = _shell_ecp->getIndex();

                        const int _ecp_l = _shell_ecp->getOffset(); //  angular momentum l is stored in offset for ECP

                        // only do the non-local parts
                        if (_ecp_l < _shell_ecp->getNumFunc()) {
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
                                    ub::matrix<double> VNL_ECP = calcVNLmatrix(_pos_row, _pos_col, _ecp_eval_pos, _decay_row, _decay_col,  _decay_matrix,_power_matrix, _coef_matrix); ////////////////

                                    // consider contractions
                                    // cut out block that is needed. sum
                                    //                                   cout << "_matrix.size1,2()   " << _matrix.size1() << "    " << _matrix.size2() << endl;
                                    for (unsigned i = 0; i < _matrix.size1(); i++) {
                                        for (unsigned j = 0; j < _matrix.size2(); j++) {
                                            _matrix(i, j) += VNL_ECP(i + _shell_row->getOffset(), j + _shell_col->getOffset()) * _contractions_row_full[i + _shell_row->getOffset()] * _contractions_col_full[j + _shell_col->getOffset()];
                                        }
                                    }


                                    // reset atom ECP containers
                                    _power_matrix = ub::zero_matrix<int>(5, 3); ///// // max 12 fit components, max non-local ECP l=0,1,2 /////////////////////
                                    _decay_matrix = ub::zero_matrix<double>(5, 3); ///// // max 12 fit components, max non-local ECP l=0,1,2
                                    _coef_matrix = ub::zero_matrix<double>(5, 3); ///// // max 12 fit components, max non-local ECP l=0,1,2
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
   
    
    
        ub::matrix<double> AOECP::calcVNLmatrix(const vec& posA, const vec& posB, const vec& posC, const double& alpha, const double& beta,const ub::matrix<double>& _gamma_ecp,const ub::matrix<int>& _power_ecp,const ub::matrix<double>& _pref_ecp) {


            const double pi = boost::math::constants::pi<double>();
            const double conv = 1.e-8;
            /* calculate the contribution of the nonlocal 
             *     ECP of atom at posC with 
             *       decay constants in _gamma_ecp
             *       coefficients in    _pref_ecp
             *       with angular momentum of max 2
             * 
             * to DFT basis shell pair 
             *    with decay alpha at posA
             *         decay beta  at posB
  
             */

            ub::matrix<double> matrix = ub::zero_matrix<double>(10, 10);

            const int nnonsep = _gamma_ecp.size1();
            const int nmax = 90;
            ub::matrix<double> XI(3, nmax);

            double f_even_r0 = .5 * sqrt(2.) * sqrt(.5 * pi);
            double f_even_r1 = .5;
            double f_even_r2 = .5 * f_even_r0;
            double f_odd_r0 = .5;
            double f_odd_r1 = .25 * sqrt(2.) * sqrt(.5 * pi);
            double f_odd_r2 = 1.;
            double DGAMAF_r0;
            double DGAMAF_r1;
            double DGAMAF_r2;

            for (int N = 0; N <= nmax - 1; N++) {

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

                for (int L = 0; L < 3; L++) {

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
            COEF.resize(extents[ range(0, 3) ][ range(0, 3) ][ range(1, 6)][range(0,42)]);

            // init it all to 0
            for ( index i1 = 0; i1 <=2; i1++ ){ //////////
                 for ( index i2 = 0; i2 <=2; i2++ ){ ///////////////
                     for ( index i3 = 1; i3 <=5; i3++ ){
                         for ( index i4 = 0; i4 <=41; i4++ ){
                           COEF[i1][i2][i3][i4] = 0.0;
                         }
                     }
                 }
             }
            for ( index i4 = 0; i4 <=41; i4++ ){  ////// for ( index i4 = 1; i4 <=42; i4++ )
            /********** ORIGINAL CKOEF SUBROUTINE *************************/
                int NU = i4 % 2; ///
                int NG = (i4 + 1) % 2;
                double FN1 = double(i4 + 1);
                double FN2 = double(i4 + 2);
                double FN3 = double(i4 + 3);
                double FN4 = double(i4 + 4);
                double FN5 = double(i4 + 5);
 
                COEF[0][0][3][i4] = NG/FN1;   /////////   M0(x)
                COEF[0][1][3][i4] = NU/FN2*sqrt(3.0);   ////////  SQ(3) * M1(x)
                COEF[0][2][3][i4] = NG/2.0*sqrt(5.0)*(3.0/FN3-1.0/FN1);   //////   SQ(5) * M2(x)

                COEF[1][0][3][i4] = COEF[0][1][3][i4];
                COEF[1][1][3][i4] = NG*3.0/FN3;                        ///////    M0(x) + 2 * M2(x)
                COEF[1][1][4][i4] = 3.0/2.0*NG*(1.0/FN1-1.0/FN3);     ////////    M0(x) - M2(x) 
                COEF[1][1][2][i4] = COEF[1][1][4][i4];
                COEF[1][2][3][i4] = sqrt(15.0)/2.0*NU*(3.0/FN4-1.0/FN2);     ///////   (2/5) * SQ(15) * ( M1(x) + (3/2) * M3(x) )
                COEF[1][2][4][i4] = sqrt(45.0)/2.0*NU*(1.0/FN2-1.0/FN4);     ///////   (SQ(45)/5) * ( M1(x) - M3(x) )
                COEF[1][2][2][i4] = COEF[1][2][4][i4];

                COEF[2][0][3][i4] = COEF[0][2][3][i4];
                COEF[2][1][3][i4] = COEF[1][2][3][i4];
                COEF[2][1][4][i4] = COEF[1][2][4][i4];
                COEF[2][1][2][i4] = COEF[1][2][4][i4];
                COEF[2][2][3][i4] = 5.0/4.0*NG*(9.0/FN5-6.0/FN3+1.0/FN1);   ///////  M0(x) + (10/7)*M2(x) + (18/7)*M4(x)
                COEF[2][2][4][i4] = NG*15.0/2.0*(1.0/FN3-1.0/FN5);          ///////  M0(x) + (5/7)*M2(x) - (12/7)*M4(x)    
                COEF[2][2][5][i4] = 15.0/8.0*NG*(1.0/FN1-2.0/FN3+1.0/FN5);  ///////  M0(x) - (10/7)*M2(x) + (3/7)*M4(x) 
                COEF[2][2][1][i4] = COEF[2][2][5][i4];
                COEF[2][2][2][i4] = COEF[2][2][4][i4];
             
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

            if (AVSSQ <= 1.0e-1) {

                NMAX1 = 1;

            } else {

                double AMAX = 0.0;
                double fak = 2.0 * alpha*AVSSQ;
                double Pow = 1;
                double factorialNN = 1;
                for (int N = 1; N <= 43; N++) {

                    int NN = N - 1;


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

                    if (NMAX1 == 0 && AMAX <= conv) NMAX1 = N;
                    if (NMAX1 != 0 && AMAX > conv) NMAX1 = 0;

                }
                if (NMAX1 == 0 && AMAX > conv) NMAX1 = 42;
            }

            // same story for B
            if (BVSSQ <= 1.0e-1) {

                NMAX2 = 1;

            } else {

                double BMAX = 0.0;
                double fak = 2.0 * beta*BVSSQ;
                double Pow = 1;
                double factorialNN = 1;
                for (int N = 1; N <= 42; N++) {

                    int NN = N - 1;

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

                    if (NMAX2 == 0 && BMAX <= conv) NMAX2 = N;
                    if (NMAX2 != 0 && BMAX > conv) NMAX2 = 0;

                }

                if (NMAX2 == 0 && BMAX > conv) NMAX2 = 42;
            }

            // something
            int INULL = 1;
            if (AVSSQ <= 1.e-1) INULL = 3;
            if (BVSSQ <= 1.e-1) INULL = 4;
            if (AVSSQ <= 1.e-1 && BVSSQ <= 1.e-1) INULL = 2;

            type_3D BLMA;
            type_3D CA;
            getBLMCOF(AVS, BLMA, CA);


            type_3D BLMB;
            type_3D CB;
            getBLMCOF(BVS, BLMB, CB);
  
            typedef boost::multi_array_types::extent_range range;
            typedef type_3D::index index;
            type_3D::extent_gen extents3D;

            type_3D CC;
            CC.resize(extents3D[ range(0,3)][range(1,6)][range(1,6)]); ///////
            for ( index L = 0; L<=2; L++){ ////
                  for ( index M1 = 1; M1<=5; M1++){
                      for ( index M2 = 1; M2<=5; M2++){

                           CC[L][M1][M2]=0.0;
                           for ( index M = 1; M<=5; M++){

                             CC[L][M1][M2] += CA[L][M][M1]*CB[L][M][M2]; /////

                           }
                      }
                  }
            }

            typedef boost::multi_array<double, 5> type_5D;
            type_5D::extent_gen extents5D;
            type_5D SUMCI;
            SUMCI.resize(extents5D[range(0,3)][range(0,3)][ range(0,3)][range(1,6)][range(1,6)]); ////
            type_3D SUMCI3;
            SUMCI3.resize(extents3D[range(0,3)][range(0,3)][range(1,6)]); ////


            switch (INULL) {

          case 1:
          {

              for ( index L = 0; L <= 2; L++  ){ /////
                  for ( index L1 = 0; L1 <= 2; L1++  ){ /////                
                      for ( index L2 = 0; L2 <= 2; L2++  ){ /////          
                          for ( index M1 = 3-L; M1 <= 3+L; M1++ ){
                              for ( index M2 = 3-L; M2 <= 3+L; M2++ ){

                                  SUMCI[L][L1][L2][M1][M2]  = 0.0; ////////////

                                  double fak1=2.0*alpha*AVSSQ;
                                  double pow1=1;
                                  double factorialN=1;

                                  for ( int N = 0; N <= NMAX1-1; N++ ) { ///

                                      if (N!=0) {
                                          pow1=pow1*fak1;
                                          factorialN=factorialN*N;
                                      }
                                            
                                      double VAR1 = COEF[L][L1][M1][N]*pow1/factorialN; ///////////
                                      double VAR2 = 0.0;
                                      double fak2=2.0*beta*BVSSQ;
                                      double pow2=1;
                                      double factorialNN=1;

                                      for ( int NN = 0; NN <= NMAX2-1; NN++ ) {

                                          if (NN!=0) {
                                              pow2=pow2*fak2;
                                              factorialNN=factorialNN*NN;
                                          }
                                          double XDUM = COEF[L][L2][M2][NN]*pow2/factorialNN;  //////
                                          VAR2  += XDUM*XI(L,N+NN+L1+L2); // L index of XI starts with 0 !! /////////

                                      }
      
                                      SUMCI[L][L1][L2][M1][M2]  += VAR1*VAR2; /////

                                  }

                              } // end M2
                          } // end M1
                      } // end L2
                  } // end L1 
              } // end L

              // now finally calculate matrix

              for (unsigned i = 0; i < matrix.size1(); i++) {   // matrix.size1() = 10
                  for (unsigned j = 0; j < matrix.size2(); j++) {   // matrix.size2() = 10

                      for (index L = 0; L <= 2; L++) { ////
                          for (index L1 = 0; L1 <= 2; L1++) { ////
                              for (index L2 = 0; L2 <= 2; L2++) { /////

                                  for (index M1 = 3-L; M1 <= 3+L; M1++) {
                                      for (index M2 = 3-L; M2 <= 3+L; M2++) {

                                          matrix(i,j) += BLMA[i][L1][M1] * BLMB[j][L2][M2] * SUMCI[L][L1][L2][M1][M2]*CC[L][M1][M2]; ////////////////

                                      }
                                  }
                              }
                          }
                      }
                  }
              }

              break;
          }


          case 2:  //  AVSSQ<=1.e-1 && BVSSQ <= 1.e-1
          {

              for (unsigned i = 0; i < matrix.size1(); i++) {
                  for (unsigned j = 0; j < matrix.size2(); j++) {

                      for (index L = 0; L <= 2; L++) {
                          double XI_L = XI(L,L+L);
                          for (index M = 3-L; M <= 3+L; M++) {

                              matrix(i,j) += BLMA[i][L][M] * BLMB[j][L][M] * XI_L;

                          }
                      }

                  }
              }


              break;
          }


          case 3:  //  AVSSQ <= 1.e-1
          {

              for ( index L = 0; L <= 2; L++  ) {
                  for ( index L2 = 0; L2 <= 2; L2++  ) {
                      for ( index M2 = 3-L; M2 <= 3+L; M2++ ) {

                          double VAR2 = 0.0;
                          double fak=2.0 * beta*BVSSQ;
                          double pow=1;
                          double factorialNN=1;

                          for (int NN = 0; NN <= NMAX2-1; NN++) {

                              if(NN!=0) {
                                  pow=pow*fak;
                                  factorialNN=factorialNN*NN;
                              }
                            
                              double XDUM = COEF[L][L2][M2][NN] * pow / factorialNN;
                              VAR2 += XDUM * XI(L, NN + L + L2);

                          }

                          SUMCI3[L][L2][M2] = VAR2;

                      } // end M2
                  } // end L2
              } // end L


              for (unsigned i = 0; i < matrix.size1(); i++) {
                  for (unsigned j = 0; j < matrix.size2(); j++) {

                      for (index L = 0; L <= 2; L++) {
                          for (index L2 = 0; L2 <= 2; L2++) {

                              for (index M1 = 3-L; M1 <= 3+L; M1++) {
                                  for (index M2 = 3-L; M2 <= 3+L; M2++) {

                                      matrix(i,j) += BLMA[i][L][M1] * BLMB[j][L2][M2] * SUMCI3[L][L2][M2] * CC[L][M1][M2];

                                  }
                              }
                          }
                      }
                  }
              }

              break;
          }


          case 4:  //  BVSSQ <= 1.e-1
          {

              for ( index L = 0; L <= 2; L++  ) {
                  for ( index L1 = 0; L1 <= 2; L1++  ) {
                      for ( index M1 = 3-L; M1 <= 3+L; M1++ ) {

                          double VAR1 = 0.0;
                          double fak=2.0 * alpha*AVSSQ;
                          double pow=1;
                          double factorialN=1;

                          for (int N = 0; N <= NMAX1-1; N++) {

                              if(N!=0) {
                                  pow=pow*fak;
                                  factorialN=factorialN*N;
                              }

                              double XDUM = COEF[L][L1][M1][N] * pow / factorialN;
                              VAR1 += XDUM * XI(L, N + L1 + L);

                          }

                          SUMCI3[L][L1][M1] = VAR1;

                      } // end M1
                  } // end L1
              } // end L


              for (unsigned i = 0; i < matrix.size1(); i++) {
                  for (unsigned j = 0; j < matrix.size2(); j++) {

                      for (index L = 0; L <= 2; L++) {
                          for (index L1 = 0; L1 <= 2; L1++) {

                              for (index M1 = 3-L; M1 <= 3+L; M1++) {
                                  for (index M2 = 3-L; M2 <= 3+L; M2++) {

                                      matrix(i,j) += BLMA[i][L1][M1] * BLMB[j][L][M2] * SUMCI3[L][L1][M1] * CC[L][M1][M2];

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
            std::vector<double> NormA(10);
            std::vector<double> NormB(10);
            getNorms(NormA, alpha);
            getNorms(NormB, beta);


            for (unsigned i = 0; i < matrix.size1(); i++) {
                for (unsigned j = 0; j < matrix.size2(); j++) {

                    matrix(i,j) = matrix(i,j) * GAUSS * NormA[i] * NormB[j];

                }
        }
     
      
        return matrix;
        }
    
    
    void AOECP::getNorms(std::vector<double>& Norms, const double decay){

            const double PI = boost::math::constants::pi<double>();
            Norms[0] = pow(2.0 * decay / PI, 0.75);   ///  Y 00
            Norms[1] = 2.0 * sqrt(decay) * Norms[0];  ///  Y 10
            Norms[2] = Norms[1];                      ///  Y 1 -1
            Norms[3] = Norms[1];                      ///  Y 11
            Norms[4] = 2.00 * decay * Norms[0] / sqrt(3.0); ///  Y 20
            Norms[5] = 4.00 * decay * Norms[0];             ///  Y 2 -1
            Norms[6] = Norms[5];                            ///  Y 21
            Norms[7] = Norms[5];                            ///  Y 2 -2
            Norms[8] = 2.00 * decay * Norms[0];             ///  Y 22
            Norms[9] = Norms[5] / sqrt(15.0);
        
            return;
    }

        void AOECP::getBLMCOF(const vec& pos, type_3D& BLC, type_3D& C) {

            typedef boost::multi_array_types::extent_range range;
            typedef type_3D::index index;
            type_3D::extent_gen extents;

            BLC.resize(extents[ range(0, 10) ][ range(0, 3) ][ range(1, 6)]); /////
            C.resize(extents[ range(0, 3) ][ range(1, 6) ][ range(1, 6)]); /////////

            type_3D BLM;
            BLM.resize(extents[ range(0, 10) ][ range(0, 3) ][ range(1, 6)]); ////

            const double PI = boost::math::constants::pi<double>();
            double SPI = sqrt(PI);
            double XS = 2.0 * SPI; /////////////          2 * SQ(pi)         
            double XP = XS / sqrt(3.0); /////////////      2 * SQ(pi/3)
            double XD1 = XP / sqrt(5.0); /////////////  2 * SQ(pi/15)
            double XD2 = 4.0 * SPI / sqrt(5.0); ////////  4 * SQ(pi/5)         Y 20
            double XD3 = XD2 / sqrt(3.0); ////////      4 * SQ(pi/15)        Y 22

            for (index I = 0; I <= 9; I++) { /////
                for (index L = 0; L <= 2; L++) { ///
                    for (index M = 1; M <= 5; M++) {

                        BLM[I][L][M] = 0.0;

                    }
                }
            }

            const tools::vec& BVS = pos;
            double XY = BVS[0] * BVS[0] + BVS[1] * BVS[1];
            double XYZ = XY + BVS[2] * BVS[2];
            double SXY = sqrt(XY); //// SXY = r * sin(theta)
            double SXYZ = sqrt(XYZ); //// SXYZ = r
            double CP = 1.0;
            double SP = 0.0;

            BLM[0][0][3] = XS; ///  Y 00

            BLM[1][0][3] = -BVS[2] * XS; ///  Y 10
            BLM[1][1][3] = XP;

            BLM[2][0][3] = -BVS[1] * XS; ///  Y 1 -1
            BLM[2][1][2] = XP;

            BLM[3][0][3] = -BVS[0] * XS; ///  Y 11
            BLM[3][1][4] = XP;

            BLM[4][0][3] = (2.0 * BVS[2] * BVS[2] - XY) * XS; ///  Y 20
            BLM[4][1][4] = 2.0 * BVS[0] * XP;
            BLM[4][1][2] = 2.0 * BVS[1] * XP;
            BLM[4][1][3] = -4.0 * BVS[2] * XP;
            BLM[4][2][3] = XD2;

            BLM[5][0][3] = BVS[1] * BVS[2] * XS; ///  Y 2 -1
            BLM[5][1][2] = -BVS[2] * XP;
            BLM[5][1][3] = -BVS[1] * XP;
            BLM[5][2][2] = XD1;

            BLM[6][0][3] = BVS[0] * BVS[2] * XS; ///  Y 21
            BLM[6][1][4] = -BVS[2] * XP;
            BLM[6][1][3] = -BVS[0] * XP;
            BLM[6][2][4] = XD1;

            BLM[7][0][3] = BVS[0] * BVS[1] * XS; ///  Y 2 -2
            BLM[7][1][4] = -BVS[1] * XP;
            BLM[7][1][2] = -BVS[0] * XP;
            BLM[7][2][1] = XD1;

            BLM[8][0][3] = (BVS[0] * BVS[0] - BVS[1] * BVS[1]) * XS; ///  Y 22
            BLM[8][1][4] = -2.0 * BVS[0] * XP;
            BLM[8][1][2] = 2.0 * BVS[1] * XP;
            BLM[8][2][5] = XD3;

            BLM[9][0][3] = (XYZ) * XS;
            BLM[9][1][4] = -2.0 * BVS[0] * XP;
            BLM[9][1][2] = -2.0 * BVS[1] * XP;
            BLM[9][1][3] = -2.0 * BVS[2] * XP;


            for (index L = 0; L <= 2; L++) { ///////
                for (index M = 1; M <= 5; M++) {
                    for (index MM = 1; MM <= 5; MM++) {

                        C[L][M][MM] = 0.0;

                    }
                }
            }


            if (SXY > 1.e-4) {
                CP = BVS[0] / SXY; //// CP = cos(phi)
                SP = BVS[1] / SXY; //// SP = sin(phi)
            }

            if (SXYZ > 1.e-4) {

                double CT = BVS[2] / SXYZ; /// CT = cos(theta)
                double ST = SXY / SXYZ; /// ST = sin(theta)

                C[0][3][3] = 1.0; //                2*SQ(pi) * Y 00   ############################

                C[1][2][2] = CP;
                C[1][2][3] = ST *SP; //             2*SQ(pi/3) * Y 1 -1  ##########################
                C[1][2][4] = CT * SP;

                C[1][3][2] = 0.0;
                C[1][3][3] = CT; //               2*SQ(pi/3) * Y 10  #############################
                C[1][3][4] = -ST;

                C[1][4][2] = -SP;
                C[1][4][3] = CP * ST; //            2*SQ(pi/3) * Y 11  ###############################
                C[1][4][4] = CT * CP;

                C[2][1][1] = CT * (2.0 * CP * CP - 1.0);
                C[2][1][2] = ST * (2.0 * CP * CP - 1.0);
                double SQ3 = sqrt(3.0);
                C[2][1][3] = SQ3 * CP * SP * ST * ST; //    2*SQ(pi/5) * Y 2 -2  #################################
                C[2][1][4] = 2.0 * CT * CP * SP * ST;
                C[2][1][5] = CP * SP * (1.0 + CT * CT);

                C[2][2][1] = -CP*ST;
                C[2][2][2] = CT*CP;
                C[2][2][3] = SQ3 * CT * ST* SP; //          2*SQ(pi/5) * Y 2 -1  ################################
                C[2][2][4] = SP * (2.0 * CT * CT - 1.0);
                C[2][2][5] = -CT * ST*SP;

                C[2][3][1] = 0.0;
                C[2][3][2] = 0.0;
                C[2][3][3] = 1.5 * CT * CT - 0.5; //          2*SQ(pi/5) * Y 20  ################################
                C[2][3][4] = -SQ3 * CT*ST;
                C[2][3][5] = .5 * SQ3 * ST * ST; /// .5 * SQ3 * (1.0 - CT * CT)

                C[2][4][1] = ST * SP;
                C[2][4][2] = -CT * SP;
                C[2][4][3] = SQ3 * CT * CP * ST; //           2*SQ(pi/5) * Y 21  #######################
                C[2][4][4] = CP * (2.0 * CT * CT - 1.0);
                C[2][4][5] = -CT * CP * ST;

                C[2][5][1] = -2.0 * CT * CP * SP;
                C[2][5][2] = -2.0 * CP * ST * SP;
                C[2][5][3] = SQ3 * ST * ST * (CP * CP - .5); //////   2*SQ(pi/5) * Y 22 ##########   0.5 * SQ3 * (CT * CT * (1.0 - 2.0 * CP * CP) + 2.0 * CP * CP - 1.0)
                C[2][5][4] = CT * ST * (2.0 * CP * CP - 1.);
                C[2][5][5] = (CP * CP - .5) * (1. + CT * CT); /// 0.5 * (CP * CP * (2.0 + 2.0 * CT * CT) - CT * CT - 1.0)
            } else {
                C[0][3][3] = 1.0; //// CT = CP = 1,        ST = SP = 0
                C[1][2][2] = 1.0;
                C[1][3][3] = 1.0;
                C[1][4][4] = 1.0;


                for (index M = 1; M <= 5; M++) {
                    C[2][M][M] = 1.0;
                }
            }
            for (index I = 0; I <= 9; I++) { ////
                for (index L = 0; L <= 2; L++) { ////
                    for (index M = 1; M <= 5; M++) {

                        BLC[I][L][M] = 0.0;
                        for (index M1 = 1; M1 <= 5; M1++) {

                            BLC[I][L][M] += BLM[I][L][M1] * C[L][M1][M]; ///

                        }

                    }
                }

            }

            return;
        } // getBLMCOF
    
    
    
}}

