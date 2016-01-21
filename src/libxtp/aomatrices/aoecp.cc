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
#include <boost/math/special_functions/factorials.hpp>
#include <boost/multi_array.hpp>
#include <votca/xtp/logger.h>
#include <votca/tools/linalg.h>
#include <votca/xtp/elements.h>
//#include <boost/timer/timer.hpp>

using namespace std;
using namespace votca::tools;



namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
    

    
    void AOECP::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col , AOBasis* ecp) {
        /*cout << "\nAO block: "<< endl;
        cout << "\t row: " << _shell_row->getType() << " at " << _shell_row->getPos() << endl;
        cout << "\t col: " << _shell_col->getType() << " at " << _shell_col->getPos() << endl;*/
        //const double pi = boost::math::constants::pi<double>();
       
        
        // cout << _gridpoint << endl;
        // shell info, only lmax tells how far to go
        //int _lmax_row = _shell_row->getLmax();
        //int _lmax_col = _shell_col->getLmax();

        // set size of internal block for recursion
        //int _nrows = this->getBlockSize( _lmax_row ); 
        //int _ncols = this->getBlockSize( _lmax_col ); 
        
        // get shell positions
        const vec& _pos_row = _shell_row->getPos();
        const vec& _pos_col = _shell_col->getPos();
        const vec  _diff    = _pos_row - _pos_col;
        // initialize some helper
        vector<double> PmA (3,0.0);
        vector<double> PmB (3,0.0);
        vector<double> PmC (3,0.0);
        double _distsq = (_diff.getX()*_diff.getX()) + (_diff.getY()*_diff.getY()) + (_diff.getZ()*_diff.getZ()); 
        
         typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
        // iterate over Gaussians in this _shell_row
        for ( GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
            // get decay constant
            const double& _decay_row = (*itr)->decay;
            
            //if ( _decay_row > 0.08 ) continue;
            
            std::vector<double> _contractions_row = (*itr)->contraction;
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
                

                for ( GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
                //get decay constant
                const double& _decay_col = (*itc)->decay;
        // if (_decay_col > 0.16) continue;
                const double _fak  = 0.5/(_decay_row + _decay_col);
                const double _fak2 = 2.0 * _fak;
                
                
                double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
        
       // check if distance between postions is big, then skip step   
       
                 if ( _exparg > 30.0 ) { continue; }
        
                
                
                
                
                std::vector<double> _contractions_col = (*itc)->contraction;
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
                ub::matrix<double> _decay_matrix = ub::zero_matrix<double>(1,3); // max 12 fit components, max non-local ECP l=0,1,2
                ub::matrix<double> _coef_matrix  = ub::zero_matrix<double>(1,3); // max 12 fit components, max non-local ECP l=0,1,2
                
                 vector< AOShell* >::iterator final_iter = ecp->lastShell();
                 --final_iter;
                vec _ecp_eval_pos;
                for (vector< AOShell* >::iterator _ecp = ecp->firstShell(); _ecp != ecp->lastShell() ; _ecp++ ) {
            
                   AOShell* _shell_ecp = ecp->getShell( _ecp );
                   const vec& _ecp_pos = _shell_ecp->getPos();
                   
                   int this_atom = _shell_ecp->getIndex();

                   const int _ecp_l = _shell_ecp->getOffset(); //  angular momentum l is stored in offset for ECP
                   
                   // only do the non-local parts
                        if (_ecp_l < _shell_ecp->getNumFunc()) {
                            int i_fit = -1;
                            for (GaussianIterator itecp = _shell_ecp->firstGaussian(); itecp != _shell_ecp->lastGaussian(); ++itecp) {
                                i_fit++;
                                
                                // get info for this angular momentum shell
                                const double& _decay_ecp = (*itecp)->decay;
                                const double& _contraction_ecp = (*itecp)->contraction[0];
                                //const int& _power_ecp = (*itecp)->power;


                                

                                // collect atom ECP
                                if ( this_atom == _atomidx ){
                                    _ecp_eval_pos = _ecp_pos;
                                    _decay_matrix(i_fit, _ecp_l ) = _decay_ecp;
                                    _coef_matrix(i_fit, _ecp_l )  = _contraction_ecp;
                                    
                                } 
                                
                                if ( (this_atom != _atomidx ) || ( _ecp == final_iter) ){
                                    
                                   // evaluate collected data, returns a (10x10) matrix of already normalized matrix elements 
                                   ub::matrix<double> VNL_ECP =  calcVNLmatrix(_pos_row,_pos_col,_ecp_eval_pos,_decay_row,_decay_col,_decay_matrix,_coef_matrix);
                                    
                                   
                                    // cout << _decay_row << " " << _decay_col << " " << VNL_ECP;
                                   
                                   
                                   // consider contractions
                                   // cut out block that is needed. sum
                                           for ( unsigned i = 0; i< _matrix.size1(); i++ ) {
                                               for (unsigned j = 0; j < _matrix.size2(); j++){
                                                _matrix(i,j) += VNL_ECP(i+_shell_row->getOffset(),j+_shell_col->getOffset()) * _contractions_row_full[i+_shell_row->getOffset()]* _contractions_col_full[j+_shell_col->getOffset()];
                                               }
                                              }
        
                                   
                                   
                                   
                                   // reset atom ECP containers
                                   _decay_matrix = ub::zero_matrix<double>(1,3); // max 12 fit components, max non-local ECP l=0,1,2
                                   _coef_matrix  = ub::zero_matrix<double>(1,3); // max 12 fit components, max non-local ECP l=0,1,2
                                   _atomidx++;
                                   i_fit = 0;
                                   //cout << "setting new matrix " << i_fit << " l " << _ecp_l << " alpha  " << _decay_ecp <<  " pref " << _contraction_ecp << endl;
                                    _decay_matrix(i_fit, _ecp_l ) = _decay_ecp;
                                    _coef_matrix(i_fit, _ecp_l )  = _contraction_ecp;
                                } // evaluate if new atom is found




                            } // all Gaussians in ecp_shell
                        } // only for non local parts
                   
                    } // all ecp_shells


                   // evaluate collected data 
                         //          ub::matrix<double> VNL_ECP =  calcVNLmatrix(_pos_row,_pos_col,_ecp_pos,_decay_row,_decay_col,_decay_matrix,_coef_matrix);

                //exit(0);




             
        
       // boost::timer::cpu_times t11 = cpu_t.elapsed();
        
        //cout << "Done with unnormalized matrix " << endl;
        
      
        //}
        //nuc.clear();
            }// _shell_col Gaussians
        }// _shell_row Gaussians
    }
   
    
    
    
    ub::matrix<double> AOECP::calcVNLmatrix(const vec& posA, const vec& posB, const vec& posC, const double& alpha, const double& beta, ub::matrix<double>& _gamma_ecp, ub::matrix<double>& _pref_ecp){
        
        
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
        
        
        
        
        ub::matrix<double> matrix = ub::zero_matrix<double>(10,10);


            const int nnonsep = _gamma_ecp.size1();
            const int nmax = 90;
            ub::matrix<double> XI(3, nmax);

            /****** ORIGINAL XINPN2 SUBROUTINE    ****************/
            for (int np = 1; np <= nmax; np++) {

                int N = np - 1;
                int NN = np + 2;
                double DGAMAF = 1.0;
                int NNG = NN % 2;

                if (NNG == 0) {
                    // NNG even
                    int NM = NN / 2 - 1;
                    for (int I = 1; I <= NM; I++) {
                        DGAMAF = DGAMAF * double(I);
                    }

                } else {
                    // NNG odd
                    int NF = (NN - 3) / 2;
                    if (NF > 0) {
                        for (int I = 1; I <= NF; I++) {
                            DGAMAF = DGAMAF * (1.5 + 1.0 * double(I - 1));
                        }


                    }
                    DGAMAF = DGAMAF * sqrt(0.25 * pi);
                }
                
                

                double DFAK = 0.5 * double(N + 3);
                //double DFAKP1 = DFAK + 1.0;

                 //cout << np << " " << DFAK << " " << DFAKP1 << " " <<  DGAMAF << endl;
                
                 
                 //cout << " XI pref_ecp " << _pref_ecp << endl;
                 //cout << " XI alpha_ecp " << _gamma_ecp << endl;
                
                for (int L = 0; L < 3; L++) {

                    
                    /***** ONLY r^2 powers in ECP ********************************/
                    XI(L, np-1) = 0.0;
                    for (int I = 0; I < nnonsep; I++) {
                        double DLI = (alpha + beta + _gamma_ecp(I, L)); // r^0 terms
                        //double DLJ = (beta + alpha + _gamma_ecp(I , L)); // r^2 terms (here!)
                        XI(L, np-1) += _pref_ecp(I, L) / pow(DLI, DFAK);// + _pref_ecp(I + nnonsep, L) * DFAK / pow(DLJ, DFAKP1);
                        //XI(L, np-1) +=  _pref_ecp(I , L) * DFAK / pow(DLJ, DFAKP1);
                        
                    }

                    XI(L,np-1) = XI(L,np-1)*DGAMAF*0.5;
                    
                }


            }
            
            //cout << "XI(0,0) " << XI(0,0) << " " << XI(1,4) << endl;
            //exit(0);
            
         /****** ORIGINAL CKO SUBROUTINE **********/   
         // get a multi dimensional array
         typedef boost::multi_array<double, 4> ma_type;
         typedef boost::multi_array_types::extent_range range;
         typedef ma_type::index index;
         ma_type::extent_gen extents;
         ma_type COEF;
         COEF.resize(extents[ range(1, 4) ][ range(1, 4) ][ range(1, 6)][range(1,43)]);
         
         for ( index i4 = 1; i4 <=42; i4++ ){

             // init it all to 0
             for ( index i1 = 1; i1 <=3; i1++ ){
                 for ( index i2 = 1; i2 <=3; i2++ ){
                     for ( index i3 = 1; i3 <=5; i3++ ){
                         
                         COEF[i1][i2][i3][i4] = 0.0;
                         
                         
                     }
                 }
                 
             } 
             
            /********** ORIGINAL CKOEF SUBROUTINE *************************/
                int N = i4 - 1;
                int NU = N % 2;
                int NN = i4;
                int NG = NN % 2;
                double FN1 = double(N + 1);
                double FN2 = double(N + 2);
                double FN3 = double(N + 3);
                double FN4 = double(N + 4);
                double FN5 = double(N + 5);

                
                COEF[1][1][3][i4] = NG/FN1;
                COEF[1][2][3][i4] = NU/FN2*sqrt(3.0);
                COEF[1][3][3][i4] = NG/2.0*sqrt(5.0)*(3.0/FN3-1.0/FN1);
                COEF[2][1][3][i4] = COEF[1][2][3][i4];
                COEF[2][2][3][i4] = NG*3.0/FN3;
                COEF[2][2][4][i4] = 3.0/2.0*NG*(1.0/FN1-1.0/FN3);
                COEF[2][2][2][i4] = COEF[2][2][4][i4];
                COEF[2][3][3][i4] = sqrt(15.0)/2.0*NU*(3.0/FN4-1.0/FN2);
                COEF[2][3][4][i4] = sqrt(45.0)/2.0*NU*(1.0/FN2-1.0/FN4);
                COEF[2][3][2][i4] = COEF[2][3][4][i4];
                COEF[3][1][3][i4] = COEF[1][3][3][i4];
                COEF[3][2][3][i4] = COEF[2][3][3][i4];
                COEF[3][2][4][i4] = COEF[2][3][4][i4];
                COEF[3][2][2][i4] = COEF[2][3][4][i4];
                COEF[3][3][3][i4] = 5.0/4.0*NG*(9.0/FN5-6.0/FN3+1.0/FN1);   
                COEF[3][3][4][i4] = NG*15.0/2.0*(1.0/FN3-1.0/FN5);        
                COEF[3][3][5][i4] = 15.0/8.0*NG*(1.0/FN1-2.0/FN3+1.0/FN5); 
                COEF[3][3][1][i4] = COEF[3][3][5][i4];
                COEF[3][3][2][i4] = COEF[3][3][4][i4];
             
         } // i4 loop (== CKO )
 
         // cout << COEF[1][1][1][1] << "  " << COEF[1][1][3][1] << " " << COEF[3][3][2][42] << endl;
         
         
         /**** PREPARATIONS DONE, NOW START ******/
         vec AVS = posA - posC;
         vec BVS = posB - posC;
         
         double AVS2 = AVS.getX()*AVS.getX() + AVS.getY()*AVS.getY() + AVS.getZ()*AVS.getZ();
         double BVS2 = BVS.getX()*BVS.getX() + BVS.getY()*BVS.getY() + BVS.getZ()*BVS.getZ();
         
         double AVSSQ = sqrt(AVS2);
         double BVSSQ = sqrt(BVS2);
         double GAUSS = exp(-alpha*AVS2-beta*BVS2);
         
         
         //cout << AVS2 << " " << BVS2 << " " << GAUSS << endl;
         //exit(0);
         
         
         // some limit determinations
         double G1 = exp(-alpha*AVS2);
         double G2 = exp(-beta*BVS2);
         
         int NMAX1=0;
         int NMAX2=0;
         
         if ( AVSSQ <= 1.0e-1 ){
             
             NMAX1=1;
             
         } else {
             
             double AMAX = 0.0;  
             for ( int N=1; N<=43; N++ ){
                 
                 int NN = N-1;
                 double factorialNN = boost::math::factorial<double>(double(NN));
                 double AF=pow(2.0*alpha*AVSSQ,NN/factorialNN)*G1;
                 double AF1 = std::abs(AF * XI(0,NN + 4));
                 double AF2 = std::abs(AF * XI(1,NN + 4));
                 double AF3 = std::abs(AF * XI(2,NN + 4));
                 AMAX = std::max(AF1,AF2);
                 AMAX = std::max(AMAX,AF3);
                 
                 if ( NMAX1 == 0 && AMAX <= conv ) NMAX1 = N;
                 if ( NMAX1 != 0 && AMAX >  conv ) NMAX1 = 0;
             
         }
             
             if (NMAX1==0 && AMAX > conv ) NMAX1 = 42;
             
         }
         
         // same story for B
         if ( BVSSQ <= 1.0e-1 ){
             
             NMAX2=1;
             
         } else {
             
             double BMAX = 0.0;  
             for ( int N=1; N<=42; N++ ){
                 
                 int NN = N-1;
                 double factorialNN = boost::math::factorial<double>(double(NN));
                 double BF=pow(2.0*beta*BVSSQ,NN)/factorialNN*G2;
                 double BF1 = std::abs(BF * XI(0,NN + 4));
                 double BF2 = std::abs(BF * XI(1,NN + 4));
                 double BF3 = std::abs(BF * XI(2,NN + 4));
                 BMAX = std::max(BF1,BF2);
                 BMAX = std::max(BMAX,BF3);
                 
                 if ( NMAX2 == 0 && BMAX <= conv ) NMAX2 = N;
                 if ( NMAX2 != 0 && BMAX >  conv ) NMAX2 = 0;
             
         }
             
             if (NMAX2 == 0 && BMAX > conv ) NMAX2 = 42;
             
         }
         
         
         //cout << "NMAX1 " << NMAX1 << " NMAX2 " << NMAX2 << endl;
         
         // something
         int INULL=1;
         if (AVSSQ <= 1.e-1 ) INULL=3;
         if (BVSSQ <= 1.e-1 ) INULL=4;   
         if (AVSSQ<=1.e-1 && BVSSQ <= 1.e-1 ) INULL=2;   
           
         
         //cout << INULL << endl;
         
         
         type_3D BLMA;
         type_3D CA;
         getBLMCOF( AVS, BLMA, CA );

          // cout << BLMA[1][1][3] << " " << BLMA[10][2][3] << " " << CA[1][3][3] << " " <<  CA[3][5][5] <<  endl;
         
         type_3D BLMB;
         type_3D CB;
         getBLMCOF( BVS, BLMB, CB );
         //         cout << BLMB[1][1][3] << " " << BLMB[10][2][3] << " " << CB[1][3][3] << " " <<  CB[3][5][5] <<  endl;
         
         
         typedef boost::multi_array_types::extent_range range;
         typedef type_3D::index index;
         type_3D::extent_gen extents3D;
         type_3D BLMAS;
         type_3D BLMBS;
         BLMAS.resize(extents3D[ range(1, 11) ][ range(1, 4) ][ range(1, 6)]);
         BLMBS.resize(extents3D[ range(1, 11) ][ range(1, 4) ][ range(1, 6)]);
         for ( index I = 1; I<=10; I++){
             for ( index L = 1; L<=3; L++){
                 for ( index M = 1; M<=5; M++){
                     BLMAS[I][L][M] = BLMA[I][L][M];
                     BLMBS[I][L][M] = BLMB[I][L][M];
                 }
             }
         } 



      BLMAS[10][1][3]=BLMA[1][1][3];
      BLMBS[10][1][3]=BLMB[1][1][3];
      BLMAS[10][2][2]=0.0;
      BLMAS[10][2][3]=0.0;
      BLMAS[10][2][4]=0.0;
      BLMBS[10][2][2]=0.0;
      BLMBS[10][2][3]=0.0;
      BLMBS[10][2][4]=0.0;
      
      
      
      
      
      
      type_3D CC;
      CC.resize(extents3D[ range(1,4)][range(1,6)][range(1,6)]);
       for ( index L = 1; L<=3; L++){
             for ( index M1 = 1; M1<=5; M1++){
                 for ( index M2 = 1; M2<=5; M2++){

                      CC[L][M1][M2]=0.0;
                      for ( index M = 1; M<=5; M++){
      
                        CC[L][M1][M2] += CA[L][M][M1]*CB[L][M][M2];
                      }
                 }
             }
       }
      
      
      //cout << CC[1][1][1] << " " << CC[3][5][5] << endl;
      
      
      
        
      
      typedef boost::multi_array<double, 5> type_5D;
      type_5D::extent_gen extents5D;
      type_5D SUMCI;
      SUMCI.resize(extents5D[ range(1,4)][ range(1,4)][ range(1,4)][range(1,6)][range(1,6)]);
      type_5D SUMCIN;
      SUMCIN.resize(extents5D[ range(1,4)][ range(1,4)][ range(1,4)][range(1,6)][range(1,6)]);
      type_5D SUMCID;
      SUMCID.resize(extents5D[ range(1,4)][ range(1,4)][ range(1,4)][range(1,6)][range(1,6)]);
      
      
      // awesome summations
      for ( index L = 1; L <= 3; L++  ){
          for ( index L1 = 1; L1 <= 3; L1++  ){                
              for ( index L2 = 1; L2 <= 3; L2++  ){
          
                  for ( index M1 = 1; M1 <=5; M1++ ){
                      for ( index M2 = 1; M2 <=5; M2++ ){
                          
                         
                                SUMCI[L][L1][L2][M1][M2]  = 0.0;
                                SUMCIN[L][L1][L2][M1][M2] = 0.0;
                                SUMCID[L][L1][L2][M1][M2] = 0.0;
                          

                                switch (INULL) {
                                    case 1:
                                    {
                                        for ( int N1 = 1; N1 <= NMAX1; N1++ ){
                                            int N = N1-1;
                                            double factorialN = boost::math::factorial<double>(double(N));
                                            double VAR1 = COEF[L][L1][M1][N1]*pow(2.0*alpha*AVSSQ,N)/factorialN;
                                            double VAR2 = 0.0;
                                            double VAR2D = 0.0;
                                            double VAR2N = 0.0;
                                            for ( int N2 = 1; N2 <= NMAX1; N2++ ){
                                                
                                                int NN   = N2 -1;
                                                int NPNS = N+NN-1+L1+L2;
                                                int NSN  = NPNS+2;
                                                int NSD = NPNS+4;
                                                double factorialNN = boost::math::factorial<double>(double(NN));
                                                double XDUM = COEF[L][L2][M2][N2]*pow(2.0*beta*BVSSQ,NN)/factorialNN;
                                                VAR2  += XDUM*XI(L-1,NPNS-1); // L index of XI starts with 0 !!
                                                VAR2N += XDUM*XI(L-1,NSN-1);  // N index of XI starts with 0 !!
                                                VAR2D += XDUM*XI(L-1,NSD-1);
                                                
                                                
                                            }
                                            
                                            SUMCI[L][L1][L2][M1][M2]  += VAR1*VAR2;
                                            SUMCIN[L][L1][L2][M1][M2] += VAR1*VAR2N;
                                            SUMCID[L][L1][L2][M1][M2] += VAR1*VAR2D;
                                            
                                            
                                        }
    
                                        
                                        break;
                                    }
                                    case 2:
                                    {

                                        int NL=L1+L2-1;
                                        SUMCI[L][L1][L2][M1][M2] = COEF[L][L1][M1][1]*COEF[L][L2][M2][1]*XI(L-1,NL-1);
                                        int NLN = NL+2;
                                        int NLD = NL+4;
                                        SUMCIN[L][L1][L2][M1][M2] = COEF[L][L1][M1][1]*COEF[L][L2][M2][1]*XI(L-1,NLN-1);
                                        SUMCID[L][L1][L2][M1][M2] = COEF[L][L1][M1][1]*COEF[L][L2][M2][1]*XI(L-1,NLD-1);

                                        break;
                                    }
                                    case 3:
                                    {
                                        double VAR2 = 0.0;
                                        double VAR2N = 0.0;
                                        double VAR2D = 0.0;
                                        for (int N2 = 1; N2 <= NMAX2; N2++) {
                                            int NN = N2 - 1;
                                            int NL = N2 + L2 + L1 - 2;
                                            int NLN = NL + 2;
                                            int NLD = NL + 4;
                                            double factorialNN = boost::math::factorial<double>(double(NN));
                                            double XDUM = COEF[L][L2][M2][N2] * pow(2.0 * beta*BVSSQ, NN) / factorialNN;
                                            VAR2 += XDUM * XI(L - 1, NL - 1);
                                            VAR2N += XDUM * XI(L - 1, NLN - 1);
                                            VAR2D += XDUM * XI(L - 1, NLD - 1);
                                        }

                                        SUMCI[L][L1][L2][M1][M2] = COEF[L][L1][M1][1] * VAR2;
                                        SUMCIN[L][L1][L2][M1][M2] = COEF[L][L1][M1][1] * VAR2N;
                                        SUMCID[L][L1][L2][M1][M2] = COEF[L][L1][M1][1] * VAR2D;
                                        
                                        break;
                                    }
                                    case 4:
                                    {
                                        double VAR1 = 0.0;
                                        double VAR1N = 0.0;
                                        double VAR1D = 0.0;

                                        for (int N1 = 1; N1 <= NMAX1; N1++) {
                                            int N = N1 - 1;
                                            int NL = N1 + L1 + L2 - 2;
                                            int NLN = NL + 2;
                                            int NLD = NL + 4;
                                            double factorialN = boost::math::factorial<double>(double(N));
                                            double XDUM = COEF[L][L1][M1][N1] * pow(2.0 * alpha*AVSSQ, N) / factorialN;
                                            VAR1 += XDUM * XI(L - 1, NL - 1);
                                            VAR1N += XDUM * XI(L - 1, NLN - 1);
                                            VAR1D += XDUM * XI(L - 1, NLD - 1);
                                        }
                                        SUMCI[L][L1][L2][M1][M2]  = COEF[L][L2][M2][1] * VAR1;
                                        SUMCIN[L][L1][L2][M1][M2] = COEF[L][L2][M2][1] * VAR1N;
                                        SUMCID[L][L1][L2][M1][M2] = COEF[L][L2][M2][1] * VAR1D;

                                        break;
                                    }

                                    default:
                                        cout << "Wrong ECP summation mode";
                                        exit(1);
                                } // switch
                      } // end M2
                  } // end M1
              } // end L2
          } // end L1 
      } // end L 
      
      
      
      //cout << SUMCI[1][1][1][1][1] << " " << SUMCI[3][3][3][5][5] << endl;
      //cout << SUMCIN[1][1][1][1][1] << " " << SUMCIN[3][3][3][5][5] << endl;
      //cout << SUMCID[1][1][1][1][1] << " " << SUMCID[3][3][3][5][5] << endl;
      
      
      // GET TRAFO HERE ALREADY
      std::vector<double> NormA(10);
      std::vector<double> NormB(10);
      getNorms(NormA,alpha);
      getNorms(NormB,beta);

      
      // now finally calculate matrix

            for (unsigned i = 0; i < matrix.size1(); i++) {
                for (unsigned j = 0; j < matrix.size2(); j++) {

                    for (index L = 1; L <= 3; L++) {
                        for (index L1 = 1; L1 <= 3; L1++) {
                            for (index L2 = 1; L2 <= 3; L2++) {

                                for (index M1 = 1; M1 <= 5; M1++) {
                                    for (index M2 = 1; M2 <= 5; M2++) {

                                        
                                        matrix(i,j) += BLMA[i+1][L1][M1] * BLMB[j+1][L2][M2] * SUMCI[L][L1][L2][M1][M2]*CC[L][M1][M2];

                                    }
                                }
                            }
                        }
                    }
                }


            }

          
      //cout << matrix(0,0) <<  " " << matrix(9,9) << endl;
      
      
      
        for (unsigned i = 0; i < matrix.size1(); i++) {
                for (unsigned j = 0; j < matrix.size2(); j++) {

                    matrix(i,j) = matrix(i,j) * GAUSS * NormA[i] * NormB[j];

                }
        }
        //cout << matrix(0,0) <<  " " << matrix(9,9) << endl;
      //exit(0);
      
        return matrix;
        
        
    }
    
    
    void AOECP::getNorms(std::vector<double>& Norms, const double& decay){

            const double PI = boost::math::constants::pi<double>();
            Norms[0] = pow(2.0 * decay / PI, 0.75);
            Norms[1] = 2.0 * sqrt(decay) * Norms[0];
            Norms[2] = Norms[1];
            Norms[3] = Norms[1];
            Norms[4] = 4.00 * decay * Norms[0];
            Norms[5] = Norms[4];
            Norms[6] = Norms[4];
            Norms[7] = 2.00 * decay * Norms[0] / sqrt(3.0);
            Norms[8] = 2.00 * decay * Norms[0];
            Norms[9] = Norms[4] / sqrt(15.0);
        
        
    }
    
    void AOECP::getBLMCOF(const vec& pos, type_3D& BLC, type_3D& C){
        
        typedef boost::multi_array_types::extent_range range;
        typedef type_3D::index index;
        type_3D::extent_gen extents;

        BLC.resize(extents[ range(1, 11) ][ range(1, 4) ][ range(1, 6)]);
        C.resize(extents[ range(1, 4) ][ range(1, 6) ][ range(1, 6)]);

        type_3D BLM;
        BLM.resize(extents[ range(1, 11) ][ range(1, 4) ][ range(1, 6)]);
          
        const double PI = boost::math::constants::pi<double>();
        
               
      double SPI=sqrt(PI);
      double XS=2.0*SPI;          
      double XP=XS/sqrt(3.0);
      double XD1 = XP/sqrt(5.0); 
      double XD2 = 4.0*SPI/sqrt(5.0);
      double XD3 = XD2/sqrt(3.0);
      
      for ( index I = 1; I <=10; I++ ){
          for ( index L = 1; L <=3; L++ ){
              for ( index M = 1; M <=5; M++ ){
                  
                  BLM[I][L][M] = 0.0;
                  
              }
          }
          
      }
      
      std::vector<double> BVS(4);
      BVS[0] = 0.0;
      BVS[1] = pos.getX(); 
      BVS[2] = pos.getY(); 
      BVS[3] = pos.getZ(); 
      
      BLM[1][1][3] = XS;
      BLM[2][1][3] = -BVS[1]*XS;
      BLM[2][2][4] = XP;
      BLM[3][1][3] = -BVS[2]*XS;
      BLM[3][2][2] = XP;
      BLM[4][1][3] = -BVS[3]*XS;
      BLM[4][2][3] = XP;
      BLM[5][1][3] = BVS[1]*BVS[3]*XS;
      BLM[5][2][4] = -BVS[3]*XP;
      BLM[5][2][3] = -BVS[1]*XP;
      BLM[5][3][4] = XD1;
      BLM[6][1][3] = BVS[2]*BVS[3]*XS;
      BLM[6][2][2] = -BVS[3]*XP;
      BLM[6][2][3] = -BVS[2]*XP;
      BLM[6][3][2] = XD1;
      BLM[7][1][3] = BVS[1]*BVS[2]*XS;
      BLM[7][2][4] = -BVS[2]*XP;
      BLM[7][2][2] = -BVS[1]*XP;
      BLM[7][3][1] = XD1;
      BLM[8][1][3] = (2.0*BVS[3]*BVS[3]-BVS[1]*BVS[1]-BVS[2]*BVS[2])*XS;
      BLM[8][2][4] = 2.0*BVS[1]*XP;
      BLM[8][2][2] = 2.0*BVS[2]*XP;
      BLM[8][2][3] = -4.0*BVS[3]*XP;
      BLM[8][3][3] = XD2;
      BLM[9][1][3] = (BVS[1]*BVS[1]-BVS[2]*BVS[2])*XS;
      BLM[9][2][4] = -2.0*BVS[1]*XP;
      BLM[9][2][2] = 2.0*BVS[2]*XP;
      BLM[9][3][5] = XD3;
      BLM[10][1][3] = (BVS[1]*BVS[1]+BVS[2]*BVS[2]+BVS[3]*BVS[3])*XS;
      BLM[10][2][4] = -2.0*BVS[1]*XP;
      BLM[10][2][2] = -2.0*BVS[2]*XP;
      BLM[10][2][3] = -2.0*BVS[3]*XP;



            for (index L = 1; L <= 3; L++) {
                for (index M = 1; M <= 5; M++) {
                    for (index MM = 1; MM <= 5; MM++) {

                        C[L][M][MM] = 0.0;

                    }
                }

            }
      

      double XY = BVS[1]*BVS[1]+BVS[2]*BVS[2];
      double XYZ=XY+BVS[3]*BVS[3];
      double SXY=sqrt(XY);
      double SXYZ=sqrt(XYZ);
      double CM=1.0;
      double CP=0.0;

      if ( SXY > 1.e-4 ) {
          
          CM = BVS[1]/SXY;
          CP = BVS[2]/SXY;
          
      }

            if (SXYZ > 1.e-4) {

                double CL = BVS[3] / SXYZ;
                double CN = SXY / SXYZ;

                C[1][3][3] = 1.0;
                C[2][2][2] = CM;
                C[2][2][3] = CN*CP;
                C[2][2][4] = CL*CP;
                C[2][3][2] = 0.0;
                C[2][3][3] = CL;
                C[2][3][4] = -CN;
                C[2][4][2] = -CP;
                C[2][4][3] = CM*CN;
                C[2][4][4] = CL*CM;
                C[3][1][1] = CL * (2.0 * CM * CM - 1.0);
                C[3][1][2] = CN * (2.0 * CM * CM - 1.0);
                double SQ3 = sqrt(3.0);
                C[3][1][3] = SQ3 * CM * CP * CN*CN;
                C[3][1][4] = 2.0 * CL * CM * CP*CN;
                C[3][1][5] = CM * CP * (1.0 + CL * CL);
                C[3][2][1] = -CM*CN;
                C[3][2][2] = CL*CM;
                C[3][2][3] = SQ3 * CL * CN*CP;
                C[3][2][4] = CP * (2.0 * CL * CL - 1.0);
                C[3][2][5] = -CL * CN*CP;
                C[3][3][1] = 0.0;
                C[3][3][2] = 0.0;
                C[3][3][3] = 1.5 * CL * CL - 0.5;
                C[3][3][4] = -SQ3 * CL*CN;
                C[3][3][5] = .5 * SQ3 * (1.0 - CL * CL);
                C[3][4][1] = CN*CP;
                C[3][4][2] = -CL*CP;
                C[3][4][3] = SQ3 * CL * CM*CN;
                C[3][4][4] = CM * (2.0 * CL * CL - 1.0);
                C[3][4][5] = -CL * CM*CN;
                C[3][5][1] = -2.0 * CL * CM*CP;
                C[3][5][2] = -2.0 * CM * CN*CP;
                C[3][5][3] = 0.5 * SQ3 * (CL * CL * (1.0 - 2.0 * CM * CM) + 2.0 * CM * CM - 1.0);
                C[3][5][4] = -CL * CN * (1.0 - 2.0 * CM * CM);
                C[3][5][5] = 0.5 * (CM * CM * (2.0 + 2.0 * CL * CL) - CL * CL - 1.0);



            } else {


                C[1][3][3] = 1.0;
                C[2][2][2] = 1.0;
                C[2][3][3] = 1.0;
                C[2][4][4] = 1.0;


                for (index M = 1; M <= 5; M++) {

                    C[3][M][M] = 1.0;

                }

            }
      
      
         for ( index I = 1; I <=10; I++ ){
          for ( index L = 1; L <=3; L++ ){
              for ( index M = 1; M <=5; M++ ){
                  
                  BLC[I][L][M] = 0.0;
                  for ( index M1 = 1; M1 <=5; M1++ ){
                      
                      BLC[I][L][M] += BLM[I][L][M1] * C[L][M1][M];
                      
                  }
                  
              }
          }
          
      }

      
       
        
        
    } // getBLMCOF
    
    
    
}}

