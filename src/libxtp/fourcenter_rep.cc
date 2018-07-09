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

// Overload of uBLAS prod function with MKL/GSL implementations


#include <votca/xtp/fourcenter.h>



namespace votca {
    namespace xtp {


 
        
        /*
         * Calculate 4-center electron repulsion integrals 
         *    R_{abcd} = int{ phi_a(r) phi_b(r) phi_c(r') phi_d(r') /(r-r') d3rd3r' }
         * for a given set of a b c d as in 
         *    Obara, Saika, J. Chem. Phys. 84, 3963 (1986)
         * section II.B for cartesian Gaussians, then transforming
         * to spherical (angular momentum) Gaussians ("complete" shells
         * from S to Lmax, and finally cutting out those angular momentum 
         * components actually present in shell-shell-shell combination.
         * Currently supported for 
         *      S,P,D,F,G  functions in DFT basis
         * 
         */

      
        bool FCMatrix::FillFourCenterRepBlock(tensor4d& block,
                const AOShell* _shell_1, const AOShell* _shell_2, const AOShell* _shell_3, const AOShell* _shell_4) {

            const double pi = boost::math::constants::pi<double>();
            
            bool _does_contribute=true;
            

            // shell info, only lmax tells how far to go
            
            int _lmax_1 = _shell_1->getLmax();
            int _lmax_2 = _shell_2->getLmax();
            int _lmax_3 = _shell_3->getLmax();
            int _lmax_4 = _shell_4->getLmax();
            
            int _mmax = _lmax_1 + _lmax_2 + _lmax_3 + _lmax_4;


            // set size of internal block for recursion
           
            const AOShell* _shell_alpha;
            const AOShell* _shell_beta;
            const AOShell* _shell_gamma;
            const AOShell* _shell_delta;
            bool alphabetaswitch = false;
            bool gammadeltaswitch = false;
            bool ab_cd_switch = false;

            // We need lmax_alpha > lmax_beta, so instead of calculating (sp,s) we calculate (ps,s), due to symmetry they are the same. 

            if (_lmax_1 < _lmax_2) {
              alphabetaswitch = true;               
            }
            if (_lmax_3 < _lmax_4) {
              gammadeltaswitch = true;               
            }
            if ( (_lmax_1 + _lmax_2) < (_lmax_3 + _lmax_4) ) {
              ab_cd_switch = true;               
            }

            if (ab_cd_switch==true) {

              if (alphabetaswitch== true) {
                _shell_alpha = _shell_4;
                _shell_beta = _shell_3;              
              } else {
                _shell_alpha = _shell_3;
                _shell_beta = _shell_4;
              }

              if (gammadeltaswitch==true) {
                _shell_gamma = _shell_2;
                _shell_delta = _shell_1;               
              } else {
                _shell_gamma = _shell_1;
                _shell_delta = _shell_2;
              }

            } else {

              if (alphabetaswitch== true) {
                _shell_alpha = _shell_2;
                _shell_beta = _shell_1;              
              } else {
                _shell_alpha = _shell_1;
                _shell_beta = _shell_2;
              }

              if (gammadeltaswitch==true) {
                _shell_gamma = _shell_4;
                _shell_delta = _shell_3;               
              } else {
                _shell_gamma = _shell_3;
                _shell_delta = _shell_4;
              }

            } 




            const tools::vec& _pos_alpha = _shell_alpha->getPos();
            const tools::vec& _pos_beta = _shell_beta->getPos();
            const tools::vec& _pos_gamma = _shell_gamma->getPos();
            const tools::vec& _pos_delta = _shell_delta->getPos();
            
            int _lmax_alpha = _shell_alpha->getLmax();
            int _lmax_beta  = _shell_beta->getLmax();
            int _lmax_gamma = _shell_gamma->getLmax();
            int _lmax_delta = _shell_delta->getLmax();

            int n_orbitals[] = { 1, 4, 10, 20, 35, 56, 84, 120, 165 };     //   n_orbitals[n] = ( (n + 1) * (n + 2) * (n + 3) ) / 6


 int nx[] = { 0,
              1, 0, 0,
              2, 1, 1, 0, 0, 0,
              3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
              4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0,
              5, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
              6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
              7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
              8, 7, 7, 6, 6, 6, 5, 5, 5, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

 int ny[] = { 0,
              0, 1, 0,
              0, 1, 0, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2, 1, 0, 8, 7, 6, 5, 4, 3, 2, 1, 0 };

 int nz[] = { 0,
              0, 0, 1,
              0, 0, 1, 0, 1, 2,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8 };


 int i_less_x[] = {
  0,
  0,  0,  0,
  1,  2,  3,  0,  0,  0,
  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,
 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  0,  0,  0,  0,  0,
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,  0,  0,  0,  0,  0,  0,
 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,  0,  0,  0,  0,  0,  0,  0,
 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83,  0,  0,  0,  0,  0,  0,  0,  0,
 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,  0,  0,  0,  0,  0,  0,  0,  0,  0 };

 int i_less_y[] = {
  0,
  0,  0,  0,
  0,  1,  0,  2,  3,  0,
  0,  4,  0,  5,  6,  0,  7,  8,  9,  0,
  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,  0,
  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,  0,
  0, 35,  0, 36, 37,  0, 38, 39, 40,  0, 41, 42, 43, 44,  0, 45, 46, 47, 48, 49,  0, 50, 51, 52, 53, 54, 55,  0,
  0, 56,  0, 57, 58,  0, 59, 60, 61,  0, 62, 63, 64, 65,  0, 66, 67, 68, 69, 70,  0, 71, 72, 73, 74, 75, 76,  0, 77, 78, 79, 80, 81, 82, 83,  0,
  0, 84,  0, 85, 86,  0, 87, 88, 89,  0, 90, 91, 92, 93,  0, 94, 95, 96, 97, 98,  0, 99,100,101,102,103,104,  0,105,106,107,108,109,110,111,  0,112,113,114,115,116,117,118,119,  0 };

 int i_less_z[] = {
  0,
  0,  0,  0,
  0,  0,  1,  0,  2,  3,
  0,  0,  4,  0,  5,  6,  0,  7,  8,  9,
  0,  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,
  0,  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,
  0,  0, 35,  0, 36, 37,  0, 38, 39, 40,  0, 41, 42, 43, 44,  0, 45, 46, 47, 48, 49,  0, 50, 51, 52, 53, 54, 55,
  0,  0, 56,  0, 57, 58,  0, 59, 60, 61,  0, 62, 63, 64, 65,  0, 66, 67, 68, 69, 70,  0, 71, 72, 73, 74, 75, 76,  0, 77, 78, 79, 80, 81, 82, 83,
  0,  0, 84,  0, 85, 86,  0, 87, 88, 89,  0, 90, 91, 92, 93,  0, 94, 95, 96, 97, 98,  0, 99,100,101,102,103,104,  0,105,106,107,108,109,110,111,  0,112,113,114,115,116,117,118,119 };


 int i_more_x[] = {  1,
                     4,  5,  6,
                    10, 11, 12, 13, 14, 15,
                    20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
                    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
                    56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76,
                    84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99,100,101,102,103,104,105,106,107,108,109,110,111,
                   120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155 };

 int i_more_y[] = {  2,
                     5,  7,  8,
                    11, 13, 14, 16, 17, 18,
                    21, 23, 24, 26, 27, 28, 30, 31, 32, 33,
                    36, 38, 39, 41, 42, 43, 45, 46, 47, 48, 50, 51, 52, 53, 54,
                    57, 59, 60, 62, 63, 64, 66, 67, 68, 69, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82,
                    85, 87, 88, 90, 91, 92, 94, 95, 96, 97, 99,100,101,102,103,105,106,107,108,109,110,112,113,114,115,116,117,118,
                   121,123,124,126,127,128,130,131,132,133,135,136,137,138,139,141,142,143,144,145,146,148,149,150,151,152,153,154,156,157,158,159,160,161,162,163 };

 int i_more_z[] = {  3,
                     6,  8,  9,
                    12, 14, 15, 17, 18, 19,
                    22, 24, 25, 27, 28, 29, 31, 32, 33, 34,
                    37, 39, 40, 42, 43, 44, 46, 47, 48, 49, 51, 52, 53, 54, 55,
                    58, 60, 61, 63, 64, 65, 67, 68, 69, 70, 72, 73, 74, 75, 76, 78, 79, 80, 81, 82, 83,
                    86, 88, 89, 91, 92, 93, 95, 96, 97, 98,100,101,102,103,104,106,107,108,109,110,111,113,114,115,116,117,118,119,
                   122,124,125,127,128,129,131,132,133,134,136,137,138,139,140,142,143,144,145,146,147,149,150,151,152,153,154,155,157,158,159,160,161,162,163,164 };



            int _nbeta = AOSuperMatrix::getBlockSize(_lmax_beta);
            int _ndelta = AOSuperMatrix::getBlockSize(_lmax_delta);
            int _ncombined_ab = AOSuperMatrix::getBlockSize(_lmax_alpha+_lmax_beta);
            int _ncombined_cd = AOSuperMatrix::getBlockSize(_lmax_gamma+_lmax_delta);
            
            tensor3d::extent_gen extents;
            tensor4d::extent_gen extents4;

            double _dist_AB = (_pos_alpha - _pos_beta) * (_pos_alpha - _pos_beta);
            double _dist_CD = (_pos_gamma - _pos_delta) * (_pos_gamma - _pos_delta);
            
            tools::vec amb = _pos_alpha - _pos_beta;
            double amb0 = 0.0;
            double amb1 = 0.0;
            double amb2 = 0.0;
            if (_dist_AB > 0.03) {
              amb0 = amb.getX();
              amb1 = amb.getY();
              amb2 = amb.getZ();
            }

            tools::vec cmd = _pos_gamma - _pos_delta;
            double cmd0 = 0.0;
            double cmd1 = 0.0;
            double cmd2 = 0.0;
            if (_dist_CD > 0.03) {
              cmd0 = cmd.getX();
              cmd1 = cmd.getY();
              cmd2 = cmd.getZ();
            }




            for ( AOShell::GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha) {
                const double _decay_alpha = italpha->getDecay();
            
              for ( AOShell::GaussianIterator itbeta = _shell_beta->firstGaussian(); itbeta != _shell_beta->lastGaussian(); ++itbeta) {
                  const double _decay_beta = itbeta->getDecay();
                    
                for ( AOShell::GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma) {
                    const double _decay_gamma = itgamma->getDecay();

                  for ( AOShell::GaussianIterator itdelta = _shell_delta->firstGaussian(); itdelta != _shell_delta->lastGaussian(); ++itdelta) {
                      const double _decay_delta = itdelta->getDecay();


          
            double zeta = _decay_alpha + _decay_beta;    
            double eta = _decay_gamma + _decay_delta;     
            double _decay = zeta + eta;
            double rho = (zeta*eta)/_decay;
            double rzeta = 0.5/zeta;
            double reta = 0.5/eta;
            double rdecay = 0.5/_decay;
            double gfak = eta/_decay;
            double cfak = zeta/_decay;         
            tools::vec _P = (_decay_alpha*_pos_alpha + _decay_beta*_pos_beta)/zeta;
            tools::vec _Q = (_decay_gamma*_pos_gamma + _decay_delta*_pos_delta)/eta;
            tools::vec _W = (zeta*_P + eta*_Q)/_decay;
            double _T = rho*(_P-_Q)*(_P-_Q);


            tools::vec pma = _P - _pos_alpha;
            double pma0 = 0.0;
            double pma1 = 0.0;
            double pma2 = 0.0;
            if (_dist_AB > 0.03) {
              pma0 = pma.getX();
              pma1 = pma.getY();
              pma2 = pma.getZ();
            }

            tools::vec qmc = _Q - _pos_gamma;
            double qmc0 = 0.0;
            double qmc1 = 0.0;
            double qmc2 = 0.0;
            if (_dist_CD > 0.03) {
              qmc0 = qmc.getX();
              qmc1 = qmc.getY();
              qmc2 = qmc.getZ();
            }

            tools::vec wmp = _W - _P;
            double wmp0 = wmp.getX();
            double wmp1 = wmp.getY();
            double wmp2 = wmp.getZ();

            tools::vec wmq = _W - _Q;
            double wmq0 = wmq.getX();
            double wmq1 = wmq.getY();
            double wmq2 = wmq.getZ();



            tensor3d R_temp;
            R_temp.resize(extents[ range(0, _ncombined_ab ) ][ range(0, _ncombined_cd ) ][ range(0, _mmax+1)]);
            //initialize to zero
            for (index3d i = 0; i != _ncombined_ab; ++i) {
              for (index3d j = 0; j != _ncombined_cd; ++j) {
                for (index3d k = 0; k != _mmax+1; ++k) {

                  R_temp[i][j][k] = 0.0;

                }
              }
             }


            tensor3d R;
            R.resize(extents[ range(0, _ncombined_ab ) ][ range(0, _nbeta ) ][ range(0, _ncombined_cd)]);
            //initialize to zero
            for (index3d i = 0; i != _ncombined_ab; ++i) {
              for (index3d j = 0; j != _nbeta; ++j) {
                for (index3d k = 0; k != _ncombined_cd; ++k) {

                  R[i][j][k] = 0.0;

                }
              }
            }
            

            const std::vector<double> _FmT=AOMatrix<double>::XIntegrate(_mmax+1, _T);

            double exp_AB = exp( -2. * _decay_alpha * _decay_beta * rzeta * _dist_AB );
            double exp_CD = exp( -2.* _decay_gamma * _decay_delta * reta * _dist_CD );
            double ssss = ( 16. * pow( _decay_alpha * _decay_beta * _decay_gamma * _decay_delta, .75 ) * exp_AB * exp_CD )
                         / ( zeta * eta * sqrt(pi * _decay) );

            //ss integrals
            for (int _i=0; _i < _mmax+1; _i++){
              R_temp[0][0][_i] = ssss * _FmT[_i];
            }

int _lmax_alpha_beta = _lmax_alpha + _lmax_beta;
int _lmax_gamma_delta = _lmax_gamma + _lmax_delta;


//Integrals     p-s - s-s
if (_lmax_alpha_beta > 0) {
  for (int m = 0; m < _mmax; m++) {
    R_temp[Cart::x][0][m] = pma0*R_temp[0][0][m] + wmp0*R_temp[0][0][m+1];
    R_temp[Cart::y][0][m] = pma1*R_temp[0][0][m] + wmp1*R_temp[0][0][m+1];
    R_temp[Cart::z][0][m] = pma2*R_temp[0][0][m] + wmp2*R_temp[0][0][m+1];
  }
}
//------------------------------------------------------

//Integrals     d-s - s-s
if (_lmax_alpha_beta > 1) {
  for (int m = 0; m < _mmax-1; m++) {
    double term = rzeta*(R_temp[0][0][m]-gfak*R_temp[0][0][m+1]);
    R_temp[Cart::xx][0][m] = pma0*R_temp[Cart::x][0][m] + wmp0*R_temp[Cart::x][0][m+1] + term;
    R_temp[Cart::xy][0][m] = pma0*R_temp[Cart::y][0][m] + wmp0*R_temp[Cart::y][0][m+1];
    R_temp[Cart::xz][0][m] = pma0*R_temp[Cart::z][0][m] + wmp0*R_temp[Cart::z][0][m+1];
    R_temp[Cart::yy][0][m] = pma1*R_temp[Cart::y][0][m] + wmp1*R_temp[Cart::y][0][m+1] + term;
    R_temp[Cart::yz][0][m] = pma1*R_temp[Cart::z][0][m] + wmp1*R_temp[Cart::z][0][m+1];
    R_temp[Cart::zz][0][m] = pma2*R_temp[Cart::z][0][m] + wmp2*R_temp[Cart::z][0][m+1] + term;
  }
}
//------------------------------------------------------

//Integrals     f-s - s-s
if (_lmax_alpha_beta > 2) {
  for (int m = 0; m < _mmax-2; m++) {
    R_temp[Cart::xxx][0][m] = pma0*R_temp[Cart::xx][0][m] + wmp0*R_temp[Cart::xx][0][m+1] + 2*rzeta*(R_temp[Cart::x][0][m]-gfak*R_temp[Cart::x][0][m+1]);
    R_temp[Cart::xxy][0][m] = pma1*R_temp[Cart::xx][0][m] + wmp1*R_temp[Cart::xx][0][m+1];
    R_temp[Cart::xxz][0][m] = pma2*R_temp[Cart::xx][0][m] + wmp2*R_temp[Cart::xx][0][m+1];
    R_temp[Cart::xyy][0][m] = pma0*R_temp[Cart::yy][0][m] + wmp0*R_temp[Cart::yy][0][m+1];
    R_temp[Cart::xyz][0][m] = pma0*R_temp[Cart::yz][0][m] + wmp0*R_temp[Cart::yz][0][m+1];
    R_temp[Cart::xzz][0][m] = pma0*R_temp[Cart::zz][0][m] + wmp0*R_temp[Cart::zz][0][m+1];
    R_temp[Cart::yyy][0][m] = pma1*R_temp[Cart::yy][0][m] + wmp1*R_temp[Cart::yy][0][m+1] + 2*rzeta*(R_temp[Cart::y][0][m]-gfak*R_temp[Cart::y][0][m+1]);
    R_temp[Cart::yyz][0][m] = pma2*R_temp[Cart::yy][0][m] + wmp2*R_temp[Cart::yy][0][m+1];
    R_temp[Cart::yzz][0][m] = pma1*R_temp[Cart::zz][0][m] + wmp1*R_temp[Cart::zz][0][m+1];
    R_temp[Cart::zzz][0][m] = pma2*R_temp[Cart::zz][0][m] + wmp2*R_temp[Cart::zz][0][m+1] + 2*rzeta*(R_temp[Cart::z][0][m]-gfak*R_temp[Cart::z][0][m+1]);
  }
}
//------------------------------------------------------

//Integrals     g-s - s-s     h-s - s-s     i-s - s-s     j-s - s-s     k-s - s-s     . . .
for (int l = 4; l < _lmax_alpha_beta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int norb_3 = n_orbitals[l-3];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  int ncart_3 = ncart_2 + 1 - l;
  for (int m = 0; m < _mmax+1-l; m++) {
    R_temp[norb_1][0][m] = pma0*R_temp[norb_2][0][m] + wmp0*R_temp[norb_2][0][m+1] + (l-1)*rzeta*(R_temp[norb_3][0][m]-gfak*R_temp[norb_3][0][m+1]);
    R_temp[norb_1 + 1][0][m] = pma1*R_temp[norb_2][0][m] + wmp1*R_temp[norb_2][0][m+1];
    R_temp[norb_1 + 2][0][m] = pma2*R_temp[norb_2][0][m] + wmp2*R_temp[norb_2][0][m+1];
    int ntimes = 3;
    int itimes = 3;
    for (int k = 3; k < ncart_2; k++) {
      R_temp[norb_1 + k][0][m] = pma0*R_temp[norb_2 + k][0][m] + wmp0*R_temp[norb_2 + k][0][m+1] + (l-ntimes)*rzeta*(R_temp[norb_3 + k][0][m]-gfak*R_temp[norb_3 + k][0][m+1]);
      itimes--;
      if (itimes == 0) {
        ntimes++;
        itimes = ntimes;
      }
    }
    for (int k = 0; k < l-1; k++) {
      int k2 = norb_2 + ncart_2 + k;
      R_temp[norb_1 + ncart_2 + k][0][m] = pma0*R_temp[k2][0][m] + wmp0*R_temp[k2][0][m+1];
      R_temp[norb_1 + ncart_1 + k][0][m] = pma1*R_temp[k2][0][m] + wmp1*R_temp[k2][0][m+1]
                                           + (l-1-k)*rzeta*(R_temp[norb_3 + ncart_3 + k][0][m]-gfak*R_temp[norb_3 + ncart_3 + k][0][m+1]);
    }
    R_temp[norb_1 + ncart_2 + l - 1][0][m] = pma0*R_temp[norb_2 + ncart_2 + l - 1][0][m] + wmp0*R_temp[norb_2 + ncart_2 + l - 1][0][m+1];
    R_temp[norb-2][0][m] = pma1*R_temp[norb_1-1][0][m] + wmp1*R_temp[norb_1-1][0][m+1];
    R_temp[norb-1][0][m] = pma2*R_temp[norb_1-1][0][m] + wmp2*R_temp[norb_1-1][0][m+1] + (l-1)*rzeta*(R_temp[norb_2-1][0][m]-gfak*R_temp[norb_2-1][0][m+1]);
  }
}
//------------------------------------------------------




if (_lmax_gamma_delta > 0) {

  //Integrals     s-s - p-s
  for (int m = 0; m < _lmax_gamma_delta; m++) {
    R_temp[0][Cart::x][m] = qmc0*R_temp[0][0][m] + wmq0*R_temp[0][0][m+1];
    R_temp[0][Cart::y][m] = qmc1*R_temp[0][0][m] + wmq1*R_temp[0][0][m+1];
    R_temp[0][Cart::z][m] = qmc2*R_temp[0][0][m] + wmq2*R_temp[0][0][m+1];
  }
  //------------------------------------------------------

  //Integrals     p-s - p-s
  if (_lmax_alpha_beta > 0) {
    for (int m = 0; m < _lmax_gamma_delta; m++) {
      double term = rdecay*R_temp[0][0][m+1];
      for (int _i =  1; _i < 4; _i++) {
        R_temp[_i][Cart::x][m] = qmc0*R_temp[_i][0][m] + wmq0*R_temp[_i][0][m+1] + nx[_i]*term;
        R_temp[_i][Cart::y][m] = qmc1*R_temp[_i][0][m] + wmq1*R_temp[_i][0][m+1] + ny[_i]*term;
        R_temp[_i][Cart::z][m] = qmc2*R_temp[_i][0][m] + wmq2*R_temp[_i][0][m+1] + nz[_i]*term;
      }
    }
  }
  //------------------------------------------------------

  //Integrals     d-s - p-s     f-s - p-s     g-s - p-s     h-s - p-s     i-s - p-s     j-s - p-s     k-s - p-s     . . .
  for (int m = 0; m < _lmax_gamma_delta; m++) {
    for (int _i = 4; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      R_temp[_i][Cart::x][m] = qmc0*R_temp[_i][0][m] + wmq0*R_temp[_i][0][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][0][m+1];
      R_temp[_i][Cart::y][m] = qmc1*R_temp[_i][0][m] + wmq1*R_temp[_i][0][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][0][m+1];
      R_temp[_i][Cart::z][m] = qmc2*R_temp[_i][0][m] + wmq2*R_temp[_i][0][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][0][m+1];
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma_delta > 0)


if (_lmax_gamma_delta > 1) {

  //Integrals     s-s - d-s
  for (int m = 0; m < _lmax_gamma_delta-1; m++) {
    double term = reta*(R_temp[0][0][m]-cfak*R_temp[0][0][m+1]);
    R_temp[0][Cart::xx][m] = qmc0*R_temp[0][Cart::x][m] + wmq0*R_temp[0][Cart::x][m+1] + term;
    R_temp[0][Cart::xy][m] = qmc0*R_temp[0][Cart::y][m] + wmq0*R_temp[0][Cart::y][m+1];
    R_temp[0][Cart::xz][m] = qmc0*R_temp[0][Cart::z][m] + wmq0*R_temp[0][Cart::z][m+1];
    R_temp[0][Cart::yy][m] = qmc1*R_temp[0][Cart::y][m] + wmq1*R_temp[0][Cart::y][m+1] + term;
    R_temp[0][Cart::yz][m] = qmc1*R_temp[0][Cart::z][m] + wmq1*R_temp[0][Cart::z][m+1];
    R_temp[0][Cart::zz][m] = qmc2*R_temp[0][Cart::z][m] + wmq2*R_temp[0][Cart::z][m+1] + term;
  }
  //------------------------------------------------------

  //Integrals     p-s - d-s     d-s - d-s     f-s - d-s     g-s - d-s     h-s - d-s     i-s - d-s     j-s - d-s     k-s - d-s     . . .
  for (int m = 0; m < _lmax_gamma_delta-1; m++) {
    for (int _i = 1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
        double term = reta*(R_temp[_i][0][m]-cfak*R_temp[_i][0][m+1]);
        R_temp[_i][Cart::xx][m] = qmc0*R_temp[_i][Cart::x][m] + wmq0*R_temp[_i][Cart::x][m+1] +nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::x][m+1] + term;
        R_temp[_i][Cart::xy][m] = qmc0*R_temp[_i][Cart::y][m] + wmq0*R_temp[_i][Cart::y][m+1] +nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::y][m+1];
        R_temp[_i][Cart::xz][m] = qmc0*R_temp[_i][Cart::z][m] + wmq0*R_temp[_i][Cart::z][m+1] +nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::z][m+1];
        R_temp[_i][Cart::yy][m] = qmc1*R_temp[_i][Cart::y][m] + wmq1*R_temp[_i][Cart::y][m+1] +ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::y][m+1] + term;
        R_temp[_i][Cart::yz][m] = qmc1*R_temp[_i][Cart::z][m] + wmq1*R_temp[_i][Cart::z][m+1] +ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::z][m+1];
        R_temp[_i][Cart::zz][m] = qmc2*R_temp[_i][Cart::z][m] + wmq2*R_temp[_i][Cart::z][m+1] +nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::z][m+1] + term;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma_delta > 1)


if (_lmax_gamma_delta > 2) {

  //Integrals     s-s - f-s
  for (int m = 0; m < _lmax_gamma_delta-2; m++) {
    R_temp[0][Cart::xxx][m] = qmc0*R_temp[0][Cart::xx][m] + wmq0*R_temp[0][Cart::xx][m+1] + 2*reta*(R_temp[0][Cart::x][m]-cfak*R_temp[0][Cart::x][m+1]);
    R_temp[0][Cart::xxy][m] = qmc1*R_temp[0][Cart::xx][m] + wmq1*R_temp[0][Cart::xx][m+1];
    R_temp[0][Cart::xxz][m] = qmc2*R_temp[0][Cart::xx][m] + wmq2*R_temp[0][Cart::xx][m+1];
    R_temp[0][Cart::xyy][m] = qmc0*R_temp[0][Cart::yy][m] + wmq0*R_temp[0][Cart::yy][m+1];
    R_temp[0][Cart::xyz][m] = qmc0*R_temp[0][Cart::yz][m] + wmq0*R_temp[0][Cart::yz][m+1];
    R_temp[0][Cart::xzz][m] = qmc0*R_temp[0][Cart::zz][m] + wmq0*R_temp[0][Cart::zz][m+1];
    R_temp[0][Cart::yyy][m] = qmc1*R_temp[0][Cart::yy][m] + wmq1*R_temp[0][Cart::yy][m+1] + 2*reta*(R_temp[0][Cart::y][m]-cfak*R_temp[0][Cart::y][m+1]);
    R_temp[0][Cart::yyz][m] = qmc2*R_temp[0][Cart::yy][m] + wmq2*R_temp[0][Cart::yy][m+1];
    R_temp[0][Cart::yzz][m] = qmc1*R_temp[0][Cart::zz][m] + wmq1*R_temp[0][Cart::zz][m+1];
    R_temp[0][Cart::zzz][m] = qmc2*R_temp[0][Cart::zz][m] + wmq2*R_temp[0][Cart::zz][m+1] + 2*reta*(R_temp[0][Cart::z][m]-cfak*R_temp[0][Cart::z][m+1]);
  }
  //------------------------------------------------------

  //Integrals     p-s - f-s     d-s - f-s     f-s - f-s     g-s - f-s     h-s - f-s     i-s - f-s     j-s - f-s     k-s - f-s     . . .
  for (int m = 0; m < _lmax_gamma_delta-2; m++) {
    for (int _i = 1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      double term_x = 2*reta*(R_temp[_i][Cart::x][m]-cfak*R_temp[_i][Cart::x][m+1]);
      double term_y = 2*reta*(R_temp[_i][Cart::y][m]-cfak*R_temp[_i][Cart::y][m+1]);
      double term_z = 2*reta*(R_temp[_i][Cart::z][m]-cfak*R_temp[_i][Cart::z][m+1]);
      R_temp[_i][Cart::xxx][m] = qmc0*R_temp[_i][Cart::xx][m] + wmq0*R_temp[_i][Cart::xx][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xx][m+1] + term_x;
      R_temp[_i][Cart::xxy][m] = qmc1*R_temp[_i][Cart::xx][m] + wmq1*R_temp[_i][Cart::xx][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xx][m+1];
      R_temp[_i][Cart::xxz][m] = qmc2*R_temp[_i][Cart::xx][m] + wmq2*R_temp[_i][Cart::xx][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xx][m+1];
      R_temp[_i][Cart::xyy][m] = qmc0*R_temp[_i][Cart::yy][m] + wmq0*R_temp[_i][Cart::yy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yy][m+1];
      R_temp[_i][Cart::xyz][m] = qmc0*R_temp[_i][Cart::yz][m] + wmq0*R_temp[_i][Cart::yz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yz][m+1];
      R_temp[_i][Cart::xzz][m] = qmc0*R_temp[_i][Cart::zz][m] + wmq0*R_temp[_i][Cart::zz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::zz][m+1];
      R_temp[_i][Cart::yyy][m] = qmc1*R_temp[_i][Cart::yy][m] + wmq1*R_temp[_i][Cart::yy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yy][m+1] + term_y;
      R_temp[_i][Cart::yyz][m] = qmc2*R_temp[_i][Cart::yy][m] + wmq2*R_temp[_i][Cart::yy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yy][m+1];
      R_temp[_i][Cart::yzz][m] = qmc1*R_temp[_i][Cart::zz][m] + wmq1*R_temp[_i][Cart::zz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::zz][m+1];
      R_temp[_i][Cart::zzz][m] = qmc2*R_temp[_i][Cart::zz][m] + wmq2*R_temp[_i][Cart::zz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::zz][m+1] + term_z;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma_delta > 2)


//Integrals     s-s - g-s     p-s - g-s     d-s - g-s     f-s - g-s     g-s - g-s     h-s - g-s     i-s - g-s     j-s - g-s     k-s - g-s     . . .
//              s-s - h-s     p-s - h-s     d-s - h-s     f-s - h-s     g-s - h-s     h-s - h-s     i-s - h-s     j-s - h-s     k-s - h-s     . . .
//              s-s - i-s     p-s - i-s     d-s - i-s     f-s - i-s     g-s - i-s     h-s - i-s     i-s - i-s     j-s - i-s     k-s - i-s     . . .
//                    j             j             j             j             j             j             j             j             j       . . .
//                    .             .             .             .             .             .             .             .             .       . . .
//                    .             .             .             .             .             .             .             .             .       . . .
for (int l = 4; l < _lmax_gamma_delta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int norb_3 = n_orbitals[l-3];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  int ncart_3 = ncart_2 + 1 - l;

  for (int m = 0; m < _lmax_gamma_delta+1-l; m++) {
    R_temp[0][norb_1][m] = qmc0*R_temp[0][norb_2][m] + wmq0*R_temp[0][norb_2][m+1] + (l-1)*reta*(R_temp[0][norb_3][m]-cfak*R_temp[0][norb_3][m+1]);
    R_temp[0][norb_1 + 1][m] = qmc1*R_temp[0][norb_2][m] + wmq1*R_temp[0][norb_2][m+1];
    R_temp[0][norb_1 + 2][m] = qmc2*R_temp[0][norb_2][m] + wmq2*R_temp[0][norb_2][m+1];
    int ntimes = 3;
    int itimes = 3;
    for (int k = 3; k < ncart_2; k++) {
      R_temp[0][norb_1 + k][m] = qmc0*R_temp[0][norb_2 + k][m] + wmq0*R_temp[0][norb_2 + k][m+1] + (l-ntimes)*reta*(R_temp[0][norb_3 + k][m]-cfak*R_temp[0][norb_3 + k][m+1]);
      itimes--;
      if (itimes == 0) {
        ntimes++;
        itimes = ntimes;
      }
    }
    for (int k = 0; k < l-1; k++) {
      R_temp[0][norb_1 + ncart_2 + k][m] = qmc0*R_temp[0][norb_2 + ncart_2 + k][m] + wmq0*R_temp[0][norb_2 + ncart_2 + k][m+1];
      R_temp[0][norb_1 + ncart_1 + k][m] = qmc1*R_temp[0][norb_2 + ncart_2 + k][m] + wmq1*R_temp[0][norb_2 + ncart_2 + k][m+1]
                                           + (l-1-k)*reta*(R_temp[0][norb_3 + ncart_3 + k][m]-cfak*R_temp[0][norb_3 + ncart_3 + k][m+1]);
    }
    R_temp[0][norb_1 + ncart_2 + l -1][m] = qmc0*R_temp[0][norb_2 + ncart_2 + l - 1][m] + wmq0*R_temp[0][norb_2 + ncart_2 + l - 1][m+1];
    R_temp[0][norb-2][m] = qmc1*R_temp[0][norb_1-1][m] + wmq1*R_temp[0][norb_1-1][m+1];
    R_temp[0][norb-1][m] = qmc2*R_temp[0][norb_1-1][m] + wmq2*R_temp[0][norb_1-1][m+1] + (l-1)*reta*(R_temp[0][norb_2-1][m]-cfak*R_temp[0][norb_2-1][m+1]);
  }

  for (int m = 0; m < _lmax_gamma_delta+1-l; m++) {
    for (int _i = 1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];

      R_temp[_i][norb_1][m] = qmc0*R_temp[_i][norb_2][m] + wmq0*R_temp[_i][norb_2][m+1] + nx_i*rdecay*R_temp[ilx_i][norb_2][m+1]
                              + (l-1)*reta*(R_temp[_i][norb_3][m]-cfak*R_temp[_i][norb_3][m+1]);
      R_temp[_i][norb_1 + 1][m] = qmc1*R_temp[_i][norb_2][m] + wmq1*R_temp[_i][norb_2][m+1] + ny_i*rdecay*R_temp[ily_i][norb_2][m+1];
      R_temp[_i][norb_1 + 2][m] = qmc2*R_temp[_i][norb_2][m] + wmq2*R_temp[_i][norb_2][m+1] + nz_i*rdecay*R_temp[ilz_i][norb_2][m+1];
      int ntimes = 3;
      int itimes = 3;
      for (int k = 3; k < ncart_2; k++) {
        R_temp[_i][norb_1 + k][m] = qmc0*R_temp[_i][norb_2 + k][m] + wmq0*R_temp[_i][norb_2 + k][m+1] + nx_i*rdecay*R_temp[ilx_i][norb_2 + k][m+1]
                                    + (l-ntimes)*reta*(R_temp[_i][norb_3 + k][m]-cfak*R_temp[_i][norb_3 + k][m+1]);
        itimes--;
        if (itimes == 0) {
          ntimes++;
          itimes = ntimes;
        }
      }
      for (int k = 0; k < l-1; k++) {
        int k2 = norb_2 + ncart_2 + k;
        R_temp[_i][norb_1 + ncart_2 + k][m] = qmc0*R_temp[_i][k2][m] + wmq0*R_temp[_i][k2][m+1] + nx_i*rdecay*R_temp[ilx_i][norb_2 + ncart_2 + k][m+1];
        R_temp[_i][norb_1 + ncart_1 + k][m] = qmc1*R_temp[_i][k2][m] + wmq1*R_temp[_i][k2][m+1] + ny_i*rdecay*R_temp[ily_i][norb_2 + ncart_2 + k][m+1]
                                              + (l-1-k)*reta*(R_temp[_i][norb_3 + ncart_3 + k][m]-cfak*R_temp[_i][norb_3 + ncart_3 + k][m+1]);
      }
      R_temp[_i][norb_1 + ncart_2 + l - 1][m] = qmc0*R_temp[_i][norb_2 + ncart_2 + l - 1][m] + wmq0*R_temp[_i][norb_2 + ncart_2 + l - 1][m+1]
                                                + nx_i*rdecay*R_temp[ilx_i][norb_2 + ncart_2 + l - 1][m+1];
      R_temp[_i][norb-2][m] = qmc1*R_temp[_i][norb_1-1][m] + wmq1*R_temp[_i][norb_1-1][m+1] + ny_i*rdecay*R_temp[ily_i][norb_1-1][m+1];
      R_temp[_i][norb-1][m] = qmc2*R_temp[_i][norb_1-1][m] + wmq2*R_temp[_i][norb_1-1][m+1] + nz_i*rdecay*R_temp[ilz_i][norb_1-1][m+1]
                              + (l-1)*reta*(R_temp[_i][norb_2-1][m]-cfak*R_temp[_i][norb_2-1][m+1]);
    }
  }

}
//------------------------------------------------------







//copy into new array for 3D use.

for (index3d i = 0; i < n_orbitals[_lmax_alpha_beta]; ++i) {
  for (index3d k = 0; k < n_orbitals[_lmax_gamma_delta]; ++k) {

    R[i][0][k] = R_temp[i][k][0];

  }
}



if (_lmax_beta > 0) {
  //Integrals     s-p - *-s     p-p - *-s     d-p - *-s     f-p - *-s     g-p - *-s     h-p - *-s     i-p - *-s     j-p - *-s     . . .
  for (int _i = 0; _i < n_orbitals[_lmax_alpha_beta-1]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < n_orbitals[_lmax_gamma_delta]; _j++) {
      R[_i][Cart::x][_j] = R[imx_i][0][_j] + amb0*R[_i][0][_j];
      R[_i][Cart::y][_j] = R[imy_i][0][_j] + amb1*R[_i][0][_j];
      R[_i][Cart::z][_j] = R[imz_i][0][_j] + amb2*R[_i][0][_j];
    }
  }
  //------------------------------------------------------
}

if (_lmax_beta > 1) {
  //Integrals     s-d - *-s     p-d - *-s     d-d - *-s     f-d - *-s     g-d - *-s     h-d - *-s     i-d - *-s     . . .
  for (int _i = 0; _i < n_orbitals[_lmax_alpha_beta-2]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < n_orbitals[_lmax_gamma_delta]; _j++) {
      R[_i][Cart::xx][_j] = R[imx_i][Cart::x][_j] + amb0*R[_i][Cart::x][_j];
      R[_i][Cart::xy][_j] = R[imx_i][Cart::y][_j] + amb0*R[_i][Cart::y][_j];
      R[_i][Cart::xz][_j] = R[imx_i][Cart::z][_j] + amb0*R[_i][Cart::z][_j];
      R[_i][Cart::yy][_j] = R[imy_i][Cart::y][_j] + amb1*R[_i][Cart::y][_j];
      R[_i][Cart::yz][_j] = R[imy_i][Cart::z][_j] + amb1*R[_i][Cart::z][_j];
      R[_i][Cart::zz][_j] = R[imz_i][Cart::z][_j] + amb2*R[_i][Cart::z][_j];
    }
  }
  //------------------------------------------------------
}

if (_lmax_beta > 2) {
  //Integrals     s-f - *-s     p-f - *-s     d-f - *-s     f-f - *-s     g-f - *-s     h-f - *-s     . . .
  for (int _i = 0; _i < n_orbitals[_lmax_alpha_beta-3]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < n_orbitals[_lmax_gamma_delta]; _j++) {
      R[_i][Cart::xxx][_j] = R[imx_i][Cart::xx][_j] + amb0*R[_i][Cart::xx][_j];
      R[_i][Cart::xxy][_j] = R[imx_i][Cart::xy][_j] + amb0*R[_i][Cart::xy][_j];
      R[_i][Cart::xxz][_j] = R[imx_i][Cart::xz][_j] + amb0*R[_i][Cart::xz][_j];
      R[_i][Cart::xyy][_j] = R[imx_i][Cart::yy][_j] + amb0*R[_i][Cart::yy][_j];
      R[_i][Cart::xyz][_j] = R[imx_i][Cart::yz][_j] + amb0*R[_i][Cart::yz][_j];
      R[_i][Cart::xzz][_j] = R[imx_i][Cart::zz][_j] + amb0*R[_i][Cart::zz][_j];
      R[_i][Cart::yyy][_j] = R[imy_i][Cart::yy][_j] + amb1*R[_i][Cart::yy][_j];
      R[_i][Cart::yyz][_j] = R[imy_i][Cart::yz][_j] + amb1*R[_i][Cart::yz][_j];
      R[_i][Cart::yzz][_j] = R[imy_i][Cart::zz][_j] + amb1*R[_i][Cart::zz][_j];
      R[_i][Cart::zzz][_j] = R[imz_i][Cart::zz][_j] + amb2*R[_i][Cart::zz][_j];
    }
  }
  //------------------------------------------------------
}




//Integrals     s-g - *-s     p-g - *-s     d-g - *-s     f-g - *-s     g-g - *-s     . . .
//              s-h - *-s     p-h - *-s     d-h - *-s     f-h - *-s     . . .
//              s-i - *-s     p-i - *-s     d-i - *-s     . . .
for (int l = 4; l < _lmax_beta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  for (int _i = 0; _i < n_orbitals[_lmax_alpha_beta-l]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < n_orbitals[_lmax_gamma_delta]; _j++) {
      for (int k = 0; k < ncart_2; k++) {
        R[_i][norb_1 + k][_j] = R[imx_i][norb_2 + k][_j] + amb0*R[_i][norb_2 + k][_j];
      }
      for (int k = 0; k < l; k++) {
        int k2 = norb_2 + ncart_2 + k;
        R[_i][norb_1 + ncart_2 + k][_j] = R[imx_i][k2][_j] + amb0*R[_i][k2][_j];
        R[_i][norb_1 + ncart_1 + k][_j] = R[imy_i][k2][_j] + amb1*R[_i][k2][_j];
      }
      R[_i][norb-1][_j] = R[imz_i][norb_1-1][_j] + amb2*R[_i][norb_1-1][_j];
    }
  }
}
//------------------------------------------------------






// Transforming alpha and beta functions to sphericals

            int istart[] = {0, 1, 1, 1, 4, 4, 4, 4, 4, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20 };
            int istop[] =  {0, 3, 3, 3, 9, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 34, 34, 34, 34, 34, 34, 34, 34, 34 };

            int _offset_alpha = _shell_alpha->getOffset();
            int _offset_beta = _shell_beta->getOffset();

            // prepare transformation matrices
            int _ntrafo_alpha = _shell_alpha->getNumFunc() + _offset_alpha;
            int _ntrafo_beta = _shell_beta->getNumFunc() + _offset_beta;

            
            // get transformation matrices
            const Eigen::MatrixXd _trafo_alpha = AOSuperMatrix::getTrafo(*italpha);
            const Eigen::MatrixXd _trafo_beta = AOSuperMatrix::getTrafo(*itbeta);

            tensor3d R3_ab_sph;
            R3_ab_sph.resize(extents[ _ntrafo_alpha ][ _ntrafo_beta ][ _ncombined_cd ]);

            for (int _i_beta = 0; _i_beta < _ntrafo_beta; _i_beta++) {
              for (int _i_alpha = 0; _i_alpha < _ntrafo_alpha; _i_alpha++) {

                for (int _j = 0; _j < n_orbitals[_lmax_gamma_delta]; _j++) {

                  R3_ab_sph[ _i_alpha ][ _i_beta ][ _j ] = 0.0;

                  for (int _i_beta_t = istart[ _i_beta ]; _i_beta_t <= istop[ _i_beta ]; _i_beta_t++) {
                    for (int _i_alpha_t = istart[ _i_alpha ]; _i_alpha_t <= istop[ _i_alpha ]; _i_alpha_t++) {

                      R3_ab_sph[ _i_alpha ][ _i_beta ][ _j ] += R[ _i_alpha_t ][ _i_beta_t][ _j]
                                                                * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_beta(_i_beta, _i_beta_t);


                    }
                  }

                }

              }
            }







//copy into new 4D array.
tensor4d R4_ab_sph;
R4_ab_sph.resize(extents4[ _ntrafo_alpha ][ _ntrafo_beta ][ _ncombined_cd ][ _ndelta ]);
//ma4_type R4_ab_sph(boost::extents[ _ntrafo_alpha ][ _ntrafo_beta ][ _ncombined_cd ][ _ndelta ]);

for (index3d j = 0; j < _ntrafo_alpha; ++j) {
  for (index3d k = 0; k < _ntrafo_beta; ++k) {
    for (index3d i = 0; i < _ncombined_cd; ++i) {

      R4_ab_sph[j][k][i][0] = R3_ab_sph[j][k][i];

    }
  }
}




if (_lmax_delta > 0) {
  //Integrals     *-* - s-p     *-* - p-p     *-* - d-p     *-* - f-p     *-* - g-p     *-* - h-p     *-* - i-p     *-* - j-p     . . .
  for (int _i = 0; _i < n_orbitals[_lmax_gamma_delta-1]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < _ntrafo_alpha; _j++) {
      for (int _k= 0; _k < _ntrafo_beta; _k++) {
        R4_ab_sph[_j][_k][_i][Cart::x] = R4_ab_sph[_j][_k][imx_i][0] + cmd0*R4_ab_sph[_j][_k][_i][0];
        R4_ab_sph[_j][_k][_i][Cart::y] = R4_ab_sph[_j][_k][imy_i][0] + cmd1*R4_ab_sph[_j][_k][_i][0];
        R4_ab_sph[_j][_k][_i][Cart::z] = R4_ab_sph[_j][_k][imz_i][0] + cmd2*R4_ab_sph[_j][_k][_i][0];
      }
    }
  }
  //------------------------------------------------------
}

if (_lmax_delta > 1) {
  //Integrals     *-* - s-d     *-* - p-d     *-* - d-d     *-* - f-d     *-* - g-d     *-* - h-d     *-* - i-d     . . .
  for (int _i = 0; _i < n_orbitals[_lmax_gamma_delta-2]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < _ntrafo_alpha; _j++) {
      for (int _k= 0; _k < _ntrafo_beta; _k++) {
        R4_ab_sph[_j][_k][_i][Cart::xx] = R4_ab_sph[_j][_k][imx_i][Cart::x] + cmd0*R4_ab_sph[_j][_k][_i][Cart::x];
        R4_ab_sph[_j][_k][_i][Cart::xy] = R4_ab_sph[_j][_k][imx_i][Cart::y] + cmd0*R4_ab_sph[_j][_k][_i][Cart::y];
        R4_ab_sph[_j][_k][_i][Cart::xz] = R4_ab_sph[_j][_k][imx_i][Cart::z] + cmd0*R4_ab_sph[_j][_k][_i][Cart::z];
        R4_ab_sph[_j][_k][_i][Cart::yy] = R4_ab_sph[_j][_k][imy_i][Cart::y] + cmd1*R4_ab_sph[_j][_k][_i][Cart::y];
        R4_ab_sph[_j][_k][_i][Cart::yz] = R4_ab_sph[_j][_k][imy_i][Cart::z] + cmd1*R4_ab_sph[_j][_k][_i][Cart::z];
        R4_ab_sph[_j][_k][_i][Cart::zz] = R4_ab_sph[_j][_k][imz_i][Cart::z] + cmd2*R4_ab_sph[_j][_k][_i][Cart::z];
      }
    }
  }
  //------------------------------------------------------
}

if (_lmax_delta > 2) {
  //Integrals     *-* - s-f     *-* - p-f     *-* - d-f     *-* - f-f     *-* - g-f     *-* - h-f     . . .
  for (int _i = 0; _i < n_orbitals[_lmax_gamma_delta-3]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < _ntrafo_alpha; _j++) {
      for (int _k= 0; _k < _ntrafo_beta; _k++) {
        R4_ab_sph[_j][_k][_i][Cart::xxx] = R4_ab_sph[_j][_k][imx_i][Cart::xx] + cmd0*R4_ab_sph[_j][_k][_i][Cart::xx];
        R4_ab_sph[_j][_k][_i][Cart::xxy] = R4_ab_sph[_j][_k][imx_i][Cart::xy] + cmd0*R4_ab_sph[_j][_k][_i][Cart::xy];
        R4_ab_sph[_j][_k][_i][Cart::xxz] = R4_ab_sph[_j][_k][imx_i][Cart::xz] + cmd0*R4_ab_sph[_j][_k][_i][Cart::xz];
        R4_ab_sph[_j][_k][_i][Cart::xyy] = R4_ab_sph[_j][_k][imx_i][Cart::yy] + cmd0*R4_ab_sph[_j][_k][_i][Cart::yy];
        R4_ab_sph[_j][_k][_i][Cart::xyz] = R4_ab_sph[_j][_k][imx_i][Cart::yz] + cmd0*R4_ab_sph[_j][_k][_i][Cart::yz];
        R4_ab_sph[_j][_k][_i][Cart::xzz] = R4_ab_sph[_j][_k][imx_i][Cart::zz] + cmd0*R4_ab_sph[_j][_k][_i][Cart::zz];
        R4_ab_sph[_j][_k][_i][Cart::yyy] = R4_ab_sph[_j][_k][imy_i][Cart::yy] + cmd1*R4_ab_sph[_j][_k][_i][Cart::yy];
        R4_ab_sph[_j][_k][_i][Cart::yyz] = R4_ab_sph[_j][_k][imy_i][Cart::yz] + cmd1*R4_ab_sph[_j][_k][_i][Cart::yz];
        R4_ab_sph[_j][_k][_i][Cart::yzz] = R4_ab_sph[_j][_k][imy_i][Cart::zz] + cmd1*R4_ab_sph[_j][_k][_i][Cart::zz];
        R4_ab_sph[_j][_k][_i][Cart::zzz] = R4_ab_sph[_j][_k][imz_i][Cart::zz] + cmd2*R4_ab_sph[_j][_k][_i][Cart::zz];
      }
    }
  }
  //------------------------------------------------------
}

//Integrals     *-* - s-g     *-* - p-g     *-* - d-g     *-* - f-g     *-* - g-g     . . .
//              *-* - s-h     *-* - p-h     *-* - d-h     *-* - f-h     . . .
//              *-* - s-i     *-* - p-i     *-* - d-i     . . .
for (int l = 4; l < _lmax_delta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  for (int _i = 0; _i < n_orbitals[_lmax_gamma_delta-l]; _i++) {
    int imx_i = i_more_x[_i];
    int imy_i = i_more_y[_i];
    int imz_i = i_more_z[_i];
    for (int _j = 0; _j < _ntrafo_alpha; _j++) {
      for (int _k= 0; _k < _ntrafo_beta; _k++) {
        for (int k = 0; k < ncart_2; k++) {
          R4_ab_sph[_j][_k][_i][norb_1 + k] = R4_ab_sph[_j][_k][imx_i][norb_2 + k] + cmd0*R4_ab_sph[_j][_k][_i][norb_2 + k];
        }
        for (int k = 0; k < l; k++) {
          int k2 = norb_2 + ncart_2 + k;
          R4_ab_sph[_j][_k][_i][norb_1 + ncart_2 + k] = R4_ab_sph[_j][_k][imx_i][k2] + cmd0*R4_ab_sph[_j][_k][_i][k2];
          R4_ab_sph[_j][_k][_i][norb_1 + ncart_1 + k] = R4_ab_sph[_j][_k][imy_i][k2] + cmd1*R4_ab_sph[_j][_k][_i][k2];
        }
        R4_ab_sph[_j][_k][_i][norb-1] = R4_ab_sph[_j][_k][imz_i][norb_1-1] + cmd2*R4_ab_sph[_j][_k][_i][norb_1-1];
      }
    }
  }
}
//------------------------------------------------------




// Transforming gamma and delta functions to sphericals

            int _offset_gamma = _shell_gamma->getOffset();
            int _offset_delta = _shell_delta->getOffset();

            // prepare transformation matrices
            int _ntrafo_gamma = _shell_gamma->getNumFunc() + _offset_gamma;
            int _ntrafo_delta = _shell_delta->getNumFunc() + _offset_delta;

            const Eigen::MatrixXd _trafo_gamma = AOSuperMatrix::getTrafo(*itgamma);
            const Eigen::MatrixXd _trafo_delta = AOSuperMatrix::getTrafo(*itdelta);


            tensor4d R4_sph;
            R4_sph.resize(extents4[ _ntrafo_alpha ][ _ntrafo_beta ][ _ntrafo_gamma ][ _ntrafo_delta ]);
            
            for (int _j = 0; _j < _ntrafo_alpha; _j++) {
                  for (int _k = 0; _k < _ntrafo_beta; _k++) {
                        for (int _i_gamma = 0; _i_gamma < _ntrafo_gamma; _i_gamma++) {
                            for (int _i_delta = 0; _i_delta < _ntrafo_delta; _i_delta++) {
             

                

                    R4_sph[ _j ][ _k ][ _i_gamma ][ _i_delta ] = 0.0;

                    for (int _i_delta_t = istart[ _i_delta ]; _i_delta_t <= istop[ _i_delta ]; _i_delta_t++) {
                      for (int _i_gamma_t = istart[ _i_gamma ]; _i_gamma_t <= istop[ _i_gamma ]; _i_gamma_t++) {

                        R4_sph[ _j ][ _k ][ _i_gamma ][ _i_delta ] += R4_ab_sph[_j][_k][ _i_gamma_t ][ _i_delta_t]
                                                                      * _trafo_gamma(_i_gamma, _i_gamma_t) * _trafo_delta(_i_delta, _i_delta_t);

                      }
                    }

                  }
                }

              }
            }

            
            

            int NumFunc_alpha = _shell_alpha->getNumFunc();
            int NumFunc_beta = _shell_beta->getNumFunc();
            int NumFunc_gamma = _shell_gamma->getNumFunc();
            int NumFunc_delta = _shell_delta->getNumFunc();
            
            
            for (int i_alpha = 0; i_alpha < NumFunc_alpha; i_alpha++) {
                for (int i_beta = 0; i_beta < NumFunc_beta; i_beta++) {
                  int a=i_alpha;
                  int b=i_beta;
                  if (alphabetaswitch) {
                    a=i_beta;
                    b=i_alpha;
                  }

                  for (int i_gamma = 0; i_gamma < NumFunc_gamma; i_gamma++) {
                    for (int i_delta = 0; i_delta < NumFunc_delta; i_delta++) {
                      int c=i_gamma;
                      int d=i_delta;
                      if (gammadeltaswitch) {
                        c=i_delta;
                        d=i_gamma;
                      }
                      if (ab_cd_switch) {
                        block[c][d][a][b] += R4_sph[_offset_alpha + i_alpha][_offset_beta + i_beta][_offset_gamma + i_gamma][_offset_delta + i_delta];
                      } else {
                        block[a][b][c][d] += R4_sph[_offset_alpha + i_alpha][_offset_beta + i_beta][_offset_gamma + i_gamma][_offset_delta + i_delta];
                      }
                    }
                  }
                }
              }

           
                 } // GaussianIterator itdelta
              } // GaussianIterator itgamma
           } // GaussianIterator itbeta
        } // GaussianIterator italpha

                        
    
       return _does_contribute;     
    } // TCrawMatrix::FillFourCenterRepBlock




    }}
