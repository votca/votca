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
 * distributed under the License is distributed on an "AS iS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/fourcenter.h>
#include <votca/xtp/aomatrix.h>


namespace votca {
    namespace xtp {

        /*
         * Calculate 4-center electron repulsion integrals 
         *    R_{abcd} = int{ phi_a(r) phi_b(r) phi_c(r') phi_d(r') /(r-r') d3rd3r' }
         * for a given set of a b c d as in 
         *    Obara, Saika, J. Chem. Phys. 84, 3963 (1986)
         * section iI.B for cartesian Gaussians, then transforming
         * to spherical (angular momentum) Gaussians ("complete" shells
         * from S to Lmax, and finally cutting out those angular momentum 
         * components actually present in shell-shell-shell combination.
         * Currently supported for 
         *      S,P,D,F,G  functions in DFT basis
         * 
         */

      
        bool FCMatrix::FillFourCenterRepBlock(tensor4d& block,
                const AOShell& shell_1, const AOShell& shell_2, const AOShell& shell_3, const AOShell& shell_4) {

            const double pi = boost::math::constants::pi<double>();
            
            bool does_contribute=true;
            

            // shell info, only lmax tells how far to go
            
            int lmax_1 = shell_1.getLmax();
            int lmax_2 = shell_2.getLmax();
            int lmax_3 = shell_3.getLmax();
            int lmax_4 = shell_4.getLmax();
            
            int mmax = lmax_1 + lmax_2 + lmax_3 + lmax_4;


            // set size of internal block for recursion
           
            const AOShell* shell_alpha;
            const AOShell* shell_beta;
            const AOShell* shell_gamma;
            const AOShell* shell_delta;
            bool alphabetaswitch = false;
            bool gammadeltaswitch = false;
            bool ab_cd_switch = false;

            // We need lmax_alpha > lmax_beta, so instead of calculating (sp,s) we calculate (ps,s), due to symmetry they are the same. 

            if (lmax_1 < lmax_2) {
              alphabetaswitch = true;               
            }
            if (lmax_3 < lmax_4) {
              gammadeltaswitch = true;               
            }
            if ( (lmax_1 + lmax_2) < (lmax_3 + lmax_4) ) {
              ab_cd_switch = true;               
            }

            if (ab_cd_switch==true) {

              if (alphabetaswitch== true) {
                shell_alpha = &shell_4;
                shell_beta = &shell_3;
              } else {
                shell_alpha = &shell_3;
                shell_beta = &shell_4;
              }

              if (gammadeltaswitch==true) {
                shell_gamma = &shell_2;
                shell_delta = &shell_1;
              } else {
                shell_gamma = &shell_1;
                shell_delta = &shell_2;
              }

            } else {

              if (alphabetaswitch== true) {
                shell_alpha = &shell_2;
                shell_beta = &shell_1;
              } else {
                shell_alpha = &shell_1;
                shell_beta = &shell_2;
              }

              if (gammadeltaswitch==true) {
                shell_gamma = &shell_4;
                shell_delta = &shell_3;
              } else {
                shell_gamma = &shell_3;
                shell_delta = &shell_4;
              }

            } 

            const Eigen::Vector3d& pos_alpha = shell_alpha->getPos();
            const Eigen::Vector3d& pos_beta = shell_beta->getPos();
            const Eigen::Vector3d& pos_gamma = shell_gamma->getPos();
            const Eigen::Vector3d& pos_delta = shell_delta->getPos();
            
            int lmax_alpha = shell_alpha->getLmax();
            int lmax_beta  = shell_beta->getLmax();
            int lmax_gamma = shell_gamma->getLmax();
            int lmax_delta = shell_delta->getLmax();

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

            int nbeta = AOSuperMatrix::getBlockSize(lmax_beta);
            int ndelta = AOSuperMatrix::getBlockSize(lmax_delta);
            int ncombined_ab = AOSuperMatrix::getBlockSize(lmax_alpha+lmax_beta);
            int ncombined_cd = AOSuperMatrix::getBlockSize(lmax_gamma+lmax_delta);
            
            tensor3d::extent_gen extents;
            tensor4d::extent_gen extents4;

            double dist_AB = (pos_alpha - pos_beta).squaredNorm();
            double dist_CD = (pos_gamma - pos_delta).squaredNorm();
            
            Eigen::Vector3d amb = pos_alpha - pos_beta;

            Eigen::Vector3d cmd = pos_gamma - pos_delta;
            
            
            for ( const auto& gaussian_alpha:*shell_alpha) {
                const double decay_alpha = gaussian_alpha.getDecay();
            
              for ( const auto& gaussian_beta:*shell_beta) {
                  const double decay_beta = gaussian_beta.getDecay();
                    
                for ( const auto& gaussian_gamma:*shell_gamma) {
                    const double decay_gamma = gaussian_gamma.getDecay();

                  for ( const auto& gaussian_delta:*shell_delta) {
                      const double decay_delta = gaussian_delta.getDecay();

            double zeta = decay_alpha + decay_beta;    
            double eta = decay_gamma + decay_delta;     
            double decay = zeta + eta;
            double rho = (zeta*eta)/decay;
            double rzeta = 0.5/zeta;
            double reta = 0.5/eta;
            double rdecay = 0.5/decay;
            double gfak = eta/decay;
            double cfak = zeta/decay;         
            Eigen::Vector3d P = (decay_alpha*pos_alpha + decay_beta*pos_beta)/zeta;
            Eigen::Vector3d Q = (decay_gamma*pos_gamma + decay_delta*pos_delta)/eta;
            Eigen::Vector3d W = (zeta*P + eta*Q)/decay;
            double U = rho*(P-Q).squaredNorm();

            Eigen::Vector3d pma = P - pos_alpha;
            Eigen::Vector3d qmc = Q - pos_gamma;
            Eigen::Vector3d wmp = W - P;
            Eigen::Vector3d wmq = W - Q;

            tensor3d R_temp;
            R_temp.resize(extents[ range(0, ncombined_ab ) ][ range(0, ncombined_cd ) ][ range(0, mmax+1)]);
            //initialize to zero
            for (index3d i = 0; i != ncombined_ab; ++i) {
              for (index3d j = 0; j != ncombined_cd; ++j) {
                for (index3d k = 0; k != mmax+1; ++k) {
                  R_temp[i][j][k] = 0.0;
                }
              }
             }

            tensor3d R;
            R.resize(extents[ range(0, ncombined_ab ) ][ range(0, nbeta ) ][ range(0, ncombined_cd)]);
            //initialize to zero
            for (index3d i = 0; i != ncombined_ab; ++i) {
              for (index3d j = 0; j != nbeta; ++j) {
                for (index3d k = 0; k != ncombined_cd; ++k) {
                  R[i][j][k] = 0.0;
                }
              }
            }
            
            const std::vector<double> FmT=AOMatrix<double>::XIntegrate(mmax+1, U);

            double exp_AB = exp( -2. * decay_alpha * decay_beta * rzeta * dist_AB );
            double exp_CD = exp( -2.* decay_gamma * decay_delta * reta * dist_CD );
            double ssss = ( 16. * pow( decay_alpha * decay_beta * decay_gamma * decay_delta, .75 ) * exp_AB * exp_CD )
                         / ( zeta * eta * sqrt(pi * decay) );

            //ss integrals
            for (int i=0; i < mmax+1; i++){
              R_temp[0][0][i] = ssss * FmT[i];
            }

int lmax_alpha_beta = lmax_alpha + lmax_beta;
int lmax_gamma_delta = lmax_gamma + lmax_delta;

//Integrals     p-s - s-s
if (lmax_alpha_beta > 0) {
  for (int m = 0; m < mmax; m++) {
    R_temp[Cart::x][0][m] = pma(0)*R_temp[0][0][m] + wmp(0)*R_temp[0][0][m+1];
    R_temp[Cart::y][0][m] = pma(1)*R_temp[0][0][m] + wmp(1)*R_temp[0][0][m+1];
    R_temp[Cart::z][0][m] = pma(2)*R_temp[0][0][m] + wmp(2)*R_temp[0][0][m+1];
  }
}
//------------------------------------------------------

//Integrals     d-s - s-s
if (lmax_alpha_beta > 1) {
  for (int m = 0; m < mmax-1; m++) {
    double term = rzeta*(R_temp[0][0][m]-gfak*R_temp[0][0][m+1]);
    R_temp[Cart::xx][0][m] = pma(0)*R_temp[Cart::x][0][m] + wmp(0)*R_temp[Cart::x][0][m+1] + term;
    R_temp[Cart::xy][0][m] = pma(0)*R_temp[Cart::y][0][m] + wmp(0)*R_temp[Cart::y][0][m+1];
    R_temp[Cart::xz][0][m] = pma(0)*R_temp[Cart::z][0][m] + wmp(0)*R_temp[Cart::z][0][m+1];
    R_temp[Cart::yy][0][m] = pma(1)*R_temp[Cart::y][0][m] + wmp(1)*R_temp[Cart::y][0][m+1] + term;
    R_temp[Cart::yz][0][m] = pma(1)*R_temp[Cart::z][0][m] + wmp(1)*R_temp[Cart::z][0][m+1];
    R_temp[Cart::zz][0][m] = pma(2)*R_temp[Cart::z][0][m] + wmp(2)*R_temp[Cart::z][0][m+1] + term;
  }
}
//------------------------------------------------------

//Integrals     f-s - s-s
if (lmax_alpha_beta > 2) {
  for (int m = 0; m < mmax-2; m++) {
    R_temp[Cart::xxx][0][m] = pma(0)*R_temp[Cart::xx][0][m] + wmp(0)*R_temp[Cart::xx][0][m+1] + 2*rzeta*(R_temp[Cart::x][0][m]-gfak*R_temp[Cart::x][0][m+1]);
    R_temp[Cart::xxy][0][m] = pma(1)*R_temp[Cart::xx][0][m] + wmp(1)*R_temp[Cart::xx][0][m+1];
    R_temp[Cart::xxz][0][m] = pma(2)*R_temp[Cart::xx][0][m] + wmp(2)*R_temp[Cart::xx][0][m+1];
    R_temp[Cart::xyy][0][m] = pma(0)*R_temp[Cart::yy][0][m] + wmp(0)*R_temp[Cart::yy][0][m+1];
    R_temp[Cart::xyz][0][m] = pma(0)*R_temp[Cart::yz][0][m] + wmp(0)*R_temp[Cart::yz][0][m+1];
    R_temp[Cart::xzz][0][m] = pma(0)*R_temp[Cart::zz][0][m] + wmp(0)*R_temp[Cart::zz][0][m+1];
    R_temp[Cart::yyy][0][m] = pma(1)*R_temp[Cart::yy][0][m] + wmp(1)*R_temp[Cart::yy][0][m+1] + 2*rzeta*(R_temp[Cart::y][0][m]-gfak*R_temp[Cart::y][0][m+1]);
    R_temp[Cart::yyz][0][m] = pma(2)*R_temp[Cart::yy][0][m] + wmp(2)*R_temp[Cart::yy][0][m+1];
    R_temp[Cart::yzz][0][m] = pma(1)*R_temp[Cart::zz][0][m] + wmp(1)*R_temp[Cart::zz][0][m+1];
    R_temp[Cart::zzz][0][m] = pma(2)*R_temp[Cart::zz][0][m] + wmp(2)*R_temp[Cart::zz][0][m+1] + 2*rzeta*(R_temp[Cart::z][0][m]-gfak*R_temp[Cart::z][0][m+1]);
  }
}
//------------------------------------------------------

//Integrals     g-s - s-s     h-s - s-s     i-s - s-s     j-s - s-s     k-s - s-s     . . .
for (int l = 4; l < lmax_alpha_beta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int norb_3 = n_orbitals[l-3];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  int ncart_3 = ncart_2 + 1 - l;
  for (int m = 0; m < mmax+1-l; m++) {
    R_temp[norb_1][0][m] = pma(0)*R_temp[norb_2][0][m] + wmp(0)*R_temp[norb_2][0][m+1] + (l-1)*rzeta*(R_temp[norb_3][0][m]-gfak*R_temp[norb_3][0][m+1]);
    R_temp[norb_1 + 1][0][m] = pma(1)*R_temp[norb_2][0][m] + wmp(1)*R_temp[norb_2][0][m+1];
    R_temp[norb_1 + 2][0][m] = pma(2)*R_temp[norb_2][0][m] + wmp(2)*R_temp[norb_2][0][m+1];
    int ntimes = 3;
    int itimes = 3;
    for (int k = 3; k < ncart_2; k++) {
      R_temp[norb_1 + k][0][m] = pma(0)*R_temp[norb_2 + k][0][m] + wmp(0)*R_temp[norb_2 + k][0][m+1] + (l-ntimes)*rzeta*(R_temp[norb_3 + k][0][m]-gfak*R_temp[norb_3 + k][0][m+1]);
      itimes--;
      if (itimes == 0) {
        ntimes++;
        itimes = ntimes;
      }
    }
    for (int k = 0; k < l-1; k++) {
      int k2 = norb_2 + ncart_2 + k;
      R_temp[norb_1 + ncart_2 + k][0][m] = pma(0)*R_temp[k2][0][m] + wmp(0)*R_temp[k2][0][m+1];
      R_temp[norb_1 + ncart_1 + k][0][m] = pma(1)*R_temp[k2][0][m] + wmp(1)*R_temp[k2][0][m+1]
                                           + (l-1-k)*rzeta*(R_temp[norb_3 + ncart_3 + k][0][m]-gfak*R_temp[norb_3 + ncart_3 + k][0][m+1]);
    }
    R_temp[norb_1 + ncart_2 + l - 1][0][m] = pma(0)*R_temp[norb_2 + ncart_2 + l - 1][0][m] + wmp(0)*R_temp[norb_2 + ncart_2 + l - 1][0][m+1];
    R_temp[norb-2][0][m] = pma(1)*R_temp[norb_1-1][0][m] + wmp(1)*R_temp[norb_1-1][0][m+1];
    R_temp[norb-1][0][m] = pma(2)*R_temp[norb_1-1][0][m] + wmp(2)*R_temp[norb_1-1][0][m+1] + (l-1)*rzeta*(R_temp[norb_2-1][0][m]-gfak*R_temp[norb_2-1][0][m+1]);
  }
}
//------------------------------------------------------

if (lmax_gamma_delta > 0) {

  //Integrals     s-s - p-s
  for (int m = 0; m < lmax_gamma_delta; m++) {
    R_temp[0][Cart::x][m] = qmc(0)*R_temp[0][0][m] + wmq(0)*R_temp[0][0][m+1];
    R_temp[0][Cart::y][m] = qmc(1)*R_temp[0][0][m] + wmq(1)*R_temp[0][0][m+1];
    R_temp[0][Cart::z][m] = qmc(2)*R_temp[0][0][m] + wmq(2)*R_temp[0][0][m+1];
  }
  //------------------------------------------------------

  //Integrals     p-s - p-s
  if (lmax_alpha_beta > 0) {
    for (int m = 0; m < lmax_gamma_delta; m++) {
      double term = rdecay*R_temp[0][0][m+1];
      for (int i =  1; i < 4; i++) {
        R_temp[i][Cart::x][m] = qmc(0)*R_temp[i][0][m] + wmq(0)*R_temp[i][0][m+1] + nx[i]*term;
        R_temp[i][Cart::y][m] = qmc(1)*R_temp[i][0][m] + wmq(1)*R_temp[i][0][m+1] + ny[i]*term;
        R_temp[i][Cart::z][m] = qmc(2)*R_temp[i][0][m] + wmq(2)*R_temp[i][0][m+1] + nz[i]*term;
      }
    }
  }
  //------------------------------------------------------

  //Integrals     d-s - p-s     f-s - p-s     g-s - p-s     h-s - p-s     i-s - p-s     j-s - p-s     k-s - p-s     . . .
  for (int m = 0; m < lmax_gamma_delta; m++) {
    for (int i = 4; i < n_orbitals[lmax_alpha_beta]; i++) {
      R_temp[i][Cart::x][m] = qmc(0)*R_temp[i][0][m] + wmq(0)*R_temp[i][0][m+1] + nx[i]*rdecay*R_temp[i_less_x[i]][0][m+1];
      R_temp[i][Cart::y][m] = qmc(1)*R_temp[i][0][m] + wmq(1)*R_temp[i][0][m+1] + ny[i]*rdecay*R_temp[i_less_y[i]][0][m+1];
      R_temp[i][Cart::z][m] = qmc(2)*R_temp[i][0][m] + wmq(2)*R_temp[i][0][m+1] + nz[i]*rdecay*R_temp[i_less_z[i]][0][m+1];
    }
  }
  //------------------------------------------------------

} // end if (lmax_gamma_delta > 0)

if (lmax_gamma_delta > 1) {

  //Integrals     s-s - d-s
  for (int m = 0; m < lmax_gamma_delta-1; m++) {
    double term = reta*(R_temp[0][0][m]-cfak*R_temp[0][0][m+1]);
    R_temp[0][Cart::xx][m] = qmc(0)*R_temp[0][Cart::x][m] + wmq(0)*R_temp[0][Cart::x][m+1] + term;
    R_temp[0][Cart::xy][m] = qmc(0)*R_temp[0][Cart::y][m] + wmq(0)*R_temp[0][Cart::y][m+1];
    R_temp[0][Cart::xz][m] = qmc(0)*R_temp[0][Cart::z][m] + wmq(0)*R_temp[0][Cart::z][m+1];
    R_temp[0][Cart::yy][m] = qmc(1)*R_temp[0][Cart::y][m] + wmq(1)*R_temp[0][Cart::y][m+1] + term;
    R_temp[0][Cart::yz][m] = qmc(1)*R_temp[0][Cart::z][m] + wmq(1)*R_temp[0][Cart::z][m+1];
    R_temp[0][Cart::zz][m] = qmc(2)*R_temp[0][Cart::z][m] + wmq(2)*R_temp[0][Cart::z][m+1] + term;
  }
  //------------------------------------------------------

  //Integrals     p-s - d-s     d-s - d-s     f-s - d-s     g-s - d-s     h-s - d-s     i-s - d-s     j-s - d-s     k-s - d-s     . . .
  for (int m = 0; m < lmax_gamma_delta-1; m++) {
    for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
        double term = reta*(R_temp[i][0][m]-cfak*R_temp[i][0][m+1]);
        R_temp[i][Cart::xx][m] = qmc(0)*R_temp[i][Cart::x][m] + wmq(0)*R_temp[i][Cart::x][m+1] +nx[i]*rdecay*R_temp[i_less_x[i]][Cart::x][m+1] + term;
        R_temp[i][Cart::xy][m] = qmc(0)*R_temp[i][Cart::y][m] + wmq(0)*R_temp[i][Cart::y][m+1] +nx[i]*rdecay*R_temp[i_less_x[i]][Cart::y][m+1];
        R_temp[i][Cart::xz][m] = qmc(0)*R_temp[i][Cart::z][m] + wmq(0)*R_temp[i][Cart::z][m+1] +nx[i]*rdecay*R_temp[i_less_x[i]][Cart::z][m+1];
        R_temp[i][Cart::yy][m] = qmc(1)*R_temp[i][Cart::y][m] + wmq(1)*R_temp[i][Cart::y][m+1] +ny[i]*rdecay*R_temp[i_less_y[i]][Cart::y][m+1] + term;
        R_temp[i][Cart::yz][m] = qmc(1)*R_temp[i][Cart::z][m] + wmq(1)*R_temp[i][Cart::z][m+1] +ny[i]*rdecay*R_temp[i_less_y[i]][Cart::z][m+1];
        R_temp[i][Cart::zz][m] = qmc(2)*R_temp[i][Cart::z][m] + wmq(2)*R_temp[i][Cart::z][m+1] +nz[i]*rdecay*R_temp[i_less_z[i]][Cart::z][m+1] + term;
    }
  }
  //------------------------------------------------------
} // end if (lmax_gamma_delta > 1)

if (lmax_gamma_delta > 2) {

  //Integrals     s-s - f-s
  for (int m = 0; m < lmax_gamma_delta-2; m++) {
    R_temp[0][Cart::xxx][m] = qmc(0)*R_temp[0][Cart::xx][m] + wmq(0)*R_temp[0][Cart::xx][m+1] + 2*reta*(R_temp[0][Cart::x][m]-cfak*R_temp[0][Cart::x][m+1]);
    R_temp[0][Cart::xxy][m] = qmc(1)*R_temp[0][Cart::xx][m] + wmq(1)*R_temp[0][Cart::xx][m+1];
    R_temp[0][Cart::xxz][m] = qmc(2)*R_temp[0][Cart::xx][m] + wmq(2)*R_temp[0][Cart::xx][m+1];
    R_temp[0][Cart::xyy][m] = qmc(0)*R_temp[0][Cart::yy][m] + wmq(0)*R_temp[0][Cart::yy][m+1];
    R_temp[0][Cart::xyz][m] = qmc(0)*R_temp[0][Cart::yz][m] + wmq(0)*R_temp[0][Cart::yz][m+1];
    R_temp[0][Cart::xzz][m] = qmc(0)*R_temp[0][Cart::zz][m] + wmq(0)*R_temp[0][Cart::zz][m+1];
    R_temp[0][Cart::yyy][m] = qmc(1)*R_temp[0][Cart::yy][m] + wmq(1)*R_temp[0][Cart::yy][m+1] + 2*reta*(R_temp[0][Cart::y][m]-cfak*R_temp[0][Cart::y][m+1]);
    R_temp[0][Cart::yyz][m] = qmc(2)*R_temp[0][Cart::yy][m] + wmq(2)*R_temp[0][Cart::yy][m+1];
    R_temp[0][Cart::yzz][m] = qmc(1)*R_temp[0][Cart::zz][m] + wmq(1)*R_temp[0][Cart::zz][m+1];
    R_temp[0][Cart::zzz][m] = qmc(2)*R_temp[0][Cart::zz][m] + wmq(2)*R_temp[0][Cart::zz][m+1] + 2*reta*(R_temp[0][Cart::z][m]-cfak*R_temp[0][Cart::z][m+1]);
  }
  //------------------------------------------------------

  //Integrals     p-s - f-s     d-s - f-s     f-s - f-s     g-s - f-s     h-s - f-s     i-s - f-s     j-s - f-s     k-s - f-s     . . .
  for (int m = 0; m < lmax_gamma_delta-2; m++) {
    for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
      double term_x = 2*reta*(R_temp[i][Cart::x][m]-cfak*R_temp[i][Cart::x][m+1]);
      double term_y = 2*reta*(R_temp[i][Cart::y][m]-cfak*R_temp[i][Cart::y][m+1]);
      double term_z = 2*reta*(R_temp[i][Cart::z][m]-cfak*R_temp[i][Cart::z][m+1]);
      R_temp[i][Cart::xxx][m] = qmc(0)*R_temp[i][Cart::xx][m] + wmq(0)*R_temp[i][Cart::xx][m+1] + nx[i]*rdecay*R_temp[i_less_x[i]][Cart::xx][m+1] + term_x;
      R_temp[i][Cart::xxy][m] = qmc(1)*R_temp[i][Cart::xx][m] + wmq(1)*R_temp[i][Cart::xx][m+1] + ny[i]*rdecay*R_temp[i_less_y[i]][Cart::xx][m+1];
      R_temp[i][Cart::xxz][m] = qmc(2)*R_temp[i][Cart::xx][m] + wmq(2)*R_temp[i][Cart::xx][m+1] + nz[i]*rdecay*R_temp[i_less_z[i]][Cart::xx][m+1];
      R_temp[i][Cart::xyy][m] = qmc(0)*R_temp[i][Cart::yy][m] + wmq(0)*R_temp[i][Cart::yy][m+1] + nx[i]*rdecay*R_temp[i_less_x[i]][Cart::yy][m+1];
      R_temp[i][Cart::xyz][m] = qmc(0)*R_temp[i][Cart::yz][m] + wmq(0)*R_temp[i][Cart::yz][m+1] + nx[i]*rdecay*R_temp[i_less_x[i]][Cart::yz][m+1];
      R_temp[i][Cart::xzz][m] = qmc(0)*R_temp[i][Cart::zz][m] + wmq(0)*R_temp[i][Cart::zz][m+1] + nx[i]*rdecay*R_temp[i_less_x[i]][Cart::zz][m+1];
      R_temp[i][Cart::yyy][m] = qmc(1)*R_temp[i][Cart::yy][m] + wmq(1)*R_temp[i][Cart::yy][m+1] + ny[i]*rdecay*R_temp[i_less_y[i]][Cart::yy][m+1] + term_y;
      R_temp[i][Cart::yyz][m] = qmc(2)*R_temp[i][Cart::yy][m] + wmq(2)*R_temp[i][Cart::yy][m+1] + nz[i]*rdecay*R_temp[i_less_z[i]][Cart::yy][m+1];
      R_temp[i][Cart::yzz][m] = qmc(1)*R_temp[i][Cart::zz][m] + wmq(1)*R_temp[i][Cart::zz][m+1] + ny[i]*rdecay*R_temp[i_less_y[i]][Cart::zz][m+1];
      R_temp[i][Cart::zzz][m] = qmc(2)*R_temp[i][Cart::zz][m] + wmq(2)*R_temp[i][Cart::zz][m+1] + nz[i]*rdecay*R_temp[i_less_z[i]][Cart::zz][m+1] + term_z;
    }
  }
  //------------------------------------------------------
} // end if (lmax_gamma_delta > 2)

//Integrals     s-s - g-s     p-s - g-s     d-s - g-s     f-s - g-s     g-s - g-s     h-s - g-s     i-s - g-s     j-s - g-s     k-s - g-s     . . .
//              s-s - h-s     p-s - h-s     d-s - h-s     f-s - h-s     g-s - h-s     h-s - h-s     i-s - h-s     j-s - h-s     k-s - h-s     . . .
//              s-s - i-s     p-s - i-s     d-s - i-s     f-s - i-s     g-s - i-s     h-s - i-s     i-s - i-s     j-s - i-s     k-s - i-s     . . .
//                    j             j             j             j             j             j             j             j             j       . . .
//                    .             .             .             .             .             .             .             .             .       . . .
//                    .             .             .             .             .             .             .             .             .       . . .
for (int l = 4; l < lmax_gamma_delta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int norb_3 = n_orbitals[l-3];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  int ncart_3 = ncart_2 + 1 - l;

  for (int m = 0; m < lmax_gamma_delta+1-l; m++) {
    R_temp[0][norb_1][m] = qmc(0)*R_temp[0][norb_2][m] + wmq(0)*R_temp[0][norb_2][m+1] + (l-1)*reta*(R_temp[0][norb_3][m]-cfak*R_temp[0][norb_3][m+1]);
    R_temp[0][norb_1 + 1][m] = qmc(1)*R_temp[0][norb_2][m] + wmq(1)*R_temp[0][norb_2][m+1];
    R_temp[0][norb_1 + 2][m] = qmc(2)*R_temp[0][norb_2][m] + wmq(2)*R_temp[0][norb_2][m+1];
    int ntimes = 3;
    int itimes = 3;
    for (int k = 3; k < ncart_2; k++) {
      R_temp[0][norb_1 + k][m] = qmc(0)*R_temp[0][norb_2 + k][m] + wmq(0)*R_temp[0][norb_2 + k][m+1] + (l-ntimes)*reta*(R_temp[0][norb_3 + k][m]-cfak*R_temp[0][norb_3 + k][m+1]);
      itimes--;
      if (itimes == 0) {
        ntimes++;
        itimes = ntimes;
      }
    }
    for (int k = 0; k < l-1; k++) {
      R_temp[0][norb_1 + ncart_2 + k][m] = qmc(0)*R_temp[0][norb_2 + ncart_2 + k][m] + wmq(0)*R_temp[0][norb_2 + ncart_2 + k][m+1];
      R_temp[0][norb_1 + ncart_1 + k][m] = qmc(1)*R_temp[0][norb_2 + ncart_2 + k][m] + wmq(1)*R_temp[0][norb_2 + ncart_2 + k][m+1]
                                           + (l-1-k)*reta*(R_temp[0][norb_3 + ncart_3 + k][m]-cfak*R_temp[0][norb_3 + ncart_3 + k][m+1]);
    }
    R_temp[0][norb_1 + ncart_2 + l -1][m] = qmc(0)*R_temp[0][norb_2 + ncart_2 + l - 1][m] + wmq(0)*R_temp[0][norb_2 + ncart_2 + l - 1][m+1];
    R_temp[0][norb-2][m] = qmc(1)*R_temp[0][norb_1-1][m] + wmq(1)*R_temp[0][norb_1-1][m+1];
    R_temp[0][norb-1][m] = qmc(2)*R_temp[0][norb_1-1][m] + wmq(2)*R_temp[0][norb_1-1][m+1] + (l-1)*reta*(R_temp[0][norb_2-1][m]-cfak*R_temp[0][norb_2-1][m+1]);
  }

  for (int m = 0; m < lmax_gamma_delta+1-l; m++) {
    for (int i = 1; i < n_orbitals[lmax_alpha_beta]; i++) {
      int nx_i = nx[i];
      int ny_i = ny[i];
      int nz_i = nz[i];
      int ilx_i = i_less_x[i];
      int ily_i = i_less_y[i];
      int ilz_i = i_less_z[i];

      R_temp[i][norb_1][m] = qmc(0)*R_temp[i][norb_2][m] + wmq(0)*R_temp[i][norb_2][m+1] + nx_i*rdecay*R_temp[ilx_i][norb_2][m+1]
                              + (l-1)*reta*(R_temp[i][norb_3][m]-cfak*R_temp[i][norb_3][m+1]);
      R_temp[i][norb_1 + 1][m] = qmc(1)*R_temp[i][norb_2][m] + wmq(1)*R_temp[i][norb_2][m+1] + ny_i*rdecay*R_temp[ily_i][norb_2][m+1];
      R_temp[i][norb_1 + 2][m] = qmc(2)*R_temp[i][norb_2][m] + wmq(2)*R_temp[i][norb_2][m+1] + nz_i*rdecay*R_temp[ilz_i][norb_2][m+1];
      int ntimes = 3;
      int itimes = 3;
      for (int k = 3; k < ncart_2; k++) {
        R_temp[i][norb_1 + k][m] = qmc(0)*R_temp[i][norb_2 + k][m] + wmq(0)*R_temp[i][norb_2 + k][m+1] + nx_i*rdecay*R_temp[ilx_i][norb_2 + k][m+1]
                                    + (l-ntimes)*reta*(R_temp[i][norb_3 + k][m]-cfak*R_temp[i][norb_3 + k][m+1]);
        itimes--;
        if (itimes == 0) {
          ntimes++;
          itimes = ntimes;
        }
      }
      for (int k = 0; k < l-1; k++) {
        int k2 = norb_2 + ncart_2 + k;
        R_temp[i][norb_1 + ncart_2 + k][m] = qmc(0)*R_temp[i][k2][m] + wmq(0)*R_temp[i][k2][m+1] + nx_i*rdecay*R_temp[ilx_i][norb_2 + ncart_2 + k][m+1];
        R_temp[i][norb_1 + ncart_1 + k][m] = qmc(1)*R_temp[i][k2][m] + wmq(1)*R_temp[i][k2][m+1] + ny_i*rdecay*R_temp[ily_i][norb_2 + ncart_2 + k][m+1]
                                              + (l-1-k)*reta*(R_temp[i][norb_3 + ncart_3 + k][m]-cfak*R_temp[i][norb_3 + ncart_3 + k][m+1]);
      }
      R_temp[i][norb_1 + ncart_2 + l - 1][m] = qmc(0)*R_temp[i][norb_2 + ncart_2 + l - 1][m] + wmq(0)*R_temp[i][norb_2 + ncart_2 + l - 1][m+1]
                                                + nx_i*rdecay*R_temp[ilx_i][norb_2 + ncart_2 + l - 1][m+1];
      R_temp[i][norb-2][m] = qmc(1)*R_temp[i][norb_1-1][m] + wmq(1)*R_temp[i][norb_1-1][m+1] + ny_i*rdecay*R_temp[ily_i][norb_1-1][m+1];
      R_temp[i][norb-1][m] = qmc(2)*R_temp[i][norb_1-1][m] + wmq(2)*R_temp[i][norb_1-1][m+1] + nz_i*rdecay*R_temp[ilz_i][norb_1-1][m+1]
                              + (l-1)*reta*(R_temp[i][norb_2-1][m]-cfak*R_temp[i][norb_2-1][m+1]);
    }
  }

}
//------------------------------------------------------

//copy into new array for 3D use.

for (index3d i = 0; i < n_orbitals[lmax_alpha_beta]; ++i) {
  for (index3d k = 0; k < n_orbitals[lmax_gamma_delta]; ++k) {
    R[i][0][k] = R_temp[i][k][0];
  }
}

if (lmax_beta > 0) {
  //Integrals     s-p - *-s     p-p - *-s     d-p - *-s     f-p - *-s     g-p - *-s     h-p - *-s     i-p - *-s     j-p - *-s     . . .
  for (int i = 0; i < n_orbitals[lmax_alpha_beta-1]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < n_orbitals[lmax_gamma_delta]; j++) {
      R[i][Cart::x][j] = R[imx_i][0][j] + amb(0)*R[i][0][j];
      R[i][Cart::y][j] = R[imy_i][0][j] + amb(1)*R[i][0][j];
      R[i][Cart::z][j] = R[imz_i][0][j] + amb(2)*R[i][0][j];
    }
  }
  //------------------------------------------------------
}

if (lmax_beta > 1) {
  //Integrals     s-d - *-s     p-d - *-s     d-d - *-s     f-d - *-s     g-d - *-s     h-d - *-s     i-d - *-s     . . .
  for (int i = 0; i < n_orbitals[lmax_alpha_beta-2]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < n_orbitals[lmax_gamma_delta]; j++) {
      R[i][Cart::xx][j] = R[imx_i][Cart::x][j] + amb(0)*R[i][Cart::x][j];
      R[i][Cart::xy][j] = R[imx_i][Cart::y][j] + amb(0)*R[i][Cart::y][j];
      R[i][Cart::xz][j] = R[imx_i][Cart::z][j] + amb(0)*R[i][Cart::z][j];
      R[i][Cart::yy][j] = R[imy_i][Cart::y][j] + amb(1)*R[i][Cart::y][j];
      R[i][Cart::yz][j] = R[imy_i][Cart::z][j] + amb(1)*R[i][Cart::z][j];
      R[i][Cart::zz][j] = R[imz_i][Cart::z][j] + amb(2)*R[i][Cart::z][j];
    }
  }
  //------------------------------------------------------
}

if (lmax_beta > 2) {
  //Integrals     s-f - *-s     p-f - *-s     d-f - *-s     f-f - *-s     g-f - *-s     h-f - *-s     . . .
  for (int i = 0; i < n_orbitals[lmax_alpha_beta-3]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < n_orbitals[lmax_gamma_delta]; j++) {
      R[i][Cart::xxx][j] = R[imx_i][Cart::xx][j] + amb(0)*R[i][Cart::xx][j];
      R[i][Cart::xxy][j] = R[imx_i][Cart::xy][j] + amb(0)*R[i][Cart::xy][j];
      R[i][Cart::xxz][j] = R[imx_i][Cart::xz][j] + amb(0)*R[i][Cart::xz][j];
      R[i][Cart::xyy][j] = R[imx_i][Cart::yy][j] + amb(0)*R[i][Cart::yy][j];
      R[i][Cart::xyz][j] = R[imx_i][Cart::yz][j] + amb(0)*R[i][Cart::yz][j];
      R[i][Cart::xzz][j] = R[imx_i][Cart::zz][j] + amb(0)*R[i][Cart::zz][j];
      R[i][Cart::yyy][j] = R[imy_i][Cart::yy][j] + amb(1)*R[i][Cart::yy][j];
      R[i][Cart::yyz][j] = R[imy_i][Cart::yz][j] + amb(1)*R[i][Cart::yz][j];
      R[i][Cart::yzz][j] = R[imy_i][Cart::zz][j] + amb(1)*R[i][Cart::zz][j];
      R[i][Cart::zzz][j] = R[imz_i][Cart::zz][j] + amb(2)*R[i][Cart::zz][j];
    }
  }
  //------------------------------------------------------
}

//Integrals     s-g - *-s     p-g - *-s     d-g - *-s     f-g - *-s     g-g - *-s     . . .
//              s-h - *-s     p-h - *-s     d-h - *-s     f-h - *-s     . . .
//              s-i - *-s     p-i - *-s     d-i - *-s     . . .
for (int l = 4; l < lmax_beta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  for (int i = 0; i < n_orbitals[lmax_alpha_beta-l]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < n_orbitals[lmax_gamma_delta]; j++) {
      for (int k = 0; k < ncart_2; k++) {
        R[i][norb_1 + k][j] = R[imx_i][norb_2 + k][j] + amb(0)*R[i][norb_2 + k][j];
      }
      for (int k = 0; k < l; k++) {
        int k2 = norb_2 + ncart_2 + k;
        R[i][norb_1 + ncart_2 + k][j] = R[imx_i][k2][j] + amb(0)*R[i][k2][j];
        R[i][norb_1 + ncart_1 + k][j] = R[imy_i][k2][j] + amb(1)*R[i][k2][j];
      }
      R[i][norb-1][j] = R[imz_i][norb_1-1][j] + amb(2)*R[i][norb_1-1][j];
    }
  }
}
//------------------------------------------------------

// Transforming alpha and beta functions to sphericals

            int istart[] = {0, 1, 1, 1, 4, 4, 4, 4, 4, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20 };
            int istop[] =  {0, 3, 3, 3, 9, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 34, 34, 34, 34, 34, 34, 34, 34, 34 };

            int offset_alpha = shell_alpha->getOffset();
            int offset_beta = shell_beta->getOffset();

            // prepare transformation matrices
            int ntrafo_alpha = shell_alpha->getNumFunc() + offset_alpha;
            int ntrafo_beta = shell_beta->getNumFunc() + offset_beta;

            
            // get transformation matrices
            const Eigen::MatrixXd trafo_alpha = AOSuperMatrix::getTrafo(gaussian_alpha);
            const Eigen::MatrixXd trafo_beta = AOSuperMatrix::getTrafo(gaussian_beta);

            tensor3d R3_ab_sph;
            R3_ab_sph.resize(extents[ ntrafo_alpha ][ ntrafo_beta ][ ncombined_cd ]);

            for (int i_beta = 0; i_beta < ntrafo_beta; i_beta++) {
              for (int i_alpha = 0; i_alpha < ntrafo_alpha; i_alpha++) {

                for (int j = 0; j < n_orbitals[lmax_gamma_delta]; j++) {

                  R3_ab_sph[ i_alpha ][ i_beta ][ j ] = 0.0;

                  for (int i_beta_t = istart[ i_beta ]; i_beta_t <= istop[ i_beta ]; i_beta_t++) {
                    for (int i_alpha_t = istart[ i_alpha ]; i_alpha_t <= istop[ i_alpha ]; i_alpha_t++) {

                      R3_ab_sph[ i_alpha ][ i_beta ][ j ] += R[ i_alpha_t ][ i_beta_t][ j]
                                                                * trafo_alpha(i_alpha_t, i_alpha) * trafo_beta(i_beta_t, i_beta);
                    }
                  }
                }
              }
            }

//copy into new 4D array.
tensor4d R4_ab_sph;
R4_ab_sph.resize(extents4[ ntrafo_alpha ][ ntrafo_beta ][ ncombined_cd ][ ndelta ]);
//ma4_type R4_ab_sph(boost::extents[ _ntrafo_alpha ][ _ntrafo_beta ][ _ncombined_cd ][ _ndelta ]);

for (index3d j = 0; j < ntrafo_alpha; ++j) {
  for (index3d k = 0; k < ntrafo_beta; ++k) {
    for (index3d i = 0; i < ncombined_cd; ++i) {

      R4_ab_sph[j][k][i][0] = R3_ab_sph[j][k][i];

    }
  }
}

if (lmax_delta > 0) {
  //Integrals     *-* - s-p     *-* - p-p     *-* - d-p     *-* - f-p     *-* - g-p     *-* - h-p     *-* - i-p     *-* - j-p     . . .
  for (int i = 0; i < n_orbitals[lmax_gamma_delta-1]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < ntrafo_alpha; j++) {
      for (int k= 0; k < ntrafo_beta; k++) {
        R4_ab_sph[j][k][i][Cart::x] = R4_ab_sph[j][k][imx_i][0] + cmd(0)*R4_ab_sph[j][k][i][0];
        R4_ab_sph[j][k][i][Cart::y] = R4_ab_sph[j][k][imy_i][0] + cmd(1)*R4_ab_sph[j][k][i][0];
        R4_ab_sph[j][k][i][Cart::z] = R4_ab_sph[j][k][imz_i][0] + cmd(2)*R4_ab_sph[j][k][i][0];
      }
    }
  }
}

if (lmax_delta > 1) {
  //Integrals     *-* - s-d     *-* - p-d     *-* - d-d     *-* - f-d     *-* - g-d     *-* - h-d     *-* - i-d     . . .
  for (int i = 0; i < n_orbitals[lmax_gamma_delta-2]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < ntrafo_alpha; j++) {
      for (int k= 0; k < ntrafo_beta; k++) {
        R4_ab_sph[j][k][i][Cart::xx] = R4_ab_sph[j][k][imx_i][Cart::x] + cmd(0)*R4_ab_sph[j][k][i][Cart::x];
        R4_ab_sph[j][k][i][Cart::xy] = R4_ab_sph[j][k][imx_i][Cart::y] + cmd(0)*R4_ab_sph[j][k][i][Cart::y];
        R4_ab_sph[j][k][i][Cart::xz] = R4_ab_sph[j][k][imx_i][Cart::z] + cmd(0)*R4_ab_sph[j][k][i][Cart::z];
        R4_ab_sph[j][k][i][Cart::yy] = R4_ab_sph[j][k][imy_i][Cart::y] + cmd(1)*R4_ab_sph[j][k][i][Cart::y];
        R4_ab_sph[j][k][i][Cart::yz] = R4_ab_sph[j][k][imy_i][Cart::z] + cmd(1)*R4_ab_sph[j][k][i][Cart::z];
        R4_ab_sph[j][k][i][Cart::zz] = R4_ab_sph[j][k][imz_i][Cart::z] + cmd(2)*R4_ab_sph[j][k][i][Cart::z];
      }
    }
  }
}

if (lmax_delta > 2) {
  //Integrals     *-* - s-f     *-* - p-f     *-* - d-f     *-* - f-f     *-* - g-f     *-* - h-f     . . .
  for (int i = 0; i < n_orbitals[lmax_gamma_delta-3]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < ntrafo_alpha; j++) {
      for (int k= 0; k < ntrafo_beta; k++) {
        R4_ab_sph[j][k][i][Cart::xxx] = R4_ab_sph[j][k][imx_i][Cart::xx] + cmd(0)*R4_ab_sph[j][k][i][Cart::xx];
        R4_ab_sph[j][k][i][Cart::xxy] = R4_ab_sph[j][k][imx_i][Cart::xy] + cmd(0)*R4_ab_sph[j][k][i][Cart::xy];
        R4_ab_sph[j][k][i][Cart::xxz] = R4_ab_sph[j][k][imx_i][Cart::xz] + cmd(0)*R4_ab_sph[j][k][i][Cart::xz];
        R4_ab_sph[j][k][i][Cart::xyy] = R4_ab_sph[j][k][imx_i][Cart::yy] + cmd(0)*R4_ab_sph[j][k][i][Cart::yy];
        R4_ab_sph[j][k][i][Cart::xyz] = R4_ab_sph[j][k][imx_i][Cart::yz] + cmd(0)*R4_ab_sph[j][k][i][Cart::yz];
        R4_ab_sph[j][k][i][Cart::xzz] = R4_ab_sph[j][k][imx_i][Cart::zz] + cmd(0)*R4_ab_sph[j][k][i][Cart::zz];
        R4_ab_sph[j][k][i][Cart::yyy] = R4_ab_sph[j][k][imy_i][Cart::yy] + cmd(1)*R4_ab_sph[j][k][i][Cart::yy];
        R4_ab_sph[j][k][i][Cart::yyz] = R4_ab_sph[j][k][imy_i][Cart::yz] + cmd(1)*R4_ab_sph[j][k][i][Cart::yz];
        R4_ab_sph[j][k][i][Cart::yzz] = R4_ab_sph[j][k][imy_i][Cart::zz] + cmd(1)*R4_ab_sph[j][k][i][Cart::zz];
        R4_ab_sph[j][k][i][Cart::zzz] = R4_ab_sph[j][k][imz_i][Cart::zz] + cmd(2)*R4_ab_sph[j][k][i][Cart::zz];
      }
    }
  }
}

//Integrals     *-* - s-g     *-* - p-g     *-* - d-g     *-* - f-g     *-* - g-g     . . .
//              *-* - s-h     *-* - p-h     *-* - d-h     *-* - f-h     . . .
//              *-* - s-i     *-* - p-i     *-* - d-i     . . .
for (int l = 4; l < lmax_delta+1; l++) {
  int norb = n_orbitals[l];
  int norb_1 = n_orbitals[l-1];
  int norb_2 = n_orbitals[l-2];
  int ncart_1 = (l*(l+1))/2;
  int ncart_2 = ncart_1 - l;
  for (int i = 0; i < n_orbitals[lmax_gamma_delta-l]; i++) {
    int imx_i = i_more_x[i];
    int imy_i = i_more_y[i];
    int imz_i = i_more_z[i];
    for (int j = 0; j < ntrafo_alpha; j++) {
      for (int k= 0; k < ntrafo_beta; k++) {
        for (int m = 0; m < ncart_2; m++) {
          R4_ab_sph[j][k][i][norb_1 + m] = R4_ab_sph[j][k][imx_i][norb_2 + m] + cmd(0)*R4_ab_sph[j][k][i][norb_2 + m];
        }
        for (int m = 0; m < l; m++) {
          int n = norb_2 + ncart_2 + m;
          R4_ab_sph[j][k][i][norb_1 + ncart_2 + m] = R4_ab_sph[j][k][imx_i][n] + cmd(0)*R4_ab_sph[j][k][i][n];
          R4_ab_sph[j][k][i][norb_1 + ncart_1 + m] = R4_ab_sph[j][k][imy_i][n] + cmd(1)*R4_ab_sph[j][k][i][n];
        }
        R4_ab_sph[j][k][i][norb-1] = R4_ab_sph[j][k][imz_i][norb_1-1] + cmd(2)*R4_ab_sph[j][k][i][norb_1-1];
      }
    }
  }
}

// Transforming gamma and delta functions to sphericals

            int offset_gamma = shell_gamma->getOffset();
            int offset_delta = shell_delta->getOffset();

            // prepare transformation matrices
            int ntrafo_gamma = shell_gamma->getNumFunc() + offset_gamma;
            int ntrafo_delta = shell_delta->getNumFunc() + offset_delta;

            const Eigen::MatrixXd trafo_gamma = AOSuperMatrix::getTrafo(gaussian_gamma);
            const Eigen::MatrixXd trafo_delta = AOSuperMatrix::getTrafo(gaussian_delta);

            tensor4d R4_sph;
            R4_sph.resize(extents4[ ntrafo_alpha ][ ntrafo_beta ][ ntrafo_gamma ][ ntrafo_delta ]);
            
            for (int j = 0; j < ntrafo_alpha; j++) {
                  for (int k = 0; k < ntrafo_beta; k++) {
                        for (int i_gamma = 0; i_gamma < ntrafo_gamma; i_gamma++) {
                            for (int i_delta = 0; i_delta < ntrafo_delta; i_delta++) {

                    R4_sph[ j ][ k ][ i_gamma ][ i_delta ] = 0.0;

                    for (int i_delta_t = istart[ i_delta ]; i_delta_t <= istop[ i_delta ]; i_delta_t++) {
                      for (int i_gamma_t = istart[ i_gamma ]; i_gamma_t <= istop[ i_gamma ]; i_gamma_t++) {

                        R4_sph[ j ][ k ][ i_gamma ][ i_delta ] += R4_ab_sph[j][k][ i_gamma_t ][ i_delta_t]
                                                                      * trafo_gamma(i_gamma_t, i_gamma) * trafo_delta(i_delta_t, i_delta);
                      }
                    }
                  }
                }
              }
            }

            int NumFunc_alpha = shell_alpha->getNumFunc();
            int NumFunc_beta = shell_beta->getNumFunc();
            int NumFunc_gamma = shell_gamma->getNumFunc();
            int NumFunc_delta = shell_delta->getNumFunc();
            
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
                        block[c][d][a][b] += R4_sph[offset_alpha + i_alpha][offset_beta + i_beta][offset_gamma + i_gamma][offset_delta + i_delta];
                      } else {
                        block[a][b][c][d] += R4_sph[offset_alpha + i_alpha][offset_beta + i_beta][offset_gamma + i_gamma][offset_delta + i_delta];
                      }
                    }
                  }
                }
              }

                 } // GaussianIterator itdelta
              } // GaussianIterator itgamma
           } // GaussianIterator itbeta
        } // GaussianIterator italpha

       return does_contribute;     
    } // TCrawMatrix::FillFourCenterRepBlock




    }}
