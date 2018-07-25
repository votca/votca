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

#include <votca/xtp/threecenter.h>

using namespace std;


namespace votca {
    namespace xtp {
 
        
        /*
         * Calculate 3-center electron repulsion integrals 
         *    R_{abc} = int{ phi_a(r)^DFT phi_b(r)^DFT phi_c(r')^AUX/(r-r') d3rd3r' }
         * for a given set of a b c as in 
         *    Obara, Saika, J. Chem. Phys. 84, 3963 (1986)
         * section II.B for cartesian Gaussians, then transforming
         * to spherical (angular momentum) Gaussians ("complete" shells
         * from S to Lmax, and finally cutting out those angular momentum 
         * components actually present in shell-shell-shell combination.
         * Currently supported for 
         *      S,P,D   functions in DFT basis and 
         *      S,P,D,F functions in AUX  basis
         * 
         */
        
      
        bool TCMatrix::FillThreeCenterRepBlock(tensor3d& threec_block, const AOShell* _shell_3, const AOShell* _shell_1, const AOShell* _shell_2) {

            const double pi = boost::math::constants::pi<double>();
            const double gwaccuracy = 1.e-11;
            
            bool _does_contribute=false;
            

            // shell info, only lmax tells how far to go
            
            int _lmax_1 = _shell_1->getLmax();
            int _lmax_2 = _shell_2->getLmax();
            int _lmax_3 = _shell_3->getLmax();
            
            int _mmax = _lmax_1+_lmax_2+_lmax_3;


            // set size of internal block for recursion
           
            const AOShell* _shell_alpha;
            const AOShell* _shell_beta;
            const AOShell* _shell_gamma;
            bool alphabetaswitch=false;

            // We need lmax_alpha > lmax_beta, so instead of calculating (sp,s) we calculate (ps,s), due to symmetry they are the same. 
            
            if (_lmax_1 < _lmax_2){
                _shell_alpha=_shell_2;
                _shell_beta =_shell_1;
                alphabetaswitch=true;
//               cout << "switched" << endl;   

                
            }
            else{
                _shell_alpha=_shell_1;
                _shell_beta =_shell_2;
            }
            _shell_gamma=_shell_3;
            
            const tools::vec& _pos_alpha = _shell_alpha->getPos();
            const tools::vec& _pos_beta = _shell_beta->getPos();
            const tools::vec& _pos_gamma = _shell_gamma->getPos();
            
            int _lmax_alpha = _shell_alpha->getLmax();
            int _lmax_beta  = _shell_beta->getLmax();
            int _lmax_gamma = _shell_gamma->getLmax();
            
            int _ngamma = AOSuperMatrix::getBlockSize(_lmax_gamma);
            int _nbeta = AOSuperMatrix::getBlockSize(_lmax_beta);
            int _ncombined =AOSuperMatrix::getBlockSize(_lmax_alpha+_lmax_beta);
            
                       

            //int n_orb = ((_lmax_gamma + 1)*(_lmax_gamma + 2)*(_lmax_gamma + 3))/6;
            int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84, 120, 165};
            
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


            
                    
            

            
            tools::vec amb=_pos_alpha-_pos_beta;
            double amb0=amb.getX();
            double amb1=amb.getY();
            double amb2=amb.getZ();
            double _dist3 = amb * amb;
         


            

            for ( AOShell::GaussianIterator italpha = _shell_alpha->begin(); italpha != _shell_alpha->end(); ++italpha){
                const double _decay_alpha = italpha->getDecay();
            
                for ( AOShell::GaussianIterator itbeta = _shell_beta->begin(); itbeta != _shell_beta->end(); ++itbeta){
                    const double _decay_beta = itbeta->getDecay();
                    double rzeta = 0.5 / (_decay_alpha+_decay_beta);
                    tools::vec _P = 2.0 * (_decay_alpha*_pos_alpha+_decay_beta*_pos_beta) * rzeta;
                    tools::vec pma = _P - _pos_alpha;
                    double pma0 = pma.getX();
                    double pma1 = pma.getY();
                    double pma2 = pma.getZ();
                    double xi = 2.0 * _decay_alpha * _decay_beta * rzeta;
                    double fact_alpha_beta = 16.0 * xi * pow(pi / (_decay_alpha * _decay_beta), 0.25) * exp(-xi * _dist3);
                    
                    for ( AOShell::GaussianIterator itgamma = _shell_gamma->begin(); itgamma != _shell_gamma->end(); ++itgamma){
                        const double _decay_gamma = itgamma->getDecay();
            
    


      
            
            double _decay=_decay_alpha + _decay_beta + _decay_gamma;
            double rgamma = 0.5/_decay_gamma; 
            double rdecay = 0.5/_decay; 

            double sss = fact_alpha_beta * pow(rdecay * rdecay * rgamma, 0.25);

            if (sss < gwaccuracy) { continue; }

            _does_contribute = true;

            double gfak=_decay_gamma/_decay;
            double cfak= (_decay_alpha + _decay_beta)/_decay;         
            tools::vec _W=(_decay_alpha*_pos_alpha+_decay_beta*_pos_beta+_decay_gamma*_pos_gamma)/_decay;
            double _T = (_decay_alpha+_decay_beta)*_decay_gamma/_decay*(_P-_pos_gamma)*(_P-_pos_gamma);

            tools::vec wmp = _W - _P; 
            tools::vec wmc = _W - _pos_gamma;

            double wmp0 = wmp.getX();
            double wmc0 = wmc.getX();

            double wmp1 = wmp.getY();
            double wmc1 = wmc.getY();

            double wmp2 = wmp.getZ();
            double wmc2 = wmc.getZ();
            
          
            tensor3d::extent_gen extents;
            tensor3d R_temp;
            R_temp.resize(extents[ range(0, _ncombined ) ][ range(0, _ngamma ) ][ range(0, max(2,_mmax+1))]);
            //initialize to zero
            for (index3d i = 0; i != _ncombined; ++i) {
                for (index3d j = 0; j != _ngamma; ++j) { 
                    for (index3d k = 0; k != _mmax+1; ++k) { 
                                       R_temp[i][j][k] = 0.0;
                                   }
                               }
                           }
            
            tensor3d R;
            R.resize(extents[ range(0, _ncombined ) ][ range(0, _nbeta ) ][ range(0, _ngamma)]);
            //initialize to zero
            for (index3d i = 0; i != _ncombined; ++i) {
                for (index3d j = 0; j != _nbeta; ++j) {
                    for (index3d k = 0; k != _ngamma; ++k) {

                                       R[i][j][k] = 0.0;
                                   }
                               }
                           }


            const std::vector<double> _FmT=AOMatrix<double>::XIntegrate(_mmax+1, _T);

            //ss integrals

            for (int _i=0;_i<_mmax+1;_i++){
                R_temp[0][0][_i]=sss*_FmT[_i];
            }

int _lmax_alpha_beta = _lmax_alpha + _lmax_beta;


//Integral  p - s - s
if (_lmax_alpha_beta > 0) {
  for (int m = 0; m < _mmax; m++) {
    R_temp[Cart::x][0][m] = pma0*R_temp[0][0][m] + wmp0*R_temp[0][0][m+1];
    R_temp[Cart::y][0][m] = pma1*R_temp[0][0][m] + wmp1*R_temp[0][0][m+1];
    R_temp[Cart::z][0][m] = pma2*R_temp[0][0][m] + wmp2*R_temp[0][0][m+1];
  }
}
//------------------------------------------------------

//Integral  d - s - s
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

//Integral  f - s - s
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

//Integral  g - s - s
if (_lmax_alpha_beta > 3) {
  for (int m = 0; m < _mmax-3; m++) {
    double term_xx = rzeta*(R_temp[Cart::xx][0][m]-gfak*R_temp[Cart::xx][0][m+1]);
    double term_yy = rzeta*(R_temp[Cart::yy][0][m]-gfak*R_temp[Cart::yy][0][m+1]);
    double term_zz = rzeta*(R_temp[Cart::zz][0][m]-gfak*R_temp[Cart::zz][0][m+1]);
    R_temp[Cart::xxxx][0][m] = pma0*R_temp[Cart::xxx][0][m] + wmp0*R_temp[Cart::xxx][0][m+1] + 3*term_xx;
    R_temp[Cart::xxxy][0][m] = pma1*R_temp[Cart::xxx][0][m] + wmp1*R_temp[Cart::xxx][0][m+1];
    R_temp[Cart::xxxz][0][m] = pma2*R_temp[Cart::xxx][0][m] + wmp2*R_temp[Cart::xxx][0][m+1];
    R_temp[Cart::xxyy][0][m] = pma0*R_temp[Cart::xyy][0][m] + wmp0*R_temp[Cart::xyy][0][m+1] + term_yy;
    R_temp[Cart::xxyz][0][m] = pma1*R_temp[Cart::xxz][0][m] + wmp1*R_temp[Cart::xxz][0][m+1];
    R_temp[Cart::xxzz][0][m] = pma0*R_temp[Cart::xzz][0][m] + wmp0*R_temp[Cart::xzz][0][m+1] + term_zz;
    R_temp[Cart::xyyy][0][m] = pma0*R_temp[Cart::yyy][0][m] + wmp0*R_temp[Cart::yyy][0][m+1];
    R_temp[Cart::xyyz][0][m] = pma0*R_temp[Cart::yyz][0][m] + wmp0*R_temp[Cart::yyz][0][m+1];
    R_temp[Cart::xyzz][0][m] = pma0*R_temp[Cart::yzz][0][m] + wmp0*R_temp[Cart::yzz][0][m+1];
    R_temp[Cart::xzzz][0][m] = pma0*R_temp[Cart::zzz][0][m] + wmp0*R_temp[Cart::zzz][0][m+1];
    R_temp[Cart::yyyy][0][m] = pma1*R_temp[Cart::yyy][0][m] + wmp1*R_temp[Cart::yyy][0][m+1] + 3*term_yy;
    R_temp[Cart::yyyz][0][m] = pma2*R_temp[Cart::yyy][0][m] + wmp2*R_temp[Cart::yyy][0][m+1];
    R_temp[Cart::yyzz][0][m] = pma1*R_temp[Cart::yzz][0][m] + wmp1*R_temp[Cart::yzz][0][m+1] + term_zz;
    R_temp[Cart::yzzz][0][m] = pma1*R_temp[Cart::zzz][0][m] + wmp1*R_temp[Cart::zzz][0][m+1];
    R_temp[Cart::zzzz][0][m] = pma2*R_temp[Cart::zzz][0][m] + wmp2*R_temp[Cart::zzz][0][m+1] + 3*term_zz;
  }
}
//------------------------------------------------------

//Integral  h - s - s
if (_lmax_alpha_beta > 4) {
  for (int m = 0; m < _mmax-4; m++) {
    double term_xxx = rzeta*(R_temp[Cart::xxx][0][m]-gfak*R_temp[Cart::xxx][0][m+1]);
    double term_yyy = rzeta*(R_temp[Cart::yyy][0][m]-gfak*R_temp[Cart::yyy][0][m+1]);
    double term_zzz = rzeta*(R_temp[Cart::zzz][0][m]-gfak*R_temp[Cart::zzz][0][m+1]);
    R_temp[Cart::xxxxx][0][m] = pma0*R_temp[Cart::xxxx][0][m] + wmp0*R_temp[Cart::xxxx][0][m+1] + 4*term_xxx;
    R_temp[Cart::xxxxy][0][m] = pma1*R_temp[Cart::xxxx][0][m] + wmp1*R_temp[Cart::xxxx][0][m+1];
    R_temp[Cart::xxxxz][0][m] = pma2*R_temp[Cart::xxxx][0][m] + wmp2*R_temp[Cart::xxxx][0][m+1];
    R_temp[Cart::xxxyy][0][m] = pma1*R_temp[Cart::xxxy][0][m] + wmp1*R_temp[Cart::xxxy][0][m+1] + term_xxx;
    R_temp[Cart::xxxyz][0][m] = pma1*R_temp[Cart::xxxz][0][m] + wmp1*R_temp[Cart::xxxz][0][m+1];
    R_temp[Cart::xxxzz][0][m] = pma2*R_temp[Cart::xxxz][0][m] + wmp2*R_temp[Cart::xxxz][0][m+1] + term_xxx;
    R_temp[Cart::xxyyy][0][m] = pma0*R_temp[Cart::xyyy][0][m] + wmp0*R_temp[Cart::xyyy][0][m+1] + term_yyy;
    R_temp[Cart::xxyyz][0][m] = pma2*R_temp[Cart::xxyy][0][m] + wmp2*R_temp[Cart::xxyy][0][m+1];
    R_temp[Cart::xxyzz][0][m] = pma1*R_temp[Cart::xxzz][0][m] + wmp1*R_temp[Cart::xxzz][0][m+1];
    R_temp[Cart::xxzzz][0][m] = pma0*R_temp[Cart::xzzz][0][m] + wmp0*R_temp[Cart::xzzz][0][m+1] + term_zzz;
    R_temp[Cart::xyyyy][0][m] = pma0*R_temp[Cart::yyyy][0][m] + wmp0*R_temp[Cart::yyyy][0][m+1];
    R_temp[Cart::xyyyz][0][m] = pma0*R_temp[Cart::yyyz][0][m] + wmp0*R_temp[Cart::yyyz][0][m+1];
    R_temp[Cart::xyyzz][0][m] = pma0*R_temp[Cart::yyzz][0][m] + wmp0*R_temp[Cart::yyzz][0][m+1];
    R_temp[Cart::xyzzz][0][m] = pma0*R_temp[Cart::yzzz][0][m] + wmp0*R_temp[Cart::yzzz][0][m+1];
    R_temp[Cart::xzzzz][0][m] = pma0*R_temp[Cart::zzzz][0][m] + wmp0*R_temp[Cart::zzzz][0][m+1];
    R_temp[Cart::yyyyy][0][m] = pma1*R_temp[Cart::yyyy][0][m] + wmp1*R_temp[Cart::yyyy][0][m+1] + 4*term_yyy;
    R_temp[Cart::yyyyz][0][m] = pma2*R_temp[Cart::yyyy][0][m] + wmp2*R_temp[Cart::yyyy][0][m+1];
    R_temp[Cart::yyyzz][0][m] = pma2*R_temp[Cart::yyyz][0][m] + wmp2*R_temp[Cart::yyyz][0][m+1] + term_yyy;
    R_temp[Cart::yyzzz][0][m] = pma1*R_temp[Cart::yzzz][0][m] + wmp1*R_temp[Cart::yzzz][0][m+1] + term_zzz;
    R_temp[Cart::yzzzz][0][m] = pma1*R_temp[Cart::zzzz][0][m] + wmp1*R_temp[Cart::zzzz][0][m+1];
    R_temp[Cart::zzzzz][0][m] = pma2*R_temp[Cart::zzzz][0][m] + wmp2*R_temp[Cart::zzzz][0][m+1] + 4*term_zzz;
  }
}
//------------------------------------------------------

//Integral  i - s - s
if (_lmax_alpha_beta > 5) {
  for (int m = 0; m < _mmax-5; m++) {
    double term_xxxx = rzeta*(R_temp[Cart::xxxx][0][m]-gfak*R_temp[Cart::xxxx][0][m+1]);
    double term_xyyy = rzeta*(R_temp[Cart::xyyy][0][m]-gfak*R_temp[Cart::xyyy][0][m+1]);
    double term_xzzz = rzeta*(R_temp[Cart::xzzz][0][m]-gfak*R_temp[Cart::xzzz][0][m+1]);
    double term_yyyy = rzeta*(R_temp[Cart::yyyy][0][m]-gfak*R_temp[Cart::yyyy][0][m+1]);
    double term_yyzz = rzeta*(R_temp[Cart::yyzz][0][m]-gfak*R_temp[Cart::yyzz][0][m+1]);
    double term_yzzz = rzeta*(R_temp[Cart::yzzz][0][m]-gfak*R_temp[Cart::yzzz][0][m+1]);
    double term_zzzz = rzeta*(R_temp[Cart::zzzz][0][m]-gfak*R_temp[Cart::zzzz][0][m+1]);
    R_temp[Cart::xxxxxx][0][m] = pma0*R_temp[Cart::xxxxx][0][m] + wmp0*R_temp[Cart::xxxxx][0][m+1] + 5*term_xxxx;
    R_temp[Cart::xxxxxy][0][m] = pma1*R_temp[Cart::xxxxx][0][m] + wmp1*R_temp[Cart::xxxxx][0][m+1];
    R_temp[Cart::xxxxxz][0][m] = pma2*R_temp[Cart::xxxxx][0][m] + wmp2*R_temp[Cart::xxxxx][0][m+1];
    R_temp[Cart::xxxxyy][0][m] = pma1*R_temp[Cart::xxxxy][0][m] + wmp1*R_temp[Cart::xxxxy][0][m+1] + term_xxxx;
    R_temp[Cart::xxxxyz][0][m] = pma1*R_temp[Cart::xxxxz][0][m] + wmp1*R_temp[Cart::xxxxz][0][m+1];
    R_temp[Cart::xxxxzz][0][m] = pma2*R_temp[Cart::xxxxz][0][m] + wmp2*R_temp[Cart::xxxxz][0][m+1] + term_xxxx;
    R_temp[Cart::xxxyyy][0][m] = pma0*R_temp[Cart::xxyyy][0][m] + wmp0*R_temp[Cart::xxyyy][0][m+1] + 2*term_xyyy;
    R_temp[Cart::xxxyyz][0][m] = pma2*R_temp[Cart::xxxyy][0][m] + wmp2*R_temp[Cart::xxxyy][0][m+1];
    R_temp[Cart::xxxyzz][0][m] = pma1*R_temp[Cart::xxxzz][0][m] + wmp1*R_temp[Cart::xxxzz][0][m+1];
    R_temp[Cart::xxxzzz][0][m] = pma0*R_temp[Cart::xxzzz][0][m] + wmp0*R_temp[Cart::xxzzz][0][m+1] + 2*term_xzzz;
    R_temp[Cart::xxyyyy][0][m] = pma0*R_temp[Cart::xyyyy][0][m] + wmp0*R_temp[Cart::xyyyy][0][m+1] + term_yyyy;
    R_temp[Cart::xxyyyz][0][m] = pma2*R_temp[Cart::xxyyy][0][m] + wmp2*R_temp[Cart::xxyyy][0][m+1];
    R_temp[Cart::xxyyzz][0][m] = pma0*R_temp[Cart::xyyzz][0][m] + wmp0*R_temp[Cart::xyyzz][0][m+1] + term_yyzz;
    R_temp[Cart::xxyzzz][0][m] = pma1*R_temp[Cart::xxzzz][0][m] + wmp1*R_temp[Cart::xxzzz][0][m+1];
    R_temp[Cart::xxzzzz][0][m] = pma0*R_temp[Cart::xzzzz][0][m] + wmp0*R_temp[Cart::xzzzz][0][m+1] + term_zzzz;
    R_temp[Cart::xyyyyy][0][m] = pma0*R_temp[Cart::yyyyy][0][m] + wmp0*R_temp[Cart::yyyyy][0][m+1];
    R_temp[Cart::xyyyyz][0][m] = pma0*R_temp[Cart::yyyyz][0][m] + wmp0*R_temp[Cart::yyyyz][0][m+1];
    R_temp[Cart::xyyyzz][0][m] = pma0*R_temp[Cart::yyyzz][0][m] + wmp0*R_temp[Cart::yyyzz][0][m+1];
    R_temp[Cart::xyyzzz][0][m] = pma0*R_temp[Cart::yyzzz][0][m] + wmp0*R_temp[Cart::yyzzz][0][m+1];
    R_temp[Cart::xyzzzz][0][m] = pma0*R_temp[Cart::yzzzz][0][m] + wmp0*R_temp[Cart::yzzzz][0][m+1];
    R_temp[Cart::xzzzzz][0][m] = pma0*R_temp[Cart::zzzzz][0][m] + wmp0*R_temp[Cart::zzzzz][0][m+1];
    R_temp[Cart::yyyyyy][0][m] = pma1*R_temp[Cart::yyyyy][0][m] + wmp1*R_temp[Cart::yyyyy][0][m+1] + 5*term_yyyy;
    R_temp[Cart::yyyyyz][0][m] = pma2*R_temp[Cart::yyyyy][0][m] + wmp2*R_temp[Cart::yyyyy][0][m+1];
    R_temp[Cart::yyyyzz][0][m] = pma2*R_temp[Cart::yyyyz][0][m] + wmp2*R_temp[Cart::yyyyz][0][m+1] + term_yyyy;
    R_temp[Cart::yyyzzz][0][m] = pma1*R_temp[Cart::yyzzz][0][m] + wmp1*R_temp[Cart::yyzzz][0][m+1] + 2*term_yzzz;
    R_temp[Cart::yyzzzz][0][m] = pma1*R_temp[Cart::yzzzz][0][m] + wmp1*R_temp[Cart::yzzzz][0][m+1] + term_zzzz;
    R_temp[Cart::yzzzzz][0][m] = pma1*R_temp[Cart::zzzzz][0][m] + wmp1*R_temp[Cart::zzzzz][0][m+1];
    R_temp[Cart::zzzzzz][0][m] = pma2*R_temp[Cart::zzzzz][0][m] + wmp2*R_temp[Cart::zzzzz][0][m+1] + 5*term_zzzz;
  }
}
//------------------------------------------------------

//Integral  j - s - s
if (_lmax_alpha_beta > 6) {
  for (int m = 0; m < _mmax-6; m++) {
    double term_xxxxx = rzeta*(R_temp[Cart::xxxxx][0][m]-gfak*R_temp[Cart::xxxxx][0][m+1]);
    double term_xxxxy = rzeta*(R_temp[Cart::xxxxy][0][m]-gfak*R_temp[Cart::xxxxy][0][m+1]);
    double term_xxxxz = rzeta*(R_temp[Cart::xxxxz][0][m]-gfak*R_temp[Cart::xxxxz][0][m+1]);
    double term_xxxzz = rzeta*(R_temp[Cart::xxxzz][0][m]-gfak*R_temp[Cart::xxxzz][0][m+1]);
    double term_xyyyy = rzeta*(R_temp[Cart::xyyyy][0][m]-gfak*R_temp[Cart::xyyyy][0][m+1]);
    double term_xzzzz = rzeta*(R_temp[Cart::xzzzz][0][m]-gfak*R_temp[Cart::xzzzz][0][m+1]);
    double term_yyyyy = rzeta*(R_temp[Cart::yyyyy][0][m]-gfak*R_temp[Cart::yyyyy][0][m+1]);
    double term_yyyyz = rzeta*(R_temp[Cart::yyyyz][0][m]-gfak*R_temp[Cart::yyyyz][0][m+1]);
    double term_yyyzz = rzeta*(R_temp[Cart::yyyzz][0][m]-gfak*R_temp[Cart::yyyzz][0][m+1]);
    double term_yyzzz = rzeta*(R_temp[Cart::yyzzz][0][m]-gfak*R_temp[Cart::yyzzz][0][m+1]);
    double term_yzzzz = rzeta*(R_temp[Cart::yzzzz][0][m]-gfak*R_temp[Cart::yzzzz][0][m+1]);
    double term_zzzzz = rzeta*(R_temp[Cart::zzzzz][0][m]-gfak*R_temp[Cart::zzzzz][0][m+1]);
    R_temp[Cart::xxxxxxx][0][m] = pma0*R_temp[Cart::xxxxxx][0][m] + wmp0*R_temp[Cart::xxxxxx][0][m+1] + 6*term_xxxxx;
    R_temp[Cart::xxxxxxy][0][m] = pma1*R_temp[Cart::xxxxxx][0][m] + wmp1*R_temp[Cart::xxxxxx][0][m+1];
    R_temp[Cart::xxxxxxz][0][m] = pma2*R_temp[Cart::xxxxxx][0][m] + wmp2*R_temp[Cart::xxxxxx][0][m+1];
    R_temp[Cart::xxxxxyy][0][m] = pma1*R_temp[Cart::xxxxxy][0][m] + wmp1*R_temp[Cart::xxxxxy][0][m+1] + term_xxxxx;
    R_temp[Cart::xxxxxyz][0][m] = pma1*R_temp[Cart::xxxxxz][0][m] + wmp1*R_temp[Cart::xxxxxz][0][m+1];
    R_temp[Cart::xxxxxzz][0][m] = pma2*R_temp[Cart::xxxxxz][0][m] + wmp2*R_temp[Cart::xxxxxz][0][m+1] + term_xxxxx;
    R_temp[Cart::xxxxyyy][0][m] = pma1*R_temp[Cart::xxxxyy][0][m] + wmp1*R_temp[Cart::xxxxyy][0][m+1] + 2*term_xxxxy;
    R_temp[Cart::xxxxyyz][0][m] = pma2*R_temp[Cart::xxxxyy][0][m] + wmp2*R_temp[Cart::xxxxyy][0][m+1];
    R_temp[Cart::xxxxyzz][0][m] = pma1*R_temp[Cart::xxxxzz][0][m] + wmp1*R_temp[Cart::xxxxzz][0][m+1];
    R_temp[Cart::xxxxzzz][0][m] = pma2*R_temp[Cart::xxxxzz][0][m] + wmp2*R_temp[Cart::xxxxzz][0][m+1] + 2*term_xxxxz;
    R_temp[Cart::xxxyyyy][0][m] = pma0*R_temp[Cart::xxyyyy][0][m] + wmp0*R_temp[Cart::xxyyyy][0][m+1] + 2*term_xyyyy;
    R_temp[Cart::xxxyyyz][0][m] = pma2*R_temp[Cart::xxxyyy][0][m] + wmp2*R_temp[Cart::xxxyyy][0][m+1];
    R_temp[Cart::xxxyyzz][0][m] = pma1*R_temp[Cart::xxxyzz][0][m] + wmp1*R_temp[Cart::xxxyzz][0][m+1] + term_xxxzz;
    R_temp[Cart::xxxyzzz][0][m] = pma1*R_temp[Cart::xxxzzz][0][m] + wmp1*R_temp[Cart::xxxzzz][0][m+1];
    R_temp[Cart::xxxzzzz][0][m] = pma0*R_temp[Cart::xxzzzz][0][m] + wmp0*R_temp[Cart::xxzzzz][0][m+1] + 2*term_xzzzz;
    R_temp[Cart::xxyyyyy][0][m] = pma0*R_temp[Cart::xyyyyy][0][m] + wmp0*R_temp[Cart::xyyyyy][0][m+1] + term_yyyyy;
    R_temp[Cart::xxyyyyz][0][m] = pma2*R_temp[Cart::xxyyyy][0][m] + wmp2*R_temp[Cart::xxyyyy][0][m+1];
    R_temp[Cart::xxyyyzz][0][m] = pma0*R_temp[Cart::xyyyzz][0][m] + wmp0*R_temp[Cart::xyyyzz][0][m+1] + term_yyyzz;
    R_temp[Cart::xxyyzzz][0][m] = pma0*R_temp[Cart::xyyzzz][0][m] + wmp0*R_temp[Cart::xyyzzz][0][m+1] + term_yyzzz;
    R_temp[Cart::xxyzzzz][0][m] = pma1*R_temp[Cart::xxzzzz][0][m] + wmp1*R_temp[Cart::xxzzzz][0][m+1];
    R_temp[Cart::xxzzzzz][0][m] = pma0*R_temp[Cart::xzzzzz][0][m] + wmp0*R_temp[Cart::xzzzzz][0][m+1] + term_zzzzz;
    R_temp[Cart::xyyyyyy][0][m] = pma0*R_temp[Cart::yyyyyy][0][m] + wmp0*R_temp[Cart::yyyyyy][0][m+1];
    R_temp[Cart::xyyyyyz][0][m] = pma0*R_temp[Cart::yyyyyz][0][m] + wmp0*R_temp[Cart::yyyyyz][0][m+1];
    R_temp[Cart::xyyyyzz][0][m] = pma0*R_temp[Cart::yyyyzz][0][m] + wmp0*R_temp[Cart::yyyyzz][0][m+1];
    R_temp[Cart::xyyyzzz][0][m] = pma0*R_temp[Cart::yyyzzz][0][m] + wmp0*R_temp[Cart::yyyzzz][0][m+1];
    R_temp[Cart::xyyzzzz][0][m] = pma0*R_temp[Cart::yyzzzz][0][m] + wmp0*R_temp[Cart::yyzzzz][0][m+1];
    R_temp[Cart::xyzzzzz][0][m] = pma0*R_temp[Cart::yzzzzz][0][m] + wmp0*R_temp[Cart::yzzzzz][0][m+1];
    R_temp[Cart::xzzzzzz][0][m] = pma0*R_temp[Cart::zzzzzz][0][m] + wmp0*R_temp[Cart::zzzzzz][0][m+1];
    R_temp[Cart::yyyyyyy][0][m] = pma1*R_temp[Cart::yyyyyy][0][m] + wmp1*R_temp[Cart::yyyyyy][0][m+1] + 6*term_yyyyy;
    R_temp[Cart::yyyyyyz][0][m] = pma2*R_temp[Cart::yyyyyy][0][m] + wmp2*R_temp[Cart::yyyyyy][0][m+1];
    R_temp[Cart::yyyyyzz][0][m] = pma2*R_temp[Cart::yyyyyz][0][m] + wmp2*R_temp[Cart::yyyyyz][0][m+1] + term_yyyyy;
    R_temp[Cart::yyyyzzz][0][m] = pma2*R_temp[Cart::yyyyzz][0][m] + wmp2*R_temp[Cart::yyyyzz][0][m+1] + 2*term_yyyyz;
    R_temp[Cart::yyyzzzz][0][m] = pma1*R_temp[Cart::yyzzzz][0][m] + wmp1*R_temp[Cart::yyzzzz][0][m+1] + 2*term_yzzzz;
    R_temp[Cart::yyzzzzz][0][m] = pma1*R_temp[Cart::yzzzzz][0][m] + wmp1*R_temp[Cart::yzzzzz][0][m+1] + term_zzzzz;
    R_temp[Cart::yzzzzzz][0][m] = pma1*R_temp[Cart::zzzzzz][0][m] + wmp1*R_temp[Cart::zzzzzz][0][m+1];
    R_temp[Cart::zzzzzzz][0][m] = pma2*R_temp[Cart::zzzzzz][0][m] + wmp2*R_temp[Cart::zzzzzz][0][m+1] + 6*term_zzzzz;
  }
}
//------------------------------------------------------

//Integral  k - s - s
if (_lmax_alpha_beta > 7) {
  for (int m = 0; m < _mmax-7; m++) {
    double term_xxxxxx = rzeta*(R_temp[Cart::xxxxxx][0][m]-gfak*R_temp[Cart::xxxxxx][0][m+1]);
    double term_xxxxxy = rzeta*(R_temp[Cart::xxxxxy][0][m]-gfak*R_temp[Cart::xxxxxy][0][m+1]);
    double term_xxxxxz = rzeta*(R_temp[Cart::xxxxxz][0][m]-gfak*R_temp[Cart::xxxxxz][0][m+1]);
    double term_xxxxzz = rzeta*(R_temp[Cart::xxxxzz][0][m]-gfak*R_temp[Cart::xxxxzz][0][m+1]);
    double term_xxxyyy = rzeta*(R_temp[Cart::xxxyyy][0][m]-gfak*R_temp[Cart::xxxyyy][0][m+1]);
    double term_xxxzzz = rzeta*(R_temp[Cart::xxxzzz][0][m]-gfak*R_temp[Cart::xxxzzz][0][m+1]);
    double term_xxyyyy = rzeta*(R_temp[Cart::xxyyyy][0][m]-gfak*R_temp[Cart::xxyyyy][0][m+1]);
    double term_xxzzzz = rzeta*(R_temp[Cart::xxzzzz][0][m]-gfak*R_temp[Cart::xxzzzz][0][m+1]);
    double term_xyyyyy = rzeta*(R_temp[Cart::xyyyyy][0][m]-gfak*R_temp[Cart::xyyyyy][0][m+1]);
    double term_xzzzzz = rzeta*(R_temp[Cart::xzzzzz][0][m]-gfak*R_temp[Cart::xzzzzz][0][m+1]);
    double term_yyyyyy = rzeta*(R_temp[Cart::yyyyyy][0][m]-gfak*R_temp[Cart::yyyyyy][0][m+1]);
    double term_yyyyyz = rzeta*(R_temp[Cart::yyyyyz][0][m]-gfak*R_temp[Cart::yyyyyz][0][m+1]);
    double term_yyyyzz = rzeta*(R_temp[Cart::yyyyzz][0][m]-gfak*R_temp[Cart::yyyyzz][0][m+1]);
    double term_yyyzzz = rzeta*(R_temp[Cart::yyyzzz][0][m]-gfak*R_temp[Cart::yyyzzz][0][m+1]);
    double term_yyzzzz = rzeta*(R_temp[Cart::yyzzzz][0][m]-gfak*R_temp[Cart::yyzzzz][0][m+1]);
    double term_yzzzzz = rzeta*(R_temp[Cart::yzzzzz][0][m]-gfak*R_temp[Cart::yzzzzz][0][m+1]);
    double term_zzzzzz = rzeta*(R_temp[Cart::zzzzzz][0][m]-gfak*R_temp[Cart::zzzzzz][0][m+1]);
    R_temp[Cart::xxxxxxxx][0][m] = pma0*R_temp[Cart::xxxxxxx][0][m] + wmp0*R_temp[Cart::xxxxxxx][0][m+1] + 7*term_xxxxxx;
    R_temp[Cart::xxxxxxxy][0][m] = pma1*R_temp[Cart::xxxxxxx][0][m] + wmp1*R_temp[Cart::xxxxxxx][0][m+1];
    R_temp[Cart::xxxxxxxz][0][m] = pma2*R_temp[Cart::xxxxxxx][0][m] + wmp2*R_temp[Cart::xxxxxxx][0][m+1];
    R_temp[Cart::xxxxxxyy][0][m] = pma1*R_temp[Cart::xxxxxxy][0][m] + wmp1*R_temp[Cart::xxxxxxy][0][m+1] + term_xxxxxx;
    R_temp[Cart::xxxxxxyz][0][m] = pma1*R_temp[Cart::xxxxxxz][0][m] + wmp1*R_temp[Cart::xxxxxxz][0][m+1];
    R_temp[Cart::xxxxxxzz][0][m] = pma2*R_temp[Cart::xxxxxxz][0][m] + wmp2*R_temp[Cart::xxxxxxz][0][m+1] + term_xxxxxx;
    R_temp[Cart::xxxxxyyy][0][m] = pma1*R_temp[Cart::xxxxxyy][0][m] + wmp1*R_temp[Cart::xxxxxyy][0][m+1] + 2*term_xxxxxy;
    R_temp[Cart::xxxxxyyz][0][m] = pma2*R_temp[Cart::xxxxxyy][0][m] + wmp2*R_temp[Cart::xxxxxyy][0][m+1];
    R_temp[Cart::xxxxxyzz][0][m] = pma1*R_temp[Cart::xxxxxzz][0][m] + wmp1*R_temp[Cart::xxxxxzz][0][m+1];
    R_temp[Cart::xxxxxzzz][0][m] = pma2*R_temp[Cart::xxxxxzz][0][m] + wmp2*R_temp[Cart::xxxxxzz][0][m+1] + 2*term_xxxxxz;
    R_temp[Cart::xxxxyyyy][0][m] = pma0*R_temp[Cart::xxxyyyy][0][m] + wmp0*R_temp[Cart::xxxyyyy][0][m+1] + 3*term_xxyyyy;
    R_temp[Cart::xxxxyyyz][0][m] = pma2*R_temp[Cart::xxxxyyy][0][m] + wmp2*R_temp[Cart::xxxxyyy][0][m+1];
    R_temp[Cart::xxxxyyzz][0][m] = pma1*R_temp[Cart::xxxxyzz][0][m] + wmp1*R_temp[Cart::xxxxyzz][0][m+1] + term_xxxxzz;
    R_temp[Cart::xxxxyzzz][0][m] = pma1*R_temp[Cart::xxxxzzz][0][m] + wmp1*R_temp[Cart::xxxxzzz][0][m+1];
    R_temp[Cart::xxxxzzzz][0][m] = pma0*R_temp[Cart::xxxzzzz][0][m] + wmp0*R_temp[Cart::xxxzzzz][0][m+1] + 3*term_xxzzzz;
    R_temp[Cart::xxxyyyyy][0][m] = pma0*R_temp[Cart::xxyyyyy][0][m] + wmp0*R_temp[Cart::xxyyyyy][0][m+1] + 2*term_xyyyyy;
    R_temp[Cart::xxxyyyyz][0][m] = pma2*R_temp[Cart::xxxyyyy][0][m] + wmp2*R_temp[Cart::xxxyyyy][0][m+1];
    R_temp[Cart::xxxyyyzz][0][m] = pma2*R_temp[Cart::xxxyyyz][0][m] + wmp2*R_temp[Cart::xxxyyyz][0][m+1] + term_xxxyyy;
    R_temp[Cart::xxxyyzzz][0][m] = pma1*R_temp[Cart::xxxyzzz][0][m] + wmp1*R_temp[Cart::xxxyzzz][0][m+1] + term_xxxzzz;
    R_temp[Cart::xxxyzzzz][0][m] = pma1*R_temp[Cart::xxxzzzz][0][m] + wmp1*R_temp[Cart::xxxzzzz][0][m+1];
    R_temp[Cart::xxxzzzzz][0][m] = pma0*R_temp[Cart::xxzzzzz][0][m] + wmp0*R_temp[Cart::xxzzzzz][0][m+1] + 2*term_xzzzzz;
    R_temp[Cart::xxyyyyyy][0][m] = pma0*R_temp[Cart::xyyyyyy][0][m] + wmp0*R_temp[Cart::xyyyyyy][0][m+1] + term_yyyyyy;
    R_temp[Cart::xxyyyyyz][0][m] = pma2*R_temp[Cart::xxyyyyy][0][m] + wmp2*R_temp[Cart::xxyyyyy][0][m+1];
    R_temp[Cart::xxyyyyzz][0][m] = pma0*R_temp[Cart::xyyyyzz][0][m] + wmp0*R_temp[Cart::xyyyyzz][0][m+1] + term_yyyyzz;
    R_temp[Cart::xxyyyzzz][0][m] = pma0*R_temp[Cart::xyyyzzz][0][m] + wmp0*R_temp[Cart::xyyyzzz][0][m+1] + term_yyyzzz;
    R_temp[Cart::xxyyzzzz][0][m] = pma0*R_temp[Cart::xyyzzzz][0][m] + wmp0*R_temp[Cart::xyyzzzz][0][m+1] + term_yyzzzz;
    R_temp[Cart::xxyzzzzz][0][m] = pma1*R_temp[Cart::xxzzzzz][0][m] + wmp1*R_temp[Cart::xxzzzzz][0][m+1];
    R_temp[Cart::xxzzzzzz][0][m] = pma0*R_temp[Cart::xzzzzzz][0][m] + wmp0*R_temp[Cart::xzzzzzz][0][m+1] + term_zzzzzz;
    R_temp[Cart::xyyyyyyy][0][m] = pma0*R_temp[Cart::yyyyyyy][0][m] + wmp0*R_temp[Cart::yyyyyyy][0][m+1];
    R_temp[Cart::xyyyyyyz][0][m] = pma0*R_temp[Cart::yyyyyyz][0][m] + wmp0*R_temp[Cart::yyyyyyz][0][m+1];
    R_temp[Cart::xyyyyyzz][0][m] = pma0*R_temp[Cart::yyyyyzz][0][m] + wmp0*R_temp[Cart::yyyyyzz][0][m+1];
    R_temp[Cart::xyyyyzzz][0][m] = pma0*R_temp[Cart::yyyyzzz][0][m] + wmp0*R_temp[Cart::yyyyzzz][0][m+1];
    R_temp[Cart::xyyyzzzz][0][m] = pma0*R_temp[Cart::yyyzzzz][0][m] + wmp0*R_temp[Cart::yyyzzzz][0][m+1];
    R_temp[Cart::xyyzzzzz][0][m] = pma0*R_temp[Cart::yyzzzzz][0][m] + wmp0*R_temp[Cart::yyzzzzz][0][m+1];
    R_temp[Cart::xyzzzzzz][0][m] = pma0*R_temp[Cart::yzzzzzz][0][m] + wmp0*R_temp[Cart::yzzzzzz][0][m+1];
    R_temp[Cart::xzzzzzzz][0][m] = pma0*R_temp[Cart::zzzzzzz][0][m] + wmp0*R_temp[Cart::zzzzzzz][0][m+1];
    R_temp[Cart::yyyyyyyy][0][m] = pma1*R_temp[Cart::yyyyyyy][0][m] + wmp1*R_temp[Cart::yyyyyyy][0][m+1] + 7*term_yyyyyy;
    R_temp[Cart::yyyyyyyz][0][m] = pma2*R_temp[Cart::yyyyyyy][0][m] + wmp2*R_temp[Cart::yyyyyyy][0][m+1];
    R_temp[Cart::yyyyyyzz][0][m] = pma2*R_temp[Cart::yyyyyyz][0][m] + wmp2*R_temp[Cart::yyyyyyz][0][m+1] + term_yyyyyy;
    R_temp[Cart::yyyyyzzz][0][m] = pma2*R_temp[Cart::yyyyyzz][0][m] + wmp2*R_temp[Cart::yyyyyzz][0][m+1] + 2*term_yyyyyz;
    R_temp[Cart::yyyyzzzz][0][m] = pma1*R_temp[Cart::yyyzzzz][0][m] + wmp1*R_temp[Cart::yyyzzzz][0][m+1] + 3*term_yyzzzz;
    R_temp[Cart::yyyzzzzz][0][m] = pma1*R_temp[Cart::yyzzzzz][0][m] + wmp1*R_temp[Cart::yyzzzzz][0][m+1] + 2*term_yzzzzz;
    R_temp[Cart::yyzzzzzz][0][m] = pma1*R_temp[Cart::yzzzzzz][0][m] + wmp1*R_temp[Cart::yzzzzzz][0][m+1] + term_zzzzzz;
    R_temp[Cart::yzzzzzzz][0][m] = pma1*R_temp[Cart::zzzzzzz][0][m] + wmp1*R_temp[Cart::zzzzzzz][0][m+1];
    R_temp[Cart::zzzzzzzz][0][m] = pma2*R_temp[Cart::zzzzzzz][0][m] + wmp2*R_temp[Cart::zzzzzzz][0][m+1] + 7*term_zzzzzz;
  }
}
//------------------------------------------------------







if (_lmax_gamma > 0) {

  //Integral  s - s - p
  for (int m = 0; m < _lmax_gamma; m++) {
    R_temp[0][Cart::x][m] = wmc0*R_temp[0][0][m+1];
    R_temp[0][Cart::y][m] = wmc1*R_temp[0][0][m+1];
    R_temp[0][Cart::z][m] = wmc2*R_temp[0][0][m+1];
  }
  //------------------------------------------------------

  //Integral  p - s - p
  if (_lmax_alpha_beta > 0) {
    for (int m = 0; m < _lmax_gamma; m++) {
      double term = rdecay*R_temp[0][0][m+1];
      for (int _i =  1; _i < 4; _i++) {
        R_temp[_i][Cart::x][m] = wmc0*R_temp[_i][0][m+1] + nx[_i]*term;
        R_temp[_i][Cart::y][m] = wmc1*R_temp[_i][0][m+1] + ny[_i]*term;
        R_temp[_i][Cart::z][m] = wmc2*R_temp[_i][0][m+1] + nz[_i]*term;
      }
    }
  }
  //------------------------------------------------------

  //Integrals     d - s - p     f - s - p     g - s - p     h - s - p     i - s - p     j - s - p     k - s - p
  for (int m = 0; m < _lmax_gamma; m++) {
    for (int _i =  4; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      R_temp[_i][Cart::x][m] = wmc0*R_temp[_i][0][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][0][m+1];
      R_temp[_i][Cart::y][m] = wmc1*R_temp[_i][0][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][0][m+1];
      R_temp[_i][Cart::z][m] = wmc2*R_temp[_i][0][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][0][m+1];
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 0)


if (_lmax_gamma > 1) {

  //Integral  s - s - d
  for (int m = 0; m < _lmax_gamma-1; m++) {
    double term = rgamma*(R_temp[0][0][m]-cfak*R_temp[0][0][m+1]);
    R_temp[0][Cart::xx][m] = wmc0*R_temp[0][Cart::x][m+1] + term;
    R_temp[0][Cart::xy][m] = wmc0*R_temp[0][Cart::y][m+1];
    R_temp[0][Cart::xz][m] = wmc0*R_temp[0][Cart::z][m+1];
    R_temp[0][Cart::yy][m] = wmc1*R_temp[0][Cart::y][m+1] + term;
    R_temp[0][Cart::yz][m] = wmc1*R_temp[0][Cart::z][m+1];
    R_temp[0][Cart::zz][m] = wmc2*R_temp[0][Cart::z][m+1] + term;
  }
  //------------------------------------------------------

  //Integrals     p - s - d     d - s - d     f - s - d     g - s - d     h - s - d     i - s - d     j - s - d     k - s - d
  for (int m = 0; m < _lmax_gamma-1; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      double term = rgamma*(R_temp[_i][0][m]-cfak*R_temp[_i][0][m+1]);
      R_temp[_i][Cart::xx][m] = wmc0*R_temp[_i][Cart::x][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::x][m+1] + term;
      R_temp[_i][Cart::xy][m] = wmc0*R_temp[_i][Cart::y][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::y][m+1];
      R_temp[_i][Cart::xz][m] = wmc0*R_temp[_i][Cart::z][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::z][m+1];
      R_temp[_i][Cart::yy][m] = wmc1*R_temp[_i][Cart::y][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::y][m+1] + term;
      R_temp[_i][Cart::yz][m] = wmc1*R_temp[_i][Cart::z][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::z][m+1];
      R_temp[_i][Cart::zz][m] = wmc2*R_temp[_i][Cart::z][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::z][m+1] + term;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 1)


if (_lmax_gamma > 2) {

  //Integral  s - s - f
  for (int m = 0; m < _lmax_gamma-2; m++) {
    R_temp[0][Cart::xxx][m] = wmc0*R_temp[0][Cart::xx][m+1] + 2*rgamma*(R_temp[0][Cart::x][m]-cfak*R_temp[0][Cart::x][m+1]);
    R_temp[0][Cart::xxy][m] = wmc1*R_temp[0][Cart::xx][m+1];
    R_temp[0][Cart::xxz][m] = wmc2*R_temp[0][Cart::xx][m+1];
    R_temp[0][Cart::xyy][m] = wmc0*R_temp[0][Cart::yy][m+1];
    R_temp[0][Cart::xyz][m] = wmc0*R_temp[0][Cart::yz][m+1];
    R_temp[0][Cart::xzz][m] = wmc0*R_temp[0][Cart::zz][m+1];
    R_temp[0][Cart::yyy][m] = wmc1*R_temp[0][Cart::yy][m+1] + 2*rgamma*(R_temp[0][Cart::y][m]-cfak*R_temp[0][Cart::y][m+1]);
    R_temp[0][Cart::yyz][m] = wmc2*R_temp[0][Cart::yy][m+1];
    R_temp[0][Cart::yzz][m] = wmc1*R_temp[0][Cart::zz][m+1];
    R_temp[0][Cart::zzz][m] = wmc2*R_temp[0][Cart::zz][m+1] + 2*rgamma*(R_temp[0][Cart::z][m]-cfak*R_temp[0][Cart::z][m+1]);
  }
  //------------------------------------------------------

  //Integrals     p - s - f     d - s - f     f - s - f     g - s - f     h - s - f     i - s - f     j - s - f     k - s - f
  for (int m = 0; m < _lmax_gamma-2; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      double term_x = 2*rgamma*(R_temp[_i][Cart::x][m]-cfak*R_temp[_i][Cart::x][m+1]);
      double term_y = 2*rgamma*(R_temp[_i][Cart::y][m]-cfak*R_temp[_i][Cart::y][m+1]);
      double term_z = 2*rgamma*(R_temp[_i][Cart::z][m]-cfak*R_temp[_i][Cart::z][m+1]);
      R_temp[_i][Cart::xxx][m] = wmc0*R_temp[_i][Cart::xx][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xx][m+1] + term_x;
      R_temp[_i][Cart::xxy][m] = wmc1*R_temp[_i][Cart::xx][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xx][m+1];
      R_temp[_i][Cart::xxz][m] = wmc2*R_temp[_i][Cart::xx][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xx][m+1];
      R_temp[_i][Cart::xyy][m] = wmc0*R_temp[_i][Cart::yy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yy][m+1];
      R_temp[_i][Cart::xyz][m] = wmc0*R_temp[_i][Cart::yz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yz][m+1];
      R_temp[_i][Cart::xzz][m] = wmc0*R_temp[_i][Cart::zz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::zz][m+1];
      R_temp[_i][Cart::yyy][m] = wmc1*R_temp[_i][Cart::yy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yy][m+1] + term_y;
      R_temp[_i][Cart::yyz][m] = wmc2*R_temp[_i][Cart::yy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yy][m+1];
      R_temp[_i][Cart::yzz][m] = wmc1*R_temp[_i][Cart::zz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::zz][m+1];
      R_temp[_i][Cart::zzz][m] = wmc2*R_temp[_i][Cart::zz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::zz][m+1] + term_z;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 2)


if (_lmax_gamma > 3) {

  //Integral  s - s - g
  for (int m = 0; m < _lmax_gamma-3; m++) {
    double term_xx = rgamma*(R_temp[0][Cart::xx][m]-cfak*R_temp[0][Cart::xx][m+1]);
    double term_yy = rgamma*(R_temp[0][Cart::yy][m]-cfak*R_temp[0][Cart::yy][m+1]);
    double term_zz = rgamma*(R_temp[0][Cart::zz][m]-cfak*R_temp[0][Cart::zz][m+1]);
    R_temp[0][Cart::xxxx][m] = wmc0*R_temp[0][Cart::xxx][m+1] + 3*term_xx;
    R_temp[0][Cart::xxxy][m] = wmc1*R_temp[0][Cart::xxx][m+1];
    R_temp[0][Cart::xxxz][m] = wmc2*R_temp[0][Cart::xxx][m+1];
    R_temp[0][Cart::xxyy][m] = wmc0*R_temp[0][Cart::xyy][m+1] + term_yy;
    R_temp[0][Cart::xxyz][m] = wmc1*R_temp[0][Cart::xxz][m+1];
    R_temp[0][Cart::xxzz][m] = wmc0*R_temp[0][Cart::xzz][m+1] + term_zz;
    R_temp[0][Cart::xyyy][m] = wmc0*R_temp[0][Cart::yyy][m+1];
    R_temp[0][Cart::xyyz][m] = wmc0*R_temp[0][Cart::yyz][m+1];
    R_temp[0][Cart::xyzz][m] = wmc0*R_temp[0][Cart::yzz][m+1];
    R_temp[0][Cart::xzzz][m] = wmc0*R_temp[0][Cart::zzz][m+1];
    R_temp[0][Cart::yyyy][m] = wmc1*R_temp[0][Cart::yyy][m+1] + 3*term_yy;
    R_temp[0][Cart::yyyz][m] = wmc2*R_temp[0][Cart::yyy][m+1];
    R_temp[0][Cart::yyzz][m] = wmc1*R_temp[0][Cart::yzz][m+1] + term_zz;
    R_temp[0][Cart::yzzz][m] = wmc1*R_temp[0][Cart::zzz][m+1];
    R_temp[0][Cart::zzzz][m] = wmc2*R_temp[0][Cart::zzz][m+1] + 3*term_zz;
  }
  //------------------------------------------------------

  //Integrals     p - s - g     d - s - g     f - s - g     g - s - g     h - s - g     i - s - g     j - s - g     k - s - g
  for (int m = 0; m < _lmax_gamma-3; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      double term_xx = rgamma*(R_temp[_i][Cart::xx][m]-cfak*R_temp[_i][Cart::xx][m+1]);
      double term_yy = rgamma*(R_temp[_i][Cart::yy][m]-cfak*R_temp[_i][Cart::yy][m+1]);
      double term_zz = rgamma*(R_temp[_i][Cart::zz][m]-cfak*R_temp[_i][Cart::zz][m+1]);
      R_temp[_i][Cart::xxxx][m] = wmc0*R_temp[_i][Cart::xxx][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xxx][m+1] + 3*term_xx;
      R_temp[_i][Cart::xxxy][m] = wmc1*R_temp[_i][Cart::xxx][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxx][m+1];
      R_temp[_i][Cart::xxxz][m] = wmc2*R_temp[_i][Cart::xxx][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxx][m+1];
      R_temp[_i][Cart::xxyy][m] = wmc0*R_temp[_i][Cart::xyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xyy][m+1] + term_yy;
      R_temp[_i][Cart::xxyz][m] = wmc1*R_temp[_i][Cart::xxz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxz][m+1];
      R_temp[_i][Cart::xxzz][m] = wmc0*R_temp[_i][Cart::xzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xzz][m+1] + term_zz;
      R_temp[_i][Cart::xyyy][m] = wmc0*R_temp[_i][Cart::yyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyy][m+1];
      R_temp[_i][Cart::xyyz][m] = wmc0*R_temp[_i][Cart::yyz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyz][m+1];
      R_temp[_i][Cart::xyzz][m] = wmc0*R_temp[_i][Cart::yzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yzz][m+1];
      R_temp[_i][Cart::xzzz][m] = wmc0*R_temp[_i][Cart::zzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::zzz][m+1];
      R_temp[_i][Cart::yyyy][m] = wmc1*R_temp[_i][Cart::yyy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yyy][m+1] + 3*term_yy;
      R_temp[_i][Cart::yyyz][m] = wmc2*R_temp[_i][Cart::yyy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yyy][m+1];
      R_temp[_i][Cart::yyzz][m] = wmc1*R_temp[_i][Cart::yzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yzz][m+1] + term_zz;
      R_temp[_i][Cart::yzzz][m] = wmc1*R_temp[_i][Cart::zzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::zzz][m+1];
      R_temp[_i][Cart::zzzz][m] = wmc2*R_temp[_i][Cart::zzz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::zzz][m+1] + 3*term_zz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 3)


if (_lmax_gamma > 4) {

  //Integral  s - s - h
  for (int m = 0; m < _lmax_gamma-4; m++) {
    double term_xxx = rgamma*(R_temp[0][Cart::xxx][m]-cfak*R_temp[0][Cart::xxx][m+1]);
    double term_yyy = rgamma*(R_temp[0][Cart::yyy][m]-cfak*R_temp[0][Cart::yyy][m+1]);
    double term_zzz = rgamma*(R_temp[0][Cart::zzz][m]-cfak*R_temp[0][Cart::zzz][m+1]);
    R_temp[0][Cart::xxxxx][m] = wmc0*R_temp[0][Cart::xxxx][m+1] + 4*term_xxx;
    R_temp[0][Cart::xxxxy][m] = wmc1*R_temp[0][Cart::xxxx][m+1];
    R_temp[0][Cart::xxxxz][m] = wmc2*R_temp[0][Cart::xxxx][m+1];
    R_temp[0][Cart::xxxyy][m] = wmc1*R_temp[0][Cart::xxxy][m+1] + term_xxx;
    R_temp[0][Cart::xxxyz][m] = wmc1*R_temp[0][Cart::xxxz][m+1];
    R_temp[0][Cart::xxxzz][m] = wmc2*R_temp[0][Cart::xxxz][m+1] + term_xxx;
    R_temp[0][Cart::xxyyy][m] = wmc0*R_temp[0][Cart::xyyy][m+1] + term_yyy;
    R_temp[0][Cart::xxyyz][m] = wmc2*R_temp[0][Cart::xxyy][m+1];
    R_temp[0][Cart::xxyzz][m] = wmc1*R_temp[0][Cart::xxzz][m+1];
    R_temp[0][Cart::xxzzz][m] = wmc0*R_temp[0][Cart::xzzz][m+1] + term_zzz;
    R_temp[0][Cart::xyyyy][m] = wmc0*R_temp[0][Cart::yyyy][m+1];
    R_temp[0][Cart::xyyyz][m] = wmc0*R_temp[0][Cart::yyyz][m+1];
    R_temp[0][Cart::xyyzz][m] = wmc0*R_temp[0][Cart::yyzz][m+1];
    R_temp[0][Cart::xyzzz][m] = wmc0*R_temp[0][Cart::yzzz][m+1];
    R_temp[0][Cart::xzzzz][m] = wmc0*R_temp[0][Cart::zzzz][m+1];
    R_temp[0][Cart::yyyyy][m] = wmc1*R_temp[0][Cart::yyyy][m+1] + 4*term_yyy;
    R_temp[0][Cart::yyyyz][m] = wmc2*R_temp[0][Cart::yyyy][m+1];
    R_temp[0][Cart::yyyzz][m] = wmc2*R_temp[0][Cart::yyyz][m+1] + term_yyy;
    R_temp[0][Cart::yyzzz][m] = wmc1*R_temp[0][Cart::yzzz][m+1] + term_zzz;
    R_temp[0][Cart::yzzzz][m] = wmc1*R_temp[0][Cart::zzzz][m+1];
    R_temp[0][Cart::zzzzz][m] = wmc2*R_temp[0][Cart::zzzz][m+1] + 4*term_zzz;
  }
  //------------------------------------------------------

  //Integrals     p - s - h     d - s - h     f - s - h     g - s - h     h - s - h     i - s - h     j - s - h     k - s - h
  for (int m = 0; m < _lmax_gamma-4; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      double term_xxx = rgamma*(R_temp[_i][Cart::xxx][m]-cfak*R_temp[_i][Cart::xxx][m+1]);
      double term_yyy = rgamma*(R_temp[_i][Cart::yyy][m]-cfak*R_temp[_i][Cart::yyy][m+1]);
      double term_zzz = rgamma*(R_temp[_i][Cart::zzz][m]-cfak*R_temp[_i][Cart::zzz][m+1]);
      R_temp[_i][Cart::xxxxx][m] = wmc0*R_temp[_i][Cart::xxxx][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xxxx][m+1] + 4*term_xxx;
      R_temp[_i][Cart::xxxxy][m] = wmc1*R_temp[_i][Cart::xxxx][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxx][m+1];
      R_temp[_i][Cart::xxxxz][m] = wmc2*R_temp[_i][Cart::xxxx][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxxx][m+1];
      R_temp[_i][Cart::xxxyy][m] = wmc1*R_temp[_i][Cart::xxxy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxy][m+1] + term_xxx;
      R_temp[_i][Cart::xxxyz][m] = wmc1*R_temp[_i][Cart::xxxz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxz][m+1];
      R_temp[_i][Cart::xxxzz][m] = wmc2*R_temp[_i][Cart::xxxz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxxz][m+1] + term_xxx;
      R_temp[_i][Cart::xxyyy][m] = wmc0*R_temp[_i][Cart::xyyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xyyy][m+1] + term_yyy;
      R_temp[_i][Cart::xxyyz][m] = wmc2*R_temp[_i][Cart::xxyy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxyy][m+1];
      R_temp[_i][Cart::xxyzz][m] = wmc1*R_temp[_i][Cart::xxzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxzz][m+1];
      R_temp[_i][Cart::xxzzz][m] = wmc0*R_temp[_i][Cart::xzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xzzz][m+1] + term_zzz;
      R_temp[_i][Cart::xyyyy][m] = wmc0*R_temp[_i][Cart::yyyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyyy][m+1];
      R_temp[_i][Cart::xyyyz][m] = wmc0*R_temp[_i][Cart::yyyz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyyz][m+1];
      R_temp[_i][Cart::xyyzz][m] = wmc0*R_temp[_i][Cart::yyzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyzz][m+1];
      R_temp[_i][Cart::xyzzz][m] = wmc0*R_temp[_i][Cart::yzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yzzz][m+1];
      R_temp[_i][Cart::xzzzz][m] = wmc0*R_temp[_i][Cart::zzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::zzzz][m+1];
      R_temp[_i][Cart::yyyyy][m] = wmc1*R_temp[_i][Cart::yyyy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yyyy][m+1] + 4*term_yyy;
      R_temp[_i][Cart::yyyyz][m] = wmc2*R_temp[_i][Cart::yyyy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yyyy][m+1];
      R_temp[_i][Cart::yyyzz][m] = wmc2*R_temp[_i][Cart::yyyz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yyyz][m+1] + term_yyy;
      R_temp[_i][Cart::yyzzz][m] = wmc1*R_temp[_i][Cart::yzzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yzzz][m+1] + term_zzz;
      R_temp[_i][Cart::yzzzz][m] = wmc1*R_temp[_i][Cart::zzzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::zzzz][m+1];
      R_temp[_i][Cart::zzzzz][m] = wmc2*R_temp[_i][Cart::zzzz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::zzzz][m+1] + 4*term_zzz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 4)


if (_lmax_gamma > 5) {

  //Integral  s - s - i
  for (int m = 0; m < _lmax_gamma-5; m++) {
    double term_xxxx = rgamma*(R_temp[0][Cart::xxxx][m]-cfak*R_temp[0][Cart::xxxx][m+1]);
    double term_xyyy = rgamma*(R_temp[0][Cart::xyyy][m]-cfak*R_temp[0][Cart::xyyy][m+1]);
    double term_xzzz = rgamma*(R_temp[0][Cart::xzzz][m]-cfak*R_temp[0][Cart::xzzz][m+1]);
    double term_yyyy = rgamma*(R_temp[0][Cart::yyyy][m]-cfak*R_temp[0][Cart::yyyy][m+1]);
    double term_yyzz = rgamma*(R_temp[0][Cart::yyzz][m]-cfak*R_temp[0][Cart::yyzz][m+1]);
    double term_yzzz = rgamma*(R_temp[0][Cart::yzzz][m]-cfak*R_temp[0][Cart::yzzz][m+1]);
    double term_zzzz = rgamma*(R_temp[0][Cart::zzzz][m]-cfak*R_temp[0][Cart::zzzz][m+1]);
    R_temp[0][Cart::xxxxxx][m] = wmc0*R_temp[0][Cart::xxxxx][m+1] + 5*term_xxxx;
    R_temp[0][Cart::xxxxxy][m] = wmc1*R_temp[0][Cart::xxxxx][m+1];
    R_temp[0][Cart::xxxxxz][m] = wmc2*R_temp[0][Cart::xxxxx][m+1];
    R_temp[0][Cart::xxxxyy][m] = wmc1*R_temp[0][Cart::xxxxy][m+1] + term_xxxx;
    R_temp[0][Cart::xxxxyz][m] = wmc1*R_temp[0][Cart::xxxxz][m+1];
    R_temp[0][Cart::xxxxzz][m] = wmc2*R_temp[0][Cart::xxxxz][m+1] + term_xxxx;
    R_temp[0][Cart::xxxyyy][m] = wmc0*R_temp[0][Cart::xxyyy][m+1] + 2*term_xyyy;
    R_temp[0][Cart::xxxyyz][m] = wmc2*R_temp[0][Cart::xxxyy][m+1];
    R_temp[0][Cart::xxxyzz][m] = wmc1*R_temp[0][Cart::xxxzz][m+1];
    R_temp[0][Cart::xxxzzz][m] = wmc0*R_temp[0][Cart::xxzzz][m+1] + 2*term_xzzz;
    R_temp[0][Cart::xxyyyy][m] = wmc0*R_temp[0][Cart::xyyyy][m+1] + term_yyyy;
    R_temp[0][Cart::xxyyyz][m] = wmc2*R_temp[0][Cart::xxyyy][m+1];
    R_temp[0][Cart::xxyyzz][m] = wmc0*R_temp[0][Cart::xyyzz][m+1] + term_yyzz;
    R_temp[0][Cart::xxyzzz][m] = wmc1*R_temp[0][Cart::xxzzz][m+1];
    R_temp[0][Cart::xxzzzz][m] = wmc0*R_temp[0][Cart::xzzzz][m+1] + term_zzzz;
    R_temp[0][Cart::xyyyyy][m] = wmc0*R_temp[0][Cart::yyyyy][m+1];
    R_temp[0][Cart::xyyyyz][m] = wmc0*R_temp[0][Cart::yyyyz][m+1];
    R_temp[0][Cart::xyyyzz][m] = wmc0*R_temp[0][Cart::yyyzz][m+1];
    R_temp[0][Cart::xyyzzz][m] = wmc0*R_temp[0][Cart::yyzzz][m+1];
    R_temp[0][Cart::xyzzzz][m] = wmc0*R_temp[0][Cart::yzzzz][m+1];
    R_temp[0][Cart::xzzzzz][m] = wmc0*R_temp[0][Cart::zzzzz][m+1];
    R_temp[0][Cart::yyyyyy][m] = wmc1*R_temp[0][Cart::yyyyy][m+1] + 5*term_yyyy;
    R_temp[0][Cart::yyyyyz][m] = wmc2*R_temp[0][Cart::yyyyy][m+1];
    R_temp[0][Cart::yyyyzz][m] = wmc2*R_temp[0][Cart::yyyyz][m+1] + term_yyyy;
    R_temp[0][Cart::yyyzzz][m] = wmc1*R_temp[0][Cart::yyzzz][m+1] + 2*term_yzzz;
    R_temp[0][Cart::yyzzzz][m] = wmc1*R_temp[0][Cart::yzzzz][m+1] + term_zzzz;
    R_temp[0][Cart::yzzzzz][m] = wmc1*R_temp[0][Cart::zzzzz][m+1];
    R_temp[0][Cart::zzzzzz][m] = wmc2*R_temp[0][Cart::zzzzz][m+1] + 5*term_zzzz;
  }
  //------------------------------------------------------

  //Integrals     p - s - i     d - s - i     f - s - i     g - s - i     h - s - i     i - s - i     j - s - i     k - s - i
  for (int m = 0; m < _lmax_gamma-5; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
      double term_xxxx = rgamma*(R_temp[_i][Cart::xxxx][m]-cfak*R_temp[_i][Cart::xxxx][m+1]);
      double term_xyyy = rgamma*(R_temp[_i][Cart::xyyy][m]-cfak*R_temp[_i][Cart::xyyy][m+1]);
      double term_xzzz = rgamma*(R_temp[_i][Cart::xzzz][m]-cfak*R_temp[_i][Cart::xzzz][m+1]);
      double term_yyyy = rgamma*(R_temp[_i][Cart::yyyy][m]-cfak*R_temp[_i][Cart::yyyy][m+1]);
      double term_yyzz = rgamma*(R_temp[_i][Cart::yyzz][m]-cfak*R_temp[_i][Cart::yyzz][m+1]);
      double term_yzzz = rgamma*(R_temp[_i][Cart::yzzz][m]-cfak*R_temp[_i][Cart::yzzz][m+1]);
      double term_zzzz = rgamma*(R_temp[_i][Cart::zzzz][m]-cfak*R_temp[_i][Cart::zzzz][m+1]);
      R_temp[_i][Cart::xxxxxx][m] = wmc0*R_temp[_i][Cart::xxxxx][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xxxxx][m+1] + 5*term_xxxx;
      R_temp[_i][Cart::xxxxxy][m] = wmc1*R_temp[_i][Cart::xxxxx][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxxx][m+1];
      R_temp[_i][Cart::xxxxxz][m] = wmc2*R_temp[_i][Cart::xxxxx][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxxxx][m+1];
      R_temp[_i][Cart::xxxxyy][m] = wmc1*R_temp[_i][Cart::xxxxy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxxy][m+1] + term_xxxx;
      R_temp[_i][Cart::xxxxyz][m] = wmc1*R_temp[_i][Cart::xxxxz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxxz][m+1];
      R_temp[_i][Cart::xxxxzz][m] = wmc2*R_temp[_i][Cart::xxxxz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxxxz][m+1] + term_xxxx;
      R_temp[_i][Cart::xxxyyy][m] = wmc0*R_temp[_i][Cart::xxyyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xxyyy][m+1] + 2*term_xyyy;
      R_temp[_i][Cart::xxxyyz][m] = wmc2*R_temp[_i][Cart::xxxyy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxxyy][m+1];
      R_temp[_i][Cart::xxxyzz][m] = wmc1*R_temp[_i][Cart::xxxzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxxzz][m+1];
      R_temp[_i][Cart::xxxzzz][m] = wmc0*R_temp[_i][Cart::xxzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xxzzz][m+1] + 2*term_xzzz;
      R_temp[_i][Cart::xxyyyy][m] = wmc0*R_temp[_i][Cart::xyyyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xyyyy][m+1] + term_yyyy;
      R_temp[_i][Cart::xxyyyz][m] = wmc2*R_temp[_i][Cart::xxyyy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::xxyyy][m+1];
      R_temp[_i][Cart::xxyyzz][m] = wmc0*R_temp[_i][Cart::xyyzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xyyzz][m+1] + term_yyzz;
      R_temp[_i][Cart::xxyzzz][m] = wmc1*R_temp[_i][Cart::xxzzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::xxzzz][m+1];
      R_temp[_i][Cart::xxzzzz][m] = wmc0*R_temp[_i][Cart::xzzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::xzzzz][m+1] + term_zzzz;
      R_temp[_i][Cart::xyyyyy][m] = wmc0*R_temp[_i][Cart::yyyyy][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyyyy][m+1];
      R_temp[_i][Cart::xyyyyz][m] = wmc0*R_temp[_i][Cart::yyyyz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyyyz][m+1];
      R_temp[_i][Cart::xyyyzz][m] = wmc0*R_temp[_i][Cart::yyyzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyyzz][m+1];
      R_temp[_i][Cart::xyyzzz][m] = wmc0*R_temp[_i][Cart::yyzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yyzzz][m+1];
      R_temp[_i][Cart::xyzzzz][m] = wmc0*R_temp[_i][Cart::yzzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::yzzzz][m+1];
      R_temp[_i][Cart::xzzzzz][m] = wmc0*R_temp[_i][Cart::zzzzz][m+1] + nx[_i]*rdecay*R_temp[i_less_x[_i]][Cart::zzzzz][m+1];
      R_temp[_i][Cart::yyyyyy][m] = wmc1*R_temp[_i][Cart::yyyyy][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yyyyy][m+1] + 5*term_yyyy;
      R_temp[_i][Cart::yyyyyz][m] = wmc2*R_temp[_i][Cart::yyyyy][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yyyyy][m+1];
      R_temp[_i][Cart::yyyyzz][m] = wmc2*R_temp[_i][Cart::yyyyz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::yyyyz][m+1] + term_yyyy;
      R_temp[_i][Cart::yyyzzz][m] = wmc1*R_temp[_i][Cart::yyzzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yyzzz][m+1] + 2*term_yzzz;
      R_temp[_i][Cart::yyzzzz][m] = wmc1*R_temp[_i][Cart::yzzzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::yzzzz][m+1] + term_zzzz;
      R_temp[_i][Cart::yzzzzz][m] = wmc1*R_temp[_i][Cart::zzzzz][m+1] + ny[_i]*rdecay*R_temp[i_less_y[_i]][Cart::zzzzz][m+1];
      R_temp[_i][Cart::zzzzzz][m] = wmc2*R_temp[_i][Cart::zzzzz][m+1] + nz[_i]*rdecay*R_temp[i_less_z[_i]][Cart::zzzzz][m+1] + 5*term_zzzz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 5)




const std::vector<double>& _contractions_gamma = itgamma->getContraction();

  // s-functions
double factor = _contractions_gamma[0];
for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
  R_temp[_i][0][1] = factor * R_temp[_i][0][0]; /// Y 0,0
}

if (_lmax_gamma > 0) {
  // p-functions
  factor = 2.*sqrt(_decay_gamma)*_contractions_gamma[1];
  for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
    R_temp[_i][1][1] = factor * R_temp[_i][3][0]; /// Y 1,0
    R_temp[_i][2][1] = factor * R_temp[_i][2][0]; /// Y 1,-1
    R_temp[_i][3][1] = factor * R_temp[_i][1][0]; /// Y 1,1
  }
}

if (_lmax_gamma > 1) {
  // d-functions
  factor = 2.0*_decay_gamma*_contractions_gamma[2];
  double factor_1 =  factor/sqrt(3.0);
  for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
    R_temp[_i][4][1] = factor_1 * ( 2.0*R_temp[_i][Cart::zz][0] - R_temp[_i][Cart::xx][0] - R_temp[_i][Cart::yy][0] );  /// d3z2-r2  Y 2,0
    R_temp[_i][5][1] = 2.*factor * R_temp[_i][Cart::yz][0];  /// dyz  Y 2,-1
    R_temp[_i][6][1] = 2.*factor * R_temp[_i][Cart::xz][0];  /// dxz  Y 2,1
    R_temp[_i][7][1] = 2.*factor * R_temp[_i][Cart::xy][0];  /// dxy  Y 2,-2
    R_temp[_i][8][1] = factor * ( R_temp[_i][Cart::xx][0] - R_temp[_i][Cart::yy][0] );  /// dx2-y2  Y 2,2
  }
}

if (_lmax_gamma > 2) {
  // f-functions
  factor = 2.0*pow(_decay_gamma,1.5)*_contractions_gamma[3];
  double factor_1 = factor*2./sqrt(15.);
  double factor_2 = factor*sqrt(2.)/sqrt(5.);
  double factor_3 = factor*sqrt(2.)/sqrt(3.);
  for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
    R_temp[_i][9][1] = factor_1 * ( 2.*R_temp[_i][Cart::zzz][0] - 3.*R_temp[_i][Cart::xxz][0] - 3.* R_temp[_i][Cart::yyz][0] ); /// Y 3,0
    R_temp[_i][10][1] = factor_2 * ( 4.*R_temp[_i][Cart::yzz][0] - R_temp[_i][Cart::xxy][0] - R_temp[_i][Cart::yyy][0] ); /// Y 3,-1
    R_temp[_i][11][1] = factor_2 * ( 4.*R_temp[_i][Cart::xzz][0] - R_temp[_i][Cart::xxx][0] - R_temp[_i][Cart::xyy][0] ); /// Y 3,1
    R_temp[_i][12][1] = 4.*factor * R_temp[_i][Cart::xyz][0]; /// Y 3,-2
    R_temp[_i][13][1] = 2.*factor * ( R_temp[_i][Cart::xxz][0] - R_temp[_i][Cart::yyz][0] ); /// Y 3,2
    R_temp[_i][14][1] = factor_3 * ( 3.*R_temp[_i][Cart::xxy][0] - R_temp[_i][Cart::yyy][0] ); /// Y 3,-3
    R_temp[_i][15][1] = factor_3 * ( R_temp[_i][Cart::xxx][0] - 3.*R_temp[_i][Cart::xyy][0] ); /// Y 3,3
  }
}

if (_lmax_gamma > 3) {
  // g-functions
  factor = 2./sqrt(3.)*_decay_gamma*_decay_gamma*_contractions_gamma[4];
  double factor_1 = factor/sqrt(35.);
  double factor_2 = factor*4./sqrt(14.);
  double factor_3 = factor*2./sqrt(7.);
  double factor_4 = factor*2.*sqrt(2.);
  for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
    R_temp[_i][16][1] = factor_1 * (    3.*(R_temp[_i][Cart::xxxx][0] + R_temp[_i][Cart::yyyy][0])
                                       + 6.*R_temp[_i][Cart::xxyy][0]
                                     - 24.*(R_temp[_i][Cart::xxzz][0] + R_temp[_i][Cart::yyzz][0])
                                       + 8.*R_temp[_i][Cart::zzzz][0] );                               /// Y 4,0
    R_temp[_i][17][1] = factor_2 * ( -3.*(R_temp[_i][Cart::xxyz][0] + R_temp[_i][Cart::yyyz][0])
                                     + 4.*R_temp[_i][Cart::yzzz][0] );                                 /// Y 4,-1
    R_temp[_i][18][1] = factor_2 * ( -3.*(R_temp[_i][Cart::xxxz][0] + R_temp[_i][Cart::xyyz][0])
                                     + 4.*R_temp[_i][Cart::xzzz][0] );                                 /// Y 4,1
    R_temp[_i][19][1] = 2.*factor_3 * (    -R_temp[_i][Cart::xxxy][0]
                                           - R_temp[_i][Cart::xyyy][0]
                                        + 6.*R_temp[_i][Cart::xyzz][0] );                              /// Y 4,-2
    R_temp[_i][20][1] = factor_3 * (      -R_temp[_i][Cart::xxxx][0]
                                     + 6.*(R_temp[_i][Cart::xxzz][0] - R_temp[_i][Cart::yyzz][0])
                                        + R_temp[_i][Cart::yyyy][0] );                                 /// Y 4,2
    R_temp[_i][21][1] = factor_4 * ( 3.*R_temp[_i][Cart::xxyz][0] 
                                      - R_temp[_i][Cart::yyyz][0] );                                   /// Y 4,-3
    R_temp[_i][22][1] = factor_4 * (      R_temp[_i][Cart::xxxz][0] 
                                     - 3.*R_temp[_i][Cart::xyyz][0] );                                 /// Y 4,3
    R_temp[_i][23][1] = 4.*factor * (   R_temp[_i][Cart::xxxy][0]
                                      - R_temp[_i][Cart::xyyy][0] );                                   /// Y 4,-4
    R_temp[_i][24][1] = factor * (      R_temp[_i][Cart::xxxx][0] 
                                   - 6.*R_temp[_i][Cart::xxyy][0]
                                      + R_temp[_i][Cart::yyyy][0] );                                   /// Y 4,4
  }
}

if (_lmax_gamma > 4) {
  // h-functions
  factor = (2./3.)*pow(_decay_gamma,2.5)*_contractions_gamma[5];
  double factor_1 = factor*2./sqrt(105.);
  double factor_2 = factor*2./sqrt(7.);
  double factor_3 = factor*sqrt(6.)/3.;
  double factor_4 = factor*2.*sqrt(3.);
  double factor_5 = factor*.2*sqrt(30.);
  for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
    R_temp[_i][25][1] = factor_1 * (   15.*(R_temp[_i][Cart::xxxxz][0] + R_temp[_i][Cart::yyyyz][0])
                                      + 30.*R_temp[_i][Cart::xxyyz][0]
                                     - 40.*(R_temp[_i][Cart::xxzzz][0] + R_temp[_i][Cart::yyzzz][0])
                                       + 8.*R_temp[_i][Cart::zzzzz][0] );                              /// Y 5,0

    R_temp[_i][26][1] = factor_2 * (        R_temp[_i][Cart::xxxxy][0]
                                       + 2.*R_temp[_i][Cart::xxyyy][0]
                                     - 12.*(R_temp[_i][Cart::xxyzz][0] + R_temp[_i][Cart::yyyzz][0])
                                          + R_temp[_i][Cart::yyyyy][0]
                                       + 8.*R_temp[_i][Cart::yzzzz][0] );                              /// Y 5,-1

    R_temp[_i][27][1] = factor_2 * (        R_temp[_i][Cart::xxxxx][0]
                                       + 2.*R_temp[_i][Cart::xxxyy][0]
                                     - 12.*(R_temp[_i][Cart::xxxzz][0] + R_temp[_i][Cart::xyyzz][0])
                                          + R_temp[_i][Cart::xyyyy][0]
                                       + 8.*R_temp[_i][Cart::xzzzz][0] );                              /// Y 5,1

    R_temp[_i][28][1] = 8.*factor * (     -R_temp[_i][Cart::xxxyz][0]
                                         - R_temp[_i][Cart::xyyyz][0]
                                      + 2.*R_temp[_i][Cart::xyzzz][0] );                               /// Y 5,-2

    R_temp[_i][29][1] = 4.*factor * (      -R_temp[_i][Cart::xxxxz][0]
                                      + 2.*(R_temp[_i][Cart::xxzzz][0] - R_temp[_i][Cart::yyzzz][0])
                                          + R_temp[_i][Cart::yyyyz][0] );                              /// Y 5,2

    R_temp[_i][30][1] = factor_3 * (   -3.*R_temp[_i][Cart::xxxxy][0]
                                      - 2.*R_temp[_i][Cart::xxyyy][0]
                                     + 24.*R_temp[_i][Cart::xxyzz][0]
                                         + R_temp[_i][Cart::yyyyy][0]
                                      - 8.*R_temp[_i][Cart::yyyzz][0] );                               /// Y 5,-3

    R_temp[_i][31][1] = factor_3 * (      -R_temp[_i][Cart::xxxxx][0]
                                      + 2.*R_temp[_i][Cart::xxxyy][0]
                                      + 8.*R_temp[_i][Cart::xxxzz][0]
                                      + 3.*R_temp[_i][Cart::xyyyy][0]
                                     - 24.*R_temp[_i][Cart::xyyzz][0] );                               /// Y 5,3

    R_temp[_i][32][1] = 4.*factor_4 * (   R_temp[_i][Cart::xxxyz][0]
                                        - R_temp[_i][Cart::xyyyz][0] );                                /// Y 5,-4

    R_temp[_i][33][1] = factor_4 * (      R_temp[_i][Cart::xxxxz][0]
                                     - 6.*R_temp[_i][Cart::xxyyz][0]
                                        + R_temp[_i][Cart::yyyyz][0] );                                /// Y 5,4

    R_temp[_i][34][1] = factor_5 * (    5.*R_temp[_i][Cart::xxxxy][0]
                                     - 10.*R_temp[_i][Cart::xxyyy][0]
                                         + R_temp[_i][Cart::yyyyy][0] );                               /// Y 5,-5

    R_temp[_i][35][1] = factor_5 * (       R_temp[_i][Cart::xxxxx][0]
                                     - 10.*R_temp[_i][Cart::xxxyy][0]
                                      + 5.*R_temp[_i][Cart::xyyyy][0] );                               /// Y 5,5
  }
}


if (_lmax_gamma > 5) {
  // i-functions
  factor = (2./3.)*_decay_gamma*_decay_gamma*_decay_gamma*_contractions_gamma[6];
  double factor_1 = factor*2./sqrt(1155.);
  double factor_2 = factor*4./sqrt(55.);
  double factor_3 = factor*sqrt(22.)/11.;
  double factor_4 = factor*2.*sqrt(165.)/55.;
  double factor_5 = factor*.4*sqrt(30.);
  double factor_6 = factor*.2*sqrt(10.);
  for (int _i =  0; _i < n_orbitals[_lmax_alpha_beta]; _i++) {
    R_temp[_i][36][1] = factor_1 * (    -5.*(R_temp[_i][Cart::xxxxxx][0] + R_temp[_i][Cart::yyyyyy][0])
                                      - 15.*(R_temp[_i][Cart::xxxxyy][0] + R_temp[_i][Cart::xxyyyy][0])
                                      + 90.*(R_temp[_i][Cart::xxxxzz][0] + R_temp[_i][Cart::yyyyzz][0])
                                      + 180.*R_temp[_i][Cart::xxyyzz][0]
                                     - 120.*(R_temp[_i][Cart::xxzzzz][0] + R_temp[_i][Cart::yyzzzz][0])
                                       + 16.*R_temp[_i][Cart::zzzzzz][0] );                                /// Y 6,0

    R_temp[_i][37][1] = factor_2 * (    5.*(R_temp[_i][Cart::xxxxyz][0] + R_temp[_i][Cart::yyyyyz][0])
                                      + 10.*R_temp[_i][Cart::xxyyyz][0]
                                     - 20.*(R_temp[_i][Cart::xxyzzz][0] + R_temp[_i][Cart::yyyzzz][0])
                                       + 8.*R_temp[_i][Cart::yzzzzz][0] );                                 /// Y 6,-1

    R_temp[_i][38][1] = factor_2 * (    5.*(R_temp[_i][Cart::xxxxxz][0] + R_temp[_i][Cart::xyyyyz][0])
                                      + 10.*R_temp[_i][Cart::xxxyyz][0]
                                     - 20.*(R_temp[_i][Cart::xxxzzz][0] + R_temp[_i][Cart::xyyzzz][0])
                                       + 8.*R_temp[_i][Cart::xzzzzz][0] );                                 /// Y 6,1

    R_temp[_i][39][1] = 2.*factor_3 * (        R_temp[_i][Cart::xxxxxy][0]
                                          + 2.*R_temp[_i][Cart::xxxyyy][0]
                                        - 16.*(R_temp[_i][Cart::xxxyzz][0] + R_temp[_i][Cart::xyyyzz][0] - R_temp[_i][Cart::xyzzzz][0])
                                             + R_temp[_i][Cart::xyyyyy][0] );                              /// Y 6,-2

    R_temp[_i][40][1] = factor_3 * (        R_temp[_i][Cart::xxxxxy][0]
                                          + R_temp[_i][Cart::xxxxyy][0]
                                     - 16.*(R_temp[_i][Cart::xxxxzz][0] - R_temp[_i][Cart::xxzzzz][0]
                                                                        - R_temp[_i][Cart::yyyyzz][0] + R_temp[_i][Cart::yyzzzz][0])
                                          - R_temp[_i][Cart::xxyyyy][0]
                                          - R_temp[_i][Cart::yyyyyy][0] );                                 /// Y 6,2

    R_temp[_i][41][1] = 2.*factor_3 * (   -9.*R_temp[_i][Cart::xxxxyz][0]
                                         - 6.*R_temp[_i][Cart::xxyyyz][0]
                                        + 24.*R_temp[_i][Cart::xxyzzz][0]
                                         + 3.*R_temp[_i][Cart::yyyyyz][0]
                                         - 8.*R_temp[_i][Cart::yyyzzz][0] );                               /// Y 6,-3

    R_temp[_i][42][1] = 2.*factor_3 * (   -3.*R_temp[_i][Cart::xxxxxz][0]
                                         + 6.*R_temp[_i][Cart::xxxyyz][0]
                                         + 8.*R_temp[_i][Cart::xxxzzz][0]
                                         + 9.*R_temp[_i][Cart::xyyyyz][0]
                                        - 24.*R_temp[_i][Cart::xyyzzz][0] );                               /// Y 6,3

    R_temp[_i][43][1] = 4.*factor_4 * (       -R_temp[_i][Cart::xxxxxy][0]
                                        + 10.*(R_temp[_i][Cart::xxxyzz][0] - R_temp[_i][Cart::xyyyzz][0])
                                             + R_temp[_i][Cart::xyyyyy][0] );                              /// Y 6,-4

    R_temp[_i][44][1] = factor_4 * (       -R_temp[_i][Cart::xxxxxx][0]
                                      + 5.*(R_temp[_i][Cart::xxxxyy][0] + R_temp[_i][Cart::xxyyyy][0])
                                     + 10.*(R_temp[_i][Cart::xxxxzz][0] + R_temp[_i][Cart::yyyyzz][0])
                                      - 60.*R_temp[_i][Cart::xxyyzz][0]
                                         -  R_temp[_i][Cart::yyyyyy][0] );                                 /// Y 6,4

    R_temp[_i][45][1] = factor_5 * (    5.*R_temp[_i][Cart::xxxxyz][0]
                                     - 10.*R_temp[_i][Cart::xxyyyz][0]
                                         + R_temp[_i][Cart::yyyyyz][0] );                                  /// Y 6,-5

    R_temp[_i][46][1] = factor_5 * (       R_temp[_i][Cart::xxxxxz][0]
                                     - 10.*R_temp[_i][Cart::xxxyyz][0]
                                      + 5.*R_temp[_i][Cart::xyyyyz][0] );                                  /// Y 6,5

    R_temp[_i][47][1] = 2.*factor_6 * (    3.*R_temp[_i][Cart::xxxxxy][0]
                                        - 10.*R_temp[_i][Cart::xxxyyy][0]
                                         + 3.*R_temp[_i][Cart::xyyyyy][0] );                               /// Y 6,-6

    R_temp[_i][48][1] = factor_6 * (        R_temp[_i][Cart::xxxxxx][0]
                                     - 15.*(R_temp[_i][Cart::xxxxyy][0] - R_temp[_i][Cart::xxyyyy][0])
                                          - R_temp[_i][Cart::yyyyyy][0] );                                 /// Y 6,6

  }
}



//copy into new array for 3D use.

for (index3d i = 0; i < n_orbitals[_lmax_alpha_beta]; ++i) {

         for (index3d k = 0; k < (_lmax_gamma+1)*(_lmax_gamma+1); ++k) {

                            R[i][0][k] = R_temp[i][k][1];
                        }

                }

            



if (_lmax_beta > 0) {
  //Integrals    s - p - *    p - p - *    d - p - *    f - p - *    g - p - *    h - p - *    i - p - *    j - p - *
  for (int _i = 0; _i < (_lmax_gamma+1)*(_lmax_gamma+1); _i++) {
    for (int _j = 0; _j < n_orbitals[_lmax_alpha_beta-1]; _j++) {
      R[_j][Cart::x][_i] = R[i_more_x[_j]][0][_i] + amb0*R[_j][0][_i];
      R[_j][Cart::y][_i] = R[i_more_y[_j]][0][_i] + amb1*R[_j][0][_i];
      R[_j][Cart::z][_i] = R[i_more_z[_j]][0][_i] + amb2*R[_j][0][_i];
    }
  }
  //------------------------------------------------------
}

if (_lmax_beta > 1) {
  //Integrals    s - d - *    p - d - *    d - d - *    f - d - *    g - d - *    h - d - *    i - d - *
  for (int _i = 0; _i < (_lmax_gamma+1)*(_lmax_gamma+1); _i++) {
    for (int _j = 0; _j < n_orbitals[_lmax_alpha_beta-2]; _j++) {
      R[_j][Cart::xx][_i] = R[i_more_x[_j]][Cart::x][_i] + amb0*R[_j][Cart::x][_i];
      R[_j][Cart::xy][_i] = R[i_more_x[_j]][Cart::y][_i] + amb0*R[_j][Cart::y][_i];
      R[_j][Cart::xz][_i] = R[i_more_x[_j]][Cart::z][_i] + amb0*R[_j][Cart::z][_i];
      R[_j][Cart::yy][_i] = R[i_more_y[_j]][Cart::y][_i] + amb1*R[_j][Cart::y][_i];
      R[_j][Cart::yz][_i] = R[i_more_y[_j]][Cart::z][_i] + amb1*R[_j][Cart::z][_i];
      R[_j][Cart::zz][_i] = R[i_more_z[_j]][Cart::z][_i] + amb2*R[_j][Cart::z][_i];
    }
  }
  //------------------------------------------------------
}

if (_lmax_beta > 2) {
  //Integrals    s - f - *    p - f - *    d - f - *    f - f - *    g - f - *    h - f - *
  for (int _i = 0; _i < (_lmax_gamma+1)*(_lmax_gamma+1); _i++) {
    for (int _j = 0; _j < n_orbitals[_lmax_alpha_beta-3]; _j++) {
      R[_j][Cart::xxx][_i] = R[i_more_x[_j]][Cart::xx][_i] + amb0*R[_j][Cart::xx][_i];
      R[_j][Cart::xxy][_i] = R[i_more_x[_j]][Cart::xy][_i] + amb0*R[_j][Cart::xy][_i];
      R[_j][Cart::xxz][_i] = R[i_more_x[_j]][Cart::xz][_i] + amb0*R[_j][Cart::xz][_i];
      R[_j][Cart::xyy][_i] = R[i_more_x[_j]][Cart::yy][_i] + amb0*R[_j][Cart::yy][_i];
      R[_j][Cart::xyz][_i] = R[i_more_x[_j]][Cart::yz][_i] + amb0*R[_j][Cart::yz][_i];
      R[_j][Cart::xzz][_i] = R[i_more_x[_j]][Cart::zz][_i] + amb0*R[_j][Cart::zz][_i];
      R[_j][Cart::yyy][_i] = R[i_more_y[_j]][Cart::yy][_i] + amb1*R[_j][Cart::yy][_i];
      R[_j][Cart::yyz][_i] = R[i_more_y[_j]][Cart::yz][_i] + amb1*R[_j][Cart::yz][_i];
      R[_j][Cart::yzz][_i] = R[i_more_y[_j]][Cart::zz][_i] + amb1*R[_j][Cart::zz][_i];
      R[_j][Cart::zzz][_i] = R[i_more_z[_j]][Cart::zz][_i] + amb2*R[_j][Cart::zz][_i];
    }
  }
  //------------------------------------------------------
}

if (_lmax_beta > 3) {
  //Integrals    s - g - *    p - g - *    d - g - *    f - g - *    g - g - *
  for (int _i = 0; _i < (_lmax_gamma+1)*(_lmax_gamma+1); _i++) {
    for (int _j = 0; _j < n_orbitals[_lmax_alpha_beta-4]; _j++) {
      R[_j][Cart::xxxx][_i] = R[i_more_x[_j]][Cart::xxx][_i] + amb0*R[_j][Cart::xxx][_i];
      R[_j][Cart::xxxy][_i] = R[i_more_x[_j]][Cart::xxy][_i] + amb0*R[_j][Cart::xxy][_i];
      R[_j][Cart::xxxz][_i] = R[i_more_x[_j]][Cart::xxz][_i] + amb0*R[_j][Cart::xxz][_i];
      R[_j][Cart::xxyy][_i] = R[i_more_x[_j]][Cart::xyy][_i] + amb0*R[_j][Cart::xyy][_i];
      R[_j][Cart::xxyz][_i] = R[i_more_x[_j]][Cart::xyz][_i] + amb0*R[_j][Cart::xyz][_i];
      R[_j][Cart::xxzz][_i] = R[i_more_x[_j]][Cart::xzz][_i] + amb0*R[_j][Cart::xzz][_i];
      R[_j][Cart::xyyy][_i] = R[i_more_x[_j]][Cart::yyy][_i] + amb0*R[_j][Cart::yyy][_i];
      R[_j][Cart::xyyz][_i] = R[i_more_x[_j]][Cart::yyz][_i] + amb0*R[_j][Cart::yyz][_i];
      R[_j][Cart::xyzz][_i] = R[i_more_x[_j]][Cart::yzz][_i] + amb0*R[_j][Cart::yzz][_i];
      R[_j][Cart::xzzz][_i] = R[i_more_x[_j]][Cart::zzz][_i] + amb0*R[_j][Cart::zzz][_i];
      R[_j][Cart::yyyy][_i] = R[i_more_y[_j]][Cart::yyy][_i] + amb1*R[_j][Cart::yyy][_i];
      R[_j][Cart::yyyz][_i] = R[i_more_y[_j]][Cart::yyz][_i] + amb1*R[_j][Cart::yyz][_i];
      R[_j][Cart::yyzz][_i] = R[i_more_y[_j]][Cart::yzz][_i] + amb1*R[_j][Cart::yzz][_i];
      R[_j][Cart::yzzz][_i] = R[i_more_y[_j]][Cart::zzz][_i] + amb1*R[_j][Cart::zzz][_i];
      R[_j][Cart::zzzz][_i] = R[i_more_z[_j]][Cart::zzz][_i] + amb2*R[_j][Cart::zzz][_i];
    }
  }
  //------------------------------------------------------
}




            int istart[] = {0, 1, 1, 1, 4, 4, 4, 4, 4, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20 };
            int istop[] =  {0, 3, 3, 3, 9, 9, 9, 9, 9, 19, 19, 19, 19, 19, 19, 19, 34, 34, 34, 34, 34, 34, 34, 34, 34 };

        
            // which ones do we want to store
            int _offset_beta = _shell_beta->getOffset();
            int _offset_alpha = _shell_alpha->getOffset();
            int _offset_gamma = _shell_gamma->getOffset();

            const Eigen::MatrixXd _trafo_beta=AOSuperMatrix::getTrafo(*itbeta);
            const Eigen::MatrixXd _trafo_alpha=AOSuperMatrix::getTrafo(*italpha);
            
       

            if (alphabetaswitch == true) {

              for (int _i_alpha = 0; _i_alpha < _shell_alpha->getNumFunc(); _i_alpha++) {
                int _i_alpha_off = _i_alpha + _offset_alpha;
                for (int _i_beta = 0; _i_beta < _shell_beta->getNumFunc(); _i_beta++) {
                  int _i_beta_off = _i_beta + _offset_beta;
                  for (int _i_gamma = 0; _i_gamma < _shell_gamma->getNumFunc(); _i_gamma++) {
                    int _i_gamma_off = _i_gamma + _offset_gamma;
                    for (int _i_beta_t = istart[ _i_beta_off ]; _i_beta_t <= istop[ _i_beta_off ]; _i_beta_t++) {
                      for (int _i_alpha_t = istart[ _i_alpha_off ]; _i_alpha_t <= istop[_i_alpha_off ]; _i_alpha_t++) {
                        threec_block[_i_gamma][_i_beta][_i_alpha] += R[ _i_alpha_t ][ _i_beta_t][ _i_gamma_off] * _trafo_alpha(_i_alpha_t, _i_alpha_off) * _trafo_beta(_i_beta_t, _i_beta_off);
                      }
                    }
                  }
                }
              }
            }
            else {
              for (int _i_alpha = 0; _i_alpha < _shell_alpha->getNumFunc(); _i_alpha++) {
                int _i_alpha_off = _i_alpha + _offset_alpha;
                for (int _i_beta = 0; _i_beta < _shell_beta->getNumFunc(); _i_beta++) {
                  int _i_beta_off = _i_beta + _offset_beta;
                  for (int _i_gamma = 0; _i_gamma < _shell_gamma->getNumFunc(); _i_gamma++) {
                    int _i_gamma_off = _i_gamma + _offset_gamma;
                    for (int _i_beta_t = istart[ _i_beta_off ]; _i_beta_t <= istop[ _i_beta_off ]; _i_beta_t++) {
                      for (int _i_alpha_t = istart[ _i_alpha_off ]; _i_alpha_t <= istop[_i_alpha_off ]; _i_alpha_t++) {
                        threec_block[_i_gamma][_i_alpha][_i_beta] += R[ _i_alpha_t ][ _i_beta_t][ _i_gamma_off] * _trafo_alpha(_i_alpha_t, _i_alpha_off) * _trafo_beta(_i_beta_t, _i_beta_off);
                      }
                    }
                  }
                }
              }


            }

                }
            }
        }

 
    
       return _does_contribute;     
    }  
        
        
        
        
        
       
    }}
