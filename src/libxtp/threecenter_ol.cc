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

namespace votca {
    namespace xtp {
                
        /*
         * Calculate 3-center overlap integrals 
         *    S_{abc} = int{ phi_a^DFT phi_b^GW phi_c^DFT d3r }
         * for a given set of a b c as in 
         *    Obara, Saika, J. Chem. Phys. 84, 3963 (1986)
         * section II.B for cartesian Gaussians, then transforming
         * to spherical (angular momentum) Gaussians ("complete" shells
         * from S to Lmax, and finally cutting out those angular momentum 
         * components actually present in shell-shell-shell combination.
         * Currently supported for 
         *      S,P,D,F,G,H,I   functions in DFT basis and 
         *      S,P,D,F,G,H,I   functions in GW  basis
         * 
         */
        bool TCMatrix::FillThreeCenterOLBlock(Eigen::MatrixXd& _subvector,const AOShell* _shell_gw,
                const AOShell* _shell_alpha,const AOShell* _shell_gamma) {
	  
            const double pi = boost::math::constants::pi<double>();
              // get shell positions
            const tools::vec& _pos_gw = _shell_gw->getPos();
            const tools::vec& _pos_alpha = _shell_alpha->getPos();
            const tools::vec& _pos_gamma = _shell_gamma->getPos();

            // shell info, only lmax tells how far to go
            int _lmax_gw = _shell_gw->getLmax();
            int _lmax_alpha = _shell_alpha->getLmax();
            int _lmax_gamma = _shell_gamma->getLmax();

            // set size of internal block for recursion
            int _ngw = AOSuperMatrix::getBlockSize(_lmax_gw);
            int _nalpha = AOSuperMatrix::getBlockSize(_lmax_alpha);
            int _ngamma = AOSuperMatrix::getBlockSize(_lmax_gamma);
            // definition of cutoff for contribution
            
            const double gwaccuracy = 1.e-9; // should become an OPTION
            
            bool _does_contribute = false;



 int nx[] = { 0,
              1, 0, 0,
              2, 1, 1, 0, 0, 0,
              3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
              4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0,
              5, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
              6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 };

 int ny[] = { 0,
              0, 1, 0,
              0, 1, 0, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0 };

 int nz[] = { 0,
              0, 0, 1,
              0, 0, 1, 0, 1, 2,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6 };


 int i_less_x[] = {  0,
                     0,  0,  0,
                     1,  2,  3,  0,  0,  0,
                     4,  5,  6,  7,  8,  9,  0,  0,  0,  0,
                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  0,  0,  0,  0,  0,
                    20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,  0,  0,  0,  0,  0,  0,
                    35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,  0,  0,  0,  0,  0,  0,  0 };

 int i_less_y[] = { 0,
                    0,  0,  0,
                    0,  1,  0,  2,  3,  0,
                    0,  4,  0,  5,  6,  0,  7,  8,  9,  0,
                    0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,  0,
                    0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,  0,
                    0, 35,  0, 36, 37,  0, 38, 39, 40,  0, 41, 42, 43, 44,  0, 45, 46, 47, 48, 49,  0, 50, 51, 52, 53, 54, 55,  0 };

 int i_less_z[] = { 0,
                    0,  0,  0,
                    0,  0,  1,  0,  2,  3,
                    0,  0,  4,  0,  5,  6,  0,  7,  8,  9,
                    0,  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,
                    0,  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,
                    0,  0, 35,  0, 36, 37,  0, 38, 39, 40,  0, 41, 42, 43, 44,  0, 45, 46, 47, 48, 49,  0, 50, 51, 52, 53, 54, 55 };



         
                // iterate over Gaussians in this _shell_row
            for ( AOShell::GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha){
            // iterate over Gaussians in this _shell_col
                const double _decay_alpha = italpha->getDecay();
            
                for ( AOShell::GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma){
                    const double _decay_gamma = itgamma->getDecay();
                    // check third threshold
                    tools::vec _diff = _pos_alpha - _pos_gamma;
                    double test = _decay_alpha * _decay_gamma * _diff*_diff;
                    
                    for ( AOShell::GaussianIterator itgw = _shell_gw->firstGaussian(); itgw != _shell_gw->lastGaussian(); ++itgw){
            // get decay constants (this all is still valid only for uncontracted functions)
                        const double _decay_gw = itgw->getDecay();
            
 
            double threshold = -(_decay_alpha + _decay_gamma + _decay_gw) * log(gwaccuracy);
            // check third threshold
            if (test > threshold) { continue; }
            // check first threshold
            _diff = _pos_alpha - _pos_gw;
            
            test += _decay_alpha * _decay_gw * _diff*_diff;
            if (test > threshold) { continue; }

            // check second threshold
            _diff = _pos_gamma - _pos_gw;
            test += _decay_gamma * _decay_gw * _diff*_diff;
            if (test > threshold) { continue; }

            
            
            

            // if all threshold test are passed, start evaluating

            // some helpers
            double fak = 0.5 / (_decay_alpha + _decay_gw + _decay_gamma);
            double fak2 = 2.0 * fak;
            
            double expo = _decay_alpha * _decay_gamma * (_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma)
                    + _decay_gamma * _decay_gw * (_pos_gamma - _pos_gw) * (_pos_gamma - _pos_gw)
                    + _decay_alpha * _decay_gw * (_pos_alpha - _pos_gw) * (_pos_alpha - _pos_gw);

            
            double prefak = pow(8.0 * _decay_alpha * _decay_gamma * _decay_gw / pi, 0.75) * pow(fak2, 1.5);

            double value = prefak * exp(-fak2 * expo);

            // check if it contributes
            if (value < gwaccuracy) { continue; }

            //double _dist1=(_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma);
            //double _dist2=(_pos_gamma - _pos_gw) * (_pos_gamma - _pos_gw);
            //double _dist3=(_pos_alpha - _pos_gw) * (_pos_alpha - _pos_gw);
            
                  
            tools::vec gvv = fak2 * (_decay_alpha * _pos_alpha + _decay_gw * _pos_gw + _decay_gamma * _pos_gamma);
            tools::vec gma = gvv - _pos_alpha;
            tools::vec gmb = gvv - _pos_gamma;
            tools::vec gmc = gvv - _pos_gw; 
            
            double gma0 = gma.getX();
            double gmb0 = gmb.getX();
            double gmc0 = gmc.getX();

            double gma1 = gma.getY();
            double gmb1 = gmb.getY();
            double gmc1 = gmc.getY();

            double gma2 = gma.getZ();
            double gmb2 = gmb.getZ();
            double gmc2 = gmc.getZ();

            
            // get s-s-s element
            

            _does_contribute = true;
            // if it does, go on and create multiarray
           
            tensor3d::extent_gen extents;
            tensor3d S;
////////            S.resize(extents[ range(1, _nalpha + 1) ][ range(1, _ngw + 1) ][ range(1, _ngamma + 1)]);
            S.resize(extents[ range(0, _nalpha) ][ range(0, _ngw) ][ range(0, _ngamma) ]); /////////////////

            //cout << S.shape()[0]<< " : "<< S.shape()[1]<< " : "<< S.shape()[2]<<endl;


//Integral     s - s - s
S[0][0][0] = value;
//------------------------------------------------------

//Integrals     s - p - s
if (_lmax_gw > 0) {
  S[0][Cart::x][0] = gmc0*S[0][0][0];
  S[0][Cart::y][0] = gmc1*S[0][0][0];
  S[0][Cart::z][0] = gmc2*S[0][0][0];
}
//------------------------------------------------------

//Integrals     s - d - s
if (_lmax_gw > 1) {
  double term = fak*S[0][0][0];
  S[0][Cart::xx][0] = gmc0*S[0][Cart::x][0] + term;
  S[0][Cart::xy][0] = gmc0*S[0][Cart::y][0];
  S[0][Cart::xz][0] = gmc0*S[0][Cart::z][0];
  S[0][Cart::yy][0] = gmc1*S[0][Cart::y][0] + term;
  S[0][Cart::yz][0] = gmc1*S[0][Cart::z][0];
  S[0][Cart::zz][0] = gmc2*S[0][Cart::z][0] + term;
}
//------------------------------------------------------

//Integrals     s - f - s
if (_lmax_gw > 2) {
  S[0][Cart::xxx][0] = gmc0*S[0][Cart::xx][0] + 2*fak*S[0][Cart::x][0];
  S[0][Cart::xxy][0] = gmc1*S[0][Cart::xx][0];
  S[0][Cart::xxz][0] = gmc2*S[0][Cart::xx][0];
  S[0][Cart::xyy][0] = gmc0*S[0][Cart::yy][0];
  S[0][Cart::xyz][0] = gmc0*S[0][Cart::yz][0];
  S[0][Cart::xzz][0] = gmc0*S[0][Cart::zz][0];
  S[0][Cart::yyy][0] = gmc1*S[0][Cart::yy][0] + 2*fak*S[0][Cart::y][0];
  S[0][Cart::yyz][0] = gmc2*S[0][Cart::yy][0];
  S[0][Cart::yzz][0] = gmc1*S[0][Cart::zz][0];
  S[0][Cart::zzz][0] = gmc2*S[0][Cart::zz][0] + 2*fak*S[0][Cart::z][0];
}
//------------------------------------------------------

//Integrals     s - g - s
if (_lmax_gw > 3) {
  double term_xx = fak*S[0][Cart::xx][0];
  double term_yy = fak*S[0][Cart::yy][0];
  double term_zz = fak*S[0][Cart::zz][0];
  S[0][Cart::xxxx][0] = gmc0*S[0][Cart::xxx][0] + 3*term_xx;
  S[0][Cart::xxxy][0] = gmc1*S[0][Cart::xxx][0];
  S[0][Cart::xxxz][0] = gmc2*S[0][Cart::xxx][0];
  S[0][Cart::xxyy][0] = gmc0*S[0][Cart::xyy][0] + term_yy;
  S[0][Cart::xxyz][0] = gmc1*S[0][Cart::xxz][0];
  S[0][Cart::xxzz][0] = gmc0*S[0][Cart::xzz][0] + term_zz;
  S[0][Cart::xyyy][0] = gmc0*S[0][Cart::yyy][0];
  S[0][Cart::xyyz][0] = gmc0*S[0][Cart::yyz][0];
  S[0][Cart::xyzz][0] = gmc0*S[0][Cart::yzz][0];
  S[0][Cart::xzzz][0] = gmc0*S[0][Cart::zzz][0];
  S[0][Cart::yyyy][0] = gmc1*S[0][Cart::yyy][0] + 3*term_yy;
  S[0][Cart::yyyz][0] = gmc2*S[0][Cart::yyy][0];
  S[0][Cart::yyzz][0] = gmc1*S[0][Cart::yzz][0] + term_zz;
  S[0][Cart::yzzz][0] = gmc1*S[0][Cart::zzz][0];
  S[0][Cart::zzzz][0] = gmc2*S[0][Cart::zzz][0] + 3*term_zz;
}
//------------------------------------------------------

//Integrals     s - h - s
if (_lmax_gw > 4) {
  double term_xxx = fak*S[0][Cart::xxx][0];
  double term_yyy = fak*S[0][Cart::yyy][0];
  double term_zzz = fak*S[0][Cart::zzz][0];
  S[0][Cart::xxxxx][0] = gmc0*S[0][Cart::xxxx][0] + 4*term_xxx;
  S[0][Cart::xxxxy][0] = gmc1*S[0][Cart::xxxx][0];
  S[0][Cart::xxxxz][0] = gmc2*S[0][Cart::xxxx][0];
  S[0][Cart::xxxyy][0] = gmc1*S[0][Cart::xxxy][0] + term_xxx;
  S[0][Cart::xxxyz][0] = gmc1*S[0][Cart::xxxz][0];
  S[0][Cart::xxxzz][0] = gmc2*S[0][Cart::xxxz][0] + term_xxx;
  S[0][Cart::xxyyy][0] = gmc0*S[0][Cart::xyyy][0] + term_yyy;
  S[0][Cart::xxyyz][0] = gmc2*S[0][Cart::xxyy][0];
  S[0][Cart::xxyzz][0] = gmc1*S[0][Cart::xxzz][0];
  S[0][Cart::xxzzz][0] = gmc0*S[0][Cart::xzzz][0] + term_zzz;
  S[0][Cart::xyyyy][0] = gmc0*S[0][Cart::yyyy][0];
  S[0][Cart::xyyyz][0] = gmc0*S[0][Cart::yyyz][0];
  S[0][Cart::xyyzz][0] = gmc0*S[0][Cart::yyzz][0];
  S[0][Cart::xyzzz][0] = gmc0*S[0][Cart::yzzz][0];
  S[0][Cart::xzzzz][0] = gmc0*S[0][Cart::zzzz][0];
  S[0][Cart::yyyyy][0] = gmc1*S[0][Cart::yyyy][0] + 4*term_yyy;
  S[0][Cart::yyyyz][0] = gmc2*S[0][Cart::yyyy][0];
  S[0][Cart::yyyzz][0] = gmc2*S[0][Cart::yyyz][0] + term_yyy;
  S[0][Cart::yyzzz][0] = gmc1*S[0][Cart::yzzz][0] + term_zzz;
  S[0][Cart::yzzzz][0] = gmc1*S[0][Cart::zzzz][0];
  S[0][Cart::zzzzz][0] = gmc2*S[0][Cart::zzzz][0] + 4*term_zzz;
}
//------------------------------------------------------

//Integrals     s - i - s
if (_lmax_gw > 5) {
  double term_xxxx = fak*S[0][Cart::xxxx][0];
  double term_xyyy = fak*S[0][Cart::xyyy][0];
  double term_xzzz = fak*S[0][Cart::xzzz][0];
  double term_yyyy = fak*S[0][Cart::yyyy][0];
  double term_yyzz = fak*S[0][Cart::yyzz][0];
  double term_yzzz = fak*S[0][Cart::yzzz][0];
  double term_zzzz = fak*S[0][Cart::zzzz][0];
  S[0][Cart::xxxxxx][0] = gmc0*S[0][Cart::xxxxx][0] + 5*term_xxxx;
  S[0][Cart::xxxxxy][0] = gmc1*S[0][Cart::xxxxx][0];
  S[0][Cart::xxxxxz][0] = gmc2*S[0][Cart::xxxxx][0];
  S[0][Cart::xxxxyy][0] = gmc1*S[0][Cart::xxxxy][0] + term_xxxx;
  S[0][Cart::xxxxyz][0] = gmc1*S[0][Cart::xxxxz][0];
  S[0][Cart::xxxxzz][0] = gmc2*S[0][Cart::xxxxz][0] + term_xxxx;
  S[0][Cart::xxxyyy][0] = gmc0*S[0][Cart::xxyyy][0] + 2*term_xyyy;
  S[0][Cart::xxxyyz][0] = gmc2*S[0][Cart::xxxyy][0];
  S[0][Cart::xxxyzz][0] = gmc1*S[0][Cart::xxxzz][0];
  S[0][Cart::xxxzzz][0] = gmc0*S[0][Cart::xxzzz][0] + 2*term_xzzz;
  S[0][Cart::xxyyyy][0] = gmc0*S[0][Cart::xyyyy][0] + term_yyyy;
  S[0][Cart::xxyyyz][0] = gmc2*S[0][Cart::xxyyy][0];
  S[0][Cart::xxyyzz][0] = gmc0*S[0][Cart::xyyzz][0] + term_yyzz;
  S[0][Cart::xxyzzz][0] = gmc1*S[0][Cart::xxzzz][0];
  S[0][Cart::xxzzzz][0] = gmc0*S[0][Cart::xzzzz][0] + term_zzzz;
  S[0][Cart::xyyyyy][0] = gmc0*S[0][Cart::yyyyy][0];
  S[0][Cart::xyyyyz][0] = gmc0*S[0][Cart::yyyyz][0];
  S[0][Cart::xyyyzz][0] = gmc0*S[0][Cart::yyyzz][0];
  S[0][Cart::xyyzzz][0] = gmc0*S[0][Cart::yyzzz][0];
  S[0][Cart::xyzzzz][0] = gmc0*S[0][Cart::yzzzz][0];
  S[0][Cart::xzzzzz][0] = gmc0*S[0][Cart::zzzzz][0];
  S[0][Cart::yyyyyy][0] = gmc1*S[0][Cart::yyyyy][0] + 5*term_yyyy;
  S[0][Cart::yyyyyz][0] = gmc2*S[0][Cart::yyyyy][0];
  S[0][Cart::yyyyzz][0] = gmc2*S[0][Cart::yyyyz][0] + term_yyyy;
  S[0][Cart::yyyzzz][0] = gmc1*S[0][Cart::yyzzz][0] + 2*term_yzzz;
  S[0][Cart::yyzzzz][0] = gmc1*S[0][Cart::yzzzz][0] + term_zzzz;
  S[0][Cart::yzzzzz][0] = gmc1*S[0][Cart::zzzzz][0];
  S[0][Cart::zzzzzz][0] = gmc2*S[0][Cart::zzzzz][0] + 5*term_zzzz;
}
//------------------------------------------------------







if (_lmax_alpha > 0) {

  //Integrals     p - s - s
  S[Cart::x][0][0] = gma0*S[0][0][0];
  S[Cart::y][0][0] = gma1*S[0][0][0];
  S[Cart::z][0][0] = gma2*S[0][0][0];
  //------------------------------------------------------

  //Integrals     p - p - s
  if (_lmax_gw > 0) {
    double term = fak*S[0][0][0];
    for (int _i =  1; _i < 4; _i++) {
      S[Cart::x][_i][0] = gma0*S[0][_i][0] + nx[_i]*term;
      S[Cart::y][_i][0] = gma1*S[0][_i][0] + ny[_i]*term;
      S[Cart::z][_i][0] = gma2*S[0][_i][0] + nz[_i]*term;
    }
  }
  //------------------------------------------------------

  //Integrals     p - d - s     p - f - s     p - g - s     p - h - s     p - i - s
  for (int _i =  4; _i < _ngw; _i++) {
    S[Cart::x][_i][0] = gma0*S[0][_i][0] + nx[_i]*fak*S[0][i_less_x[_i]][0];
    S[Cart::y][_i][0] = gma1*S[0][_i][0] + ny[_i]*fak*S[0][i_less_y[_i]][0];
    S[Cart::z][_i][0] = gma2*S[0][_i][0] + nz[_i]*fak*S[0][i_less_z[_i]][0];
  }
  //------------------------------------------------------

} // end if (_lmax_alpha > 0)


if (_lmax_alpha > 1) {

  //Integrals     d - s - s
  double term = fak*S[0][0][0];
  S[Cart::xx][0][0] = gma0*S[Cart::x][0][0] + term;
  S[Cart::xy][0][0] = gma0*S[Cart::y][0][0];
  S[Cart::xz][0][0] = gma0*S[Cart::z][0][0];
  S[Cart::yy][0][0] = gma1*S[Cart::y][0][0] + term;
  S[Cart::yz][0][0] = gma1*S[Cart::z][0][0];
  S[Cart::zz][0][0] = gma2*S[Cart::z][0][0] + term;
  //------------------------------------------------------

  //Integrals     d - p - s     d - d - s     d - f - s     d - g - s     d - h - s     d - i - s
  for (int _i =  1; _i < _ngw; _i++) {
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term = fak*S[0][_i][0];
    S[Cart::xx][_i][0] = gma0*S[Cart::x][_i][0] + nx_i*fak*S[Cart::x][ilx_i][0] + term;
    S[Cart::xy][_i][0] = gma0*S[Cart::y][_i][0] + nx_i*fak*S[Cart::y][ilx_i][0];
    S[Cart::xz][_i][0] = gma0*S[Cart::z][_i][0] + nx_i*fak*S[Cart::z][ilx_i][0];
    S[Cart::yy][_i][0] = gma1*S[Cart::y][_i][0] + ny_i*fak*S[Cart::y][ily_i][0] + term;
    S[Cart::yz][_i][0] = gma1*S[Cart::z][_i][0] + ny_i*fak*S[Cart::z][ily_i][0];
    S[Cart::zz][_i][0] = gma2*S[Cart::z][_i][0] + nz_i*fak*S[Cart::z][ilz_i][0] + term;
  }
  //------------------------------------------------------

} // end if (_lmax_alpha > 1)


if (_lmax_alpha > 2) {

  //Integrals     f - s - s
  S[Cart::xxx][0][0] = gma0*S[Cart::xx][0][0] + 2*fak*S[Cart::x][0][0];
  S[Cart::xxy][0][0] = gma1*S[Cart::xx][0][0];
  S[Cart::xxz][0][0] = gma2*S[Cart::xx][0][0];
  S[Cart::xyy][0][0] = gma0*S[Cart::yy][0][0];
  S[Cart::xyz][0][0] = gma0*S[Cart::yz][0][0];
  S[Cart::xzz][0][0] = gma0*S[Cart::zz][0][0];
  S[Cart::yyy][0][0] = gma1*S[Cart::yy][0][0] + 2*fak*S[Cart::y][0][0];
  S[Cart::yyz][0][0] = gma2*S[Cart::yy][0][0];
  S[Cart::yzz][0][0] = gma1*S[Cart::zz][0][0];
  S[Cart::zzz][0][0] = gma2*S[Cart::zz][0][0] + 2*fak*S[Cart::z][0][0];
  //------------------------------------------------------

  //Integrals     f - p - s     f - d - s     f - f - s     f - g - s     f - h - s     f - i - s
  for (int _i =  1; _i < _ngw; _i++) {
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_x = 2*fak*S[Cart::x][_i][0];
    double term_y = 2*fak*S[Cart::y][_i][0];
    double term_z = 2*fak*S[Cart::z][_i][0];
    S[Cart::xxx][_i][0] = gma0*S[Cart::xx][_i][0] + nx_i*fak*S[Cart::xx][ilx_i][0] + term_x;
    S[Cart::xxy][_i][0] = gma1*S[Cart::xx][_i][0] + ny_i*fak*S[Cart::xx][ily_i][0];
    S[Cart::xxz][_i][0] = gma2*S[Cart::xx][_i][0] + nz_i*fak*S[Cart::xx][ilz_i][0];
    S[Cart::xyy][_i][0] = gma0*S[Cart::yy][_i][0] + nx_i*fak*S[Cart::yy][ilx_i][0];
    S[Cart::xyz][_i][0] = gma0*S[Cart::yz][_i][0] + nx_i*fak*S[Cart::yz][ilx_i][0];
    S[Cart::xzz][_i][0] = gma0*S[Cart::zz][_i][0] + nx_i*fak*S[Cart::zz][ilx_i][0];
    S[Cart::yyy][_i][0] = gma1*S[Cart::yy][_i][0] + ny_i*fak*S[Cart::yy][ily_i][0] + term_y;
    S[Cart::yyz][_i][0] = gma2*S[Cart::yy][_i][0] + nz_i*fak*S[Cart::yy][ilz_i][0];
    S[Cart::yzz][_i][0] = gma1*S[Cart::zz][_i][0] + ny_i*fak*S[Cart::zz][ily_i][0];
    S[Cart::zzz][_i][0] = gma2*S[Cart::zz][_i][0] + nz_i*fak*S[Cart::zz][ilz_i][0] + term_z;
  }
  //------------------------------------------------------

} // end if (_lmax_alpha > 2)


if (_lmax_alpha > 3) {

  //Integrals     g - s - s
  double term_xx = fak*S[Cart::xx][0][0];
  double term_yy = fak*S[Cart::yy][0][0];
  double term_zz = fak*S[Cart::zz][0][0];
  S[Cart::xxxx][0][0] = gma0*S[Cart::xxx][0][0] + 3*term_xx;
  S[Cart::xxxy][0][0] = gma1*S[Cart::xxx][0][0];
  S[Cart::xxxz][0][0] = gma2*S[Cart::xxx][0][0];
  S[Cart::xxyy][0][0] = gma0*S[Cart::xyy][0][0] + term_yy;
  S[Cart::xxyz][0][0] = gma1*S[Cart::xxz][0][0];
  S[Cart::xxzz][0][0] = gma0*S[Cart::xzz][0][0] + term_zz;
  S[Cart::xyyy][0][0] = gma0*S[Cart::yyy][0][0];
  S[Cart::xyyz][0][0] = gma0*S[Cart::yyz][0][0];
  S[Cart::xyzz][0][0] = gma0*S[Cart::yzz][0][0];
  S[Cart::xzzz][0][0] = gma0*S[Cart::zzz][0][0];
  S[Cart::yyyy][0][0] = gma1*S[Cart::yyy][0][0] + 3*term_yy;
  S[Cart::yyyz][0][0] = gma2*S[Cart::yyy][0][0];
  S[Cart::yyzz][0][0] = gma1*S[Cart::yzz][0][0] + term_zz;
  S[Cart::yzzz][0][0] = gma1*S[Cart::zzz][0][0];
  S[Cart::zzzz][0][0] = gma2*S[Cart::zzz][0][0] + 3*term_zz;
  //------------------------------------------------------

  //Integrals     g - p - s     g - d - s     g - f - s     g - g - s     g - h - s     g - i - s
  for (int _i =  1; _i < _ngw; _i++) {
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_xx = fak*S[Cart::xx][_i][0];
    double term_yy = fak*S[Cart::yy][_i][0];
    double term_zz = fak*S[Cart::zz][_i][0];
    S[Cart::xxxx][_i][0] = gma0*S[Cart::xxx][_i][0] + nx_i*fak*S[Cart::xxx][ilx_i][0] + 3*term_xx;
    S[Cart::xxxy][_i][0] = gma1*S[Cart::xxx][_i][0] + ny_i*fak*S[Cart::xxx][ily_i][0];
    S[Cart::xxxz][_i][0] = gma2*S[Cart::xxx][_i][0] + nz_i*fak*S[Cart::xxx][ilz_i][0];
    S[Cart::xxyy][_i][0] = gma0*S[Cart::xyy][_i][0] + nx_i*fak*S[Cart::xyy][ilx_i][0] + term_yy;
    S[Cart::xxyz][_i][0] = gma1*S[Cart::xxz][_i][0] + ny_i*fak*S[Cart::xxz][ily_i][0];
    S[Cart::xxzz][_i][0] = gma0*S[Cart::xzz][_i][0] + nx_i*fak*S[Cart::xzz][ilx_i][0] + term_zz;
    S[Cart::xyyy][_i][0] = gma0*S[Cart::yyy][_i][0] + nx_i*fak*S[Cart::yyy][ilx_i][0];
    S[Cart::xyyz][_i][0] = gma0*S[Cart::yyz][_i][0] + nx_i*fak*S[Cart::yyz][ilx_i][0];
    S[Cart::xyzz][_i][0] = gma0*S[Cart::yzz][_i][0] + nx_i*fak*S[Cart::yzz][ilx_i][0];
    S[Cart::xzzz][_i][0] = gma0*S[Cart::zzz][_i][0] + nx_i*fak*S[Cart::zzz][ilx_i][0];
    S[Cart::yyyy][_i][0] = gma1*S[Cart::yyy][_i][0] + ny_i*fak*S[Cart::yyy][ily_i][0] + 3*term_yy;
    S[Cart::yyyz][_i][0] = gma2*S[Cart::yyy][_i][0] + nz_i*fak*S[Cart::yyy][ilz_i][0];
    S[Cart::yyzz][_i][0] = gma1*S[Cart::yzz][_i][0] + ny_i*fak*S[Cart::yzz][ily_i][0] + term_zz;
    S[Cart::yzzz][_i][0] = gma1*S[Cart::zzz][_i][0] + ny_i*fak*S[Cart::zzz][ily_i][0];
    S[Cart::zzzz][_i][0] = gma2*S[Cart::zzz][_i][0] + nz_i*fak*S[Cart::zzz][ilz_i][0] + 3*term_zz;
    }
  //------------------------------------------------------

} // end if (_lmax_alpha > 3)


if (_lmax_alpha > 4) {

  //Integrals     h - s - s
  double term_xxx = fak*S[Cart::xxx][0][0];
  double term_yyy = fak*S[Cart::yyy][0][0];
  double term_zzz = fak*S[Cart::zzz][0][0];
  S[Cart::xxxxx][0][0] = gma0*S[Cart::xxxx][0][0] + 4*term_xxx;
  S[Cart::xxxxy][0][0] = gma1*S[Cart::xxxx][0][0];
  S[Cart::xxxxz][0][0] = gma2*S[Cart::xxxx][0][0];
  S[Cart::xxxyy][0][0] = gma1*S[Cart::xxxy][0][0] + term_xxx;
  S[Cart::xxxyz][0][0] = gma1*S[Cart::xxxz][0][0];
  S[Cart::xxxzz][0][0] = gma2*S[Cart::xxxz][0][0] + term_xxx;
  S[Cart::xxyyy][0][0] = gma0*S[Cart::xyyy][0][0] + term_yyy;
  S[Cart::xxyyz][0][0] = gma2*S[Cart::xxyy][0][0];
  S[Cart::xxyzz][0][0] = gma1*S[Cart::xxzz][0][0];
  S[Cart::xxzzz][0][0] = gma0*S[Cart::xzzz][0][0] + term_zzz;
  S[Cart::xyyyy][0][0] = gma0*S[Cart::yyyy][0][0];
  S[Cart::xyyyz][0][0] = gma0*S[Cart::yyyz][0][0];
  S[Cart::xyyzz][0][0] = gma0*S[Cart::yyzz][0][0];
  S[Cart::xyzzz][0][0] = gma0*S[Cart::yzzz][0][0];
  S[Cart::xzzzz][0][0] = gma0*S[Cart::zzzz][0][0];
  S[Cart::yyyyy][0][0] = gma1*S[Cart::yyyy][0][0] + 4*term_yyy;
  S[Cart::yyyyz][0][0] = gma2*S[Cart::yyyy][0][0];
  S[Cart::yyyzz][0][0] = gma2*S[Cart::yyyz][0][0] + term_yyy;
  S[Cart::yyzzz][0][0] = gma1*S[Cart::yzzz][0][0] + term_zzz;
  S[Cart::yzzzz][0][0] = gma1*S[Cart::zzzz][0][0];
  S[Cart::zzzzz][0][0] = gma2*S[Cart::zzzz][0][0] + 4*term_zzz;
  //------------------------------------------------------

  //Integrals     h - p - s     h - d - s     h - f - s     h - g - s     h - h - s     h - i - s
  for (int _i =  1; _i < _ngw; _i++) {
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_xxx = fak*S[Cart::xxx][_i][0];
    double term_yyy = fak*S[Cart::yyy][_i][0];
    double term_zzz = fak*S[Cart::zzz][_i][0];
    S[Cart::xxxxx][_i][0] = gma0*S[Cart::xxxx][_i][0] + nx_i*fak*S[Cart::xxxx][ilx_i][0] + 4*term_xxx;
    S[Cart::xxxxy][_i][0] = gma1*S[Cart::xxxx][_i][0] + ny_i*fak*S[Cart::xxxx][ily_i][0];
    S[Cart::xxxxz][_i][0] = gma2*S[Cart::xxxx][_i][0] + nz_i*fak*S[Cart::xxxx][ilz_i][0];
    S[Cart::xxxyy][_i][0] = gma1*S[Cart::xxxy][_i][0] + ny_i*fak*S[Cart::xxxy][ily_i][0] + term_xxx;
    S[Cart::xxxyz][_i][0] = gma1*S[Cart::xxxz][_i][0] + ny_i*fak*S[Cart::xxxz][ily_i][0];
    S[Cart::xxxzz][_i][0] = gma2*S[Cart::xxxz][_i][0] + nz_i*fak*S[Cart::xxxz][ilz_i][0] + term_xxx;
    S[Cart::xxyyy][_i][0] = gma0*S[Cart::xyyy][_i][0] + nx_i*fak*S[Cart::xyyy][ilx_i][0] + term_yyy;
    S[Cart::xxyyz][_i][0] = gma2*S[Cart::xxyy][_i][0] + nz_i*fak*S[Cart::xxyy][ilz_i][0];
    S[Cart::xxyzz][_i][0] = gma1*S[Cart::xxzz][_i][0] + ny_i*fak*S[Cart::xxzz][ily_i][0];
    S[Cart::xxzzz][_i][0] = gma0*S[Cart::xzzz][_i][0] + nx_i*fak*S[Cart::xzzz][ilx_i][0] + term_zzz;
    S[Cart::xyyyy][_i][0] = gma0*S[Cart::yyyy][_i][0] + nx_i*fak*S[Cart::yyyy][ilx_i][0];
    S[Cart::xyyyz][_i][0] = gma0*S[Cart::yyyz][_i][0] + nx_i*fak*S[Cart::yyyz][ilx_i][0];
    S[Cart::xyyzz][_i][0] = gma0*S[Cart::yyzz][_i][0] + nx_i*fak*S[Cart::yyzz][ilx_i][0];
    S[Cart::xyzzz][_i][0] = gma0*S[Cart::yzzz][_i][0] + nx_i*fak*S[Cart::yzzz][ilx_i][0];
    S[Cart::xzzzz][_i][0] = gma0*S[Cart::zzzz][_i][0] + nx_i*fak*S[Cart::zzzz][ilx_i][0];
    S[Cart::yyyyy][_i][0] = gma1*S[Cart::yyyy][_i][0] + ny_i*fak*S[Cart::yyyy][ily_i][0] + 4*term_yyy;
    S[Cart::yyyyz][_i][0] = gma2*S[Cart::yyyy][_i][0] + nz_i*fak*S[Cart::yyyy][ilz_i][0];
    S[Cart::yyyzz][_i][0] = gma2*S[Cart::yyyz][_i][0] + nz_i*fak*S[Cart::yyyz][ilz_i][0] + term_yyy;
    S[Cart::yyzzz][_i][0] = gma1*S[Cart::yzzz][_i][0] + ny_i*fak*S[Cart::yzzz][ily_i][0] + term_zzz;
    S[Cart::yzzzz][_i][0] = gma1*S[Cart::zzzz][_i][0] + ny_i*fak*S[Cart::zzzz][ily_i][0];
    S[Cart::zzzzz][_i][0] = gma2*S[Cart::zzzz][_i][0] + nz_i*fak*S[Cart::zzzz][ilz_i][0] + 4*term_zzz;
  }
  //------------------------------------------------------

} // end if (_lmax_alpha > 4)


if (_lmax_alpha > 5) {

  //Integrals     i - s - s
  double term_xxxx = fak*S[Cart::xxxx][0][0];
  double term_xyyy = fak*S[Cart::xyyy][0][0];
  double term_xzzz = fak*S[Cart::xzzz][0][0];
  double term_yyyy = fak*S[Cart::yyyy][0][0];
  double term_yyzz = fak*S[Cart::yyzz][0][0];
  double term_yzzz = fak*S[Cart::yzzz][0][0];
  double term_zzzz = fak*S[Cart::zzzz][0][0];
  S[Cart::xxxxxx][0][0] = gma0*S[Cart::xxxxx][0][0] + 5*term_xxxx;
  S[Cart::xxxxxy][0][0] = gma1*S[Cart::xxxxx][0][0];
  S[Cart::xxxxxz][0][0] = gma2*S[Cart::xxxxx][0][0];
  S[Cart::xxxxyy][0][0] = gma1*S[Cart::xxxxy][0][0] + term_xxxx;
  S[Cart::xxxxyz][0][0] = gma1*S[Cart::xxxxz][0][0];
  S[Cart::xxxxzz][0][0] = gma2*S[Cart::xxxxz][0][0] + term_xxxx;
  S[Cart::xxxyyy][0][0] = gma0*S[Cart::xxyyy][0][0] + 2*term_xyyy;
  S[Cart::xxxyyz][0][0] = gma2*S[Cart::xxxyy][0][0];
  S[Cart::xxxyzz][0][0] = gma1*S[Cart::xxxzz][0][0];
  S[Cart::xxxzzz][0][0] = gma0*S[Cart::xxzzz][0][0] + 2*term_xzzz;
  S[Cart::xxyyyy][0][0] = gma0*S[Cart::xyyyy][0][0] + term_yyyy;
  S[Cart::xxyyyz][0][0] = gma2*S[Cart::xxyyy][0][0];
  S[Cart::xxyyzz][0][0] = gma0*S[Cart::xyyzz][0][0] + term_yyzz;
  S[Cart::xxyzzz][0][0] = gma1*S[Cart::xxzzz][0][0];
  S[Cart::xxzzzz][0][0] = gma0*S[Cart::xzzzz][0][0] + term_zzzz;
  S[Cart::xyyyyy][0][0] = gma0*S[Cart::yyyyy][0][0];
  S[Cart::xyyyyz][0][0] = gma0*S[Cart::yyyyz][0][0];
  S[Cart::xyyyzz][0][0] = gma0*S[Cart::yyyzz][0][0];
  S[Cart::xyyzzz][0][0] = gma0*S[Cart::yyzzz][0][0];
  S[Cart::xyzzzz][0][0] = gma0*S[Cart::yzzzz][0][0];
  S[Cart::xzzzzz][0][0] = gma0*S[Cart::zzzzz][0][0];
  S[Cart::yyyyyy][0][0] = gma1*S[Cart::yyyyy][0][0] + 5*term_yyyy;
  S[Cart::yyyyyz][0][0] = gma2*S[Cart::yyyyy][0][0];
  S[Cart::yyyyzz][0][0] = gma2*S[Cart::yyyyz][0][0] + term_yyyy;
  S[Cart::yyyzzz][0][0] = gma1*S[Cart::yyzzz][0][0] + 2*term_yzzz;
  S[Cart::yyzzzz][0][0] = gma1*S[Cart::yzzzz][0][0] + term_zzzz;
  S[Cart::yzzzzz][0][0] = gma1*S[Cart::zzzzz][0][0];
  S[Cart::zzzzzz][0][0] = gma2*S[Cart::zzzzz][0][0] + 5*term_zzzz;
  //------------------------------------------------------

  //Integrals     i - p - s     i - d - s     i - f - s     i - g - s     i - h - s     i - i - s
  for (int _i =  1; _i < _ngw; _i++) {
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_xxxx = fak*S[Cart::xxxx][_i][0];
    double term_xyyy = fak*S[Cart::xyyy][_i][0];
    double term_xzzz = fak*S[Cart::xzzz][_i][0];
    double term_yyyy = fak*S[Cart::yyyy][_i][0];
    double term_yyzz = fak*S[Cart::yyzz][_i][0];
    double term_yzzz = fak*S[Cart::yzzz][_i][0];
    double term_zzzz = fak*S[Cart::zzzz][_i][0];
    S[Cart::xxxxxx][_i][0] = gma0*S[Cart::xxxxx][_i][0] + nx_i*fak*S[Cart::xxxxx][ilx_i][0] + 5*term_xxxx;
    S[Cart::xxxxxy][_i][0] = gma1*S[Cart::xxxxx][_i][0] + ny_i*fak*S[Cart::xxxxx][ily_i][0];
    S[Cart::xxxxxz][_i][0] = gma2*S[Cart::xxxxx][_i][0] + nz_i*fak*S[Cart::xxxxx][ilz_i][0];
    S[Cart::xxxxyy][_i][0] = gma1*S[Cart::xxxxy][_i][0] + ny_i*fak*S[Cart::xxxxy][ily_i][0] + term_xxxx;
    S[Cart::xxxxyz][_i][0] = gma1*S[Cart::xxxxz][_i][0] + ny_i*fak*S[Cart::xxxxz][ily_i][0];
    S[Cart::xxxxzz][_i][0] = gma2*S[Cart::xxxxz][_i][0] + nz_i*fak*S[Cart::xxxxz][ilz_i][0] + term_xxxx;
    S[Cart::xxxyyy][_i][0] = gma0*S[Cart::xxyyy][_i][0] + nx_i*fak*S[Cart::xxyyy][ilx_i][0] + 2*term_xyyy;
    S[Cart::xxxyyz][_i][0] = gma2*S[Cart::xxxyy][_i][0] + nz_i*fak*S[Cart::xxxyy][ilz_i][0];
    S[Cart::xxxyzz][_i][0] = gma1*S[Cart::xxxzz][_i][0] + ny_i*fak*S[Cart::xxxzz][ily_i][0];
    S[Cart::xxxzzz][_i][0] = gma0*S[Cart::xxzzz][_i][0] + nx_i*fak*S[Cart::xxzzz][ilx_i][0] + 2*term_xzzz;
    S[Cart::xxyyyy][_i][0] = gma0*S[Cart::xyyyy][_i][0] + nx_i*fak*S[Cart::xyyyy][ilx_i][0] + term_yyyy;
    S[Cart::xxyyyz][_i][0] = gma2*S[Cart::xxyyy][_i][0] + nz_i*fak*S[Cart::xxyyy][ilz_i][0];
    S[Cart::xxyyzz][_i][0] = gma0*S[Cart::xyyzz][_i][0] + nx_i*fak*S[Cart::xyyzz][ilx_i][0] + term_yyzz;
    S[Cart::xxyzzz][_i][0] = gma1*S[Cart::xxzzz][_i][0] + ny_i*fak*S[Cart::xxzzz][ily_i][0];
    S[Cart::xxzzzz][_i][0] = gma0*S[Cart::xzzzz][_i][0] + nx_i*fak*S[Cart::xzzzz][ilx_i][0] + term_zzzz;
    S[Cart::xyyyyy][_i][0] = gma0*S[Cart::yyyyy][_i][0] + nx_i*fak*S[Cart::yyyyy][ilx_i][0];
    S[Cart::xyyyyz][_i][0] = gma0*S[Cart::yyyyz][_i][0] + nx_i*fak*S[Cart::yyyyz][ilx_i][0];
    S[Cart::xyyyzz][_i][0] = gma0*S[Cart::yyyzz][_i][0] + nx_i*fak*S[Cart::yyyzz][ilx_i][0];
    S[Cart::xyyzzz][_i][0] = gma0*S[Cart::yyzzz][_i][0] + nx_i*fak*S[Cart::yyzzz][ilx_i][0];
    S[Cart::xyzzzz][_i][0] = gma0*S[Cart::yzzzz][_i][0] + nx_i*fak*S[Cart::yzzzz][ilx_i][0];
    S[Cart::xzzzzz][_i][0] = gma0*S[Cart::zzzzz][_i][0] + nx_i*fak*S[Cart::zzzzz][ilx_i][0];
    S[Cart::yyyyyy][_i][0] = gma1*S[Cart::yyyyy][_i][0] + ny_i*fak*S[Cart::yyyyy][ily_i][0] + 5*term_yyyy;
    S[Cart::yyyyyz][_i][0] = gma2*S[Cart::yyyyy][_i][0] + nz_i*fak*S[Cart::yyyyy][ilz_i][0];
    S[Cart::yyyyzz][_i][0] = gma2*S[Cart::yyyyz][_i][0] + nz_i*fak*S[Cart::yyyyz][ilz_i][0] + term_yyyy;
    S[Cart::yyyzzz][_i][0] = gma1*S[Cart::yyzzz][_i][0] + ny_i*fak*S[Cart::yyzzz][ily_i][0] + 2*term_yzzz;
    S[Cart::yyzzzz][_i][0] = gma1*S[Cart::yzzzz][_i][0] + ny_i*fak*S[Cart::yzzzz][ily_i][0] + term_zzzz;
    S[Cart::yzzzzz][_i][0] = gma1*S[Cart::zzzzz][_i][0] + ny_i*fak*S[Cart::zzzzz][ily_i][0];
    S[Cart::zzzzzz][_i][0] = gma2*S[Cart::zzzzz][_i][0] + nz_i*fak*S[Cart::zzzzz][ilz_i][0] + 5*term_zzzz;
  }
  //------------------------------------------------------

} // end if (_lmax_alpha > 5)







if (_lmax_gamma > 0) {

  //Integrals     * - * - p
  for (int _j =  0; _j < _nalpha; _j++) {
    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    for (int _i =  0; _i < _ngw; _i++) {
      S[_j][_i][Cart::x] = gmb0*S[_j][_i][0] + fak*(nx[_i]*S[_j][i_less_x[_i]][0]+nx_j*S[ilx_j][_i][0]);
      S[_j][_i][Cart::y] = gmb1*S[_j][_i][0] + fak*(ny[_i]*S[_j][i_less_y[_i]][0]+ny_j*S[ily_j][_i][0]);
      S[_j][_i][Cart::z] = gmb2*S[_j][_i][0] + fak*(nz[_i]*S[_j][i_less_z[_i]][0]+nz_j*S[ilz_j][_i][0]);
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 0)


if (_lmax_gamma > 1) {

  //Integrals     * - * - d
  for (int _j =  0; _j < _nalpha; _j++) {
    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    for (int _i =  0; _i < _ngw; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      double term = fak*S[_j][_i][0];
      S[_j][_i][Cart::xx] = gmb0*S[_j][_i][Cart::x] + fak*(nx_i*S[_j][ilx_i][Cart::x]+nx_j*S[ilx_j][_i][Cart::x]) + term;
      S[_j][_i][Cart::xy] = gmb0*S[_j][_i][Cart::y] + fak*(nx_i*S[_j][ilx_i][Cart::y]+nx_j*S[ilx_j][_i][Cart::y]);
      S[_j][_i][Cart::xz] = gmb0*S[_j][_i][Cart::z] + fak*(nx_i*S[_j][ilx_i][Cart::z]+nx_j*S[ilx_j][_i][Cart::z]);
      S[_j][_i][Cart::yy] = gmb1*S[_j][_i][Cart::y] + fak*(ny_i*S[_j][ily_i][Cart::y]+ny_j*S[ily_j][_i][Cart::y]) + term;
      S[_j][_i][Cart::yz] = gmb1*S[_j][_i][Cart::z] + fak*(ny_i*S[_j][ily_i][Cart::z]+ny_j*S[ily_j][_i][Cart::z]);
      S[_j][_i][Cart::zz] = gmb2*S[_j][_i][Cart::z] + fak*(nz_i*S[_j][ilz_i][Cart::z]+nz_j*S[ilz_j][_i][Cart::z]) + term;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 1)


if (_lmax_gamma > 2) {

  //Integrals     * - * - f
  for (int _j =  0; _j < _nalpha; _j++) {
    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    for (int _i =  0; _i < _ngw; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      double term_x = 2*fak*S[_j][_i][Cart::x];
      double term_y = 2*fak*S[_j][_i][Cart::y];
      double term_z = 2*fak*S[_j][_i][Cart::z];
      S[_j][_i][Cart::xxx] = gmb0*S[_j][_i][Cart::xx] + fak*(nx_i*S[_j][ilx_i][Cart::xx]+nx_j*S[ilx_j][_i][Cart::xx]) + term_x;
      S[_j][_i][Cart::xxy] = gmb1*S[_j][_i][Cart::xx] + fak*(ny_i*S[_j][ily_i][Cart::xx]+ny_j*S[ily_j][_i][Cart::xx]);
      S[_j][_i][Cart::xxz] = gmb2*S[_j][_i][Cart::xx] + fak*(nz_i*S[_j][ilz_i][Cart::xx]+nz_j*S[ilz_j][_i][Cart::xx]);
      S[_j][_i][Cart::xyy] = gmb0*S[_j][_i][Cart::yy] + fak*(nx_i*S[_j][ilx_i][Cart::yy]+nx_j*S[ilx_j][_i][Cart::yy]);
      S[_j][_i][Cart::xyz] = gmb0*S[_j][_i][Cart::yz] + fak*(nx_i*S[_j][ilx_i][Cart::yz]+nx_j*S[ilx_j][_i][Cart::yz]);
      S[_j][_i][Cart::xzz] = gmb0*S[_j][_i][Cart::zz] + fak*(nx_i*S[_j][ilx_i][Cart::zz]+nx_j*S[ilx_j][_i][Cart::zz]);
      S[_j][_i][Cart::yyy] = gmb1*S[_j][_i][Cart::yy] + fak*(ny_i*S[_j][ily_i][Cart::yy]+ny_j*S[ily_j][_i][Cart::yy]) + term_y;
      S[_j][_i][Cart::yyz] = gmb2*S[_j][_i][Cart::yy] + fak*(nz_i*S[_j][ilz_i][Cart::yy]+nz_j*S[ilz_j][_i][Cart::yy]);
      S[_j][_i][Cart::yzz] = gmb1*S[_j][_i][Cart::zz] + fak*(ny_i*S[_j][ily_i][Cart::zz]+ny_j*S[ily_j][_i][Cart::zz]);
      S[_j][_i][Cart::zzz] = gmb2*S[_j][_i][Cart::zz] + fak*(nz_i*S[_j][ilz_i][Cart::zz]+nz_j*S[ilz_j][_i][Cart::zz]) + term_z;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 2)


if (_lmax_gamma > 3) {

  //Integrals     * - * - g
  for (int _j =  0; _j < _nalpha; _j++) {
    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    for (int _i =  0; _i < _ngw; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      double term_xx = fak*S[_j][_i][Cart::xx];
      double term_yy = fak*S[_j][_i][Cart::yy];
      double term_zz = fak*S[_j][_i][Cart::zz];
      S[_j][_i][Cart::xxxx] = gmb0*S[_j][_i][Cart::xxx] + fak*(nx_i*S[_j][ilx_i][Cart::xxx]+nx_j*S[ilx_j][_i][Cart::xxx]) + 3*term_xx;
      S[_j][_i][Cart::xxxy] = gmb1*S[_j][_i][Cart::xxx] + fak*(ny_i*S[_j][ily_i][Cart::xxx]+ny_j*S[ily_j][_i][Cart::xxx]);
      S[_j][_i][Cart::xxxz] = gmb2*S[_j][_i][Cart::xxx] + fak*(nz_i*S[_j][ilz_i][Cart::xxx]+nz_j*S[ilz_j][_i][Cart::xxx]);
      S[_j][_i][Cart::xxyy] = gmb0*S[_j][_i][Cart::xyy] + fak*(nx_i*S[_j][ilx_i][Cart::xyy]+nx_j*S[ilx_j][_i][Cart::xyy]) + term_yy;
      S[_j][_i][Cart::xxyz] = gmb1*S[_j][_i][Cart::xxz] + fak*(ny_i*S[_j][ily_i][Cart::xxz]+ny_j*S[ily_j][_i][Cart::xxz]);
      S[_j][_i][Cart::xxzz] = gmb0*S[_j][_i][Cart::xzz] + fak*(nx_i*S[_j][ilx_i][Cart::xzz]+nx_j*S[ilx_j][_i][Cart::xzz]) + term_zz;
      S[_j][_i][Cart::xyyy] = gmb0*S[_j][_i][Cart::yyy] + fak*(nx_i*S[_j][ilx_i][Cart::yyy]+nx_j*S[ilx_j][_i][Cart::yyy]);
      S[_j][_i][Cart::xyyz] = gmb0*S[_j][_i][Cart::yyz] + fak*(nx_i*S[_j][ilx_i][Cart::yyz]+nx_j*S[ilx_j][_i][Cart::yyz]);
      S[_j][_i][Cart::xyzz] = gmb0*S[_j][_i][Cart::yzz] + fak*(nx_i*S[_j][ilx_i][Cart::yzz]+nx_j*S[ilx_j][_i][Cart::yzz]);
      S[_j][_i][Cart::xzzz] = gmb0*S[_j][_i][Cart::zzz] + fak*(nx_i*S[_j][ilx_i][Cart::zzz]+nx_j*S[ilx_j][_i][Cart::zzz]);
      S[_j][_i][Cart::yyyy] = gmb1*S[_j][_i][Cart::yyy] + fak*(ny_i*S[_j][ily_i][Cart::yyy]+ny_j*S[ily_j][_i][Cart::yyy]) + 3*term_yy;
      S[_j][_i][Cart::yyyz] = gmb2*S[_j][_i][Cart::yyy] + fak*(nz_i*S[_j][ilz_i][Cart::yyy]+nz_j*S[ilz_j][_i][Cart::yyy]);
      S[_j][_i][Cart::yyzz] = gmb1*S[_j][_i][Cart::yzz] + fak*(ny_i*S[_j][ily_i][Cart::yzz]+ny_j*S[ily_j][_i][Cart::yzz]) + term_zz;
      S[_j][_i][Cart::yzzz] = gmb1*S[_j][_i][Cart::zzz] + fak*(ny_i*S[_j][ily_i][Cart::zzz]+ny_j*S[ily_j][_i][Cart::zzz]);
      S[_j][_i][Cart::zzzz] = gmb2*S[_j][_i][Cart::zzz] + fak*(nz_i*S[_j][ilz_i][Cart::zzz]+nz_j*S[ilz_j][_i][Cart::zzz]) + 3*term_zz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 3)


if (_lmax_gamma > 4) {

  //Integrals     * - * - h
  for (int _j =  0; _j < _nalpha; _j++) {
    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    for (int _i =  0; _i < _ngw; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      double term_xxx = fak*S[_j][_i][Cart::xxx];
      double term_yyy = fak*S[_j][_i][Cart::yyy];
      double term_zzz = fak*S[_j][_i][Cart::zzz];
      S[_j][_i][Cart::xxxxx] = gmb0*S[_j][_i][Cart::xxxx] + fak*(nx_i*S[_j][ilx_i][Cart::xxxx]+nx_j*S[ilx_j][_i][Cart::xxxx]) + 4*term_xxx;
      S[_j][_i][Cart::xxxxy] = gmb1*S[_j][_i][Cart::xxxx] + fak*(ny_i*S[_j][ily_i][Cart::xxxx]+ny_j*S[ily_j][_i][Cart::xxxx]);
      S[_j][_i][Cart::xxxxz] = gmb2*S[_j][_i][Cart::xxxx] + fak*(nz_i*S[_j][ilz_i][Cart::xxxx]+nz_j*S[ilz_j][_i][Cart::xxxx]);
      S[_j][_i][Cart::xxxyy] = gmb1*S[_j][_i][Cart::xxxy] + fak*(ny_i*S[_j][ily_i][Cart::xxxy]+ny_j*S[ily_j][_i][Cart::xxxy]) + term_xxx;
      S[_j][_i][Cart::xxxyz] = gmb1*S[_j][_i][Cart::xxxz] + fak*(ny_i*S[_j][ily_i][Cart::xxxz]+ny_j*S[ily_j][_i][Cart::xxxz]);
      S[_j][_i][Cart::xxxzz] = gmb2*S[_j][_i][Cart::xxxz] + fak*(nz_i*S[_j][ilz_i][Cart::xxxz]+nz_j*S[ilz_j][_i][Cart::xxxz]) + term_xxx;
      S[_j][_i][Cart::xxyyy] = gmb0*S[_j][_i][Cart::xyyy] + fak*(nx_i*S[_j][ilx_i][Cart::xyyy]+nx_j*S[ilx_j][_i][Cart::xyyy]) + term_yyy;
      S[_j][_i][Cart::xxyyz] = gmb2*S[_j][_i][Cart::xxyy] + fak*(nz_i*S[_j][ilz_i][Cart::xxyy]+nz_j*S[ilz_j][_i][Cart::xxyy]);
      S[_j][_i][Cart::xxyzz] = gmb1*S[_j][_i][Cart::xxzz] + fak*(ny_i*S[_j][ily_i][Cart::xxzz]+ny_j*S[ily_j][_i][Cart::xxzz]);
      S[_j][_i][Cart::xxzzz] = gmb0*S[_j][_i][Cart::xzzz] + fak*(nx_i*S[_j][ilx_i][Cart::xzzz]+nx_j*S[ilx_j][_i][Cart::xzzz]) + term_zzz;
      S[_j][_i][Cart::xyyyy] = gmb0*S[_j][_i][Cart::yyyy] + fak*(nx_i*S[_j][ilx_i][Cart::yyyy]+nx_j*S[ilx_j][_i][Cart::yyyy]);
      S[_j][_i][Cart::xyyyz] = gmb0*S[_j][_i][Cart::yyyz] + fak*(nx_i*S[_j][ilx_i][Cart::yyyz]+nx_j*S[ilx_j][_i][Cart::yyyz]);
      S[_j][_i][Cart::xyyzz] = gmb0*S[_j][_i][Cart::yyzz] + fak*(nx_i*S[_j][ilx_i][Cart::yyzz]+nx_j*S[ilx_j][_i][Cart::yyzz]);
      S[_j][_i][Cart::xyzzz] = gmb0*S[_j][_i][Cart::yzzz] + fak*(nx_i*S[_j][ilx_i][Cart::yzzz]+nx_j*S[ilx_j][_i][Cart::yzzz]);
      S[_j][_i][Cart::xzzzz] = gmb0*S[_j][_i][Cart::zzzz] + fak*(nx_i*S[_j][ilx_i][Cart::zzzz]+nx_j*S[ilx_j][_i][Cart::zzzz]);
      S[_j][_i][Cart::yyyyy] = gmb1*S[_j][_i][Cart::yyyy] + fak*(ny_i*S[_j][ily_i][Cart::yyyy]+ny_j*S[ily_j][_i][Cart::yyyy]) + 4*term_yyy;
      S[_j][_i][Cart::yyyyz] = gmb2*S[_j][_i][Cart::yyyy] + fak*(nz_i*S[_j][ilz_i][Cart::yyyy]+nz_j*S[ilz_j][_i][Cart::yyyy]);
      S[_j][_i][Cart::yyyzz] = gmb2*S[_j][_i][Cart::yyyz] + fak*(nz_i*S[_j][ilz_i][Cart::yyyz]+nz_j*S[ilz_j][_i][Cart::yyyz]) + term_yyy;
      S[_j][_i][Cart::yyzzz] = gmb1*S[_j][_i][Cart::yzzz] + fak*(ny_i*S[_j][ily_i][Cart::yzzz]+ny_j*S[ily_j][_i][Cart::yzzz]) + term_zzz;
      S[_j][_i][Cart::yzzzz] = gmb1*S[_j][_i][Cart::zzzz] + fak*(ny_i*S[_j][ily_i][Cart::zzzz]+ny_j*S[ily_j][_i][Cart::zzzz]);
      S[_j][_i][Cart::zzzzz] = gmb2*S[_j][_i][Cart::zzzz] + fak*(nz_i*S[_j][ilz_i][Cart::zzzz]+nz_j*S[ilz_j][_i][Cart::zzzz]) + 4*term_zzz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 4)


if (_lmax_gamma > 5) {

  //Integrals     * - * - i
  for (int _j =  0; _j < _nalpha; _j++) {
    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    for (int _i =  0; _i < _ngw; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      double term_xxxx = fak*S[_j][_i][Cart::xxxx];
      double term_xyyy = fak*S[_j][_i][Cart::xyyy];
      double term_xzzz = fak*S[_j][_i][Cart::xzzz];
      double term_yyyy = fak*S[_j][_i][Cart::yyyy];
      double term_yyzz = fak*S[_j][_i][Cart::yyzz];
      double term_yzzz = fak*S[_j][_i][Cart::yzzz];
      double term_zzzz = fak*S[_j][_i][Cart::zzzz];
      S[_j][_i][Cart::xxxxxx] = gmb0*S[_j][_i][Cart::xxxxx] + fak*(nx_i*S[_j][ilx_i][Cart::xxxxx]+nx_j*S[ilx_j][_i][Cart::xxxxx]) + 5*term_xxxx;
      S[_j][_i][Cart::xxxxxy] = gmb1*S[_j][_i][Cart::xxxxx] + fak*(ny_i*S[_j][ily_i][Cart::xxxxx]+ny_j*S[ily_j][_i][Cart::xxxxx]);
      S[_j][_i][Cart::xxxxxz] = gmb2*S[_j][_i][Cart::xxxxx] + fak*(nz_i*S[_j][ilz_i][Cart::xxxxx]+nz_j*S[ilz_j][_i][Cart::xxxxx]);
      S[_j][_i][Cart::xxxxyy] = gmb1*S[_j][_i][Cart::xxxxy] + fak*(ny_i*S[_j][ily_i][Cart::xxxxy]+ny_j*S[ily_j][_i][Cart::xxxxy]) + term_xxxx;
      S[_j][_i][Cart::xxxxyz] = gmb1*S[_j][_i][Cart::xxxxz] + fak*(ny_i*S[_j][ily_i][Cart::xxxxz]+ny_j*S[ily_j][_i][Cart::xxxxz]);
      S[_j][_i][Cart::xxxxzz] = gmb2*S[_j][_i][Cart::xxxxz] + fak*(nz_i*S[_j][ilz_i][Cart::xxxxz]+nz_j*S[ilz_j][_i][Cart::xxxxz]) + term_xxxx;
      S[_j][_i][Cart::xxxyyy] = gmb0*S[_j][_i][Cart::xxyyy] + fak*(nx_i*S[_j][ilx_i][Cart::xxyyy]+nx_j*S[ilx_j][_i][Cart::xxyyy]) + 2*term_xyyy;
      S[_j][_i][Cart::xxxyyz] = gmb2*S[_j][_i][Cart::xxxyy] + fak*(nz_i*S[_j][ilz_i][Cart::xxxyy]+nz_j*S[ilz_j][_i][Cart::xxxyy]);
      S[_j][_i][Cart::xxxyzz] = gmb1*S[_j][_i][Cart::xxxzz] + fak*(ny_i*S[_j][ily_i][Cart::xxxzz]+ny_j*S[ily_j][_i][Cart::xxxzz]);
      S[_j][_i][Cart::xxxzzz] = gmb0*S[_j][_i][Cart::xxzzz] + fak*(nx_i*S[_j][ilx_i][Cart::xxzzz]+nx_j*S[ilx_j][_i][Cart::xxzzz]) + 2*term_xzzz;
      S[_j][_i][Cart::xxyyyy] = gmb0*S[_j][_i][Cart::xyyyy] + fak*(nx_i*S[_j][ilx_i][Cart::xyyyy]+nx_j*S[ilx_j][_i][Cart::xyyyy]) + term_yyyy;
      S[_j][_i][Cart::xxyyyz] = gmb2*S[_j][_i][Cart::xxyyy] + fak*(nz_i*S[_j][ilz_i][Cart::xxyyy]+nz_j*S[ilz_j][_i][Cart::xxyyy]);
      S[_j][_i][Cart::xxyyzz] = gmb0*S[_j][_i][Cart::xyyzz] + fak*(nx_i*S[_j][ilx_i][Cart::xyyzz]+nx_j*S[ilx_j][_i][Cart::xyyzz]) + term_yyzz;
      S[_j][_i][Cart::xxyzzz] = gmb1*S[_j][_i][Cart::xxzzz] + fak*(ny_i*S[_j][ily_i][Cart::xxzzz]+ny_j*S[ily_j][_i][Cart::xxzzz]);
      S[_j][_i][Cart::xxzzzz] = gmb0*S[_j][_i][Cart::xzzzz] + fak*(nx_i*S[_j][ilx_i][Cart::xzzzz]+nx_j*S[ilx_j][_i][Cart::xzzzz]) + term_zzzz;
      S[_j][_i][Cart::xyyyyy] = gmb0*S[_j][_i][Cart::yyyyy] + fak*(nx_i*S[_j][ilx_i][Cart::yyyyy]+nx_j*S[ilx_j][_i][Cart::yyyyy]);
      S[_j][_i][Cart::xyyyyz] = gmb0*S[_j][_i][Cart::yyyyz] + fak*(nx_i*S[_j][ilx_i][Cart::yyyyz]+nx_j*S[ilx_j][_i][Cart::yyyyz]);
      S[_j][_i][Cart::xyyyzz] = gmb0*S[_j][_i][Cart::yyyzz] + fak*(nx_i*S[_j][ilx_i][Cart::yyyzz]+nx_j*S[ilx_j][_i][Cart::yyyzz]);
      S[_j][_i][Cart::xyyzzz] = gmb0*S[_j][_i][Cart::yyzzz] + fak*(nx_i*S[_j][ilx_i][Cart::yyzzz]+nx_j*S[ilx_j][_i][Cart::yyzzz]);
      S[_j][_i][Cart::xyzzzz] = gmb0*S[_j][_i][Cart::yzzzz] + fak*(nx_i*S[_j][ilx_i][Cart::yzzzz]+nx_j*S[ilx_j][_i][Cart::yzzzz]);
      S[_j][_i][Cart::xzzzzz] = gmb0*S[_j][_i][Cart::zzzzz] + fak*(nx_i*S[_j][ilx_i][Cart::zzzzz]+nx_j*S[ilx_j][_i][Cart::zzzzz]);
      S[_j][_i][Cart::yyyyyy] = gmb1*S[_j][_i][Cart::yyyyy] + fak*(ny_i*S[_j][ily_i][Cart::yyyyy]+ny_j*S[ily_j][_i][Cart::yyyyy]) + 5*term_yyyy;
      S[_j][_i][Cart::yyyyyz] = gmb2*S[_j][_i][Cart::yyyyy] + fak*(nz_i*S[_j][ilz_i][Cart::yyyyy]+nz_j*S[ilz_j][_i][Cart::yyyyy]);
      S[_j][_i][Cart::yyyyzz] = gmb2*S[_j][_i][Cart::yyyyz] + fak*(nz_i*S[_j][ilz_i][Cart::yyyyz]+nz_j*S[ilz_j][_i][Cart::yyyyz]) + term_yyyy;
      S[_j][_i][Cart::yyyzzz] = gmb1*S[_j][_i][Cart::yyzzz] + fak*(ny_i*S[_j][ily_i][Cart::yyzzz]+ny_j*S[ily_j][_i][Cart::yyzzz]) + 2*term_yzzz;
      S[_j][_i][Cart::yyzzzz] = gmb1*S[_j][_i][Cart::yzzzz] + fak*(ny_i*S[_j][ily_i][Cart::yzzzz]+ny_j*S[ily_j][_i][Cart::yzzzz]) + term_zzzz;
      S[_j][_i][Cart::yzzzzz] = gmb1*S[_j][_i][Cart::zzzzz] + fak*(ny_i*S[_j][ily_i][Cart::zzzzz]+ny_j*S[ily_j][_i][Cart::zzzzz]);
      S[_j][_i][Cart::zzzzzz] = gmb2*S[_j][_i][Cart::zzzzz] + fak*(nz_i*S[_j][ilz_i][Cart::zzzzz]+nz_j*S[ilz_j][_i][Cart::zzzzz]) + 5*term_zzzz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_gamma > 5)






            // data is now stored in unnormalized cartesian Gaussians in the multiarray
            // Now, weird-looking construction since multiarray is not accessible for product
            //              s px py pz dxz dyz dxy d3z2-r2 dx2-y2  f1  f2  f3  f4  f5  f6  f7  g1  g2  g3  g4  g5  g6  g7  g8  g9 
/////            int istart[] = {0, 1, 2, 3, 5,  6,  4,   7,     7,    12,  10, 11, 11, 10, 19, 15,  5, 25, 27, 23, 20, 25, 27, 23, 20}; //extend for g
/////            int istop[] =  {0, 1, 2, 3, 5,  6,  4,   9,     8,    17,  16, 18, 13, 14, 19, 17, 31, 33, 32 ,34, 30, 33, 32, 24, 31}; // extend for g

            int istart[] = { 0,  1,  1,  1,  4,  4,  4,  4,  4, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, //////////////
                            35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56 };   //////////////
            int istop[] =  { 0,  3,  3,  3,  9,  9,  9,  9,  9, 19, 19, 19, 19, 19, 19, 19, 34, 34, 34, 34, 34, 34, 34, 34, 34, //////////////
                            55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83 };   //////////////

           
            // which ones do we want to store
            int _offset_gw = _shell_gw->getOffset();
            int _offset_alpha = _shell_alpha->getOffset();
            int _offset_gamma = _shell_gamma->getOffset();

            // prepare transformation matrices
            int _ntrafo_gw = _shell_gw->getNumFunc() + _offset_gw;
            int _ntrafo_alpha = _shell_alpha->getNumFunc() + _offset_alpha;
            int _ntrafo_gamma = _shell_gamma->getNumFunc() + _offset_gamma;
            
            
            

            const Eigen::MatrixXd _trafo_gw = AOSuperMatrix::getTrafo(*itgw);
            const Eigen::MatrixXd _trafo_alpha = AOSuperMatrix::getTrafo(*italpha);
            const Eigen::MatrixXd _trafo_gamma = AOSuperMatrix::getTrafo(*itgamma);

       
            
            // transform from unnormalized cartesians to normalized sphericals
            // container with indices starting at zero
           
            double S_sph;
            
            for (int _i_alpha = _offset_alpha; _i_alpha < _ntrafo_alpha; _i_alpha++) {
                int alpha=_i_alpha-_offset_alpha;
                for (int _i_gw =  _offset_gw; _i_gw < _ntrafo_gw; _i_gw++) {
                    int g_w=_i_gw-_offset_gw;
                    for (int _i_gamma = _offset_gamma; _i_gamma < _ntrafo_gamma; _i_gamma++) {

                       
                        S_sph=0.0;
                        
                        for (int _i_alpha_t = istart[ _i_alpha ]; _i_alpha_t <= istop[ _i_alpha ]; _i_alpha_t++) {
                            for (int _i_gw_t = istart[ _i_gw ]; _i_gw_t <= istop[ _i_gw ]; _i_gw_t++) {
                                for (int _i_gamma_t = istart[ _i_gamma ]; _i_gamma_t <= istop[ _i_gamma ]; _i_gamma_t++) {

                                  
                                    S_sph+= S[ _i_alpha_t ][ _i_gw_t ][ _i_gamma_t ]  /////////////////
                                            * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_gw(_i_gw, _i_gw_t) * _trafo_gamma(_i_gamma, _i_gamma_t);



                                }
                            }
                        }
                        int _i_index = _shell_gamma->getNumFunc() * g_w + _i_gamma-_offset_gamma;

                        _subvector(alpha, _i_index) += S_sph;
                        
                    }
                }
            }

                    }
                }
            }

            return _does_contribute;
        }


    }
}
