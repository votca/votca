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
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Overload of uBLAS prod function with MKL/GSL implementations


#include <votca/xtp/threecenters.h>


using namespace votca::tools;

namespace votca {
    namespace xtp {
        namespace ub = boost::numeric::ublas;

 
        
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
         *      S,P,D,F,G   functions in DFT basis and 
         *      S,P,D,F,G,H,I   functions in GW  basis
         * 
         */
        bool TCrawMatrix::FillThreeCenterOLBlock(ub::matrix<double>& _subvector, AOShell* _shell_gw, AOShell* _shell_alpha, AOShell* _shell_gamma) {
	  //bool TCrawMatrix::FillThreeCenterOLBlock(ub::matrix<float>& _subvector, AOShell* _shell_gw, AOShell* _shell_alpha, AOShell* _shell_gamma) {
            const double pi = boost::math::constants::pi<double>();
              // get shell positions
            const vec& _pos_gw = _shell_gw->getPos();
            const vec& _pos_alpha = _shell_alpha->getPos();
            const vec& _pos_gamma = _shell_gamma->getPos();

            // shell info, only lmax tells how far to go
            int _lmax_gw = _shell_gw->getLmax();
            int _lmax_alpha = _shell_alpha->getLmax();
            int _lmax_gamma = _shell_gamma->getLmax();

            // set size of internal block for recursion
            int _ngw = this->getBlockSize(_lmax_gw);
            int _nalpha = this->getBlockSize(_lmax_alpha);
            int _ngamma = this->getBlockSize(_lmax_gamma);
            // definition of cutoff for contribution
            
            const double gwaccuracy = 1.e-9; // should become an OPTION
            
            bool _does_contribute = false;



 int nx[] = {
 0,
 1, 0, 0,
 2, 1, 1, 0, 0, 0,
 3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0,
 5, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
 6, 5, 5, 4, 4, 4, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0
 };

 int ny[] = {
 0,
 0, 1, 0,
 0, 1, 0, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0
 };

 int nz[] = {
 0,
 0, 0, 1,
 0, 0, 1, 0, 1, 2,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6
 };


 int i_less_x[] = {
  0,
  0,  0,  0,
  1,  2,  3,  0,  0,  0,
  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,
 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  0,  0,  0,  0,  0,
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,  0,  0,  0,  0,  0,  0,
 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,  0,  0,  0,  0,  0,  0,  0
 };

 int i_less_y[] = {
  0,
  0,  0,  0,
  0,  1,  0,  2,  3,  0,
  0,  4,  0,  5,  6,  0,  7,  8,  9,  0,
  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,  0,
  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,  0,
  0, 35,  0, 36, 37,  0, 38, 39, 40,  0, 41, 42, 43, 44,  0, 45, 46, 47, 48, 49,  0, 50, 51, 52, 53, 54, 55,  0
 };

 int i_less_z[] = {
  0,
  0,  0,  0,
  0,  0,  1,  0,  2,  3,
  0,  0,  4,  0,  5,  6,  0,  7,  8,  9,
  0,  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,
  0,  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,
  0,  0, 35,  0, 36, 37,  0, 38, 39, 40,  0, 41, 42, 43, 44,  0, 45, 46, 47, 48, 49,  0, 50, 51, 52, 53, 54, 55
 };



            typedef std::vector< AOGaussianPrimitive* >::iterator GaussianIterator;
                // iterate over Gaussians in this _shell_row
            for ( GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha){
            // iterate over Gaussians in this _shell_col
                const double& _decay_alpha = (*italpha)->decay;
            
                for ( GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma){
                    const double& _decay_gamma = (*itgamma)->decay;
                    // check third threshold
                    vec _diff = _pos_alpha - _pos_gamma;
                    double test = _decay_alpha * _decay_gamma * _diff*_diff;
                    
                    for ( GaussianIterator itgw = _shell_gw->firstGaussian(); itgw != _shell_gw->lastGaussian(); ++itgw){
            // get decay constants (this all is still valid only for uncontracted functions)
                        const double& _decay_gw = (*itgw)->decay;
            
 
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
            
                  
            vec gvv = fak2 * (_decay_alpha * _pos_alpha + _decay_gw * _pos_gw + _decay_gamma * _pos_gamma);
            vec gma = gvv - _pos_alpha;
            vec gmb = gvv - _pos_gamma;
            vec gmc = gvv - _pos_gw; 
            
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
            typedef boost::multi_array<double, 3> ma_type;
            typedef boost::multi_array_types::extent_range range;

            ma_type::extent_gen extents;
            ma_type S;
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

} 
            

            // data is now stored in unnormalized cartesian Gaussians in the multiarray
            // Now, weird-looking construction since multiarray is not accessible for ub::prod
            //              s px py pz dxz dyz dxy d3z2-r2 dx2-y2  f1  f2  f3  f4  f5  f6  f7  g1  g2  g3  g4  g5  g6  g7  g8  g9 
/////            int istart[] = {0, 1, 2, 3, 5,  6,  4,   7,     7,    12,  10, 11, 11, 10, 19, 15,  5, 25, 27, 23, 20, 25, 27, 23, 20}; //extend for g
/////            int istop[] =  {0, 1, 2, 3, 5,  6,  4,   9,     8,    17,  16, 18, 13, 14, 19, 17, 31, 33, 32 ,34, 30, 33, 32, 24, 31}; // extend for g

            int istart[] = { 0,  1,  1,  1,  4,  4,  4,  4,  4, 10, 10, 10, 10, 10, 10, 10, 20, 20, 20, 20, 20, 20, 20, 20, 20, //////////////
                            35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56 };   //////////////
            int istop[] =  { 0,  3,  3,  3,  9,  9,  9,  9,  9, 19, 19, 19, 19, 19, 19, 19, 34, 34, 34, 34, 34, 34, 34, 34, 34, //////////////
                            55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83, 83 };   //////////////

            // ub::vector<ub::matrix<double> >& _subvector
            // which ones do we want to store
            int _offset_gw = _shell_gw->getOffset();
            int _offset_alpha = _shell_alpha->getOffset();
            int _offset_gamma = _shell_gamma->getOffset();

            // prepare transformation matrices
            int _ntrafo_gw = _shell_gw->getNumFunc() + _offset_gw;
            int _ntrafo_alpha = _shell_alpha->getNumFunc() + _offset_alpha;
            int _ntrafo_gamma = _shell_gamma->getNumFunc() + _offset_gamma;

            ub::matrix<double> _trafo_gw = ub::zero_matrix<double>(_ntrafo_gw, _ngw);
            ub::matrix<double> _trafo_alpha = ub::zero_matrix<double>(_ntrafo_alpha, _nalpha);
            ub::matrix<double> _trafo_gamma = ub::zero_matrix<double>(_ntrafo_gamma, _ngamma);

            
            //std::vector<double> _contractions_alpha = (*italpha)->contraction;
            //std::vector<double> _contractions_gamma = (*itgamma)->contraction;
            //std::vector<double> _contractions_gw    = (*itgw)->contraction;
            
            // get transformation matrices
            this->getTrafo(_trafo_gw, _lmax_gw, _decay_gw, (*itgw)->contraction);
            this->getTrafo(_trafo_alpha, _lmax_alpha, _decay_alpha, (*italpha)->contraction);
            this->getTrafo(_trafo_gamma, _lmax_gamma, _decay_gamma, (*itgamma)->contraction);

            // transform from unnormalized cartesians to normalized sphericals
            // container with indices starting at zero
            //ma_type S_sph;
            //S_sph.resize(extents[range(_offset_alpha, _ntrafo_alpha) ][range(_offset_gw, _ntrafo_gw) ][range(_offset_gamma,_ntrafo_gamma )]);
            double S_sph;
            
            for (int _i_alpha = _offset_alpha; _i_alpha < _ntrafo_alpha; _i_alpha++) {
                int alpha=_i_alpha-_offset_alpha;
                for (int _i_gw =  _offset_gw; _i_gw < _ntrafo_gw; _i_gw++) {
                    int g_w=_i_gw-_offset_gw;
                    for (int _i_gamma = _offset_gamma; _i_gamma < _ntrafo_gamma; _i_gamma++) {

                        //S_sph[ _i_alpha ][  _i_gw ][ _i_gamma ] = 0.0;
                        S_sph=0.0;
                        
                        for (int _i_alpha_t = istart[ _i_alpha ]; _i_alpha_t <= istop[ _i_alpha ]; _i_alpha_t++) {
                            for (int _i_gw_t = istart[ _i_gw ]; _i_gw_t <= istop[ _i_gw ]; _i_gw_t++) {
                                for (int _i_gamma_t = istart[ _i_gamma ]; _i_gamma_t <= istop[ _i_gamma ]; _i_gamma_t++) {

                                    //S_sph[_i_alpha ][ _i_gw ][  _i_gamma ] += S[ _i_alpha_t + 1 ][ _i_gw_t + 1 ][ _i_gamma_t + 1]
                                    //        * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_gw(_i_gw, _i_gw_t) * _trafo_gamma(_i_gamma, _i_gamma_t);
///////////                                    S_sph+= S[ _i_alpha_t + 1 ][ _i_gw_t + 1 ][ _i_gamma_t + 1]
                                    S_sph+= S[ _i_alpha_t ][ _i_gw_t ][ _i_gamma_t ]  /////////////////
                                            * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_gw(_i_gw, _i_gw_t) * _trafo_gamma(_i_gamma, _i_gamma_t);



                                }
                            }
                        }
                        int _i_index = _shell_gamma->getNumFunc() * g_w + _i_gamma-_offset_gamma;

                        _subvector(alpha, _i_index) += S_sph;//[ _i_alpha ][ _i_gw ][ _i_gamma ];
                        
                    }
                }
            }

          

                    }
                }
            }

            return _does_contribute;


        }

     

        void TCrawMatrix::getTrafo(ub::matrix<double>& _trafo,const int _lmax, const double& _decay,const std::vector<double>& contractions) {
        // s-functions
        _trafo(0,0) = contractions[0]; // s
        ///         0    1  2  3    4  5  6  7  8  9   10  11  12  13  14  15  16  17  18  19       20    21    22    23    24    25    26    27    28    29    30    31    32    33    34 
        ///         s,   x, y, z,   xy xz yz xx yy zz, xxy xyy xyz xxz xzz yyz yzz xxx yyy zzz,    xxxy, xxxz, xxyy, xxyz, xxzz, xyyy, xyyz, xyzz, xzzz, yyyz, yyzz, yzzz, xxxx, yyyy, zzzz,
        // p-functions
        if (_lmax > 0) {
          //cout << _trafo_row.size1() << ":" << _trafo_row.size2() << endl;
          double factor = 2.*sqrt(_decay)*contractions[1];
          _trafo(1,3) = factor;  // Y 1,0
          _trafo(2,2) = factor;  // Y 1,-1
          _trafo(3,1) = factor;  // Y 1,1
        }

        // d-functions
        if (_lmax > 1) { // order of functions changed
          double factor = 2.*_decay*contractions[2];
          double factor_1 =  factor/sqrt(3.);
          _trafo(4,Cart::xx) = -factor_1;    // d3z2-r2 (dxx)
          _trafo(4,Cart::yy) = -factor_1;    // d3z2-r2 (dyy)  Y 2,0
          _trafo(4,Cart::zz) = 2.*factor_1;  // d3z2-r2 (dzz)

          _trafo(5,Cart::yz) = 2.*factor;     // dyz           Y 2,-1

          _trafo(6,Cart::xz) = 2.*factor;     // dxz           Y 2,1

          _trafo(7,Cart::xy) = 2.*factor;     // dxy           Y 2,-2

          _trafo(8,Cart::xx) = factor;       // dx2-y2 (dxx)   Y 2,2
          _trafo(8,Cart::yy) = -factor;      // dx2-y2 (dzz)
        }
       
        // f-functions
        if (_lmax > 2) { // order of functions changed
          double factor = 2.*pow(_decay,1.5)*contractions[3];
          double factor_1 = factor*2./sqrt(15.);
          double factor_2 = factor*sqrt(2.)/sqrt(5.);
          double factor_3 = factor*sqrt(2.)/sqrt(3.);

          _trafo(9,Cart::xxz) = -3.*factor_1;        // f1 (f??) xxz 13
          _trafo(9,Cart::yyz) = -3.*factor_1;        // f1 (f??) yyz 15        Y 3,0
          _trafo(9,Cart::zzz) = 2.*factor_1;         // f1 (f??) zzz 19

          _trafo(10,Cart::xxy) = -factor_2;          // f3 xxy 10
          _trafo(10,Cart::yyy) = -factor_2;          // f3 yyy 18   Y 3,-1
          _trafo(10,Cart::yzz) = 4.*factor_2;        // f3 yzz 16

          _trafo(11,Cart::xxx) = -factor_2;          // f2 xxx 17
          _trafo(11,Cart::xyy) = -factor_2;          // f2 xyy 11   Y 3,1
          _trafo(11,Cart::xzz) = 4.*factor_2;        // f2 xzz 14

          _trafo(12,Cart::xyz) = 4.*factor;          // f6 xyz 12     Y 3,-2

          _trafo(13,Cart::xxz) = 2.*factor;          // f7 (f??)   xxz   13
          _trafo(13,Cart::yyz) = -2.*factor;         // f7 (f??)   yyz   15   Y 3,2

          _trafo(14,Cart::xxy) = 3.*factor_3;        // f4 xxy 10
          _trafo(14,Cart::yyy) = -factor_3;          // f4 yyy 18   Y 3,-3

          _trafo(15,Cart::xxx) = factor_3;           // f5 (f??) xxx 17
          _trafo(15,Cart::xyy) = -3.*factor_3;       // f5 (f??) xyy 11     Y 3,3
        }

        // g-functions
        if (_lmax > 3) {
          double factor = 2./sqrt(3.)*_decay*_decay*contractions[4];
          double factor_1 = factor/sqrt(35.);
          double factor_2 = factor*4./sqrt(14.);
          double factor_3 = factor*2./sqrt(7.);
          double factor_4 = factor*2.*sqrt(2.);

          _trafo(16,Cart::xxxx) = 3.*factor_1;   /// Y 4,0
          _trafo(16,Cart::xxyy) = 6.*factor_1;
          _trafo(16,Cart::xxzz) = -24.*factor_1;
          _trafo(16,Cart::yyyy) = 3.*factor_1;
          _trafo(16,Cart::yyzz) = -24.*factor_1;
          _trafo(16,Cart::zzzz) = 8.*factor_1;

          _trafo(17,Cart::xxyz) = -3.*factor_2;  /// Y 4,-1
          _trafo(17,Cart::yyyz) = -3.*factor_2;
          _trafo(17,Cart::yzzz) = 4.*factor_2;

          _trafo(18,Cart::xxxz) = -3.*factor_2;  /// Y 4,1
          _trafo(18,Cart::xyyz) = -3.*factor_2;
          _trafo(18,Cart::xzzz) = 4.*factor_2;

          _trafo(19,Cart::xxxy) = -2.*factor_3;  /// Y 4,-2
          _trafo(19,Cart::xyyy) = -2.*factor_3;
          _trafo(19,Cart::xyzz) = 12.*factor_3;

          _trafo(20,Cart::xxxx) = -factor_3;     /// Y 4,2
          _trafo(20,Cart::xxzz) = 6.*factor_3;
          _trafo(20,Cart::yyyy) = factor_3;
          _trafo(20,Cart::yyzz) = -6.*factor_3;

          _trafo(21,Cart::xxyz) = 3.*factor_4;   /// Y 4,-3
          _trafo(21,Cart::yyyz) = -factor_4;

          _trafo(22,Cart::xxxz) = factor_4;      /// Y 4,3
          _trafo(22,Cart::xyyz) = -3.*factor_4;

          _trafo(23,Cart::xxxy) = 4.*factor;     /// Y 4,-4
          _trafo(23,Cart::xyyy) = -4.*factor;

          _trafo(24,Cart::xxxx) = factor;        /// Y 4,4
          _trafo(24,Cart::xxyy) = -6.*factor;
          _trafo(24,Cart::yyyy) = factor;
        }

        // h-functions
        if (_lmax > 4) {
          double factor = (2./3.)*pow(_decay,2.5)*contractions[5];
          double factor_1 = factor*2./sqrt(105.);
          double factor_2 = factor*2./sqrt(7.);
          double factor_3 = factor*sqrt(6.)/3.;
          double factor_4 = factor*2.*sqrt(3.);
          double factor_5 = factor*.2*sqrt(30.);

          _trafo(25,Cart::xxxxz) = 15.*factor_1;      /// Y 5,0
          _trafo(25,Cart::xxyyz) = 30.*factor_1;
          _trafo(25,Cart::xxzzz) = -40.*factor_1;
          _trafo(25,Cart::yyyyz) = 15.*factor_1;
          _trafo(25,Cart::yyzzz) = -40.*factor_1;
          _trafo(25,Cart::zzzzz) = 8.*factor_1;

          _trafo(26,Cart::xxxxy) = factor_2;          /// Y 5,-1
          _trafo(26,Cart::xxyyy) = 2.*factor_2;
          _trafo(26,Cart::xxyzz) = -12.*factor_2;
          _trafo(26,Cart::yyyyy) = factor_2;
          _trafo(26,Cart::yyyzz) = -12.*factor_2;
          _trafo(26,Cart::yzzzz) = 8.*factor_2;

          _trafo(27,Cart::xxxxx) = factor_2;          /// Y 5,1
          _trafo(27,Cart::xxxyy) = 2.*factor_2;
          _trafo(27,Cart::xxxzz) = -12.*factor_2;
          _trafo(27,Cart::xyyyy) = factor_2;
          _trafo(27,Cart::xyyzz) = -12.*factor_2;
          _trafo(27,Cart::xzzzz) = 8.*factor_2;

          _trafo(28,Cart::xxxyz) = -8.*factor;        /// Y 5,-2
          _trafo(28,Cart::xyyyz) = -8.*factor;
          _trafo(28,Cart::xyzzz) = 16.*factor;

          _trafo(29,Cart::xxxxz) = -4.*factor;        /// Y 5,2
          _trafo(29,Cart::xxzzz) = 8.*factor;
          _trafo(29,Cart::yyyyz) = 4.*factor;
          _trafo(29,Cart::yyzzz) = -8.*factor;

          _trafo(30,Cart::xxxxy) = -3.*factor_3;      /// Y 5,-3
          _trafo(30,Cart::xxyyy) = -2.*factor_3;
          _trafo(30,Cart::xxyzz) = 24.*factor_3;
          _trafo(30,Cart::yyyyy) = factor_3;
          _trafo(30,Cart::yyyzz) = -8.*factor_3;

          _trafo(31,Cart::xxxxx) = -factor_3;         /// Y 5,3
          _trafo(31,Cart::xxxyy) = 2.*factor_3;
          _trafo(31,Cart::xxxzz) = 8.*factor_3;
          _trafo(31,Cart::xyyyy) = 3.*factor_3;
          _trafo(31,Cart::xyyzz) = -24.*factor_3;

          _trafo(32,Cart::xxxyz) = 4.*factor_4;       /// Y 5,-4
          _trafo(32,Cart::xyyyz) = -4.*factor_4;

          _trafo(33,Cart::xxxxz) = factor_4;          /// Y 5,4
          _trafo(33,Cart::xxyyz) = -6.*factor_4;
          _trafo(33,Cart::yyyyz) = factor_4;

          _trafo(34,Cart::xxxxy) = 5.*factor_5;       /// Y 5,-5
          _trafo(34,Cart::xxyyy) = -10.*factor_5;
          _trafo(34,Cart::yyyyy) = factor_5;

          _trafo(35,Cart::xxxxx) = factor_5;          /// Y 5,5
          _trafo(35,Cart::xxxyy) = -10.*factor_5;
          _trafo(35,Cart::xyyyy) = 5.*factor_5;
        }

        // i-functions
        if (_lmax > 5) {
          double factor = (2./3.)*_decay*_decay*_decay*contractions[6];
          double factor_1 = factor*2./sqrt(1155.);
          double factor_2 = factor*4./sqrt(55.);
          double factor_3 = factor*sqrt(22.)/11.;
          double factor_4 = factor*2.*sqrt(165.)/55.;
          double factor_5 = factor*.4*sqrt(30.);
          double factor_6 = factor*.2*sqrt(10.);

          _trafo(36,Cart::xxxxxx) = -5.*factor_1;     /// Y 6,0
          _trafo(36,Cart::xxxxyy) = -15.*factor_1;
          _trafo(36,Cart::xxxxzz) = 90.*factor_1;
          _trafo(36,Cart::xxyyyy) = -15.*factor_1;
          _trafo(36,Cart::xxyyzz) = 180.*factor_1;
          _trafo(36,Cart::xxzzzz) = -120.*factor_1;
          _trafo(36,Cart::yyyyyy) = -5.*factor_1;
          _trafo(36,Cart::yyyyzz) = 90.*factor_1;
          _trafo(36,Cart::yyzzzz) = -120.*factor_1;
          _trafo(36,Cart::zzzzzz) = 16.*factor_1;

          _trafo(37,Cart::xxxxyz) = 5.*factor_2;      /// Y 6,-1
          _trafo(37,Cart::xxyyyz) = 10.*factor_2;
          _trafo(37,Cart::xxyzzz) = -20.*factor_2;
          _trafo(37,Cart::yyyyyz) = 5.*factor_2;
          _trafo(37,Cart::yyyzzz) = -20.*factor_2;
          _trafo(37,Cart::yzzzzz) = 8.*factor_2;

          _trafo(38,Cart::xxxxxz) = 5.*factor_2;      /// Y 6,1
          _trafo(38,Cart::xxxyyz) = 10.*factor_2;
          _trafo(38,Cart::xxxzzz) = -20.*factor_2;
          _trafo(38,Cart::xyyyyz) = 5.*factor_2;
          _trafo(38,Cart::xyyzzz) = -20.*factor_2;
          _trafo(38,Cart::xzzzzz) = 8.*factor_2;

          _trafo(39,Cart::xxxxxy) = 2.*factor_3;      /// Y 6,-2
          _trafo(39,Cart::xxxyyy) = 4.*factor_3;
          _trafo(39,Cart::xxxyzz) = -32.*factor_3;
          _trafo(39,Cart::xyyyyy) = 2.*factor_3;
          _trafo(39,Cart::xyyyzz) = -32.*factor_3;
          _trafo(39,Cart::xyzzzz) = 32.*factor_3;

          _trafo(40,Cart::xxxxxy) = factor_3;         /// Y 6,2
          _trafo(40,Cart::xxxxyy) = factor_3;
          _trafo(40,Cart::xxxxzz) = -16.*factor_3;
          _trafo(40,Cart::xxyyyy) = -factor_3;
          _trafo(40,Cart::xxzzzz) = 16.*factor_3;
          _trafo(40,Cart::yyyyyy) = -factor_3;
          _trafo(40,Cart::yyyyzz) = 16.*factor_3;
          _trafo(40,Cart::yyzzzz) = -16.*factor_3;

          _trafo(41,Cart::xxxxyz) = -18.*factor_3;    /// Y 6,-3
          _trafo(41,Cart::xxyyyz) = -12.*factor_3;
          _trafo(41,Cart::xxyzzz) = 48.*factor_3;
          _trafo(41,Cart::yyyyyz) = 6.*factor_3;
          _trafo(41,Cart::yyyzzz) = -16.*factor_3;

          _trafo(42,Cart::xxxxxz) = -6.*factor_3;     /// Y 6,3
          _trafo(42,Cart::xxxyyz) = 12.*factor_3;
          _trafo(42,Cart::xxxzzz) = 16.*factor_3;
          _trafo(42,Cart::xyyyyz) = 18.*factor_3;
          _trafo(42,Cart::xyyzzz) = -48.*factor_3;

          _trafo(43,Cart::xxxxxy) = -4.*factor_4;     /// Y 6,-4
          _trafo(43,Cart::xxxyzz) = 40.*factor_4;
          _trafo(43,Cart::xyyyyy) = 4.*factor_4;
          _trafo(43,Cart::xyyyzz) = -40.*factor_4;

          _trafo(44,Cart::xxxxxx) = -factor_4;        /// Y 6,4
          _trafo(44,Cart::xxxxyy) = 5.*factor_4;
          _trafo(44,Cart::xxxxzz) = 10.*factor_4;
          _trafo(44,Cart::xxyyyy) = 5.*factor_4;
          _trafo(44,Cart::xxyyzz) = -60.*factor_4;
          _trafo(44,Cart::yyyyyy) = -factor_4;
          _trafo(44,Cart::yyyyzz) = 10.*factor_4;

          _trafo(45,Cart::xxxxyz) = 5.*factor_5;      /// Y 6,-5
          _trafo(45,Cart::xxyyyz) = -10.*factor_5;
          _trafo(45,Cart::yyyyyz) = factor_5;

          _trafo(46,Cart::xxxxxz) = factor_5;         /// Y 6,5
          _trafo(46,Cart::xxxyyz) = -10.*factor_5;
          _trafo(46,Cart::xyyyyz) = 5.*factor_5;

          _trafo(47,Cart::xxxxxy) = 6.*factor_6;      /// Y 6,-6
          _trafo(47,Cart::xxxyyy) = -20.*factor_6;
          _trafo(47,Cart::xyyyyy) = 6.*factor_6;

          _trafo(48,Cart::xxxxxx) = factor_6;         /// Y 6,6
          _trafo(48,Cart::xxxxyy) = -15.*factor_6;
          _trafo(48,Cart::xxyyyy) = 15.*factor_6;
          _trafo(48,Cart::yyyyyy) = -factor_6;
        }

        return;
        }
        

        int TCrawMatrix::getBlockSize(const int _lmax) {
            int _block_size=-1;
            if (_lmax == 0) {
                _block_size = 1;
            } // s
            else if (_lmax == 1) {
                _block_size = 4;
            } // p
            else if (_lmax == 2) {
                _block_size = 10;
            } // d
            else if (_lmax == 3) {
                _block_size = 20;
            } // f
            else if (_lmax == 4) {
                _block_size = 35;
            } // g
            else if (_lmax == 5) { ////
                _block_size = 56; /////
            } // h
            else if (_lmax == 6) { /////
                _block_size = 84; /////
            } // i
            else if (_lmax == 7) { /////
                _block_size = 120; /////
            } // j
            else if (_lmax == 8) { /////
                _block_size = 165; /////
            } // k
            else{
                throw runtime_error("lmax for getBlocksize not known.");
            }

            return _block_size;
        }


    }
}

