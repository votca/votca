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
         *      S,P,D   functions in DFT basis and 
         *      S,P,D,F functions in GW  basis
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
            
            
            int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};

 // for alphabetical order

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
            
            const double gwaccuracy = 1.e-9; // should become an OPTION
            
            bool _does_contribute = false;
            
             typedef std::vector< AOGaussianPrimitive* >::iterator GaussianIterator;
                // iterate over Gaussians in this _shell_row
            for ( GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha){
            // iterate over Gaussians in this _shell_col
                const double _decay_alpha = (*italpha)->decay;
            
                for ( GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma){
                    const double& _decay_gamma = (*itgamma)->decay;
                    // check third threshold
                    vec _diff = _pos_alpha - _pos_gamma;
                    double test = _decay_alpha * _decay_gamma * _diff*_diff;
                    
                    for ( GaussianIterator itgw = _shell_gw->firstGaussian(); itgw != _shell_gw->lastGaussian(); ++itgw){
            // get decay constants (this all is still valid only for uncontracted functions)
                        const double _decay_gw = (*itgw)->decay;
            
 
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
            
            //double fak4=  4.0 * fak;
            
            double expo = _decay_alpha * _decay_gamma * (_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma)
                    + _decay_gamma * _decay_gw * (_pos_gamma - _pos_gw) * (_pos_gamma - _pos_gw)
                    + _decay_alpha * _decay_gw * (_pos_alpha - _pos_gw) * (_pos_alpha - _pos_gw);

            
            double prefak = pow(8.0 * _decay_alpha * _decay_gamma * _decay_gw / pi, 0.75) * pow(fak2, 1.5);

            double value = prefak * exp(-fak2 * expo);
            double fak3 = 3.0 * fak;
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
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type S;
            S.resize(extents[ range(0, _nalpha) ][ range(0, _ngw) ][ range(0, _ngamma)]);

            //cout << S.shape()[0]<< " : "<< S.shape()[1]<< " : "<< S.shape()[2]<<endl;
            
            S[0][0][0] = value; ////////////////

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




  if (_lmax_gamma > 0) {

    //Integrals     s - s - p
    S[0][0][Cart::x] = gmb0*S[0][0][0];
    S[0][0][Cart::y] = gmb1*S[0][0][0];
    S[0][0][Cart::z] = gmb2*S[0][0][0];
    //------------------------------------------------------

    //Integrals     s - p - p
    if (_lmax_gw > 0) {
      double term = fak*S[0][0][0];
      for (int _i =  1; _i < 4; _i++) {
        S[0][_i][Cart::x] = gmb0*S[0][_i][0] + nx[_i]*term;
        S[0][_i][Cart::y] = gmb1*S[0][_i][0] + ny[_i]*term;
        S[0][_i][Cart::z] = gmb2*S[0][_i][0] + nz[_i]*term;
      }
    }
    //------------------------------------------------------

    //Integrals     s - d - p     s - f - p     s - g - p     s - h - p     s - i - p
    for (int _i_gw = 2; _i_gw < _lmax_gw+1; _i_gw++) {
      for (int _i =  n_orbitals[_i_gw-1]; _i < n_orbitals[_i_gw]; _i++) {
        S[0][_i][Cart::x] = gmb0*S[0][_i][0] + nx[_i]*fak*S[0][i_less_x[_i]][0];
        S[0][_i][Cart::y] = gmb1*S[0][_i][0] + ny[_i]*fak*S[0][i_less_y[_i]][0];
        S[0][_i][Cart::z] = gmb2*S[0][_i][0] + nz[_i]*fak*S[0][i_less_z[_i]][0];
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_gamma > 0)


  if (_lmax_gamma > 1) {

    //Integrals     s - s - d
    double term = fak*S[0][0][0];
    S[0][0][Cart::xx] = gmb0*S[0][0][Cart::x] + term;
    S[0][0][Cart::xy] = gmb0*S[0][0][Cart::y];
    S[0][0][Cart::xz] = gmb0*S[0][0][Cart::z];
    S[0][0][Cart::yy] = gmb1*S[0][0][Cart::y] + term;
    S[0][0][Cart::yz] = gmb1*S[0][0][Cart::z];
    S[0][0][Cart::zz] = gmb2*S[0][0][Cart::z] + term;
    //------------------------------------------------------

    //Integrals     s - p - d     s - d - d     s - f - d     s - g - d     s - h - d     s - i - d
    for (int _i_gw = 1; _i_gw < _lmax_gw+1; _i_gw++) {
      for (int _i =  n_orbitals[_i_gw-1]; _i < n_orbitals[_i_gw]; _i++) {
        double term = fak*S[0][_i][0];
        S[0][_i][Cart::xx] = gmb0*S[0][_i][Cart::x] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::x] + term;
        S[0][_i][Cart::xy] = gmb0*S[0][_i][Cart::y] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::y];
        S[0][_i][Cart::xz] = gmb0*S[0][_i][Cart::z] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::z];
        S[0][_i][Cart::yy] = gmb1*S[0][_i][Cart::y] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::y] + term;
        S[0][_i][Cart::yz] = gmb1*S[0][_i][Cart::z] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::z];
        S[0][_i][Cart::zz] = gmb2*S[0][_i][Cart::z] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::z] + term;
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_gamma > 1)


  if (_lmax_gamma > 2) {

    //Integrals     s - s - f
    S[0][0][Cart::xxx] = gmb0*S[0][0][Cart::xx] + 2*fak*S[0][0][Cart::x];
    S[0][0][Cart::xxy] = gmb1*S[0][0][Cart::xx];
    S[0][0][Cart::xxz] = gmb2*S[0][0][Cart::xx];
    S[0][0][Cart::xyy] = gmb0*S[0][0][Cart::yy];
    S[0][0][Cart::xyz] = gmb0*S[0][0][Cart::yz];
    S[0][0][Cart::xzz] = gmb0*S[0][0][Cart::zz];
    S[0][0][Cart::yyy] = gmb1*S[0][0][Cart::yy] + 2*fak*S[0][0][Cart::y];
    S[0][0][Cart::yyz] = gmb2*S[0][0][Cart::yy];
    S[0][0][Cart::yzz] = gmb1*S[0][0][Cart::zz];
    S[0][0][Cart::zzz] = gmb2*S[0][0][Cart::zz] + 2*fak*S[0][0][Cart::z];
    //------------------------------------------------------

    //Integrals     s - p - f     s - d - f     s - f - f     s - g - f     s - h - f     s - i - f
    for (int _i_gw = 1; _i_gw < _lmax_gw+1; _i_gw++) {
      for (int _i =  n_orbitals[_i_gw-1]; _i < n_orbitals[_i_gw]; _i++) {
        double term_x = 2*fak*S[0][_i][Cart::x];
        double term_y = 2*fak*S[0][_i][Cart::y];
        double term_z = 2*fak*S[0][_i][Cart::z];
        S[0][_i][Cart::xxx] = gmb0*S[0][_i][Cart::xx] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::xx] + term_x;
        S[0][_i][Cart::xxy] = gmb1*S[0][_i][Cart::xx] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::xx];
        S[0][_i][Cart::xxz] = gmb2*S[0][_i][Cart::xx] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::xx];
        S[0][_i][Cart::xyy] = gmb0*S[0][_i][Cart::yy] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::yy];
        S[0][_i][Cart::xyz] = gmb0*S[0][_i][Cart::yz] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::yz];
        S[0][_i][Cart::xzz] = gmb0*S[0][_i][Cart::zz] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::zz];
        S[0][_i][Cart::yyy] = gmb1*S[0][_i][Cart::yy] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::yy] + term_y;
        S[0][_i][Cart::yyz] = gmb2*S[0][_i][Cart::yy] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::yy];
        S[0][_i][Cart::yzz] = gmb1*S[0][_i][Cart::zz] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::zz];
        S[0][_i][Cart::zzz] = gmb2*S[0][_i][Cart::zz] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::zz] + term_z;
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_gamma > 2)


  if (_lmax_gamma > 3) {

    //Integrals     s - s - g
    double term_xx = fak*S[0][0][Cart::xx];
    double term_yy = fak*S[0][0][Cart::yy];
    double term_zz = fak*S[0][0][Cart::zz];
    S[0][0][Cart::xxxx] = gmb0*S[0][0][Cart::xxx] + 3*term_xx;
    S[0][0][Cart::xxxy] = gmb1*S[0][0][Cart::xxx];
    S[0][0][Cart::xxxz] = gmb2*S[0][0][Cart::xxx];
    S[0][0][Cart::xxyy] = gmb0*S[0][0][Cart::xyy] + term_yy;
    S[0][0][Cart::xxyz] = gmb1*S[0][0][Cart::xxz];
    S[0][0][Cart::xxzz] = gmb0*S[0][0][Cart::xzz] + term_zz;
    S[0][0][Cart::xyyy] = gmb0*S[0][0][Cart::yyy];
    S[0][0][Cart::xyyz] = gmb0*S[0][0][Cart::yyz];
    S[0][0][Cart::xyzz] = gmb0*S[0][0][Cart::yzz];
    S[0][0][Cart::xzzz] = gmb0*S[0][0][Cart::zzz];
    S[0][0][Cart::yyyy] = gmb1*S[0][0][Cart::yyy] + 3*term_yy;
    S[0][0][Cart::yyyz] = gmb2*S[0][0][Cart::yyy];
    S[0][0][Cart::yyzz] = gmb1*S[0][0][Cart::yzz] + term_zz;
    S[0][0][Cart::yzzz] = gmb1*S[0][0][Cart::zzz];
    S[0][0][Cart::zzzz] = gmb2*S[0][0][Cart::zzz] + 3*term_zz;
    //------------------------------------------------------

    //Integrals     s - p - g     s - d - g     s - f - g     s - g - g     s - h - g     s - i - g
    for (int _i_gw = 1; _i_gw < _lmax_gw+1; _i_gw++) {
      for (int _i =  n_orbitals[_i_gw-1]; _i < n_orbitals[_i_gw]; _i++) {
        double term_xx = fak*S[0][_i][Cart::xx];
        double term_yy = fak*S[0][_i][Cart::yy];
        double term_zz = fak*S[0][_i][Cart::zz];
        S[0][_i][Cart::xxxx] = gmb0*S[0][_i][Cart::xxx] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::xxx] + 3*term_xx;
        S[0][_i][Cart::xxxy] = gmb1*S[0][_i][Cart::xxx] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::xxx];
        S[0][_i][Cart::xxxz] = gmb2*S[0][_i][Cart::xxx] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::xxx];
        S[0][_i][Cart::xxyy] = gmb0*S[0][_i][Cart::xyy] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::xyy] + term_yy;
        S[0][_i][Cart::xxyz] = gmb1*S[0][_i][Cart::xxz] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::xxz];
        S[0][_i][Cart::xxzz] = gmb0*S[0][_i][Cart::xzz] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::xzz] + term_zz;
        S[0][_i][Cart::xyyy] = gmb0*S[0][_i][Cart::yyy] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::yyy];
        S[0][_i][Cart::xyyz] = gmb0*S[0][_i][Cart::yyz] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::yyz];
        S[0][_i][Cart::xyzz] = gmb0*S[0][_i][Cart::yzz] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::yzz];
        S[0][_i][Cart::xzzz] = gmb0*S[0][_i][Cart::zzz] + nx[_i]*fak*S[0][i_less_x[_i]][Cart::zzz];
        S[0][_i][Cart::yyyy] = gmb1*S[0][_i][Cart::yyy] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::yyy] + 3*term_yy;
        S[0][_i][Cart::yyyz] = gmb2*S[0][_i][Cart::yyy] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::yyy];
        S[0][_i][Cart::yyzz] = gmb1*S[0][_i][Cart::yzz] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::yzz] + term_zz;
        S[0][_i][Cart::yzzz] = gmb1*S[0][_i][Cart::zzz] + ny[_i]*fak*S[0][i_less_y[_i]][Cart::zzz];
        S[0][_i][Cart::zzzz] = gmb2*S[0][_i][Cart::zzz] + nz[_i]*fak*S[0][i_less_z[_i]][Cart::zzz] + 3*term_zz;
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_gamma > 3)





  if (_lmax_alpha > 0) {

    //Integrals     p - s - s     p - p - s     p - d - s     p - f - s     p - g - s     p - h - s     p - i - s
    //              p - s - p     p - p - p     p - d - p     p - f - p     p - g - p     p - h - p     p - i - p
    //              p - s - d     p - p - d     p - d - d     p - f - d     p - g - d     p - h - d     p - i - d
    //              p - s - f     p - p - f     p - d - f     p - f - f     p - g - f     p - h - f     p - i - f
    //              p - s - g     p - p - g     p - d - g     p - f - g     p - g - g     p - h - g     p - i - g
    for (int _i = 0; _i < _ngamma; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      for (int _j = 0; _j < _ngw; _j++) {
        S[Cart::x][_j][_i] = gma0*S[0][_j][_i] + fak*(nx[_j]*S[0][i_less_x[_j]][_i]+nx_i*S[0][_j][ilx_i]);
        S[Cart::y][_j][_i] = gma1*S[0][_j][_i] + fak*(ny[_j]*S[0][i_less_y[_j]][_i]+ny_i*S[0][_j][ily_i]);
        S[Cart::z][_j][_i] = gma2*S[0][_j][_i] + fak*(nz[_j]*S[0][i_less_z[_j]][_i]+nz_i*S[0][_j][ilz_i]);
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_alpha > 0)


  if (_lmax_alpha > 1) {

    //Integrals     d - s - s     d - p - s     d - d - s     d - f - s     d - g - s     d - h - s     d - i - s
    //              d - s - p     d - p - p     d - d - p     d - f - p     d - g - p     d - h - p     d - i - p
    //              d - s - d     d - p - d     d - d - d     d - f - d     d - g - d     d - h - d     d - i - d
    //              d - s - f     d - p - f     d - d - f     d - f - f     d - g - f     d - h - f     d - i - f
    //              d - s - g     d - p - g     d - d - g     d - f - g     d - g - g     d - h - g     d - i - g
    for (int _i = 0; _i < _ngamma; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      for (int _j = 0; _j < _ngw; _j++) {
        double term = fak*S[0][_j][_i];
        S[Cart::xx][_j][_i] = gma0*S[Cart::x][_j][_i] + fak*(nx[_j]*S[Cart::x][i_less_x[_j]][_i]+nx_i*S[Cart::x][_j][ilx_i]) + term;
        S[Cart::xy][_j][_i] = gma0*S[Cart::y][_j][_i] + fak*(nx[_j]*S[Cart::y][i_less_x[_j]][_i]+nx_i*S[Cart::y][_j][ilx_i]);
        S[Cart::xz][_j][_i] = gma0*S[Cart::z][_j][_i] + fak*(nx[_j]*S[Cart::z][i_less_x[_j]][_i]+nx_i*S[Cart::z][_j][ilx_i]);
        S[Cart::yy][_j][_i] = gma1*S[Cart::y][_j][_i] + fak*(ny[_j]*S[Cart::y][i_less_y[_j]][_i]+ny_i*S[Cart::y][_j][ily_i]) + term;
        S[Cart::yz][_j][_i] = gma1*S[Cart::z][_j][_i] + fak*(ny[_j]*S[Cart::z][i_less_y[_j]][_i]+ny_i*S[Cart::z][_j][ily_i]);
        S[Cart::zz][_j][_i] = gma2*S[Cart::z][_j][_i] + fak*(nz[_j]*S[Cart::z][i_less_z[_j]][_i]+nz_i*S[Cart::z][_j][ilz_i]) + term;
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_alpha > 1)


  if (_lmax_alpha > 2) {

    //Integrals     f - s - s     f - p - s     f - d - s     f - f - s     f - g - s     f - h - s     f - i - s
    //              f - s - p     f - p - p     f - d - p     f - f - p     f - g - p     f - h - p     f - i - p
    //              f - s - d     f - p - d     f - d - d     f - f - d     f - g - d     f - h - d     f - i - d
    //              f - s - f     f - p - f     f - d - f     f - f - f     f - g - f     f - h - f     f - i - f
    //              f - s - g     f - p - g     f - d - g     f - f - g     f - g - g     f - h - g     f - i - g
    for (int _i = 0; _i < _ngamma; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      for (int _j = 0; _j < _ngw; _j++) {
        double term_x = 2*fak*S[Cart::x][_j][_i];
        double term_y = 2*fak*S[Cart::y][_j][_i];
        double term_z = 2*fak*S[Cart::z][_j][_i];
        S[Cart::xxx][_j][_i] = gma0*S[Cart::xx][_j][_i] + fak*(nx[_j]*S[Cart::xx][i_less_x[_j]][_i]+nx_i*S[Cart::xx][_j][ilx_i]) + term_x;
        S[Cart::xxy][_j][_i] = gma1*S[Cart::xx][_j][_i] + fak*(ny[_j]*S[Cart::xx][i_less_y[_j]][_i]+ny_i*S[Cart::xx][_j][ily_i]);
        S[Cart::xxz][_j][_i] = gma2*S[Cart::xx][_j][_i] + fak*(nz[_j]*S[Cart::xx][i_less_z[_j]][_i]+nz_i*S[Cart::xx][_j][ilz_i]);
        S[Cart::xyy][_j][_i] = gma0*S[Cart::yy][_j][_i] + fak*(nx[_j]*S[Cart::yy][i_less_x[_j]][_i]+nx_i*S[Cart::yy][_j][ilx_i]);
        S[Cart::xyz][_j][_i] = gma0*S[Cart::yz][_j][_i] + fak*(nx[_j]*S[Cart::yz][i_less_x[_j]][_i]+nx_i*S[Cart::yz][_j][ilx_i]);
        S[Cart::xzz][_j][_i] = gma0*S[Cart::zz][_j][_i] + fak*(nx[_j]*S[Cart::zz][i_less_x[_j]][_i]+nx_i*S[Cart::zz][_j][ilx_i]);
        S[Cart::yyy][_j][_i] = gma1*S[Cart::yy][_j][_i] + fak*(ny[_j]*S[Cart::yy][i_less_y[_j]][_i]+ny_i*S[Cart::yy][_j][ily_i]) + term_y;
        S[Cart::yyz][_j][_i] = gma2*S[Cart::yy][_j][_i] + fak*(nz[_j]*S[Cart::yy][i_less_z[_j]][_i]+nz_i*S[Cart::yy][_j][ilz_i]);
        S[Cart::yzz][_j][_i] = gma1*S[Cart::zz][_j][_i] + fak*(ny[_j]*S[Cart::zz][i_less_y[_j]][_i]+ny_i*S[Cart::zz][_j][ily_i]);
        S[Cart::zzz][_j][_i] = gma2*S[Cart::zz][_j][_i] + fak*(nz[_j]*S[Cart::zz][i_less_z[_j]][_i]+nz_i*S[Cart::zz][_j][ilz_i]) + term_z;
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_alpha > 2)


  if (_lmax_alpha > 3) {

    //Integrals     g - s - s     g - p - s     g - d - s     g - f - s     g - g - s     g - h - s     g - i - s
    //              g - s - p     g - p - p     g - d - p     g - f - p     g - g - p     g - h - p     g - i - p
    //              g - s - d     g - p - d     g - d - d     g - f - d     g - g - d     g - h - d     g - i - d
    //              g - s - f     g - p - f     g - d - f     g - f - f     g - g - f     g - h - f     g - i - f
    //              g - s - g     g - p - g     g - d - g     g - f - g     g - g - g     g - h - g     g - i - g
    for (int _i = 0; _i < _ngamma; _i++) {
      int nx_i = nx[_i];
      int ny_i = ny[_i];
      int nz_i = nz[_i];
      int ilx_i = i_less_x[_i];
      int ily_i = i_less_y[_i];
      int ilz_i = i_less_z[_i];
      for (int _j = 0; _j < _ngw; _j++) {
        int nx_j = nx[_j];
        int ny_j = ny[_j];
        int nz_j = nz[_j];
        int ilx_j = i_less_x[_j];
        int ily_j = i_less_y[_j];
        int ilz_j = i_less_z[_j];
        double term_xx = fak*S[Cart::xx][_j][_i];
        double term_yy = fak*S[Cart::yy][_j][_i];
        double term_zz = fak*S[Cart::zz][_j][_i];
        S[Cart::xxxx][_j][_i] = gma0*S[Cart::xxx][_j][_i] + fak*(nx_j*S[Cart::xxx][ilx_j][_i]+nx_i*S[Cart::xxx][_j][ilx_i]) + 3*term_xx;
        S[Cart::xxxy][_j][_i] = gma1*S[Cart::xxx][_j][_i] + fak*(ny_j*S[Cart::xxx][ily_j][_i]+ny_i*S[Cart::xxx][_j][ily_i]);
        S[Cart::xxxz][_j][_i] = gma2*S[Cart::xxx][_j][_i] + fak*(nz_j*S[Cart::xxx][ilz_j][_i]+nz_i*S[Cart::xxx][_j][ilz_i]);
        S[Cart::xxyy][_j][_i] = gma0*S[Cart::xyy][_j][_i] + fak*(nx_j*S[Cart::xyy][ilx_j][_i]+nx_i*S[Cart::xyy][_j][ilx_i]) + term_yy;
        S[Cart::xxyz][_j][_i] = gma1*S[Cart::xxz][_j][_i] + fak*(ny_j*S[Cart::xxz][ily_j][_i]+ny_i*S[Cart::xxz][_j][ily_i]);
        S[Cart::xxzz][_j][_i] = gma0*S[Cart::xzz][_j][_i] + fak*(nx_j*S[Cart::xzz][ilx_j][_i]+nx_i*S[Cart::xzz][_j][ilx_i]) + term_zz;
        S[Cart::xyyy][_j][_i] = gma0*S[Cart::yyy][_j][_i] + fak*(nx_j*S[Cart::yyy][ilx_j][_i]+nx_i*S[Cart::yyy][_j][ilx_i]);
        S[Cart::xyyz][_j][_i] = gma0*S[Cart::yyz][_j][_i] + fak*(nx_j*S[Cart::yyz][ilx_j][_i]+nx_i*S[Cart::yyz][_j][ilx_i]);
        S[Cart::xyzz][_j][_i] = gma0*S[Cart::yzz][_j][_i] + fak*(nx_j*S[Cart::yzz][ilx_j][_i]+nx_i*S[Cart::yzz][_j][ilx_i]);
        S[Cart::xzzz][_j][_i] = gma0*S[Cart::zzz][_j][_i] + fak*(nx_j*S[Cart::zzz][ilx_j][_i]+nx_i*S[Cart::zzz][_j][ilx_i]);
        S[Cart::yyyy][_j][_i] = gma1*S[Cart::yyy][_j][_i] + fak*(ny_j*S[Cart::yyy][ily_j][_i]+ny_i*S[Cart::yyy][_j][ily_i]) + 3*term_yy;
        S[Cart::yyyz][_j][_i] = gma2*S[Cart::yyy][_j][_i] + fak*(nz_j*S[Cart::yyy][ilz_j][_i]+nz_i*S[Cart::yyy][_j][ilz_i]);
        S[Cart::yyzz][_j][_i] = gma1*S[Cart::yzz][_j][_i] + fak*(ny_j*S[Cart::yzz][ily_j][_i]+ny_i*S[Cart::yzz][_j][ily_i]) + term_zz;
        S[Cart::yzzz][_j][_i] = gma1*S[Cart::zzz][_j][_i] + fak*(ny_j*S[Cart::zzz][ily_j][_i]+ny_i*S[Cart::zzz][_j][ily_i]);
        S[Cart::zzzz][_j][_i] = gma2*S[Cart::zzz][_j][_i] + fak*(nz_j*S[Cart::zzz][ilz_j][_i]+nz_i*S[Cart::zzz][_j][ilz_i]) + 3*term_zz;
      }
    }
    //------------------------------------------------------

  } // end if (_lmax_alpha > 3)
            
            

            // data is now stored in unnormalized cartesian Gaussians in the multiarray
            // Now, weird-looking construction since multiarray is not accessible for ub::prod
                //              s   pz        py      px      d3z2-r2   dyz        dxz       dxy dx2-y2  Y3,0  Y3,-1  Y3,1  Y3,-2  Y3,2  Y3,-3  Y3,3    Y4,0      Y4,-1    Y4,1         Y4,-2       Y4,2          Y4,-3     Y4,3        Y4,-4 Y4,4 
            int istart[] = {0, Cart::z, Cart::y, Cart::x, Cart::xx,  Cart::yz, Cart::xz, Cart::xx,   13,     10,   11,    12,   13,     10,  11,  Cart::xxxx,  Cart::xxyz, Cart::xxxz, Cart::xxxy, Cart::xxxx, Cart::xxyz, Cart::xxxz,  Cart::xxxy, Cart::xxxx}; //extend for g
            int istop[] =  {0, Cart::z, Cart::y, Cart::x, Cart::zz,  Cart::yz, Cart::xz, Cart::yy,   19,     18,   17,    12,   15,     18,  17,  Cart::zzzz , Cart::yzzz, Cart::xzzz, Cart::xyzz ,Cart::yyzz, Cart::yyyz, Cart::xyyz,  Cart::xyyy, Cart::yyyy}; // extend for g

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
                                    S_sph+= S[ _i_alpha_t  ][ _i_gw_t  ][ _i_gamma_t ]
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
      
        // p-functions
        if ( _lmax > 0 ){
            //cout << _trafo_row.size1() << ":" << _trafo_row.size2() << endl;
            double factor = 2.*sqrt(_decay)*contractions[1];
            _trafo(1,Cart::z) = factor;  // Y 1,0
            _trafo(2,Cart::y) = factor;  // Y 1,-1
            _trafo(3,Cart::x) = factor;  // Y 1,1
        }

         // d-functions
        if ( _lmax > 1 ) { // order of functions changed
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
        if ( _lmax > 2 ) { // order of functions changed
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
        if ( _lmax > 3 ) {
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

