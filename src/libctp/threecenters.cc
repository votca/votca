/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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


#include <votca/ctp/threecenters.h>

using namespace std;
using namespace votca::tools;

namespace votca {
    namespace ctp {
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
            
             typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
                // iterate over Gaussians in this _shell_row
            for ( GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha){
            // iterate over Gaussians in this _shell_col
                const double& _decay_alpha = (*italpha)->decay;
            
                for ( GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma){
                    const double& _decay_gamma = (*itgamma)->decay;
                    
                    for ( GaussianIterator itgw = _shell_gw->firstGaussian(); itgw != _shell_gw->lastGaussian(); ++itgw){
            // get decay constants (this all is still valid only for uncontracted functions)
                        const double& _decay_gw = (*itgw)->decay;
            
 
            double threshold = -(_decay_alpha + _decay_gamma + _decay_gw) * log(gwaccuracy);

            // check first threshold
            vec _diff = _pos_alpha - _pos_gw;
            
            double test = _decay_alpha * _decay_gw * _diff*_diff;
            if (test > threshold) { continue; }

            // check second threshold
            _diff = _pos_gamma - _pos_gw;
            test += _decay_gamma * _decay_gw * _diff*_diff;
            if (test > threshold) { continue; }

            // check third threshold
            _diff = _pos_alpha - _pos_gamma;
            test += _decay_alpha * _decay_gamma * _diff*_diff;
            if (test > threshold) { continue; }

            // if all threshold test are passed, start evaluating

            // some helpers
            double fak = 0.5 / (_decay_alpha + _decay_gw + _decay_gamma);
            double fak2 = 2.0 * fak;
            double fak3 = 3.0 * fak;
            double fak4=  4.0 * fak;

            double _dist1=(_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma);
            double _dist2=(_pos_gamma - _pos_gw) * (_pos_gamma - _pos_gw);
            double _dist3=(_pos_alpha - _pos_gw) * (_pos_alpha - _pos_gw);
            
                  
            vec gvv = fak2 * (_decay_alpha * _pos_alpha + _decay_gw * _pos_gw + _decay_gamma * _pos_gamma);
            vec gma = gvv - _pos_alpha;
            vec gmb = gvv - _pos_gamma;
            vec gmc = gvv - _pos_gw; 
            
            double gma0 = 0.0;
            double gmb0 = 0.0;
            double gmc0 = 0.0;

            double gma1 = 0.0;
            double gmb1 = 0.0;
            double gmc1 = 0.0;

            double gma2 = 0.0;
            double gmb2 = 0.0;
            double gmc2 = 0.0;

            if ((_dist1 + _dist2 + _dist3)>0.01){
          
           
            gma0 = gma.getX();
            gmb0 = gmb.getX();
            gmc0 = gmc.getX();

            gma1 = gma.getY();
            gmb1 = gmb.getY();
            gmc1 = gmc.getY();

            gma2 = gma.getZ();
            gmb2 = gmb.getZ();
            gmc2 = gmc.getZ();

            }
            // get s-s-s element
            double expo = _decay_alpha * _decay_gamma * (_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma)
                    + _decay_gamma * _decay_gw * (_pos_gamma - _pos_gw) * (_pos_gamma - _pos_gw)
                    + _decay_alpha * _decay_gw * (_pos_alpha - _pos_gw) * (_pos_alpha - _pos_gw);

            const double pi = boost::math::constants::pi<double>();
            double prefak = pow(8.0 * _decay_alpha * _decay_gamma * _decay_gw / pi, 0.75) * pow(fak2, 1.5);

            double value = prefak * exp(-fak2 * expo);

            // check if it contributes
            if (value < gwaccuracy) { continue; }

            _does_contribute = true;
            // if it does, go on and create multiarray
            typedef boost::multi_array<double, 3> ma_type;
            typedef boost::multi_array_types::extent_range range;
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type S;
            S.resize(extents[ range(1, _nalpha + 1) ][ range(1, _ngw + 1) ][ range(1, _ngamma + 1)]);

            // now fill s-s-s element
            S[1][1][1] = value;

            // s-p-s elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 1 && _lmax_gamma >= 0) {
                S[1][2][1] = gmc0 * S[1][1][1];
                S[1][3][1] = gmc1 * S[1][1][1];
                S[1][4][1] = gmc2 * S[1][1][1];
            }

            // s-d-s elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 2 && _lmax_gamma >= 0) {
                S[1][8][1] = gmc0 * S[1][2][1] + fak * S[1][1][1];
                S[1][5][1] = gmc1 * S[1][2][1];
                S[1][6][1] = gmc2 * S[1][2][1];
                S[1][9][1] = gmc1 * S[1][3][1] + fak * S[1][1][1];
                S[1][7][1] = gmc2 * S[1][3][1];
                S[1][10][1] = gmc2 * S[1][4][1] + fak * S[1][1][1];
            }

            // p-s-s elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 0 && _lmax_gamma >= 0) {
                S[2][1][1] = gma0 * S[1][1][1];
                S[3][1][1] = gma1 * S[1][1][1];
                S[4][1][1] = gma2 * S[1][1][1];
            }

            // p-p-s elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 1 && _lmax_gamma >= 0) {
                S[2][2][1] = gma0 * S[1][2][1] + fak * S[1][1][1];
                S[3][2][1] = gma1 * S[1][2][1];
                S[4][2][1] = gma2 * S[1][2][1];
                S[2][3][1] = gma0 * S[1][3][1];
                S[3][3][1] = gma1 * S[1][3][1] + fak * S[1][1][1];
                S[4][3][1] = gma2 * S[1][3][1];
                S[2][4][1] = gma0 * S[1][4][1];
                S[3][4][1] = gma1 * S[1][4][1];
                S[4][4][1] = gma2 * S[1][4][1] + fak * S[1][1][1];
            }

            // p-d-s elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 2 && _lmax_gamma >= 0) {
                S[2][5][1] = gma0 * S[1][5][1] + fak * S[1][3][1];
                S[3][5][1] = gma1 * S[1][5][1] + fak * S[1][2][1];
                S[4][5][1] = gma2 * S[1][5][1];
                S[2][6][1] = gma0 * S[1][6][1] + fak * S[1][4][1];
                S[3][6][1] = gma1 * S[1][6][1];
                S[4][6][1] = gma2 * S[1][6][1] + fak * S[1][2][1];
                S[2][7][1] = gma0 * S[1][7][1];
                S[3][7][1] = gma1 * S[1][7][1] + fak * S[1][4][1];
                S[4][7][1] = gma2 * S[1][7][1] + fak * S[1][3][1];
                S[2][8][1] = gma0 * S[1][8][1] + fak2 * S[1][2][1];
                S[3][8][1] = gma1 * S[1][8][1];
                S[4][8][1] = gma2 * S[1][8][1];
                S[2][9][1] = gma0 * S[1][9][1];
                S[3][9][1] = gma1 * S[1][9][1] + fak2 * S[1][3][1];
                S[4][9][1] = gma2 * S[1][9][1];
                S[2][10][1] = gma0 * S[1][10][1];
                S[3][10][1] = gma1 * S[1][10][1];
                S[4][10][1] = gma2 * S[1][10][1] + fak2 * S[1][4][1];
            }

            // s-s-p elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 0 && _lmax_gamma >= 1) {
                S[1][1][2] = gmb0 * S[1][1][1];
                S[1][1][3] = gmb1 * S[1][1][1];
                S[1][1][4] = gmb2 * S[1][1][1];
            }

            // s-p-p elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 1 && _lmax_gamma >= 1) {
                S[1][2][2] = gmc0 * S[1][1][2] + fak * S[1][1][1];
                S[1][3][2] = gmc1 * S[1][1][2];
                S[1][4][2] = gmc2 * S[1][1][2];
                S[1][2][3] = gmc0 * S[1][1][3];
                S[1][3][3] = gmc1 * S[1][1][3] + fak * S[1][1][1];
                S[1][4][3] = gmc2 * S[1][1][3];
                S[1][2][4] = gmc0 * S[1][1][4];
                S[1][3][4] = gmc1 * S[1][1][4];
                S[1][4][4] = gmc2 * S[1][1][4] + fak * S[1][1][1];
            }

            // s-d-p elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 2 && _lmax_gamma >= 1) {
                S[1][8][2] = gmc0 * S[1][2][2] + fak * (S[1][1][2] + S[1][2][1]);
                S[1][5][2] = gmc1 * S[1][2][2];
                S[1][6][2] = gmc2 * S[1][2][2];
                S[1][8][3] = gmc0 * S[1][2][3] + fak * S[1][1][3];
                S[1][5][3] = gmc1 * S[1][2][3] + fak * S[1][2][1];
                S[1][6][3] = gmc2 * S[1][2][3];
                S[1][8][4] = gmc0 * S[1][2][4] + fak * S[1][1][4];
                S[1][5][4] = gmc1 * S[1][2][4];
                S[1][6][4] = gmc2 * S[1][2][4] + fak * S[1][2][1];
                S[1][9][2] = gmc1 * S[1][3][2] + fak * S[1][1][2];
                S[1][7][2] = gmc2 * S[1][3][2];
                S[1][9][3] = gmc1 * S[1][3][3] + fak * (S[1][1][3] + S[1][3][1]);
                S[1][7][3] = gmc2 * S[1][3][3];
                S[1][9][4] = gmc1 * S[1][3][4] + fak * S[1][1][4];
                S[1][7][4] = gmc2 * S[1][3][4] + fak * S[1][3][1];
                S[1][10][2] = gmc2 * S[1][4][2] + fak * S[1][1][2];
                S[1][10][3] = gmc2 * S[1][4][3] + fak * S[1][1][3];
                S[1][10][4] = gmc2 * S[1][4][4] + fak * (S[1][1][4] + S[1][4][1]);
            }

            // p-s-p elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 0 && _lmax_gamma >= 1) {
                S[2][1][2] = gma0 * S[1][1][2] + fak * S[1][1][1];
                S[3][1][2] = gma1 * S[1][1][2];
                S[4][1][2] = gma2 * S[1][1][2];
                S[2][1][3] = gma0 * S[1][1][3];
                S[3][1][3] = gma1 * S[1][1][3] + fak * S[1][1][1];
                S[4][1][3] = gma2 * S[1][1][3];
                S[2][1][4] = gma0 * S[1][1][4];
                S[3][1][4] = gma1 * S[1][1][4];
                S[4][1][4] = gma2 * S[1][1][4] + fak * S[1][1][1];
            }

            // p-p-p elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 1 && _lmax_gamma >= 1) {
                S[2][2][2] = gma0 * S[1][2][2] + fak * (S[1][1][2] + S[1][2][1]);
                S[3][2][2] = gma1 * S[1][2][2];
                S[4][2][2] = gma2 * S[1][2][2];
                S[2][2][3] = gma0 * S[1][2][3] + fak * S[1][1][3];
                S[3][2][3] = gma1 * S[1][2][3] + fak * S[1][2][1];
                S[4][2][3] = gma2 * S[1][2][3];
                S[2][2][4] = gma0 * S[1][2][4] + fak * S[1][1][4];
                S[3][2][4] = gma1 * S[1][2][4];
                S[4][2][4] = gma2 * S[1][2][4] + fak * S[1][2][1];
                S[2][3][2] = gma0 * S[1][3][2] + fak * S[1][3][1];
                S[3][3][2] = gma1 * S[1][3][2] + fak * S[1][1][2];
                S[4][3][2] = gma2 * S[1][3][2];
                S[2][3][3] = gma0 * S[1][3][3];
                S[3][3][3] = gma1 * S[1][3][3] + fak * (S[1][1][3] + S[1][3][1]);
                S[4][3][3] = gma2 * S[1][3][3];
                S[2][3][4] = gma0 * S[1][3][4];
                S[3][3][4] = gma1 * S[1][3][4] + fak * S[1][1][4];
                S[4][3][4] = gma2 * S[1][3][4] + fak * S[1][3][1];
                S[2][4][2] = gma0 * S[1][4][2] + fak * S[1][4][1];
                S[3][4][2] = gma1 * S[1][4][2];
                S[4][4][2] = gma2 * S[1][4][2] + fak * S[1][1][2];
                S[2][4][3] = gma0 * S[1][4][3];
                S[3][4][3] = gma1 * S[1][4][3] + fak * S[1][4][1];
                S[4][4][3] = gma2 * S[1][4][3] + fak * S[1][1][3];
                S[2][4][4] = gma0 * S[1][4][4];
                S[3][4][4] = gma1 * S[1][4][4];
                S[4][4][4] = gma2 * S[1][4][4] + fak * (S[1][1][4] + S[1][4][1]);
            }

            // p-d-p elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 2 && _lmax_gamma >= 1) {
                S[2][5][2] = gma0 * S[1][5][2] + fak * (S[1][3][2] + S[1][5][1]);
                S[3][5][2] = gma1 * S[1][5][2] + fak * S[1][2][2];
                S[4][5][2] = gma2 * S[1][5][2];
                S[2][5][3] = gma0 * S[1][5][3] + fak * S[1][3][3];
                S[3][5][3] = gma1 * S[1][5][3] + fak * (S[1][2][3] + S[1][5][1]);
                S[4][5][3] = gma2 * S[1][5][3];
                S[2][5][4] = gma0 * S[1][5][4] + fak * S[1][3][4];
                S[3][5][4] = gma1 * S[1][5][4] + fak * S[1][2][4];
                S[4][5][4] = gma2 * S[1][5][4] + fak * S[1][5][1];
                S[2][6][2] = gma0 * S[1][6][2] + fak * (S[1][4][2] + S[1][6][1]);
                S[3][6][2] = gma1 * S[1][6][2];
                S[4][6][2] = gma2 * S[1][6][2] + fak * S[1][2][2];
                S[2][6][3] = gma0 * S[1][6][3] + fak * S[1][4][3];
                S[3][6][3] = gma1 * S[1][6][3] + fak * S[1][6][1];
                S[4][6][3] = gma2 * S[1][6][3] + fak * S[1][2][3];
                S[2][6][4] = gma0 * S[1][6][4] + fak * S[1][4][4];
                S[3][6][4] = gma1 * S[1][6][4];
                S[4][6][4] = gma2 * S[1][6][4] + fak * (S[1][2][4] + S[1][6][1]);
                S[2][7][2] = gma0 * S[1][7][2] + fak * S[1][7][1];
                S[3][7][2] = gma1 * S[1][7][2] + fak * S[1][4][2];
                S[4][7][2] = gma2 * S[1][7][2] + fak * S[1][3][2];
                S[2][7][3] = gma0 * S[1][7][3];
                S[3][7][3] = gma1 * S[1][7][3] + fak * (S[1][4][3] + S[1][7][1]);
                S[4][7][3] = gma2 * S[1][7][3] + fak * S[1][3][3];
                S[2][7][4] = gma0 * S[1][7][4];
                S[3][7][4] = gma1 * S[1][7][4] + fak * S[1][4][4];
                S[4][7][4] = gma2 * S[1][7][4] + fak * (S[1][3][4] + S[1][7][1]);
                S[2][8][2] = gma0 * S[1][8][2] + fak * (S[1][2][2] + S[1][2][2] + S[1][8][1]);
                S[3][8][2] = gma1 * S[1][8][2];
                S[4][8][2] = gma2 * S[1][8][2];
                S[2][8][3] = gma0 * S[1][8][3] + fak2 * S[1][2][3];
                S[3][8][3] = gma1 * S[1][8][3] + fak * S[1][8][1];
                S[4][8][3] = gma2 * S[1][8][3];
                S[2][8][4] = gma0 * S[1][8][4] + fak * (2 * S[1][2][4]);
                S[3][8][4] = gma1 * S[1][8][4];
                S[4][8][4] = gma2 * S[1][8][4] + fak * S[1][8][1];
                S[2][9][2] = gma0 * S[1][9][2] + fak * S[1][9][1];
                S[3][9][2] = gma1 * S[1][9][2] + fak2 * S[1][3][2];
                S[4][9][2] = gma2 * S[1][9][2];
                S[2][9][3] = gma0 * S[1][9][3];
                S[3][9][3] = gma1 * S[1][9][3] + fak * (S[1][3][3] + S[1][3][3] + S[1][9][1]);
                S[4][9][3] = gma2 * S[1][9][3];
                S[2][9][4] = gma0 * S[1][9][4];
                S[3][9][4] = gma1 * S[1][9][4] + fak2 * S[1][3][4];
                S[4][9][4] = gma2 * S[1][9][4] + fak * S[1][9][1];
                S[2][10][2] = gma0 * S[1][10][2] + fak * S[1][10][1];
                S[3][10][2] = gma1 * S[1][10][2];
                S[4][10][2] = gma2 * S[1][10][2] + fak2 * S[1][4][2];
                S[2][10][3] = gma0 * S[1][10][3];
                S[3][10][3] = gma1 * S[1][10][3] + fak * S[1][10][1];
                S[4][10][3] = gma2 * S[1][10][3] + fak2 * S[1][4][3];
                S[2][10][4] = gma0 * S[1][10][4];
                S[3][10][4] = gma1 * S[1][10][4];
                S[4][10][4] = gma2 * S[1][10][4] + fak * (S[1][4][4] + S[1][4][4] + S[1][10][1]);
            }

            // s-s-d elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 0 && _lmax_gamma >= 2) {
                S[1][1][8] = gmb0 * S[1][1][2] + fak * S[1][1][1];
                S[1][1][5] = gmb1 * S[1][1][2];
                S[1][1][6] = gmb2 * S[1][1][2];
                S[1][1][9] = gmb1 * S[1][1][3] + fak * S[1][1][1];
                S[1][1][7] = gmb2 * S[1][1][3];
                S[1][1][10] = gmb2 * S[1][1][4] + fak * S[1][1][1];
            }

            // s-p-d elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 1 && _lmax_gamma >= 2) {
                S[1][2][5] = gmc0 * S[1][1][5] + fak * S[1][1][3];
                S[1][3][5] = gmc1 * S[1][1][5] + fak * S[1][1][2];
                S[1][4][5] = gmc2 * S[1][1][5];
                S[1][2][6] = gmc0 * S[1][1][6] + fak * S[1][1][4];
                S[1][3][6] = gmc1 * S[1][1][6];
                S[1][4][6] = gmc2 * S[1][1][6] + fak * S[1][1][2];
                S[1][2][7] = gmc0 * S[1][1][7];
                S[1][3][7] = gmc1 * S[1][1][7] + fak * S[1][1][4];
                S[1][4][7] = gmc2 * S[1][1][7] + fak * S[1][1][3];
                S[1][2][8] = gmc0 * S[1][1][8] + fak2 * S[1][1][2];
                S[1][3][8] = gmc1 * S[1][1][8];
                S[1][4][8] = gmc2 * S[1][1][8];
                S[1][2][9] = gmc0 * S[1][1][9];
                S[1][3][9] = gmc1 * S[1][1][9] + fak2 * S[1][1][3];
                S[1][4][9] = gmc2 * S[1][1][9];
                S[1][2][10] = gmc0 * S[1][1][10];
                S[1][3][10] = gmc1 * S[1][1][10];
                S[1][4][10] = gmc2 * S[1][1][10] + fak2 * S[1][1][4];
            }

            // s-d-d elements
            if (_lmax_alpha >= 0 && _lmax_gw >= 2 && _lmax_gamma >= 2) {
                S[1][8][5] = gmc0 * S[1][2][5] + fak * (S[1][1][5] + S[1][2][3]);
                S[1][5][5] = gmc1 * S[1][2][5] + fak * S[1][2][2];
                S[1][6][5] = gmc2 * S[1][2][5];
                S[1][8][6] = gmc0 * S[1][2][6] + fak * (S[1][1][6] + S[1][2][4]);
                S[1][5][6] = gmc1 * S[1][2][6];
                S[1][6][6] = gmc2 * S[1][2][6] + fak * S[1][2][2];
                S[1][8][7] = gmc0 * S[1][2][7] + fak * S[1][1][7];
                S[1][5][7] = gmc1 * S[1][2][7] + fak * S[1][2][4];
                S[1][6][7] = gmc2 * S[1][2][7] + fak * S[1][2][3];
                S[1][8][8] = gmc0 * S[1][2][8] + fak * (S[1][1][8] + S[1][2][2] + S[1][2][2]);
                S[1][5][8] = gmc1 * S[1][2][8];
                S[1][6][8] = gmc2 * S[1][2][8];
                S[1][8][9] = gmc0 * S[1][2][9] + fak * S[1][1][9];
                S[1][5][9] = gmc1 * S[1][2][9] + fak2 * S[1][2][3];
                S[1][6][9] = gmc2 * S[1][2][9];
                S[1][8][10] = gmc0 * S[1][2][10] + fak * S[1][1][10];
                S[1][5][10] = gmc1 * S[1][2][10];
                S[1][6][10] = gmc2 * S[1][2][10] + fak2 * S[1][2][4];
                S[1][9][5] = gmc1 * S[1][3][5] + fak * (S[1][1][5] + S[1][3][2]);
                S[1][7][5] = gmc2 * S[1][3][5];
                S[1][9][6] = gmc1 * S[1][3][6] + fak * S[1][1][6];
                S[1][7][6] = gmc2 * S[1][3][6] + fak * S[1][3][2];
                S[1][9][7] = gmc1 * S[1][3][7] + fak * (S[1][1][7] + S[1][3][4]);
                S[1][7][7] = gmc2 * S[1][3][7] + fak * S[1][3][3];
                S[1][9][8] = gmc1 * S[1][3][8] + fak * S[1][1][8];
                S[1][7][8] = gmc2 * S[1][3][8];
                S[1][9][9] = gmc1 * S[1][3][9] + fak * (S[1][1][9] + S[1][3][3] + S[1][3][3]);
                S[1][7][9] = gmc2 * S[1][3][9];
                S[1][9][10] = gmc1 * S[1][3][10] + fak * S[1][1][10];
                S[1][7][10] = gmc2 * S[1][3][10] + fak2 * S[1][3][4];
                S[1][10][5] = gmc2 * S[1][4][5] + fak * S[1][1][5];
                S[1][10][6] = gmc2 * S[1][4][6] + fak * (S[1][1][6] + S[1][4][2]);
                S[1][10][7] = gmc2 * S[1][4][7] + fak * (S[1][1][7] + S[1][4][3]);
                S[1][10][8] = gmc2 * S[1][4][8] + fak * S[1][1][8];
                S[1][10][9] = gmc2 * S[1][4][9] + fak * S[1][1][9];
                S[1][10][10] = gmc2 * S[1][4][10] + fak * (S[1][1][10] + S[1][4][4] + S[1][4][4]);
            }

            // d-s-s elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 0 && _lmax_gamma >= 0) {
                S[8][1][1] = gma0 * S[2][1][1] + fak * S[1][1][1];
                S[5][1][1] = gma1 * S[2][1][1];
                S[6][1][1] = gma2 * S[2][1][1];
                S[9][1][1] = gma1 * S[3][1][1] + fak * S[1][1][1];
                S[7][1][1] = gma2 * S[3][1][1];
                S[10][1][1] = gma2 * S[4][1][1] + fak * S[1][1][1];
            }

            // d-p-s elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 1 && _lmax_gamma >= 0) {
                S[8][2][1] = gma0 * S[2][2][1] + fak * (S[1][2][1] + S[2][1][1]);
                S[5][2][1] = gma1 * S[2][2][1];
                S[6][2][1] = gma2 * S[2][2][1];
                S[8][3][1] = gma0 * S[2][3][1] + fak * S[1][3][1];
                S[5][3][1] = gma1 * S[2][3][1] + fak * S[2][1][1];
                S[6][3][1] = gma2 * S[2][3][1];
                S[8][4][1] = gma0 * S[2][4][1] + fak * S[1][4][1];
                S[5][4][1] = gma1 * S[2][4][1];
                S[6][4][1] = gma2 * S[2][4][1] + fak * S[2][1][1];
                S[9][2][1] = gma1 * S[3][2][1] + fak * S[1][2][1];
                S[7][2][1] = gma2 * S[3][2][1];
                S[9][3][1] = gma1 * S[3][3][1] + fak * (S[1][3][1] + S[3][1][1]);
                S[7][3][1] = gma2 * S[3][3][1];
                S[9][4][1] = gma1 * S[3][4][1] + fak * S[1][4][1];
                S[7][4][1] = gma2 * S[3][4][1] + fak * S[3][1][1];
                S[10][2][1] = gma2 * S[4][2][1] + fak * S[1][2][1];
                S[10][3][1] = gma2 * S[4][3][1] + fak * S[1][3][1];
                S[10][4][1] = gma2 * S[4][4][1] + fak * (S[1][4][1] + S[4][1][1]);
            }

            // d-d-s elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 2 && _lmax_gamma >= 0) {
                S[8][5][1] = gma0 * S[2][5][1] + fak * (S[1][5][1] + S[2][3][1]);
                S[5][5][1] = gma1 * S[2][5][1] + fak * S[2][2][1];
                S[6][5][1] = gma2 * S[2][5][1];
                S[8][6][1] = gma0 * S[2][6][1] + fak * (S[1][6][1] + S[2][4][1]);
                S[5][6][1] = gma1 * S[2][6][1];
                S[6][6][1] = gma2 * S[2][6][1] + fak * S[2][2][1];
                S[8][7][1] = gma0 * S[2][7][1] + fak * S[1][7][1];
                S[5][7][1] = gma1 * S[2][7][1] + fak * S[2][4][1];
                S[6][7][1] = gma2 * S[2][7][1] + fak * S[2][3][1];
                S[8][8][1] = gma0 * S[2][8][1] + fak * (S[1][8][1] + S[2][2][1] + S[2][2][1]);
                S[5][8][1] = gma1 * S[2][8][1];
                S[6][8][1] = gma2 * S[2][8][1];
                S[8][9][1] = gma0 * S[2][9][1] + fak * S[1][9][1];
                S[5][9][1] = gma1 * S[2][9][1] + fak2 * S[2][3][1];
                S[6][9][1] = gma2 * S[2][9][1];
                S[8][10][1] = gma0 * S[2][10][1] + fak * S[1][10][1];
                S[5][10][1] = gma1 * S[2][10][1];
                S[6][10][1] = gma2 * S[2][10][1] + fak2 * S[2][4][1];
                S[9][5][1] = gma1 * S[3][5][1] + fak * (S[1][5][1] + S[3][2][1]);
                S[7][5][1] = gma2 * S[3][5][1];
                S[9][6][1] = gma1 * S[3][6][1] + fak * S[1][6][1];
                S[7][6][1] = gma2 * S[3][6][1] + fak * S[3][2][1];
                S[9][7][1] = gma1 * S[3][7][1] + fak * (S[1][7][1] + S[3][4][1]);
                S[7][7][1] = gma2 * S[3][7][1] + fak * S[3][3][1];
                S[9][8][1] = gma1 * S[3][8][1] + fak * S[1][8][1];
                S[7][8][1] = gma2 * S[3][8][1];
                S[9][9][1] = gma1 * S[3][9][1] + fak * (S[1][9][1] + S[3][3][1] + S[3][3][1]);
                S[7][9][1] = gma2 * S[3][9][1];
                S[9][10][1] = gma1 * S[3][10][1] + fak * S[1][10][1];
                S[7][10][1] = gma2 * S[3][10][1] + fak2 * S[3][4][1];
                S[10][5][1] = gma2 * S[4][5][1] + fak * S[1][5][1];
                S[10][6][1] = gma2 * S[4][6][1] + fak * (S[1][6][1] + S[4][2][1]);
                S[10][7][1] = gma2 * S[4][7][1] + fak * (S[1][7][1] + S[4][3][1]);
                S[10][8][1] = gma2 * S[4][8][1] + fak * S[1][8][1];
                S[10][9][1] = gma2 * S[4][9][1] + fak * S[1][9][1];
                S[10][10][1] = gma2 * S[4][10][1] + fak * (S[1][10][1] + S[4][4][1] + S[4][4][1]);
            }

            // p-s-d elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 0 && _lmax_gamma >= 2) {
                S[2][1][5] = gma0 * S[1][1][5] + fak * S[1][1][3];
                S[3][1][5] = gma1 * S[1][1][5] + fak * S[1][1][2];
                S[4][1][5] = gma2 * S[1][1][5];
                S[2][1][6] = gma0 * S[1][1][6] + fak * S[1][1][4];
                S[3][1][6] = gma1 * S[1][1][6];
                S[4][1][6] = gma2 * S[1][1][6] + fak * S[1][1][2];
                S[2][1][7] = gma0 * S[1][1][7];
                S[3][1][7] = gma1 * S[1][1][7] + fak * S[1][1][4];
                S[4][1][7] = gma2 * S[1][1][7] + fak * S[1][1][3];
                S[2][1][8] = gma0 * S[1][1][8] + fak2 * S[1][1][2];
                S[3][1][8] = gma1 * S[1][1][8];
                S[4][1][8] = gma2 * S[1][1][8];
                S[2][1][9] = gma0 * S[1][1][9];
                S[3][1][9] = gma1 * S[1][1][9] + fak2 * S[1][1][3];
                S[4][1][9] = gma2 * S[1][1][9];
                S[2][1][10] = gma0 * S[1][1][10];
                S[3][1][10] = gma1 * S[1][1][10];
                S[4][1][10] = gma2 * S[1][1][10] + fak2 * S[1][1][4];
            }

            // p-p-d elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 1 && _lmax_gamma >= 2) {
                S[2][2][5] = gma0 * S[1][2][5] + fak * (S[1][1][5] + S[1][2][3]);
                S[3][2][5] = gma1 * S[1][2][5] + fak * S[1][2][2];
                S[4][2][5] = gma2 * S[1][2][5];
                S[2][2][6] = gma0 * S[1][2][6] + fak * (S[1][1][6] + S[1][2][4]);
                S[3][2][6] = gma1 * S[1][2][6];
                S[4][2][6] = gma2 * S[1][2][6] + fak * S[1][2][2];
                S[2][2][7] = gma0 * S[1][2][7] + fak * S[1][1][7];
                S[3][2][7] = gma1 * S[1][2][7] + fak * S[1][2][4];
                S[4][2][7] = gma2 * S[1][2][7] + fak * S[1][2][3];
                S[2][2][8] = gma0 * S[1][2][8] + fak * (S[1][1][8] + S[1][2][2] + S[1][2][2]);
                S[3][2][8] = gma1 * S[1][2][8];
                S[4][2][8] = gma2 * S[1][2][8];
                S[2][2][9] = gma0 * S[1][2][9] + fak * S[1][1][9];
                S[3][2][9] = gma1 * S[1][2][9] + fak2 * S[1][2][3];
                S[4][2][9] = gma2 * S[1][2][9];
                S[2][2][10] = gma0 * S[1][2][10] + fak * S[1][1][10];
                S[3][2][10] = gma1 * S[1][2][10];
                S[4][2][10] = gma2 * S[1][2][10] + fak2 * S[1][2][4];
                S[2][3][5] = gma0 * S[1][3][5] + fak * S[1][3][3];
                S[3][3][5] = gma1 * S[1][3][5] + fak * (S[1][1][5] + S[1][3][2]);
                S[4][3][5] = gma2 * S[1][3][5];
                S[2][3][6] = gma0 * S[1][3][6] + fak * S[1][3][4];
                S[3][3][6] = gma1 * S[1][3][6] + fak * S[1][1][6];
                S[4][3][6] = gma2 * S[1][3][6] + fak * S[1][3][2];
                S[2][3][7] = gma0 * S[1][3][7];
                S[3][3][7] = gma1 * S[1][3][7] + fak * (S[1][1][7] + S[1][3][4]);
                S[4][3][7] = gma2 * S[1][3][7] + fak * S[1][3][3];
                S[2][3][8] = gma0 * S[1][3][8] + fak2 * S[1][3][2];
                S[3][3][8] = gma1 * S[1][3][8] + fak * S[1][1][8];
                S[4][3][8] = gma2 * S[1][3][8];
                S[2][3][9] = gma0 * S[1][3][9];
                S[3][3][9] = gma1 * S[1][3][9] + fak * (S[1][1][9] + S[1][3][3] + S[1][3][3]);
                S[4][3][9] = gma2 * S[1][3][9];
                S[2][3][10] = gma0 * S[1][3][10];
                S[3][3][10] = gma1 * S[1][3][10] + fak * S[1][1][10];
                S[4][3][10] = gma2 * S[1][3][10] + fak2 * S[1][3][4];
                S[2][4][5] = gma0 * S[1][4][5] + fak * S[1][4][3];
                S[3][4][5] = gma1 * S[1][4][5] + fak * S[1][4][2];
                S[4][4][5] = gma2 * S[1][4][5] + fak * S[1][1][5];
                S[2][4][6] = gma0 * S[1][4][6] + fak * S[1][4][4];
                S[3][4][6] = gma1 * S[1][4][6];
                S[4][4][6] = gma2 * S[1][4][6] + fak * (S[1][1][6] + S[1][4][2]);
                S[2][4][7] = gma0 * S[1][4][7];
                S[3][4][7] = gma1 * S[1][4][7] + fak * S[1][4][4];
                S[4][4][7] = gma2 * S[1][4][7] + fak * (S[1][1][7] + S[1][4][3]);
                S[2][4][8] = gma0 * S[1][4][8] + fak2 * S[1][4][2];
                S[3][4][8] = gma1 * S[1][4][8];
                S[4][4][8] = gma2 * S[1][4][8] + fak * S[1][1][8];
                S[2][4][9] = gma0 * S[1][4][9];
                S[3][4][9] = gma1 * S[1][4][9] + fak2 * S[1][4][3];
                S[4][4][9] = gma2 * S[1][4][9] + fak * S[1][1][9];
                S[2][4][10] = gma0 * S[1][4][10];
                S[3][4][10] = gma1 * S[1][4][10];
                S[4][4][10] = gma2 * S[1][4][10] + fak * (S[1][1][10] + S[1][4][4] + S[1][4][4]);
            }

            // p-d-d elements
            if (_lmax_alpha >= 1 && _lmax_gw >= 2 && _lmax_gamma >= 2) {
                S[2][5][5] = gma0 * S[1][5][5] + fak * (S[1][3][5] + S[1][5][3]);
                S[3][5][5] = gma1 * S[1][5][5] + fak * (S[1][2][5] + S[1][5][2]);
                S[4][5][5] = gma2 * S[1][5][5];
                S[2][5][6] = gma0 * S[1][5][6] + fak * (S[1][3][6] + S[1][5][4]);
                S[3][5][6] = gma1 * S[1][5][6] + fak * S[1][2][6];
                S[4][5][6] = gma2 * S[1][5][6] + fak * S[1][5][2];
                S[2][5][7] = gma0 * S[1][5][7] + fak * S[1][3][7];
                S[3][5][7] = gma1 * S[1][5][7] + fak * (S[1][2][7] + S[1][5][4]);
                S[4][5][7] = gma2 * S[1][5][7] + fak * S[1][5][3];
                S[2][5][8] = gma0 * S[1][5][8] + fak * (S[1][3][8] + S[1][5][2] + S[1][5][2]);
                S[3][5][8] = gma1 * S[1][5][8] + fak * S[1][2][8];
                S[4][5][8] = gma2 * S[1][5][8];
                S[2][5][9] = gma0 * S[1][5][9] + fak * S[1][3][9];
                S[3][5][9] = gma1 * S[1][5][9] + fak * (S[1][2][9] + S[1][5][3] + S[1][5][3]);
                S[4][5][9] = gma2 * S[1][5][9];
                S[2][5][10] = gma0 * S[1][5][10] + fak * S[1][3][10];
                S[3][5][10] = gma1 * S[1][5][10] + fak * S[1][2][10];
                S[4][5][10] = gma2 * S[1][5][10] + fak2 * S[1][5][4];
                S[2][6][5] = gma0 * S[1][6][5] + fak * (S[1][4][5] + S[1][6][3]);
                S[3][6][5] = gma1 * S[1][6][5] + fak * S[1][6][2];
                S[4][6][5] = gma2 * S[1][6][5] + fak * S[1][2][5];
                S[2][6][6] = gma0 * S[1][6][6] + fak * (S[1][4][6] + S[1][6][4]);
                S[3][6][6] = gma1 * S[1][6][6];
                S[4][6][6] = gma2 * S[1][6][6] + fak * (S[1][2][6] + S[1][6][2]);
                S[2][6][7] = gma0 * S[1][6][7] + fak * S[1][4][7];
                S[3][6][7] = gma1 * S[1][6][7] + fak * S[1][6][4];
                S[4][6][7] = gma2 * S[1][6][7] + fak * (S[1][2][7] + S[1][6][3]);
                S[2][6][8] = gma0 * S[1][6][8] + fak * (S[1][4][8] + S[1][6][2] + S[1][6][2]);
                S[3][6][8] = gma1 * S[1][6][8];
                S[4][6][8] = gma2 * S[1][6][8] + fak * S[1][2][8];
                S[2][6][9] = gma0 * S[1][6][9] + fak * S[1][4][9];
                S[3][6][9] = gma1 * S[1][6][9] + fak2 * S[1][6][3];
                S[4][6][9] = gma2 * S[1][6][9] + fak * S[1][2][9];
                S[2][6][10] = gma0 * S[1][6][10] + fak * S[1][4][10];
                S[3][6][10] = gma1 * S[1][6][10];
                S[4][6][10] = gma2 * S[1][6][10] + fak * (S[1][2][10] + S[1][6][4] + S[1][6][4]);
                S[2][7][5] = gma0 * S[1][7][5] + fak * S[1][7][3];
                S[3][7][5] = gma1 * S[1][7][5] + fak * (S[1][4][5] + S[1][7][2]);
                S[4][7][5] = gma2 * S[1][7][5] + fak * S[1][3][5];
                S[2][7][6] = gma0 * S[1][7][6] + fak * S[1][7][4];
                S[3][7][6] = gma1 * S[1][7][6] + fak * S[1][4][6];
                S[4][7][6] = gma2 * S[1][7][6] + fak * (S[1][3][6] + S[1][7][2]);
                S[2][7][7] = gma0 * S[1][7][7];
                S[3][7][7] = gma1 * S[1][7][7] + fak * (S[1][4][7] + S[1][7][4]);
                S[4][7][7] = gma2 * S[1][7][7] + fak * (S[1][3][7] + S[1][7][3]);
                S[2][7][8] = gma0 * S[1][7][8] + fak * (2 * S[1][7][2]);
                S[3][7][8] = gma1 * S[1][7][8] + fak * S[1][4][8];
                S[4][7][8] = gma2 * S[1][7][8] + fak * S[1][3][8];
                S[2][7][9] = gma0 * S[1][7][9];
                S[3][7][9] = gma1 * S[1][7][9] + fak * (S[1][4][9] + S[1][7][3] + S[1][7][3]);
                S[4][7][9] = gma2 * S[1][7][9] + fak * S[1][3][9];
                S[2][7][10] = gma0 * S[1][7][10];
                S[3][7][10] = gma1 * S[1][7][10] + fak * S[1][4][10];
                S[4][7][10] = gma2 * S[1][7][10] + fak * (S[1][3][10] + S[1][7][4] + S[1][7][4]);
                S[2][8][5] = gma0 * S[1][8][5] + fak * (S[1][2][5] + S[1][2][5] + S[1][8][3]);
                S[3][8][5] = gma1 * S[1][8][5] + fak * S[1][8][2];
                S[4][8][5] = gma2 * S[1][8][5];
                S[2][8][6] = gma0 * S[1][8][6] + fak * (S[1][2][6] + S[1][2][6] + S[1][8][4]);
                S[3][8][6] = gma1 * S[1][8][6];
                S[4][8][6] = gma2 * S[1][8][6] + fak * S[1][8][2];
                S[2][8][7] = gma0 * S[1][8][7] + fak2 * S[1][2][7];
                S[3][8][7] = gma1 * S[1][8][7] + fak * S[1][8][4];
                S[4][8][7] = gma2 * S[1][8][7] + fak * S[1][8][3];
                S[2][8][8] = gma0 * S[1][8][8] + fak2 * (S[1][2][8] + S[1][8][2]);
                S[3][8][8] = gma1 * S[1][8][8];
                S[4][8][8] = gma2 * S[1][8][8];
                S[2][8][9] = gma0 * S[1][8][9] + fak2 * S[1][2][9];
                S[3][8][9] = gma1 * S[1][8][9] + fak2 * S[1][8][3];
                S[4][8][9] = gma2 * S[1][8][9];
                S[2][8][10] = gma0 * S[1][8][10] + fak2 * S[1][2][10];
                S[3][8][10] = gma1 * S[1][8][10];
                S[4][8][10] = gma2 * S[1][8][10] + fak2 * S[1][8][4];
                S[2][9][5] = gma0 * S[1][9][5] + fak * S[1][9][3];
                S[3][9][5] = gma1 * S[1][9][5] + fak * (S[1][3][5] + S[1][3][5] + S[1][9][2]);
                S[4][9][5] = gma2 * S[1][9][5];
                S[2][9][6] = gma0 * S[1][9][6] + fak * S[1][9][4];
                S[3][9][6] = gma1 * S[1][9][6] + fak2 * S[1][3][6];
                S[4][9][6] = gma2 * S[1][9][6] + fak * S[1][9][2];
                S[2][9][7] = gma0 * S[1][9][7];
                S[3][9][7] = gma1 * S[1][9][7] + fak * (S[1][3][7] + S[1][3][7] + S[1][9][4]);
                S[4][9][7] = gma2 * S[1][9][7] + fak * S[1][9][3];
                S[2][9][8] = gma0 * S[1][9][8] + fak2 * S[1][9][2];
                S[3][9][8] = gma1 * S[1][9][8] + fak2 * S[1][3][8];
                S[4][9][8] = gma2 * S[1][9][8];
                S[2][9][9] = gma0 * S[1][9][9];
                S[3][9][9] = gma1 * S[1][9][9] + fak2 * (S[1][3][9] + S[1][9][3]);
                S[4][9][9] = gma2 * S[1][9][9];
                S[2][9][10] = gma0 * S[1][9][10];
                S[3][9][10] = gma1 * S[1][9][10] + fak2 * S[1][3][10];
                S[4][9][10] = gma2 * S[1][9][10] + fak2 * S[1][9][4];
                S[2][10][5] = gma0 * S[1][10][5] + fak * S[1][10][3];
                S[3][10][5] = gma1 * S[1][10][5] + fak * S[1][10][2];
                S[4][10][5] = gma2 * S[1][10][5] + fak2 * S[1][4][5];
                S[2][10][6] = gma0 * S[1][10][6] + fak * S[1][10][4];
                S[3][10][6] = gma1 * S[1][10][6];
                S[4][10][6] = gma2 * S[1][10][6] + fak * (S[1][4][6] + S[1][4][6] + S[1][10][2]);
                S[2][10][7] = gma0 * S[1][10][7];
                S[3][10][7] = gma1 * S[1][10][7] + fak * S[1][10][4];
                S[4][10][7] = gma2 * S[1][10][7] + fak * (S[1][4][7] + S[1][4][7] + S[1][10][3]);
                S[2][10][8] = gma0 * S[1][10][8] + fak2 * S[1][10][2];
                S[3][10][8] = gma1 * S[1][10][8];
                S[4][10][8] = gma2 * S[1][10][8] + fak2 * S[1][4][8];
                S[2][10][9] = gma0 * S[1][10][9];
                S[3][10][9] = gma1 * S[1][10][9] + fak2 * S[1][10][3];
                S[4][10][9] = gma2 * S[1][10][9] + fak2 * S[1][4][9];
                S[2][10][10] = gma0 * S[1][10][10];
                S[3][10][10] = gma1 * S[1][10][10];
                S[4][10][10] = gma2 * S[1][10][10] + fak2 * (S[1][4][10] + S[1][10][4]);
            }

            // d-s-p elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 0 && _lmax_gamma >= 1) {
                S[8][1][2] = gma0 * S[2][1][2] + fak * (S[1][1][2] + S[2][1][1]);
                S[5][1][2] = gma1 * S[2][1][2];
                S[6][1][2] = gma2 * S[2][1][2];
                S[8][1][3] = gma0 * S[2][1][3] + fak * S[1][1][3];
                S[5][1][3] = gma1 * S[2][1][3] + fak * S[2][1][1];
                S[6][1][3] = gma2 * S[2][1][3];
                S[8][1][4] = gma0 * S[2][1][4] + fak * S[1][1][4];
                S[5][1][4] = gma1 * S[2][1][4];
                S[6][1][4] = gma2 * S[2][1][4] + fak * S[2][1][1];
                S[9][1][2] = gma1 * S[3][1][2] + fak * S[1][1][2];
                S[7][1][2] = gma2 * S[3][1][2];
                S[9][1][3] = gma1 * S[3][1][3] + fak * (S[1][1][3] + S[3][1][1]);
                S[7][1][3] = gma2 * S[3][1][3];
                S[9][1][4] = gma1 * S[3][1][4] + fak * S[1][1][4];
                S[7][1][4] = gma2 * S[3][1][4] + fak * S[3][1][1];
                S[10][1][2] = gma2 * S[4][1][2] + fak * S[1][1][2];
                S[10][1][3] = gma2 * S[4][1][3] + fak * S[1][1][3];
                S[10][1][4] = gma2 * S[4][1][4] + fak * (S[1][1][4] + S[4][1][1]);
            }

            // d-p-d elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 1 && _lmax_gamma >= 1) {
                S[8][2][2] = gma0 * S[2][2][2] + fak * (S[1][2][2] + S[2][1][2] + S[2][2][1]);
                S[5][2][2] = gma1 * S[2][2][2];
                S[6][2][2] = gma2 * S[2][2][2];
                S[8][2][3] = gma0 * S[2][2][3] + fak * (S[1][2][3] + S[2][1][3]);
                S[5][2][3] = gma1 * S[2][2][3] + fak * S[2][2][1];
                S[6][2][3] = gma2 * S[2][2][3];
                S[8][2][4] = gma0 * S[2][2][4] + fak * (S[1][2][4] + S[2][1][4]);
                S[5][2][4] = gma1 * S[2][2][4];
                S[6][2][4] = gma2 * S[2][2][4] + fak * S[2][2][1];
                S[8][3][2] = gma0 * S[2][3][2] + fak * (S[1][3][2] + S[2][3][1]);
                S[5][3][2] = gma1 * S[2][3][2] + fak * S[2][1][2];
                S[6][3][2] = gma2 * S[2][3][2];
                S[8][3][3] = gma0 * S[2][3][3] + fak * S[1][3][3];
                S[5][3][3] = gma1 * S[2][3][3] + fak * (S[2][1][3] + S[2][3][1]);
                S[6][3][3] = gma2 * S[2][3][3];
                S[8][3][4] = gma0 * S[2][3][4] + fak * S[1][3][4];
                S[5][3][4] = gma1 * S[2][3][4] + fak * S[2][1][4];
                S[6][3][4] = gma2 * S[2][3][4] + fak * S[2][3][1];
                S[8][4][2] = gma0 * S[2][4][2] + fak * (S[1][4][2] + S[2][4][1]);
                S[5][4][2] = gma1 * S[2][4][2];
                S[6][4][2] = gma2 * S[2][4][2] + fak * S[2][1][2];
                S[8][4][3] = gma0 * S[2][4][3] + fak * S[1][4][3];
                S[5][4][3] = gma1 * S[2][4][3] + fak * S[2][4][1];
                S[6][4][3] = gma2 * S[2][4][3] + fak * S[2][1][3];
                S[8][4][4] = gma0 * S[2][4][4] + fak * S[1][4][4];
                S[5][4][4] = gma1 * S[2][4][4];
                S[6][4][4] = gma2 * S[2][4][4] + fak * (S[2][1][4] + S[2][4][1]);
                S[9][2][2] = gma1 * S[3][2][2] + fak * S[1][2][2];
                S[7][2][2] = gma2 * S[3][2][2];
                S[9][2][3] = gma1 * S[3][2][3] + fak * (S[1][2][3] + S[3][2][1]);
                S[7][2][3] = gma2 * S[3][2][3];
                S[9][2][4] = gma1 * S[3][2][4] + fak * S[1][2][4];
                S[7][2][4] = gma2 * S[3][2][4] + fak * S[3][2][1];
                S[9][3][2] = gma1 * S[3][3][2] + fak * (S[1][3][2] + S[3][1][2]);
                S[7][3][2] = gma2 * S[3][3][2];
                S[9][3][3] = gma1 * S[3][3][3] + fak * (S[1][3][3] + S[3][1][3] + S[3][3][1]);
                S[7][3][3] = gma2 * S[3][3][3];
                S[9][3][4] = gma1 * S[3][3][4] + fak * (S[1][3][4] + S[3][1][4]);
                S[7][3][4] = gma2 * S[3][3][4] + fak * S[3][3][1];
                S[9][4][2] = gma1 * S[3][4][2] + fak * S[1][4][2];
                S[7][4][2] = gma2 * S[3][4][2] + fak * S[3][1][2];
                S[9][4][3] = gma1 * S[3][4][3] + fak * (S[1][4][3] + S[3][4][1]);
                S[7][4][3] = gma2 * S[3][4][3] + fak * S[3][1][3];
                S[9][4][4] = gma1 * S[3][4][4] + fak * S[1][4][4];
                S[7][4][4] = gma2 * S[3][4][4] + fak * (S[3][1][4] + S[3][4][1]);
                S[10][2][2] = gma2 * S[4][2][2] + fak * S[1][2][2];
                S[10][2][3] = gma2 * S[4][2][3] + fak * S[1][2][3];
                S[10][2][4] = gma2 * S[4][2][4] + fak * (S[1][2][4] + S[4][2][1]);
                S[10][3][2] = gma2 * S[4][3][2] + fak * S[1][3][2];
                S[10][3][3] = gma2 * S[4][3][3] + fak * S[1][3][3];
                S[10][3][4] = gma2 * S[4][3][4] + fak * (S[1][3][4] + S[4][3][1]);
                S[10][4][2] = gma2 * S[4][4][2] + fak * (S[1][4][2] + S[4][1][2]);
                S[10][4][3] = gma2 * S[4][4][3] + fak * (S[1][4][3] + S[4][1][3]);
                S[10][4][4] = gma2 * S[4][4][4] + fak * (S[1][4][4] + S[4][1][4] + S[4][4][1]);
            }

            // d-d-p elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 2 && _lmax_gamma >= 1) {
                S[8][5][2] = gma0 * S[2][5][2] + fak * (S[1][5][2] + S[2][3][2] + S[2][5][1]);
                S[5][5][2] = gma1 * S[2][5][2] + fak * S[2][2][2];
                S[6][5][2] = gma2 * S[2][5][2];
                S[8][5][3] = gma0 * S[2][5][3] + fak * (S[1][5][3] + S[2][3][3]);
                S[5][5][3] = gma1 * S[2][5][3] + fak * (S[2][2][3] + S[2][5][1]);
                S[6][5][3] = gma2 * S[2][5][3];
                S[8][5][4] = gma0 * S[2][5][4] + fak * (S[1][5][4] + S[2][3][4]);
                S[5][5][4] = gma1 * S[2][5][4] + fak * S[2][2][4];
                S[6][5][4] = gma2 * S[2][5][4] + fak * S[2][5][1];
                S[8][6][2] = gma0 * S[2][6][2] + fak * (S[1][6][2] + S[2][4][2] + S[2][6][1]);
                S[5][6][2] = gma1 * S[2][6][2];
                S[6][6][2] = gma2 * S[2][6][2] + fak * S[2][2][2];
                S[8][6][3] = gma0 * S[2][6][3] + fak * (S[1][6][3] + S[2][4][3]);
                S[5][6][3] = gma1 * S[2][6][3] + fak * S[2][6][1];
                S[6][6][3] = gma2 * S[2][6][3] + fak * S[2][2][3];
                S[8][6][4] = gma0 * S[2][6][4] + fak * (S[1][6][4] + S[2][4][4]);
                S[5][6][4] = gma1 * S[2][6][4];
                S[6][6][4] = gma2 * S[2][6][4] + fak * (S[2][2][4] + S[2][6][1]);
                S[8][7][2] = gma0 * S[2][7][2] + fak * (S[1][7][2] + S[2][7][1]);
                S[5][7][2] = gma1 * S[2][7][2] + fak * S[2][4][2];
                S[6][7][2] = gma2 * S[2][7][2] + fak * S[2][3][2];
                S[8][7][3] = gma0 * S[2][7][3] + fak * S[1][7][3];
                S[5][7][3] = gma1 * S[2][7][3] + fak * (S[2][4][3] + S[2][7][1]);
                S[6][7][3] = gma2 * S[2][7][3] + fak * S[2][3][3];
                S[8][7][4] = gma0 * S[2][7][4] + fak * S[1][7][4];
                S[5][7][4] = gma1 * S[2][7][4] + fak * S[2][4][4];
                S[6][7][4] = gma2 * S[2][7][4] + fak * (S[2][3][4] + S[2][7][1]);
                S[8][8][2] = gma0 * S[2][8][2] + fak * (S[1][8][2] + S[2][2][2] + S[2][2][2] + S[2][8][1]);
                S[5][8][2] = gma1 * S[2][8][2];
                S[6][8][2] = gma2 * S[2][8][2];
                S[8][8][3] = gma0 * S[2][8][3] + fak * (S[1][8][3] + S[2][2][3] + S[2][2][3]);
                S[5][8][3] = gma1 * S[2][8][3] + fak * S[2][8][1];
                S[6][8][3] = gma2 * S[2][8][3];
                S[8][8][4] = gma0 * S[2][8][4] + fak * (S[1][8][4] + S[2][2][4] + S[2][2][4]);
                S[5][8][4] = gma1 * S[2][8][4];
                S[6][8][4] = gma2 * S[2][8][4] + fak * S[2][8][1];
                S[8][9][2] = gma0 * S[2][9][2] + fak * (S[1][9][2] + S[2][9][1]);
                S[5][9][2] = gma1 * S[2][9][2] + fak2 * S[2][3][2];
                S[6][9][2] = gma2 * S[2][9][2];
                S[8][9][3] = gma0 * S[2][9][3] + fak * S[1][9][3];
                S[5][9][3] = gma1 * S[2][9][3] + fak * (S[2][3][3] + S[2][3][3] + S[2][9][1]);
                S[6][9][3] = gma2 * S[2][9][3];
                S[8][9][4] = gma0 * S[2][9][4] + fak * S[1][9][4];
                S[5][9][4] = gma1 * S[2][9][4] + fak * (2 * S[2][3][4]);
                S[6][9][4] = gma2 * S[2][9][4] + fak * S[2][9][1];
                S[8][10][2] = gma0 * S[2][10][2] + fak * (S[1][10][2] + S[2][10][1]);
                S[5][10][2] = gma1 * S[2][10][2];
                S[6][10][2] = gma2 * S[2][10][2] + fak * (2 * S[2][4][2]);
                S[8][10][3] = gma0 * S[2][10][3] + fak * S[1][10][3];
                S[5][10][3] = gma1 * S[2][10][3] + fak * S[2][10][1];
                S[6][10][3] = gma2 * S[2][10][3] + fak * (2 * S[2][4][3]);
                S[8][10][4] = gma0 * S[2][10][4] + fak * S[1][10][4];
                S[5][10][4] = gma1 * S[2][10][4];
                S[6][10][4] = gma2 * S[2][10][4] + fak * (S[2][4][4] + S[2][4][4] + S[2][10][1]);
                S[9][5][2] = gma1 * S[3][5][2] + fak * (S[1][5][2] + S[3][2][2]);
                S[7][5][2] = gma2 * S[3][5][2];
                S[9][5][3] = gma1 * S[3][5][3] + fak * (S[1][5][3] + S[3][2][3] + S[3][5][1]);
                S[7][5][3] = gma2 * S[3][5][3];
                S[9][5][4] = gma1 * S[3][5][4] + fak * (S[1][5][4] + S[3][2][4]);
                S[7][5][4] = gma2 * S[3][5][4] + fak * S[3][5][1];
                S[9][6][2] = gma1 * S[3][6][2] + fak * S[1][6][2];
                S[7][6][2] = gma2 * S[3][6][2] + fak * S[3][2][2];
                S[9][6][3] = gma1 * S[3][6][3] + fak * (S[1][6][3] + S[3][6][1]);
                S[7][6][3] = gma2 * S[3][6][3] + fak * S[3][2][3];
                S[9][6][4] = gma1 * S[3][6][4] + fak * S[1][6][4];
                S[7][6][4] = gma2 * S[3][6][4] + fak * (S[3][2][4] + S[3][6][1]);
                S[9][7][2] = gma1 * S[3][7][2] + fak * (S[1][7][2] + S[3][4][2]);
                S[7][7][2] = gma2 * S[3][7][2] + fak * S[3][3][2];
                S[9][7][3] = gma1 * S[3][7][3] + fak * (S[1][7][3] + S[3][4][3] + S[3][7][1]);
                S[7][7][3] = gma2 * S[3][7][3] + fak * S[3][3][3];
                S[9][7][4] = gma1 * S[3][7][4] + fak * (S[1][7][4] + S[3][4][4]);
                S[7][7][4] = gma2 * S[3][7][4] + fak * (S[3][3][4] + S[3][7][1]);
                S[9][8][2] = gma1 * S[3][8][2] + fak * S[1][8][2];
                S[7][8][2] = gma2 * S[3][8][2];
                S[9][8][3] = gma1 * S[3][8][3] + fak * (S[1][8][3] + S[3][8][1]);
                S[7][8][3] = gma2 * S[3][8][3];
                S[9][8][4] = gma1 * S[3][8][4] + fak * S[1][8][4];
                S[7][8][4] = gma2 * S[3][8][4] + fak * S[3][8][1];
                S[9][9][2] = gma1 * S[3][9][2] + fak * (S[1][9][2] + S[3][3][2] + S[3][3][2]);
                S[7][9][2] = gma2 * S[3][9][2];
                S[9][9][3] = gma1 * S[3][9][3] + fak * (S[1][9][3] + S[3][3][3] + S[3][3][3] + S[3][9][1]);
                S[7][9][3] = gma2 * S[3][9][3];
                S[9][9][4] = gma1 * S[3][9][4] + fak * (S[1][9][4] + S[3][3][4] + S[3][3][4]);
                S[7][9][4] = gma2 * S[3][9][4] + fak * S[3][9][1];
                S[9][10][2] = gma1 * S[3][10][2] + fak * S[1][10][2];
                S[7][10][2] = gma2 * S[3][10][2] + fak2 * S[3][4][2];
                S[9][10][3] = gma1 * S[3][10][3] + fak * (S[1][10][3] + S[3][10][1]);
                S[7][10][3] = gma2 * S[3][10][3] + fak2 * S[3][4][3];
                S[9][10][4] = gma1 * S[3][10][4] + fak * S[1][10][4];
                S[7][10][4] = gma2 * S[3][10][4] + fak * (S[3][4][4] + S[3][4][4] + S[3][10][1]);
                S[10][5][2] = gma2 * S[4][5][2] + fak * S[1][5][2];
                S[10][5][3] = gma2 * S[4][5][3] + fak * S[1][5][3];
                S[10][5][4] = gma2 * S[4][5][4] + fak * (S[1][5][4] + S[4][5][1]);
                S[10][6][2] = gma2 * S[4][6][2] + fak * (S[1][6][2] + S[4][2][2]);
                S[10][6][3] = gma2 * S[4][6][3] + fak * (S[1][6][3] + S[4][2][3]);
                S[10][6][4] = gma2 * S[4][6][4] + fak * (S[1][6][4] + S[4][2][4] + S[4][6][1]);
                S[10][7][2] = gma2 * S[4][7][2] + fak * (S[1][7][2] + S[4][3][2]);
                S[10][7][3] = gma2 * S[4][7][3] + fak * (S[1][7][3] + S[4][3][3]);
                S[10][7][4] = gma2 * S[4][7][4] + fak * (S[1][7][4] + S[4][3][4] + S[4][7][1]);
                S[10][8][2] = gma2 * S[4][8][2] + fak * S[1][8][2];
                S[10][8][3] = gma2 * S[4][8][3] + fak * S[1][8][3];
                S[10][8][4] = gma2 * S[4][8][4] + fak * (S[1][8][4] + S[4][8][1]);
                S[10][9][2] = gma2 * S[4][9][2] + fak * S[1][9][2];
                S[10][9][3] = gma2 * S[4][9][3] + fak * S[1][9][3];
                S[10][9][4] = gma2 * S[4][9][4] + fak * (S[1][9][4] + S[4][9][1]);
                S[10][10][2] = gma2 * S[4][10][2] + fak * (S[1][10][2] + S[4][4][2] + S[4][4][2]);
                S[10][10][3] = gma2 * S[4][10][3] + fak * (S[1][10][3] + S[4][4][3] + S[4][4][3]);
                S[10][10][4] = gma2 * S[4][10][4] + fak * (S[1][10][4] + S[4][4][4] + S[4][4][4] + S[4][10][1]);
            }

            // d-s-d elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 0 && _lmax_gamma >= 2) {
                S[8][1][5] = gma0 * S[2][1][5] + fak * (S[1][1][5] + S[2][1][3]);
                S[5][1][5] = gma1 * S[2][1][5] + fak * S[2][1][2];
                S[6][1][5] = gma2 * S[2][1][5];
                S[8][1][6] = gma0 * S[2][1][6] + fak * (S[1][1][6] + S[2][1][4]);
                S[5][1][6] = gma1 * S[2][1][6];
                S[6][1][6] = gma2 * S[2][1][6] + fak * S[2][1][2];
                S[8][1][7] = gma0 * S[2][1][7] + fak * S[1][1][7];
                S[5][1][7] = gma1 * S[2][1][7] + fak * S[2][1][4];
                S[6][1][7] = gma2 * S[2][1][7] + fak * S[2][1][3];
                S[8][1][8] = gma0 * S[2][1][8] + fak * (S[1][1][8] + S[2][1][2] + S[2][1][2]);
                S[5][1][8] = gma1 * S[2][1][8];
                S[6][1][8] = gma2 * S[2][1][8];
                S[8][1][9] = gma0 * S[2][1][9] + fak * S[1][1][9];
                S[5][1][9] = gma1 * S[2][1][9] + fak2 * S[2][1][3];
                S[6][1][9] = gma2 * S[2][1][9];
                S[8][1][10] = gma0 * S[2][1][10] + fak * S[1][1][10];
                S[5][1][10] = gma1 * S[2][1][10];
                S[6][1][10] = gma2 * S[2][1][10] + fak2 * S[2][1][4];
                S[9][1][5] = gma1 * S[3][1][5] + fak * (S[1][1][5] + S[3][1][2]);
                S[7][1][5] = gma2 * S[3][1][5];
                S[9][1][6] = gma1 * S[3][1][6] + fak * S[1][1][6];
                S[7][1][6] = gma2 * S[3][1][6] + fak * S[3][1][2];
                S[9][1][7] = gma1 * S[3][1][7] + fak * (S[1][1][7] + S[3][1][4]);
                S[7][1][7] = gma2 * S[3][1][7] + fak * S[3][1][3];
                S[9][1][8] = gma1 * S[3][1][8] + fak * S[1][1][8];
                S[7][1][8] = gma2 * S[3][1][8];
                S[9][1][9] = gma1 * S[3][1][9] + fak * (S[1][1][9] + S[3][1][3] + S[3][1][3]);
                S[7][1][9] = gma2 * S[3][1][9];
                S[9][1][10] = gma1 * S[3][1][10] + fak * S[1][1][10];
                S[7][1][10] = gma2 * S[3][1][10] + fak2 * S[3][1][4];
                S[10][1][5] = gma2 * S[4][1][5] + fak * S[1][1][5];
                S[10][1][6] = gma2 * S[4][1][6] + fak * (S[1][1][6] + S[4][1][2]);
                S[10][1][7] = gma2 * S[4][1][7] + fak * (S[1][1][7] + S[4][1][3]);
                S[10][1][8] = gma2 * S[4][1][8] + fak * S[1][1][8];
                S[10][1][9] = gma2 * S[4][1][9] + fak * S[1][1][9];
                S[10][1][10] = gma2 * S[4][1][10] + fak * (S[1][1][10] + S[4][1][4] + S[4][1][4]);
            }

            // d-p-d elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 1 && _lmax_gamma >= 2) {
                S[8][2][5] = gma0 * S[2][2][5] + fak * (S[1][2][5] + S[2][1][5] + S[2][2][3]);
                S[5][2][5] = gma1 * S[2][2][5] + fak * S[2][2][2];
                S[6][2][5] = gma2 * S[2][2][5];
                S[8][2][6] = gma0 * S[2][2][6] + fak * (S[1][2][6] + S[2][1][6] + S[2][2][4]);
                S[5][2][6] = gma1 * S[2][2][6];
                S[6][2][6] = gma2 * S[2][2][6] + fak * S[2][2][2];
                S[8][2][7] = gma0 * S[2][2][7] + fak * (S[1][2][7] + S[2][1][7]);
                S[5][2][7] = gma1 * S[2][2][7] + fak * S[2][2][4];
                S[6][2][7] = gma2 * S[2][2][7] + fak * S[2][2][3];
                S[8][2][8] = gma0 * S[2][2][8] + fak * (S[1][2][8] + S[2][1][8] + S[2][2][2] + S[2][2][2]);
                S[5][2][8] = gma1 * S[2][2][8];
                S[6][2][8] = gma2 * S[2][2][8];
                S[8][2][9] = gma0 * S[2][2][9] + fak * (S[1][2][9] + S[2][1][9]);
                S[5][2][9] = gma1 * S[2][2][9] + fak2 * S[2][2][3];
                S[6][2][9] = gma2 * S[2][2][9];
                S[8][2][10] = gma0 * S[2][2][10] + fak * (S[1][2][10] + S[2][1][10]);
                S[5][2][10] = gma1 * S[2][2][10];
                S[6][2][10] = gma2 * S[2][2][10] + fak2 * S[2][2][4];
                S[8][3][5] = gma0 * S[2][3][5] + fak * (S[1][3][5] + S[2][3][3]);
                S[5][3][5] = gma1 * S[2][3][5] + fak * (S[2][1][5] + S[2][3][2]);
                S[6][3][5] = gma2 * S[2][3][5];
                S[8][3][6] = gma0 * S[2][3][6] + fak * (S[1][3][6] + S[2][3][4]);
                S[5][3][6] = gma1 * S[2][3][6] + fak * S[2][1][6];
                S[6][3][6] = gma2 * S[2][3][6] + fak * S[2][3][2];
                S[8][3][7] = gma0 * S[2][3][7] + fak * S[1][3][7];
                S[5][3][7] = gma1 * S[2][3][7] + fak * (S[2][1][7] + S[2][3][4]);
                S[6][3][7] = gma2 * S[2][3][7] + fak * S[2][3][3];
                S[8][3][8] = gma0 * S[2][3][8] + fak * (S[1][3][8] + S[2][3][2] + S[2][3][2]);
                S[5][3][8] = gma1 * S[2][3][8] + fak * S[2][1][8];
                S[6][3][8] = gma2 * S[2][3][8];
                S[8][3][9] = gma0 * S[2][3][9] + fak * S[1][3][9];
                S[5][3][9] = gma1 * S[2][3][9] + fak * (S[2][1][9] + S[2][3][3] + S[2][3][3]);
                S[6][3][9] = gma2 * S[2][3][9];
                S[8][3][10] = gma0 * S[2][3][10] + fak * S[1][3][10];
                S[5][3][10] = gma1 * S[2][3][10] + fak * S[2][1][10];
                S[6][3][10] = gma2 * S[2][3][10] + fak2 * S[2][3][4];
                S[8][4][5] = gma0 * S[2][4][5] + fak * (S[1][4][5] + S[2][4][3]);
                S[5][4][5] = gma1 * S[2][4][5] + fak * S[2][4][2];
                S[6][4][5] = gma2 * S[2][4][5] + fak * S[2][1][5];
                S[8][4][6] = gma0 * S[2][4][6] + fak * (S[1][4][6] + S[2][4][4]);
                S[5][4][6] = gma1 * S[2][4][6];
                S[6][4][6] = gma2 * S[2][4][6] + fak * (S[2][1][6] + S[2][4][2]);
                S[8][4][7] = gma0 * S[2][4][7] + fak * S[1][4][7];
                S[5][4][7] = gma1 * S[2][4][7] + fak * S[2][4][4];
                S[6][4][7] = gma2 * S[2][4][7] + fak * (S[2][1][7] + S[2][4][3]);
                S[8][4][8] = gma0 * S[2][4][8] + fak * (S[1][4][8] + S[2][4][2] + S[2][4][2]);
                S[5][4][8] = gma1 * S[2][4][8];
                S[6][4][8] = gma2 * S[2][4][8] + fak * S[2][1][8];
                S[8][4][9] = gma0 * S[2][4][9] + fak * S[1][4][9];
                S[5][4][9] = gma1 * S[2][4][9] + fak2 * S[2][4][3];
                S[6][4][9] = gma2 * S[2][4][9] + fak * S[2][1][9];
                S[8][4][10] = gma0 * S[2][4][10] + fak * S[1][4][10];
                S[5][4][10] = gma1 * S[2][4][10];
                S[6][4][10] = gma2 * S[2][4][10] + fak * (S[2][1][10] + S[2][4][4] + S[2][4][4]);
                S[9][2][5] = gma1 * S[3][2][5] + fak * (S[1][2][5] + S[3][2][2]);
                S[7][2][5] = gma2 * S[3][2][5];
                S[9][2][6] = gma1 * S[3][2][6] + fak * S[1][2][6];
                S[7][2][6] = gma2 * S[3][2][6] + fak * S[3][2][2];
                S[9][2][7] = gma1 * S[3][2][7] + fak * (S[1][2][7] + S[3][2][4]);
                S[7][2][7] = gma2 * S[3][2][7] + fak * S[3][2][3];
                S[9][2][8] = gma1 * S[3][2][8] + fak * S[1][2][8];
                S[7][2][8] = gma2 * S[3][2][8];
                S[9][2][9] = gma1 * S[3][2][9] + fak * (S[1][2][9] + S[3][2][3] + S[3][2][3]);
                S[7][2][9] = gma2 * S[3][2][9];
                S[9][2][10] = gma1 * S[3][2][10] + fak * S[1][2][10];
                S[7][2][10] = gma2 * S[3][2][10] + fak2 * S[3][2][4];
                S[9][3][5] = gma1 * S[3][3][5] + fak * (S[1][3][5] + S[3][1][5] + S[3][3][2]);
                S[7][3][5] = gma2 * S[3][3][5];
                S[9][3][6] = gma1 * S[3][3][6] + fak * (S[1][3][6] + S[3][1][6]);
                S[7][3][6] = gma2 * S[3][3][6] + fak * S[3][3][2];
                S[9][3][7] = gma1 * S[3][3][7] + fak * (S[1][3][7] + S[3][1][7] + S[3][3][4]);
                S[7][3][7] = gma2 * S[3][3][7] + fak * S[3][3][3];
                S[9][3][8] = gma1 * S[3][3][8] + fak * (S[1][3][8] + S[3][1][8]);
                S[7][3][8] = gma2 * S[3][3][8];
                S[9][3][9] = gma1 * S[3][3][9] + fak * (S[1][3][9] + S[3][1][9] + S[3][3][3] + S[3][3][3]);
                S[7][3][9] = gma2 * S[3][3][9];
                S[9][3][10] = gma1 * S[3][3][10] + fak * (S[1][3][10] + S[3][1][10]);
                S[7][3][10] = gma2 * S[3][3][10] + fak2 * S[3][3][4];
                S[9][4][5] = gma1 * S[3][4][5] + fak * (S[1][4][5] + S[3][4][2]);
                S[7][4][5] = gma2 * S[3][4][5] + fak * S[3][1][5];
                S[9][4][6] = gma1 * S[3][4][6] + fak * S[1][4][6];
                S[7][4][6] = gma2 * S[3][4][6] + fak * (S[3][1][6] + S[3][4][2]);
                S[9][4][7] = gma1 * S[3][4][7] + fak * (S[1][4][7] + S[3][4][4]);
                S[7][4][7] = gma2 * S[3][4][7] + fak * (S[3][1][7] + S[3][4][3]);
                S[9][4][8] = gma1 * S[3][4][8] + fak * S[1][4][8];
                S[7][4][8] = gma2 * S[3][4][8] + fak * S[3][1][8];
                S[9][4][9] = gma1 * S[3][4][9] + fak * (S[1][4][9] + S[3][4][3] + S[3][4][3]);
                S[7][4][9] = gma2 * S[3][4][9] + fak * S[3][1][9];
                S[9][4][10] = gma1 * S[3][4][10] + fak * S[1][4][10];
                S[7][4][10] = gma2 * S[3][4][10] + fak * (S[3][1][10] + S[3][4][4] + S[3][4][4]);
                S[10][2][5] = gma2 * S[4][2][5] + fak * S[1][2][5];
                S[10][2][6] = gma2 * S[4][2][6] + fak * (S[1][2][6] + S[4][2][2]);
                S[10][2][7] = gma2 * S[4][2][7] + fak * (S[1][2][7] + S[4][2][3]);
                S[10][2][8] = gma2 * S[4][2][8] + fak * S[1][2][8];
                S[10][2][9] = gma2 * S[4][2][9] + fak * S[1][2][9];
                S[10][2][10] = gma2 * S[4][2][10] + fak * (S[1][2][10] + S[4][2][4] + S[4][2][4]);
                S[10][3][5] = gma2 * S[4][3][5] + fak * S[1][3][5];
                S[10][3][6] = gma2 * S[4][3][6] + fak * (S[1][3][6] + S[4][3][2]);
                S[10][3][7] = gma2 * S[4][3][7] + fak * (S[1][3][7] + S[4][3][3]);
                S[10][3][8] = gma2 * S[4][3][8] + fak * S[1][3][8];
                S[10][3][9] = gma2 * S[4][3][9] + fak * S[1][3][9];
                S[10][3][10] = gma2 * S[4][3][10] + fak * (S[1][3][10] + S[4][3][4] + S[4][3][4]);
                S[10][4][5] = gma2 * S[4][4][5] + fak * (S[1][4][5] + S[4][1][5]);
                S[10][4][6] = gma2 * S[4][4][6] + fak * (S[1][4][6] + S[4][1][6] + S[4][4][2]);
                S[10][4][7] = gma2 * S[4][4][7] + fak * (S[1][4][7] + S[4][1][7] + S[4][4][3]);
                S[10][4][8] = gma2 * S[4][4][8] + fak * (S[1][4][8] + S[4][1][8]);
                S[10][4][9] = gma2 * S[4][4][9] + fak * (S[1][4][9] + S[4][1][9]);
                S[10][4][10] = gma2 * S[4][4][10] + fak * (S[1][4][10] + S[4][1][10] + S[4][4][4] + S[4][4][4]);
            }

            // d-d-d elements
            if (_lmax_alpha >= 2 && _lmax_gw >= 2 && _lmax_gamma >= 2) {
                S[8][5][5] = gma0 * S[2][5][5] + fak * (S[1][5][5] + S[2][3][5] + S[2][5][3]);
                S[5][5][5] = gma1 * S[2][5][5] + fak * (S[2][2][5] + S[2][5][2]);
                S[6][5][5] = gma2 * S[2][5][5];
                S[8][5][6] = gma0 * S[2][5][6] + fak * (S[1][5][6] + S[2][3][6] + S[2][5][4]);
                S[5][5][6] = gma1 * S[2][5][6] + fak * S[2][2][6];
                S[6][5][6] = gma2 * S[2][5][6] + fak * S[2][5][2];
                S[8][5][7] = gma0 * S[2][5][7] + fak * (S[1][5][7] + S[2][3][7]);
                S[5][5][7] = gma1 * S[2][5][7] + fak * (S[2][2][7] + S[2][5][4]);
                S[6][5][7] = gma2 * S[2][5][7] + fak * S[2][5][3];
                S[8][5][8] = gma0 * S[2][5][8] + fak * (S[1][5][8] + S[2][3][8] + S[2][5][2] + S[2][5][2]);
                S[5][5][8] = gma1 * S[2][5][8] + fak * S[2][2][8];
                S[6][5][8] = gma2 * S[2][5][8];
                S[8][5][9] = gma0 * S[2][5][9] + fak * (S[1][5][9] + S[2][3][9]);
                S[5][5][9] = gma1 * S[2][5][9] + fak * (S[2][2][9] + S[2][5][3] + S[2][5][3]);
                S[6][5][9] = gma2 * S[2][5][9];
                S[8][5][10] = gma0 * S[2][5][10] + fak * (S[1][5][10] + S[2][3][10]);
                S[5][5][10] = gma1 * S[2][5][10] + fak * S[2][2][10];
                S[6][5][10] = gma2 * S[2][5][10] + fak2 * S[2][5][4];
                S[8][6][5] = gma0 * S[2][6][5] + fak * (S[1][6][5] + S[2][4][5] + S[2][6][3]);
                S[5][6][5] = gma1 * S[2][6][5] + fak * S[2][6][2];
                S[6][6][5] = gma2 * S[2][6][5] + fak * S[2][2][5];
                S[8][6][6] = gma0 * S[2][6][6] + fak * (S[1][6][6] + S[2][4][6] + S[2][6][4]);
                S[5][6][6] = gma1 * S[2][6][6];
                S[6][6][6] = gma2 * S[2][6][6] + fak * (S[2][2][6] + S[2][6][2]);
                S[8][6][7] = gma0 * S[2][6][7] + fak * (S[1][6][7] + S[2][4][7]);
                S[5][6][7] = gma1 * S[2][6][7] + fak * S[2][6][4];
                S[6][6][7] = gma2 * S[2][6][7] + fak * (S[2][2][7] + S[2][6][3]);
                S[8][6][8] = gma0 * S[2][6][8] + fak * (S[1][6][8] + S[2][4][8] + S[2][6][2] + S[2][6][2]);
                S[5][6][8] = gma1 * S[2][6][8];
                S[6][6][8] = gma2 * S[2][6][8] + fak * S[2][2][8];
                S[8][6][9] = gma0 * S[2][6][9] + fak * (S[1][6][9] + S[2][4][9]);
                S[5][6][9] = gma1 * S[2][6][9] + fak2 * S[2][6][3];
                S[6][6][9] = gma2 * S[2][6][9] + fak * S[2][2][9];
                S[8][6][10] = gma0 * S[2][6][10] + fak * (S[1][6][10] + S[2][4][10]);
                S[5][6][10] = gma1 * S[2][6][10];
                S[6][6][10] = gma2 * S[2][6][10] + fak * (S[2][2][10] + S[2][6][4] + S[2][6][4]);
                S[8][7][5] = gma0 * S[2][7][5] + fak * (S[1][7][5] + S[2][7][3]);
                S[5][7][5] = gma1 * S[2][7][5] + fak * (S[2][4][5] + S[2][7][2]);
                S[6][7][5] = gma2 * S[2][7][5] + fak * S[2][3][5];
                S[8][7][6] = gma0 * S[2][7][6] + fak * (S[1][7][6] + S[2][7][4]);
                S[5][7][6] = gma1 * S[2][7][6] + fak * S[2][4][6];
                S[6][7][6] = gma2 * S[2][7][6] + fak * (S[2][3][6] + S[2][7][2]);
                S[8][7][7] = gma0 * S[2][7][7] + fak * S[1][7][7];
                S[5][7][7] = gma1 * S[2][7][7] + fak * (S[2][4][7] + S[2][7][4]);
                S[6][7][7] = gma2 * S[2][7][7] + fak * (S[2][3][7] + S[2][7][3]);
                S[8][7][8] = gma0 * S[2][7][8] + fak * (S[1][7][8] + S[2][7][2] + S[2][7][2]);
                S[5][7][8] = gma1 * S[2][7][8] + fak * S[2][4][8];
                S[6][7][8] = gma2 * S[2][7][8] + fak * S[2][3][8];
                S[8][7][9] = gma0 * S[2][7][9] + fak * S[1][7][9];
                S[5][7][9] = gma1 * S[2][7][9] + fak * (S[2][4][9] + S[2][7][3] + S[2][7][3]);
                S[6][7][9] = gma2 * S[2][7][9] + fak * S[2][3][9];
                S[8][7][10] = gma0 * S[2][7][10] + fak * S[1][7][10];
                S[5][7][10] = gma1 * S[2][7][10] + fak * S[2][4][10];
                S[6][7][10] = gma2 * S[2][7][10] + fak * (S[2][3][10] + S[2][7][4] + S[2][7][4]);
                S[8][8][5] = gma0 * S[2][8][5] + fak * (S[1][8][5] + S[2][2][5] + S[2][2][5] + S[2][8][3]);
                S[5][8][5] = gma1 * S[2][8][5] + fak * S[2][8][2];
                S[6][8][5] = gma2 * S[2][8][5];
                S[8][8][6] = gma0 * S[2][8][6] + fak * (S[1][8][6] + S[2][2][6] + S[2][2][6] + S[2][8][4]);
                S[5][8][6] = gma1 * S[2][8][6];
                S[6][8][6] = gma2 * S[2][8][6] + fak * S[2][8][2];
                S[8][8][7] = gma0 * S[2][8][7] + fak * (S[1][8][7] + S[2][2][7] + S[2][2][7]);
                S[5][8][7] = gma1 * S[2][8][7] + fak * S[2][8][4];
                S[6][8][7] = gma2 * S[2][8][7] + fak * S[2][8][3];
                S[8][8][8] = gma0 * S[2][8][8] + fak * (S[1][8][8] + S[2][2][8] + S[2][2][8] + S[2][8][2] + S[2][8][2]);
                S[5][8][8] = gma1 * S[2][8][8];
                S[6][8][8] = gma2 * S[2][8][8];
                S[8][8][9] = gma0 * S[2][8][9] + fak * (S[1][8][9] + S[2][2][9] + S[2][2][9]);
                S[5][8][9] = gma1 * S[2][8][9] + fak2 * S[2][8][3];
                S[6][8][9] = gma2 * S[2][8][9];
                S[8][8][10] = gma0 * S[2][8][10] + fak * (S[1][8][10] + S[2][2][10] + S[2][2][10]);
                S[5][8][10] = gma1 * S[2][8][10];
                S[6][8][10] = gma2 * S[2][8][10] + fak2 * S[2][8][4];
                S[8][9][5] = gma0 * S[2][9][5] + fak * (S[1][9][5] + S[2][9][3]);
                S[5][9][5] = gma1 * S[2][9][5] + fak * (S[2][3][5] + S[2][3][5] + S[2][9][2]);
                S[6][9][5] = gma2 * S[2][9][5];
                S[8][9][6] = gma0 * S[2][9][6] + fak * (S[1][9][6] + S[2][9][4]);
                S[5][9][6] = gma1 * S[2][9][6] + fak2 * S[2][3][6];
                S[6][9][6] = gma2 * S[2][9][6] + fak * S[2][9][2];
                S[8][9][7] = gma0 * S[2][9][7] + fak * S[1][9][7];
                S[5][9][7] = gma1 * S[2][9][7] + fak * (S[2][3][7] + S[2][3][7] + S[2][9][4]);
                S[6][9][7] = gma2 * S[2][9][7] + fak * S[2][9][3];
                S[8][9][8] = gma0 * S[2][9][8] + fak * (S[1][9][8] + S[2][9][2] + S[2][9][2]);
                S[5][9][8] = gma1 * S[2][9][8] + fak2 * S[2][3][8];
                S[6][9][8] = gma2 * S[2][9][8];
                S[8][9][9] = gma0 * S[2][9][9] + fak * S[1][9][9];
                S[5][9][9] = gma1 * S[2][9][9] + fak2 * (S[2][3][9] + S[2][9][3]);
                S[6][9][9] = gma2 * S[2][9][9];
                S[8][9][10] = gma0 * S[2][9][10] + fak * S[1][9][10];
                S[5][9][10] = gma1 * S[2][9][10] + fak2 * S[2][3][10];
                S[6][9][10] = gma2 * S[2][9][10] + fak2 * S[2][9][4];
                S[8][10][5] = gma0 * S[2][10][5] + fak * (S[1][10][5] + S[2][10][3]);
                S[5][10][5] = gma1 * S[2][10][5] + fak * S[2][10][2];
                S[6][10][5] = gma2 * S[2][10][5] + fak2 * S[2][4][5];
                S[8][10][6] = gma0 * S[2][10][6] + fak * (S[1][10][6] + S[2][10][4]);
                S[5][10][6] = gma1 * S[2][10][6];
                S[6][10][6] = gma2 * S[2][10][6] + fak * (S[2][4][6] + S[2][4][6] + S[2][10][2]);
                S[8][10][7] = gma0 * S[2][10][7] + fak * S[1][10][7];
                S[5][10][7] = gma1 * S[2][10][7] + fak * S[2][10][4];
                S[6][10][7] = gma2 * S[2][10][7] + fak * (S[2][4][7] + S[2][4][7] + S[2][10][3]);
                S[8][10][8] = gma0 * S[2][10][8] + fak * (S[1][10][8] + S[2][10][2] + S[2][10][2]);
                S[5][10][8] = gma1 * S[2][10][8];
                S[6][10][8] = gma2 * S[2][10][8] + fak2 * S[2][4][8];
                S[8][10][9] = gma0 * S[2][10][9] + fak * S[1][10][9];
                S[5][10][9] = gma1 * S[2][10][9] + fak2 * S[2][10][3];
                S[6][10][9] = gma2 * S[2][10][9] + fak2 * S[2][4][9];
                S[8][10][10] = gma0 * S[2][10][10] + fak * S[1][10][10];
                S[5][10][10] = gma1 * S[2][10][10];
                S[6][10][10] = gma2 * S[2][10][10] + fak2 * (S[2][4][10] + S[2][10][4]);
                S[9][5][5] = gma1 * S[3][5][5] + fak * (S[1][5][5] + S[3][2][5] + S[3][5][2]);
                S[7][5][5] = gma2 * S[3][5][5];
                S[9][5][6] = gma1 * S[3][5][6] + fak * (S[1][5][6] + S[3][2][6]);
                S[7][5][6] = gma2 * S[3][5][6] + fak * S[3][5][2];
                S[9][5][7] = gma1 * S[3][5][7] + fak * (S[1][5][7] + S[3][2][7] + S[3][5][4]);
                S[7][5][7] = gma2 * S[3][5][7] + fak * S[3][5][3];
                S[9][5][8] = gma1 * S[3][5][8] + fak * (S[1][5][8] + S[3][2][8]);
                S[7][5][8] = gma2 * S[3][5][8];
                S[9][5][9] = gma1 * S[3][5][9] + fak * (S[1][5][9] + S[3][2][9] + S[3][5][3] + S[3][5][3]);
                S[7][5][9] = gma2 * S[3][5][9];
                S[9][5][10] = gma1 * S[3][5][10] + fak * (S[1][5][10] + S[3][2][10]);
                S[7][5][10] = gma2 * S[3][5][10] + fak2 * S[3][5][4];
                S[9][6][5] = gma1 * S[3][6][5] + fak * (S[1][6][5] + S[3][6][2]);
                S[7][6][5] = gma2 * S[3][6][5] + fak * S[3][2][5];
                S[9][6][6] = gma1 * S[3][6][6] + fak * S[1][6][6];
                S[7][6][6] = gma2 * S[3][6][6] + fak * (S[3][2][6] + S[3][6][2]);
                S[9][6][7] = gma1 * S[3][6][7] + fak * (S[1][6][7] + S[3][6][4]);
                S[7][6][7] = gma2 * S[3][6][7] + fak * (S[3][2][7] + S[3][6][3]);
                S[9][6][8] = gma1 * S[3][6][8] + fak * S[1][6][8];
                S[7][6][8] = gma2 * S[3][6][8] + fak * S[3][2][8];
                S[9][6][9] = gma1 * S[3][6][9] + fak * (S[1][6][9] + S[3][6][3] + S[3][6][3]);
                S[7][6][9] = gma2 * S[3][6][9] + fak * S[3][2][9];
                S[9][6][10] = gma1 * S[3][6][10] + fak * S[1][6][10];
                S[7][6][10] = gma2 * S[3][6][10] + fak * (S[3][2][10] + S[3][6][4] + S[3][6][4]);
                S[9][7][5] = gma1 * S[3][7][5] + fak * (S[1][7][5] + S[3][4][5] + S[3][7][2]);
                S[7][7][5] = gma2 * S[3][7][5] + fak * S[3][3][5];
                S[9][7][6] = gma1 * S[3][7][6] + fak * (S[1][7][6] + S[3][4][6]);
                S[7][7][6] = gma2 * S[3][7][6] + fak * (S[3][3][6] + S[3][7][2]);
                S[9][7][7] = gma1 * S[3][7][7] + fak * (S[1][7][7] + S[3][4][7] + S[3][7][4]);
                S[7][7][7] = gma2 * S[3][7][7] + fak * (S[3][3][7] + S[3][7][3]);
                S[9][7][8] = gma1 * S[3][7][8] + fak * (S[1][7][8] + S[3][4][8]);
                S[7][7][8] = gma2 * S[3][7][8] + fak * S[3][3][8];
                S[9][7][9] = gma1 * S[3][7][9] + fak * (S[1][7][9] + S[3][4][9] + S[3][7][3] + S[3][7][3]);
                S[7][7][9] = gma2 * S[3][7][9] + fak * S[3][3][9];
                S[9][7][10] = gma1 * S[3][7][10] + fak * (S[1][7][10] + S[3][4][10]);
                S[7][7][10] = gma2 * S[3][7][10] + fak * (S[3][3][10] + S[3][7][4] + S[3][7][4]);
                S[9][8][5] = gma1 * S[3][8][5] + fak * (S[1][8][5] + S[3][8][2]);
                S[7][8][5] = gma2 * S[3][8][5];
                S[9][8][6] = gma1 * S[3][8][6] + fak * S[1][8][6];
                S[7][8][6] = gma2 * S[3][8][6] + fak * S[3][8][2];
                S[9][8][7] = gma1 * S[3][8][7] + fak * (S[1][8][7] + S[3][8][4]);
                S[7][8][7] = gma2 * S[3][8][7] + fak * S[3][8][3];
                S[9][8][8] = gma1 * S[3][8][8] + fak * S[1][8][8];
                S[7][8][8] = gma2 * S[3][8][8];
                S[9][8][9] = gma1 * S[3][8][9] + fak * (S[1][8][9] + S[3][8][3] + S[3][8][3]);
                S[7][8][9] = gma2 * S[3][8][9];
                S[9][8][10] = gma1 * S[3][8][10] + fak * S[1][8][10];
                S[7][8][10] = gma2 * S[3][8][10] + fak2 * S[3][8][4];
                S[9][9][5] = gma1 * S[3][9][5] + fak * (S[1][9][5] + S[3][3][5] + S[3][3][5] + S[3][9][2]);
                S[7][9][5] = gma2 * S[3][9][5];
                S[9][9][6] = gma1 * S[3][9][6] + fak * (S[1][9][6] + S[3][3][6] + S[3][3][6]);
                S[7][9][6] = gma2 * S[3][9][6] + fak * S[3][9][2];
                S[9][9][7] = gma1 * S[3][9][7] + fak * (S[1][9][7] + S[3][3][7] + S[3][3][7] + S[3][9][4]);
                S[7][9][7] = gma2 * S[3][9][7] + fak * S[3][9][3];
                S[9][9][8] = gma1 * S[3][9][8] + fak * (S[1][9][8] + S[3][3][8] + S[3][3][8]);
                S[7][9][8] = gma2 * S[3][9][8];
                S[9][9][9] = gma1 * S[3][9][9] + fak * (S[1][9][9] + S[3][3][9] + S[3][3][9] + S[3][9][3] + S[3][9][3]);
                S[7][9][9] = gma2 * S[3][9][9];
                S[9][9][10] = gma1 * S[3][9][10] + fak * (S[1][9][10] + S[3][3][10] + S[3][3][10]);
                S[7][9][10] = gma2 * S[3][9][10] + fak2 * S[3][9][4];
                S[9][10][5] = gma1 * S[3][10][5] + fak * (S[1][10][5] + S[3][10][2]);
                S[7][10][5] = gma2 * S[3][10][5] + fak2 * S[3][4][5];
                S[9][10][6] = gma1 * S[3][10][6] + fak * S[1][10][6];
                S[7][10][6] = gma2 * S[3][10][6] + fak * (S[3][4][6] + S[3][4][6] + S[3][10][2]);
                S[9][10][7] = gma1 * S[3][10][7] + fak * (S[1][10][7] + S[3][10][4]);
                S[7][10][7] = gma2 * S[3][10][7] + fak * (S[3][4][7] + S[3][4][7] + S[3][10][3]);
                S[9][10][8] = gma1 * S[3][10][8] + fak * S[1][10][8];
                S[7][10][8] = gma2 * S[3][10][8] + fak2 * S[3][4][8];
                S[9][10][9] = gma1 * S[3][10][9] + fak * (S[1][10][9] + S[3][10][3] + S[3][10][3]);
                S[7][10][9] = gma2 * S[3][10][9] + fak2 * S[3][4][9];
                S[9][10][10] = gma1 * S[3][10][10] + fak * S[1][10][10];
                S[7][10][10] = gma2 * S[3][10][10] + fak2 * (S[3][4][10] + S[3][10][4]);
                S[10][5][5] = gma2 * S[4][5][5] + fak * S[1][5][5];
                S[10][5][6] = gma2 * S[4][5][6] + fak * (S[1][5][6] + S[4][5][2]);
                S[10][5][7] = gma2 * S[4][5][7] + fak * (S[1][5][7] + S[4][5][3]);
                S[10][5][8] = gma2 * S[4][5][8] + fak * S[1][5][8];
                S[10][5][9] = gma2 * S[4][5][9] + fak * S[1][5][9];
                S[10][5][10] = gma2 * S[4][5][10] + fak * (S[1][5][10] + S[4][5][4] + S[4][5][4]);
                S[10][6][5] = gma2 * S[4][6][5] + fak * (S[1][6][5] + S[4][2][5]);
                S[10][6][6] = gma2 * S[4][6][6] + fak * (S[1][6][6] + S[4][2][6] + S[4][6][2]);
                S[10][6][7] = gma2 * S[4][6][7] + fak * (S[1][6][7] + S[4][2][7] + S[4][6][3]);
                S[10][6][8] = gma2 * S[4][6][8] + fak * (S[1][6][8] + S[4][2][8]);
                S[10][6][9] = gma2 * S[4][6][9] + fak * (S[1][6][9] + S[4][2][9]);
                S[10][6][10] = gma2 * S[4][6][10] + fak * (S[1][6][10] + S[4][2][10] + S[4][6][4] + S[4][6][4]);
                S[10][7][5] = gma2 * S[4][7][5] + fak * (S[1][7][5] + S[4][3][5]);
                S[10][7][6] = gma2 * S[4][7][6] + fak * (S[1][7][6] + S[4][3][6] + S[4][7][2]);
                S[10][7][7] = gma2 * S[4][7][7] + fak * (S[1][7][7] + S[4][3][7] + S[4][7][3]);
                S[10][7][8] = gma2 * S[4][7][8] + fak * (S[1][7][8] + S[4][3][8]);
                S[10][7][9] = gma2 * S[4][7][9] + fak * (S[1][7][9] + S[4][3][9]);
                S[10][7][10] = gma2 * S[4][7][10] + fak * (S[1][7][10] + S[4][3][10] + S[4][7][4] + S[4][7][4]);
                S[10][8][5] = gma2 * S[4][8][5] + fak * (1 * S[1][8][5]);
                S[10][8][6] = gma2 * S[4][8][6] + fak * (S[1][8][6] + S[4][8][2]);
                S[10][8][7] = gma2 * S[4][8][7] + fak * (S[1][8][7] + S[4][8][3]);
                S[10][8][8] = gma2 * S[4][8][8] + fak * S[1][8][8];
                S[10][8][9] = gma2 * S[4][8][9] + fak * S[1][8][9];
                S[10][8][10] = gma2 * S[4][8][10] + fak * (S[1][8][10] + S[4][8][4] + S[4][8][4]);
                S[10][9][5] = gma2 * S[4][9][5] + fak * S[1][9][5];
                S[10][9][6] = gma2 * S[4][9][6] + fak * (S[1][9][6] + S[4][9][2]);
                S[10][9][7] = gma2 * S[4][9][7] + fak * (S[1][9][7] + S[4][9][3]);
                S[10][9][8] = gma2 * S[4][9][8] + fak * S[1][9][8];
                S[10][9][9] = gma2 * S[4][9][9] + fak * S[1][9][9];
                S[10][9][10] = gma2 * S[4][9][10] + fak * (S[1][9][10] + S[4][9][4] + S[4][9][4]);
                S[10][10][5] = gma2 * S[4][10][5] + fak * (S[1][10][5] + S[4][4][5] + S[4][4][5]);
                S[10][10][6] = gma2 * S[4][10][6] + fak * (S[1][10][6] + S[4][4][6] + S[4][4][6] + S[4][10][2]);
                S[10][10][7] = gma2 * S[4][10][7] + fak * (S[1][10][7] + S[4][4][7] + S[4][4][7] + S[4][10][3]);
                S[10][10][8] = gma2 * S[4][10][8] + fak * (S[1][10][8] + S[4][4][8] + S[4][4][8]);
                S[10][10][9] = gma2 * S[4][10][9] + fak * (S[1][10][9] + S[4][4][9] + S[4][4][9]);
                S[10][10][10] = gma2 * S[4][10][10] + fak * (S[1][10][10] + S[4][4][10] + S[4][4][10] + S[4][10][4] + S[4][10][4]);
            }
            
            // s-f-s elements
             if ( _lmax_alpha >= 0 && _lmax_gw >= 3 && _lmax_gamma >= 0) {
                S[1][14][1] = gmc0*S[1][5][1] + fak * S[1][3][1];
                S[1][15][1] = gmc1*S[1][5][1] + fak * S[1][2][1];
                S[1][20][1] = gmc2*S[1][5][1];
                S[1][16][1] = gmc0*S[1][6][1] + fak * S[1][4][1];
                S[1][17][1] = gmc2*S[1][6][1] + fak * S[1][2][1];
                S[1][18][1] = gmc1*S[1][7][1] + fak * S[1][4][1];
                S[1][19][1] = gmc2*S[1][7][1] + fak * S[1][3][1];
                S[1][11][1] = gmc0*S[1][8][1] + fak2* S[1][2][1];
                S[1][12][1] = gmc1*S[1][9][1] + fak2* S[1][3][1];
                S[1][13][1] = gmc2*S[1][10][1] + fak2* S[1][4][1];
             }

            // s-f-p elements 
             if ( _lmax_alpha >= 0 && _lmax_gw >= 3 && _lmax_gamma >= 1) {
                S[1][14][2] = gmc0*S[1][5][2] + fak * (S[1][3][2] +S[1][5][1] );
                S[1][15][2] = gmc1*S[1][5][2] + fak * S[1][2][2];
                S[1][20][2] = gmc2*S[1][5][2];
                S[1][14][3] = gmc0*S[1][5][3] + fak * S[1][3][3];
                S[1][15][3] = gmc1*S[1][5][3] + fak * (S[1][2][3] +S[1][5][1] );
                S[1][20][3] = gmc2*S[1][5][3];
                S[1][14][4] = gmc0*S[1][5][4] + fak * S[1][3][4];
                S[1][15][4] = gmc1*S[1][5][4] + fak * S[1][2][4];
                S[1][20][4] = gmc2*S[1][5][4] + fak * S[1][5][1];
                S[1][16][2] = gmc0*S[1][6][2] + fak * (S[1][4][2] +S[1][6][1] );
                S[1][17][2] = gmc2*S[1][6][2] + fak * S[1][2][2];
                S[1][16][3] = gmc0*S[1][6][3] + fak * S[1][4][3];
                S[1][17][3] = gmc2*S[1][6][3] + fak * S[1][2][3];
                S[1][16][4] = gmc0*S[1][6][4] + fak * S[1][4][4];
                S[1][17][4] = gmc2*S[1][6][4] + fak * (S[1][2][4] +S[1][6][1] );
                S[1][18][2] = gmc1*S[1][7][2] + fak * S[1][4][2] ;
                S[1][19][2] = gmc2*S[1][7][2] + fak * S[1][3][2];
                S[1][18][3] = gmc1*S[1][7][3] + fak * (S[1][4][3] +S[1][7][1] );
                S[1][19][3] = gmc2*S[1][7][3] + fak * S[1][3][3];
                S[1][18][4] = gmc1*S[1][7][4] + fak * S[1][4][4];
                S[1][19][4] = gmc2*S[1][7][4] + fak * (S[1][3][4] +S[1][7][1] );
                S[1][11][2] = gmc0*S[1][8][2] + fak * (2.0*S[1][2][2] +S[1][8][1] );
                S[1][11][3] = gmc0*S[1][8][3] + fak2* S[1][2][3];
                S[1][11][4] = gmc0*S[1][8][4] + fak2* S[1][2][4];
                S[1][12][2] = gmc1*S[1][9][2] + fak2* S[1][3][2];
                S[1][12][3] = gmc1*S[1][9][3] + fak * (2.0*S[1][3][3] +S[1][9][1] );
                S[1][12][4] = gmc1*S[1][9][4] + fak2* S[1][3][4];
                S[1][13][2] = gmc2*S[1][10][2] + fak2* S[1][4][2];
                S[1][13][3] = gmc2*S[1][10][3] + fak2* S[1][4][3];
                S[1][13][4] = gmc2*S[1][10][4] + fak * (2.0*S[1][4][4] +S[1][10][1] );
             }

            // p-f-s elements
             if ( _lmax_alpha >= 1 && _lmax_gw >= 3 && _lmax_gamma >= 0) {
                S[2][11][1] = gma0*S[1][11][1] + fak3* S[1][8][1];
                S[3][11][1] = gma1*S[1][11][1];
                S[4][11][1] = gma2*S[1][11][1];
                S[2][12][1] = gma0*S[1][12][1];
                S[3][12][1] = gma1*S[1][12][1] + fak3* S[1][9][1];
                S[4][12][1] = gma2*S[1][12][1];
                S[2][13][1] = gma0*S[1][13][1];
                S[3][13][1] = gma1*S[1][13][1];
                S[4][13][1] = gma2*S[1][13][1] + fak3* S[1][10][1];
                S[2][14][1] = gma0*S[1][14][1] + fak2* S[1][5][1];
                S[3][14][1] = gma1*S[1][14][1] + fak * S[1][8][1];
                S[4][14][1] = gma2*S[1][14][1];
                S[2][15][1] = gma0*S[1][15][1] + fak * S[1][9][1];
                S[3][15][1] = gma1*S[1][15][1] + fak2* S[1][5][1];
                S[4][15][1] = gma2*S[1][15][1];
                S[2][16][1] = gma0*S[1][16][1] + fak2* S[1][6][1];
                S[3][16][1] = gma1*S[1][16][1];
                S[4][16][1] = gma2*S[1][16][1] + fak * S[1][8][1];
                S[2][17][1] = gma0*S[1][17][1] + fak * S[1][10][1];
                S[3][17][1] = gma1*S[1][17][1];
                S[4][17][1] = gma2*S[1][17][1] + fak2* S[1][6][1];
                S[2][18][1] = gma0*S[1][18][1];
                S[3][18][1] = gma1*S[1][18][1] + fak2* S[1][7][1];
                S[4][18][1] = gma2*S[1][18][1] + fak * S[1][9][1];
                S[2][19][1] = gma0*S[1][19][1];
                S[3][19][1] = gma1*S[1][19][1] + fak * S[1][10][1];
                S[4][19][1] = gma2*S[1][19][1] + fak2* S[1][7][1];
                S[2][20][1] = gma0*S[1][20][1] + fak * S[1][7][1];
                S[3][20][1] = gma1*S[1][20][1] + fak * S[1][6][1];
                S[4][20][1] = gma2*S[1][20][1] + fak * S[1][5][1];
             }

            // s-f-d elements
             if ( _lmax_alpha >= 0 && _lmax_gw >= 3 && _lmax_gamma >= 2) {
                S[1][14][5] = gmc0*S[1][5][5] + fak * (S[1][3][5] +S[1][5][3] );
                S[1][15][5] = gmc1*S[1][5][5] + fak * (S[1][2][5] +S[1][5][2] );
                S[1][20][5] = gmc2*S[1][5][5];
                S[1][14][6] = gmc0*S[1][5][6] + fak * (S[1][3][6] +S[1][5][4] );
                S[1][15][6] = gmc1*S[1][5][6] + fak * S[1][2][6];
                S[1][20][6] = gmc2*S[1][5][6] + fak * S[1][5][2];
                S[1][14][7] = gmc0*S[1][5][7] + fak * S[1][3][7];
                S[1][15][7] = gmc1*S[1][5][7] + fak * (S[1][2][7] +S[1][5][4] );
                S[1][20][7] = gmc2*S[1][5][7] + fak * S[1][5][3];
                S[1][14][8] = gmc0*S[1][5][8] + fak * (S[1][3][8] +2.0*S[1][5][2] );
                S[1][15][8] = gmc1*S[1][5][8] + fak * S[1][2][8];
                S[1][20][8] = gmc2*S[1][5][8];
                S[1][14][9] = gmc0*S[1][5][9] + fak * S[1][3][9];
                S[1][15][9] = gmc1*S[1][5][9] + fak * (S[1][2][9] +2.0*S[1][5][3] );
                S[1][20][9] = gmc2*S[1][5][9];
                S[1][14][10] = gmc0*S[1][5][10] + fak * S[1][3][10];
                S[1][15][10] = gmc1*S[1][5][10] + fak * S[1][2][10];
                S[1][20][10] = gmc2*S[1][5][10] + fak2* S[1][5][4];
                S[1][16][5] = gmc0*S[1][6][5] + fak * (S[1][4][5] +S[1][6][3] );
                S[1][17][5] = gmc2*S[1][6][5] + fak * S[1][2][5];
                S[1][16][6] = gmc0*S[1][6][6] + fak * (S[1][4][6] +S[1][6][4] );
                S[1][17][6] = gmc2*S[1][6][6] + fak * (S[1][2][6] +S[1][6][2] );
                S[1][16][7] = gmc0*S[1][6][7] + fak * S[1][4][7];
                S[1][17][7] = gmc2*S[1][6][7] + fak * (S[1][2][7] +S[1][6][3] );
                S[1][16][8] = gmc0*S[1][6][8] + fak * (S[1][4][8] +2.0*S[1][6][2] );
                S[1][17][8] = gmc2*S[1][6][8] + fak * S[1][2][8];
                S[1][16][9] = gmc0*S[1][6][9] + fak * S[1][4][9];
                S[1][17][9] = gmc2*S[1][6][9] + fak * S[1][2][9];
                S[1][16][10] = gmc0*S[1][6][10] + fak * S[1][4][10];
                S[1][17][10] = gmc2*S[1][6][10] + fak * (S[1][2][10] +2.0*S[1][6][4] );
                S[1][18][5] = gmc1*S[1][7][5] + fak * (S[1][4][5] +S[1][7][2] );
                S[1][19][5] = gmc2*S[1][7][5] + fak * S[1][3][5];
                S[1][18][6] = gmc1*S[1][7][6] + fak * S[1][4][6];
                S[1][19][6] = gmc2*S[1][7][6] + fak * (S[1][3][6] +S[1][7][2] );
                S[1][18][7] = gmc1*S[1][7][7] + fak * (S[1][4][7] +S[1][7][4] );
                S[1][19][7] = gmc2*S[1][7][7] + fak * (S[1][3][7] +S[1][7][3] );
                S[1][18][8] = gmc1*S[1][7][8] + fak * S[1][4][8];
                S[1][19][8] = gmc2*S[1][7][8] + fak * S[1][3][8];
                S[1][18][9] = gmc1*S[1][7][9] + fak * (S[1][4][9] +2.0*S[1][7][3] );
                S[1][19][9] = gmc2*S[1][7][9] + fak * S[1][3][9];
                S[1][18][10] = gmc1*S[1][7][10] + fak * S[1][4][10];
                S[1][19][10] = gmc2*S[1][7][10] + fak * (S[1][3][10] +2.0*S[1][7][4] );
                S[1][11][5] = gmc0*S[1][8][5] + fak * (2.0*S[1][2][5] +S[1][8][3] );
                S[1][11][6] = gmc0*S[1][8][6] + fak * (2.0*S[1][2][6] +S[1][8][4] );
                S[1][11][7] = gmc0*S[1][8][7] + fak2* S[1][2][7];
                S[1][11][8] = gmc0*S[1][8][8] + fak2* (S[1][2][8] +S[1][8][2] );
                S[1][11][9] = gmc0*S[1][8][9] + fak2* S[1][2][9];
                S[1][11][10] = gmc0*S[1][8][10] + fak2* S[1][2][10];
                S[1][12][5] = gmc1*S[1][9][5] + fak * (2.0*S[1][3][5] +S[1][9][2] );
                S[1][12][6] = gmc1*S[1][9][6] + fak2* S[1][3][6];
                S[1][12][7] = gmc1*S[1][9][7] + fak * (2.0*S[1][3][7] +S[1][9][4] );
                S[1][12][8] = gmc1*S[1][9][8] + fak2* S[1][3][8];
                S[1][12][9] = gmc1*S[1][9][9] + fak2* (S[1][3][9] +S[1][9][3] );
                S[1][12][10] = gmc1*S[1][9][10] + fak2* S[1][3][10];
                S[1][13][5] = gmc2*S[1][10][5] + fak2* S[1][4][5];
                S[1][13][6] = gmc2*S[1][10][6] + fak * (2.0*S[1][4][6] +S[1][10][2] );
                S[1][13][7] = gmc2*S[1][10][7] + fak * (2.0*S[1][4][7] +S[1][10][3] );
                S[1][13][8] = gmc2*S[1][10][8] + fak2* S[1][4][8];
                S[1][13][9] = gmc2*S[1][10][9] + fak2* S[1][4][9];
                S[1][13][10] = gmc2*S[1][10][10] + fak2* (S[1][4][10] +S[1][10][4] );
             }

            // p-f-p elements
             if ( _lmax_alpha >= 1 && _lmax_gw >= 3 && _lmax_gamma >= 1) {
                S[2][11][2] = gma0*S[1][11][2] + fak * (3.0*S[1][8][2] +S[1][11][1] );
                S[3][11][2] = gma1*S[1][11][2];
                S[4][11][2] = gma2*S[1][11][2];
                S[2][11][3] = gma0*S[1][11][3] + fak3* S[1][8][3];
                S[3][11][3] = gma1*S[1][11][3] + fak * S[1][11][1];
                S[4][11][3] = gma2*S[1][11][3];
                S[2][11][4] = gma0*S[1][11][4] + fak3* S[1][8][4];
                S[3][11][4] = gma1*S[1][11][4];
                S[4][11][4] = gma2*S[1][11][4] + fak * S[1][11][1];
                S[2][12][2] = gma0*S[1][12][2] + fak * S[1][12][1];
                S[3][12][2] = gma1*S[1][12][2] + fak3* S[1][9][2];
                S[4][12][2] = gma2*S[1][12][2];
                S[2][12][3] = gma0*S[1][12][3];
                S[3][12][3] = gma1*S[1][12][3] + fak * (3.0*S[1][9][3] +S[1][12][1] );
                S[4][12][3] = gma2*S[1][12][3];
                S[2][12][4] = gma0*S[1][12][4];
                S[3][12][4] = gma1*S[1][12][4] + fak3* S[1][9][4];
                S[4][12][4] = gma2*S[1][12][4] + fak * S[1][12][1];
                S[2][13][2] = gma0*S[1][13][2] + fak * S[1][13][1];
                S[3][13][2] = gma1*S[1][13][2];
                S[4][13][2] = gma2*S[1][13][2] + fak3* S[1][10][2];
                S[2][13][3] = gma0*S[1][13][3];
                S[3][13][3] = gma1*S[1][13][3] + fak * S[1][13][1];
                S[4][13][3] = gma2*S[1][13][3] + fak3* S[1][10][3];
                S[2][13][4] = gma0*S[1][13][4];
                S[3][13][4] = gma1*S[1][13][4];
                S[4][13][4] = gma2*S[1][13][4] + fak * (3.0*S[1][10][4] +S[1][13][1] );
                S[2][14][2] = gma0*S[1][14][2] + fak * (2.0*S[1][5][2] +S[1][14][1] );
                S[3][14][2] = gma1*S[1][14][2] + fak * S[1][8][2];
                S[4][14][2] = gma2*S[1][14][2];
                S[2][14][3] = gma0*S[1][14][3] + fak2* S[1][5][3];
                S[3][14][3] = gma1*S[1][14][3] + fak * (S[1][8][3] +S[1][14][1] );
                S[4][14][3] = gma2*S[1][14][3];
                S[2][14][4] = gma0*S[1][14][4] + fak2* S[1][5][4];
                S[3][14][4] = gma1*S[1][14][4] + fak * S[1][8][4];
                S[4][14][4] = gma2*S[1][14][4] + fak * S[1][14][1];
                S[2][15][2] = gma0*S[1][15][2] + fak * (S[1][9][2] +S[1][15][1] );
                S[3][15][2] = gma1*S[1][15][2] + fak2* S[1][5][2];
                S[4][15][2] = gma2*S[1][15][2];
                S[2][15][3] = gma0*S[1][15][3] + fak * S[1][9][3];
                S[3][15][3] = gma1*S[1][15][3] + fak * (2.0*S[1][5][3] +S[1][15][1] );
                S[4][15][3] = gma2*S[1][15][3];
                S[2][15][4] = gma0*S[1][15][4] + fak * S[1][9][4];
                S[3][15][4] = gma1*S[1][15][4] + fak2* S[1][5][4];
                S[4][15][4] = gma2*S[1][15][4] + fak * S[1][15][1];
                S[2][16][2] = gma0*S[1][16][2] + fak * (2.0*S[1][6][2] +S[1][16][1] );
                S[3][16][2] = gma1*S[1][16][2];
                S[4][16][2] = gma2*S[1][16][2] + fak * S[1][8][2];
                S[2][16][3] = gma0*S[1][16][3] + fak2* S[1][6][3];
                S[3][16][3] = gma1*S[1][16][3] + fak * S[1][16][1];
                S[4][16][3] = gma2*S[1][16][3] + fak * S[1][8][3];
                S[2][16][4] = gma0*S[1][16][4] + fak2* S[1][6][4];
                S[3][16][4] = gma1*S[1][16][4];
                S[4][16][4] = gma2*S[1][16][4] + fak * (S[1][8][4] +S[1][16][1] );
                S[2][17][2] = gma0*S[1][17][2] + fak * (S[1][10][2] +S[1][17][1] );
                S[3][17][2] = gma1*S[1][17][2];
                S[4][17][2] = gma2*S[1][17][2] + fak2* S[1][6][2];
                S[2][17][3] = gma0*S[1][17][3] + fak * S[1][10][3];
                S[3][17][3] = gma1*S[1][17][3] + fak * S[1][17][1];
                S[4][17][3] = gma2*S[1][17][3] + fak2* S[1][6][3];
                S[2][17][4] = gma0*S[1][17][4] + fak * S[1][10][4];
                S[3][17][4] = gma1*S[1][17][4];
                S[4][17][4] = gma2*S[1][17][4] + fak * (2.0*S[1][6][4] +S[1][17][1] );
                S[2][18][2] = gma0*S[1][18][2] + fak * S[1][18][1];
                S[3][18][2] = gma1*S[1][18][2] + fak2* S[1][7][2];
                S[4][18][2] = gma2*S[1][18][2] + fak * S[1][9][2];
                S[2][18][3] = gma0*S[1][18][3];
                S[3][18][3] = gma1*S[1][18][3] + fak * (2.0*S[1][7][3] +S[1][18][1] );
                S[4][18][3] = gma2*S[1][18][3] + fak * S[1][9][3];
                S[2][18][4] = gma0*S[1][18][4];
                S[3][18][4] = gma1*S[1][18][4] + fak2* S[1][7][4];
                S[4][18][4] = gma2*S[1][18][4] + fak * (S[1][9][4] +S[1][18][1] );
                S[2][19][2] = gma0*S[1][19][2] + fak * S[1][19][1];
                S[3][19][2] = gma1*S[1][19][2] + fak * S[1][10][2];
                S[4][19][2] = gma2*S[1][19][2] + fak2* S[1][7][2];
                S[2][19][3] = gma0*S[1][19][3];
                S[3][19][3] = gma1*S[1][19][3] + fak * (S[1][10][3] +S[1][19][1] );
                S[4][19][3] = gma2*S[1][19][3] + fak2* S[1][7][3];
                S[2][19][4] = gma0*S[1][19][4];
                S[3][19][4] = gma1*S[1][19][4] + fak * S[1][10][4];
                S[4][19][4] = gma2*S[1][19][4] + fak * (2.0*S[1][7][4] +S[1][19][1] );
                S[2][20][2] = gma0*S[1][20][2] + fak * (S[1][7][2] +S[1][20][1] );
                S[3][20][2] = gma1*S[1][20][2] + fak * S[1][6][2];
                S[4][20][2] = gma2*S[1][20][2] + fak * S[1][5][2];
                S[2][20][3] = gma0*S[1][20][3] + fak * S[1][7][3];
                S[3][20][3] = gma1*S[1][20][3] + fak * (S[1][6][3] +S[1][20][1] );
                S[4][20][3] = gma2*S[1][20][3] + fak * S[1][5][3];
                S[2][20][4] = gma0*S[1][20][4] + fak * S[1][7][4];
                S[3][20][4] = gma1*S[1][20][4] + fak * S[1][6][4];
                S[4][20][4] = gma2*S[1][20][4] + fak * (S[1][5][4] +S[1][20][1] );
             }

            // d-f-s
             if ( _lmax_alpha >= 2 && _lmax_gw >= 3 && _lmax_gamma >= 0) {
                S[8][11][1] = gma0*S[2][11][1] + fak * (S[1][11][1] +3.0*S[2][8][1] );
                S[5][11][1] = gma1*S[2][11][1];
                S[6][11][1] = gma2*S[2][11][1];
                S[8][12][1] = gma0*S[2][12][1] + fak * S[1][12][1];
                S[5][12][1] = gma1*S[2][12][1] + fak3* S[2][9][1];
                S[6][12][1] = gma2*S[2][12][1];
                S[8][13][1] = gma0*S[2][13][1] + fak * S[1][13][1];
                S[5][13][1] = gma1*S[2][13][1];
                S[6][13][1] = gma2*S[2][13][1] + fak3* S[2][10][1];
                S[8][14][1] = gma0*S[2][14][1] + fak * (S[1][14][1] +2.0*S[2][5][1] );
                S[5][14][1] = gma1*S[2][14][1] + fak * S[2][8][1];
                S[6][14][1] = gma2*S[2][14][1];
                S[8][15][1] = gma0*S[2][15][1] + fak * (S[1][15][1] +S[2][9][1] );
                S[5][15][1] = gma1*S[2][15][1] + fak2* S[2][5][1];
                S[6][15][1] = gma2*S[2][15][1];
                S[8][16][1] = gma0*S[2][16][1] + fak * (S[1][16][1] +2.0*S[2][6][1] );
                S[5][16][1] = gma1*S[2][16][1];
                S[6][16][1] = gma2*S[2][16][1] + fak * S[2][8][1];
                S[8][17][1] = gma0*S[2][17][1] + fak * (S[1][17][1] +S[2][10][1] );
                S[5][17][1] = gma1*S[2][17][1];
                S[6][17][1] = gma2*S[2][17][1] + fak2* S[2][6][1];
                S[8][18][1] = gma0*S[2][18][1] + fak * S[1][18][1];
                S[5][18][1] = gma1*S[2][18][1] + fak2* S[2][7][1];
                S[6][18][1] = gma2*S[2][18][1] + fak * S[2][9][1];
                S[8][19][1] = gma0*S[2][19][1] + fak * S[1][19][1];
                S[5][19][1] = gma1*S[2][19][1] + fak * S[2][10][1];
                S[6][19][1] = gma2*S[2][19][1] + fak2* S[2][7][1];
                S[8][20][1] = gma0*S[2][20][1] + fak * (S[1][20][1] +S[2][7][1] );
                S[5][20][1] = gma1*S[2][20][1] + fak * S[2][6][1];
                S[6][20][1] = gma2*S[2][20][1] + fak * S[2][5][1];
                S[9][11][1] = gma1*S[3][11][1] + fak * S[1][11][1];
                S[7][11][1] = gma2*S[3][11][1];
                S[9][12][1] = gma1*S[3][12][1] + fak * (S[1][12][1] +3.0*S[3][9][1] );
                S[7][12][1] = gma2*S[3][12][1];
                S[9][13][1] = gma1*S[3][13][1] + fak * S[1][13][1];
                S[7][13][1] = gma2*S[3][13][1] + fak3* S[3][10][1];
                S[9][14][1] = gma1*S[3][14][1] + fak * (S[1][14][1] +S[3][8][1] );
                S[7][14][1] = gma2*S[3][14][1];
                S[9][15][1] = gma1*S[3][15][1] + fak * (S[1][15][1] +2.0*S[3][5][1] );
                S[7][15][1] = gma2*S[3][15][1];
                S[9][16][1] = gma1*S[3][16][1] + fak * S[1][16][1];
                S[7][16][1] = gma2*S[3][16][1] + fak * S[3][8][1];
                S[9][17][1] = gma1*S[3][17][1] + fak * S[1][17][1];
                S[7][17][1] = gma2*S[3][17][1] + fak2* S[3][6][1];
                S[9][18][1] = gma1*S[3][18][1] + fak * (S[1][18][1] +2.0*S[3][7][1] );
                S[7][18][1] = gma2*S[3][18][1] + fak * S[3][9][1];
                S[9][19][1] = gma1*S[3][19][1] + fak * (S[1][19][1] +S[3][10][1] );
                S[7][19][1] = gma2*S[3][19][1] + fak2* S[3][7][1];
                S[9][20][1] = gma1*S[3][20][1] + fak * (S[1][20][1] +S[3][6][1] );
                S[7][20][1] = gma2*S[3][20][1] + fak * S[3][5][1];
                S[10][11][1] = gma2*S[4][11][1] + fak * S[1][11][1];
                S[10][12][1] = gma2*S[4][12][1] + fak * S[1][12][1];
                S[10][13][1] = gma2*S[4][13][1] + fak * (S[1][13][1] +3.0*S[4][10][1] );
                S[10][14][1] = gma2*S[4][14][1] + fak * S[1][14][1];
                S[10][15][1] = gma2*S[4][15][1] + fak * S[1][15][1];
                S[10][16][1] = gma2*S[4][16][1] + fak * (S[1][16][1] +S[4][8][1] );
                S[10][17][1] = gma2*S[4][17][1] + fak * (S[1][17][1] +2.0*S[4][6][1] );
                S[10][18][1] = gma2*S[4][18][1] + fak * (S[1][18][1] +S[4][9][1] );
                S[10][19][1] = gma2*S[4][19][1] + fak * (S[1][19][1] +2.0*S[4][7][1] );
                S[10][20][1] = gma2*S[4][20][1] + fak * (S[1][20][1] +S[4][5][1] );
             }

            // p-f-d elements
             if ( _lmax_alpha >= 1 && _lmax_gw >= 3 && _lmax_gamma >= 2) {
                S[2][11][5] = gma0*S[1][11][5] + fak * (3.0*S[1][8][5] +S[1][11][3] );
                S[3][11][5] = gma1*S[1][11][5] + fak * S[1][11][2];
                S[4][11][5] = gma2*S[1][11][5];
                S[2][11][6] = gma0*S[1][11][6] + fak * (3.0*S[1][8][6] +S[1][11][4] );
                S[3][11][6] = gma1*S[1][11][6];
                S[4][11][6] = gma2*S[1][11][6] + fak * S[1][11][2];
                S[2][11][7] = gma0*S[1][11][7] + fak3* S[1][8][7];
                S[3][11][7] = gma1*S[1][11][7] + fak * S[1][11][4];
                S[4][11][7] = gma2*S[1][11][7] + fak * S[1][11][3];
                S[2][11][8] = gma0*S[1][11][8] + fak *(3.0*S[1][8][8]+2.0*S[1][11][2]);
                S[3][11][8] = gma1*S[1][11][8];
                S[4][11][8] = gma2*S[1][11][8];
                S[2][11][9] = gma0*S[1][11][9] + fak3* S[1][8][9];
                S[3][11][9] = gma1*S[1][11][9] + fak2* S[1][11][3];
                S[4][11][9] = gma2*S[1][11][9];
                S[2][11][10] = gma0*S[1][11][10] + fak3* S[1][8][10];
                S[3][11][10] = gma1*S[1][11][10];
                S[4][11][10] = gma2*S[1][11][10] + fak2* S[1][11][4];
                S[2][12][5] = gma0*S[1][12][5] + fak * S[1][12][3];
                S[3][12][5] = gma1*S[1][12][5] + fak * (3.0*S[1][9][5] +S[1][12][2] );
                S[4][12][5] = gma2*S[1][12][5];
                S[2][12][6] = gma0*S[1][12][6] + fak * S[1][12][4];
                S[3][12][6] = gma1*S[1][12][6] + fak3* S[1][9][6];
                S[4][12][6] = gma2*S[1][12][6] + fak * S[1][12][2];
                S[2][12][7] = gma0*S[1][12][7];
                S[3][12][7] = gma1*S[1][12][7] + fak * (3.0*S[1][9][7] +S[1][12][4] );
                S[4][12][7] = gma2*S[1][12][7] + fak * S[1][12][3];
                S[2][12][8] = gma0*S[1][12][8] + fak2* S[1][12][2];
                S[3][12][8] = gma1*S[1][12][8] + fak3* S[1][9][8];
                S[4][12][8] = gma2*S[1][12][8];
                S[2][12][9] = gma0*S[1][12][9];
                S[3][12][9] = gma1*S[1][12][9] + fak *(3.0*S[1][9][9]+2.0*S[1][12][3]);
                S[4][12][9] = gma2*S[1][12][9];
                S[2][12][10] = gma0*S[1][12][10];
                S[3][12][10] = gma1*S[1][12][10] + fak3* S[1][9][10];
                S[4][12][10] = gma2*S[1][12][10] + fak2* S[1][12][4];
                S[2][13][5] = gma0*S[1][13][5] + fak * S[1][13][3];
                S[3][13][5] = gma1*S[1][13][5] + fak * S[1][13][2];
                S[4][13][5] = gma2*S[1][13][5] + fak3* S[1][10][5];
                S[2][13][6] = gma0*S[1][13][6] + fak * S[1][13][4];
                S[3][13][6] = gma1*S[1][13][6];
                S[4][13][6] = gma2*S[1][13][6] + fak * (3.0*S[1][10][6] +S[1][13][2] );
                S[2][13][7] = gma0*S[1][13][7];
                S[3][13][7] = gma1*S[1][13][7] + fak * S[1][13][4];
                S[4][13][7] = gma2*S[1][13][7] + fak * (3.0*S[1][10][7] +S[1][13][3] );
                S[2][13][8] = gma0*S[1][13][8] + fak2* S[1][13][2];
                S[3][13][8] = gma1*S[1][13][8];
                S[4][13][8] = gma2*S[1][13][8] + fak3* S[1][10][8];
                S[2][13][9] = gma0*S[1][13][9];
                S[3][13][9] = gma1*S[1][13][9] + fak2* S[1][13][3];
                S[4][13][9] = gma2*S[1][13][9] + fak3* S[1][10][9];
                S[2][13][10] = gma0*S[1][13][10];
                S[3][13][10] = gma1*S[1][13][10];
                S[4][13][10] = gma2*S[1][13][10] +    fak * (3.0*S[1][10][10] +2.0*S[1][13][4] );
                S[2][14][5] = gma0*S[1][14][5] + fak * (2.0*S[1][5][5] +S[1][14][3] );
                S[3][14][5] = gma1*S[1][14][5] + fak * (S[1][8][5] +S[1][14][2] );
                S[4][14][5] = gma2*S[1][14][5];
                S[2][14][6] = gma0*S[1][14][6] + fak * (2.0*S[1][5][6] +S[1][14][4] );
                S[3][14][6] = gma1*S[1][14][6] + fak * S[1][8][6];
                S[4][14][6] = gma2*S[1][14][6] + fak * S[1][14][2];
                S[2][14][7] = gma0*S[1][14][7] + fak2* S[1][5][7];
                S[3][14][7] = gma1*S[1][14][7] + fak * (S[1][8][7] +S[1][14][4] );
                S[4][14][7] = gma2*S[1][14][7] + fak * S[1][14][3];
                S[2][14][8] = gma0*S[1][14][8] + fak2* (S[1][5][8] +S[1][14][2] );
                S[3][14][8] = gma1*S[1][14][8] + fak * S[1][8][8];
                S[4][14][8] = gma2*S[1][14][8];
                S[2][14][9] = gma0*S[1][14][9] + fak2* S[1][5][9];
                S[3][14][9] = gma1*S[1][14][9] + fak * (S[1][8][9] +2.0*S[1][14][3] );
                S[4][14][9] = gma2*S[1][14][9];
                S[2][14][10] = gma0*S[1][14][10] + fak2* S[1][5][10];
                S[3][14][10] = gma1*S[1][14][10] + fak * S[1][8][10];
                S[4][14][10] = gma2*S[1][14][10] + fak2* S[1][14][4];
                S[2][15][5] = gma0*S[1][15][5] + fak * (S[1][9][5] +S[1][15][3] );
                S[3][15][5] = gma1*S[1][15][5] + fak * (2.0*S[1][5][5] +S[1][15][2] );
                S[4][15][5] = gma2*S[1][15][5];
                S[2][15][6] = gma0*S[1][15][6] + fak * (S[1][9][6] +S[1][15][4] );
                S[3][15][6] = gma1*S[1][15][6] + fak2* S[1][5][6];
                S[4][15][6] = gma2*S[1][15][6] + fak * S[1][15][2];
                S[2][15][7] = gma0*S[1][15][7] + fak * S[1][9][7];
                S[3][15][7] = gma1*S[1][15][7] + fak * (2.0*S[1][5][7] +S[1][15][4] );
                S[4][15][7] = gma2*S[1][15][7] + fak * S[1][15][3];
                S[2][15][8] = gma0*S[1][15][8] + fak * (S[1][9][8] +2.0*S[1][15][2] );
                S[3][15][8] = gma1*S[1][15][8] + fak2* S[1][5][8];
                S[4][15][8] = gma2*S[1][15][8];
                S[2][15][9] = gma0*S[1][15][9] + fak * S[1][9][9];
                S[3][15][9] = gma1*S[1][15][9] + fak2* (S[1][5][9] +S[1][15][3] );
                S[4][15][9] = gma2*S[1][15][9];
                S[2][15][10] = gma0*S[1][15][10] + fak * S[1][9][10];
                S[3][15][10] = gma1*S[1][15][10] + fak2* S[1][5][10];
                S[4][15][10] = gma2*S[1][15][10] + fak2* S[1][15][4];
                S[2][16][5] = gma0*S[1][16][5] + fak * (2.0*S[1][6][5] +S[1][16][3] );
                S[3][16][5] = gma1*S[1][16][5] + fak * S[1][16][2];
                S[4][16][5] = gma2*S[1][16][5] + fak * S[1][8][5];
                S[2][16][6] = gma0*S[1][16][6] + fak * (2.0*S[1][6][6] +S[1][16][4] );
                S[3][16][6] = gma1*S[1][16][6];
                S[4][16][6] = gma2*S[1][16][6] + fak * (S[1][8][6] +S[1][16][2] );
                S[2][16][7] = gma0*S[1][16][7] + fak2* S[1][6][7];
                S[3][16][7] = gma1*S[1][16][7] + fak * S[1][16][4];
                S[4][16][7] = gma2*S[1][16][7] + fak * (S[1][8][7] +S[1][16][3] );
                S[2][16][8] = gma0*S[1][16][8] + fak2* (S[1][6][8] +S[1][16][2] );
                S[3][16][8] = gma1*S[1][16][8];
                S[4][16][8] = gma2*S[1][16][8] + fak * S[1][8][8];
                S[2][16][9] = gma0*S[1][16][9] + fak2* S[1][6][9];
                S[3][16][9] = gma1*S[1][16][9] + fak2* S[1][16][3];
                S[4][16][9] = gma2*S[1][16][9] + fak * S[1][8][9];
                S[2][16][10] = gma0*S[1][16][10] + fak2* S[1][6][10];
                S[3][16][10] = gma1*S[1][16][10];
                S[4][16][10] = gma2*S[1][16][10] + fak * (S[1][8][10] +2.0*S[1][16][4]);
                S[2][17][5] = gma0*S[1][17][5] + fak * (S[1][10][5] +S[1][17][3] );
                S[3][17][5] = gma1*S[1][17][5] + fak * S[1][17][2];
                S[4][17][5] = gma2*S[1][17][5] + fak2* S[1][6][5];
                S[2][17][6] = gma0*S[1][17][6] + fak * (S[1][10][6] +S[1][17][4] );
                S[3][17][6] = gma1*S[1][17][6];
                S[4][17][6] = gma2*S[1][17][6] + fak * (2.0*S[1][6][6] +S[1][17][2] );
                S[2][17][7] = gma0*S[1][17][7] + fak * S[1][10][7];
                S[3][17][7] = gma1*S[1][17][7] + fak * S[1][17][4];
                S[4][17][7] = gma2*S[1][17][7] + fak * (2.0*S[1][6][7] +S[1][17][3] );
                S[2][17][8] = gma0*S[1][17][8] + fak * (S[1][10][8] +2.0*S[1][17][2] );
                S[3][17][8] = gma1*S[1][17][8];
                S[4][17][8] = gma2*S[1][17][8] + fak2* S[1][6][8];
                S[2][17][9] = gma0*S[1][17][9] + fak * S[1][10][9];
                S[3][17][9] = gma1*S[1][17][9] + fak2* S[1][17][3];
                S[4][17][9] = gma2*S[1][17][9] + fak2* S[1][6][9];
                S[2][17][10] = gma0*S[1][17][10] + fak * S[1][10][10];
                S[3][17][10] = gma1*S[1][17][10];
                S[4][17][10] = gma2*S[1][17][10] + fak2* (S[1][6][10] +S[1][17][4] );
                S[2][18][5] = gma0*S[1][18][5] + fak * S[1][18][3];
                S[3][18][5] = gma1*S[1][18][5] + fak * (2.0*S[1][7][5] +S[1][18][2] );
                S[4][18][5] = gma2*S[1][18][5] + fak * S[1][9][5];
                S[2][18][6] = gma0*S[1][18][6] + fak * S[1][18][4];
                S[3][18][6] = gma1*S[1][18][6] + fak2* S[1][7][6];
                S[4][18][6] = gma2*S[1][18][6] + fak * (S[1][9][6] +S[1][18][2] );
                S[2][18][7] = gma0*S[1][18][7];
                S[3][18][7] = gma1*S[1][18][7] + fak * (2.0*S[1][7][7] +S[1][18][4] );
                S[4][18][7] = gma2*S[1][18][7] + fak * (S[1][9][7] +S[1][18][3] );
                S[2][18][8] = gma0*S[1][18][8] + fak2* S[1][18][2];
                S[3][18][8] = gma1*S[1][18][8] + fak2* S[1][7][8];
                S[4][18][8] = gma2*S[1][18][8] + fak * S[1][9][8];
                S[2][18][9] = gma0*S[1][18][9];
                S[3][18][9] = gma1*S[1][18][9] + fak2* (S[1][7][9] +S[1][18][3] );
                S[4][18][9] = gma2*S[1][18][9] + fak * S[1][9][9];
                S[2][18][10] = gma0*S[1][18][10];
                S[3][18][10] = gma1*S[1][18][10] + fak2* S[1][7][10];
                S[4][18][10] = gma2*S[1][18][10] + fak * (S[1][9][10] +2.0*S[1][18][4]);
                S[2][19][5] = gma0*S[1][19][5] + fak * S[1][19][3];
                S[3][19][5] = gma1*S[1][19][5] + fak * (S[1][10][5] +S[1][19][2] );
                S[4][19][5] = gma2*S[1][19][5] + fak2* S[1][7][5];
                S[2][19][6] = gma0*S[1][19][6] + fak * S[1][19][4];
                S[3][19][6] = gma1*S[1][19][6] + fak * S[1][10][6];
                S[4][19][6] = gma2*S[1][19][6] + fak * (2.0*S[1][7][6] +S[1][19][2] );
                S[2][19][7] = gma0*S[1][19][7];
                S[3][19][7] = gma1*S[1][19][7] + fak * (S[1][10][7] +S[1][19][4] );
                S[4][19][7] = gma2*S[1][19][7] + fak * (2.0*S[1][7][7] +S[1][19][3] );
                S[2][19][8] = gma0*S[1][19][8] + fak2* S[1][19][2];
                S[3][19][8] = gma1*S[1][19][8] + fak * S[1][10][8];
                S[4][19][8] = gma2*S[1][19][8] + fak2* S[1][7][8];
                S[2][19][9] = gma0*S[1][19][9];
                S[3][19][9] = gma1*S[1][19][9] + fak * (S[1][10][9] +2.0*S[1][19][3] );
                S[4][19][9] = gma2*S[1][19][9] + fak2* S[1][7][9];
                S[2][19][10] = gma0*S[1][19][10];
                S[3][19][10] = gma1*S[1][19][10] + fak * S[1][10][10];
                S[4][19][10] = gma2*S[1][19][10] + fak2* (S[1][7][10] +S[1][19][4] );
                S[2][20][5] = gma0*S[1][20][5] + fak * (S[1][7][5] +S[1][20][3] );
                S[3][20][5] = gma1*S[1][20][5] + fak * (S[1][6][5] +S[1][20][2] );
                S[4][20][5] = gma2*S[1][20][5] + fak * S[1][5][5];
                S[2][20][6] = gma0*S[1][20][6] + fak * (S[1][7][6] +S[1][20][4] );
                S[3][20][6] = gma1*S[1][20][6] + fak * S[1][6][6];
                S[4][20][6] = gma2*S[1][20][6] + fak * (S[1][5][6] +S[1][20][2] );
                S[2][20][7] = gma0*S[1][20][7] + fak * S[1][7][7];
                S[3][20][7] = gma1*S[1][20][7] + fak * (S[1][6][7] +S[1][20][4] );
                S[4][20][7] = gma2*S[1][20][7] + fak * (S[1][5][7] +S[1][20][3] );
                S[2][20][8] = gma0*S[1][20][8] + fak * (S[1][7][8] +2.0*S[1][20][2] );
                S[3][20][8] = gma1*S[1][20][8] + fak * S[1][6][8];
                S[4][20][8] = gma2*S[1][20][8] + fak * S[1][5][8];
                S[2][20][9] = gma0*S[1][20][9] + fak * S[1][7][9];
                S[3][20][9] = gma1*S[1][20][9] + fak * (S[1][6][9] +2.0*S[1][20][3] );
                S[4][20][9] = gma2*S[1][20][9] + fak * S[1][5][9];
                S[2][20][10] = gma0*S[1][20][10] + fak * S[1][7][10];
                S[3][20][10] = gma1*S[1][20][10] + fak * S[1][6][10];
                S[4][20][10] = gma2*S[1][20][10] + fak * (S[1][5][10] +2.0*S[1][20][4]);
             }

            // d-f-p elements
             if ( _lmax_alpha >= 2 && _lmax_gw >= 3 && _lmax_gamma >= 1) {
                S[8][11][2] = gma0*S[2][11][2] + fak * (S[1][11][2] +3.0*S[2][8][2] +S[2][11][1] );
                S[5][11][2] = gma1*S[2][11][2];
                S[6][11][2] = gma2*S[2][11][2];
                S[8][11][3] = gma0*S[2][11][3] + fak * (S[1][11][3] +3.0*S[2][8][3] );
                S[5][11][3] = gma1*S[2][11][3] + fak * S[2][11][1];
                S[6][11][3] = gma2*S[2][11][3];
                S[8][11][4] = gma0*S[2][11][4] + fak * (S[1][11][4] +3.0*S[2][8][4] );
                S[5][11][4] = gma1*S[2][11][4];
                S[6][11][4] = gma2*S[2][11][4] + fak * S[2][11][1];
                S[8][12][2] = gma0*S[2][12][2] + fak * (S[1][12][2] +S[2][12][1] );
                S[5][12][2] = gma1*S[2][12][2] + fak3* S[2][9][2];
                S[6][12][2] = gma2*S[2][12][2];
                S[8][12][3] = gma0*S[2][12][3] + fak * S[1][12][3];
                S[5][12][3] = gma1*S[2][12][3] + fak * (3.0*S[2][9][3] +S[2][12][1] );
                S[6][12][3] = gma2*S[2][12][3];
                S[8][12][4] = gma0*S[2][12][4] + fak * S[1][12][4];
                S[5][12][4] = gma1*S[2][12][4] + fak3* S[2][9][4];
                S[6][12][4] = gma2*S[2][12][4] + fak * S[2][12][1];
                S[8][13][2] = gma0*S[2][13][2] + fak * (S[1][13][2] +S[2][13][1] );
                S[5][13][2] = gma1*S[2][13][2];
                S[6][13][2] = gma2*S[2][13][2] + fak3* S[2][10][2];
                S[8][13][3] = gma0*S[2][13][3] + fak * S[1][13][3];
                S[5][13][3] = gma1*S[2][13][3] + fak * S[2][13][1];
                S[6][13][3] = gma2*S[2][13][3] + fak3* S[2][10][3];
                S[8][13][4] = gma0*S[2][13][4] + fak * S[1][13][4];
                S[5][13][4] = gma1*S[2][13][4];
                S[6][13][4] = gma2*S[2][13][4] + fak * (3.0*S[2][10][4] +S[2][13][1] );
                S[8][14][2] = gma0*S[2][14][2] + fak * (S[1][14][2] +2.0*S[2][5][2] +S[2][14][1] );
                S[5][14][2] = gma1*S[2][14][2] + fak * S[2][8][2];
                S[6][14][2] = gma2*S[2][14][2];
                S[8][14][3] = gma0*S[2][14][3] + fak * (S[1][14][3] +2.0*S[2][5][3] );
                S[5][14][3] = gma1*S[2][14][3] + fak * (S[2][8][3] +S[2][14][1] );
                S[6][14][3] = gma2*S[2][14][3];
                S[8][14][4] = gma0*S[2][14][4] + fak * (S[1][14][4] +2.0*S[2][5][4] );
                S[5][14][4] = gma1*S[2][14][4] + fak * S[2][8][4];
                S[6][14][4] = gma2*S[2][14][4] + fak * S[2][14][1];
                S[8][15][2] = gma0*S[2][15][2] + fak *(S[1][15][2]+S[2][9][2]+S[2][15][1]);
                S[5][15][2] = gma1*S[2][15][2] + fak2* S[2][5][2];
                S[6][15][2] = gma2*S[2][15][2];
                S[8][15][3] = gma0*S[2][15][3] + fak * (S[1][15][3] +S[2][9][3] );
                S[5][15][3] = gma1*S[2][15][3] + fak * (2.0*S[2][5][3] +S[2][15][1] );
                S[6][15][3] = gma2*S[2][15][3];
                S[8][15][4] = gma0*S[2][15][4] + fak * (S[1][15][4] +S[2][9][4] );
                S[5][15][4] = gma1*S[2][15][4] + fak2* S[2][5][4];
                S[6][15][4] = gma2*S[2][15][4] + fak * S[2][15][1];
                S[8][16][2] = gma0*S[2][16][2] + fak * (S[1][16][2] +2.0*S[2][6][2] +S[2][16][1] );
                S[5][16][2] = gma1*S[2][16][2];
                S[6][16][2] = gma2*S[2][16][2] + fak * S[2][8][2];
                S[8][16][3] = gma0*S[2][16][3] + fak * (S[1][16][3] +2.0*S[2][6][3] );
                S[5][16][3] = gma1*S[2][16][3] + fak * S[2][16][1];
                S[6][16][3] = gma2*S[2][16][3] + fak * S[2][8][3];
                S[8][16][4] = gma0*S[2][16][4] + fak * (S[1][16][4] +2.0*S[2][6][4] );
                S[5][16][4] = gma1*S[2][16][4];
                S[6][16][4] = gma2*S[2][16][4] + fak * (S[2][8][4] +S[2][16][1] );
                S[8][17][2] = gma0*S[2][17][2] + fak*(S[1][17][2]+S[2][10][2]+S[2][17][1]);
                S[5][17][2] = gma1*S[2][17][2];
                S[6][17][2] = gma2*S[2][17][2] + fak2* S[2][6][2];
                S[8][17][3] = gma0*S[2][17][3] + fak * (S[1][17][3] +S[2][10][3] );
                S[5][17][3] = gma1*S[2][17][3] + fak * S[2][17][1];
                S[6][17][3] = gma2*S[2][17][3] + fak2* S[2][6][3];
                S[8][17][4] = gma0*S[2][17][4] + fak * (S[1][17][4] +S[2][10][4] );
                S[5][17][4] = gma1*S[2][17][4];
                S[6][17][4] = gma2*S[2][17][4] + fak * (2.0*S[2][6][4] +S[2][17][1] );
                S[8][18][2] = gma0*S[2][18][2] + fak * (S[1][18][2] +S[2][18][1] );
                S[5][18][2] = gma1*S[2][18][2] + fak2* S[2][7][2];
                S[6][18][2] = gma2*S[2][18][2] + fak * S[2][9][2];
                S[8][18][3] = gma0*S[2][18][3] + fak * S[1][18][3];
                S[5][18][3] = gma1*S[2][18][3] + fak * (2.0*S[2][7][3] +S[2][18][1] );
                S[6][18][3] = gma2*S[2][18][3] + fak * S[2][9][3];
                S[8][18][4] = gma0*S[2][18][4] + fak * S[1][18][4];
                S[5][18][4] = gma1*S[2][18][4] + fak2* S[2][7][4];
                S[6][18][4] = gma2*S[2][18][4] + fak * (S[2][9][4] +S[2][18][1] );
                S[8][19][2] = gma0*S[2][19][2] + fak * (S[1][19][2] +S[2][19][1] );
                S[5][19][2] = gma1*S[2][19][2] + fak * S[2][10][2];
                S[6][19][2] = gma2*S[2][19][2] + fak2* S[2][7][2];
                S[8][19][3] = gma0*S[2][19][3] + fak * S[1][19][3];
                S[5][19][3] = gma1*S[2][19][3] + fak * (S[2][10][3] +S[2][19][1] );
                S[6][19][3] = gma2*S[2][19][3] + fak2* S[2][7][3];
                S[8][19][4] = gma0*S[2][19][4] + fak * S[1][19][4];
                S[5][19][4] = gma1*S[2][19][4] + fak * S[2][10][4];
                S[6][19][4] = gma2*S[2][19][4] + fak * (2.0*S[2][7][4] +S[2][19][1] );
                S[8][20][2] = gma0*S[2][20][2] + fak *(S[1][20][2]+S[2][7][2]+S[2][20][1]);
                S[5][20][2] = gma1*S[2][20][2] + fak * S[2][6][2];
                S[6][20][2] = gma2*S[2][20][2] + fak * S[2][5][2];
                S[8][20][3] = gma0*S[2][20][3] + fak * (S[1][20][3] +S[2][7][3] );
                S[5][20][3] = gma1*S[2][20][3] + fak * (S[2][6][3] +S[2][20][1] );
                S[6][20][3] = gma2*S[2][20][3] + fak * S[2][5][3];
                S[8][20][4] = gma0*S[2][20][4] + fak * (S[1][20][4] +S[2][7][4] );
                S[5][20][4] = gma1*S[2][20][4] + fak * S[2][6][4];
                S[6][20][4] = gma2*S[2][20][4] + fak * (S[2][5][4] +S[2][20][1] );
                S[9][11][2] = gma1*S[3][11][2] + fak * S[1][11][2];
                S[7][11][2] = gma2*S[3][11][2];
                S[9][11][3] = gma1*S[3][11][3] + fak * (S[1][11][3] +S[3][11][1] );
                S[7][11][3] = gma2*S[3][11][3];
                S[9][11][4] = gma1*S[3][11][4] + fak * S[1][11][4];
                S[7][11][4] = gma2*S[3][11][4] + fak * S[3][11][1];
                S[9][12][2] = gma1*S[3][12][2] + fak * (S[1][12][2] +3.0*S[3][9][2] );
                S[7][12][2] = gma2*S[3][12][2];
                S[9][12][3] = gma1*S[3][12][3] + fak * (S[1][12][3] +3.0*S[3][9][3] +S[3][12][1] );
                S[7][12][3] = gma2*S[3][12][3];
                S[9][12][4] = gma1*S[3][12][4] + fak * (S[1][12][4] +3.0*S[3][9][4] );
                S[7][12][4] = gma2*S[3][12][4] + fak * S[3][12][1];
                S[9][13][2] = gma1*S[3][13][2] + fak * S[1][13][2];
                S[7][13][2] = gma2*S[3][13][2] + fak3* S[3][10][2];
                S[9][13][3] = gma1*S[3][13][3] + fak * (S[1][13][3] +S[3][13][1] );
                S[7][13][3] = gma2*S[3][13][3] + fak3* S[3][10][3];
                S[9][13][4] = gma1*S[3][13][4] + fak * S[1][13][4];
                S[7][13][4] = gma2*S[3][13][4] + fak * (3.0*S[3][10][4] +S[3][13][1] );
                S[9][14][2] = gma1*S[3][14][2] + fak * (S[1][14][2] +S[3][8][2] );
                S[7][14][2] = gma2*S[3][14][2];
                S[9][14][3] = gma1*S[3][14][3] + fak *(S[1][14][3]+S[3][8][3]+S[3][14][1]);
                S[7][14][3] = gma2*S[3][14][3];
                S[9][14][4] = gma1*S[3][14][4] + fak * (S[1][14][4] +S[3][8][4] );
                S[7][14][4] = gma2*S[3][14][4] + fak * S[3][14][1];
                S[9][15][2] = gma1*S[3][15][2] + fak * (S[1][15][2] +2.0*S[3][5][2] );
                S[7][15][2] = gma2*S[3][15][2];
                S[9][15][3] = gma1*S[3][15][3] + fak * (S[1][15][3] +2.0*S[3][5][3] +S[3][15][1] );
                S[7][15][3] = gma2*S[3][15][3];
                S[9][15][4] = gma1*S[3][15][4] + fak * (S[1][15][4] +2.0*S[3][5][4] );
                S[7][15][4] = gma2*S[3][15][4] + fak * S[3][15][1];
                S[9][16][2] = gma1*S[3][16][2] + fak * S[1][16][2];
                S[7][16][2] = gma2*S[3][16][2] + fak * S[3][8][2];
                S[9][16][3] = gma1*S[3][16][3] + fak * (S[1][16][3] +S[3][16][1] );
                S[7][16][3] = gma2*S[3][16][3] + fak * S[3][8][3];
                S[9][16][4] = gma1*S[3][16][4] + fak * S[1][16][4];
                S[7][16][4] = gma2*S[3][16][4] + fak * (S[3][8][4] +S[3][16][1] );
                S[9][17][2] = gma1*S[3][17][2] + fak * S[1][17][2];
                S[7][17][2] = gma2*S[3][17][2] + fak2* S[3][6][2];
                S[9][17][3] = gma1*S[3][17][3] + fak * (S[1][17][3] +S[3][17][1] );
                S[7][17][3] = gma2*S[3][17][3] + fak2* S[3][6][3];
                S[9][17][4] = gma1*S[3][17][4] + fak * S[1][17][4];
                S[7][17][4] = gma2*S[3][17][4] + fak * (2.0*S[3][6][4] +S[3][17][1] );
                S[9][18][2] = gma1*S[3][18][2] + fak * (S[1][18][2] +2.0*S[3][7][2] );
                S[7][18][2] = gma2*S[3][18][2] + fak * S[3][9][2];
                S[9][18][3] = gma1*S[3][18][3] + fak * (S[1][18][3] +2.0*S[3][7][3] +S[3][18][1] );
                S[7][18][3] = gma2*S[3][18][3] + fak * S[3][9][3];
                S[9][18][4] = gma1*S[3][18][4] + fak * (S[1][18][4] +2.0*S[3][7][4] );
                S[7][18][4] = gma2*S[3][18][4] + fak * (S[3][9][4] +S[3][18][1] );
                S[9][19][2] = gma1*S[3][19][2] + fak * (S[1][19][2] +S[3][10][2] );
                S[7][19][2] = gma2*S[3][19][2] + fak2* S[3][7][2];
                S[9][19][3] = gma1*S[3][19][3] + fak*(S[1][19][3]+S[3][10][3]+S[3][19][1]);
                S[7][19][3] = gma2*S[3][19][3] + fak2* S[3][7][3];
                S[9][19][4] = gma1*S[3][19][4] + fak * (S[1][19][4] +S[3][10][4] );
                S[7][19][4] = gma2*S[3][19][4] + fak * (2.0*S[3][7][4] +S[3][19][1] );
                S[9][20][2] = gma1*S[3][20][2] + fak * (S[1][20][2] +S[3][6][2] );
                S[7][20][2] = gma2*S[3][20][2] + fak * S[3][5][2];
                S[9][20][3] = gma1*S[3][20][3] + fak *(S[1][20][3]+S[3][6][3]+S[3][20][1]);
                S[7][20][3] = gma2*S[3][20][3] + fak * S[3][5][3];
                S[9][20][4] = gma1*S[3][20][4] + fak * (S[1][20][4] +S[3][6][4] );
                S[7][20][4] = gma2*S[3][20][4] + fak * (S[3][5][4] +S[3][20][1] );
                S[10][11][2] = gma2*S[4][11][2] + fak * S[1][11][2];
                S[10][11][3] = gma2*S[4][11][3] + fak * S[1][11][3];
                S[10][11][4] = gma2*S[4][11][4] + fak * (S[1][11][4] +S[4][11][1] );
                S[10][12][2] = gma2*S[4][12][2] + fak * S[1][12][2];
                S[10][12][3] = gma2*S[4][12][3] + fak * S[1][12][3];
                S[10][12][4] = gma2*S[4][12][4] + fak * (S[1][12][4] +S[4][12][1] );
                S[10][13][2] = gma2*S[4][13][2] + fak * (S[1][13][2] +3.0*S[4][10][2] );
                S[10][13][3] = gma2*S[4][13][3] + fak * (S[1][13][3] +3.0*S[4][10][3] );
                S[10][13][4] = gma2*S[4][13][4] + fak * (S[1][13][4] +3.0*S[4][10][4] +S[4][13][1] );
                S[10][14][2] = gma2*S[4][14][2] + fak * S[1][14][2];
                S[10][14][3] = gma2*S[4][14][3] + fak * S[1][14][3];
                S[10][14][4] = gma2*S[4][14][4] + fak * (S[1][14][4] +S[4][14][1] );
                S[10][15][2] = gma2*S[4][15][2] + fak * S[1][15][2];
                S[10][15][3] = gma2*S[4][15][3] + fak * S[1][15][3];
                S[10][15][4] = gma2*S[4][15][4] + fak * (S[1][15][4] +S[4][15][1] );
                S[10][16][2] = gma2*S[4][16][2] + fak * (S[1][16][2] +S[4][8][2] );
                S[10][16][3] = gma2*S[4][16][3] + fak * (S[1][16][3] +S[4][8][3] );
                S[10][16][4] = gma2*S[4][16][4] + fak*(S[1][16][4]+S[4][8][4]+S[4][16][1]);
                S[10][17][2] = gma2*S[4][17][2] + fak * (S[1][17][2] +2.0*S[4][6][2] );
                S[10][17][3] = gma2*S[4][17][3] + fak * (S[1][17][3] +2.0*S[4][6][3] );
                S[10][17][4] = gma2*S[4][17][4] + fak * (S[1][17][4] +2.0*S[4][6][4] +S[4][17][1] );
                S[10][18][2] = gma2*S[4][18][2] + fak * (S[1][18][2] +S[4][9][2] );
                S[10][18][3] = gma2*S[4][18][3] + fak * (S[1][18][3] +S[4][9][3] );
                S[10][18][4] = gma2*S[4][18][4] + fak*(S[1][18][4]+S[4][9][4]+S[4][18][1]);
                S[10][19][2] = gma2*S[4][19][2] + fak * (S[1][19][2] +2.0*S[4][7][2] );
                S[10][19][3] = gma2*S[4][19][3] + fak * (S[1][19][3] +2.0*S[4][7][3] );
                S[10][19][4] = gma2*S[4][19][4] + fak * (S[1][19][4] +2.0*S[4][7][4] +S[4][19][1] );
                S[10][20][2] = gma2*S[4][20][2] + fak * (S[1][20][2] +S[4][5][2] );
                S[10][20][3] = gma2*S[4][20][3] + fak * (S[1][20][3] +S[4][5][3] );
                S[10][20][4] = gma2*S[4][20][4] + fak*(S[1][20][4]+S[4][5][4]+S[4][20][1]);
             }

            // d-f-d elements
             if ( _lmax_alpha >= 2 && _lmax_gw >= 3 && _lmax_gamma >= 2) {
                S[8][11][5] = gma0*S[2][11][5] + fak * (S[1][11][5] +3.0*S[2][8][5] +S[2][11][3] );
                S[5][11][5] = gma1*S[2][11][5] + fak * S[2][11][2];
                S[6][11][5] = gma2*S[2][11][5];
                S[8][11][6] = gma0*S[2][11][6] + fak * (S[1][11][6] +3.0*S[2][8][6] +S[2][11][4] );
                S[5][11][6] = gma1*S[2][11][6];
                S[6][11][6] = gma2*S[2][11][6] + fak * S[2][11][2];
                S[8][11][7] = gma0*S[2][11][7] + fak * (S[1][11][7] +3.0*S[2][8][7] );
                S[5][11][7] = gma1*S[2][11][7] + fak * S[2][11][4];
                S[6][11][7] = gma2*S[2][11][7] + fak * S[2][11][3];
                S[8][11][8] = gma0*S[2][11][8] + fak * (S[1][11][8] +3.0*S[2][8][8] +2.0*S[2][11][2] );
                S[5][11][8] = gma1*S[2][11][8];
                S[6][11][8] = gma2*S[2][11][8];
                S[8][11][9] = gma0*S[2][11][9] + fak * (S[1][11][9] +3.0*S[2][8][9] );
                S[5][11][9] = gma1*S[2][11][9] + fak2* S[2][11][3];
                S[6][11][9] = gma2*S[2][11][9];
                S[8][11][10] = gma0*S[2][11][10] + fak * (S[1][11][10]+3.0*S[2][8][10]);
                S[5][11][10] = gma1*S[2][11][10];
                S[6][11][10] = gma2*S[2][11][10] + fak2* S[2][11][4];
                S[8][12][5] = gma0*S[2][12][5] + fak * (S[1][12][5] +S[2][12][3] );
                S[5][12][5] = gma1*S[2][12][5] + fak * (3.0*S[2][9][5] +S[2][12][2] );
                S[6][12][5] = gma2*S[2][12][5];
                S[8][12][6] = gma0*S[2][12][6] + fak * (S[1][12][6] +S[2][12][4] );
                S[5][12][6] = gma1*S[2][12][6] + fak3* S[2][9][6];
                S[6][12][6] = gma2*S[2][12][6] + fak * S[2][12][2];
                S[8][12][7] = gma0*S[2][12][7] + fak * S[1][12][7];
                S[5][12][7] = gma1*S[2][12][7] + fak * (3.0*S[2][9][7] +S[2][12][4] );
                S[6][12][7] = gma2*S[2][12][7] + fak * S[2][12][3];
                S[8][12][8] = gma0*S[2][12][8] + fak * (S[1][12][8] +2.0*S[2][12][2] );
                S[5][12][8] = gma1*S[2][12][8] + fak3* S[2][9][8];
                S[6][12][8] = gma2*S[2][12][8];
                S[8][12][9] = gma0*S[2][12][9] + fak * S[1][12][9];
                S[5][12][9] = gma1*S[2][12][9] + fak *(3.0*S[2][9][9]+2.0*S[2][12][3]);
                S[6][12][9] = gma2*S[2][12][9];
                S[8][12][10] = gma0*S[2][12][10] + fak * S[1][12][10];
                S[5][12][10] = gma1*S[2][12][10] + fak3* S[2][9][10];
                S[6][12][10] = gma2*S[2][12][10] + fak2* S[2][12][4];
                S[8][13][5] = gma0*S[2][13][5] + fak * (S[1][13][5] +S[2][13][3] );
                S[5][13][5] = gma1*S[2][13][5] + fak * S[2][13][2];
                S[6][13][5] = gma2*S[2][13][5] + fak3* S[2][10][5];
                S[8][13][6] = gma0*S[2][13][6] + fak * (S[1][13][6] +S[2][13][4] );
                S[5][13][6] = gma1*S[2][13][6];
                S[6][13][6] = gma2*S[2][13][6] + fak * (3.0*S[2][10][6] +S[2][13][2] );
                S[8][13][7] = gma0*S[2][13][7] + fak * S[1][13][7];
                S[5][13][7] = gma1*S[2][13][7] + fak * S[2][13][4];
                S[6][13][7] = gma2*S[2][13][7] + fak * (3.0*S[2][10][7] +S[2][13][3] );
                S[8][13][8] = gma0*S[2][13][8] + fak * (S[1][13][8] +2.0*S[2][13][2] );
                S[5][13][8] = gma1*S[2][13][8];
                S[6][13][8] = gma2*S[2][13][8] + fak3* S[2][10][8];
                S[8][13][9] = gma0*S[2][13][9] + fak * S[1][13][9];
                S[5][13][9] = gma1*S[2][13][9] + fak2* S[2][13][3];
                S[6][13][9] = gma2*S[2][13][9] + fak3* S[2][10][9];
                S[8][13][10] = gma0*S[2][13][10] + fak * S[1][13][10];
                S[5][13][10] = gma1*S[2][13][10];
                S[6][13][10] = gma2*S[2][13][10] + fak * (3.0*S[2][10][10] +2.0*S[2][13][4] );
                S[8][14][5] = gma0*S[2][14][5] +  fak * (S[1][14][5] +2.0*S[2][5][5] +S[2][14][3] );
                S[5][14][5] = gma1*S[2][14][5] + fak * (S[2][8][5] +S[2][14][2] );
                S[6][14][5] = gma2*S[2][14][5];
                S[8][14][6] = gma0*S[2][14][6] + fak * (S[1][14][6] +2.0*S[2][5][6] +S[2][14][4] );
                S[5][14][6] = gma1*S[2][14][6] + fak * S[2][8][6];
                S[6][14][6] = gma2*S[2][14][6] + fak * S[2][14][2];
                S[8][14][7] = gma0*S[2][14][7] + fak * (S[1][14][7] +2.0*S[2][5][7] );
                S[5][14][7] = gma1*S[2][14][7] + fak * (S[2][8][7] +S[2][14][4] );
                S[6][14][7] = gma2*S[2][14][7] + fak * S[2][14][3];
                S[8][14][8] = gma0*S[2][14][8] + fak * (S[1][14][8] +2.0*S[2][5][8] +2.0*S[2][14][2] );
                S[5][14][8] = gma1*S[2][14][8] + fak * S[2][8][8];
                S[6][14][8] = gma2*S[2][14][8];
                S[8][14][9] = gma0*S[2][14][9] + fak * (S[1][14][9] +2.0*S[2][5][9] );
                S[5][14][9] = gma1*S[2][14][9] + fak * (S[2][8][9] +2.0*S[2][14][3] );
                S[6][14][9] = gma2*S[2][14][9];
                S[8][14][10] = gma0*S[2][14][10] + fak * (S[1][14][10]+2.0*S[2][5][10]);
                S[5][14][10] = gma1*S[2][14][10] + fak * S[2][8][10];
                S[6][14][10] = gma2*S[2][14][10] + fak2* S[2][14][4];
                S[8][15][5] = gma0*S[2][15][5] + fak *(S[1][15][5]+S[2][9][5]+S[2][15][3]);
                S[5][15][5] = gma1*S[2][15][5] + fak * (2.0*S[2][5][5] +S[2][15][2] );
                S[6][15][5] = gma2*S[2][15][5];
                S[8][15][6] = gma0*S[2][15][6] + fak *(S[1][15][6]+S[2][9][6]+S[2][15][4]);
                S[5][15][6] = gma1*S[2][15][6] + fak2* S[2][5][6];
                S[6][15][6] = gma2*S[2][15][6] + fak * S[2][15][2];
                S[8][15][7] = gma0*S[2][15][7] + fak * (S[1][15][7] +S[2][9][7] );
                S[5][15][7] = gma1*S[2][15][7] + fak * (2.0*S[2][5][7] +S[2][15][4] );
                S[6][15][7] = gma2*S[2][15][7] + fak * S[2][15][3];
                S[8][15][8] = gma0*S[2][15][8] + fak * (S[1][15][8] +S[2][9][8] +2.0*S[2][15][2] );
                S[5][15][8] = gma1*S[2][15][8] + fak2* S[2][5][8];
                S[6][15][8] = gma2*S[2][15][8];
                S[8][15][9] = gma0*S[2][15][9] + fak * (S[1][15][9] +S[2][9][9] );
                S[5][15][9] = gma1*S[2][15][9] + fak2* (S[2][5][9] +S[2][15][3] );
                S[6][15][9] = gma2*S[2][15][9];
                S[8][15][10] = gma0*S[2][15][10] + fak * (S[1][15][10] +S[2][9][10] );
                S[5][15][10] = gma1*S[2][15][10] + fak2* S[2][5][10];
                S[6][15][10] = gma2*S[2][15][10] + fak2* S[2][15][4];
                S[8][16][5] = gma0*S[2][16][5] + fak * (S[1][16][5] +2.0*S[2][6][5] +S[2][16][3] );
                S[5][16][5] = gma1*S[2][16][5] + fak * S[2][16][2];
                S[6][16][5] = gma2*S[2][16][5] + fak * S[2][8][5];
                S[8][16][6] = gma0*S[2][16][6] + fak * (S[1][16][6] +2.0*S[2][6][6] +S[2][16][4] );
                S[5][16][6] = gma1*S[2][16][6];
                S[6][16][6] = gma2*S[2][16][6] + fak * (S[2][8][6] +S[2][16][2] );
                S[8][16][7] = gma0*S[2][16][7] + fak * (S[1][16][7] +2.0*S[2][6][7] );
                S[5][16][7] = gma1*S[2][16][7] + fak * S[2][16][4];
                S[6][16][7] = gma2*S[2][16][7] + fak * (S[2][8][7] +S[2][16][3] );
                S[8][16][8] = gma0*S[2][16][8] + fak * (S[1][16][8] +2.0*S[2][6][8] +2.0*S[2][16][2] );
                S[5][16][8] = gma1*S[2][16][8];
                S[6][16][8] = gma2*S[2][16][8] + fak * S[2][8][8];
                S[8][16][9] = gma0*S[2][16][9] + fak * (S[1][16][9] +2.0*S[2][6][9] );
                S[5][16][9] = gma1*S[2][16][9] + fak2* S[2][16][3];
                S[6][16][9] = gma2*S[2][16][9] + fak * S[2][8][9];
                S[8][16][10] = gma0*S[2][16][10] + fak * (S[1][16][10]+2.0*S[2][6][10]);
                S[5][16][10] = gma1*S[2][16][10];
                S[6][16][10] = gma2*S[2][16][10] + fak * (S[2][8][10] +2.0*S[2][16][4]);
                S[8][17][5] = gma0*S[2][17][5] + fak*(S[1][17][5]+S[2][10][5]+S[2][17][3]);
                S[5][17][5] = gma1*S[2][17][5] + fak * S[2][17][2];
                S[6][17][5] = gma2*S[2][17][5] + fak2* S[2][6][5];
                S[8][17][6] = gma0*S[2][17][6] + fak*(S[1][17][6]+S[2][10][6]+S[2][17][4]);
                S[5][17][6] = gma1*S[2][17][6];
                S[6][17][6] = gma2*S[2][17][6] + fak * (2.0*S[2][6][6] +S[2][17][2] );
                S[8][17][7] = gma0*S[2][17][7] + fak * (S[1][17][7] +S[2][10][7] );
                S[5][17][7] = gma1*S[2][17][7] + fak * S[2][17][4];
                S[6][17][7] = gma2*S[2][17][7] + fak * (2.0*S[2][6][7] +S[2][17][3] );
                S[8][17][8] = gma0*S[2][17][8] + fak * (S[1][17][8] +S[2][10][8] +2.0*S[2][17][2] );
                S[5][17][8] = gma1*S[2][17][8];
                S[6][17][8] = gma2*S[2][17][8] + fak2* S[2][6][8];
                S[8][17][9] = gma0*S[2][17][9] + fak * (S[1][17][9] +S[2][10][9] );
                S[5][17][9] = gma1*S[2][17][9] + fak2* S[2][17][3];
                S[6][17][9] = gma2*S[2][17][9] + fak2* S[2][6][9];
                S[8][17][10] = gma0*S[2][17][10] + fak * (S[1][17][10] +S[2][10][10] );
                S[5][17][10] = gma1*S[2][17][10];
                S[6][17][10] = gma2*S[2][17][10] + fak2* (S[2][6][10] +S[2][17][4] );
                S[8][18][5] = gma0*S[2][18][5] + fak * (S[1][18][5] +S[2][18][3] );
                S[5][18][5] = gma1*S[2][18][5] + fak * (2.0*S[2][7][5] +S[2][18][2] );
                S[6][18][5] = gma2*S[2][18][5] + fak * S[2][9][5];
                S[8][18][6] = gma0*S[2][18][6] + fak * (S[1][18][6] +S[2][18][4] );
                S[5][18][6] = gma1*S[2][18][6] + fak2* S[2][7][6];
                S[6][18][6] = gma2*S[2][18][6] + fak * (S[2][9][6] +S[2][18][2] );
                S[8][18][7] = gma0*S[2][18][7] + fak * S[1][18][7];
                S[5][18][7] = gma1*S[2][18][7] + fak * (2.0*S[2][7][7] +S[2][18][4] );
                S[6][18][7] = gma2*S[2][18][7] + fak * (S[2][9][7] +S[2][18][3] );
                S[8][18][8] = gma0*S[2][18][8] + fak * (S[1][18][8] +2.0*S[2][18][2] );
                S[5][18][8] = gma1*S[2][18][8] + fak2* S[2][7][8];
                S[6][18][8] = gma2*S[2][18][8] + fak * S[2][9][8];
                S[8][18][9] = gma0*S[2][18][9] + fak * S[1][18][9];
                S[5][18][9] = gma1*S[2][18][9] + fak2* (S[2][7][9] +S[2][18][3] );
                S[6][18][9] = gma2*S[2][18][9] + fak * S[2][9][9];
                S[8][18][10] = gma0*S[2][18][10] + fak * S[1][18][10];
                S[5][18][10] = gma1*S[2][18][10] + fak2* S[2][7][10];
                S[6][18][10] = gma2*S[2][18][10] + fak * (S[2][9][10] +2.0*S[2][18][4]);
                S[8][19][5] = gma0*S[2][19][5] + fak * (S[1][19][5] +S[2][19][3] );
                S[5][19][5] = gma1*S[2][19][5] + fak * (S[2][10][5] +S[2][19][2] );
                S[6][19][5] = gma2*S[2][19][5] + fak2* S[2][7][5];
                S[8][19][6] = gma0*S[2][19][6] + fak * (S[1][19][6] +S[2][19][4] );
                S[5][19][6] = gma1*S[2][19][6] + fak * S[2][10][6];
                S[6][19][6] = gma2*S[2][19][6] + fak * (2.0*S[2][7][6] +S[2][19][2] );
                S[8][19][7] = gma0*S[2][19][7] + fak * S[1][19][7];
                S[5][19][7] = gma1*S[2][19][7] + fak * (S[2][10][7] +S[2][19][4] );
                S[6][19][7] = gma2*S[2][19][7] + fak * (2.0*S[2][7][7] +S[2][19][3] );
                S[8][19][8] = gma0*S[2][19][8] + fak * (S[1][19][8] +2.0*S[2][19][2] );
                S[5][19][8] = gma1*S[2][19][8] + fak * S[2][10][8];
                S[6][19][8] = gma2*S[2][19][8] + fak2* S[2][7][8];
                S[8][19][9] = gma0*S[2][19][9] + fak * S[1][19][9];
                S[5][19][9] = gma1*S[2][19][9] + fak * (S[2][10][9] +2.0*S[2][19][3] );
                S[6][19][9] = gma2*S[2][19][9] + fak2* S[2][7][9];
                S[8][19][10] = gma0*S[2][19][10] + fak * S[1][19][10];
                S[5][19][10] = gma1*S[2][19][10] + fak * S[2][10][10];
                S[6][19][10] = gma2*S[2][19][10] + fak2* (S[2][7][10] +S[2][19][4] );
                S[8][20][5] = gma0*S[2][20][5] + fak *(S[1][20][5]+S[2][7][5]+S[2][20][3]);
                S[5][20][5] = gma1*S[2][20][5] + fak * (S[2][6][5] +S[2][20][2] );
                S[6][20][5] = gma2*S[2][20][5] + fak * S[2][5][5];
                S[8][20][6] = gma0*S[2][20][6] + fak *(S[1][20][6]+S[2][7][6]+S[2][20][4]);
                S[5][20][6] = gma1*S[2][20][6] + fak * S[2][6][6];
                S[6][20][6] = gma2*S[2][20][6] + fak * (S[2][5][6] +S[2][20][2] );
                S[8][20][7] = gma0*S[2][20][7] + fak * (S[1][20][7] +S[2][7][7] );
                S[5][20][7] = gma1*S[2][20][7] + fak * (S[2][6][7] +S[2][20][4] );
                S[6][20][7] = gma2*S[2][20][7] + fak * (S[2][5][7] +S[2][20][3] );
                S[8][20][8] = gma0*S[2][20][8] + fak * (S[1][20][8] +S[2][7][8] +2.0*S[2][20][2] );
                S[5][20][8] = gma1*S[2][20][8] + fak * S[2][6][8];
                S[6][20][8] = gma2*S[2][20][8] + fak * S[2][5][8];
                S[8][20][9] = gma0*S[2][20][9] + fak * (S[1][20][9] +S[2][7][9] );
                S[5][20][9] = gma1*S[2][20][9] + fak * (S[2][6][9] +2.0*S[2][20][3] );
                S[6][20][9] = gma2*S[2][20][9] + fak * S[2][5][9];
                S[8][20][10] = gma0*S[2][20][10] + fak * (S[1][20][10] +S[2][7][10] );
                S[5][20][10] = gma1*S[2][20][10] + fak * S[2][6][10];
                S[6][20][10] = gma2*S[2][20][10] + fak * (S[2][5][10] +2.0*S[2][20][4]);
                S[9][11][5] = gma1*S[3][11][5] + fak * (S[1][11][5] +S[3][11][2] );
                S[7][11][5] = gma2*S[3][11][5];
                S[9][11][6] = gma1*S[3][11][6] + fak * S[1][11][6];
                S[7][11][6] = gma2*S[3][11][6] + fak * S[3][11][2];
                S[9][11][7] = gma1*S[3][11][7] + fak * (S[1][11][7] +S[3][11][4] );
                S[7][11][7] = gma2*S[3][11][7] + fak * S[3][11][3];
                S[9][11][8] = gma1*S[3][11][8] + fak * S[1][11][8];
                S[7][11][8] = gma2*S[3][11][8];
                S[9][11][9] = gma1*S[3][11][9] + fak * (S[1][11][9] +2.0*S[3][11][3] );
                S[7][11][9] = gma2*S[3][11][9];
                S[9][11][10] = gma1*S[3][11][10] + fak * S[1][11][10];
                S[7][11][10] = gma2*S[3][11][10] + fak2* S[3][11][4];
                S[9][12][5] = gma1*S[3][12][5]   + fak * (S[1][12][5] +3.0*S[3][9][5] +S[3][12][2] );
                S[7][12][5] = gma2*S[3][12][5];
                S[9][12][6] = gma1*S[3][12][6] + fak * (S[1][12][6] +3.0*S[3][9][6] );
                S[7][12][6] = gma2*S[3][12][6] + fak * S[3][12][2];
                S[9][12][7] = gma1*S[3][12][7] + fak * (S[1][12][7] +3.0*S[3][9][7] +S[3][12][4] );
                S[7][12][7] = gma2*S[3][12][7] + fak * S[3][12][3];
                S[9][12][8] = gma1*S[3][12][8] + fak * (S[1][12][8] +3.0*S[3][9][8] );
                S[7][12][8] = gma2*S[3][12][8];
                S[9][12][9] = gma1*S[3][12][9] + fak * (S[1][12][9] +3.0*S[3][9][9] +2.0*S[3][12][3] );
                S[7][12][9] = gma2*S[3][12][9];
                S[9][12][10] = gma1*S[3][12][10] + fak * (S[1][12][10]+3.0*S[3][9][10]);
                S[7][12][10] = gma2*S[3][12][10] + fak2* S[3][12][4];
                S[9][13][5] = gma1*S[3][13][5] + fak * (S[1][13][5] +S[3][13][2] );
                S[7][13][5] = gma2*S[3][13][5] + fak3* S[3][10][5];
                S[9][13][6] = gma1*S[3][13][6] + fak * S[1][13][6];
                S[7][13][6] = gma2*S[3][13][6] + fak * (3.0*S[3][10][6] +S[3][13][2] );
                S[9][13][7] = gma1*S[3][13][7] + fak * (S[1][13][7] +S[3][13][4] );
                S[7][13][7] = gma2*S[3][13][7] + fak * (3.0*S[3][10][7] +S[3][13][3] );
                S[9][13][8] = gma1*S[3][13][8] + fak * S[1][13][8];
                S[7][13][8] = gma2*S[3][13][8] + fak3* S[3][10][8];
                S[9][13][9] = gma1*S[3][13][9] + fak * (S[1][13][9] +2.0*S[3][13][3] );
                S[7][13][9] = gma2*S[3][13][9] + fak3* S[3][10][9];
                S[9][13][10] = gma1*S[3][13][10] + fak * S[1][13][10];
                S[7][13][10] = gma2*S[3][13][10] + fak * (3.0*S[3][10][10] +2.0*S[3][13][4] );
                S[9][14][5] = gma1*S[3][14][5] + fak *(S[1][14][5]+S[3][8][5]+S[3][14][2]);
                S[7][14][5] = gma2*S[3][14][5];
                S[9][14][6] = gma1*S[3][14][6] + fak * (S[1][14][6] +S[3][8][6] );
                S[7][14][6] = gma2*S[3][14][6] + fak * S[3][14][2];
                S[9][14][7] = gma1*S[3][14][7] + fak *(S[1][14][7]+S[3][8][7]+S[3][14][4]);
                S[7][14][7] = gma2*S[3][14][7] + fak * S[3][14][3];
                S[9][14][8] = gma1*S[3][14][8] + fak * (S[1][14][8] +S[3][8][8] );
                S[7][14][8] = gma2*S[3][14][8];
                S[9][14][9] = gma1*S[3][14][9] + fak * (S[1][14][9] +S[3][8][9] +2.0*S[3][14][3] );
                S[7][14][9] = gma2*S[3][14][9];
                S[9][14][10] = gma1*S[3][14][10] + fak * (S[1][14][10] +S[3][8][10] );
                S[7][14][10] = gma2*S[3][14][10] + fak2* S[3][14][4];
                S[9][15][5] = gma1*S[3][15][5] + fak * (S[1][15][5] +2.0*S[3][5][5] +S[3][15][2] );
                S[7][15][5] = gma2*S[3][15][5];
                S[9][15][6] = gma1*S[3][15][6] + fak * (S[1][15][6] +2.0*S[3][5][6] );
                S[7][15][6] = gma2*S[3][15][6] + fak * S[3][15][2];
                S[9][15][7] = gma1*S[3][15][7] + fak * (S[1][15][7] +2.0*S[3][5][7] +S[3][15][4] );
                S[7][15][7] = gma2*S[3][15][7] + fak * S[3][15][3];
                S[9][15][8] = gma1*S[3][15][8] + fak * (S[1][15][8] +2.0*S[3][5][8] );
                S[7][15][8] = gma2*S[3][15][8];
                S[9][15][9] = gma1*S[3][15][9] + fak * (S[1][15][9] +2.0*S[3][5][9] +2.0*S[3][15][3] );
                S[7][15][9] = gma2*S[3][15][9];
                S[9][15][10] = gma1*S[3][15][10] + fak * (S[1][15][10]+2.0*S[3][5][10]);
                S[7][15][10] = gma2*S[3][15][10] + fak2* S[3][15][4];
                S[9][16][5] = gma1*S[3][16][5] + fak * (S[1][16][5] +S[3][16][2] );
                S[7][16][5] = gma2*S[3][16][5] + fak * S[3][8][5];
                S[9][16][6] = gma1*S[3][16][6] + fak * S[1][16][6];
                S[7][16][6] = gma2*S[3][16][6] + fak * (S[3][8][6] +S[3][16][2] );
                S[9][16][7] = gma1*S[3][16][7] + fak * (S[1][16][7] +S[3][16][4] );
                S[7][16][7] = gma2*S[3][16][7] + fak * (S[3][8][7] +S[3][16][3] );
                S[9][16][8] = gma1*S[3][16][8] + fak * S[1][16][8];
                S[7][16][8] = gma2*S[3][16][8] + fak * S[3][8][8];
                S[9][16][9] = gma1*S[3][16][9] + fak * (S[1][16][9] +2.0*S[3][16][3] );
                S[7][16][9] = gma2*S[3][16][9] + fak * S[3][8][9];
                S[9][16][10] = gma1*S[3][16][10] + fak * S[1][16][10];
                S[7][16][10] = gma2*S[3][16][10] + fak * (S[3][8][10] +2.0*S[3][16][4]);
                S[9][17][5] = gma1*S[3][17][5] + fak * (S[1][17][5] +S[3][17][2] );
                S[7][17][5] = gma2*S[3][17][5] + fak2* S[3][6][5];
                S[9][17][6] = gma1*S[3][17][6] + fak * S[1][17][6];
                S[7][17][6] = gma2*S[3][17][6] + fak * (2.0*S[3][6][6] +S[3][17][2] );
                S[9][17][7] = gma1*S[3][17][7] + fak * (S[1][17][7] +S[3][17][4] );
                S[7][17][7] = gma2*S[3][17][7] + fak * (2.0*S[3][6][7] +S[3][17][3] );
                S[9][17][8] = gma1*S[3][17][8] + fak * S[1][17][8];
                S[7][17][8] = gma2*S[3][17][8] + fak2* S[3][6][8];
                S[9][17][9] = gma1*S[3][17][9] + fak * (S[1][17][9] +2.0*S[3][17][3] );
                S[7][17][9] = gma2*S[3][17][9] + fak2* S[3][6][9];
                S[9][17][10] = gma1*S[3][17][10] + fak * S[1][17][10];
                S[7][17][10] = gma2*S[3][17][10] + fak2* (S[3][6][10] +S[3][17][4] );
                S[9][18][5] = gma1*S[3][18][5]  + fak * (S[1][18][5] +2.0*S[3][7][5] +S[3][18][2] );
                S[7][18][5] = gma2*S[3][18][5] + fak * S[3][9][5];
                S[9][18][6] = gma1*S[3][18][6] + fak * (S[1][18][6] +2.0*S[3][7][6] );
                S[7][18][6] = gma2*S[3][18][6] + fak * (S[3][9][6] +S[3][18][2] );
                S[9][18][7] = gma1*S[3][18][7] + fak * (S[1][18][7] +2.0*S[3][7][7] +S[3][18][4] );
                S[7][18][7] = gma2*S[3][18][7] + fak * (S[3][9][7] +S[3][18][3] );
                S[9][18][8] = gma1*S[3][18][8] + fak * (S[1][18][8] +2.0*S[3][7][8] );
                S[7][18][8] = gma2*S[3][18][8] + fak * S[3][9][8];
                S[9][18][9] = gma1*S[3][18][9] + fak * (S[1][18][9] +2.0*S[3][7][9] +2.0*S[3][18][3] );
                S[7][18][9] = gma2*S[3][18][9] + fak * S[3][9][9];
                S[9][18][10] = gma1*S[3][18][10] + fak * (S[1][18][10]+2.0*S[3][7][10]);
                S[7][18][10] = gma2*S[3][18][10] + fak * (S[3][9][10] +2.0*S[3][18][4]);
                S[9][19][5] = gma1*S[3][19][5] + fak*(S[1][19][5]+S[3][10][5]+S[3][19][2]);
                S[7][19][5] = gma2*S[3][19][5] + fak2* S[3][7][5];
                S[9][19][6] = gma1*S[3][19][6] + fak * (S[1][19][6] +S[3][10][6] );
                S[7][19][6] = gma2*S[3][19][6] + fak * (2.0*S[3][7][6] +S[3][19][2] );
                S[9][19][7] = gma1*S[3][19][7] + fak*(S[1][19][7]+S[3][10][7]+S[3][19][4]);
                S[7][19][7] = gma2*S[3][19][7] + fak * (2.0*S[3][7][7] +S[3][19][3] );
                S[9][19][8] = gma1*S[3][19][8] + fak * (S[1][19][8] +S[3][10][8] );
                S[7][19][8] = gma2*S[3][19][8] + fak2* S[3][7][8];
                S[9][19][9] = gma1*S[3][19][9] + fak * (S[1][19][9] +S[3][10][9] +2.0*S[3][19][3] );
                S[7][19][9] = gma2*S[3][19][9] + fak2* S[3][7][9];
                S[9][19][10] = gma1*S[3][19][10] + fak * (S[1][19][10] +S[3][10][10] );
                S[7][19][10] = gma2*S[3][19][10] + fak2* (S[3][7][10] +S[3][19][4] );
                S[9][20][5] = gma1*S[3][20][5] + fak * (S[1][20][5] +S[3][6][5] + S[3][20][2] );
                S[7][20][5] = gma2*S[3][20][5] + fak * S[3][5][5];
                S[9][20][6] = gma1*S[3][20][6] + fak * (S[1][20][6] +S[3][6][6] );
                S[7][20][6] = gma2*S[3][20][6] + fak * (S[3][5][6] +S[3][20][2] );
                S[9][20][7] = gma1*S[3][20][7] + fak *(S[1][20][7]+S[3][6][7]+S[3][20][4]);
                S[7][20][7] = gma2*S[3][20][7] + fak * (S[3][5][7] +S[3][20][3] );
                S[9][20][8] = gma1*S[3][20][8] + fak * (S[1][20][8] +S[3][6][8] );
                S[7][20][8] = gma2*S[3][20][8] + fak * S[3][5][8];
                S[9][20][9] = gma1*S[3][20][9] + fak * (S[1][20][9] +S[3][6][9] +2.0*S[3][20][3] );
                S[7][20][9] = gma2*S[3][20][9] + fak * S[3][5][9];
                S[9][20][10] = gma1*S[3][20][10] + fak * (S[1][20][10] +S[3][6][10] );
                S[7][20][10] = gma2*S[3][20][10] + fak * (S[3][5][10] +2.0*S[3][20][4]);
                S[10][11][5] = gma2*S[4][11][5] + fak * S[1][11][5];
                S[10][11][6] = gma2*S[4][11][6] + fak * (S[1][11][6] +S[4][11][2] );
                S[10][11][7] = gma2*S[4][11][7] + fak * (S[1][11][7] +S[4][11][3] );
                S[10][11][8] = gma2*S[4][11][8] + fak * S[1][11][8];
                S[10][11][9] = gma2*S[4][11][9] + fak * S[1][11][9];
                S[10][11][10] = gma2*S[4][11][10] + fak *(S[1][11][10]+2.0*S[4][11][4]);
                S[10][12][5] = gma2*S[4][12][5] + fak * S[1][12][5];
                S[10][12][6] = gma2*S[4][12][6] + fak * (S[1][12][6] +S[4][12][2] );
                S[10][12][7] = gma2*S[4][12][7] + fak * (S[1][12][7] +S[4][12][3] );
                S[10][12][8] = gma2*S[4][12][8] + fak * S[1][12][8];
                S[10][12][9] = gma2*S[4][12][9] + fak * S[1][12][9];
                S[10][12][10] = gma2*S[4][12][10] + fak *(S[1][12][10]+2.0*S[4][12][4]);
                S[10][13][5] = gma2*S[4][13][5] + fak * (S[1][13][5] +3.0*S[4][10][5] );
                S[10][13][6] = gma2*S[4][13][6] + fak * (S[1][13][6] +3.0*S[4][10][6] +S[4][13][2] );
                S[10][13][7] = gma2*S[4][13][7] + fak * (S[1][13][7] +3.0*S[4][10][7] +S[4][13][3] );
                S[10][13][8] = gma2*S[4][13][8] + fak * (S[1][13][8] +3.0*S[4][10][8] );
                S[10][13][9] = gma2*S[4][13][9] + fak * (S[1][13][9] +3.0*S[4][10][9] );
                S[10][13][10] = gma2*S[4][13][10] + fak * (S[1][13][10] +3.0*S[4][10][10] +2.0*S[4][13][4] );
                S[10][14][5] = gma2*S[4][14][5] + fak * S[1][14][5];
                S[10][14][6] = gma2*S[4][14][6] + fak * (S[1][14][6] +S[4][14][2] );
                S[10][14][7] = gma2*S[4][14][7] + fak * (S[1][14][7] +S[4][14][3] );
                S[10][14][8] = gma2*S[4][14][8] + fak * S[1][14][8];
                S[10][14][9] = gma2*S[4][14][9] + fak * S[1][14][9];
                S[10][14][10] = gma2*S[4][14][10] + fak *(S[1][14][10]+2.0*S[4][14][4]);
                S[10][15][5] = gma2*S[4][15][5] + fak * S[1][15][5];
                S[10][15][6] = gma2*S[4][15][6] + fak * (S[1][15][6] +S[4][15][2] );
                S[10][15][7] = gma2*S[4][15][7] + fak * (S[1][15][7] +S[4][15][3] );
                S[10][15][8] = gma2*S[4][15][8] + fak * S[1][15][8];
                S[10][15][9] = gma2*S[4][15][9] + fak * S[1][15][9];
                S[10][15][10] = gma2*S[4][15][10] + fak *(S[1][15][10]+2.0*S[4][15][4]);
                S[10][16][5] = gma2*S[4][16][5] + fak * (S[1][16][5] +S[4][8][5] );
                S[10][16][6] = gma2*S[4][16][6] + fak*(S[1][16][6]+S[4][8][6]+S[4][16][2]);
                S[10][16][7] = gma2*S[4][16][7] + fak*(S[1][16][7]+S[4][8][7]+S[4][16][3]);
                S[10][16][8] = gma2*S[4][16][8] + fak * (S[1][16][8] +S[4][8][8] );
                S[10][16][9] = gma2*S[4][16][9] + fak * (S[1][16][9] +S[4][8][9] );
                S[10][16][10] = gma2*S[4][16][10] + fak * (S[1][16][10] +S[4][8][10] +2.0*S[4][16][4] );
                S[10][17][5] = gma2*S[4][17][5] + fak * (S[1][17][5] +2.0*S[4][6][5] );
                S[10][17][6] = gma2*S[4][17][6] + fak * (S[1][17][6] +2.0*S[4][6][6] +S[4][17][2] );
                S[10][17][7] = gma2*S[4][17][7] + fak * (S[1][17][7] +2.0*S[4][6][7] +S[4][17][3] );
                S[10][17][8] = gma2*S[4][17][8] + fak * (S[1][17][8] +2.0*S[4][6][8] );
                S[10][17][9] = gma2*S[4][17][9] + fak * (S[1][17][9] +2.0*S[4][6][9] );
                S[10][17][10] = gma2*S[4][17][10] + fak * (S[1][17][10] +2.0*S[4][6][10] +2.0*S[4][17][4] );
                S[10][18][5] = gma2*S[4][18][5] + fak * (S[1][18][5] +S[4][9][5] );
                S[10][18][6] = gma2*S[4][18][6] + fak*(S[1][18][6]+S[4][9][6]+S[4][18][2]);
                S[10][18][7] = gma2*S[4][18][7] + fak*(S[1][18][7]+S[4][9][7]+S[4][18][3]);
                S[10][18][8] = gma2*S[4][18][8] + fak * (S[1][18][8] +S[4][9][8] );
                S[10][18][9] = gma2*S[4][18][9] + fak * (S[1][18][9] +S[4][9][9] );
                S[10][18][10] = gma2*S[4][18][10] + fak * (S[1][18][10] +S[4][9][10] +2.0*S[4][18][4] );
                S[10][19][5] = gma2*S[4][19][5] + fak * (S[1][19][5] +2.0*S[4][7][5] );
                S[10][19][6] = gma2*S[4][19][6] + fak * (S[1][19][6] +2.0*S[4][7][6] +S[4][19][2] );
                S[10][19][7] = gma2*S[4][19][7] + fak * (S[1][19][7] +2.0*S[4][7][7] +S[4][19][3] );
                S[10][19][8] = gma2*S[4][19][8] + fak * (S[1][19][8] +2.0*S[4][7][8] );
                S[10][19][9] = gma2*S[4][19][9] + fak * (S[1][19][9] +2.0*S[4][7][9] );
                S[10][19][10] = gma2*S[4][19][10] + fak * (S[1][19][10] +2.0*S[4][7][10] +2.0*S[4][19][4] );
                S[10][20][5] = gma2*S[4][20][5] + fak * (S[1][20][5] +S[4][5][5] );
                S[10][20][6] = gma2*S[4][20][6] + fak*(S[1][20][6]+S[4][5][6]+S[4][20][2]);
                S[10][20][7] = gma2*S[4][20][7] + fak*(S[1][20][7]+S[4][5][7]+S[4][20][3]);
                S[10][20][8] = gma2*S[4][20][8] + fak * (S[1][20][8] +S[4][5][8] );
                S[10][20][9] = gma2*S[4][20][9] + fak * (S[1][20][9] +S[4][5][9] );
                S[10][20][10] = gma2*S[4][20][10] + fak * (S[1][20][10] +S[4][5][10] +2.0*S[4][20][4] );
             } 
            
            if ( _lmax_gw > 3) {
            
            if ( _lmax_alpha>0  &&  _lmax_gamma>0 ) {
     S[0][20][0] = gmc0*S[0][10][0]+     fak3* S[0][7][0];
     S[0][23][0] = gmc1*S[0][10][0];
     S[0][25][0] = gmc2*S[0][10][0];
     S[0][24][0] = gmc0*S[0][11][0];
     S[0][21][0] = gmc1*S[0][11][0]+     fak3* S[0][8][0];
     S[0][27][0] = gmc2*S[0][11][0];
     S[0][26][0] = gmc0*S[0][12][0];
     S[0][28][0] = gmc1*S[0][12][0];
     S[0][22][0] = gmc2*S[0][12][0]+     fak3* S[0][9][0];
     S[0][31][0] = gmc1*S[0][13][0]+     fak * S[0][7][0];
     S[0][32][0] = gmc2*S[0][13][0];
     S[0][33][0] = gmc2*S[0][14][0];
     S[0][29][0] = gmc2*S[0][15][0]+     fak * S[0][7][0];
     S[0][34][0] = gmc1*S[0][16][0];
     S[0][30][0] = gmc2*S[0][17][0]+     fak * S[0][8][0];
  }
  if ( _lmax_alpha>0  &&  _lmax_gamma>1 ) {
     S[0][20][1] = gmc0*S[0][10][1]+     fak * (3e0*S[0][7][1] + S[0][10][0] );
     S[0][23][1] = gmc1*S[0][10][1];
     S[0][25][1] = gmc2*S[0][10][1];
     S[0][20][2] = gmc0*S[0][10][2]+     fak3* S[0][7][2];
     S[0][23][2] = gmc1*S[0][10][2]+     fak * S[0][10][0];
     S[0][25][2] = gmc2*S[0][10][2];
     S[0][20][3] = gmc0*S[0][10][3]+     fak3* S[0][7][3];
     S[0][23][3] = gmc1*S[0][10][3];
     S[0][25][3] = gmc2*S[0][10][3]+     fak * S[0][10][0];
     S[0][24][1] = gmc0*S[0][11][1]+     fak * S[0][11][0];
     S[0][21][1] = gmc1*S[0][11][1]+     fak3* S[0][8][1];
     S[0][27][1] = gmc2*S[0][11][1];
     S[0][24][2] = gmc0*S[0][11][2];
     S[0][21][2] = gmc1*S[0][11][2]+     fak * (3e0*S[0][8][2] + S[0][11][0] );
     S[0][27][2] = gmc2*S[0][11][2];
     S[0][24][3] = gmc0*S[0][11][3];
     S[0][21][3] = gmc1*S[0][11][3]+     fak3* S[0][8][3];
     S[0][27][3] = gmc2*S[0][11][3]+     fak * S[0][11][0];
     S[0][26][1] = gmc0*S[0][12][1]+     fak * S[0][12][0];
     S[0][28][1] = gmc1*S[0][12][1];
     S[0][22][1] = gmc2*S[0][12][1]+     fak3* S[0][9][1];
     S[0][26][2] = gmc0*S[0][12][2];
     S[0][28][2] = gmc1*S[0][12][2]+     fak * S[0][12][0];
     S[0][22][2] = gmc2*S[0][12][2]+     fak3* S[0][9][2];
     S[0][26][3] = gmc0*S[0][12][3];
     S[0][28][3] = gmc1*S[0][12][3];
     S[0][22][3] = gmc2*S[0][12][3]+     fak * (3e0*S[0][9][3] + S[0][12][0] );
     S[0][31][1] = gmc1*S[0][13][1]+     fak * S[0][7][1];
     S[0][32][1] = gmc2*S[0][13][1];
     S[0][31][2] = gmc1*S[0][13][2]+     fak * (S[0][7][2] + S[0][13][0] );
     S[0][32][2] = gmc2*S[0][13][2];
     S[0][31][3] = gmc1*S[0][13][3]+     fak * S[0][7][3];
     S[0][32][3] = gmc2*S[0][13][3]+     fak * S[0][13][0];
     S[0][33][1] = gmc2*S[0][14][1];
     S[0][33][2] = gmc2*S[0][14][2];
     S[0][33][3] = gmc2*S[0][14][3]+     fak * S[0][14][0];
     S[0][29][1] = gmc2*S[0][15][1]+     fak * S[0][7][1];
     S[0][29][2] = gmc2*S[0][15][2]+     fak * S[0][7][2];
     S[0][29][3] = gmc2*S[0][15][3]+     fak * (S[0][7][3] + S[0][15][0] );
     S[0][34][1] = gmc1*S[0][16][1];
     S[0][34][2] = gmc1*S[0][16][2]+     fak * S[0][16][0];
     S[0][34][3] = gmc1*S[0][16][3];
     S[0][30][1] = gmc2*S[0][17][1]+     fak * S[0][8][1];
     S[0][30][2] = gmc2*S[0][17][2]+     fak * S[0][8][2];
     S[0][30][3] = gmc2*S[0][17][3]+     fak * (S[0][8][3] + S[0][17][0] );
  }
  if ( _lmax_alpha>1  &&  _lmax_gamma>0 ) {
     S[1][20][0] = gma0*S[0][20][0]+     fak4* S[0][10][0];
     S[2][20][0] = gma1*S[0][20][0];
     S[3][20][0] = gma2*S[0][20][0];
     S[1][21][0] = gma0*S[0][21][0];
     S[2][21][0] = gma1*S[0][21][0]+     fak4* S[0][11][0];
     S[3][21][0] = gma2*S[0][21][0];
     S[1][22][0] = gma0*S[0][22][0];
     S[2][22][0] = gma1*S[0][22][0];
     S[3][22][0] = gma2*S[0][22][0]+     fak4* S[0][12][0];
     S[1][23][0] = gma0*S[0][23][0]+     fak3* S[0][13][0];
     S[2][23][0] = gma1*S[0][23][0]+     fak * S[0][10][0];
     S[3][23][0] = gma2*S[0][23][0];
     S[1][24][0] = gma0*S[0][24][0]+     fak * S[0][11][0];
     S[2][24][0] = gma1*S[0][24][0]+     fak3* S[0][14][0];
     S[3][24][0] = gma2*S[0][24][0];
     S[1][25][0] = gma0*S[0][25][0]+     fak3* S[0][15][0];
     S[2][25][0] = gma1*S[0][25][0];
     S[3][25][0] = gma2*S[0][25][0]+     fak * S[0][10][0];
     S[1][26][0] = gma0*S[0][26][0]+     fak * S[0][12][0];
     S[2][26][0] = gma1*S[0][26][0];
     S[3][26][0] = gma2*S[0][26][0]+     fak3* S[0][16][0];
     S[1][27][0] = gma0*S[0][27][0];
     S[2][27][0] = gma1*S[0][27][0]+     fak3* S[0][17][0];
     S[3][27][0] = gma2*S[0][27][0]+     fak * S[0][11][0];
     S[1][28][0] = gma0*S[0][28][0];
     S[2][28][0] = gma1*S[0][28][0]+     fak * S[0][12][0];
     S[3][28][0] = gma2*S[0][28][0]+     fak3* S[0][18][0];
     S[1][29][0] = gma0*S[0][29][0]+     fak2* S[0][16][0];
     S[2][29][0] = gma1*S[0][29][0];
     S[3][29][0] = gma2*S[0][29][0]+     fak2* S[0][15][0];
     S[1][30][0] = gma0*S[0][30][0];
     S[2][30][0] = gma1*S[0][30][0]+     fak2* S[0][18][0];
     S[3][30][0] = gma2*S[0][30][0]+     fak2* S[0][17][0];
     S[1][31][0] = gma0*S[0][31][0]+     fak2* S[0][14][0];
     S[2][31][0] = gma1*S[0][31][0]+     fak2* S[0][13][0];
     S[3][31][0] = gma2*S[0][31][0];
     S[1][32][0] = gma0*S[0][32][0]+     fak2* S[0][19][0];
     S[2][32][0] = gma1*S[0][32][0]+     fak * S[0][15][0];
     S[3][32][0] = gma2*S[0][32][0]+     fak * S[0][13][0];
     S[1][33][0] = gma0*S[0][33][0]+     fak * S[0][17][0];
     S[2][33][0] = gma1*S[0][33][0]+     fak2* S[0][19][0];
     S[3][33][0] = gma2*S[0][33][0]+     fak * S[0][14][0];
     S[1][34][0] = gma0*S[0][34][0]+     fak * S[0][18][0];
     S[2][34][0] = gma1*S[0][34][0]+     fak * S[0][16][0];
     S[3][34][0] = gma2*S[0][34][0]+     fak2* S[0][19][0];
  }
  if ( _lmax_alpha>0  &&  _lmax_gamma>2 ) {
     S[0][20][4] = gmc0*S[0][10][4]+     fak * (3.e0*S[0][7][4] + S[0][10][2] );
     S[0][23][4] = gmc1*S[0][10][4]+     fak * S[0][10][1];
     S[0][25][4] = gmc2*S[0][10][4];
     S[0][20][5] = gmc0*S[0][10][5]+     fak * (3.e0*S[0][7][5] + S[0][10][3] );
     S[0][23][5] = gmc1*S[0][10][5];
     S[0][25][5] = gmc2*S[0][10][5]+     fak * S[0][10][1];
     S[0][20][6] = gmc0*S[0][10][6]+     fak3* S[0][7][6];
     S[0][23][6] = gmc1*S[0][10][6]+     fak * S[0][10][3];
     S[0][25][6] = gmc2*S[0][10][6]+     fak * S[0][10][2];
     S[0][20][7] = gmc0*S[0][10][7]+     fak * (3.e0*S[0][7][7] + 2.e0*S[0][10][1] );
     S[0][23][7] = gmc1*S[0][10][7];
     S[0][25][7] = gmc2*S[0][10][7];
     S[0][20][8] = gmc0*S[0][10][8]+     fak3* S[0][7][8];
     S[0][23][8] = gmc1*S[0][10][8]+     fak2* S[0][10][2];
     S[0][25][8] = gmc2*S[0][10][8];
     S[0][20][9] = gmc0*S[0][10][9]+     fak3* S[0][7][9];
     S[0][23][9] = gmc1*S[0][10][9];
     S[0][25][9] = gmc2*S[0][10][9]+     fak2* S[0][10][3];
     S[0][24][4] = gmc0*S[0][11][4]+     fak * S[0][11][2];
     S[0][21][4] = gmc1*S[0][11][4]+     fak * (3.e0*S[0][8][4] + S[0][11][1] );
     S[0][27][4] = gmc2*S[0][11][4];
     S[0][24][5] = gmc0*S[0][11][5]+     fak * S[0][11][3];
     S[0][21][5] = gmc1*S[0][11][5]+     fak3* S[0][8][5];
     S[0][27][5] = gmc2*S[0][11][5]+     fak * S[0][11][1];
     S[0][24][6] = gmc0*S[0][11][6];
     S[0][21][6] = gmc1*S[0][11][6]+     fak * (3.e0*S[0][8][6] + S[0][11][3] );
     S[0][27][6] = gmc2*S[0][11][6]+     fak * S[0][11][2];
     S[0][24][7] = gmc0*S[0][11][7]+     fak2* S[0][11][1];
     S[0][21][7] = gmc1*S[0][11][7]+     fak3* S[0][8][7];
     S[0][27][7] = gmc2*S[0][11][7];
     S[0][24][8] = gmc0*S[0][11][8];
     S[0][21][8] = gmc1*S[0][11][8]+     fak * (3.e0*S[0][8][8] + 2.e0*S[0][11][2] );
     S[0][27][8] = gmc2*S[0][11][8];
     S[0][24][9] = gmc0*S[0][11][9];
     S[0][21][9] = gmc1*S[0][11][9]+     fak3* S[0][8][9];
     S[0][27][9] = gmc2*S[0][11][9]+     fak2* S[0][11][3];
     S[0][26][4] = gmc0*S[0][12][4]+     fak * S[0][12][2];
     S[0][28][4] = gmc1*S[0][12][4]+     fak * S[0][12][1];
     S[0][22][4] = gmc2*S[0][12][4]+     fak3* S[0][9][4];
     S[0][26][5] = gmc0*S[0][12][5]+     fak * S[0][12][3];
     S[0][28][5] = gmc1*S[0][12][5];
     S[0][22][5] = gmc2*S[0][12][5]+     fak * (3.e0*S[0][9][5] + S[0][12][1] );
     S[0][26][6] = gmc0*S[0][12][6];
     S[0][28][6] = gmc1*S[0][12][6]+     fak * S[0][12][3];
     S[0][22][6] = gmc2*S[0][12][6]+     fak * (3.e0*S[0][9][6] + S[0][12][2] );
     S[0][26][7] = gmc0*S[0][12][7]+     fak2* S[0][12][1];
     S[0][28][7] = gmc1*S[0][12][7];
     S[0][22][7] = gmc2*S[0][12][7]+     fak3* S[0][9][7];
     S[0][26][8] = gmc0*S[0][12][8];
     S[0][28][8] = gmc1*S[0][12][8]+     fak2* S[0][12][2];
     S[0][22][8] = gmc2*S[0][12][8]+     fak3* S[0][9][8];
     S[0][26][9] = gmc0*S[0][12][9];
     S[0][28][9] = gmc1*S[0][12][9];
     S[0][22][9] = gmc2*S[0][12][9]+     fak * (3.e0*S[0][9][9] + 2.e0*S[0][12][3] );
     S[0][31][4] = gmc1*S[0][13][4]+     fak * (S[0][7][4] + S[0][13][1] );
     S[0][32][4] = gmc2*S[0][13][4];
     S[0][31][5] = gmc1*S[0][13][5]+     fak * S[0][7][5];
     S[0][32][5] = gmc2*S[0][13][5]+     fak * S[0][13][1];
     S[0][31][6] = gmc1*S[0][13][6]+     fak * (S[0][7][6] + S[0][13][3] );
     S[0][32][6] = gmc2*S[0][13][6]+     fak * S[0][13][2];
     S[0][31][7] = gmc1*S[0][13][7]+     fak * S[0][7][7];
     S[0][32][7] = gmc2*S[0][13][7];
     S[0][31][8] = gmc1*S[0][13][8]+     fak * (S[0][7][8] + 2.e0*S[0][13][2] );
     S[0][32][8] = gmc2*S[0][13][8];
     S[0][31][9] = gmc1*S[0][13][9]+     fak * S[0][7][9];
     S[0][32][9] = gmc2*S[0][13][9]+     fak2* S[0][13][3];
     S[0][33][4] = gmc2*S[0][14][4];
     S[0][33][5] = gmc2*S[0][14][5]+     fak * S[0][14][1];
     S[0][33][6] = gmc2*S[0][14][6]+     fak * S[0][14][2];
     S[0][33][7] = gmc2*S[0][14][7];
     S[0][33][8] = gmc2*S[0][14][8];
     S[0][33][9] = gmc2*S[0][14][9]+     fak2* S[0][14][3];
     S[0][29][4] = gmc2*S[0][15][4]+     fak * S[0][7][4];
     S[0][29][5] = gmc2*S[0][15][5]+     fak * (S[0][7][5] + S[0][15][1] );
     S[0][29][6] = gmc2*S[0][15][6]+     fak * (S[0][7][6] + S[0][15][2] );
     S[0][29][7] = gmc2*S[0][15][7]+     fak * S[0][7][7];
     S[0][29][8] = gmc2*S[0][15][8]+     fak * S[0][7][8];
     S[0][29][9] = gmc2*S[0][15][9]+     fak * (S[0][7][9] + 2.e0*S[0][15][3] );
     S[0][34][4] = gmc1*S[0][16][4]+     fak * S[0][16][1];
     S[0][34][5] = gmc1*S[0][16][5];
     S[0][34][6] = gmc1*S[0][16][6]+     fak * S[0][16][3];
     S[0][34][7] = gmc1*S[0][16][7];
     S[0][34][8] = gmc1*S[0][16][8]+     fak2* S[0][16][2];
     S[0][34][9] = gmc1*S[0][16][9];
     S[0][30][4] = gmc2*S[0][17][4]+     fak * S[0][8][4];
     S[0][30][5] = gmc2*S[0][17][5]+     fak * (S[0][8][5] + S[0][17][1] );
     S[0][30][6] = gmc2*S[0][17][6]+     fak * (S[0][8][6] + S[0][17][2] );
     S[0][30][7] = gmc2*S[0][17][7]+     fak * S[0][8][7];
     S[0][30][8] = gmc2*S[0][17][8]+     fak * S[0][8][8];
     S[0][30][9] = gmc2*S[0][17][9]+     fak * (S[0][8][9] + 2.e0*S[0][17][3] );
  }
  if ( _lmax_alpha>1  &&  _lmax_gamma>1 ) {
     S[1][20][1] = gma0*S[0][20][1]+     fak * (4.e0*S[0][10][1] + S[0][20][0] );
     S[2][20][1] = gma1*S[0][20][1];
     S[3][20][1] = gma2*S[0][20][1];
     S[1][20][2] = gma0*S[0][20][2]+     fak4* S[0][10][2];
     S[2][20][2] = gma1*S[0][20][2]+     fak * S[0][20][0];
     S[3][20][2] = gma2*S[0][20][2];
     S[1][20][3] = gma0*S[0][20][3]+     fak4* S[0][10][3];
     S[2][20][3] = gma1*S[0][20][3];
     S[3][20][3] = gma2*S[0][20][3]+     fak * S[0][20][0];
     S[1][21][1] = gma0*S[0][21][1]+     fak * S[0][21][0];
     S[2][21][1] = gma1*S[0][21][1]+     fak4* S[0][11][1];
     S[3][21][1] = gma2*S[0][21][1];
     S[1][21][2] = gma0*S[0][21][2];
     S[2][21][2] = gma1*S[0][21][2]+     fak * (4.e0*S[0][11][2] + S[0][21][0] );
     S[3][21][2] = gma2*S[0][21][2];
     S[1][21][3] = gma0*S[0][21][3];
     S[2][21][3] = gma1*S[0][21][3]+     fak4* S[0][11][3];
     S[3][21][3] = gma2*S[0][21][3]+     fak * S[0][21][0];
     S[1][22][1] = gma0*S[0][22][1]+     fak * S[0][22][0];
     S[2][22][1] = gma1*S[0][22][1];
     S[3][22][1] = gma2*S[0][22][1]+     fak4* S[0][12][1];
     S[1][22][2] = gma0*S[0][22][2];
     S[2][22][2] = gma1*S[0][22][2]+     fak * S[0][22][0];
     S[3][22][2] = gma2*S[0][22][2]+     fak4* S[0][12][2];
     S[1][22][3] = gma0*S[0][22][3];
     S[2][22][3] = gma1*S[0][22][3];
     S[3][22][3] = gma2*S[0][22][3]+     fak * (4.e0*S[0][12][3] + S[0][22][0] );
     S[1][23][1] = gma0*S[0][23][1]+     fak * (3.e0*S[0][13][1] + S[0][23][0] );
     S[2][23][1] = gma1*S[0][23][1]+     fak * S[0][10][1];
     S[3][23][1] = gma2*S[0][23][1];
     S[1][23][2] = gma0*S[0][23][2]+     fak3* S[0][13][2];
     S[2][23][2] = gma1*S[0][23][2]+     fak * (S[0][10][2] + S[0][23][0] );
     S[3][23][2] = gma2*S[0][23][2];
     S[1][23][3] = gma0*S[0][23][3]+     fak3* S[0][13][3];
     S[2][23][3] = gma1*S[0][23][3]+     fak * S[0][10][3];
     S[3][23][3] = gma2*S[0][23][3]+     fak * S[0][23][0];
     S[1][24][1] = gma0*S[0][24][1]+     fak * (S[0][11][1] + S[0][24][0] );
     S[2][24][1] = gma1*S[0][24][1]+     fak3* S[0][14][1];
     S[3][24][1] = gma2*S[0][24][1];
     S[1][24][2] = gma0*S[0][24][2]+     fak * S[0][11][2];
     S[2][24][2] = gma1*S[0][24][2]+     fak * (3.e0*S[0][14][2] + S[0][24][0] );
     S[3][24][2] = gma2*S[0][24][2];
     S[1][24][3] = gma0*S[0][24][3]+     fak * S[0][11][3];
     S[2][24][3] = gma1*S[0][24][3]+     fak3* S[0][14][3];
     S[3][24][3] = gma2*S[0][24][3]+     fak * S[0][24][0];
     S[1][25][1] = gma0*S[0][25][1]+     fak * (3.e0*S[0][15][1] + S[0][25][0] );
     S[2][25][1] = gma1*S[0][25][1];
     S[3][25][1] = gma2*S[0][25][1]+     fak * S[0][10][1];
     S[1][25][2] = gma0*S[0][25][2]+     fak3* S[0][15][2];
     S[2][25][2] = gma1*S[0][25][2]+     fak * S[0][25][0];
     S[3][25][2] = gma2*S[0][25][2]+     fak * S[0][10][2];
     S[1][25][3] = gma0*S[0][25][3]+     fak3* S[0][15][3];
     S[2][25][3] = gma1*S[0][25][3];
     S[3][25][3] = gma2*S[0][25][3]+     fak * (S[0][10][3] + S[0][25][0] );
     S[1][26][1] = gma0*S[0][26][1]+     fak * (S[0][12][1] + S[0][26][0] );
     S[2][26][1] = gma1*S[0][26][1];
     S[3][26][1] = gma2*S[0][26][1]+     fak3* S[0][16][1];
     S[1][26][2] = gma0*S[0][26][2]+     fak * S[0][12][2];
     S[2][26][2] = gma1*S[0][26][2]+     fak * S[0][26][0];
     S[3][26][2] = gma2*S[0][26][2]+     fak3* S[0][16][2];
     S[1][26][3] = gma0*S[0][26][3]+     fak * S[0][12][3];
     S[2][26][3] = gma1*S[0][26][3];
     S[3][26][3] = gma2*S[0][26][3]+     fak * (3.e0*S[0][16][3] + S[0][26][0] );
     S[1][27][1] = gma0*S[0][27][1]+     fak * S[0][27][0];
     S[2][27][1] = gma1*S[0][27][1]+     fak3* S[0][17][1];
     S[3][27][1] = gma2*S[0][27][1]+     fak * S[0][11][1];
     S[1][27][2] = gma0*S[0][27][2];
     S[2][27][2] = gma1*S[0][27][2]+     fak * (3.e0*S[0][17][2] + S[0][27][0] );
     S[3][27][2] = gma2*S[0][27][2]+     fak * S[0][11][2];
     S[1][27][3] = gma0*S[0][27][3];
     S[2][27][3] = gma1*S[0][27][3]+     fak3* S[0][17][3];
     S[3][27][3] = gma2*S[0][27][3]+     fak * (S[0][11][3] + S[0][27][0] );
     S[1][28][1] = gma0*S[0][28][1]+     fak * S[0][28][0];
     S[2][28][1] = gma1*S[0][28][1]+     fak * S[0][12][1];
     S[3][28][1] = gma2*S[0][28][1]+     fak3* S[0][18][1];
     S[1][28][2] = gma0*S[0][28][2];
     S[2][28][2] = gma1*S[0][28][2]+     fak * (S[0][12][2] + S[0][28][0] );
     S[3][28][2] = gma2*S[0][28][2]+     fak3* S[0][18][2];
     S[1][28][3] = gma0*S[0][28][3];
     S[2][28][3] = gma1*S[0][28][3]+     fak * S[0][12][3];
     S[3][28][3] = gma2*S[0][28][3]+     fak * (3.e0*S[0][18][3] + S[0][28][0] );
     S[1][29][1] = gma0*S[0][29][1]+     fak * (2.e0*S[0][16][1] + S[0][29][0] );
     S[2][29][1] = gma1*S[0][29][1];
     S[3][29][1] = gma2*S[0][29][1]+     fak2* S[0][15][1];
     S[1][29][2] = gma0*S[0][29][2]+     fak2* S[0][16][2];
     S[2][29][2] = gma1*S[0][29][2]+     fak * S[0][29][0];
     S[3][29][2] = gma2*S[0][29][2]+     fak2* S[0][15][2];
     S[1][29][3] = gma0*S[0][29][3]+     fak2* S[0][16][3];
     S[2][29][3] = gma1*S[0][29][3];
     S[3][29][3] = gma2*S[0][29][3]+     fak * (2.e0*S[0][15][3] + S[0][29][0] );
     S[1][30][1] = gma0*S[0][30][1]+     fak * S[0][30][0];
     S[2][30][1] = gma1*S[0][30][1]+     fak2* S[0][18][1];
     S[3][30][1] = gma2*S[0][30][1]+     fak2* S[0][17][1];
     S[1][30][2] = gma0*S[0][30][2];
     S[2][30][2] = gma1*S[0][30][2]+     fak * (2.e0*S[0][18][2] + S[0][30][0] );
     S[3][30][2] = gma2*S[0][30][2]+     fak2* S[0][17][2];
     S[1][30][3] = gma0*S[0][30][3];
     S[2][30][3] = gma1*S[0][30][3]+     fak2* S[0][18][3];
     S[3][30][3] = gma2*S[0][30][3]+     fak * (2.e0*S[0][17][3] + S[0][30][0] );
     S[1][31][1] = gma0*S[0][31][1]+     fak * (2.e0*S[0][14][1] + S[0][31][0] );
     S[2][31][1] = gma1*S[0][31][1]+     fak2* S[0][13][1];
     S[3][31][1] = gma2*S[0][31][1];
     S[1][31][2] = gma0*S[0][31][2]+     fak2* S[0][14][2];
     S[2][31][2] = gma1*S[0][31][2]+     fak * (2.e0*S[0][13][2] + S[0][31][0] );
     S[3][31][2] = gma2*S[0][31][2];
     S[1][31][3] = gma0*S[0][31][3]+     fak2* S[0][14][3];
     S[2][31][3] = gma1*S[0][31][3]+     fak2* S[0][13][3];
     S[3][31][3] = gma2*S[0][31][3]+     fak * S[0][31][0];
     S[1][32][1] = gma0*S[0][32][1]+     fak * (2.e0*S[0][19][1] + S[0][32][0] );
     S[2][32][1] = gma1*S[0][32][1]+     fak * S[0][15][1];
     S[3][32][1] = gma2*S[0][32][1]+     fak * S[0][13][1];
     S[1][32][2] = gma0*S[0][32][2]+     fak2* S[0][19][2];
     S[2][32][2] = gma1*S[0][32][2]+     fak * (S[0][15][2] + S[0][32][0] );
     S[3][32][2] = gma2*S[0][32][2]+     fak * S[0][13][2];
     S[1][32][3] = gma0*S[0][32][3]+     fak2* S[0][19][3];
     S[2][32][3] = gma1*S[0][32][3]+     fak * S[0][15][3];
     S[3][32][3] = gma2*S[0][32][3]+     fak * (S[0][13][3] + S[0][32][0] );
     S[1][33][1] = gma0*S[0][33][1]+     fak * (S[0][17][1] + S[0][33][0] );
     S[2][33][1] = gma1*S[0][33][1]+     fak2* S[0][19][1];
     S[3][33][1] = gma2*S[0][33][1]+     fak * S[0][14][1];
     S[1][33][2] = gma0*S[0][33][2]+     fak * S[0][17][2];
     S[2][33][2] = gma1*S[0][33][2]+     fak * (2.e0*S[0][19][2] + S[0][33][0] );
     S[3][33][2] = gma2*S[0][33][2]+     fak * S[0][14][2];
     S[1][33][3] = gma0*S[0][33][3]+     fak * S[0][17][3];
     S[2][33][3] = gma1*S[0][33][3]+     fak2* S[0][19][3];
     S[3][33][3] = gma2*S[0][33][3]+     fak * (S[0][14][3] + S[0][33][0] );
     S[1][34][1] = gma0*S[0][34][1]+     fak * (S[0][18][1] + S[0][34][0] );
     S[2][34][1] = gma1*S[0][34][1]+     fak * S[0][16][1];
     S[3][34][1] = gma2*S[0][34][1]+     fak2* S[0][19][1];
     S[1][34][2] = gma0*S[0][34][2]+     fak * S[0][18][2];
     S[2][34][2] = gma1*S[0][34][2]+     fak * (S[0][16][2] + S[0][34][0] );
     S[3][34][2] = gma2*S[0][34][2]+     fak2* S[0][19][2];
     S[1][34][3] = gma0*S[0][34][3]+     fak * S[0][18][3];
     S[2][34][3] = gma1*S[0][34][3]+     fak * S[0][16][3];
     S[3][34][3] = gma2*S[0][34][3]+     fak * (2.e0*S[0][19][3] + S[0][34][0] );
  }
  if ( _lmax_alpha>2  &&  _lmax_gamma>0 ) {
     S[7][20][0] = gma0*S[1][20][0]+     fak * (S[0][20][0] + 4.e0*S[1][10][0] ) ;
     S[4][20][0] = gma1*S[1][20][0];
     S[5][20][0] = gma2*S[1][20][0];
     S[7][21][0] = gma0*S[1][21][0]+     fak * S[0][21][0];
     S[4][21][0] = gma1*S[1][21][0]+     fak4* S[1][11][0];
     S[5][21][0] = gma2*S[1][21][0];
     S[7][22][0] = gma0*S[1][22][0]+     fak * S[0][22][0];
     S[4][22][0] = gma1*S[1][22][0];
     S[5][22][0] = gma2*S[1][22][0]+     fak4* S[1][12][0];
     S[7][23][0] = gma0*S[1][23][0]+     fak * (S[0][23][0] + 3.e0*S[1][13][0] ) ;
     S[4][23][0] = gma1*S[1][23][0]+     fak * S[1][10][0];
     S[5][23][0] = gma2*S[1][23][0];
     S[7][24][0] = gma0*S[1][24][0]+     fak * (S[0][24][0] + S[1][11][0] ) ;
     S[4][24][0] = gma1*S[1][24][0]+     fak3* S[1][14][0];
     S[5][24][0] = gma2*S[1][24][0];
     S[7][25][0] = gma0*S[1][25][0]+     fak * (S[0][25][0] + 3.e0*S[1][15][0] ) ;
     S[4][25][0] = gma1*S[1][25][0];
     S[5][25][0] = gma2*S[1][25][0]+     fak * S[1][10][0];
     S[7][26][0] = gma0*S[1][26][0]+     fak * (S[0][26][0] + S[1][12][0] ) ;
     S[4][26][0] = gma1*S[1][26][0];
     S[5][26][0] = gma2*S[1][26][0]+     fak3* S[1][16][0];
     S[7][27][0] = gma0*S[1][27][0]+     fak * S[0][27][0];
     S[4][27][0] = gma1*S[1][27][0]+     fak3* S[1][17][0];
     S[5][27][0] = gma2*S[1][27][0]+     fak * S[1][11][0];
     S[7][28][0] = gma0*S[1][28][0]+     fak * S[0][28][0];
     S[4][28][0] = gma1*S[1][28][0]+     fak * S[1][12][0];
     S[5][28][0] = gma2*S[1][28][0]+     fak3* S[1][18][0];
     S[7][29][0] = gma0*S[1][29][0]+     fak * (S[0][29][0] + 2.e0*S[1][16][0] ) ;
     S[4][29][0] = gma1*S[1][29][0];
     S[5][29][0] = gma2*S[1][29][0]+     fak2* S[1][15][0];
     S[7][30][0] = gma0*S[1][30][0]+     fak * S[0][30][0];
     S[4][30][0] = gma1*S[1][30][0]+     fak2* S[1][18][0];
     S[5][30][0] = gma2*S[1][30][0]+     fak2* S[1][17][0];
     S[7][31][0] = gma0*S[1][31][0]+     fak * (S[0][31][0] + 2.e0*S[1][14][0] ) ;
     S[4][31][0] = gma1*S[1][31][0]+     fak2* S[1][13][0];
     S[5][31][0] = gma2*S[1][31][0];
     S[7][32][0] = gma0*S[1][32][0]+     fak * (S[0][32][0] + 2.e0*S[1][19][0] ) ;
     S[4][32][0] = gma1*S[1][32][0]+     fak * S[1][15][0];
     S[5][32][0] = gma2*S[1][32][0]+     fak * S[1][13][0];
     S[7][33][0] = gma0*S[1][33][0]+     fak * (S[0][33][0] + S[1][17][0] ) ;
     S[4][33][0] = gma1*S[1][33][0]+     fak2* S[1][19][0];
     S[5][33][0] = gma2*S[1][33][0]+     fak * S[1][14][0];
     S[7][34][0] = gma0*S[1][34][0]+     fak * (S[0][34][0] + S[1][18][0] ) ;
     S[4][34][0] = gma1*S[1][34][0]+     fak * S[1][16][0];
     S[5][34][0] = gma2*S[1][34][0]+     fak2* S[1][19][0];
     S[8][20][0] = gma1*S[2][20][0]+     fak * S[0][20][0];
     S[6][20][0] = gma2*S[2][20][0];
     S[8][21][0] = gma1*S[2][21][0]+     fak * (S[0][21][0] + 4.e0*S[2][11][0] ) ;
     S[6][21][0] = gma2*S[2][21][0];
     S[8][22][0] = gma1*S[2][22][0]+     fak * S[0][22][0];
     S[6][22][0] = gma2*S[2][22][0]+     fak4* S[2][12][0];
     S[8][23][0] = gma1*S[2][23][0]+     fak * (S[0][23][0] + S[2][10][0] ) ;
     S[6][23][0] = gma2*S[2][23][0];
     S[8][24][0] = gma1*S[2][24][0]+     fak * (S[0][24][0] + 3.e0*S[2][14][0] ) ;
     S[6][24][0] = gma2*S[2][24][0];
     S[8][25][0] = gma1*S[2][25][0]+     fak * S[0][25][0];
     S[6][25][0] = gma2*S[2][25][0]+     fak * S[2][10][0];
     S[8][26][0] = gma1*S[2][26][0]+     fak * S[0][26][0];
     S[6][26][0] = gma2*S[2][26][0]+     fak3* S[2][16][0];
     S[8][27][0] = gma1*S[2][27][0]+     fak * (S[0][27][0] + 3.e0*S[2][17][0] ) ;
     S[6][27][0] = gma2*S[2][27][0]+     fak * S[2][11][0];
     S[8][28][0] = gma1*S[2][28][0]+     fak * (S[0][28][0] + S[2][12][0] ) ;
     S[6][28][0] = gma2*S[2][28][0]+     fak3* S[2][18][0];
     S[8][29][0] = gma1*S[2][29][0]+     fak * S[0][29][0];
     S[6][29][0] = gma2*S[2][29][0]+     fak2* S[2][15][0];
     S[8][30][0] = gma1*S[2][30][0]+     fak * (S[0][30][0] + 2.e0*S[2][18][0] ) ;
     S[6][30][0] = gma2*S[2][30][0]+     fak2* S[2][17][0];
     S[8][31][0] = gma1*S[2][31][0]+     fak * (S[0][31][0] + 2.e0*S[2][13][0] ) ;
     S[6][31][0] = gma2*S[2][31][0];
     S[8][32][0] = gma1*S[2][32][0]+     fak * (S[0][32][0] + S[2][15][0] ) ;
     S[6][32][0] = gma2*S[2][32][0]+     fak * S[2][13][0];
     S[8][33][0] = gma1*S[2][33][0]+     fak * (S[0][33][0] + 2.e0*S[2][19][0] ) ;
     S[6][33][0] = gma2*S[2][33][0]+     fak * S[2][14][0];
     S[8][34][0] = gma1*S[2][34][0]+     fak * (S[0][34][0] + S[2][16][0] ) ;
     S[6][34][0] = gma2*S[2][34][0]+     fak2* S[2][19][0];
     S[9][20][0] = gma2*S[3][20][0]+     fak * S[0][20][0];
     S[9][21][0] = gma2*S[3][21][0]+     fak * S[0][21][0];
     S[9][22][0] = gma2*S[3][22][0]+     fak * (S[0][22][0] + 4.e0*S[3][12][0] ) ;
     S[9][23][0] = gma2*S[3][23][0]+     fak * S[0][23][0];
     S[9][24][0] = gma2*S[3][24][0]+     fak * S[0][24][0];
     S[9][25][0] = gma2*S[3][25][0]+     fak * (S[0][25][0] + S[3][10][0] ) ;
     S[9][26][0] = gma2*S[3][26][0]+     fak * (S[0][26][0] + 3.e0*S[3][16][0] ) ;
     S[9][27][0] = gma2*S[3][27][0]+     fak * (S[0][27][0] + S[3][11][0] ) ;
     S[9][28][0] = gma2*S[3][28][0]+     fak * (S[0][28][0] + 3.e0*S[3][18][0] ) ;
     S[9][29][0] = gma2*S[3][29][0]+     fak * (S[0][29][0] + 2.e0*S[3][15][0] ) ;
     S[9][30][0] = gma2*S[3][30][0]+     fak * (S[0][30][0] + 2.e0*S[3][17][0] ) ;
     S[9][31][0] = gma2*S[3][31][0]+     fak * S[0][31][0];
     S[9][32][0] = gma2*S[3][32][0]+     fak * (S[0][32][0] + S[3][13][0] ) ;
     S[9][33][0] = gma2*S[3][33][0]+     fak * (S[0][33][0] + S[3][14][0] ) ;
     S[9][34][0] = gma2*S[3][34][0]+     fak * (S[0][34][0] + 2.e0*S[3][19][0] ) ;
  }
  if ( _lmax_alpha>1  &&  _lmax_gamma>2 ) {
     S[1][20][4] = gma0*S[0][20][4]+     fak * (4.e0*S[0][10][4] + S[0][20][2] );
     S[2][20][4] = gma1*S[0][20][4]+     fak * S[0][20][1];
     S[3][20][4] = gma2*S[0][20][4];
     S[1][20][5] = gma0*S[0][20][5]+     fak * (4.e0*S[0][10][5] + S[0][20][3] );
     S[2][20][5] = gma1*S[0][20][5];
     S[3][20][5] = gma2*S[0][20][5]+     fak * S[0][20][1];
     S[1][20][6] = gma0*S[0][20][6]+     fak4* S[0][10][6];
     S[2][20][6] = gma1*S[0][20][6]+     fak * S[0][20][3];
     S[3][20][6] = gma2*S[0][20][6]+     fak * S[0][20][2];
     S[1][20][7] = gma0*S[0][20][7]+     fak2* (2.e0*S[0][10][7] + S[0][20][1] );
     S[2][20][7] = gma1*S[0][20][7];
     S[3][20][7] = gma2*S[0][20][7];
     S[1][20][8] = gma0*S[0][20][8]+     fak4* S[0][10][8];
     S[2][20][8] = gma1*S[0][20][8]+     fak2* S[0][20][2];
     S[3][20][8] = gma2*S[0][20][8];
     S[1][20][9] = gma0*S[0][20][9]+     fak4* S[0][10][9];
     S[2][20][9] = gma1*S[0][20][9];
     S[3][20][9] = gma2*S[0][20][9]+     fak2* S[0][20][3];
     S[1][21][4] = gma0*S[0][21][4]+     fak * S[0][21][2];
     S[2][21][4] = gma1*S[0][21][4]+     fak * (4.e0*S[0][11][4] + S[0][21][1] );
     S[3][21][4] = gma2*S[0][21][4];
     S[1][21][5] = gma0*S[0][21][5]+     fak * S[0][21][3];
     S[2][21][5] = gma1*S[0][21][5]+     fak4* S[0][11][5];
     S[3][21][5] = gma2*S[0][21][5]+     fak * S[0][21][1];
     S[1][21][6] = gma0*S[0][21][6];
     S[2][21][6] = gma1*S[0][21][6]+     fak * (4.e0*S[0][11][6] + S[0][21][3] );
     S[3][21][6] = gma2*S[0][21][6]+     fak * S[0][21][2];
     S[1][21][7] = gma0*S[0][21][7]+     fak2* S[0][21][1];
     S[2][21][7] = gma1*S[0][21][7]+     fak4* S[0][11][7];
     S[3][21][7] = gma2*S[0][21][7];
     S[1][21][8] = gma0*S[0][21][8];
     S[2][21][8] = gma1*S[0][21][8]+     fak2* (2.e0*S[0][11][8] + S[0][21][2] );
     S[3][21][8] = gma2*S[0][21][8];
     S[1][21][9] = gma0*S[0][21][9];
     S[2][21][9] = gma1*S[0][21][9]+     fak4* S[0][11][9];
     S[3][21][9] = gma2*S[0][21][9]+     fak2* S[0][21][3];
     S[1][22][4] = gma0*S[0][22][4]+     fak * S[0][22][2];
     S[2][22][4] = gma1*S[0][22][4]+     fak * S[0][22][1];
     S[3][22][4] = gma2*S[0][22][4]+     fak4* S[0][12][4];
     S[1][22][5] = gma0*S[0][22][5]+     fak * S[0][22][3];
     S[2][22][5] = gma1*S[0][22][5];
     S[3][22][5] = gma2*S[0][22][5]+     fak * (4.e0*S[0][12][5] + S[0][22][1] );
     S[1][22][6] = gma0*S[0][22][6];
     S[2][22][6] = gma1*S[0][22][6]+     fak * S[0][22][3];
     S[3][22][6] = gma2*S[0][22][6]+     fak * (4.e0*S[0][12][6] + S[0][22][2] );
     S[1][22][7] = gma0*S[0][22][7]+     fak2* S[0][22][1];
     S[2][22][7] = gma1*S[0][22][7];
     S[3][22][7] = gma2*S[0][22][7]+     fak4* S[0][12][7];
     S[1][22][8] = gma0*S[0][22][8];
     S[2][22][8] = gma1*S[0][22][8]+     fak2* S[0][22][2];
     S[3][22][8] = gma2*S[0][22][8]+     fak4* S[0][12][8];
     S[1][22][9] = gma0*S[0][22][9];
     S[2][22][9] = gma1*S[0][22][9];
     S[3][22][9] = gma2*S[0][22][9]+     fak2* (2.e0*S[0][12][9] + S[0][22][3] );
     S[1][23][4] = gma0*S[0][23][4]+     fak * (3.e0*S[0][13][4] + S[0][23][2] );
     S[2][23][4] = gma1*S[0][23][4]+     fak * (S[0][10][4] + S[0][23][1] );
     S[3][23][4] = gma2*S[0][23][4];
     S[1][23][5] = gma0*S[0][23][5]+     fak * (3.e0*S[0][13][5] + S[0][23][3] );
     S[2][23][5] = gma1*S[0][23][5]+     fak * S[0][10][5];
     S[3][23][5] = gma2*S[0][23][5]+     fak * S[0][23][1];
     S[1][23][6] = gma0*S[0][23][6]+     fak3* S[0][13][6];
     S[2][23][6] = gma1*S[0][23][6]+     fak * (S[0][10][6] + S[0][23][3] );
     S[3][23][6] = gma2*S[0][23][6]+     fak * S[0][23][2];
     S[1][23][7] = gma0*S[0][23][7]+     fak * (3.e0*S[0][13][7] + 2.e0*S[0][23][1] );
     S[2][23][7] = gma1*S[0][23][7]+     fak * S[0][10][7];
     S[3][23][7] = gma2*S[0][23][7];
     S[1][23][8] = gma0*S[0][23][8]+     fak3* S[0][13][8];
     S[2][23][8] = gma1*S[0][23][8]+     fak * (S[0][10][8] + 2.e0*S[0][23][2] );
     S[3][23][8] = gma2*S[0][23][8];
     S[1][23][9] = gma0*S[0][23][9]+     fak3* S[0][13][9];
     S[2][23][9] = gma1*S[0][23][9]+     fak * S[0][10][9];
     S[3][23][9] = gma2*S[0][23][9]+     fak2* S[0][23][3];
     S[1][24][4] = gma0*S[0][24][4]+     fak * (S[0][11][4] + S[0][24][2] );
     S[2][24][4] = gma1*S[0][24][4]+     fak * (3.e0*S[0][14][4] + S[0][24][1] );
     S[3][24][4] = gma2*S[0][24][4];
     S[1][24][5] = gma0*S[0][24][5]+     fak * (S[0][11][5] + S[0][24][3] );
     S[2][24][5] = gma1*S[0][24][5]+     fak3* S[0][14][5];
     S[3][24][5] = gma2*S[0][24][5]+     fak * S[0][24][1];
     S[1][24][6] = gma0*S[0][24][6]+     fak * S[0][11][6];
     S[2][24][6] = gma1*S[0][24][6]+     fak * (3.e0*S[0][14][6] + S[0][24][3] );
     S[3][24][6] = gma2*S[0][24][6]+     fak * S[0][24][2];
     S[1][24][7] = gma0*S[0][24][7]+     fak * (S[0][11][7] + 2.e0*S[0][24][1] );
     S[2][24][7] = gma1*S[0][24][7]+     fak3* S[0][14][7];
     S[3][24][7] = gma2*S[0][24][7];
     S[1][24][8] = gma0*S[0][24][8]+     fak * S[0][11][8];
     S[2][24][8] = gma1*S[0][24][8]+     fak * (3.e0*S[0][14][8] + 2.e0*S[0][24][2] );
     S[3][24][8] = gma2*S[0][24][8];
     S[1][24][9] = gma0*S[0][24][9]+     fak * S[0][11][9];
     S[2][24][9] = gma1*S[0][24][9]+     fak3* S[0][14][9];
     S[3][24][9] = gma2*S[0][24][9]+     fak2* S[0][24][3];
     S[1][25][4] = gma0*S[0][25][4]+     fak * (3.e0*S[0][15][4] + S[0][25][2] );
     S[2][25][4] = gma1*S[0][25][4]+     fak * S[0][25][1];
     S[3][25][4] = gma2*S[0][25][4]+     fak * S[0][10][4];
     S[1][25][5] = gma0*S[0][25][5]+     fak * (3.e0*S[0][15][5] + S[0][25][3] );
     S[2][25][5] = gma1*S[0][25][5];
     S[3][25][5] = gma2*S[0][25][5]+     fak * (S[0][10][5] + S[0][25][1] );
     S[1][25][6] = gma0*S[0][25][6]+     fak3* S[0][15][6];
     S[2][25][6] = gma1*S[0][25][6]+     fak * S[0][25][3];
     S[3][25][6] = gma2*S[0][25][6]+     fak * (S[0][10][6] + S[0][25][2] );
     S[1][25][7] = gma0*S[0][25][7]+     fak * (3.e0*S[0][15][7] + 2.e0*S[0][25][1] );
     S[2][25][7] = gma1*S[0][25][7];
     S[3][25][7] = gma2*S[0][25][7]+     fak * S[0][10][7];
     S[1][25][8] = gma0*S[0][25][8]+     fak3* S[0][15][8];
     S[2][25][8] = gma1*S[0][25][8]+     fak2* S[0][25][2];
     S[3][25][8] = gma2*S[0][25][8]+     fak * S[0][10][8];
     S[1][25][9] = gma0*S[0][25][9]+     fak3* S[0][15][9];
     S[2][25][9] = gma1*S[0][25][9];
     S[3][25][9] = gma2*S[0][25][9]+     fak * (S[0][10][9] + 2.e0*S[0][25][3] );
     S[1][26][4] = gma0*S[0][26][4]+     fak * (S[0][12][4] + S[0][26][2] );
     S[2][26][4] = gma1*S[0][26][4]+     fak * S[0][26][1];
     S[3][26][4] = gma2*S[0][26][4]+     fak3* S[0][16][4];
     S[1][26][5] = gma0*S[0][26][5]+     fak * (S[0][12][5] + S[0][26][3] );
     S[2][26][5] = gma1*S[0][26][5];
     S[3][26][5] = gma2*S[0][26][5]+     fak * (3.e0*S[0][16][5] + S[0][26][1] );
     S[1][26][6] = gma0*S[0][26][6]+     fak * S[0][12][6];
     S[2][26][6] = gma1*S[0][26][6]+     fak * S[0][26][3];
     S[3][26][6] = gma2*S[0][26][6]+     fak * (3.e0*S[0][16][6] + S[0][26][2] );
     S[1][26][7] = gma0*S[0][26][7]+     fak * (S[0][12][7] + 2.e0*S[0][26][1] );
     S[2][26][7] = gma1*S[0][26][7];
     S[3][26][7] = gma2*S[0][26][7]+     fak3* S[0][16][7];
     S[1][26][8] = gma0*S[0][26][8]+     fak * S[0][12][8];
     S[2][26][8] = gma1*S[0][26][8]+     fak2* S[0][26][2];
     S[3][26][8] = gma2*S[0][26][8]+     fak3* S[0][16][8];
     S[1][26][9] = gma0*S[0][26][9]+     fak * S[0][12][9];
     S[2][26][9] = gma1*S[0][26][9];
     S[3][26][9] = gma2*S[0][26][9]+     fak * (3.e0*S[0][16][9] + 2.e0*S[0][26][3] );
     S[1][27][4] = gma0*S[0][27][4]+     fak * S[0][27][2];
     S[2][27][4] = gma1*S[0][27][4]+     fak * (3.e0*S[0][17][4] + S[0][27][1] );
     S[3][27][4] = gma2*S[0][27][4]+     fak * S[0][11][4];
     S[1][27][5] = gma0*S[0][27][5]+     fak * S[0][27][3];
     S[2][27][5] = gma1*S[0][27][5]+     fak3* S[0][17][5];
     S[3][27][5] = gma2*S[0][27][5]+     fak * (S[0][11][5] + S[0][27][1] );
     S[1][27][6] = gma0*S[0][27][6];
     S[2][27][6] = gma1*S[0][27][6]+     fak * (3.e0*S[0][17][6] + S[0][27][3] );
     S[3][27][6] = gma2*S[0][27][6]+     fak * (S[0][11][6] + S[0][27][2] );
     S[1][27][7] = gma0*S[0][27][7]+     fak2* S[0][27][1];
     S[2][27][7] = gma1*S[0][27][7]+     fak3* S[0][17][7];
     S[3][27][7] = gma2*S[0][27][7]+     fak * S[0][11][7];
     S[1][27][8] = gma0*S[0][27][8];
     S[2][27][8] = gma1*S[0][27][8]+     fak * (3.e0*S[0][17][8] + 2.e0*S[0][27][2] );
     S[3][27][8] = gma2*S[0][27][8]+     fak * S[0][11][8];
     S[1][27][9] = gma0*S[0][27][9];
     S[2][27][9] = gma1*S[0][27][9]+     fak3* S[0][17][9];
     S[3][27][9] = gma2*S[0][27][9]+     fak * (S[0][11][9] + 2.e0*S[0][27][3] );
     S[1][28][4] = gma0*S[0][28][4]+     fak * S[0][28][2];
     S[2][28][4] = gma1*S[0][28][4]+     fak * (S[0][12][4] + S[0][28][1] );
     S[3][28][4] = gma2*S[0][28][4]+     fak3* S[0][18][4];
     S[1][28][5] = gma0*S[0][28][5]+     fak * S[0][28][3];
     S[2][28][5] = gma1*S[0][28][5]+     fak * S[0][12][5];
     S[3][28][5] = gma2*S[0][28][5]+     fak * (3.e0*S[0][18][5] + S[0][28][1] );
     S[1][28][6] = gma0*S[0][28][6];
     S[2][28][6] = gma1*S[0][28][6]+     fak * (S[0][12][6] + S[0][28][3] );
     S[3][28][6] = gma2*S[0][28][6]+     fak * (3.e0*S[0][18][6] + S[0][28][2] );
     S[1][28][7] = gma0*S[0][28][7]+     fak2* S[0][28][1];
     S[2][28][7] = gma1*S[0][28][7]+     fak * S[0][12][7];
     S[3][28][7] = gma2*S[0][28][7]+     fak3* S[0][18][7];
     S[1][28][8] = gma0*S[0][28][8];
     S[2][28][8] = gma1*S[0][28][8]+     fak * (S[0][12][8] + 2.e0*S[0][28][2] );
     S[3][28][8] = gma2*S[0][28][8]+     fak3* S[0][18][8];
     S[1][28][9] = gma0*S[0][28][9];
     S[2][28][9] = gma1*S[0][28][9]+     fak * S[0][12][9];
     S[3][28][9] = gma2*S[0][28][9]+     fak * (3.e0*S[0][18][9] + 2.e0*S[0][28][3] );
     S[1][29][4] = gma0*S[0][29][4]+     fak * (2.e0*S[0][16][4] + S[0][29][2] );
     S[2][29][4] = gma1*S[0][29][4]+     fak * S[0][29][1];
     S[3][29][4] = gma2*S[0][29][4]+     fak2* S[0][15][4];
     S[1][29][5] = gma0*S[0][29][5]+     fak * (2.e0*S[0][16][5] + S[0][29][3] );
     S[2][29][5] = gma1*S[0][29][5];
     S[3][29][5] = gma2*S[0][29][5]+     fak * (2.e0*S[0][15][5] + S[0][29][1] );
     S[1][29][6] = gma0*S[0][29][6]+     fak2* S[0][16][6];
     S[2][29][6] = gma1*S[0][29][6]+     fak * S[0][29][3];
     S[3][29][6] = gma2*S[0][29][6]+     fak * (2.e0*S[0][15][6] + S[0][29][2] );
     S[1][29][7] = gma0*S[0][29][7]+     fak2* (S[0][16][7] + S[0][29][1] );
     S[2][29][7] = gma1*S[0][29][7];
     S[3][29][7] = gma2*S[0][29][7]+     fak2* S[0][15][7];
     S[1][29][8] = gma0*S[0][29][8]+     fak2* S[0][16][8];
     S[2][29][8] = gma1*S[0][29][8]+     fak2* S[0][29][2];
     S[3][29][8] = gma2*S[0][29][8]+     fak2* S[0][15][8];
     S[1][29][9] = gma0*S[0][29][9]+     fak2* S[0][16][9];
     S[2][29][9] = gma1*S[0][29][9];
     S[3][29][9] = gma2*S[0][29][9]+     fak2* (S[0][15][9] + S[0][29][3] );
     S[1][30][4] = gma0*S[0][30][4]+     fak * S[0][30][2];
     S[2][30][4] = gma1*S[0][30][4]+     fak * (2.e0*S[0][18][4] + S[0][30][1] );
     S[3][30][4] = gma2*S[0][30][4]+     fak2* S[0][17][4];
     S[1][30][5] = gma0*S[0][30][5]+     fak * S[0][30][3];
     S[2][30][5] = gma1*S[0][30][5]+     fak2* S[0][18][5];
     S[3][30][5] = gma2*S[0][30][5]+     fak * (2.e0*S[0][17][5] + S[0][30][1] );
     S[1][30][6] = gma0*S[0][30][6];
     S[2][30][6] = gma1*S[0][30][6]+     fak * (2.e0*S[0][18][6] + S[0][30][3] );
     S[3][30][6] = gma2*S[0][30][6]+     fak * (2.e0*S[0][17][6] + S[0][30][2] );
     S[1][30][7] = gma0*S[0][30][7]+     fak2* S[0][30][1];
     S[2][30][7] = gma1*S[0][30][7]+     fak2* S[0][18][7];
     S[3][30][7] = gma2*S[0][30][7]+     fak2* S[0][17][7];
     S[1][30][8] = gma0*S[0][30][8];
     S[2][30][8] = gma1*S[0][30][8]+     fak2* (S[0][18][8] + S[0][30][2] );
     S[3][30][8] = gma2*S[0][30][8]+     fak2* S[0][17][8];
     S[1][30][9] = gma0*S[0][30][9];
     S[2][30][9] = gma1*S[0][30][9]+     fak2* S[0][18][9];
     S[3][30][9] = gma2*S[0][30][9]+     fak2* (S[0][17][9] + S[0][30][3] );
     S[1][31][4] = gma0*S[0][31][4]+     fak * (2.e0*S[0][14][4] + S[0][31][2] );
     S[2][31][4] = gma1*S[0][31][4]+     fak * (2.e0*S[0][13][4] + S[0][31][1] );
     S[3][31][4] = gma2*S[0][31][4];
     S[1][31][5] = gma0*S[0][31][5]+     fak * (2.e0*S[0][14][5] + S[0][31][3] );
     S[2][31][5] = gma1*S[0][31][5]+     fak2* S[0][13][5];
     S[3][31][5] = gma2*S[0][31][5]+     fak * S[0][31][1];
     S[1][31][6] = gma0*S[0][31][6]+     fak2* S[0][14][6];
     S[2][31][6] = gma1*S[0][31][6]+     fak * (2.e0*S[0][13][6] + S[0][31][3] );
     S[3][31][6] = gma2*S[0][31][6]+     fak * S[0][31][2];
     S[1][31][7] = gma0*S[0][31][7]+     fak2* (S[0][14][7] + S[0][31][1] );
     S[2][31][7] = gma1*S[0][31][7]+     fak2* S[0][13][7];
     S[3][31][7] = gma2*S[0][31][7];
     S[1][31][8] = gma0*S[0][31][8]+     fak2* S[0][14][8];
     S[2][31][8] = gma1*S[0][31][8]+     fak2* (S[0][13][8] + S[0][31][2] );
     S[3][31][8] = gma2*S[0][31][8];
     S[1][31][9] = gma0*S[0][31][9]+     fak2* S[0][14][9];
     S[2][31][9] = gma1*S[0][31][9]+     fak2* S[0][13][9];
     S[3][31][9] = gma2*S[0][31][9]+     fak2* S[0][31][3];
     S[1][32][4] = gma0*S[0][32][4]+     fak * (2.e0*S[0][19][4] + S[0][32][2] );
     S[2][32][4] = gma1*S[0][32][4]+     fak * (S[0][15][4] + S[0][32][1] );
     S[3][32][4] = gma2*S[0][32][4]+     fak * S[0][13][4];
     S[1][32][5] = gma0*S[0][32][5]+     fak * (2.e0*S[0][19][5] + S[0][32][3] );
     S[2][32][5] = gma1*S[0][32][5]+     fak * S[0][15][5];
     S[3][32][5] = gma2*S[0][32][5]+     fak * (S[0][13][5] + S[0][32][1] );
     S[1][32][6] = gma0*S[0][32][6]+     fak2* S[0][19][6];
     S[2][32][6] = gma1*S[0][32][6]+     fak * (S[0][15][6] + S[0][32][3] );
     S[3][32][6] = gma2*S[0][32][6]+     fak * (S[0][13][6] + S[0][32][2] );
     S[1][32][7] = gma0*S[0][32][7]+     fak2* (S[0][19][7] + S[0][32][1] );
     S[2][32][7] = gma1*S[0][32][7]+     fak * S[0][15][7];
     S[3][32][7] = gma2*S[0][32][7]+     fak * S[0][13][7];
     S[1][32][8] = gma0*S[0][32][8]+     fak2* S[0][19][8];
     S[2][32][8] = gma1*S[0][32][8]+     fak * (S[0][15][8] + 2.e0*S[0][32][2] );
     S[3][32][8] = gma2*S[0][32][8]+     fak * S[0][13][8];
     S[1][32][9] = gma0*S[0][32][9]+     fak2* S[0][19][9];
     S[2][32][9] = gma1*S[0][32][9]+     fak * S[0][15][9];
     S[3][32][9] = gma2*S[0][32][9]+     fak * (S[0][13][9] + 2.e0*S[0][32][3] );
     S[1][33][4] = gma0*S[0][33][4]+     fak * (S[0][17][4] + S[0][33][2] );
     S[2][33][4] = gma1*S[0][33][4]+     fak * (2.e0*S[0][19][4] + S[0][33][1] );
     S[3][33][4] = gma2*S[0][33][4]+     fak * S[0][14][4];
     S[1][33][5] = gma0*S[0][33][5]+     fak * (S[0][17][5] + S[0][33][3] );
     S[2][33][5] = gma1*S[0][33][5]+     fak2* S[0][19][5];
     S[3][33][5] = gma2*S[0][33][5]+     fak * (S[0][14][5] + S[0][33][1] );
     S[1][33][6] = gma0*S[0][33][6]+     fak * S[0][17][6];
     S[2][33][6] = gma1*S[0][33][6]+     fak * (2.e0*S[0][19][6] + S[0][33][3] );
     S[3][33][6] = gma2*S[0][33][6]+     fak * (S[0][14][6] + S[0][33][2] );
     S[1][33][7] = gma0*S[0][33][7]+     fak * (S[0][17][7] + 2.e0*S[0][33][1] );
     S[2][33][7] = gma1*S[0][33][7]+     fak2* S[0][19][7];
     S[3][33][7] = gma2*S[0][33][7]+     fak * S[0][14][7];
     S[1][33][8] = gma0*S[0][33][8]+     fak * S[0][17][8];
     S[2][33][8] = gma1*S[0][33][8]+     fak2* (S[0][19][8] + S[0][33][2] );
     S[3][33][8] = gma2*S[0][33][8]+     fak * S[0][14][8];
     S[1][33][9] = gma0*S[0][33][9]+     fak * S[0][17][9];
     S[2][33][9] = gma1*S[0][33][9]+     fak2* S[0][19][9];
     S[3][33][9] = gma2*S[0][33][9]+     fak * (S[0][14][9] + 2.e0*S[0][33][3] );
     S[1][34][4] = gma0*S[0][34][4]+     fak * (S[0][18][4] + S[0][34][2] );
     S[2][34][4] = gma1*S[0][34][4]+     fak * (S[0][16][4] + S[0][34][1] );
     S[3][34][4] = gma2*S[0][34][4]+     fak2* S[0][19][4];
     S[1][34][5] = gma0*S[0][34][5]+     fak * (S[0][18][5] + S[0][34][3] );
     S[2][34][5] = gma1*S[0][34][5]+     fak * S[0][16][5];
     S[3][34][5] = gma2*S[0][34][5]+     fak * (2.e0*S[0][19][5] + S[0][34][1] );
     S[1][34][6] = gma0*S[0][34][6]+     fak * S[0][18][6];
     S[2][34][6] = gma1*S[0][34][6]+     fak * (S[0][16][6] + S[0][34][3] );
     S[3][34][6] = gma2*S[0][34][6]+     fak * (2.e0*S[0][19][6] + S[0][34][2] );
     S[1][34][7] = gma0*S[0][34][7]+     fak * (S[0][18][7] + 2.e0*S[0][34][1] );
     S[2][34][7] = gma1*S[0][34][7]+     fak * S[0][16][7];
     S[3][34][7] = gma2*S[0][34][7]+     fak2* S[0][19][7];
     S[1][34][8] = gma0*S[0][34][8]+     fak * S[0][18][8];
     S[2][34][8] = gma1*S[0][34][8]+     fak * (S[0][16][8] + 2.e0*S[0][34][2] );
     S[3][34][8] = gma2*S[0][34][8]+     fak2* S[0][19][8];
     S[1][34][9] = gma0*S[0][34][9]+     fak * S[0][18][9];
     S[2][34][9] = gma1*S[0][34][9]+     fak * S[0][16][9];
     S[3][34][9] = gma2*S[0][34][9]+     fak2* (S[0][19][9] + S[0][34][3] );
  }
  if ( _lmax_alpha>2  &&  _lmax_gamma>1 ) {
     S[7][20][1] = gma0*S[1][20][1]+     fak * (S[0][20][1] + 4.e0*S[1][10][1] + S[1][20][0] );
     S[4][20][1] = gma1*S[1][20][1];
     S[5][20][1] = gma2*S[1][20][1];
     S[7][20][2] = gma0*S[1][20][2]+     fak * (S[0][20][2] + 4.e0*S[1][10][2] );
     S[4][20][2] = gma1*S[1][20][2]+     fak * S[1][20][0];
     S[5][20][2] = gma2*S[1][20][2];
     S[7][20][3] = gma0*S[1][20][3]+     fak * (S[0][20][3] + 4.e0*S[1][10][3] );
     S[4][20][3] = gma1*S[1][20][3];
     S[5][20][3] = gma2*S[1][20][3]+     fak * S[1][20][0];
     S[7][21][1] = gma0*S[1][21][1]+     fak * (S[0][21][1] + S[1][21][0] );
     S[4][21][1] = gma1*S[1][21][1]+     fak4* S[1][11][1];
     S[5][21][1] = gma2*S[1][21][1];
     S[7][21][2] = gma0*S[1][21][2]+     fak * S[0][21][2];
     S[4][21][2] = gma1*S[1][21][2]+     fak * (4.e0*S[1][11][2] + S[1][21][0] );
     S[5][21][2] = gma2*S[1][21][2];
     S[7][21][3] = gma0*S[1][21][3]+     fak * S[0][21][3];
     S[4][21][3] = gma1*S[1][21][3]+     fak4* S[1][11][3];
     S[5][21][3] = gma2*S[1][21][3]+     fak * S[1][21][0];
     S[7][22][1] = gma0*S[1][22][1]+     fak * (S[0][22][1] + S[1][22][0] );
     S[4][22][1] = gma1*S[1][22][1];
     S[5][22][1] = gma2*S[1][22][1]+     fak4* S[1][12][1];
     S[7][22][2] = gma0*S[1][22][2]+     fak * S[0][22][2];
     S[4][22][2] = gma1*S[1][22][2]+     fak * S[1][22][0];
     S[5][22][2] = gma2*S[1][22][2]+     fak4* S[1][12][2];
     S[7][22][3] = gma0*S[1][22][3]+     fak * S[0][22][3];
     S[4][22][3] = gma1*S[1][22][3];
     S[5][22][3] = gma2*S[1][22][3]+     fak * (4.e0*S[1][12][3] + S[1][22][0] );
     S[7][23][1] = gma0*S[1][23][1]+     fak * (S[0][23][1] + 3.e0*S[1][13][1] + S[1][23][0] );
     S[4][23][1] = gma1*S[1][23][1]+     fak * S[1][10][1];
     S[5][23][1] = gma2*S[1][23][1];
     S[7][23][2] = gma0*S[1][23][2]+     fak * (S[0][23][2] + 3.e0*S[1][13][2] );
     S[4][23][2] = gma1*S[1][23][2]+     fak * (S[1][10][2] + S[1][23][0] );
     S[5][23][2] = gma2*S[1][23][2];
     S[7][23][3] = gma0*S[1][23][3]+     fak * (S[0][23][3] + 3.e0*S[1][13][3] );
     S[4][23][3] = gma1*S[1][23][3]+     fak * S[1][10][3];
     S[5][23][3] = gma2*S[1][23][3]+     fak * S[1][23][0];
     S[7][24][1] = gma0*S[1][24][1]+     fak * (S[0][24][1] + S[1][11][1] + S[1][24][0] );
     S[4][24][1] = gma1*S[1][24][1]+     fak3* S[1][14][1];
     S[5][24][1] = gma2*S[1][24][1];
     S[7][24][2] = gma0*S[1][24][2]+     fak * (S[0][24][2] + S[1][11][2] );
     S[4][24][2] = gma1*S[1][24][2]+     fak * (3.e0*S[1][14][2] + S[1][24][0] );
     S[5][24][2] = gma2*S[1][24][2];
     S[7][24][3] = gma0*S[1][24][3]+     fak * (S[0][24][3] + S[1][11][3] );
     S[4][24][3] = gma1*S[1][24][3]+     fak3* S[1][14][3];
     S[5][24][3] = gma2*S[1][24][3]+     fak * S[1][24][0];
     S[7][25][1] = gma0*S[1][25][1]+     fak * (S[0][25][1] + 3.e0*S[1][15][1] + S[1][25][0] );
     S[4][25][1] = gma1*S[1][25][1];
     S[5][25][1] = gma2*S[1][25][1]+     fak * S[1][10][1];
     S[7][25][2] = gma0*S[1][25][2]+     fak * (S[0][25][2] + 3.e0*S[1][15][2] );
     S[4][25][2] = gma1*S[1][25][2]+     fak * S[1][25][0];
     S[5][25][2] = gma2*S[1][25][2]+     fak * S[1][10][2];
     S[7][25][3] = gma0*S[1][25][3]+     fak * (S[0][25][3] + 3.e0*S[1][15][3] );
     S[4][25][3] = gma1*S[1][25][3];
     S[5][25][3] = gma2*S[1][25][3]+     fak * (S[1][10][3] + S[1][25][0] );
     S[7][26][1] = gma0*S[1][26][1]+     fak * (S[0][26][1] + S[1][12][1] + S[1][26][0] );
     S[4][26][1] = gma1*S[1][26][1];
     S[5][26][1] = gma2*S[1][26][1]+     fak3* S[1][16][1];
     S[7][26][2] = gma0*S[1][26][2]+     fak * (S[0][26][2] + S[1][12][2] );
     S[4][26][2] = gma1*S[1][26][2]+     fak * S[1][26][0];
     S[5][26][2] = gma2*S[1][26][2]+     fak3* S[1][16][2];
     S[7][26][3] = gma0*S[1][26][3]+     fak * (S[0][26][3] + S[1][12][3] );
     S[4][26][3] = gma1*S[1][26][3];
     S[5][26][3] = gma2*S[1][26][3]+     fak * (3.e0*S[1][16][3] + S[1][26][0] );
     S[7][27][1] = gma0*S[1][27][1]+     fak * (S[0][27][1] + S[1][27][0] );
     S[4][27][1] = gma1*S[1][27][1]+     fak3* S[1][17][1];
     S[5][27][1] = gma2*S[1][27][1]+     fak * S[1][11][1];
     S[7][27][2] = gma0*S[1][27][2]+     fak * S[0][27][2];
     S[4][27][2] = gma1*S[1][27][2]+     fak * (3.e0*S[1][17][2] + S[1][27][0] );
     S[5][27][2] = gma2*S[1][27][2]+     fak * S[1][11][2];
     S[7][27][3] = gma0*S[1][27][3]+     fak * S[0][27][3];
     S[4][27][3] = gma1*S[1][27][3]+     fak3* S[1][17][3];
     S[5][27][3] = gma2*S[1][27][3]+     fak * (S[1][11][3] + S[1][27][0] );
     S[7][28][1] = gma0*S[1][28][1]+     fak * (S[0][28][1] + S[1][28][0] );
     S[4][28][1] = gma1*S[1][28][1]+     fak * S[1][12][1];
     S[5][28][1] = gma2*S[1][28][1]+     fak3* S[1][18][1];
     S[7][28][2] = gma0*S[1][28][2]+     fak * S[0][28][2];
     S[4][28][2] = gma1*S[1][28][2]+     fak * (S[1][12][2] + S[1][28][0] );
     S[5][28][2] = gma2*S[1][28][2]+     fak3* S[1][18][2];
     S[7][28][3] = gma0*S[1][28][3]+     fak * S[0][28][3];
     S[4][28][3] = gma1*S[1][28][3]+     fak * S[1][12][3];
     S[5][28][3] = gma2*S[1][28][3]+     fak * (3.e0*S[1][18][3] + S[1][28][0] );
     S[7][29][1] = gma0*S[1][29][1]+     fak * (S[0][29][1] + 2.e0*S[1][16][1] + S[1][29][0] );
     S[4][29][1] = gma1*S[1][29][1];
     S[5][29][1] = gma2*S[1][29][1]+     fak2* S[1][15][1];
     S[7][29][2] = gma0*S[1][29][2]+     fak * (S[0][29][2] + 2.e0*S[1][16][2] );
     S[4][29][2] = gma1*S[1][29][2]+     fak * S[1][29][0];
     S[5][29][2] = gma2*S[1][29][2]+     fak2* S[1][15][2];
     S[7][29][3] = gma0*S[1][29][3]+     fak * (S[0][29][3] + 2.e0*S[1][16][3] );
     S[4][29][3] = gma1*S[1][29][3];
     S[5][29][3] = gma2*S[1][29][3]+     fak * (2.e0*S[1][15][3] + S[1][29][0] );
     S[7][30][1] = gma0*S[1][30][1]+     fak * (S[0][30][1] + S[1][30][0] );
     S[4][30][1] = gma1*S[1][30][1]+     fak2* S[1][18][1];
     S[5][30][1] = gma2*S[1][30][1]+     fak2* S[1][17][1];
     S[7][30][2] = gma0*S[1][30][2]+     fak * S[0][30][2];
     S[4][30][2] = gma1*S[1][30][2]+     fak * (2.e0*S[1][18][2] + S[1][30][0] );
     S[5][30][2] = gma2*S[1][30][2]+     fak2* S[1][17][2];
     S[7][30][3] = gma0*S[1][30][3]+     fak * S[0][30][3];
     S[4][30][3] = gma1*S[1][30][3]+     fak2* S[1][18][3];
     S[5][30][3] = gma2*S[1][30][3]+     fak * (2.e0*S[1][17][3] + S[1][30][0] );
     S[7][31][1] = gma0*S[1][31][1]+     fak * (S[0][31][1] + 2.e0*S[1][14][1] + S[1][31][0] );
     S[4][31][1] = gma1*S[1][31][1]+     fak2* S[1][13][1];
     S[5][31][1] = gma2*S[1][31][1];
     S[7][31][2] = gma0*S[1][31][2]+     fak * (S[0][31][2] + 2.e0*S[1][14][2] );
     S[4][31][2] = gma1*S[1][31][2]+     fak * (2.e0*S[1][13][2] + S[1][31][0] );
     S[5][31][2] = gma2*S[1][31][2];
     S[7][31][3] = gma0*S[1][31][3]+     fak * (S[0][31][3] + 2.e0*S[1][14][3] );
     S[4][31][3] = gma1*S[1][31][3]+     fak2* S[1][13][3];
     S[5][31][3] = gma2*S[1][31][3]+     fak * S[1][31][0];
     S[7][32][1] = gma0*S[1][32][1]+     fak * (S[0][32][1] + 2.e0*S[1][19][1] + S[1][32][0] );
     S[4][32][1] = gma1*S[1][32][1]+     fak * S[1][15][1];
     S[5][32][1] = gma2*S[1][32][1]+     fak * S[1][13][1];
     S[7][32][2] = gma0*S[1][32][2]+     fak * (S[0][32][2] + 2.e0*S[1][19][2] );
     S[4][32][2] = gma1*S[1][32][2]+     fak * (S[1][15][2] + S[1][32][0] );
     S[5][32][2] = gma2*S[1][32][2]+     fak * S[1][13][2];
     S[7][32][3] = gma0*S[1][32][3]+     fak * (S[0][32][3] + 2.e0*S[1][19][3] );
     S[4][32][3] = gma1*S[1][32][3]+     fak * S[1][15][3];
     S[5][32][3] = gma2*S[1][32][3]+     fak * (S[1][13][3] + S[1][32][0] );
     S[7][33][1] = gma0*S[1][33][1]+     fak * (S[0][33][1] + S[1][17][1] + S[1][33][0] );
     S[4][33][1] = gma1*S[1][33][1]+     fak2* S[1][19][1];
     S[5][33][1] = gma2*S[1][33][1]+     fak * S[1][14][1];
     S[7][33][2] = gma0*S[1][33][2]+     fak * (S[0][33][2] + S[1][17][2] );
     S[4][33][2] = gma1*S[1][33][2]+     fak * (2.e0*S[1][19][2] + S[1][33][0] );
     S[5][33][2] = gma2*S[1][33][2]+     fak * S[1][14][2];
     S[7][33][3] = gma0*S[1][33][3]+     fak * (S[0][33][3] + S[1][17][3] );
     S[4][33][3] = gma1*S[1][33][3]+     fak2* S[1][19][3];
     S[5][33][3] = gma2*S[1][33][3]+     fak * (S[1][14][3] + S[1][33][0] );
     S[7][34][1] = gma0*S[1][34][1]+     fak * (S[0][34][1] + S[1][18][1] + S[1][34][0] );
     S[4][34][1] = gma1*S[1][34][1]+     fak * S[1][16][1];
     S[5][34][1] = gma2*S[1][34][1]+     fak2* S[1][19][1];
     S[7][34][2] = gma0*S[1][34][2]+     fak * (S[0][34][2] + S[1][18][2] );
     S[4][34][2] = gma1*S[1][34][2]+     fak * (S[1][16][2] + S[1][34][0] );
     S[5][34][2] = gma2*S[1][34][2]+     fak2* S[1][19][2];
     S[7][34][3] = gma0*S[1][34][3]+     fak * (S[0][34][3] + S[1][18][3] );
     S[4][34][3] = gma1*S[1][34][3]+     fak * S[1][16][3];
     S[5][34][3] = gma2*S[1][34][3]+     fak * (2.e0*S[1][19][3] + S[1][34][0] );
     S[8][20][1] = gma1*S[2][20][1]+     fak * S[0][20][1];
     S[6][20][1] = gma2*S[2][20][1];
     S[8][20][2] = gma1*S[2][20][2]+     fak * (S[0][20][2] + S[2][20][0] );
     S[6][20][2] = gma2*S[2][20][2];
     S[8][20][3] = gma1*S[2][20][3]+     fak * S[0][20][3];
     S[6][20][3] = gma2*S[2][20][3]+     fak * S[2][20][0];
     S[8][21][1] = gma1*S[2][21][1]+     fak * (S[0][21][1] + 4.e0*S[2][11][1] );
     S[6][21][1] = gma2*S[2][21][1];
     S[8][21][2] = gma1*S[2][21][2]+     fak * (S[0][21][2] + 4.e0*S[2][11][2] + S[2][21][0] );
     S[6][21][2] = gma2*S[2][21][2];
     S[8][21][3] = gma1*S[2][21][3]+     fak * (S[0][21][3] + 4.e0*S[2][11][3] );
     S[6][21][3] = gma2*S[2][21][3]+     fak * S[2][21][0];
     S[8][22][1] = gma1*S[2][22][1]+     fak * S[0][22][1];
     S[6][22][1] = gma2*S[2][22][1]+     fak4* S[2][12][1];
     S[8][22][2] = gma1*S[2][22][2]+     fak * (S[0][22][2] + S[2][22][0] );
     S[6][22][2] = gma2*S[2][22][2]+     fak4* S[2][12][2];
     S[8][22][3] = gma1*S[2][22][3]+     fak * S[0][22][3];
     S[6][22][3] = gma2*S[2][22][3]+     fak * (4.e0*S[2][12][3] + S[2][22][0] );
     S[8][23][1] = gma1*S[2][23][1]+     fak * (S[0][23][1] + S[2][10][1] );
     S[6][23][1] = gma2*S[2][23][1];
     S[8][23][2] = gma1*S[2][23][2]+     fak * (S[0][23][2] + S[2][10][2] + S[2][23][0] );
     S[6][23][2] = gma2*S[2][23][2];
     S[8][23][3] = gma1*S[2][23][3]+     fak * (S[0][23][3] + S[2][10][3] );
     S[6][23][3] = gma2*S[2][23][3]+     fak * S[2][23][0];
     S[8][24][1] = gma1*S[2][24][1]+     fak * (S[0][24][1] + 3.e0*S[2][14][1] );
     S[6][24][1] = gma2*S[2][24][1];
     S[8][24][2] = gma1*S[2][24][2]+     fak * (S[0][24][2] + 3.e0*S[2][14][2] + S[2][24][0] );
     S[6][24][2] = gma2*S[2][24][2];
     S[8][24][3] = gma1*S[2][24][3]+     fak * (S[0][24][3] + 3.e0*S[2][14][3] );
     S[6][24][3] = gma2*S[2][24][3]+     fak * S[2][24][0];
     S[8][25][1] = gma1*S[2][25][1]+     fak * S[0][25][1];
     S[6][25][1] = gma2*S[2][25][1]+     fak * S[2][10][1];
     S[8][25][2] = gma1*S[2][25][2]+     fak * (S[0][25][2] + S[2][25][0] );
     S[6][25][2] = gma2*S[2][25][2]+     fak * S[2][10][2];
     S[8][25][3] = gma1*S[2][25][3]+     fak * S[0][25][3];
     S[6][25][3] = gma2*S[2][25][3]+     fak * (S[2][10][3] + S[2][25][0] );
     S[8][26][1] = gma1*S[2][26][1]+     fak * S[0][26][1];
     S[6][26][1] = gma2*S[2][26][1]+     fak3* S[2][16][1];
     S[8][26][2] = gma1*S[2][26][2]+     fak * (S[0][26][2] + S[2][26][0] );
     S[6][26][2] = gma2*S[2][26][2]+     fak3* S[2][16][2];
     S[8][26][3] = gma1*S[2][26][3]+     fak * S[0][26][3];
     S[6][26][3] = gma2*S[2][26][3]+     fak * (3.e0*S[2][16][3] + S[2][26][0] );
     S[8][27][1] = gma1*S[2][27][1]+     fak * (S[0][27][1] + 3.e0*S[2][17][1] );
     S[6][27][1] = gma2*S[2][27][1]+     fak * S[2][11][1];
     S[8][27][2] = gma1*S[2][27][2]+     fak * (S[0][27][2] + 3.e0*S[2][17][2] + S[2][27][0] );
     S[6][27][2] = gma2*S[2][27][2]+     fak * S[2][11][2];
     S[8][27][3] = gma1*S[2][27][3]+     fak * (S[0][27][3] + 3.e0*S[2][17][3] );
     S[6][27][3] = gma2*S[2][27][3]+     fak * (S[2][11][3] + S[2][27][0] );
     S[8][28][1] = gma1*S[2][28][1]+     fak * (S[0][28][1] + S[2][12][1] );
     S[6][28][1] = gma2*S[2][28][1]+     fak3* S[2][18][1];
     S[8][28][2] = gma1*S[2][28][2]+     fak * (S[0][28][2] + S[2][12][2] + S[2][28][0] );
     S[6][28][2] = gma2*S[2][28][2]+     fak3* S[2][18][2];
     S[8][28][3] = gma1*S[2][28][3]+     fak * (S[0][28][3] + S[2][12][3] );
     S[6][28][3] = gma2*S[2][28][3]+     fak * (3.e0*S[2][18][3] + S[2][28][0] );
     S[8][29][1] = gma1*S[2][29][1]+     fak * S[0][29][1];
     S[6][29][1] = gma2*S[2][29][1]+     fak2* S[2][15][1];
     S[8][29][2] = gma1*S[2][29][2]+     fak * (S[0][29][2] + S[2][29][0] );
     S[6][29][2] = gma2*S[2][29][2]+     fak2* S[2][15][2];
     S[8][29][3] = gma1*S[2][29][3]+     fak * S[0][29][3];
     S[6][29][3] = gma2*S[2][29][3]+     fak * (2.e0*S[2][15][3] + S[2][29][0] );
     S[8][30][1] = gma1*S[2][30][1]+     fak * (S[0][30][1] + 2.e0*S[2][18][1] );
     S[6][30][1] = gma2*S[2][30][1]+     fak2* S[2][17][1];
     S[8][30][2] = gma1*S[2][30][2]+     fak * (S[0][30][2] + 2.e0*S[2][18][2] + S[2][30][0] );
     S[6][30][2] = gma2*S[2][30][2]+     fak2* S[2][17][2];
     S[8][30][3] = gma1*S[2][30][3]+     fak * (S[0][30][3] + 2.e0*S[2][18][3] );
     S[6][30][3] = gma2*S[2][30][3]+     fak * (2.e0*S[2][17][3] + S[2][30][0] );
     S[8][31][1] = gma1*S[2][31][1]+     fak * (S[0][31][1] + 2.e0*S[2][13][1] );
     S[6][31][1] = gma2*S[2][31][1];
     S[8][31][2] = gma1*S[2][31][2]+     fak * (S[0][31][2] + 2.e0*S[2][13][2] + S[2][31][0] );
     S[6][31][2] = gma2*S[2][31][2];
     S[8][31][3] = gma1*S[2][31][3]+     fak * (S[0][31][3] + 2.e0*S[2][13][3] );
     S[6][31][3] = gma2*S[2][31][3]+     fak * S[2][31][0];
     S[8][32][1] = gma1*S[2][32][1]+     fak * (S[0][32][1] + S[2][15][1] );
     S[6][32][1] = gma2*S[2][32][1]+     fak * S[2][13][1];
     S[8][32][2] = gma1*S[2][32][2]+     fak * (S[0][32][2] + S[2][15][2] + S[2][32][0] );
     S[6][32][2] = gma2*S[2][32][2]+     fak * S[2][13][2];
     S[8][32][3] = gma1*S[2][32][3]+     fak * (S[0][32][3] + S[2][15][3] );
     S[6][32][3] = gma2*S[2][32][3]+     fak * (S[2][13][3] + S[2][32][0] );
     S[8][33][1] = gma1*S[2][33][1]+     fak * (S[0][33][1] + 2.e0*S[2][19][1] );
     S[6][33][1] = gma2*S[2][33][1]+     fak * S[2][14][1];
     S[8][33][2] = gma1*S[2][33][2]+     fak * (S[0][33][2] + 2.e0*S[2][19][2] + S[2][33][0] );
     S[6][33][2] = gma2*S[2][33][2]+     fak * S[2][14][2];
     S[8][33][3] = gma1*S[2][33][3]+     fak * (S[0][33][3] + 2.e0*S[2][19][3] );
     S[6][33][3] = gma2*S[2][33][3]+     fak * (S[2][14][3] + S[2][33][0] );
     S[8][34][1] = gma1*S[2][34][1]+     fak * (S[0][34][1] + S[2][16][1] );
     S[6][34][1] = gma2*S[2][34][1]+     fak2* S[2][19][1];
     S[8][34][2] = gma1*S[2][34][2]+     fak * (S[0][34][2] + S[2][16][2] + S[2][34][0] );
     S[6][34][2] = gma2*S[2][34][2]+     fak2* S[2][19][2];
     S[8][34][3] = gma1*S[2][34][3]+     fak * (S[0][34][3] + S[2][16][3] );
     S[6][34][3] = gma2*S[2][34][3]+     fak * (2.e0*S[2][19][3] + S[2][34][0] );
     S[9][20][1] = gma2*S[3][20][1]+     fak * S[0][20][1];
     S[9][20][2] = gma2*S[3][20][2]+     fak * S[0][20][2];
     S[9][20][3] = gma2*S[3][20][3]+     fak * (S[0][20][3] + S[3][20][0] );
     S[9][21][1] = gma2*S[3][21][1]+     fak * S[0][21][1];
     S[9][21][2] = gma2*S[3][21][2]+     fak * S[0][21][2];
     S[9][21][3] = gma2*S[3][21][3]+     fak * (S[0][21][3] + S[3][21][0] );
     S[9][22][1] = gma2*S[3][22][1]+     fak * (S[0][22][1] + 4.e0*S[3][12][1] );
     S[9][22][2] = gma2*S[3][22][2]+     fak * (S[0][22][2] + 4.e0*S[3][12][2] );
     S[9][22][3] = gma2*S[3][22][3]+     fak * (S[0][22][3] + 4.e0*S[3][12][3] + S[3][22][0] );
     S[9][23][1] = gma2*S[3][23][1]+     fak * S[0][23][1];
     S[9][23][2] = gma2*S[3][23][2]+     fak * S[0][23][2];
     S[9][23][3] = gma2*S[3][23][3]+     fak * (S[0][23][3] + S[3][23][0] );
     S[9][24][1] = gma2*S[3][24][1]+     fak * S[0][24][1];
     S[9][24][2] = gma2*S[3][24][2]+     fak * S[0][24][2];
     S[9][24][3] = gma2*S[3][24][3]+     fak * (S[0][24][3] + S[3][24][0] );
     S[9][25][1] = gma2*S[3][25][1]+     fak * (S[0][25][1] + S[3][10][1] );
     S[9][25][2] = gma2*S[3][25][2]+     fak * (S[0][25][2] + S[3][10][2] );
     S[9][25][3] = gma2*S[3][25][3]+     fak * (S[0][25][3] + S[3][10][3] + S[3][25][0] );
     S[9][26][1] = gma2*S[3][26][1]+     fak * (S[0][26][1] + 3.e0*S[3][16][1] );
     S[9][26][2] = gma2*S[3][26][2]+     fak * (S[0][26][2] + 3.e0*S[3][16][2] );
     S[9][26][3] = gma2*S[3][26][3]+     fak * (S[0][26][3] + 3.e0*S[3][16][3] + S[3][26][0] );
     S[9][27][1] = gma2*S[3][27][1]+     fak * (S[0][27][1] + S[3][11][1] );
     S[9][27][2] = gma2*S[3][27][2]+     fak * (S[0][27][2] + S[3][11][2] );
     S[9][27][3] = gma2*S[3][27][3]+     fak * (S[0][27][3] + S[3][11][3] + S[3][27][0] );
     S[9][28][1] = gma2*S[3][28][1]+     fak * (S[0][28][1] + 3.e0*S[3][18][1] );
     S[9][28][2] = gma2*S[3][28][2]+     fak * (S[0][28][2] + 3.e0*S[3][18][2] );
     S[9][28][3] = gma2*S[3][28][3]+     fak * (S[0][28][3] + 3.e0*S[3][18][3] + S[3][28][0] );
     S[9][29][1] = gma2*S[3][29][1]+     fak * (S[0][29][1] + 2.e0*S[3][15][1] );
     S[9][29][2] = gma2*S[3][29][2]+     fak * (S[0][29][2] + 2.e0*S[3][15][2] );
     S[9][29][3] = gma2*S[3][29][3]+     fak * (S[0][29][3] + 2.e0*S[3][15][3] + S[3][29][0] );
     S[9][30][1] = gma2*S[3][30][1]+     fak * (S[0][30][1] + 2.e0*S[3][17][1] );
     S[9][30][2] = gma2*S[3][30][2]+     fak * (S[0][30][2] + 2.e0*S[3][17][2] );
     S[9][30][3] = gma2*S[3][30][3]+     fak * (S[0][30][3] + 2.e0*S[3][17][3] + S[3][30][0] );
     S[9][31][1] = gma2*S[3][31][1]+     fak * S[0][31][1];
     S[9][31][2] = gma2*S[3][31][2]+     fak * S[0][31][2];
     S[9][31][3] = gma2*S[3][31][3]+     fak * (S[0][31][3] + S[3][31][0] );
     S[9][32][1] = gma2*S[3][32][1]+     fak * (S[0][32][1] + S[3][13][1] );
     S[9][32][2] = gma2*S[3][32][2]+     fak * (S[0][32][2] + S[3][13][2] );
     S[9][32][3] = gma2*S[3][32][3]+     fak * (S[0][32][3] + S[3][13][3] + S[3][32][0] );
     S[9][33][1] = gma2*S[3][33][1]+     fak * (S[0][33][1] + S[3][14][1] );
     S[9][33][2] = gma2*S[3][33][2]+     fak * (S[0][33][2] + S[3][14][2] );
     S[9][33][3] = gma2*S[3][33][3]+     fak * (S[0][33][3] + S[3][14][3] + S[3][33][0] );
     S[9][34][1] = gma2*S[3][34][1]+     fak * (S[0][34][1] + 2.e0*S[3][19][1] );
     S[9][34][2] = gma2*S[3][34][2]+     fak * (S[0][34][2] + 2.e0*S[3][19][2] );
     S[9][34][3] = gma2*S[3][34][3]+     fak * (S[0][34][3] + 2.e0*S[3][19][3] + S[3][34][0] );
  }
  if ( _lmax_alpha>2  &&  _lmax_gamma>2 ) {
     S[7][20][4] = gma0*S[1][20][4]+     fak * (S[0][20][4] + 4.e0*S[1][10][4] + S[1][20][2] );
     S[4][20][4] = gma1*S[1][20][4]+     fak * S[1][20][1];
     S[5][20][4] = gma2*S[1][20][4];
     S[7][20][5] = gma0*S[1][20][5]+     fak * (S[0][20][5] + 4.e0*S[1][10][5] + S[1][20][3] );
     S[4][20][5] = gma1*S[1][20][5];
     S[5][20][5] = gma2*S[1][20][5]+     fak * S[1][20][1];
     S[7][20][6] = gma0*S[1][20][6]+     fak * (S[0][20][6] + 4.e0*S[1][10][6] );
     S[4][20][6] = gma1*S[1][20][6]+     fak * S[1][20][3];
     S[5][20][6] = gma2*S[1][20][6]+     fak * S[1][20][2];
     S[7][20][7] = gma0*S[1][20][7]+     fak * (S[0][20][7] + 4.e0*S[1][10][7] + 2.e0*S[1][20][1] );
     S[4][20][7] = gma1*S[1][20][7];
     S[5][20][7] = gma2*S[1][20][7];
     S[7][20][8] = gma0*S[1][20][8]+     fak * (S[0][20][8] + 4.e0*S[1][10][8] );
     S[4][20][8] = gma1*S[1][20][8]+     fak2* S[1][20][2];
     S[5][20][8] = gma2*S[1][20][8];
     S[7][20][9] = gma0*S[1][20][9]+     fak * (S[0][20][9] + 4.e0*S[1][10][9] );
     S[4][20][9] = gma1*S[1][20][9];
     S[5][20][9] = gma2*S[1][20][9]+     fak2* S[1][20][3];
     S[7][21][4] = gma0*S[1][21][4]+     fak * (S[0][21][4] + S[1][21][2] );
     S[4][21][4] = gma1*S[1][21][4]+     fak * (4.e0*S[1][11][4] + S[1][21][1] );
     S[5][21][4] = gma2*S[1][21][4];
     S[7][21][5] = gma0*S[1][21][5]+     fak * (S[0][21][5] + S[1][21][3] );
     S[4][21][5] = gma1*S[1][21][5]+     fak4* S[1][11][5];
     S[5][21][5] = gma2*S[1][21][5]+     fak * S[1][21][1];
     S[7][21][6] = gma0*S[1][21][6]+     fak * S[0][21][6];
     S[4][21][6] = gma1*S[1][21][6]+     fak * (4.e0*S[1][11][6] + S[1][21][3] );
     S[5][21][6] = gma2*S[1][21][6]+     fak * S[1][21][2];
     S[7][21][7] = gma0*S[1][21][7]+     fak * (S[0][21][7] + 2.e0*S[1][21][1] );
     S[4][21][7] = gma1*S[1][21][7]+     fak4* S[1][11][7];
     S[5][21][7] = gma2*S[1][21][7];
     S[7][21][8] = gma0*S[1][21][8]+     fak * S[0][21][8];
     S[4][21][8] = gma1*S[1][21][8]+     fak2* (2.e0*S[1][11][8] + S[1][21][2] );
     S[5][21][8] = gma2*S[1][21][8];
     S[7][21][9] = gma0*S[1][21][9]+     fak * S[0][21][9];
     S[4][21][9] = gma1*S[1][21][9]+     fak4* S[1][11][9];
     S[5][21][9] = gma2*S[1][21][9]+     fak2* S[1][21][3];
     S[7][22][4] = gma0*S[1][22][4]+     fak * (S[0][22][4] + S[1][22][2] );
     S[4][22][4] = gma1*S[1][22][4]+     fak * S[1][22][1];
     S[5][22][4] = gma2*S[1][22][4]+     fak4* S[1][12][4];
     S[7][22][5] = gma0*S[1][22][5]+     fak * (S[0][22][5] + S[1][22][3] );
     S[4][22][5] = gma1*S[1][22][5];
     S[5][22][5] = gma2*S[1][22][5]+     fak * (4.e0*S[1][12][5] + S[1][22][1] );
     S[7][22][6] = gma0*S[1][22][6]+     fak * S[0][22][6];
     S[4][22][6] = gma1*S[1][22][6]+     fak * S[1][22][3];
     S[5][22][6] = gma2*S[1][22][6]+     fak * (4.e0*S[1][12][6] + S[1][22][2] );
     S[7][22][7] = gma0*S[1][22][7]+     fak * (S[0][22][7] + 2.e0*S[1][22][1] );
     S[4][22][7] = gma1*S[1][22][7];
     S[5][22][7] = gma2*S[1][22][7]+     fak4* S[1][12][7];
     S[7][22][8] = gma0*S[1][22][8]+     fak * S[0][22][8];
     S[4][22][8] = gma1*S[1][22][8]+     fak2* S[1][22][2];
     S[5][22][8] = gma2*S[1][22][8]+     fak4* S[1][12][8];
     S[7][22][9] = gma0*S[1][22][9]+     fak * S[0][22][9];
     S[4][22][9] = gma1*S[1][22][9];
     S[5][22][9] = gma2*S[1][22][9]+     fak2* (2.e0*S[1][12][9] + S[1][22][3] );
     S[7][23][4] = gma0*S[1][23][4]+     fak * (S[0][23][4] + 3.e0*S[1][13][4] + S[1][23][2] );
     S[4][23][4] = gma1*S[1][23][4]+     fak * (S[1][10][4] + S[1][23][1] );
     S[5][23][4] = gma2*S[1][23][4];
     S[7][23][5] = gma0*S[1][23][5]+     fak * (S[0][23][5] + 3.e0*S[1][13][5] + S[1][23][3] );
     S[4][23][5] = gma1*S[1][23][5]+     fak * S[1][10][5];
     S[5][23][5] = gma2*S[1][23][5]+     fak * S[1][23][1];
     S[7][23][6] = gma0*S[1][23][6]+     fak * (S[0][23][6] + 3.e0*S[1][13][6] );
     S[4][23][6] = gma1*S[1][23][6]+     fak * (S[1][10][6] + S[1][23][3] );
     S[5][23][6] = gma2*S[1][23][6]+     fak * S[1][23][2];
     S[7][23][7] = gma0*S[1][23][7]+     fak * (S[0][23][7] + 3.e0*S[1][13][7] + 2.e0*S[1][23][1] );
     S[4][23][7] = gma1*S[1][23][7]+     fak * S[1][10][7];
     S[5][23][7] = gma2*S[1][23][7];
     S[7][23][8] = gma0*S[1][23][8]+     fak * (S[0][23][8] + 3.e0*S[1][13][8] );
     S[4][23][8] = gma1*S[1][23][8]+     fak * (S[1][10][8] + 2.e0*S[1][23][2] );
     S[5][23][8] = gma2*S[1][23][8];
     S[7][23][9] = gma0*S[1][23][9]+     fak * (S[0][23][9] + 3.e0*S[1][13][9] );
     S[4][23][9] = gma1*S[1][23][9]+     fak * S[1][10][9];
     S[5][23][9] = gma2*S[1][23][9]+     fak2* S[1][23][3];
     S[7][24][4] = gma0*S[1][24][4]+     fak * (S[0][24][4] + S[1][11][4] + S[1][24][2] );
     S[4][24][4] = gma1*S[1][24][4]+     fak * (3.e0*S[1][14][4] + S[1][24][1] );
     S[5][24][4] = gma2*S[1][24][4];
     S[7][24][5] = gma0*S[1][24][5]+     fak * (S[0][24][5] + S[1][11][5] + S[1][24][3] );
     S[4][24][5] = gma1*S[1][24][5]+     fak3* S[1][14][5];
     S[5][24][5] = gma2*S[1][24][5]+     fak * S[1][24][1];
     S[7][24][6] = gma0*S[1][24][6]+     fak * (S[0][24][6] + S[1][11][6] );
     S[4][24][6] = gma1*S[1][24][6]+     fak * (3.e0*S[1][14][6] + S[1][24][3] );
     S[5][24][6] = gma2*S[1][24][6]+     fak * S[1][24][2];
     S[7][24][7] = gma0*S[1][24][7]+     fak * (S[0][24][7] + S[1][11][7] + 2.e0*S[1][24][1] );
     S[4][24][7] = gma1*S[1][24][7]+     fak3* S[1][14][7];
     S[5][24][7] = gma2*S[1][24][7];
     S[7][24][8] = gma0*S[1][24][8]+     fak * (S[0][24][8] + S[1][11][8] );
     S[4][24][8] = gma1*S[1][24][8]+     fak * (3.e0*S[1][14][8] + 2.e0*S[1][24][2] );
     S[5][24][8] = gma2*S[1][24][8];
     S[7][24][9] = gma0*S[1][24][9]+     fak * (S[0][24][9] + S[1][11][9] );
     S[4][24][9] = gma1*S[1][24][9]+     fak3* S[1][14][9];
     S[5][24][9] = gma2*S[1][24][9]+     fak2* S[1][24][3];
     S[7][25][4] = gma0*S[1][25][4]+     fak * (S[0][25][4] + 3.e0*S[1][15][4] + S[1][25][2] );
     S[4][25][4] = gma1*S[1][25][4]+     fak * S[1][25][1];
     S[5][25][4] = gma2*S[1][25][4]+     fak * S[1][10][4];
     S[7][25][5] = gma0*S[1][25][5]+     fak * (S[0][25][5] + 3.e0*S[1][15][5] + S[1][25][3] );
     S[4][25][5] = gma1*S[1][25][5];
     S[5][25][5] = gma2*S[1][25][5]+     fak * (S[1][10][5] + S[1][25][1] );
     S[7][25][6] = gma0*S[1][25][6]+     fak * (S[0][25][6] + 3.e0*S[1][15][6] );
     S[4][25][6] = gma1*S[1][25][6]+     fak * S[1][25][3];
     S[5][25][6] = gma2*S[1][25][6]+     fak * (S[1][10][6] + S[1][25][2] );
     S[7][25][7] = gma0*S[1][25][7]+     fak * (S[0][25][7] + 3.e0*S[1][15][7] + 2.e0*S[1][25][1] );
     S[4][25][7] = gma1*S[1][25][7];
     S[5][25][7] = gma2*S[1][25][7]+     fak * S[1][10][7];
     S[7][25][8] = gma0*S[1][25][8]+     fak * (S[0][25][8] + 3.e0*S[1][15][8] );
     S[4][25][8] = gma1*S[1][25][8]+     fak2* S[1][25][2];
     S[5][25][8] = gma2*S[1][25][8]+     fak * S[1][10][8];
     S[7][25][9] = gma0*S[1][25][9]+     fak * (S[0][25][9] + 3.e0*S[1][15][9] );
     S[4][25][9] = gma1*S[1][25][9];
     S[5][25][9] = gma2*S[1][25][9]+     fak * (S[1][10][9] + 2.e0*S[1][25][3] );
     S[7][26][4] = gma0*S[1][26][4]+     fak * (S[0][26][4] + S[1][12][4] + S[1][26][2] );
     S[4][26][4] = gma1*S[1][26][4]+     fak * S[1][26][1];
     S[5][26][4] = gma2*S[1][26][4]+     fak3* S[1][16][4];
     S[7][26][5] = gma0*S[1][26][5]+     fak * (S[0][26][5] + S[1][12][5] + S[1][26][3] );
     S[4][26][5] = gma1*S[1][26][5];
     S[5][26][5] = gma2*S[1][26][5]+     fak * (3.e0*S[1][16][5] + S[1][26][1] );
     S[7][26][6] = gma0*S[1][26][6]+     fak * (S[0][26][6] + S[1][12][6] );
     S[4][26][6] = gma1*S[1][26][6]+     fak * S[1][26][3];
     S[5][26][6] = gma2*S[1][26][6]+     fak * (3.e0*S[1][16][6] + S[1][26][2] );
     S[7][26][7] = gma0*S[1][26][7]+     fak * (S[0][26][7] + S[1][12][7] + 2.e0*S[1][26][1] );
     S[4][26][7] = gma1*S[1][26][7];
     S[5][26][7] = gma2*S[1][26][7]+     fak3* S[1][16][7];
     S[7][26][8] = gma0*S[1][26][8]+     fak * (S[0][26][8] + S[1][12][8] );
     S[4][26][8] = gma1*S[1][26][8]+     fak2* S[1][26][2];
     S[5][26][8] = gma2*S[1][26][8]+     fak3* S[1][16][8];
     S[7][26][9] = gma0*S[1][26][9]+     fak * (S[0][26][9] + S[1][12][9] );
     S[4][26][9] = gma1*S[1][26][9];
     S[5][26][9] = gma2*S[1][26][9]+     fak * (3.e0*S[1][16][9] + 2.e0*S[1][26][3] );
     S[7][27][4] = gma0*S[1][27][4]+     fak * (S[0][27][4] + S[1][27][2] );
     S[4][27][4] = gma1*S[1][27][4]+     fak * (3.e0*S[1][17][4] + S[1][27][1] );
     S[5][27][4] = gma2*S[1][27][4]+     fak * S[1][11][4];
     S[7][27][5] = gma0*S[1][27][5]+     fak * (S[0][27][5] + S[1][27][3] );
     S[4][27][5] = gma1*S[1][27][5]+     fak3* S[1][17][5];
     S[5][27][5] = gma2*S[1][27][5]+     fak * (S[1][11][5] + S[1][27][1] );
     S[7][27][6] = gma0*S[1][27][6]+     fak * S[0][27][6];
     S[4][27][6] = gma1*S[1][27][6]+     fak * (3.e0*S[1][17][6] + S[1][27][3] );
     S[5][27][6] = gma2*S[1][27][6]+     fak * (S[1][11][6] + S[1][27][2] );
     S[7][27][7] = gma0*S[1][27][7]+     fak * (S[0][27][7] + 2.e0*S[1][27][1] );
     S[4][27][7] = gma1*S[1][27][7]+     fak3* S[1][17][7];
     S[5][27][7] = gma2*S[1][27][7]+     fak * S[1][11][7];
     S[7][27][8] = gma0*S[1][27][8]+     fak * S[0][27][8];
     S[4][27][8] = gma1*S[1][27][8]+     fak * (3.e0*S[1][17][8] + 2.e0*S[1][27][2] );
     S[5][27][8] = gma2*S[1][27][8]+     fak * S[1][11][8];
     S[7][27][9] = gma0*S[1][27][9]+     fak * S[0][27][9];
     S[4][27][9] = gma1*S[1][27][9]+     fak3* S[1][17][9];
     S[5][27][9] = gma2*S[1][27][9]+     fak * (S[1][11][9] + 2.e0*S[1][27][3] );
     S[7][28][4] = gma0*S[1][28][4]+     fak * (S[0][28][4] + S[1][28][2] );
     S[4][28][4] = gma1*S[1][28][4]+     fak * (S[1][12][4] + S[1][28][1] );
     S[5][28][4] = gma2*S[1][28][4]+     fak3* S[1][18][4];
     S[7][28][5] = gma0*S[1][28][5]+     fak * (S[0][28][5] + S[1][28][3] );
     S[4][28][5] = gma1*S[1][28][5]+     fak * S[1][12][5];
     S[5][28][5] = gma2*S[1][28][5]+     fak * (3.e0*S[1][18][5] + S[1][28][1] );
     S[7][28][6] = gma0*S[1][28][6]+     fak * S[0][28][6];
     S[4][28][6] = gma1*S[1][28][6]+     fak * (S[1][12][6] + S[1][28][3] );
     S[5][28][6] = gma2*S[1][28][6]+     fak * (3.e0*S[1][18][6] + S[1][28][2] );
     S[7][28][7] = gma0*S[1][28][7]+     fak * (S[0][28][7] + 2.e0*S[1][28][1] );
     S[4][28][7] = gma1*S[1][28][7]+     fak * S[1][12][7];
     S[5][28][7] = gma2*S[1][28][7]+     fak3* S[1][18][7];
     S[7][28][8] = gma0*S[1][28][8]+     fak * S[0][28][8];
     S[4][28][8] = gma1*S[1][28][8]+     fak * (S[1][12][8] + 2.e0*S[1][28][2] );
     S[5][28][8] = gma2*S[1][28][8]+     fak3* S[1][18][8];
     S[7][28][9] = gma0*S[1][28][9]+     fak * S[0][28][9];
     S[4][28][9] = gma1*S[1][28][9]+     fak * S[1][12][9];
     S[5][28][9] = gma2*S[1][28][9]+     fak * (3.e0*S[1][18][9] + 2.e0*S[1][28][3] );
     S[7][29][4] = gma0*S[1][29][4]+     fak * (S[0][29][4] + 2.e0*S[1][16][4] + S[1][29][2] );
     S[4][29][4] = gma1*S[1][29][4]+     fak * S[1][29][1];
     S[5][29][4] = gma2*S[1][29][4]+     fak2* S[1][15][4];
     S[7][29][5] = gma0*S[1][29][5]+     fak * (S[0][29][5] + 2.e0*S[1][16][5] + S[1][29][3] );
     S[4][29][5] = gma1*S[1][29][5];
     S[5][29][5] = gma2*S[1][29][5]+     fak * (2.e0*S[1][15][5] + S[1][29][1] );
     S[7][29][6] = gma0*S[1][29][6]+     fak * (S[0][29][6] + 2.e0*S[1][16][6] );
     S[4][29][6] = gma1*S[1][29][6]+     fak * S[1][29][3];
     S[5][29][6] = gma2*S[1][29][6]+     fak * (2.e0*S[1][15][6] + S[1][29][2] );
     S[7][29][7] = gma0*S[1][29][7]+     fak * (S[0][29][7] + 2.e0*S[1][16][7] + 2.e0*S[1][29][1] );
     S[4][29][7] = gma1*S[1][29][7];
     S[5][29][7] = gma2*S[1][29][7]+     fak2* S[1][15][7];
     S[7][29][8] = gma0*S[1][29][8]+     fak * (S[0][29][8] + 2.e0*S[1][16][8] );
     S[4][29][8] = gma1*S[1][29][8]+     fak2* S[1][29][2];
     S[5][29][8] = gma2*S[1][29][8]+     fak2* S[1][15][8];
     S[7][29][9] = gma0*S[1][29][9]+     fak * (S[0][29][9] + 2.e0*S[1][16][9] );
     S[4][29][9] = gma1*S[1][29][9];
     S[5][29][9] = gma2*S[1][29][9]+     fak2* (S[1][15][9] + S[1][29][3] );
     S[7][30][4] = gma0*S[1][30][4]+     fak * (S[0][30][4] + S[1][30][2] );
     S[4][30][4] = gma1*S[1][30][4]+     fak * (2.e0*S[1][18][4] + S[1][30][1] );
     S[5][30][4] = gma2*S[1][30][4]+     fak2* S[1][17][4];
     S[7][30][5] = gma0*S[1][30][5]+     fak * (S[0][30][5] + S[1][30][3] );
     S[4][30][5] = gma1*S[1][30][5]+     fak2* S[1][18][5];
     S[5][30][5] = gma2*S[1][30][5]+     fak * (2.e0*S[1][17][5] + S[1][30][1] );
     S[7][30][6] = gma0*S[1][30][6]+     fak * S[0][30][6];
     S[4][30][6] = gma1*S[1][30][6]+     fak * (2.e0*S[1][18][6] + S[1][30][3] );
     S[5][30][6] = gma2*S[1][30][6]+     fak * (2.e0*S[1][17][6] + S[1][30][2] );
     S[7][30][7] = gma0*S[1][30][7]+     fak * (S[0][30][7] + 2.e0*S[1][30][1] );
     S[4][30][7] = gma1*S[1][30][7]+     fak2* S[1][18][7];
     S[5][30][7] = gma2*S[1][30][7]+     fak2* S[1][17][7];
     S[7][30][8] = gma0*S[1][30][8]+     fak * S[0][30][8];
     S[4][30][8] = gma1*S[1][30][8]+     fak2* (S[1][18][8] + S[1][30][2] );
     S[5][30][8] = gma2*S[1][30][8]+     fak2* S[1][17][8];
     S[7][30][9] = gma0*S[1][30][9]+     fak * S[0][30][9];
     S[4][30][9] = gma1*S[1][30][9]+     fak2* S[1][18][9];
     S[5][30][9] = gma2*S[1][30][9]+     fak2* (S[1][17][9] + S[1][30][3] );
     S[7][31][4] = gma0*S[1][31][4]+     fak * (S[0][31][4] + 2.e0*S[1][14][4] + S[1][31][2] );
     S[4][31][4] = gma1*S[1][31][4]+     fak * (2.e0*S[1][13][4] + S[1][31][1] );
     S[5][31][4] = gma2*S[1][31][4];
     S[7][31][5] = gma0*S[1][31][5]+     fak * (S[0][31][5] + 2.e0*S[1][14][5] + S[1][31][3] );
     S[4][31][5] = gma1*S[1][31][5]+     fak2* S[1][13][5];
     S[5][31][5] = gma2*S[1][31][5]+     fak * S[1][31][1];
     S[7][31][6] = gma0*S[1][31][6]+     fak * (S[0][31][6] + 2.e0*S[1][14][6] );
     S[4][31][6] = gma1*S[1][31][6]+     fak * (2.e0*S[1][13][6] + S[1][31][3] );
     S[5][31][6] = gma2*S[1][31][6]+     fak * S[1][31][2];
     S[7][31][7] = gma0*S[1][31][7]+     fak * (S[0][31][7] + 2.e0*S[1][14][7] + 2.e0*S[1][31][1] );
     S[4][31][7] = gma1*S[1][31][7]+     fak2* S[1][13][7];
     S[5][31][7] = gma2*S[1][31][7];
     S[7][31][8] = gma0*S[1][31][8]+     fak * (S[0][31][8] + 2.e0*S[1][14][8] );
     S[4][31][8] = gma1*S[1][31][8]+     fak2* (S[1][13][8] + S[1][31][2] );
     S[5][31][8] = gma2*S[1][31][8];
     S[7][31][9] = gma0*S[1][31][9]+     fak * (S[0][31][9] + 2.e0*S[1][14][9] );
     S[4][31][9] = gma1*S[1][31][9]+     fak2* S[1][13][9];
     S[5][31][9] = gma2*S[1][31][9]+     fak2* S[1][31][3];
     S[7][32][4] = gma0*S[1][32][4]+     fak * (S[0][32][4] + 2.e0*S[1][19][4] + S[1][32][2] );
     S[4][32][4] = gma1*S[1][32][4]+     fak * (S[1][15][4] + S[1][32][1] );
     S[5][32][4] = gma2*S[1][32][4]+     fak * S[1][13][4];
     S[7][32][5] = gma0*S[1][32][5]+     fak * (S[0][32][5] + 2.e0*S[1][19][5] + S[1][32][3] );
     S[4][32][5] = gma1*S[1][32][5]+     fak * S[1][15][5];
     S[5][32][5] = gma2*S[1][32][5]+     fak * (S[1][13][5] + S[1][32][1] );
     S[7][32][6] = gma0*S[1][32][6]+     fak * (S[0][32][6] + 2.e0*S[1][19][6] );
     S[4][32][6] = gma1*S[1][32][6]+     fak * (S[1][15][6] + S[1][32][3] );
     S[5][32][6] = gma2*S[1][32][6]+     fak * (S[1][13][6] + S[1][32][2] );
     S[7][32][7] = gma0*S[1][32][7]+     fak * (S[0][32][7] + 2.e0*S[1][19][7] + 2.e0*S[1][32][1] );
     S[4][32][7] = gma1*S[1][32][7]+     fak * S[1][15][7];
     S[5][32][7] = gma2*S[1][32][7]+     fak * S[1][13][7];
     S[7][32][8] = gma0*S[1][32][8]+     fak * (S[0][32][8] + 2.e0*S[1][19][8] );
     S[4][32][8] = gma1*S[1][32][8]+     fak * (S[1][15][8] + 2.e0*S[1][32][2] );
     S[5][32][8] = gma2*S[1][32][8]+     fak * S[1][13][8];
     S[7][32][9] = gma0*S[1][32][9]+     fak * (S[0][32][9] + 2.e0*S[1][19][9] );
     S[4][32][9] = gma1*S[1][32][9]+     fak * S[1][15][9];
     S[5][32][9] = gma2*S[1][32][9]+     fak * (S[1][13][9] + 2.e0*S[1][32][3] );
     S[7][33][4] = gma0*S[1][33][4]+     fak * (S[0][33][4] + S[1][17][4] + S[1][33][2] );
     S[4][33][4] = gma1*S[1][33][4]+     fak * (2.e0*S[1][19][4] + S[1][33][1] );
     S[5][33][4] = gma2*S[1][33][4]+     fak * S[1][14][4];
     S[7][33][5] = gma0*S[1][33][5]+     fak * (S[0][33][5] + S[1][17][5] + S[1][33][3] );
     S[4][33][5] = gma1*S[1][33][5]+     fak2* S[1][19][5];
     S[5][33][5] = gma2*S[1][33][5]+     fak * (S[1][14][5] + S[1][33][1] );
     S[7][33][6] = gma0*S[1][33][6]+     fak * (S[0][33][6] + S[1][17][6] );
     S[4][33][6] = gma1*S[1][33][6]+     fak * (2.e0*S[1][19][6] + S[1][33][3] );
     S[5][33][6] = gma2*S[1][33][6]+     fak * (S[1][14][6] + S[1][33][2] );
     S[7][33][7] = gma0*S[1][33][7]+     fak * (S[0][33][7] + S[1][17][7] + 2.e0*S[1][33][1] );
     S[4][33][7] = gma1*S[1][33][7]+     fak2* S[1][19][7];
     S[5][33][7] = gma2*S[1][33][7]+     fak * S[1][14][7];
     S[7][33][8] = gma0*S[1][33][8]+     fak * (S[0][33][8] + S[1][17][8] );
     S[4][33][8] = gma1*S[1][33][8]+     fak2* (S[1][19][8] + S[1][33][2] );
     S[5][33][8] = gma2*S[1][33][8]+     fak * S[1][14][8];
     S[7][33][9] = gma0*S[1][33][9]+     fak * (S[0][33][9] + S[1][17][9] );
     S[4][33][9] = gma1*S[1][33][9]+     fak2* S[1][19][9];
     S[5][33][9] = gma2*S[1][33][9]+     fak * (S[1][14][9] + 2.e0*S[1][33][3] );
     S[7][34][4] = gma0*S[1][34][4]+     fak * (S[0][34][4] + S[1][18][4] + S[1][34][2] );
     S[4][34][4] = gma1*S[1][34][4]+     fak * (S[1][16][4] + S[1][34][1] );
     S[5][34][4] = gma2*S[1][34][4]+     fak2* S[1][19][4];
     S[7][34][5] = gma0*S[1][34][5]+     fak * (S[0][34][5] + S[1][18][5] + S[1][34][3] );
     S[4][34][5] = gma1*S[1][34][5]+     fak * S[1][16][5];
     S[5][34][5] = gma2*S[1][34][5]+     fak * (2.e0*S[1][19][5] + S[1][34][1] );
     S[7][34][6] = gma0*S[1][34][6]+     fak * (S[0][34][6] + S[1][18][6] );
     S[4][34][6] = gma1*S[1][34][6]+     fak * (S[1][16][6] + S[1][34][3] );
     S[5][34][6] = gma2*S[1][34][6]+     fak * (2.e0*S[1][19][6] + S[1][34][2] );
     S[7][34][7] = gma0*S[1][34][7]+     fak * (S[0][34][7] + S[1][18][7] + 2.e0*S[1][34][1] );
     S[4][34][7] = gma1*S[1][34][7]+     fak * S[1][16][7];
     S[5][34][7] = gma2*S[1][34][7]+     fak2* S[1][19][7];
     S[7][34][8] = gma0*S[1][34][8]+     fak * (S[0][34][8] + S[1][18][8] );
     S[4][34][8] = gma1*S[1][34][8]+     fak * (S[1][16][8] + 2.e0*S[1][34][2] );
     S[5][34][8] = gma2*S[1][34][8]+     fak2* S[1][19][8];
     S[7][34][9] = gma0*S[1][34][9]+     fak * (S[0][34][9] + S[1][18][9] );
     S[4][34][9] = gma1*S[1][34][9]+     fak * S[1][16][9];
     S[5][34][9] = gma2*S[1][34][9]+     fak2* (S[1][19][9] + S[1][34][3] );
     S[8][20][4] = gma1*S[2][20][4]+     fak * (S[0][20][4] + S[2][20][1] );
     S[6][20][4] = gma2*S[2][20][4];
     S[8][20][5] = gma1*S[2][20][5]+     fak * S[0][20][5];
     S[6][20][5] = gma2*S[2][20][5]+     fak * S[2][20][1];
     S[8][20][6] = gma1*S[2][20][6]+     fak * (S[0][20][6] + S[2][20][3] );
     S[6][20][6] = gma2*S[2][20][6]+     fak * S[2][20][2];
     S[8][20][7] = gma1*S[2][20][7]+     fak * S[0][20][7];
     S[6][20][7] = gma2*S[2][20][7];
     S[8][20][8] = gma1*S[2][20][8]+     fak * (S[0][20][8] + 2.e0*S[2][20][2] );
     S[6][20][8] = gma2*S[2][20][8];
     S[8][20][9] = gma1*S[2][20][9]+     fak * S[0][20][9];
     S[6][20][9] = gma2*S[2][20][9]+     fak2* S[2][20][3];
     S[8][21][4] = gma1*S[2][21][4]+     fak * (S[0][21][4] + 4.e0*S[2][11][4] + S[2][21][1] );
     S[6][21][4] = gma2*S[2][21][4];
     S[8][21][5] = gma1*S[2][21][5]+     fak * (S[0][21][5] + 4.e0*S[2][11][5] );
     S[6][21][5] = gma2*S[2][21][5]+     fak * S[2][21][1];
     S[8][21][6] = gma1*S[2][21][6]+     fak * (S[0][21][6] + 4.e0*S[2][11][6] + S[2][21][3] );
     S[6][21][6] = gma2*S[2][21][6]+     fak * S[2][21][2];
     S[8][21][7] = gma1*S[2][21][7]+     fak * (S[0][21][7] + 4.e0*S[2][11][7] );
     S[6][21][7] = gma2*S[2][21][7];
     S[8][21][8] = gma1*S[2][21][8]+     fak * (S[0][21][8] + 4.e0*S[2][11][8] + 2.e0*S[2][21][2] );
     S[6][21][8] = gma2*S[2][21][8];
     S[8][21][9] = gma1*S[2][21][9]+     fak * (S[0][21][9] + 4.e0*S[2][11][9] );
     S[6][21][9] = gma2*S[2][21][9]+     fak2* S[2][21][3];
     S[8][22][4] = gma1*S[2][22][4]+     fak * (S[0][22][4] + S[2][22][1] );
     S[6][22][4] = gma2*S[2][22][4]+     fak4* S[2][12][4];
     S[8][22][5] = gma1*S[2][22][5]+     fak * S[0][22][5];
     S[6][22][5] = gma2*S[2][22][5]+     fak * (4.e0*S[2][12][5] + S[2][22][1] );
     S[8][22][6] = gma1*S[2][22][6]+     fak * (S[0][22][6] + S[2][22][3] );
     S[6][22][6] = gma2*S[2][22][6]+     fak * (4.e0*S[2][12][6] + S[2][22][2] );
     S[8][22][7] = gma1*S[2][22][7]+     fak * S[0][22][7];
     S[6][22][7] = gma2*S[2][22][7]+     fak4* S[2][12][7];
     S[8][22][8] = gma1*S[2][22][8]+     fak * (S[0][22][8] + 2.e0*S[2][22][2] );
     S[6][22][8] = gma2*S[2][22][8]+     fak4* S[2][12][8];
     S[8][22][9] = gma1*S[2][22][9]+     fak * S[0][22][9];
     S[6][22][9] = gma2*S[2][22][9]+     fak2* (2.e0*S[2][12][9] + S[2][22][3] );
     S[8][23][4] = gma1*S[2][23][4]+     fak * (S[0][23][4] + S[2][10][4] + S[2][23][1] );
     S[6][23][4] = gma2*S[2][23][4];
     S[8][23][5] = gma1*S[2][23][5]+     fak * (S[0][23][5] + S[2][10][5] );
     S[6][23][5] = gma2*S[2][23][5]+     fak * S[2][23][1];
     S[8][23][6] = gma1*S[2][23][6]+     fak * (S[0][23][6] + S[2][10][6] + S[2][23][3] );
     S[6][23][6] = gma2*S[2][23][6]+     fak * S[2][23][2];
     S[8][23][7] = gma1*S[2][23][7]+     fak * (S[0][23][7] + S[2][10][7] );
     S[6][23][7] = gma2*S[2][23][7];
     S[8][23][8] = gma1*S[2][23][8]+     fak * (S[0][23][8] + S[2][10][8] + 2.e0*S[2][23][2] );
     S[6][23][8] = gma2*S[2][23][8];
     S[8][23][9] = gma1*S[2][23][9]+     fak * (S[0][23][9] + S[2][10][9] );
     S[6][23][9] = gma2*S[2][23][9]+     fak2* S[2][23][3];
     S[8][24][4] = gma1*S[2][24][4]+     fak * (S[0][24][4] + 3.e0*S[2][14][4] + S[2][24][1] );
     S[6][24][4] = gma2*S[2][24][4];
     S[8][24][5] = gma1*S[2][24][5]+     fak * (S[0][24][5] + 3.e0*S[2][14][5] );
     S[6][24][5] = gma2*S[2][24][5]+     fak * S[2][24][1];
     S[8][24][6] = gma1*S[2][24][6]+     fak * (S[0][24][6] + 3.e0*S[2][14][6] + S[2][24][3] );
     S[6][24][6] = gma2*S[2][24][6]+     fak * S[2][24][2];
     S[8][24][7] = gma1*S[2][24][7]+     fak * (S[0][24][7] + 3.e0*S[2][14][7] );
     S[6][24][7] = gma2*S[2][24][7];
     S[8][24][8] = gma1*S[2][24][8]+     fak * (S[0][24][8] + 3.e0*S[2][14][8] + 2.e0*S[2][24][2] );
     S[6][24][8] = gma2*S[2][24][8];
     S[8][24][9] = gma1*S[2][24][9]+     fak * (S[0][24][9] + 3.e0*S[2][14][9] );
     S[6][24][9] = gma2*S[2][24][9]+     fak2* S[2][24][3];
     S[8][25][4] = gma1*S[2][25][4]+     fak * (S[0][25][4] + S[2][25][1] );
     S[6][25][4] = gma2*S[2][25][4]+     fak * S[2][10][4];
     S[8][25][5] = gma1*S[2][25][5]+     fak * S[0][25][5];
     S[6][25][5] = gma2*S[2][25][5]+     fak * (S[2][10][5] + S[2][25][1] );
     S[8][25][6] = gma1*S[2][25][6]+     fak * (S[0][25][6] + S[2][25][3] );
     S[6][25][6] = gma2*S[2][25][6]+     fak * (S[2][10][6] + S[2][25][2] );
     S[8][25][7] = gma1*S[2][25][7]+     fak * S[0][25][7];
     S[6][25][7] = gma2*S[2][25][7]+     fak * S[2][10][7];
     S[8][25][8] = gma1*S[2][25][8]+     fak * (S[0][25][8] + 2.e0*S[2][25][2] );
     S[6][25][8] = gma2*S[2][25][8]+     fak * S[2][10][8];
     S[8][25][9] = gma1*S[2][25][9]+     fak * S[0][25][9];
     S[6][25][9] = gma2*S[2][25][9]+     fak * (S[2][10][9] + 2.e0*S[2][25][3] );
     S[8][26][4] = gma1*S[2][26][4]+     fak * (S[0][26][4] + S[2][26][1] );
     S[6][26][4] = gma2*S[2][26][4]+     fak3* S[2][16][4];
     S[8][26][5] = gma1*S[2][26][5]+     fak * S[0][26][5];
     S[6][26][5] = gma2*S[2][26][5]+     fak * (3.e0*S[2][16][5] + S[2][26][1] );
     S[8][26][6] = gma1*S[2][26][6]+     fak * (S[0][26][6] + S[2][26][3] );
     S[6][26][6] = gma2*S[2][26][6]+     fak * (3.e0*S[2][16][6] + S[2][26][2] );
     S[8][26][7] = gma1*S[2][26][7]+     fak * S[0][26][7];
     S[6][26][7] = gma2*S[2][26][7]+     fak3* S[2][16][7];
     S[8][26][8] = gma1*S[2][26][8]+     fak * (S[0][26][8] + 2.e0*S[2][26][2] );
     S[6][26][8] = gma2*S[2][26][8]+     fak3* S[2][16][8];
     S[8][26][9] = gma1*S[2][26][9]+     fak * S[0][26][9];
     S[6][26][9] = gma2*S[2][26][9]+     fak * (3.e0*S[2][16][9] + 2.e0*S[2][26][3] );
     S[8][27][4] = gma1*S[2][27][4]+     fak * (S[0][27][4] + 3.e0*S[2][17][4] + S[2][27][1] );
     S[6][27][4] = gma2*S[2][27][4]+     fak * S[2][11][4];
     S[8][27][5] = gma1*S[2][27][5]+     fak * (S[0][27][5] + 3.e0*S[2][17][5] );
     S[6][27][5] = gma2*S[2][27][5]+     fak * (S[2][11][5] + S[2][27][1] );
     S[8][27][6] = gma1*S[2][27][6]+     fak * (S[0][27][6] + 3.e0*S[2][17][6] + S[2][27][3] );
     S[6][27][6] = gma2*S[2][27][6]+     fak * (S[2][11][6] + S[2][27][2] );
     S[8][27][7] = gma1*S[2][27][7]+     fak * (S[0][27][7] + 3.e0*S[2][17][7] );
     S[6][27][7] = gma2*S[2][27][7]+     fak * S[2][11][7];
     S[8][27][8] = gma1*S[2][27][8]+     fak * (S[0][27][8] + 3.e0*S[2][17][8] + 2.e0*S[2][27][2] );
     S[6][27][8] = gma2*S[2][27][8]+     fak * S[2][11][8];
     S[8][27][9] = gma1*S[2][27][9]+     fak * (S[0][27][9] + 3.e0*S[2][17][9] );
     S[6][27][9] = gma2*S[2][27][9]+     fak * (S[2][11][9] + 2.e0*S[2][27][3] );
     S[8][28][4] = gma1*S[2][28][4]+     fak * (S[0][28][4] + S[2][12][4] + S[2][28][1] );
     S[6][28][4] = gma2*S[2][28][4]+     fak3* S[2][18][4];
     S[8][28][5] = gma1*S[2][28][5]+     fak * (S[0][28][5] + S[2][12][5] );
     S[6][28][5] = gma2*S[2][28][5]+     fak * (3.e0*S[2][18][5] + S[2][28][1] );
     S[8][28][6] = gma1*S[2][28][6]+     fak * (S[0][28][6] + S[2][12][6] + S[2][28][3] );
     S[6][28][6] = gma2*S[2][28][6]+     fak * (3.e0*S[2][18][6] + S[2][28][2] );
     S[8][28][7] = gma1*S[2][28][7]+     fak * (S[0][28][7] + S[2][12][7] );
     S[6][28][7] = gma2*S[2][28][7]+     fak3* S[2][18][7];
     S[8][28][8] = gma1*S[2][28][8]+     fak * (S[0][28][8] + S[2][12][8] + 2.e0*S[2][28][2] );
     S[6][28][8] = gma2*S[2][28][8]+     fak3* S[2][18][8];
     S[8][28][9] = gma1*S[2][28][9]+     fak * (S[0][28][9] + S[2][12][9] );
     S[6][28][9] = gma2*S[2][28][9]+     fak * (3.e0*S[2][18][9] + 2.e0*S[2][28][3] );
     S[8][29][4] = gma1*S[2][29][4]+     fak * (S[0][29][4] + S[2][29][1] );
     S[6][29][4] = gma2*S[2][29][4]+     fak2* S[2][15][4];
     S[8][29][5] = gma1*S[2][29][5]+     fak * S[0][29][5];
     S[6][29][5] = gma2*S[2][29][5]+     fak * (2.e0*S[2][15][5] + S[2][29][1] );
     S[8][29][6] = gma1*S[2][29][6]+     fak * (S[0][29][6] + S[2][29][3] );
     S[6][29][6] = gma2*S[2][29][6]+     fak * (2.e0*S[2][15][6] + S[2][29][2] );
     S[8][29][7] = gma1*S[2][29][7]+     fak * S[0][29][7];
     S[6][29][7] = gma2*S[2][29][7]+     fak2* S[2][15][7];
     S[8][29][8] = gma1*S[2][29][8]+     fak * (S[0][29][8] + 2.e0*S[2][29][2] );
     S[6][29][8] = gma2*S[2][29][8]+     fak2* S[2][15][8];
     S[8][29][9] = gma1*S[2][29][9]+     fak * S[0][29][9];
     S[6][29][9] = gma2*S[2][29][9]+     fak2* (S[2][15][9] + S[2][29][3] );
     S[8][30][4] = gma1*S[2][30][4]+     fak * (S[0][30][4] + 2.e0*S[2][18][4] + S[2][30][1] );
     S[6][30][4] = gma2*S[2][30][4]+     fak2* S[2][17][4];
     S[8][30][5] = gma1*S[2][30][5]+     fak * (S[0][30][5] + 2.e0*S[2][18][5] );
     S[6][30][5] = gma2*S[2][30][5]+     fak * (2.e0*S[2][17][5] + S[2][30][1] );
     S[8][30][6] = gma1*S[2][30][6]+     fak * (S[0][30][6] + 2.e0*S[2][18][6] + S[2][30][3] );
     S[6][30][6] = gma2*S[2][30][6]+     fak * (2.e0*S[2][17][6] + S[2][30][2] );
     S[8][30][7] = gma1*S[2][30][7]+     fak * (S[0][30][7] + 2.e0*S[2][18][7] );
     S[6][30][7] = gma2*S[2][30][7]+     fak2* S[2][17][7];
     S[8][30][8] = gma1*S[2][30][8]+     fak * (S[0][30][8] + 2.e0*S[2][18][8] + 2.e0*S[2][30][2] );
     S[6][30][8] = gma2*S[2][30][8]+     fak2* S[2][17][8];
     S[8][30][9] = gma1*S[2][30][9]+     fak * (S[0][30][9] + 2.e0*S[2][18][9] );
     S[6][30][9] = gma2*S[2][30][9]+     fak2* (S[2][17][9] + S[2][30][3] );
     S[8][31][4] = gma1*S[2][31][4]+     fak * (S[0][31][4] + 2.e0*S[2][13][4] + S[2][31][1] );
     S[6][31][4] = gma2*S[2][31][4];
     S[8][31][5] = gma1*S[2][31][5]+     fak * (S[0][31][5] + 2.e0*S[2][13][5] );
     S[6][31][5] = gma2*S[2][31][5]+     fak * S[2][31][1];
     S[8][31][6] = gma1*S[2][31][6]+     fak * (S[0][31][6] + 2.e0*S[2][13][6] + S[2][31][3] );
     S[6][31][6] = gma2*S[2][31][6]+     fak * S[2][31][2];
     S[8][31][7] = gma1*S[2][31][7]+     fak * (S[0][31][7] + 2.e0*S[2][13][7] );
     S[6][31][7] = gma2*S[2][31][7];
     S[8][31][8] = gma1*S[2][31][8]+     fak * (S[0][31][8] + 2.e0*S[2][13][8] + 2.e0*S[2][31][2] );
     S[6][31][8] = gma2*S[2][31][8];
     S[8][31][9] = gma1*S[2][31][9]+     fak * (S[0][31][9] + 2.e0*S[2][13][9] );
     S[6][31][9] = gma2*S[2][31][9]+     fak2* S[2][31][3];
     S[8][32][4] = gma1*S[2][32][4]+     fak * (S[0][32][4] + S[2][15][4] + S[2][32][1] );
     S[6][32][4] = gma2*S[2][32][4]+     fak * S[2][13][4];
     S[8][32][5] = gma1*S[2][32][5]+     fak * (S[0][32][5] + S[2][15][5] );
     S[6][32][5] = gma2*S[2][32][5]+     fak * (S[2][13][5] + S[2][32][1] );
     S[8][32][6] = gma1*S[2][32][6]+     fak * (S[0][32][6] + S[2][15][6] + S[2][32][3] );
     S[6][32][6] = gma2*S[2][32][6]+     fak * (S[2][13][6] + S[2][32][2] );
     S[8][32][7] = gma1*S[2][32][7]+     fak * (S[0][32][7] + S[2][15][7] );
     S[6][32][7] = gma2*S[2][32][7]+     fak * S[2][13][7];
     S[8][32][8] = gma1*S[2][32][8]+     fak * (S[0][32][8] + S[2][15][8] + 2.e0*S[2][32][2] );
     S[6][32][8] = gma2*S[2][32][8]+     fak * S[2][13][8];
     S[8][32][9] = gma1*S[2][32][9]+     fak * (S[0][32][9] + S[2][15][9] );
     S[6][32][9] = gma2*S[2][32][9]+     fak * (S[2][13][9] + 2.e0*S[2][32][3] );
     S[8][33][4] = gma1*S[2][33][4]+     fak * (S[0][33][4] + 2.e0*S[2][19][4] + S[2][33][1] );
     S[6][33][4] = gma2*S[2][33][4]+     fak * S[2][14][4];
     S[8][33][5] = gma1*S[2][33][5]+     fak * (S[0][33][5] + 2.e0*S[2][19][5] );
     S[6][33][5] = gma2*S[2][33][5]+     fak * (S[2][14][5] + S[2][33][1] );
     S[8][33][6] = gma1*S[2][33][6]+     fak * (S[0][33][6] + 2.e0*S[2][19][6] + S[2][33][3] );
     S[6][33][6] = gma2*S[2][33][6]+     fak * (S[2][14][6] + S[2][33][2] );
     S[8][33][7] = gma1*S[2][33][7]+     fak * (S[0][33][7] + 2.e0*S[2][19][7] );
     S[6][33][7] = gma2*S[2][33][7]+     fak * S[2][14][7];
     S[8][33][8] = gma1*S[2][33][8]+     fak * (S[0][33][8] + 2.e0*S[2][19][8] + 2.e0*S[2][33][2] );
     S[6][33][8] = gma2*S[2][33][8]+     fak * S[2][14][8];
     S[8][33][9] = gma1*S[2][33][9]+     fak * (S[0][33][9] + 2.e0*S[2][19][9] );
     S[6][33][9] = gma2*S[2][33][9]+     fak * (S[2][14][9] + 2.e0*S[2][33][3] );
     S[8][34][4] = gma1*S[2][34][4]+     fak * (S[0][34][4] + S[2][16][4] + S[2][34][1] );
     S[6][34][4] = gma2*S[2][34][4]+     fak2* S[2][19][4];
     S[8][34][5] = gma1*S[2][34][5]+     fak * (S[0][34][5] + S[2][16][5] );
     S[6][34][5] = gma2*S[2][34][5]+     fak * (2.e0*S[2][19][5] + S[2][34][1] );
     S[8][34][6] = gma1*S[2][34][6]+     fak * (S[0][34][6] + S[2][16][6] + S[2][34][3] );
     S[6][34][6] = gma2*S[2][34][6]+     fak * (2.e0*S[2][19][6] + S[2][34][2] );
     S[8][34][7] = gma1*S[2][34][7]+     fak * (S[0][34][7] + S[2][16][7] );
     S[6][34][7] = gma2*S[2][34][7]+     fak2* S[2][19][7];
     S[8][34][8] = gma1*S[2][34][8]+     fak * (S[0][34][8] + S[2][16][8] + 2.e0*S[2][34][2] );
     S[6][34][8] = gma2*S[2][34][8]+     fak2* S[2][19][8];
     S[8][34][9] = gma1*S[2][34][9]+     fak * (S[0][34][9] + S[2][16][9] );
     S[6][34][9] = gma2*S[2][34][9]+     fak2* (S[2][19][9] + S[2][34][3] );
     S[9][20][4] = gma2*S[3][20][4]+     fak * S[0][20][4];
     S[9][20][5] = gma2*S[3][20][5]+     fak * (S[0][20][5] + S[3][20][1] );
     S[9][20][6] = gma2*S[3][20][6]+     fak * (S[0][20][6] + S[3][20][2] );
     S[9][20][7] = gma2*S[3][20][7]+     fak * S[0][20][7];
     S[9][20][8] = gma2*S[3][20][8]+     fak * S[0][20][8];
     S[9][20][9] = gma2*S[3][20][9]+     fak * (S[0][20][9] + 2.e0*S[3][20][3] );
     S[9][21][4] = gma2*S[3][21][4]+     fak * S[0][21][4];
     S[9][21][5] = gma2*S[3][21][5]+     fak * (S[0][21][5] + S[3][21][1] );
     S[9][21][6] = gma2*S[3][21][6]+     fak * (S[0][21][6] + S[3][21][2] );
     S[9][21][7] = gma2*S[3][21][7]+     fak * S[0][21][7];
     S[9][21][8] = gma2*S[3][21][8]+     fak * S[0][21][8];
     S[9][21][9] = gma2*S[3][21][9]+     fak * (S[0][21][9] + 2.e0*S[3][21][3] );
     S[9][22][4] = gma2*S[3][22][4]+     fak * (S[0][22][4] + 4.e0*S[3][12][4] );
     S[9][22][5] = gma2*S[3][22][5]+     fak * (S[0][22][5] + 4.e0*S[3][12][5] + S[3][22][1] );
     S[9][22][6] = gma2*S[3][22][6]+     fak * (S[0][22][6] + 4.e0*S[3][12][6] + S[3][22][2] );
     S[9][22][7] = gma2*S[3][22][7]+     fak * (S[0][22][7] + 4.e0*S[3][12][7] );
     S[9][22][8] = gma2*S[3][22][8]+     fak * (S[0][22][8] + 4.e0*S[3][12][8] );
     S[9][22][9] = gma2*S[3][22][9]+     fak * (S[0][22][9] + 4.e0*S[3][12][9] + 2.e0*S[3][22][3] );
     S[9][23][4] = gma2*S[3][23][4]+     fak * S[0][23][4];
     S[9][23][5] = gma2*S[3][23][5]+     fak * (S[0][23][5] + S[3][23][1] );
     S[9][23][6] = gma2*S[3][23][6]+     fak * (S[0][23][6] + S[3][23][2] );
     S[9][23][7] = gma2*S[3][23][7]+     fak * S[0][23][7];
     S[9][23][8] = gma2*S[3][23][8]+     fak * S[0][23][8];
     S[9][23][9] = gma2*S[3][23][9]+     fak * (S[0][23][9] + 2.e0*S[3][23][3] );
     S[9][24][4] = gma2*S[3][24][4]+     fak * S[0][24][4];
     S[9][24][5] = gma2*S[3][24][5]+     fak * (S[0][24][5] + S[3][24][1] );
     S[9][24][6] = gma2*S[3][24][6]+     fak * (S[0][24][6] + S[3][24][2] );
     S[9][24][7] = gma2*S[3][24][7]+     fak * S[0][24][7];
     S[9][24][8] = gma2*S[3][24][8]+     fak * S[0][24][8];
     S[9][24][9] = gma2*S[3][24][9]+     fak * (S[0][24][9] + 2.e0*S[3][24][3] );
     S[9][25][4] = gma2*S[3][25][4]+     fak * (S[0][25][4] + S[3][10][4] );
     S[9][25][5] = gma2*S[3][25][5]+     fak * (S[0][25][5] + S[3][10][5] + S[3][25][1] );
     S[9][25][6] = gma2*S[3][25][6]+     fak * (S[0][25][6] + S[3][10][6] + S[3][25][2] );
     S[9][25][7] = gma2*S[3][25][7]+     fak * (S[0][25][7] + S[3][10][7] );
     S[9][25][8] = gma2*S[3][25][8]+     fak * (S[0][25][8] + S[3][10][8] );
     S[9][25][9] = gma2*S[3][25][9]+     fak * (S[0][25][9] + S[3][10][9] + 2.e0*S[3][25][3] );
     S[9][26][4] = gma2*S[3][26][4]+     fak * (S[0][26][4] + 3.e0*S[3][16][4] );
     S[9][26][5] = gma2*S[3][26][5]+     fak * (S[0][26][5] + 3.e0*S[3][16][5] + S[3][26][1] );
     S[9][26][6] = gma2*S[3][26][6]+     fak * (S[0][26][6] + 3.e0*S[3][16][6] + S[3][26][2] );
     S[9][26][7] = gma2*S[3][26][7]+     fak * (S[0][26][7] + 3.e0*S[3][16][7] );
     S[9][26][8] = gma2*S[3][26][8]+     fak * (S[0][26][8] + 3.e0*S[3][16][8] );
     S[9][26][9] = gma2*S[3][26][9]+     fak * (S[0][26][9] + 3.e0*S[3][16][9] + 2.e0*S[3][26][3] );
     S[9][27][4] = gma2*S[3][27][4]+     fak * (S[0][27][4] + S[3][11][4] );
     S[9][27][5] = gma2*S[3][27][5]+     fak * (S[0][27][5] + S[3][11][5] + S[3][27][1] );
     S[9][27][6] = gma2*S[3][27][6]+     fak * (S[0][27][6] + S[3][11][6] + S[3][27][2] );
     S[9][27][7] = gma2*S[3][27][7]+     fak * (S[0][27][7] + S[3][11][7] );
     S[9][27][8] = gma2*S[3][27][8]+     fak * (S[0][27][8] + S[3][11][8] );
     S[9][27][9] = gma2*S[3][27][9]+     fak * (S[0][27][9] + S[3][11][9] + 2.e0*S[3][27][3] );
     S[9][28][4] = gma2*S[3][28][4]+     fak * (S[0][28][4] + 3.e0*S[3][18][4] );
     S[9][28][5] = gma2*S[3][28][5]+     fak * (S[0][28][5] + 3.e0*S[3][18][5] + S[3][28][1] );
     S[9][28][6] = gma2*S[3][28][6]+     fak * (S[0][28][6] + 3.e0*S[3][18][6] + S[3][28][2] );
     S[9][28][7] = gma2*S[3][28][7]+     fak * (S[0][28][7] + 3.e0*S[3][18][7] );
     S[9][28][8] = gma2*S[3][28][8]+     fak * (S[0][28][8] + 3.e0*S[3][18][8] );
     S[9][28][9] = gma2*S[3][28][9]+     fak * (S[0][28][9] + 3.e0*S[3][18][9] + 2.e0*S[3][28][3] );
     S[9][29][4] = gma2*S[3][29][4]+     fak * (S[0][29][4] + 2.e0*S[3][15][4] );
     S[9][29][5] = gma2*S[3][29][5]+     fak * (S[0][29][5] + 2.e0*S[3][15][5] + S[3][29][1] );
     S[9][29][6] = gma2*S[3][29][6]+     fak * (S[0][29][6] + 2.e0*S[3][15][6] + S[3][29][2] );
     S[9][29][7] = gma2*S[3][29][7]+     fak * (S[0][29][7] + 2.e0*S[3][15][7] );
     S[9][29][8] = gma2*S[3][29][8]+     fak * (S[0][29][8] + 2.e0*S[3][15][8] );
     S[9][29][9] = gma2*S[3][29][9]+     fak * (S[0][29][9] + 2.e0*S[3][15][9] + 2.e0*S[3][29][3] );
     S[9][30][4] = gma2*S[3][30][4]+     fak * (S[0][30][4] + 2.e0*S[3][17][4] );
     S[9][30][5] = gma2*S[3][30][5]+     fak * (S[0][30][5] + 2.e0*S[3][17][5] + S[3][30][1] );
     S[9][30][6] = gma2*S[3][30][6]+     fak * (S[0][30][6] + 2.e0*S[3][17][6] + S[3][30][2] );
     S[9][30][7] = gma2*S[3][30][7]+     fak * (S[0][30][7] + 2.e0*S[3][17][7] );
     S[9][30][8] = gma2*S[3][30][8]+     fak * (S[0][30][8] + 2.e0*S[3][17][8] );
     S[9][30][9] = gma2*S[3][30][9]+     fak * (S[0][30][9] + 2.e0*S[3][17][9] + 2.e0*S[3][30][3] );
     S[9][31][4] = gma2*S[3][31][4]+     fak * S[0][31][4];
     S[9][31][5] = gma2*S[3][31][5]+     fak * (S[0][31][5] + S[3][31][1] );
     S[9][31][6] = gma2*S[3][31][6]+     fak * (S[0][31][6] + S[3][31][2] );
     S[9][31][7] = gma2*S[3][31][7]+     fak * S[0][31][7];
     S[9][31][8] = gma2*S[3][31][8]+     fak * S[0][31][8];
     S[9][31][9] = gma2*S[3][31][9]+     fak * (S[0][31][9] + 2.e0*S[3][31][3] );
     S[9][32][4] = gma2*S[3][32][4]+     fak * (S[0][32][4] + S[3][13][4] );
     S[9][32][5] = gma2*S[3][32][5]+     fak * (S[0][32][5] + S[3][13][5] + S[3][32][1] );
     S[9][32][6] = gma2*S[3][32][6]+     fak * (S[0][32][6] + S[3][13][6] + S[3][32][2] );
     S[9][32][7] = gma2*S[3][32][7]+     fak * (S[0][32][7] + S[3][13][7] );
     S[9][32][8] = gma2*S[3][32][8]+     fak * (S[0][32][8] + S[3][13][8] );
     S[9][32][9] = gma2*S[3][32][9]+     fak * (S[0][32][9] + S[3][13][9] + 2.e0*S[3][32][3] );
     S[9][33][4] = gma2*S[3][33][4]+     fak * (S[0][33][4] + S[3][14][4] );
     S[9][33][5] = gma2*S[3][33][5]+     fak * (S[0][33][5] + S[3][14][5] + S[3][33][1] );
     S[9][33][6] = gma2*S[3][33][6]+     fak * (S[0][33][6] + S[3][14][6] + S[3][33][2] );
     S[9][33][7] = gma2*S[3][33][7]+     fak * (S[0][33][7] + S[3][14][7] );
     S[9][33][8] = gma2*S[3][33][8]+     fak * (S[0][33][8] + S[3][14][8] );
     S[9][33][9] = gma2*S[3][33][9]+     fak * (S[0][33][9] + S[3][14][9] + 2.e0*S[3][33][3] );
     S[9][34][4] = gma2*S[3][34][4]+     fak * (S[0][34][4] + 2.e0*S[3][19][4] );
     S[9][34][5] = gma2*S[3][34][5]+     fak * (S[0][34][5] + 2.e0*S[3][19][5] + S[3][34][1] );
     S[9][34][6] = gma2*S[3][34][6]+     fak * (S[0][34][6] + 2.e0*S[3][19][6] + S[3][34][2] );
     S[9][34][7] = gma2*S[3][34][7]+     fak * (S[0][34][7] + 2.e0*S[3][19][7] );
     S[9][34][8] = gma2*S[3][34][8]+     fak * (S[0][34][8] + 2.e0*S[3][19][8] );
     S[9][34][9] = gma2*S[3][34][9]+     fak * (S[0][34][9] + 2.e0*S[3][19][9] + 2.e0*S[3][34][3] );
  }
            
            
            
            
            }
            
            

            // data is now stored in unnormalized cartesian Gaussians in the multiarray
            // Now, weird-looking construction since multiarray is not accessible for ub::prod
            //              s px py pz dxz dyz dxy d3z2-r2 dx2-y2  f1  f2  f3  f4  f5  f6  f7  g1  g2  g3  g4  g5  g6  g7  g8  g9 
            int istart[] = {0, 1, 2, 3, 5,  6,  4,   7,     7,    12,  10, 11, 11, 10, 19, 15,  5, 25, 27, 23, 20, 25, 27, 23, 20}; //extend for g
            int istop[] =  {0, 1, 2, 3, 5,  6,  4,   9,     8,    17,  16, 18, 13, 14, 19, 17, 31, 33, 32 ,34, 30, 33, 32, 24, 31}; // extend for g

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

            
            std::vector<double> _contractions_alpha = (*italpha)->contraction;
            std::vector<double> _contractions_gamma = (*itgamma)->contraction;
            std::vector<double> _contractions_gw    = (*itgw)->contraction;
            
            // get transformation matrices
            this->getTrafo(_trafo_gw, _lmax_gw, _decay_gw, _contractions_gw);
            this->getTrafo(_trafo_alpha, _lmax_alpha, _decay_alpha, _contractions_alpha);
            this->getTrafo(_trafo_gamma, _lmax_gamma, _decay_gamma, _contractions_gamma);

            // transform from unnormalized cartesians to normalized sphericals
            // container with indices starting at zero
            ma_type S_sph;
            S_sph.resize(extents[ _ntrafo_alpha ][ _ntrafo_gw ][ _ntrafo_gamma ]);

            for (int _i_gw = 0; _i_gw < _ntrafo_gw; _i_gw++) {
                for (int _i_alpha = 0; _i_alpha < _ntrafo_alpha; _i_alpha++) {
                    for (int _i_gamma = 0; _i_gamma < _ntrafo_gamma; _i_gamma++) {

                        S_sph[ _i_alpha ][ _i_gw ][ _i_gamma ] = 0.0;

                        for (int _i_gw_t = istart[ _i_gw ]; _i_gw_t <= istop[ _i_gw ]; _i_gw_t++) {
                            for (int _i_alpha_t = istart[ _i_alpha ]; _i_alpha_t <= istop[ _i_alpha ]; _i_alpha_t++) {
                                for (int _i_gamma_t = istart[ _i_gamma ]; _i_gamma_t <= istop[ _i_gamma ]; _i_gamma_t++) {

                                    S_sph[ _i_alpha ][ _i_gw ][ _i_gamma ] += S[ _i_alpha_t + 1 ][ _i_gw_t + 1 ][ _i_gamma_t + 1]
                                            * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_gw(_i_gw, _i_gw_t) * _trafo_gamma(_i_gamma, _i_gamma_t);


                                }
                            }
                        }
                    }
                }
            }

            // only store the parts, we need
            for (int _i_gw = 0; _i_gw < _shell_gw->getNumFunc(); _i_gw++) {
                for (int _i_alpha = 0; _i_alpha < _shell_alpha->getNumFunc(); _i_alpha++) {
                    for (int _i_gamma = 0; _i_gamma < _shell_gamma->getNumFunc(); _i_gamma++) {


                        int _i_index = _shell_gamma->getNumFunc() * _i_gw + _i_gamma;

                        _subvector(_i_alpha, _i_index) += S_sph[ _offset_alpha + _i_alpha ][ _offset_gw + _i_gw ][ _offset_gamma + _i_gamma ];

                    }
                }
            }

                    }
                }
            }

            return _does_contribute;


        }

     

        void TCrawMatrix::getTrafo(ub::matrix<double>& _trafo, int _lmax, const double& _decay,std::vector<double> contractions) {
        // s-functions
        _trafo(0,0) = 1.0*contractions[0]; // s
       
        // p-functions
        if ( _lmax > 0 ){
            //cout << _trafo_row.size1() << ":" << _trafo_row.size2() << endl;
            _trafo(1,1) = 2.0*sqrt(_decay)*contractions[1];
            _trafo(2,2) = 2.0*sqrt(_decay)*contractions[1];
            _trafo(3,3) = 2.0*sqrt(_decay)*contractions[1];
        }

        // d-functions
        if ( _lmax > 1 ){
            _trafo(4,5) = 4.0*_decay*contractions[2];             // dxz
            _trafo(5,6) = _trafo(4,5);            // dyz
            _trafo(6,4) = _trafo(4,5);            // dxy
            _trafo(7,7) = -2.0*_decay/sqrt(3.0)*contractions[2];  // d3z2-r2 (dxx)
            _trafo(7,8) = _trafo(7,7);            // d3z2-r2 (dyy)
            _trafo(7,9) = -2.0*_trafo(7,7);       // d3z2-r2 (dzz)
            _trafo(8,7) = 2.0*_decay*contractions[2];             // dx2-y2 (dxx)
            _trafo(8,8) = -_trafo(8,7);           // dx2-y2 (dzz)
        }
        
        // f-functions
        if ( _lmax > 2 ){
            _trafo(9,12) = 4.0 * 2.0 *pow(_decay,1.5)/sqrt(15.)*contractions[3]; // f1 (f??)
            _trafo(9,15) = -1.5 * _trafo(9,12);        // f1 (f??)
            _trafo(9,17) = _trafo(9,15);               // f1 (f??)
            
            _trafo(10,16) = 4.0 * 2.0 * sqrt(2.0)/sqrt(5.0) * pow(_decay,1.5)*contractions[3]; // f2 (f??)
            _trafo(10,10) = -0.25 * _trafo(10,16);                             // f2 f(??)
            _trafo(10,14) = _trafo(10,10);                                     // f2 f(??)
            
            _trafo(11,18) = _trafo(10,16);                                     // f3 (f??)
            _trafo(11,13) = -0.25 * _trafo(11,18);                             // f3 f(??)
            _trafo(11,11) = _trafo(11,13);                                     // f3 f(??)            
                   
            _trafo(12,13) = 3.0 * 2.0 * sqrt(2.0)/sqrt(3.0) * pow(_decay,1.5)*contractions[3]; // f4 (f??)
            _trafo(12,11) = -_trafo(12,13)/3.0;                                // f4 (f??)
            
            _trafo(13,10) = -_trafo(12,11);                                    // f5 (f??)
            _trafo(13,14) = -_trafo(12,13);                                    // f5 (f??)
            
            _trafo(14,19) = 8.0 * pow(_decay,1.5)*contractions[3];                             // f6 (f??)
            
            _trafo(15,15) = 0.5 * _trafo(14,19);                               // f7 (f??)
            _trafo(15,17) = -_trafo(15,15);                                    // f7 (f??)
        }
        
        // g-functions
        if ( _lmax > 3 ){
            //g1
            _trafo(16,22) = 8.0 * 2.0/sqrt(105.0) * pow(_decay,2.0)*contractions[4];
            _trafo(16,21) = 3.0 * 2.0/sqrt(105.0) * pow(_decay,2.0)*contractions[4];
            _trafo(16,20) = _trafo(16,21);
            _trafo(16,29) = -3.0 * _trafo(16,22);
            _trafo(16,31) = 2.0 * _trafo(16,21);
            _trafo(16,30) = _trafo(16,29);
            _trafo(16,5)  = _trafo(16,31);
            
             /* vv(17,:) =  (/   23,  22, 21, 30, 32, 31,   6 /) ! g
                cc(17,:) =  (/    8,  3, 3, -24, 6, -24,    6 /)
                normConst(17,:) = (/ 2.d0/sqrt(105.d0) ,2.d0  /)
              */
            //g2
            _trafo(17,26) = 4.0 * 4.0*sqrt(2.0)/sqrt(21.0) * pow(_decay,2.0)*contractions[4];
            _trafo(17,25) = -0.75 * _trafo(17,26);
            _trafo(17,33) = _trafo(17,25);
             
             /* vv(18,:) =  (/   27,  26, 34,  0,  0,  0,   3 /) ! g
                cc(18,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
                normConst(18,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
              */
            //g3
            _trafo(18,28) = _trafo(17,26);
            _trafo(18,32) = _trafo(17,25);
            _trafo(18,27) = _trafo(17,25);
             
            /* vv(19,:) =  (/   29,  33, 28,  0,  0,  0,   3 /) ! g 
               cc(19,:) =  (/    4,  -3, -3,  0,  0,  0,   3 /)
               normConst(19,:) = (/ 4.d0*sqrt(2.d0)/sqrt(21.d0) ,2.d0  /)
             */
            //g4
            _trafo(19,34) = 6.0 * 8.0/sqrt(21.0) * pow(_decay,2.0)*contractions[4];
            _trafo(19,23) = -_trafo(19,34)/6.0;
            _trafo(19,24) = _trafo(19,23);
             
            /* vv(20,:) =  (/   35,  24, 25,  0,  0,  0,   3 /) ! g
               cc(20,:) =  (/    6,  -1, -1,  0,  0,  0,   3 /)
               normConst(20,:) = (/ 8.d0/sqrt(21.d0) ,2.d0  /)
             */
            //g5
            _trafo(20,29) = 6.0 * 4.0/sqrt(21.0) * pow(_decay,2.0)*contractions[4];
            _trafo(20,20) = -_trafo(20,29)/6.0;
            _trafo(20,30) = -_trafo(20,29);
            _trafo(20,21) = -_trafo(20,20);

            /* vv(21,:) =  (/   30,  21, 31, 22,  0,  0,   4 /) ! g
               cc(21,:) =  (/    6,  -1, -6, 1,  0,  0,    4 /)
               normConst(21,:) = (/ 4.d0/sqrt(21.d0) ,2.d0  /)
             */
            //g6
            _trafo(21,25) = 4.0 * sqrt(2.0)/sqrt(3.0) * pow(_decay,2.0)*contractions[4];
            _trafo(21,33) = -3.0 * _trafo(21,25);
             
            /* vv(22,:) =  (/   26,  34,  0,  0,  0,  0,   2 /) ! g
               cc(22,:) =  (/    1,  -3,  0,  0,  0,  0,   2 /)
               normConst(22,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
             */
            //g7
            _trafo(22,32) = -_trafo(21,33);
            _trafo(22,27) = -_trafo(21,25);
            
            /* vv(23,:) =  (/   33,  28,  0,  0,  0,  0,   2 /) ! g
               cc(23,:) =  (/    3,  -1,  0,  0,  0,  0,   2 /)
               normConst(23,:) = (/ 4.d0*sqrt(2.d0)/sqrt(3.d0) ,2.d0  /)
             */
            //g8
            _trafo(23,23) = 8.0/sqrt(3.0) * pow(_decay,2.0)*contractions[4];
            _trafo(23,24) = -_trafo(23,23);
             
            /* vv(24,:) =  (/   24,  25,  0,  0,  0,  0,   2 /) ! g 
               cc(24,:) =  (/    1,  -1,  0,  0,  0,  0,   2 /)
               normConst(24,:) = (/ 8.d0/sqrt(3.d0) ,2.d0  /)
             */
            //g9
            _trafo(24,20) = 2.0/sqrt(3.0) * pow(_decay,2.0)*contractions[4];
            _trafo(24,21) = _trafo(24,20);
            _trafo(24,31) = -6.0 * _trafo(24,20);
             
            /* vv(25,:) =  (/   21,  22, 32,  0,  0,  0,   3 /) ! g
               cc(25,:) =  (/    1,  1, -6,  0,  0,  0,   3  /)
               normConst(25,:) = (/ 2.d0/sqrt(3.d0) ,2.d0  /)
             */
           
       
       }

        }

        int TCrawMatrix::getBlockSize(int _lmax) {
            int _block_size;
            if (_lmax == 0) {
                _block_size = 1;
            } // s
            if (_lmax == 1) {
                _block_size = 4;
            } // p
            if (_lmax == 2) {
                _block_size = 10;
            } // d
            if (_lmax == 3) {
                _block_size = 20;
            } // f
            if (_lmax == 4) {
                _block_size = 35;
            } // g

            return _block_size;
        }


    }
}

