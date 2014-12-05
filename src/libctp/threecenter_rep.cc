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
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/threecenters.h>

#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/linalg.h>

using namespace std;
using namespace votca::tools;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;

 
        
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
        
      
        bool TCrawMatrix::FillThreeCenterRepBlock(ub::matrix<double>& _subvector,  AOShell* _shell_1, AOShell* _shell_2, AOShell* _shell_3) {
            
            const double pi = boost::math::constants::pi<double>();
            
          
            
            
            // shell info, only lmax tells how far to go
            
            int _lmax_1 = _shell_1->getLmax();
            int _lmax_2 = _shell_2->getLmax();
            int _lmax_3 = _shell_3->getLmax();
            
            int _mmax = _lmax_1+_lmax_2+_lmax_3;

            // set size of internal block for recursion
           
            AOShell* _shell_alpha;
            AOShell* _shell_beta;
            AOShell* _shell_gamma;
            bool alphabetaswitch=false;

            // We need lmax_alpha > lmax_beta, so instead of calculating (sp,s) we calculate (ps,s), due to symmetry they are the same. 
            
            if (_lmax_1 < _lmax_2){
                _shell_alpha=_shell_2;
                _shell_beta =_shell_1;
                alphabetaswitch=true;
                
            }
            else{
                _shell_alpha=_shell_1;
                _shell_beta =_shell_2;
            }
            _shell_gamma=_shell_3;
            
            const vec& _pos_alpha = _shell_alpha->getPos();
            const vec& _pos_beta = _shell_beta->getPos();
            const vec& _pos_gamma = _shell_gamma->getPos();
            
            int _lmax_alpha = _shell_alpha->getLmax();
            int _lmax_beta  = _shell_beta->getLmax();
            int _lmax_gamma = _shell_gamma->getLmax();
            
            int _ngamma = this->getBlockSize(_lmax_gamma);
            int _nalpha = this->getBlockSize(_lmax_alpha);
            int _nbeta = this->getBlockSize(_lmax_beta);
            
            //start vertical recurrence
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;

            for ( GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha){

                const double& _decay_alpha = (*italpha)->decay;
            
                for ( GaussianIterator itbeta = _shell_beta->firstGaussian(); itbeta != _shell_beta->lastGaussian(); ++itbeta){
                    const double& _decay_beta = (*itbeta)->decay;
                    
                    for ( GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma){

                        const double& _decay_gamma = (*itgamma)->decay;
            
            
            typedef boost::multi_array<double, 4> ma_type;
            typedef boost::multi_array_types::extent_range range;
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type R;
            R.resize(extents[ range(0, _nalpha+_nbeta ) ][ range(0, _nbeta ) ][ range(0, _ngamma)][range(0,_mmax)]);
            //initialize to zero
            for (index i = 0; i != _nalpha+_nbeta; ++i) {
                for (index j = 0; j != _nbeta; ++j) {
                    for (index k = 0; k != _ngamma; ++k) {
                        for (index l = 0; l != _mmax; ++l) {
                                       R[i][j][k][l] = 0.0;
                                   }
                               }
                           }
            
            
            double _decay=_decay_alpha + _decay_beta + _decay_gamma;
            double rzeta=0.5/(_decay_alpha+_decay_beta);
            double fak = 0.5 / (_decay);
            double gfak=_decay_gamma/_decay;
            double cfak= (_decay_alpha + _decay_beta)/_decay;
            double _T = (_decay_alpha+_decay_beta)*_decay_gamma/_decay*(_P-_pos_gamma)*(_P-_pos_gamma);
            vec _P=(_decay_alpha*_pos_alpha+_decay_beta*_pos_beta)/(_decay_alpha+_decay_beta);
            vec _W=(_decay_alpha*_pos_alpha+_decay_beta*_pos_beta+_decay*gamma*_pos_gamma)/_decay;
            double _dist1=(_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma);
            double _dist2=(_pos_gamma - _pos_beta) * (_pos_gamma - _pos_beta);
            double _dist3=(_pos_alpha - _pos_beta) * (_pos_alpha - _pos_beta);
            
            
            
            vec pma = _P - _pos_alpha;
            //vec pmb = _P - _pos_beta;
            vec wmp = _W - _P; 
            vec wmc = _W - _pos_gamma;
            
            
            double pma0 = 0.0;
            //double pmb0 = 0.0;
            double wmp0 = 0.0;
            double wmc0 = 0.0;

            double pma1 = 0.0;
            //double pmb1 = 0.0;
            double wmp1 = 0.0;
            double wmc1 = 0.0;

            double pma2 = 0.0;
            //double pmb2 = 0.0;
            double wmp2 = 0.0;
            double wmc2 = 0.0;
            
            if ((_dist1 + _dist2 + _dist3)<0.01){
          
            pma0 = pma.getX();
            //pmb0 = pmb.getX();
            wmp0 = wmp.getX();
            wmc0 = wmc.getX();
            pma1 = pma.getY();
            //pmb1 = pmb.getY();
            wmp1 = wmp.getY();
            wmc1 = wmc.getY():
            pma2 = pma.getZ();
            //pmb2 = pmb.getZ();
            wmp2 = wmp.getZ();
            wmc2 = wmc.getZ();
            }
            
            
            vector<double> _FmT(_mmax, 0.0); // that size needs to be checked!
            // call xint01(FmT,8,T,u_lower)
            XIntegrate(_FmT, _T);
            
            double sss = 8*pow(2*pi,0.25)*pow(decay_alpha*decay_beta*decay_gamma,0.75)/((decay_alpha+decay_beta)*decay_gamma);
            
            for (int _i=0;_i<_mmax;_i++){
                R[Cart::s][Cart::ss][Cart::s][_i]=sss*_FmT[_i];
            }
            
            
//omitting s-s-s
//Integral s - s - p - m5
if (_mmax >1 ){
if (_lmax_gamma>0){

R[Cart::s][Cart::s][Cart::y][5]+=wmc1*R[Cart::s][Cart::s][Cart::s][6]
R[Cart::s][Cart::s][Cart::x][5]+=wmc0*R[Cart::s][Cart::s][Cart::s][6]
R[Cart::s][Cart::s][Cart::z][5]+=wmc2*R[Cart::s][Cart::s][Cart::s][6]
}}
//------------------------------------------------------

//Integral s - s - d - m4
if (_mmax >2 ){
if (_lmax_gamma>1){

R[Cart::s][Cart::s][Cart::yy][4]+=wmc1*R[Cart::s][Cart::s][Cart::y][5]
R[Cart::s][Cart::s][Cart::xy][4]+=wmc0*R[Cart::s][Cart::s][Cart::y][5]
R[Cart::s][Cart::s][Cart::yz][4]+=wmc1*R[Cart::s][Cart::s][Cart::z][5]
R[Cart::s][Cart::s][Cart::xx][4]+=wmc0*R[Cart::s][Cart::s][Cart::x][5]
R[Cart::s][Cart::s][Cart::xz][4]+=wmc0*R[Cart::s][Cart::s][Cart::z][5]
R[Cart::s][Cart::s][Cart::zz][4]+=wmc2*R[Cart::s][Cart::s][Cart::z][5]
}}
//------------------------------------------------------

//Integral s - s - f - m3
if (_mmax >3 ){
if (_lmax_gamma>2){

R[Cart::s][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::s][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::s][Cart::s][Cart::y][3]-cfak*R[Cart::s][Cart::s][Cart::y][4])
R[Cart::s][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::s][Cart::s][Cart::yy][4]
R[Cart::s][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::s][Cart::s][Cart::yy][4]
R[Cart::s][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::s][Cart::s][Cart::xx][4]
R[Cart::s][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::s][Cart::s][Cart::yz][4]
R[Cart::s][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::s][Cart::s][Cart::zz][4]
R[Cart::s][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::s][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::s][Cart::s][Cart::x][3]-cfak*R[Cart::s][Cart::s][Cart::x][4])
R[Cart::s][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::s][Cart::s][Cart::xx][4]
R[Cart::s][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::s][Cart::s][Cart::zz][4]
R[Cart::s][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::s][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::s][Cart::s][Cart::z][3]-cfak*R[Cart::s][Cart::s][Cart::z][4])
}}
//------------------------------------------------------

//Integral p - s - s - m5
if (_mmax >1 ){
if (_lmax_alpha>0){

R[Cart::y][Cart::s][Cart::s][5]+=pma1*R[Cart::s][Cart::s][Cart::s][5]+wmp1*R[Cart::s][Cart::s][Cart::s][6]
R[Cart::x][Cart::s][Cart::s][5]+=pma0*R[Cart::s][Cart::s][Cart::s][5]+wmp0*R[Cart::s][Cart::s][Cart::s][6]
R[Cart::z][Cart::s][Cart::s][5]+=pma2*R[Cart::s][Cart::s][Cart::s][5]+wmp2*R[Cart::s][Cart::s][Cart::s][6]
}}
//------------------------------------------------------

//Integral p - s - p - m5
if (_mmax >2 ){
if (_lmax_alpha>0 && _lmax_gamma>0){

R[Cart::y][Cart::s][Cart::y][5]+=pma1*R[Cart::s][Cart::s][Cart::y][5]+wmp1*R[Cart::s][Cart::s][Cart::y][6]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::s][6]
R[Cart::y][Cart::s][Cart::x][5]+=pma1*R[Cart::s][Cart::s][Cart::x][5]+wmp1*R[Cart::s][Cart::s][Cart::x][6]
R[Cart::y][Cart::s][Cart::z][5]+=pma1*R[Cart::s][Cart::s][Cart::z][5]+wmp1*R[Cart::s][Cart::s][Cart::z][6]
R[Cart::x][Cart::s][Cart::y][5]+=pma0*R[Cart::s][Cart::s][Cart::y][5]+wmp0*R[Cart::s][Cart::s][Cart::y][6]
R[Cart::x][Cart::s][Cart::x][5]+=pma0*R[Cart::s][Cart::s][Cart::x][5]+wmp0*R[Cart::s][Cart::s][Cart::x][6]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::s][6]
R[Cart::x][Cart::s][Cart::z][5]+=pma0*R[Cart::s][Cart::s][Cart::z][5]+wmp0*R[Cart::s][Cart::s][Cart::z][6]
R[Cart::z][Cart::s][Cart::y][5]+=pma2*R[Cart::s][Cart::s][Cart::y][5]+wmp2*R[Cart::s][Cart::s][Cart::y][6]
R[Cart::z][Cart::s][Cart::x][5]+=pma2*R[Cart::s][Cart::s][Cart::x][5]+wmp2*R[Cart::s][Cart::s][Cart::x][6]
R[Cart::z][Cart::s][Cart::z][5]+=pma2*R[Cart::s][Cart::s][Cart::z][5]+wmp2*R[Cart::s][Cart::s][Cart::z][6]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::s][6]
}}
//------------------------------------------------------

//Integral p - s - d - m4
if (_mmax >3 ){
if (_lmax_alpha>0 && _lmax_gamma>1){

R[Cart::y][Cart::s][Cart::yy][4]+=wmc1*R[Cart::y][Cart::s][Cart::y][5]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::y][5]
R[Cart::y][Cart::s][Cart::xy][4]+=wmc0*R[Cart::y][Cart::s][Cart::y][5]
R[Cart::y][Cart::s][Cart::yz][4]+=wmc1*R[Cart::y][Cart::s][Cart::z][5]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::z][5]
R[Cart::y][Cart::s][Cart::xx][4]+=wmc0*R[Cart::y][Cart::s][Cart::x][5]
R[Cart::y][Cart::s][Cart::xz][4]+=wmc0*R[Cart::y][Cart::s][Cart::z][5]
R[Cart::y][Cart::s][Cart::zz][4]+=wmc2*R[Cart::y][Cart::s][Cart::z][5]
R[Cart::x][Cart::s][Cart::yy][4]+=wmc1*R[Cart::x][Cart::s][Cart::y][5]
R[Cart::x][Cart::s][Cart::xy][4]+=wmc0*R[Cart::x][Cart::s][Cart::y][5]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::y][5]
R[Cart::x][Cart::s][Cart::yz][4]+=wmc1*R[Cart::x][Cart::s][Cart::z][5]
R[Cart::x][Cart::s][Cart::xx][4]+=wmc0*R[Cart::x][Cart::s][Cart::x][5]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::x][5]
R[Cart::x][Cart::s][Cart::xz][4]+=wmc0*R[Cart::x][Cart::s][Cart::z][5]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::z][5]
R[Cart::x][Cart::s][Cart::zz][4]+=wmc2*R[Cart::x][Cart::s][Cart::z][5]
R[Cart::z][Cart::s][Cart::yy][4]+=wmc1*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::z][Cart::s][Cart::xy][4]+=wmc0*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::z][Cart::s][Cart::yz][4]+=wmc1*R[Cart::z][Cart::s][Cart::z][5]
R[Cart::z][Cart::s][Cart::xx][4]+=wmc0*R[Cart::z][Cart::s][Cart::x][5]
R[Cart::z][Cart::s][Cart::xz][4]+=wmc0*R[Cart::z][Cart::s][Cart::z][5]
R[Cart::z][Cart::s][Cart::zz][4]+=wmc2*R[Cart::z][Cart::s][Cart::z][5]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::z][5]
}}
//------------------------------------------------------

//Integral p - s - f - m3
if (_mmax >4 ){
if (_lmax_alpha>0 && _lmax_gamma>2){

R[Cart::y][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::y][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::y][Cart::s][Cart::y][3]-cfak*R[Cart::y][Cart::s][Cart::y][4])+0.5/_decay*1*R[Cart::s][Cart::s][Cart::yy][4]
R[Cart::y][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::y][Cart::s][Cart::yy][4]
R[Cart::y][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::y][Cart::s][Cart::yy][4]
R[Cart::y][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::y][Cart::s][Cart::xx][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::xx][4]
R[Cart::y][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::y][Cart::s][Cart::yz][4]
R[Cart::y][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::y][Cart::s][Cart::zz][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::zz][4]
R[Cart::y][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::y][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::y][Cart::s][Cart::x][3]-cfak*R[Cart::y][Cart::s][Cart::x][4])
R[Cart::y][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::y][Cart::s][Cart::xx][4]
R[Cart::y][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::y][Cart::s][Cart::zz][4]
R[Cart::y][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::y][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::y][Cart::s][Cart::z][3]-cfak*R[Cart::y][Cart::s][Cart::z][4])
R[Cart::x][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::x][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::x][Cart::s][Cart::y][3]-cfak*R[Cart::x][Cart::s][Cart::y][4])
R[Cart::x][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::x][Cart::s][Cart::yy][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::yy][4]
R[Cart::x][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::x][Cart::s][Cart::yy][4]
R[Cart::x][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::x][Cart::s][Cart::xx][4]
R[Cart::x][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::x][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::yz][4]
R[Cart::x][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::x][Cart::s][Cart::zz][4]
R[Cart::x][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::x][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::x][Cart::s][Cart::x][3]-cfak*R[Cart::x][Cart::s][Cart::x][4])+0.5/_decay*1*R[Cart::s][Cart::s][Cart::xx][4]
R[Cart::x][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::x][Cart::s][Cart::xx][4]
R[Cart::x][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::x][Cart::s][Cart::zz][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::zz][4]
R[Cart::x][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::x][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::x][Cart::s][Cart::z][3]-cfak*R[Cart::x][Cart::s][Cart::z][4])
R[Cart::z][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::z][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::z][Cart::s][Cart::y][3]-cfak*R[Cart::z][Cart::s][Cart::y][4])
R[Cart::z][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::z][Cart::s][Cart::yy][4]
R[Cart::z][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::z][Cart::s][Cart::yy][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::yy][4]
R[Cart::z][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::z][Cart::s][Cart::xx][4]
R[Cart::z][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::z][Cart::s][Cart::yz][4]
R[Cart::z][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::z][Cart::s][Cart::zz][4]
R[Cart::z][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::z][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::z][Cart::s][Cart::x][3]-cfak*R[Cart::z][Cart::s][Cart::x][4])
R[Cart::z][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::z][Cart::s][Cart::xx][4]+0.5/_decay*1*R[Cart::s][Cart::s][Cart::xx][4]
R[Cart::z][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::z][Cart::s][Cart::zz][4]
R[Cart::z][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::z][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::z][Cart::s][Cart::z][3]-cfak*R[Cart::z][Cart::s][Cart::z][4])+0.5/_decay*1*R[Cart::s][Cart::s][Cart::zz][4]
}}
//------------------------------------------------------

//Integral d - s - s - m4
if (_mmax >2 ){
if (_lmax_alpha>1){

R[Cart::yy][Cart::s][Cart::s][4]+=pma1*R[Cart::y][Cart::s][Cart::s][4]+wmp1*R[Cart::y][Cart::s][Cart::s][5]
R[Cart::xy][Cart::s][Cart::s][4]+=pma0*R[Cart::y][Cart::s][Cart::s][4]+wmp0*R[Cart::y][Cart::s][Cart::s][5]
R[Cart::yz][Cart::s][Cart::s][4]+=pma1*R[Cart::z][Cart::s][Cart::s][4]+wmp1*R[Cart::z][Cart::s][Cart::s][5]
R[Cart::xx][Cart::s][Cart::s][4]+=pma0*R[Cart::x][Cart::s][Cart::s][4]+wmp0*R[Cart::x][Cart::s][Cart::s][5]
R[Cart::xz][Cart::s][Cart::s][4]+=pma0*R[Cart::z][Cart::s][Cart::s][4]+wmp0*R[Cart::z][Cart::s][Cart::s][5]
R[Cart::zz][Cart::s][Cart::s][4]+=pma2*R[Cart::z][Cart::s][Cart::s][4]+wmp2*R[Cart::z][Cart::s][Cart::s][5]
}}
//------------------------------------------------------

//Integral d - s - p - m4
if (_mmax >3 ){
if (_lmax_alpha>1 && _lmax_gamma>0){

R[Cart::yy][Cart::s][Cart::y][4]+=pma1*R[Cart::y][Cart::s][Cart::y][4]+wmp1*R[Cart::y][Cart::s][Cart::y][5]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::s][5]
R[Cart::yy][Cart::s][Cart::x][4]+=pma1*R[Cart::y][Cart::s][Cart::x][4]+wmp1*R[Cart::y][Cart::s][Cart::x][5]
R[Cart::yy][Cart::s][Cart::z][4]+=pma1*R[Cart::y][Cart::s][Cart::z][4]+wmp1*R[Cart::y][Cart::s][Cart::z][5]
R[Cart::xy][Cart::s][Cart::y][4]+=pma0*R[Cart::y][Cart::s][Cart::y][4]+wmp0*R[Cart::y][Cart::s][Cart::y][5]
R[Cart::xy][Cart::s][Cart::x][4]+=pma0*R[Cart::y][Cart::s][Cart::x][4]+wmp0*R[Cart::y][Cart::s][Cart::x][5]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::s][5]
R[Cart::xy][Cart::s][Cart::z][4]+=pma0*R[Cart::y][Cart::s][Cart::z][4]+wmp0*R[Cart::y][Cart::s][Cart::z][5]
R[Cart::yz][Cart::s][Cart::y][4]+=pma1*R[Cart::z][Cart::s][Cart::y][4]+wmp1*R[Cart::z][Cart::s][Cart::y][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::s][5]
R[Cart::yz][Cart::s][Cart::x][4]+=pma1*R[Cart::z][Cart::s][Cart::x][4]+wmp1*R[Cart::z][Cart::s][Cart::x][5]
R[Cart::yz][Cart::s][Cart::z][4]+=pma1*R[Cart::z][Cart::s][Cart::z][4]+wmp1*R[Cart::z][Cart::s][Cart::z][5]
R[Cart::xx][Cart::s][Cart::y][4]+=pma0*R[Cart::x][Cart::s][Cart::y][4]+wmp0*R[Cart::x][Cart::s][Cart::y][5]
R[Cart::xx][Cart::s][Cart::x][4]+=pma0*R[Cart::x][Cart::s][Cart::x][4]+wmp0*R[Cart::x][Cart::s][Cart::x][5]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::s][5]
R[Cart::xx][Cart::s][Cart::z][4]+=pma0*R[Cart::x][Cart::s][Cart::z][4]+wmp0*R[Cart::x][Cart::s][Cart::z][5]
R[Cart::xz][Cart::s][Cart::y][4]+=pma0*R[Cart::z][Cart::s][Cart::y][4]+wmp0*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::xz][Cart::s][Cart::x][4]+=pma0*R[Cart::z][Cart::s][Cart::x][4]+wmp0*R[Cart::z][Cart::s][Cart::x][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::s][5]
R[Cart::xz][Cart::s][Cart::z][4]+=pma0*R[Cart::z][Cart::s][Cart::z][4]+wmp0*R[Cart::z][Cart::s][Cart::z][5]
R[Cart::zz][Cart::s][Cart::y][4]+=pma2*R[Cart::z][Cart::s][Cart::y][4]+wmp2*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::zz][Cart::s][Cart::x][4]+=pma2*R[Cart::z][Cart::s][Cart::x][4]+wmp2*R[Cart::z][Cart::s][Cart::x][5]
R[Cart::zz][Cart::s][Cart::z][4]+=pma2*R[Cart::z][Cart::s][Cart::z][4]+wmp2*R[Cart::z][Cart::s][Cart::z][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::s][5]
}}
//------------------------------------------------------

//Integral d - s - d - m4
if (_mmax >4 ){
if (_lmax_alpha>1 && _lmax_gamma>1){

R[Cart::yy][Cart::s][Cart::yy][4]+=pma1*R[Cart::y][Cart::s][Cart::yy][4]+wmp1*R[Cart::y][Cart::s][Cart::yy][5]+0.5/_decay*2*R[Cart::y][Cart::s][Cart::y][5]
R[Cart::yy][Cart::s][Cart::xy][4]+=pma1*R[Cart::y][Cart::s][Cart::xy][4]+wmp1*R[Cart::y][Cart::s][Cart::xy][5]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::x][5]
R[Cart::yy][Cart::s][Cart::yz][4]+=pma1*R[Cart::y][Cart::s][Cart::yz][4]+wmp1*R[Cart::y][Cart::s][Cart::yz][5]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::z][5]
R[Cart::yy][Cart::s][Cart::xx][4]+=pma1*R[Cart::y][Cart::s][Cart::xx][4]+wmp1*R[Cart::y][Cart::s][Cart::xx][5]
R[Cart::yy][Cart::s][Cart::xz][4]+=pma1*R[Cart::y][Cart::s][Cart::xz][4]+wmp1*R[Cart::y][Cart::s][Cart::xz][5]
R[Cart::yy][Cart::s][Cart::zz][4]+=pma1*R[Cart::y][Cart::s][Cart::zz][4]+wmp1*R[Cart::y][Cart::s][Cart::zz][5]
R[Cart::xy][Cart::s][Cart::yy][4]+=pma0*R[Cart::y][Cart::s][Cart::yy][4]+wmp0*R[Cart::y][Cart::s][Cart::yy][5]
R[Cart::xy][Cart::s][Cart::xy][4]+=pma0*R[Cart::y][Cart::s][Cart::xy][4]+wmp0*R[Cart::y][Cart::s][Cart::xy][5]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::y][5]
R[Cart::xy][Cart::s][Cart::yz][4]+=pma0*R[Cart::y][Cart::s][Cart::yz][4]+wmp0*R[Cart::y][Cart::s][Cart::yz][5]
R[Cart::xy][Cart::s][Cart::xx][4]+=pma0*R[Cart::y][Cart::s][Cart::xx][4]+wmp0*R[Cart::y][Cart::s][Cart::xx][5]+0.5/_decay*2*R[Cart::y][Cart::s][Cart::x][5]
R[Cart::xy][Cart::s][Cart::xz][4]+=pma0*R[Cart::y][Cart::s][Cart::xz][4]+wmp0*R[Cart::y][Cart::s][Cart::xz][5]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::z][5]
R[Cart::xy][Cart::s][Cart::zz][4]+=pma0*R[Cart::y][Cart::s][Cart::zz][4]+wmp0*R[Cart::y][Cart::s][Cart::zz][5]
R[Cart::yz][Cart::s][Cart::yy][4]+=pma1*R[Cart::z][Cart::s][Cart::yy][4]+wmp1*R[Cart::z][Cart::s][Cart::yy][5]+0.5/_decay*2*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::yz][Cart::s][Cart::xy][4]+=pma1*R[Cart::z][Cart::s][Cart::xy][4]+wmp1*R[Cart::z][Cart::s][Cart::xy][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::x][5]
R[Cart::yz][Cart::s][Cart::yz][4]+=pma1*R[Cart::z][Cart::s][Cart::yz][4]+wmp1*R[Cart::z][Cart::s][Cart::yz][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::z][5]
R[Cart::yz][Cart::s][Cart::xx][4]+=pma1*R[Cart::z][Cart::s][Cart::xx][4]+wmp1*R[Cart::z][Cart::s][Cart::xx][5]
R[Cart::yz][Cart::s][Cart::xz][4]+=pma1*R[Cart::z][Cart::s][Cart::xz][4]+wmp1*R[Cart::z][Cart::s][Cart::xz][5]
R[Cart::yz][Cart::s][Cart::zz][4]+=pma1*R[Cart::z][Cart::s][Cart::zz][4]+wmp1*R[Cart::z][Cart::s][Cart::zz][5]
R[Cart::xx][Cart::s][Cart::yy][4]+=pma0*R[Cart::x][Cart::s][Cart::yy][4]+wmp0*R[Cart::x][Cart::s][Cart::yy][5]
R[Cart::xx][Cart::s][Cart::xy][4]+=pma0*R[Cart::x][Cart::s][Cart::xy][4]+wmp0*R[Cart::x][Cart::s][Cart::xy][5]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::y][5]
R[Cart::xx][Cart::s][Cart::yz][4]+=pma0*R[Cart::x][Cart::s][Cart::yz][4]+wmp0*R[Cart::x][Cart::s][Cart::yz][5]
R[Cart::xx][Cart::s][Cart::xx][4]+=pma0*R[Cart::x][Cart::s][Cart::xx][4]+wmp0*R[Cart::x][Cart::s][Cart::xx][5]+0.5/_decay*2*R[Cart::x][Cart::s][Cart::x][5]
R[Cart::xx][Cart::s][Cart::xz][4]+=pma0*R[Cart::x][Cart::s][Cart::xz][4]+wmp0*R[Cart::x][Cart::s][Cart::xz][5]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::z][5]
R[Cart::xx][Cart::s][Cart::zz][4]+=pma0*R[Cart::x][Cart::s][Cart::zz][4]+wmp0*R[Cart::x][Cart::s][Cart::zz][5]
R[Cart::xz][Cart::s][Cart::yy][4]+=pma0*R[Cart::z][Cart::s][Cart::yy][4]+wmp0*R[Cart::z][Cart::s][Cart::yy][5]
R[Cart::xz][Cart::s][Cart::xy][4]+=pma0*R[Cart::z][Cart::s][Cart::xy][4]+wmp0*R[Cart::z][Cart::s][Cart::xy][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::xz][Cart::s][Cart::yz][4]+=pma0*R[Cart::z][Cart::s][Cart::yz][4]+wmp0*R[Cart::z][Cart::s][Cart::yz][5]
R[Cart::xz][Cart::s][Cart::xx][4]+=pma0*R[Cart::z][Cart::s][Cart::xx][4]+wmp0*R[Cart::z][Cart::s][Cart::xx][5]+0.5/_decay*2*R[Cart::z][Cart::s][Cart::x][5]
R[Cart::xz][Cart::s][Cart::xz][4]+=pma0*R[Cart::z][Cart::s][Cart::xz][4]+wmp0*R[Cart::z][Cart::s][Cart::xz][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::z][5]
R[Cart::xz][Cart::s][Cart::zz][4]+=pma0*R[Cart::z][Cart::s][Cart::zz][4]+wmp0*R[Cart::z][Cart::s][Cart::zz][5]
R[Cart::zz][Cart::s][Cart::yy][4]+=pma2*R[Cart::z][Cart::s][Cart::yy][4]+wmp2*R[Cart::z][Cart::s][Cart::yy][5]
R[Cart::zz][Cart::s][Cart::xy][4]+=pma2*R[Cart::z][Cart::s][Cart::xy][4]+wmp2*R[Cart::z][Cart::s][Cart::xy][5]
R[Cart::zz][Cart::s][Cart::yz][4]+=pma2*R[Cart::z][Cart::s][Cart::yz][4]+wmp2*R[Cart::z][Cart::s][Cart::yz][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::y][5]
R[Cart::zz][Cart::s][Cart::xx][4]+=pma2*R[Cart::z][Cart::s][Cart::xx][4]+wmp2*R[Cart::z][Cart::s][Cart::xx][5]
R[Cart::zz][Cart::s][Cart::xz][4]+=pma2*R[Cart::z][Cart::s][Cart::xz][4]+wmp2*R[Cart::z][Cart::s][Cart::xz][5]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::x][5]
R[Cart::zz][Cart::s][Cart::zz][4]+=pma2*R[Cart::z][Cart::s][Cart::zz][4]+wmp2*R[Cart::z][Cart::s][Cart::zz][5]+0.5/_decay*2*R[Cart::z][Cart::s][Cart::z][5]
}}
//------------------------------------------------------

//Integral d - s - f - m3
if (_mmax >5 ){
if (_lmax_alpha>1 && _lmax_gamma>2){

R[Cart::yy][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::yy][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::yy][Cart::s][Cart::y][3]-cfak*R[Cart::yy][Cart::s][Cart::y][4])+0.5/_decay*2*R[Cart::y][Cart::s][Cart::yy][4]
R[Cart::yy][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::yy][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::yy][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::yy][Cart::s][Cart::xx][4]+0.5/_decay*2*R[Cart::y][Cart::s][Cart::xx][4]
R[Cart::yy][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::yy][Cart::s][Cart::yz][4]
R[Cart::yy][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::yy][Cart::s][Cart::zz][4]+0.5/_decay*2*R[Cart::y][Cart::s][Cart::zz][4]
R[Cart::yy][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::yy][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::yy][Cart::s][Cart::x][3]-cfak*R[Cart::yy][Cart::s][Cart::x][4])
R[Cart::yy][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::yy][Cart::s][Cart::xx][4]
R[Cart::yy][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::yy][Cart::s][Cart::zz][4]
R[Cart::yy][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::yy][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::yy][Cart::s][Cart::z][3]-cfak*R[Cart::yy][Cart::s][Cart::z][4])
R[Cart::xy][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::xy][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::xy][Cart::s][Cart::y][3]-cfak*R[Cart::xy][Cart::s][Cart::y][4])+0.5/_decay*1*R[Cart::x][Cart::s][Cart::yy][4]
R[Cart::xy][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::xy][Cart::s][Cart::yy][4]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::yy][4]
R[Cart::xy][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::xy][Cart::s][Cart::yy][4]
R[Cart::xy][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::xy][Cart::s][Cart::xx][4]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::xx][4]
R[Cart::xy][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::xy][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::yz][4]
R[Cart::xy][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::xy][Cart::s][Cart::zz][4]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::zz][4]
R[Cart::xy][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::xy][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::xy][Cart::s][Cart::x][3]-cfak*R[Cart::xy][Cart::s][Cart::x][4])+0.5/_decay*1*R[Cart::y][Cart::s][Cart::xx][4]
R[Cart::xy][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::xy][Cart::s][Cart::xx][4]
R[Cart::xy][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::xy][Cart::s][Cart::zz][4]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::zz][4]
R[Cart::xy][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::xy][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::xy][Cart::s][Cart::z][3]-cfak*R[Cart::xy][Cart::s][Cart::z][4])
R[Cart::yz][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::yz][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::yz][Cart::s][Cart::y][3]-cfak*R[Cart::yz][Cart::s][Cart::y][4])+0.5/_decay*1*R[Cart::z][Cart::s][Cart::yy][4]
R[Cart::yz][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::yz][Cart::s][Cart::yy][4]
R[Cart::yz][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::yz][Cart::s][Cart::yy][4]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::yy][4]
R[Cart::yz][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::yz][Cart::s][Cart::xx][4]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::xx][4]
R[Cart::yz][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::yz][Cart::s][Cart::yz][4]
R[Cart::yz][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::yz][Cart::s][Cart::zz][4]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::zz][4]
R[Cart::yz][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::yz][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::yz][Cart::s][Cart::x][3]-cfak*R[Cart::yz][Cart::s][Cart::x][4])
R[Cart::yz][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::yz][Cart::s][Cart::xx][4]+0.5/_decay*1*R[Cart::y][Cart::s][Cart::xx][4]
R[Cart::yz][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::yz][Cart::s][Cart::zz][4]
R[Cart::yz][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::yz][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::yz][Cart::s][Cart::z][3]-cfak*R[Cart::yz][Cart::s][Cart::z][4])+0.5/_decay*1*R[Cart::y][Cart::s][Cart::zz][4]
R[Cart::xx][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::xx][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::xx][Cart::s][Cart::y][3]-cfak*R[Cart::xx][Cart::s][Cart::y][4])
R[Cart::xx][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::xx][Cart::s][Cart::yy][4]+0.5/_decay*2*R[Cart::x][Cart::s][Cart::yy][4]
R[Cart::xx][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::xx][Cart::s][Cart::yy][4]
R[Cart::xx][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xx][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::xx][Cart::s][Cart::yz][4]+0.5/_decay*2*R[Cart::x][Cart::s][Cart::yz][4]
R[Cart::xx][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::xx][Cart::s][Cart::zz][4]
R[Cart::xx][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::xx][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::xx][Cart::s][Cart::x][3]-cfak*R[Cart::xx][Cart::s][Cart::x][4])+0.5/_decay*2*R[Cart::x][Cart::s][Cart::xx][4]
R[Cart::xx][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xx][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::xx][Cart::s][Cart::zz][4]+0.5/_decay*2*R[Cart::x][Cart::s][Cart::zz][4]
R[Cart::xx][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::xx][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::xx][Cart::s][Cart::z][3]-cfak*R[Cart::xx][Cart::s][Cart::z][4])
R[Cart::xz][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::xz][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::xz][Cart::s][Cart::y][3]-cfak*R[Cart::xz][Cart::s][Cart::y][4])
R[Cart::xz][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::xz][Cart::s][Cart::yy][4]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::yy][4]
R[Cart::xz][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::xz][Cart::s][Cart::yy][4]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::yy][4]
R[Cart::xz][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::xz][Cart::s][Cart::xx][4]
R[Cart::xz][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::xz][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::yz][4]
R[Cart::xz][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::xz][Cart::s][Cart::zz][4]
R[Cart::xz][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::xz][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::xz][Cart::s][Cart::x][3]-cfak*R[Cart::xz][Cart::s][Cart::x][4])+0.5/_decay*1*R[Cart::z][Cart::s][Cart::xx][4]
R[Cart::xz][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::xz][Cart::s][Cart::xx][4]+0.5/_decay*1*R[Cart::x][Cart::s][Cart::xx][4]
R[Cart::xz][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::xz][Cart::s][Cart::zz][4]+0.5/_decay*1*R[Cart::z][Cart::s][Cart::zz][4]
R[Cart::xz][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::xz][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::xz][Cart::s][Cart::z][3]-cfak*R[Cart::xz][Cart::s][Cart::z][4])+0.5/_decay*1*R[Cart::x][Cart::s][Cart::zz][4]
R[Cart::zz][Cart::s][Cart::yyy][3]+=wmc1*R[Cart::zz][Cart::s][Cart::yy][4]+1/_decay_gamma*(R[Cart::zz][Cart::s][Cart::y][3]-cfak*R[Cart::zz][Cart::s][Cart::y][4])
R[Cart::zz][Cart::s][Cart::xyy][3]+=wmc0*R[Cart::zz][Cart::s][Cart::yy][4]
R[Cart::zz][Cart::s][Cart::yyz][3]+=wmc2*R[Cart::zz][Cart::s][Cart::yy][4]+0.5/_decay*2*R[Cart::z][Cart::s][Cart::yy][4]
R[Cart::zz][Cart::s][Cart::xxy][3]+=wmc1*R[Cart::zz][Cart::s][Cart::xx][4]
R[Cart::zz][Cart::s][Cart::xyz][3]+=wmc0*R[Cart::zz][Cart::s][Cart::yz][4]
R[Cart::zz][Cart::s][Cart::yzz][3]+=wmc1*R[Cart::zz][Cart::s][Cart::zz][4]
R[Cart::zz][Cart::s][Cart::xxx][3]+=wmc0*R[Cart::zz][Cart::s][Cart::xx][4]+1/_decay_gamma*(R[Cart::zz][Cart::s][Cart::x][3]-cfak*R[Cart::zz][Cart::s][Cart::x][4])
R[Cart::zz][Cart::s][Cart::xxz][3]+=wmc2*R[Cart::zz][Cart::s][Cart::xx][4]+0.5/_decay*2*R[Cart::z][Cart::s][Cart::xx][4]
R[Cart::zz][Cart::s][Cart::xzz][3]+=wmc0*R[Cart::zz][Cart::s][Cart::zz][4]
R[Cart::zz][Cart::s][Cart::zzz][3]+=wmc2*R[Cart::zz][Cart::s][Cart::zz][4]+1/_decay_gamma*(R[Cart::zz][Cart::s][Cart::z][3]-cfak*R[Cart::zz][Cart::s][Cart::z][4])+0.5/_decay*2*R[Cart::z][Cart::s][Cart::zz][4]
}}
//------------------------------------------------------

//Integral f - s - s - m3
if (_mmax >3 ){
if (_lmax_alpha>2){

R[Cart::yyy][Cart::s][Cart::s][3]+=pma1*R[Cart::yy][Cart::s][Cart::s][3]+wmp1*R[Cart::yy][Cart::s][Cart::s][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::s][3]-gfak*R[Cart::y][Cart::s][Cart::s][4])
R[Cart::xyy][Cart::s][Cart::s][3]+=pma0*R[Cart::yy][Cart::s][Cart::s][3]+wmp0*R[Cart::yy][Cart::s][Cart::s][4]
R[Cart::yyz][Cart::s][Cart::s][3]+=pma2*R[Cart::yy][Cart::s][Cart::s][3]+wmp2*R[Cart::yy][Cart::s][Cart::s][4]
R[Cart::xxy][Cart::s][Cart::s][3]+=pma1*R[Cart::xx][Cart::s][Cart::s][3]+wmp1*R[Cart::xx][Cart::s][Cart::s][4]
R[Cart::xyz][Cart::s][Cart::s][3]+=pma0*R[Cart::yz][Cart::s][Cart::s][3]+wmp0*R[Cart::yz][Cart::s][Cart::s][4]
R[Cart::yzz][Cart::s][Cart::s][3]+=pma1*R[Cart::zz][Cart::s][Cart::s][3]+wmp1*R[Cart::zz][Cart::s][Cart::s][4]
R[Cart::xxx][Cart::s][Cart::s][3]+=pma0*R[Cart::xx][Cart::s][Cart::s][3]+wmp0*R[Cart::xx][Cart::s][Cart::s][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::s][3]-gfak*R[Cart::x][Cart::s][Cart::s][4])
R[Cart::xxz][Cart::s][Cart::s][3]+=pma2*R[Cart::xx][Cart::s][Cart::s][3]+wmp2*R[Cart::xx][Cart::s][Cart::s][4]
R[Cart::xzz][Cart::s][Cart::s][3]+=pma0*R[Cart::zz][Cart::s][Cart::s][3]+wmp0*R[Cart::zz][Cart::s][Cart::s][4]
R[Cart::zzz][Cart::s][Cart::s][3]+=pma2*R[Cart::zz][Cart::s][Cart::s][3]+wmp2*R[Cart::zz][Cart::s][Cart::s][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::s][3]-gfak*R[Cart::z][Cart::s][Cart::s][4])
}}
//------------------------------------------------------

//Integral f - s - p - m3
if (_mmax >4 ){
if (_lmax_alpha>2 && _lmax_gamma>0){

R[Cart::yyy][Cart::s][Cart::y][3]+=pma1*R[Cart::yy][Cart::s][Cart::y][3]+wmp1*R[Cart::yy][Cart::s][Cart::y][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::y][3]-gfak*R[Cart::y][Cart::s][Cart::y][4])+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::s][4]
R[Cart::yyy][Cart::s][Cart::x][3]+=pma1*R[Cart::yy][Cart::s][Cart::x][3]+wmp1*R[Cart::yy][Cart::s][Cart::x][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::x][3]-gfak*R[Cart::y][Cart::s][Cart::x][4])
R[Cart::yyy][Cart::s][Cart::z][3]+=pma1*R[Cart::yy][Cart::s][Cart::z][3]+wmp1*R[Cart::yy][Cart::s][Cart::z][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::z][3]-gfak*R[Cart::y][Cart::s][Cart::z][4])
R[Cart::xyy][Cart::s][Cart::y][3]+=pma0*R[Cart::yy][Cart::s][Cart::y][3]+wmp0*R[Cart::yy][Cart::s][Cart::y][4]
R[Cart::xyy][Cart::s][Cart::x][3]+=pma0*R[Cart::yy][Cart::s][Cart::x][3]+wmp0*R[Cart::yy][Cart::s][Cart::x][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::s][4]
R[Cart::xyy][Cart::s][Cart::z][3]+=pma0*R[Cart::yy][Cart::s][Cart::z][3]+wmp0*R[Cart::yy][Cart::s][Cart::z][4]
R[Cart::yyz][Cart::s][Cart::y][3]+=pma2*R[Cart::yy][Cart::s][Cart::y][3]+wmp2*R[Cart::yy][Cart::s][Cart::y][4]
R[Cart::yyz][Cart::s][Cart::x][3]+=pma2*R[Cart::yy][Cart::s][Cart::x][3]+wmp2*R[Cart::yy][Cart::s][Cart::x][4]
R[Cart::yyz][Cart::s][Cart::z][3]+=pma2*R[Cart::yy][Cart::s][Cart::z][3]+wmp2*R[Cart::yy][Cart::s][Cart::z][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::s][4]
R[Cart::xxy][Cart::s][Cart::y][3]+=pma1*R[Cart::xx][Cart::s][Cart::y][3]+wmp1*R[Cart::xx][Cart::s][Cart::y][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::s][4]
R[Cart::xxy][Cart::s][Cart::x][3]+=pma1*R[Cart::xx][Cart::s][Cart::x][3]+wmp1*R[Cart::xx][Cart::s][Cart::x][4]
R[Cart::xxy][Cart::s][Cart::z][3]+=pma1*R[Cart::xx][Cart::s][Cart::z][3]+wmp1*R[Cart::xx][Cart::s][Cart::z][4]
R[Cart::xyz][Cart::s][Cart::y][3]+=pma0*R[Cart::yz][Cart::s][Cart::y][3]+wmp0*R[Cart::yz][Cart::s][Cart::y][4]
R[Cart::xyz][Cart::s][Cart::x][3]+=pma0*R[Cart::yz][Cart::s][Cart::x][3]+wmp0*R[Cart::yz][Cart::s][Cart::x][4]+0.5/_decay*1*R[Cart::yz][Cart::s][Cart::s][4]
R[Cart::xyz][Cart::s][Cart::z][3]+=pma0*R[Cart::yz][Cart::s][Cart::z][3]+wmp0*R[Cart::yz][Cart::s][Cart::z][4]
R[Cart::yzz][Cart::s][Cart::y][3]+=pma1*R[Cart::zz][Cart::s][Cart::y][3]+wmp1*R[Cart::zz][Cart::s][Cart::y][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::s][4]
R[Cart::yzz][Cart::s][Cart::x][3]+=pma1*R[Cart::zz][Cart::s][Cart::x][3]+wmp1*R[Cart::zz][Cart::s][Cart::x][4]
R[Cart::yzz][Cart::s][Cart::z][3]+=pma1*R[Cart::zz][Cart::s][Cart::z][3]+wmp1*R[Cart::zz][Cart::s][Cart::z][4]
R[Cart::xxx][Cart::s][Cart::y][3]+=pma0*R[Cart::xx][Cart::s][Cart::y][3]+wmp0*R[Cart::xx][Cart::s][Cart::y][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::y][3]-gfak*R[Cart::x][Cart::s][Cart::y][4])
R[Cart::xxx][Cart::s][Cart::x][3]+=pma0*R[Cart::xx][Cart::s][Cart::x][3]+wmp0*R[Cart::xx][Cart::s][Cart::x][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::x][3]-gfak*R[Cart::x][Cart::s][Cart::x][4])+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::s][4]
R[Cart::xxx][Cart::s][Cart::z][3]+=pma0*R[Cart::xx][Cart::s][Cart::z][3]+wmp0*R[Cart::xx][Cart::s][Cart::z][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::z][3]-gfak*R[Cart::x][Cart::s][Cart::z][4])
R[Cart::xxz][Cart::s][Cart::y][3]+=pma2*R[Cart::xx][Cart::s][Cart::y][3]+wmp2*R[Cart::xx][Cart::s][Cart::y][4]
R[Cart::xxz][Cart::s][Cart::x][3]+=pma2*R[Cart::xx][Cart::s][Cart::x][3]+wmp2*R[Cart::xx][Cart::s][Cart::x][4]
R[Cart::xxz][Cart::s][Cart::z][3]+=pma2*R[Cart::xx][Cart::s][Cart::z][3]+wmp2*R[Cart::xx][Cart::s][Cart::z][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::s][4]
R[Cart::xzz][Cart::s][Cart::y][3]+=pma0*R[Cart::zz][Cart::s][Cart::y][3]+wmp0*R[Cart::zz][Cart::s][Cart::y][4]
R[Cart::xzz][Cart::s][Cart::x][3]+=pma0*R[Cart::zz][Cart::s][Cart::x][3]+wmp0*R[Cart::zz][Cart::s][Cart::x][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::s][4]
R[Cart::xzz][Cart::s][Cart::z][3]+=pma0*R[Cart::zz][Cart::s][Cart::z][3]+wmp0*R[Cart::zz][Cart::s][Cart::z][4]
R[Cart::zzz][Cart::s][Cart::y][3]+=pma2*R[Cart::zz][Cart::s][Cart::y][3]+wmp2*R[Cart::zz][Cart::s][Cart::y][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::y][3]-gfak*R[Cart::z][Cart::s][Cart::y][4])
R[Cart::zzz][Cart::s][Cart::x][3]+=pma2*R[Cart::zz][Cart::s][Cart::x][3]+wmp2*R[Cart::zz][Cart::s][Cart::x][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::x][3]-gfak*R[Cart::z][Cart::s][Cart::x][4])
R[Cart::zzz][Cart::s][Cart::z][3]+=pma2*R[Cart::zz][Cart::s][Cart::z][3]+wmp2*R[Cart::zz][Cart::s][Cart::z][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::z][3]-gfak*R[Cart::z][Cart::s][Cart::z][4])+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::s][4]
}}
//------------------------------------------------------

//Integral f - s - d - m3
if (_mmax >5 ){
if (_lmax_alpha>2 && _lmax_gamma>1){

R[Cart::yyy][Cart::s][Cart::yy][3]+=pma1*R[Cart::yy][Cart::s][Cart::yy][3]+wmp1*R[Cart::yy][Cart::s][Cart::yy][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::yy][3]-gfak*R[Cart::y][Cart::s][Cart::yy][4])+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::y][4]
R[Cart::yyy][Cart::s][Cart::xy][3]+=pma1*R[Cart::yy][Cart::s][Cart::xy][3]+wmp1*R[Cart::yy][Cart::s][Cart::xy][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xy][3]-gfak*R[Cart::y][Cart::s][Cart::xy][4])+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::x][4]
R[Cart::yyy][Cart::s][Cart::yz][3]+=pma1*R[Cart::yy][Cart::s][Cart::yz][3]+wmp1*R[Cart::yy][Cart::s][Cart::yz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::yz][3]-gfak*R[Cart::y][Cart::s][Cart::yz][4])+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::z][4]
R[Cart::yyy][Cart::s][Cart::xx][3]+=pma1*R[Cart::yy][Cart::s][Cart::xx][3]+wmp1*R[Cart::yy][Cart::s][Cart::xx][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xx][3]-gfak*R[Cart::y][Cart::s][Cart::xx][4])
R[Cart::yyy][Cart::s][Cart::xz][3]+=pma1*R[Cart::yy][Cart::s][Cart::xz][3]+wmp1*R[Cart::yy][Cart::s][Cart::xz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xz][3]-gfak*R[Cart::y][Cart::s][Cart::xz][4])
R[Cart::yyy][Cart::s][Cart::zz][3]+=pma1*R[Cart::yy][Cart::s][Cart::zz][3]+wmp1*R[Cart::yy][Cart::s][Cart::zz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::zz][3]-gfak*R[Cart::y][Cart::s][Cart::zz][4])
R[Cart::xyy][Cart::s][Cart::yy][3]+=pma0*R[Cart::yy][Cart::s][Cart::yy][3]+wmp0*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::xyy][Cart::s][Cart::xy][3]+=pma0*R[Cart::yy][Cart::s][Cart::xy][3]+wmp0*R[Cart::yy][Cart::s][Cart::xy][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::y][4]
R[Cart::xyy][Cart::s][Cart::yz][3]+=pma0*R[Cart::yy][Cart::s][Cart::yz][3]+wmp0*R[Cart::yy][Cart::s][Cart::yz][4]
R[Cart::xyy][Cart::s][Cart::xx][3]+=pma0*R[Cart::yy][Cart::s][Cart::xx][3]+wmp0*R[Cart::yy][Cart::s][Cart::xx][4]+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::x][4]
R[Cart::xyy][Cart::s][Cart::xz][3]+=pma0*R[Cart::yy][Cart::s][Cart::xz][3]+wmp0*R[Cart::yy][Cart::s][Cart::xz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::z][4]
R[Cart::xyy][Cart::s][Cart::zz][3]+=pma0*R[Cart::yy][Cart::s][Cart::zz][3]+wmp0*R[Cart::yy][Cart::s][Cart::zz][4]
R[Cart::yyz][Cart::s][Cart::yy][3]+=pma2*R[Cart::yy][Cart::s][Cart::yy][3]+wmp2*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::yyz][Cart::s][Cart::xy][3]+=pma2*R[Cart::yy][Cart::s][Cart::xy][3]+wmp2*R[Cart::yy][Cart::s][Cart::xy][4]
R[Cart::yyz][Cart::s][Cart::yz][3]+=pma2*R[Cart::yy][Cart::s][Cart::yz][3]+wmp2*R[Cart::yy][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::y][4]
R[Cart::yyz][Cart::s][Cart::xx][3]+=pma2*R[Cart::yy][Cart::s][Cart::xx][3]+wmp2*R[Cart::yy][Cart::s][Cart::xx][4]
R[Cart::yyz][Cart::s][Cart::xz][3]+=pma2*R[Cart::yy][Cart::s][Cart::xz][3]+wmp2*R[Cart::yy][Cart::s][Cart::xz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::x][4]
R[Cart::yyz][Cart::s][Cart::zz][3]+=pma2*R[Cart::yy][Cart::s][Cart::zz][3]+wmp2*R[Cart::yy][Cart::s][Cart::zz][4]+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::z][4]
R[Cart::xxy][Cart::s][Cart::yy][3]+=pma1*R[Cart::xx][Cart::s][Cart::yy][3]+wmp1*R[Cart::xx][Cart::s][Cart::yy][4]+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::y][4]
R[Cart::xxy][Cart::s][Cart::xy][3]+=pma1*R[Cart::xx][Cart::s][Cart::xy][3]+wmp1*R[Cart::xx][Cart::s][Cart::xy][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::x][4]
R[Cart::xxy][Cart::s][Cart::yz][3]+=pma1*R[Cart::xx][Cart::s][Cart::yz][3]+wmp1*R[Cart::xx][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::z][4]
R[Cart::xxy][Cart::s][Cart::xx][3]+=pma1*R[Cart::xx][Cart::s][Cart::xx][3]+wmp1*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xxy][Cart::s][Cart::xz][3]+=pma1*R[Cart::xx][Cart::s][Cart::xz][3]+wmp1*R[Cart::xx][Cart::s][Cart::xz][4]
R[Cart::xxy][Cart::s][Cart::zz][3]+=pma1*R[Cart::xx][Cart::s][Cart::zz][3]+wmp1*R[Cart::xx][Cart::s][Cart::zz][4]
R[Cart::xyz][Cart::s][Cart::yy][3]+=pma0*R[Cart::yz][Cart::s][Cart::yy][3]+wmp0*R[Cart::yz][Cart::s][Cart::yy][4]
R[Cart::xyz][Cart::s][Cart::xy][3]+=pma0*R[Cart::yz][Cart::s][Cart::xy][3]+wmp0*R[Cart::yz][Cart::s][Cart::xy][4]+0.5/_decay*1*R[Cart::yz][Cart::s][Cart::y][4]
R[Cart::xyz][Cart::s][Cart::yz][3]+=pma0*R[Cart::yz][Cart::s][Cart::yz][3]+wmp0*R[Cart::yz][Cart::s][Cart::yz][4]
R[Cart::xyz][Cart::s][Cart::xx][3]+=pma0*R[Cart::yz][Cart::s][Cart::xx][3]+wmp0*R[Cart::yz][Cart::s][Cart::xx][4]+0.5/_decay*2*R[Cart::yz][Cart::s][Cart::x][4]
R[Cart::xyz][Cart::s][Cart::xz][3]+=pma0*R[Cart::yz][Cart::s][Cart::xz][3]+wmp0*R[Cart::yz][Cart::s][Cart::xz][4]+0.5/_decay*1*R[Cart::yz][Cart::s][Cart::z][4]
R[Cart::xyz][Cart::s][Cart::zz][3]+=pma0*R[Cart::yz][Cart::s][Cart::zz][3]+wmp0*R[Cart::yz][Cart::s][Cart::zz][4]
R[Cart::yzz][Cart::s][Cart::yy][3]+=pma1*R[Cart::zz][Cart::s][Cart::yy][3]+wmp1*R[Cart::zz][Cart::s][Cart::yy][4]+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::y][4]
R[Cart::yzz][Cart::s][Cart::xy][3]+=pma1*R[Cart::zz][Cart::s][Cart::xy][3]+wmp1*R[Cart::zz][Cart::s][Cart::xy][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::x][4]
R[Cart::yzz][Cart::s][Cart::yz][3]+=pma1*R[Cart::zz][Cart::s][Cart::yz][3]+wmp1*R[Cart::zz][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::z][4]
R[Cart::yzz][Cart::s][Cart::xx][3]+=pma1*R[Cart::zz][Cart::s][Cart::xx][3]+wmp1*R[Cart::zz][Cart::s][Cart::xx][4]
R[Cart::yzz][Cart::s][Cart::xz][3]+=pma1*R[Cart::zz][Cart::s][Cart::xz][3]+wmp1*R[Cart::zz][Cart::s][Cart::xz][4]
R[Cart::yzz][Cart::s][Cart::zz][3]+=pma1*R[Cart::zz][Cart::s][Cart::zz][3]+wmp1*R[Cart::zz][Cart::s][Cart::zz][4]
R[Cart::xxx][Cart::s][Cart::yy][3]+=pma0*R[Cart::xx][Cart::s][Cart::yy][3]+wmp0*R[Cart::xx][Cart::s][Cart::yy][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::yy][3]-gfak*R[Cart::x][Cart::s][Cart::yy][4])
R[Cart::xxx][Cart::s][Cart::xy][3]+=pma0*R[Cart::xx][Cart::s][Cart::xy][3]+wmp0*R[Cart::xx][Cart::s][Cart::xy][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xy][3]-gfak*R[Cart::x][Cart::s][Cart::xy][4])+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::y][4]
R[Cart::xxx][Cart::s][Cart::yz][3]+=pma0*R[Cart::xx][Cart::s][Cart::yz][3]+wmp0*R[Cart::xx][Cart::s][Cart::yz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::yz][3]-gfak*R[Cart::x][Cart::s][Cart::yz][4])
R[Cart::xxx][Cart::s][Cart::xx][3]+=pma0*R[Cart::xx][Cart::s][Cart::xx][3]+wmp0*R[Cart::xx][Cart::s][Cart::xx][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xx][3]-gfak*R[Cart::x][Cart::s][Cart::xx][4])+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::x][4]
R[Cart::xxx][Cart::s][Cart::xz][3]+=pma0*R[Cart::xx][Cart::s][Cart::xz][3]+wmp0*R[Cart::xx][Cart::s][Cart::xz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xz][3]-gfak*R[Cart::x][Cart::s][Cart::xz][4])+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::z][4]
R[Cart::xxx][Cart::s][Cart::zz][3]+=pma0*R[Cart::xx][Cart::s][Cart::zz][3]+wmp0*R[Cart::xx][Cart::s][Cart::zz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::zz][3]-gfak*R[Cart::x][Cart::s][Cart::zz][4])
R[Cart::xxz][Cart::s][Cart::yy][3]+=pma2*R[Cart::xx][Cart::s][Cart::yy][3]+wmp2*R[Cart::xx][Cart::s][Cart::yy][4]
R[Cart::xxz][Cart::s][Cart::xy][3]+=pma2*R[Cart::xx][Cart::s][Cart::xy][3]+wmp2*R[Cart::xx][Cart::s][Cart::xy][4]
R[Cart::xxz][Cart::s][Cart::yz][3]+=pma2*R[Cart::xx][Cart::s][Cart::yz][3]+wmp2*R[Cart::xx][Cart::s][Cart::yz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::y][4]
R[Cart::xxz][Cart::s][Cart::xx][3]+=pma2*R[Cart::xx][Cart::s][Cart::xx][3]+wmp2*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xxz][Cart::s][Cart::xz][3]+=pma2*R[Cart::xx][Cart::s][Cart::xz][3]+wmp2*R[Cart::xx][Cart::s][Cart::xz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::x][4]
R[Cart::xxz][Cart::s][Cart::zz][3]+=pma2*R[Cart::xx][Cart::s][Cart::zz][3]+wmp2*R[Cart::xx][Cart::s][Cart::zz][4]+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::z][4]
R[Cart::xzz][Cart::s][Cart::yy][3]+=pma0*R[Cart::zz][Cart::s][Cart::yy][3]+wmp0*R[Cart::zz][Cart::s][Cart::yy][4]
R[Cart::xzz][Cart::s][Cart::xy][3]+=pma0*R[Cart::zz][Cart::s][Cart::xy][3]+wmp0*R[Cart::zz][Cart::s][Cart::xy][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::y][4]
R[Cart::xzz][Cart::s][Cart::yz][3]+=pma0*R[Cart::zz][Cart::s][Cart::yz][3]+wmp0*R[Cart::zz][Cart::s][Cart::yz][4]
R[Cart::xzz][Cart::s][Cart::xx][3]+=pma0*R[Cart::zz][Cart::s][Cart::xx][3]+wmp0*R[Cart::zz][Cart::s][Cart::xx][4]+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::x][4]
R[Cart::xzz][Cart::s][Cart::xz][3]+=pma0*R[Cart::zz][Cart::s][Cart::xz][3]+wmp0*R[Cart::zz][Cart::s][Cart::xz][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::z][4]
R[Cart::xzz][Cart::s][Cart::zz][3]+=pma0*R[Cart::zz][Cart::s][Cart::zz][3]+wmp0*R[Cart::zz][Cart::s][Cart::zz][4]
R[Cart::zzz][Cart::s][Cart::yy][3]+=pma2*R[Cart::zz][Cart::s][Cart::yy][3]+wmp2*R[Cart::zz][Cart::s][Cart::yy][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::yy][3]-gfak*R[Cart::z][Cart::s][Cart::yy][4])
R[Cart::zzz][Cart::s][Cart::xy][3]+=pma2*R[Cart::zz][Cart::s][Cart::xy][3]+wmp2*R[Cart::zz][Cart::s][Cart::xy][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xy][3]-gfak*R[Cart::z][Cart::s][Cart::xy][4])
R[Cart::zzz][Cart::s][Cart::yz][3]+=pma2*R[Cart::zz][Cart::s][Cart::yz][3]+wmp2*R[Cart::zz][Cart::s][Cart::yz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::yz][3]-gfak*R[Cart::z][Cart::s][Cart::yz][4])+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::y][4]
R[Cart::zzz][Cart::s][Cart::xx][3]+=pma2*R[Cart::zz][Cart::s][Cart::xx][3]+wmp2*R[Cart::zz][Cart::s][Cart::xx][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xx][3]-gfak*R[Cart::z][Cart::s][Cart::xx][4])
R[Cart::zzz][Cart::s][Cart::xz][3]+=pma2*R[Cart::zz][Cart::s][Cart::xz][3]+wmp2*R[Cart::zz][Cart::s][Cart::xz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xz][3]-gfak*R[Cart::z][Cart::s][Cart::xz][4])+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::x][4]
R[Cart::zzz][Cart::s][Cart::zz][3]+=pma2*R[Cart::zz][Cart::s][Cart::zz][3]+wmp2*R[Cart::zz][Cart::s][Cart::zz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::zz][3]-gfak*R[Cart::z][Cart::s][Cart::zz][4])+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::z][4]
}}
//------------------------------------------------------

//Integral f - s - f - m3
if (_mmax >6 ){
if (_lmax_alpha>2 && _lmax_gamma>2){

R[Cart::yyy][Cart::s][Cart::yyy][3]+=pma1*R[Cart::yy][Cart::s][Cart::yyy][3]+wmp1*R[Cart::yy][Cart::s][Cart::yyy][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::yyy][3]-gfak*R[Cart::y][Cart::s][Cart::yyy][4])+0.5/_decay*3*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::yyy][Cart::s][Cart::xyy][3]+=pma1*R[Cart::yy][Cart::s][Cart::xyy][3]+wmp1*R[Cart::yy][Cart::s][Cart::xyy][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xyy][3]-gfak*R[Cart::y][Cart::s][Cart::xyy][4])+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::xy][4]
R[Cart::yyy][Cart::s][Cart::yyz][3]+=pma1*R[Cart::yy][Cart::s][Cart::yyz][3]+wmp1*R[Cart::yy][Cart::s][Cart::yyz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::yyz][3]-gfak*R[Cart::y][Cart::s][Cart::yyz][4])+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::yz][4]
R[Cart::yyy][Cart::s][Cart::xxy][3]+=pma1*R[Cart::yy][Cart::s][Cart::xxy][3]+wmp1*R[Cart::yy][Cart::s][Cart::xxy][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xxy][3]-gfak*R[Cart::y][Cart::s][Cart::xxy][4])+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::xx][4]
R[Cart::yyy][Cart::s][Cart::xyz][3]+=pma1*R[Cart::yy][Cart::s][Cart::xyz][3]+wmp1*R[Cart::yy][Cart::s][Cart::xyz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xyz][3]-gfak*R[Cart::y][Cart::s][Cart::xyz][4])+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::xz][4]
R[Cart::yyy][Cart::s][Cart::yzz][3]+=pma1*R[Cart::yy][Cart::s][Cart::yzz][3]+wmp1*R[Cart::yy][Cart::s][Cart::yzz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::yzz][3]-gfak*R[Cart::y][Cart::s][Cart::yzz][4])+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::zz][4]
R[Cart::yyy][Cart::s][Cart::xxx][3]+=pma1*R[Cart::yy][Cart::s][Cart::xxx][3]+wmp1*R[Cart::yy][Cart::s][Cart::xxx][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xxx][3]-gfak*R[Cart::y][Cart::s][Cart::xxx][4])
R[Cart::yyy][Cart::s][Cart::xxz][3]+=pma1*R[Cart::yy][Cart::s][Cart::xxz][3]+wmp1*R[Cart::yy][Cart::s][Cart::xxz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xxz][3]-gfak*R[Cart::y][Cart::s][Cart::xxz][4])
R[Cart::yyy][Cart::s][Cart::xzz][3]+=pma1*R[Cart::yy][Cart::s][Cart::xzz][3]+wmp1*R[Cart::yy][Cart::s][Cart::xzz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::xzz][3]-gfak*R[Cart::y][Cart::s][Cart::xzz][4])
R[Cart::yyy][Cart::s][Cart::zzz][3]+=pma1*R[Cart::yy][Cart::s][Cart::zzz][3]+wmp1*R[Cart::yy][Cart::s][Cart::zzz][4]+1*rzeta*(R[Cart::y][Cart::s][Cart::zzz][3]-gfak*R[Cart::y][Cart::s][Cart::zzz][4])
R[Cart::xyy][Cart::s][Cart::yyy][3]+=pma0*R[Cart::yy][Cart::s][Cart::yyy][3]+wmp0*R[Cart::yy][Cart::s][Cart::yyy][4]
R[Cart::xyy][Cart::s][Cart::xyy][3]+=pma0*R[Cart::yy][Cart::s][Cart::xyy][3]+wmp0*R[Cart::yy][Cart::s][Cart::xyy][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::xyy][Cart::s][Cart::yyz][3]+=pma0*R[Cart::yy][Cart::s][Cart::yyz][3]+wmp0*R[Cart::yy][Cart::s][Cart::yyz][4]
R[Cart::xyy][Cart::s][Cart::xxy][3]+=pma0*R[Cart::yy][Cart::s][Cart::xxy][3]+wmp0*R[Cart::yy][Cart::s][Cart::xxy][4]+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::xy][4]
R[Cart::xyy][Cart::s][Cart::xyz][3]+=pma0*R[Cart::yy][Cart::s][Cart::xyz][3]+wmp0*R[Cart::yy][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::yz][4]
R[Cart::xyy][Cart::s][Cart::yzz][3]+=pma0*R[Cart::yy][Cart::s][Cart::yzz][3]+wmp0*R[Cart::yy][Cart::s][Cart::yzz][4]
R[Cart::xyy][Cart::s][Cart::xxx][3]+=pma0*R[Cart::yy][Cart::s][Cart::xxx][3]+wmp0*R[Cart::yy][Cart::s][Cart::xxx][4]+0.5/_decay*3*R[Cart::yy][Cart::s][Cart::xx][4]
R[Cart::xyy][Cart::s][Cart::xxz][3]+=pma0*R[Cart::yy][Cart::s][Cart::xxz][3]+wmp0*R[Cart::yy][Cart::s][Cart::xxz][4]+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::xz][4]
R[Cart::xyy][Cart::s][Cart::xzz][3]+=pma0*R[Cart::yy][Cart::s][Cart::xzz][3]+wmp0*R[Cart::yy][Cart::s][Cart::xzz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::zz][4]
R[Cart::xyy][Cart::s][Cart::zzz][3]+=pma0*R[Cart::yy][Cart::s][Cart::zzz][3]+wmp0*R[Cart::yy][Cart::s][Cart::zzz][4]
R[Cart::yyz][Cart::s][Cart::yyy][3]+=pma2*R[Cart::yy][Cart::s][Cart::yyy][3]+wmp2*R[Cart::yy][Cart::s][Cart::yyy][4]
R[Cart::yyz][Cart::s][Cart::xyy][3]+=pma2*R[Cart::yy][Cart::s][Cart::xyy][3]+wmp2*R[Cart::yy][Cart::s][Cart::xyy][4]
R[Cart::yyz][Cart::s][Cart::yyz][3]+=pma2*R[Cart::yy][Cart::s][Cart::yyz][3]+wmp2*R[Cart::yy][Cart::s][Cart::yyz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::yy][4]
R[Cart::yyz][Cart::s][Cart::xxy][3]+=pma2*R[Cart::yy][Cart::s][Cart::xxy][3]+wmp2*R[Cart::yy][Cart::s][Cart::xxy][4]
R[Cart::yyz][Cart::s][Cart::xyz][3]+=pma2*R[Cart::yy][Cart::s][Cart::xyz][3]+wmp2*R[Cart::yy][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::xy][4]
R[Cart::yyz][Cart::s][Cart::yzz][3]+=pma2*R[Cart::yy][Cart::s][Cart::yzz][3]+wmp2*R[Cart::yy][Cart::s][Cart::yzz][4]+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::yz][4]
R[Cart::yyz][Cart::s][Cart::xxx][3]+=pma2*R[Cart::yy][Cart::s][Cart::xxx][3]+wmp2*R[Cart::yy][Cart::s][Cart::xxx][4]
R[Cart::yyz][Cart::s][Cart::xxz][3]+=pma2*R[Cart::yy][Cart::s][Cart::xxz][3]+wmp2*R[Cart::yy][Cart::s][Cart::xxz][4]+0.5/_decay*1*R[Cart::yy][Cart::s][Cart::xx][4]
R[Cart::yyz][Cart::s][Cart::xzz][3]+=pma2*R[Cart::yy][Cart::s][Cart::xzz][3]+wmp2*R[Cart::yy][Cart::s][Cart::xzz][4]+0.5/_decay*2*R[Cart::yy][Cart::s][Cart::xz][4]
R[Cart::yyz][Cart::s][Cart::zzz][3]+=pma2*R[Cart::yy][Cart::s][Cart::zzz][3]+wmp2*R[Cart::yy][Cart::s][Cart::zzz][4]+0.5/_decay*3*R[Cart::yy][Cart::s][Cart::zz][4]
R[Cart::xxy][Cart::s][Cart::yyy][3]+=pma1*R[Cart::xx][Cart::s][Cart::yyy][3]+wmp1*R[Cart::xx][Cart::s][Cart::yyy][4]+0.5/_decay*3*R[Cart::xx][Cart::s][Cart::yy][4]
R[Cart::xxy][Cart::s][Cart::xyy][3]+=pma1*R[Cart::xx][Cart::s][Cart::xyy][3]+wmp1*R[Cart::xx][Cart::s][Cart::xyy][4]+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::xy][4]
R[Cart::xxy][Cart::s][Cart::yyz][3]+=pma1*R[Cart::xx][Cart::s][Cart::yyz][3]+wmp1*R[Cart::xx][Cart::s][Cart::yyz][4]+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::yz][4]
R[Cart::xxy][Cart::s][Cart::xxy][3]+=pma1*R[Cart::xx][Cart::s][Cart::xxy][3]+wmp1*R[Cart::xx][Cart::s][Cart::xxy][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xxy][Cart::s][Cart::xyz][3]+=pma1*R[Cart::xx][Cart::s][Cart::xyz][3]+wmp1*R[Cart::xx][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::xz][4]
R[Cart::xxy][Cart::s][Cart::yzz][3]+=pma1*R[Cart::xx][Cart::s][Cart::yzz][3]+wmp1*R[Cart::xx][Cart::s][Cart::yzz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::zz][4]
R[Cart::xxy][Cart::s][Cart::xxx][3]+=pma1*R[Cart::xx][Cart::s][Cart::xxx][3]+wmp1*R[Cart::xx][Cart::s][Cart::xxx][4]
R[Cart::xxy][Cart::s][Cart::xxz][3]+=pma1*R[Cart::xx][Cart::s][Cart::xxz][3]+wmp1*R[Cart::xx][Cart::s][Cart::xxz][4]
R[Cart::xxy][Cart::s][Cart::xzz][3]+=pma1*R[Cart::xx][Cart::s][Cart::xzz][3]+wmp1*R[Cart::xx][Cart::s][Cart::xzz][4]
R[Cart::xxy][Cart::s][Cart::zzz][3]+=pma1*R[Cart::xx][Cart::s][Cart::zzz][3]+wmp1*R[Cart::xx][Cart::s][Cart::zzz][4]
R[Cart::xyz][Cart::s][Cart::yyy][3]+=pma0*R[Cart::yz][Cart::s][Cart::yyy][3]+wmp0*R[Cart::yz][Cart::s][Cart::yyy][4]
R[Cart::xyz][Cart::s][Cart::xyy][3]+=pma0*R[Cart::yz][Cart::s][Cart::xyy][3]+wmp0*R[Cart::yz][Cart::s][Cart::xyy][4]+0.5/_decay*1*R[Cart::yz][Cart::s][Cart::yy][4]
R[Cart::xyz][Cart::s][Cart::yyz][3]+=pma0*R[Cart::yz][Cart::s][Cart::yyz][3]+wmp0*R[Cart::yz][Cart::s][Cart::yyz][4]
R[Cart::xyz][Cart::s][Cart::xxy][3]+=pma0*R[Cart::yz][Cart::s][Cart::xxy][3]+wmp0*R[Cart::yz][Cart::s][Cart::xxy][4]+0.5/_decay*2*R[Cart::yz][Cart::s][Cart::xy][4]
R[Cart::xyz][Cart::s][Cart::xyz][3]+=pma0*R[Cart::yz][Cart::s][Cart::xyz][3]+wmp0*R[Cart::yz][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::yz][Cart::s][Cart::yz][4]
R[Cart::xyz][Cart::s][Cart::yzz][3]+=pma0*R[Cart::yz][Cart::s][Cart::yzz][3]+wmp0*R[Cart::yz][Cart::s][Cart::yzz][4]
R[Cart::xyz][Cart::s][Cart::xxx][3]+=pma0*R[Cart::yz][Cart::s][Cart::xxx][3]+wmp0*R[Cart::yz][Cart::s][Cart::xxx][4]+0.5/_decay*3*R[Cart::yz][Cart::s][Cart::xx][4]
R[Cart::xyz][Cart::s][Cart::xxz][3]+=pma0*R[Cart::yz][Cart::s][Cart::xxz][3]+wmp0*R[Cart::yz][Cart::s][Cart::xxz][4]+0.5/_decay*2*R[Cart::yz][Cart::s][Cart::xz][4]
R[Cart::xyz][Cart::s][Cart::xzz][3]+=pma0*R[Cart::yz][Cart::s][Cart::xzz][3]+wmp0*R[Cart::yz][Cart::s][Cart::xzz][4]+0.5/_decay*1*R[Cart::yz][Cart::s][Cart::zz][4]
R[Cart::xyz][Cart::s][Cart::zzz][3]+=pma0*R[Cart::yz][Cart::s][Cart::zzz][3]+wmp0*R[Cart::yz][Cart::s][Cart::zzz][4]
R[Cart::yzz][Cart::s][Cart::yyy][3]+=pma1*R[Cart::zz][Cart::s][Cart::yyy][3]+wmp1*R[Cart::zz][Cart::s][Cart::yyy][4]+0.5/_decay*3*R[Cart::zz][Cart::s][Cart::yy][4]
R[Cart::yzz][Cart::s][Cart::xyy][3]+=pma1*R[Cart::zz][Cart::s][Cart::xyy][3]+wmp1*R[Cart::zz][Cart::s][Cart::xyy][4]+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::xy][4]
R[Cart::yzz][Cart::s][Cart::yyz][3]+=pma1*R[Cart::zz][Cart::s][Cart::yyz][3]+wmp1*R[Cart::zz][Cart::s][Cart::yyz][4]+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::yz][4]
R[Cart::yzz][Cart::s][Cart::xxy][3]+=pma1*R[Cart::zz][Cart::s][Cart::xxy][3]+wmp1*R[Cart::zz][Cart::s][Cart::xxy][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::xx][4]
R[Cart::yzz][Cart::s][Cart::xyz][3]+=pma1*R[Cart::zz][Cart::s][Cart::xyz][3]+wmp1*R[Cart::zz][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::xz][4]
R[Cart::yzz][Cart::s][Cart::yzz][3]+=pma1*R[Cart::zz][Cart::s][Cart::yzz][3]+wmp1*R[Cart::zz][Cart::s][Cart::yzz][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::zz][4]
R[Cart::yzz][Cart::s][Cart::xxx][3]+=pma1*R[Cart::zz][Cart::s][Cart::xxx][3]+wmp1*R[Cart::zz][Cart::s][Cart::xxx][4]
R[Cart::yzz][Cart::s][Cart::xxz][3]+=pma1*R[Cart::zz][Cart::s][Cart::xxz][3]+wmp1*R[Cart::zz][Cart::s][Cart::xxz][4]
R[Cart::yzz][Cart::s][Cart::xzz][3]+=pma1*R[Cart::zz][Cart::s][Cart::xzz][3]+wmp1*R[Cart::zz][Cart::s][Cart::xzz][4]
R[Cart::yzz][Cart::s][Cart::zzz][3]+=pma1*R[Cart::zz][Cart::s][Cart::zzz][3]+wmp1*R[Cart::zz][Cart::s][Cart::zzz][4]
R[Cart::xxx][Cart::s][Cart::yyy][3]+=pma0*R[Cart::xx][Cart::s][Cart::yyy][3]+wmp0*R[Cart::xx][Cart::s][Cart::yyy][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::yyy][3]-gfak*R[Cart::x][Cart::s][Cart::yyy][4])
R[Cart::xxx][Cart::s][Cart::xyy][3]+=pma0*R[Cart::xx][Cart::s][Cart::xyy][3]+wmp0*R[Cart::xx][Cart::s][Cart::xyy][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xyy][3]-gfak*R[Cart::x][Cart::s][Cart::xyy][4])+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::yy][4]
R[Cart::xxx][Cart::s][Cart::yyz][3]+=pma0*R[Cart::xx][Cart::s][Cart::yyz][3]+wmp0*R[Cart::xx][Cart::s][Cart::yyz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::yyz][3]-gfak*R[Cart::x][Cart::s][Cart::yyz][4])
R[Cart::xxx][Cart::s][Cart::xxy][3]+=pma0*R[Cart::xx][Cart::s][Cart::xxy][3]+wmp0*R[Cart::xx][Cart::s][Cart::xxy][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xxy][3]-gfak*R[Cart::x][Cart::s][Cart::xxy][4])+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::xy][4]
R[Cart::xxx][Cart::s][Cart::xyz][3]+=pma0*R[Cart::xx][Cart::s][Cart::xyz][3]+wmp0*R[Cart::xx][Cart::s][Cart::xyz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xyz][3]-gfak*R[Cart::x][Cart::s][Cart::xyz][4])+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::yz][4]
R[Cart::xxx][Cart::s][Cart::yzz][3]+=pma0*R[Cart::xx][Cart::s][Cart::yzz][3]+wmp0*R[Cart::xx][Cart::s][Cart::yzz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::yzz][3]-gfak*R[Cart::x][Cart::s][Cart::yzz][4])
R[Cart::xxx][Cart::s][Cart::xxx][3]+=pma0*R[Cart::xx][Cart::s][Cart::xxx][3]+wmp0*R[Cart::xx][Cart::s][Cart::xxx][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xxx][3]-gfak*R[Cart::x][Cart::s][Cart::xxx][4])+0.5/_decay*3*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xxx][Cart::s][Cart::xxz][3]+=pma0*R[Cart::xx][Cart::s][Cart::xxz][3]+wmp0*R[Cart::xx][Cart::s][Cart::xxz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xxz][3]-gfak*R[Cart::x][Cart::s][Cart::xxz][4])+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::xz][4]
R[Cart::xxx][Cart::s][Cart::xzz][3]+=pma0*R[Cart::xx][Cart::s][Cart::xzz][3]+wmp0*R[Cart::xx][Cart::s][Cart::xzz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::xzz][3]-gfak*R[Cart::x][Cart::s][Cart::xzz][4])+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::zz][4]
R[Cart::xxx][Cart::s][Cart::zzz][3]+=pma0*R[Cart::xx][Cart::s][Cart::zzz][3]+wmp0*R[Cart::xx][Cart::s][Cart::zzz][4]+1*rzeta*(R[Cart::x][Cart::s][Cart::zzz][3]-gfak*R[Cart::x][Cart::s][Cart::zzz][4])
R[Cart::xxz][Cart::s][Cart::yyy][3]+=pma2*R[Cart::xx][Cart::s][Cart::yyy][3]+wmp2*R[Cart::xx][Cart::s][Cart::yyy][4]
R[Cart::xxz][Cart::s][Cart::xyy][3]+=pma2*R[Cart::xx][Cart::s][Cart::xyy][3]+wmp2*R[Cart::xx][Cart::s][Cart::xyy][4]
R[Cart::xxz][Cart::s][Cart::yyz][3]+=pma2*R[Cart::xx][Cart::s][Cart::yyz][3]+wmp2*R[Cart::xx][Cart::s][Cart::yyz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::yy][4]
R[Cart::xxz][Cart::s][Cart::xxy][3]+=pma2*R[Cart::xx][Cart::s][Cart::xxy][3]+wmp2*R[Cart::xx][Cart::s][Cart::xxy][4]
R[Cart::xxz][Cart::s][Cart::xyz][3]+=pma2*R[Cart::xx][Cart::s][Cart::xyz][3]+wmp2*R[Cart::xx][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::xy][4]
R[Cart::xxz][Cart::s][Cart::yzz][3]+=pma2*R[Cart::xx][Cart::s][Cart::yzz][3]+wmp2*R[Cart::xx][Cart::s][Cart::yzz][4]+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::yz][4]
R[Cart::xxz][Cart::s][Cart::xxx][3]+=pma2*R[Cart::xx][Cart::s][Cart::xxx][3]+wmp2*R[Cart::xx][Cart::s][Cart::xxx][4]
R[Cart::xxz][Cart::s][Cart::xxz][3]+=pma2*R[Cart::xx][Cart::s][Cart::xxz][3]+wmp2*R[Cart::xx][Cart::s][Cart::xxz][4]+0.5/_decay*1*R[Cart::xx][Cart::s][Cart::xx][4]
R[Cart::xxz][Cart::s][Cart::xzz][3]+=pma2*R[Cart::xx][Cart::s][Cart::xzz][3]+wmp2*R[Cart::xx][Cart::s][Cart::xzz][4]+0.5/_decay*2*R[Cart::xx][Cart::s][Cart::xz][4]
R[Cart::xxz][Cart::s][Cart::zzz][3]+=pma2*R[Cart::xx][Cart::s][Cart::zzz][3]+wmp2*R[Cart::xx][Cart::s][Cart::zzz][4]+0.5/_decay*3*R[Cart::xx][Cart::s][Cart::zz][4]
R[Cart::xzz][Cart::s][Cart::yyy][3]+=pma0*R[Cart::zz][Cart::s][Cart::yyy][3]+wmp0*R[Cart::zz][Cart::s][Cart::yyy][4]
R[Cart::xzz][Cart::s][Cart::xyy][3]+=pma0*R[Cart::zz][Cart::s][Cart::xyy][3]+wmp0*R[Cart::zz][Cart::s][Cart::xyy][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::yy][4]
R[Cart::xzz][Cart::s][Cart::yyz][3]+=pma0*R[Cart::zz][Cart::s][Cart::yyz][3]+wmp0*R[Cart::zz][Cart::s][Cart::yyz][4]
R[Cart::xzz][Cart::s][Cart::xxy][3]+=pma0*R[Cart::zz][Cart::s][Cart::xxy][3]+wmp0*R[Cart::zz][Cart::s][Cart::xxy][4]+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::xy][4]
R[Cart::xzz][Cart::s][Cart::xyz][3]+=pma0*R[Cart::zz][Cart::s][Cart::xyz][3]+wmp0*R[Cart::zz][Cart::s][Cart::xyz][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::yz][4]
R[Cart::xzz][Cart::s][Cart::yzz][3]+=pma0*R[Cart::zz][Cart::s][Cart::yzz][3]+wmp0*R[Cart::zz][Cart::s][Cart::yzz][4]
R[Cart::xzz][Cart::s][Cart::xxx][3]+=pma0*R[Cart::zz][Cart::s][Cart::xxx][3]+wmp0*R[Cart::zz][Cart::s][Cart::xxx][4]+0.5/_decay*3*R[Cart::zz][Cart::s][Cart::xx][4]
R[Cart::xzz][Cart::s][Cart::xxz][3]+=pma0*R[Cart::zz][Cart::s][Cart::xxz][3]+wmp0*R[Cart::zz][Cart::s][Cart::xxz][4]+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::xz][4]
R[Cart::xzz][Cart::s][Cart::xzz][3]+=pma0*R[Cart::zz][Cart::s][Cart::xzz][3]+wmp0*R[Cart::zz][Cart::s][Cart::xzz][4]+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::zz][4]
R[Cart::xzz][Cart::s][Cart::zzz][3]+=pma0*R[Cart::zz][Cart::s][Cart::zzz][3]+wmp0*R[Cart::zz][Cart::s][Cart::zzz][4]
R[Cart::zzz][Cart::s][Cart::yyy][3]+=pma2*R[Cart::zz][Cart::s][Cart::yyy][3]+wmp2*R[Cart::zz][Cart::s][Cart::yyy][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::yyy][3]-gfak*R[Cart::z][Cart::s][Cart::yyy][4])
R[Cart::zzz][Cart::s][Cart::xyy][3]+=pma2*R[Cart::zz][Cart::s][Cart::xyy][3]+wmp2*R[Cart::zz][Cart::s][Cart::xyy][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xyy][3]-gfak*R[Cart::z][Cart::s][Cart::xyy][4])
R[Cart::zzz][Cart::s][Cart::yyz][3]+=pma2*R[Cart::zz][Cart::s][Cart::yyz][3]+wmp2*R[Cart::zz][Cart::s][Cart::yyz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::yyz][3]-gfak*R[Cart::z][Cart::s][Cart::yyz][4])+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::yy][4]
R[Cart::zzz][Cart::s][Cart::xxy][3]+=pma2*R[Cart::zz][Cart::s][Cart::xxy][3]+wmp2*R[Cart::zz][Cart::s][Cart::xxy][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xxy][3]-gfak*R[Cart::z][Cart::s][Cart::xxy][4])
R[Cart::zzz][Cart::s][Cart::xyz][3]+=pma2*R[Cart::zz][Cart::s][Cart::xyz][3]+wmp2*R[Cart::zz][Cart::s][Cart::xyz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xyz][3]-gfak*R[Cart::z][Cart::s][Cart::xyz][4])+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::xy][4]
R[Cart::zzz][Cart::s][Cart::yzz][3]+=pma2*R[Cart::zz][Cart::s][Cart::yzz][3]+wmp2*R[Cart::zz][Cart::s][Cart::yzz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::yzz][3]-gfak*R[Cart::z][Cart::s][Cart::yzz][4])+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::yz][4]
R[Cart::zzz][Cart::s][Cart::xxx][3]+=pma2*R[Cart::zz][Cart::s][Cart::xxx][3]+wmp2*R[Cart::zz][Cart::s][Cart::xxx][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xxx][3]-gfak*R[Cart::z][Cart::s][Cart::xxx][4])
R[Cart::zzz][Cart::s][Cart::xxz][3]+=pma2*R[Cart::zz][Cart::s][Cart::xxz][3]+wmp2*R[Cart::zz][Cart::s][Cart::xxz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xxz][3]-gfak*R[Cart::z][Cart::s][Cart::xxz][4])+0.5/_decay*1*R[Cart::zz][Cart::s][Cart::xx][4]
R[Cart::zzz][Cart::s][Cart::xzz][3]+=pma2*R[Cart::zz][Cart::s][Cart::xzz][3]+wmp2*R[Cart::zz][Cart::s][Cart::xzz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::xzz][3]-gfak*R[Cart::z][Cart::s][Cart::xzz][4])+0.5/_decay*2*R[Cart::zz][Cart::s][Cart::xz][4]
R[Cart::zzz][Cart::s][Cart::zzz][3]+=pma2*R[Cart::zz][Cart::s][Cart::zzz][3]+wmp2*R[Cart::zz][Cart::s][Cart::zzz][4]+1*rzeta*(R[Cart::z][Cart::s][Cart::zzz][3]-gfak*R[Cart::z][Cart::s][Cart::zzz][4])+0.5/_decay*3*R[Cart::zz][Cart::s][Cart::zz][4]
}}
//------------------------------------------------------

//Integral g - s - s - m2
if (_mmax >4 ){
if (_lmax_alpha>3){

R[Cart::yyyy][Cart::s][Cart::s][2]+=pma1*R[Cart::yyy][Cart::s][Cart::s][2]+wmp1*R[Cart::yyy][Cart::s][Cart::s][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::s][2]-gfak*R[Cart::yy][Cart::s][Cart::s][3])
R[Cart::xyyy][Cart::s][Cart::s][2]+=pma0*R[Cart::yyy][Cart::s][Cart::s][2]+wmp0*R[Cart::yyy][Cart::s][Cart::s][3]
R[Cart::yyyz][Cart::s][Cart::s][2]+=pma2*R[Cart::yyy][Cart::s][Cart::s][2]+wmp2*R[Cart::yyy][Cart::s][Cart::s][3]
R[Cart::xxyy][Cart::s][Cart::s][2]+=pma0*R[Cart::xyy][Cart::s][Cart::s][2]+wmp0*R[Cart::xyy][Cart::s][Cart::s][3]
R[Cart::xyyz][Cart::s][Cart::s][2]+=pma0*R[Cart::yyz][Cart::s][Cart::s][2]+wmp0*R[Cart::yyz][Cart::s][Cart::s][3]
R[Cart::yyzz][Cart::s][Cart::s][2]+=pma1*R[Cart::yzz][Cart::s][Cart::s][2]+wmp1*R[Cart::yzz][Cart::s][Cart::s][3]
R[Cart::xxxy][Cart::s][Cart::s][2]+=pma1*R[Cart::xxx][Cart::s][Cart::s][2]+wmp1*R[Cart::xxx][Cart::s][Cart::s][3]
R[Cart::xxyz][Cart::s][Cart::s][2]+=pma1*R[Cart::xxz][Cart::s][Cart::s][2]+wmp1*R[Cart::xxz][Cart::s][Cart::s][3]
R[Cart::xyzz][Cart::s][Cart::s][2]+=pma0*R[Cart::yzz][Cart::s][Cart::s][2]+wmp0*R[Cart::yzz][Cart::s][Cart::s][3]
R[Cart::yzzz][Cart::s][Cart::s][2]+=pma1*R[Cart::zzz][Cart::s][Cart::s][2]+wmp1*R[Cart::zzz][Cart::s][Cart::s][3]
R[Cart::xxxx][Cart::s][Cart::s][2]+=pma0*R[Cart::xxx][Cart::s][Cart::s][2]+wmp0*R[Cart::xxx][Cart::s][Cart::s][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::s][2]-gfak*R[Cart::xx][Cart::s][Cart::s][3])
R[Cart::xxxz][Cart::s][Cart::s][2]+=pma2*R[Cart::xxx][Cart::s][Cart::s][2]+wmp2*R[Cart::xxx][Cart::s][Cart::s][3]
R[Cart::xxzz][Cart::s][Cart::s][2]+=pma0*R[Cart::xzz][Cart::s][Cart::s][2]+wmp0*R[Cart::xzz][Cart::s][Cart::s][3]
R[Cart::xzzz][Cart::s][Cart::s][2]+=pma0*R[Cart::zzz][Cart::s][Cart::s][2]+wmp0*R[Cart::zzz][Cart::s][Cart::s][3]
R[Cart::zzzz][Cart::s][Cart::s][2]+=pma2*R[Cart::zzz][Cart::s][Cart::s][2]+wmp2*R[Cart::zzz][Cart::s][Cart::s][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::s][2]-gfak*R[Cart::zz][Cart::s][Cart::s][3])
}}
//------------------------------------------------------

//Integral g - s - p - m2
if (_mmax >5 ){
if (_lmax_alpha>3 && _lmax_gamma>0){

R[Cart::yyyy][Cart::s][Cart::y][2]+=pma1*R[Cart::yyy][Cart::s][Cart::y][2]+wmp1*R[Cart::yyy][Cart::s][Cart::y][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::y][2]-gfak*R[Cart::yy][Cart::s][Cart::y][3])+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::s][3]
R[Cart::yyyy][Cart::s][Cart::x][2]+=pma1*R[Cart::yyy][Cart::s][Cart::x][2]+wmp1*R[Cart::yyy][Cart::s][Cart::x][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::x][2]-gfak*R[Cart::yy][Cart::s][Cart::x][3])
R[Cart::yyyy][Cart::s][Cart::z][2]+=pma1*R[Cart::yyy][Cart::s][Cart::z][2]+wmp1*R[Cart::yyy][Cart::s][Cart::z][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::z][2]-gfak*R[Cart::yy][Cart::s][Cart::z][3])
R[Cart::xyyy][Cart::s][Cart::y][2]+=pma0*R[Cart::yyy][Cart::s][Cart::y][2]+wmp0*R[Cart::yyy][Cart::s][Cart::y][3]
R[Cart::xyyy][Cart::s][Cart::x][2]+=pma0*R[Cart::yyy][Cart::s][Cart::x][2]+wmp0*R[Cart::yyy][Cart::s][Cart::x][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::s][3]
R[Cart::xyyy][Cart::s][Cart::z][2]+=pma0*R[Cart::yyy][Cart::s][Cart::z][2]+wmp0*R[Cart::yyy][Cart::s][Cart::z][3]
R[Cart::yyyz][Cart::s][Cart::y][2]+=pma2*R[Cart::yyy][Cart::s][Cart::y][2]+wmp2*R[Cart::yyy][Cart::s][Cart::y][3]
R[Cart::yyyz][Cart::s][Cart::x][2]+=pma2*R[Cart::yyy][Cart::s][Cart::x][2]+wmp2*R[Cart::yyy][Cart::s][Cart::x][3]
R[Cart::yyyz][Cart::s][Cart::z][2]+=pma2*R[Cart::yyy][Cart::s][Cart::z][2]+wmp2*R[Cart::yyy][Cart::s][Cart::z][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::s][3]
R[Cart::xxyy][Cart::s][Cart::y][2]+=pma0*R[Cart::xyy][Cart::s][Cart::y][2]+wmp0*R[Cart::xyy][Cart::s][Cart::y][3]
R[Cart::xxyy][Cart::s][Cart::x][2]+=pma0*R[Cart::xyy][Cart::s][Cart::x][2]+wmp0*R[Cart::xyy][Cart::s][Cart::x][3]+0.5/_decay*1*R[Cart::xyy][Cart::s][Cart::s][3]
R[Cart::xxyy][Cart::s][Cart::z][2]+=pma0*R[Cart::xyy][Cart::s][Cart::z][2]+wmp0*R[Cart::xyy][Cart::s][Cart::z][3]
R[Cart::xyyz][Cart::s][Cart::y][2]+=pma0*R[Cart::yyz][Cart::s][Cart::y][2]+wmp0*R[Cart::yyz][Cart::s][Cart::y][3]
R[Cart::xyyz][Cart::s][Cart::x][2]+=pma0*R[Cart::yyz][Cart::s][Cart::x][2]+wmp0*R[Cart::yyz][Cart::s][Cart::x][3]+0.5/_decay*1*R[Cart::yyz][Cart::s][Cart::s][3]
R[Cart::xyyz][Cart::s][Cart::z][2]+=pma0*R[Cart::yyz][Cart::s][Cart::z][2]+wmp0*R[Cart::yyz][Cart::s][Cart::z][3]
R[Cart::yyzz][Cart::s][Cart::y][2]+=pma1*R[Cart::yzz][Cart::s][Cart::y][2]+wmp1*R[Cart::yzz][Cart::s][Cart::y][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::s][3]
R[Cart::yyzz][Cart::s][Cart::x][2]+=pma1*R[Cart::yzz][Cart::s][Cart::x][2]+wmp1*R[Cart::yzz][Cart::s][Cart::x][3]
R[Cart::yyzz][Cart::s][Cart::z][2]+=pma1*R[Cart::yzz][Cart::s][Cart::z][2]+wmp1*R[Cart::yzz][Cart::s][Cart::z][3]
R[Cart::xxxy][Cart::s][Cart::y][2]+=pma1*R[Cart::xxx][Cart::s][Cart::y][2]+wmp1*R[Cart::xxx][Cart::s][Cart::y][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::s][3]
R[Cart::xxxy][Cart::s][Cart::x][2]+=pma1*R[Cart::xxx][Cart::s][Cart::x][2]+wmp1*R[Cart::xxx][Cart::s][Cart::x][3]
R[Cart::xxxy][Cart::s][Cart::z][2]+=pma1*R[Cart::xxx][Cart::s][Cart::z][2]+wmp1*R[Cart::xxx][Cart::s][Cart::z][3]
R[Cart::xxyz][Cart::s][Cart::y][2]+=pma1*R[Cart::xxz][Cart::s][Cart::y][2]+wmp1*R[Cart::xxz][Cart::s][Cart::y][3]+0.5/_decay*1*R[Cart::xxz][Cart::s][Cart::s][3]
R[Cart::xxyz][Cart::s][Cart::x][2]+=pma1*R[Cart::xxz][Cart::s][Cart::x][2]+wmp1*R[Cart::xxz][Cart::s][Cart::x][3]
R[Cart::xxyz][Cart::s][Cart::z][2]+=pma1*R[Cart::xxz][Cart::s][Cart::z][2]+wmp1*R[Cart::xxz][Cart::s][Cart::z][3]
R[Cart::xyzz][Cart::s][Cart::y][2]+=pma0*R[Cart::yzz][Cart::s][Cart::y][2]+wmp0*R[Cart::yzz][Cart::s][Cart::y][3]
R[Cart::xyzz][Cart::s][Cart::x][2]+=pma0*R[Cart::yzz][Cart::s][Cart::x][2]+wmp0*R[Cart::yzz][Cart::s][Cart::x][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::s][3]
R[Cart::xyzz][Cart::s][Cart::z][2]+=pma0*R[Cart::yzz][Cart::s][Cart::z][2]+wmp0*R[Cart::yzz][Cart::s][Cart::z][3]
R[Cart::yzzz][Cart::s][Cart::y][2]+=pma1*R[Cart::zzz][Cart::s][Cart::y][2]+wmp1*R[Cart::zzz][Cart::s][Cart::y][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::s][3]
R[Cart::yzzz][Cart::s][Cart::x][2]+=pma1*R[Cart::zzz][Cart::s][Cart::x][2]+wmp1*R[Cart::zzz][Cart::s][Cart::x][3]
R[Cart::yzzz][Cart::s][Cart::z][2]+=pma1*R[Cart::zzz][Cart::s][Cart::z][2]+wmp1*R[Cart::zzz][Cart::s][Cart::z][3]
R[Cart::xxxx][Cart::s][Cart::y][2]+=pma0*R[Cart::xxx][Cart::s][Cart::y][2]+wmp0*R[Cart::xxx][Cart::s][Cart::y][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::y][2]-gfak*R[Cart::xx][Cart::s][Cart::y][3])
R[Cart::xxxx][Cart::s][Cart::x][2]+=pma0*R[Cart::xxx][Cart::s][Cart::x][2]+wmp0*R[Cart::xxx][Cart::s][Cart::x][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::x][2]-gfak*R[Cart::xx][Cart::s][Cart::x][3])+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::s][3]
R[Cart::xxxx][Cart::s][Cart::z][2]+=pma0*R[Cart::xxx][Cart::s][Cart::z][2]+wmp0*R[Cart::xxx][Cart::s][Cart::z][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::z][2]-gfak*R[Cart::xx][Cart::s][Cart::z][3])
R[Cart::xxxz][Cart::s][Cart::y][2]+=pma2*R[Cart::xxx][Cart::s][Cart::y][2]+wmp2*R[Cart::xxx][Cart::s][Cart::y][3]
R[Cart::xxxz][Cart::s][Cart::x][2]+=pma2*R[Cart::xxx][Cart::s][Cart::x][2]+wmp2*R[Cart::xxx][Cart::s][Cart::x][3]
R[Cart::xxxz][Cart::s][Cart::z][2]+=pma2*R[Cart::xxx][Cart::s][Cart::z][2]+wmp2*R[Cart::xxx][Cart::s][Cart::z][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::s][3]
R[Cart::xxzz][Cart::s][Cart::y][2]+=pma0*R[Cart::xzz][Cart::s][Cart::y][2]+wmp0*R[Cart::xzz][Cart::s][Cart::y][3]
R[Cart::xxzz][Cart::s][Cart::x][2]+=pma0*R[Cart::xzz][Cart::s][Cart::x][2]+wmp0*R[Cart::xzz][Cart::s][Cart::x][3]+0.5/_decay*1*R[Cart::xzz][Cart::s][Cart::s][3]
R[Cart::xxzz][Cart::s][Cart::z][2]+=pma0*R[Cart::xzz][Cart::s][Cart::z][2]+wmp0*R[Cart::xzz][Cart::s][Cart::z][3]
R[Cart::xzzz][Cart::s][Cart::y][2]+=pma0*R[Cart::zzz][Cart::s][Cart::y][2]+wmp0*R[Cart::zzz][Cart::s][Cart::y][3]
R[Cart::xzzz][Cart::s][Cart::x][2]+=pma0*R[Cart::zzz][Cart::s][Cart::x][2]+wmp0*R[Cart::zzz][Cart::s][Cart::x][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::s][3]
R[Cart::xzzz][Cart::s][Cart::z][2]+=pma0*R[Cart::zzz][Cart::s][Cart::z][2]+wmp0*R[Cart::zzz][Cart::s][Cart::z][3]
R[Cart::zzzz][Cart::s][Cart::y][2]+=pma2*R[Cart::zzz][Cart::s][Cart::y][2]+wmp2*R[Cart::zzz][Cart::s][Cart::y][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::y][2]-gfak*R[Cart::zz][Cart::s][Cart::y][3])
R[Cart::zzzz][Cart::s][Cart::x][2]+=pma2*R[Cart::zzz][Cart::s][Cart::x][2]+wmp2*R[Cart::zzz][Cart::s][Cart::x][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::x][2]-gfak*R[Cart::zz][Cart::s][Cart::x][3])
R[Cart::zzzz][Cart::s][Cart::z][2]+=pma2*R[Cart::zzz][Cart::s][Cart::z][2]+wmp2*R[Cart::zzz][Cart::s][Cart::z][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::z][2]-gfak*R[Cart::zz][Cart::s][Cart::z][3])+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::s][3]
}}
//------------------------------------------------------

//Integral g - s - d - m2
if (_mmax >6 ){
if (_lmax_alpha>3 && _lmax_gamma>1){

R[Cart::yyyy][Cart::s][Cart::yy][2]+=pma1*R[Cart::yyy][Cart::s][Cart::yy][2]+wmp1*R[Cart::yyy][Cart::s][Cart::yy][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::yy][2]-gfak*R[Cart::yy][Cart::s][Cart::yy][3])+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::y][3]
R[Cart::yyyy][Cart::s][Cart::xy][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xy][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xy][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xy][2]-gfak*R[Cart::yy][Cart::s][Cart::xy][3])+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::x][3]
R[Cart::yyyy][Cart::s][Cart::yz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::yz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::yz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::yz][2]-gfak*R[Cart::yy][Cart::s][Cart::yz][3])+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::z][3]
R[Cart::yyyy][Cart::s][Cart::xx][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xx][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xx][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xx][2]-gfak*R[Cart::yy][Cart::s][Cart::xx][3])
R[Cart::yyyy][Cart::s][Cart::xz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xz][2]-gfak*R[Cart::yy][Cart::s][Cart::xz][3])
R[Cart::yyyy][Cart::s][Cart::zz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::zz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::zz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::zz][2]-gfak*R[Cart::yy][Cart::s][Cart::zz][3])
R[Cart::xyyy][Cart::s][Cart::yy][2]+=pma0*R[Cart::yyy][Cart::s][Cart::yy][2]+wmp0*R[Cart::yyy][Cart::s][Cart::yy][3]
R[Cart::xyyy][Cart::s][Cart::xy][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xy][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::y][3]
R[Cart::xyyy][Cart::s][Cart::yz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::yz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::yz][3]
R[Cart::xyyy][Cart::s][Cart::xx][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xx][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xx][3]+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::x][3]
R[Cart::xyyy][Cart::s][Cart::xz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::z][3]
R[Cart::xyyy][Cart::s][Cart::zz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::zz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::zz][3]
R[Cart::yyyz][Cart::s][Cart::yy][2]+=pma2*R[Cart::yyy][Cart::s][Cart::yy][2]+wmp2*R[Cart::yyy][Cart::s][Cart::yy][3]
R[Cart::yyyz][Cart::s][Cart::xy][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xy][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xy][3]
R[Cart::yyyz][Cart::s][Cart::yz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::yz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::yz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::y][3]
R[Cart::yyyz][Cart::s][Cart::xx][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xx][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xx][3]
R[Cart::yyyz][Cart::s][Cart::xz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::x][3]
R[Cart::yyyz][Cart::s][Cart::zz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::zz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::zz][3]+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::z][3]
R[Cart::xxyy][Cart::s][Cart::yy][2]+=pma0*R[Cart::xyy][Cart::s][Cart::yy][2]+wmp0*R[Cart::xyy][Cart::s][Cart::yy][3]
R[Cart::xxyy][Cart::s][Cart::xy][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xy][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::xyy][Cart::s][Cart::y][3]
R[Cart::xxyy][Cart::s][Cart::yz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::yz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::yz][3]
R[Cart::xxyy][Cart::s][Cart::xx][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xx][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xx][3]+0.5/_decay*2*R[Cart::xyy][Cart::s][Cart::x][3]
R[Cart::xxyy][Cart::s][Cart::xz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::xyy][Cart::s][Cart::z][3]
R[Cart::xxyy][Cart::s][Cart::zz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::zz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::zz][3]
R[Cart::xyyz][Cart::s][Cart::yy][2]+=pma0*R[Cart::yyz][Cart::s][Cart::yy][2]+wmp0*R[Cart::yyz][Cart::s][Cart::yy][3]
R[Cart::xyyz][Cart::s][Cart::xy][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xy][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::yyz][Cart::s][Cart::y][3]
R[Cart::xyyz][Cart::s][Cart::yz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::yz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::yz][3]
R[Cart::xyyz][Cart::s][Cart::xx][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xx][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xx][3]+0.5/_decay*2*R[Cart::yyz][Cart::s][Cart::x][3]
R[Cart::xyyz][Cart::s][Cart::xz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::yyz][Cart::s][Cart::z][3]
R[Cart::xyyz][Cart::s][Cart::zz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::zz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::zz][3]
R[Cart::yyzz][Cart::s][Cart::yy][2]+=pma1*R[Cart::yzz][Cart::s][Cart::yy][2]+wmp1*R[Cart::yzz][Cart::s][Cart::yy][3]+0.5/_decay*2*R[Cart::yzz][Cart::s][Cart::y][3]
R[Cart::yyzz][Cart::s][Cart::xy][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xy][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::x][3]
R[Cart::yyzz][Cart::s][Cart::yz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::yz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::yz][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::z][3]
R[Cart::yyzz][Cart::s][Cart::xx][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xx][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xx][3]
R[Cart::yyzz][Cart::s][Cart::xz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xz][3]
R[Cart::yyzz][Cart::s][Cart::zz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::zz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::zz][3]
R[Cart::xxxy][Cart::s][Cart::yy][2]+=pma1*R[Cart::xxx][Cart::s][Cart::yy][2]+wmp1*R[Cart::xxx][Cart::s][Cart::yy][3]+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::y][3]
R[Cart::xxxy][Cart::s][Cart::xy][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xy][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::x][3]
R[Cart::xxxy][Cart::s][Cart::yz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::yz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::yz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::z][3]
R[Cart::xxxy][Cart::s][Cart::xx][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xx][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xx][3]
R[Cart::xxxy][Cart::s][Cart::xz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xz][3]
R[Cart::xxxy][Cart::s][Cart::zz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::zz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::zz][3]
R[Cart::xxyz][Cart::s][Cart::yy][2]+=pma1*R[Cart::xxz][Cart::s][Cart::yy][2]+wmp1*R[Cart::xxz][Cart::s][Cart::yy][3]+0.5/_decay*2*R[Cart::xxz][Cart::s][Cart::y][3]
R[Cart::xxyz][Cart::s][Cart::xy][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xy][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::xxz][Cart::s][Cart::x][3]
R[Cart::xxyz][Cart::s][Cart::yz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::yz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::yz][3]+0.5/_decay*1*R[Cart::xxz][Cart::s][Cart::z][3]
R[Cart::xxyz][Cart::s][Cart::xx][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xx][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xx][3]
R[Cart::xxyz][Cart::s][Cart::xz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xz][3]
R[Cart::xxyz][Cart::s][Cart::zz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::zz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::zz][3]
R[Cart::xyzz][Cart::s][Cart::yy][2]+=pma0*R[Cart::yzz][Cart::s][Cart::yy][2]+wmp0*R[Cart::yzz][Cart::s][Cart::yy][3]
R[Cart::xyzz][Cart::s][Cart::xy][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xy][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::y][3]
R[Cart::xyzz][Cart::s][Cart::yz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::yz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::yz][3]
R[Cart::xyzz][Cart::s][Cart::xx][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xx][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xx][3]+0.5/_decay*2*R[Cart::yzz][Cart::s][Cart::x][3]
R[Cart::xyzz][Cart::s][Cart::xz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::z][3]
R[Cart::xyzz][Cart::s][Cart::zz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::zz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::zz][3]
R[Cart::yzzz][Cart::s][Cart::yy][2]+=pma1*R[Cart::zzz][Cart::s][Cart::yy][2]+wmp1*R[Cart::zzz][Cart::s][Cart::yy][3]+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::y][3]
R[Cart::yzzz][Cart::s][Cart::xy][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xy][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::x][3]
R[Cart::yzzz][Cart::s][Cart::yz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::yz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::yz][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::z][3]
R[Cart::yzzz][Cart::s][Cart::xx][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xx][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xx][3]
R[Cart::yzzz][Cart::s][Cart::xz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xz][3]
R[Cart::yzzz][Cart::s][Cart::zz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::zz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::zz][3]
R[Cart::xxxx][Cart::s][Cart::yy][2]+=pma0*R[Cart::xxx][Cart::s][Cart::yy][2]+wmp0*R[Cart::xxx][Cart::s][Cart::yy][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::yy][2]-gfak*R[Cart::xx][Cart::s][Cart::yy][3])
R[Cart::xxxx][Cart::s][Cart::xy][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xy][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xy][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xy][2]-gfak*R[Cart::xx][Cart::s][Cart::xy][3])+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::y][3]
R[Cart::xxxx][Cart::s][Cart::yz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::yz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::yz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::yz][2]-gfak*R[Cart::xx][Cart::s][Cart::yz][3])
R[Cart::xxxx][Cart::s][Cart::xx][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xx][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xx][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xx][2]-gfak*R[Cart::xx][Cart::s][Cart::xx][3])+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::x][3]
R[Cart::xxxx][Cart::s][Cart::xz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xz][2]-gfak*R[Cart::xx][Cart::s][Cart::xz][3])+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::z][3]
R[Cart::xxxx][Cart::s][Cart::zz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::zz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::zz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::zz][2]-gfak*R[Cart::xx][Cart::s][Cart::zz][3])
R[Cart::xxxz][Cart::s][Cart::yy][2]+=pma2*R[Cart::xxx][Cart::s][Cart::yy][2]+wmp2*R[Cart::xxx][Cart::s][Cart::yy][3]
R[Cart::xxxz][Cart::s][Cart::xy][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xy][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xy][3]
R[Cart::xxxz][Cart::s][Cart::yz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::yz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::yz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::y][3]
R[Cart::xxxz][Cart::s][Cart::xx][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xx][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xx][3]
R[Cart::xxxz][Cart::s][Cart::xz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::x][3]
R[Cart::xxxz][Cart::s][Cart::zz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::zz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::zz][3]+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::z][3]
R[Cart::xxzz][Cart::s][Cart::yy][2]+=pma0*R[Cart::xzz][Cart::s][Cart::yy][2]+wmp0*R[Cart::xzz][Cart::s][Cart::yy][3]
R[Cart::xxzz][Cart::s][Cart::xy][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xy][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::xzz][Cart::s][Cart::y][3]
R[Cart::xxzz][Cart::s][Cart::yz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::yz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::yz][3]
R[Cart::xxzz][Cart::s][Cart::xx][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xx][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xx][3]+0.5/_decay*2*R[Cart::xzz][Cart::s][Cart::x][3]
R[Cart::xxzz][Cart::s][Cart::xz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::xzz][Cart::s][Cart::z][3]
R[Cart::xxzz][Cart::s][Cart::zz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::zz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::zz][3]
R[Cart::xzzz][Cart::s][Cart::yy][2]+=pma0*R[Cart::zzz][Cart::s][Cart::yy][2]+wmp0*R[Cart::zzz][Cart::s][Cart::yy][3]
R[Cart::xzzz][Cart::s][Cart::xy][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xy][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xy][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::y][3]
R[Cart::xzzz][Cart::s][Cart::yz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::yz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::yz][3]
R[Cart::xzzz][Cart::s][Cart::xx][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xx][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xx][3]+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::x][3]
R[Cart::xzzz][Cart::s][Cart::xz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xz][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::z][3]
R[Cart::xzzz][Cart::s][Cart::zz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::zz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::zz][3]
R[Cart::zzzz][Cart::s][Cart::yy][2]+=pma2*R[Cart::zzz][Cart::s][Cart::yy][2]+wmp2*R[Cart::zzz][Cart::s][Cart::yy][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::yy][2]-gfak*R[Cart::zz][Cart::s][Cart::yy][3])
R[Cart::zzzz][Cart::s][Cart::xy][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xy][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xy][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xy][2]-gfak*R[Cart::zz][Cart::s][Cart::xy][3])
R[Cart::zzzz][Cart::s][Cart::yz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::yz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::yz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::yz][2]-gfak*R[Cart::zz][Cart::s][Cart::yz][3])+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::y][3]
R[Cart::zzzz][Cart::s][Cart::xx][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xx][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xx][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xx][2]-gfak*R[Cart::zz][Cart::s][Cart::xx][3])
R[Cart::zzzz][Cart::s][Cart::xz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xz][2]-gfak*R[Cart::zz][Cart::s][Cart::xz][3])+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::x][3]
R[Cart::zzzz][Cart::s][Cart::zz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::zz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::zz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::zz][2]-gfak*R[Cart::zz][Cart::s][Cart::zz][3])+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::z][3]
}}
//------------------------------------------------------

//Integral g - s - f - m2
if (_mmax >7 ){
if (_lmax_alpha>3 && _lmax_gamma>2){

R[Cart::yyyy][Cart::s][Cart::yyy][2]+=pma1*R[Cart::yyy][Cart::s][Cart::yyy][2]+wmp1*R[Cart::yyy][Cart::s][Cart::yyy][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::yyy][2]-gfak*R[Cart::yy][Cart::s][Cart::yyy][3])+0.5/_decay*3*R[Cart::yyy][Cart::s][Cart::yy][3]
R[Cart::yyyy][Cart::s][Cart::xyy][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xyy][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xyy][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xyy][2]-gfak*R[Cart::yy][Cart::s][Cart::xyy][3])+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::xy][3]
R[Cart::yyyy][Cart::s][Cart::yyz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::yyz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::yyz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::yyz][2]-gfak*R[Cart::yy][Cart::s][Cart::yyz][3])+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::yz][3]
R[Cart::yyyy][Cart::s][Cart::xxy][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xxy][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xxy][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xxy][2]-gfak*R[Cart::yy][Cart::s][Cart::xxy][3])+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::xx][3]
R[Cart::yyyy][Cart::s][Cart::xyz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xyz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xyz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xyz][2]-gfak*R[Cart::yy][Cart::s][Cart::xyz][3])+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::xz][3]
R[Cart::yyyy][Cart::s][Cart::yzz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::yzz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::yzz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::yzz][2]-gfak*R[Cart::yy][Cart::s][Cart::yzz][3])+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::zz][3]
R[Cart::yyyy][Cart::s][Cart::xxx][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xxx][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xxx][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xxx][2]-gfak*R[Cart::yy][Cart::s][Cart::xxx][3])
R[Cart::yyyy][Cart::s][Cart::xxz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xxz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xxz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xxz][2]-gfak*R[Cart::yy][Cart::s][Cart::xxz][3])
R[Cart::yyyy][Cart::s][Cart::xzz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::xzz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::xzz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::xzz][2]-gfak*R[Cart::yy][Cart::s][Cart::xzz][3])
R[Cart::yyyy][Cart::s][Cart::zzz][2]+=pma1*R[Cart::yyy][Cart::s][Cart::zzz][2]+wmp1*R[Cart::yyy][Cart::s][Cart::zzz][3]+2*rzeta*(R[Cart::yy][Cart::s][Cart::zzz][2]-gfak*R[Cart::yy][Cart::s][Cart::zzz][3])
R[Cart::xyyy][Cart::s][Cart::yyy][2]+=pma0*R[Cart::yyy][Cart::s][Cart::yyy][2]+wmp0*R[Cart::yyy][Cart::s][Cart::yyy][3]
R[Cart::xyyy][Cart::s][Cart::xyy][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xyy][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xyy][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::yy][3]
R[Cart::xyyy][Cart::s][Cart::yyz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::yyz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::yyz][3]
R[Cart::xyyy][Cart::s][Cart::xxy][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xxy][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xxy][3]+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::xy][3]
R[Cart::xyyy][Cart::s][Cart::xyz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xyz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::yz][3]
R[Cart::xyyy][Cart::s][Cart::yzz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::yzz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::yzz][3]
R[Cart::xyyy][Cart::s][Cart::xxx][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xxx][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xxx][3]+0.5/_decay*3*R[Cart::yyy][Cart::s][Cart::xx][3]
R[Cart::xyyy][Cart::s][Cart::xxz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xxz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xxz][3]+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::xz][3]
R[Cart::xyyy][Cart::s][Cart::xzz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::xzz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::xzz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::zz][3]
R[Cart::xyyy][Cart::s][Cart::zzz][2]+=pma0*R[Cart::yyy][Cart::s][Cart::zzz][2]+wmp0*R[Cart::yyy][Cart::s][Cart::zzz][3]
R[Cart::yyyz][Cart::s][Cart::yyy][2]+=pma2*R[Cart::yyy][Cart::s][Cart::yyy][2]+wmp2*R[Cart::yyy][Cart::s][Cart::yyy][3]
R[Cart::yyyz][Cart::s][Cart::xyy][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xyy][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xyy][3]
R[Cart::yyyz][Cart::s][Cart::yyz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::yyz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::yyz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::yy][3]
R[Cart::yyyz][Cart::s][Cart::xxy][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xxy][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xxy][3]
R[Cart::yyyz][Cart::s][Cart::xyz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xyz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::xy][3]
R[Cart::yyyz][Cart::s][Cart::yzz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::yzz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::yzz][3]+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::yz][3]
R[Cart::yyyz][Cart::s][Cart::xxx][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xxx][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xxx][3]
R[Cart::yyyz][Cart::s][Cart::xxz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xxz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xxz][3]+0.5/_decay*1*R[Cart::yyy][Cart::s][Cart::xx][3]
R[Cart::yyyz][Cart::s][Cart::xzz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::xzz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::xzz][3]+0.5/_decay*2*R[Cart::yyy][Cart::s][Cart::xz][3]
R[Cart::yyyz][Cart::s][Cart::zzz][2]+=pma2*R[Cart::yyy][Cart::s][Cart::zzz][2]+wmp2*R[Cart::yyy][Cart::s][Cart::zzz][3]+0.5/_decay*3*R[Cart::yyy][Cart::s][Cart::zz][3]
R[Cart::xxyy][Cart::s][Cart::yyy][2]+=pma0*R[Cart::xyy][Cart::s][Cart::yyy][2]+wmp0*R[Cart::xyy][Cart::s][Cart::yyy][3]
R[Cart::xxyy][Cart::s][Cart::xyy][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xyy][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xyy][3]+0.5/_decay*1*R[Cart::xyy][Cart::s][Cart::yy][3]
R[Cart::xxyy][Cart::s][Cart::yyz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::yyz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::yyz][3]
R[Cart::xxyy][Cart::s][Cart::xxy][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xxy][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xxy][3]+0.5/_decay*2*R[Cart::xyy][Cart::s][Cart::xy][3]
R[Cart::xxyy][Cart::s][Cart::xyz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xyz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::xyy][Cart::s][Cart::yz][3]
R[Cart::xxyy][Cart::s][Cart::yzz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::yzz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::yzz][3]
R[Cart::xxyy][Cart::s][Cart::xxx][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xxx][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xxx][3]+0.5/_decay*3*R[Cart::xyy][Cart::s][Cart::xx][3]
R[Cart::xxyy][Cart::s][Cart::xxz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xxz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xxz][3]+0.5/_decay*2*R[Cart::xyy][Cart::s][Cart::xz][3]
R[Cart::xxyy][Cart::s][Cart::xzz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::xzz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::xzz][3]+0.5/_decay*1*R[Cart::xyy][Cart::s][Cart::zz][3]
R[Cart::xxyy][Cart::s][Cart::zzz][2]+=pma0*R[Cart::xyy][Cart::s][Cart::zzz][2]+wmp0*R[Cart::xyy][Cart::s][Cart::zzz][3]
R[Cart::xyyz][Cart::s][Cart::yyy][2]+=pma0*R[Cart::yyz][Cart::s][Cart::yyy][2]+wmp0*R[Cart::yyz][Cart::s][Cart::yyy][3]
R[Cart::xyyz][Cart::s][Cart::xyy][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xyy][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xyy][3]+0.5/_decay*1*R[Cart::yyz][Cart::s][Cart::yy][3]
R[Cart::xyyz][Cart::s][Cart::yyz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::yyz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::yyz][3]
R[Cart::xyyz][Cart::s][Cart::xxy][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xxy][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xxy][3]+0.5/_decay*2*R[Cart::yyz][Cart::s][Cart::xy][3]
R[Cart::xyyz][Cart::s][Cart::xyz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xyz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::yyz][Cart::s][Cart::yz][3]
R[Cart::xyyz][Cart::s][Cart::yzz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::yzz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::yzz][3]
R[Cart::xyyz][Cart::s][Cart::xxx][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xxx][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xxx][3]+0.5/_decay*3*R[Cart::yyz][Cart::s][Cart::xx][3]
R[Cart::xyyz][Cart::s][Cart::xxz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xxz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xxz][3]+0.5/_decay*2*R[Cart::yyz][Cart::s][Cart::xz][3]
R[Cart::xyyz][Cart::s][Cart::xzz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::xzz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::xzz][3]+0.5/_decay*1*R[Cart::yyz][Cart::s][Cart::zz][3]
R[Cart::xyyz][Cart::s][Cart::zzz][2]+=pma0*R[Cart::yyz][Cart::s][Cart::zzz][2]+wmp0*R[Cart::yyz][Cart::s][Cart::zzz][3]
R[Cart::yyzz][Cart::s][Cart::yyy][2]+=pma1*R[Cart::yzz][Cart::s][Cart::yyy][2]+wmp1*R[Cart::yzz][Cart::s][Cart::yyy][3]+0.5/_decay*3*R[Cart::yzz][Cart::s][Cart::yy][3]
R[Cart::yyzz][Cart::s][Cart::xyy][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xyy][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xyy][3]+0.5/_decay*2*R[Cart::yzz][Cart::s][Cart::xy][3]
R[Cart::yyzz][Cart::s][Cart::yyz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::yyz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::yyz][3]+0.5/_decay*2*R[Cart::yzz][Cart::s][Cart::yz][3]
R[Cart::yyzz][Cart::s][Cart::xxy][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xxy][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xxy][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::xx][3]
R[Cart::yyzz][Cart::s][Cart::xyz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xyz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::xz][3]
R[Cart::yyzz][Cart::s][Cart::yzz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::yzz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::yzz][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::zz][3]
R[Cart::yyzz][Cart::s][Cart::xxx][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xxx][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xxx][3]
R[Cart::yyzz][Cart::s][Cart::xxz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xxz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xxz][3]
R[Cart::yyzz][Cart::s][Cart::xzz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::xzz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::xzz][3]
R[Cart::yyzz][Cart::s][Cart::zzz][2]+=pma1*R[Cart::yzz][Cart::s][Cart::zzz][2]+wmp1*R[Cart::yzz][Cart::s][Cart::zzz][3]
R[Cart::xxxy][Cart::s][Cart::yyy][2]+=pma1*R[Cart::xxx][Cart::s][Cart::yyy][2]+wmp1*R[Cart::xxx][Cart::s][Cart::yyy][3]+0.5/_decay*3*R[Cart::xxx][Cart::s][Cart::yy][3]
R[Cart::xxxy][Cart::s][Cart::xyy][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xyy][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xyy][3]+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::xy][3]
R[Cart::xxxy][Cart::s][Cart::yyz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::yyz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::yyz][3]+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::yz][3]
R[Cart::xxxy][Cart::s][Cart::xxy][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xxy][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xxy][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::xx][3]
R[Cart::xxxy][Cart::s][Cart::xyz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xyz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::xz][3]
R[Cart::xxxy][Cart::s][Cart::yzz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::yzz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::yzz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::zz][3]
R[Cart::xxxy][Cart::s][Cart::xxx][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xxx][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xxx][3]
R[Cart::xxxy][Cart::s][Cart::xxz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xxz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xxz][3]
R[Cart::xxxy][Cart::s][Cart::xzz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::xzz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::xzz][3]
R[Cart::xxxy][Cart::s][Cart::zzz][2]+=pma1*R[Cart::xxx][Cart::s][Cart::zzz][2]+wmp1*R[Cart::xxx][Cart::s][Cart::zzz][3]
R[Cart::xxyz][Cart::s][Cart::yyy][2]+=pma1*R[Cart::xxz][Cart::s][Cart::yyy][2]+wmp1*R[Cart::xxz][Cart::s][Cart::yyy][3]+0.5/_decay*3*R[Cart::xxz][Cart::s][Cart::yy][3]
R[Cart::xxyz][Cart::s][Cart::xyy][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xyy][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xyy][3]+0.5/_decay*2*R[Cart::xxz][Cart::s][Cart::xy][3]
R[Cart::xxyz][Cart::s][Cart::yyz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::yyz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::yyz][3]+0.5/_decay*2*R[Cart::xxz][Cart::s][Cart::yz][3]
R[Cart::xxyz][Cart::s][Cart::xxy][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xxy][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xxy][3]+0.5/_decay*1*R[Cart::xxz][Cart::s][Cart::xx][3]
R[Cart::xxyz][Cart::s][Cart::xyz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xyz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::xxz][Cart::s][Cart::xz][3]
R[Cart::xxyz][Cart::s][Cart::yzz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::yzz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::yzz][3]+0.5/_decay*1*R[Cart::xxz][Cart::s][Cart::zz][3]
R[Cart::xxyz][Cart::s][Cart::xxx][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xxx][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xxx][3]
R[Cart::xxyz][Cart::s][Cart::xxz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xxz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xxz][3]
R[Cart::xxyz][Cart::s][Cart::xzz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::xzz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::xzz][3]
R[Cart::xxyz][Cart::s][Cart::zzz][2]+=pma1*R[Cart::xxz][Cart::s][Cart::zzz][2]+wmp1*R[Cart::xxz][Cart::s][Cart::zzz][3]
R[Cart::xyzz][Cart::s][Cart::yyy][2]+=pma0*R[Cart::yzz][Cart::s][Cart::yyy][2]+wmp0*R[Cart::yzz][Cart::s][Cart::yyy][3]
R[Cart::xyzz][Cart::s][Cart::xyy][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xyy][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xyy][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::yy][3]
R[Cart::xyzz][Cart::s][Cart::yyz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::yyz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::yyz][3]
R[Cart::xyzz][Cart::s][Cart::xxy][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xxy][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xxy][3]+0.5/_decay*2*R[Cart::yzz][Cart::s][Cart::xy][3]
R[Cart::xyzz][Cart::s][Cart::xyz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xyz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::yz][3]
R[Cart::xyzz][Cart::s][Cart::yzz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::yzz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::yzz][3]
R[Cart::xyzz][Cart::s][Cart::xxx][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xxx][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xxx][3]+0.5/_decay*3*R[Cart::yzz][Cart::s][Cart::xx][3]
R[Cart::xyzz][Cart::s][Cart::xxz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xxz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xxz][3]+0.5/_decay*2*R[Cart::yzz][Cart::s][Cart::xz][3]
R[Cart::xyzz][Cart::s][Cart::xzz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::xzz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::xzz][3]+0.5/_decay*1*R[Cart::yzz][Cart::s][Cart::zz][3]
R[Cart::xyzz][Cart::s][Cart::zzz][2]+=pma0*R[Cart::yzz][Cart::s][Cart::zzz][2]+wmp0*R[Cart::yzz][Cart::s][Cart::zzz][3]
R[Cart::yzzz][Cart::s][Cart::yyy][2]+=pma1*R[Cart::zzz][Cart::s][Cart::yyy][2]+wmp1*R[Cart::zzz][Cart::s][Cart::yyy][3]+0.5/_decay*3*R[Cart::zzz][Cart::s][Cart::yy][3]
R[Cart::yzzz][Cart::s][Cart::xyy][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xyy][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xyy][3]+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::xy][3]
R[Cart::yzzz][Cart::s][Cart::yyz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::yyz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::yyz][3]+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::yz][3]
R[Cart::yzzz][Cart::s][Cart::xxy][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xxy][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xxy][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::xx][3]
R[Cart::yzzz][Cart::s][Cart::xyz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xyz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::xz][3]
R[Cart::yzzz][Cart::s][Cart::yzz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::yzz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::yzz][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::zz][3]
R[Cart::yzzz][Cart::s][Cart::xxx][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xxx][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xxx][3]
R[Cart::yzzz][Cart::s][Cart::xxz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xxz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xxz][3]
R[Cart::yzzz][Cart::s][Cart::xzz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::xzz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::xzz][3]
R[Cart::yzzz][Cart::s][Cart::zzz][2]+=pma1*R[Cart::zzz][Cart::s][Cart::zzz][2]+wmp1*R[Cart::zzz][Cart::s][Cart::zzz][3]
R[Cart::xxxx][Cart::s][Cart::yyy][2]+=pma0*R[Cart::xxx][Cart::s][Cart::yyy][2]+wmp0*R[Cart::xxx][Cart::s][Cart::yyy][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::yyy][2]-gfak*R[Cart::xx][Cart::s][Cart::yyy][3])
R[Cart::xxxx][Cart::s][Cart::xyy][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xyy][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xyy][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xyy][2]-gfak*R[Cart::xx][Cart::s][Cart::xyy][3])+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::yy][3]
R[Cart::xxxx][Cart::s][Cart::yyz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::yyz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::yyz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::yyz][2]-gfak*R[Cart::xx][Cart::s][Cart::yyz][3])
R[Cart::xxxx][Cart::s][Cart::xxy][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xxy][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xxy][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xxy][2]-gfak*R[Cart::xx][Cart::s][Cart::xxy][3])+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::xy][3]
R[Cart::xxxx][Cart::s][Cart::xyz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xyz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xyz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xyz][2]-gfak*R[Cart::xx][Cart::s][Cart::xyz][3])+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::yz][3]
R[Cart::xxxx][Cart::s][Cart::yzz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::yzz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::yzz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::yzz][2]-gfak*R[Cart::xx][Cart::s][Cart::yzz][3])
R[Cart::xxxx][Cart::s][Cart::xxx][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xxx][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xxx][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xxx][2]-gfak*R[Cart::xx][Cart::s][Cart::xxx][3])+0.5/_decay*3*R[Cart::xxx][Cart::s][Cart::xx][3]
R[Cart::xxxx][Cart::s][Cart::xxz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xxz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xxz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xxz][2]-gfak*R[Cart::xx][Cart::s][Cart::xxz][3])+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::xz][3]
R[Cart::xxxx][Cart::s][Cart::xzz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::xzz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::xzz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::xzz][2]-gfak*R[Cart::xx][Cart::s][Cart::xzz][3])+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::zz][3]
R[Cart::xxxx][Cart::s][Cart::zzz][2]+=pma0*R[Cart::xxx][Cart::s][Cart::zzz][2]+wmp0*R[Cart::xxx][Cart::s][Cart::zzz][3]+2*rzeta*(R[Cart::xx][Cart::s][Cart::zzz][2]-gfak*R[Cart::xx][Cart::s][Cart::zzz][3])
R[Cart::xxxz][Cart::s][Cart::yyy][2]+=pma2*R[Cart::xxx][Cart::s][Cart::yyy][2]+wmp2*R[Cart::xxx][Cart::s][Cart::yyy][3]
R[Cart::xxxz][Cart::s][Cart::xyy][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xyy][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xyy][3]
R[Cart::xxxz][Cart::s][Cart::yyz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::yyz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::yyz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::yy][3]
R[Cart::xxxz][Cart::s][Cart::xxy][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xxy][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xxy][3]
R[Cart::xxxz][Cart::s][Cart::xyz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xyz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::xy][3]
R[Cart::xxxz][Cart::s][Cart::yzz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::yzz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::yzz][3]+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::yz][3]
R[Cart::xxxz][Cart::s][Cart::xxx][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xxx][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xxx][3]
R[Cart::xxxz][Cart::s][Cart::xxz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xxz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xxz][3]+0.5/_decay*1*R[Cart::xxx][Cart::s][Cart::xx][3]
R[Cart::xxxz][Cart::s][Cart::xzz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::xzz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::xzz][3]+0.5/_decay*2*R[Cart::xxx][Cart::s][Cart::xz][3]
R[Cart::xxxz][Cart::s][Cart::zzz][2]+=pma2*R[Cart::xxx][Cart::s][Cart::zzz][2]+wmp2*R[Cart::xxx][Cart::s][Cart::zzz][3]+0.5/_decay*3*R[Cart::xxx][Cart::s][Cart::zz][3]
R[Cart::xxzz][Cart::s][Cart::yyy][2]+=pma0*R[Cart::xzz][Cart::s][Cart::yyy][2]+wmp0*R[Cart::xzz][Cart::s][Cart::yyy][3]
R[Cart::xxzz][Cart::s][Cart::xyy][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xyy][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xyy][3]+0.5/_decay*1*R[Cart::xzz][Cart::s][Cart::yy][3]
R[Cart::xxzz][Cart::s][Cart::yyz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::yyz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::yyz][3]
R[Cart::xxzz][Cart::s][Cart::xxy][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xxy][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xxy][3]+0.5/_decay*2*R[Cart::xzz][Cart::s][Cart::xy][3]
R[Cart::xxzz][Cart::s][Cart::xyz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xyz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::xzz][Cart::s][Cart::yz][3]
R[Cart::xxzz][Cart::s][Cart::yzz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::yzz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::yzz][3]
R[Cart::xxzz][Cart::s][Cart::xxx][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xxx][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xxx][3]+0.5/_decay*3*R[Cart::xzz][Cart::s][Cart::xx][3]
R[Cart::xxzz][Cart::s][Cart::xxz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xxz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xxz][3]+0.5/_decay*2*R[Cart::xzz][Cart::s][Cart::xz][3]
R[Cart::xxzz][Cart::s][Cart::xzz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::xzz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::xzz][3]+0.5/_decay*1*R[Cart::xzz][Cart::s][Cart::zz][3]
R[Cart::xxzz][Cart::s][Cart::zzz][2]+=pma0*R[Cart::xzz][Cart::s][Cart::zzz][2]+wmp0*R[Cart::xzz][Cart::s][Cart::zzz][3]
R[Cart::xzzz][Cart::s][Cart::yyy][2]+=pma0*R[Cart::zzz][Cart::s][Cart::yyy][2]+wmp0*R[Cart::zzz][Cart::s][Cart::yyy][3]
R[Cart::xzzz][Cart::s][Cart::xyy][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xyy][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xyy][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::yy][3]
R[Cart::xzzz][Cart::s][Cart::yyz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::yyz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::yyz][3]
R[Cart::xzzz][Cart::s][Cart::xxy][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xxy][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xxy][3]+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::xy][3]
R[Cart::xzzz][Cart::s][Cart::xyz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xyz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xyz][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::yz][3]
R[Cart::xzzz][Cart::s][Cart::yzz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::yzz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::yzz][3]
R[Cart::xzzz][Cart::s][Cart::xxx][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xxx][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xxx][3]+0.5/_decay*3*R[Cart::zzz][Cart::s][Cart::xx][3]
R[Cart::xzzz][Cart::s][Cart::xxz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xxz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xxz][3]+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::xz][3]
R[Cart::xzzz][Cart::s][Cart::xzz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::xzz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::xzz][3]+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::zz][3]
R[Cart::xzzz][Cart::s][Cart::zzz][2]+=pma0*R[Cart::zzz][Cart::s][Cart::zzz][2]+wmp0*R[Cart::zzz][Cart::s][Cart::zzz][3]
R[Cart::zzzz][Cart::s][Cart::yyy][2]+=pma2*R[Cart::zzz][Cart::s][Cart::yyy][2]+wmp2*R[Cart::zzz][Cart::s][Cart::yyy][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::yyy][2]-gfak*R[Cart::zz][Cart::s][Cart::yyy][3])
R[Cart::zzzz][Cart::s][Cart::xyy][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xyy][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xyy][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xyy][2]-gfak*R[Cart::zz][Cart::s][Cart::xyy][3])
R[Cart::zzzz][Cart::s][Cart::yyz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::yyz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::yyz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::yyz][2]-gfak*R[Cart::zz][Cart::s][Cart::yyz][3])+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::yy][3]
R[Cart::zzzz][Cart::s][Cart::xxy][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xxy][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xxy][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xxy][2]-gfak*R[Cart::zz][Cart::s][Cart::xxy][3])
R[Cart::zzzz][Cart::s][Cart::xyz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xyz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xyz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xyz][2]-gfak*R[Cart::zz][Cart::s][Cart::xyz][3])+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::xy][3]
R[Cart::zzzz][Cart::s][Cart::yzz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::yzz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::yzz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::yzz][2]-gfak*R[Cart::zz][Cart::s][Cart::yzz][3])+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::yz][3]
R[Cart::zzzz][Cart::s][Cart::xxx][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xxx][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xxx][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xxx][2]-gfak*R[Cart::zz][Cart::s][Cart::xxx][3])
R[Cart::zzzz][Cart::s][Cart::xxz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xxz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xxz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xxz][2]-gfak*R[Cart::zz][Cart::s][Cart::xxz][3])+0.5/_decay*1*R[Cart::zzz][Cart::s][Cart::xx][3]
R[Cart::zzzz][Cart::s][Cart::xzz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::xzz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::xzz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::xzz][2]-gfak*R[Cart::zz][Cart::s][Cart::xzz][3])+0.5/_decay*2*R[Cart::zzz][Cart::s][Cart::xz][3]
R[Cart::zzzz][Cart::s][Cart::zzz][2]+=pma2*R[Cart::zzz][Cart::s][Cart::zzz][2]+wmp2*R[Cart::zzz][Cart::s][Cart::zzz][3]+2*rzeta*(R[Cart::zz][Cart::s][Cart::zzz][2]-gfak*R[Cart::zz][Cart::s][Cart::zzz][3])+0.5/_decay*3*R[Cart::zzz][Cart::s][Cart::zz][3]
}}
//------------------------------------------------------

//Integral h - s - s - m1
if (_mmax >5 ){
if (_lmax_alpha>4){

R[Cart::yyyyy][Cart::s][Cart::s][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::s][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::s][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::s][1]-gfak*R[Cart::yyy][Cart::s][Cart::s][2])
R[Cart::xyyyy][Cart::s][Cart::s][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::s][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::s][2]
R[Cart::yyyyz][Cart::s][Cart::s][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::s][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::s][2]
R[Cart::xxyyy][Cart::s][Cart::s][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::s][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::s][2]
R[Cart::xyyyz][Cart::s][Cart::s][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::s][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::s][2]
R[Cart::yyyzz][Cart::s][Cart::s][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::s][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::s][2]
R[Cart::xxxyy][Cart::s][Cart::s][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::s][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::s][2]
R[Cart::xxyyz][Cart::s][Cart::s][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::s][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::s][2]
R[Cart::xyyzz][Cart::s][Cart::s][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::s][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::s][2]
R[Cart::yyzzz][Cart::s][Cart::s][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::s][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::s][2]
R[Cart::xxxxy][Cart::s][Cart::s][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::s][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::s][2]
R[Cart::xxxyz][Cart::s][Cart::s][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::s][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::s][2]
R[Cart::xxyzz][Cart::s][Cart::s][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::s][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::s][2]
R[Cart::xyzzz][Cart::s][Cart::s][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::s][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::s][2]
R[Cart::yzzzz][Cart::s][Cart::s][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::s][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::s][2]
R[Cart::xxxxx][Cart::s][Cart::s][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::s][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::s][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::s][1]-gfak*R[Cart::xxx][Cart::s][Cart::s][2])
R[Cart::xxxxz][Cart::s][Cart::s][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::s][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::s][2]
R[Cart::xxxzz][Cart::s][Cart::s][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::s][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::s][2]
R[Cart::xxzzz][Cart::s][Cart::s][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::s][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::s][2]
R[Cart::xzzzz][Cart::s][Cart::s][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::s][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::s][2]
R[Cart::zzzzz][Cart::s][Cart::s][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::s][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::s][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::s][1]-gfak*R[Cart::zzz][Cart::s][Cart::s][2])
}}
//------------------------------------------------------

//Integral h - s - p - m1
if (_mmax >6 ){
if (_lmax_alpha>4 && _lmax_gamma>0){

R[Cart::yyyyy][Cart::s][Cart::y][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::y][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::y][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::y][1]-gfak*R[Cart::yyy][Cart::s][Cart::y][2])+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::s][2]
R[Cart::yyyyy][Cart::s][Cart::x][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::x][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::x][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::x][1]-gfak*R[Cart::yyy][Cart::s][Cart::x][2])
R[Cart::yyyyy][Cart::s][Cart::z][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::z][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::z][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::z][1]-gfak*R[Cart::yyy][Cart::s][Cart::z][2])
R[Cart::xyyyy][Cart::s][Cart::y][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::y][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::y][2]
R[Cart::xyyyy][Cart::s][Cart::x][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::x][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::s][2]
R[Cart::xyyyy][Cart::s][Cart::z][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::z][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::z][2]
R[Cart::yyyyz][Cart::s][Cart::y][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::y][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::y][2]
R[Cart::yyyyz][Cart::s][Cart::x][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::x][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::x][2]
R[Cart::yyyyz][Cart::s][Cart::z][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::z][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::z][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::s][2]
R[Cart::xxyyy][Cart::s][Cart::y][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::y][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::y][2]
R[Cart::xxyyy][Cart::s][Cart::x][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::x][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::xyyy][Cart::s][Cart::s][2]
R[Cart::xxyyy][Cart::s][Cart::z][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::z][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::z][2]
R[Cart::xyyyz][Cart::s][Cart::y][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::y][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::y][2]
R[Cart::xyyyz][Cart::s][Cart::x][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::x][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::s][2]
R[Cart::xyyyz][Cart::s][Cart::z][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::z][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::z][2]
R[Cart::yyyzz][Cart::s][Cart::y][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::y][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::y][2]
R[Cart::yyyzz][Cart::s][Cart::x][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::x][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::x][2]
R[Cart::yyyzz][Cart::s][Cart::z][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::z][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::z][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::s][2]
R[Cart::xxxyy][Cart::s][Cart::y][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::y][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::y][2]+0.5/_decay*1*R[Cart::xxxy][Cart::s][Cart::s][2]
R[Cart::xxxyy][Cart::s][Cart::x][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::x][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::x][2]
R[Cart::xxxyy][Cart::s][Cart::z][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::z][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::z][2]
R[Cart::xxyyz][Cart::s][Cart::y][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::y][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::y][2]
R[Cart::xxyyz][Cart::s][Cart::x][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::x][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::x][2]
R[Cart::xxyyz][Cart::s][Cart::z][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::z][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::z][2]+0.5/_decay*1*R[Cart::xxyy][Cart::s][Cart::s][2]
R[Cart::xyyzz][Cart::s][Cart::y][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::y][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::y][2]
R[Cart::xyyzz][Cart::s][Cart::x][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::x][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::yyzz][Cart::s][Cart::s][2]
R[Cart::xyyzz][Cart::s][Cart::z][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::z][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::z][2]
R[Cart::yyzzz][Cart::s][Cart::y][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::y][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::y][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::s][2]
R[Cart::yyzzz][Cart::s][Cart::x][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::x][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::x][2]
R[Cart::yyzzz][Cart::s][Cart::z][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::z][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::z][2]
R[Cart::xxxxy][Cart::s][Cart::y][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::y][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::y][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::s][2]
R[Cart::xxxxy][Cart::s][Cart::x][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::x][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::x][2]
R[Cart::xxxxy][Cart::s][Cart::z][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::z][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::z][2]
R[Cart::xxxyz][Cart::s][Cart::y][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::y][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::y][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::s][2]
R[Cart::xxxyz][Cart::s][Cart::x][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::x][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::x][2]
R[Cart::xxxyz][Cart::s][Cart::z][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::z][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::z][2]
R[Cart::xxyzz][Cart::s][Cart::y][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::y][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::y][2]+0.5/_decay*1*R[Cart::xxzz][Cart::s][Cart::s][2]
R[Cart::xxyzz][Cart::s][Cart::x][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::x][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::x][2]
R[Cart::xxyzz][Cart::s][Cart::z][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::z][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::z][2]
R[Cart::xyzzz][Cart::s][Cart::y][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::y][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::y][2]
R[Cart::xyzzz][Cart::s][Cart::x][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::x][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::s][2]
R[Cart::xyzzz][Cart::s][Cart::z][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::z][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::z][2]
R[Cart::yzzzz][Cart::s][Cart::y][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::y][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::y][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::s][2]
R[Cart::yzzzz][Cart::s][Cart::x][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::x][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::x][2]
R[Cart::yzzzz][Cart::s][Cart::z][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::z][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::z][2]
R[Cart::xxxxx][Cart::s][Cart::y][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::y][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::y][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::y][1]-gfak*R[Cart::xxx][Cart::s][Cart::y][2])
R[Cart::xxxxx][Cart::s][Cart::x][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::x][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::x][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::x][1]-gfak*R[Cart::xxx][Cart::s][Cart::x][2])+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::s][2]
R[Cart::xxxxx][Cart::s][Cart::z][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::z][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::z][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::z][1]-gfak*R[Cart::xxx][Cart::s][Cart::z][2])
R[Cart::xxxxz][Cart::s][Cart::y][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::y][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::y][2]
R[Cart::xxxxz][Cart::s][Cart::x][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::x][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::x][2]
R[Cart::xxxxz][Cart::s][Cart::z][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::z][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::z][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::s][2]
R[Cart::xxxzz][Cart::s][Cart::y][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::y][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::y][2]
R[Cart::xxxzz][Cart::s][Cart::x][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::x][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::x][2]
R[Cart::xxxzz][Cart::s][Cart::z][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::z][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::z][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::s][2]
R[Cart::xxzzz][Cart::s][Cart::y][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::y][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::y][2]
R[Cart::xxzzz][Cart::s][Cart::x][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::x][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::xzzz][Cart::s][Cart::s][2]
R[Cart::xxzzz][Cart::s][Cart::z][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::z][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::z][2]
R[Cart::xzzzz][Cart::s][Cart::y][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::y][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::y][2]
R[Cart::xzzzz][Cart::s][Cart::x][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::x][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::x][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::s][2]
R[Cart::xzzzz][Cart::s][Cart::z][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::z][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::z][2]
R[Cart::zzzzz][Cart::s][Cart::y][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::y][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::y][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::y][1]-gfak*R[Cart::zzz][Cart::s][Cart::y][2])
R[Cart::zzzzz][Cart::s][Cart::x][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::x][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::x][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::x][1]-gfak*R[Cart::zzz][Cart::s][Cart::x][2])
R[Cart::zzzzz][Cart::s][Cart::z][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::z][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::z][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::z][1]-gfak*R[Cart::zzz][Cart::s][Cart::z][2])+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::s][2]
}}
//------------------------------------------------------

//Integral h - s - d - m1
if (_mmax >7 ){
if (_lmax_alpha>4 && _lmax_gamma>1){

R[Cart::yyyyy][Cart::s][Cart::yy][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::yy][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::yy][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::yy][1]-gfak*R[Cart::yyy][Cart::s][Cart::yy][2])+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::y][2]
R[Cart::yyyyy][Cart::s][Cart::xy][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xy][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xy][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xy][1]-gfak*R[Cart::yyy][Cart::s][Cart::xy][2])+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::x][2]
R[Cart::yyyyy][Cart::s][Cart::yz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::yz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::yz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::yz][1]-gfak*R[Cart::yyy][Cart::s][Cart::yz][2])+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::z][2]
R[Cart::yyyyy][Cart::s][Cart::xx][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xx][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xx][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xx][1]-gfak*R[Cart::yyy][Cart::s][Cart::xx][2])
R[Cart::yyyyy][Cart::s][Cart::xz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xz][1]-gfak*R[Cart::yyy][Cart::s][Cart::xz][2])
R[Cart::yyyyy][Cart::s][Cart::zz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::zz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::zz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::zz][1]-gfak*R[Cart::yyy][Cart::s][Cart::zz][2])
R[Cart::xyyyy][Cart::s][Cart::yy][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::yy][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::yy][2]
R[Cart::xyyyy][Cart::s][Cart::xy][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xy][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::y][2]
R[Cart::xyyyy][Cart::s][Cart::yz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::yz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::yz][2]
R[Cart::xyyyy][Cart::s][Cart::xx][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xx][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::x][2]
R[Cart::xyyyy][Cart::s][Cart::xz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::z][2]
R[Cart::xyyyy][Cart::s][Cart::zz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::zz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::zz][2]
R[Cart::yyyyz][Cart::s][Cart::yy][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::yy][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::yy][2]
R[Cart::yyyyz][Cart::s][Cart::xy][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xy][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xy][2]
R[Cart::yyyyz][Cart::s][Cart::yz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::yz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::y][2]
R[Cart::yyyyz][Cart::s][Cart::xx][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xx][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xx][2]
R[Cart::yyyyz][Cart::s][Cart::xz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::x][2]
R[Cart::yyyyz][Cart::s][Cart::zz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::zz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::zz][2]+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::z][2]
R[Cart::xxyyy][Cart::s][Cart::yy][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::yy][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::yy][2]
R[Cart::xxyyy][Cart::s][Cart::xy][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xy][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::xyyy][Cart::s][Cart::y][2]
R[Cart::xxyyy][Cart::s][Cart::yz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::yz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::yz][2]
R[Cart::xxyyy][Cart::s][Cart::xx][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xx][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::xyyy][Cart::s][Cart::x][2]
R[Cart::xxyyy][Cart::s][Cart::xz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::xyyy][Cart::s][Cart::z][2]
R[Cart::xxyyy][Cart::s][Cart::zz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::zz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::zz][2]
R[Cart::xyyyz][Cart::s][Cart::yy][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::yy][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::yy][2]
R[Cart::xyyyz][Cart::s][Cart::xy][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xy][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::y][2]
R[Cart::xyyyz][Cart::s][Cart::yz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::yz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::yz][2]
R[Cart::xyyyz][Cart::s][Cart::xx][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xx][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::yyyz][Cart::s][Cart::x][2]
R[Cart::xyyyz][Cart::s][Cart::xz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::z][2]
R[Cart::xyyyz][Cart::s][Cart::zz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::zz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::zz][2]
R[Cart::yyyzz][Cart::s][Cart::yy][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::yy][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::yy][2]
R[Cart::yyyzz][Cart::s][Cart::xy][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xy][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xy][2]
R[Cart::yyyzz][Cart::s][Cart::yz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::yz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::y][2]
R[Cart::yyyzz][Cart::s][Cart::xx][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xx][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xx][2]
R[Cart::yyyzz][Cart::s][Cart::xz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::x][2]
R[Cart::yyyzz][Cart::s][Cart::zz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::zz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::zz][2]+0.5/_decay*2*R[Cart::yyyz][Cart::s][Cart::z][2]
R[Cart::xxxyy][Cart::s][Cart::yy][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::yy][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::yy][2]+0.5/_decay*2*R[Cart::xxxy][Cart::s][Cart::y][2]
R[Cart::xxxyy][Cart::s][Cart::xy][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xy][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::xxxy][Cart::s][Cart::x][2]
R[Cart::xxxyy][Cart::s][Cart::yz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::yz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxxy][Cart::s][Cart::z][2]
R[Cart::xxxyy][Cart::s][Cart::xx][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xx][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xx][2]
R[Cart::xxxyy][Cart::s][Cart::xz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xz][2]
R[Cart::xxxyy][Cart::s][Cart::zz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::zz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::zz][2]
R[Cart::xxyyz][Cart::s][Cart::yy][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::yy][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::yy][2]
R[Cart::xxyyz][Cart::s][Cart::xy][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xy][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xy][2]
R[Cart::xxyyz][Cart::s][Cart::yz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::yz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxyy][Cart::s][Cart::y][2]
R[Cart::xxyyz][Cart::s][Cart::xx][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xx][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xx][2]
R[Cart::xxyyz][Cart::s][Cart::xz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::xxyy][Cart::s][Cart::x][2]
R[Cart::xxyyz][Cart::s][Cart::zz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::zz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::zz][2]+0.5/_decay*2*R[Cart::xxyy][Cart::s][Cart::z][2]
R[Cart::xyyzz][Cart::s][Cart::yy][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::yy][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::yy][2]
R[Cart::xyyzz][Cart::s][Cart::xy][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xy][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::yyzz][Cart::s][Cart::y][2]
R[Cart::xyyzz][Cart::s][Cart::yz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::yz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::yz][2]
R[Cart::xyyzz][Cart::s][Cart::xx][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xx][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::yyzz][Cart::s][Cart::x][2]
R[Cart::xyyzz][Cart::s][Cart::xz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::yyzz][Cart::s][Cart::z][2]
R[Cart::xyyzz][Cart::s][Cart::zz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::zz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::zz][2]
R[Cart::yyzzz][Cart::s][Cart::yy][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::yy][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::yy][2]+0.5/_decay*2*R[Cart::yzzz][Cart::s][Cart::y][2]
R[Cart::yyzzz][Cart::s][Cart::xy][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xy][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::x][2]
R[Cart::yyzzz][Cart::s][Cart::yz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::yz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::z][2]
R[Cart::yyzzz][Cart::s][Cart::xx][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xx][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xx][2]
R[Cart::yyzzz][Cart::s][Cart::xz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xz][2]
R[Cart::yyzzz][Cart::s][Cart::zz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::zz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::zz][2]
R[Cart::xxxxy][Cart::s][Cart::yy][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::yy][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::yy][2]+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::y][2]
R[Cart::xxxxy][Cart::s][Cart::xy][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xy][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::x][2]
R[Cart::xxxxy][Cart::s][Cart::yz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::yz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::z][2]
R[Cart::xxxxy][Cart::s][Cart::xx][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xx][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xx][2]
R[Cart::xxxxy][Cart::s][Cart::xz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xz][2]
R[Cart::xxxxy][Cart::s][Cart::zz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::zz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::zz][2]
R[Cart::xxxyz][Cart::s][Cart::yy][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::yy][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::yy][2]+0.5/_decay*2*R[Cart::xxxz][Cart::s][Cart::y][2]
R[Cart::xxxyz][Cart::s][Cart::xy][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xy][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::x][2]
R[Cart::xxxyz][Cart::s][Cart::yz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::yz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::z][2]
R[Cart::xxxyz][Cart::s][Cart::xx][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xx][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xx][2]
R[Cart::xxxyz][Cart::s][Cart::xz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xz][2]
R[Cart::xxxyz][Cart::s][Cart::zz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::zz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::zz][2]
R[Cart::xxyzz][Cart::s][Cart::yy][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::yy][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::yy][2]+0.5/_decay*2*R[Cart::xxzz][Cart::s][Cart::y][2]
R[Cart::xxyzz][Cart::s][Cart::xy][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xy][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::xxzz][Cart::s][Cart::x][2]
R[Cart::xxyzz][Cart::s][Cart::yz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::yz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxzz][Cart::s][Cart::z][2]
R[Cart::xxyzz][Cart::s][Cart::xx][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xx][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xx][2]
R[Cart::xxyzz][Cart::s][Cart::xz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xz][2]
R[Cart::xxyzz][Cart::s][Cart::zz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::zz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::zz][2]
R[Cart::xyzzz][Cart::s][Cart::yy][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::yy][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::yy][2]
R[Cart::xyzzz][Cart::s][Cart::xy][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xy][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::y][2]
R[Cart::xyzzz][Cart::s][Cart::yz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::yz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::yz][2]
R[Cart::xyzzz][Cart::s][Cart::xx][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xx][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::yzzz][Cart::s][Cart::x][2]
R[Cart::xyzzz][Cart::s][Cart::xz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::z][2]
R[Cart::xyzzz][Cart::s][Cart::zz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::zz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::zz][2]
R[Cart::yzzzz][Cart::s][Cart::yy][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::yy][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::yy][2]+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::y][2]
R[Cart::yzzzz][Cart::s][Cart::xy][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xy][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::x][2]
R[Cart::yzzzz][Cart::s][Cart::yz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::yz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::z][2]
R[Cart::yzzzz][Cart::s][Cart::xx][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xx][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xx][2]
R[Cart::yzzzz][Cart::s][Cart::xz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xz][2]
R[Cart::yzzzz][Cart::s][Cart::zz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::zz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::zz][2]
R[Cart::xxxxx][Cart::s][Cart::yy][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::yy][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::yy][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::yy][1]-gfak*R[Cart::xxx][Cart::s][Cart::yy][2])
R[Cart::xxxxx][Cart::s][Cart::xy][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xy][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xy][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xy][1]-gfak*R[Cart::xxx][Cart::s][Cart::xy][2])+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::y][2]
R[Cart::xxxxx][Cart::s][Cart::yz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::yz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::yz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::yz][1]-gfak*R[Cart::xxx][Cart::s][Cart::yz][2])
R[Cart::xxxxx][Cart::s][Cart::xx][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xx][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xx][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xx][1]-gfak*R[Cart::xxx][Cart::s][Cart::xx][2])+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::x][2]
R[Cart::xxxxx][Cart::s][Cart::xz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xz][1]-gfak*R[Cart::xxx][Cart::s][Cart::xz][2])+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::z][2]
R[Cart::xxxxx][Cart::s][Cart::zz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::zz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::zz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::zz][1]-gfak*R[Cart::xxx][Cart::s][Cart::zz][2])
R[Cart::xxxxz][Cart::s][Cart::yy][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::yy][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::yy][2]
R[Cart::xxxxz][Cart::s][Cart::xy][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xy][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xy][2]
R[Cart::xxxxz][Cart::s][Cart::yz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::yz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::y][2]
R[Cart::xxxxz][Cart::s][Cart::xx][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xx][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xx][2]
R[Cart::xxxxz][Cart::s][Cart::xz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::x][2]
R[Cart::xxxxz][Cart::s][Cart::zz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::zz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::zz][2]+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::z][2]
R[Cart::xxxzz][Cart::s][Cart::yy][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::yy][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::yy][2]
R[Cart::xxxzz][Cart::s][Cart::xy][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xy][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xy][2]
R[Cart::xxxzz][Cart::s][Cart::yz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::yz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::yz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::y][2]
R[Cart::xxxzz][Cart::s][Cart::xx][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xx][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xx][2]
R[Cart::xxxzz][Cart::s][Cart::xz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::x][2]
R[Cart::xxxzz][Cart::s][Cart::zz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::zz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::zz][2]+0.5/_decay*2*R[Cart::xxxz][Cart::s][Cart::z][2]
R[Cart::xxzzz][Cart::s][Cart::yy][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::yy][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::yy][2]
R[Cart::xxzzz][Cart::s][Cart::xy][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xy][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::xzzz][Cart::s][Cart::y][2]
R[Cart::xxzzz][Cart::s][Cart::yz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::yz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::yz][2]
R[Cart::xxzzz][Cart::s][Cart::xx][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xx][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::xzzz][Cart::s][Cart::x][2]
R[Cart::xxzzz][Cart::s][Cart::xz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::xzzz][Cart::s][Cart::z][2]
R[Cart::xxzzz][Cart::s][Cart::zz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::zz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::zz][2]
R[Cart::xzzzz][Cart::s][Cart::yy][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::yy][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::yy][2]
R[Cart::xzzzz][Cart::s][Cart::xy][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xy][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xy][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::y][2]
R[Cart::xzzzz][Cart::s][Cart::yz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::yz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::yz][2]
R[Cart::xzzzz][Cart::s][Cart::xx][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xx][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xx][2]+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::x][2]
R[Cart::xzzzz][Cart::s][Cart::xz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xz][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::z][2]
R[Cart::xzzzz][Cart::s][Cart::zz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::zz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::zz][2]
R[Cart::zzzzz][Cart::s][Cart::yy][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::yy][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::yy][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::yy][1]-gfak*R[Cart::zzz][Cart::s][Cart::yy][2])
R[Cart::zzzzz][Cart::s][Cart::xy][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xy][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xy][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xy][1]-gfak*R[Cart::zzz][Cart::s][Cart::xy][2])
R[Cart::zzzzz][Cart::s][Cart::yz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::yz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::yz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::yz][1]-gfak*R[Cart::zzz][Cart::s][Cart::yz][2])+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::y][2]
R[Cart::zzzzz][Cart::s][Cart::xx][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xx][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xx][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xx][1]-gfak*R[Cart::zzz][Cart::s][Cart::xx][2])
R[Cart::zzzzz][Cart::s][Cart::xz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xz][1]-gfak*R[Cart::zzz][Cart::s][Cart::xz][2])+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::x][2]
R[Cart::zzzzz][Cart::s][Cart::zz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::zz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::zz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::zz][1]-gfak*R[Cart::zzz][Cart::s][Cart::zz][2])+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::z][2]
}}
//------------------------------------------------------

//Integral h - s - f - m1
if (_mmax >8 ){
if (_lmax_alpha>4 && _lmax_gamma>2){

R[Cart::yyyyy][Cart::s][Cart::yyy][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::yyy][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::yyy][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::yyy][1]-gfak*R[Cart::yyy][Cart::s][Cart::yyy][2])+0.5/_decay*3*R[Cart::yyyy][Cart::s][Cart::yy][2]
R[Cart::yyyyy][Cart::s][Cart::xyy][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xyy][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xyy][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xyy][1]-gfak*R[Cart::yyy][Cart::s][Cart::xyy][2])+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::xy][2]
R[Cart::yyyyy][Cart::s][Cart::yyz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::yyz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::yyz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::yyz][1]-gfak*R[Cart::yyy][Cart::s][Cart::yyz][2])+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::yz][2]
R[Cart::yyyyy][Cart::s][Cart::xxy][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xxy][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xxy][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xxy][1]-gfak*R[Cart::yyy][Cart::s][Cart::xxy][2])+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::xx][2]
R[Cart::yyyyy][Cart::s][Cart::xyz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xyz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xyz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xyz][1]-gfak*R[Cart::yyy][Cart::s][Cart::xyz][2])+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::xz][2]
R[Cart::yyyyy][Cart::s][Cart::yzz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::yzz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::yzz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::yzz][1]-gfak*R[Cart::yyy][Cart::s][Cart::yzz][2])+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::zz][2]
R[Cart::yyyyy][Cart::s][Cart::xxx][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xxx][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xxx][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xxx][1]-gfak*R[Cart::yyy][Cart::s][Cart::xxx][2])
R[Cart::yyyyy][Cart::s][Cart::xxz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xxz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xxz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xxz][1]-gfak*R[Cart::yyy][Cart::s][Cart::xxz][2])
R[Cart::yyyyy][Cart::s][Cart::xzz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::xzz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::xzz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::xzz][1]-gfak*R[Cart::yyy][Cart::s][Cart::xzz][2])
R[Cart::yyyyy][Cart::s][Cart::zzz][1]+=pma1*R[Cart::yyyy][Cart::s][Cart::zzz][1]+wmp1*R[Cart::yyyy][Cart::s][Cart::zzz][2]+3*rzeta*(R[Cart::yyy][Cart::s][Cart::zzz][1]-gfak*R[Cart::yyy][Cart::s][Cart::zzz][2])
R[Cart::xyyyy][Cart::s][Cart::yyy][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::yyy][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::yyy][2]
R[Cart::xyyyy][Cart::s][Cart::xyy][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xyy][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::yy][2]
R[Cart::xyyyy][Cart::s][Cart::yyz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::yyz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::yyz][2]
R[Cart::xyyyy][Cart::s][Cart::xxy][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xxy][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::xy][2]
R[Cart::xyyyy][Cart::s][Cart::xyz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xyz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::yz][2]
R[Cart::xyyyy][Cart::s][Cart::yzz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::yzz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::yzz][2]
R[Cart::xyyyy][Cart::s][Cart::xxx][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xxx][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::yyyy][Cart::s][Cart::xx][2]
R[Cart::xyyyy][Cart::s][Cart::xxz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xxz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::xz][2]
R[Cart::xyyyy][Cart::s][Cart::xzz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::xzz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::zz][2]
R[Cart::xyyyy][Cart::s][Cart::zzz][1]+=pma0*R[Cart::yyyy][Cart::s][Cart::zzz][1]+wmp0*R[Cart::yyyy][Cart::s][Cart::zzz][2]
R[Cart::yyyyz][Cart::s][Cart::yyy][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::yyy][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::yyy][2]
R[Cart::yyyyz][Cart::s][Cart::xyy][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xyy][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xyy][2]
R[Cart::yyyyz][Cart::s][Cart::yyz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::yyz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::yyz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::yy][2]
R[Cart::yyyyz][Cart::s][Cart::xxy][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xxy][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xxy][2]
R[Cart::yyyyz][Cart::s][Cart::xyz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xyz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::xy][2]
R[Cart::yyyyz][Cart::s][Cart::yzz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::yzz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::yzz][2]+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::yz][2]
R[Cart::yyyyz][Cart::s][Cart::xxx][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xxx][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xxx][2]
R[Cart::yyyyz][Cart::s][Cart::xxz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xxz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xxz][2]+0.5/_decay*1*R[Cart::yyyy][Cart::s][Cart::xx][2]
R[Cart::yyyyz][Cart::s][Cart::xzz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::xzz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::xzz][2]+0.5/_decay*2*R[Cart::yyyy][Cart::s][Cart::xz][2]
R[Cart::yyyyz][Cart::s][Cart::zzz][1]+=pma2*R[Cart::yyyy][Cart::s][Cart::zzz][1]+wmp2*R[Cart::yyyy][Cart::s][Cart::zzz][2]+0.5/_decay*3*R[Cart::yyyy][Cart::s][Cart::zz][2]
R[Cart::xxyyy][Cart::s][Cart::yyy][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::yyy][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::yyy][2]
R[Cart::xxyyy][Cart::s][Cart::xyy][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xyy][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::xyyy][Cart::s][Cart::yy][2]
R[Cart::xxyyy][Cart::s][Cart::yyz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::yyz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::yyz][2]
R[Cart::xxyyy][Cart::s][Cart::xxy][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xxy][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::xyyy][Cart::s][Cart::xy][2]
R[Cart::xxyyy][Cart::s][Cart::xyz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xyz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xyyy][Cart::s][Cart::yz][2]
R[Cart::xxyyy][Cart::s][Cart::yzz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::yzz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::yzz][2]
R[Cart::xxyyy][Cart::s][Cart::xxx][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xxx][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::xyyy][Cart::s][Cart::xx][2]
R[Cart::xxyyy][Cart::s][Cart::xxz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xxz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::xyyy][Cart::s][Cart::xz][2]
R[Cart::xxyyy][Cart::s][Cart::xzz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::xzz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::xyyy][Cart::s][Cart::zz][2]
R[Cart::xxyyy][Cart::s][Cart::zzz][1]+=pma0*R[Cart::xyyy][Cart::s][Cart::zzz][1]+wmp0*R[Cart::xyyy][Cart::s][Cart::zzz][2]
R[Cart::xyyyz][Cart::s][Cart::yyy][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::yyy][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::yyy][2]
R[Cart::xyyyz][Cart::s][Cart::xyy][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xyy][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::yy][2]
R[Cart::xyyyz][Cart::s][Cart::yyz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::yyz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::yyz][2]
R[Cart::xyyyz][Cart::s][Cart::xxy][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xxy][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::yyyz][Cart::s][Cart::xy][2]
R[Cart::xyyyz][Cart::s][Cart::xyz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xyz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::yz][2]
R[Cart::xyyyz][Cart::s][Cart::yzz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::yzz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::yzz][2]
R[Cart::xyyyz][Cart::s][Cart::xxx][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xxx][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::yyyz][Cart::s][Cart::xx][2]
R[Cart::xyyyz][Cart::s][Cart::xxz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xxz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::yyyz][Cart::s][Cart::xz][2]
R[Cart::xyyyz][Cart::s][Cart::xzz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::xzz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::zz][2]
R[Cart::xyyyz][Cart::s][Cart::zzz][1]+=pma0*R[Cart::yyyz][Cart::s][Cart::zzz][1]+wmp0*R[Cart::yyyz][Cart::s][Cart::zzz][2]
R[Cart::yyyzz][Cart::s][Cart::yyy][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::yyy][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::yyy][2]
R[Cart::yyyzz][Cart::s][Cart::xyy][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xyy][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xyy][2]
R[Cart::yyyzz][Cart::s][Cart::yyz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::yyz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::yyz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::yy][2]
R[Cart::yyyzz][Cart::s][Cart::xxy][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xxy][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xxy][2]
R[Cart::yyyzz][Cart::s][Cart::xyz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xyz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::xy][2]
R[Cart::yyyzz][Cart::s][Cart::yzz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::yzz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::yzz][2]+0.5/_decay*2*R[Cart::yyyz][Cart::s][Cart::yz][2]
R[Cart::yyyzz][Cart::s][Cart::xxx][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xxx][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xxx][2]
R[Cart::yyyzz][Cart::s][Cart::xxz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xxz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xxz][2]+0.5/_decay*1*R[Cart::yyyz][Cart::s][Cart::xx][2]
R[Cart::yyyzz][Cart::s][Cart::xzz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::xzz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::xzz][2]+0.5/_decay*2*R[Cart::yyyz][Cart::s][Cart::xz][2]
R[Cart::yyyzz][Cart::s][Cart::zzz][1]+=pma2*R[Cart::yyyz][Cart::s][Cart::zzz][1]+wmp2*R[Cart::yyyz][Cart::s][Cart::zzz][2]+0.5/_decay*3*R[Cart::yyyz][Cart::s][Cart::zz][2]
R[Cart::xxxyy][Cart::s][Cart::yyy][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::yyy][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::yyy][2]+0.5/_decay*3*R[Cart::xxxy][Cart::s][Cart::yy][2]
R[Cart::xxxyy][Cart::s][Cart::xyy][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xyy][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xyy][2]+0.5/_decay*2*R[Cart::xxxy][Cart::s][Cart::xy][2]
R[Cart::xxxyy][Cart::s][Cart::yyz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::yyz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::yyz][2]+0.5/_decay*2*R[Cart::xxxy][Cart::s][Cart::yz][2]
R[Cart::xxxyy][Cart::s][Cart::xxy][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xxy][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xxy][2]+0.5/_decay*1*R[Cart::xxxy][Cart::s][Cart::xx][2]
R[Cart::xxxyy][Cart::s][Cart::xyz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xyz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxxy][Cart::s][Cart::xz][2]
R[Cart::xxxyy][Cart::s][Cart::yzz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::yzz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::yzz][2]+0.5/_decay*1*R[Cart::xxxy][Cart::s][Cart::zz][2]
R[Cart::xxxyy][Cart::s][Cart::xxx][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xxx][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xxx][2]
R[Cart::xxxyy][Cart::s][Cart::xxz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xxz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xxz][2]
R[Cart::xxxyy][Cart::s][Cart::xzz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::xzz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::xzz][2]
R[Cart::xxxyy][Cart::s][Cart::zzz][1]+=pma1*R[Cart::xxxy][Cart::s][Cart::zzz][1]+wmp1*R[Cart::xxxy][Cart::s][Cart::zzz][2]
R[Cart::xxyyz][Cart::s][Cart::yyy][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::yyy][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::yyy][2]
R[Cart::xxyyz][Cart::s][Cart::xyy][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xyy][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xyy][2]
R[Cart::xxyyz][Cart::s][Cart::yyz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::yyz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::yyz][2]+0.5/_decay*1*R[Cart::xxyy][Cart::s][Cart::yy][2]
R[Cart::xxyyz][Cart::s][Cart::xxy][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xxy][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xxy][2]
R[Cart::xxyyz][Cart::s][Cart::xyz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xyz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxyy][Cart::s][Cart::xy][2]
R[Cart::xxyyz][Cart::s][Cart::yzz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::yzz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::yzz][2]+0.5/_decay*2*R[Cart::xxyy][Cart::s][Cart::yz][2]
R[Cart::xxyyz][Cart::s][Cart::xxx][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xxx][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xxx][2]
R[Cart::xxyyz][Cart::s][Cart::xxz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xxz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xxz][2]+0.5/_decay*1*R[Cart::xxyy][Cart::s][Cart::xx][2]
R[Cart::xxyyz][Cart::s][Cart::xzz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::xzz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::xzz][2]+0.5/_decay*2*R[Cart::xxyy][Cart::s][Cart::xz][2]
R[Cart::xxyyz][Cart::s][Cart::zzz][1]+=pma2*R[Cart::xxyy][Cart::s][Cart::zzz][1]+wmp2*R[Cart::xxyy][Cart::s][Cart::zzz][2]+0.5/_decay*3*R[Cart::xxyy][Cart::s][Cart::zz][2]
R[Cart::xyyzz][Cart::s][Cart::yyy][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::yyy][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::yyy][2]
R[Cart::xyyzz][Cart::s][Cart::xyy][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xyy][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::yyzz][Cart::s][Cart::yy][2]
R[Cart::xyyzz][Cart::s][Cart::yyz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::yyz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::yyz][2]
R[Cart::xyyzz][Cart::s][Cart::xxy][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xxy][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::yyzz][Cart::s][Cart::xy][2]
R[Cart::xyyzz][Cart::s][Cart::xyz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xyz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yyzz][Cart::s][Cart::yz][2]
R[Cart::xyyzz][Cart::s][Cart::yzz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::yzz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::yzz][2]
R[Cart::xyyzz][Cart::s][Cart::xxx][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xxx][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::yyzz][Cart::s][Cart::xx][2]
R[Cart::xyyzz][Cart::s][Cart::xxz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xxz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::yyzz][Cart::s][Cart::xz][2]
R[Cart::xyyzz][Cart::s][Cart::xzz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::xzz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::yyzz][Cart::s][Cart::zz][2]
R[Cart::xyyzz][Cart::s][Cart::zzz][1]+=pma0*R[Cart::yyzz][Cart::s][Cart::zzz][1]+wmp0*R[Cart::yyzz][Cart::s][Cart::zzz][2]
R[Cart::yyzzz][Cart::s][Cart::yyy][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::yyy][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::yyy][2]+0.5/_decay*3*R[Cart::yzzz][Cart::s][Cart::yy][2]
R[Cart::yyzzz][Cart::s][Cart::xyy][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xyy][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xyy][2]+0.5/_decay*2*R[Cart::yzzz][Cart::s][Cart::xy][2]
R[Cart::yyzzz][Cart::s][Cart::yyz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::yyz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::yyz][2]+0.5/_decay*2*R[Cart::yzzz][Cart::s][Cart::yz][2]
R[Cart::yyzzz][Cart::s][Cart::xxy][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xxy][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xxy][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::xx][2]
R[Cart::yyzzz][Cart::s][Cart::xyz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xyz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::xz][2]
R[Cart::yyzzz][Cart::s][Cart::yzz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::yzz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::yzz][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::zz][2]
R[Cart::yyzzz][Cart::s][Cart::xxx][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xxx][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xxx][2]
R[Cart::yyzzz][Cart::s][Cart::xxz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xxz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xxz][2]
R[Cart::yyzzz][Cart::s][Cart::xzz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::xzz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::xzz][2]
R[Cart::yyzzz][Cart::s][Cart::zzz][1]+=pma1*R[Cart::yzzz][Cart::s][Cart::zzz][1]+wmp1*R[Cart::yzzz][Cart::s][Cart::zzz][2]
R[Cart::xxxxy][Cart::s][Cart::yyy][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::yyy][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::yyy][2]+0.5/_decay*3*R[Cart::xxxx][Cart::s][Cart::yy][2]
R[Cart::xxxxy][Cart::s][Cart::xyy][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xyy][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xyy][2]+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::xy][2]
R[Cart::xxxxy][Cart::s][Cart::yyz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::yyz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::yyz][2]+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::yz][2]
R[Cart::xxxxy][Cart::s][Cart::xxy][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xxy][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xxy][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::xx][2]
R[Cart::xxxxy][Cart::s][Cart::xyz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xyz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::xz][2]
R[Cart::xxxxy][Cart::s][Cart::yzz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::yzz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::yzz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::zz][2]
R[Cart::xxxxy][Cart::s][Cart::xxx][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xxx][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xxx][2]
R[Cart::xxxxy][Cart::s][Cart::xxz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xxz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xxz][2]
R[Cart::xxxxy][Cart::s][Cart::xzz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::xzz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::xzz][2]
R[Cart::xxxxy][Cart::s][Cart::zzz][1]+=pma1*R[Cart::xxxx][Cart::s][Cart::zzz][1]+wmp1*R[Cart::xxxx][Cart::s][Cart::zzz][2]
R[Cart::xxxyz][Cart::s][Cart::yyy][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::yyy][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::yyy][2]+0.5/_decay*3*R[Cart::xxxz][Cart::s][Cart::yy][2]
R[Cart::xxxyz][Cart::s][Cart::xyy][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xyy][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xyy][2]+0.5/_decay*2*R[Cart::xxxz][Cart::s][Cart::xy][2]
R[Cart::xxxyz][Cart::s][Cart::yyz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::yyz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::yyz][2]+0.5/_decay*2*R[Cart::xxxz][Cart::s][Cart::yz][2]
R[Cart::xxxyz][Cart::s][Cart::xxy][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xxy][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xxy][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::xx][2]
R[Cart::xxxyz][Cart::s][Cart::xyz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xyz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::xz][2]
R[Cart::xxxyz][Cart::s][Cart::yzz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::yzz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::yzz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::zz][2]
R[Cart::xxxyz][Cart::s][Cart::xxx][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xxx][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xxx][2]
R[Cart::xxxyz][Cart::s][Cart::xxz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xxz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xxz][2]
R[Cart::xxxyz][Cart::s][Cart::xzz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::xzz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::xzz][2]
R[Cart::xxxyz][Cart::s][Cart::zzz][1]+=pma1*R[Cart::xxxz][Cart::s][Cart::zzz][1]+wmp1*R[Cart::xxxz][Cart::s][Cart::zzz][2]
R[Cart::xxyzz][Cart::s][Cart::yyy][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::yyy][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::yyy][2]+0.5/_decay*3*R[Cart::xxzz][Cart::s][Cart::yy][2]
R[Cart::xxyzz][Cart::s][Cart::xyy][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xyy][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xyy][2]+0.5/_decay*2*R[Cart::xxzz][Cart::s][Cart::xy][2]
R[Cart::xxyzz][Cart::s][Cart::yyz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::yyz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::yyz][2]+0.5/_decay*2*R[Cart::xxzz][Cart::s][Cart::yz][2]
R[Cart::xxyzz][Cart::s][Cart::xxy][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xxy][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xxy][2]+0.5/_decay*1*R[Cart::xxzz][Cart::s][Cart::xx][2]
R[Cart::xxyzz][Cart::s][Cart::xyz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xyz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxzz][Cart::s][Cart::xz][2]
R[Cart::xxyzz][Cart::s][Cart::yzz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::yzz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::yzz][2]+0.5/_decay*1*R[Cart::xxzz][Cart::s][Cart::zz][2]
R[Cart::xxyzz][Cart::s][Cart::xxx][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xxx][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xxx][2]
R[Cart::xxyzz][Cart::s][Cart::xxz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xxz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xxz][2]
R[Cart::xxyzz][Cart::s][Cart::xzz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::xzz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::xzz][2]
R[Cart::xxyzz][Cart::s][Cart::zzz][1]+=pma1*R[Cart::xxzz][Cart::s][Cart::zzz][1]+wmp1*R[Cart::xxzz][Cart::s][Cart::zzz][2]
R[Cart::xyzzz][Cart::s][Cart::yyy][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::yyy][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::yyy][2]
R[Cart::xyzzz][Cart::s][Cart::xyy][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xyy][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::yy][2]
R[Cart::xyzzz][Cart::s][Cart::yyz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::yyz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::yyz][2]
R[Cart::xyzzz][Cart::s][Cart::xxy][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xxy][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::yzzz][Cart::s][Cart::xy][2]
R[Cart::xyzzz][Cart::s][Cart::xyz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xyz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::yz][2]
R[Cart::xyzzz][Cart::s][Cart::yzz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::yzz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::yzz][2]
R[Cart::xyzzz][Cart::s][Cart::xxx][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xxx][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::yzzz][Cart::s][Cart::xx][2]
R[Cart::xyzzz][Cart::s][Cart::xxz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xxz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::yzzz][Cart::s][Cart::xz][2]
R[Cart::xyzzz][Cart::s][Cart::xzz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::xzz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::yzzz][Cart::s][Cart::zz][2]
R[Cart::xyzzz][Cart::s][Cart::zzz][1]+=pma0*R[Cart::yzzz][Cart::s][Cart::zzz][1]+wmp0*R[Cart::yzzz][Cart::s][Cart::zzz][2]
R[Cart::yzzzz][Cart::s][Cart::yyy][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::yyy][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::yyy][2]+0.5/_decay*3*R[Cart::zzzz][Cart::s][Cart::yy][2]
R[Cart::yzzzz][Cart::s][Cart::xyy][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xyy][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xyy][2]+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::xy][2]
R[Cart::yzzzz][Cart::s][Cart::yyz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::yyz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::yyz][2]+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::yz][2]
R[Cart::yzzzz][Cart::s][Cart::xxy][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xxy][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xxy][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::xx][2]
R[Cart::yzzzz][Cart::s][Cart::xyz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xyz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::xz][2]
R[Cart::yzzzz][Cart::s][Cart::yzz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::yzz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::yzz][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::zz][2]
R[Cart::yzzzz][Cart::s][Cart::xxx][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xxx][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xxx][2]
R[Cart::yzzzz][Cart::s][Cart::xxz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xxz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xxz][2]
R[Cart::yzzzz][Cart::s][Cart::xzz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::xzz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::xzz][2]
R[Cart::yzzzz][Cart::s][Cart::zzz][1]+=pma1*R[Cart::zzzz][Cart::s][Cart::zzz][1]+wmp1*R[Cart::zzzz][Cart::s][Cart::zzz][2]
R[Cart::xxxxx][Cart::s][Cart::yyy][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::yyy][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::yyy][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::yyy][1]-gfak*R[Cart::xxx][Cart::s][Cart::yyy][2])
R[Cart::xxxxx][Cart::s][Cart::xyy][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xyy][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xyy][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xyy][1]-gfak*R[Cart::xxx][Cart::s][Cart::xyy][2])+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::yy][2]
R[Cart::xxxxx][Cart::s][Cart::yyz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::yyz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::yyz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::yyz][1]-gfak*R[Cart::xxx][Cart::s][Cart::yyz][2])
R[Cart::xxxxx][Cart::s][Cart::xxy][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xxy][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xxy][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xxy][1]-gfak*R[Cart::xxx][Cart::s][Cart::xxy][2])+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::xy][2]
R[Cart::xxxxx][Cart::s][Cart::xyz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xyz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xyz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xyz][1]-gfak*R[Cart::xxx][Cart::s][Cart::xyz][2])+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::yz][2]
R[Cart::xxxxx][Cart::s][Cart::yzz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::yzz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::yzz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::yzz][1]-gfak*R[Cart::xxx][Cart::s][Cart::yzz][2])
R[Cart::xxxxx][Cart::s][Cart::xxx][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xxx][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xxx][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xxx][1]-gfak*R[Cart::xxx][Cart::s][Cart::xxx][2])+0.5/_decay*3*R[Cart::xxxx][Cart::s][Cart::xx][2]
R[Cart::xxxxx][Cart::s][Cart::xxz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xxz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xxz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xxz][1]-gfak*R[Cart::xxx][Cart::s][Cart::xxz][2])+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::xz][2]
R[Cart::xxxxx][Cart::s][Cart::xzz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::xzz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::xzz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::xzz][1]-gfak*R[Cart::xxx][Cart::s][Cart::xzz][2])+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::zz][2]
R[Cart::xxxxx][Cart::s][Cart::zzz][1]+=pma0*R[Cart::xxxx][Cart::s][Cart::zzz][1]+wmp0*R[Cart::xxxx][Cart::s][Cart::zzz][2]+3*rzeta*(R[Cart::xxx][Cart::s][Cart::zzz][1]-gfak*R[Cart::xxx][Cart::s][Cart::zzz][2])
R[Cart::xxxxz][Cart::s][Cart::yyy][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::yyy][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::yyy][2]
R[Cart::xxxxz][Cart::s][Cart::xyy][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xyy][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xyy][2]
R[Cart::xxxxz][Cart::s][Cart::yyz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::yyz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::yyz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::yy][2]
R[Cart::xxxxz][Cart::s][Cart::xxy][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xxy][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xxy][2]
R[Cart::xxxxz][Cart::s][Cart::xyz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xyz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::xy][2]
R[Cart::xxxxz][Cart::s][Cart::yzz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::yzz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::yzz][2]+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::yz][2]
R[Cart::xxxxz][Cart::s][Cart::xxx][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xxx][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xxx][2]
R[Cart::xxxxz][Cart::s][Cart::xxz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xxz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xxz][2]+0.5/_decay*1*R[Cart::xxxx][Cart::s][Cart::xx][2]
R[Cart::xxxxz][Cart::s][Cart::xzz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::xzz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::xzz][2]+0.5/_decay*2*R[Cart::xxxx][Cart::s][Cart::xz][2]
R[Cart::xxxxz][Cart::s][Cart::zzz][1]+=pma2*R[Cart::xxxx][Cart::s][Cart::zzz][1]+wmp2*R[Cart::xxxx][Cart::s][Cart::zzz][2]+0.5/_decay*3*R[Cart::xxxx][Cart::s][Cart::zz][2]
R[Cart::xxxzz][Cart::s][Cart::yyy][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::yyy][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::yyy][2]
R[Cart::xxxzz][Cart::s][Cart::xyy][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xyy][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xyy][2]
R[Cart::xxxzz][Cart::s][Cart::yyz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::yyz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::yyz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::yy][2]
R[Cart::xxxzz][Cart::s][Cart::xxy][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xxy][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xxy][2]
R[Cart::xxxzz][Cart::s][Cart::xyz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xyz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::xy][2]
R[Cart::xxxzz][Cart::s][Cart::yzz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::yzz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::yzz][2]+0.5/_decay*2*R[Cart::xxxz][Cart::s][Cart::yz][2]
R[Cart::xxxzz][Cart::s][Cart::xxx][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xxx][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xxx][2]
R[Cart::xxxzz][Cart::s][Cart::xxz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xxz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xxz][2]+0.5/_decay*1*R[Cart::xxxz][Cart::s][Cart::xx][2]
R[Cart::xxxzz][Cart::s][Cart::xzz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::xzz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::xzz][2]+0.5/_decay*2*R[Cart::xxxz][Cart::s][Cart::xz][2]
R[Cart::xxxzz][Cart::s][Cart::zzz][1]+=pma2*R[Cart::xxxz][Cart::s][Cart::zzz][1]+wmp2*R[Cart::xxxz][Cart::s][Cart::zzz][2]+0.5/_decay*3*R[Cart::xxxz][Cart::s][Cart::zz][2]
R[Cart::xxzzz][Cart::s][Cart::yyy][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::yyy][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::yyy][2]
R[Cart::xxzzz][Cart::s][Cart::xyy][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xyy][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::xzzz][Cart::s][Cart::yy][2]
R[Cart::xxzzz][Cart::s][Cart::yyz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::yyz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::yyz][2]
R[Cart::xxzzz][Cart::s][Cart::xxy][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xxy][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::xzzz][Cart::s][Cart::xy][2]
R[Cart::xxzzz][Cart::s][Cart::xyz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xyz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::xzzz][Cart::s][Cart::yz][2]
R[Cart::xxzzz][Cart::s][Cart::yzz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::yzz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::yzz][2]
R[Cart::xxzzz][Cart::s][Cart::xxx][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xxx][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::xzzz][Cart::s][Cart::xx][2]
R[Cart::xxzzz][Cart::s][Cart::xxz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xxz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::xzzz][Cart::s][Cart::xz][2]
R[Cart::xxzzz][Cart::s][Cart::xzz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::xzz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::xzzz][Cart::s][Cart::zz][2]
R[Cart::xxzzz][Cart::s][Cart::zzz][1]+=pma0*R[Cart::xzzz][Cart::s][Cart::zzz][1]+wmp0*R[Cart::xzzz][Cart::s][Cart::zzz][2]
R[Cart::xzzzz][Cart::s][Cart::yyy][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::yyy][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::yyy][2]
R[Cart::xzzzz][Cart::s][Cart::xyy][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xyy][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xyy][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::yy][2]
R[Cart::xzzzz][Cart::s][Cart::yyz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::yyz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::yyz][2]
R[Cart::xzzzz][Cart::s][Cart::xxy][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xxy][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xxy][2]+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::xy][2]
R[Cart::xzzzz][Cart::s][Cart::xyz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xyz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xyz][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::yz][2]
R[Cart::xzzzz][Cart::s][Cart::yzz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::yzz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::yzz][2]
R[Cart::xzzzz][Cart::s][Cart::xxx][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xxx][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xxx][2]+0.5/_decay*3*R[Cart::zzzz][Cart::s][Cart::xx][2]
R[Cart::xzzzz][Cart::s][Cart::xxz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xxz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xxz][2]+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::xz][2]
R[Cart::xzzzz][Cart::s][Cart::xzz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::xzz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::xzz][2]+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::zz][2]
R[Cart::xzzzz][Cart::s][Cart::zzz][1]+=pma0*R[Cart::zzzz][Cart::s][Cart::zzz][1]+wmp0*R[Cart::zzzz][Cart::s][Cart::zzz][2]
R[Cart::zzzzz][Cart::s][Cart::yyy][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::yyy][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::yyy][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::yyy][1]-gfak*R[Cart::zzz][Cart::s][Cart::yyy][2])
R[Cart::zzzzz][Cart::s][Cart::xyy][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xyy][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xyy][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xyy][1]-gfak*R[Cart::zzz][Cart::s][Cart::xyy][2])
R[Cart::zzzzz][Cart::s][Cart::yyz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::yyz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::yyz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::yyz][1]-gfak*R[Cart::zzz][Cart::s][Cart::yyz][2])+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::yy][2]
R[Cart::zzzzz][Cart::s][Cart::xxy][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xxy][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xxy][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xxy][1]-gfak*R[Cart::zzz][Cart::s][Cart::xxy][2])
R[Cart::zzzzz][Cart::s][Cart::xyz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xyz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xyz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xyz][1]-gfak*R[Cart::zzz][Cart::s][Cart::xyz][2])+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::xy][2]
R[Cart::zzzzz][Cart::s][Cart::yzz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::yzz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::yzz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::yzz][1]-gfak*R[Cart::zzz][Cart::s][Cart::yzz][2])+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::yz][2]
R[Cart::zzzzz][Cart::s][Cart::xxx][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xxx][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xxx][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xxx][1]-gfak*R[Cart::zzz][Cart::s][Cart::xxx][2])
R[Cart::zzzzz][Cart::s][Cart::xxz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xxz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xxz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xxz][1]-gfak*R[Cart::zzz][Cart::s][Cart::xxz][2])+0.5/_decay*1*R[Cart::zzzz][Cart::s][Cart::xx][2]
R[Cart::zzzzz][Cart::s][Cart::xzz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::xzz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::xzz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::xzz][1]-gfak*R[Cart::zzz][Cart::s][Cart::xzz][2])+0.5/_decay*2*R[Cart::zzzz][Cart::s][Cart::xz][2]
R[Cart::zzzzz][Cart::s][Cart::zzz][1]+=pma2*R[Cart::zzzz][Cart::s][Cart::zzz][1]+wmp2*R[Cart::zzzz][Cart::s][Cart::zzz][2]+3*rzeta*(R[Cart::zzz][Cart::s][Cart::zzz][1]-gfak*R[Cart::zzz][Cart::s][Cart::zzz][2])+0.5/_decay*3*R[Cart::zzzz][Cart::s][Cart::zz][2]
}}
//------------------------------------------------------

//Integral i - s - s - m0
if (_mmax >6 ){
if (_lmax_alpha>5){

R[Cart::yyyyyy][Cart::s][Cart::s][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::s][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::s][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::s][0]-gfak*R[Cart::yyyy][Cart::s][Cart::s][1])
R[Cart::xyyyyy][Cart::s][Cart::s][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::s][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::s][1]
R[Cart::yyyyyz][Cart::s][Cart::s][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::s][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::s][1]
R[Cart::xxyyyy][Cart::s][Cart::s][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::s][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::s][1]
R[Cart::xyyyyz][Cart::s][Cart::s][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::s][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::s][1]
R[Cart::yyyyzz][Cart::s][Cart::s][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::s][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::s][1]
R[Cart::xxxyyy][Cart::s][Cart::s][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::s][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::s][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::s][0]-gfak*R[Cart::xyyy][Cart::s][Cart::s][1])
R[Cart::xxyyyz][Cart::s][Cart::s][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::s][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::s][1]
R[Cart::xyyyzz][Cart::s][Cart::s][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::s][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::s][1]
R[Cart::yyyzzz][Cart::s][Cart::s][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::s][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::s][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::s][0]-gfak*R[Cart::yzzz][Cart::s][Cart::s][1])
R[Cart::xxxxyy][Cart::s][Cart::s][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::s][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::s][1]
R[Cart::xxxyyz][Cart::s][Cart::s][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::s][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::s][1]
R[Cart::xxyyzz][Cart::s][Cart::s][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::s][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::s][1]
R[Cart::xyyzzz][Cart::s][Cart::s][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::s][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::s][1]
R[Cart::yyzzzz][Cart::s][Cart::s][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::s][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::s][1]
R[Cart::xxxxxy][Cart::s][Cart::s][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::s][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::s][1]
R[Cart::xxxxyz][Cart::s][Cart::s][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::s][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::s][1]
R[Cart::xxxyzz][Cart::s][Cart::s][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::s][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::s][1]
R[Cart::xxyzzz][Cart::s][Cart::s][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::s][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::s][1]
R[Cart::xyzzzz][Cart::s][Cart::s][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::s][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::s][1]
R[Cart::yzzzzz][Cart::s][Cart::s][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::s][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::s][1]
R[Cart::xxxxxx][Cart::s][Cart::s][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::s][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::s][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::s][0]-gfak*R[Cart::xxxx][Cart::s][Cart::s][1])
R[Cart::xxxxxz][Cart::s][Cart::s][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::s][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::s][1]
R[Cart::xxxxzz][Cart::s][Cart::s][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::s][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::s][1]
R[Cart::xxxzzz][Cart::s][Cart::s][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::s][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::s][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::s][0]-gfak*R[Cart::xzzz][Cart::s][Cart::s][1])
R[Cart::xxzzzz][Cart::s][Cart::s][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::s][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::s][1]
R[Cart::xzzzzz][Cart::s][Cart::s][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::s][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::s][1]
R[Cart::zzzzzz][Cart::s][Cart::s][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::s][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::s][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::s][0]-gfak*R[Cart::zzzz][Cart::s][Cart::s][1])
}}
//------------------------------------------------------

//Integral i - s - p - m0
if (_mmax >7 ){
if (_lmax_alpha>5 && _lmax_gamma>0){

R[Cart::yyyyyy][Cart::s][Cart::y][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::y][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::y][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::y][0]-gfak*R[Cart::yyyy][Cart::s][Cart::y][1])+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::s][1]
R[Cart::yyyyyy][Cart::s][Cart::x][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::x][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::x][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::x][0]-gfak*R[Cart::yyyy][Cart::s][Cart::x][1])
R[Cart::yyyyyy][Cart::s][Cart::z][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::z][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::z][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::z][0]-gfak*R[Cart::yyyy][Cart::s][Cart::z][1])
R[Cart::xyyyyy][Cart::s][Cart::y][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::y][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::y][1]
R[Cart::xyyyyy][Cart::s][Cart::x][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::x][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::s][1]
R[Cart::xyyyyy][Cart::s][Cart::z][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::z][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::z][1]
R[Cart::yyyyyz][Cart::s][Cart::y][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::y][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::y][1]
R[Cart::yyyyyz][Cart::s][Cart::x][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::x][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::x][1]
R[Cart::yyyyyz][Cart::s][Cart::z][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::z][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::z][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::s][1]
R[Cart::xxyyyy][Cart::s][Cart::y][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::y][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::y][1]
R[Cart::xxyyyy][Cart::s][Cart::x][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::x][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::xyyyy][Cart::s][Cart::s][1]
R[Cart::xxyyyy][Cart::s][Cart::z][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::z][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::z][1]
R[Cart::xyyyyz][Cart::s][Cart::y][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::y][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::y][1]
R[Cart::xyyyyz][Cart::s][Cart::x][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::x][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::s][1]
R[Cart::xyyyyz][Cart::s][Cart::z][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::z][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::z][1]
R[Cart::yyyyzz][Cart::s][Cart::y][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::y][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::y][1]
R[Cart::yyyyzz][Cart::s][Cart::x][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::x][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::x][1]
R[Cart::yyyyzz][Cart::s][Cart::z][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::z][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::z][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::s][1]
R[Cart::xxxyyy][Cart::s][Cart::y][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::y][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::y][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::y][0]-gfak*R[Cart::xyyy][Cart::s][Cart::y][1])
R[Cart::xxxyyy][Cart::s][Cart::x][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::x][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::x][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::x][0]-gfak*R[Cart::xyyy][Cart::s][Cart::x][1])+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::s][1]
R[Cart::xxxyyy][Cart::s][Cart::z][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::z][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::z][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::z][0]-gfak*R[Cart::xyyy][Cart::s][Cart::z][1])
R[Cart::xxyyyz][Cart::s][Cart::y][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::y][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::y][1]
R[Cart::xxyyyz][Cart::s][Cart::x][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::x][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::x][1]
R[Cart::xxyyyz][Cart::s][Cart::z][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::z][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::z][1]+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::s][1]
R[Cart::xyyyzz][Cart::s][Cart::y][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::y][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::y][1]
R[Cart::xyyyzz][Cart::s][Cart::x][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::x][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::yyyzz][Cart::s][Cart::s][1]
R[Cart::xyyyzz][Cart::s][Cart::z][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::z][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::z][1]
R[Cart::yyyzzz][Cart::s][Cart::y][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::y][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::y][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::y][0]-gfak*R[Cart::yzzz][Cart::s][Cart::y][1])+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::s][1]
R[Cart::yyyzzz][Cart::s][Cart::x][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::x][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::x][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::x][0]-gfak*R[Cart::yzzz][Cart::s][Cart::x][1])
R[Cart::yyyzzz][Cart::s][Cart::z][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::z][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::z][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::z][0]-gfak*R[Cart::yzzz][Cart::s][Cart::z][1])
R[Cart::xxxxyy][Cart::s][Cart::y][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::y][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::xxxxy][Cart::s][Cart::s][1]
R[Cart::xxxxyy][Cart::s][Cart::x][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::x][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::x][1]
R[Cart::xxxxyy][Cart::s][Cart::z][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::z][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::z][1]
R[Cart::xxxyyz][Cart::s][Cart::y][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::y][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::y][1]
R[Cart::xxxyyz][Cart::s][Cart::x][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::x][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::x][1]
R[Cart::xxxyyz][Cart::s][Cart::z][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::z][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::z][1]+0.5/_decay*1*R[Cart::xxxyy][Cart::s][Cart::s][1]
R[Cart::xxyyzz][Cart::s][Cart::y][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::y][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::y][1]
R[Cart::xxyyzz][Cart::s][Cart::x][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::x][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::xyyzz][Cart::s][Cart::s][1]
R[Cart::xxyyzz][Cart::s][Cart::z][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::z][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::z][1]
R[Cart::xyyzzz][Cart::s][Cart::y][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::y][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::y][1]
R[Cart::xyyzzz][Cart::s][Cart::x][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::x][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::s][1]
R[Cart::xyyzzz][Cart::s][Cart::z][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::z][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::z][1]
R[Cart::yyzzzz][Cart::s][Cart::y][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::y][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::s][1]
R[Cart::yyzzzz][Cart::s][Cart::x][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::x][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::x][1]
R[Cart::yyzzzz][Cart::s][Cart::z][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::z][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::z][1]
R[Cart::xxxxxy][Cart::s][Cart::y][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::y][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::s][1]
R[Cart::xxxxxy][Cart::s][Cart::x][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::x][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::x][1]
R[Cart::xxxxxy][Cart::s][Cart::z][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::z][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::z][1]
R[Cart::xxxxyz][Cart::s][Cart::y][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::y][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::s][1]
R[Cart::xxxxyz][Cart::s][Cart::x][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::x][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::x][1]
R[Cart::xxxxyz][Cart::s][Cart::z][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::z][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::z][1]
R[Cart::xxxyzz][Cart::s][Cart::y][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::y][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::xxxzz][Cart::s][Cart::s][1]
R[Cart::xxxyzz][Cart::s][Cart::x][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::x][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::x][1]
R[Cart::xxxyzz][Cart::s][Cart::z][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::z][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::z][1]
R[Cart::xxyzzz][Cart::s][Cart::y][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::y][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::s][1]
R[Cart::xxyzzz][Cart::s][Cart::x][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::x][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::x][1]
R[Cart::xxyzzz][Cart::s][Cart::z][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::z][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::z][1]
R[Cart::xyzzzz][Cart::s][Cart::y][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::y][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::y][1]
R[Cart::xyzzzz][Cart::s][Cart::x][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::x][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::s][1]
R[Cart::xyzzzz][Cart::s][Cart::z][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::z][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::z][1]
R[Cart::yzzzzz][Cart::s][Cart::y][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::y][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::y][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::s][1]
R[Cart::yzzzzz][Cart::s][Cart::x][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::x][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::x][1]
R[Cart::yzzzzz][Cart::s][Cart::z][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::z][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::z][1]
R[Cart::xxxxxx][Cart::s][Cart::y][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::y][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::y][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::y][0]-gfak*R[Cart::xxxx][Cart::s][Cart::y][1])
R[Cart::xxxxxx][Cart::s][Cart::x][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::x][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::x][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::x][0]-gfak*R[Cart::xxxx][Cart::s][Cart::x][1])+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::s][1]
R[Cart::xxxxxx][Cart::s][Cart::z][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::z][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::z][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::z][0]-gfak*R[Cart::xxxx][Cart::s][Cart::z][1])
R[Cart::xxxxxz][Cart::s][Cart::y][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::y][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::y][1]
R[Cart::xxxxxz][Cart::s][Cart::x][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::x][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::x][1]
R[Cart::xxxxxz][Cart::s][Cart::z][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::z][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::z][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::s][1]
R[Cart::xxxxzz][Cart::s][Cart::y][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::y][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::y][1]
R[Cart::xxxxzz][Cart::s][Cart::x][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::x][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::x][1]
R[Cart::xxxxzz][Cart::s][Cart::z][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::z][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::z][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::s][1]
R[Cart::xxxzzz][Cart::s][Cart::y][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::y][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::y][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::y][0]-gfak*R[Cart::xzzz][Cart::s][Cart::y][1])
R[Cart::xxxzzz][Cart::s][Cart::x][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::x][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::x][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::x][0]-gfak*R[Cart::xzzz][Cart::s][Cart::x][1])+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::s][1]
R[Cart::xxxzzz][Cart::s][Cart::z][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::z][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::z][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::z][0]-gfak*R[Cart::xzzz][Cart::s][Cart::z][1])
R[Cart::xxzzzz][Cart::s][Cart::y][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::y][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::y][1]
R[Cart::xxzzzz][Cart::s][Cart::x][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::x][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::xzzzz][Cart::s][Cart::s][1]
R[Cart::xxzzzz][Cart::s][Cart::z][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::z][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::z][1]
R[Cart::xzzzzz][Cart::s][Cart::y][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::y][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::y][1]
R[Cart::xzzzzz][Cart::s][Cart::x][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::x][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::x][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::s][1]
R[Cart::xzzzzz][Cart::s][Cart::z][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::z][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::z][1]
R[Cart::zzzzzz][Cart::s][Cart::y][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::y][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::y][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::y][0]-gfak*R[Cart::zzzz][Cart::s][Cart::y][1])
R[Cart::zzzzzz][Cart::s][Cart::x][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::x][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::x][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::x][0]-gfak*R[Cart::zzzz][Cart::s][Cart::x][1])
R[Cart::zzzzzz][Cart::s][Cart::z][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::z][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::z][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::z][0]-gfak*R[Cart::zzzz][Cart::s][Cart::z][1])+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::s][1]
}}
//------------------------------------------------------

//Integral i - s - d - m0
if (_mmax >8 ){
if (_lmax_alpha>5 && _lmax_gamma>1){

R[Cart::yyyyyy][Cart::s][Cart::yy][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::yy][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::yy][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::yy][0]-gfak*R[Cart::yyyy][Cart::s][Cart::yy][1])+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::y][1]
R[Cart::yyyyyy][Cart::s][Cart::xy][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xy][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xy][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xy][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xy][1])+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::x][1]
R[Cart::yyyyyy][Cart::s][Cart::yz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::yz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::yz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::yz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::yz][1])+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::z][1]
R[Cart::yyyyyy][Cart::s][Cart::xx][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xx][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xx][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xx][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xx][1])
R[Cart::yyyyyy][Cart::s][Cart::xz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xz][1])
R[Cart::yyyyyy][Cart::s][Cart::zz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::zz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::zz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::zz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::zz][1])
R[Cart::xyyyyy][Cart::s][Cart::yy][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::yy][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::yy][1]
R[Cart::xyyyyy][Cart::s][Cart::xy][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xy][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::y][1]
R[Cart::xyyyyy][Cart::s][Cart::yz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::yz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::yz][1]
R[Cart::xyyyyy][Cart::s][Cart::xx][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xx][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::x][1]
R[Cart::xyyyyy][Cart::s][Cart::xz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::z][1]
R[Cart::xyyyyy][Cart::s][Cart::zz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::zz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::zz][1]
R[Cart::yyyyyz][Cart::s][Cart::yy][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::yy][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::yy][1]
R[Cart::yyyyyz][Cart::s][Cart::xy][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xy][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xy][1]
R[Cart::yyyyyz][Cart::s][Cart::yz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::yz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::y][1]
R[Cart::yyyyyz][Cart::s][Cart::xx][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xx][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xx][1]
R[Cart::yyyyyz][Cart::s][Cart::xz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::x][1]
R[Cart::yyyyyz][Cart::s][Cart::zz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::zz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::zz][1]+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::z][1]
R[Cart::xxyyyy][Cart::s][Cart::yy][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::yy][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::yy][1]
R[Cart::xxyyyy][Cart::s][Cart::xy][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xy][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xyyyy][Cart::s][Cart::y][1]
R[Cart::xxyyyy][Cart::s][Cart::yz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::yz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::yz][1]
R[Cart::xxyyyy][Cart::s][Cart::xx][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xx][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::xyyyy][Cart::s][Cart::x][1]
R[Cart::xxyyyy][Cart::s][Cart::xz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xyyyy][Cart::s][Cart::z][1]
R[Cart::xxyyyy][Cart::s][Cart::zz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::zz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::zz][1]
R[Cart::xyyyyz][Cart::s][Cart::yy][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::yy][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::yy][1]
R[Cart::xyyyyz][Cart::s][Cart::xy][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xy][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::y][1]
R[Cart::xyyyyz][Cart::s][Cart::yz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::yz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::yz][1]
R[Cart::xyyyyz][Cart::s][Cart::xx][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xx][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::yyyyz][Cart::s][Cart::x][1]
R[Cart::xyyyyz][Cart::s][Cart::xz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::z][1]
R[Cart::xyyyyz][Cart::s][Cart::zz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::zz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::zz][1]
R[Cart::yyyyzz][Cart::s][Cart::yy][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::yy][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::yy][1]
R[Cart::yyyyzz][Cart::s][Cart::xy][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xy][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xy][1]
R[Cart::yyyyzz][Cart::s][Cart::yz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::yz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::y][1]
R[Cart::yyyyzz][Cart::s][Cart::xx][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xx][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xx][1]
R[Cart::yyyyzz][Cart::s][Cart::xz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::x][1]
R[Cart::yyyyzz][Cart::s][Cart::zz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::zz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::zz][1]+0.5/_decay*2*R[Cart::yyyyz][Cart::s][Cart::z][1]
R[Cart::xxxyyy][Cart::s][Cart::yy][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::yy][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::yy][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::yy][0]-gfak*R[Cart::xyyy][Cart::s][Cart::yy][1])
R[Cart::xxxyyy][Cart::s][Cart::xy][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xy][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xy][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xy][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xy][1])+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::y][1]
R[Cart::xxxyyy][Cart::s][Cart::yz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::yz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::yz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::yz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::yz][1])
R[Cart::xxxyyy][Cart::s][Cart::xx][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xx][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xx][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xx][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xx][1])+0.5/_decay*2*R[Cart::xxyyy][Cart::s][Cart::x][1]
R[Cart::xxxyyy][Cart::s][Cart::xz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xz][1])+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::z][1]
R[Cart::xxxyyy][Cart::s][Cart::zz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::zz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::zz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::zz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::zz][1])
R[Cart::xxyyyz][Cart::s][Cart::yy][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::yy][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::yy][1]
R[Cart::xxyyyz][Cart::s][Cart::xy][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xy][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xy][1]
R[Cart::xxyyyz][Cart::s][Cart::yz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::yz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::y][1]
R[Cart::xxyyyz][Cart::s][Cart::xx][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xx][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xx][1]
R[Cart::xxyyyz][Cart::s][Cart::xz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::x][1]
R[Cart::xxyyyz][Cart::s][Cart::zz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::zz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::zz][1]+0.5/_decay*2*R[Cart::xxyyy][Cart::s][Cart::z][1]
R[Cart::xyyyzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::yy][1]
R[Cart::xyyyzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::yyyzz][Cart::s][Cart::y][1]
R[Cart::xyyyzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::yz][1]
R[Cart::xyyyzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::yyyzz][Cart::s][Cart::x][1]
R[Cart::xyyyzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yyyzz][Cart::s][Cart::z][1]
R[Cart::xyyyzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::zz][1]
R[Cart::yyyzzz][Cart::s][Cart::yy][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::yy][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::yy][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::yy][0]-gfak*R[Cart::yzzz][Cart::s][Cart::yy][1])+0.5/_decay*2*R[Cart::yyzzz][Cart::s][Cart::y][1]
R[Cart::yyyzzz][Cart::s][Cart::xy][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xy][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xy][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xy][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xy][1])+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::x][1]
R[Cart::yyyzzz][Cart::s][Cart::yz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::yz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::yz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::yz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::yz][1])+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::z][1]
R[Cart::yyyzzz][Cart::s][Cart::xx][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xx][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xx][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xx][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xx][1])
R[Cart::yyyzzz][Cart::s][Cart::xz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xz][1])
R[Cart::yyyzzz][Cart::s][Cart::zz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::zz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::zz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::zz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::zz][1])
R[Cart::xxxxyy][Cart::s][Cart::yy][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::yy][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::xxxxy][Cart::s][Cart::y][1]
R[Cart::xxxxyy][Cart::s][Cart::xy][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xy][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xxxxy][Cart::s][Cart::x][1]
R[Cart::xxxxyy][Cart::s][Cart::yz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::yz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxxy][Cart::s][Cart::z][1]
R[Cart::xxxxyy][Cart::s][Cart::xx][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xx][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xx][1]
R[Cart::xxxxyy][Cart::s][Cart::xz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xz][1]
R[Cart::xxxxyy][Cart::s][Cart::zz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::zz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::zz][1]
R[Cart::xxxyyz][Cart::s][Cart::yy][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::yy][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::yy][1]
R[Cart::xxxyyz][Cart::s][Cart::xy][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xy][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xy][1]
R[Cart::xxxyyz][Cart::s][Cart::yz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::yz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxyy][Cart::s][Cart::y][1]
R[Cart::xxxyyz][Cart::s][Cart::xx][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xx][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xx][1]
R[Cart::xxxyyz][Cart::s][Cart::xz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xxxyy][Cart::s][Cart::x][1]
R[Cart::xxxyyz][Cart::s][Cart::zz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::zz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::zz][1]+0.5/_decay*2*R[Cart::xxxyy][Cart::s][Cart::z][1]
R[Cart::xxyyzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::yy][1]
R[Cart::xxyyzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xyyzz][Cart::s][Cart::y][1]
R[Cart::xxyyzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::yz][1]
R[Cart::xxyyzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::xyyzz][Cart::s][Cart::x][1]
R[Cart::xxyyzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xyyzz][Cart::s][Cart::z][1]
R[Cart::xxyyzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::zz][1]
R[Cart::xyyzzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::yy][1]
R[Cart::xyyzzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::y][1]
R[Cart::xyyzzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::yz][1]
R[Cart::xyyzzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::yyzzz][Cart::s][Cart::x][1]
R[Cart::xyyzzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::z][1]
R[Cart::xyyzzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::zz][1]
R[Cart::yyzzzz][Cart::s][Cart::yy][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::yy][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::yzzzz][Cart::s][Cart::y][1]
R[Cart::yyzzzz][Cart::s][Cart::xy][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xy][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::x][1]
R[Cart::yyzzzz][Cart::s][Cart::yz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::yz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::z][1]
R[Cart::yyzzzz][Cart::s][Cart::xx][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xx][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xx][1]
R[Cart::yyzzzz][Cart::s][Cart::xz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xz][1]
R[Cart::yyzzzz][Cart::s][Cart::zz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::zz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::zz][1]
R[Cart::xxxxxy][Cart::s][Cart::yy][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::yy][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::y][1]
R[Cart::xxxxxy][Cart::s][Cart::xy][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xy][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::x][1]
R[Cart::xxxxxy][Cart::s][Cart::yz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::yz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::z][1]
R[Cart::xxxxxy][Cart::s][Cart::xx][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xx][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xx][1]
R[Cart::xxxxxy][Cart::s][Cart::xz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xz][1]
R[Cart::xxxxxy][Cart::s][Cart::zz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::zz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::zz][1]
R[Cart::xxxxyz][Cart::s][Cart::yy][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::yy][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::xxxxz][Cart::s][Cart::y][1]
R[Cart::xxxxyz][Cart::s][Cart::xy][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xy][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::x][1]
R[Cart::xxxxyz][Cart::s][Cart::yz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::yz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::z][1]
R[Cart::xxxxyz][Cart::s][Cart::xx][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xx][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xx][1]
R[Cart::xxxxyz][Cart::s][Cart::xz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xz][1]
R[Cart::xxxxyz][Cart::s][Cart::zz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::zz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::zz][1]
R[Cart::xxxyzz][Cart::s][Cart::yy][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::yy][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::xxxzz][Cart::s][Cart::y][1]
R[Cart::xxxyzz][Cart::s][Cart::xy][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xy][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xxxzz][Cart::s][Cart::x][1]
R[Cart::xxxyzz][Cart::s][Cart::yz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::yz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxzz][Cart::s][Cart::z][1]
R[Cart::xxxyzz][Cart::s][Cart::xx][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xx][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xx][1]
R[Cart::xxxyzz][Cart::s][Cart::xz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xz][1]
R[Cart::xxxyzz][Cart::s][Cart::zz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::zz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::zz][1]
R[Cart::xxyzzz][Cart::s][Cart::yy][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::yy][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::xxzzz][Cart::s][Cart::y][1]
R[Cart::xxyzzz][Cart::s][Cart::xy][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xy][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::x][1]
R[Cart::xxyzzz][Cart::s][Cart::yz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::yz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::z][1]
R[Cart::xxyzzz][Cart::s][Cart::xx][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xx][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xx][1]
R[Cart::xxyzzz][Cart::s][Cart::xz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xz][1]
R[Cart::xxyzzz][Cart::s][Cart::zz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::zz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::zz][1]
R[Cart::xyzzzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::yy][1]
R[Cart::xyzzzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::y][1]
R[Cart::xyzzzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::yz][1]
R[Cart::xyzzzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::yzzzz][Cart::s][Cart::x][1]
R[Cart::xyzzzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::z][1]
R[Cart::xyzzzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::zz][1]
R[Cart::yzzzzz][Cart::s][Cart::yy][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::yy][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::yy][1]+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::y][1]
R[Cart::yzzzzz][Cart::s][Cart::xy][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xy][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::x][1]
R[Cart::yzzzzz][Cart::s][Cart::yz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::yz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::z][1]
R[Cart::yzzzzz][Cart::s][Cart::xx][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xx][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xx][1]
R[Cart::yzzzzz][Cart::s][Cart::xz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xz][1]
R[Cart::yzzzzz][Cart::s][Cart::zz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::zz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::zz][1]
R[Cart::xxxxxx][Cart::s][Cart::yy][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::yy][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::yy][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::yy][0]-gfak*R[Cart::xxxx][Cart::s][Cart::yy][1])
R[Cart::xxxxxx][Cart::s][Cart::xy][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xy][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xy][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xy][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xy][1])+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::y][1]
R[Cart::xxxxxx][Cart::s][Cart::yz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::yz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::yz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::yz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::yz][1])
R[Cart::xxxxxx][Cart::s][Cart::xx][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xx][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xx][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xx][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xx][1])+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::x][1]
R[Cart::xxxxxx][Cart::s][Cart::xz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xz][1])+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::z][1]
R[Cart::xxxxxx][Cart::s][Cart::zz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::zz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::zz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::zz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::zz][1])
R[Cart::xxxxxz][Cart::s][Cart::yy][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::yy][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::yy][1]
R[Cart::xxxxxz][Cart::s][Cart::xy][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xy][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xy][1]
R[Cart::xxxxxz][Cart::s][Cart::yz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::yz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::y][1]
R[Cart::xxxxxz][Cart::s][Cart::xx][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xx][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xx][1]
R[Cart::xxxxxz][Cart::s][Cart::xz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::x][1]
R[Cart::xxxxxz][Cart::s][Cart::zz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::zz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::zz][1]+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::z][1]
R[Cart::xxxxzz][Cart::s][Cart::yy][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::yy][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::yy][1]
R[Cart::xxxxzz][Cart::s][Cart::xy][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xy][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xy][1]
R[Cart::xxxxzz][Cart::s][Cart::yz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::yz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::yz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::y][1]
R[Cart::xxxxzz][Cart::s][Cart::xx][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xx][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xx][1]
R[Cart::xxxxzz][Cart::s][Cart::xz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::x][1]
R[Cart::xxxxzz][Cart::s][Cart::zz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::zz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::zz][1]+0.5/_decay*2*R[Cart::xxxxz][Cart::s][Cart::z][1]
R[Cart::xxxzzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::yy][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::yy][0]-gfak*R[Cart::xzzz][Cart::s][Cart::yy][1])
R[Cart::xxxzzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xy][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xy][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xy][1])+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::y][1]
R[Cart::xxxzzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::yz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::yz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::yz][1])
R[Cart::xxxzzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xx][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xx][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xx][1])+0.5/_decay*2*R[Cart::xxzzz][Cart::s][Cart::x][1]
R[Cart::xxxzzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xz][1])+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::z][1]
R[Cart::xxxzzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::zz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::zz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::zz][1])
R[Cart::xxzzzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::yy][1]
R[Cart::xxzzzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::xzzzz][Cart::s][Cart::y][1]
R[Cart::xxzzzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::yz][1]
R[Cart::xxzzzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::xzzzz][Cart::s][Cart::x][1]
R[Cart::xxzzzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::xzzzz][Cart::s][Cart::z][1]
R[Cart::xxzzzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::zz][1]
R[Cart::xzzzzz][Cart::s][Cart::yy][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::yy][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::yy][1]
R[Cart::xzzzzz][Cart::s][Cart::xy][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xy][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xy][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::y][1]
R[Cart::xzzzzz][Cart::s][Cart::yz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::yz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::yz][1]
R[Cart::xzzzzz][Cart::s][Cart::xx][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xx][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xx][1]+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::x][1]
R[Cart::xzzzzz][Cart::s][Cart::xz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xz][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::z][1]
R[Cart::xzzzzz][Cart::s][Cart::zz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::zz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::zz][1]
R[Cart::zzzzzz][Cart::s][Cart::yy][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::yy][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::yy][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::yy][0]-gfak*R[Cart::zzzz][Cart::s][Cart::yy][1])
R[Cart::zzzzzz][Cart::s][Cart::xy][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xy][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xy][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xy][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xy][1])
R[Cart::zzzzzz][Cart::s][Cart::yz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::yz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::yz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::yz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::yz][1])+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::y][1]
R[Cart::zzzzzz][Cart::s][Cart::xx][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xx][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xx][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xx][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xx][1])
R[Cart::zzzzzz][Cart::s][Cart::xz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xz][1])+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::x][1]
R[Cart::zzzzzz][Cart::s][Cart::zz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::zz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::zz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::zz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::zz][1])+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::z][1]
}}
//------------------------------------------------------

//Integral i - s - f - m0
if (_mmax >9 ){
if (_lmax_alpha>5 && _lmax_gamma>2){

R[Cart::yyyyyy][Cart::s][Cart::yyy][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::yyy][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::yyy][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::yyy][0]-gfak*R[Cart::yyyy][Cart::s][Cart::yyy][1])+0.5/_decay*3*R[Cart::yyyyy][Cart::s][Cart::yy][1]
R[Cart::yyyyyy][Cart::s][Cart::xyy][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xyy][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xyy][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xyy][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xyy][1])+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::xy][1]
R[Cart::yyyyyy][Cart::s][Cart::yyz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::yyz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::yyz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::yyz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::yyz][1])+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::yz][1]
R[Cart::yyyyyy][Cart::s][Cart::xxy][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xxy][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xxy][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xxy][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xxy][1])+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::xx][1]
R[Cart::yyyyyy][Cart::s][Cart::xyz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xyz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xyz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xyz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xyz][1])+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::xz][1]
R[Cart::yyyyyy][Cart::s][Cart::yzz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::yzz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::yzz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::yzz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::yzz][1])+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::zz][1]
R[Cart::yyyyyy][Cart::s][Cart::xxx][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xxx][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xxx][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xxx][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xxx][1])
R[Cart::yyyyyy][Cart::s][Cart::xxz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xxz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xxz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xxz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xxz][1])
R[Cart::yyyyyy][Cart::s][Cart::xzz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::xzz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::xzz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::xzz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::xzz][1])
R[Cart::yyyyyy][Cart::s][Cart::zzz][0]+=pma1*R[Cart::yyyyy][Cart::s][Cart::zzz][0]+wmp1*R[Cart::yyyyy][Cart::s][Cart::zzz][1]+4*rzeta*(R[Cart::yyyy][Cart::s][Cart::zzz][0]-gfak*R[Cart::yyyy][Cart::s][Cart::zzz][1])
R[Cart::xyyyyy][Cart::s][Cart::yyy][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::yyy][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::yyy][1]
R[Cart::xyyyyy][Cart::s][Cart::xyy][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xyy][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::yy][1]
R[Cart::xyyyyy][Cart::s][Cart::yyz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::yyz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::yyz][1]
R[Cart::xyyyyy][Cart::s][Cart::xxy][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xxy][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::xy][1]
R[Cart::xyyyyy][Cart::s][Cart::xyz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xyz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::yz][1]
R[Cart::xyyyyy][Cart::s][Cart::yzz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::yzz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::yzz][1]
R[Cart::xyyyyy][Cart::s][Cart::xxx][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xxx][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::yyyyy][Cart::s][Cart::xx][1]
R[Cart::xyyyyy][Cart::s][Cart::xxz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xxz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::xz][1]
R[Cart::xyyyyy][Cart::s][Cart::xzz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::xzz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::zz][1]
R[Cart::xyyyyy][Cart::s][Cart::zzz][0]+=pma0*R[Cart::yyyyy][Cart::s][Cart::zzz][0]+wmp0*R[Cart::yyyyy][Cart::s][Cart::zzz][1]
R[Cart::yyyyyz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::yyy][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::yyy][1]
R[Cart::yyyyyz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xyy][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xyy][1]
R[Cart::yyyyyz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::yyz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::yyz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::yy][1]
R[Cart::yyyyyz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xxy][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xxy][1]
R[Cart::yyyyyz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xyz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::xy][1]
R[Cart::yyyyyz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::yzz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::yzz][1]+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::yz][1]
R[Cart::yyyyyz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xxx][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xxx][1]
R[Cart::yyyyyz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xxz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xxz][1]+0.5/_decay*1*R[Cart::yyyyy][Cart::s][Cart::xx][1]
R[Cart::yyyyyz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::xzz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::xzz][1]+0.5/_decay*2*R[Cart::yyyyy][Cart::s][Cart::xz][1]
R[Cart::yyyyyz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::yyyyy][Cart::s][Cart::zzz][0]+wmp2*R[Cart::yyyyy][Cart::s][Cart::zzz][1]+0.5/_decay*3*R[Cart::yyyyy][Cart::s][Cart::zz][1]
R[Cart::xxyyyy][Cart::s][Cart::yyy][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::yyy][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::yyy][1]
R[Cart::xxyyyy][Cart::s][Cart::xyy][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xyy][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::xyyyy][Cart::s][Cart::yy][1]
R[Cart::xxyyyy][Cart::s][Cart::yyz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::yyz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::yyz][1]
R[Cart::xxyyyy][Cart::s][Cart::xxy][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xxy][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::xyyyy][Cart::s][Cart::xy][1]
R[Cart::xxyyyy][Cart::s][Cart::xyz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xyz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xyyyy][Cart::s][Cart::yz][1]
R[Cart::xxyyyy][Cart::s][Cart::yzz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::yzz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::yzz][1]
R[Cart::xxyyyy][Cart::s][Cart::xxx][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xxx][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::xyyyy][Cart::s][Cart::xx][1]
R[Cart::xxyyyy][Cart::s][Cart::xxz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xxz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::xyyyy][Cart::s][Cart::xz][1]
R[Cart::xxyyyy][Cart::s][Cart::xzz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::xzz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::xyyyy][Cart::s][Cart::zz][1]
R[Cart::xxyyyy][Cart::s][Cart::zzz][0]+=pma0*R[Cart::xyyyy][Cart::s][Cart::zzz][0]+wmp0*R[Cart::xyyyy][Cart::s][Cart::zzz][1]
R[Cart::xyyyyz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::yyy][1]
R[Cart::xyyyyz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::yy][1]
R[Cart::xyyyyz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::yyz][1]
R[Cart::xyyyyz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::yyyyz][Cart::s][Cart::xy][1]
R[Cart::xyyyyz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::yz][1]
R[Cart::xyyyyz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::yzz][1]
R[Cart::xyyyyz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::yyyyz][Cart::s][Cart::xx][1]
R[Cart::xyyyyz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::yyyyz][Cart::s][Cart::xz][1]
R[Cart::xyyyyz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::zz][1]
R[Cart::xyyyyz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::yyyyz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::yyyyz][Cart::s][Cart::zzz][1]
R[Cart::yyyyzz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::yyy][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::yyy][1]
R[Cart::yyyyzz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xyy][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xyy][1]
R[Cart::yyyyzz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::yyz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::yyz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::yy][1]
R[Cart::yyyyzz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xxy][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xxy][1]
R[Cart::yyyyzz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xyz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::xy][1]
R[Cart::yyyyzz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::yzz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::yzz][1]+0.5/_decay*2*R[Cart::yyyyz][Cart::s][Cart::yz][1]
R[Cart::yyyyzz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xxx][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xxx][1]
R[Cart::yyyyzz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xxz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xxz][1]+0.5/_decay*1*R[Cart::yyyyz][Cart::s][Cart::xx][1]
R[Cart::yyyyzz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::xzz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::xzz][1]+0.5/_decay*2*R[Cart::yyyyz][Cart::s][Cart::xz][1]
R[Cart::yyyyzz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::yyyyz][Cart::s][Cart::zzz][0]+wmp2*R[Cart::yyyyz][Cart::s][Cart::zzz][1]+0.5/_decay*3*R[Cart::yyyyz][Cart::s][Cart::zz][1]
R[Cart::xxxyyy][Cart::s][Cart::yyy][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::yyy][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::yyy][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::yyy][0]-gfak*R[Cart::xyyy][Cart::s][Cart::yyy][1])
R[Cart::xxxyyy][Cart::s][Cart::xyy][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xyy][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xyy][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xyy][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xyy][1])+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::yy][1]
R[Cart::xxxyyy][Cart::s][Cart::yyz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::yyz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::yyz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::yyz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::yyz][1])
R[Cart::xxxyyy][Cart::s][Cart::xxy][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xxy][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xxy][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xxy][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xxy][1])+0.5/_decay*2*R[Cart::xxyyy][Cart::s][Cart::xy][1]
R[Cart::xxxyyy][Cart::s][Cart::xyz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xyz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xyz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xyz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xyz][1])+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::yz][1]
R[Cart::xxxyyy][Cart::s][Cart::yzz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::yzz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::yzz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::yzz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::yzz][1])
R[Cart::xxxyyy][Cart::s][Cart::xxx][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xxx][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xxx][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xxx][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xxx][1])+0.5/_decay*3*R[Cart::xxyyy][Cart::s][Cart::xx][1]
R[Cart::xxxyyy][Cart::s][Cart::xxz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xxz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xxz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xxz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xxz][1])+0.5/_decay*2*R[Cart::xxyyy][Cart::s][Cart::xz][1]
R[Cart::xxxyyy][Cart::s][Cart::xzz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::xzz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::xzz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::xzz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::xzz][1])+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::zz][1]
R[Cart::xxxyyy][Cart::s][Cart::zzz][0]+=pma0*R[Cart::xxyyy][Cart::s][Cart::zzz][0]+wmp0*R[Cart::xxyyy][Cart::s][Cart::zzz][1]+1*rzeta*(R[Cart::xyyy][Cart::s][Cart::zzz][0]-gfak*R[Cart::xyyy][Cart::s][Cart::zzz][1])
R[Cart::xxyyyz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::yyy][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::yyy][1]
R[Cart::xxyyyz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xyy][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xyy][1]
R[Cart::xxyyyz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::yyz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::yyz][1]+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::yy][1]
R[Cart::xxyyyz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xxy][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xxy][1]
R[Cart::xxyyyz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xyz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::xy][1]
R[Cart::xxyyyz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::yzz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::yzz][1]+0.5/_decay*2*R[Cart::xxyyy][Cart::s][Cart::yz][1]
R[Cart::xxyyyz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xxx][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xxx][1]
R[Cart::xxyyyz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xxz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xxz][1]+0.5/_decay*1*R[Cart::xxyyy][Cart::s][Cart::xx][1]
R[Cart::xxyyyz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::xzz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::xzz][1]+0.5/_decay*2*R[Cart::xxyyy][Cart::s][Cart::xz][1]
R[Cart::xxyyyz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::xxyyy][Cart::s][Cart::zzz][0]+wmp2*R[Cart::xxyyy][Cart::s][Cart::zzz][1]+0.5/_decay*3*R[Cart::xxyyy][Cart::s][Cart::zz][1]
R[Cart::xyyyzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::yyy][1]
R[Cart::xyyyzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::yyyzz][Cart::s][Cart::yy][1]
R[Cart::xyyyzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::yyz][1]
R[Cart::xyyyzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::yyyzz][Cart::s][Cart::xy][1]
R[Cart::xyyyzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yyyzz][Cart::s][Cart::yz][1]
R[Cart::xyyyzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::yzz][1]
R[Cart::xyyyzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::yyyzz][Cart::s][Cart::xx][1]
R[Cart::xyyyzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::yyyzz][Cart::s][Cart::xz][1]
R[Cart::xyyyzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::yyyzz][Cart::s][Cart::zz][1]
R[Cart::xyyyzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::yyyzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::yyyzz][Cart::s][Cart::zzz][1]
R[Cart::yyyzzz][Cart::s][Cart::yyy][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::yyy][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::yyy][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::yyy][0]-gfak*R[Cart::yzzz][Cart::s][Cart::yyy][1])+0.5/_decay*3*R[Cart::yyzzz][Cart::s][Cart::yy][1]
R[Cart::yyyzzz][Cart::s][Cart::xyy][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xyy][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xyy][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xyy][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xyy][1])+0.5/_decay*2*R[Cart::yyzzz][Cart::s][Cart::xy][1]
R[Cart::yyyzzz][Cart::s][Cart::yyz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::yyz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::yyz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::yyz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::yyz][1])+0.5/_decay*2*R[Cart::yyzzz][Cart::s][Cart::yz][1]
R[Cart::yyyzzz][Cart::s][Cart::xxy][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xxy][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xxy][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xxy][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xxy][1])+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::xx][1]
R[Cart::yyyzzz][Cart::s][Cart::xyz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xyz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xyz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xyz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xyz][1])+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::xz][1]
R[Cart::yyyzzz][Cart::s][Cart::yzz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::yzz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::yzz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::yzz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::yzz][1])+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::zz][1]
R[Cart::yyyzzz][Cart::s][Cart::xxx][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xxx][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xxx][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xxx][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xxx][1])
R[Cart::yyyzzz][Cart::s][Cart::xxz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xxz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xxz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xxz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xxz][1])
R[Cart::yyyzzz][Cart::s][Cart::xzz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::xzz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::xzz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::xzz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::xzz][1])
R[Cart::yyyzzz][Cart::s][Cart::zzz][0]+=pma1*R[Cart::yyzzz][Cart::s][Cart::zzz][0]+wmp1*R[Cart::yyzzz][Cart::s][Cart::zzz][1]+1*rzeta*(R[Cart::yzzz][Cart::s][Cart::zzz][0]-gfak*R[Cart::yzzz][Cart::s][Cart::zzz][1])
R[Cart::xxxxyy][Cart::s][Cart::yyy][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::yyy][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::xxxxy][Cart::s][Cart::yy][1]
R[Cart::xxxxyy][Cart::s][Cart::xyy][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xyy][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::xxxxy][Cart::s][Cart::xy][1]
R[Cart::xxxxyy][Cart::s][Cart::yyz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::yyz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::xxxxy][Cart::s][Cart::yz][1]
R[Cart::xxxxyy][Cart::s][Cart::xxy][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xxy][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::xxxxy][Cart::s][Cart::xx][1]
R[Cart::xxxxyy][Cart::s][Cart::xyz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xyz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxxy][Cart::s][Cart::xz][1]
R[Cart::xxxxyy][Cart::s][Cart::yzz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::yzz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::xxxxy][Cart::s][Cart::zz][1]
R[Cart::xxxxyy][Cart::s][Cart::xxx][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xxx][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xxx][1]
R[Cart::xxxxyy][Cart::s][Cart::xxz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xxz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xxz][1]
R[Cart::xxxxyy][Cart::s][Cart::xzz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::xzz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::xzz][1]
R[Cart::xxxxyy][Cart::s][Cart::zzz][0]+=pma1*R[Cart::xxxxy][Cart::s][Cart::zzz][0]+wmp1*R[Cart::xxxxy][Cart::s][Cart::zzz][1]
R[Cart::xxxyyz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::yyy][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::yyy][1]
R[Cart::xxxyyz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xyy][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xyy][1]
R[Cart::xxxyyz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::yyz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::yyz][1]+0.5/_decay*1*R[Cart::xxxyy][Cart::s][Cart::yy][1]
R[Cart::xxxyyz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xxy][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xxy][1]
R[Cart::xxxyyz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xyz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxyy][Cart::s][Cart::xy][1]
R[Cart::xxxyyz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::yzz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::yzz][1]+0.5/_decay*2*R[Cart::xxxyy][Cart::s][Cart::yz][1]
R[Cart::xxxyyz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xxx][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xxx][1]
R[Cart::xxxyyz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xxz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xxz][1]+0.5/_decay*1*R[Cart::xxxyy][Cart::s][Cart::xx][1]
R[Cart::xxxyyz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::xzz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::xzz][1]+0.5/_decay*2*R[Cart::xxxyy][Cart::s][Cart::xz][1]
R[Cart::xxxyyz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::xxxyy][Cart::s][Cart::zzz][0]+wmp2*R[Cart::xxxyy][Cart::s][Cart::zzz][1]+0.5/_decay*3*R[Cart::xxxyy][Cart::s][Cart::zz][1]
R[Cart::xxyyzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::yyy][1]
R[Cart::xxyyzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::xyyzz][Cart::s][Cart::yy][1]
R[Cart::xxyyzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::yyz][1]
R[Cart::xxyyzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::xyyzz][Cart::s][Cart::xy][1]
R[Cart::xxyyzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xyyzz][Cart::s][Cart::yz][1]
R[Cart::xxyyzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::yzz][1]
R[Cart::xxyyzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::xyyzz][Cart::s][Cart::xx][1]
R[Cart::xxyyzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::xyyzz][Cart::s][Cart::xz][1]
R[Cart::xxyyzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::xyyzz][Cart::s][Cart::zz][1]
R[Cart::xxyyzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::xyyzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::xyyzz][Cart::s][Cart::zzz][1]
R[Cart::xyyzzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::yyy][1]
R[Cart::xyyzzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::yy][1]
R[Cart::xyyzzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::yyz][1]
R[Cart::xyyzzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::yyzzz][Cart::s][Cart::xy][1]
R[Cart::xyyzzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::yz][1]
R[Cart::xyyzzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::yzz][1]
R[Cart::xyyzzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::yyzzz][Cart::s][Cart::xx][1]
R[Cart::xyyzzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::yyzzz][Cart::s][Cart::xz][1]
R[Cart::xyyzzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::yyzzz][Cart::s][Cart::zz][1]
R[Cart::xyyzzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::yyzzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::yyzzz][Cart::s][Cart::zzz][1]
R[Cart::yyzzzz][Cart::s][Cart::yyy][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::yyy][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::yzzzz][Cart::s][Cart::yy][1]
R[Cart::yyzzzz][Cart::s][Cart::xyy][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xyy][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::yzzzz][Cart::s][Cart::xy][1]
R[Cart::yyzzzz][Cart::s][Cart::yyz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::yyz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::yzzzz][Cart::s][Cart::yz][1]
R[Cart::yyzzzz][Cart::s][Cart::xxy][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xxy][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::xx][1]
R[Cart::yyzzzz][Cart::s][Cart::xyz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xyz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::xz][1]
R[Cart::yyzzzz][Cart::s][Cart::yzz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::yzz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::zz][1]
R[Cart::yyzzzz][Cart::s][Cart::xxx][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xxx][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xxx][1]
R[Cart::yyzzzz][Cart::s][Cart::xxz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xxz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xxz][1]
R[Cart::yyzzzz][Cart::s][Cart::xzz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::xzz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::xzz][1]
R[Cart::yyzzzz][Cart::s][Cart::zzz][0]+=pma1*R[Cart::yzzzz][Cart::s][Cart::zzz][0]+wmp1*R[Cart::yzzzz][Cart::s][Cart::zzz][1]
R[Cart::xxxxxy][Cart::s][Cart::yyy][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::yyy][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::xxxxx][Cart::s][Cart::yy][1]
R[Cart::xxxxxy][Cart::s][Cart::xyy][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xyy][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::xy][1]
R[Cart::xxxxxy][Cart::s][Cart::yyz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::yyz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::yz][1]
R[Cart::xxxxxy][Cart::s][Cart::xxy][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xxy][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::xx][1]
R[Cart::xxxxxy][Cart::s][Cart::xyz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xyz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::xz][1]
R[Cart::xxxxxy][Cart::s][Cart::yzz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::yzz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::zz][1]
R[Cart::xxxxxy][Cart::s][Cart::xxx][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xxx][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xxx][1]
R[Cart::xxxxxy][Cart::s][Cart::xxz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xxz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xxz][1]
R[Cart::xxxxxy][Cart::s][Cart::xzz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::xzz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::xzz][1]
R[Cart::xxxxxy][Cart::s][Cart::zzz][0]+=pma1*R[Cart::xxxxx][Cart::s][Cart::zzz][0]+wmp1*R[Cart::xxxxx][Cart::s][Cart::zzz][1]
R[Cart::xxxxyz][Cart::s][Cart::yyy][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::yyy][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::xxxxz][Cart::s][Cart::yy][1]
R[Cart::xxxxyz][Cart::s][Cart::xyy][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xyy][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::xxxxz][Cart::s][Cart::xy][1]
R[Cart::xxxxyz][Cart::s][Cart::yyz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::yyz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::xxxxz][Cart::s][Cart::yz][1]
R[Cart::xxxxyz][Cart::s][Cart::xxy][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xxy][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::xx][1]
R[Cart::xxxxyz][Cart::s][Cart::xyz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xyz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::xz][1]
R[Cart::xxxxyz][Cart::s][Cart::yzz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::yzz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::zz][1]
R[Cart::xxxxyz][Cart::s][Cart::xxx][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xxx][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xxx][1]
R[Cart::xxxxyz][Cart::s][Cart::xxz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xxz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xxz][1]
R[Cart::xxxxyz][Cart::s][Cart::xzz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::xzz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::xzz][1]
R[Cart::xxxxyz][Cart::s][Cart::zzz][0]+=pma1*R[Cart::xxxxz][Cart::s][Cart::zzz][0]+wmp1*R[Cart::xxxxz][Cart::s][Cart::zzz][1]
R[Cart::xxxyzz][Cart::s][Cart::yyy][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::yyy][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::xxxzz][Cart::s][Cart::yy][1]
R[Cart::xxxyzz][Cart::s][Cart::xyy][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xyy][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::xxxzz][Cart::s][Cart::xy][1]
R[Cart::xxxyzz][Cart::s][Cart::yyz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::yyz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::xxxzz][Cart::s][Cart::yz][1]
R[Cart::xxxyzz][Cart::s][Cart::xxy][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xxy][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::xxxzz][Cart::s][Cart::xx][1]
R[Cart::xxxyzz][Cart::s][Cart::xyz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xyz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxzz][Cart::s][Cart::xz][1]
R[Cart::xxxyzz][Cart::s][Cart::yzz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::yzz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::xxxzz][Cart::s][Cart::zz][1]
R[Cart::xxxyzz][Cart::s][Cart::xxx][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xxx][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xxx][1]
R[Cart::xxxyzz][Cart::s][Cart::xxz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xxz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xxz][1]
R[Cart::xxxyzz][Cart::s][Cart::xzz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::xzz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::xzz][1]
R[Cart::xxxyzz][Cart::s][Cart::zzz][0]+=pma1*R[Cart::xxxzz][Cart::s][Cart::zzz][0]+wmp1*R[Cart::xxxzz][Cart::s][Cart::zzz][1]
R[Cart::xxyzzz][Cart::s][Cart::yyy][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::yyy][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::xxzzz][Cart::s][Cart::yy][1]
R[Cart::xxyzzz][Cart::s][Cart::xyy][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xyy][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::xxzzz][Cart::s][Cart::xy][1]
R[Cart::xxyzzz][Cart::s][Cart::yyz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::yyz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::xxzzz][Cart::s][Cart::yz][1]
R[Cart::xxyzzz][Cart::s][Cart::xxy][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xxy][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::xx][1]
R[Cart::xxyzzz][Cart::s][Cart::xyz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xyz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::xz][1]
R[Cart::xxyzzz][Cart::s][Cart::yzz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::yzz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::zz][1]
R[Cart::xxyzzz][Cart::s][Cart::xxx][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xxx][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xxx][1]
R[Cart::xxyzzz][Cart::s][Cart::xxz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xxz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xxz][1]
R[Cart::xxyzzz][Cart::s][Cart::xzz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::xzz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::xzz][1]
R[Cart::xxyzzz][Cart::s][Cart::zzz][0]+=pma1*R[Cart::xxzzz][Cart::s][Cart::zzz][0]+wmp1*R[Cart::xxzzz][Cart::s][Cart::zzz][1]
R[Cart::xyzzzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::yyy][1]
R[Cart::xyzzzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::yy][1]
R[Cart::xyzzzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::yyz][1]
R[Cart::xyzzzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::yzzzz][Cart::s][Cart::xy][1]
R[Cart::xyzzzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::yz][1]
R[Cart::xyzzzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::yzz][1]
R[Cart::xyzzzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::yzzzz][Cart::s][Cart::xx][1]
R[Cart::xyzzzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::yzzzz][Cart::s][Cart::xz][1]
R[Cart::xyzzzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::yzzzz][Cart::s][Cart::zz][1]
R[Cart::xyzzzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::yzzzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::yzzzz][Cart::s][Cart::zzz][1]
R[Cart::yzzzzz][Cart::s][Cart::yyy][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::yyy][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::yyy][1]+0.5/_decay*3*R[Cart::zzzzz][Cart::s][Cart::yy][1]
R[Cart::yzzzzz][Cart::s][Cart::xyy][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xyy][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xyy][1]+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::xy][1]
R[Cart::yzzzzz][Cart::s][Cart::yyz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::yyz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::yyz][1]+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::yz][1]
R[Cart::yzzzzz][Cart::s][Cart::xxy][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xxy][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xxy][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::xx][1]
R[Cart::yzzzzz][Cart::s][Cart::xyz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xyz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::xz][1]
R[Cart::yzzzzz][Cart::s][Cart::yzz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::yzz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::yzz][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::zz][1]
R[Cart::yzzzzz][Cart::s][Cart::xxx][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xxx][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xxx][1]
R[Cart::yzzzzz][Cart::s][Cart::xxz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xxz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xxz][1]
R[Cart::yzzzzz][Cart::s][Cart::xzz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::xzz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::xzz][1]
R[Cart::yzzzzz][Cart::s][Cart::zzz][0]+=pma1*R[Cart::zzzzz][Cart::s][Cart::zzz][0]+wmp1*R[Cart::zzzzz][Cart::s][Cart::zzz][1]
R[Cart::xxxxxx][Cart::s][Cart::yyy][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::yyy][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::yyy][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::yyy][0]-gfak*R[Cart::xxxx][Cart::s][Cart::yyy][1])
R[Cart::xxxxxx][Cart::s][Cart::xyy][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xyy][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xyy][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xyy][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xyy][1])+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::yy][1]
R[Cart::xxxxxx][Cart::s][Cart::yyz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::yyz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::yyz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::yyz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::yyz][1])
R[Cart::xxxxxx][Cart::s][Cart::xxy][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xxy][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xxy][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xxy][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xxy][1])+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::xy][1]
R[Cart::xxxxxx][Cart::s][Cart::xyz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xyz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xyz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xyz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xyz][1])+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::yz][1]
R[Cart::xxxxxx][Cart::s][Cart::yzz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::yzz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::yzz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::yzz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::yzz][1])
R[Cart::xxxxxx][Cart::s][Cart::xxx][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xxx][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xxx][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xxx][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xxx][1])+0.5/_decay*3*R[Cart::xxxxx][Cart::s][Cart::xx][1]
R[Cart::xxxxxx][Cart::s][Cart::xxz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xxz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xxz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xxz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xxz][1])+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::xz][1]
R[Cart::xxxxxx][Cart::s][Cart::xzz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::xzz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::xzz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::xzz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::xzz][1])+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::zz][1]
R[Cart::xxxxxx][Cart::s][Cart::zzz][0]+=pma0*R[Cart::xxxxx][Cart::s][Cart::zzz][0]+wmp0*R[Cart::xxxxx][Cart::s][Cart::zzz][1]+4*rzeta*(R[Cart::xxxx][Cart::s][Cart::zzz][0]-gfak*R[Cart::xxxx][Cart::s][Cart::zzz][1])
R[Cart::xxxxxz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::yyy][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::yyy][1]
R[Cart::xxxxxz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xyy][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xyy][1]
R[Cart::xxxxxz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::yyz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::yyz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::yy][1]
R[Cart::xxxxxz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xxy][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xxy][1]
R[Cart::xxxxxz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xyz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::xy][1]
R[Cart::xxxxxz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::yzz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::yzz][1]+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::yz][1]
R[Cart::xxxxxz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xxx][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xxx][1]
R[Cart::xxxxxz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xxz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xxz][1]+0.5/_decay*1*R[Cart::xxxxx][Cart::s][Cart::xx][1]
R[Cart::xxxxxz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::xzz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::xzz][1]+0.5/_decay*2*R[Cart::xxxxx][Cart::s][Cart::xz][1]
R[Cart::xxxxxz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::xxxxx][Cart::s][Cart::zzz][0]+wmp2*R[Cart::xxxxx][Cart::s][Cart::zzz][1]+0.5/_decay*3*R[Cart::xxxxx][Cart::s][Cart::zz][1]
R[Cart::xxxxzz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::yyy][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::yyy][1]
R[Cart::xxxxzz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xyy][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xyy][1]
R[Cart::xxxxzz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::yyz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::yyz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::yy][1]
R[Cart::xxxxzz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xxy][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xxy][1]
R[Cart::xxxxzz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xyz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::xy][1]
R[Cart::xxxxzz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::yzz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::yzz][1]+0.5/_decay*2*R[Cart::xxxxz][Cart::s][Cart::yz][1]
R[Cart::xxxxzz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xxx][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xxx][1]
R[Cart::xxxxzz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xxz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xxz][1]+0.5/_decay*1*R[Cart::xxxxz][Cart::s][Cart::xx][1]
R[Cart::xxxxzz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::xzz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::xzz][1]+0.5/_decay*2*R[Cart::xxxxz][Cart::s][Cart::xz][1]
R[Cart::xxxxzz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::xxxxz][Cart::s][Cart::zzz][0]+wmp2*R[Cart::xxxxz][Cart::s][Cart::zzz][1]+0.5/_decay*3*R[Cart::xxxxz][Cart::s][Cart::zz][1]
R[Cart::xxxzzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::yyy][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::yyy][0]-gfak*R[Cart::xzzz][Cart::s][Cart::yyy][1])
R[Cart::xxxzzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xyy][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xyy][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xyy][1])+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::yy][1]
R[Cart::xxxzzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::yyz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::yyz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::yyz][1])
R[Cart::xxxzzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xxy][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xxy][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xxy][1])+0.5/_decay*2*R[Cart::xxzzz][Cart::s][Cart::xy][1]
R[Cart::xxxzzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xyz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xyz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xyz][1])+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::yz][1]
R[Cart::xxxzzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::yzz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::yzz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::yzz][1])
R[Cart::xxxzzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xxx][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xxx][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xxx][1])+0.5/_decay*3*R[Cart::xxzzz][Cart::s][Cart::xx][1]
R[Cart::xxxzzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xxz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xxz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xxz][1])+0.5/_decay*2*R[Cart::xxzzz][Cart::s][Cart::xz][1]
R[Cart::xxxzzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::xzz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::xzz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::xzz][1])+0.5/_decay*1*R[Cart::xxzzz][Cart::s][Cart::zz][1]
R[Cart::xxxzzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::xxzzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::xxzzz][Cart::s][Cart::zzz][1]+1*rzeta*(R[Cart::xzzz][Cart::s][Cart::zzz][0]-gfak*R[Cart::xzzz][Cart::s][Cart::zzz][1])
R[Cart::xxzzzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::yyy][1]
R[Cart::xxzzzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::xzzzz][Cart::s][Cart::yy][1]
R[Cart::xxzzzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::yyz][1]
R[Cart::xxzzzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::xzzzz][Cart::s][Cart::xy][1]
R[Cart::xxzzzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::xzzzz][Cart::s][Cart::yz][1]
R[Cart::xxzzzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::yzz][1]
R[Cart::xxzzzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::xzzzz][Cart::s][Cart::xx][1]
R[Cart::xxzzzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::xzzzz][Cart::s][Cart::xz][1]
R[Cart::xxzzzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::xzzzz][Cart::s][Cart::zz][1]
R[Cart::xxzzzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::xzzzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::xzzzz][Cart::s][Cart::zzz][1]
R[Cart::xzzzzz][Cart::s][Cart::yyy][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::yyy][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::yyy][1]
R[Cart::xzzzzz][Cart::s][Cart::xyy][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xyy][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xyy][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::yy][1]
R[Cart::xzzzzz][Cart::s][Cart::yyz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::yyz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::yyz][1]
R[Cart::xzzzzz][Cart::s][Cart::xxy][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xxy][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xxy][1]+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::xy][1]
R[Cart::xzzzzz][Cart::s][Cart::xyz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xyz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xyz][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::yz][1]
R[Cart::xzzzzz][Cart::s][Cart::yzz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::yzz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::yzz][1]
R[Cart::xzzzzz][Cart::s][Cart::xxx][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xxx][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xxx][1]+0.5/_decay*3*R[Cart::zzzzz][Cart::s][Cart::xx][1]
R[Cart::xzzzzz][Cart::s][Cart::xxz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xxz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xxz][1]+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::xz][1]
R[Cart::xzzzzz][Cart::s][Cart::xzz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::xzz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::xzz][1]+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::zz][1]
R[Cart::xzzzzz][Cart::s][Cart::zzz][0]+=pma0*R[Cart::zzzzz][Cart::s][Cart::zzz][0]+wmp0*R[Cart::zzzzz][Cart::s][Cart::zzz][1]
R[Cart::zzzzzz][Cart::s][Cart::yyy][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::yyy][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::yyy][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::yyy][0]-gfak*R[Cart::zzzz][Cart::s][Cart::yyy][1])
R[Cart::zzzzzz][Cart::s][Cart::xyy][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xyy][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xyy][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xyy][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xyy][1])
R[Cart::zzzzzz][Cart::s][Cart::yyz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::yyz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::yyz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::yyz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::yyz][1])+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::yy][1]
R[Cart::zzzzzz][Cart::s][Cart::xxy][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xxy][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xxy][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xxy][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xxy][1])
R[Cart::zzzzzz][Cart::s][Cart::xyz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xyz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xyz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xyz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xyz][1])+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::xy][1]
R[Cart::zzzzzz][Cart::s][Cart::yzz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::yzz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::yzz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::yzz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::yzz][1])+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::yz][1]
R[Cart::zzzzzz][Cart::s][Cart::xxx][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xxx][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xxx][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xxx][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xxx][1])
R[Cart::zzzzzz][Cart::s][Cart::xxz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xxz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xxz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xxz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xxz][1])+0.5/_decay*1*R[Cart::zzzzz][Cart::s][Cart::xx][1]
R[Cart::zzzzzz][Cart::s][Cart::xzz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::xzz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::xzz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::xzz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::xzz][1])+0.5/_decay*2*R[Cart::zzzzz][Cart::s][Cart::xz][1]
R[Cart::zzzzzz][Cart::s][Cart::zzz][0]+=pma2*R[Cart::zzzzz][Cart::s][Cart::zzz][0]+wmp2*R[Cart::zzzzz][Cart::s][Cart::zzz][1]+4*rzeta*(R[Cart::zzzz][Cart::s][Cart::zzz][0]-gfak*R[Cart::zzzz][Cart::s][Cart::zzz][1])+0.5/_decay*3*R[Cart::zzzzz][Cart::s][Cart::zz][1]
}}
//------------------------------------------------------


            
            
            
            
            
            
     
            
            
            
            
            
            
            
            


                }
            }
        }

    }  
          
        
        
        
        
        
        void TCrawMatrix::XIntegrate(vector<double>& _FmT, const double& _T  ){
        
        const int _mm = _FmT.size() - 1;
        const double pi = boost::math::constants::pi<double>();
        if ( _mm < 0 || _mm > 10){
            cerr << "mm is: " << _mm << " This should not have happened!" << flush;
            exit(1);
        }
        
        if ( _T < 0.0 ) {
            cerr << "T is: " << _T << " This should not have happened!" << flush;
            exit(1);
        }
  
        if ( _T >= 10.0 ) {
            // forward iteration
            _FmT[0]=0.50*sqrt(pi/_T)* erf(sqrt(_T));

            for (int m = 1; m < _FmT.size(); m++ ){
                _FmT[m] = (2*m-1) * _FmT[m-1]/(2.0*_T) - exp(-_T)/(2.0*_T) ;
            }
        }

        if ( _T < 1e-10 ){
           for ( int m=0; m < _FmT.size(); m++){
               _FmT[m] = 1.0/(2.0*m+1.0) - _T/(2.0*m+3.0); 
           }
        }

        
        if ( _T >= 1e-10 && _T < 10.0 ){
            // backward iteration
            double fm = 0.0;
            for ( int m = 60; m >= _mm; m--){
                fm = (2.0*_T)/(2.0*m+1.0) * ( fm + exp(-_T)/(2.0*_T));
            } 
            _FmT[_mm] = fm;
            for (int m = _mm-1 ; m >= 0; m--){
                _FmT[m] = (2.0*_T)/(2.0*m+1.0) * (_FmT[m+1] + exp(-_T)/(2.0*_T));
            }
        }
        

    }
    }}
