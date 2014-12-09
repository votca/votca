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
            int _ncombined =this->getBlockSize(_lmax_alpha+_lmax_beta)
            typedef boost::multi_array<double, 3> ma_type;
            typedef boost::multi_array_types::extent_range range;
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            ma_type R;
            R.resize(extents[ range(0, _ncombined ) ][ range(0, _nbeta ) ][ range(0, _ngamma)]);
            //initialize to zero
            for (index i = 0; i != _ncombined; ++i) {
                for (index j = 0; j != _nbeta; ++j) {
                    for (index k = 0; k != _ngamma; ++k) {

                                       R[i][j][k] = 0.0;
                                   }
                               }
                           }
            
            //start vertical recurrence
            typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;

            for ( GaussianIterator italpha = _shell_alpha->firstGaussian(); italpha != _shell_alpha->lastGaussian(); ++italpha){

                const double& _decay_alpha = (*italpha)->decay;
            
                for ( GaussianIterator itbeta = _shell_beta->firstGaussian(); itbeta != _shell_beta->lastGaussian(); ++itbeta){
                    const double& _decay_beta = (*itbeta)->decay;
                    
                    for ( GaussianIterator itgamma = _shell_gamma->firstGaussian(); itgamma != _shell_gamma->lastGaussian(); ++itgamma){

                        const double& _decay_gamma = (*itgamma)->decay;
            

            
            
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
            

            ma_type R_temp;
            R_temp.resize(extents[ range(0, _nalpha+_nbeta ) ][ range(0, _ngamma ) ][ range(0, _mmax)]);
            //initialize to zero
            for (index i = 0; i != _nalpha+_nbeta; ++i) {
                for (index j = 0; j != _nbeta; ++j) {
                    for (index k = 0; k != _mmax; ++k) {

                                       R_temp[i][j][k] = 0.0;
                                   }
                               }
                           }
            
            
            vector<double> _FmT(_mmax, 0.0); 
           
            XIntegrate(_FmT, _T);
            
            double sss = 8*pow(2*pi,0.25)*pow(decay_alpha*decay_beta*decay_gamma,0.75)/((decay_alpha+decay_beta)*decay_gamma);
            
            for (int _i=0;_i<_mmax;_i++){
                R_temp[Cart::s][Cart::s][_i]=sss*_FmT[_i];
            }
            

//Integral s - s - p - m5
if (_mmax >1 ){
if (_lmax_gamma>0){

R_temp[Cart::s][Cart::y][5]+=wmc1*R_temp[Cart::s][Cart::s][6]
R_temp[Cart::s][Cart::x][5]+=wmc0*R_temp[Cart::s][Cart::s][6]
R_temp[Cart::s][Cart::z][5]+=wmc2*R_temp[Cart::s][Cart::s][6]
}}
//------------------------------------------------------

//Integral s - s - d - m4
if (_mmax >2 ){
if (_lmax_gamma>1){

R_temp[Cart::s][Cart::yy][4]+=wmc1*R_temp[Cart::s][Cart::y][5]
R_temp[Cart::s][Cart::xy][4]+=wmc0*R_temp[Cart::s][Cart::y][5]
R_temp[Cart::s][Cart::yz][4]+=wmc1*R_temp[Cart::s][Cart::z][5]
R_temp[Cart::s][Cart::xx][4]+=wmc0*R_temp[Cart::s][Cart::x][5]
R_temp[Cart::s][Cart::xz][4]+=wmc0*R_temp[Cart::s][Cart::z][5]
R_temp[Cart::s][Cart::zz][4]+=wmc2*R_temp[Cart::s][Cart::z][5]
}}
//------------------------------------------------------

//Integral s - s - f - m3
if (_mmax >3 ){
if (_lmax_gamma>2){

R_temp[Cart::s][Cart::yyy][3]+=wmc1*R_temp[Cart::s][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::s][Cart::y][3]-cfak*R_temp[Cart::s][Cart::y][4])
R_temp[Cart::s][Cart::xyy][3]+=wmc0*R_temp[Cart::s][Cart::yy][4]
R_temp[Cart::s][Cart::yyz][3]+=wmc2*R_temp[Cart::s][Cart::yy][4]
R_temp[Cart::s][Cart::xxy][3]+=wmc1*R_temp[Cart::s][Cart::xx][4]
R_temp[Cart::s][Cart::xyz][3]+=wmc0*R_temp[Cart::s][Cart::yz][4]
R_temp[Cart::s][Cart::yzz][3]+=wmc1*R_temp[Cart::s][Cart::zz][4]
R_temp[Cart::s][Cart::xxx][3]+=wmc0*R_temp[Cart::s][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::s][Cart::x][3]-cfak*R_temp[Cart::s][Cart::x][4])
R_temp[Cart::s][Cart::xxz][3]+=wmc2*R_temp[Cart::s][Cart::xx][4]
R_temp[Cart::s][Cart::xzz][3]+=wmc0*R_temp[Cart::s][Cart::zz][4]
R_temp[Cart::s][Cart::zzz][3]+=wmc2*R_temp[Cart::s][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::s][Cart::z][3]-cfak*R_temp[Cart::s][Cart::z][4])
}}
//------------------------------------------------------

//Integral p - s - s - m5
if (_mmax >1 ){
if (_lmax_alpha>0){

R_temp[Cart::y][Cart::s][5]+=pma1*R_temp[Cart::s][Cart::s][5]+wmp1*R_temp[Cart::s][Cart::s][6]
R_temp[Cart::x][Cart::s][5]+=pma0*R_temp[Cart::s][Cart::s][5]+wmp0*R_temp[Cart::s][Cart::s][6]
R_temp[Cart::z][Cart::s][5]+=pma2*R_temp[Cart::s][Cart::s][5]+wmp2*R_temp[Cart::s][Cart::s][6]
}}
//------------------------------------------------------

//Integral p - s - p - m5
if (_mmax >2 ){
if (_lmax_alpha>0 && _lmax_gamma>0){

R_temp[Cart::y][Cart::y][5]+=pma1*R_temp[Cart::s][Cart::y][5]+wmp1*R_temp[Cart::s][Cart::y][6]+0.5/_decay*1*R_temp[Cart::s][Cart::s][6]
R_temp[Cart::y][Cart::x][5]+=pma1*R_temp[Cart::s][Cart::x][5]+wmp1*R_temp[Cart::s][Cart::x][6]
R_temp[Cart::y][Cart::z][5]+=pma1*R_temp[Cart::s][Cart::z][5]+wmp1*R_temp[Cart::s][Cart::z][6]
R_temp[Cart::x][Cart::y][5]+=pma0*R_temp[Cart::s][Cart::y][5]+wmp0*R_temp[Cart::s][Cart::y][6]
R_temp[Cart::x][Cart::x][5]+=pma0*R_temp[Cart::s][Cart::x][5]+wmp0*R_temp[Cart::s][Cart::x][6]+0.5/_decay*1*R_temp[Cart::s][Cart::s][6]
R_temp[Cart::x][Cart::z][5]+=pma0*R_temp[Cart::s][Cart::z][5]+wmp0*R_temp[Cart::s][Cart::z][6]
R_temp[Cart::z][Cart::y][5]+=pma2*R_temp[Cart::s][Cart::y][5]+wmp2*R_temp[Cart::s][Cart::y][6]
R_temp[Cart::z][Cart::x][5]+=pma2*R_temp[Cart::s][Cart::x][5]+wmp2*R_temp[Cart::s][Cart::x][6]
R_temp[Cart::z][Cart::z][5]+=pma2*R_temp[Cart::s][Cart::z][5]+wmp2*R_temp[Cart::s][Cart::z][6]+0.5/_decay*1*R_temp[Cart::s][Cart::s][6]
}}
//------------------------------------------------------

//Integral p - s - d - m4
if (_mmax >3 ){
if (_lmax_alpha>0 && _lmax_gamma>1){

R_temp[Cart::y][Cart::yy][4]+=wmc1*R_temp[Cart::y][Cart::y][5]+0.5/_decay*1*R_temp[Cart::s][Cart::y][5]
R_temp[Cart::y][Cart::xy][4]+=wmc0*R_temp[Cart::y][Cart::y][5]
R_temp[Cart::y][Cart::yz][4]+=wmc1*R_temp[Cart::y][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::z][5]
R_temp[Cart::y][Cart::xx][4]+=wmc0*R_temp[Cart::y][Cart::x][5]
R_temp[Cart::y][Cart::xz][4]+=wmc0*R_temp[Cart::y][Cart::z][5]
R_temp[Cart::y][Cart::zz][4]+=wmc2*R_temp[Cart::y][Cart::z][5]
R_temp[Cart::x][Cart::yy][4]+=wmc1*R_temp[Cart::x][Cart::y][5]
R_temp[Cart::x][Cart::xy][4]+=wmc0*R_temp[Cart::x][Cart::y][5]+0.5/_decay*1*R_temp[Cart::s][Cart::y][5]
R_temp[Cart::x][Cart::yz][4]+=wmc1*R_temp[Cart::x][Cart::z][5]
R_temp[Cart::x][Cart::xx][4]+=wmc0*R_temp[Cart::x][Cart::x][5]+0.5/_decay*1*R_temp[Cart::s][Cart::x][5]
R_temp[Cart::x][Cart::xz][4]+=wmc0*R_temp[Cart::x][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::z][5]
R_temp[Cart::x][Cart::zz][4]+=wmc2*R_temp[Cart::x][Cart::z][5]
R_temp[Cart::z][Cart::yy][4]+=wmc1*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::z][Cart::xy][4]+=wmc0*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::z][Cart::yz][4]+=wmc1*R_temp[Cart::z][Cart::z][5]
R_temp[Cart::z][Cart::xx][4]+=wmc0*R_temp[Cart::z][Cart::x][5]
R_temp[Cart::z][Cart::xz][4]+=wmc0*R_temp[Cart::z][Cart::z][5]
R_temp[Cart::z][Cart::zz][4]+=wmc2*R_temp[Cart::z][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::z][5]
}}
//------------------------------------------------------

//Integral p - s - f - m3
if (_mmax >4 ){
if (_lmax_alpha>0 && _lmax_gamma>2){

R_temp[Cart::y][Cart::yyy][3]+=wmc1*R_temp[Cart::y][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::y][Cart::y][3]-cfak*R_temp[Cart::y][Cart::y][4])+0.5/_decay*1*R_temp[Cart::s][Cart::yy][4]
R_temp[Cart::y][Cart::xyy][3]+=wmc0*R_temp[Cart::y][Cart::yy][4]
R_temp[Cart::y][Cart::yyz][3]+=wmc2*R_temp[Cart::y][Cart::yy][4]
R_temp[Cart::y][Cart::xxy][3]+=wmc1*R_temp[Cart::y][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][4]
R_temp[Cart::y][Cart::xyz][3]+=wmc0*R_temp[Cart::y][Cart::yz][4]
R_temp[Cart::y][Cart::yzz][3]+=wmc1*R_temp[Cart::y][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][4]
R_temp[Cart::y][Cart::xxx][3]+=wmc0*R_temp[Cart::y][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::y][Cart::x][3]-cfak*R_temp[Cart::y][Cart::x][4])
R_temp[Cart::y][Cart::xxz][3]+=wmc2*R_temp[Cart::y][Cart::xx][4]
R_temp[Cart::y][Cart::xzz][3]+=wmc0*R_temp[Cart::y][Cart::zz][4]
R_temp[Cart::y][Cart::zzz][3]+=wmc2*R_temp[Cart::y][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::y][Cart::z][3]-cfak*R_temp[Cart::y][Cart::z][4])
R_temp[Cart::x][Cart::yyy][3]+=wmc1*R_temp[Cart::x][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::x][Cart::y][3]-cfak*R_temp[Cart::x][Cart::y][4])
R_temp[Cart::x][Cart::xyy][3]+=wmc0*R_temp[Cart::x][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][4]
R_temp[Cart::x][Cart::yyz][3]+=wmc2*R_temp[Cart::x][Cart::yy][4]
R_temp[Cart::x][Cart::xxy][3]+=wmc1*R_temp[Cart::x][Cart::xx][4]
R_temp[Cart::x][Cart::xyz][3]+=wmc0*R_temp[Cart::x][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::s][Cart::yz][4]
R_temp[Cart::x][Cart::yzz][3]+=wmc1*R_temp[Cart::x][Cart::zz][4]
R_temp[Cart::x][Cart::xxx][3]+=wmc0*R_temp[Cart::x][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::x][Cart::x][3]-cfak*R_temp[Cart::x][Cart::x][4])+0.5/_decay*1*R_temp[Cart::s][Cart::xx][4]
R_temp[Cart::x][Cart::xxz][3]+=wmc2*R_temp[Cart::x][Cart::xx][4]
R_temp[Cart::x][Cart::xzz][3]+=wmc0*R_temp[Cart::x][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][4]
R_temp[Cart::x][Cart::zzz][3]+=wmc2*R_temp[Cart::x][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::x][Cart::z][3]-cfak*R_temp[Cart::x][Cart::z][4])
R_temp[Cart::z][Cart::yyy][3]+=wmc1*R_temp[Cart::z][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::z][Cart::y][3]-cfak*R_temp[Cart::z][Cart::y][4])
R_temp[Cart::z][Cart::xyy][3]+=wmc0*R_temp[Cart::z][Cart::yy][4]
R_temp[Cart::z][Cart::yyz][3]+=wmc2*R_temp[Cart::z][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][4]
R_temp[Cart::z][Cart::xxy][3]+=wmc1*R_temp[Cart::z][Cart::xx][4]
R_temp[Cart::z][Cart::xyz][3]+=wmc0*R_temp[Cart::z][Cart::yz][4]
R_temp[Cart::z][Cart::yzz][3]+=wmc1*R_temp[Cart::z][Cart::zz][4]
R_temp[Cart::z][Cart::xxx][3]+=wmc0*R_temp[Cart::z][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::z][Cart::x][3]-cfak*R_temp[Cart::z][Cart::x][4])
R_temp[Cart::z][Cart::xxz][3]+=wmc2*R_temp[Cart::z][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][4]
R_temp[Cart::z][Cart::xzz][3]+=wmc0*R_temp[Cart::z][Cart::zz][4]
R_temp[Cart::z][Cart::zzz][3]+=wmc2*R_temp[Cart::z][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::z][Cart::z][3]-cfak*R_temp[Cart::z][Cart::z][4])+0.5/_decay*1*R_temp[Cart::s][Cart::zz][4]
}}
//------------------------------------------------------

//Integral d - s - s - m4
if (_mmax >2 ){
if (_lmax_alpha>1){

R_temp[Cart::yy][Cart::s][4]+=pma1*R_temp[Cart::y][Cart::s][4]+wmp1*R_temp[Cart::y][Cart::s][5]
R_temp[Cart::xy][Cart::s][4]+=pma0*R_temp[Cart::y][Cart::s][4]+wmp0*R_temp[Cart::y][Cart::s][5]
R_temp[Cart::yz][Cart::s][4]+=pma1*R_temp[Cart::z][Cart::s][4]+wmp1*R_temp[Cart::z][Cart::s][5]
R_temp[Cart::xx][Cart::s][4]+=pma0*R_temp[Cart::x][Cart::s][4]+wmp0*R_temp[Cart::x][Cart::s][5]
R_temp[Cart::xz][Cart::s][4]+=pma0*R_temp[Cart::z][Cart::s][4]+wmp0*R_temp[Cart::z][Cart::s][5]
R_temp[Cart::zz][Cart::s][4]+=pma2*R_temp[Cart::z][Cart::s][4]+wmp2*R_temp[Cart::z][Cart::s][5]
}}
//------------------------------------------------------

//Integral d - s - p - m4
if (_mmax >3 ){
if (_lmax_alpha>1 && _lmax_gamma>0){

R_temp[Cart::yy][Cart::y][4]+=pma1*R_temp[Cart::y][Cart::y][4]+wmp1*R_temp[Cart::y][Cart::y][5]+0.5/_decay*1*R_temp[Cart::y][Cart::s][5]
R_temp[Cart::yy][Cart::x][4]+=pma1*R_temp[Cart::y][Cart::x][4]+wmp1*R_temp[Cart::y][Cart::x][5]
R_temp[Cart::yy][Cart::z][4]+=pma1*R_temp[Cart::y][Cart::z][4]+wmp1*R_temp[Cart::y][Cart::z][5]
R_temp[Cart::xy][Cart::y][4]+=pma0*R_temp[Cart::y][Cart::y][4]+wmp0*R_temp[Cart::y][Cart::y][5]
R_temp[Cart::xy][Cart::x][4]+=pma0*R_temp[Cart::y][Cart::x][4]+wmp0*R_temp[Cart::y][Cart::x][5]+0.5/_decay*1*R_temp[Cart::y][Cart::s][5]
R_temp[Cart::xy][Cart::z][4]+=pma0*R_temp[Cart::y][Cart::z][4]+wmp0*R_temp[Cart::y][Cart::z][5]
R_temp[Cart::yz][Cart::y][4]+=pma1*R_temp[Cart::z][Cart::y][4]+wmp1*R_temp[Cart::z][Cart::y][5]+0.5/_decay*1*R_temp[Cart::z][Cart::s][5]
R_temp[Cart::yz][Cart::x][4]+=pma1*R_temp[Cart::z][Cart::x][4]+wmp1*R_temp[Cart::z][Cart::x][5]
R_temp[Cart::yz][Cart::z][4]+=pma1*R_temp[Cart::z][Cart::z][4]+wmp1*R_temp[Cart::z][Cart::z][5]
R_temp[Cart::xx][Cart::y][4]+=pma0*R_temp[Cart::x][Cart::y][4]+wmp0*R_temp[Cart::x][Cart::y][5]
R_temp[Cart::xx][Cart::x][4]+=pma0*R_temp[Cart::x][Cart::x][4]+wmp0*R_temp[Cart::x][Cart::x][5]+0.5/_decay*1*R_temp[Cart::x][Cart::s][5]
R_temp[Cart::xx][Cart::z][4]+=pma0*R_temp[Cart::x][Cart::z][4]+wmp0*R_temp[Cart::x][Cart::z][5]
R_temp[Cart::xz][Cart::y][4]+=pma0*R_temp[Cart::z][Cart::y][4]+wmp0*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::xz][Cart::x][4]+=pma0*R_temp[Cart::z][Cart::x][4]+wmp0*R_temp[Cart::z][Cart::x][5]+0.5/_decay*1*R_temp[Cart::z][Cart::s][5]
R_temp[Cart::xz][Cart::z][4]+=pma0*R_temp[Cart::z][Cart::z][4]+wmp0*R_temp[Cart::z][Cart::z][5]
R_temp[Cart::zz][Cart::y][4]+=pma2*R_temp[Cart::z][Cart::y][4]+wmp2*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::zz][Cart::x][4]+=pma2*R_temp[Cart::z][Cart::x][4]+wmp2*R_temp[Cart::z][Cart::x][5]
R_temp[Cart::zz][Cart::z][4]+=pma2*R_temp[Cart::z][Cart::z][4]+wmp2*R_temp[Cart::z][Cart::z][5]+0.5/_decay*1*R_temp[Cart::z][Cart::s][5]
}}
//------------------------------------------------------

//Integral d - s - d - m4
if (_mmax >4 ){
if (_lmax_alpha>1 && _lmax_gamma>1){

R_temp[Cart::yy][Cart::yy][4]+=pma1*R_temp[Cart::y][Cart::yy][4]+wmp1*R_temp[Cart::y][Cart::yy][5]+0.5/_decay*2*R_temp[Cart::y][Cart::y][5]
R_temp[Cart::yy][Cart::xy][4]+=pma1*R_temp[Cart::y][Cart::xy][4]+wmp1*R_temp[Cart::y][Cart::xy][5]+0.5/_decay*1*R_temp[Cart::y][Cart::x][5]
R_temp[Cart::yy][Cart::yz][4]+=pma1*R_temp[Cart::y][Cart::yz][4]+wmp1*R_temp[Cart::y][Cart::yz][5]+0.5/_decay*1*R_temp[Cart::y][Cart::z][5]
R_temp[Cart::yy][Cart::xx][4]+=pma1*R_temp[Cart::y][Cart::xx][4]+wmp1*R_temp[Cart::y][Cart::xx][5]
R_temp[Cart::yy][Cart::xz][4]+=pma1*R_temp[Cart::y][Cart::xz][4]+wmp1*R_temp[Cart::y][Cart::xz][5]
R_temp[Cart::yy][Cart::zz][4]+=pma1*R_temp[Cart::y][Cart::zz][4]+wmp1*R_temp[Cart::y][Cart::zz][5]
R_temp[Cart::xy][Cart::yy][4]+=pma0*R_temp[Cart::y][Cart::yy][4]+wmp0*R_temp[Cart::y][Cart::yy][5]
R_temp[Cart::xy][Cart::xy][4]+=pma0*R_temp[Cart::y][Cart::xy][4]+wmp0*R_temp[Cart::y][Cart::xy][5]+0.5/_decay*1*R_temp[Cart::y][Cart::y][5]
R_temp[Cart::xy][Cart::yz][4]+=pma0*R_temp[Cart::y][Cart::yz][4]+wmp0*R_temp[Cart::y][Cart::yz][5]
R_temp[Cart::xy][Cart::xx][4]+=pma0*R_temp[Cart::y][Cart::xx][4]+wmp0*R_temp[Cart::y][Cart::xx][5]+0.5/_decay*2*R_temp[Cart::y][Cart::x][5]
R_temp[Cart::xy][Cart::xz][4]+=pma0*R_temp[Cart::y][Cart::xz][4]+wmp0*R_temp[Cart::y][Cart::xz][5]+0.5/_decay*1*R_temp[Cart::y][Cart::z][5]
R_temp[Cart::xy][Cart::zz][4]+=pma0*R_temp[Cart::y][Cart::zz][4]+wmp0*R_temp[Cart::y][Cart::zz][5]
R_temp[Cart::yz][Cart::yy][4]+=pma1*R_temp[Cart::z][Cart::yy][4]+wmp1*R_temp[Cart::z][Cart::yy][5]+0.5/_decay*2*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::yz][Cart::xy][4]+=pma1*R_temp[Cart::z][Cart::xy][4]+wmp1*R_temp[Cart::z][Cart::xy][5]+0.5/_decay*1*R_temp[Cart::z][Cart::x][5]
R_temp[Cart::yz][Cart::yz][4]+=pma1*R_temp[Cart::z][Cart::yz][4]+wmp1*R_temp[Cart::z][Cart::yz][5]+0.5/_decay*1*R_temp[Cart::z][Cart::z][5]
R_temp[Cart::yz][Cart::xx][4]+=pma1*R_temp[Cart::z][Cart::xx][4]+wmp1*R_temp[Cart::z][Cart::xx][5]
R_temp[Cart::yz][Cart::xz][4]+=pma1*R_temp[Cart::z][Cart::xz][4]+wmp1*R_temp[Cart::z][Cart::xz][5]
R_temp[Cart::yz][Cart::zz][4]+=pma1*R_temp[Cart::z][Cart::zz][4]+wmp1*R_temp[Cart::z][Cart::zz][5]
R_temp[Cart::xx][Cart::yy][4]+=pma0*R_temp[Cart::x][Cart::yy][4]+wmp0*R_temp[Cart::x][Cart::yy][5]
R_temp[Cart::xx][Cart::xy][4]+=pma0*R_temp[Cart::x][Cart::xy][4]+wmp0*R_temp[Cart::x][Cart::xy][5]+0.5/_decay*1*R_temp[Cart::x][Cart::y][5]
R_temp[Cart::xx][Cart::yz][4]+=pma0*R_temp[Cart::x][Cart::yz][4]+wmp0*R_temp[Cart::x][Cart::yz][5]
R_temp[Cart::xx][Cart::xx][4]+=pma0*R_temp[Cart::x][Cart::xx][4]+wmp0*R_temp[Cart::x][Cart::xx][5]+0.5/_decay*2*R_temp[Cart::x][Cart::x][5]
R_temp[Cart::xx][Cart::xz][4]+=pma0*R_temp[Cart::x][Cart::xz][4]+wmp0*R_temp[Cart::x][Cart::xz][5]+0.5/_decay*1*R_temp[Cart::x][Cart::z][5]
R_temp[Cart::xx][Cart::zz][4]+=pma0*R_temp[Cart::x][Cart::zz][4]+wmp0*R_temp[Cart::x][Cart::zz][5]
R_temp[Cart::xz][Cart::yy][4]+=pma0*R_temp[Cart::z][Cart::yy][4]+wmp0*R_temp[Cart::z][Cart::yy][5]
R_temp[Cart::xz][Cart::xy][4]+=pma0*R_temp[Cart::z][Cart::xy][4]+wmp0*R_temp[Cart::z][Cart::xy][5]+0.5/_decay*1*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::xz][Cart::yz][4]+=pma0*R_temp[Cart::z][Cart::yz][4]+wmp0*R_temp[Cart::z][Cart::yz][5]
R_temp[Cart::xz][Cart::xx][4]+=pma0*R_temp[Cart::z][Cart::xx][4]+wmp0*R_temp[Cart::z][Cart::xx][5]+0.5/_decay*2*R_temp[Cart::z][Cart::x][5]
R_temp[Cart::xz][Cart::xz][4]+=pma0*R_temp[Cart::z][Cart::xz][4]+wmp0*R_temp[Cart::z][Cart::xz][5]+0.5/_decay*1*R_temp[Cart::z][Cart::z][5]
R_temp[Cart::xz][Cart::zz][4]+=pma0*R_temp[Cart::z][Cart::zz][4]+wmp0*R_temp[Cart::z][Cart::zz][5]
R_temp[Cart::zz][Cart::yy][4]+=pma2*R_temp[Cart::z][Cart::yy][4]+wmp2*R_temp[Cart::z][Cart::yy][5]
R_temp[Cart::zz][Cart::xy][4]+=pma2*R_temp[Cart::z][Cart::xy][4]+wmp2*R_temp[Cart::z][Cart::xy][5]
R_temp[Cart::zz][Cart::yz][4]+=pma2*R_temp[Cart::z][Cart::yz][4]+wmp2*R_temp[Cart::z][Cart::yz][5]+0.5/_decay*1*R_temp[Cart::z][Cart::y][5]
R_temp[Cart::zz][Cart::xx][4]+=pma2*R_temp[Cart::z][Cart::xx][4]+wmp2*R_temp[Cart::z][Cart::xx][5]
R_temp[Cart::zz][Cart::xz][4]+=pma2*R_temp[Cart::z][Cart::xz][4]+wmp2*R_temp[Cart::z][Cart::xz][5]+0.5/_decay*1*R_temp[Cart::z][Cart::x][5]
R_temp[Cart::zz][Cart::zz][4]+=pma2*R_temp[Cart::z][Cart::zz][4]+wmp2*R_temp[Cart::z][Cart::zz][5]+0.5/_decay*2*R_temp[Cart::z][Cart::z][5]
}}
//------------------------------------------------------

//Integral d - s - f - m3
if (_mmax >5 ){
if (_lmax_alpha>1 && _lmax_gamma>2){

R_temp[Cart::yy][Cart::yyy][3]+=wmc1*R_temp[Cart::yy][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::yy][Cart::y][3]-cfak*R_temp[Cart::yy][Cart::y][4])+0.5/_decay*2*R_temp[Cart::y][Cart::yy][4]
R_temp[Cart::yy][Cart::xyy][3]+=wmc0*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::yy][Cart::yyz][3]+=wmc2*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::yy][Cart::xxy][3]+=wmc1*R_temp[Cart::yy][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::y][Cart::xx][4]
R_temp[Cart::yy][Cart::xyz][3]+=wmc0*R_temp[Cart::yy][Cart::yz][4]
R_temp[Cart::yy][Cart::yzz][3]+=wmc1*R_temp[Cart::yy][Cart::zz][4]+0.5/_decay*2*R_temp[Cart::y][Cart::zz][4]
R_temp[Cart::yy][Cart::xxx][3]+=wmc0*R_temp[Cart::yy][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::yy][Cart::x][3]-cfak*R_temp[Cart::yy][Cart::x][4])
R_temp[Cart::yy][Cart::xxz][3]+=wmc2*R_temp[Cart::yy][Cart::xx][4]
R_temp[Cart::yy][Cart::xzz][3]+=wmc0*R_temp[Cart::yy][Cart::zz][4]
R_temp[Cart::yy][Cart::zzz][3]+=wmc2*R_temp[Cart::yy][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::yy][Cart::z][3]-cfak*R_temp[Cart::yy][Cart::z][4])
R_temp[Cart::xy][Cart::yyy][3]+=wmc1*R_temp[Cart::xy][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::xy][Cart::y][3]-cfak*R_temp[Cart::xy][Cart::y][4])+0.5/_decay*1*R_temp[Cart::x][Cart::yy][4]
R_temp[Cart::xy][Cart::xyy][3]+=wmc0*R_temp[Cart::xy][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][4]
R_temp[Cart::xy][Cart::yyz][3]+=wmc2*R_temp[Cart::xy][Cart::yy][4]
R_temp[Cart::xy][Cart::xxy][3]+=wmc1*R_temp[Cart::xy][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][4]
R_temp[Cart::xy][Cart::xyz][3]+=wmc0*R_temp[Cart::xy][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::y][Cart::yz][4]
R_temp[Cart::xy][Cart::yzz][3]+=wmc1*R_temp[Cart::xy][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::x][Cart::zz][4]
R_temp[Cart::xy][Cart::xxx][3]+=wmc0*R_temp[Cart::xy][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::xy][Cart::x][3]-cfak*R_temp[Cart::xy][Cart::x][4])+0.5/_decay*1*R_temp[Cart::y][Cart::xx][4]
R_temp[Cart::xy][Cart::xxz][3]+=wmc2*R_temp[Cart::xy][Cart::xx][4]
R_temp[Cart::xy][Cart::xzz][3]+=wmc0*R_temp[Cart::xy][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::y][Cart::zz][4]
R_temp[Cart::xy][Cart::zzz][3]+=wmc2*R_temp[Cart::xy][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::xy][Cart::z][3]-cfak*R_temp[Cart::xy][Cart::z][4])
R_temp[Cart::yz][Cart::yyy][3]+=wmc1*R_temp[Cart::yz][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::yz][Cart::y][3]-cfak*R_temp[Cart::yz][Cart::y][4])+0.5/_decay*1*R_temp[Cart::z][Cart::yy][4]
R_temp[Cart::yz][Cart::xyy][3]+=wmc0*R_temp[Cart::yz][Cart::yy][4]
R_temp[Cart::yz][Cart::yyz][3]+=wmc2*R_temp[Cart::yz][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][4]
R_temp[Cart::yz][Cart::xxy][3]+=wmc1*R_temp[Cart::yz][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::z][Cart::xx][4]
R_temp[Cart::yz][Cart::xyz][3]+=wmc0*R_temp[Cart::yz][Cart::yz][4]
R_temp[Cart::yz][Cart::yzz][3]+=wmc1*R_temp[Cart::yz][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][4]
R_temp[Cart::yz][Cart::xxx][3]+=wmc0*R_temp[Cart::yz][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::yz][Cart::x][3]-cfak*R_temp[Cart::yz][Cart::x][4])
R_temp[Cart::yz][Cart::xxz][3]+=wmc2*R_temp[Cart::yz][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::y][Cart::xx][4]
R_temp[Cart::yz][Cart::xzz][3]+=wmc0*R_temp[Cart::yz][Cart::zz][4]
R_temp[Cart::yz][Cart::zzz][3]+=wmc2*R_temp[Cart::yz][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::yz][Cart::z][3]-cfak*R_temp[Cart::yz][Cart::z][4])+0.5/_decay*1*R_temp[Cart::y][Cart::zz][4]
R_temp[Cart::xx][Cart::yyy][3]+=wmc1*R_temp[Cart::xx][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::xx][Cart::y][3]-cfak*R_temp[Cart::xx][Cart::y][4])
R_temp[Cart::xx][Cart::xyy][3]+=wmc0*R_temp[Cart::xx][Cart::yy][4]+0.5/_decay*2*R_temp[Cart::x][Cart::yy][4]
R_temp[Cart::xx][Cart::yyz][3]+=wmc2*R_temp[Cart::xx][Cart::yy][4]
R_temp[Cart::xx][Cart::xxy][3]+=wmc1*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xx][Cart::xyz][3]+=wmc0*R_temp[Cart::xx][Cart::yz][4]+0.5/_decay*2*R_temp[Cart::x][Cart::yz][4]
R_temp[Cart::xx][Cart::yzz][3]+=wmc1*R_temp[Cart::xx][Cart::zz][4]
R_temp[Cart::xx][Cart::xxx][3]+=wmc0*R_temp[Cart::xx][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::xx][Cart::x][3]-cfak*R_temp[Cart::xx][Cart::x][4])+0.5/_decay*2*R_temp[Cart::x][Cart::xx][4]
R_temp[Cart::xx][Cart::xxz][3]+=wmc2*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xx][Cart::xzz][3]+=wmc0*R_temp[Cart::xx][Cart::zz][4]+0.5/_decay*2*R_temp[Cart::x][Cart::zz][4]
R_temp[Cart::xx][Cart::zzz][3]+=wmc2*R_temp[Cart::xx][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::xx][Cart::z][3]-cfak*R_temp[Cart::xx][Cart::z][4])
R_temp[Cart::xz][Cart::yyy][3]+=wmc1*R_temp[Cart::xz][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::xz][Cart::y][3]-cfak*R_temp[Cart::xz][Cart::y][4])
R_temp[Cart::xz][Cart::xyy][3]+=wmc0*R_temp[Cart::xz][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::z][Cart::yy][4]
R_temp[Cart::xz][Cart::yyz][3]+=wmc2*R_temp[Cart::xz][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::x][Cart::yy][4]
R_temp[Cart::xz][Cart::xxy][3]+=wmc1*R_temp[Cart::xz][Cart::xx][4]
R_temp[Cart::xz][Cart::xyz][3]+=wmc0*R_temp[Cart::xz][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::yz][4]
R_temp[Cart::xz][Cart::yzz][3]+=wmc1*R_temp[Cart::xz][Cart::zz][4]
R_temp[Cart::xz][Cart::xxx][3]+=wmc0*R_temp[Cart::xz][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::xz][Cart::x][3]-cfak*R_temp[Cart::xz][Cart::x][4])+0.5/_decay*1*R_temp[Cart::z][Cart::xx][4]
R_temp[Cart::xz][Cart::xxz][3]+=wmc2*R_temp[Cart::xz][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][4]
R_temp[Cart::xz][Cart::xzz][3]+=wmc0*R_temp[Cart::xz][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][4]
R_temp[Cart::xz][Cart::zzz][3]+=wmc2*R_temp[Cart::xz][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::xz][Cart::z][3]-cfak*R_temp[Cart::xz][Cart::z][4])+0.5/_decay*1*R_temp[Cart::x][Cart::zz][4]
R_temp[Cart::zz][Cart::yyy][3]+=wmc1*R_temp[Cart::zz][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::zz][Cart::y][3]-cfak*R_temp[Cart::zz][Cart::y][4])
R_temp[Cart::zz][Cart::xyy][3]+=wmc0*R_temp[Cart::zz][Cart::yy][4]
R_temp[Cart::zz][Cart::yyz][3]+=wmc2*R_temp[Cart::zz][Cart::yy][4]+0.5/_decay*2*R_temp[Cart::z][Cart::yy][4]
R_temp[Cart::zz][Cart::xxy][3]+=wmc1*R_temp[Cart::zz][Cart::xx][4]
R_temp[Cart::zz][Cart::xyz][3]+=wmc0*R_temp[Cart::zz][Cart::yz][4]
R_temp[Cart::zz][Cart::yzz][3]+=wmc1*R_temp[Cart::zz][Cart::zz][4]
R_temp[Cart::zz][Cart::xxx][3]+=wmc0*R_temp[Cart::zz][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::zz][Cart::x][3]-cfak*R_temp[Cart::zz][Cart::x][4])
R_temp[Cart::zz][Cart::xxz][3]+=wmc2*R_temp[Cart::zz][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::z][Cart::xx][4]
R_temp[Cart::zz][Cart::xzz][3]+=wmc0*R_temp[Cart::zz][Cart::zz][4]
R_temp[Cart::zz][Cart::zzz][3]+=wmc2*R_temp[Cart::zz][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::zz][Cart::z][3]-cfak*R_temp[Cart::zz][Cart::z][4])+0.5/_decay*2*R_temp[Cart::z][Cart::zz][4]
}}
//------------------------------------------------------

//Integral f - s - s - m3
if (_mmax >3 ){
if (_lmax_alpha>2){

R_temp[Cart::yyy][Cart::s][3]+=pma1*R_temp[Cart::yy][Cart::s][3]+wmp1*R_temp[Cart::yy][Cart::s][4]+1*rzeta*(R_temp[Cart::y][Cart::s][3]-gfak*R_temp[Cart::y][Cart::s][4])
R_temp[Cart::xyy][Cart::s][3]+=pma0*R_temp[Cart::yy][Cart::s][3]+wmp0*R_temp[Cart::yy][Cart::s][4]
R_temp[Cart::yyz][Cart::s][3]+=pma2*R_temp[Cart::yy][Cart::s][3]+wmp2*R_temp[Cart::yy][Cart::s][4]
R_temp[Cart::xxy][Cart::s][3]+=pma1*R_temp[Cart::xx][Cart::s][3]+wmp1*R_temp[Cart::xx][Cart::s][4]
R_temp[Cart::xyz][Cart::s][3]+=pma0*R_temp[Cart::yz][Cart::s][3]+wmp0*R_temp[Cart::yz][Cart::s][4]
R_temp[Cart::yzz][Cart::s][3]+=pma1*R_temp[Cart::zz][Cart::s][3]+wmp1*R_temp[Cart::zz][Cart::s][4]
R_temp[Cart::xxx][Cart::s][3]+=pma0*R_temp[Cart::xx][Cart::s][3]+wmp0*R_temp[Cart::xx][Cart::s][4]+1*rzeta*(R_temp[Cart::x][Cart::s][3]-gfak*R_temp[Cart::x][Cart::s][4])
R_temp[Cart::xxz][Cart::s][3]+=pma2*R_temp[Cart::xx][Cart::s][3]+wmp2*R_temp[Cart::xx][Cart::s][4]
R_temp[Cart::xzz][Cart::s][3]+=pma0*R_temp[Cart::zz][Cart::s][3]+wmp0*R_temp[Cart::zz][Cart::s][4]
R_temp[Cart::zzz][Cart::s][3]+=pma2*R_temp[Cart::zz][Cart::s][3]+wmp2*R_temp[Cart::zz][Cart::s][4]+1*rzeta*(R_temp[Cart::z][Cart::s][3]-gfak*R_temp[Cart::z][Cart::s][4])
}}
//------------------------------------------------------

//Integral f - s - p - m3
if (_mmax >4 ){
if (_lmax_alpha>2 && _lmax_gamma>0){

R_temp[Cart::yyy][Cart::y][3]+=pma1*R_temp[Cart::yy][Cart::y][3]+wmp1*R_temp[Cart::yy][Cart::y][4]+1*rzeta*(R_temp[Cart::y][Cart::y][3]-gfak*R_temp[Cart::y][Cart::y][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::s][4]
R_temp[Cart::yyy][Cart::x][3]+=pma1*R_temp[Cart::yy][Cart::x][3]+wmp1*R_temp[Cart::yy][Cart::x][4]+1*rzeta*(R_temp[Cart::y][Cart::x][3]-gfak*R_temp[Cart::y][Cart::x][4])
R_temp[Cart::yyy][Cart::z][3]+=pma1*R_temp[Cart::yy][Cart::z][3]+wmp1*R_temp[Cart::yy][Cart::z][4]+1*rzeta*(R_temp[Cart::y][Cart::z][3]-gfak*R_temp[Cart::y][Cart::z][4])
R_temp[Cart::xyy][Cart::y][3]+=pma0*R_temp[Cart::yy][Cart::y][3]+wmp0*R_temp[Cart::yy][Cart::y][4]
R_temp[Cart::xyy][Cart::x][3]+=pma0*R_temp[Cart::yy][Cart::x][3]+wmp0*R_temp[Cart::yy][Cart::x][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][4]
R_temp[Cart::xyy][Cart::z][3]+=pma0*R_temp[Cart::yy][Cart::z][3]+wmp0*R_temp[Cart::yy][Cart::z][4]
R_temp[Cart::yyz][Cart::y][3]+=pma2*R_temp[Cart::yy][Cart::y][3]+wmp2*R_temp[Cart::yy][Cart::y][4]
R_temp[Cart::yyz][Cart::x][3]+=pma2*R_temp[Cart::yy][Cart::x][3]+wmp2*R_temp[Cart::yy][Cart::x][4]
R_temp[Cart::yyz][Cart::z][3]+=pma2*R_temp[Cart::yy][Cart::z][3]+wmp2*R_temp[Cart::yy][Cart::z][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][4]
R_temp[Cart::xxy][Cart::y][3]+=pma1*R_temp[Cart::xx][Cart::y][3]+wmp1*R_temp[Cart::xx][Cart::y][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][4]
R_temp[Cart::xxy][Cart::x][3]+=pma1*R_temp[Cart::xx][Cart::x][3]+wmp1*R_temp[Cart::xx][Cart::x][4]
R_temp[Cart::xxy][Cart::z][3]+=pma1*R_temp[Cart::xx][Cart::z][3]+wmp1*R_temp[Cart::xx][Cart::z][4]
R_temp[Cart::xyz][Cart::y][3]+=pma0*R_temp[Cart::yz][Cart::y][3]+wmp0*R_temp[Cart::yz][Cart::y][4]
R_temp[Cart::xyz][Cart::x][3]+=pma0*R_temp[Cart::yz][Cart::x][3]+wmp0*R_temp[Cart::yz][Cart::x][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::s][4]
R_temp[Cart::xyz][Cart::z][3]+=pma0*R_temp[Cart::yz][Cart::z][3]+wmp0*R_temp[Cart::yz][Cart::z][4]
R_temp[Cart::yzz][Cart::y][3]+=pma1*R_temp[Cart::zz][Cart::y][3]+wmp1*R_temp[Cart::zz][Cart::y][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][4]
R_temp[Cart::yzz][Cart::x][3]+=pma1*R_temp[Cart::zz][Cart::x][3]+wmp1*R_temp[Cart::zz][Cart::x][4]
R_temp[Cart::yzz][Cart::z][3]+=pma1*R_temp[Cart::zz][Cart::z][3]+wmp1*R_temp[Cart::zz][Cart::z][4]
R_temp[Cart::xxx][Cart::y][3]+=pma0*R_temp[Cart::xx][Cart::y][3]+wmp0*R_temp[Cart::xx][Cart::y][4]+1*rzeta*(R_temp[Cart::x][Cart::y][3]-gfak*R_temp[Cart::x][Cart::y][4])
R_temp[Cart::xxx][Cart::x][3]+=pma0*R_temp[Cart::xx][Cart::x][3]+wmp0*R_temp[Cart::xx][Cart::x][4]+1*rzeta*(R_temp[Cart::x][Cart::x][3]-gfak*R_temp[Cart::x][Cart::x][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::s][4]
R_temp[Cart::xxx][Cart::z][3]+=pma0*R_temp[Cart::xx][Cart::z][3]+wmp0*R_temp[Cart::xx][Cart::z][4]+1*rzeta*(R_temp[Cart::x][Cart::z][3]-gfak*R_temp[Cart::x][Cart::z][4])
R_temp[Cart::xxz][Cart::y][3]+=pma2*R_temp[Cart::xx][Cart::y][3]+wmp2*R_temp[Cart::xx][Cart::y][4]
R_temp[Cart::xxz][Cart::x][3]+=pma2*R_temp[Cart::xx][Cart::x][3]+wmp2*R_temp[Cart::xx][Cart::x][4]
R_temp[Cart::xxz][Cart::z][3]+=pma2*R_temp[Cart::xx][Cart::z][3]+wmp2*R_temp[Cart::xx][Cart::z][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][4]
R_temp[Cart::xzz][Cart::y][3]+=pma0*R_temp[Cart::zz][Cart::y][3]+wmp0*R_temp[Cart::zz][Cart::y][4]
R_temp[Cart::xzz][Cart::x][3]+=pma0*R_temp[Cart::zz][Cart::x][3]+wmp0*R_temp[Cart::zz][Cart::x][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][4]
R_temp[Cart::xzz][Cart::z][3]+=pma0*R_temp[Cart::zz][Cart::z][3]+wmp0*R_temp[Cart::zz][Cart::z][4]
R_temp[Cart::zzz][Cart::y][3]+=pma2*R_temp[Cart::zz][Cart::y][3]+wmp2*R_temp[Cart::zz][Cart::y][4]+1*rzeta*(R_temp[Cart::z][Cart::y][3]-gfak*R_temp[Cart::z][Cart::y][4])
R_temp[Cart::zzz][Cart::x][3]+=pma2*R_temp[Cart::zz][Cart::x][3]+wmp2*R_temp[Cart::zz][Cart::x][4]+1*rzeta*(R_temp[Cart::z][Cart::x][3]-gfak*R_temp[Cart::z][Cart::x][4])
R_temp[Cart::zzz][Cart::z][3]+=pma2*R_temp[Cart::zz][Cart::z][3]+wmp2*R_temp[Cart::zz][Cart::z][4]+1*rzeta*(R_temp[Cart::z][Cart::z][3]-gfak*R_temp[Cart::z][Cart::z][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::s][4]
}}
//------------------------------------------------------

//Integral f - s - d - m3
if (_mmax >5 ){
if (_lmax_alpha>2 && _lmax_gamma>1){

R_temp[Cart::yyy][Cart::yy][3]+=pma1*R_temp[Cart::yy][Cart::yy][3]+wmp1*R_temp[Cart::yy][Cart::yy][4]+1*rzeta*(R_temp[Cart::y][Cart::yy][3]-gfak*R_temp[Cart::y][Cart::yy][4])+0.5/_decay*2*R_temp[Cart::yy][Cart::y][4]
R_temp[Cart::yyy][Cart::xy][3]+=pma1*R_temp[Cart::yy][Cart::xy][3]+wmp1*R_temp[Cart::yy][Cart::xy][4]+1*rzeta*(R_temp[Cart::y][Cart::xy][3]-gfak*R_temp[Cart::y][Cart::xy][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::x][4]
R_temp[Cart::yyy][Cart::yz][3]+=pma1*R_temp[Cart::yy][Cart::yz][3]+wmp1*R_temp[Cart::yy][Cart::yz][4]+1*rzeta*(R_temp[Cart::y][Cart::yz][3]-gfak*R_temp[Cart::y][Cart::yz][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::z][4]
R_temp[Cart::yyy][Cart::xx][3]+=pma1*R_temp[Cart::yy][Cart::xx][3]+wmp1*R_temp[Cart::yy][Cart::xx][4]+1*rzeta*(R_temp[Cart::y][Cart::xx][3]-gfak*R_temp[Cart::y][Cart::xx][4])
R_temp[Cart::yyy][Cart::xz][3]+=pma1*R_temp[Cart::yy][Cart::xz][3]+wmp1*R_temp[Cart::yy][Cart::xz][4]+1*rzeta*(R_temp[Cart::y][Cart::xz][3]-gfak*R_temp[Cart::y][Cart::xz][4])
R_temp[Cart::yyy][Cart::zz][3]+=pma1*R_temp[Cart::yy][Cart::zz][3]+wmp1*R_temp[Cart::yy][Cart::zz][4]+1*rzeta*(R_temp[Cart::y][Cart::zz][3]-gfak*R_temp[Cart::y][Cart::zz][4])
R_temp[Cart::xyy][Cart::yy][3]+=pma0*R_temp[Cart::yy][Cart::yy][3]+wmp0*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::xyy][Cart::xy][3]+=pma0*R_temp[Cart::yy][Cart::xy][3]+wmp0*R_temp[Cart::yy][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][4]
R_temp[Cart::xyy][Cart::yz][3]+=pma0*R_temp[Cart::yy][Cart::yz][3]+wmp0*R_temp[Cart::yy][Cart::yz][4]
R_temp[Cart::xyy][Cart::xx][3]+=pma0*R_temp[Cart::yy][Cart::xx][3]+wmp0*R_temp[Cart::yy][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::yy][Cart::x][4]
R_temp[Cart::xyy][Cart::xz][3]+=pma0*R_temp[Cart::yy][Cart::xz][3]+wmp0*R_temp[Cart::yy][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::z][4]
R_temp[Cart::xyy][Cart::zz][3]+=pma0*R_temp[Cart::yy][Cart::zz][3]+wmp0*R_temp[Cart::yy][Cart::zz][4]
R_temp[Cart::yyz][Cart::yy][3]+=pma2*R_temp[Cart::yy][Cart::yy][3]+wmp2*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::yyz][Cart::xy][3]+=pma2*R_temp[Cart::yy][Cart::xy][3]+wmp2*R_temp[Cart::yy][Cart::xy][4]
R_temp[Cart::yyz][Cart::yz][3]+=pma2*R_temp[Cart::yy][Cart::yz][3]+wmp2*R_temp[Cart::yy][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][4]
R_temp[Cart::yyz][Cart::xx][3]+=pma2*R_temp[Cart::yy][Cart::xx][3]+wmp2*R_temp[Cart::yy][Cart::xx][4]
R_temp[Cart::yyz][Cart::xz][3]+=pma2*R_temp[Cart::yy][Cart::xz][3]+wmp2*R_temp[Cart::yy][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::x][4]
R_temp[Cart::yyz][Cart::zz][3]+=pma2*R_temp[Cart::yy][Cart::zz][3]+wmp2*R_temp[Cart::yy][Cart::zz][4]+0.5/_decay*2*R_temp[Cart::yy][Cart::z][4]
R_temp[Cart::xxy][Cart::yy][3]+=pma1*R_temp[Cart::xx][Cart::yy][3]+wmp1*R_temp[Cart::xx][Cart::yy][4]+0.5/_decay*2*R_temp[Cart::xx][Cart::y][4]
R_temp[Cart::xxy][Cart::xy][3]+=pma1*R_temp[Cart::xx][Cart::xy][3]+wmp1*R_temp[Cart::xx][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][4]
R_temp[Cart::xxy][Cart::yz][3]+=pma1*R_temp[Cart::xx][Cart::yz][3]+wmp1*R_temp[Cart::xx][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::z][4]
R_temp[Cart::xxy][Cart::xx][3]+=pma1*R_temp[Cart::xx][Cart::xx][3]+wmp1*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xxy][Cart::xz][3]+=pma1*R_temp[Cart::xx][Cart::xz][3]+wmp1*R_temp[Cart::xx][Cart::xz][4]
R_temp[Cart::xxy][Cart::zz][3]+=pma1*R_temp[Cart::xx][Cart::zz][3]+wmp1*R_temp[Cart::xx][Cart::zz][4]
R_temp[Cart::xyz][Cart::yy][3]+=pma0*R_temp[Cart::yz][Cart::yy][3]+wmp0*R_temp[Cart::yz][Cart::yy][4]
R_temp[Cart::xyz][Cart::xy][3]+=pma0*R_temp[Cart::yz][Cart::xy][3]+wmp0*R_temp[Cart::yz][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::y][4]
R_temp[Cart::xyz][Cart::yz][3]+=pma0*R_temp[Cart::yz][Cart::yz][3]+wmp0*R_temp[Cart::yz][Cart::yz][4]
R_temp[Cart::xyz][Cart::xx][3]+=pma0*R_temp[Cart::yz][Cart::xx][3]+wmp0*R_temp[Cart::yz][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::yz][Cart::x][4]
R_temp[Cart::xyz][Cart::xz][3]+=pma0*R_temp[Cart::yz][Cart::xz][3]+wmp0*R_temp[Cart::yz][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::z][4]
R_temp[Cart::xyz][Cart::zz][3]+=pma0*R_temp[Cart::yz][Cart::zz][3]+wmp0*R_temp[Cart::yz][Cart::zz][4]
R_temp[Cart::yzz][Cart::yy][3]+=pma1*R_temp[Cart::zz][Cart::yy][3]+wmp1*R_temp[Cart::zz][Cart::yy][4]+0.5/_decay*2*R_temp[Cart::zz][Cart::y][4]
R_temp[Cart::yzz][Cart::xy][3]+=pma1*R_temp[Cart::zz][Cart::xy][3]+wmp1*R_temp[Cart::zz][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::x][4]
R_temp[Cart::yzz][Cart::yz][3]+=pma1*R_temp[Cart::zz][Cart::yz][3]+wmp1*R_temp[Cart::zz][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][4]
R_temp[Cart::yzz][Cart::xx][3]+=pma1*R_temp[Cart::zz][Cart::xx][3]+wmp1*R_temp[Cart::zz][Cart::xx][4]
R_temp[Cart::yzz][Cart::xz][3]+=pma1*R_temp[Cart::zz][Cart::xz][3]+wmp1*R_temp[Cart::zz][Cart::xz][4]
R_temp[Cart::yzz][Cart::zz][3]+=pma1*R_temp[Cart::zz][Cart::zz][3]+wmp1*R_temp[Cart::zz][Cart::zz][4]
R_temp[Cart::xxx][Cart::yy][3]+=pma0*R_temp[Cart::xx][Cart::yy][3]+wmp0*R_temp[Cart::xx][Cart::yy][4]+1*rzeta*(R_temp[Cart::x][Cart::yy][3]-gfak*R_temp[Cart::x][Cart::yy][4])
R_temp[Cart::xxx][Cart::xy][3]+=pma0*R_temp[Cart::xx][Cart::xy][3]+wmp0*R_temp[Cart::xx][Cart::xy][4]+1*rzeta*(R_temp[Cart::x][Cart::xy][3]-gfak*R_temp[Cart::x][Cart::xy][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::y][4]
R_temp[Cart::xxx][Cart::yz][3]+=pma0*R_temp[Cart::xx][Cart::yz][3]+wmp0*R_temp[Cart::xx][Cart::yz][4]+1*rzeta*(R_temp[Cart::x][Cart::yz][3]-gfak*R_temp[Cart::x][Cart::yz][4])
R_temp[Cart::xxx][Cart::xx][3]+=pma0*R_temp[Cart::xx][Cart::xx][3]+wmp0*R_temp[Cart::xx][Cart::xx][4]+1*rzeta*(R_temp[Cart::x][Cart::xx][3]-gfak*R_temp[Cart::x][Cart::xx][4])+0.5/_decay*2*R_temp[Cart::xx][Cart::x][4]
R_temp[Cart::xxx][Cart::xz][3]+=pma0*R_temp[Cart::xx][Cart::xz][3]+wmp0*R_temp[Cart::xx][Cart::xz][4]+1*rzeta*(R_temp[Cart::x][Cart::xz][3]-gfak*R_temp[Cart::x][Cart::xz][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::z][4]
R_temp[Cart::xxx][Cart::zz][3]+=pma0*R_temp[Cart::xx][Cart::zz][3]+wmp0*R_temp[Cart::xx][Cart::zz][4]+1*rzeta*(R_temp[Cart::x][Cart::zz][3]-gfak*R_temp[Cart::x][Cart::zz][4])
R_temp[Cart::xxz][Cart::yy][3]+=pma2*R_temp[Cart::xx][Cart::yy][3]+wmp2*R_temp[Cart::xx][Cart::yy][4]
R_temp[Cart::xxz][Cart::xy][3]+=pma2*R_temp[Cart::xx][Cart::xy][3]+wmp2*R_temp[Cart::xx][Cart::xy][4]
R_temp[Cart::xxz][Cart::yz][3]+=pma2*R_temp[Cart::xx][Cart::yz][3]+wmp2*R_temp[Cart::xx][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::y][4]
R_temp[Cart::xxz][Cart::xx][3]+=pma2*R_temp[Cart::xx][Cart::xx][3]+wmp2*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xxz][Cart::xz][3]+=pma2*R_temp[Cart::xx][Cart::xz][3]+wmp2*R_temp[Cart::xx][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][4]
R_temp[Cart::xxz][Cart::zz][3]+=pma2*R_temp[Cart::xx][Cart::zz][3]+wmp2*R_temp[Cart::xx][Cart::zz][4]+0.5/_decay*2*R_temp[Cart::xx][Cart::z][4]
R_temp[Cart::xzz][Cart::yy][3]+=pma0*R_temp[Cart::zz][Cart::yy][3]+wmp0*R_temp[Cart::zz][Cart::yy][4]
R_temp[Cart::xzz][Cart::xy][3]+=pma0*R_temp[Cart::zz][Cart::xy][3]+wmp0*R_temp[Cart::zz][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::y][4]
R_temp[Cart::xzz][Cart::yz][3]+=pma0*R_temp[Cart::zz][Cart::yz][3]+wmp0*R_temp[Cart::zz][Cart::yz][4]
R_temp[Cart::xzz][Cart::xx][3]+=pma0*R_temp[Cart::zz][Cart::xx][3]+wmp0*R_temp[Cart::zz][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::zz][Cart::x][4]
R_temp[Cart::xzz][Cart::xz][3]+=pma0*R_temp[Cart::zz][Cart::xz][3]+wmp0*R_temp[Cart::zz][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][4]
R_temp[Cart::xzz][Cart::zz][3]+=pma0*R_temp[Cart::zz][Cart::zz][3]+wmp0*R_temp[Cart::zz][Cart::zz][4]
R_temp[Cart::zzz][Cart::yy][3]+=pma2*R_temp[Cart::zz][Cart::yy][3]+wmp2*R_temp[Cart::zz][Cart::yy][4]+1*rzeta*(R_temp[Cart::z][Cart::yy][3]-gfak*R_temp[Cart::z][Cart::yy][4])
R_temp[Cart::zzz][Cart::xy][3]+=pma2*R_temp[Cart::zz][Cart::xy][3]+wmp2*R_temp[Cart::zz][Cart::xy][4]+1*rzeta*(R_temp[Cart::z][Cart::xy][3]-gfak*R_temp[Cart::z][Cart::xy][4])
R_temp[Cart::zzz][Cart::yz][3]+=pma2*R_temp[Cart::zz][Cart::yz][3]+wmp2*R_temp[Cart::zz][Cart::yz][4]+1*rzeta*(R_temp[Cart::z][Cart::yz][3]-gfak*R_temp[Cart::z][Cart::yz][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::y][4]
R_temp[Cart::zzz][Cart::xx][3]+=pma2*R_temp[Cart::zz][Cart::xx][3]+wmp2*R_temp[Cart::zz][Cart::xx][4]+1*rzeta*(R_temp[Cart::z][Cart::xx][3]-gfak*R_temp[Cart::z][Cart::xx][4])
R_temp[Cart::zzz][Cart::xz][3]+=pma2*R_temp[Cart::zz][Cart::xz][3]+wmp2*R_temp[Cart::zz][Cart::xz][4]+1*rzeta*(R_temp[Cart::z][Cart::xz][3]-gfak*R_temp[Cart::z][Cart::xz][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::x][4]
R_temp[Cart::zzz][Cart::zz][3]+=pma2*R_temp[Cart::zz][Cart::zz][3]+wmp2*R_temp[Cart::zz][Cart::zz][4]+1*rzeta*(R_temp[Cart::z][Cart::zz][3]-gfak*R_temp[Cart::z][Cart::zz][4])+0.5/_decay*2*R_temp[Cart::zz][Cart::z][4]
}}
//------------------------------------------------------

//Integral f - s - f - m3
if (_mmax >6 ){
if (_lmax_alpha>2 && _lmax_gamma>2){

R_temp[Cart::yyy][Cart::yyy][3]+=pma1*R_temp[Cart::yy][Cart::yyy][3]+wmp1*R_temp[Cart::yy][Cart::yyy][4]+1*rzeta*(R_temp[Cart::y][Cart::yyy][3]-gfak*R_temp[Cart::y][Cart::yyy][4])+0.5/_decay*3*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::yyy][Cart::xyy][3]+=pma1*R_temp[Cart::yy][Cart::xyy][3]+wmp1*R_temp[Cart::yy][Cart::xyy][4]+1*rzeta*(R_temp[Cart::y][Cart::xyy][3]-gfak*R_temp[Cart::y][Cart::xyy][4])+0.5/_decay*2*R_temp[Cart::yy][Cart::xy][4]
R_temp[Cart::yyy][Cart::yyz][3]+=pma1*R_temp[Cart::yy][Cart::yyz][3]+wmp1*R_temp[Cart::yy][Cart::yyz][4]+1*rzeta*(R_temp[Cart::y][Cart::yyz][3]-gfak*R_temp[Cart::y][Cart::yyz][4])+0.5/_decay*2*R_temp[Cart::yy][Cart::yz][4]
R_temp[Cart::yyy][Cart::xxy][3]+=pma1*R_temp[Cart::yy][Cart::xxy][3]+wmp1*R_temp[Cart::yy][Cart::xxy][4]+1*rzeta*(R_temp[Cart::y][Cart::xxy][3]-gfak*R_temp[Cart::y][Cart::xxy][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::xx][4]
R_temp[Cart::yyy][Cart::xyz][3]+=pma1*R_temp[Cart::yy][Cart::xyz][3]+wmp1*R_temp[Cart::yy][Cart::xyz][4]+1*rzeta*(R_temp[Cart::y][Cart::xyz][3]-gfak*R_temp[Cart::y][Cart::xyz][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::xz][4]
R_temp[Cart::yyy][Cart::yzz][3]+=pma1*R_temp[Cart::yy][Cart::yzz][3]+wmp1*R_temp[Cart::yy][Cart::yzz][4]+1*rzeta*(R_temp[Cart::y][Cart::yzz][3]-gfak*R_temp[Cart::y][Cart::yzz][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::zz][4]
R_temp[Cart::yyy][Cart::xxx][3]+=pma1*R_temp[Cart::yy][Cart::xxx][3]+wmp1*R_temp[Cart::yy][Cart::xxx][4]+1*rzeta*(R_temp[Cart::y][Cart::xxx][3]-gfak*R_temp[Cart::y][Cart::xxx][4])
R_temp[Cart::yyy][Cart::xxz][3]+=pma1*R_temp[Cart::yy][Cart::xxz][3]+wmp1*R_temp[Cart::yy][Cart::xxz][4]+1*rzeta*(R_temp[Cart::y][Cart::xxz][3]-gfak*R_temp[Cart::y][Cart::xxz][4])
R_temp[Cart::yyy][Cart::xzz][3]+=pma1*R_temp[Cart::yy][Cart::xzz][3]+wmp1*R_temp[Cart::yy][Cart::xzz][4]+1*rzeta*(R_temp[Cart::y][Cart::xzz][3]-gfak*R_temp[Cart::y][Cart::xzz][4])
R_temp[Cart::yyy][Cart::zzz][3]+=pma1*R_temp[Cart::yy][Cart::zzz][3]+wmp1*R_temp[Cart::yy][Cart::zzz][4]+1*rzeta*(R_temp[Cart::y][Cart::zzz][3]-gfak*R_temp[Cart::y][Cart::zzz][4])
R_temp[Cart::xyy][Cart::yyy][3]+=pma0*R_temp[Cart::yy][Cart::yyy][3]+wmp0*R_temp[Cart::yy][Cart::yyy][4]
R_temp[Cart::xyy][Cart::xyy][3]+=pma0*R_temp[Cart::yy][Cart::xyy][3]+wmp0*R_temp[Cart::yy][Cart::xyy][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::xyy][Cart::yyz][3]+=pma0*R_temp[Cart::yy][Cart::yyz][3]+wmp0*R_temp[Cart::yy][Cart::yyz][4]
R_temp[Cart::xyy][Cart::xxy][3]+=pma0*R_temp[Cart::yy][Cart::xxy][3]+wmp0*R_temp[Cart::yy][Cart::xxy][4]+0.5/_decay*2*R_temp[Cart::yy][Cart::xy][4]
R_temp[Cart::xyy][Cart::xyz][3]+=pma0*R_temp[Cart::yy][Cart::xyz][3]+wmp0*R_temp[Cart::yy][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::yz][4]
R_temp[Cart::xyy][Cart::yzz][3]+=pma0*R_temp[Cart::yy][Cart::yzz][3]+wmp0*R_temp[Cart::yy][Cart::yzz][4]
R_temp[Cart::xyy][Cart::xxx][3]+=pma0*R_temp[Cart::yy][Cart::xxx][3]+wmp0*R_temp[Cart::yy][Cart::xxx][4]+0.5/_decay*3*R_temp[Cart::yy][Cart::xx][4]
R_temp[Cart::xyy][Cart::xxz][3]+=pma0*R_temp[Cart::yy][Cart::xxz][3]+wmp0*R_temp[Cart::yy][Cart::xxz][4]+0.5/_decay*2*R_temp[Cart::yy][Cart::xz][4]
R_temp[Cart::xyy][Cart::xzz][3]+=pma0*R_temp[Cart::yy][Cart::xzz][3]+wmp0*R_temp[Cart::yy][Cart::xzz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::zz][4]
R_temp[Cart::xyy][Cart::zzz][3]+=pma0*R_temp[Cart::yy][Cart::zzz][3]+wmp0*R_temp[Cart::yy][Cart::zzz][4]
R_temp[Cart::yyz][Cart::yyy][3]+=pma2*R_temp[Cart::yy][Cart::yyy][3]+wmp2*R_temp[Cart::yy][Cart::yyy][4]
R_temp[Cart::yyz][Cart::xyy][3]+=pma2*R_temp[Cart::yy][Cart::xyy][3]+wmp2*R_temp[Cart::yy][Cart::xyy][4]
R_temp[Cart::yyz][Cart::yyz][3]+=pma2*R_temp[Cart::yy][Cart::yyz][3]+wmp2*R_temp[Cart::yy][Cart::yyz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::yy][4]
R_temp[Cart::yyz][Cart::xxy][3]+=pma2*R_temp[Cart::yy][Cart::xxy][3]+wmp2*R_temp[Cart::yy][Cart::xxy][4]
R_temp[Cart::yyz][Cart::xyz][3]+=pma2*R_temp[Cart::yy][Cart::xyz][3]+wmp2*R_temp[Cart::yy][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::xy][4]
R_temp[Cart::yyz][Cart::yzz][3]+=pma2*R_temp[Cart::yy][Cart::yzz][3]+wmp2*R_temp[Cart::yy][Cart::yzz][4]+0.5/_decay*2*R_temp[Cart::yy][Cart::yz][4]
R_temp[Cart::yyz][Cart::xxx][3]+=pma2*R_temp[Cart::yy][Cart::xxx][3]+wmp2*R_temp[Cart::yy][Cart::xxx][4]
R_temp[Cart::yyz][Cart::xxz][3]+=pma2*R_temp[Cart::yy][Cart::xxz][3]+wmp2*R_temp[Cart::yy][Cart::xxz][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::xx][4]
R_temp[Cart::yyz][Cart::xzz][3]+=pma2*R_temp[Cart::yy][Cart::xzz][3]+wmp2*R_temp[Cart::yy][Cart::xzz][4]+0.5/_decay*2*R_temp[Cart::yy][Cart::xz][4]
R_temp[Cart::yyz][Cart::zzz][3]+=pma2*R_temp[Cart::yy][Cart::zzz][3]+wmp2*R_temp[Cart::yy][Cart::zzz][4]+0.5/_decay*3*R_temp[Cart::yy][Cart::zz][4]
R_temp[Cart::xxy][Cart::yyy][3]+=pma1*R_temp[Cart::xx][Cart::yyy][3]+wmp1*R_temp[Cart::xx][Cart::yyy][4]+0.5/_decay*3*R_temp[Cart::xx][Cart::yy][4]
R_temp[Cart::xxy][Cart::xyy][3]+=pma1*R_temp[Cart::xx][Cart::xyy][3]+wmp1*R_temp[Cart::xx][Cart::xyy][4]+0.5/_decay*2*R_temp[Cart::xx][Cart::xy][4]
R_temp[Cart::xxy][Cart::yyz][3]+=pma1*R_temp[Cart::xx][Cart::yyz][3]+wmp1*R_temp[Cart::xx][Cart::yyz][4]+0.5/_decay*2*R_temp[Cart::xx][Cart::yz][4]
R_temp[Cart::xxy][Cart::xxy][3]+=pma1*R_temp[Cart::xx][Cart::xxy][3]+wmp1*R_temp[Cart::xx][Cart::xxy][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xxy][Cart::xyz][3]+=pma1*R_temp[Cart::xx][Cart::xyz][3]+wmp1*R_temp[Cart::xx][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::xz][4]
R_temp[Cart::xxy][Cart::yzz][3]+=pma1*R_temp[Cart::xx][Cart::yzz][3]+wmp1*R_temp[Cart::xx][Cart::yzz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::zz][4]
R_temp[Cart::xxy][Cart::xxx][3]+=pma1*R_temp[Cart::xx][Cart::xxx][3]+wmp1*R_temp[Cart::xx][Cart::xxx][4]
R_temp[Cart::xxy][Cart::xxz][3]+=pma1*R_temp[Cart::xx][Cart::xxz][3]+wmp1*R_temp[Cart::xx][Cart::xxz][4]
R_temp[Cart::xxy][Cart::xzz][3]+=pma1*R_temp[Cart::xx][Cart::xzz][3]+wmp1*R_temp[Cart::xx][Cart::xzz][4]
R_temp[Cart::xxy][Cart::zzz][3]+=pma1*R_temp[Cart::xx][Cart::zzz][3]+wmp1*R_temp[Cart::xx][Cart::zzz][4]
R_temp[Cart::xyz][Cart::yyy][3]+=pma0*R_temp[Cart::yz][Cart::yyy][3]+wmp0*R_temp[Cart::yz][Cart::yyy][4]
R_temp[Cart::xyz][Cart::xyy][3]+=pma0*R_temp[Cart::yz][Cart::xyy][3]+wmp0*R_temp[Cart::yz][Cart::xyy][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::yy][4]
R_temp[Cart::xyz][Cart::yyz][3]+=pma0*R_temp[Cart::yz][Cart::yyz][3]+wmp0*R_temp[Cart::yz][Cart::yyz][4]
R_temp[Cart::xyz][Cart::xxy][3]+=pma0*R_temp[Cart::yz][Cart::xxy][3]+wmp0*R_temp[Cart::yz][Cart::xxy][4]+0.5/_decay*2*R_temp[Cart::yz][Cart::xy][4]
R_temp[Cart::xyz][Cart::xyz][3]+=pma0*R_temp[Cart::yz][Cart::xyz][3]+wmp0*R_temp[Cart::yz][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::yz][4]
R_temp[Cart::xyz][Cart::yzz][3]+=pma0*R_temp[Cart::yz][Cart::yzz][3]+wmp0*R_temp[Cart::yz][Cart::yzz][4]
R_temp[Cart::xyz][Cart::xxx][3]+=pma0*R_temp[Cart::yz][Cart::xxx][3]+wmp0*R_temp[Cart::yz][Cart::xxx][4]+0.5/_decay*3*R_temp[Cart::yz][Cart::xx][4]
R_temp[Cart::xyz][Cart::xxz][3]+=pma0*R_temp[Cart::yz][Cart::xxz][3]+wmp0*R_temp[Cart::yz][Cart::xxz][4]+0.5/_decay*2*R_temp[Cart::yz][Cart::xz][4]
R_temp[Cart::xyz][Cart::xzz][3]+=pma0*R_temp[Cart::yz][Cart::xzz][3]+wmp0*R_temp[Cart::yz][Cart::xzz][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::zz][4]
R_temp[Cart::xyz][Cart::zzz][3]+=pma0*R_temp[Cart::yz][Cart::zzz][3]+wmp0*R_temp[Cart::yz][Cart::zzz][4]
R_temp[Cart::yzz][Cart::yyy][3]+=pma1*R_temp[Cart::zz][Cart::yyy][3]+wmp1*R_temp[Cart::zz][Cart::yyy][4]+0.5/_decay*3*R_temp[Cart::zz][Cart::yy][4]
R_temp[Cart::yzz][Cart::xyy][3]+=pma1*R_temp[Cart::zz][Cart::xyy][3]+wmp1*R_temp[Cart::zz][Cart::xyy][4]+0.5/_decay*2*R_temp[Cart::zz][Cart::xy][4]
R_temp[Cart::yzz][Cart::yyz][3]+=pma1*R_temp[Cart::zz][Cart::yyz][3]+wmp1*R_temp[Cart::zz][Cart::yyz][4]+0.5/_decay*2*R_temp[Cart::zz][Cart::yz][4]
R_temp[Cart::yzz][Cart::xxy][3]+=pma1*R_temp[Cart::zz][Cart::xxy][3]+wmp1*R_temp[Cart::zz][Cart::xxy][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::xx][4]
R_temp[Cart::yzz][Cart::xyz][3]+=pma1*R_temp[Cart::zz][Cart::xyz][3]+wmp1*R_temp[Cart::zz][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::xz][4]
R_temp[Cart::yzz][Cart::yzz][3]+=pma1*R_temp[Cart::zz][Cart::yzz][3]+wmp1*R_temp[Cart::zz][Cart::yzz][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::zz][4]
R_temp[Cart::yzz][Cart::xxx][3]+=pma1*R_temp[Cart::zz][Cart::xxx][3]+wmp1*R_temp[Cart::zz][Cart::xxx][4]
R_temp[Cart::yzz][Cart::xxz][3]+=pma1*R_temp[Cart::zz][Cart::xxz][3]+wmp1*R_temp[Cart::zz][Cart::xxz][4]
R_temp[Cart::yzz][Cart::xzz][3]+=pma1*R_temp[Cart::zz][Cart::xzz][3]+wmp1*R_temp[Cart::zz][Cart::xzz][4]
R_temp[Cart::yzz][Cart::zzz][3]+=pma1*R_temp[Cart::zz][Cart::zzz][3]+wmp1*R_temp[Cart::zz][Cart::zzz][4]
R_temp[Cart::xxx][Cart::yyy][3]+=pma0*R_temp[Cart::xx][Cart::yyy][3]+wmp0*R_temp[Cart::xx][Cart::yyy][4]+1*rzeta*(R_temp[Cart::x][Cart::yyy][3]-gfak*R_temp[Cart::x][Cart::yyy][4])
R_temp[Cart::xxx][Cart::xyy][3]+=pma0*R_temp[Cart::xx][Cart::xyy][3]+wmp0*R_temp[Cart::xx][Cart::xyy][4]+1*rzeta*(R_temp[Cart::x][Cart::xyy][3]-gfak*R_temp[Cart::x][Cart::xyy][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::yy][4]
R_temp[Cart::xxx][Cart::yyz][3]+=pma0*R_temp[Cart::xx][Cart::yyz][3]+wmp0*R_temp[Cart::xx][Cart::yyz][4]+1*rzeta*(R_temp[Cart::x][Cart::yyz][3]-gfak*R_temp[Cart::x][Cart::yyz][4])
R_temp[Cart::xxx][Cart::xxy][3]+=pma0*R_temp[Cart::xx][Cart::xxy][3]+wmp0*R_temp[Cart::xx][Cart::xxy][4]+1*rzeta*(R_temp[Cart::x][Cart::xxy][3]-gfak*R_temp[Cart::x][Cart::xxy][4])+0.5/_decay*2*R_temp[Cart::xx][Cart::xy][4]
R_temp[Cart::xxx][Cart::xyz][3]+=pma0*R_temp[Cart::xx][Cart::xyz][3]+wmp0*R_temp[Cart::xx][Cart::xyz][4]+1*rzeta*(R_temp[Cart::x][Cart::xyz][3]-gfak*R_temp[Cart::x][Cart::xyz][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::yz][4]
R_temp[Cart::xxx][Cart::yzz][3]+=pma0*R_temp[Cart::xx][Cart::yzz][3]+wmp0*R_temp[Cart::xx][Cart::yzz][4]+1*rzeta*(R_temp[Cart::x][Cart::yzz][3]-gfak*R_temp[Cart::x][Cart::yzz][4])
R_temp[Cart::xxx][Cart::xxx][3]+=pma0*R_temp[Cart::xx][Cart::xxx][3]+wmp0*R_temp[Cart::xx][Cart::xxx][4]+1*rzeta*(R_temp[Cart::x][Cart::xxx][3]-gfak*R_temp[Cart::x][Cart::xxx][4])+0.5/_decay*3*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xxx][Cart::xxz][3]+=pma0*R_temp[Cart::xx][Cart::xxz][3]+wmp0*R_temp[Cart::xx][Cart::xxz][4]+1*rzeta*(R_temp[Cart::x][Cart::xxz][3]-gfak*R_temp[Cart::x][Cart::xxz][4])+0.5/_decay*2*R_temp[Cart::xx][Cart::xz][4]
R_temp[Cart::xxx][Cart::xzz][3]+=pma0*R_temp[Cart::xx][Cart::xzz][3]+wmp0*R_temp[Cart::xx][Cart::xzz][4]+1*rzeta*(R_temp[Cart::x][Cart::xzz][3]-gfak*R_temp[Cart::x][Cart::xzz][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::zz][4]
R_temp[Cart::xxx][Cart::zzz][3]+=pma0*R_temp[Cart::xx][Cart::zzz][3]+wmp0*R_temp[Cart::xx][Cart::zzz][4]+1*rzeta*(R_temp[Cart::x][Cart::zzz][3]-gfak*R_temp[Cart::x][Cart::zzz][4])
R_temp[Cart::xxz][Cart::yyy][3]+=pma2*R_temp[Cart::xx][Cart::yyy][3]+wmp2*R_temp[Cart::xx][Cart::yyy][4]
R_temp[Cart::xxz][Cart::xyy][3]+=pma2*R_temp[Cart::xx][Cart::xyy][3]+wmp2*R_temp[Cart::xx][Cart::xyy][4]
R_temp[Cart::xxz][Cart::yyz][3]+=pma2*R_temp[Cart::xx][Cart::yyz][3]+wmp2*R_temp[Cart::xx][Cart::yyz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::yy][4]
R_temp[Cart::xxz][Cart::xxy][3]+=pma2*R_temp[Cart::xx][Cart::xxy][3]+wmp2*R_temp[Cart::xx][Cart::xxy][4]
R_temp[Cart::xxz][Cart::xyz][3]+=pma2*R_temp[Cart::xx][Cart::xyz][3]+wmp2*R_temp[Cart::xx][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::xy][4]
R_temp[Cart::xxz][Cart::yzz][3]+=pma2*R_temp[Cart::xx][Cart::yzz][3]+wmp2*R_temp[Cart::xx][Cart::yzz][4]+0.5/_decay*2*R_temp[Cart::xx][Cart::yz][4]
R_temp[Cart::xxz][Cart::xxx][3]+=pma2*R_temp[Cart::xx][Cart::xxx][3]+wmp2*R_temp[Cart::xx][Cart::xxx][4]
R_temp[Cart::xxz][Cart::xxz][3]+=pma2*R_temp[Cart::xx][Cart::xxz][3]+wmp2*R_temp[Cart::xx][Cart::xxz][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::xx][4]
R_temp[Cart::xxz][Cart::xzz][3]+=pma2*R_temp[Cart::xx][Cart::xzz][3]+wmp2*R_temp[Cart::xx][Cart::xzz][4]+0.5/_decay*2*R_temp[Cart::xx][Cart::xz][4]
R_temp[Cart::xxz][Cart::zzz][3]+=pma2*R_temp[Cart::xx][Cart::zzz][3]+wmp2*R_temp[Cart::xx][Cart::zzz][4]+0.5/_decay*3*R_temp[Cart::xx][Cart::zz][4]
R_temp[Cart::xzz][Cart::yyy][3]+=pma0*R_temp[Cart::zz][Cart::yyy][3]+wmp0*R_temp[Cart::zz][Cart::yyy][4]
R_temp[Cart::xzz][Cart::xyy][3]+=pma0*R_temp[Cart::zz][Cart::xyy][3]+wmp0*R_temp[Cart::zz][Cart::xyy][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::yy][4]
R_temp[Cart::xzz][Cart::yyz][3]+=pma0*R_temp[Cart::zz][Cart::yyz][3]+wmp0*R_temp[Cart::zz][Cart::yyz][4]
R_temp[Cart::xzz][Cart::xxy][3]+=pma0*R_temp[Cart::zz][Cart::xxy][3]+wmp0*R_temp[Cart::zz][Cart::xxy][4]+0.5/_decay*2*R_temp[Cart::zz][Cart::xy][4]
R_temp[Cart::xzz][Cart::xyz][3]+=pma0*R_temp[Cart::zz][Cart::xyz][3]+wmp0*R_temp[Cart::zz][Cart::xyz][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::yz][4]
R_temp[Cart::xzz][Cart::yzz][3]+=pma0*R_temp[Cart::zz][Cart::yzz][3]+wmp0*R_temp[Cart::zz][Cart::yzz][4]
R_temp[Cart::xzz][Cart::xxx][3]+=pma0*R_temp[Cart::zz][Cart::xxx][3]+wmp0*R_temp[Cart::zz][Cart::xxx][4]+0.5/_decay*3*R_temp[Cart::zz][Cart::xx][4]
R_temp[Cart::xzz][Cart::xxz][3]+=pma0*R_temp[Cart::zz][Cart::xxz][3]+wmp0*R_temp[Cart::zz][Cart::xxz][4]+0.5/_decay*2*R_temp[Cart::zz][Cart::xz][4]
R_temp[Cart::xzz][Cart::xzz][3]+=pma0*R_temp[Cart::zz][Cart::xzz][3]+wmp0*R_temp[Cart::zz][Cart::xzz][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::zz][4]
R_temp[Cart::xzz][Cart::zzz][3]+=pma0*R_temp[Cart::zz][Cart::zzz][3]+wmp0*R_temp[Cart::zz][Cart::zzz][4]
R_temp[Cart::zzz][Cart::yyy][3]+=pma2*R_temp[Cart::zz][Cart::yyy][3]+wmp2*R_temp[Cart::zz][Cart::yyy][4]+1*rzeta*(R_temp[Cart::z][Cart::yyy][3]-gfak*R_temp[Cart::z][Cart::yyy][4])
R_temp[Cart::zzz][Cart::xyy][3]+=pma2*R_temp[Cart::zz][Cart::xyy][3]+wmp2*R_temp[Cart::zz][Cart::xyy][4]+1*rzeta*(R_temp[Cart::z][Cart::xyy][3]-gfak*R_temp[Cart::z][Cart::xyy][4])
R_temp[Cart::zzz][Cart::yyz][3]+=pma2*R_temp[Cart::zz][Cart::yyz][3]+wmp2*R_temp[Cart::zz][Cart::yyz][4]+1*rzeta*(R_temp[Cart::z][Cart::yyz][3]-gfak*R_temp[Cart::z][Cart::yyz][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::yy][4]
R_temp[Cart::zzz][Cart::xxy][3]+=pma2*R_temp[Cart::zz][Cart::xxy][3]+wmp2*R_temp[Cart::zz][Cart::xxy][4]+1*rzeta*(R_temp[Cart::z][Cart::xxy][3]-gfak*R_temp[Cart::z][Cart::xxy][4])
R_temp[Cart::zzz][Cart::xyz][3]+=pma2*R_temp[Cart::zz][Cart::xyz][3]+wmp2*R_temp[Cart::zz][Cart::xyz][4]+1*rzeta*(R_temp[Cart::z][Cart::xyz][3]-gfak*R_temp[Cart::z][Cart::xyz][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::xy][4]
R_temp[Cart::zzz][Cart::yzz][3]+=pma2*R_temp[Cart::zz][Cart::yzz][3]+wmp2*R_temp[Cart::zz][Cart::yzz][4]+1*rzeta*(R_temp[Cart::z][Cart::yzz][3]-gfak*R_temp[Cart::z][Cart::yzz][4])+0.5/_decay*2*R_temp[Cart::zz][Cart::yz][4]
R_temp[Cart::zzz][Cart::xxx][3]+=pma2*R_temp[Cart::zz][Cart::xxx][3]+wmp2*R_temp[Cart::zz][Cart::xxx][4]+1*rzeta*(R_temp[Cart::z][Cart::xxx][3]-gfak*R_temp[Cart::z][Cart::xxx][4])
R_temp[Cart::zzz][Cart::xxz][3]+=pma2*R_temp[Cart::zz][Cart::xxz][3]+wmp2*R_temp[Cart::zz][Cart::xxz][4]+1*rzeta*(R_temp[Cart::z][Cart::xxz][3]-gfak*R_temp[Cart::z][Cart::xxz][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::xx][4]
R_temp[Cart::zzz][Cart::xzz][3]+=pma2*R_temp[Cart::zz][Cart::xzz][3]+wmp2*R_temp[Cart::zz][Cart::xzz][4]+1*rzeta*(R_temp[Cart::z][Cart::xzz][3]-gfak*R_temp[Cart::z][Cart::xzz][4])+0.5/_decay*2*R_temp[Cart::zz][Cart::xz][4]
R_temp[Cart::zzz][Cart::zzz][3]+=pma2*R_temp[Cart::zz][Cart::zzz][3]+wmp2*R_temp[Cart::zz][Cart::zzz][4]+1*rzeta*(R_temp[Cart::z][Cart::zzz][3]-gfak*R_temp[Cart::z][Cart::zzz][4])+0.5/_decay*3*R_temp[Cart::zz][Cart::zz][4]
}}
//------------------------------------------------------

//Integral g - s - s - m2
if (_mmax >4 ){
if (_lmax_alpha>3){

R_temp[Cart::yyyy][Cart::s][2]+=pma1*R_temp[Cart::yyy][Cart::s][2]+wmp1*R_temp[Cart::yyy][Cart::s][3]+2*rzeta*(R_temp[Cart::yy][Cart::s][2]-gfak*R_temp[Cart::yy][Cart::s][3])
R_temp[Cart::xyyy][Cart::s][2]+=pma0*R_temp[Cart::yyy][Cart::s][2]+wmp0*R_temp[Cart::yyy][Cart::s][3]
R_temp[Cart::yyyz][Cart::s][2]+=pma2*R_temp[Cart::yyy][Cart::s][2]+wmp2*R_temp[Cart::yyy][Cart::s][3]
R_temp[Cart::xxyy][Cart::s][2]+=pma0*R_temp[Cart::xyy][Cart::s][2]+wmp0*R_temp[Cart::xyy][Cart::s][3]
R_temp[Cart::xyyz][Cart::s][2]+=pma0*R_temp[Cart::yyz][Cart::s][2]+wmp0*R_temp[Cart::yyz][Cart::s][3]
R_temp[Cart::yyzz][Cart::s][2]+=pma1*R_temp[Cart::yzz][Cart::s][2]+wmp1*R_temp[Cart::yzz][Cart::s][3]
R_temp[Cart::xxxy][Cart::s][2]+=pma1*R_temp[Cart::xxx][Cart::s][2]+wmp1*R_temp[Cart::xxx][Cart::s][3]
R_temp[Cart::xxyz][Cart::s][2]+=pma1*R_temp[Cart::xxz][Cart::s][2]+wmp1*R_temp[Cart::xxz][Cart::s][3]
R_temp[Cart::xyzz][Cart::s][2]+=pma0*R_temp[Cart::yzz][Cart::s][2]+wmp0*R_temp[Cart::yzz][Cart::s][3]
R_temp[Cart::yzzz][Cart::s][2]+=pma1*R_temp[Cart::zzz][Cart::s][2]+wmp1*R_temp[Cart::zzz][Cart::s][3]
R_temp[Cart::xxxx][Cart::s][2]+=pma0*R_temp[Cart::xxx][Cart::s][2]+wmp0*R_temp[Cart::xxx][Cart::s][3]+2*rzeta*(R_temp[Cart::xx][Cart::s][2]-gfak*R_temp[Cart::xx][Cart::s][3])
R_temp[Cart::xxxz][Cart::s][2]+=pma2*R_temp[Cart::xxx][Cart::s][2]+wmp2*R_temp[Cart::xxx][Cart::s][3]
R_temp[Cart::xxzz][Cart::s][2]+=pma0*R_temp[Cart::xzz][Cart::s][2]+wmp0*R_temp[Cart::xzz][Cart::s][3]
R_temp[Cart::xzzz][Cart::s][2]+=pma0*R_temp[Cart::zzz][Cart::s][2]+wmp0*R_temp[Cart::zzz][Cart::s][3]
R_temp[Cart::zzzz][Cart::s][2]+=pma2*R_temp[Cart::zzz][Cart::s][2]+wmp2*R_temp[Cart::zzz][Cart::s][3]+2*rzeta*(R_temp[Cart::zz][Cart::s][2]-gfak*R_temp[Cart::zz][Cart::s][3])
}}
//------------------------------------------------------

//Integral g - s - p - m2
if (_mmax >5 ){
if (_lmax_alpha>3 && _lmax_gamma>0){

R_temp[Cart::yyyy][Cart::y][2]+=pma1*R_temp[Cart::yyy][Cart::y][2]+wmp1*R_temp[Cart::yyy][Cart::y][3]+2*rzeta*(R_temp[Cart::yy][Cart::y][2]-gfak*R_temp[Cart::yy][Cart::y][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][3]
R_temp[Cart::yyyy][Cart::x][2]+=pma1*R_temp[Cart::yyy][Cart::x][2]+wmp1*R_temp[Cart::yyy][Cart::x][3]+2*rzeta*(R_temp[Cart::yy][Cart::x][2]-gfak*R_temp[Cart::yy][Cart::x][3])
R_temp[Cart::yyyy][Cart::z][2]+=pma1*R_temp[Cart::yyy][Cart::z][2]+wmp1*R_temp[Cart::yyy][Cart::z][3]+2*rzeta*(R_temp[Cart::yy][Cart::z][2]-gfak*R_temp[Cart::yy][Cart::z][3])
R_temp[Cart::xyyy][Cart::y][2]+=pma0*R_temp[Cart::yyy][Cart::y][2]+wmp0*R_temp[Cart::yyy][Cart::y][3]
R_temp[Cart::xyyy][Cart::x][2]+=pma0*R_temp[Cart::yyy][Cart::x][2]+wmp0*R_temp[Cart::yyy][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][3]
R_temp[Cart::xyyy][Cart::z][2]+=pma0*R_temp[Cart::yyy][Cart::z][2]+wmp0*R_temp[Cart::yyy][Cart::z][3]
R_temp[Cart::yyyz][Cart::y][2]+=pma2*R_temp[Cart::yyy][Cart::y][2]+wmp2*R_temp[Cart::yyy][Cart::y][3]
R_temp[Cart::yyyz][Cart::x][2]+=pma2*R_temp[Cart::yyy][Cart::x][2]+wmp2*R_temp[Cart::yyy][Cart::x][3]
R_temp[Cart::yyyz][Cart::z][2]+=pma2*R_temp[Cart::yyy][Cart::z][2]+wmp2*R_temp[Cart::yyy][Cart::z][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][3]
R_temp[Cart::xxyy][Cart::y][2]+=pma0*R_temp[Cart::xyy][Cart::y][2]+wmp0*R_temp[Cart::xyy][Cart::y][3]
R_temp[Cart::xxyy][Cart::x][2]+=pma0*R_temp[Cart::xyy][Cart::x][2]+wmp0*R_temp[Cart::xyy][Cart::x][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::s][3]
R_temp[Cart::xxyy][Cart::z][2]+=pma0*R_temp[Cart::xyy][Cart::z][2]+wmp0*R_temp[Cart::xyy][Cart::z][3]
R_temp[Cart::xyyz][Cart::y][2]+=pma0*R_temp[Cart::yyz][Cart::y][2]+wmp0*R_temp[Cart::yyz][Cart::y][3]
R_temp[Cart::xyyz][Cart::x][2]+=pma0*R_temp[Cart::yyz][Cart::x][2]+wmp0*R_temp[Cart::yyz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::s][3]
R_temp[Cart::xyyz][Cart::z][2]+=pma0*R_temp[Cart::yyz][Cart::z][2]+wmp0*R_temp[Cart::yyz][Cart::z][3]
R_temp[Cart::yyzz][Cart::y][2]+=pma1*R_temp[Cart::yzz][Cart::y][2]+wmp1*R_temp[Cart::yzz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][3]
R_temp[Cart::yyzz][Cart::x][2]+=pma1*R_temp[Cart::yzz][Cart::x][2]+wmp1*R_temp[Cart::yzz][Cart::x][3]
R_temp[Cart::yyzz][Cart::z][2]+=pma1*R_temp[Cart::yzz][Cart::z][2]+wmp1*R_temp[Cart::yzz][Cart::z][3]
R_temp[Cart::xxxy][Cart::y][2]+=pma1*R_temp[Cart::xxx][Cart::y][2]+wmp1*R_temp[Cart::xxx][Cart::y][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][3]
R_temp[Cart::xxxy][Cart::x][2]+=pma1*R_temp[Cart::xxx][Cart::x][2]+wmp1*R_temp[Cart::xxx][Cart::x][3]
R_temp[Cart::xxxy][Cart::z][2]+=pma1*R_temp[Cart::xxx][Cart::z][2]+wmp1*R_temp[Cart::xxx][Cart::z][3]
R_temp[Cart::xxyz][Cart::y][2]+=pma1*R_temp[Cart::xxz][Cart::y][2]+wmp1*R_temp[Cart::xxz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::s][3]
R_temp[Cart::xxyz][Cart::x][2]+=pma1*R_temp[Cart::xxz][Cart::x][2]+wmp1*R_temp[Cart::xxz][Cart::x][3]
R_temp[Cart::xxyz][Cart::z][2]+=pma1*R_temp[Cart::xxz][Cart::z][2]+wmp1*R_temp[Cart::xxz][Cart::z][3]
R_temp[Cart::xyzz][Cart::y][2]+=pma0*R_temp[Cart::yzz][Cart::y][2]+wmp0*R_temp[Cart::yzz][Cart::y][3]
R_temp[Cart::xyzz][Cart::x][2]+=pma0*R_temp[Cart::yzz][Cart::x][2]+wmp0*R_temp[Cart::yzz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][3]
R_temp[Cart::xyzz][Cart::z][2]+=pma0*R_temp[Cart::yzz][Cart::z][2]+wmp0*R_temp[Cart::yzz][Cart::z][3]
R_temp[Cart::yzzz][Cart::y][2]+=pma1*R_temp[Cart::zzz][Cart::y][2]+wmp1*R_temp[Cart::zzz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][3]
R_temp[Cart::yzzz][Cart::x][2]+=pma1*R_temp[Cart::zzz][Cart::x][2]+wmp1*R_temp[Cart::zzz][Cart::x][3]
R_temp[Cart::yzzz][Cart::z][2]+=pma1*R_temp[Cart::zzz][Cart::z][2]+wmp1*R_temp[Cart::zzz][Cart::z][3]
R_temp[Cart::xxxx][Cart::y][2]+=pma0*R_temp[Cart::xxx][Cart::y][2]+wmp0*R_temp[Cart::xxx][Cart::y][3]+2*rzeta*(R_temp[Cart::xx][Cart::y][2]-gfak*R_temp[Cart::xx][Cart::y][3])
R_temp[Cart::xxxx][Cart::x][2]+=pma0*R_temp[Cart::xxx][Cart::x][2]+wmp0*R_temp[Cart::xxx][Cart::x][3]+2*rzeta*(R_temp[Cart::xx][Cart::x][2]-gfak*R_temp[Cart::xx][Cart::x][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][3]
R_temp[Cart::xxxx][Cart::z][2]+=pma0*R_temp[Cart::xxx][Cart::z][2]+wmp0*R_temp[Cart::xxx][Cart::z][3]+2*rzeta*(R_temp[Cart::xx][Cart::z][2]-gfak*R_temp[Cart::xx][Cart::z][3])
R_temp[Cart::xxxz][Cart::y][2]+=pma2*R_temp[Cart::xxx][Cart::y][2]+wmp2*R_temp[Cart::xxx][Cart::y][3]
R_temp[Cart::xxxz][Cart::x][2]+=pma2*R_temp[Cart::xxx][Cart::x][2]+wmp2*R_temp[Cart::xxx][Cart::x][3]
R_temp[Cart::xxxz][Cart::z][2]+=pma2*R_temp[Cart::xxx][Cart::z][2]+wmp2*R_temp[Cart::xxx][Cart::z][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][3]
R_temp[Cart::xxzz][Cart::y][2]+=pma0*R_temp[Cart::xzz][Cart::y][2]+wmp0*R_temp[Cart::xzz][Cart::y][3]
R_temp[Cart::xxzz][Cart::x][2]+=pma0*R_temp[Cart::xzz][Cart::x][2]+wmp0*R_temp[Cart::xzz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::s][3]
R_temp[Cart::xxzz][Cart::z][2]+=pma0*R_temp[Cart::xzz][Cart::z][2]+wmp0*R_temp[Cart::xzz][Cart::z][3]
R_temp[Cart::xzzz][Cart::y][2]+=pma0*R_temp[Cart::zzz][Cart::y][2]+wmp0*R_temp[Cart::zzz][Cart::y][3]
R_temp[Cart::xzzz][Cart::x][2]+=pma0*R_temp[Cart::zzz][Cart::x][2]+wmp0*R_temp[Cart::zzz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][3]
R_temp[Cart::xzzz][Cart::z][2]+=pma0*R_temp[Cart::zzz][Cart::z][2]+wmp0*R_temp[Cart::zzz][Cart::z][3]
R_temp[Cart::zzzz][Cart::y][2]+=pma2*R_temp[Cart::zzz][Cart::y][2]+wmp2*R_temp[Cart::zzz][Cart::y][3]+2*rzeta*(R_temp[Cart::zz][Cart::y][2]-gfak*R_temp[Cart::zz][Cart::y][3])
R_temp[Cart::zzzz][Cart::x][2]+=pma2*R_temp[Cart::zzz][Cart::x][2]+wmp2*R_temp[Cart::zzz][Cart::x][3]+2*rzeta*(R_temp[Cart::zz][Cart::x][2]-gfak*R_temp[Cart::zz][Cart::x][3])
R_temp[Cart::zzzz][Cart::z][2]+=pma2*R_temp[Cart::zzz][Cart::z][2]+wmp2*R_temp[Cart::zzz][Cart::z][3]+2*rzeta*(R_temp[Cart::zz][Cart::z][2]-gfak*R_temp[Cart::zz][Cart::z][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][3]
}}
//------------------------------------------------------

//Integral g - s - d - m2
if (_mmax >6 ){
if (_lmax_alpha>3 && _lmax_gamma>1){

R_temp[Cart::yyyy][Cart::yy][2]+=pma1*R_temp[Cart::yyy][Cart::yy][2]+wmp1*R_temp[Cart::yyy][Cart::yy][3]+2*rzeta*(R_temp[Cart::yy][Cart::yy][2]-gfak*R_temp[Cart::yy][Cart::yy][3])+0.5/_decay*2*R_temp[Cart::yyy][Cart::y][3]
R_temp[Cart::yyyy][Cart::xy][2]+=pma1*R_temp[Cart::yyy][Cart::xy][2]+wmp1*R_temp[Cart::yyy][Cart::xy][3]+2*rzeta*(R_temp[Cart::yy][Cart::xy][2]-gfak*R_temp[Cart::yy][Cart::xy][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::x][3]
R_temp[Cart::yyyy][Cart::yz][2]+=pma1*R_temp[Cart::yyy][Cart::yz][2]+wmp1*R_temp[Cart::yyy][Cart::yz][3]+2*rzeta*(R_temp[Cart::yy][Cart::yz][2]-gfak*R_temp[Cart::yy][Cart::yz][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::z][3]
R_temp[Cart::yyyy][Cart::xx][2]+=pma1*R_temp[Cart::yyy][Cart::xx][2]+wmp1*R_temp[Cart::yyy][Cart::xx][3]+2*rzeta*(R_temp[Cart::yy][Cart::xx][2]-gfak*R_temp[Cart::yy][Cart::xx][3])
R_temp[Cart::yyyy][Cart::xz][2]+=pma1*R_temp[Cart::yyy][Cart::xz][2]+wmp1*R_temp[Cart::yyy][Cart::xz][3]+2*rzeta*(R_temp[Cart::yy][Cart::xz][2]-gfak*R_temp[Cart::yy][Cart::xz][3])
R_temp[Cart::yyyy][Cart::zz][2]+=pma1*R_temp[Cart::yyy][Cart::zz][2]+wmp1*R_temp[Cart::yyy][Cart::zz][3]+2*rzeta*(R_temp[Cart::yy][Cart::zz][2]-gfak*R_temp[Cart::yy][Cart::zz][3])
R_temp[Cart::xyyy][Cart::yy][2]+=pma0*R_temp[Cart::yyy][Cart::yy][2]+wmp0*R_temp[Cart::yyy][Cart::yy][3]
R_temp[Cart::xyyy][Cart::xy][2]+=pma0*R_temp[Cart::yyy][Cart::xy][2]+wmp0*R_temp[Cart::yyy][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::y][3]
R_temp[Cart::xyyy][Cart::yz][2]+=pma0*R_temp[Cart::yyy][Cart::yz][2]+wmp0*R_temp[Cart::yyy][Cart::yz][3]
R_temp[Cart::xyyy][Cart::xx][2]+=pma0*R_temp[Cart::yyy][Cart::xx][2]+wmp0*R_temp[Cart::yyy][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::yyy][Cart::x][3]
R_temp[Cart::xyyy][Cart::xz][2]+=pma0*R_temp[Cart::yyy][Cart::xz][2]+wmp0*R_temp[Cart::yyy][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::z][3]
R_temp[Cart::xyyy][Cart::zz][2]+=pma0*R_temp[Cart::yyy][Cart::zz][2]+wmp0*R_temp[Cart::yyy][Cart::zz][3]
R_temp[Cart::yyyz][Cart::yy][2]+=pma2*R_temp[Cart::yyy][Cart::yy][2]+wmp2*R_temp[Cart::yyy][Cart::yy][3]
R_temp[Cart::yyyz][Cart::xy][2]+=pma2*R_temp[Cart::yyy][Cart::xy][2]+wmp2*R_temp[Cart::yyy][Cart::xy][3]
R_temp[Cart::yyyz][Cart::yz][2]+=pma2*R_temp[Cart::yyy][Cart::yz][2]+wmp2*R_temp[Cart::yyy][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::y][3]
R_temp[Cart::yyyz][Cart::xx][2]+=pma2*R_temp[Cart::yyy][Cart::xx][2]+wmp2*R_temp[Cart::yyy][Cart::xx][3]
R_temp[Cart::yyyz][Cart::xz][2]+=pma2*R_temp[Cart::yyy][Cart::xz][2]+wmp2*R_temp[Cart::yyy][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::x][3]
R_temp[Cart::yyyz][Cart::zz][2]+=pma2*R_temp[Cart::yyy][Cart::zz][2]+wmp2*R_temp[Cart::yyy][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::yyy][Cart::z][3]
R_temp[Cart::xxyy][Cart::yy][2]+=pma0*R_temp[Cart::xyy][Cart::yy][2]+wmp0*R_temp[Cart::xyy][Cart::yy][3]
R_temp[Cart::xxyy][Cart::xy][2]+=pma0*R_temp[Cart::xyy][Cart::xy][2]+wmp0*R_temp[Cart::xyy][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::y][3]
R_temp[Cart::xxyy][Cart::yz][2]+=pma0*R_temp[Cart::xyy][Cart::yz][2]+wmp0*R_temp[Cart::xyy][Cart::yz][3]
R_temp[Cart::xxyy][Cart::xx][2]+=pma0*R_temp[Cart::xyy][Cart::xx][2]+wmp0*R_temp[Cart::xyy][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::xyy][Cart::x][3]
R_temp[Cart::xxyy][Cart::xz][2]+=pma0*R_temp[Cart::xyy][Cart::xz][2]+wmp0*R_temp[Cart::xyy][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::z][3]
R_temp[Cart::xxyy][Cart::zz][2]+=pma0*R_temp[Cart::xyy][Cart::zz][2]+wmp0*R_temp[Cart::xyy][Cart::zz][3]
R_temp[Cart::xyyz][Cart::yy][2]+=pma0*R_temp[Cart::yyz][Cart::yy][2]+wmp0*R_temp[Cart::yyz][Cart::yy][3]
R_temp[Cart::xyyz][Cart::xy][2]+=pma0*R_temp[Cart::yyz][Cart::xy][2]+wmp0*R_temp[Cart::yyz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::y][3]
R_temp[Cart::xyyz][Cart::yz][2]+=pma0*R_temp[Cart::yyz][Cart::yz][2]+wmp0*R_temp[Cart::yyz][Cart::yz][3]
R_temp[Cart::xyyz][Cart::xx][2]+=pma0*R_temp[Cart::yyz][Cart::xx][2]+wmp0*R_temp[Cart::yyz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::yyz][Cart::x][3]
R_temp[Cart::xyyz][Cart::xz][2]+=pma0*R_temp[Cart::yyz][Cart::xz][2]+wmp0*R_temp[Cart::yyz][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::z][3]
R_temp[Cart::xyyz][Cart::zz][2]+=pma0*R_temp[Cart::yyz][Cart::zz][2]+wmp0*R_temp[Cart::yyz][Cart::zz][3]
R_temp[Cart::yyzz][Cart::yy][2]+=pma1*R_temp[Cart::yzz][Cart::yy][2]+wmp1*R_temp[Cart::yzz][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::yzz][Cart::y][3]
R_temp[Cart::yyzz][Cart::xy][2]+=pma1*R_temp[Cart::yzz][Cart::xy][2]+wmp1*R_temp[Cart::yzz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::x][3]
R_temp[Cart::yyzz][Cart::yz][2]+=pma1*R_temp[Cart::yzz][Cart::yz][2]+wmp1*R_temp[Cart::yzz][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::z][3]
R_temp[Cart::yyzz][Cart::xx][2]+=pma1*R_temp[Cart::yzz][Cart::xx][2]+wmp1*R_temp[Cart::yzz][Cart::xx][3]
R_temp[Cart::yyzz][Cart::xz][2]+=pma1*R_temp[Cart::yzz][Cart::xz][2]+wmp1*R_temp[Cart::yzz][Cart::xz][3]
R_temp[Cart::yyzz][Cart::zz][2]+=pma1*R_temp[Cart::yzz][Cart::zz][2]+wmp1*R_temp[Cart::yzz][Cart::zz][3]
R_temp[Cart::xxxy][Cart::yy][2]+=pma1*R_temp[Cart::xxx][Cart::yy][2]+wmp1*R_temp[Cart::xxx][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::xxx][Cart::y][3]
R_temp[Cart::xxxy][Cart::xy][2]+=pma1*R_temp[Cart::xxx][Cart::xy][2]+wmp1*R_temp[Cart::xxx][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::x][3]
R_temp[Cart::xxxy][Cart::yz][2]+=pma1*R_temp[Cart::xxx][Cart::yz][2]+wmp1*R_temp[Cart::xxx][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::z][3]
R_temp[Cart::xxxy][Cart::xx][2]+=pma1*R_temp[Cart::xxx][Cart::xx][2]+wmp1*R_temp[Cart::xxx][Cart::xx][3]
R_temp[Cart::xxxy][Cart::xz][2]+=pma1*R_temp[Cart::xxx][Cart::xz][2]+wmp1*R_temp[Cart::xxx][Cart::xz][3]
R_temp[Cart::xxxy][Cart::zz][2]+=pma1*R_temp[Cart::xxx][Cart::zz][2]+wmp1*R_temp[Cart::xxx][Cart::zz][3]
R_temp[Cart::xxyz][Cart::yy][2]+=pma1*R_temp[Cart::xxz][Cart::yy][2]+wmp1*R_temp[Cart::xxz][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::xxz][Cart::y][3]
R_temp[Cart::xxyz][Cart::xy][2]+=pma1*R_temp[Cart::xxz][Cart::xy][2]+wmp1*R_temp[Cart::xxz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::x][3]
R_temp[Cart::xxyz][Cart::yz][2]+=pma1*R_temp[Cart::xxz][Cart::yz][2]+wmp1*R_temp[Cart::xxz][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::z][3]
R_temp[Cart::xxyz][Cart::xx][2]+=pma1*R_temp[Cart::xxz][Cart::xx][2]+wmp1*R_temp[Cart::xxz][Cart::xx][3]
R_temp[Cart::xxyz][Cart::xz][2]+=pma1*R_temp[Cart::xxz][Cart::xz][2]+wmp1*R_temp[Cart::xxz][Cart::xz][3]
R_temp[Cart::xxyz][Cart::zz][2]+=pma1*R_temp[Cart::xxz][Cart::zz][2]+wmp1*R_temp[Cart::xxz][Cart::zz][3]
R_temp[Cart::xyzz][Cart::yy][2]+=pma0*R_temp[Cart::yzz][Cart::yy][2]+wmp0*R_temp[Cart::yzz][Cart::yy][3]
R_temp[Cart::xyzz][Cart::xy][2]+=pma0*R_temp[Cart::yzz][Cart::xy][2]+wmp0*R_temp[Cart::yzz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::y][3]
R_temp[Cart::xyzz][Cart::yz][2]+=pma0*R_temp[Cart::yzz][Cart::yz][2]+wmp0*R_temp[Cart::yzz][Cart::yz][3]
R_temp[Cart::xyzz][Cart::xx][2]+=pma0*R_temp[Cart::yzz][Cart::xx][2]+wmp0*R_temp[Cart::yzz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::yzz][Cart::x][3]
R_temp[Cart::xyzz][Cart::xz][2]+=pma0*R_temp[Cart::yzz][Cart::xz][2]+wmp0*R_temp[Cart::yzz][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::z][3]
R_temp[Cart::xyzz][Cart::zz][2]+=pma0*R_temp[Cart::yzz][Cart::zz][2]+wmp0*R_temp[Cart::yzz][Cart::zz][3]
R_temp[Cart::yzzz][Cart::yy][2]+=pma1*R_temp[Cart::zzz][Cart::yy][2]+wmp1*R_temp[Cart::zzz][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::zzz][Cart::y][3]
R_temp[Cart::yzzz][Cart::xy][2]+=pma1*R_temp[Cart::zzz][Cart::xy][2]+wmp1*R_temp[Cart::zzz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::x][3]
R_temp[Cart::yzzz][Cart::yz][2]+=pma1*R_temp[Cart::zzz][Cart::yz][2]+wmp1*R_temp[Cart::zzz][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::z][3]
R_temp[Cart::yzzz][Cart::xx][2]+=pma1*R_temp[Cart::zzz][Cart::xx][2]+wmp1*R_temp[Cart::zzz][Cart::xx][3]
R_temp[Cart::yzzz][Cart::xz][2]+=pma1*R_temp[Cart::zzz][Cart::xz][2]+wmp1*R_temp[Cart::zzz][Cart::xz][3]
R_temp[Cart::yzzz][Cart::zz][2]+=pma1*R_temp[Cart::zzz][Cart::zz][2]+wmp1*R_temp[Cart::zzz][Cart::zz][3]
R_temp[Cart::xxxx][Cart::yy][2]+=pma0*R_temp[Cart::xxx][Cart::yy][2]+wmp0*R_temp[Cart::xxx][Cart::yy][3]+2*rzeta*(R_temp[Cart::xx][Cart::yy][2]-gfak*R_temp[Cart::xx][Cart::yy][3])
R_temp[Cart::xxxx][Cart::xy][2]+=pma0*R_temp[Cart::xxx][Cart::xy][2]+wmp0*R_temp[Cart::xxx][Cart::xy][3]+2*rzeta*(R_temp[Cart::xx][Cart::xy][2]-gfak*R_temp[Cart::xx][Cart::xy][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::y][3]
R_temp[Cart::xxxx][Cart::yz][2]+=pma0*R_temp[Cart::xxx][Cart::yz][2]+wmp0*R_temp[Cart::xxx][Cart::yz][3]+2*rzeta*(R_temp[Cart::xx][Cart::yz][2]-gfak*R_temp[Cart::xx][Cart::yz][3])
R_temp[Cart::xxxx][Cart::xx][2]+=pma0*R_temp[Cart::xxx][Cart::xx][2]+wmp0*R_temp[Cart::xxx][Cart::xx][3]+2*rzeta*(R_temp[Cart::xx][Cart::xx][2]-gfak*R_temp[Cart::xx][Cart::xx][3])+0.5/_decay*2*R_temp[Cart::xxx][Cart::x][3]
R_temp[Cart::xxxx][Cart::xz][2]+=pma0*R_temp[Cart::xxx][Cart::xz][2]+wmp0*R_temp[Cart::xxx][Cart::xz][3]+2*rzeta*(R_temp[Cart::xx][Cart::xz][2]-gfak*R_temp[Cart::xx][Cart::xz][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::z][3]
R_temp[Cart::xxxx][Cart::zz][2]+=pma0*R_temp[Cart::xxx][Cart::zz][2]+wmp0*R_temp[Cart::xxx][Cart::zz][3]+2*rzeta*(R_temp[Cart::xx][Cart::zz][2]-gfak*R_temp[Cart::xx][Cart::zz][3])
R_temp[Cart::xxxz][Cart::yy][2]+=pma2*R_temp[Cart::xxx][Cart::yy][2]+wmp2*R_temp[Cart::xxx][Cart::yy][3]
R_temp[Cart::xxxz][Cart::xy][2]+=pma2*R_temp[Cart::xxx][Cart::xy][2]+wmp2*R_temp[Cart::xxx][Cart::xy][3]
R_temp[Cart::xxxz][Cart::yz][2]+=pma2*R_temp[Cart::xxx][Cart::yz][2]+wmp2*R_temp[Cart::xxx][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::y][3]
R_temp[Cart::xxxz][Cart::xx][2]+=pma2*R_temp[Cart::xxx][Cart::xx][2]+wmp2*R_temp[Cart::xxx][Cart::xx][3]
R_temp[Cart::xxxz][Cart::xz][2]+=pma2*R_temp[Cart::xxx][Cart::xz][2]+wmp2*R_temp[Cart::xxx][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::x][3]
R_temp[Cart::xxxz][Cart::zz][2]+=pma2*R_temp[Cart::xxx][Cart::zz][2]+wmp2*R_temp[Cart::xxx][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::xxx][Cart::z][3]
R_temp[Cart::xxzz][Cart::yy][2]+=pma0*R_temp[Cart::xzz][Cart::yy][2]+wmp0*R_temp[Cart::xzz][Cart::yy][3]
R_temp[Cart::xxzz][Cart::xy][2]+=pma0*R_temp[Cart::xzz][Cart::xy][2]+wmp0*R_temp[Cart::xzz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::y][3]
R_temp[Cart::xxzz][Cart::yz][2]+=pma0*R_temp[Cart::xzz][Cart::yz][2]+wmp0*R_temp[Cart::xzz][Cart::yz][3]
R_temp[Cart::xxzz][Cart::xx][2]+=pma0*R_temp[Cart::xzz][Cart::xx][2]+wmp0*R_temp[Cart::xzz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::xzz][Cart::x][3]
R_temp[Cart::xxzz][Cart::xz][2]+=pma0*R_temp[Cart::xzz][Cart::xz][2]+wmp0*R_temp[Cart::xzz][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::z][3]
R_temp[Cart::xxzz][Cart::zz][2]+=pma0*R_temp[Cart::xzz][Cart::zz][2]+wmp0*R_temp[Cart::xzz][Cart::zz][3]
R_temp[Cart::xzzz][Cart::yy][2]+=pma0*R_temp[Cart::zzz][Cart::yy][2]+wmp0*R_temp[Cart::zzz][Cart::yy][3]
R_temp[Cart::xzzz][Cart::xy][2]+=pma0*R_temp[Cart::zzz][Cart::xy][2]+wmp0*R_temp[Cart::zzz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::y][3]
R_temp[Cart::xzzz][Cart::yz][2]+=pma0*R_temp[Cart::zzz][Cart::yz][2]+wmp0*R_temp[Cart::zzz][Cart::yz][3]
R_temp[Cart::xzzz][Cart::xx][2]+=pma0*R_temp[Cart::zzz][Cart::xx][2]+wmp0*R_temp[Cart::zzz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::zzz][Cart::x][3]
R_temp[Cart::xzzz][Cart::xz][2]+=pma0*R_temp[Cart::zzz][Cart::xz][2]+wmp0*R_temp[Cart::zzz][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::z][3]
R_temp[Cart::xzzz][Cart::zz][2]+=pma0*R_temp[Cart::zzz][Cart::zz][2]+wmp0*R_temp[Cart::zzz][Cart::zz][3]
R_temp[Cart::zzzz][Cart::yy][2]+=pma2*R_temp[Cart::zzz][Cart::yy][2]+wmp2*R_temp[Cart::zzz][Cart::yy][3]+2*rzeta*(R_temp[Cart::zz][Cart::yy][2]-gfak*R_temp[Cart::zz][Cart::yy][3])
R_temp[Cart::zzzz][Cart::xy][2]+=pma2*R_temp[Cart::zzz][Cart::xy][2]+wmp2*R_temp[Cart::zzz][Cart::xy][3]+2*rzeta*(R_temp[Cart::zz][Cart::xy][2]-gfak*R_temp[Cart::zz][Cart::xy][3])
R_temp[Cart::zzzz][Cart::yz][2]+=pma2*R_temp[Cart::zzz][Cart::yz][2]+wmp2*R_temp[Cart::zzz][Cart::yz][3]+2*rzeta*(R_temp[Cart::zz][Cart::yz][2]-gfak*R_temp[Cart::zz][Cart::yz][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::y][3]
R_temp[Cart::zzzz][Cart::xx][2]+=pma2*R_temp[Cart::zzz][Cart::xx][2]+wmp2*R_temp[Cart::zzz][Cart::xx][3]+2*rzeta*(R_temp[Cart::zz][Cart::xx][2]-gfak*R_temp[Cart::zz][Cart::xx][3])
R_temp[Cart::zzzz][Cart::xz][2]+=pma2*R_temp[Cart::zzz][Cart::xz][2]+wmp2*R_temp[Cart::zzz][Cart::xz][3]+2*rzeta*(R_temp[Cart::zz][Cart::xz][2]-gfak*R_temp[Cart::zz][Cart::xz][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::x][3]
R_temp[Cart::zzzz][Cart::zz][2]+=pma2*R_temp[Cart::zzz][Cart::zz][2]+wmp2*R_temp[Cart::zzz][Cart::zz][3]+2*rzeta*(R_temp[Cart::zz][Cart::zz][2]-gfak*R_temp[Cart::zz][Cart::zz][3])+0.5/_decay*2*R_temp[Cart::zzz][Cart::z][3]
}}
//------------------------------------------------------

//Integral g - s - f - m2
if (_mmax >7 ){
if (_lmax_alpha>3 && _lmax_gamma>2){

R_temp[Cart::yyyy][Cart::yyy][2]+=pma1*R_temp[Cart::yyy][Cart::yyy][2]+wmp1*R_temp[Cart::yyy][Cart::yyy][3]+2*rzeta*(R_temp[Cart::yy][Cart::yyy][2]-gfak*R_temp[Cart::yy][Cart::yyy][3])+0.5/_decay*3*R_temp[Cart::yyy][Cart::yy][3]
R_temp[Cart::yyyy][Cart::xyy][2]+=pma1*R_temp[Cart::yyy][Cart::xyy][2]+wmp1*R_temp[Cart::yyy][Cart::xyy][3]+2*rzeta*(R_temp[Cart::yy][Cart::xyy][2]-gfak*R_temp[Cart::yy][Cart::xyy][3])+0.5/_decay*2*R_temp[Cart::yyy][Cart::xy][3]
R_temp[Cart::yyyy][Cart::yyz][2]+=pma1*R_temp[Cart::yyy][Cart::yyz][2]+wmp1*R_temp[Cart::yyy][Cart::yyz][3]+2*rzeta*(R_temp[Cart::yy][Cart::yyz][2]-gfak*R_temp[Cart::yy][Cart::yyz][3])+0.5/_decay*2*R_temp[Cart::yyy][Cart::yz][3]
R_temp[Cart::yyyy][Cart::xxy][2]+=pma1*R_temp[Cart::yyy][Cart::xxy][2]+wmp1*R_temp[Cart::yyy][Cart::xxy][3]+2*rzeta*(R_temp[Cart::yy][Cart::xxy][2]-gfak*R_temp[Cart::yy][Cart::xxy][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::xx][3]
R_temp[Cart::yyyy][Cart::xyz][2]+=pma1*R_temp[Cart::yyy][Cart::xyz][2]+wmp1*R_temp[Cart::yyy][Cart::xyz][3]+2*rzeta*(R_temp[Cart::yy][Cart::xyz][2]-gfak*R_temp[Cart::yy][Cart::xyz][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::xz][3]
R_temp[Cart::yyyy][Cart::yzz][2]+=pma1*R_temp[Cart::yyy][Cart::yzz][2]+wmp1*R_temp[Cart::yyy][Cart::yzz][3]+2*rzeta*(R_temp[Cart::yy][Cart::yzz][2]-gfak*R_temp[Cart::yy][Cart::yzz][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::zz][3]
R_temp[Cart::yyyy][Cart::xxx][2]+=pma1*R_temp[Cart::yyy][Cart::xxx][2]+wmp1*R_temp[Cart::yyy][Cart::xxx][3]+2*rzeta*(R_temp[Cart::yy][Cart::xxx][2]-gfak*R_temp[Cart::yy][Cart::xxx][3])
R_temp[Cart::yyyy][Cart::xxz][2]+=pma1*R_temp[Cart::yyy][Cart::xxz][2]+wmp1*R_temp[Cart::yyy][Cart::xxz][3]+2*rzeta*(R_temp[Cart::yy][Cart::xxz][2]-gfak*R_temp[Cart::yy][Cart::xxz][3])
R_temp[Cart::yyyy][Cart::xzz][2]+=pma1*R_temp[Cart::yyy][Cart::xzz][2]+wmp1*R_temp[Cart::yyy][Cart::xzz][3]+2*rzeta*(R_temp[Cart::yy][Cart::xzz][2]-gfak*R_temp[Cart::yy][Cart::xzz][3])
R_temp[Cart::yyyy][Cart::zzz][2]+=pma1*R_temp[Cart::yyy][Cart::zzz][2]+wmp1*R_temp[Cart::yyy][Cart::zzz][3]+2*rzeta*(R_temp[Cart::yy][Cart::zzz][2]-gfak*R_temp[Cart::yy][Cart::zzz][3])
R_temp[Cart::xyyy][Cart::yyy][2]+=pma0*R_temp[Cart::yyy][Cart::yyy][2]+wmp0*R_temp[Cart::yyy][Cart::yyy][3]
R_temp[Cart::xyyy][Cart::xyy][2]+=pma0*R_temp[Cart::yyy][Cart::xyy][2]+wmp0*R_temp[Cart::yyy][Cart::xyy][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::yy][3]
R_temp[Cart::xyyy][Cart::yyz][2]+=pma0*R_temp[Cart::yyy][Cart::yyz][2]+wmp0*R_temp[Cart::yyy][Cart::yyz][3]
R_temp[Cart::xyyy][Cart::xxy][2]+=pma0*R_temp[Cart::yyy][Cart::xxy][2]+wmp0*R_temp[Cart::yyy][Cart::xxy][3]+0.5/_decay*2*R_temp[Cart::yyy][Cart::xy][3]
R_temp[Cart::xyyy][Cart::xyz][2]+=pma0*R_temp[Cart::yyy][Cart::xyz][2]+wmp0*R_temp[Cart::yyy][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::yz][3]
R_temp[Cart::xyyy][Cart::yzz][2]+=pma0*R_temp[Cart::yyy][Cart::yzz][2]+wmp0*R_temp[Cart::yyy][Cart::yzz][3]
R_temp[Cart::xyyy][Cart::xxx][2]+=pma0*R_temp[Cart::yyy][Cart::xxx][2]+wmp0*R_temp[Cart::yyy][Cart::xxx][3]+0.5/_decay*3*R_temp[Cart::yyy][Cart::xx][3]
R_temp[Cart::xyyy][Cart::xxz][2]+=pma0*R_temp[Cart::yyy][Cart::xxz][2]+wmp0*R_temp[Cart::yyy][Cart::xxz][3]+0.5/_decay*2*R_temp[Cart::yyy][Cart::xz][3]
R_temp[Cart::xyyy][Cart::xzz][2]+=pma0*R_temp[Cart::yyy][Cart::xzz][2]+wmp0*R_temp[Cart::yyy][Cart::xzz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::zz][3]
R_temp[Cart::xyyy][Cart::zzz][2]+=pma0*R_temp[Cart::yyy][Cart::zzz][2]+wmp0*R_temp[Cart::yyy][Cart::zzz][3]
R_temp[Cart::yyyz][Cart::yyy][2]+=pma2*R_temp[Cart::yyy][Cart::yyy][2]+wmp2*R_temp[Cart::yyy][Cart::yyy][3]
R_temp[Cart::yyyz][Cart::xyy][2]+=pma2*R_temp[Cart::yyy][Cart::xyy][2]+wmp2*R_temp[Cart::yyy][Cart::xyy][3]
R_temp[Cart::yyyz][Cart::yyz][2]+=pma2*R_temp[Cart::yyy][Cart::yyz][2]+wmp2*R_temp[Cart::yyy][Cart::yyz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::yy][3]
R_temp[Cart::yyyz][Cart::xxy][2]+=pma2*R_temp[Cart::yyy][Cart::xxy][2]+wmp2*R_temp[Cart::yyy][Cart::xxy][3]
R_temp[Cart::yyyz][Cart::xyz][2]+=pma2*R_temp[Cart::yyy][Cart::xyz][2]+wmp2*R_temp[Cart::yyy][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::xy][3]
R_temp[Cart::yyyz][Cart::yzz][2]+=pma2*R_temp[Cart::yyy][Cart::yzz][2]+wmp2*R_temp[Cart::yyy][Cart::yzz][3]+0.5/_decay*2*R_temp[Cart::yyy][Cart::yz][3]
R_temp[Cart::yyyz][Cart::xxx][2]+=pma2*R_temp[Cart::yyy][Cart::xxx][2]+wmp2*R_temp[Cart::yyy][Cart::xxx][3]
R_temp[Cart::yyyz][Cart::xxz][2]+=pma2*R_temp[Cart::yyy][Cart::xxz][2]+wmp2*R_temp[Cart::yyy][Cart::xxz][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::xx][3]
R_temp[Cart::yyyz][Cart::xzz][2]+=pma2*R_temp[Cart::yyy][Cart::xzz][2]+wmp2*R_temp[Cart::yyy][Cart::xzz][3]+0.5/_decay*2*R_temp[Cart::yyy][Cart::xz][3]
R_temp[Cart::yyyz][Cart::zzz][2]+=pma2*R_temp[Cart::yyy][Cart::zzz][2]+wmp2*R_temp[Cart::yyy][Cart::zzz][3]+0.5/_decay*3*R_temp[Cart::yyy][Cart::zz][3]
R_temp[Cart::xxyy][Cart::yyy][2]+=pma0*R_temp[Cart::xyy][Cart::yyy][2]+wmp0*R_temp[Cart::xyy][Cart::yyy][3]
R_temp[Cart::xxyy][Cart::xyy][2]+=pma0*R_temp[Cart::xyy][Cart::xyy][2]+wmp0*R_temp[Cart::xyy][Cart::xyy][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::yy][3]
R_temp[Cart::xxyy][Cart::yyz][2]+=pma0*R_temp[Cart::xyy][Cart::yyz][2]+wmp0*R_temp[Cart::xyy][Cart::yyz][3]
R_temp[Cart::xxyy][Cart::xxy][2]+=pma0*R_temp[Cart::xyy][Cart::xxy][2]+wmp0*R_temp[Cart::xyy][Cart::xxy][3]+0.5/_decay*2*R_temp[Cart::xyy][Cart::xy][3]
R_temp[Cart::xxyy][Cart::xyz][2]+=pma0*R_temp[Cart::xyy][Cart::xyz][2]+wmp0*R_temp[Cart::xyy][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::yz][3]
R_temp[Cart::xxyy][Cart::yzz][2]+=pma0*R_temp[Cart::xyy][Cart::yzz][2]+wmp0*R_temp[Cart::xyy][Cart::yzz][3]
R_temp[Cart::xxyy][Cart::xxx][2]+=pma0*R_temp[Cart::xyy][Cart::xxx][2]+wmp0*R_temp[Cart::xyy][Cart::xxx][3]+0.5/_decay*3*R_temp[Cart::xyy][Cart::xx][3]
R_temp[Cart::xxyy][Cart::xxz][2]+=pma0*R_temp[Cart::xyy][Cart::xxz][2]+wmp0*R_temp[Cart::xyy][Cart::xxz][3]+0.5/_decay*2*R_temp[Cart::xyy][Cart::xz][3]
R_temp[Cart::xxyy][Cart::xzz][2]+=pma0*R_temp[Cart::xyy][Cart::xzz][2]+wmp0*R_temp[Cart::xyy][Cart::xzz][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::zz][3]
R_temp[Cart::xxyy][Cart::zzz][2]+=pma0*R_temp[Cart::xyy][Cart::zzz][2]+wmp0*R_temp[Cart::xyy][Cart::zzz][3]
R_temp[Cart::xyyz][Cart::yyy][2]+=pma0*R_temp[Cart::yyz][Cart::yyy][2]+wmp0*R_temp[Cart::yyz][Cart::yyy][3]
R_temp[Cart::xyyz][Cart::xyy][2]+=pma0*R_temp[Cart::yyz][Cart::xyy][2]+wmp0*R_temp[Cart::yyz][Cart::xyy][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::yy][3]
R_temp[Cart::xyyz][Cart::yyz][2]+=pma0*R_temp[Cart::yyz][Cart::yyz][2]+wmp0*R_temp[Cart::yyz][Cart::yyz][3]
R_temp[Cart::xyyz][Cart::xxy][2]+=pma0*R_temp[Cart::yyz][Cart::xxy][2]+wmp0*R_temp[Cart::yyz][Cart::xxy][3]+0.5/_decay*2*R_temp[Cart::yyz][Cart::xy][3]
R_temp[Cart::xyyz][Cart::xyz][2]+=pma0*R_temp[Cart::yyz][Cart::xyz][2]+wmp0*R_temp[Cart::yyz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::yz][3]
R_temp[Cart::xyyz][Cart::yzz][2]+=pma0*R_temp[Cart::yyz][Cart::yzz][2]+wmp0*R_temp[Cart::yyz][Cart::yzz][3]
R_temp[Cart::xyyz][Cart::xxx][2]+=pma0*R_temp[Cart::yyz][Cart::xxx][2]+wmp0*R_temp[Cart::yyz][Cart::xxx][3]+0.5/_decay*3*R_temp[Cart::yyz][Cart::xx][3]
R_temp[Cart::xyyz][Cart::xxz][2]+=pma0*R_temp[Cart::yyz][Cart::xxz][2]+wmp0*R_temp[Cart::yyz][Cart::xxz][3]+0.5/_decay*2*R_temp[Cart::yyz][Cart::xz][3]
R_temp[Cart::xyyz][Cart::xzz][2]+=pma0*R_temp[Cart::yyz][Cart::xzz][2]+wmp0*R_temp[Cart::yyz][Cart::xzz][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::zz][3]
R_temp[Cart::xyyz][Cart::zzz][2]+=pma0*R_temp[Cart::yyz][Cart::zzz][2]+wmp0*R_temp[Cart::yyz][Cart::zzz][3]
R_temp[Cart::yyzz][Cart::yyy][2]+=pma1*R_temp[Cart::yzz][Cart::yyy][2]+wmp1*R_temp[Cart::yzz][Cart::yyy][3]+0.5/_decay*3*R_temp[Cart::yzz][Cart::yy][3]
R_temp[Cart::yyzz][Cart::xyy][2]+=pma1*R_temp[Cart::yzz][Cart::xyy][2]+wmp1*R_temp[Cart::yzz][Cart::xyy][3]+0.5/_decay*2*R_temp[Cart::yzz][Cart::xy][3]
R_temp[Cart::yyzz][Cart::yyz][2]+=pma1*R_temp[Cart::yzz][Cart::yyz][2]+wmp1*R_temp[Cart::yzz][Cart::yyz][3]+0.5/_decay*2*R_temp[Cart::yzz][Cart::yz][3]
R_temp[Cart::yyzz][Cart::xxy][2]+=pma1*R_temp[Cart::yzz][Cart::xxy][2]+wmp1*R_temp[Cart::yzz][Cart::xxy][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::xx][3]
R_temp[Cart::yyzz][Cart::xyz][2]+=pma1*R_temp[Cart::yzz][Cart::xyz][2]+wmp1*R_temp[Cart::yzz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::xz][3]
R_temp[Cart::yyzz][Cart::yzz][2]+=pma1*R_temp[Cart::yzz][Cart::yzz][2]+wmp1*R_temp[Cart::yzz][Cart::yzz][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::zz][3]
R_temp[Cart::yyzz][Cart::xxx][2]+=pma1*R_temp[Cart::yzz][Cart::xxx][2]+wmp1*R_temp[Cart::yzz][Cart::xxx][3]
R_temp[Cart::yyzz][Cart::xxz][2]+=pma1*R_temp[Cart::yzz][Cart::xxz][2]+wmp1*R_temp[Cart::yzz][Cart::xxz][3]
R_temp[Cart::yyzz][Cart::xzz][2]+=pma1*R_temp[Cart::yzz][Cart::xzz][2]+wmp1*R_temp[Cart::yzz][Cart::xzz][3]
R_temp[Cart::yyzz][Cart::zzz][2]+=pma1*R_temp[Cart::yzz][Cart::zzz][2]+wmp1*R_temp[Cart::yzz][Cart::zzz][3]
R_temp[Cart::xxxy][Cart::yyy][2]+=pma1*R_temp[Cart::xxx][Cart::yyy][2]+wmp1*R_temp[Cart::xxx][Cart::yyy][3]+0.5/_decay*3*R_temp[Cart::xxx][Cart::yy][3]
R_temp[Cart::xxxy][Cart::xyy][2]+=pma1*R_temp[Cart::xxx][Cart::xyy][2]+wmp1*R_temp[Cart::xxx][Cart::xyy][3]+0.5/_decay*2*R_temp[Cart::xxx][Cart::xy][3]
R_temp[Cart::xxxy][Cart::yyz][2]+=pma1*R_temp[Cart::xxx][Cart::yyz][2]+wmp1*R_temp[Cart::xxx][Cart::yyz][3]+0.5/_decay*2*R_temp[Cart::xxx][Cart::yz][3]
R_temp[Cart::xxxy][Cart::xxy][2]+=pma1*R_temp[Cart::xxx][Cart::xxy][2]+wmp1*R_temp[Cart::xxx][Cart::xxy][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xx][3]
R_temp[Cart::xxxy][Cart::xyz][2]+=pma1*R_temp[Cart::xxx][Cart::xyz][2]+wmp1*R_temp[Cart::xxx][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xz][3]
R_temp[Cart::xxxy][Cart::yzz][2]+=pma1*R_temp[Cart::xxx][Cart::yzz][2]+wmp1*R_temp[Cart::xxx][Cart::yzz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::zz][3]
R_temp[Cart::xxxy][Cart::xxx][2]+=pma1*R_temp[Cart::xxx][Cart::xxx][2]+wmp1*R_temp[Cart::xxx][Cart::xxx][3]
R_temp[Cart::xxxy][Cart::xxz][2]+=pma1*R_temp[Cart::xxx][Cart::xxz][2]+wmp1*R_temp[Cart::xxx][Cart::xxz][3]
R_temp[Cart::xxxy][Cart::xzz][2]+=pma1*R_temp[Cart::xxx][Cart::xzz][2]+wmp1*R_temp[Cart::xxx][Cart::xzz][3]
R_temp[Cart::xxxy][Cart::zzz][2]+=pma1*R_temp[Cart::xxx][Cart::zzz][2]+wmp1*R_temp[Cart::xxx][Cart::zzz][3]
R_temp[Cart::xxyz][Cart::yyy][2]+=pma1*R_temp[Cart::xxz][Cart::yyy][2]+wmp1*R_temp[Cart::xxz][Cart::yyy][3]+0.5/_decay*3*R_temp[Cart::xxz][Cart::yy][3]
R_temp[Cart::xxyz][Cart::xyy][2]+=pma1*R_temp[Cart::xxz][Cart::xyy][2]+wmp1*R_temp[Cart::xxz][Cart::xyy][3]+0.5/_decay*2*R_temp[Cart::xxz][Cart::xy][3]
R_temp[Cart::xxyz][Cart::yyz][2]+=pma1*R_temp[Cart::xxz][Cart::yyz][2]+wmp1*R_temp[Cart::xxz][Cart::yyz][3]+0.5/_decay*2*R_temp[Cart::xxz][Cart::yz][3]
R_temp[Cart::xxyz][Cart::xxy][2]+=pma1*R_temp[Cart::xxz][Cart::xxy][2]+wmp1*R_temp[Cart::xxz][Cart::xxy][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::xx][3]
R_temp[Cart::xxyz][Cart::xyz][2]+=pma1*R_temp[Cart::xxz][Cart::xyz][2]+wmp1*R_temp[Cart::xxz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::xz][3]
R_temp[Cart::xxyz][Cart::yzz][2]+=pma1*R_temp[Cart::xxz][Cart::yzz][2]+wmp1*R_temp[Cart::xxz][Cart::yzz][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::zz][3]
R_temp[Cart::xxyz][Cart::xxx][2]+=pma1*R_temp[Cart::xxz][Cart::xxx][2]+wmp1*R_temp[Cart::xxz][Cart::xxx][3]
R_temp[Cart::xxyz][Cart::xxz][2]+=pma1*R_temp[Cart::xxz][Cart::xxz][2]+wmp1*R_temp[Cart::xxz][Cart::xxz][3]
R_temp[Cart::xxyz][Cart::xzz][2]+=pma1*R_temp[Cart::xxz][Cart::xzz][2]+wmp1*R_temp[Cart::xxz][Cart::xzz][3]
R_temp[Cart::xxyz][Cart::zzz][2]+=pma1*R_temp[Cart::xxz][Cart::zzz][2]+wmp1*R_temp[Cart::xxz][Cart::zzz][3]
R_temp[Cart::xyzz][Cart::yyy][2]+=pma0*R_temp[Cart::yzz][Cart::yyy][2]+wmp0*R_temp[Cart::yzz][Cart::yyy][3]
R_temp[Cart::xyzz][Cart::xyy][2]+=pma0*R_temp[Cart::yzz][Cart::xyy][2]+wmp0*R_temp[Cart::yzz][Cart::xyy][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::yy][3]
R_temp[Cart::xyzz][Cart::yyz][2]+=pma0*R_temp[Cart::yzz][Cart::yyz][2]+wmp0*R_temp[Cart::yzz][Cart::yyz][3]
R_temp[Cart::xyzz][Cart::xxy][2]+=pma0*R_temp[Cart::yzz][Cart::xxy][2]+wmp0*R_temp[Cart::yzz][Cart::xxy][3]+0.5/_decay*2*R_temp[Cart::yzz][Cart::xy][3]
R_temp[Cart::xyzz][Cart::xyz][2]+=pma0*R_temp[Cart::yzz][Cart::xyz][2]+wmp0*R_temp[Cart::yzz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::yz][3]
R_temp[Cart::xyzz][Cart::yzz][2]+=pma0*R_temp[Cart::yzz][Cart::yzz][2]+wmp0*R_temp[Cart::yzz][Cart::yzz][3]
R_temp[Cart::xyzz][Cart::xxx][2]+=pma0*R_temp[Cart::yzz][Cart::xxx][2]+wmp0*R_temp[Cart::yzz][Cart::xxx][3]+0.5/_decay*3*R_temp[Cart::yzz][Cart::xx][3]
R_temp[Cart::xyzz][Cart::xxz][2]+=pma0*R_temp[Cart::yzz][Cart::xxz][2]+wmp0*R_temp[Cart::yzz][Cart::xxz][3]+0.5/_decay*2*R_temp[Cart::yzz][Cart::xz][3]
R_temp[Cart::xyzz][Cart::xzz][2]+=pma0*R_temp[Cart::yzz][Cart::xzz][2]+wmp0*R_temp[Cart::yzz][Cart::xzz][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::zz][3]
R_temp[Cart::xyzz][Cart::zzz][2]+=pma0*R_temp[Cart::yzz][Cart::zzz][2]+wmp0*R_temp[Cart::yzz][Cart::zzz][3]
R_temp[Cart::yzzz][Cart::yyy][2]+=pma1*R_temp[Cart::zzz][Cart::yyy][2]+wmp1*R_temp[Cart::zzz][Cart::yyy][3]+0.5/_decay*3*R_temp[Cart::zzz][Cart::yy][3]
R_temp[Cart::yzzz][Cart::xyy][2]+=pma1*R_temp[Cart::zzz][Cart::xyy][2]+wmp1*R_temp[Cart::zzz][Cart::xyy][3]+0.5/_decay*2*R_temp[Cart::zzz][Cart::xy][3]
R_temp[Cart::yzzz][Cart::yyz][2]+=pma1*R_temp[Cart::zzz][Cart::yyz][2]+wmp1*R_temp[Cart::zzz][Cart::yyz][3]+0.5/_decay*2*R_temp[Cart::zzz][Cart::yz][3]
R_temp[Cart::yzzz][Cart::xxy][2]+=pma1*R_temp[Cart::zzz][Cart::xxy][2]+wmp1*R_temp[Cart::zzz][Cart::xxy][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::xx][3]
R_temp[Cart::yzzz][Cart::xyz][2]+=pma1*R_temp[Cart::zzz][Cart::xyz][2]+wmp1*R_temp[Cart::zzz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::xz][3]
R_temp[Cart::yzzz][Cart::yzz][2]+=pma1*R_temp[Cart::zzz][Cart::yzz][2]+wmp1*R_temp[Cart::zzz][Cart::yzz][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::zz][3]
R_temp[Cart::yzzz][Cart::xxx][2]+=pma1*R_temp[Cart::zzz][Cart::xxx][2]+wmp1*R_temp[Cart::zzz][Cart::xxx][3]
R_temp[Cart::yzzz][Cart::xxz][2]+=pma1*R_temp[Cart::zzz][Cart::xxz][2]+wmp1*R_temp[Cart::zzz][Cart::xxz][3]
R_temp[Cart::yzzz][Cart::xzz][2]+=pma1*R_temp[Cart::zzz][Cart::xzz][2]+wmp1*R_temp[Cart::zzz][Cart::xzz][3]
R_temp[Cart::yzzz][Cart::zzz][2]+=pma1*R_temp[Cart::zzz][Cart::zzz][2]+wmp1*R_temp[Cart::zzz][Cart::zzz][3]
R_temp[Cart::xxxx][Cart::yyy][2]+=pma0*R_temp[Cart::xxx][Cart::yyy][2]+wmp0*R_temp[Cart::xxx][Cart::yyy][3]+2*rzeta*(R_temp[Cart::xx][Cart::yyy][2]-gfak*R_temp[Cart::xx][Cart::yyy][3])
R_temp[Cart::xxxx][Cart::xyy][2]+=pma0*R_temp[Cart::xxx][Cart::xyy][2]+wmp0*R_temp[Cart::xxx][Cart::xyy][3]+2*rzeta*(R_temp[Cart::xx][Cart::xyy][2]-gfak*R_temp[Cart::xx][Cart::xyy][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::yy][3]
R_temp[Cart::xxxx][Cart::yyz][2]+=pma0*R_temp[Cart::xxx][Cart::yyz][2]+wmp0*R_temp[Cart::xxx][Cart::yyz][3]+2*rzeta*(R_temp[Cart::xx][Cart::yyz][2]-gfak*R_temp[Cart::xx][Cart::yyz][3])
R_temp[Cart::xxxx][Cart::xxy][2]+=pma0*R_temp[Cart::xxx][Cart::xxy][2]+wmp0*R_temp[Cart::xxx][Cart::xxy][3]+2*rzeta*(R_temp[Cart::xx][Cart::xxy][2]-gfak*R_temp[Cart::xx][Cart::xxy][3])+0.5/_decay*2*R_temp[Cart::xxx][Cart::xy][3]
R_temp[Cart::xxxx][Cart::xyz][2]+=pma0*R_temp[Cart::xxx][Cart::xyz][2]+wmp0*R_temp[Cart::xxx][Cart::xyz][3]+2*rzeta*(R_temp[Cart::xx][Cart::xyz][2]-gfak*R_temp[Cart::xx][Cart::xyz][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::yz][3]
R_temp[Cart::xxxx][Cart::yzz][2]+=pma0*R_temp[Cart::xxx][Cart::yzz][2]+wmp0*R_temp[Cart::xxx][Cart::yzz][3]+2*rzeta*(R_temp[Cart::xx][Cart::yzz][2]-gfak*R_temp[Cart::xx][Cart::yzz][3])
R_temp[Cart::xxxx][Cart::xxx][2]+=pma0*R_temp[Cart::xxx][Cart::xxx][2]+wmp0*R_temp[Cart::xxx][Cart::xxx][3]+2*rzeta*(R_temp[Cart::xx][Cart::xxx][2]-gfak*R_temp[Cart::xx][Cart::xxx][3])+0.5/_decay*3*R_temp[Cart::xxx][Cart::xx][3]
R_temp[Cart::xxxx][Cart::xxz][2]+=pma0*R_temp[Cart::xxx][Cart::xxz][2]+wmp0*R_temp[Cart::xxx][Cart::xxz][3]+2*rzeta*(R_temp[Cart::xx][Cart::xxz][2]-gfak*R_temp[Cart::xx][Cart::xxz][3])+0.5/_decay*2*R_temp[Cart::xxx][Cart::xz][3]
R_temp[Cart::xxxx][Cart::xzz][2]+=pma0*R_temp[Cart::xxx][Cart::xzz][2]+wmp0*R_temp[Cart::xxx][Cart::xzz][3]+2*rzeta*(R_temp[Cart::xx][Cart::xzz][2]-gfak*R_temp[Cart::xx][Cart::xzz][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::zz][3]
R_temp[Cart::xxxx][Cart::zzz][2]+=pma0*R_temp[Cart::xxx][Cart::zzz][2]+wmp0*R_temp[Cart::xxx][Cart::zzz][3]+2*rzeta*(R_temp[Cart::xx][Cart::zzz][2]-gfak*R_temp[Cart::xx][Cart::zzz][3])
R_temp[Cart::xxxz][Cart::yyy][2]+=pma2*R_temp[Cart::xxx][Cart::yyy][2]+wmp2*R_temp[Cart::xxx][Cart::yyy][3]
R_temp[Cart::xxxz][Cart::xyy][2]+=pma2*R_temp[Cart::xxx][Cart::xyy][2]+wmp2*R_temp[Cart::xxx][Cart::xyy][3]
R_temp[Cart::xxxz][Cart::yyz][2]+=pma2*R_temp[Cart::xxx][Cart::yyz][2]+wmp2*R_temp[Cart::xxx][Cart::yyz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::yy][3]
R_temp[Cart::xxxz][Cart::xxy][2]+=pma2*R_temp[Cart::xxx][Cart::xxy][2]+wmp2*R_temp[Cart::xxx][Cart::xxy][3]
R_temp[Cart::xxxz][Cart::xyz][2]+=pma2*R_temp[Cart::xxx][Cart::xyz][2]+wmp2*R_temp[Cart::xxx][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xy][3]
R_temp[Cart::xxxz][Cart::yzz][2]+=pma2*R_temp[Cart::xxx][Cart::yzz][2]+wmp2*R_temp[Cart::xxx][Cart::yzz][3]+0.5/_decay*2*R_temp[Cart::xxx][Cart::yz][3]
R_temp[Cart::xxxz][Cart::xxx][2]+=pma2*R_temp[Cart::xxx][Cart::xxx][2]+wmp2*R_temp[Cart::xxx][Cart::xxx][3]
R_temp[Cart::xxxz][Cart::xxz][2]+=pma2*R_temp[Cart::xxx][Cart::xxz][2]+wmp2*R_temp[Cart::xxx][Cart::xxz][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xx][3]
R_temp[Cart::xxxz][Cart::xzz][2]+=pma2*R_temp[Cart::xxx][Cart::xzz][2]+wmp2*R_temp[Cart::xxx][Cart::xzz][3]+0.5/_decay*2*R_temp[Cart::xxx][Cart::xz][3]
R_temp[Cart::xxxz][Cart::zzz][2]+=pma2*R_temp[Cart::xxx][Cart::zzz][2]+wmp2*R_temp[Cart::xxx][Cart::zzz][3]+0.5/_decay*3*R_temp[Cart::xxx][Cart::zz][3]
R_temp[Cart::xxzz][Cart::yyy][2]+=pma0*R_temp[Cart::xzz][Cart::yyy][2]+wmp0*R_temp[Cart::xzz][Cart::yyy][3]
R_temp[Cart::xxzz][Cart::xyy][2]+=pma0*R_temp[Cart::xzz][Cart::xyy][2]+wmp0*R_temp[Cart::xzz][Cart::xyy][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::yy][3]
R_temp[Cart::xxzz][Cart::yyz][2]+=pma0*R_temp[Cart::xzz][Cart::yyz][2]+wmp0*R_temp[Cart::xzz][Cart::yyz][3]
R_temp[Cart::xxzz][Cart::xxy][2]+=pma0*R_temp[Cart::xzz][Cart::xxy][2]+wmp0*R_temp[Cart::xzz][Cart::xxy][3]+0.5/_decay*2*R_temp[Cart::xzz][Cart::xy][3]
R_temp[Cart::xxzz][Cart::xyz][2]+=pma0*R_temp[Cart::xzz][Cart::xyz][2]+wmp0*R_temp[Cart::xzz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::yz][3]
R_temp[Cart::xxzz][Cart::yzz][2]+=pma0*R_temp[Cart::xzz][Cart::yzz][2]+wmp0*R_temp[Cart::xzz][Cart::yzz][3]
R_temp[Cart::xxzz][Cart::xxx][2]+=pma0*R_temp[Cart::xzz][Cart::xxx][2]+wmp0*R_temp[Cart::xzz][Cart::xxx][3]+0.5/_decay*3*R_temp[Cart::xzz][Cart::xx][3]
R_temp[Cart::xxzz][Cart::xxz][2]+=pma0*R_temp[Cart::xzz][Cart::xxz][2]+wmp0*R_temp[Cart::xzz][Cart::xxz][3]+0.5/_decay*2*R_temp[Cart::xzz][Cart::xz][3]
R_temp[Cart::xxzz][Cart::xzz][2]+=pma0*R_temp[Cart::xzz][Cart::xzz][2]+wmp0*R_temp[Cart::xzz][Cart::xzz][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::zz][3]
R_temp[Cart::xxzz][Cart::zzz][2]+=pma0*R_temp[Cart::xzz][Cart::zzz][2]+wmp0*R_temp[Cart::xzz][Cart::zzz][3]
R_temp[Cart::xzzz][Cart::yyy][2]+=pma0*R_temp[Cart::zzz][Cart::yyy][2]+wmp0*R_temp[Cart::zzz][Cart::yyy][3]
R_temp[Cart::xzzz][Cart::xyy][2]+=pma0*R_temp[Cart::zzz][Cart::xyy][2]+wmp0*R_temp[Cart::zzz][Cart::xyy][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::yy][3]
R_temp[Cart::xzzz][Cart::yyz][2]+=pma0*R_temp[Cart::zzz][Cart::yyz][2]+wmp0*R_temp[Cart::zzz][Cart::yyz][3]
R_temp[Cart::xzzz][Cart::xxy][2]+=pma0*R_temp[Cart::zzz][Cart::xxy][2]+wmp0*R_temp[Cart::zzz][Cart::xxy][3]+0.5/_decay*2*R_temp[Cart::zzz][Cart::xy][3]
R_temp[Cart::xzzz][Cart::xyz][2]+=pma0*R_temp[Cart::zzz][Cart::xyz][2]+wmp0*R_temp[Cart::zzz][Cart::xyz][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::yz][3]
R_temp[Cart::xzzz][Cart::yzz][2]+=pma0*R_temp[Cart::zzz][Cart::yzz][2]+wmp0*R_temp[Cart::zzz][Cart::yzz][3]
R_temp[Cart::xzzz][Cart::xxx][2]+=pma0*R_temp[Cart::zzz][Cart::xxx][2]+wmp0*R_temp[Cart::zzz][Cart::xxx][3]+0.5/_decay*3*R_temp[Cart::zzz][Cart::xx][3]
R_temp[Cart::xzzz][Cart::xxz][2]+=pma0*R_temp[Cart::zzz][Cart::xxz][2]+wmp0*R_temp[Cart::zzz][Cart::xxz][3]+0.5/_decay*2*R_temp[Cart::zzz][Cart::xz][3]
R_temp[Cart::xzzz][Cart::xzz][2]+=pma0*R_temp[Cart::zzz][Cart::xzz][2]+wmp0*R_temp[Cart::zzz][Cart::xzz][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::zz][3]
R_temp[Cart::xzzz][Cart::zzz][2]+=pma0*R_temp[Cart::zzz][Cart::zzz][2]+wmp0*R_temp[Cart::zzz][Cart::zzz][3]
R_temp[Cart::zzzz][Cart::yyy][2]+=pma2*R_temp[Cart::zzz][Cart::yyy][2]+wmp2*R_temp[Cart::zzz][Cart::yyy][3]+2*rzeta*(R_temp[Cart::zz][Cart::yyy][2]-gfak*R_temp[Cart::zz][Cart::yyy][3])
R_temp[Cart::zzzz][Cart::xyy][2]+=pma2*R_temp[Cart::zzz][Cart::xyy][2]+wmp2*R_temp[Cart::zzz][Cart::xyy][3]+2*rzeta*(R_temp[Cart::zz][Cart::xyy][2]-gfak*R_temp[Cart::zz][Cart::xyy][3])
R_temp[Cart::zzzz][Cart::yyz][2]+=pma2*R_temp[Cart::zzz][Cart::yyz][2]+wmp2*R_temp[Cart::zzz][Cart::yyz][3]+2*rzeta*(R_temp[Cart::zz][Cart::yyz][2]-gfak*R_temp[Cart::zz][Cart::yyz][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::yy][3]
R_temp[Cart::zzzz][Cart::xxy][2]+=pma2*R_temp[Cart::zzz][Cart::xxy][2]+wmp2*R_temp[Cart::zzz][Cart::xxy][3]+2*rzeta*(R_temp[Cart::zz][Cart::xxy][2]-gfak*R_temp[Cart::zz][Cart::xxy][3])
R_temp[Cart::zzzz][Cart::xyz][2]+=pma2*R_temp[Cart::zzz][Cart::xyz][2]+wmp2*R_temp[Cart::zzz][Cart::xyz][3]+2*rzeta*(R_temp[Cart::zz][Cart::xyz][2]-gfak*R_temp[Cart::zz][Cart::xyz][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::xy][3]
R_temp[Cart::zzzz][Cart::yzz][2]+=pma2*R_temp[Cart::zzz][Cart::yzz][2]+wmp2*R_temp[Cart::zzz][Cart::yzz][3]+2*rzeta*(R_temp[Cart::zz][Cart::yzz][2]-gfak*R_temp[Cart::zz][Cart::yzz][3])+0.5/_decay*2*R_temp[Cart::zzz][Cart::yz][3]
R_temp[Cart::zzzz][Cart::xxx][2]+=pma2*R_temp[Cart::zzz][Cart::xxx][2]+wmp2*R_temp[Cart::zzz][Cart::xxx][3]+2*rzeta*(R_temp[Cart::zz][Cart::xxx][2]-gfak*R_temp[Cart::zz][Cart::xxx][3])
R_temp[Cart::zzzz][Cart::xxz][2]+=pma2*R_temp[Cart::zzz][Cart::xxz][2]+wmp2*R_temp[Cart::zzz][Cart::xxz][3]+2*rzeta*(R_temp[Cart::zz][Cart::xxz][2]-gfak*R_temp[Cart::zz][Cart::xxz][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::xx][3]
R_temp[Cart::zzzz][Cart::xzz][2]+=pma2*R_temp[Cart::zzz][Cart::xzz][2]+wmp2*R_temp[Cart::zzz][Cart::xzz][3]+2*rzeta*(R_temp[Cart::zz][Cart::xzz][2]-gfak*R_temp[Cart::zz][Cart::xzz][3])+0.5/_decay*2*R_temp[Cart::zzz][Cart::xz][3]
R_temp[Cart::zzzz][Cart::zzz][2]+=pma2*R_temp[Cart::zzz][Cart::zzz][2]+wmp2*R_temp[Cart::zzz][Cart::zzz][3]+2*rzeta*(R_temp[Cart::zz][Cart::zzz][2]-gfak*R_temp[Cart::zz][Cart::zzz][3])+0.5/_decay*3*R_temp[Cart::zzz][Cart::zz][3]
}}
//------------------------------------------------------

//Integral h - s - s - m1
if (_mmax >5 ){
if (_lmax_alpha>4){

R_temp[Cart::yyyyy][Cart::s][1]+=pma1*R_temp[Cart::yyyy][Cart::s][1]+wmp1*R_temp[Cart::yyyy][Cart::s][2]+3*rzeta*(R_temp[Cart::yyy][Cart::s][1]-gfak*R_temp[Cart::yyy][Cart::s][2])
R_temp[Cart::xyyyy][Cart::s][1]+=pma0*R_temp[Cart::yyyy][Cart::s][1]+wmp0*R_temp[Cart::yyyy][Cart::s][2]
R_temp[Cart::yyyyz][Cart::s][1]+=pma2*R_temp[Cart::yyyy][Cart::s][1]+wmp2*R_temp[Cart::yyyy][Cart::s][2]
R_temp[Cart::xxyyy][Cart::s][1]+=pma0*R_temp[Cart::xyyy][Cart::s][1]+wmp0*R_temp[Cart::xyyy][Cart::s][2]
R_temp[Cart::xyyyz][Cart::s][1]+=pma0*R_temp[Cart::yyyz][Cart::s][1]+wmp0*R_temp[Cart::yyyz][Cart::s][2]
R_temp[Cart::yyyzz][Cart::s][1]+=pma2*R_temp[Cart::yyyz][Cart::s][1]+wmp2*R_temp[Cart::yyyz][Cart::s][2]
R_temp[Cart::xxxyy][Cart::s][1]+=pma1*R_temp[Cart::xxxy][Cart::s][1]+wmp1*R_temp[Cart::xxxy][Cart::s][2]
R_temp[Cart::xxyyz][Cart::s][1]+=pma2*R_temp[Cart::xxyy][Cart::s][1]+wmp2*R_temp[Cart::xxyy][Cart::s][2]
R_temp[Cart::xyyzz][Cart::s][1]+=pma0*R_temp[Cart::yyzz][Cart::s][1]+wmp0*R_temp[Cart::yyzz][Cart::s][2]
R_temp[Cart::yyzzz][Cart::s][1]+=pma1*R_temp[Cart::yzzz][Cart::s][1]+wmp1*R_temp[Cart::yzzz][Cart::s][2]
R_temp[Cart::xxxxy][Cart::s][1]+=pma1*R_temp[Cart::xxxx][Cart::s][1]+wmp1*R_temp[Cart::xxxx][Cart::s][2]
R_temp[Cart::xxxyz][Cart::s][1]+=pma1*R_temp[Cart::xxxz][Cart::s][1]+wmp1*R_temp[Cart::xxxz][Cart::s][2]
R_temp[Cart::xxyzz][Cart::s][1]+=pma1*R_temp[Cart::xxzz][Cart::s][1]+wmp1*R_temp[Cart::xxzz][Cart::s][2]
R_temp[Cart::xyzzz][Cart::s][1]+=pma0*R_temp[Cart::yzzz][Cart::s][1]+wmp0*R_temp[Cart::yzzz][Cart::s][2]
R_temp[Cart::yzzzz][Cart::s][1]+=pma1*R_temp[Cart::zzzz][Cart::s][1]+wmp1*R_temp[Cart::zzzz][Cart::s][2]
R_temp[Cart::xxxxx][Cart::s][1]+=pma0*R_temp[Cart::xxxx][Cart::s][1]+wmp0*R_temp[Cart::xxxx][Cart::s][2]+3*rzeta*(R_temp[Cart::xxx][Cart::s][1]-gfak*R_temp[Cart::xxx][Cart::s][2])
R_temp[Cart::xxxxz][Cart::s][1]+=pma2*R_temp[Cart::xxxx][Cart::s][1]+wmp2*R_temp[Cart::xxxx][Cart::s][2]
R_temp[Cart::xxxzz][Cart::s][1]+=pma2*R_temp[Cart::xxxz][Cart::s][1]+wmp2*R_temp[Cart::xxxz][Cart::s][2]
R_temp[Cart::xxzzz][Cart::s][1]+=pma0*R_temp[Cart::xzzz][Cart::s][1]+wmp0*R_temp[Cart::xzzz][Cart::s][2]
R_temp[Cart::xzzzz][Cart::s][1]+=pma0*R_temp[Cart::zzzz][Cart::s][1]+wmp0*R_temp[Cart::zzzz][Cart::s][2]
R_temp[Cart::zzzzz][Cart::s][1]+=pma2*R_temp[Cart::zzzz][Cart::s][1]+wmp2*R_temp[Cart::zzzz][Cart::s][2]+3*rzeta*(R_temp[Cart::zzz][Cart::s][1]-gfak*R_temp[Cart::zzz][Cart::s][2])
}}
//------------------------------------------------------

//Integral h - s - p - m1
if (_mmax >6 ){
if (_lmax_alpha>4 && _lmax_gamma>0){

R_temp[Cart::yyyyy][Cart::y][1]+=pma1*R_temp[Cart::yyyy][Cart::y][1]+wmp1*R_temp[Cart::yyyy][Cart::y][2]+3*rzeta*(R_temp[Cart::yyy][Cart::y][1]-gfak*R_temp[Cart::yyy][Cart::y][2])+0.5/_decay*1*R_temp[Cart::yyyy][Cart::s][2]
R_temp[Cart::yyyyy][Cart::x][1]+=pma1*R_temp[Cart::yyyy][Cart::x][1]+wmp1*R_temp[Cart::yyyy][Cart::x][2]+3*rzeta*(R_temp[Cart::yyy][Cart::x][1]-gfak*R_temp[Cart::yyy][Cart::x][2])
R_temp[Cart::yyyyy][Cart::z][1]+=pma1*R_temp[Cart::yyyy][Cart::z][1]+wmp1*R_temp[Cart::yyyy][Cart::z][2]+3*rzeta*(R_temp[Cart::yyy][Cart::z][1]-gfak*R_temp[Cart::yyy][Cart::z][2])
R_temp[Cart::xyyyy][Cart::y][1]+=pma0*R_temp[Cart::yyyy][Cart::y][1]+wmp0*R_temp[Cart::yyyy][Cart::y][2]
R_temp[Cart::xyyyy][Cart::x][1]+=pma0*R_temp[Cart::yyyy][Cart::x][1]+wmp0*R_temp[Cart::yyyy][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::s][2]
R_temp[Cart::xyyyy][Cart::z][1]+=pma0*R_temp[Cart::yyyy][Cart::z][1]+wmp0*R_temp[Cart::yyyy][Cart::z][2]
R_temp[Cart::yyyyz][Cart::y][1]+=pma2*R_temp[Cart::yyyy][Cart::y][1]+wmp2*R_temp[Cart::yyyy][Cart::y][2]
R_temp[Cart::yyyyz][Cart::x][1]+=pma2*R_temp[Cart::yyyy][Cart::x][1]+wmp2*R_temp[Cart::yyyy][Cart::x][2]
R_temp[Cart::yyyyz][Cart::z][1]+=pma2*R_temp[Cart::yyyy][Cart::z][1]+wmp2*R_temp[Cart::yyyy][Cart::z][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::s][2]
R_temp[Cart::xxyyy][Cart::y][1]+=pma0*R_temp[Cart::xyyy][Cart::y][1]+wmp0*R_temp[Cart::xyyy][Cart::y][2]
R_temp[Cart::xxyyy][Cart::x][1]+=pma0*R_temp[Cart::xyyy][Cart::x][1]+wmp0*R_temp[Cart::xyyy][Cart::x][2]+0.5/_decay*1*R_temp[Cart::xyyy][Cart::s][2]
R_temp[Cart::xxyyy][Cart::z][1]+=pma0*R_temp[Cart::xyyy][Cart::z][1]+wmp0*R_temp[Cart::xyyy][Cart::z][2]
R_temp[Cart::xyyyz][Cart::y][1]+=pma0*R_temp[Cart::yyyz][Cart::y][1]+wmp0*R_temp[Cart::yyyz][Cart::y][2]
R_temp[Cart::xyyyz][Cart::x][1]+=pma0*R_temp[Cart::yyyz][Cart::x][1]+wmp0*R_temp[Cart::yyyz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::s][2]
R_temp[Cart::xyyyz][Cart::z][1]+=pma0*R_temp[Cart::yyyz][Cart::z][1]+wmp0*R_temp[Cart::yyyz][Cart::z][2]
R_temp[Cart::yyyzz][Cart::y][1]+=pma2*R_temp[Cart::yyyz][Cart::y][1]+wmp2*R_temp[Cart::yyyz][Cart::y][2]
R_temp[Cart::yyyzz][Cart::x][1]+=pma2*R_temp[Cart::yyyz][Cart::x][1]+wmp2*R_temp[Cart::yyyz][Cart::x][2]
R_temp[Cart::yyyzz][Cart::z][1]+=pma2*R_temp[Cart::yyyz][Cart::z][1]+wmp2*R_temp[Cart::yyyz][Cart::z][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::s][2]
R_temp[Cart::xxxyy][Cart::y][1]+=pma1*R_temp[Cart::xxxy][Cart::y][1]+wmp1*R_temp[Cart::xxxy][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xxxy][Cart::s][2]
R_temp[Cart::xxxyy][Cart::x][1]+=pma1*R_temp[Cart::xxxy][Cart::x][1]+wmp1*R_temp[Cart::xxxy][Cart::x][2]
R_temp[Cart::xxxyy][Cart::z][1]+=pma1*R_temp[Cart::xxxy][Cart::z][1]+wmp1*R_temp[Cart::xxxy][Cart::z][2]
R_temp[Cart::xxyyz][Cart::y][1]+=pma2*R_temp[Cart::xxyy][Cart::y][1]+wmp2*R_temp[Cart::xxyy][Cart::y][2]
R_temp[Cart::xxyyz][Cart::x][1]+=pma2*R_temp[Cart::xxyy][Cart::x][1]+wmp2*R_temp[Cart::xxyy][Cart::x][2]
R_temp[Cart::xxyyz][Cart::z][1]+=pma2*R_temp[Cart::xxyy][Cart::z][1]+wmp2*R_temp[Cart::xxyy][Cart::z][2]+0.5/_decay*1*R_temp[Cart::xxyy][Cart::s][2]
R_temp[Cart::xyyzz][Cart::y][1]+=pma0*R_temp[Cart::yyzz][Cart::y][1]+wmp0*R_temp[Cart::yyzz][Cart::y][2]
R_temp[Cart::xyyzz][Cart::x][1]+=pma0*R_temp[Cart::yyzz][Cart::x][1]+wmp0*R_temp[Cart::yyzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yyzz][Cart::s][2]
R_temp[Cart::xyyzz][Cart::z][1]+=pma0*R_temp[Cart::yyzz][Cart::z][1]+wmp0*R_temp[Cart::yyzz][Cart::z][2]
R_temp[Cart::yyzzz][Cart::y][1]+=pma1*R_temp[Cart::yzzz][Cart::y][1]+wmp1*R_temp[Cart::yzzz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::s][2]
R_temp[Cart::yyzzz][Cart::x][1]+=pma1*R_temp[Cart::yzzz][Cart::x][1]+wmp1*R_temp[Cart::yzzz][Cart::x][2]
R_temp[Cart::yyzzz][Cart::z][1]+=pma1*R_temp[Cart::yzzz][Cart::z][1]+wmp1*R_temp[Cart::yzzz][Cart::z][2]
R_temp[Cart::xxxxy][Cart::y][1]+=pma1*R_temp[Cart::xxxx][Cart::y][1]+wmp1*R_temp[Cart::xxxx][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::s][2]
R_temp[Cart::xxxxy][Cart::x][1]+=pma1*R_temp[Cart::xxxx][Cart::x][1]+wmp1*R_temp[Cart::xxxx][Cart::x][2]
R_temp[Cart::xxxxy][Cart::z][1]+=pma1*R_temp[Cart::xxxx][Cart::z][1]+wmp1*R_temp[Cart::xxxx][Cart::z][2]
R_temp[Cart::xxxyz][Cart::y][1]+=pma1*R_temp[Cart::xxxz][Cart::y][1]+wmp1*R_temp[Cart::xxxz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::s][2]
R_temp[Cart::xxxyz][Cart::x][1]+=pma1*R_temp[Cart::xxxz][Cart::x][1]+wmp1*R_temp[Cart::xxxz][Cart::x][2]
R_temp[Cart::xxxyz][Cart::z][1]+=pma1*R_temp[Cart::xxxz][Cart::z][1]+wmp1*R_temp[Cart::xxxz][Cart::z][2]
R_temp[Cart::xxyzz][Cart::y][1]+=pma1*R_temp[Cart::xxzz][Cart::y][1]+wmp1*R_temp[Cart::xxzz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xxzz][Cart::s][2]
R_temp[Cart::xxyzz][Cart::x][1]+=pma1*R_temp[Cart::xxzz][Cart::x][1]+wmp1*R_temp[Cart::xxzz][Cart::x][2]
R_temp[Cart::xxyzz][Cart::z][1]+=pma1*R_temp[Cart::xxzz][Cart::z][1]+wmp1*R_temp[Cart::xxzz][Cart::z][2]
R_temp[Cart::xyzzz][Cart::y][1]+=pma0*R_temp[Cart::yzzz][Cart::y][1]+wmp0*R_temp[Cart::yzzz][Cart::y][2]
R_temp[Cart::xyzzz][Cart::x][1]+=pma0*R_temp[Cart::yzzz][Cart::x][1]+wmp0*R_temp[Cart::yzzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::s][2]
R_temp[Cart::xyzzz][Cart::z][1]+=pma0*R_temp[Cart::yzzz][Cart::z][1]+wmp0*R_temp[Cart::yzzz][Cart::z][2]
R_temp[Cart::yzzzz][Cart::y][1]+=pma1*R_temp[Cart::zzzz][Cart::y][1]+wmp1*R_temp[Cart::zzzz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::s][2]
R_temp[Cart::yzzzz][Cart::x][1]+=pma1*R_temp[Cart::zzzz][Cart::x][1]+wmp1*R_temp[Cart::zzzz][Cart::x][2]
R_temp[Cart::yzzzz][Cart::z][1]+=pma1*R_temp[Cart::zzzz][Cart::z][1]+wmp1*R_temp[Cart::zzzz][Cart::z][2]
R_temp[Cart::xxxxx][Cart::y][1]+=pma0*R_temp[Cart::xxxx][Cart::y][1]+wmp0*R_temp[Cart::xxxx][Cart::y][2]+3*rzeta*(R_temp[Cart::xxx][Cart::y][1]-gfak*R_temp[Cart::xxx][Cart::y][2])
R_temp[Cart::xxxxx][Cart::x][1]+=pma0*R_temp[Cart::xxxx][Cart::x][1]+wmp0*R_temp[Cart::xxxx][Cart::x][2]+3*rzeta*(R_temp[Cart::xxx][Cart::x][1]-gfak*R_temp[Cart::xxx][Cart::x][2])+0.5/_decay*1*R_temp[Cart::xxxx][Cart::s][2]
R_temp[Cart::xxxxx][Cart::z][1]+=pma0*R_temp[Cart::xxxx][Cart::z][1]+wmp0*R_temp[Cart::xxxx][Cart::z][2]+3*rzeta*(R_temp[Cart::xxx][Cart::z][1]-gfak*R_temp[Cart::xxx][Cart::z][2])
R_temp[Cart::xxxxz][Cart::y][1]+=pma2*R_temp[Cart::xxxx][Cart::y][1]+wmp2*R_temp[Cart::xxxx][Cart::y][2]
R_temp[Cart::xxxxz][Cart::x][1]+=pma2*R_temp[Cart::xxxx][Cart::x][1]+wmp2*R_temp[Cart::xxxx][Cart::x][2]
R_temp[Cart::xxxxz][Cart::z][1]+=pma2*R_temp[Cart::xxxx][Cart::z][1]+wmp2*R_temp[Cart::xxxx][Cart::z][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::s][2]
R_temp[Cart::xxxzz][Cart::y][1]+=pma2*R_temp[Cart::xxxz][Cart::y][1]+wmp2*R_temp[Cart::xxxz][Cart::y][2]
R_temp[Cart::xxxzz][Cart::x][1]+=pma2*R_temp[Cart::xxxz][Cart::x][1]+wmp2*R_temp[Cart::xxxz][Cart::x][2]
R_temp[Cart::xxxzz][Cart::z][1]+=pma2*R_temp[Cart::xxxz][Cart::z][1]+wmp2*R_temp[Cart::xxxz][Cart::z][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::s][2]
R_temp[Cart::xxzzz][Cart::y][1]+=pma0*R_temp[Cart::xzzz][Cart::y][1]+wmp0*R_temp[Cart::xzzz][Cart::y][2]
R_temp[Cart::xxzzz][Cart::x][1]+=pma0*R_temp[Cart::xzzz][Cart::x][1]+wmp0*R_temp[Cart::xzzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::xzzz][Cart::s][2]
R_temp[Cart::xxzzz][Cart::z][1]+=pma0*R_temp[Cart::xzzz][Cart::z][1]+wmp0*R_temp[Cart::xzzz][Cart::z][2]
R_temp[Cart::xzzzz][Cart::y][1]+=pma0*R_temp[Cart::zzzz][Cart::y][1]+wmp0*R_temp[Cart::zzzz][Cart::y][2]
R_temp[Cart::xzzzz][Cart::x][1]+=pma0*R_temp[Cart::zzzz][Cart::x][1]+wmp0*R_temp[Cart::zzzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::s][2]
R_temp[Cart::xzzzz][Cart::z][1]+=pma0*R_temp[Cart::zzzz][Cart::z][1]+wmp0*R_temp[Cart::zzzz][Cart::z][2]
R_temp[Cart::zzzzz][Cart::y][1]+=pma2*R_temp[Cart::zzzz][Cart::y][1]+wmp2*R_temp[Cart::zzzz][Cart::y][2]+3*rzeta*(R_temp[Cart::zzz][Cart::y][1]-gfak*R_temp[Cart::zzz][Cart::y][2])
R_temp[Cart::zzzzz][Cart::x][1]+=pma2*R_temp[Cart::zzzz][Cart::x][1]+wmp2*R_temp[Cart::zzzz][Cart::x][2]+3*rzeta*(R_temp[Cart::zzz][Cart::x][1]-gfak*R_temp[Cart::zzz][Cart::x][2])
R_temp[Cart::zzzzz][Cart::z][1]+=pma2*R_temp[Cart::zzzz][Cart::z][1]+wmp2*R_temp[Cart::zzzz][Cart::z][2]+3*rzeta*(R_temp[Cart::zzz][Cart::z][1]-gfak*R_temp[Cart::zzz][Cart::z][2])+0.5/_decay*1*R_temp[Cart::zzzz][Cart::s][2]
}}
//------------------------------------------------------

//Integral h - s - d - m1
if (_mmax >7 ){
if (_lmax_alpha>4 && _lmax_gamma>1){

R_temp[Cart::yyyyy][Cart::yy][1]+=pma1*R_temp[Cart::yyyy][Cart::yy][1]+wmp1*R_temp[Cart::yyyy][Cart::yy][2]+3*rzeta*(R_temp[Cart::yyy][Cart::yy][1]-gfak*R_temp[Cart::yyy][Cart::yy][2])+0.5/_decay*2*R_temp[Cart::yyyy][Cart::y][2]
R_temp[Cart::yyyyy][Cart::xy][1]+=pma1*R_temp[Cart::yyyy][Cart::xy][1]+wmp1*R_temp[Cart::yyyy][Cart::xy][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xy][1]-gfak*R_temp[Cart::yyy][Cart::xy][2])+0.5/_decay*1*R_temp[Cart::yyyy][Cart::x][2]
R_temp[Cart::yyyyy][Cart::yz][1]+=pma1*R_temp[Cart::yyyy][Cart::yz][1]+wmp1*R_temp[Cart::yyyy][Cart::yz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::yz][1]-gfak*R_temp[Cart::yyy][Cart::yz][2])+0.5/_decay*1*R_temp[Cart::yyyy][Cart::z][2]
R_temp[Cart::yyyyy][Cart::xx][1]+=pma1*R_temp[Cart::yyyy][Cart::xx][1]+wmp1*R_temp[Cart::yyyy][Cart::xx][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xx][1]-gfak*R_temp[Cart::yyy][Cart::xx][2])
R_temp[Cart::yyyyy][Cart::xz][1]+=pma1*R_temp[Cart::yyyy][Cart::xz][1]+wmp1*R_temp[Cart::yyyy][Cart::xz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xz][1]-gfak*R_temp[Cart::yyy][Cart::xz][2])
R_temp[Cart::yyyyy][Cart::zz][1]+=pma1*R_temp[Cart::yyyy][Cart::zz][1]+wmp1*R_temp[Cart::yyyy][Cart::zz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::zz][1]-gfak*R_temp[Cart::yyy][Cart::zz][2])
R_temp[Cart::xyyyy][Cart::yy][1]+=pma0*R_temp[Cart::yyyy][Cart::yy][1]+wmp0*R_temp[Cart::yyyy][Cart::yy][2]
R_temp[Cart::xyyyy][Cart::xy][1]+=pma0*R_temp[Cart::yyyy][Cart::xy][1]+wmp0*R_temp[Cart::yyyy][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::y][2]
R_temp[Cart::xyyyy][Cart::yz][1]+=pma0*R_temp[Cart::yyyy][Cart::yz][1]+wmp0*R_temp[Cart::yyyy][Cart::yz][2]
R_temp[Cart::xyyyy][Cart::xx][1]+=pma0*R_temp[Cart::yyyy][Cart::xx][1]+wmp0*R_temp[Cart::yyyy][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yyyy][Cart::x][2]
R_temp[Cart::xyyyy][Cart::xz][1]+=pma0*R_temp[Cart::yyyy][Cart::xz][1]+wmp0*R_temp[Cart::yyyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::z][2]
R_temp[Cart::xyyyy][Cart::zz][1]+=pma0*R_temp[Cart::yyyy][Cart::zz][1]+wmp0*R_temp[Cart::yyyy][Cart::zz][2]
R_temp[Cart::yyyyz][Cart::yy][1]+=pma2*R_temp[Cart::yyyy][Cart::yy][1]+wmp2*R_temp[Cart::yyyy][Cart::yy][2]
R_temp[Cart::yyyyz][Cart::xy][1]+=pma2*R_temp[Cart::yyyy][Cart::xy][1]+wmp2*R_temp[Cart::yyyy][Cart::xy][2]
R_temp[Cart::yyyyz][Cart::yz][1]+=pma2*R_temp[Cart::yyyy][Cart::yz][1]+wmp2*R_temp[Cart::yyyy][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::y][2]
R_temp[Cart::yyyyz][Cart::xx][1]+=pma2*R_temp[Cart::yyyy][Cart::xx][1]+wmp2*R_temp[Cart::yyyy][Cart::xx][2]
R_temp[Cart::yyyyz][Cart::xz][1]+=pma2*R_temp[Cart::yyyy][Cart::xz][1]+wmp2*R_temp[Cart::yyyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::x][2]
R_temp[Cart::yyyyz][Cart::zz][1]+=pma2*R_temp[Cart::yyyy][Cart::zz][1]+wmp2*R_temp[Cart::yyyy][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::yyyy][Cart::z][2]
R_temp[Cart::xxyyy][Cart::yy][1]+=pma0*R_temp[Cart::xyyy][Cart::yy][1]+wmp0*R_temp[Cart::xyyy][Cart::yy][2]
R_temp[Cart::xxyyy][Cart::xy][1]+=pma0*R_temp[Cart::xyyy][Cart::xy][1]+wmp0*R_temp[Cart::xyyy][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xyyy][Cart::y][2]
R_temp[Cart::xxyyy][Cart::yz][1]+=pma0*R_temp[Cart::xyyy][Cart::yz][1]+wmp0*R_temp[Cart::xyyy][Cart::yz][2]
R_temp[Cart::xxyyy][Cart::xx][1]+=pma0*R_temp[Cart::xyyy][Cart::xx][1]+wmp0*R_temp[Cart::xyyy][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::xyyy][Cart::x][2]
R_temp[Cart::xxyyy][Cart::xz][1]+=pma0*R_temp[Cart::xyyy][Cart::xz][1]+wmp0*R_temp[Cart::xyyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xyyy][Cart::z][2]
R_temp[Cart::xxyyy][Cart::zz][1]+=pma0*R_temp[Cart::xyyy][Cart::zz][1]+wmp0*R_temp[Cart::xyyy][Cart::zz][2]
R_temp[Cart::xyyyz][Cart::yy][1]+=pma0*R_temp[Cart::yyyz][Cart::yy][1]+wmp0*R_temp[Cart::yyyz][Cart::yy][2]
R_temp[Cart::xyyyz][Cart::xy][1]+=pma0*R_temp[Cart::yyyz][Cart::xy][1]+wmp0*R_temp[Cart::yyyz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::y][2]
R_temp[Cart::xyyyz][Cart::yz][1]+=pma0*R_temp[Cart::yyyz][Cart::yz][1]+wmp0*R_temp[Cart::yyyz][Cart::yz][2]
R_temp[Cart::xyyyz][Cart::xx][1]+=pma0*R_temp[Cart::yyyz][Cart::xx][1]+wmp0*R_temp[Cart::yyyz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yyyz][Cart::x][2]
R_temp[Cart::xyyyz][Cart::xz][1]+=pma0*R_temp[Cart::yyyz][Cart::xz][1]+wmp0*R_temp[Cart::yyyz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::z][2]
R_temp[Cart::xyyyz][Cart::zz][1]+=pma0*R_temp[Cart::yyyz][Cart::zz][1]+wmp0*R_temp[Cart::yyyz][Cart::zz][2]
R_temp[Cart::yyyzz][Cart::yy][1]+=pma2*R_temp[Cart::yyyz][Cart::yy][1]+wmp2*R_temp[Cart::yyyz][Cart::yy][2]
R_temp[Cart::yyyzz][Cart::xy][1]+=pma2*R_temp[Cart::yyyz][Cart::xy][1]+wmp2*R_temp[Cart::yyyz][Cart::xy][2]
R_temp[Cart::yyyzz][Cart::yz][1]+=pma2*R_temp[Cart::yyyz][Cart::yz][1]+wmp2*R_temp[Cart::yyyz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::y][2]
R_temp[Cart::yyyzz][Cart::xx][1]+=pma2*R_temp[Cart::yyyz][Cart::xx][1]+wmp2*R_temp[Cart::yyyz][Cart::xx][2]
R_temp[Cart::yyyzz][Cart::xz][1]+=pma2*R_temp[Cart::yyyz][Cart::xz][1]+wmp2*R_temp[Cart::yyyz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::x][2]
R_temp[Cart::yyyzz][Cart::zz][1]+=pma2*R_temp[Cart::yyyz][Cart::zz][1]+wmp2*R_temp[Cart::yyyz][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::yyyz][Cart::z][2]
R_temp[Cart::xxxyy][Cart::yy][1]+=pma1*R_temp[Cart::xxxy][Cart::yy][1]+wmp1*R_temp[Cart::xxxy][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xxxy][Cart::y][2]
R_temp[Cart::xxxyy][Cart::xy][1]+=pma1*R_temp[Cart::xxxy][Cart::xy][1]+wmp1*R_temp[Cart::xxxy][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xxxy][Cart::x][2]
R_temp[Cart::xxxyy][Cart::yz][1]+=pma1*R_temp[Cart::xxxy][Cart::yz][1]+wmp1*R_temp[Cart::xxxy][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxxy][Cart::z][2]
R_temp[Cart::xxxyy][Cart::xx][1]+=pma1*R_temp[Cart::xxxy][Cart::xx][1]+wmp1*R_temp[Cart::xxxy][Cart::xx][2]
R_temp[Cart::xxxyy][Cart::xz][1]+=pma1*R_temp[Cart::xxxy][Cart::xz][1]+wmp1*R_temp[Cart::xxxy][Cart::xz][2]
R_temp[Cart::xxxyy][Cart::zz][1]+=pma1*R_temp[Cart::xxxy][Cart::zz][1]+wmp1*R_temp[Cart::xxxy][Cart::zz][2]
R_temp[Cart::xxyyz][Cart::yy][1]+=pma2*R_temp[Cart::xxyy][Cart::yy][1]+wmp2*R_temp[Cart::xxyy][Cart::yy][2]
R_temp[Cart::xxyyz][Cart::xy][1]+=pma2*R_temp[Cart::xxyy][Cart::xy][1]+wmp2*R_temp[Cart::xxyy][Cart::xy][2]
R_temp[Cart::xxyyz][Cart::yz][1]+=pma2*R_temp[Cart::xxyy][Cart::yz][1]+wmp2*R_temp[Cart::xxyy][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxyy][Cart::y][2]
R_temp[Cart::xxyyz][Cart::xx][1]+=pma2*R_temp[Cart::xxyy][Cart::xx][1]+wmp2*R_temp[Cart::xxyy][Cart::xx][2]
R_temp[Cart::xxyyz][Cart::xz][1]+=pma2*R_temp[Cart::xxyy][Cart::xz][1]+wmp2*R_temp[Cart::xxyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xxyy][Cart::x][2]
R_temp[Cart::xxyyz][Cart::zz][1]+=pma2*R_temp[Cart::xxyy][Cart::zz][1]+wmp2*R_temp[Cart::xxyy][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::xxyy][Cart::z][2]
R_temp[Cart::xyyzz][Cart::yy][1]+=pma0*R_temp[Cart::yyzz][Cart::yy][1]+wmp0*R_temp[Cart::yyzz][Cart::yy][2]
R_temp[Cart::xyyzz][Cart::xy][1]+=pma0*R_temp[Cart::yyzz][Cart::xy][1]+wmp0*R_temp[Cart::yyzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yyzz][Cart::y][2]
R_temp[Cart::xyyzz][Cart::yz][1]+=pma0*R_temp[Cart::yyzz][Cart::yz][1]+wmp0*R_temp[Cart::yyzz][Cart::yz][2]
R_temp[Cart::xyyzz][Cart::xx][1]+=pma0*R_temp[Cart::yyzz][Cart::xx][1]+wmp0*R_temp[Cart::yyzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yyzz][Cart::x][2]
R_temp[Cart::xyyzz][Cart::xz][1]+=pma0*R_temp[Cart::yyzz][Cart::xz][1]+wmp0*R_temp[Cart::yyzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyzz][Cart::z][2]
R_temp[Cart::xyyzz][Cart::zz][1]+=pma0*R_temp[Cart::yyzz][Cart::zz][1]+wmp0*R_temp[Cart::yyzz][Cart::zz][2]
R_temp[Cart::yyzzz][Cart::yy][1]+=pma1*R_temp[Cart::yzzz][Cart::yy][1]+wmp1*R_temp[Cart::yzzz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::yzzz][Cart::y][2]
R_temp[Cart::yyzzz][Cart::xy][1]+=pma1*R_temp[Cart::yzzz][Cart::xy][1]+wmp1*R_temp[Cart::yzzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::x][2]
R_temp[Cart::yyzzz][Cart::yz][1]+=pma1*R_temp[Cart::yzzz][Cart::yz][1]+wmp1*R_temp[Cart::yzzz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::z][2]
R_temp[Cart::yyzzz][Cart::xx][1]+=pma1*R_temp[Cart::yzzz][Cart::xx][1]+wmp1*R_temp[Cart::yzzz][Cart::xx][2]
R_temp[Cart::yyzzz][Cart::xz][1]+=pma1*R_temp[Cart::yzzz][Cart::xz][1]+wmp1*R_temp[Cart::yzzz][Cart::xz][2]
R_temp[Cart::yyzzz][Cart::zz][1]+=pma1*R_temp[Cart::yzzz][Cart::zz][1]+wmp1*R_temp[Cart::yzzz][Cart::zz][2]
R_temp[Cart::xxxxy][Cart::yy][1]+=pma1*R_temp[Cart::xxxx][Cart::yy][1]+wmp1*R_temp[Cart::xxxx][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xxxx][Cart::y][2]
R_temp[Cart::xxxxy][Cart::xy][1]+=pma1*R_temp[Cart::xxxx][Cart::xy][1]+wmp1*R_temp[Cart::xxxx][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::x][2]
R_temp[Cart::xxxxy][Cart::yz][1]+=pma1*R_temp[Cart::xxxx][Cart::yz][1]+wmp1*R_temp[Cart::xxxx][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::z][2]
R_temp[Cart::xxxxy][Cart::xx][1]+=pma1*R_temp[Cart::xxxx][Cart::xx][1]+wmp1*R_temp[Cart::xxxx][Cart::xx][2]
R_temp[Cart::xxxxy][Cart::xz][1]+=pma1*R_temp[Cart::xxxx][Cart::xz][1]+wmp1*R_temp[Cart::xxxx][Cart::xz][2]
R_temp[Cart::xxxxy][Cart::zz][1]+=pma1*R_temp[Cart::xxxx][Cart::zz][1]+wmp1*R_temp[Cart::xxxx][Cart::zz][2]
R_temp[Cart::xxxyz][Cart::yy][1]+=pma1*R_temp[Cart::xxxz][Cart::yy][1]+wmp1*R_temp[Cart::xxxz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xxxz][Cart::y][2]
R_temp[Cart::xxxyz][Cart::xy][1]+=pma1*R_temp[Cart::xxxz][Cart::xy][1]+wmp1*R_temp[Cart::xxxz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::x][2]
R_temp[Cart::xxxyz][Cart::yz][1]+=pma1*R_temp[Cart::xxxz][Cart::yz][1]+wmp1*R_temp[Cart::xxxz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::z][2]
R_temp[Cart::xxxyz][Cart::xx][1]+=pma1*R_temp[Cart::xxxz][Cart::xx][1]+wmp1*R_temp[Cart::xxxz][Cart::xx][2]
R_temp[Cart::xxxyz][Cart::xz][1]+=pma1*R_temp[Cart::xxxz][Cart::xz][1]+wmp1*R_temp[Cart::xxxz][Cart::xz][2]
R_temp[Cart::xxxyz][Cart::zz][1]+=pma1*R_temp[Cart::xxxz][Cart::zz][1]+wmp1*R_temp[Cart::xxxz][Cart::zz][2]
R_temp[Cart::xxyzz][Cart::yy][1]+=pma1*R_temp[Cart::xxzz][Cart::yy][1]+wmp1*R_temp[Cart::xxzz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xxzz][Cart::y][2]
R_temp[Cart::xxyzz][Cart::xy][1]+=pma1*R_temp[Cart::xxzz][Cart::xy][1]+wmp1*R_temp[Cart::xxzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xxzz][Cart::x][2]
R_temp[Cart::xxyzz][Cart::yz][1]+=pma1*R_temp[Cart::xxzz][Cart::yz][1]+wmp1*R_temp[Cart::xxzz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxzz][Cart::z][2]
R_temp[Cart::xxyzz][Cart::xx][1]+=pma1*R_temp[Cart::xxzz][Cart::xx][1]+wmp1*R_temp[Cart::xxzz][Cart::xx][2]
R_temp[Cart::xxyzz][Cart::xz][1]+=pma1*R_temp[Cart::xxzz][Cart::xz][1]+wmp1*R_temp[Cart::xxzz][Cart::xz][2]
R_temp[Cart::xxyzz][Cart::zz][1]+=pma1*R_temp[Cart::xxzz][Cart::zz][1]+wmp1*R_temp[Cart::xxzz][Cart::zz][2]
R_temp[Cart::xyzzz][Cart::yy][1]+=pma0*R_temp[Cart::yzzz][Cart::yy][1]+wmp0*R_temp[Cart::yzzz][Cart::yy][2]
R_temp[Cart::xyzzz][Cart::xy][1]+=pma0*R_temp[Cart::yzzz][Cart::xy][1]+wmp0*R_temp[Cart::yzzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::y][2]
R_temp[Cart::xyzzz][Cart::yz][1]+=pma0*R_temp[Cart::yzzz][Cart::yz][1]+wmp0*R_temp[Cart::yzzz][Cart::yz][2]
R_temp[Cart::xyzzz][Cart::xx][1]+=pma0*R_temp[Cart::yzzz][Cart::xx][1]+wmp0*R_temp[Cart::yzzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yzzz][Cart::x][2]
R_temp[Cart::xyzzz][Cart::xz][1]+=pma0*R_temp[Cart::yzzz][Cart::xz][1]+wmp0*R_temp[Cart::yzzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::z][2]
R_temp[Cart::xyzzz][Cart::zz][1]+=pma0*R_temp[Cart::yzzz][Cart::zz][1]+wmp0*R_temp[Cart::yzzz][Cart::zz][2]
R_temp[Cart::yzzzz][Cart::yy][1]+=pma1*R_temp[Cart::zzzz][Cart::yy][1]+wmp1*R_temp[Cart::zzzz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::zzzz][Cart::y][2]
R_temp[Cart::yzzzz][Cart::xy][1]+=pma1*R_temp[Cart::zzzz][Cart::xy][1]+wmp1*R_temp[Cart::zzzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::x][2]
R_temp[Cart::yzzzz][Cart::yz][1]+=pma1*R_temp[Cart::zzzz][Cart::yz][1]+wmp1*R_temp[Cart::zzzz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::z][2]
R_temp[Cart::yzzzz][Cart::xx][1]+=pma1*R_temp[Cart::zzzz][Cart::xx][1]+wmp1*R_temp[Cart::zzzz][Cart::xx][2]
R_temp[Cart::yzzzz][Cart::xz][1]+=pma1*R_temp[Cart::zzzz][Cart::xz][1]+wmp1*R_temp[Cart::zzzz][Cart::xz][2]
R_temp[Cart::yzzzz][Cart::zz][1]+=pma1*R_temp[Cart::zzzz][Cart::zz][1]+wmp1*R_temp[Cart::zzzz][Cart::zz][2]
R_temp[Cart::xxxxx][Cart::yy][1]+=pma0*R_temp[Cart::xxxx][Cart::yy][1]+wmp0*R_temp[Cart::xxxx][Cart::yy][2]+3*rzeta*(R_temp[Cart::xxx][Cart::yy][1]-gfak*R_temp[Cart::xxx][Cart::yy][2])
R_temp[Cart::xxxxx][Cart::xy][1]+=pma0*R_temp[Cart::xxxx][Cart::xy][1]+wmp0*R_temp[Cart::xxxx][Cart::xy][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xy][1]-gfak*R_temp[Cart::xxx][Cart::xy][2])+0.5/_decay*1*R_temp[Cart::xxxx][Cart::y][2]
R_temp[Cart::xxxxx][Cart::yz][1]+=pma0*R_temp[Cart::xxxx][Cart::yz][1]+wmp0*R_temp[Cart::xxxx][Cart::yz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::yz][1]-gfak*R_temp[Cart::xxx][Cart::yz][2])
R_temp[Cart::xxxxx][Cart::xx][1]+=pma0*R_temp[Cart::xxxx][Cart::xx][1]+wmp0*R_temp[Cart::xxxx][Cart::xx][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xx][1]-gfak*R_temp[Cart::xxx][Cart::xx][2])+0.5/_decay*2*R_temp[Cart::xxxx][Cart::x][2]
R_temp[Cart::xxxxx][Cart::xz][1]+=pma0*R_temp[Cart::xxxx][Cart::xz][1]+wmp0*R_temp[Cart::xxxx][Cart::xz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xz][1]-gfak*R_temp[Cart::xxx][Cart::xz][2])+0.5/_decay*1*R_temp[Cart::xxxx][Cart::z][2]
R_temp[Cart::xxxxx][Cart::zz][1]+=pma0*R_temp[Cart::xxxx][Cart::zz][1]+wmp0*R_temp[Cart::xxxx][Cart::zz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::zz][1]-gfak*R_temp[Cart::xxx][Cart::zz][2])
R_temp[Cart::xxxxz][Cart::yy][1]+=pma2*R_temp[Cart::xxxx][Cart::yy][1]+wmp2*R_temp[Cart::xxxx][Cart::yy][2]
R_temp[Cart::xxxxz][Cart::xy][1]+=pma2*R_temp[Cart::xxxx][Cart::xy][1]+wmp2*R_temp[Cart::xxxx][Cart::xy][2]
R_temp[Cart::xxxxz][Cart::yz][1]+=pma2*R_temp[Cart::xxxx][Cart::yz][1]+wmp2*R_temp[Cart::xxxx][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::y][2]
R_temp[Cart::xxxxz][Cart::xx][1]+=pma2*R_temp[Cart::xxxx][Cart::xx][1]+wmp2*R_temp[Cart::xxxx][Cart::xx][2]
R_temp[Cart::xxxxz][Cart::xz][1]+=pma2*R_temp[Cart::xxxx][Cart::xz][1]+wmp2*R_temp[Cart::xxxx][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::x][2]
R_temp[Cart::xxxxz][Cart::zz][1]+=pma2*R_temp[Cart::xxxx][Cart::zz][1]+wmp2*R_temp[Cart::xxxx][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::xxxx][Cart::z][2]
R_temp[Cart::xxxzz][Cart::yy][1]+=pma2*R_temp[Cart::xxxz][Cart::yy][1]+wmp2*R_temp[Cart::xxxz][Cart::yy][2]
R_temp[Cart::xxxzz][Cart::xy][1]+=pma2*R_temp[Cart::xxxz][Cart::xy][1]+wmp2*R_temp[Cart::xxxz][Cart::xy][2]
R_temp[Cart::xxxzz][Cart::yz][1]+=pma2*R_temp[Cart::xxxz][Cart::yz][1]+wmp2*R_temp[Cart::xxxz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::y][2]
R_temp[Cart::xxxzz][Cart::xx][1]+=pma2*R_temp[Cart::xxxz][Cart::xx][1]+wmp2*R_temp[Cart::xxxz][Cart::xx][2]
R_temp[Cart::xxxzz][Cart::xz][1]+=pma2*R_temp[Cart::xxxz][Cart::xz][1]+wmp2*R_temp[Cart::xxxz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::x][2]
R_temp[Cart::xxxzz][Cart::zz][1]+=pma2*R_temp[Cart::xxxz][Cart::zz][1]+wmp2*R_temp[Cart::xxxz][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::xxxz][Cart::z][2]
R_temp[Cart::xxzzz][Cart::yy][1]+=pma0*R_temp[Cart::xzzz][Cart::yy][1]+wmp0*R_temp[Cart::xzzz][Cart::yy][2]
R_temp[Cart::xxzzz][Cart::xy][1]+=pma0*R_temp[Cart::xzzz][Cart::xy][1]+wmp0*R_temp[Cart::xzzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xzzz][Cart::y][2]
R_temp[Cart::xxzzz][Cart::yz][1]+=pma0*R_temp[Cart::xzzz][Cart::yz][1]+wmp0*R_temp[Cart::xzzz][Cart::yz][2]
R_temp[Cart::xxzzz][Cart::xx][1]+=pma0*R_temp[Cart::xzzz][Cart::xx][1]+wmp0*R_temp[Cart::xzzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::xzzz][Cart::x][2]
R_temp[Cart::xxzzz][Cart::xz][1]+=pma0*R_temp[Cart::xzzz][Cart::xz][1]+wmp0*R_temp[Cart::xzzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xzzz][Cart::z][2]
R_temp[Cart::xxzzz][Cart::zz][1]+=pma0*R_temp[Cart::xzzz][Cart::zz][1]+wmp0*R_temp[Cart::xzzz][Cart::zz][2]
R_temp[Cart::xzzzz][Cart::yy][1]+=pma0*R_temp[Cart::zzzz][Cart::yy][1]+wmp0*R_temp[Cart::zzzz][Cart::yy][2]
R_temp[Cart::xzzzz][Cart::xy][1]+=pma0*R_temp[Cart::zzzz][Cart::xy][1]+wmp0*R_temp[Cart::zzzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::y][2]
R_temp[Cart::xzzzz][Cart::yz][1]+=pma0*R_temp[Cart::zzzz][Cart::yz][1]+wmp0*R_temp[Cart::zzzz][Cart::yz][2]
R_temp[Cart::xzzzz][Cart::xx][1]+=pma0*R_temp[Cart::zzzz][Cart::xx][1]+wmp0*R_temp[Cart::zzzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::zzzz][Cart::x][2]
R_temp[Cart::xzzzz][Cart::xz][1]+=pma0*R_temp[Cart::zzzz][Cart::xz][1]+wmp0*R_temp[Cart::zzzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::z][2]
R_temp[Cart::xzzzz][Cart::zz][1]+=pma0*R_temp[Cart::zzzz][Cart::zz][1]+wmp0*R_temp[Cart::zzzz][Cart::zz][2]
R_temp[Cart::zzzzz][Cart::yy][1]+=pma2*R_temp[Cart::zzzz][Cart::yy][1]+wmp2*R_temp[Cart::zzzz][Cart::yy][2]+3*rzeta*(R_temp[Cart::zzz][Cart::yy][1]-gfak*R_temp[Cart::zzz][Cart::yy][2])
R_temp[Cart::zzzzz][Cart::xy][1]+=pma2*R_temp[Cart::zzzz][Cart::xy][1]+wmp2*R_temp[Cart::zzzz][Cart::xy][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xy][1]-gfak*R_temp[Cart::zzz][Cart::xy][2])
R_temp[Cart::zzzzz][Cart::yz][1]+=pma2*R_temp[Cart::zzzz][Cart::yz][1]+wmp2*R_temp[Cart::zzzz][Cart::yz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::yz][1]-gfak*R_temp[Cart::zzz][Cart::yz][2])+0.5/_decay*1*R_temp[Cart::zzzz][Cart::y][2]
R_temp[Cart::zzzzz][Cart::xx][1]+=pma2*R_temp[Cart::zzzz][Cart::xx][1]+wmp2*R_temp[Cart::zzzz][Cart::xx][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xx][1]-gfak*R_temp[Cart::zzz][Cart::xx][2])
R_temp[Cart::zzzzz][Cart::xz][1]+=pma2*R_temp[Cart::zzzz][Cart::xz][1]+wmp2*R_temp[Cart::zzzz][Cart::xz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xz][1]-gfak*R_temp[Cart::zzz][Cart::xz][2])+0.5/_decay*1*R_temp[Cart::zzzz][Cart::x][2]
R_temp[Cart::zzzzz][Cart::zz][1]+=pma2*R_temp[Cart::zzzz][Cart::zz][1]+wmp2*R_temp[Cart::zzzz][Cart::zz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::zz][1]-gfak*R_temp[Cart::zzz][Cart::zz][2])+0.5/_decay*2*R_temp[Cart::zzzz][Cart::z][2]
}}
//------------------------------------------------------

//Integral h - s - f - m1
if (_mmax >8 ){
if (_lmax_alpha>4 && _lmax_gamma>2){

R_temp[Cart::yyyyy][Cart::yyy][1]+=pma1*R_temp[Cart::yyyy][Cart::yyy][1]+wmp1*R_temp[Cart::yyyy][Cart::yyy][2]+3*rzeta*(R_temp[Cart::yyy][Cart::yyy][1]-gfak*R_temp[Cart::yyy][Cart::yyy][2])+0.5/_decay*3*R_temp[Cart::yyyy][Cart::yy][2]
R_temp[Cart::yyyyy][Cart::xyy][1]+=pma1*R_temp[Cart::yyyy][Cart::xyy][1]+wmp1*R_temp[Cart::yyyy][Cart::xyy][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xyy][1]-gfak*R_temp[Cart::yyy][Cart::xyy][2])+0.5/_decay*2*R_temp[Cart::yyyy][Cart::xy][2]
R_temp[Cart::yyyyy][Cart::yyz][1]+=pma1*R_temp[Cart::yyyy][Cart::yyz][1]+wmp1*R_temp[Cart::yyyy][Cart::yyz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::yyz][1]-gfak*R_temp[Cart::yyy][Cart::yyz][2])+0.5/_decay*2*R_temp[Cart::yyyy][Cart::yz][2]
R_temp[Cart::yyyyy][Cart::xxy][1]+=pma1*R_temp[Cart::yyyy][Cart::xxy][1]+wmp1*R_temp[Cart::yyyy][Cart::xxy][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xxy][1]-gfak*R_temp[Cart::yyy][Cart::xxy][2])+0.5/_decay*1*R_temp[Cart::yyyy][Cart::xx][2]
R_temp[Cart::yyyyy][Cart::xyz][1]+=pma1*R_temp[Cart::yyyy][Cart::xyz][1]+wmp1*R_temp[Cart::yyyy][Cart::xyz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xyz][1]-gfak*R_temp[Cart::yyy][Cart::xyz][2])+0.5/_decay*1*R_temp[Cart::yyyy][Cart::xz][2]
R_temp[Cart::yyyyy][Cart::yzz][1]+=pma1*R_temp[Cart::yyyy][Cart::yzz][1]+wmp1*R_temp[Cart::yyyy][Cart::yzz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::yzz][1]-gfak*R_temp[Cart::yyy][Cart::yzz][2])+0.5/_decay*1*R_temp[Cart::yyyy][Cart::zz][2]
R_temp[Cart::yyyyy][Cart::xxx][1]+=pma1*R_temp[Cart::yyyy][Cart::xxx][1]+wmp1*R_temp[Cart::yyyy][Cart::xxx][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xxx][1]-gfak*R_temp[Cart::yyy][Cart::xxx][2])
R_temp[Cart::yyyyy][Cart::xxz][1]+=pma1*R_temp[Cart::yyyy][Cart::xxz][1]+wmp1*R_temp[Cart::yyyy][Cart::xxz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xxz][1]-gfak*R_temp[Cart::yyy][Cart::xxz][2])
R_temp[Cart::yyyyy][Cart::xzz][1]+=pma1*R_temp[Cart::yyyy][Cart::xzz][1]+wmp1*R_temp[Cart::yyyy][Cart::xzz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::xzz][1]-gfak*R_temp[Cart::yyy][Cart::xzz][2])
R_temp[Cart::yyyyy][Cart::zzz][1]+=pma1*R_temp[Cart::yyyy][Cart::zzz][1]+wmp1*R_temp[Cart::yyyy][Cart::zzz][2]+3*rzeta*(R_temp[Cart::yyy][Cart::zzz][1]-gfak*R_temp[Cart::yyy][Cart::zzz][2])
R_temp[Cart::xyyyy][Cart::yyy][1]+=pma0*R_temp[Cart::yyyy][Cart::yyy][1]+wmp0*R_temp[Cart::yyyy][Cart::yyy][2]
R_temp[Cart::xyyyy][Cart::xyy][1]+=pma0*R_temp[Cart::yyyy][Cart::xyy][1]+wmp0*R_temp[Cart::yyyy][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::yy][2]
R_temp[Cart::xyyyy][Cart::yyz][1]+=pma0*R_temp[Cart::yyyy][Cart::yyz][1]+wmp0*R_temp[Cart::yyyy][Cart::yyz][2]
R_temp[Cart::xyyyy][Cart::xxy][1]+=pma0*R_temp[Cart::yyyy][Cart::xxy][1]+wmp0*R_temp[Cart::yyyy][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::yyyy][Cart::xy][2]
R_temp[Cart::xyyyy][Cart::xyz][1]+=pma0*R_temp[Cart::yyyy][Cart::xyz][1]+wmp0*R_temp[Cart::yyyy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::yz][2]
R_temp[Cart::xyyyy][Cart::yzz][1]+=pma0*R_temp[Cart::yyyy][Cart::yzz][1]+wmp0*R_temp[Cart::yyyy][Cart::yzz][2]
R_temp[Cart::xyyyy][Cart::xxx][1]+=pma0*R_temp[Cart::yyyy][Cart::xxx][1]+wmp0*R_temp[Cart::yyyy][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::yyyy][Cart::xx][2]
R_temp[Cart::xyyyy][Cart::xxz][1]+=pma0*R_temp[Cart::yyyy][Cart::xxz][1]+wmp0*R_temp[Cart::yyyy][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::yyyy][Cart::xz][2]
R_temp[Cart::xyyyy][Cart::xzz][1]+=pma0*R_temp[Cart::yyyy][Cart::xzz][1]+wmp0*R_temp[Cart::yyyy][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::zz][2]
R_temp[Cart::xyyyy][Cart::zzz][1]+=pma0*R_temp[Cart::yyyy][Cart::zzz][1]+wmp0*R_temp[Cart::yyyy][Cart::zzz][2]
R_temp[Cart::yyyyz][Cart::yyy][1]+=pma2*R_temp[Cart::yyyy][Cart::yyy][1]+wmp2*R_temp[Cart::yyyy][Cart::yyy][2]
R_temp[Cart::yyyyz][Cart::xyy][1]+=pma2*R_temp[Cart::yyyy][Cart::xyy][1]+wmp2*R_temp[Cart::yyyy][Cart::xyy][2]
R_temp[Cart::yyyyz][Cart::yyz][1]+=pma2*R_temp[Cart::yyyy][Cart::yyz][1]+wmp2*R_temp[Cart::yyyy][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::yy][2]
R_temp[Cart::yyyyz][Cart::xxy][1]+=pma2*R_temp[Cart::yyyy][Cart::xxy][1]+wmp2*R_temp[Cart::yyyy][Cart::xxy][2]
R_temp[Cart::yyyyz][Cart::xyz][1]+=pma2*R_temp[Cart::yyyy][Cart::xyz][1]+wmp2*R_temp[Cart::yyyy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::xy][2]
R_temp[Cart::yyyyz][Cart::yzz][1]+=pma2*R_temp[Cart::yyyy][Cart::yzz][1]+wmp2*R_temp[Cart::yyyy][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::yyyy][Cart::yz][2]
R_temp[Cart::yyyyz][Cart::xxx][1]+=pma2*R_temp[Cart::yyyy][Cart::xxx][1]+wmp2*R_temp[Cart::yyyy][Cart::xxx][2]
R_temp[Cart::yyyyz][Cart::xxz][1]+=pma2*R_temp[Cart::yyyy][Cart::xxz][1]+wmp2*R_temp[Cart::yyyy][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::yyyy][Cart::xx][2]
R_temp[Cart::yyyyz][Cart::xzz][1]+=pma2*R_temp[Cart::yyyy][Cart::xzz][1]+wmp2*R_temp[Cart::yyyy][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::yyyy][Cart::xz][2]
R_temp[Cart::yyyyz][Cart::zzz][1]+=pma2*R_temp[Cart::yyyy][Cart::zzz][1]+wmp2*R_temp[Cart::yyyy][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::yyyy][Cart::zz][2]
R_temp[Cart::xxyyy][Cart::yyy][1]+=pma0*R_temp[Cart::xyyy][Cart::yyy][1]+wmp0*R_temp[Cart::xyyy][Cart::yyy][2]
R_temp[Cart::xxyyy][Cart::xyy][1]+=pma0*R_temp[Cart::xyyy][Cart::xyy][1]+wmp0*R_temp[Cart::xyyy][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::xyyy][Cart::yy][2]
R_temp[Cart::xxyyy][Cart::yyz][1]+=pma0*R_temp[Cart::xyyy][Cart::yyz][1]+wmp0*R_temp[Cart::xyyy][Cart::yyz][2]
R_temp[Cart::xxyyy][Cart::xxy][1]+=pma0*R_temp[Cart::xyyy][Cart::xxy][1]+wmp0*R_temp[Cart::xyyy][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::xyyy][Cart::xy][2]
R_temp[Cart::xxyyy][Cart::xyz][1]+=pma0*R_temp[Cart::xyyy][Cart::xyz][1]+wmp0*R_temp[Cart::xyyy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xyyy][Cart::yz][2]
R_temp[Cart::xxyyy][Cart::yzz][1]+=pma0*R_temp[Cart::xyyy][Cart::yzz][1]+wmp0*R_temp[Cart::xyyy][Cart::yzz][2]
R_temp[Cart::xxyyy][Cart::xxx][1]+=pma0*R_temp[Cart::xyyy][Cart::xxx][1]+wmp0*R_temp[Cart::xyyy][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::xyyy][Cart::xx][2]
R_temp[Cart::xxyyy][Cart::xxz][1]+=pma0*R_temp[Cart::xyyy][Cart::xxz][1]+wmp0*R_temp[Cart::xyyy][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::xyyy][Cart::xz][2]
R_temp[Cart::xxyyy][Cart::xzz][1]+=pma0*R_temp[Cart::xyyy][Cart::xzz][1]+wmp0*R_temp[Cart::xyyy][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::xyyy][Cart::zz][2]
R_temp[Cart::xxyyy][Cart::zzz][1]+=pma0*R_temp[Cart::xyyy][Cart::zzz][1]+wmp0*R_temp[Cart::xyyy][Cart::zzz][2]
R_temp[Cart::xyyyz][Cart::yyy][1]+=pma0*R_temp[Cart::yyyz][Cart::yyy][1]+wmp0*R_temp[Cart::yyyz][Cart::yyy][2]
R_temp[Cart::xyyyz][Cart::xyy][1]+=pma0*R_temp[Cart::yyyz][Cart::xyy][1]+wmp0*R_temp[Cart::yyyz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::yy][2]
R_temp[Cart::xyyyz][Cart::yyz][1]+=pma0*R_temp[Cart::yyyz][Cart::yyz][1]+wmp0*R_temp[Cart::yyyz][Cart::yyz][2]
R_temp[Cart::xyyyz][Cart::xxy][1]+=pma0*R_temp[Cart::yyyz][Cart::xxy][1]+wmp0*R_temp[Cart::yyyz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::yyyz][Cart::xy][2]
R_temp[Cart::xyyyz][Cart::xyz][1]+=pma0*R_temp[Cart::yyyz][Cart::xyz][1]+wmp0*R_temp[Cart::yyyz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::yz][2]
R_temp[Cart::xyyyz][Cart::yzz][1]+=pma0*R_temp[Cart::yyyz][Cart::yzz][1]+wmp0*R_temp[Cart::yyyz][Cart::yzz][2]
R_temp[Cart::xyyyz][Cart::xxx][1]+=pma0*R_temp[Cart::yyyz][Cart::xxx][1]+wmp0*R_temp[Cart::yyyz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::yyyz][Cart::xx][2]
R_temp[Cart::xyyyz][Cart::xxz][1]+=pma0*R_temp[Cart::yyyz][Cart::xxz][1]+wmp0*R_temp[Cart::yyyz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::yyyz][Cart::xz][2]
R_temp[Cart::xyyyz][Cart::xzz][1]+=pma0*R_temp[Cart::yyyz][Cart::xzz][1]+wmp0*R_temp[Cart::yyyz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::zz][2]
R_temp[Cart::xyyyz][Cart::zzz][1]+=pma0*R_temp[Cart::yyyz][Cart::zzz][1]+wmp0*R_temp[Cart::yyyz][Cart::zzz][2]
R_temp[Cart::yyyzz][Cart::yyy][1]+=pma2*R_temp[Cart::yyyz][Cart::yyy][1]+wmp2*R_temp[Cart::yyyz][Cart::yyy][2]
R_temp[Cart::yyyzz][Cart::xyy][1]+=pma2*R_temp[Cart::yyyz][Cart::xyy][1]+wmp2*R_temp[Cart::yyyz][Cart::xyy][2]
R_temp[Cart::yyyzz][Cart::yyz][1]+=pma2*R_temp[Cart::yyyz][Cart::yyz][1]+wmp2*R_temp[Cart::yyyz][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::yy][2]
R_temp[Cart::yyyzz][Cart::xxy][1]+=pma2*R_temp[Cart::yyyz][Cart::xxy][1]+wmp2*R_temp[Cart::yyyz][Cart::xxy][2]
R_temp[Cart::yyyzz][Cart::xyz][1]+=pma2*R_temp[Cart::yyyz][Cart::xyz][1]+wmp2*R_temp[Cart::yyyz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::xy][2]
R_temp[Cart::yyyzz][Cart::yzz][1]+=pma2*R_temp[Cart::yyyz][Cart::yzz][1]+wmp2*R_temp[Cart::yyyz][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::yyyz][Cart::yz][2]
R_temp[Cart::yyyzz][Cart::xxx][1]+=pma2*R_temp[Cart::yyyz][Cart::xxx][1]+wmp2*R_temp[Cart::yyyz][Cart::xxx][2]
R_temp[Cart::yyyzz][Cart::xxz][1]+=pma2*R_temp[Cart::yyyz][Cart::xxz][1]+wmp2*R_temp[Cart::yyyz][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::yyyz][Cart::xx][2]
R_temp[Cart::yyyzz][Cart::xzz][1]+=pma2*R_temp[Cart::yyyz][Cart::xzz][1]+wmp2*R_temp[Cart::yyyz][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::yyyz][Cart::xz][2]
R_temp[Cart::yyyzz][Cart::zzz][1]+=pma2*R_temp[Cart::yyyz][Cart::zzz][1]+wmp2*R_temp[Cart::yyyz][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::yyyz][Cart::zz][2]
R_temp[Cart::xxxyy][Cart::yyy][1]+=pma1*R_temp[Cart::xxxy][Cart::yyy][1]+wmp1*R_temp[Cart::xxxy][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::xxxy][Cart::yy][2]
R_temp[Cart::xxxyy][Cart::xyy][1]+=pma1*R_temp[Cart::xxxy][Cart::xyy][1]+wmp1*R_temp[Cart::xxxy][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::xxxy][Cart::xy][2]
R_temp[Cart::xxxyy][Cart::yyz][1]+=pma1*R_temp[Cart::xxxy][Cart::yyz][1]+wmp1*R_temp[Cart::xxxy][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::xxxy][Cart::yz][2]
R_temp[Cart::xxxyy][Cart::xxy][1]+=pma1*R_temp[Cart::xxxy][Cart::xxy][1]+wmp1*R_temp[Cart::xxxy][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::xxxy][Cart::xx][2]
R_temp[Cart::xxxyy][Cart::xyz][1]+=pma1*R_temp[Cart::xxxy][Cart::xyz][1]+wmp1*R_temp[Cart::xxxy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxxy][Cart::xz][2]
R_temp[Cart::xxxyy][Cart::yzz][1]+=pma1*R_temp[Cart::xxxy][Cart::yzz][1]+wmp1*R_temp[Cart::xxxy][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::xxxy][Cart::zz][2]
R_temp[Cart::xxxyy][Cart::xxx][1]+=pma1*R_temp[Cart::xxxy][Cart::xxx][1]+wmp1*R_temp[Cart::xxxy][Cart::xxx][2]
R_temp[Cart::xxxyy][Cart::xxz][1]+=pma1*R_temp[Cart::xxxy][Cart::xxz][1]+wmp1*R_temp[Cart::xxxy][Cart::xxz][2]
R_temp[Cart::xxxyy][Cart::xzz][1]+=pma1*R_temp[Cart::xxxy][Cart::xzz][1]+wmp1*R_temp[Cart::xxxy][Cart::xzz][2]
R_temp[Cart::xxxyy][Cart::zzz][1]+=pma1*R_temp[Cart::xxxy][Cart::zzz][1]+wmp1*R_temp[Cart::xxxy][Cart::zzz][2]
R_temp[Cart::xxyyz][Cart::yyy][1]+=pma2*R_temp[Cart::xxyy][Cart::yyy][1]+wmp2*R_temp[Cart::xxyy][Cart::yyy][2]
R_temp[Cart::xxyyz][Cart::xyy][1]+=pma2*R_temp[Cart::xxyy][Cart::xyy][1]+wmp2*R_temp[Cart::xxyy][Cart::xyy][2]
R_temp[Cart::xxyyz][Cart::yyz][1]+=pma2*R_temp[Cart::xxyy][Cart::yyz][1]+wmp2*R_temp[Cart::xxyy][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::xxyy][Cart::yy][2]
R_temp[Cart::xxyyz][Cart::xxy][1]+=pma2*R_temp[Cart::xxyy][Cart::xxy][1]+wmp2*R_temp[Cart::xxyy][Cart::xxy][2]
R_temp[Cart::xxyyz][Cart::xyz][1]+=pma2*R_temp[Cart::xxyy][Cart::xyz][1]+wmp2*R_temp[Cart::xxyy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxyy][Cart::xy][2]
R_temp[Cart::xxyyz][Cart::yzz][1]+=pma2*R_temp[Cart::xxyy][Cart::yzz][1]+wmp2*R_temp[Cart::xxyy][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::xxyy][Cart::yz][2]
R_temp[Cart::xxyyz][Cart::xxx][1]+=pma2*R_temp[Cart::xxyy][Cart::xxx][1]+wmp2*R_temp[Cart::xxyy][Cart::xxx][2]
R_temp[Cart::xxyyz][Cart::xxz][1]+=pma2*R_temp[Cart::xxyy][Cart::xxz][1]+wmp2*R_temp[Cart::xxyy][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::xxyy][Cart::xx][2]
R_temp[Cart::xxyyz][Cart::xzz][1]+=pma2*R_temp[Cart::xxyy][Cart::xzz][1]+wmp2*R_temp[Cart::xxyy][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::xxyy][Cart::xz][2]
R_temp[Cart::xxyyz][Cart::zzz][1]+=pma2*R_temp[Cart::xxyy][Cart::zzz][1]+wmp2*R_temp[Cart::xxyy][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::xxyy][Cart::zz][2]
R_temp[Cart::xyyzz][Cart::yyy][1]+=pma0*R_temp[Cart::yyzz][Cart::yyy][1]+wmp0*R_temp[Cart::yyzz][Cart::yyy][2]
R_temp[Cart::xyyzz][Cart::xyy][1]+=pma0*R_temp[Cart::yyzz][Cart::xyy][1]+wmp0*R_temp[Cart::yyzz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::yyzz][Cart::yy][2]
R_temp[Cart::xyyzz][Cart::yyz][1]+=pma0*R_temp[Cart::yyzz][Cart::yyz][1]+wmp0*R_temp[Cart::yyzz][Cart::yyz][2]
R_temp[Cart::xyyzz][Cart::xxy][1]+=pma0*R_temp[Cart::yyzz][Cart::xxy][1]+wmp0*R_temp[Cart::yyzz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::yyzz][Cart::xy][2]
R_temp[Cart::xyyzz][Cart::xyz][1]+=pma0*R_temp[Cart::yyzz][Cart::xyz][1]+wmp0*R_temp[Cart::yyzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yyzz][Cart::yz][2]
R_temp[Cart::xyyzz][Cart::yzz][1]+=pma0*R_temp[Cart::yyzz][Cart::yzz][1]+wmp0*R_temp[Cart::yyzz][Cart::yzz][2]
R_temp[Cart::xyyzz][Cart::xxx][1]+=pma0*R_temp[Cart::yyzz][Cart::xxx][1]+wmp0*R_temp[Cart::yyzz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::yyzz][Cart::xx][2]
R_temp[Cart::xyyzz][Cart::xxz][1]+=pma0*R_temp[Cart::yyzz][Cart::xxz][1]+wmp0*R_temp[Cart::yyzz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::yyzz][Cart::xz][2]
R_temp[Cart::xyyzz][Cart::xzz][1]+=pma0*R_temp[Cart::yyzz][Cart::xzz][1]+wmp0*R_temp[Cart::yyzz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::yyzz][Cart::zz][2]
R_temp[Cart::xyyzz][Cart::zzz][1]+=pma0*R_temp[Cart::yyzz][Cart::zzz][1]+wmp0*R_temp[Cart::yyzz][Cart::zzz][2]
R_temp[Cart::yyzzz][Cart::yyy][1]+=pma1*R_temp[Cart::yzzz][Cart::yyy][1]+wmp1*R_temp[Cart::yzzz][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::yzzz][Cart::yy][2]
R_temp[Cart::yyzzz][Cart::xyy][1]+=pma1*R_temp[Cart::yzzz][Cart::xyy][1]+wmp1*R_temp[Cart::yzzz][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::yzzz][Cart::xy][2]
R_temp[Cart::yyzzz][Cart::yyz][1]+=pma1*R_temp[Cart::yzzz][Cart::yyz][1]+wmp1*R_temp[Cart::yzzz][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::yzzz][Cart::yz][2]
R_temp[Cart::yyzzz][Cart::xxy][1]+=pma1*R_temp[Cart::yzzz][Cart::xxy][1]+wmp1*R_temp[Cart::yzzz][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::xx][2]
R_temp[Cart::yyzzz][Cart::xyz][1]+=pma1*R_temp[Cart::yzzz][Cart::xyz][1]+wmp1*R_temp[Cart::yzzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::xz][2]
R_temp[Cart::yyzzz][Cart::yzz][1]+=pma1*R_temp[Cart::yzzz][Cart::yzz][1]+wmp1*R_temp[Cart::yzzz][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::zz][2]
R_temp[Cart::yyzzz][Cart::xxx][1]+=pma1*R_temp[Cart::yzzz][Cart::xxx][1]+wmp1*R_temp[Cart::yzzz][Cart::xxx][2]
R_temp[Cart::yyzzz][Cart::xxz][1]+=pma1*R_temp[Cart::yzzz][Cart::xxz][1]+wmp1*R_temp[Cart::yzzz][Cart::xxz][2]
R_temp[Cart::yyzzz][Cart::xzz][1]+=pma1*R_temp[Cart::yzzz][Cart::xzz][1]+wmp1*R_temp[Cart::yzzz][Cart::xzz][2]
R_temp[Cart::yyzzz][Cart::zzz][1]+=pma1*R_temp[Cart::yzzz][Cart::zzz][1]+wmp1*R_temp[Cart::yzzz][Cart::zzz][2]
R_temp[Cart::xxxxy][Cart::yyy][1]+=pma1*R_temp[Cart::xxxx][Cart::yyy][1]+wmp1*R_temp[Cart::xxxx][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::xxxx][Cart::yy][2]
R_temp[Cart::xxxxy][Cart::xyy][1]+=pma1*R_temp[Cart::xxxx][Cart::xyy][1]+wmp1*R_temp[Cart::xxxx][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::xxxx][Cart::xy][2]
R_temp[Cart::xxxxy][Cart::yyz][1]+=pma1*R_temp[Cart::xxxx][Cart::yyz][1]+wmp1*R_temp[Cart::xxxx][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::xxxx][Cart::yz][2]
R_temp[Cart::xxxxy][Cart::xxy][1]+=pma1*R_temp[Cart::xxxx][Cart::xxy][1]+wmp1*R_temp[Cart::xxxx][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::xx][2]
R_temp[Cart::xxxxy][Cart::xyz][1]+=pma1*R_temp[Cart::xxxx][Cart::xyz][1]+wmp1*R_temp[Cart::xxxx][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::xz][2]
R_temp[Cart::xxxxy][Cart::yzz][1]+=pma1*R_temp[Cart::xxxx][Cart::yzz][1]+wmp1*R_temp[Cart::xxxx][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::zz][2]
R_temp[Cart::xxxxy][Cart::xxx][1]+=pma1*R_temp[Cart::xxxx][Cart::xxx][1]+wmp1*R_temp[Cart::xxxx][Cart::xxx][2]
R_temp[Cart::xxxxy][Cart::xxz][1]+=pma1*R_temp[Cart::xxxx][Cart::xxz][1]+wmp1*R_temp[Cart::xxxx][Cart::xxz][2]
R_temp[Cart::xxxxy][Cart::xzz][1]+=pma1*R_temp[Cart::xxxx][Cart::xzz][1]+wmp1*R_temp[Cart::xxxx][Cart::xzz][2]
R_temp[Cart::xxxxy][Cart::zzz][1]+=pma1*R_temp[Cart::xxxx][Cart::zzz][1]+wmp1*R_temp[Cart::xxxx][Cart::zzz][2]
R_temp[Cart::xxxyz][Cart::yyy][1]+=pma1*R_temp[Cart::xxxz][Cart::yyy][1]+wmp1*R_temp[Cart::xxxz][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::xxxz][Cart::yy][2]
R_temp[Cart::xxxyz][Cart::xyy][1]+=pma1*R_temp[Cart::xxxz][Cart::xyy][1]+wmp1*R_temp[Cart::xxxz][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::xxxz][Cart::xy][2]
R_temp[Cart::xxxyz][Cart::yyz][1]+=pma1*R_temp[Cart::xxxz][Cart::yyz][1]+wmp1*R_temp[Cart::xxxz][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::xxxz][Cart::yz][2]
R_temp[Cart::xxxyz][Cart::xxy][1]+=pma1*R_temp[Cart::xxxz][Cart::xxy][1]+wmp1*R_temp[Cart::xxxz][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::xx][2]
R_temp[Cart::xxxyz][Cart::xyz][1]+=pma1*R_temp[Cart::xxxz][Cart::xyz][1]+wmp1*R_temp[Cart::xxxz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::xz][2]
R_temp[Cart::xxxyz][Cart::yzz][1]+=pma1*R_temp[Cart::xxxz][Cart::yzz][1]+wmp1*R_temp[Cart::xxxz][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::zz][2]
R_temp[Cart::xxxyz][Cart::xxx][1]+=pma1*R_temp[Cart::xxxz][Cart::xxx][1]+wmp1*R_temp[Cart::xxxz][Cart::xxx][2]
R_temp[Cart::xxxyz][Cart::xxz][1]+=pma1*R_temp[Cart::xxxz][Cart::xxz][1]+wmp1*R_temp[Cart::xxxz][Cart::xxz][2]
R_temp[Cart::xxxyz][Cart::xzz][1]+=pma1*R_temp[Cart::xxxz][Cart::xzz][1]+wmp1*R_temp[Cart::xxxz][Cart::xzz][2]
R_temp[Cart::xxxyz][Cart::zzz][1]+=pma1*R_temp[Cart::xxxz][Cart::zzz][1]+wmp1*R_temp[Cart::xxxz][Cart::zzz][2]
R_temp[Cart::xxyzz][Cart::yyy][1]+=pma1*R_temp[Cart::xxzz][Cart::yyy][1]+wmp1*R_temp[Cart::xxzz][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::xxzz][Cart::yy][2]
R_temp[Cart::xxyzz][Cart::xyy][1]+=pma1*R_temp[Cart::xxzz][Cart::xyy][1]+wmp1*R_temp[Cart::xxzz][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::xxzz][Cart::xy][2]
R_temp[Cart::xxyzz][Cart::yyz][1]+=pma1*R_temp[Cart::xxzz][Cart::yyz][1]+wmp1*R_temp[Cart::xxzz][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::xxzz][Cart::yz][2]
R_temp[Cart::xxyzz][Cart::xxy][1]+=pma1*R_temp[Cart::xxzz][Cart::xxy][1]+wmp1*R_temp[Cart::xxzz][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::xxzz][Cart::xx][2]
R_temp[Cart::xxyzz][Cart::xyz][1]+=pma1*R_temp[Cart::xxzz][Cart::xyz][1]+wmp1*R_temp[Cart::xxzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxzz][Cart::xz][2]
R_temp[Cart::xxyzz][Cart::yzz][1]+=pma1*R_temp[Cart::xxzz][Cart::yzz][1]+wmp1*R_temp[Cart::xxzz][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::xxzz][Cart::zz][2]
R_temp[Cart::xxyzz][Cart::xxx][1]+=pma1*R_temp[Cart::xxzz][Cart::xxx][1]+wmp1*R_temp[Cart::xxzz][Cart::xxx][2]
R_temp[Cart::xxyzz][Cart::xxz][1]+=pma1*R_temp[Cart::xxzz][Cart::xxz][1]+wmp1*R_temp[Cart::xxzz][Cart::xxz][2]
R_temp[Cart::xxyzz][Cart::xzz][1]+=pma1*R_temp[Cart::xxzz][Cart::xzz][1]+wmp1*R_temp[Cart::xxzz][Cart::xzz][2]
R_temp[Cart::xxyzz][Cart::zzz][1]+=pma1*R_temp[Cart::xxzz][Cart::zzz][1]+wmp1*R_temp[Cart::xxzz][Cart::zzz][2]
R_temp[Cart::xyzzz][Cart::yyy][1]+=pma0*R_temp[Cart::yzzz][Cart::yyy][1]+wmp0*R_temp[Cart::yzzz][Cart::yyy][2]
R_temp[Cart::xyzzz][Cart::xyy][1]+=pma0*R_temp[Cart::yzzz][Cart::xyy][1]+wmp0*R_temp[Cart::yzzz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::yy][2]
R_temp[Cart::xyzzz][Cart::yyz][1]+=pma0*R_temp[Cart::yzzz][Cart::yyz][1]+wmp0*R_temp[Cart::yzzz][Cart::yyz][2]
R_temp[Cart::xyzzz][Cart::xxy][1]+=pma0*R_temp[Cart::yzzz][Cart::xxy][1]+wmp0*R_temp[Cart::yzzz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::yzzz][Cart::xy][2]
R_temp[Cart::xyzzz][Cart::xyz][1]+=pma0*R_temp[Cart::yzzz][Cart::xyz][1]+wmp0*R_temp[Cart::yzzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::yz][2]
R_temp[Cart::xyzzz][Cart::yzz][1]+=pma0*R_temp[Cart::yzzz][Cart::yzz][1]+wmp0*R_temp[Cart::yzzz][Cart::yzz][2]
R_temp[Cart::xyzzz][Cart::xxx][1]+=pma0*R_temp[Cart::yzzz][Cart::xxx][1]+wmp0*R_temp[Cart::yzzz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::yzzz][Cart::xx][2]
R_temp[Cart::xyzzz][Cart::xxz][1]+=pma0*R_temp[Cart::yzzz][Cart::xxz][1]+wmp0*R_temp[Cart::yzzz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::yzzz][Cart::xz][2]
R_temp[Cart::xyzzz][Cart::xzz][1]+=pma0*R_temp[Cart::yzzz][Cart::xzz][1]+wmp0*R_temp[Cart::yzzz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::yzzz][Cart::zz][2]
R_temp[Cart::xyzzz][Cart::zzz][1]+=pma0*R_temp[Cart::yzzz][Cart::zzz][1]+wmp0*R_temp[Cart::yzzz][Cart::zzz][2]
R_temp[Cart::yzzzz][Cart::yyy][1]+=pma1*R_temp[Cart::zzzz][Cart::yyy][1]+wmp1*R_temp[Cart::zzzz][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::zzzz][Cart::yy][2]
R_temp[Cart::yzzzz][Cart::xyy][1]+=pma1*R_temp[Cart::zzzz][Cart::xyy][1]+wmp1*R_temp[Cart::zzzz][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::zzzz][Cart::xy][2]
R_temp[Cart::yzzzz][Cart::yyz][1]+=pma1*R_temp[Cart::zzzz][Cart::yyz][1]+wmp1*R_temp[Cart::zzzz][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::zzzz][Cart::yz][2]
R_temp[Cart::yzzzz][Cart::xxy][1]+=pma1*R_temp[Cart::zzzz][Cart::xxy][1]+wmp1*R_temp[Cart::zzzz][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::xx][2]
R_temp[Cart::yzzzz][Cart::xyz][1]+=pma1*R_temp[Cart::zzzz][Cart::xyz][1]+wmp1*R_temp[Cart::zzzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::xz][2]
R_temp[Cart::yzzzz][Cart::yzz][1]+=pma1*R_temp[Cart::zzzz][Cart::yzz][1]+wmp1*R_temp[Cart::zzzz][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::zz][2]
R_temp[Cart::yzzzz][Cart::xxx][1]+=pma1*R_temp[Cart::zzzz][Cart::xxx][1]+wmp1*R_temp[Cart::zzzz][Cart::xxx][2]
R_temp[Cart::yzzzz][Cart::xxz][1]+=pma1*R_temp[Cart::zzzz][Cart::xxz][1]+wmp1*R_temp[Cart::zzzz][Cart::xxz][2]
R_temp[Cart::yzzzz][Cart::xzz][1]+=pma1*R_temp[Cart::zzzz][Cart::xzz][1]+wmp1*R_temp[Cart::zzzz][Cart::xzz][2]
R_temp[Cart::yzzzz][Cart::zzz][1]+=pma1*R_temp[Cart::zzzz][Cart::zzz][1]+wmp1*R_temp[Cart::zzzz][Cart::zzz][2]
R_temp[Cart::xxxxx][Cart::yyy][1]+=pma0*R_temp[Cart::xxxx][Cart::yyy][1]+wmp0*R_temp[Cart::xxxx][Cart::yyy][2]+3*rzeta*(R_temp[Cart::xxx][Cart::yyy][1]-gfak*R_temp[Cart::xxx][Cart::yyy][2])
R_temp[Cart::xxxxx][Cart::xyy][1]+=pma0*R_temp[Cart::xxxx][Cart::xyy][1]+wmp0*R_temp[Cart::xxxx][Cart::xyy][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xyy][1]-gfak*R_temp[Cart::xxx][Cart::xyy][2])+0.5/_decay*1*R_temp[Cart::xxxx][Cart::yy][2]
R_temp[Cart::xxxxx][Cart::yyz][1]+=pma0*R_temp[Cart::xxxx][Cart::yyz][1]+wmp0*R_temp[Cart::xxxx][Cart::yyz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::yyz][1]-gfak*R_temp[Cart::xxx][Cart::yyz][2])
R_temp[Cart::xxxxx][Cart::xxy][1]+=pma0*R_temp[Cart::xxxx][Cart::xxy][1]+wmp0*R_temp[Cart::xxxx][Cart::xxy][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xxy][1]-gfak*R_temp[Cart::xxx][Cart::xxy][2])+0.5/_decay*2*R_temp[Cart::xxxx][Cart::xy][2]
R_temp[Cart::xxxxx][Cart::xyz][1]+=pma0*R_temp[Cart::xxxx][Cart::xyz][1]+wmp0*R_temp[Cart::xxxx][Cart::xyz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xyz][1]-gfak*R_temp[Cart::xxx][Cart::xyz][2])+0.5/_decay*1*R_temp[Cart::xxxx][Cart::yz][2]
R_temp[Cart::xxxxx][Cart::yzz][1]+=pma0*R_temp[Cart::xxxx][Cart::yzz][1]+wmp0*R_temp[Cart::xxxx][Cart::yzz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::yzz][1]-gfak*R_temp[Cart::xxx][Cart::yzz][2])
R_temp[Cart::xxxxx][Cart::xxx][1]+=pma0*R_temp[Cart::xxxx][Cart::xxx][1]+wmp0*R_temp[Cart::xxxx][Cart::xxx][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xxx][1]-gfak*R_temp[Cart::xxx][Cart::xxx][2])+0.5/_decay*3*R_temp[Cart::xxxx][Cart::xx][2]
R_temp[Cart::xxxxx][Cart::xxz][1]+=pma0*R_temp[Cart::xxxx][Cart::xxz][1]+wmp0*R_temp[Cart::xxxx][Cart::xxz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xxz][1]-gfak*R_temp[Cart::xxx][Cart::xxz][2])+0.5/_decay*2*R_temp[Cart::xxxx][Cart::xz][2]
R_temp[Cart::xxxxx][Cart::xzz][1]+=pma0*R_temp[Cart::xxxx][Cart::xzz][1]+wmp0*R_temp[Cart::xxxx][Cart::xzz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::xzz][1]-gfak*R_temp[Cart::xxx][Cart::xzz][2])+0.5/_decay*1*R_temp[Cart::xxxx][Cart::zz][2]
R_temp[Cart::xxxxx][Cart::zzz][1]+=pma0*R_temp[Cart::xxxx][Cart::zzz][1]+wmp0*R_temp[Cart::xxxx][Cart::zzz][2]+3*rzeta*(R_temp[Cart::xxx][Cart::zzz][1]-gfak*R_temp[Cart::xxx][Cart::zzz][2])
R_temp[Cart::xxxxz][Cart::yyy][1]+=pma2*R_temp[Cart::xxxx][Cart::yyy][1]+wmp2*R_temp[Cart::xxxx][Cart::yyy][2]
R_temp[Cart::xxxxz][Cart::xyy][1]+=pma2*R_temp[Cart::xxxx][Cart::xyy][1]+wmp2*R_temp[Cart::xxxx][Cart::xyy][2]
R_temp[Cart::xxxxz][Cart::yyz][1]+=pma2*R_temp[Cart::xxxx][Cart::yyz][1]+wmp2*R_temp[Cart::xxxx][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::yy][2]
R_temp[Cart::xxxxz][Cart::xxy][1]+=pma2*R_temp[Cart::xxxx][Cart::xxy][1]+wmp2*R_temp[Cart::xxxx][Cart::xxy][2]
R_temp[Cart::xxxxz][Cart::xyz][1]+=pma2*R_temp[Cart::xxxx][Cart::xyz][1]+wmp2*R_temp[Cart::xxxx][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::xy][2]
R_temp[Cart::xxxxz][Cart::yzz][1]+=pma2*R_temp[Cart::xxxx][Cart::yzz][1]+wmp2*R_temp[Cart::xxxx][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::xxxx][Cart::yz][2]
R_temp[Cart::xxxxz][Cart::xxx][1]+=pma2*R_temp[Cart::xxxx][Cart::xxx][1]+wmp2*R_temp[Cart::xxxx][Cart::xxx][2]
R_temp[Cart::xxxxz][Cart::xxz][1]+=pma2*R_temp[Cart::xxxx][Cart::xxz][1]+wmp2*R_temp[Cart::xxxx][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::xxxx][Cart::xx][2]
R_temp[Cart::xxxxz][Cart::xzz][1]+=pma2*R_temp[Cart::xxxx][Cart::xzz][1]+wmp2*R_temp[Cart::xxxx][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::xxxx][Cart::xz][2]
R_temp[Cart::xxxxz][Cart::zzz][1]+=pma2*R_temp[Cart::xxxx][Cart::zzz][1]+wmp2*R_temp[Cart::xxxx][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::xxxx][Cart::zz][2]
R_temp[Cart::xxxzz][Cart::yyy][1]+=pma2*R_temp[Cart::xxxz][Cart::yyy][1]+wmp2*R_temp[Cart::xxxz][Cart::yyy][2]
R_temp[Cart::xxxzz][Cart::xyy][1]+=pma2*R_temp[Cart::xxxz][Cart::xyy][1]+wmp2*R_temp[Cart::xxxz][Cart::xyy][2]
R_temp[Cart::xxxzz][Cart::yyz][1]+=pma2*R_temp[Cart::xxxz][Cart::yyz][1]+wmp2*R_temp[Cart::xxxz][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::yy][2]
R_temp[Cart::xxxzz][Cart::xxy][1]+=pma2*R_temp[Cart::xxxz][Cart::xxy][1]+wmp2*R_temp[Cart::xxxz][Cart::xxy][2]
R_temp[Cart::xxxzz][Cart::xyz][1]+=pma2*R_temp[Cart::xxxz][Cart::xyz][1]+wmp2*R_temp[Cart::xxxz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::xy][2]
R_temp[Cart::xxxzz][Cart::yzz][1]+=pma2*R_temp[Cart::xxxz][Cart::yzz][1]+wmp2*R_temp[Cart::xxxz][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::xxxz][Cart::yz][2]
R_temp[Cart::xxxzz][Cart::xxx][1]+=pma2*R_temp[Cart::xxxz][Cart::xxx][1]+wmp2*R_temp[Cart::xxxz][Cart::xxx][2]
R_temp[Cart::xxxzz][Cart::xxz][1]+=pma2*R_temp[Cart::xxxz][Cart::xxz][1]+wmp2*R_temp[Cart::xxxz][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::xxxz][Cart::xx][2]
R_temp[Cart::xxxzz][Cart::xzz][1]+=pma2*R_temp[Cart::xxxz][Cart::xzz][1]+wmp2*R_temp[Cart::xxxz][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::xxxz][Cart::xz][2]
R_temp[Cart::xxxzz][Cart::zzz][1]+=pma2*R_temp[Cart::xxxz][Cart::zzz][1]+wmp2*R_temp[Cart::xxxz][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::xxxz][Cart::zz][2]
R_temp[Cart::xxzzz][Cart::yyy][1]+=pma0*R_temp[Cart::xzzz][Cart::yyy][1]+wmp0*R_temp[Cart::xzzz][Cart::yyy][2]
R_temp[Cart::xxzzz][Cart::xyy][1]+=pma0*R_temp[Cart::xzzz][Cart::xyy][1]+wmp0*R_temp[Cart::xzzz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::xzzz][Cart::yy][2]
R_temp[Cart::xxzzz][Cart::yyz][1]+=pma0*R_temp[Cart::xzzz][Cart::yyz][1]+wmp0*R_temp[Cart::xzzz][Cart::yyz][2]
R_temp[Cart::xxzzz][Cart::xxy][1]+=pma0*R_temp[Cart::xzzz][Cart::xxy][1]+wmp0*R_temp[Cart::xzzz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::xzzz][Cart::xy][2]
R_temp[Cart::xxzzz][Cart::xyz][1]+=pma0*R_temp[Cart::xzzz][Cart::xyz][1]+wmp0*R_temp[Cart::xzzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xzzz][Cart::yz][2]
R_temp[Cart::xxzzz][Cart::yzz][1]+=pma0*R_temp[Cart::xzzz][Cart::yzz][1]+wmp0*R_temp[Cart::xzzz][Cart::yzz][2]
R_temp[Cart::xxzzz][Cart::xxx][1]+=pma0*R_temp[Cart::xzzz][Cart::xxx][1]+wmp0*R_temp[Cart::xzzz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::xzzz][Cart::xx][2]
R_temp[Cart::xxzzz][Cart::xxz][1]+=pma0*R_temp[Cart::xzzz][Cart::xxz][1]+wmp0*R_temp[Cart::xzzz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::xzzz][Cart::xz][2]
R_temp[Cart::xxzzz][Cart::xzz][1]+=pma0*R_temp[Cart::xzzz][Cart::xzz][1]+wmp0*R_temp[Cart::xzzz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::xzzz][Cart::zz][2]
R_temp[Cart::xxzzz][Cart::zzz][1]+=pma0*R_temp[Cart::xzzz][Cart::zzz][1]+wmp0*R_temp[Cart::xzzz][Cart::zzz][2]
R_temp[Cart::xzzzz][Cart::yyy][1]+=pma0*R_temp[Cart::zzzz][Cart::yyy][1]+wmp0*R_temp[Cart::zzzz][Cart::yyy][2]
R_temp[Cart::xzzzz][Cart::xyy][1]+=pma0*R_temp[Cart::zzzz][Cart::xyy][1]+wmp0*R_temp[Cart::zzzz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::yy][2]
R_temp[Cart::xzzzz][Cart::yyz][1]+=pma0*R_temp[Cart::zzzz][Cart::yyz][1]+wmp0*R_temp[Cart::zzzz][Cart::yyz][2]
R_temp[Cart::xzzzz][Cart::xxy][1]+=pma0*R_temp[Cart::zzzz][Cart::xxy][1]+wmp0*R_temp[Cart::zzzz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::zzzz][Cart::xy][2]
R_temp[Cart::xzzzz][Cart::xyz][1]+=pma0*R_temp[Cart::zzzz][Cart::xyz][1]+wmp0*R_temp[Cart::zzzz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::yz][2]
R_temp[Cart::xzzzz][Cart::yzz][1]+=pma0*R_temp[Cart::zzzz][Cart::yzz][1]+wmp0*R_temp[Cart::zzzz][Cart::yzz][2]
R_temp[Cart::xzzzz][Cart::xxx][1]+=pma0*R_temp[Cart::zzzz][Cart::xxx][1]+wmp0*R_temp[Cart::zzzz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::zzzz][Cart::xx][2]
R_temp[Cart::xzzzz][Cart::xxz][1]+=pma0*R_temp[Cart::zzzz][Cart::xxz][1]+wmp0*R_temp[Cart::zzzz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::zzzz][Cart::xz][2]
R_temp[Cart::xzzzz][Cart::xzz][1]+=pma0*R_temp[Cart::zzzz][Cart::xzz][1]+wmp0*R_temp[Cart::zzzz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::zzzz][Cart::zz][2]
R_temp[Cart::xzzzz][Cart::zzz][1]+=pma0*R_temp[Cart::zzzz][Cart::zzz][1]+wmp0*R_temp[Cart::zzzz][Cart::zzz][2]
R_temp[Cart::zzzzz][Cart::yyy][1]+=pma2*R_temp[Cart::zzzz][Cart::yyy][1]+wmp2*R_temp[Cart::zzzz][Cart::yyy][2]+3*rzeta*(R_temp[Cart::zzz][Cart::yyy][1]-gfak*R_temp[Cart::zzz][Cart::yyy][2])
R_temp[Cart::zzzzz][Cart::xyy][1]+=pma2*R_temp[Cart::zzzz][Cart::xyy][1]+wmp2*R_temp[Cart::zzzz][Cart::xyy][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xyy][1]-gfak*R_temp[Cart::zzz][Cart::xyy][2])
R_temp[Cart::zzzzz][Cart::yyz][1]+=pma2*R_temp[Cart::zzzz][Cart::yyz][1]+wmp2*R_temp[Cart::zzzz][Cart::yyz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::yyz][1]-gfak*R_temp[Cart::zzz][Cart::yyz][2])+0.5/_decay*1*R_temp[Cart::zzzz][Cart::yy][2]
R_temp[Cart::zzzzz][Cart::xxy][1]+=pma2*R_temp[Cart::zzzz][Cart::xxy][1]+wmp2*R_temp[Cart::zzzz][Cart::xxy][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xxy][1]-gfak*R_temp[Cart::zzz][Cart::xxy][2])
R_temp[Cart::zzzzz][Cart::xyz][1]+=pma2*R_temp[Cart::zzzz][Cart::xyz][1]+wmp2*R_temp[Cart::zzzz][Cart::xyz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xyz][1]-gfak*R_temp[Cart::zzz][Cart::xyz][2])+0.5/_decay*1*R_temp[Cart::zzzz][Cart::xy][2]
R_temp[Cart::zzzzz][Cart::yzz][1]+=pma2*R_temp[Cart::zzzz][Cart::yzz][1]+wmp2*R_temp[Cart::zzzz][Cart::yzz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::yzz][1]-gfak*R_temp[Cart::zzz][Cart::yzz][2])+0.5/_decay*2*R_temp[Cart::zzzz][Cart::yz][2]
R_temp[Cart::zzzzz][Cart::xxx][1]+=pma2*R_temp[Cart::zzzz][Cart::xxx][1]+wmp2*R_temp[Cart::zzzz][Cart::xxx][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xxx][1]-gfak*R_temp[Cart::zzz][Cart::xxx][2])
R_temp[Cart::zzzzz][Cart::xxz][1]+=pma2*R_temp[Cart::zzzz][Cart::xxz][1]+wmp2*R_temp[Cart::zzzz][Cart::xxz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xxz][1]-gfak*R_temp[Cart::zzz][Cart::xxz][2])+0.5/_decay*1*R_temp[Cart::zzzz][Cart::xx][2]
R_temp[Cart::zzzzz][Cart::xzz][1]+=pma2*R_temp[Cart::zzzz][Cart::xzz][1]+wmp2*R_temp[Cart::zzzz][Cart::xzz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::xzz][1]-gfak*R_temp[Cart::zzz][Cart::xzz][2])+0.5/_decay*2*R_temp[Cart::zzzz][Cart::xz][2]
R_temp[Cart::zzzzz][Cart::zzz][1]+=pma2*R_temp[Cart::zzzz][Cart::zzz][1]+wmp2*R_temp[Cart::zzzz][Cart::zzz][2]+3*rzeta*(R_temp[Cart::zzz][Cart::zzz][1]-gfak*R_temp[Cart::zzz][Cart::zzz][2])+0.5/_decay*3*R_temp[Cart::zzzz][Cart::zz][2]
}}
//------------------------------------------------------

//Integral i - s - s - m0
if (_mmax >6 ){
if (_lmax_alpha>5){

R_temp[Cart::yyyyyy][Cart::s][0]+=pma1*R_temp[Cart::yyyyy][Cart::s][0]+wmp1*R_temp[Cart::yyyyy][Cart::s][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::s][0]-gfak*R_temp[Cart::yyyy][Cart::s][1])
R_temp[Cart::xyyyyy][Cart::s][0]+=pma0*R_temp[Cart::yyyyy][Cart::s][0]+wmp0*R_temp[Cart::yyyyy][Cart::s][1]
R_temp[Cart::yyyyyz][Cart::s][0]+=pma2*R_temp[Cart::yyyyy][Cart::s][0]+wmp2*R_temp[Cart::yyyyy][Cart::s][1]
R_temp[Cart::xxyyyy][Cart::s][0]+=pma0*R_temp[Cart::xyyyy][Cart::s][0]+wmp0*R_temp[Cart::xyyyy][Cart::s][1]
R_temp[Cart::xyyyyz][Cart::s][0]+=pma0*R_temp[Cart::yyyyz][Cart::s][0]+wmp0*R_temp[Cart::yyyyz][Cart::s][1]
R_temp[Cart::yyyyzz][Cart::s][0]+=pma2*R_temp[Cart::yyyyz][Cart::s][0]+wmp2*R_temp[Cart::yyyyz][Cart::s][1]
R_temp[Cart::xxxyyy][Cart::s][0]+=pma0*R_temp[Cart::xxyyy][Cart::s][0]+wmp0*R_temp[Cart::xxyyy][Cart::s][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::s][0]-gfak*R_temp[Cart::xyyy][Cart::s][1])
R_temp[Cart::xxyyyz][Cart::s][0]+=pma2*R_temp[Cart::xxyyy][Cart::s][0]+wmp2*R_temp[Cart::xxyyy][Cart::s][1]
R_temp[Cart::xyyyzz][Cart::s][0]+=pma0*R_temp[Cart::yyyzz][Cart::s][0]+wmp0*R_temp[Cart::yyyzz][Cart::s][1]
R_temp[Cart::yyyzzz][Cart::s][0]+=pma1*R_temp[Cart::yyzzz][Cart::s][0]+wmp1*R_temp[Cart::yyzzz][Cart::s][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::s][0]-gfak*R_temp[Cart::yzzz][Cart::s][1])
R_temp[Cart::xxxxyy][Cart::s][0]+=pma1*R_temp[Cart::xxxxy][Cart::s][0]+wmp1*R_temp[Cart::xxxxy][Cart::s][1]
R_temp[Cart::xxxyyz][Cart::s][0]+=pma2*R_temp[Cart::xxxyy][Cart::s][0]+wmp2*R_temp[Cart::xxxyy][Cart::s][1]
R_temp[Cart::xxyyzz][Cart::s][0]+=pma0*R_temp[Cart::xyyzz][Cart::s][0]+wmp0*R_temp[Cart::xyyzz][Cart::s][1]
R_temp[Cart::xyyzzz][Cart::s][0]+=pma0*R_temp[Cart::yyzzz][Cart::s][0]+wmp0*R_temp[Cart::yyzzz][Cart::s][1]
R_temp[Cart::yyzzzz][Cart::s][0]+=pma1*R_temp[Cart::yzzzz][Cart::s][0]+wmp1*R_temp[Cart::yzzzz][Cart::s][1]
R_temp[Cart::xxxxxy][Cart::s][0]+=pma1*R_temp[Cart::xxxxx][Cart::s][0]+wmp1*R_temp[Cart::xxxxx][Cart::s][1]
R_temp[Cart::xxxxyz][Cart::s][0]+=pma1*R_temp[Cart::xxxxz][Cart::s][0]+wmp1*R_temp[Cart::xxxxz][Cart::s][1]
R_temp[Cart::xxxyzz][Cart::s][0]+=pma1*R_temp[Cart::xxxzz][Cart::s][0]+wmp1*R_temp[Cart::xxxzz][Cart::s][1]
R_temp[Cart::xxyzzz][Cart::s][0]+=pma1*R_temp[Cart::xxzzz][Cart::s][0]+wmp1*R_temp[Cart::xxzzz][Cart::s][1]
R_temp[Cart::xyzzzz][Cart::s][0]+=pma0*R_temp[Cart::yzzzz][Cart::s][0]+wmp0*R_temp[Cart::yzzzz][Cart::s][1]
R_temp[Cart::yzzzzz][Cart::s][0]+=pma1*R_temp[Cart::zzzzz][Cart::s][0]+wmp1*R_temp[Cart::zzzzz][Cart::s][1]
R_temp[Cart::xxxxxx][Cart::s][0]+=pma0*R_temp[Cart::xxxxx][Cart::s][0]+wmp0*R_temp[Cart::xxxxx][Cart::s][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::s][0]-gfak*R_temp[Cart::xxxx][Cart::s][1])
R_temp[Cart::xxxxxz][Cart::s][0]+=pma2*R_temp[Cart::xxxxx][Cart::s][0]+wmp2*R_temp[Cart::xxxxx][Cart::s][1]
R_temp[Cart::xxxxzz][Cart::s][0]+=pma2*R_temp[Cart::xxxxz][Cart::s][0]+wmp2*R_temp[Cart::xxxxz][Cart::s][1]
R_temp[Cart::xxxzzz][Cart::s][0]+=pma0*R_temp[Cart::xxzzz][Cart::s][0]+wmp0*R_temp[Cart::xxzzz][Cart::s][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::s][0]-gfak*R_temp[Cart::xzzz][Cart::s][1])
R_temp[Cart::xxzzzz][Cart::s][0]+=pma0*R_temp[Cart::xzzzz][Cart::s][0]+wmp0*R_temp[Cart::xzzzz][Cart::s][1]
R_temp[Cart::xzzzzz][Cart::s][0]+=pma0*R_temp[Cart::zzzzz][Cart::s][0]+wmp0*R_temp[Cart::zzzzz][Cart::s][1]
R_temp[Cart::zzzzzz][Cart::s][0]+=pma2*R_temp[Cart::zzzzz][Cart::s][0]+wmp2*R_temp[Cart::zzzzz][Cart::s][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::s][0]-gfak*R_temp[Cart::zzzz][Cart::s][1])
}}
//------------------------------------------------------

//Integral i - s - p - m0
if (_mmax >7 ){
if (_lmax_alpha>5 && _lmax_gamma>0){

R_temp[Cart::yyyyyy][Cart::y][0]+=pma1*R_temp[Cart::yyyyy][Cart::y][0]+wmp1*R_temp[Cart::yyyyy][Cart::y][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::y][0]-gfak*R_temp[Cart::yyyy][Cart::y][1])+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::s][1]
R_temp[Cart::yyyyyy][Cart::x][0]+=pma1*R_temp[Cart::yyyyy][Cart::x][0]+wmp1*R_temp[Cart::yyyyy][Cart::x][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::x][0]-gfak*R_temp[Cart::yyyy][Cart::x][1])
R_temp[Cart::yyyyyy][Cart::z][0]+=pma1*R_temp[Cart::yyyyy][Cart::z][0]+wmp1*R_temp[Cart::yyyyy][Cart::z][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::z][0]-gfak*R_temp[Cart::yyyy][Cart::z][1])
R_temp[Cart::xyyyyy][Cart::y][0]+=pma0*R_temp[Cart::yyyyy][Cart::y][0]+wmp0*R_temp[Cart::yyyyy][Cart::y][1]
R_temp[Cart::xyyyyy][Cart::x][0]+=pma0*R_temp[Cart::yyyyy][Cart::x][0]+wmp0*R_temp[Cart::yyyyy][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::s][1]
R_temp[Cart::xyyyyy][Cart::z][0]+=pma0*R_temp[Cart::yyyyy][Cart::z][0]+wmp0*R_temp[Cart::yyyyy][Cart::z][1]
R_temp[Cart::yyyyyz][Cart::y][0]+=pma2*R_temp[Cart::yyyyy][Cart::y][0]+wmp2*R_temp[Cart::yyyyy][Cart::y][1]
R_temp[Cart::yyyyyz][Cart::x][0]+=pma2*R_temp[Cart::yyyyy][Cart::x][0]+wmp2*R_temp[Cart::yyyyy][Cart::x][1]
R_temp[Cart::yyyyyz][Cart::z][0]+=pma2*R_temp[Cart::yyyyy][Cart::z][0]+wmp2*R_temp[Cart::yyyyy][Cart::z][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::s][1]
R_temp[Cart::xxyyyy][Cart::y][0]+=pma0*R_temp[Cart::xyyyy][Cart::y][0]+wmp0*R_temp[Cart::xyyyy][Cart::y][1]
R_temp[Cart::xxyyyy][Cart::x][0]+=pma0*R_temp[Cart::xyyyy][Cart::x][0]+wmp0*R_temp[Cart::xyyyy][Cart::x][1]+0.5/_decay*1*R_temp[Cart::xyyyy][Cart::s][1]
R_temp[Cart::xxyyyy][Cart::z][0]+=pma0*R_temp[Cart::xyyyy][Cart::z][0]+wmp0*R_temp[Cart::xyyyy][Cart::z][1]
R_temp[Cart::xyyyyz][Cart::y][0]+=pma0*R_temp[Cart::yyyyz][Cart::y][0]+wmp0*R_temp[Cart::yyyyz][Cart::y][1]
R_temp[Cart::xyyyyz][Cart::x][0]+=pma0*R_temp[Cart::yyyyz][Cart::x][0]+wmp0*R_temp[Cart::yyyyz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::s][1]
R_temp[Cart::xyyyyz][Cart::z][0]+=pma0*R_temp[Cart::yyyyz][Cart::z][0]+wmp0*R_temp[Cart::yyyyz][Cart::z][1]
R_temp[Cart::yyyyzz][Cart::y][0]+=pma2*R_temp[Cart::yyyyz][Cart::y][0]+wmp2*R_temp[Cart::yyyyz][Cart::y][1]
R_temp[Cart::yyyyzz][Cart::x][0]+=pma2*R_temp[Cart::yyyyz][Cart::x][0]+wmp2*R_temp[Cart::yyyyz][Cart::x][1]
R_temp[Cart::yyyyzz][Cart::z][0]+=pma2*R_temp[Cart::yyyyz][Cart::z][0]+wmp2*R_temp[Cart::yyyyz][Cart::z][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::s][1]
R_temp[Cart::xxxyyy][Cart::y][0]+=pma0*R_temp[Cart::xxyyy][Cart::y][0]+wmp0*R_temp[Cart::xxyyy][Cart::y][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::y][0]-gfak*R_temp[Cart::xyyy][Cart::y][1])
R_temp[Cart::xxxyyy][Cart::x][0]+=pma0*R_temp[Cart::xxyyy][Cart::x][0]+wmp0*R_temp[Cart::xxyyy][Cart::x][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::x][0]-gfak*R_temp[Cart::xyyy][Cart::x][1])+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::s][1]
R_temp[Cart::xxxyyy][Cart::z][0]+=pma0*R_temp[Cart::xxyyy][Cart::z][0]+wmp0*R_temp[Cart::xxyyy][Cart::z][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::z][0]-gfak*R_temp[Cart::xyyy][Cart::z][1])
R_temp[Cart::xxyyyz][Cart::y][0]+=pma2*R_temp[Cart::xxyyy][Cart::y][0]+wmp2*R_temp[Cart::xxyyy][Cart::y][1]
R_temp[Cart::xxyyyz][Cart::x][0]+=pma2*R_temp[Cart::xxyyy][Cart::x][0]+wmp2*R_temp[Cart::xxyyy][Cart::x][1]
R_temp[Cart::xxyyyz][Cart::z][0]+=pma2*R_temp[Cart::xxyyy][Cart::z][0]+wmp2*R_temp[Cart::xxyyy][Cart::z][1]+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::s][1]
R_temp[Cart::xyyyzz][Cart::y][0]+=pma0*R_temp[Cart::yyyzz][Cart::y][0]+wmp0*R_temp[Cart::yyyzz][Cart::y][1]
R_temp[Cart::xyyyzz][Cart::x][0]+=pma0*R_temp[Cart::yyyzz][Cart::x][0]+wmp0*R_temp[Cart::yyyzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yyyzz][Cart::s][1]
R_temp[Cart::xyyyzz][Cart::z][0]+=pma0*R_temp[Cart::yyyzz][Cart::z][0]+wmp0*R_temp[Cart::yyyzz][Cart::z][1]
R_temp[Cart::yyyzzz][Cart::y][0]+=pma1*R_temp[Cart::yyzzz][Cart::y][0]+wmp1*R_temp[Cart::yyzzz][Cart::y][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::y][0]-gfak*R_temp[Cart::yzzz][Cart::y][1])+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::s][1]
R_temp[Cart::yyyzzz][Cart::x][0]+=pma1*R_temp[Cart::yyzzz][Cart::x][0]+wmp1*R_temp[Cart::yyzzz][Cart::x][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::x][0]-gfak*R_temp[Cart::yzzz][Cart::x][1])
R_temp[Cart::yyyzzz][Cart::z][0]+=pma1*R_temp[Cart::yyzzz][Cart::z][0]+wmp1*R_temp[Cart::yyzzz][Cart::z][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::z][0]-gfak*R_temp[Cart::yzzz][Cart::z][1])
R_temp[Cart::xxxxyy][Cart::y][0]+=pma1*R_temp[Cart::xxxxy][Cart::y][0]+wmp1*R_temp[Cart::xxxxy][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxxxy][Cart::s][1]
R_temp[Cart::xxxxyy][Cart::x][0]+=pma1*R_temp[Cart::xxxxy][Cart::x][0]+wmp1*R_temp[Cart::xxxxy][Cart::x][1]
R_temp[Cart::xxxxyy][Cart::z][0]+=pma1*R_temp[Cart::xxxxy][Cart::z][0]+wmp1*R_temp[Cart::xxxxy][Cart::z][1]
R_temp[Cart::xxxyyz][Cart::y][0]+=pma2*R_temp[Cart::xxxyy][Cart::y][0]+wmp2*R_temp[Cart::xxxyy][Cart::y][1]
R_temp[Cart::xxxyyz][Cart::x][0]+=pma2*R_temp[Cart::xxxyy][Cart::x][0]+wmp2*R_temp[Cart::xxxyy][Cart::x][1]
R_temp[Cart::xxxyyz][Cart::z][0]+=pma2*R_temp[Cart::xxxyy][Cart::z][0]+wmp2*R_temp[Cart::xxxyy][Cart::z][1]+0.5/_decay*1*R_temp[Cart::xxxyy][Cart::s][1]
R_temp[Cart::xxyyzz][Cart::y][0]+=pma0*R_temp[Cart::xyyzz][Cart::y][0]+wmp0*R_temp[Cart::xyyzz][Cart::y][1]
R_temp[Cart::xxyyzz][Cart::x][0]+=pma0*R_temp[Cart::xyyzz][Cart::x][0]+wmp0*R_temp[Cart::xyyzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::xyyzz][Cart::s][1]
R_temp[Cart::xxyyzz][Cart::z][0]+=pma0*R_temp[Cart::xyyzz][Cart::z][0]+wmp0*R_temp[Cart::xyyzz][Cart::z][1]
R_temp[Cart::xyyzzz][Cart::y][0]+=pma0*R_temp[Cart::yyzzz][Cart::y][0]+wmp0*R_temp[Cart::yyzzz][Cart::y][1]
R_temp[Cart::xyyzzz][Cart::x][0]+=pma0*R_temp[Cart::yyzzz][Cart::x][0]+wmp0*R_temp[Cart::yyzzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::s][1]
R_temp[Cart::xyyzzz][Cart::z][0]+=pma0*R_temp[Cart::yyzzz][Cart::z][0]+wmp0*R_temp[Cart::yyzzz][Cart::z][1]
R_temp[Cart::yyzzzz][Cart::y][0]+=pma1*R_temp[Cart::yzzzz][Cart::y][0]+wmp1*R_temp[Cart::yzzzz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::s][1]
R_temp[Cart::yyzzzz][Cart::x][0]+=pma1*R_temp[Cart::yzzzz][Cart::x][0]+wmp1*R_temp[Cart::yzzzz][Cart::x][1]
R_temp[Cart::yyzzzz][Cart::z][0]+=pma1*R_temp[Cart::yzzzz][Cart::z][0]+wmp1*R_temp[Cart::yzzzz][Cart::z][1]
R_temp[Cart::xxxxxy][Cart::y][0]+=pma1*R_temp[Cart::xxxxx][Cart::y][0]+wmp1*R_temp[Cart::xxxxx][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::s][1]
R_temp[Cart::xxxxxy][Cart::x][0]+=pma1*R_temp[Cart::xxxxx][Cart::x][0]+wmp1*R_temp[Cart::xxxxx][Cart::x][1]
R_temp[Cart::xxxxxy][Cart::z][0]+=pma1*R_temp[Cart::xxxxx][Cart::z][0]+wmp1*R_temp[Cart::xxxxx][Cart::z][1]
R_temp[Cart::xxxxyz][Cart::y][0]+=pma1*R_temp[Cart::xxxxz][Cart::y][0]+wmp1*R_temp[Cart::xxxxz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::s][1]
R_temp[Cart::xxxxyz][Cart::x][0]+=pma1*R_temp[Cart::xxxxz][Cart::x][0]+wmp1*R_temp[Cart::xxxxz][Cart::x][1]
R_temp[Cart::xxxxyz][Cart::z][0]+=pma1*R_temp[Cart::xxxxz][Cart::z][0]+wmp1*R_temp[Cart::xxxxz][Cart::z][1]
R_temp[Cart::xxxyzz][Cart::y][0]+=pma1*R_temp[Cart::xxxzz][Cart::y][0]+wmp1*R_temp[Cart::xxxzz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxxzz][Cart::s][1]
R_temp[Cart::xxxyzz][Cart::x][0]+=pma1*R_temp[Cart::xxxzz][Cart::x][0]+wmp1*R_temp[Cart::xxxzz][Cart::x][1]
R_temp[Cart::xxxyzz][Cart::z][0]+=pma1*R_temp[Cart::xxxzz][Cart::z][0]+wmp1*R_temp[Cart::xxxzz][Cart::z][1]
R_temp[Cart::xxyzzz][Cart::y][0]+=pma1*R_temp[Cart::xxzzz][Cart::y][0]+wmp1*R_temp[Cart::xxzzz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::s][1]
R_temp[Cart::xxyzzz][Cart::x][0]+=pma1*R_temp[Cart::xxzzz][Cart::x][0]+wmp1*R_temp[Cart::xxzzz][Cart::x][1]
R_temp[Cart::xxyzzz][Cart::z][0]+=pma1*R_temp[Cart::xxzzz][Cart::z][0]+wmp1*R_temp[Cart::xxzzz][Cart::z][1]
R_temp[Cart::xyzzzz][Cart::y][0]+=pma0*R_temp[Cart::yzzzz][Cart::y][0]+wmp0*R_temp[Cart::yzzzz][Cart::y][1]
R_temp[Cart::xyzzzz][Cart::x][0]+=pma0*R_temp[Cart::yzzzz][Cart::x][0]+wmp0*R_temp[Cart::yzzzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::s][1]
R_temp[Cart::xyzzzz][Cart::z][0]+=pma0*R_temp[Cart::yzzzz][Cart::z][0]+wmp0*R_temp[Cart::yzzzz][Cart::z][1]
R_temp[Cart::yzzzzz][Cart::y][0]+=pma1*R_temp[Cart::zzzzz][Cart::y][0]+wmp1*R_temp[Cart::zzzzz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::s][1]
R_temp[Cart::yzzzzz][Cart::x][0]+=pma1*R_temp[Cart::zzzzz][Cart::x][0]+wmp1*R_temp[Cart::zzzzz][Cart::x][1]
R_temp[Cart::yzzzzz][Cart::z][0]+=pma1*R_temp[Cart::zzzzz][Cart::z][0]+wmp1*R_temp[Cart::zzzzz][Cart::z][1]
R_temp[Cart::xxxxxx][Cart::y][0]+=pma0*R_temp[Cart::xxxxx][Cart::y][0]+wmp0*R_temp[Cart::xxxxx][Cart::y][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::y][0]-gfak*R_temp[Cart::xxxx][Cart::y][1])
R_temp[Cart::xxxxxx][Cart::x][0]+=pma0*R_temp[Cart::xxxxx][Cart::x][0]+wmp0*R_temp[Cart::xxxxx][Cart::x][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::x][0]-gfak*R_temp[Cart::xxxx][Cart::x][1])+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::s][1]
R_temp[Cart::xxxxxx][Cart::z][0]+=pma0*R_temp[Cart::xxxxx][Cart::z][0]+wmp0*R_temp[Cart::xxxxx][Cart::z][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::z][0]-gfak*R_temp[Cart::xxxx][Cart::z][1])
R_temp[Cart::xxxxxz][Cart::y][0]+=pma2*R_temp[Cart::xxxxx][Cart::y][0]+wmp2*R_temp[Cart::xxxxx][Cart::y][1]
R_temp[Cart::xxxxxz][Cart::x][0]+=pma2*R_temp[Cart::xxxxx][Cart::x][0]+wmp2*R_temp[Cart::xxxxx][Cart::x][1]
R_temp[Cart::xxxxxz][Cart::z][0]+=pma2*R_temp[Cart::xxxxx][Cart::z][0]+wmp2*R_temp[Cart::xxxxx][Cart::z][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::s][1]
R_temp[Cart::xxxxzz][Cart::y][0]+=pma2*R_temp[Cart::xxxxz][Cart::y][0]+wmp2*R_temp[Cart::xxxxz][Cart::y][1]
R_temp[Cart::xxxxzz][Cart::x][0]+=pma2*R_temp[Cart::xxxxz][Cart::x][0]+wmp2*R_temp[Cart::xxxxz][Cart::x][1]
R_temp[Cart::xxxxzz][Cart::z][0]+=pma2*R_temp[Cart::xxxxz][Cart::z][0]+wmp2*R_temp[Cart::xxxxz][Cart::z][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::s][1]
R_temp[Cart::xxxzzz][Cart::y][0]+=pma0*R_temp[Cart::xxzzz][Cart::y][0]+wmp0*R_temp[Cart::xxzzz][Cart::y][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::y][0]-gfak*R_temp[Cart::xzzz][Cart::y][1])
R_temp[Cart::xxxzzz][Cart::x][0]+=pma0*R_temp[Cart::xxzzz][Cart::x][0]+wmp0*R_temp[Cart::xxzzz][Cart::x][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::x][0]-gfak*R_temp[Cart::xzzz][Cart::x][1])+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::s][1]
R_temp[Cart::xxxzzz][Cart::z][0]+=pma0*R_temp[Cart::xxzzz][Cart::z][0]+wmp0*R_temp[Cart::xxzzz][Cart::z][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::z][0]-gfak*R_temp[Cart::xzzz][Cart::z][1])
R_temp[Cart::xxzzzz][Cart::y][0]+=pma0*R_temp[Cart::xzzzz][Cart::y][0]+wmp0*R_temp[Cart::xzzzz][Cart::y][1]
R_temp[Cart::xxzzzz][Cart::x][0]+=pma0*R_temp[Cart::xzzzz][Cart::x][0]+wmp0*R_temp[Cart::xzzzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::xzzzz][Cart::s][1]
R_temp[Cart::xxzzzz][Cart::z][0]+=pma0*R_temp[Cart::xzzzz][Cart::z][0]+wmp0*R_temp[Cart::xzzzz][Cart::z][1]
R_temp[Cart::xzzzzz][Cart::y][0]+=pma0*R_temp[Cart::zzzzz][Cart::y][0]+wmp0*R_temp[Cart::zzzzz][Cart::y][1]
R_temp[Cart::xzzzzz][Cart::x][0]+=pma0*R_temp[Cart::zzzzz][Cart::x][0]+wmp0*R_temp[Cart::zzzzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::s][1]
R_temp[Cart::xzzzzz][Cart::z][0]+=pma0*R_temp[Cart::zzzzz][Cart::z][0]+wmp0*R_temp[Cart::zzzzz][Cart::z][1]
R_temp[Cart::zzzzzz][Cart::y][0]+=pma2*R_temp[Cart::zzzzz][Cart::y][0]+wmp2*R_temp[Cart::zzzzz][Cart::y][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::y][0]-gfak*R_temp[Cart::zzzz][Cart::y][1])
R_temp[Cart::zzzzzz][Cart::x][0]+=pma2*R_temp[Cart::zzzzz][Cart::x][0]+wmp2*R_temp[Cart::zzzzz][Cart::x][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::x][0]-gfak*R_temp[Cart::zzzz][Cart::x][1])
R_temp[Cart::zzzzzz][Cart::z][0]+=pma2*R_temp[Cart::zzzzz][Cart::z][0]+wmp2*R_temp[Cart::zzzzz][Cart::z][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::z][0]-gfak*R_temp[Cart::zzzz][Cart::z][1])+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::s][1]
}}
//------------------------------------------------------

//Integral i - s - d - m0
if (_mmax >8 ){
if (_lmax_alpha>5 && _lmax_gamma>1){

R_temp[Cart::yyyyyy][Cart::yy][0]+=pma1*R_temp[Cart::yyyyy][Cart::yy][0]+wmp1*R_temp[Cart::yyyyy][Cart::yy][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::yy][0]-gfak*R_temp[Cart::yyyy][Cart::yy][1])+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::y][1]
R_temp[Cart::yyyyyy][Cart::xy][0]+=pma1*R_temp[Cart::yyyyy][Cart::xy][0]+wmp1*R_temp[Cart::yyyyy][Cart::xy][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xy][0]-gfak*R_temp[Cart::yyyy][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::x][1]
R_temp[Cart::yyyyyy][Cart::yz][0]+=pma1*R_temp[Cart::yyyyy][Cart::yz][0]+wmp1*R_temp[Cart::yyyyy][Cart::yz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::yz][0]-gfak*R_temp[Cart::yyyy][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::z][1]
R_temp[Cart::yyyyyy][Cart::xx][0]+=pma1*R_temp[Cart::yyyyy][Cart::xx][0]+wmp1*R_temp[Cart::yyyyy][Cart::xx][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xx][0]-gfak*R_temp[Cart::yyyy][Cart::xx][1])
R_temp[Cart::yyyyyy][Cart::xz][0]+=pma1*R_temp[Cart::yyyyy][Cart::xz][0]+wmp1*R_temp[Cart::yyyyy][Cart::xz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xz][0]-gfak*R_temp[Cart::yyyy][Cart::xz][1])
R_temp[Cart::yyyyyy][Cart::zz][0]+=pma1*R_temp[Cart::yyyyy][Cart::zz][0]+wmp1*R_temp[Cart::yyyyy][Cart::zz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::zz][0]-gfak*R_temp[Cart::yyyy][Cart::zz][1])
R_temp[Cart::xyyyyy][Cart::yy][0]+=pma0*R_temp[Cart::yyyyy][Cart::yy][0]+wmp0*R_temp[Cart::yyyyy][Cart::yy][1]
R_temp[Cart::xyyyyy][Cart::xy][0]+=pma0*R_temp[Cart::yyyyy][Cart::xy][0]+wmp0*R_temp[Cart::yyyyy][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::y][1]
R_temp[Cart::xyyyyy][Cart::yz][0]+=pma0*R_temp[Cart::yyyyy][Cart::yz][0]+wmp0*R_temp[Cart::yyyyy][Cart::yz][1]
R_temp[Cart::xyyyyy][Cart::xx][0]+=pma0*R_temp[Cart::yyyyy][Cart::xx][0]+wmp0*R_temp[Cart::yyyyy][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::x][1]
R_temp[Cart::xyyyyy][Cart::xz][0]+=pma0*R_temp[Cart::yyyyy][Cart::xz][0]+wmp0*R_temp[Cart::yyyyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::z][1]
R_temp[Cart::xyyyyy][Cart::zz][0]+=pma0*R_temp[Cart::yyyyy][Cart::zz][0]+wmp0*R_temp[Cart::yyyyy][Cart::zz][1]
R_temp[Cart::yyyyyz][Cart::yy][0]+=pma2*R_temp[Cart::yyyyy][Cart::yy][0]+wmp2*R_temp[Cart::yyyyy][Cart::yy][1]
R_temp[Cart::yyyyyz][Cart::xy][0]+=pma2*R_temp[Cart::yyyyy][Cart::xy][0]+wmp2*R_temp[Cart::yyyyy][Cart::xy][1]
R_temp[Cart::yyyyyz][Cart::yz][0]+=pma2*R_temp[Cart::yyyyy][Cart::yz][0]+wmp2*R_temp[Cart::yyyyy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::y][1]
R_temp[Cart::yyyyyz][Cart::xx][0]+=pma2*R_temp[Cart::yyyyy][Cart::xx][0]+wmp2*R_temp[Cart::yyyyy][Cart::xx][1]
R_temp[Cart::yyyyyz][Cart::xz][0]+=pma2*R_temp[Cart::yyyyy][Cart::xz][0]+wmp2*R_temp[Cart::yyyyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::x][1]
R_temp[Cart::yyyyyz][Cart::zz][0]+=pma2*R_temp[Cart::yyyyy][Cart::zz][0]+wmp2*R_temp[Cart::yyyyy][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::z][1]
R_temp[Cart::xxyyyy][Cart::yy][0]+=pma0*R_temp[Cart::xyyyy][Cart::yy][0]+wmp0*R_temp[Cart::xyyyy][Cart::yy][1]
R_temp[Cart::xxyyyy][Cart::xy][0]+=pma0*R_temp[Cart::xyyyy][Cart::xy][0]+wmp0*R_temp[Cart::xyyyy][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xyyyy][Cart::y][1]
R_temp[Cart::xxyyyy][Cart::yz][0]+=pma0*R_temp[Cart::xyyyy][Cart::yz][0]+wmp0*R_temp[Cart::xyyyy][Cart::yz][1]
R_temp[Cart::xxyyyy][Cart::xx][0]+=pma0*R_temp[Cart::xyyyy][Cart::xx][0]+wmp0*R_temp[Cart::xyyyy][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::xyyyy][Cart::x][1]
R_temp[Cart::xxyyyy][Cart::xz][0]+=pma0*R_temp[Cart::xyyyy][Cart::xz][0]+wmp0*R_temp[Cart::xyyyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xyyyy][Cart::z][1]
R_temp[Cart::xxyyyy][Cart::zz][0]+=pma0*R_temp[Cart::xyyyy][Cart::zz][0]+wmp0*R_temp[Cart::xyyyy][Cart::zz][1]
R_temp[Cart::xyyyyz][Cart::yy][0]+=pma0*R_temp[Cart::yyyyz][Cart::yy][0]+wmp0*R_temp[Cart::yyyyz][Cart::yy][1]
R_temp[Cart::xyyyyz][Cart::xy][0]+=pma0*R_temp[Cart::yyyyz][Cart::xy][0]+wmp0*R_temp[Cart::yyyyz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::y][1]
R_temp[Cart::xyyyyz][Cart::yz][0]+=pma0*R_temp[Cart::yyyyz][Cart::yz][0]+wmp0*R_temp[Cart::yyyyz][Cart::yz][1]
R_temp[Cart::xyyyyz][Cart::xx][0]+=pma0*R_temp[Cart::yyyyz][Cart::xx][0]+wmp0*R_temp[Cart::yyyyz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yyyyz][Cart::x][1]
R_temp[Cart::xyyyyz][Cart::xz][0]+=pma0*R_temp[Cart::yyyyz][Cart::xz][0]+wmp0*R_temp[Cart::yyyyz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::z][1]
R_temp[Cart::xyyyyz][Cart::zz][0]+=pma0*R_temp[Cart::yyyyz][Cart::zz][0]+wmp0*R_temp[Cart::yyyyz][Cart::zz][1]
R_temp[Cart::yyyyzz][Cart::yy][0]+=pma2*R_temp[Cart::yyyyz][Cart::yy][0]+wmp2*R_temp[Cart::yyyyz][Cart::yy][1]
R_temp[Cart::yyyyzz][Cart::xy][0]+=pma2*R_temp[Cart::yyyyz][Cart::xy][0]+wmp2*R_temp[Cart::yyyyz][Cart::xy][1]
R_temp[Cart::yyyyzz][Cart::yz][0]+=pma2*R_temp[Cart::yyyyz][Cart::yz][0]+wmp2*R_temp[Cart::yyyyz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::y][1]
R_temp[Cart::yyyyzz][Cart::xx][0]+=pma2*R_temp[Cart::yyyyz][Cart::xx][0]+wmp2*R_temp[Cart::yyyyz][Cart::xx][1]
R_temp[Cart::yyyyzz][Cart::xz][0]+=pma2*R_temp[Cart::yyyyz][Cart::xz][0]+wmp2*R_temp[Cart::yyyyz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::x][1]
R_temp[Cart::yyyyzz][Cart::zz][0]+=pma2*R_temp[Cart::yyyyz][Cart::zz][0]+wmp2*R_temp[Cart::yyyyz][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::yyyyz][Cart::z][1]
R_temp[Cart::xxxyyy][Cart::yy][0]+=pma0*R_temp[Cart::xxyyy][Cart::yy][0]+wmp0*R_temp[Cart::xxyyy][Cart::yy][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::yy][0]-gfak*R_temp[Cart::xyyy][Cart::yy][1])
R_temp[Cart::xxxyyy][Cart::xy][0]+=pma0*R_temp[Cart::xxyyy][Cart::xy][0]+wmp0*R_temp[Cart::xxyyy][Cart::xy][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xy][0]-gfak*R_temp[Cart::xyyy][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::y][1]
R_temp[Cart::xxxyyy][Cart::yz][0]+=pma0*R_temp[Cart::xxyyy][Cart::yz][0]+wmp0*R_temp[Cart::xxyyy][Cart::yz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::yz][0]-gfak*R_temp[Cart::xyyy][Cart::yz][1])
R_temp[Cart::xxxyyy][Cart::xx][0]+=pma0*R_temp[Cart::xxyyy][Cart::xx][0]+wmp0*R_temp[Cart::xxyyy][Cart::xx][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xx][0]-gfak*R_temp[Cart::xyyy][Cart::xx][1])+0.5/_decay*2*R_temp[Cart::xxyyy][Cart::x][1]
R_temp[Cart::xxxyyy][Cart::xz][0]+=pma0*R_temp[Cart::xxyyy][Cart::xz][0]+wmp0*R_temp[Cart::xxyyy][Cart::xz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xz][0]-gfak*R_temp[Cart::xyyy][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::z][1]
R_temp[Cart::xxxyyy][Cart::zz][0]+=pma0*R_temp[Cart::xxyyy][Cart::zz][0]+wmp0*R_temp[Cart::xxyyy][Cart::zz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::zz][0]-gfak*R_temp[Cart::xyyy][Cart::zz][1])
R_temp[Cart::xxyyyz][Cart::yy][0]+=pma2*R_temp[Cart::xxyyy][Cart::yy][0]+wmp2*R_temp[Cart::xxyyy][Cart::yy][1]
R_temp[Cart::xxyyyz][Cart::xy][0]+=pma2*R_temp[Cart::xxyyy][Cart::xy][0]+wmp2*R_temp[Cart::xxyyy][Cart::xy][1]
R_temp[Cart::xxyyyz][Cart::yz][0]+=pma2*R_temp[Cart::xxyyy][Cart::yz][0]+wmp2*R_temp[Cart::xxyyy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::y][1]
R_temp[Cart::xxyyyz][Cart::xx][0]+=pma2*R_temp[Cart::xxyyy][Cart::xx][0]+wmp2*R_temp[Cart::xxyyy][Cart::xx][1]
R_temp[Cart::xxyyyz][Cart::xz][0]+=pma2*R_temp[Cart::xxyyy][Cart::xz][0]+wmp2*R_temp[Cart::xxyyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::x][1]
R_temp[Cart::xxyyyz][Cart::zz][0]+=pma2*R_temp[Cart::xxyyy][Cart::zz][0]+wmp2*R_temp[Cart::xxyyy][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::xxyyy][Cart::z][1]
R_temp[Cart::xyyyzz][Cart::yy][0]+=pma0*R_temp[Cart::yyyzz][Cart::yy][0]+wmp0*R_temp[Cart::yyyzz][Cart::yy][1]
R_temp[Cart::xyyyzz][Cart::xy][0]+=pma0*R_temp[Cart::yyyzz][Cart::xy][0]+wmp0*R_temp[Cart::yyyzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yyyzz][Cart::y][1]
R_temp[Cart::xyyyzz][Cart::yz][0]+=pma0*R_temp[Cart::yyyzz][Cart::yz][0]+wmp0*R_temp[Cart::yyyzz][Cart::yz][1]
R_temp[Cart::xyyyzz][Cart::xx][0]+=pma0*R_temp[Cart::yyyzz][Cart::xx][0]+wmp0*R_temp[Cart::yyyzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yyyzz][Cart::x][1]
R_temp[Cart::xyyyzz][Cart::xz][0]+=pma0*R_temp[Cart::yyyzz][Cart::xz][0]+wmp0*R_temp[Cart::yyyzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyyzz][Cart::z][1]
R_temp[Cart::xyyyzz][Cart::zz][0]+=pma0*R_temp[Cart::yyyzz][Cart::zz][0]+wmp0*R_temp[Cart::yyyzz][Cart::zz][1]
R_temp[Cart::yyyzzz][Cart::yy][0]+=pma1*R_temp[Cart::yyzzz][Cart::yy][0]+wmp1*R_temp[Cart::yyzzz][Cart::yy][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::yy][0]-gfak*R_temp[Cart::yzzz][Cart::yy][1])+0.5/_decay*2*R_temp[Cart::yyzzz][Cart::y][1]
R_temp[Cart::yyyzzz][Cart::xy][0]+=pma1*R_temp[Cart::yyzzz][Cart::xy][0]+wmp1*R_temp[Cart::yyzzz][Cart::xy][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xy][0]-gfak*R_temp[Cart::yzzz][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::x][1]
R_temp[Cart::yyyzzz][Cart::yz][0]+=pma1*R_temp[Cart::yyzzz][Cart::yz][0]+wmp1*R_temp[Cart::yyzzz][Cart::yz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::yz][0]-gfak*R_temp[Cart::yzzz][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::z][1]
R_temp[Cart::yyyzzz][Cart::xx][0]+=pma1*R_temp[Cart::yyzzz][Cart::xx][0]+wmp1*R_temp[Cart::yyzzz][Cart::xx][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xx][0]-gfak*R_temp[Cart::yzzz][Cart::xx][1])
R_temp[Cart::yyyzzz][Cart::xz][0]+=pma1*R_temp[Cart::yyzzz][Cart::xz][0]+wmp1*R_temp[Cart::yyzzz][Cart::xz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xz][0]-gfak*R_temp[Cart::yzzz][Cart::xz][1])
R_temp[Cart::yyyzzz][Cart::zz][0]+=pma1*R_temp[Cart::yyzzz][Cart::zz][0]+wmp1*R_temp[Cart::yyzzz][Cart::zz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::zz][0]-gfak*R_temp[Cart::yzzz][Cart::zz][1])
R_temp[Cart::xxxxyy][Cart::yy][0]+=pma1*R_temp[Cart::xxxxy][Cart::yy][0]+wmp1*R_temp[Cart::xxxxy][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxxxy][Cart::y][1]
R_temp[Cart::xxxxyy][Cart::xy][0]+=pma1*R_temp[Cart::xxxxy][Cart::xy][0]+wmp1*R_temp[Cart::xxxxy][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxxxy][Cart::x][1]
R_temp[Cart::xxxxyy][Cart::yz][0]+=pma1*R_temp[Cart::xxxxy][Cart::yz][0]+wmp1*R_temp[Cart::xxxxy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxxy][Cart::z][1]
R_temp[Cart::xxxxyy][Cart::xx][0]+=pma1*R_temp[Cart::xxxxy][Cart::xx][0]+wmp1*R_temp[Cart::xxxxy][Cart::xx][1]
R_temp[Cart::xxxxyy][Cart::xz][0]+=pma1*R_temp[Cart::xxxxy][Cart::xz][0]+wmp1*R_temp[Cart::xxxxy][Cart::xz][1]
R_temp[Cart::xxxxyy][Cart::zz][0]+=pma1*R_temp[Cart::xxxxy][Cart::zz][0]+wmp1*R_temp[Cart::xxxxy][Cart::zz][1]
R_temp[Cart::xxxyyz][Cart::yy][0]+=pma2*R_temp[Cart::xxxyy][Cart::yy][0]+wmp2*R_temp[Cart::xxxyy][Cart::yy][1]
R_temp[Cart::xxxyyz][Cart::xy][0]+=pma2*R_temp[Cart::xxxyy][Cart::xy][0]+wmp2*R_temp[Cart::xxxyy][Cart::xy][1]
R_temp[Cart::xxxyyz][Cart::yz][0]+=pma2*R_temp[Cart::xxxyy][Cart::yz][0]+wmp2*R_temp[Cart::xxxyy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxyy][Cart::y][1]
R_temp[Cart::xxxyyz][Cart::xx][0]+=pma2*R_temp[Cart::xxxyy][Cart::xx][0]+wmp2*R_temp[Cart::xxxyy][Cart::xx][1]
R_temp[Cart::xxxyyz][Cart::xz][0]+=pma2*R_temp[Cart::xxxyy][Cart::xz][0]+wmp2*R_temp[Cart::xxxyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xxxyy][Cart::x][1]
R_temp[Cart::xxxyyz][Cart::zz][0]+=pma2*R_temp[Cart::xxxyy][Cart::zz][0]+wmp2*R_temp[Cart::xxxyy][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::xxxyy][Cart::z][1]
R_temp[Cart::xxyyzz][Cart::yy][0]+=pma0*R_temp[Cart::xyyzz][Cart::yy][0]+wmp0*R_temp[Cart::xyyzz][Cart::yy][1]
R_temp[Cart::xxyyzz][Cart::xy][0]+=pma0*R_temp[Cart::xyyzz][Cart::xy][0]+wmp0*R_temp[Cart::xyyzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xyyzz][Cart::y][1]
R_temp[Cart::xxyyzz][Cart::yz][0]+=pma0*R_temp[Cart::xyyzz][Cart::yz][0]+wmp0*R_temp[Cart::xyyzz][Cart::yz][1]
R_temp[Cart::xxyyzz][Cart::xx][0]+=pma0*R_temp[Cart::xyyzz][Cart::xx][0]+wmp0*R_temp[Cart::xyyzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::xyyzz][Cart::x][1]
R_temp[Cart::xxyyzz][Cart::xz][0]+=pma0*R_temp[Cart::xyyzz][Cart::xz][0]+wmp0*R_temp[Cart::xyyzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xyyzz][Cart::z][1]
R_temp[Cart::xxyyzz][Cart::zz][0]+=pma0*R_temp[Cart::xyyzz][Cart::zz][0]+wmp0*R_temp[Cart::xyyzz][Cart::zz][1]
R_temp[Cart::xyyzzz][Cart::yy][0]+=pma0*R_temp[Cart::yyzzz][Cart::yy][0]+wmp0*R_temp[Cart::yyzzz][Cart::yy][1]
R_temp[Cart::xyyzzz][Cart::xy][0]+=pma0*R_temp[Cart::yyzzz][Cart::xy][0]+wmp0*R_temp[Cart::yyzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::y][1]
R_temp[Cart::xyyzzz][Cart::yz][0]+=pma0*R_temp[Cart::yyzzz][Cart::yz][0]+wmp0*R_temp[Cart::yyzzz][Cart::yz][1]
R_temp[Cart::xyyzzz][Cart::xx][0]+=pma0*R_temp[Cart::yyzzz][Cart::xx][0]+wmp0*R_temp[Cart::yyzzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yyzzz][Cart::x][1]
R_temp[Cart::xyyzzz][Cart::xz][0]+=pma0*R_temp[Cart::yyzzz][Cart::xz][0]+wmp0*R_temp[Cart::yyzzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::z][1]
R_temp[Cart::xyyzzz][Cart::zz][0]+=pma0*R_temp[Cart::yyzzz][Cart::zz][0]+wmp0*R_temp[Cart::yyzzz][Cart::zz][1]
R_temp[Cart::yyzzzz][Cart::yy][0]+=pma1*R_temp[Cart::yzzzz][Cart::yy][0]+wmp1*R_temp[Cart::yzzzz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::yzzzz][Cart::y][1]
R_temp[Cart::yyzzzz][Cart::xy][0]+=pma1*R_temp[Cart::yzzzz][Cart::xy][0]+wmp1*R_temp[Cart::yzzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::x][1]
R_temp[Cart::yyzzzz][Cart::yz][0]+=pma1*R_temp[Cart::yzzzz][Cart::yz][0]+wmp1*R_temp[Cart::yzzzz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::z][1]
R_temp[Cart::yyzzzz][Cart::xx][0]+=pma1*R_temp[Cart::yzzzz][Cart::xx][0]+wmp1*R_temp[Cart::yzzzz][Cart::xx][1]
R_temp[Cart::yyzzzz][Cart::xz][0]+=pma1*R_temp[Cart::yzzzz][Cart::xz][0]+wmp1*R_temp[Cart::yzzzz][Cart::xz][1]
R_temp[Cart::yyzzzz][Cart::zz][0]+=pma1*R_temp[Cart::yzzzz][Cart::zz][0]+wmp1*R_temp[Cart::yzzzz][Cart::zz][1]
R_temp[Cart::xxxxxy][Cart::yy][0]+=pma1*R_temp[Cart::xxxxx][Cart::yy][0]+wmp1*R_temp[Cart::xxxxx][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::y][1]
R_temp[Cart::xxxxxy][Cart::xy][0]+=pma1*R_temp[Cart::xxxxx][Cart::xy][0]+wmp1*R_temp[Cart::xxxxx][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::x][1]
R_temp[Cart::xxxxxy][Cart::yz][0]+=pma1*R_temp[Cart::xxxxx][Cart::yz][0]+wmp1*R_temp[Cart::xxxxx][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::z][1]
R_temp[Cart::xxxxxy][Cart::xx][0]+=pma1*R_temp[Cart::xxxxx][Cart::xx][0]+wmp1*R_temp[Cart::xxxxx][Cart::xx][1]
R_temp[Cart::xxxxxy][Cart::xz][0]+=pma1*R_temp[Cart::xxxxx][Cart::xz][0]+wmp1*R_temp[Cart::xxxxx][Cart::xz][1]
R_temp[Cart::xxxxxy][Cart::zz][0]+=pma1*R_temp[Cart::xxxxx][Cart::zz][0]+wmp1*R_temp[Cart::xxxxx][Cart::zz][1]
R_temp[Cart::xxxxyz][Cart::yy][0]+=pma1*R_temp[Cart::xxxxz][Cart::yy][0]+wmp1*R_temp[Cart::xxxxz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxxxz][Cart::y][1]
R_temp[Cart::xxxxyz][Cart::xy][0]+=pma1*R_temp[Cart::xxxxz][Cart::xy][0]+wmp1*R_temp[Cart::xxxxz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::x][1]
R_temp[Cart::xxxxyz][Cart::yz][0]+=pma1*R_temp[Cart::xxxxz][Cart::yz][0]+wmp1*R_temp[Cart::xxxxz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::z][1]
R_temp[Cart::xxxxyz][Cart::xx][0]+=pma1*R_temp[Cart::xxxxz][Cart::xx][0]+wmp1*R_temp[Cart::xxxxz][Cart::xx][1]
R_temp[Cart::xxxxyz][Cart::xz][0]+=pma1*R_temp[Cart::xxxxz][Cart::xz][0]+wmp1*R_temp[Cart::xxxxz][Cart::xz][1]
R_temp[Cart::xxxxyz][Cart::zz][0]+=pma1*R_temp[Cart::xxxxz][Cart::zz][0]+wmp1*R_temp[Cart::xxxxz][Cart::zz][1]
R_temp[Cart::xxxyzz][Cart::yy][0]+=pma1*R_temp[Cart::xxxzz][Cart::yy][0]+wmp1*R_temp[Cart::xxxzz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxxzz][Cart::y][1]
R_temp[Cart::xxxyzz][Cart::xy][0]+=pma1*R_temp[Cart::xxxzz][Cart::xy][0]+wmp1*R_temp[Cart::xxxzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxxzz][Cart::x][1]
R_temp[Cart::xxxyzz][Cart::yz][0]+=pma1*R_temp[Cart::xxxzz][Cart::yz][0]+wmp1*R_temp[Cart::xxxzz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxzz][Cart::z][1]
R_temp[Cart::xxxyzz][Cart::xx][0]+=pma1*R_temp[Cart::xxxzz][Cart::xx][0]+wmp1*R_temp[Cart::xxxzz][Cart::xx][1]
R_temp[Cart::xxxyzz][Cart::xz][0]+=pma1*R_temp[Cart::xxxzz][Cart::xz][0]+wmp1*R_temp[Cart::xxxzz][Cart::xz][1]
R_temp[Cart::xxxyzz][Cart::zz][0]+=pma1*R_temp[Cart::xxxzz][Cart::zz][0]+wmp1*R_temp[Cart::xxxzz][Cart::zz][1]
R_temp[Cart::xxyzzz][Cart::yy][0]+=pma1*R_temp[Cart::xxzzz][Cart::yy][0]+wmp1*R_temp[Cart::xxzzz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxzzz][Cart::y][1]
R_temp[Cart::xxyzzz][Cart::xy][0]+=pma1*R_temp[Cart::xxzzz][Cart::xy][0]+wmp1*R_temp[Cart::xxzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::x][1]
R_temp[Cart::xxyzzz][Cart::yz][0]+=pma1*R_temp[Cart::xxzzz][Cart::yz][0]+wmp1*R_temp[Cart::xxzzz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::z][1]
R_temp[Cart::xxyzzz][Cart::xx][0]+=pma1*R_temp[Cart::xxzzz][Cart::xx][0]+wmp1*R_temp[Cart::xxzzz][Cart::xx][1]
R_temp[Cart::xxyzzz][Cart::xz][0]+=pma1*R_temp[Cart::xxzzz][Cart::xz][0]+wmp1*R_temp[Cart::xxzzz][Cart::xz][1]
R_temp[Cart::xxyzzz][Cart::zz][0]+=pma1*R_temp[Cart::xxzzz][Cart::zz][0]+wmp1*R_temp[Cart::xxzzz][Cart::zz][1]
R_temp[Cart::xyzzzz][Cart::yy][0]+=pma0*R_temp[Cart::yzzzz][Cart::yy][0]+wmp0*R_temp[Cart::yzzzz][Cart::yy][1]
R_temp[Cart::xyzzzz][Cart::xy][0]+=pma0*R_temp[Cart::yzzzz][Cart::xy][0]+wmp0*R_temp[Cart::yzzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::y][1]
R_temp[Cart::xyzzzz][Cart::yz][0]+=pma0*R_temp[Cart::yzzzz][Cart::yz][0]+wmp0*R_temp[Cart::yzzzz][Cart::yz][1]
R_temp[Cart::xyzzzz][Cart::xx][0]+=pma0*R_temp[Cart::yzzzz][Cart::xx][0]+wmp0*R_temp[Cart::yzzzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yzzzz][Cart::x][1]
R_temp[Cart::xyzzzz][Cart::xz][0]+=pma0*R_temp[Cart::yzzzz][Cart::xz][0]+wmp0*R_temp[Cart::yzzzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::z][1]
R_temp[Cart::xyzzzz][Cart::zz][0]+=pma0*R_temp[Cart::yzzzz][Cart::zz][0]+wmp0*R_temp[Cart::yzzzz][Cart::zz][1]
R_temp[Cart::yzzzzz][Cart::yy][0]+=pma1*R_temp[Cart::zzzzz][Cart::yy][0]+wmp1*R_temp[Cart::zzzzz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::y][1]
R_temp[Cart::yzzzzz][Cart::xy][0]+=pma1*R_temp[Cart::zzzzz][Cart::xy][0]+wmp1*R_temp[Cart::zzzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::x][1]
R_temp[Cart::yzzzzz][Cart::yz][0]+=pma1*R_temp[Cart::zzzzz][Cart::yz][0]+wmp1*R_temp[Cart::zzzzz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::z][1]
R_temp[Cart::yzzzzz][Cart::xx][0]+=pma1*R_temp[Cart::zzzzz][Cart::xx][0]+wmp1*R_temp[Cart::zzzzz][Cart::xx][1]
R_temp[Cart::yzzzzz][Cart::xz][0]+=pma1*R_temp[Cart::zzzzz][Cart::xz][0]+wmp1*R_temp[Cart::zzzzz][Cart::xz][1]
R_temp[Cart::yzzzzz][Cart::zz][0]+=pma1*R_temp[Cart::zzzzz][Cart::zz][0]+wmp1*R_temp[Cart::zzzzz][Cart::zz][1]
R_temp[Cart::xxxxxx][Cart::yy][0]+=pma0*R_temp[Cart::xxxxx][Cart::yy][0]+wmp0*R_temp[Cart::xxxxx][Cart::yy][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::yy][0]-gfak*R_temp[Cart::xxxx][Cart::yy][1])
R_temp[Cart::xxxxxx][Cart::xy][0]+=pma0*R_temp[Cart::xxxxx][Cart::xy][0]+wmp0*R_temp[Cart::xxxxx][Cart::xy][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xy][0]-gfak*R_temp[Cart::xxxx][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::y][1]
R_temp[Cart::xxxxxx][Cart::yz][0]+=pma0*R_temp[Cart::xxxxx][Cart::yz][0]+wmp0*R_temp[Cart::xxxxx][Cart::yz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::yz][0]-gfak*R_temp[Cart::xxxx][Cart::yz][1])
R_temp[Cart::xxxxxx][Cart::xx][0]+=pma0*R_temp[Cart::xxxxx][Cart::xx][0]+wmp0*R_temp[Cart::xxxxx][Cart::xx][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xx][0]-gfak*R_temp[Cart::xxxx][Cart::xx][1])+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::x][1]
R_temp[Cart::xxxxxx][Cart::xz][0]+=pma0*R_temp[Cart::xxxxx][Cart::xz][0]+wmp0*R_temp[Cart::xxxxx][Cart::xz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xz][0]-gfak*R_temp[Cart::xxxx][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::z][1]
R_temp[Cart::xxxxxx][Cart::zz][0]+=pma0*R_temp[Cart::xxxxx][Cart::zz][0]+wmp0*R_temp[Cart::xxxxx][Cart::zz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::zz][0]-gfak*R_temp[Cart::xxxx][Cart::zz][1])
R_temp[Cart::xxxxxz][Cart::yy][0]+=pma2*R_temp[Cart::xxxxx][Cart::yy][0]+wmp2*R_temp[Cart::xxxxx][Cart::yy][1]
R_temp[Cart::xxxxxz][Cart::xy][0]+=pma2*R_temp[Cart::xxxxx][Cart::xy][0]+wmp2*R_temp[Cart::xxxxx][Cart::xy][1]
R_temp[Cart::xxxxxz][Cart::yz][0]+=pma2*R_temp[Cart::xxxxx][Cart::yz][0]+wmp2*R_temp[Cart::xxxxx][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::y][1]
R_temp[Cart::xxxxxz][Cart::xx][0]+=pma2*R_temp[Cart::xxxxx][Cart::xx][0]+wmp2*R_temp[Cart::xxxxx][Cart::xx][1]
R_temp[Cart::xxxxxz][Cart::xz][0]+=pma2*R_temp[Cart::xxxxx][Cart::xz][0]+wmp2*R_temp[Cart::xxxxx][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::x][1]
R_temp[Cart::xxxxxz][Cart::zz][0]+=pma2*R_temp[Cart::xxxxx][Cart::zz][0]+wmp2*R_temp[Cart::xxxxx][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::z][1]
R_temp[Cart::xxxxzz][Cart::yy][0]+=pma2*R_temp[Cart::xxxxz][Cart::yy][0]+wmp2*R_temp[Cart::xxxxz][Cart::yy][1]
R_temp[Cart::xxxxzz][Cart::xy][0]+=pma2*R_temp[Cart::xxxxz][Cart::xy][0]+wmp2*R_temp[Cart::xxxxz][Cart::xy][1]
R_temp[Cart::xxxxzz][Cart::yz][0]+=pma2*R_temp[Cart::xxxxz][Cart::yz][0]+wmp2*R_temp[Cart::xxxxz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::y][1]
R_temp[Cart::xxxxzz][Cart::xx][0]+=pma2*R_temp[Cart::xxxxz][Cart::xx][0]+wmp2*R_temp[Cart::xxxxz][Cart::xx][1]
R_temp[Cart::xxxxzz][Cart::xz][0]+=pma2*R_temp[Cart::xxxxz][Cart::xz][0]+wmp2*R_temp[Cart::xxxxz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::x][1]
R_temp[Cart::xxxxzz][Cart::zz][0]+=pma2*R_temp[Cart::xxxxz][Cart::zz][0]+wmp2*R_temp[Cart::xxxxz][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::xxxxz][Cart::z][1]
R_temp[Cart::xxxzzz][Cart::yy][0]+=pma0*R_temp[Cart::xxzzz][Cart::yy][0]+wmp0*R_temp[Cart::xxzzz][Cart::yy][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::yy][0]-gfak*R_temp[Cart::xzzz][Cart::yy][1])
R_temp[Cart::xxxzzz][Cart::xy][0]+=pma0*R_temp[Cart::xxzzz][Cart::xy][0]+wmp0*R_temp[Cart::xxzzz][Cart::xy][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xy][0]-gfak*R_temp[Cart::xzzz][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::y][1]
R_temp[Cart::xxxzzz][Cart::yz][0]+=pma0*R_temp[Cart::xxzzz][Cart::yz][0]+wmp0*R_temp[Cart::xxzzz][Cart::yz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::yz][0]-gfak*R_temp[Cart::xzzz][Cart::yz][1])
R_temp[Cart::xxxzzz][Cart::xx][0]+=pma0*R_temp[Cart::xxzzz][Cart::xx][0]+wmp0*R_temp[Cart::xxzzz][Cart::xx][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xx][0]-gfak*R_temp[Cart::xzzz][Cart::xx][1])+0.5/_decay*2*R_temp[Cart::xxzzz][Cart::x][1]
R_temp[Cart::xxxzzz][Cart::xz][0]+=pma0*R_temp[Cart::xxzzz][Cart::xz][0]+wmp0*R_temp[Cart::xxzzz][Cart::xz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xz][0]-gfak*R_temp[Cart::xzzz][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::z][1]
R_temp[Cart::xxxzzz][Cart::zz][0]+=pma0*R_temp[Cart::xxzzz][Cart::zz][0]+wmp0*R_temp[Cart::xxzzz][Cart::zz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::zz][0]-gfak*R_temp[Cart::xzzz][Cart::zz][1])
R_temp[Cart::xxzzzz][Cart::yy][0]+=pma0*R_temp[Cart::xzzzz][Cart::yy][0]+wmp0*R_temp[Cart::xzzzz][Cart::yy][1]
R_temp[Cart::xxzzzz][Cart::xy][0]+=pma0*R_temp[Cart::xzzzz][Cart::xy][0]+wmp0*R_temp[Cart::xzzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xzzzz][Cart::y][1]
R_temp[Cart::xxzzzz][Cart::yz][0]+=pma0*R_temp[Cart::xzzzz][Cart::yz][0]+wmp0*R_temp[Cart::xzzzz][Cart::yz][1]
R_temp[Cart::xxzzzz][Cart::xx][0]+=pma0*R_temp[Cart::xzzzz][Cart::xx][0]+wmp0*R_temp[Cart::xzzzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::xzzzz][Cart::x][1]
R_temp[Cart::xxzzzz][Cart::xz][0]+=pma0*R_temp[Cart::xzzzz][Cart::xz][0]+wmp0*R_temp[Cart::xzzzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xzzzz][Cart::z][1]
R_temp[Cart::xxzzzz][Cart::zz][0]+=pma0*R_temp[Cart::xzzzz][Cart::zz][0]+wmp0*R_temp[Cart::xzzzz][Cart::zz][1]
R_temp[Cart::xzzzzz][Cart::yy][0]+=pma0*R_temp[Cart::zzzzz][Cart::yy][0]+wmp0*R_temp[Cart::zzzzz][Cart::yy][1]
R_temp[Cart::xzzzzz][Cart::xy][0]+=pma0*R_temp[Cart::zzzzz][Cart::xy][0]+wmp0*R_temp[Cart::zzzzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::y][1]
R_temp[Cart::xzzzzz][Cart::yz][0]+=pma0*R_temp[Cart::zzzzz][Cart::yz][0]+wmp0*R_temp[Cart::zzzzz][Cart::yz][1]
R_temp[Cart::xzzzzz][Cart::xx][0]+=pma0*R_temp[Cart::zzzzz][Cart::xx][0]+wmp0*R_temp[Cart::zzzzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::x][1]
R_temp[Cart::xzzzzz][Cart::xz][0]+=pma0*R_temp[Cart::zzzzz][Cart::xz][0]+wmp0*R_temp[Cart::zzzzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::z][1]
R_temp[Cart::xzzzzz][Cart::zz][0]+=pma0*R_temp[Cart::zzzzz][Cart::zz][0]+wmp0*R_temp[Cart::zzzzz][Cart::zz][1]
R_temp[Cart::zzzzzz][Cart::yy][0]+=pma2*R_temp[Cart::zzzzz][Cart::yy][0]+wmp2*R_temp[Cart::zzzzz][Cart::yy][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::yy][0]-gfak*R_temp[Cart::zzzz][Cart::yy][1])
R_temp[Cart::zzzzzz][Cart::xy][0]+=pma2*R_temp[Cart::zzzzz][Cart::xy][0]+wmp2*R_temp[Cart::zzzzz][Cart::xy][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xy][0]-gfak*R_temp[Cart::zzzz][Cart::xy][1])
R_temp[Cart::zzzzzz][Cart::yz][0]+=pma2*R_temp[Cart::zzzzz][Cart::yz][0]+wmp2*R_temp[Cart::zzzzz][Cart::yz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::yz][0]-gfak*R_temp[Cart::zzzz][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::y][1]
R_temp[Cart::zzzzzz][Cart::xx][0]+=pma2*R_temp[Cart::zzzzz][Cart::xx][0]+wmp2*R_temp[Cart::zzzzz][Cart::xx][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xx][0]-gfak*R_temp[Cart::zzzz][Cart::xx][1])
R_temp[Cart::zzzzzz][Cart::xz][0]+=pma2*R_temp[Cart::zzzzz][Cart::xz][0]+wmp2*R_temp[Cart::zzzzz][Cart::xz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xz][0]-gfak*R_temp[Cart::zzzz][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::x][1]
R_temp[Cart::zzzzzz][Cart::zz][0]+=pma2*R_temp[Cart::zzzzz][Cart::zz][0]+wmp2*R_temp[Cart::zzzzz][Cart::zz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::zz][0]-gfak*R_temp[Cart::zzzz][Cart::zz][1])+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::z][1]
}}
//------------------------------------------------------

//Integral i - s - f - m0
if (_mmax >9 ){
if (_lmax_alpha>5 && _lmax_gamma>2){

R_temp[Cart::yyyyyy][Cart::yyy][0]+=pma1*R_temp[Cart::yyyyy][Cart::yyy][0]+wmp1*R_temp[Cart::yyyyy][Cart::yyy][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::yyy][0]-gfak*R_temp[Cart::yyyy][Cart::yyy][1])+0.5/_decay*3*R_temp[Cart::yyyyy][Cart::yy][1]
R_temp[Cart::yyyyyy][Cart::xyy][0]+=pma1*R_temp[Cart::yyyyy][Cart::xyy][0]+wmp1*R_temp[Cart::yyyyy][Cart::xyy][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xyy][0]-gfak*R_temp[Cart::yyyy][Cart::xyy][1])+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::xy][1]
R_temp[Cart::yyyyyy][Cart::yyz][0]+=pma1*R_temp[Cart::yyyyy][Cart::yyz][0]+wmp1*R_temp[Cart::yyyyy][Cart::yyz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::yyz][0]-gfak*R_temp[Cart::yyyy][Cart::yyz][1])+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::yz][1]
R_temp[Cart::yyyyyy][Cart::xxy][0]+=pma1*R_temp[Cart::yyyyy][Cart::xxy][0]+wmp1*R_temp[Cart::yyyyy][Cart::xxy][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xxy][0]-gfak*R_temp[Cart::yyyy][Cart::xxy][1])+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::xx][1]
R_temp[Cart::yyyyyy][Cart::xyz][0]+=pma1*R_temp[Cart::yyyyy][Cart::xyz][0]+wmp1*R_temp[Cart::yyyyy][Cart::xyz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xyz][0]-gfak*R_temp[Cart::yyyy][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::xz][1]
R_temp[Cart::yyyyyy][Cart::yzz][0]+=pma1*R_temp[Cart::yyyyy][Cart::yzz][0]+wmp1*R_temp[Cart::yyyyy][Cart::yzz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::yzz][0]-gfak*R_temp[Cart::yyyy][Cart::yzz][1])+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::zz][1]
R_temp[Cart::yyyyyy][Cart::xxx][0]+=pma1*R_temp[Cart::yyyyy][Cart::xxx][0]+wmp1*R_temp[Cart::yyyyy][Cart::xxx][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xxx][0]-gfak*R_temp[Cart::yyyy][Cart::xxx][1])
R_temp[Cart::yyyyyy][Cart::xxz][0]+=pma1*R_temp[Cart::yyyyy][Cart::xxz][0]+wmp1*R_temp[Cart::yyyyy][Cart::xxz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xxz][0]-gfak*R_temp[Cart::yyyy][Cart::xxz][1])
R_temp[Cart::yyyyyy][Cart::xzz][0]+=pma1*R_temp[Cart::yyyyy][Cart::xzz][0]+wmp1*R_temp[Cart::yyyyy][Cart::xzz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::xzz][0]-gfak*R_temp[Cart::yyyy][Cart::xzz][1])
R_temp[Cart::yyyyyy][Cart::zzz][0]+=pma1*R_temp[Cart::yyyyy][Cart::zzz][0]+wmp1*R_temp[Cart::yyyyy][Cart::zzz][1]+4*rzeta*(R_temp[Cart::yyyy][Cart::zzz][0]-gfak*R_temp[Cart::yyyy][Cart::zzz][1])
R_temp[Cart::xyyyyy][Cart::yyy][0]+=pma0*R_temp[Cart::yyyyy][Cart::yyy][0]+wmp0*R_temp[Cart::yyyyy][Cart::yyy][1]
R_temp[Cart::xyyyyy][Cart::xyy][0]+=pma0*R_temp[Cart::yyyyy][Cart::xyy][0]+wmp0*R_temp[Cart::yyyyy][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::yy][1]
R_temp[Cart::xyyyyy][Cart::yyz][0]+=pma0*R_temp[Cart::yyyyy][Cart::yyz][0]+wmp0*R_temp[Cart::yyyyy][Cart::yyz][1]
R_temp[Cart::xyyyyy][Cart::xxy][0]+=pma0*R_temp[Cart::yyyyy][Cart::xxy][0]+wmp0*R_temp[Cart::yyyyy][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::xy][1]
R_temp[Cart::xyyyyy][Cart::xyz][0]+=pma0*R_temp[Cart::yyyyy][Cart::xyz][0]+wmp0*R_temp[Cart::yyyyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::yz][1]
R_temp[Cart::xyyyyy][Cart::yzz][0]+=pma0*R_temp[Cart::yyyyy][Cart::yzz][0]+wmp0*R_temp[Cart::yyyyy][Cart::yzz][1]
R_temp[Cart::xyyyyy][Cart::xxx][0]+=pma0*R_temp[Cart::yyyyy][Cart::xxx][0]+wmp0*R_temp[Cart::yyyyy][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yyyyy][Cart::xx][1]
R_temp[Cart::xyyyyy][Cart::xxz][0]+=pma0*R_temp[Cart::yyyyy][Cart::xxz][0]+wmp0*R_temp[Cart::yyyyy][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::xz][1]
R_temp[Cart::xyyyyy][Cart::xzz][0]+=pma0*R_temp[Cart::yyyyy][Cart::xzz][0]+wmp0*R_temp[Cart::yyyyy][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::zz][1]
R_temp[Cart::xyyyyy][Cart::zzz][0]+=pma0*R_temp[Cart::yyyyy][Cart::zzz][0]+wmp0*R_temp[Cart::yyyyy][Cart::zzz][1]
R_temp[Cart::yyyyyz][Cart::yyy][0]+=pma2*R_temp[Cart::yyyyy][Cart::yyy][0]+wmp2*R_temp[Cart::yyyyy][Cart::yyy][1]
R_temp[Cart::yyyyyz][Cart::xyy][0]+=pma2*R_temp[Cart::yyyyy][Cart::xyy][0]+wmp2*R_temp[Cart::yyyyy][Cart::xyy][1]
R_temp[Cart::yyyyyz][Cart::yyz][0]+=pma2*R_temp[Cart::yyyyy][Cart::yyz][0]+wmp2*R_temp[Cart::yyyyy][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::yy][1]
R_temp[Cart::yyyyyz][Cart::xxy][0]+=pma2*R_temp[Cart::yyyyy][Cart::xxy][0]+wmp2*R_temp[Cart::yyyyy][Cart::xxy][1]
R_temp[Cart::yyyyyz][Cart::xyz][0]+=pma2*R_temp[Cart::yyyyy][Cart::xyz][0]+wmp2*R_temp[Cart::yyyyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::xy][1]
R_temp[Cart::yyyyyz][Cart::yzz][0]+=pma2*R_temp[Cart::yyyyy][Cart::yzz][0]+wmp2*R_temp[Cart::yyyyy][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::yz][1]
R_temp[Cart::yyyyyz][Cart::xxx][0]+=pma2*R_temp[Cart::yyyyy][Cart::xxx][0]+wmp2*R_temp[Cart::yyyyy][Cart::xxx][1]
R_temp[Cart::yyyyyz][Cart::xxz][0]+=pma2*R_temp[Cart::yyyyy][Cart::xxz][0]+wmp2*R_temp[Cart::yyyyy][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::yyyyy][Cart::xx][1]
R_temp[Cart::yyyyyz][Cart::xzz][0]+=pma2*R_temp[Cart::yyyyy][Cart::xzz][0]+wmp2*R_temp[Cart::yyyyy][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::yyyyy][Cart::xz][1]
R_temp[Cart::yyyyyz][Cart::zzz][0]+=pma2*R_temp[Cart::yyyyy][Cart::zzz][0]+wmp2*R_temp[Cart::yyyyy][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::yyyyy][Cart::zz][1]
R_temp[Cart::xxyyyy][Cart::yyy][0]+=pma0*R_temp[Cart::xyyyy][Cart::yyy][0]+wmp0*R_temp[Cart::xyyyy][Cart::yyy][1]
R_temp[Cart::xxyyyy][Cart::xyy][0]+=pma0*R_temp[Cart::xyyyy][Cart::xyy][0]+wmp0*R_temp[Cart::xyyyy][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::xyyyy][Cart::yy][1]
R_temp[Cart::xxyyyy][Cart::yyz][0]+=pma0*R_temp[Cart::xyyyy][Cart::yyz][0]+wmp0*R_temp[Cart::xyyyy][Cart::yyz][1]
R_temp[Cart::xxyyyy][Cart::xxy][0]+=pma0*R_temp[Cart::xyyyy][Cart::xxy][0]+wmp0*R_temp[Cart::xyyyy][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::xyyyy][Cart::xy][1]
R_temp[Cart::xxyyyy][Cart::xyz][0]+=pma0*R_temp[Cart::xyyyy][Cart::xyz][0]+wmp0*R_temp[Cart::xyyyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xyyyy][Cart::yz][1]
R_temp[Cart::xxyyyy][Cart::yzz][0]+=pma0*R_temp[Cart::xyyyy][Cart::yzz][0]+wmp0*R_temp[Cart::xyyyy][Cart::yzz][1]
R_temp[Cart::xxyyyy][Cart::xxx][0]+=pma0*R_temp[Cart::xyyyy][Cart::xxx][0]+wmp0*R_temp[Cart::xyyyy][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::xyyyy][Cart::xx][1]
R_temp[Cart::xxyyyy][Cart::xxz][0]+=pma0*R_temp[Cart::xyyyy][Cart::xxz][0]+wmp0*R_temp[Cart::xyyyy][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::xyyyy][Cart::xz][1]
R_temp[Cart::xxyyyy][Cart::xzz][0]+=pma0*R_temp[Cart::xyyyy][Cart::xzz][0]+wmp0*R_temp[Cart::xyyyy][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::xyyyy][Cart::zz][1]
R_temp[Cart::xxyyyy][Cart::zzz][0]+=pma0*R_temp[Cart::xyyyy][Cart::zzz][0]+wmp0*R_temp[Cart::xyyyy][Cart::zzz][1]
R_temp[Cart::xyyyyz][Cart::yyy][0]+=pma0*R_temp[Cart::yyyyz][Cart::yyy][0]+wmp0*R_temp[Cart::yyyyz][Cart::yyy][1]
R_temp[Cart::xyyyyz][Cart::xyy][0]+=pma0*R_temp[Cart::yyyyz][Cart::xyy][0]+wmp0*R_temp[Cart::yyyyz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::yy][1]
R_temp[Cart::xyyyyz][Cart::yyz][0]+=pma0*R_temp[Cart::yyyyz][Cart::yyz][0]+wmp0*R_temp[Cart::yyyyz][Cart::yyz][1]
R_temp[Cart::xyyyyz][Cart::xxy][0]+=pma0*R_temp[Cart::yyyyz][Cart::xxy][0]+wmp0*R_temp[Cart::yyyyz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yyyyz][Cart::xy][1]
R_temp[Cart::xyyyyz][Cart::xyz][0]+=pma0*R_temp[Cart::yyyyz][Cart::xyz][0]+wmp0*R_temp[Cart::yyyyz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::yz][1]
R_temp[Cart::xyyyyz][Cart::yzz][0]+=pma0*R_temp[Cart::yyyyz][Cart::yzz][0]+wmp0*R_temp[Cart::yyyyz][Cart::yzz][1]
R_temp[Cart::xyyyyz][Cart::xxx][0]+=pma0*R_temp[Cart::yyyyz][Cart::xxx][0]+wmp0*R_temp[Cart::yyyyz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yyyyz][Cart::xx][1]
R_temp[Cart::xyyyyz][Cart::xxz][0]+=pma0*R_temp[Cart::yyyyz][Cart::xxz][0]+wmp0*R_temp[Cart::yyyyz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yyyyz][Cart::xz][1]
R_temp[Cart::xyyyyz][Cart::xzz][0]+=pma0*R_temp[Cart::yyyyz][Cart::xzz][0]+wmp0*R_temp[Cart::yyyyz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::zz][1]
R_temp[Cart::xyyyyz][Cart::zzz][0]+=pma0*R_temp[Cart::yyyyz][Cart::zzz][0]+wmp0*R_temp[Cart::yyyyz][Cart::zzz][1]
R_temp[Cart::yyyyzz][Cart::yyy][0]+=pma2*R_temp[Cart::yyyyz][Cart::yyy][0]+wmp2*R_temp[Cart::yyyyz][Cart::yyy][1]
R_temp[Cart::yyyyzz][Cart::xyy][0]+=pma2*R_temp[Cart::yyyyz][Cart::xyy][0]+wmp2*R_temp[Cart::yyyyz][Cart::xyy][1]
R_temp[Cart::yyyyzz][Cart::yyz][0]+=pma2*R_temp[Cart::yyyyz][Cart::yyz][0]+wmp2*R_temp[Cart::yyyyz][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::yy][1]
R_temp[Cart::yyyyzz][Cart::xxy][0]+=pma2*R_temp[Cart::yyyyz][Cart::xxy][0]+wmp2*R_temp[Cart::yyyyz][Cart::xxy][1]
R_temp[Cart::yyyyzz][Cart::xyz][0]+=pma2*R_temp[Cart::yyyyz][Cart::xyz][0]+wmp2*R_temp[Cart::yyyyz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::xy][1]
R_temp[Cart::yyyyzz][Cart::yzz][0]+=pma2*R_temp[Cart::yyyyz][Cart::yzz][0]+wmp2*R_temp[Cart::yyyyz][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::yyyyz][Cart::yz][1]
R_temp[Cart::yyyyzz][Cart::xxx][0]+=pma2*R_temp[Cart::yyyyz][Cart::xxx][0]+wmp2*R_temp[Cart::yyyyz][Cart::xxx][1]
R_temp[Cart::yyyyzz][Cart::xxz][0]+=pma2*R_temp[Cart::yyyyz][Cart::xxz][0]+wmp2*R_temp[Cart::yyyyz][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::yyyyz][Cart::xx][1]
R_temp[Cart::yyyyzz][Cart::xzz][0]+=pma2*R_temp[Cart::yyyyz][Cart::xzz][0]+wmp2*R_temp[Cart::yyyyz][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::yyyyz][Cart::xz][1]
R_temp[Cart::yyyyzz][Cart::zzz][0]+=pma2*R_temp[Cart::yyyyz][Cart::zzz][0]+wmp2*R_temp[Cart::yyyyz][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::yyyyz][Cart::zz][1]
R_temp[Cart::xxxyyy][Cart::yyy][0]+=pma0*R_temp[Cart::xxyyy][Cart::yyy][0]+wmp0*R_temp[Cart::xxyyy][Cart::yyy][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::yyy][0]-gfak*R_temp[Cart::xyyy][Cart::yyy][1])
R_temp[Cart::xxxyyy][Cart::xyy][0]+=pma0*R_temp[Cart::xxyyy][Cart::xyy][0]+wmp0*R_temp[Cart::xxyyy][Cart::xyy][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xyy][0]-gfak*R_temp[Cart::xyyy][Cart::xyy][1])+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::yy][1]
R_temp[Cart::xxxyyy][Cart::yyz][0]+=pma0*R_temp[Cart::xxyyy][Cart::yyz][0]+wmp0*R_temp[Cart::xxyyy][Cart::yyz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::yyz][0]-gfak*R_temp[Cart::xyyy][Cart::yyz][1])
R_temp[Cart::xxxyyy][Cart::xxy][0]+=pma0*R_temp[Cart::xxyyy][Cart::xxy][0]+wmp0*R_temp[Cart::xxyyy][Cart::xxy][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xxy][0]-gfak*R_temp[Cart::xyyy][Cart::xxy][1])+0.5/_decay*2*R_temp[Cart::xxyyy][Cart::xy][1]
R_temp[Cart::xxxyyy][Cart::xyz][0]+=pma0*R_temp[Cart::xxyyy][Cart::xyz][0]+wmp0*R_temp[Cart::xxyyy][Cart::xyz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xyz][0]-gfak*R_temp[Cart::xyyy][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::yz][1]
R_temp[Cart::xxxyyy][Cart::yzz][0]+=pma0*R_temp[Cart::xxyyy][Cart::yzz][0]+wmp0*R_temp[Cart::xxyyy][Cart::yzz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::yzz][0]-gfak*R_temp[Cart::xyyy][Cart::yzz][1])
R_temp[Cart::xxxyyy][Cart::xxx][0]+=pma0*R_temp[Cart::xxyyy][Cart::xxx][0]+wmp0*R_temp[Cart::xxyyy][Cart::xxx][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xxx][0]-gfak*R_temp[Cart::xyyy][Cart::xxx][1])+0.5/_decay*3*R_temp[Cart::xxyyy][Cart::xx][1]
R_temp[Cart::xxxyyy][Cart::xxz][0]+=pma0*R_temp[Cart::xxyyy][Cart::xxz][0]+wmp0*R_temp[Cart::xxyyy][Cart::xxz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xxz][0]-gfak*R_temp[Cart::xyyy][Cart::xxz][1])+0.5/_decay*2*R_temp[Cart::xxyyy][Cart::xz][1]
R_temp[Cart::xxxyyy][Cart::xzz][0]+=pma0*R_temp[Cart::xxyyy][Cart::xzz][0]+wmp0*R_temp[Cart::xxyyy][Cart::xzz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::xzz][0]-gfak*R_temp[Cart::xyyy][Cart::xzz][1])+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::zz][1]
R_temp[Cart::xxxyyy][Cart::zzz][0]+=pma0*R_temp[Cart::xxyyy][Cart::zzz][0]+wmp0*R_temp[Cart::xxyyy][Cart::zzz][1]+1*rzeta*(R_temp[Cart::xyyy][Cart::zzz][0]-gfak*R_temp[Cart::xyyy][Cart::zzz][1])
R_temp[Cart::xxyyyz][Cart::yyy][0]+=pma2*R_temp[Cart::xxyyy][Cart::yyy][0]+wmp2*R_temp[Cart::xxyyy][Cart::yyy][1]
R_temp[Cart::xxyyyz][Cart::xyy][0]+=pma2*R_temp[Cart::xxyyy][Cart::xyy][0]+wmp2*R_temp[Cart::xxyyy][Cart::xyy][1]
R_temp[Cart::xxyyyz][Cart::yyz][0]+=pma2*R_temp[Cart::xxyyy][Cart::yyz][0]+wmp2*R_temp[Cart::xxyyy][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::yy][1]
R_temp[Cart::xxyyyz][Cart::xxy][0]+=pma2*R_temp[Cart::xxyyy][Cart::xxy][0]+wmp2*R_temp[Cart::xxyyy][Cart::xxy][1]
R_temp[Cart::xxyyyz][Cart::xyz][0]+=pma2*R_temp[Cart::xxyyy][Cart::xyz][0]+wmp2*R_temp[Cart::xxyyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::xy][1]
R_temp[Cart::xxyyyz][Cart::yzz][0]+=pma2*R_temp[Cart::xxyyy][Cart::yzz][0]+wmp2*R_temp[Cart::xxyyy][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::xxyyy][Cart::yz][1]
R_temp[Cart::xxyyyz][Cart::xxx][0]+=pma2*R_temp[Cart::xxyyy][Cart::xxx][0]+wmp2*R_temp[Cart::xxyyy][Cart::xxx][1]
R_temp[Cart::xxyyyz][Cart::xxz][0]+=pma2*R_temp[Cart::xxyyy][Cart::xxz][0]+wmp2*R_temp[Cart::xxyyy][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::xxyyy][Cart::xx][1]
R_temp[Cart::xxyyyz][Cart::xzz][0]+=pma2*R_temp[Cart::xxyyy][Cart::xzz][0]+wmp2*R_temp[Cart::xxyyy][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::xxyyy][Cart::xz][1]
R_temp[Cart::xxyyyz][Cart::zzz][0]+=pma2*R_temp[Cart::xxyyy][Cart::zzz][0]+wmp2*R_temp[Cart::xxyyy][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::xxyyy][Cart::zz][1]
R_temp[Cart::xyyyzz][Cart::yyy][0]+=pma0*R_temp[Cart::yyyzz][Cart::yyy][0]+wmp0*R_temp[Cart::yyyzz][Cart::yyy][1]
R_temp[Cart::xyyyzz][Cart::xyy][0]+=pma0*R_temp[Cart::yyyzz][Cart::xyy][0]+wmp0*R_temp[Cart::yyyzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yyyzz][Cart::yy][1]
R_temp[Cart::xyyyzz][Cart::yyz][0]+=pma0*R_temp[Cart::yyyzz][Cart::yyz][0]+wmp0*R_temp[Cart::yyyzz][Cart::yyz][1]
R_temp[Cart::xyyyzz][Cart::xxy][0]+=pma0*R_temp[Cart::yyyzz][Cart::xxy][0]+wmp0*R_temp[Cart::yyyzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yyyzz][Cart::xy][1]
R_temp[Cart::xyyyzz][Cart::xyz][0]+=pma0*R_temp[Cart::yyyzz][Cart::xyz][0]+wmp0*R_temp[Cart::yyyzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyyzz][Cart::yz][1]
R_temp[Cart::xyyyzz][Cart::yzz][0]+=pma0*R_temp[Cart::yyyzz][Cart::yzz][0]+wmp0*R_temp[Cart::yyyzz][Cart::yzz][1]
R_temp[Cart::xyyyzz][Cart::xxx][0]+=pma0*R_temp[Cart::yyyzz][Cart::xxx][0]+wmp0*R_temp[Cart::yyyzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yyyzz][Cart::xx][1]
R_temp[Cart::xyyyzz][Cart::xxz][0]+=pma0*R_temp[Cart::yyyzz][Cart::xxz][0]+wmp0*R_temp[Cart::yyyzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yyyzz][Cart::xz][1]
R_temp[Cart::xyyyzz][Cart::xzz][0]+=pma0*R_temp[Cart::yyyzz][Cart::xzz][0]+wmp0*R_temp[Cart::yyyzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yyyzz][Cart::zz][1]
R_temp[Cart::xyyyzz][Cart::zzz][0]+=pma0*R_temp[Cart::yyyzz][Cart::zzz][0]+wmp0*R_temp[Cart::yyyzz][Cart::zzz][1]
R_temp[Cart::yyyzzz][Cart::yyy][0]+=pma1*R_temp[Cart::yyzzz][Cart::yyy][0]+wmp1*R_temp[Cart::yyzzz][Cart::yyy][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::yyy][0]-gfak*R_temp[Cart::yzzz][Cart::yyy][1])+0.5/_decay*3*R_temp[Cart::yyzzz][Cart::yy][1]
R_temp[Cart::yyyzzz][Cart::xyy][0]+=pma1*R_temp[Cart::yyzzz][Cart::xyy][0]+wmp1*R_temp[Cart::yyzzz][Cart::xyy][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xyy][0]-gfak*R_temp[Cart::yzzz][Cart::xyy][1])+0.5/_decay*2*R_temp[Cart::yyzzz][Cart::xy][1]
R_temp[Cart::yyyzzz][Cart::yyz][0]+=pma1*R_temp[Cart::yyzzz][Cart::yyz][0]+wmp1*R_temp[Cart::yyzzz][Cart::yyz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::yyz][0]-gfak*R_temp[Cart::yzzz][Cart::yyz][1])+0.5/_decay*2*R_temp[Cart::yyzzz][Cart::yz][1]
R_temp[Cart::yyyzzz][Cart::xxy][0]+=pma1*R_temp[Cart::yyzzz][Cart::xxy][0]+wmp1*R_temp[Cart::yyzzz][Cart::xxy][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xxy][0]-gfak*R_temp[Cart::yzzz][Cart::xxy][1])+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::xx][1]
R_temp[Cart::yyyzzz][Cart::xyz][0]+=pma1*R_temp[Cart::yyzzz][Cart::xyz][0]+wmp1*R_temp[Cart::yyzzz][Cart::xyz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xyz][0]-gfak*R_temp[Cart::yzzz][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::xz][1]
R_temp[Cart::yyyzzz][Cart::yzz][0]+=pma1*R_temp[Cart::yyzzz][Cart::yzz][0]+wmp1*R_temp[Cart::yyzzz][Cart::yzz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::yzz][0]-gfak*R_temp[Cart::yzzz][Cart::yzz][1])+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::zz][1]
R_temp[Cart::yyyzzz][Cart::xxx][0]+=pma1*R_temp[Cart::yyzzz][Cart::xxx][0]+wmp1*R_temp[Cart::yyzzz][Cart::xxx][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xxx][0]-gfak*R_temp[Cart::yzzz][Cart::xxx][1])
R_temp[Cart::yyyzzz][Cart::xxz][0]+=pma1*R_temp[Cart::yyzzz][Cart::xxz][0]+wmp1*R_temp[Cart::yyzzz][Cart::xxz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xxz][0]-gfak*R_temp[Cart::yzzz][Cart::xxz][1])
R_temp[Cart::yyyzzz][Cart::xzz][0]+=pma1*R_temp[Cart::yyzzz][Cart::xzz][0]+wmp1*R_temp[Cart::yyzzz][Cart::xzz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::xzz][0]-gfak*R_temp[Cart::yzzz][Cart::xzz][1])
R_temp[Cart::yyyzzz][Cart::zzz][0]+=pma1*R_temp[Cart::yyzzz][Cart::zzz][0]+wmp1*R_temp[Cart::yyzzz][Cart::zzz][1]+1*rzeta*(R_temp[Cart::yzzz][Cart::zzz][0]-gfak*R_temp[Cart::yzzz][Cart::zzz][1])
R_temp[Cart::xxxxyy][Cart::yyy][0]+=pma1*R_temp[Cart::xxxxy][Cart::yyy][0]+wmp1*R_temp[Cart::xxxxy][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxxxy][Cart::yy][1]
R_temp[Cart::xxxxyy][Cart::xyy][0]+=pma1*R_temp[Cart::xxxxy][Cart::xyy][0]+wmp1*R_temp[Cart::xxxxy][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxxxy][Cart::xy][1]
R_temp[Cart::xxxxyy][Cart::yyz][0]+=pma1*R_temp[Cart::xxxxy][Cart::yyz][0]+wmp1*R_temp[Cart::xxxxy][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxxxy][Cart::yz][1]
R_temp[Cart::xxxxyy][Cart::xxy][0]+=pma1*R_temp[Cart::xxxxy][Cart::xxy][0]+wmp1*R_temp[Cart::xxxxy][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxxxy][Cart::xx][1]
R_temp[Cart::xxxxyy][Cart::xyz][0]+=pma1*R_temp[Cart::xxxxy][Cart::xyz][0]+wmp1*R_temp[Cart::xxxxy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxxy][Cart::xz][1]
R_temp[Cart::xxxxyy][Cart::yzz][0]+=pma1*R_temp[Cart::xxxxy][Cart::yzz][0]+wmp1*R_temp[Cart::xxxxy][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxxxy][Cart::zz][1]
R_temp[Cart::xxxxyy][Cart::xxx][0]+=pma1*R_temp[Cart::xxxxy][Cart::xxx][0]+wmp1*R_temp[Cart::xxxxy][Cart::xxx][1]
R_temp[Cart::xxxxyy][Cart::xxz][0]+=pma1*R_temp[Cart::xxxxy][Cart::xxz][0]+wmp1*R_temp[Cart::xxxxy][Cart::xxz][1]
R_temp[Cart::xxxxyy][Cart::xzz][0]+=pma1*R_temp[Cart::xxxxy][Cart::xzz][0]+wmp1*R_temp[Cart::xxxxy][Cart::xzz][1]
R_temp[Cart::xxxxyy][Cart::zzz][0]+=pma1*R_temp[Cart::xxxxy][Cart::zzz][0]+wmp1*R_temp[Cart::xxxxy][Cart::zzz][1]
R_temp[Cart::xxxyyz][Cart::yyy][0]+=pma2*R_temp[Cart::xxxyy][Cart::yyy][0]+wmp2*R_temp[Cart::xxxyy][Cart::yyy][1]
R_temp[Cart::xxxyyz][Cart::xyy][0]+=pma2*R_temp[Cart::xxxyy][Cart::xyy][0]+wmp2*R_temp[Cart::xxxyy][Cart::xyy][1]
R_temp[Cart::xxxyyz][Cart::yyz][0]+=pma2*R_temp[Cart::xxxyy][Cart::yyz][0]+wmp2*R_temp[Cart::xxxyy][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::xxxyy][Cart::yy][1]
R_temp[Cart::xxxyyz][Cart::xxy][0]+=pma2*R_temp[Cart::xxxyy][Cart::xxy][0]+wmp2*R_temp[Cart::xxxyy][Cart::xxy][1]
R_temp[Cart::xxxyyz][Cart::xyz][0]+=pma2*R_temp[Cart::xxxyy][Cart::xyz][0]+wmp2*R_temp[Cart::xxxyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxyy][Cart::xy][1]
R_temp[Cart::xxxyyz][Cart::yzz][0]+=pma2*R_temp[Cart::xxxyy][Cart::yzz][0]+wmp2*R_temp[Cart::xxxyy][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::xxxyy][Cart::yz][1]
R_temp[Cart::xxxyyz][Cart::xxx][0]+=pma2*R_temp[Cart::xxxyy][Cart::xxx][0]+wmp2*R_temp[Cart::xxxyy][Cart::xxx][1]
R_temp[Cart::xxxyyz][Cart::xxz][0]+=pma2*R_temp[Cart::xxxyy][Cart::xxz][0]+wmp2*R_temp[Cart::xxxyy][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::xxxyy][Cart::xx][1]
R_temp[Cart::xxxyyz][Cart::xzz][0]+=pma2*R_temp[Cart::xxxyy][Cart::xzz][0]+wmp2*R_temp[Cart::xxxyy][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::xxxyy][Cart::xz][1]
R_temp[Cart::xxxyyz][Cart::zzz][0]+=pma2*R_temp[Cart::xxxyy][Cart::zzz][0]+wmp2*R_temp[Cart::xxxyy][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::xxxyy][Cart::zz][1]
R_temp[Cart::xxyyzz][Cart::yyy][0]+=pma0*R_temp[Cart::xyyzz][Cart::yyy][0]+wmp0*R_temp[Cart::xyyzz][Cart::yyy][1]
R_temp[Cart::xxyyzz][Cart::xyy][0]+=pma0*R_temp[Cart::xyyzz][Cart::xyy][0]+wmp0*R_temp[Cart::xyyzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::xyyzz][Cart::yy][1]
R_temp[Cart::xxyyzz][Cart::yyz][0]+=pma0*R_temp[Cart::xyyzz][Cart::yyz][0]+wmp0*R_temp[Cart::xyyzz][Cart::yyz][1]
R_temp[Cart::xxyyzz][Cart::xxy][0]+=pma0*R_temp[Cart::xyyzz][Cart::xxy][0]+wmp0*R_temp[Cart::xyyzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::xyyzz][Cart::xy][1]
R_temp[Cart::xxyyzz][Cart::xyz][0]+=pma0*R_temp[Cart::xyyzz][Cart::xyz][0]+wmp0*R_temp[Cart::xyyzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xyyzz][Cart::yz][1]
R_temp[Cart::xxyyzz][Cart::yzz][0]+=pma0*R_temp[Cart::xyyzz][Cart::yzz][0]+wmp0*R_temp[Cart::xyyzz][Cart::yzz][1]
R_temp[Cart::xxyyzz][Cart::xxx][0]+=pma0*R_temp[Cart::xyyzz][Cart::xxx][0]+wmp0*R_temp[Cart::xyyzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::xyyzz][Cart::xx][1]
R_temp[Cart::xxyyzz][Cart::xxz][0]+=pma0*R_temp[Cart::xyyzz][Cart::xxz][0]+wmp0*R_temp[Cart::xyyzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::xyyzz][Cart::xz][1]
R_temp[Cart::xxyyzz][Cart::xzz][0]+=pma0*R_temp[Cart::xyyzz][Cart::xzz][0]+wmp0*R_temp[Cart::xyyzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::xyyzz][Cart::zz][1]
R_temp[Cart::xxyyzz][Cart::zzz][0]+=pma0*R_temp[Cart::xyyzz][Cart::zzz][0]+wmp0*R_temp[Cart::xyyzz][Cart::zzz][1]
R_temp[Cart::xyyzzz][Cart::yyy][0]+=pma0*R_temp[Cart::yyzzz][Cart::yyy][0]+wmp0*R_temp[Cart::yyzzz][Cart::yyy][1]
R_temp[Cart::xyyzzz][Cart::xyy][0]+=pma0*R_temp[Cart::yyzzz][Cart::xyy][0]+wmp0*R_temp[Cart::yyzzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::yy][1]
R_temp[Cart::xyyzzz][Cart::yyz][0]+=pma0*R_temp[Cart::yyzzz][Cart::yyz][0]+wmp0*R_temp[Cart::yyzzz][Cart::yyz][1]
R_temp[Cart::xyyzzz][Cart::xxy][0]+=pma0*R_temp[Cart::yyzzz][Cart::xxy][0]+wmp0*R_temp[Cart::yyzzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yyzzz][Cart::xy][1]
R_temp[Cart::xyyzzz][Cart::xyz][0]+=pma0*R_temp[Cart::yyzzz][Cart::xyz][0]+wmp0*R_temp[Cart::yyzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::yz][1]
R_temp[Cart::xyyzzz][Cart::yzz][0]+=pma0*R_temp[Cart::yyzzz][Cart::yzz][0]+wmp0*R_temp[Cart::yyzzz][Cart::yzz][1]
R_temp[Cart::xyyzzz][Cart::xxx][0]+=pma0*R_temp[Cart::yyzzz][Cart::xxx][0]+wmp0*R_temp[Cart::yyzzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yyzzz][Cart::xx][1]
R_temp[Cart::xyyzzz][Cart::xxz][0]+=pma0*R_temp[Cart::yyzzz][Cart::xxz][0]+wmp0*R_temp[Cart::yyzzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yyzzz][Cart::xz][1]
R_temp[Cart::xyyzzz][Cart::xzz][0]+=pma0*R_temp[Cart::yyzzz][Cart::xzz][0]+wmp0*R_temp[Cart::yyzzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yyzzz][Cart::zz][1]
R_temp[Cart::xyyzzz][Cart::zzz][0]+=pma0*R_temp[Cart::yyzzz][Cart::zzz][0]+wmp0*R_temp[Cart::yyzzz][Cart::zzz][1]
R_temp[Cart::yyzzzz][Cart::yyy][0]+=pma1*R_temp[Cart::yzzzz][Cart::yyy][0]+wmp1*R_temp[Cart::yzzzz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::yzzzz][Cart::yy][1]
R_temp[Cart::yyzzzz][Cart::xyy][0]+=pma1*R_temp[Cart::yzzzz][Cart::xyy][0]+wmp1*R_temp[Cart::yzzzz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::yzzzz][Cart::xy][1]
R_temp[Cart::yyzzzz][Cart::yyz][0]+=pma1*R_temp[Cart::yzzzz][Cart::yyz][0]+wmp1*R_temp[Cart::yzzzz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::yzzzz][Cart::yz][1]
R_temp[Cart::yyzzzz][Cart::xxy][0]+=pma1*R_temp[Cart::yzzzz][Cart::xxy][0]+wmp1*R_temp[Cart::yzzzz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::xx][1]
R_temp[Cart::yyzzzz][Cart::xyz][0]+=pma1*R_temp[Cart::yzzzz][Cart::xyz][0]+wmp1*R_temp[Cart::yzzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::xz][1]
R_temp[Cart::yyzzzz][Cart::yzz][0]+=pma1*R_temp[Cart::yzzzz][Cart::yzz][0]+wmp1*R_temp[Cart::yzzzz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::zz][1]
R_temp[Cart::yyzzzz][Cart::xxx][0]+=pma1*R_temp[Cart::yzzzz][Cart::xxx][0]+wmp1*R_temp[Cart::yzzzz][Cart::xxx][1]
R_temp[Cart::yyzzzz][Cart::xxz][0]+=pma1*R_temp[Cart::yzzzz][Cart::xxz][0]+wmp1*R_temp[Cart::yzzzz][Cart::xxz][1]
R_temp[Cart::yyzzzz][Cart::xzz][0]+=pma1*R_temp[Cart::yzzzz][Cart::xzz][0]+wmp1*R_temp[Cart::yzzzz][Cart::xzz][1]
R_temp[Cart::yyzzzz][Cart::zzz][0]+=pma1*R_temp[Cart::yzzzz][Cart::zzz][0]+wmp1*R_temp[Cart::yzzzz][Cart::zzz][1]
R_temp[Cart::xxxxxy][Cart::yyy][0]+=pma1*R_temp[Cart::xxxxx][Cart::yyy][0]+wmp1*R_temp[Cart::xxxxx][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxxxx][Cart::yy][1]
R_temp[Cart::xxxxxy][Cart::xyy][0]+=pma1*R_temp[Cart::xxxxx][Cart::xyy][0]+wmp1*R_temp[Cart::xxxxx][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::xy][1]
R_temp[Cart::xxxxxy][Cart::yyz][0]+=pma1*R_temp[Cart::xxxxx][Cart::yyz][0]+wmp1*R_temp[Cart::xxxxx][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::yz][1]
R_temp[Cart::xxxxxy][Cart::xxy][0]+=pma1*R_temp[Cart::xxxxx][Cart::xxy][0]+wmp1*R_temp[Cart::xxxxx][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::xx][1]
R_temp[Cart::xxxxxy][Cart::xyz][0]+=pma1*R_temp[Cart::xxxxx][Cart::xyz][0]+wmp1*R_temp[Cart::xxxxx][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::xz][1]
R_temp[Cart::xxxxxy][Cart::yzz][0]+=pma1*R_temp[Cart::xxxxx][Cart::yzz][0]+wmp1*R_temp[Cart::xxxxx][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::zz][1]
R_temp[Cart::xxxxxy][Cart::xxx][0]+=pma1*R_temp[Cart::xxxxx][Cart::xxx][0]+wmp1*R_temp[Cart::xxxxx][Cart::xxx][1]
R_temp[Cart::xxxxxy][Cart::xxz][0]+=pma1*R_temp[Cart::xxxxx][Cart::xxz][0]+wmp1*R_temp[Cart::xxxxx][Cart::xxz][1]
R_temp[Cart::xxxxxy][Cart::xzz][0]+=pma1*R_temp[Cart::xxxxx][Cart::xzz][0]+wmp1*R_temp[Cart::xxxxx][Cart::xzz][1]
R_temp[Cart::xxxxxy][Cart::zzz][0]+=pma1*R_temp[Cart::xxxxx][Cart::zzz][0]+wmp1*R_temp[Cart::xxxxx][Cart::zzz][1]
R_temp[Cart::xxxxyz][Cart::yyy][0]+=pma1*R_temp[Cart::xxxxz][Cart::yyy][0]+wmp1*R_temp[Cart::xxxxz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxxxz][Cart::yy][1]
R_temp[Cart::xxxxyz][Cart::xyy][0]+=pma1*R_temp[Cart::xxxxz][Cart::xyy][0]+wmp1*R_temp[Cart::xxxxz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxxxz][Cart::xy][1]
R_temp[Cart::xxxxyz][Cart::yyz][0]+=pma1*R_temp[Cart::xxxxz][Cart::yyz][0]+wmp1*R_temp[Cart::xxxxz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxxxz][Cart::yz][1]
R_temp[Cart::xxxxyz][Cart::xxy][0]+=pma1*R_temp[Cart::xxxxz][Cart::xxy][0]+wmp1*R_temp[Cart::xxxxz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::xx][1]
R_temp[Cart::xxxxyz][Cart::xyz][0]+=pma1*R_temp[Cart::xxxxz][Cart::xyz][0]+wmp1*R_temp[Cart::xxxxz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::xz][1]
R_temp[Cart::xxxxyz][Cart::yzz][0]+=pma1*R_temp[Cart::xxxxz][Cart::yzz][0]+wmp1*R_temp[Cart::xxxxz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::zz][1]
R_temp[Cart::xxxxyz][Cart::xxx][0]+=pma1*R_temp[Cart::xxxxz][Cart::xxx][0]+wmp1*R_temp[Cart::xxxxz][Cart::xxx][1]
R_temp[Cart::xxxxyz][Cart::xxz][0]+=pma1*R_temp[Cart::xxxxz][Cart::xxz][0]+wmp1*R_temp[Cart::xxxxz][Cart::xxz][1]
R_temp[Cart::xxxxyz][Cart::xzz][0]+=pma1*R_temp[Cart::xxxxz][Cart::xzz][0]+wmp1*R_temp[Cart::xxxxz][Cart::xzz][1]
R_temp[Cart::xxxxyz][Cart::zzz][0]+=pma1*R_temp[Cart::xxxxz][Cart::zzz][0]+wmp1*R_temp[Cart::xxxxz][Cart::zzz][1]
R_temp[Cart::xxxyzz][Cart::yyy][0]+=pma1*R_temp[Cart::xxxzz][Cart::yyy][0]+wmp1*R_temp[Cart::xxxzz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxxzz][Cart::yy][1]
R_temp[Cart::xxxyzz][Cart::xyy][0]+=pma1*R_temp[Cart::xxxzz][Cart::xyy][0]+wmp1*R_temp[Cart::xxxzz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxxzz][Cart::xy][1]
R_temp[Cart::xxxyzz][Cart::yyz][0]+=pma1*R_temp[Cart::xxxzz][Cart::yyz][0]+wmp1*R_temp[Cart::xxxzz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxxzz][Cart::yz][1]
R_temp[Cart::xxxyzz][Cart::xxy][0]+=pma1*R_temp[Cart::xxxzz][Cart::xxy][0]+wmp1*R_temp[Cart::xxxzz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxxzz][Cart::xx][1]
R_temp[Cart::xxxyzz][Cart::xyz][0]+=pma1*R_temp[Cart::xxxzz][Cart::xyz][0]+wmp1*R_temp[Cart::xxxzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxzz][Cart::xz][1]
R_temp[Cart::xxxyzz][Cart::yzz][0]+=pma1*R_temp[Cart::xxxzz][Cart::yzz][0]+wmp1*R_temp[Cart::xxxzz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxxzz][Cart::zz][1]
R_temp[Cart::xxxyzz][Cart::xxx][0]+=pma1*R_temp[Cart::xxxzz][Cart::xxx][0]+wmp1*R_temp[Cart::xxxzz][Cart::xxx][1]
R_temp[Cart::xxxyzz][Cart::xxz][0]+=pma1*R_temp[Cart::xxxzz][Cart::xxz][0]+wmp1*R_temp[Cart::xxxzz][Cart::xxz][1]
R_temp[Cart::xxxyzz][Cart::xzz][0]+=pma1*R_temp[Cart::xxxzz][Cart::xzz][0]+wmp1*R_temp[Cart::xxxzz][Cart::xzz][1]
R_temp[Cart::xxxyzz][Cart::zzz][0]+=pma1*R_temp[Cart::xxxzz][Cart::zzz][0]+wmp1*R_temp[Cart::xxxzz][Cart::zzz][1]
R_temp[Cart::xxyzzz][Cart::yyy][0]+=pma1*R_temp[Cart::xxzzz][Cart::yyy][0]+wmp1*R_temp[Cart::xxzzz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxzzz][Cart::yy][1]
R_temp[Cart::xxyzzz][Cart::xyy][0]+=pma1*R_temp[Cart::xxzzz][Cart::xyy][0]+wmp1*R_temp[Cart::xxzzz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxzzz][Cart::xy][1]
R_temp[Cart::xxyzzz][Cart::yyz][0]+=pma1*R_temp[Cart::xxzzz][Cart::yyz][0]+wmp1*R_temp[Cart::xxzzz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxzzz][Cart::yz][1]
R_temp[Cart::xxyzzz][Cart::xxy][0]+=pma1*R_temp[Cart::xxzzz][Cart::xxy][0]+wmp1*R_temp[Cart::xxzzz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::xx][1]
R_temp[Cart::xxyzzz][Cart::xyz][0]+=pma1*R_temp[Cart::xxzzz][Cart::xyz][0]+wmp1*R_temp[Cart::xxzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::xz][1]
R_temp[Cart::xxyzzz][Cart::yzz][0]+=pma1*R_temp[Cart::xxzzz][Cart::yzz][0]+wmp1*R_temp[Cart::xxzzz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::zz][1]
R_temp[Cart::xxyzzz][Cart::xxx][0]+=pma1*R_temp[Cart::xxzzz][Cart::xxx][0]+wmp1*R_temp[Cart::xxzzz][Cart::xxx][1]
R_temp[Cart::xxyzzz][Cart::xxz][0]+=pma1*R_temp[Cart::xxzzz][Cart::xxz][0]+wmp1*R_temp[Cart::xxzzz][Cart::xxz][1]
R_temp[Cart::xxyzzz][Cart::xzz][0]+=pma1*R_temp[Cart::xxzzz][Cart::xzz][0]+wmp1*R_temp[Cart::xxzzz][Cart::xzz][1]
R_temp[Cart::xxyzzz][Cart::zzz][0]+=pma1*R_temp[Cart::xxzzz][Cart::zzz][0]+wmp1*R_temp[Cart::xxzzz][Cart::zzz][1]
R_temp[Cart::xyzzzz][Cart::yyy][0]+=pma0*R_temp[Cart::yzzzz][Cart::yyy][0]+wmp0*R_temp[Cart::yzzzz][Cart::yyy][1]
R_temp[Cart::xyzzzz][Cart::xyy][0]+=pma0*R_temp[Cart::yzzzz][Cart::xyy][0]+wmp0*R_temp[Cart::yzzzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::yy][1]
R_temp[Cart::xyzzzz][Cart::yyz][0]+=pma0*R_temp[Cart::yzzzz][Cart::yyz][0]+wmp0*R_temp[Cart::yzzzz][Cart::yyz][1]
R_temp[Cart::xyzzzz][Cart::xxy][0]+=pma0*R_temp[Cart::yzzzz][Cart::xxy][0]+wmp0*R_temp[Cart::yzzzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yzzzz][Cart::xy][1]
R_temp[Cart::xyzzzz][Cart::xyz][0]+=pma0*R_temp[Cart::yzzzz][Cart::xyz][0]+wmp0*R_temp[Cart::yzzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::yz][1]
R_temp[Cart::xyzzzz][Cart::yzz][0]+=pma0*R_temp[Cart::yzzzz][Cart::yzz][0]+wmp0*R_temp[Cart::yzzzz][Cart::yzz][1]
R_temp[Cart::xyzzzz][Cart::xxx][0]+=pma0*R_temp[Cart::yzzzz][Cart::xxx][0]+wmp0*R_temp[Cart::yzzzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yzzzz][Cart::xx][1]
R_temp[Cart::xyzzzz][Cart::xxz][0]+=pma0*R_temp[Cart::yzzzz][Cart::xxz][0]+wmp0*R_temp[Cart::yzzzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yzzzz][Cart::xz][1]
R_temp[Cart::xyzzzz][Cart::xzz][0]+=pma0*R_temp[Cart::yzzzz][Cart::xzz][0]+wmp0*R_temp[Cart::yzzzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yzzzz][Cart::zz][1]
R_temp[Cart::xyzzzz][Cart::zzz][0]+=pma0*R_temp[Cart::yzzzz][Cart::zzz][0]+wmp0*R_temp[Cart::yzzzz][Cart::zzz][1]
R_temp[Cart::yzzzzz][Cart::yyy][0]+=pma1*R_temp[Cart::zzzzz][Cart::yyy][0]+wmp1*R_temp[Cart::zzzzz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::zzzzz][Cart::yy][1]
R_temp[Cart::yzzzzz][Cart::xyy][0]+=pma1*R_temp[Cart::zzzzz][Cart::xyy][0]+wmp1*R_temp[Cart::zzzzz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::xy][1]
R_temp[Cart::yzzzzz][Cart::yyz][0]+=pma1*R_temp[Cart::zzzzz][Cart::yyz][0]+wmp1*R_temp[Cart::zzzzz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::yz][1]
R_temp[Cart::yzzzzz][Cart::xxy][0]+=pma1*R_temp[Cart::zzzzz][Cart::xxy][0]+wmp1*R_temp[Cart::zzzzz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::xx][1]
R_temp[Cart::yzzzzz][Cart::xyz][0]+=pma1*R_temp[Cart::zzzzz][Cart::xyz][0]+wmp1*R_temp[Cart::zzzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::xz][1]
R_temp[Cart::yzzzzz][Cart::yzz][0]+=pma1*R_temp[Cart::zzzzz][Cart::yzz][0]+wmp1*R_temp[Cart::zzzzz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::zz][1]
R_temp[Cart::yzzzzz][Cart::xxx][0]+=pma1*R_temp[Cart::zzzzz][Cart::xxx][0]+wmp1*R_temp[Cart::zzzzz][Cart::xxx][1]
R_temp[Cart::yzzzzz][Cart::xxz][0]+=pma1*R_temp[Cart::zzzzz][Cart::xxz][0]+wmp1*R_temp[Cart::zzzzz][Cart::xxz][1]
R_temp[Cart::yzzzzz][Cart::xzz][0]+=pma1*R_temp[Cart::zzzzz][Cart::xzz][0]+wmp1*R_temp[Cart::zzzzz][Cart::xzz][1]
R_temp[Cart::yzzzzz][Cart::zzz][0]+=pma1*R_temp[Cart::zzzzz][Cart::zzz][0]+wmp1*R_temp[Cart::zzzzz][Cart::zzz][1]
R_temp[Cart::xxxxxx][Cart::yyy][0]+=pma0*R_temp[Cart::xxxxx][Cart::yyy][0]+wmp0*R_temp[Cart::xxxxx][Cart::yyy][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::yyy][0]-gfak*R_temp[Cart::xxxx][Cart::yyy][1])
R_temp[Cart::xxxxxx][Cart::xyy][0]+=pma0*R_temp[Cart::xxxxx][Cart::xyy][0]+wmp0*R_temp[Cart::xxxxx][Cart::xyy][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xyy][0]-gfak*R_temp[Cart::xxxx][Cart::xyy][1])+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::yy][1]
R_temp[Cart::xxxxxx][Cart::yyz][0]+=pma0*R_temp[Cart::xxxxx][Cart::yyz][0]+wmp0*R_temp[Cart::xxxxx][Cart::yyz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::yyz][0]-gfak*R_temp[Cart::xxxx][Cart::yyz][1])
R_temp[Cart::xxxxxx][Cart::xxy][0]+=pma0*R_temp[Cart::xxxxx][Cart::xxy][0]+wmp0*R_temp[Cart::xxxxx][Cart::xxy][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xxy][0]-gfak*R_temp[Cart::xxxx][Cart::xxy][1])+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::xy][1]
R_temp[Cart::xxxxxx][Cart::xyz][0]+=pma0*R_temp[Cart::xxxxx][Cart::xyz][0]+wmp0*R_temp[Cart::xxxxx][Cart::xyz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xyz][0]-gfak*R_temp[Cart::xxxx][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::yz][1]
R_temp[Cart::xxxxxx][Cart::yzz][0]+=pma0*R_temp[Cart::xxxxx][Cart::yzz][0]+wmp0*R_temp[Cart::xxxxx][Cart::yzz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::yzz][0]-gfak*R_temp[Cart::xxxx][Cart::yzz][1])
R_temp[Cart::xxxxxx][Cart::xxx][0]+=pma0*R_temp[Cart::xxxxx][Cart::xxx][0]+wmp0*R_temp[Cart::xxxxx][Cart::xxx][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xxx][0]-gfak*R_temp[Cart::xxxx][Cart::xxx][1])+0.5/_decay*3*R_temp[Cart::xxxxx][Cart::xx][1]
R_temp[Cart::xxxxxx][Cart::xxz][0]+=pma0*R_temp[Cart::xxxxx][Cart::xxz][0]+wmp0*R_temp[Cart::xxxxx][Cart::xxz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xxz][0]-gfak*R_temp[Cart::xxxx][Cart::xxz][1])+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::xz][1]
R_temp[Cart::xxxxxx][Cart::xzz][0]+=pma0*R_temp[Cart::xxxxx][Cart::xzz][0]+wmp0*R_temp[Cart::xxxxx][Cart::xzz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::xzz][0]-gfak*R_temp[Cart::xxxx][Cart::xzz][1])+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::zz][1]
R_temp[Cart::xxxxxx][Cart::zzz][0]+=pma0*R_temp[Cart::xxxxx][Cart::zzz][0]+wmp0*R_temp[Cart::xxxxx][Cart::zzz][1]+4*rzeta*(R_temp[Cart::xxxx][Cart::zzz][0]-gfak*R_temp[Cart::xxxx][Cart::zzz][1])
R_temp[Cart::xxxxxz][Cart::yyy][0]+=pma2*R_temp[Cart::xxxxx][Cart::yyy][0]+wmp2*R_temp[Cart::xxxxx][Cart::yyy][1]
R_temp[Cart::xxxxxz][Cart::xyy][0]+=pma2*R_temp[Cart::xxxxx][Cart::xyy][0]+wmp2*R_temp[Cart::xxxxx][Cart::xyy][1]
R_temp[Cart::xxxxxz][Cart::yyz][0]+=pma2*R_temp[Cart::xxxxx][Cart::yyz][0]+wmp2*R_temp[Cart::xxxxx][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::yy][1]
R_temp[Cart::xxxxxz][Cart::xxy][0]+=pma2*R_temp[Cart::xxxxx][Cart::xxy][0]+wmp2*R_temp[Cart::xxxxx][Cart::xxy][1]
R_temp[Cart::xxxxxz][Cart::xyz][0]+=pma2*R_temp[Cart::xxxxx][Cart::xyz][0]+wmp2*R_temp[Cart::xxxxx][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::xy][1]
R_temp[Cart::xxxxxz][Cart::yzz][0]+=pma2*R_temp[Cart::xxxxx][Cart::yzz][0]+wmp2*R_temp[Cart::xxxxx][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::yz][1]
R_temp[Cart::xxxxxz][Cart::xxx][0]+=pma2*R_temp[Cart::xxxxx][Cart::xxx][0]+wmp2*R_temp[Cart::xxxxx][Cart::xxx][1]
R_temp[Cart::xxxxxz][Cart::xxz][0]+=pma2*R_temp[Cart::xxxxx][Cart::xxz][0]+wmp2*R_temp[Cart::xxxxx][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::xxxxx][Cart::xx][1]
R_temp[Cart::xxxxxz][Cart::xzz][0]+=pma2*R_temp[Cart::xxxxx][Cart::xzz][0]+wmp2*R_temp[Cart::xxxxx][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::xxxxx][Cart::xz][1]
R_temp[Cart::xxxxxz][Cart::zzz][0]+=pma2*R_temp[Cart::xxxxx][Cart::zzz][0]+wmp2*R_temp[Cart::xxxxx][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::xxxxx][Cart::zz][1]
R_temp[Cart::xxxxzz][Cart::yyy][0]+=pma2*R_temp[Cart::xxxxz][Cart::yyy][0]+wmp2*R_temp[Cart::xxxxz][Cart::yyy][1]
R_temp[Cart::xxxxzz][Cart::xyy][0]+=pma2*R_temp[Cart::xxxxz][Cart::xyy][0]+wmp2*R_temp[Cart::xxxxz][Cart::xyy][1]
R_temp[Cart::xxxxzz][Cart::yyz][0]+=pma2*R_temp[Cart::xxxxz][Cart::yyz][0]+wmp2*R_temp[Cart::xxxxz][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::yy][1]
R_temp[Cart::xxxxzz][Cart::xxy][0]+=pma2*R_temp[Cart::xxxxz][Cart::xxy][0]+wmp2*R_temp[Cart::xxxxz][Cart::xxy][1]
R_temp[Cart::xxxxzz][Cart::xyz][0]+=pma2*R_temp[Cart::xxxxz][Cart::xyz][0]+wmp2*R_temp[Cart::xxxxz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::xy][1]
R_temp[Cart::xxxxzz][Cart::yzz][0]+=pma2*R_temp[Cart::xxxxz][Cart::yzz][0]+wmp2*R_temp[Cart::xxxxz][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::xxxxz][Cart::yz][1]
R_temp[Cart::xxxxzz][Cart::xxx][0]+=pma2*R_temp[Cart::xxxxz][Cart::xxx][0]+wmp2*R_temp[Cart::xxxxz][Cart::xxx][1]
R_temp[Cart::xxxxzz][Cart::xxz][0]+=pma2*R_temp[Cart::xxxxz][Cart::xxz][0]+wmp2*R_temp[Cart::xxxxz][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::xxxxz][Cart::xx][1]
R_temp[Cart::xxxxzz][Cart::xzz][0]+=pma2*R_temp[Cart::xxxxz][Cart::xzz][0]+wmp2*R_temp[Cart::xxxxz][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::xxxxz][Cart::xz][1]
R_temp[Cart::xxxxzz][Cart::zzz][0]+=pma2*R_temp[Cart::xxxxz][Cart::zzz][0]+wmp2*R_temp[Cart::xxxxz][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::xxxxz][Cart::zz][1]
R_temp[Cart::xxxzzz][Cart::yyy][0]+=pma0*R_temp[Cart::xxzzz][Cart::yyy][0]+wmp0*R_temp[Cart::xxzzz][Cart::yyy][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::yyy][0]-gfak*R_temp[Cart::xzzz][Cart::yyy][1])
R_temp[Cart::xxxzzz][Cart::xyy][0]+=pma0*R_temp[Cart::xxzzz][Cart::xyy][0]+wmp0*R_temp[Cart::xxzzz][Cart::xyy][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xyy][0]-gfak*R_temp[Cart::xzzz][Cart::xyy][1])+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::yy][1]
R_temp[Cart::xxxzzz][Cart::yyz][0]+=pma0*R_temp[Cart::xxzzz][Cart::yyz][0]+wmp0*R_temp[Cart::xxzzz][Cart::yyz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::yyz][0]-gfak*R_temp[Cart::xzzz][Cart::yyz][1])
R_temp[Cart::xxxzzz][Cart::xxy][0]+=pma0*R_temp[Cart::xxzzz][Cart::xxy][0]+wmp0*R_temp[Cart::xxzzz][Cart::xxy][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xxy][0]-gfak*R_temp[Cart::xzzz][Cart::xxy][1])+0.5/_decay*2*R_temp[Cart::xxzzz][Cart::xy][1]
R_temp[Cart::xxxzzz][Cart::xyz][0]+=pma0*R_temp[Cart::xxzzz][Cart::xyz][0]+wmp0*R_temp[Cart::xxzzz][Cart::xyz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xyz][0]-gfak*R_temp[Cart::xzzz][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::yz][1]
R_temp[Cart::xxxzzz][Cart::yzz][0]+=pma0*R_temp[Cart::xxzzz][Cart::yzz][0]+wmp0*R_temp[Cart::xxzzz][Cart::yzz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::yzz][0]-gfak*R_temp[Cart::xzzz][Cart::yzz][1])
R_temp[Cart::xxxzzz][Cart::xxx][0]+=pma0*R_temp[Cart::xxzzz][Cart::xxx][0]+wmp0*R_temp[Cart::xxzzz][Cart::xxx][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xxx][0]-gfak*R_temp[Cart::xzzz][Cart::xxx][1])+0.5/_decay*3*R_temp[Cart::xxzzz][Cart::xx][1]
R_temp[Cart::xxxzzz][Cart::xxz][0]+=pma0*R_temp[Cart::xxzzz][Cart::xxz][0]+wmp0*R_temp[Cart::xxzzz][Cart::xxz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xxz][0]-gfak*R_temp[Cart::xzzz][Cart::xxz][1])+0.5/_decay*2*R_temp[Cart::xxzzz][Cart::xz][1]
R_temp[Cart::xxxzzz][Cart::xzz][0]+=pma0*R_temp[Cart::xxzzz][Cart::xzz][0]+wmp0*R_temp[Cart::xxzzz][Cart::xzz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::xzz][0]-gfak*R_temp[Cart::xzzz][Cart::xzz][1])+0.5/_decay*1*R_temp[Cart::xxzzz][Cart::zz][1]
R_temp[Cart::xxxzzz][Cart::zzz][0]+=pma0*R_temp[Cart::xxzzz][Cart::zzz][0]+wmp0*R_temp[Cart::xxzzz][Cart::zzz][1]+1*rzeta*(R_temp[Cart::xzzz][Cart::zzz][0]-gfak*R_temp[Cart::xzzz][Cart::zzz][1])
R_temp[Cart::xxzzzz][Cart::yyy][0]+=pma0*R_temp[Cart::xzzzz][Cart::yyy][0]+wmp0*R_temp[Cart::xzzzz][Cart::yyy][1]
R_temp[Cart::xxzzzz][Cart::xyy][0]+=pma0*R_temp[Cart::xzzzz][Cart::xyy][0]+wmp0*R_temp[Cart::xzzzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::xzzzz][Cart::yy][1]
R_temp[Cart::xxzzzz][Cart::yyz][0]+=pma0*R_temp[Cart::xzzzz][Cart::yyz][0]+wmp0*R_temp[Cart::xzzzz][Cart::yyz][1]
R_temp[Cart::xxzzzz][Cart::xxy][0]+=pma0*R_temp[Cart::xzzzz][Cart::xxy][0]+wmp0*R_temp[Cart::xzzzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::xzzzz][Cart::xy][1]
R_temp[Cart::xxzzzz][Cart::xyz][0]+=pma0*R_temp[Cart::xzzzz][Cart::xyz][0]+wmp0*R_temp[Cart::xzzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xzzzz][Cart::yz][1]
R_temp[Cart::xxzzzz][Cart::yzz][0]+=pma0*R_temp[Cart::xzzzz][Cart::yzz][0]+wmp0*R_temp[Cart::xzzzz][Cart::yzz][1]
R_temp[Cart::xxzzzz][Cart::xxx][0]+=pma0*R_temp[Cart::xzzzz][Cart::xxx][0]+wmp0*R_temp[Cart::xzzzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::xzzzz][Cart::xx][1]
R_temp[Cart::xxzzzz][Cart::xxz][0]+=pma0*R_temp[Cart::xzzzz][Cart::xxz][0]+wmp0*R_temp[Cart::xzzzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::xzzzz][Cart::xz][1]
R_temp[Cart::xxzzzz][Cart::xzz][0]+=pma0*R_temp[Cart::xzzzz][Cart::xzz][0]+wmp0*R_temp[Cart::xzzzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::xzzzz][Cart::zz][1]
R_temp[Cart::xxzzzz][Cart::zzz][0]+=pma0*R_temp[Cart::xzzzz][Cart::zzz][0]+wmp0*R_temp[Cart::xzzzz][Cart::zzz][1]
R_temp[Cart::xzzzzz][Cart::yyy][0]+=pma0*R_temp[Cart::zzzzz][Cart::yyy][0]+wmp0*R_temp[Cart::zzzzz][Cart::yyy][1]
R_temp[Cart::xzzzzz][Cart::xyy][0]+=pma0*R_temp[Cart::zzzzz][Cart::xyy][0]+wmp0*R_temp[Cart::zzzzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::yy][1]
R_temp[Cart::xzzzzz][Cart::yyz][0]+=pma0*R_temp[Cart::zzzzz][Cart::yyz][0]+wmp0*R_temp[Cart::zzzzz][Cart::yyz][1]
R_temp[Cart::xzzzzz][Cart::xxy][0]+=pma0*R_temp[Cart::zzzzz][Cart::xxy][0]+wmp0*R_temp[Cart::zzzzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::xy][1]
R_temp[Cart::xzzzzz][Cart::xyz][0]+=pma0*R_temp[Cart::zzzzz][Cart::xyz][0]+wmp0*R_temp[Cart::zzzzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::yz][1]
R_temp[Cart::xzzzzz][Cart::yzz][0]+=pma0*R_temp[Cart::zzzzz][Cart::yzz][0]+wmp0*R_temp[Cart::zzzzz][Cart::yzz][1]
R_temp[Cart::xzzzzz][Cart::xxx][0]+=pma0*R_temp[Cart::zzzzz][Cart::xxx][0]+wmp0*R_temp[Cart::zzzzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::zzzzz][Cart::xx][1]
R_temp[Cart::xzzzzz][Cart::xxz][0]+=pma0*R_temp[Cart::zzzzz][Cart::xxz][0]+wmp0*R_temp[Cart::zzzzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::xz][1]
R_temp[Cart::xzzzzz][Cart::xzz][0]+=pma0*R_temp[Cart::zzzzz][Cart::xzz][0]+wmp0*R_temp[Cart::zzzzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::zz][1]
R_temp[Cart::xzzzzz][Cart::zzz][0]+=pma0*R_temp[Cart::zzzzz][Cart::zzz][0]+wmp0*R_temp[Cart::zzzzz][Cart::zzz][1]
R_temp[Cart::zzzzzz][Cart::yyy][0]+=pma2*R_temp[Cart::zzzzz][Cart::yyy][0]+wmp2*R_temp[Cart::zzzzz][Cart::yyy][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::yyy][0]-gfak*R_temp[Cart::zzzz][Cart::yyy][1])
R_temp[Cart::zzzzzz][Cart::xyy][0]+=pma2*R_temp[Cart::zzzzz][Cart::xyy][0]+wmp2*R_temp[Cart::zzzzz][Cart::xyy][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xyy][0]-gfak*R_temp[Cart::zzzz][Cart::xyy][1])
R_temp[Cart::zzzzzz][Cart::yyz][0]+=pma2*R_temp[Cart::zzzzz][Cart::yyz][0]+wmp2*R_temp[Cart::zzzzz][Cart::yyz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::yyz][0]-gfak*R_temp[Cart::zzzz][Cart::yyz][1])+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::yy][1]
R_temp[Cart::zzzzzz][Cart::xxy][0]+=pma2*R_temp[Cart::zzzzz][Cart::xxy][0]+wmp2*R_temp[Cart::zzzzz][Cart::xxy][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xxy][0]-gfak*R_temp[Cart::zzzz][Cart::xxy][1])
R_temp[Cart::zzzzzz][Cart::xyz][0]+=pma2*R_temp[Cart::zzzzz][Cart::xyz][0]+wmp2*R_temp[Cart::zzzzz][Cart::xyz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xyz][0]-gfak*R_temp[Cart::zzzz][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::xy][1]
R_temp[Cart::zzzzzz][Cart::yzz][0]+=pma2*R_temp[Cart::zzzzz][Cart::yzz][0]+wmp2*R_temp[Cart::zzzzz][Cart::yzz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::yzz][0]-gfak*R_temp[Cart::zzzz][Cart::yzz][1])+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::yz][1]
R_temp[Cart::zzzzzz][Cart::xxx][0]+=pma2*R_temp[Cart::zzzzz][Cart::xxx][0]+wmp2*R_temp[Cart::zzzzz][Cart::xxx][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xxx][0]-gfak*R_temp[Cart::zzzz][Cart::xxx][1])
R_temp[Cart::zzzzzz][Cart::xxz][0]+=pma2*R_temp[Cart::zzzzz][Cart::xxz][0]+wmp2*R_temp[Cart::zzzzz][Cart::xxz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xxz][0]-gfak*R_temp[Cart::zzzz][Cart::xxz][1])+0.5/_decay*1*R_temp[Cart::zzzzz][Cart::xx][1]
R_temp[Cart::zzzzzz][Cart::xzz][0]+=pma2*R_temp[Cart::zzzzz][Cart::xzz][0]+wmp2*R_temp[Cart::zzzzz][Cart::xzz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::xzz][0]-gfak*R_temp[Cart::zzzz][Cart::xzz][1])+0.5/_decay*2*R_temp[Cart::zzzzz][Cart::xz][1]
R_temp[Cart::zzzzzz][Cart::zzz][0]+=pma2*R_temp[Cart::zzzzz][Cart::zzz][0]+wmp2*R_temp[Cart::zzzzz][Cart::zzz][1]+4*rzeta*(R_temp[Cart::zzzz][Cart::zzz][0]-gfak*R_temp[Cart::zzzz][Cart::zzz][1])+0.5/_decay*3*R_temp[Cart::zzzzz][Cart::zz][1]
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
