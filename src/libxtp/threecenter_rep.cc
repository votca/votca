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




using namespace std;
using namespace votca::tools;

namespace votca {
    namespace xtp {
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
        
      
        bool TCrawMatrix::FillThreeCenterRepBlock(ub::matrix<double>& _subvector,  AOShell* _shell_3, AOShell* _shell_1, AOShell* _shell_2) {
            
            const double pi = boost::math::constants::pi<double>();
            
            bool _does_contribute=true;
            

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
                cout << "switched" << endl;    

                
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
            int _ncombined =this->getBlockSize(_lmax_alpha+_lmax_beta);
            
            typedef boost::multi_array<double, 3> ma_type;
            typedef boost::multi_array_types::extent_range range;
            typedef ma_type::index index;
            ma_type::extent_gen extents;
            
            
            
            
            double _dist1=(_pos_alpha - _pos_gamma)*(_pos_alpha - _pos_gamma);
            double _dist2=(_pos_gamma - _pos_beta) * (_pos_gamma - _pos_beta);
            double _dist3=(_pos_alpha - _pos_beta) * (_pos_alpha - _pos_beta);
            
            vec amb=_pos_alpha-_pos_beta;
            double amb0=0.0;
            double amb1=0.0;
            double amb2=0.0;
         
            if (_dist3<0.03){
                amb0=amb.getX();
                amb1=amb.getY();
                amb2=amb.getZ();
            }
            cout << "before contraction" << endl;

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
            //double fak = 0.5 / (_decay);
            double gfak=_decay_gamma/_decay;
            double cfak= (_decay_alpha + _decay_beta)/_decay;         
            vec _P=(_decay_alpha*_pos_alpha+_decay_beta*_pos_beta)/(_decay_alpha+_decay_beta);
            vec _W=(_decay_alpha*_pos_alpha+_decay_beta*_pos_beta+_decay_gamma*_pos_gamma)/_decay;
            double _T = (_decay_alpha+_decay_beta)*_decay_gamma/_decay*(_P-_pos_gamma)*(_P-_pos_gamma);
            
            
            
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
            
            if ((_dist1 + _dist2 + _dist3)>0.01){
          
            pma0 = pma.getX();
            //pmb0 = pmb.getX();
            wmp0 = wmp.getX();
            wmc0 = wmc.getX();
            pma1 = pma.getY();
            //pmb1 = pmb.getY();
            wmp1 = wmp.getY();
            wmc1 = wmc.getY();
            pma2 = pma.getZ();
            //pmb2 = pmb.getZ();
            wmp2 = wmp.getZ();
            wmc2 = wmc.getZ();
            }
            
            cout << "los" << endl;
            ma_type R_temp;
            R_temp.resize(extents[ range(0, _ncombined ) ][ range(0, _ngamma ) ][ range(0, _mmax+1)]);
            //initialize to zero
            for (index i = 0; i != _ncombined; ++i) {
                for (index j = 0; j != _nbeta; ++j) {
                    for (index k = 0; k != _mmax; ++k) {

                                       R_temp[i][j][k] = 0.0;
                                   }
                               }
                           }
            
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
            
                   
            vector<double> _FmT(_mmax, 0.0); 
           
            XIntegrate(_FmT, _T);
            
            double sss = 8*pow(2*pi,0.25)*pow(_decay_alpha*_decay_beta*_decay_gamma,0.75)/((_decay_alpha+_decay_beta)*_decay_gamma);
            //ss integrals
            for (int _i=0;_i<_mmax;_i++){
                R_temp[Cart::s][Cart::s][_i]=sss*_FmT[_i];
            }
            
            cout << _mmax<< "mmax" << endl;
          
//Integral s - s - p - m0
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][0]+=wmc1*R_temp[Cart::s][Cart::s][1];
R_temp[Cart::s][Cart::x][0]+=wmc0*R_temp[Cart::s][Cart::s][1];
R_temp[Cart::s][Cart::z][0]+=wmc2*R_temp[Cart::s][Cart::s][1];
}
//------------------------------------------------------
       

//Integral s - s - p - m1
if (_mmax >1 ){
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][1]+=wmc1*R_temp[Cart::s][Cart::s][2];
R_temp[Cart::s][Cart::x][1]+=wmc0*R_temp[Cart::s][Cart::s][2];
R_temp[Cart::s][Cart::z][1]+=wmc2*R_temp[Cart::s][Cart::s][2];
}
}
//------------------------------------------------------
           
//Integral s - s - p - m2
if (_mmax >2 ){
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][2]+=wmc1*R_temp[Cart::s][Cart::s][3];
R_temp[Cart::s][Cart::x][2]+=wmc0*R_temp[Cart::s][Cart::s][3];
R_temp[Cart::s][Cart::z][2]+=wmc2*R_temp[Cart::s][Cart::s][3];
}
}
//------------------------------------------------------

//Integral s - s - p - m3
if (_mmax >3 ){
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][3]+=wmc1*R_temp[Cart::s][Cart::s][4];
R_temp[Cart::s][Cart::x][3]+=wmc0*R_temp[Cart::s][Cart::s][4];
R_temp[Cart::s][Cart::z][3]+=wmc2*R_temp[Cart::s][Cart::s][4];
}
}
//------------------------------------------------------

//Integral s - s - p - m4
if (_mmax >4 ){
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][4]+=wmc1*R_temp[Cart::s][Cart::s][5];
R_temp[Cart::s][Cart::x][4]+=wmc0*R_temp[Cart::s][Cart::s][5];
R_temp[Cart::s][Cart::z][4]+=wmc2*R_temp[Cart::s][Cart::s][5];
}
}
//------------------------------------------------------

//Integral s - s - p - m5
if (_mmax >5 ){
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][5]+=wmc1*R_temp[Cart::s][Cart::s][6];
R_temp[Cart::s][Cart::x][5]+=wmc0*R_temp[Cart::s][Cart::s][6];
R_temp[Cart::s][Cart::z][5]+=wmc2*R_temp[Cart::s][Cart::s][6];
}
}
//------------------------------------------------------

//Integral s - s - p - m6
if (_mmax >6 ){
if (_lmax_gamma>0){
R_temp[Cart::s][Cart::y][6]+=wmc1*R_temp[Cart::s][Cart::s][7];
R_temp[Cart::s][Cart::x][6]+=wmc0*R_temp[Cart::s][Cart::s][7];
R_temp[Cart::s][Cart::z][6]+=wmc2*R_temp[Cart::s][Cart::s][7];
}
}
//------------------------------------------------------

//Integral s - s - d - m0
if (_lmax_gamma>1){
R_temp[Cart::s][Cart::yy][0]+=wmc1*R_temp[Cart::s][Cart::y][1];
R_temp[Cart::s][Cart::xy][0]+=wmc0*R_temp[Cart::s][Cart::y][1];
R_temp[Cart::s][Cart::yz][0]+=wmc1*R_temp[Cart::s][Cart::z][1];
R_temp[Cart::s][Cart::xx][0]+=wmc0*R_temp[Cart::s][Cart::x][1];
R_temp[Cart::s][Cart::xz][0]+=wmc0*R_temp[Cart::s][Cart::z][1];
R_temp[Cart::s][Cart::zz][0]+=wmc2*R_temp[Cart::s][Cart::z][1];
}
//------------------------------------------------------

//Integral s - s - d - m1
if (_mmax >1 ){
if (_lmax_gamma>1){
R_temp[Cart::s][Cart::yy][1]+=wmc1*R_temp[Cart::s][Cart::y][2];
R_temp[Cart::s][Cart::xy][1]+=wmc0*R_temp[Cart::s][Cart::y][2];
R_temp[Cart::s][Cart::yz][1]+=wmc1*R_temp[Cart::s][Cart::z][2];
R_temp[Cart::s][Cart::xx][1]+=wmc0*R_temp[Cart::s][Cart::x][2];
R_temp[Cart::s][Cart::xz][1]+=wmc0*R_temp[Cart::s][Cart::z][2];
R_temp[Cart::s][Cart::zz][1]+=wmc2*R_temp[Cart::s][Cart::z][2];
}
}
//------------------------------------------------------

//Integral s - s - d - m2
if (_mmax >2 ){
if (_lmax_gamma>1){
R_temp[Cart::s][Cart::yy][2]+=wmc1*R_temp[Cart::s][Cart::y][3];
R_temp[Cart::s][Cart::xy][2]+=wmc0*R_temp[Cart::s][Cart::y][3];
R_temp[Cart::s][Cart::yz][2]+=wmc1*R_temp[Cart::s][Cart::z][3];
R_temp[Cart::s][Cart::xx][2]+=wmc0*R_temp[Cart::s][Cart::x][3];
R_temp[Cart::s][Cart::xz][2]+=wmc0*R_temp[Cart::s][Cart::z][3];
R_temp[Cart::s][Cart::zz][2]+=wmc2*R_temp[Cart::s][Cart::z][3];
}
}
//------------------------------------------------------

//Integral s - s - d - m3
if (_mmax >3 ){
if (_lmax_gamma>1){
R_temp[Cart::s][Cart::yy][3]+=wmc1*R_temp[Cart::s][Cart::y][4];
R_temp[Cart::s][Cart::xy][3]+=wmc0*R_temp[Cart::s][Cart::y][4];
R_temp[Cart::s][Cart::yz][3]+=wmc1*R_temp[Cart::s][Cart::z][4];
R_temp[Cart::s][Cart::xx][3]+=wmc0*R_temp[Cart::s][Cart::x][4];
R_temp[Cart::s][Cart::xz][3]+=wmc0*R_temp[Cart::s][Cart::z][4];
R_temp[Cart::s][Cart::zz][3]+=wmc2*R_temp[Cart::s][Cart::z][4];
}
}
//------------------------------------------------------

//Integral s - s - d - m4
if (_mmax >4 ){
if (_lmax_gamma>1){
R_temp[Cart::s][Cart::yy][4]+=wmc1*R_temp[Cart::s][Cart::y][5];
R_temp[Cart::s][Cart::xy][4]+=wmc0*R_temp[Cart::s][Cart::y][5];
R_temp[Cart::s][Cart::yz][4]+=wmc1*R_temp[Cart::s][Cart::z][5];
R_temp[Cart::s][Cart::xx][4]+=wmc0*R_temp[Cart::s][Cart::x][5];
R_temp[Cart::s][Cart::xz][4]+=wmc0*R_temp[Cart::s][Cart::z][5];
R_temp[Cart::s][Cart::zz][4]+=wmc2*R_temp[Cart::s][Cart::z][5];
}
}
//------------------------------------------------------

//Integral s - s - d - m5
if (_mmax >5 ){
if (_lmax_gamma>1){
R_temp[Cart::s][Cart::yy][5]+=wmc1*R_temp[Cart::s][Cart::y][6];
R_temp[Cart::s][Cart::xy][5]+=wmc0*R_temp[Cart::s][Cart::y][6];
R_temp[Cart::s][Cart::yz][5]+=wmc1*R_temp[Cart::s][Cart::z][6];
R_temp[Cart::s][Cart::xx][5]+=wmc0*R_temp[Cart::s][Cart::x][6];
R_temp[Cart::s][Cart::xz][5]+=wmc0*R_temp[Cart::s][Cart::z][6];
R_temp[Cart::s][Cart::zz][5]+=wmc2*R_temp[Cart::s][Cart::z][6];
}
}
//------------------------------------------------------

//Integral s - s - f - m0
if (_lmax_gamma>2){
R_temp[Cart::s][Cart::yyy][0]+=wmc1*R_temp[Cart::s][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::s][Cart::y][0]-cfak*R_temp[Cart::s][Cart::y][1]);
R_temp[Cart::s][Cart::xyy][0]+=wmc0*R_temp[Cart::s][Cart::yy][1];
R_temp[Cart::s][Cart::yyz][0]+=wmc2*R_temp[Cart::s][Cart::yy][1];
R_temp[Cart::s][Cart::xxy][0]+=wmc1*R_temp[Cart::s][Cart::xx][1];
R_temp[Cart::s][Cart::xyz][0]+=wmc0*R_temp[Cart::s][Cart::yz][1];
R_temp[Cart::s][Cart::yzz][0]+=wmc1*R_temp[Cart::s][Cart::zz][1];
R_temp[Cart::s][Cart::xxx][0]+=wmc0*R_temp[Cart::s][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::s][Cart::x][0]-cfak*R_temp[Cart::s][Cart::x][1]);
R_temp[Cart::s][Cart::xxz][0]+=wmc2*R_temp[Cart::s][Cart::xx][1];
R_temp[Cart::s][Cart::xzz][0]+=wmc0*R_temp[Cart::s][Cart::zz][1];
R_temp[Cart::s][Cart::zzz][0]+=wmc2*R_temp[Cart::s][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::s][Cart::z][0]-cfak*R_temp[Cart::s][Cart::z][1]);
}
//------------------------------------------------------

//Integral s - s - f - m1
if (_mmax >1 ){
if (_lmax_gamma>2){
R_temp[Cart::s][Cart::yyy][1]+=wmc1*R_temp[Cart::s][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::s][Cart::y][1]-cfak*R_temp[Cart::s][Cart::y][2]);
R_temp[Cart::s][Cart::xyy][1]+=wmc0*R_temp[Cart::s][Cart::yy][2];
R_temp[Cart::s][Cart::yyz][1]+=wmc2*R_temp[Cart::s][Cart::yy][2];
R_temp[Cart::s][Cart::xxy][1]+=wmc1*R_temp[Cart::s][Cart::xx][2];
R_temp[Cart::s][Cart::xyz][1]+=wmc0*R_temp[Cart::s][Cart::yz][2];
R_temp[Cart::s][Cart::yzz][1]+=wmc1*R_temp[Cart::s][Cart::zz][2];
R_temp[Cart::s][Cart::xxx][1]+=wmc0*R_temp[Cart::s][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::s][Cart::x][1]-cfak*R_temp[Cart::s][Cart::x][2]);
R_temp[Cart::s][Cart::xxz][1]+=wmc2*R_temp[Cart::s][Cart::xx][2];
R_temp[Cart::s][Cart::xzz][1]+=wmc0*R_temp[Cart::s][Cart::zz][2];
R_temp[Cart::s][Cart::zzz][1]+=wmc2*R_temp[Cart::s][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::s][Cart::z][1]-cfak*R_temp[Cart::s][Cart::z][2]);
}
}
//------------------------------------------------------

//Integral s - s - f - m2
if (_mmax >2 ){
if (_lmax_gamma>2){
R_temp[Cart::s][Cart::yyy][2]+=wmc1*R_temp[Cart::s][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::s][Cart::y][2]-cfak*R_temp[Cart::s][Cart::y][3]);
R_temp[Cart::s][Cart::xyy][2]+=wmc0*R_temp[Cart::s][Cart::yy][3];
R_temp[Cart::s][Cart::yyz][2]+=wmc2*R_temp[Cart::s][Cart::yy][3];
R_temp[Cart::s][Cart::xxy][2]+=wmc1*R_temp[Cart::s][Cart::xx][3];
R_temp[Cart::s][Cart::xyz][2]+=wmc0*R_temp[Cart::s][Cart::yz][3];
R_temp[Cart::s][Cart::yzz][2]+=wmc1*R_temp[Cart::s][Cart::zz][3];
R_temp[Cart::s][Cart::xxx][2]+=wmc0*R_temp[Cart::s][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::s][Cart::x][2]-cfak*R_temp[Cart::s][Cart::x][3]);
R_temp[Cart::s][Cart::xxz][2]+=wmc2*R_temp[Cart::s][Cart::xx][3];
R_temp[Cart::s][Cart::xzz][2]+=wmc0*R_temp[Cart::s][Cart::zz][3];
R_temp[Cart::s][Cart::zzz][2]+=wmc2*R_temp[Cart::s][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::s][Cart::z][2]-cfak*R_temp[Cart::s][Cart::z][3]);
}
}
//------------------------------------------------------

//Integral s - s - f - m3
if (_mmax >3 ){
if (_lmax_gamma>2){
R_temp[Cart::s][Cart::yyy][3]+=wmc1*R_temp[Cart::s][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::s][Cart::y][3]-cfak*R_temp[Cart::s][Cart::y][4]);
R_temp[Cart::s][Cart::xyy][3]+=wmc0*R_temp[Cart::s][Cart::yy][4];
R_temp[Cart::s][Cart::yyz][3]+=wmc2*R_temp[Cart::s][Cart::yy][4];
R_temp[Cart::s][Cart::xxy][3]+=wmc1*R_temp[Cart::s][Cart::xx][4];
R_temp[Cart::s][Cart::xyz][3]+=wmc0*R_temp[Cart::s][Cart::yz][4];
R_temp[Cart::s][Cart::yzz][3]+=wmc1*R_temp[Cart::s][Cart::zz][4];
R_temp[Cart::s][Cart::xxx][3]+=wmc0*R_temp[Cart::s][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::s][Cart::x][3]-cfak*R_temp[Cart::s][Cart::x][4]);
R_temp[Cart::s][Cart::xxz][3]+=wmc2*R_temp[Cart::s][Cart::xx][4];
R_temp[Cart::s][Cart::xzz][3]+=wmc0*R_temp[Cart::s][Cart::zz][4];
R_temp[Cart::s][Cart::zzz][3]+=wmc2*R_temp[Cart::s][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::s][Cart::z][3]-cfak*R_temp[Cart::s][Cart::z][4]);
}
}
//------------------------------------------------------

//Integral s - s - f - m4
if (_mmax >4 ){
if (_lmax_gamma>2){
R_temp[Cart::s][Cart::yyy][4]+=wmc1*R_temp[Cart::s][Cart::yy][5]+1/_decay_gamma*(R_temp[Cart::s][Cart::y][4]-cfak*R_temp[Cart::s][Cart::y][5]);
R_temp[Cart::s][Cart::xyy][4]+=wmc0*R_temp[Cart::s][Cart::yy][5];
R_temp[Cart::s][Cart::yyz][4]+=wmc2*R_temp[Cart::s][Cart::yy][5];
R_temp[Cart::s][Cart::xxy][4]+=wmc1*R_temp[Cart::s][Cart::xx][5];
R_temp[Cart::s][Cart::xyz][4]+=wmc0*R_temp[Cart::s][Cart::yz][5];
R_temp[Cart::s][Cart::yzz][4]+=wmc1*R_temp[Cart::s][Cart::zz][5];
R_temp[Cart::s][Cart::xxx][4]+=wmc0*R_temp[Cart::s][Cart::xx][5]+1/_decay_gamma*(R_temp[Cart::s][Cart::x][4]-cfak*R_temp[Cart::s][Cart::x][5]);
R_temp[Cart::s][Cart::xxz][4]+=wmc2*R_temp[Cart::s][Cart::xx][5];
R_temp[Cart::s][Cart::xzz][4]+=wmc0*R_temp[Cart::s][Cart::zz][5];
R_temp[Cart::s][Cart::zzz][4]+=wmc2*R_temp[Cart::s][Cart::zz][5]+1/_decay_gamma*(R_temp[Cart::s][Cart::z][4]-cfak*R_temp[Cart::s][Cart::z][5]);
}
}
//------------------------------------------------------

//Integral p - s - s - m0
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][0]+=pma1*R_temp[Cart::s][Cart::s][0]+wmp1*R_temp[Cart::s][Cart::s][1];
R_temp[Cart::x][Cart::s][0]+=pma0*R_temp[Cart::s][Cart::s][0]+wmp0*R_temp[Cart::s][Cart::s][1];
R_temp[Cart::z][Cart::s][0]+=pma2*R_temp[Cart::s][Cart::s][0]+wmp2*R_temp[Cart::s][Cart::s][1];
}
//------------------------------------------------------

//Integral p - s - s - m1
if (_mmax >1 ){
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][1]+=pma1*R_temp[Cart::s][Cart::s][1]+wmp1*R_temp[Cart::s][Cart::s][2];
R_temp[Cart::x][Cart::s][1]+=pma0*R_temp[Cart::s][Cart::s][1]+wmp0*R_temp[Cart::s][Cart::s][2];
R_temp[Cart::z][Cart::s][1]+=pma2*R_temp[Cart::s][Cart::s][1]+wmp2*R_temp[Cart::s][Cart::s][2];
}
}
//------------------------------------------------------

//Integral p - s - s - m2
if (_mmax >2 ){
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][2]+=pma1*R_temp[Cart::s][Cart::s][2]+wmp1*R_temp[Cart::s][Cart::s][3];
R_temp[Cart::x][Cart::s][2]+=pma0*R_temp[Cart::s][Cart::s][2]+wmp0*R_temp[Cart::s][Cart::s][3];
R_temp[Cart::z][Cart::s][2]+=pma2*R_temp[Cart::s][Cart::s][2]+wmp2*R_temp[Cart::s][Cart::s][3];
}
}
//------------------------------------------------------

//Integral p - s - s - m3
if (_mmax >3 ){
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][3]+=pma1*R_temp[Cart::s][Cart::s][3]+wmp1*R_temp[Cart::s][Cart::s][4];
R_temp[Cart::x][Cart::s][3]+=pma0*R_temp[Cart::s][Cart::s][3]+wmp0*R_temp[Cart::s][Cart::s][4];
R_temp[Cart::z][Cart::s][3]+=pma2*R_temp[Cart::s][Cart::s][3]+wmp2*R_temp[Cart::s][Cart::s][4];
}
}
//------------------------------------------------------

//Integral p - s - s - m4
if (_mmax >4 ){
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][4]+=pma1*R_temp[Cart::s][Cart::s][4]+wmp1*R_temp[Cart::s][Cart::s][5];
R_temp[Cart::x][Cart::s][4]+=pma0*R_temp[Cart::s][Cart::s][4]+wmp0*R_temp[Cart::s][Cart::s][5];
R_temp[Cart::z][Cart::s][4]+=pma2*R_temp[Cart::s][Cart::s][4]+wmp2*R_temp[Cart::s][Cart::s][5];
}
}
//------------------------------------------------------

//Integral p - s - s - m5
if (_mmax >5 ){
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][5]+=pma1*R_temp[Cart::s][Cart::s][5]+wmp1*R_temp[Cart::s][Cart::s][6];
R_temp[Cart::x][Cart::s][5]+=pma0*R_temp[Cart::s][Cart::s][5]+wmp0*R_temp[Cart::s][Cart::s][6];
R_temp[Cart::z][Cart::s][5]+=pma2*R_temp[Cart::s][Cart::s][5]+wmp2*R_temp[Cart::s][Cart::s][6];
}
}
//------------------------------------------------------

//Integral p - s - s - m6
if (_mmax >6 ){
if (_lmax_alpha>0){
R_temp[Cart::y][Cart::s][6]+=pma1*R_temp[Cart::s][Cart::s][6]+wmp1*R_temp[Cart::s][Cart::s][7];
R_temp[Cart::x][Cart::s][6]+=pma0*R_temp[Cart::s][Cart::s][6]+wmp0*R_temp[Cart::s][Cart::s][7];
R_temp[Cart::z][Cart::s][6]+=pma2*R_temp[Cart::s][Cart::s][6]+wmp2*R_temp[Cart::s][Cart::s][7];
}
}
//------------------------------------------------------

//Integral p - s - p - m0
if (_lmax_alpha>0 && _lmax_gamma>0){
R_temp[Cart::y][Cart::y][0]+=pma1*R_temp[Cart::s][Cart::y][0]+wmp1*R_temp[Cart::s][Cart::y][1]+0.5/_decay*1*R_temp[Cart::s][Cart::s][1];
R_temp[Cart::y][Cart::x][0]+=pma1*R_temp[Cart::s][Cart::x][0]+wmp1*R_temp[Cart::s][Cart::x][1];
R_temp[Cart::y][Cart::z][0]+=pma1*R_temp[Cart::s][Cart::z][0]+wmp1*R_temp[Cart::s][Cart::z][1];
R_temp[Cart::x][Cart::y][0]+=pma0*R_temp[Cart::s][Cart::y][0]+wmp0*R_temp[Cart::s][Cart::y][1];
R_temp[Cart::x][Cart::x][0]+=pma0*R_temp[Cart::s][Cart::x][0]+wmp0*R_temp[Cart::s][Cart::x][1]+0.5/_decay*1*R_temp[Cart::s][Cart::s][1];
R_temp[Cart::x][Cart::z][0]+=pma0*R_temp[Cart::s][Cart::z][0]+wmp0*R_temp[Cart::s][Cart::z][1];
R_temp[Cart::z][Cart::y][0]+=pma2*R_temp[Cart::s][Cart::y][0]+wmp2*R_temp[Cart::s][Cart::y][1];
R_temp[Cart::z][Cart::x][0]+=pma2*R_temp[Cart::s][Cart::x][0]+wmp2*R_temp[Cart::s][Cart::x][1];
R_temp[Cart::z][Cart::z][0]+=pma2*R_temp[Cart::s][Cart::z][0]+wmp2*R_temp[Cart::s][Cart::z][1]+0.5/_decay*1*R_temp[Cart::s][Cart::s][1];
}
//------------------------------------------------------

//Integral p - s - p - m1
if (_mmax >1 ){
if (_lmax_alpha>0 && _lmax_gamma>0){
R_temp[Cart::y][Cart::y][1]+=pma1*R_temp[Cart::s][Cart::y][1]+wmp1*R_temp[Cart::s][Cart::y][2]+0.5/_decay*1*R_temp[Cart::s][Cart::s][2];
R_temp[Cart::y][Cart::x][1]+=pma1*R_temp[Cart::s][Cart::x][1]+wmp1*R_temp[Cart::s][Cart::x][2];
R_temp[Cart::y][Cart::z][1]+=pma1*R_temp[Cart::s][Cart::z][1]+wmp1*R_temp[Cart::s][Cart::z][2];
R_temp[Cart::x][Cart::y][1]+=pma0*R_temp[Cart::s][Cart::y][1]+wmp0*R_temp[Cart::s][Cart::y][2];
R_temp[Cart::x][Cart::x][1]+=pma0*R_temp[Cart::s][Cart::x][1]+wmp0*R_temp[Cart::s][Cart::x][2]+0.5/_decay*1*R_temp[Cart::s][Cart::s][2];
R_temp[Cart::x][Cart::z][1]+=pma0*R_temp[Cart::s][Cart::z][1]+wmp0*R_temp[Cart::s][Cart::z][2];
R_temp[Cart::z][Cart::y][1]+=pma2*R_temp[Cart::s][Cart::y][1]+wmp2*R_temp[Cart::s][Cart::y][2];
R_temp[Cart::z][Cart::x][1]+=pma2*R_temp[Cart::s][Cart::x][1]+wmp2*R_temp[Cart::s][Cart::x][2];
R_temp[Cart::z][Cart::z][1]+=pma2*R_temp[Cart::s][Cart::z][1]+wmp2*R_temp[Cart::s][Cart::z][2]+0.5/_decay*1*R_temp[Cart::s][Cart::s][2];
}
}
//------------------------------------------------------

//Integral p - s - p - m2
if (_mmax >2 ){
if (_lmax_alpha>0 && _lmax_gamma>0){
R_temp[Cart::y][Cart::y][2]+=pma1*R_temp[Cart::s][Cart::y][2]+wmp1*R_temp[Cart::s][Cart::y][3]+0.5/_decay*1*R_temp[Cart::s][Cart::s][3];
R_temp[Cart::y][Cart::x][2]+=pma1*R_temp[Cart::s][Cart::x][2]+wmp1*R_temp[Cart::s][Cart::x][3];
R_temp[Cart::y][Cart::z][2]+=pma1*R_temp[Cart::s][Cart::z][2]+wmp1*R_temp[Cart::s][Cart::z][3];
R_temp[Cart::x][Cart::y][2]+=pma0*R_temp[Cart::s][Cart::y][2]+wmp0*R_temp[Cart::s][Cart::y][3];
R_temp[Cart::x][Cart::x][2]+=pma0*R_temp[Cart::s][Cart::x][2]+wmp0*R_temp[Cart::s][Cart::x][3]+0.5/_decay*1*R_temp[Cart::s][Cart::s][3];
R_temp[Cart::x][Cart::z][2]+=pma0*R_temp[Cart::s][Cart::z][2]+wmp0*R_temp[Cart::s][Cart::z][3];
R_temp[Cart::z][Cart::y][2]+=pma2*R_temp[Cart::s][Cart::y][2]+wmp2*R_temp[Cart::s][Cart::y][3];
R_temp[Cart::z][Cart::x][2]+=pma2*R_temp[Cart::s][Cart::x][2]+wmp2*R_temp[Cart::s][Cart::x][3];
R_temp[Cart::z][Cart::z][2]+=pma2*R_temp[Cart::s][Cart::z][2]+wmp2*R_temp[Cart::s][Cart::z][3]+0.5/_decay*1*R_temp[Cart::s][Cart::s][3];
}
}
//------------------------------------------------------

//Integral p - s - p - m3
if (_mmax >3 ){
if (_lmax_alpha>0 && _lmax_gamma>0){
R_temp[Cart::y][Cart::y][3]+=pma1*R_temp[Cart::s][Cart::y][3]+wmp1*R_temp[Cart::s][Cart::y][4]+0.5/_decay*1*R_temp[Cart::s][Cart::s][4];
R_temp[Cart::y][Cart::x][3]+=pma1*R_temp[Cart::s][Cart::x][3]+wmp1*R_temp[Cart::s][Cart::x][4];
R_temp[Cart::y][Cart::z][3]+=pma1*R_temp[Cart::s][Cart::z][3]+wmp1*R_temp[Cart::s][Cart::z][4];
R_temp[Cart::x][Cart::y][3]+=pma0*R_temp[Cart::s][Cart::y][3]+wmp0*R_temp[Cart::s][Cart::y][4];
R_temp[Cart::x][Cart::x][3]+=pma0*R_temp[Cart::s][Cart::x][3]+wmp0*R_temp[Cart::s][Cart::x][4]+0.5/_decay*1*R_temp[Cart::s][Cart::s][4];
R_temp[Cart::x][Cart::z][3]+=pma0*R_temp[Cart::s][Cart::z][3]+wmp0*R_temp[Cart::s][Cart::z][4];
R_temp[Cart::z][Cart::y][3]+=pma2*R_temp[Cart::s][Cart::y][3]+wmp2*R_temp[Cart::s][Cart::y][4];
R_temp[Cart::z][Cart::x][3]+=pma2*R_temp[Cart::s][Cart::x][3]+wmp2*R_temp[Cart::s][Cart::x][4];
R_temp[Cart::z][Cart::z][3]+=pma2*R_temp[Cart::s][Cart::z][3]+wmp2*R_temp[Cart::s][Cart::z][4]+0.5/_decay*1*R_temp[Cart::s][Cart::s][4];
}
}
//------------------------------------------------------

//Integral p - s - p - m4
if (_mmax >4 ){
if (_lmax_alpha>0 && _lmax_gamma>0){
R_temp[Cart::y][Cart::y][4]+=pma1*R_temp[Cart::s][Cart::y][4]+wmp1*R_temp[Cart::s][Cart::y][5]+0.5/_decay*1*R_temp[Cart::s][Cart::s][5];
R_temp[Cart::y][Cart::x][4]+=pma1*R_temp[Cart::s][Cart::x][4]+wmp1*R_temp[Cart::s][Cart::x][5];
R_temp[Cart::y][Cart::z][4]+=pma1*R_temp[Cart::s][Cart::z][4]+wmp1*R_temp[Cart::s][Cart::z][5];
R_temp[Cart::x][Cart::y][4]+=pma0*R_temp[Cart::s][Cart::y][4]+wmp0*R_temp[Cart::s][Cart::y][5];
R_temp[Cart::x][Cart::x][4]+=pma0*R_temp[Cart::s][Cart::x][4]+wmp0*R_temp[Cart::s][Cart::x][5]+0.5/_decay*1*R_temp[Cart::s][Cart::s][5];
R_temp[Cart::x][Cart::z][4]+=pma0*R_temp[Cart::s][Cart::z][4]+wmp0*R_temp[Cart::s][Cart::z][5];
R_temp[Cart::z][Cart::y][4]+=pma2*R_temp[Cart::s][Cart::y][4]+wmp2*R_temp[Cart::s][Cart::y][5];
R_temp[Cart::z][Cart::x][4]+=pma2*R_temp[Cart::s][Cart::x][4]+wmp2*R_temp[Cart::s][Cart::x][5];
R_temp[Cart::z][Cart::z][4]+=pma2*R_temp[Cart::s][Cart::z][4]+wmp2*R_temp[Cart::s][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::s][5];
}
}
//------------------------------------------------------

//Integral p - s - p - m5
if (_mmax >5 ){
if (_lmax_alpha>0 && _lmax_gamma>0){
R_temp[Cart::y][Cart::y][5]+=pma1*R_temp[Cart::s][Cart::y][5]+wmp1*R_temp[Cart::s][Cart::y][6]+0.5/_decay*1*R_temp[Cart::s][Cart::s][6];
R_temp[Cart::y][Cart::x][5]+=pma1*R_temp[Cart::s][Cart::x][5]+wmp1*R_temp[Cart::s][Cart::x][6];
R_temp[Cart::y][Cart::z][5]+=pma1*R_temp[Cart::s][Cart::z][5]+wmp1*R_temp[Cart::s][Cart::z][6];
R_temp[Cart::x][Cart::y][5]+=pma0*R_temp[Cart::s][Cart::y][5]+wmp0*R_temp[Cart::s][Cart::y][6];
R_temp[Cart::x][Cart::x][5]+=pma0*R_temp[Cart::s][Cart::x][5]+wmp0*R_temp[Cart::s][Cart::x][6]+0.5/_decay*1*R_temp[Cart::s][Cart::s][6];
R_temp[Cart::x][Cart::z][5]+=pma0*R_temp[Cart::s][Cart::z][5]+wmp0*R_temp[Cart::s][Cart::z][6];
R_temp[Cart::z][Cart::y][5]+=pma2*R_temp[Cart::s][Cart::y][5]+wmp2*R_temp[Cart::s][Cart::y][6];
R_temp[Cart::z][Cart::x][5]+=pma2*R_temp[Cart::s][Cart::x][5]+wmp2*R_temp[Cart::s][Cart::x][6];
R_temp[Cart::z][Cart::z][5]+=pma2*R_temp[Cart::s][Cart::z][5]+wmp2*R_temp[Cart::s][Cart::z][6]+0.5/_decay*1*R_temp[Cart::s][Cart::s][6];
}
}
//------------------------------------------------------

//Integral p - s - d - m0
if (_lmax_alpha>0 && _lmax_gamma>1){
R_temp[Cart::y][Cart::yy][0]+=wmc1*R_temp[Cart::y][Cart::y][1]+0.5/_decay*1*R_temp[Cart::s][Cart::y][1];
R_temp[Cart::y][Cart::xy][0]+=wmc0*R_temp[Cart::y][Cart::y][1];
R_temp[Cart::y][Cart::yz][0]+=wmc1*R_temp[Cart::y][Cart::z][1]+0.5/_decay*1*R_temp[Cart::s][Cart::z][1];
R_temp[Cart::y][Cart::xx][0]+=wmc0*R_temp[Cart::y][Cart::x][1];
R_temp[Cart::y][Cart::xz][0]+=wmc0*R_temp[Cart::y][Cart::z][1];
R_temp[Cart::y][Cart::zz][0]+=wmc2*R_temp[Cart::y][Cart::z][1];
R_temp[Cart::x][Cart::yy][0]+=wmc1*R_temp[Cart::x][Cart::y][1];
R_temp[Cart::x][Cart::xy][0]+=wmc0*R_temp[Cart::x][Cart::y][1]+0.5/_decay*1*R_temp[Cart::s][Cart::y][1];
R_temp[Cart::x][Cart::yz][0]+=wmc1*R_temp[Cart::x][Cart::z][1];
R_temp[Cart::x][Cart::xx][0]+=wmc0*R_temp[Cart::x][Cart::x][1]+0.5/_decay*1*R_temp[Cart::s][Cart::x][1];
R_temp[Cart::x][Cart::xz][0]+=wmc0*R_temp[Cart::x][Cart::z][1]+0.5/_decay*1*R_temp[Cart::s][Cart::z][1];
R_temp[Cart::x][Cart::zz][0]+=wmc2*R_temp[Cart::x][Cart::z][1];
R_temp[Cart::z][Cart::yy][0]+=wmc1*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::z][Cart::xy][0]+=wmc0*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::z][Cart::yz][0]+=wmc1*R_temp[Cart::z][Cart::z][1];
R_temp[Cart::z][Cart::xx][0]+=wmc0*R_temp[Cart::z][Cart::x][1];
R_temp[Cart::z][Cart::xz][0]+=wmc0*R_temp[Cart::z][Cart::z][1];
R_temp[Cart::z][Cart::zz][0]+=wmc2*R_temp[Cart::z][Cart::z][1]+0.5/_decay*1*R_temp[Cart::s][Cart::z][1];
}
//------------------------------------------------------

//Integral p - s - d - m1
if (_mmax >1 ){
if (_lmax_alpha>0 && _lmax_gamma>1){
R_temp[Cart::y][Cart::yy][1]+=wmc1*R_temp[Cart::y][Cart::y][2]+0.5/_decay*1*R_temp[Cart::s][Cart::y][2];
R_temp[Cart::y][Cart::xy][1]+=wmc0*R_temp[Cart::y][Cart::y][2];
R_temp[Cart::y][Cart::yz][1]+=wmc1*R_temp[Cart::y][Cart::z][2]+0.5/_decay*1*R_temp[Cart::s][Cart::z][2];
R_temp[Cart::y][Cart::xx][1]+=wmc0*R_temp[Cart::y][Cart::x][2];
R_temp[Cart::y][Cart::xz][1]+=wmc0*R_temp[Cart::y][Cart::z][2];
R_temp[Cart::y][Cart::zz][1]+=wmc2*R_temp[Cart::y][Cart::z][2];
R_temp[Cart::x][Cart::yy][1]+=wmc1*R_temp[Cart::x][Cart::y][2];
R_temp[Cart::x][Cart::xy][1]+=wmc0*R_temp[Cart::x][Cart::y][2]+0.5/_decay*1*R_temp[Cart::s][Cart::y][2];
R_temp[Cart::x][Cart::yz][1]+=wmc1*R_temp[Cart::x][Cart::z][2];
R_temp[Cart::x][Cart::xx][1]+=wmc0*R_temp[Cart::x][Cart::x][2]+0.5/_decay*1*R_temp[Cart::s][Cart::x][2];
R_temp[Cart::x][Cart::xz][1]+=wmc0*R_temp[Cart::x][Cart::z][2]+0.5/_decay*1*R_temp[Cart::s][Cart::z][2];
R_temp[Cart::x][Cart::zz][1]+=wmc2*R_temp[Cart::x][Cart::z][2];
R_temp[Cart::z][Cart::yy][1]+=wmc1*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::z][Cart::xy][1]+=wmc0*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::z][Cart::yz][1]+=wmc1*R_temp[Cart::z][Cart::z][2];
R_temp[Cart::z][Cart::xx][1]+=wmc0*R_temp[Cart::z][Cart::x][2];
R_temp[Cart::z][Cart::xz][1]+=wmc0*R_temp[Cart::z][Cart::z][2];
R_temp[Cart::z][Cart::zz][1]+=wmc2*R_temp[Cart::z][Cart::z][2]+0.5/_decay*1*R_temp[Cart::s][Cart::z][2];
}
}
//------------------------------------------------------

//Integral p - s - d - m2
if (_mmax >2 ){
if (_lmax_alpha>0 && _lmax_gamma>1){
R_temp[Cart::y][Cart::yy][2]+=wmc1*R_temp[Cart::y][Cart::y][3]+0.5/_decay*1*R_temp[Cart::s][Cart::y][3];
R_temp[Cart::y][Cart::xy][2]+=wmc0*R_temp[Cart::y][Cart::y][3];
R_temp[Cart::y][Cart::yz][2]+=wmc1*R_temp[Cart::y][Cart::z][3]+0.5/_decay*1*R_temp[Cart::s][Cart::z][3];
R_temp[Cart::y][Cart::xx][2]+=wmc0*R_temp[Cart::y][Cart::x][3];
R_temp[Cart::y][Cart::xz][2]+=wmc0*R_temp[Cart::y][Cart::z][3];
R_temp[Cart::y][Cart::zz][2]+=wmc2*R_temp[Cart::y][Cart::z][3];
R_temp[Cart::x][Cart::yy][2]+=wmc1*R_temp[Cart::x][Cart::y][3];
R_temp[Cart::x][Cart::xy][2]+=wmc0*R_temp[Cart::x][Cart::y][3]+0.5/_decay*1*R_temp[Cart::s][Cart::y][3];
R_temp[Cart::x][Cart::yz][2]+=wmc1*R_temp[Cart::x][Cart::z][3];
R_temp[Cart::x][Cart::xx][2]+=wmc0*R_temp[Cart::x][Cart::x][3]+0.5/_decay*1*R_temp[Cart::s][Cart::x][3];
R_temp[Cart::x][Cart::xz][2]+=wmc0*R_temp[Cart::x][Cart::z][3]+0.5/_decay*1*R_temp[Cart::s][Cart::z][3];
R_temp[Cart::x][Cart::zz][2]+=wmc2*R_temp[Cart::x][Cart::z][3];
R_temp[Cart::z][Cart::yy][2]+=wmc1*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::z][Cart::xy][2]+=wmc0*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::z][Cart::yz][2]+=wmc1*R_temp[Cart::z][Cart::z][3];
R_temp[Cart::z][Cart::xx][2]+=wmc0*R_temp[Cart::z][Cart::x][3];
R_temp[Cart::z][Cart::xz][2]+=wmc0*R_temp[Cart::z][Cart::z][3];
R_temp[Cart::z][Cart::zz][2]+=wmc2*R_temp[Cart::z][Cart::z][3]+0.5/_decay*1*R_temp[Cart::s][Cart::z][3];
}
}
//------------------------------------------------------

//Integral p - s - d - m3
if (_mmax >3 ){
if (_lmax_alpha>0 && _lmax_gamma>1){
R_temp[Cart::y][Cart::yy][3]+=wmc1*R_temp[Cart::y][Cart::y][4]+0.5/_decay*1*R_temp[Cart::s][Cart::y][4];
R_temp[Cart::y][Cart::xy][3]+=wmc0*R_temp[Cart::y][Cart::y][4];
R_temp[Cart::y][Cart::yz][3]+=wmc1*R_temp[Cart::y][Cart::z][4]+0.5/_decay*1*R_temp[Cart::s][Cart::z][4];
R_temp[Cart::y][Cart::xx][3]+=wmc0*R_temp[Cart::y][Cart::x][4];
R_temp[Cart::y][Cart::xz][3]+=wmc0*R_temp[Cart::y][Cart::z][4];
R_temp[Cart::y][Cart::zz][3]+=wmc2*R_temp[Cart::y][Cart::z][4];
R_temp[Cart::x][Cart::yy][3]+=wmc1*R_temp[Cart::x][Cart::y][4];
R_temp[Cart::x][Cart::xy][3]+=wmc0*R_temp[Cart::x][Cart::y][4]+0.5/_decay*1*R_temp[Cart::s][Cart::y][4];
R_temp[Cart::x][Cart::yz][3]+=wmc1*R_temp[Cart::x][Cart::z][4];
R_temp[Cart::x][Cart::xx][3]+=wmc0*R_temp[Cart::x][Cart::x][4]+0.5/_decay*1*R_temp[Cart::s][Cart::x][4];
R_temp[Cart::x][Cart::xz][3]+=wmc0*R_temp[Cart::x][Cart::z][4]+0.5/_decay*1*R_temp[Cart::s][Cart::z][4];
R_temp[Cart::x][Cart::zz][3]+=wmc2*R_temp[Cart::x][Cart::z][4];
R_temp[Cart::z][Cart::yy][3]+=wmc1*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::z][Cart::xy][3]+=wmc0*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::z][Cart::yz][3]+=wmc1*R_temp[Cart::z][Cart::z][4];
R_temp[Cart::z][Cart::xx][3]+=wmc0*R_temp[Cart::z][Cart::x][4];
R_temp[Cart::z][Cart::xz][3]+=wmc0*R_temp[Cart::z][Cart::z][4];
R_temp[Cart::z][Cart::zz][3]+=wmc2*R_temp[Cart::z][Cart::z][4]+0.5/_decay*1*R_temp[Cart::s][Cart::z][4];
}
}
//------------------------------------------------------

//Integral p - s - d - m4
if (_mmax >4 ){
if (_lmax_alpha>0 && _lmax_gamma>1){
R_temp[Cart::y][Cart::yy][4]+=wmc1*R_temp[Cart::y][Cart::y][5]+0.5/_decay*1*R_temp[Cart::s][Cart::y][5];
R_temp[Cart::y][Cart::xy][4]+=wmc0*R_temp[Cart::y][Cart::y][5];
R_temp[Cart::y][Cart::yz][4]+=wmc1*R_temp[Cart::y][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::z][5];
R_temp[Cart::y][Cart::xx][4]+=wmc0*R_temp[Cart::y][Cart::x][5];
R_temp[Cart::y][Cart::xz][4]+=wmc0*R_temp[Cart::y][Cart::z][5];
R_temp[Cart::y][Cart::zz][4]+=wmc2*R_temp[Cart::y][Cart::z][5];
R_temp[Cart::x][Cart::yy][4]+=wmc1*R_temp[Cart::x][Cart::y][5];
R_temp[Cart::x][Cart::xy][4]+=wmc0*R_temp[Cart::x][Cart::y][5]+0.5/_decay*1*R_temp[Cart::s][Cart::y][5];
R_temp[Cart::x][Cart::yz][4]+=wmc1*R_temp[Cart::x][Cart::z][5];
R_temp[Cart::x][Cart::xx][4]+=wmc0*R_temp[Cart::x][Cart::x][5]+0.5/_decay*1*R_temp[Cart::s][Cart::x][5];
R_temp[Cart::x][Cart::xz][4]+=wmc0*R_temp[Cart::x][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::z][5];
R_temp[Cart::x][Cart::zz][4]+=wmc2*R_temp[Cart::x][Cart::z][5];
R_temp[Cart::z][Cart::yy][4]+=wmc1*R_temp[Cart::z][Cart::y][5];
R_temp[Cart::z][Cart::xy][4]+=wmc0*R_temp[Cart::z][Cart::y][5];
R_temp[Cart::z][Cart::yz][4]+=wmc1*R_temp[Cart::z][Cart::z][5];
R_temp[Cart::z][Cart::xx][4]+=wmc0*R_temp[Cart::z][Cart::x][5];
R_temp[Cart::z][Cart::xz][4]+=wmc0*R_temp[Cart::z][Cart::z][5];
R_temp[Cart::z][Cart::zz][4]+=wmc2*R_temp[Cart::z][Cart::z][5]+0.5/_decay*1*R_temp[Cart::s][Cart::z][5];
}
}
//------------------------------------------------------

//Integral p - s - f - m0
if (_lmax_alpha>0 && _lmax_gamma>2){
R_temp[Cart::y][Cart::yyy][0]+=wmc1*R_temp[Cart::y][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::y][Cart::y][0]-cfak*R_temp[Cart::y][Cart::y][1])+0.5/_decay*1*R_temp[Cart::s][Cart::yy][1];
R_temp[Cart::y][Cart::xyy][0]+=wmc0*R_temp[Cart::y][Cart::yy][1];
R_temp[Cart::y][Cart::yyz][0]+=wmc2*R_temp[Cart::y][Cart::yy][1];
R_temp[Cart::y][Cart::xxy][0]+=wmc1*R_temp[Cart::y][Cart::xx][1]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][1];
R_temp[Cart::y][Cart::xyz][0]+=wmc0*R_temp[Cart::y][Cart::yz][1];
R_temp[Cart::y][Cart::yzz][0]+=wmc1*R_temp[Cart::y][Cart::zz][1]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][1];
R_temp[Cart::y][Cart::xxx][0]+=wmc0*R_temp[Cart::y][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::y][Cart::x][0]-cfak*R_temp[Cart::y][Cart::x][1]);
R_temp[Cart::y][Cart::xxz][0]+=wmc2*R_temp[Cart::y][Cart::xx][1];
R_temp[Cart::y][Cart::xzz][0]+=wmc0*R_temp[Cart::y][Cart::zz][1];
R_temp[Cart::y][Cart::zzz][0]+=wmc2*R_temp[Cart::y][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::y][Cart::z][0]-cfak*R_temp[Cart::y][Cart::z][1]);
R_temp[Cart::x][Cart::yyy][0]+=wmc1*R_temp[Cart::x][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::x][Cart::y][0]-cfak*R_temp[Cart::x][Cart::y][1]);
R_temp[Cart::x][Cart::xyy][0]+=wmc0*R_temp[Cart::x][Cart::yy][1]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][1];
R_temp[Cart::x][Cart::yyz][0]+=wmc2*R_temp[Cart::x][Cart::yy][1];
R_temp[Cart::x][Cart::xxy][0]+=wmc1*R_temp[Cart::x][Cart::xx][1];
R_temp[Cart::x][Cart::xyz][0]+=wmc0*R_temp[Cart::x][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::s][Cart::yz][1];
R_temp[Cart::x][Cart::yzz][0]+=wmc1*R_temp[Cart::x][Cart::zz][1];
R_temp[Cart::x][Cart::xxx][0]+=wmc0*R_temp[Cart::x][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::x][Cart::x][0]-cfak*R_temp[Cart::x][Cart::x][1])+0.5/_decay*1*R_temp[Cart::s][Cart::xx][1];
R_temp[Cart::x][Cart::xxz][0]+=wmc2*R_temp[Cart::x][Cart::xx][1];
R_temp[Cart::x][Cart::xzz][0]+=wmc0*R_temp[Cart::x][Cart::zz][1]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][1];
R_temp[Cart::x][Cart::zzz][0]+=wmc2*R_temp[Cart::x][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::x][Cart::z][0]-cfak*R_temp[Cart::x][Cart::z][1]);
R_temp[Cart::z][Cart::yyy][0]+=wmc1*R_temp[Cart::z][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::z][Cart::y][0]-cfak*R_temp[Cart::z][Cart::y][1]);
R_temp[Cart::z][Cart::xyy][0]+=wmc0*R_temp[Cart::z][Cart::yy][1];
R_temp[Cart::z][Cart::yyz][0]+=wmc2*R_temp[Cart::z][Cart::yy][1]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][1];
R_temp[Cart::z][Cart::xxy][0]+=wmc1*R_temp[Cart::z][Cart::xx][1];
R_temp[Cart::z][Cart::xyz][0]+=wmc0*R_temp[Cart::z][Cart::yz][1];
R_temp[Cart::z][Cart::yzz][0]+=wmc1*R_temp[Cart::z][Cart::zz][1];
R_temp[Cart::z][Cart::xxx][0]+=wmc0*R_temp[Cart::z][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::z][Cart::x][0]-cfak*R_temp[Cart::z][Cart::x][1]);
R_temp[Cart::z][Cart::xxz][0]+=wmc2*R_temp[Cart::z][Cart::xx][1]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][1];
R_temp[Cart::z][Cart::xzz][0]+=wmc0*R_temp[Cart::z][Cart::zz][1];
R_temp[Cart::z][Cart::zzz][0]+=wmc2*R_temp[Cart::z][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::z][Cart::z][0]-cfak*R_temp[Cart::z][Cart::z][1])+0.5/_decay*1*R_temp[Cart::s][Cart::zz][1];
}
//------------------------------------------------------

//Integral p - s - f - m1
if (_mmax >1 ){
if (_lmax_alpha>0 && _lmax_gamma>2){
R_temp[Cart::y][Cart::yyy][1]+=wmc1*R_temp[Cart::y][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::y][Cart::y][1]-cfak*R_temp[Cart::y][Cart::y][2])+0.5/_decay*1*R_temp[Cart::s][Cart::yy][2];
R_temp[Cart::y][Cart::xyy][1]+=wmc0*R_temp[Cart::y][Cart::yy][2];
R_temp[Cart::y][Cart::yyz][1]+=wmc2*R_temp[Cart::y][Cart::yy][2];
R_temp[Cart::y][Cart::xxy][1]+=wmc1*R_temp[Cart::y][Cart::xx][2]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][2];
R_temp[Cart::y][Cart::xyz][1]+=wmc0*R_temp[Cart::y][Cart::yz][2];
R_temp[Cart::y][Cart::yzz][1]+=wmc1*R_temp[Cart::y][Cart::zz][2]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][2];
R_temp[Cart::y][Cart::xxx][1]+=wmc0*R_temp[Cart::y][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::y][Cart::x][1]-cfak*R_temp[Cart::y][Cart::x][2]);
R_temp[Cart::y][Cart::xxz][1]+=wmc2*R_temp[Cart::y][Cart::xx][2];
R_temp[Cart::y][Cart::xzz][1]+=wmc0*R_temp[Cart::y][Cart::zz][2];
R_temp[Cart::y][Cart::zzz][1]+=wmc2*R_temp[Cart::y][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::y][Cart::z][1]-cfak*R_temp[Cart::y][Cart::z][2]);
R_temp[Cart::x][Cart::yyy][1]+=wmc1*R_temp[Cart::x][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::x][Cart::y][1]-cfak*R_temp[Cart::x][Cart::y][2]);
R_temp[Cart::x][Cart::xyy][1]+=wmc0*R_temp[Cart::x][Cart::yy][2]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][2];
R_temp[Cart::x][Cart::yyz][1]+=wmc2*R_temp[Cart::x][Cart::yy][2];
R_temp[Cart::x][Cart::xxy][1]+=wmc1*R_temp[Cart::x][Cart::xx][2];
R_temp[Cart::x][Cart::xyz][1]+=wmc0*R_temp[Cart::x][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::s][Cart::yz][2];
R_temp[Cart::x][Cart::yzz][1]+=wmc1*R_temp[Cart::x][Cart::zz][2];
R_temp[Cart::x][Cart::xxx][1]+=wmc0*R_temp[Cart::x][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::x][Cart::x][1]-cfak*R_temp[Cart::x][Cart::x][2])+0.5/_decay*1*R_temp[Cart::s][Cart::xx][2];
R_temp[Cart::x][Cart::xxz][1]+=wmc2*R_temp[Cart::x][Cart::xx][2];
R_temp[Cart::x][Cart::xzz][1]+=wmc0*R_temp[Cart::x][Cart::zz][2]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][2];
R_temp[Cart::x][Cart::zzz][1]+=wmc2*R_temp[Cart::x][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::x][Cart::z][1]-cfak*R_temp[Cart::x][Cart::z][2]);
R_temp[Cart::z][Cart::yyy][1]+=wmc1*R_temp[Cart::z][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::z][Cart::y][1]-cfak*R_temp[Cart::z][Cart::y][2]);
R_temp[Cart::z][Cart::xyy][1]+=wmc0*R_temp[Cart::z][Cart::yy][2];
R_temp[Cart::z][Cart::yyz][1]+=wmc2*R_temp[Cart::z][Cart::yy][2]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][2];
R_temp[Cart::z][Cart::xxy][1]+=wmc1*R_temp[Cart::z][Cart::xx][2];
R_temp[Cart::z][Cart::xyz][1]+=wmc0*R_temp[Cart::z][Cart::yz][2];
R_temp[Cart::z][Cart::yzz][1]+=wmc1*R_temp[Cart::z][Cart::zz][2];
R_temp[Cart::z][Cart::xxx][1]+=wmc0*R_temp[Cart::z][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::z][Cart::x][1]-cfak*R_temp[Cart::z][Cart::x][2]);
R_temp[Cart::z][Cart::xxz][1]+=wmc2*R_temp[Cart::z][Cart::xx][2]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][2];
R_temp[Cart::z][Cart::xzz][1]+=wmc0*R_temp[Cart::z][Cart::zz][2];
R_temp[Cart::z][Cart::zzz][1]+=wmc2*R_temp[Cart::z][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::z][Cart::z][1]-cfak*R_temp[Cart::z][Cart::z][2])+0.5/_decay*1*R_temp[Cart::s][Cart::zz][2];
}
}
//------------------------------------------------------

//Integral p - s - f - m2
if (_mmax >2 ){
if (_lmax_alpha>0 && _lmax_gamma>2){
R_temp[Cart::y][Cart::yyy][2]+=wmc1*R_temp[Cart::y][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::y][Cart::y][2]-cfak*R_temp[Cart::y][Cart::y][3])+0.5/_decay*1*R_temp[Cart::s][Cart::yy][3];
R_temp[Cart::y][Cart::xyy][2]+=wmc0*R_temp[Cart::y][Cart::yy][3];
R_temp[Cart::y][Cart::yyz][2]+=wmc2*R_temp[Cart::y][Cart::yy][3];
R_temp[Cart::y][Cart::xxy][2]+=wmc1*R_temp[Cart::y][Cart::xx][3]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][3];
R_temp[Cart::y][Cart::xyz][2]+=wmc0*R_temp[Cart::y][Cart::yz][3];
R_temp[Cart::y][Cart::yzz][2]+=wmc1*R_temp[Cart::y][Cart::zz][3]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][3];
R_temp[Cart::y][Cart::xxx][2]+=wmc0*R_temp[Cart::y][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::y][Cart::x][2]-cfak*R_temp[Cart::y][Cart::x][3]);
R_temp[Cart::y][Cart::xxz][2]+=wmc2*R_temp[Cart::y][Cart::xx][3];
R_temp[Cart::y][Cart::xzz][2]+=wmc0*R_temp[Cart::y][Cart::zz][3];
R_temp[Cart::y][Cart::zzz][2]+=wmc2*R_temp[Cart::y][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::y][Cart::z][2]-cfak*R_temp[Cart::y][Cart::z][3]);
R_temp[Cart::x][Cart::yyy][2]+=wmc1*R_temp[Cart::x][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::x][Cart::y][2]-cfak*R_temp[Cart::x][Cart::y][3]);
R_temp[Cart::x][Cart::xyy][2]+=wmc0*R_temp[Cart::x][Cart::yy][3]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][3];
R_temp[Cart::x][Cart::yyz][2]+=wmc2*R_temp[Cart::x][Cart::yy][3];
R_temp[Cart::x][Cart::xxy][2]+=wmc1*R_temp[Cart::x][Cart::xx][3];
R_temp[Cart::x][Cart::xyz][2]+=wmc0*R_temp[Cart::x][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::s][Cart::yz][3];
R_temp[Cart::x][Cart::yzz][2]+=wmc1*R_temp[Cart::x][Cart::zz][3];
R_temp[Cart::x][Cart::xxx][2]+=wmc0*R_temp[Cart::x][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::x][Cart::x][2]-cfak*R_temp[Cart::x][Cart::x][3])+0.5/_decay*1*R_temp[Cart::s][Cart::xx][3];
R_temp[Cart::x][Cart::xxz][2]+=wmc2*R_temp[Cart::x][Cart::xx][3];
R_temp[Cart::x][Cart::xzz][2]+=wmc0*R_temp[Cart::x][Cart::zz][3]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][3];
R_temp[Cart::x][Cart::zzz][2]+=wmc2*R_temp[Cart::x][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::x][Cart::z][2]-cfak*R_temp[Cart::x][Cart::z][3]);
R_temp[Cart::z][Cart::yyy][2]+=wmc1*R_temp[Cart::z][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::z][Cart::y][2]-cfak*R_temp[Cart::z][Cart::y][3]);
R_temp[Cart::z][Cart::xyy][2]+=wmc0*R_temp[Cart::z][Cart::yy][3];
R_temp[Cart::z][Cart::yyz][2]+=wmc2*R_temp[Cart::z][Cart::yy][3]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][3];
R_temp[Cart::z][Cart::xxy][2]+=wmc1*R_temp[Cart::z][Cart::xx][3];
R_temp[Cart::z][Cart::xyz][2]+=wmc0*R_temp[Cart::z][Cart::yz][3];
R_temp[Cart::z][Cart::yzz][2]+=wmc1*R_temp[Cart::z][Cart::zz][3];
R_temp[Cart::z][Cart::xxx][2]+=wmc0*R_temp[Cart::z][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::z][Cart::x][2]-cfak*R_temp[Cart::z][Cart::x][3]);
R_temp[Cart::z][Cart::xxz][2]+=wmc2*R_temp[Cart::z][Cart::xx][3]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][3];
R_temp[Cart::z][Cart::xzz][2]+=wmc0*R_temp[Cart::z][Cart::zz][3];
R_temp[Cart::z][Cart::zzz][2]+=wmc2*R_temp[Cart::z][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::z][Cart::z][2]-cfak*R_temp[Cart::z][Cart::z][3])+0.5/_decay*1*R_temp[Cart::s][Cart::zz][3];
}
}
//------------------------------------------------------

//Integral p - s - f - m3
if (_mmax >3 ){
if (_lmax_alpha>0 && _lmax_gamma>2){
R_temp[Cart::y][Cart::yyy][3]+=wmc1*R_temp[Cart::y][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::y][Cart::y][3]-cfak*R_temp[Cart::y][Cart::y][4])+0.5/_decay*1*R_temp[Cart::s][Cart::yy][4];
R_temp[Cart::y][Cart::xyy][3]+=wmc0*R_temp[Cart::y][Cart::yy][4];
R_temp[Cart::y][Cart::yyz][3]+=wmc2*R_temp[Cart::y][Cart::yy][4];
R_temp[Cart::y][Cart::xxy][3]+=wmc1*R_temp[Cart::y][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][4];
R_temp[Cart::y][Cart::xyz][3]+=wmc0*R_temp[Cart::y][Cart::yz][4];
R_temp[Cart::y][Cart::yzz][3]+=wmc1*R_temp[Cart::y][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][4];
R_temp[Cart::y][Cart::xxx][3]+=wmc0*R_temp[Cart::y][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::y][Cart::x][3]-cfak*R_temp[Cart::y][Cart::x][4]);
R_temp[Cart::y][Cart::xxz][3]+=wmc2*R_temp[Cart::y][Cart::xx][4];
R_temp[Cart::y][Cart::xzz][3]+=wmc0*R_temp[Cart::y][Cart::zz][4];
R_temp[Cart::y][Cart::zzz][3]+=wmc2*R_temp[Cart::y][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::y][Cart::z][3]-cfak*R_temp[Cart::y][Cart::z][4]);
R_temp[Cart::x][Cart::yyy][3]+=wmc1*R_temp[Cart::x][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::x][Cart::y][3]-cfak*R_temp[Cart::x][Cart::y][4]);
R_temp[Cart::x][Cart::xyy][3]+=wmc0*R_temp[Cart::x][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][4];
R_temp[Cart::x][Cart::yyz][3]+=wmc2*R_temp[Cart::x][Cart::yy][4];
R_temp[Cart::x][Cart::xxy][3]+=wmc1*R_temp[Cart::x][Cart::xx][4];
R_temp[Cart::x][Cart::xyz][3]+=wmc0*R_temp[Cart::x][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::s][Cart::yz][4];
R_temp[Cart::x][Cart::yzz][3]+=wmc1*R_temp[Cart::x][Cart::zz][4];
R_temp[Cart::x][Cart::xxx][3]+=wmc0*R_temp[Cart::x][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::x][Cart::x][3]-cfak*R_temp[Cart::x][Cart::x][4])+0.5/_decay*1*R_temp[Cart::s][Cart::xx][4];
R_temp[Cart::x][Cart::xxz][3]+=wmc2*R_temp[Cart::x][Cart::xx][4];
R_temp[Cart::x][Cart::xzz][3]+=wmc0*R_temp[Cart::x][Cart::zz][4]+0.5/_decay*1*R_temp[Cart::s][Cart::zz][4];
R_temp[Cart::x][Cart::zzz][3]+=wmc2*R_temp[Cart::x][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::x][Cart::z][3]-cfak*R_temp[Cart::x][Cart::z][4]);
R_temp[Cart::z][Cart::yyy][3]+=wmc1*R_temp[Cart::z][Cart::yy][4]+1/_decay_gamma*(R_temp[Cart::z][Cart::y][3]-cfak*R_temp[Cart::z][Cart::y][4]);
R_temp[Cart::z][Cart::xyy][3]+=wmc0*R_temp[Cart::z][Cart::yy][4];
R_temp[Cart::z][Cart::yyz][3]+=wmc2*R_temp[Cart::z][Cart::yy][4]+0.5/_decay*1*R_temp[Cart::s][Cart::yy][4];
R_temp[Cart::z][Cart::xxy][3]+=wmc1*R_temp[Cart::z][Cart::xx][4];
R_temp[Cart::z][Cart::xyz][3]+=wmc0*R_temp[Cart::z][Cart::yz][4];
R_temp[Cart::z][Cart::yzz][3]+=wmc1*R_temp[Cart::z][Cart::zz][4];
R_temp[Cart::z][Cart::xxx][3]+=wmc0*R_temp[Cart::z][Cart::xx][4]+1/_decay_gamma*(R_temp[Cart::z][Cart::x][3]-cfak*R_temp[Cart::z][Cart::x][4]);
R_temp[Cart::z][Cart::xxz][3]+=wmc2*R_temp[Cart::z][Cart::xx][4]+0.5/_decay*1*R_temp[Cart::s][Cart::xx][4];
R_temp[Cart::z][Cart::xzz][3]+=wmc0*R_temp[Cart::z][Cart::zz][4];
R_temp[Cart::z][Cart::zzz][3]+=wmc2*R_temp[Cart::z][Cart::zz][4]+1/_decay_gamma*(R_temp[Cart::z][Cart::z][3]-cfak*R_temp[Cart::z][Cart::z][4])+0.5/_decay*1*R_temp[Cart::s][Cart::zz][4];
}
}
//------------------------------------------------------

//Integral d - s - s - m0
if (_lmax_alpha>1){
R_temp[Cart::yy][Cart::s][0]+=pma1*R_temp[Cart::y][Cart::s][0]+wmp1*R_temp[Cart::y][Cart::s][1];
R_temp[Cart::xy][Cart::s][0]+=pma0*R_temp[Cart::y][Cart::s][0]+wmp0*R_temp[Cart::y][Cart::s][1];
R_temp[Cart::yz][Cart::s][0]+=pma1*R_temp[Cart::z][Cart::s][0]+wmp1*R_temp[Cart::z][Cart::s][1];
R_temp[Cart::xx][Cart::s][0]+=pma0*R_temp[Cart::x][Cart::s][0]+wmp0*R_temp[Cart::x][Cart::s][1];
R_temp[Cart::xz][Cart::s][0]+=pma0*R_temp[Cart::z][Cart::s][0]+wmp0*R_temp[Cart::z][Cart::s][1];
R_temp[Cart::zz][Cart::s][0]+=pma2*R_temp[Cart::z][Cart::s][0]+wmp2*R_temp[Cart::z][Cart::s][1];
}
//------------------------------------------------------

//Integral d - s - s - m1
if (_mmax >1 ){
if (_lmax_alpha>1){
R_temp[Cart::yy][Cart::s][1]+=pma1*R_temp[Cart::y][Cart::s][1]+wmp1*R_temp[Cart::y][Cart::s][2];
R_temp[Cart::xy][Cart::s][1]+=pma0*R_temp[Cart::y][Cart::s][1]+wmp0*R_temp[Cart::y][Cart::s][2];
R_temp[Cart::yz][Cart::s][1]+=pma1*R_temp[Cart::z][Cart::s][1]+wmp1*R_temp[Cart::z][Cart::s][2];
R_temp[Cart::xx][Cart::s][1]+=pma0*R_temp[Cart::x][Cart::s][1]+wmp0*R_temp[Cart::x][Cart::s][2];
R_temp[Cart::xz][Cart::s][1]+=pma0*R_temp[Cart::z][Cart::s][1]+wmp0*R_temp[Cart::z][Cart::s][2];
R_temp[Cart::zz][Cart::s][1]+=pma2*R_temp[Cart::z][Cart::s][1]+wmp2*R_temp[Cart::z][Cart::s][2];
}
}
//------------------------------------------------------

//Integral d - s - s - m2
if (_mmax >2 ){
if (_lmax_alpha>1){
R_temp[Cart::yy][Cart::s][2]+=pma1*R_temp[Cart::y][Cart::s][2]+wmp1*R_temp[Cart::y][Cart::s][3];
R_temp[Cart::xy][Cart::s][2]+=pma0*R_temp[Cart::y][Cart::s][2]+wmp0*R_temp[Cart::y][Cart::s][3];
R_temp[Cart::yz][Cart::s][2]+=pma1*R_temp[Cart::z][Cart::s][2]+wmp1*R_temp[Cart::z][Cart::s][3];
R_temp[Cart::xx][Cart::s][2]+=pma0*R_temp[Cart::x][Cart::s][2]+wmp0*R_temp[Cart::x][Cart::s][3];
R_temp[Cart::xz][Cart::s][2]+=pma0*R_temp[Cart::z][Cart::s][2]+wmp0*R_temp[Cart::z][Cart::s][3];
R_temp[Cart::zz][Cart::s][2]+=pma2*R_temp[Cart::z][Cart::s][2]+wmp2*R_temp[Cart::z][Cart::s][3];
}
}
//------------------------------------------------------

//Integral d - s - s - m3
if (_mmax >3 ){
if (_lmax_alpha>1){
R_temp[Cart::yy][Cart::s][3]+=pma1*R_temp[Cart::y][Cart::s][3]+wmp1*R_temp[Cart::y][Cart::s][4];
R_temp[Cart::xy][Cart::s][3]+=pma0*R_temp[Cart::y][Cart::s][3]+wmp0*R_temp[Cart::y][Cart::s][4];
R_temp[Cart::yz][Cart::s][3]+=pma1*R_temp[Cart::z][Cart::s][3]+wmp1*R_temp[Cart::z][Cart::s][4];
R_temp[Cart::xx][Cart::s][3]+=pma0*R_temp[Cart::x][Cart::s][3]+wmp0*R_temp[Cart::x][Cart::s][4];
R_temp[Cart::xz][Cart::s][3]+=pma0*R_temp[Cart::z][Cart::s][3]+wmp0*R_temp[Cart::z][Cart::s][4];
R_temp[Cart::zz][Cart::s][3]+=pma2*R_temp[Cart::z][Cart::s][3]+wmp2*R_temp[Cart::z][Cart::s][4];
}
}
//------------------------------------------------------

//Integral d - s - s - m4
if (_mmax >4 ){
if (_lmax_alpha>1){
R_temp[Cart::yy][Cart::s][4]+=pma1*R_temp[Cart::y][Cart::s][4]+wmp1*R_temp[Cart::y][Cart::s][5];
R_temp[Cart::xy][Cart::s][4]+=pma0*R_temp[Cart::y][Cart::s][4]+wmp0*R_temp[Cart::y][Cart::s][5];
R_temp[Cart::yz][Cart::s][4]+=pma1*R_temp[Cart::z][Cart::s][4]+wmp1*R_temp[Cart::z][Cart::s][5];
R_temp[Cart::xx][Cart::s][4]+=pma0*R_temp[Cart::x][Cart::s][4]+wmp0*R_temp[Cart::x][Cart::s][5];
R_temp[Cart::xz][Cart::s][4]+=pma0*R_temp[Cart::z][Cart::s][4]+wmp0*R_temp[Cart::z][Cart::s][5];
R_temp[Cart::zz][Cart::s][4]+=pma2*R_temp[Cart::z][Cart::s][4]+wmp2*R_temp[Cart::z][Cart::s][5];
}
}
//------------------------------------------------------

//Integral d - s - s - m5
if (_mmax >5 ){
if (_lmax_alpha>1){
R_temp[Cart::yy][Cart::s][5]+=pma1*R_temp[Cart::y][Cart::s][5]+wmp1*R_temp[Cart::y][Cart::s][6];
R_temp[Cart::xy][Cart::s][5]+=pma0*R_temp[Cart::y][Cart::s][5]+wmp0*R_temp[Cart::y][Cart::s][6];
R_temp[Cart::yz][Cart::s][5]+=pma1*R_temp[Cart::z][Cart::s][5]+wmp1*R_temp[Cart::z][Cart::s][6];
R_temp[Cart::xx][Cart::s][5]+=pma0*R_temp[Cart::x][Cart::s][5]+wmp0*R_temp[Cart::x][Cart::s][6];
R_temp[Cart::xz][Cart::s][5]+=pma0*R_temp[Cart::z][Cart::s][5]+wmp0*R_temp[Cart::z][Cart::s][6];
R_temp[Cart::zz][Cart::s][5]+=pma2*R_temp[Cart::z][Cart::s][5]+wmp2*R_temp[Cart::z][Cart::s][6];
}
}
//------------------------------------------------------

//Integral d - s - p - m0
if (_lmax_alpha>1 && _lmax_gamma>0){
R_temp[Cart::yy][Cart::y][0]+=pma1*R_temp[Cart::y][Cart::y][0]+wmp1*R_temp[Cart::y][Cart::y][1]+0.5/_decay*1*R_temp[Cart::y][Cart::s][1];
R_temp[Cart::yy][Cart::x][0]+=pma1*R_temp[Cart::y][Cart::x][0]+wmp1*R_temp[Cart::y][Cart::x][1];
R_temp[Cart::yy][Cart::z][0]+=pma1*R_temp[Cart::y][Cart::z][0]+wmp1*R_temp[Cart::y][Cart::z][1];
R_temp[Cart::xy][Cart::y][0]+=pma0*R_temp[Cart::y][Cart::y][0]+wmp0*R_temp[Cart::y][Cart::y][1];
R_temp[Cart::xy][Cart::x][0]+=pma0*R_temp[Cart::y][Cart::x][0]+wmp0*R_temp[Cart::y][Cart::x][1]+0.5/_decay*1*R_temp[Cart::y][Cart::s][1];
R_temp[Cart::xy][Cart::z][0]+=pma0*R_temp[Cart::y][Cart::z][0]+wmp0*R_temp[Cart::y][Cart::z][1];
R_temp[Cart::yz][Cart::y][0]+=pma1*R_temp[Cart::z][Cart::y][0]+wmp1*R_temp[Cart::z][Cart::y][1]+0.5/_decay*1*R_temp[Cart::z][Cart::s][1];
R_temp[Cart::yz][Cart::x][0]+=pma1*R_temp[Cart::z][Cart::x][0]+wmp1*R_temp[Cart::z][Cart::x][1];
R_temp[Cart::yz][Cart::z][0]+=pma1*R_temp[Cart::z][Cart::z][0]+wmp1*R_temp[Cart::z][Cart::z][1];
R_temp[Cart::xx][Cart::y][0]+=pma0*R_temp[Cart::x][Cart::y][0]+wmp0*R_temp[Cart::x][Cart::y][1];
R_temp[Cart::xx][Cart::x][0]+=pma0*R_temp[Cart::x][Cart::x][0]+wmp0*R_temp[Cart::x][Cart::x][1]+0.5/_decay*1*R_temp[Cart::x][Cart::s][1];
R_temp[Cart::xx][Cart::z][0]+=pma0*R_temp[Cart::x][Cart::z][0]+wmp0*R_temp[Cart::x][Cart::z][1];
R_temp[Cart::xz][Cart::y][0]+=pma0*R_temp[Cart::z][Cart::y][0]+wmp0*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::xz][Cart::x][0]+=pma0*R_temp[Cart::z][Cart::x][0]+wmp0*R_temp[Cart::z][Cart::x][1]+0.5/_decay*1*R_temp[Cart::z][Cart::s][1];
R_temp[Cart::xz][Cart::z][0]+=pma0*R_temp[Cart::z][Cart::z][0]+wmp0*R_temp[Cart::z][Cart::z][1];
R_temp[Cart::zz][Cart::y][0]+=pma2*R_temp[Cart::z][Cart::y][0]+wmp2*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::zz][Cart::x][0]+=pma2*R_temp[Cart::z][Cart::x][0]+wmp2*R_temp[Cart::z][Cart::x][1];
R_temp[Cart::zz][Cart::z][0]+=pma2*R_temp[Cart::z][Cart::z][0]+wmp2*R_temp[Cart::z][Cart::z][1]+0.5/_decay*1*R_temp[Cart::z][Cart::s][1];
}
//------------------------------------------------------

//Integral d - s - p - m1
if (_mmax >1 ){
if (_lmax_alpha>1 && _lmax_gamma>0){
R_temp[Cart::yy][Cart::y][1]+=pma1*R_temp[Cart::y][Cart::y][1]+wmp1*R_temp[Cart::y][Cart::y][2]+0.5/_decay*1*R_temp[Cart::y][Cart::s][2];
R_temp[Cart::yy][Cart::x][1]+=pma1*R_temp[Cart::y][Cart::x][1]+wmp1*R_temp[Cart::y][Cart::x][2];
R_temp[Cart::yy][Cart::z][1]+=pma1*R_temp[Cart::y][Cart::z][1]+wmp1*R_temp[Cart::y][Cart::z][2];
R_temp[Cart::xy][Cart::y][1]+=pma0*R_temp[Cart::y][Cart::y][1]+wmp0*R_temp[Cart::y][Cart::y][2];
R_temp[Cart::xy][Cart::x][1]+=pma0*R_temp[Cart::y][Cart::x][1]+wmp0*R_temp[Cart::y][Cart::x][2]+0.5/_decay*1*R_temp[Cart::y][Cart::s][2];
R_temp[Cart::xy][Cart::z][1]+=pma0*R_temp[Cart::y][Cart::z][1]+wmp0*R_temp[Cart::y][Cart::z][2];
R_temp[Cart::yz][Cart::y][1]+=pma1*R_temp[Cart::z][Cart::y][1]+wmp1*R_temp[Cart::z][Cart::y][2]+0.5/_decay*1*R_temp[Cart::z][Cart::s][2];
R_temp[Cart::yz][Cart::x][1]+=pma1*R_temp[Cart::z][Cart::x][1]+wmp1*R_temp[Cart::z][Cart::x][2];
R_temp[Cart::yz][Cart::z][1]+=pma1*R_temp[Cart::z][Cart::z][1]+wmp1*R_temp[Cart::z][Cart::z][2];
R_temp[Cart::xx][Cart::y][1]+=pma0*R_temp[Cart::x][Cart::y][1]+wmp0*R_temp[Cart::x][Cart::y][2];
R_temp[Cart::xx][Cart::x][1]+=pma0*R_temp[Cart::x][Cart::x][1]+wmp0*R_temp[Cart::x][Cart::x][2]+0.5/_decay*1*R_temp[Cart::x][Cart::s][2];
R_temp[Cart::xx][Cart::z][1]+=pma0*R_temp[Cart::x][Cart::z][1]+wmp0*R_temp[Cart::x][Cart::z][2];
R_temp[Cart::xz][Cart::y][1]+=pma0*R_temp[Cart::z][Cart::y][1]+wmp0*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::xz][Cart::x][1]+=pma0*R_temp[Cart::z][Cart::x][1]+wmp0*R_temp[Cart::z][Cart::x][2]+0.5/_decay*1*R_temp[Cart::z][Cart::s][2];
R_temp[Cart::xz][Cart::z][1]+=pma0*R_temp[Cart::z][Cart::z][1]+wmp0*R_temp[Cart::z][Cart::z][2];
R_temp[Cart::zz][Cart::y][1]+=pma2*R_temp[Cart::z][Cart::y][1]+wmp2*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::zz][Cart::x][1]+=pma2*R_temp[Cart::z][Cart::x][1]+wmp2*R_temp[Cart::z][Cart::x][2];
R_temp[Cart::zz][Cart::z][1]+=pma2*R_temp[Cart::z][Cart::z][1]+wmp2*R_temp[Cart::z][Cart::z][2]+0.5/_decay*1*R_temp[Cart::z][Cart::s][2];
}
}
//------------------------------------------------------

//Integral d - s - p - m2
if (_mmax >2 ){
if (_lmax_alpha>1 && _lmax_gamma>0){
R_temp[Cart::yy][Cart::y][2]+=pma1*R_temp[Cart::y][Cart::y][2]+wmp1*R_temp[Cart::y][Cart::y][3]+0.5/_decay*1*R_temp[Cart::y][Cart::s][3];
R_temp[Cart::yy][Cart::x][2]+=pma1*R_temp[Cart::y][Cart::x][2]+wmp1*R_temp[Cart::y][Cart::x][3];
R_temp[Cart::yy][Cart::z][2]+=pma1*R_temp[Cart::y][Cart::z][2]+wmp1*R_temp[Cart::y][Cart::z][3];
R_temp[Cart::xy][Cart::y][2]+=pma0*R_temp[Cart::y][Cart::y][2]+wmp0*R_temp[Cart::y][Cart::y][3];
R_temp[Cart::xy][Cart::x][2]+=pma0*R_temp[Cart::y][Cart::x][2]+wmp0*R_temp[Cart::y][Cart::x][3]+0.5/_decay*1*R_temp[Cart::y][Cart::s][3];
R_temp[Cart::xy][Cart::z][2]+=pma0*R_temp[Cart::y][Cart::z][2]+wmp0*R_temp[Cart::y][Cart::z][3];
R_temp[Cart::yz][Cart::y][2]+=pma1*R_temp[Cart::z][Cart::y][2]+wmp1*R_temp[Cart::z][Cart::y][3]+0.5/_decay*1*R_temp[Cart::z][Cart::s][3];
R_temp[Cart::yz][Cart::x][2]+=pma1*R_temp[Cart::z][Cart::x][2]+wmp1*R_temp[Cart::z][Cart::x][3];
R_temp[Cart::yz][Cart::z][2]+=pma1*R_temp[Cart::z][Cart::z][2]+wmp1*R_temp[Cart::z][Cart::z][3];
R_temp[Cart::xx][Cart::y][2]+=pma0*R_temp[Cart::x][Cart::y][2]+wmp0*R_temp[Cart::x][Cart::y][3];
R_temp[Cart::xx][Cart::x][2]+=pma0*R_temp[Cart::x][Cart::x][2]+wmp0*R_temp[Cart::x][Cart::x][3]+0.5/_decay*1*R_temp[Cart::x][Cart::s][3];
R_temp[Cart::xx][Cart::z][2]+=pma0*R_temp[Cart::x][Cart::z][2]+wmp0*R_temp[Cart::x][Cart::z][3];
R_temp[Cart::xz][Cart::y][2]+=pma0*R_temp[Cart::z][Cart::y][2]+wmp0*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::xz][Cart::x][2]+=pma0*R_temp[Cart::z][Cart::x][2]+wmp0*R_temp[Cart::z][Cart::x][3]+0.5/_decay*1*R_temp[Cart::z][Cart::s][3];
R_temp[Cart::xz][Cart::z][2]+=pma0*R_temp[Cart::z][Cart::z][2]+wmp0*R_temp[Cart::z][Cart::z][3];
R_temp[Cart::zz][Cart::y][2]+=pma2*R_temp[Cart::z][Cart::y][2]+wmp2*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::zz][Cart::x][2]+=pma2*R_temp[Cart::z][Cart::x][2]+wmp2*R_temp[Cart::z][Cart::x][3];
R_temp[Cart::zz][Cart::z][2]+=pma2*R_temp[Cart::z][Cart::z][2]+wmp2*R_temp[Cart::z][Cart::z][3]+0.5/_decay*1*R_temp[Cart::z][Cart::s][3];
}
}
//------------------------------------------------------

//Integral d - s - p - m3
if (_mmax >3 ){
if (_lmax_alpha>1 && _lmax_gamma>0){
R_temp[Cart::yy][Cart::y][3]+=pma1*R_temp[Cart::y][Cart::y][3]+wmp1*R_temp[Cart::y][Cart::y][4]+0.5/_decay*1*R_temp[Cart::y][Cart::s][4];
R_temp[Cart::yy][Cart::x][3]+=pma1*R_temp[Cart::y][Cart::x][3]+wmp1*R_temp[Cart::y][Cart::x][4];
R_temp[Cart::yy][Cart::z][3]+=pma1*R_temp[Cart::y][Cart::z][3]+wmp1*R_temp[Cart::y][Cart::z][4];
R_temp[Cart::xy][Cart::y][3]+=pma0*R_temp[Cart::y][Cart::y][3]+wmp0*R_temp[Cart::y][Cart::y][4];
R_temp[Cart::xy][Cart::x][3]+=pma0*R_temp[Cart::y][Cart::x][3]+wmp0*R_temp[Cart::y][Cart::x][4]+0.5/_decay*1*R_temp[Cart::y][Cart::s][4];
R_temp[Cart::xy][Cart::z][3]+=pma0*R_temp[Cart::y][Cart::z][3]+wmp0*R_temp[Cart::y][Cart::z][4];
R_temp[Cart::yz][Cart::y][3]+=pma1*R_temp[Cart::z][Cart::y][3]+wmp1*R_temp[Cart::z][Cart::y][4]+0.5/_decay*1*R_temp[Cart::z][Cart::s][4];
R_temp[Cart::yz][Cart::x][3]+=pma1*R_temp[Cart::z][Cart::x][3]+wmp1*R_temp[Cart::z][Cart::x][4];
R_temp[Cart::yz][Cart::z][3]+=pma1*R_temp[Cart::z][Cart::z][3]+wmp1*R_temp[Cart::z][Cart::z][4];
R_temp[Cart::xx][Cart::y][3]+=pma0*R_temp[Cart::x][Cart::y][3]+wmp0*R_temp[Cart::x][Cart::y][4];
R_temp[Cart::xx][Cart::x][3]+=pma0*R_temp[Cart::x][Cart::x][3]+wmp0*R_temp[Cart::x][Cart::x][4]+0.5/_decay*1*R_temp[Cart::x][Cart::s][4];
R_temp[Cart::xx][Cart::z][3]+=pma0*R_temp[Cart::x][Cart::z][3]+wmp0*R_temp[Cart::x][Cart::z][4];
R_temp[Cart::xz][Cart::y][3]+=pma0*R_temp[Cart::z][Cart::y][3]+wmp0*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::xz][Cart::x][3]+=pma0*R_temp[Cart::z][Cart::x][3]+wmp0*R_temp[Cart::z][Cart::x][4]+0.5/_decay*1*R_temp[Cart::z][Cart::s][4];
R_temp[Cart::xz][Cart::z][3]+=pma0*R_temp[Cart::z][Cart::z][3]+wmp0*R_temp[Cart::z][Cart::z][4];
R_temp[Cart::zz][Cart::y][3]+=pma2*R_temp[Cart::z][Cart::y][3]+wmp2*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::zz][Cart::x][3]+=pma2*R_temp[Cart::z][Cart::x][3]+wmp2*R_temp[Cart::z][Cart::x][4];
R_temp[Cart::zz][Cart::z][3]+=pma2*R_temp[Cart::z][Cart::z][3]+wmp2*R_temp[Cart::z][Cart::z][4]+0.5/_decay*1*R_temp[Cart::z][Cart::s][4];
}
}
//------------------------------------------------------

//Integral d - s - p - m4
if (_mmax >4 ){
if (_lmax_alpha>1 && _lmax_gamma>0){
R_temp[Cart::yy][Cart::y][4]+=pma1*R_temp[Cart::y][Cart::y][4]+wmp1*R_temp[Cart::y][Cart::y][5]+0.5/_decay*1*R_temp[Cart::y][Cart::s][5];
R_temp[Cart::yy][Cart::x][4]+=pma1*R_temp[Cart::y][Cart::x][4]+wmp1*R_temp[Cart::y][Cart::x][5];
R_temp[Cart::yy][Cart::z][4]+=pma1*R_temp[Cart::y][Cart::z][4]+wmp1*R_temp[Cart::y][Cart::z][5];
R_temp[Cart::xy][Cart::y][4]+=pma0*R_temp[Cart::y][Cart::y][4]+wmp0*R_temp[Cart::y][Cart::y][5];
R_temp[Cart::xy][Cart::x][4]+=pma0*R_temp[Cart::y][Cart::x][4]+wmp0*R_temp[Cart::y][Cart::x][5]+0.5/_decay*1*R_temp[Cart::y][Cart::s][5];
R_temp[Cart::xy][Cart::z][4]+=pma0*R_temp[Cart::y][Cart::z][4]+wmp0*R_temp[Cart::y][Cart::z][5];
R_temp[Cart::yz][Cart::y][4]+=pma1*R_temp[Cart::z][Cart::y][4]+wmp1*R_temp[Cart::z][Cart::y][5]+0.5/_decay*1*R_temp[Cart::z][Cart::s][5];
R_temp[Cart::yz][Cart::x][4]+=pma1*R_temp[Cart::z][Cart::x][4]+wmp1*R_temp[Cart::z][Cart::x][5];
R_temp[Cart::yz][Cart::z][4]+=pma1*R_temp[Cart::z][Cart::z][4]+wmp1*R_temp[Cart::z][Cart::z][5];
R_temp[Cart::xx][Cart::y][4]+=pma0*R_temp[Cart::x][Cart::y][4]+wmp0*R_temp[Cart::x][Cart::y][5];
R_temp[Cart::xx][Cart::x][4]+=pma0*R_temp[Cart::x][Cart::x][4]+wmp0*R_temp[Cart::x][Cart::x][5]+0.5/_decay*1*R_temp[Cart::x][Cart::s][5];
R_temp[Cart::xx][Cart::z][4]+=pma0*R_temp[Cart::x][Cart::z][4]+wmp0*R_temp[Cart::x][Cart::z][5];
R_temp[Cart::xz][Cart::y][4]+=pma0*R_temp[Cart::z][Cart::y][4]+wmp0*R_temp[Cart::z][Cart::y][5];
R_temp[Cart::xz][Cart::x][4]+=pma0*R_temp[Cart::z][Cart::x][4]+wmp0*R_temp[Cart::z][Cart::x][5]+0.5/_decay*1*R_temp[Cart::z][Cart::s][5];
R_temp[Cart::xz][Cart::z][4]+=pma0*R_temp[Cart::z][Cart::z][4]+wmp0*R_temp[Cart::z][Cart::z][5];
R_temp[Cart::zz][Cart::y][4]+=pma2*R_temp[Cart::z][Cart::y][4]+wmp2*R_temp[Cart::z][Cart::y][5];
R_temp[Cart::zz][Cart::x][4]+=pma2*R_temp[Cart::z][Cart::x][4]+wmp2*R_temp[Cart::z][Cart::x][5];
R_temp[Cart::zz][Cart::z][4]+=pma2*R_temp[Cart::z][Cart::z][4]+wmp2*R_temp[Cart::z][Cart::z][5]+0.5/_decay*1*R_temp[Cart::z][Cart::s][5];
}
}
//------------------------------------------------------

//Integral d - s - d - m0
if (_lmax_alpha>1 && _lmax_gamma>1){
R_temp[Cart::yy][Cart::yy][0]+=pma1*R_temp[Cart::y][Cart::yy][0]+wmp1*R_temp[Cart::y][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::y][Cart::y][1];
R_temp[Cart::yy][Cart::xy][0]+=pma1*R_temp[Cart::y][Cart::xy][0]+wmp1*R_temp[Cart::y][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::y][Cart::x][1];
R_temp[Cart::yy][Cart::yz][0]+=pma1*R_temp[Cart::y][Cart::yz][0]+wmp1*R_temp[Cart::y][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::y][Cart::z][1];
R_temp[Cart::yy][Cart::xx][0]+=pma1*R_temp[Cart::y][Cart::xx][0]+wmp1*R_temp[Cart::y][Cart::xx][1];
R_temp[Cart::yy][Cart::xz][0]+=pma1*R_temp[Cart::y][Cart::xz][0]+wmp1*R_temp[Cart::y][Cart::xz][1];
R_temp[Cart::yy][Cart::zz][0]+=pma1*R_temp[Cart::y][Cart::zz][0]+wmp1*R_temp[Cart::y][Cart::zz][1];
R_temp[Cart::xy][Cart::yy][0]+=pma0*R_temp[Cart::y][Cart::yy][0]+wmp0*R_temp[Cart::y][Cart::yy][1];
R_temp[Cart::xy][Cart::xy][0]+=pma0*R_temp[Cart::y][Cart::xy][0]+wmp0*R_temp[Cart::y][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::y][Cart::y][1];
R_temp[Cart::xy][Cart::yz][0]+=pma0*R_temp[Cart::y][Cart::yz][0]+wmp0*R_temp[Cart::y][Cart::yz][1];
R_temp[Cart::xy][Cart::xx][0]+=pma0*R_temp[Cart::y][Cart::xx][0]+wmp0*R_temp[Cart::y][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::y][Cart::x][1];
R_temp[Cart::xy][Cart::xz][0]+=pma0*R_temp[Cart::y][Cart::xz][0]+wmp0*R_temp[Cart::y][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::y][Cart::z][1];
R_temp[Cart::xy][Cart::zz][0]+=pma0*R_temp[Cart::y][Cart::zz][0]+wmp0*R_temp[Cart::y][Cart::zz][1];
R_temp[Cart::yz][Cart::yy][0]+=pma1*R_temp[Cart::z][Cart::yy][0]+wmp1*R_temp[Cart::z][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::yz][Cart::xy][0]+=pma1*R_temp[Cart::z][Cart::xy][0]+wmp1*R_temp[Cart::z][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::z][Cart::x][1];
R_temp[Cart::yz][Cart::yz][0]+=pma1*R_temp[Cart::z][Cart::yz][0]+wmp1*R_temp[Cart::z][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::z][1];
R_temp[Cart::yz][Cart::xx][0]+=pma1*R_temp[Cart::z][Cart::xx][0]+wmp1*R_temp[Cart::z][Cart::xx][1];
R_temp[Cart::yz][Cart::xz][0]+=pma1*R_temp[Cart::z][Cart::xz][0]+wmp1*R_temp[Cart::z][Cart::xz][1];
R_temp[Cart::yz][Cart::zz][0]+=pma1*R_temp[Cart::z][Cart::zz][0]+wmp1*R_temp[Cart::z][Cart::zz][1];
R_temp[Cart::xx][Cart::yy][0]+=pma0*R_temp[Cart::x][Cart::yy][0]+wmp0*R_temp[Cart::x][Cart::yy][1];
R_temp[Cart::xx][Cart::xy][0]+=pma0*R_temp[Cart::x][Cart::xy][0]+wmp0*R_temp[Cart::x][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::x][Cart::y][1];
R_temp[Cart::xx][Cart::yz][0]+=pma0*R_temp[Cart::x][Cart::yz][0]+wmp0*R_temp[Cart::x][Cart::yz][1];
R_temp[Cart::xx][Cart::xx][0]+=pma0*R_temp[Cart::x][Cart::xx][0]+wmp0*R_temp[Cart::x][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::x][Cart::x][1];
R_temp[Cart::xx][Cart::xz][0]+=pma0*R_temp[Cart::x][Cart::xz][0]+wmp0*R_temp[Cart::x][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::x][Cart::z][1];
R_temp[Cart::xx][Cart::zz][0]+=pma0*R_temp[Cart::x][Cart::zz][0]+wmp0*R_temp[Cart::x][Cart::zz][1];
R_temp[Cart::xz][Cart::yy][0]+=pma0*R_temp[Cart::z][Cart::yy][0]+wmp0*R_temp[Cart::z][Cart::yy][1];
R_temp[Cart::xz][Cart::xy][0]+=pma0*R_temp[Cart::z][Cart::xy][0]+wmp0*R_temp[Cart::z][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::xz][Cart::yz][0]+=pma0*R_temp[Cart::z][Cart::yz][0]+wmp0*R_temp[Cart::z][Cart::yz][1];
R_temp[Cart::xz][Cart::xx][0]+=pma0*R_temp[Cart::z][Cart::xx][0]+wmp0*R_temp[Cart::z][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::z][Cart::x][1];
R_temp[Cart::xz][Cart::xz][0]+=pma0*R_temp[Cart::z][Cart::xz][0]+wmp0*R_temp[Cart::z][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::z][1];
R_temp[Cart::xz][Cart::zz][0]+=pma0*R_temp[Cart::z][Cart::zz][0]+wmp0*R_temp[Cart::z][Cart::zz][1];
R_temp[Cart::zz][Cart::yy][0]+=pma2*R_temp[Cart::z][Cart::yy][0]+wmp2*R_temp[Cart::z][Cart::yy][1];
R_temp[Cart::zz][Cart::xy][0]+=pma2*R_temp[Cart::z][Cart::xy][0]+wmp2*R_temp[Cart::z][Cart::xy][1];
R_temp[Cart::zz][Cart::yz][0]+=pma2*R_temp[Cart::z][Cart::yz][0]+wmp2*R_temp[Cart::z][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::y][1];
R_temp[Cart::zz][Cart::xx][0]+=pma2*R_temp[Cart::z][Cart::xx][0]+wmp2*R_temp[Cart::z][Cart::xx][1];
R_temp[Cart::zz][Cart::xz][0]+=pma2*R_temp[Cart::z][Cart::xz][0]+wmp2*R_temp[Cart::z][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::x][1];
R_temp[Cart::zz][Cart::zz][0]+=pma2*R_temp[Cart::z][Cart::zz][0]+wmp2*R_temp[Cart::z][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::z][Cart::z][1];
}
//------------------------------------------------------

//Integral d - s - d - m1
if (_mmax >1 ){
if (_lmax_alpha>1 && _lmax_gamma>1){
R_temp[Cart::yy][Cart::yy][1]+=pma1*R_temp[Cart::y][Cart::yy][1]+wmp1*R_temp[Cart::y][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::y][Cart::y][2];
R_temp[Cart::yy][Cart::xy][1]+=pma1*R_temp[Cart::y][Cart::xy][1]+wmp1*R_temp[Cart::y][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::y][Cart::x][2];
R_temp[Cart::yy][Cart::yz][1]+=pma1*R_temp[Cart::y][Cart::yz][1]+wmp1*R_temp[Cart::y][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::y][Cart::z][2];
R_temp[Cart::yy][Cart::xx][1]+=pma1*R_temp[Cart::y][Cart::xx][1]+wmp1*R_temp[Cart::y][Cart::xx][2];
R_temp[Cart::yy][Cart::xz][1]+=pma1*R_temp[Cart::y][Cart::xz][1]+wmp1*R_temp[Cart::y][Cart::xz][2];
R_temp[Cart::yy][Cart::zz][1]+=pma1*R_temp[Cart::y][Cart::zz][1]+wmp1*R_temp[Cart::y][Cart::zz][2];
R_temp[Cart::xy][Cart::yy][1]+=pma0*R_temp[Cart::y][Cart::yy][1]+wmp0*R_temp[Cart::y][Cart::yy][2];
R_temp[Cart::xy][Cart::xy][1]+=pma0*R_temp[Cart::y][Cart::xy][1]+wmp0*R_temp[Cart::y][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::y][Cart::y][2];
R_temp[Cart::xy][Cart::yz][1]+=pma0*R_temp[Cart::y][Cart::yz][1]+wmp0*R_temp[Cart::y][Cart::yz][2];
R_temp[Cart::xy][Cart::xx][1]+=pma0*R_temp[Cart::y][Cart::xx][1]+wmp0*R_temp[Cart::y][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::y][Cart::x][2];
R_temp[Cart::xy][Cart::xz][1]+=pma0*R_temp[Cart::y][Cart::xz][1]+wmp0*R_temp[Cart::y][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::y][Cart::z][2];
R_temp[Cart::xy][Cart::zz][1]+=pma0*R_temp[Cart::y][Cart::zz][1]+wmp0*R_temp[Cart::y][Cart::zz][2];
R_temp[Cart::yz][Cart::yy][1]+=pma1*R_temp[Cart::z][Cart::yy][1]+wmp1*R_temp[Cart::z][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::yz][Cart::xy][1]+=pma1*R_temp[Cart::z][Cart::xy][1]+wmp1*R_temp[Cart::z][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::z][Cart::x][2];
R_temp[Cart::yz][Cart::yz][1]+=pma1*R_temp[Cart::z][Cart::yz][1]+wmp1*R_temp[Cart::z][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::z][2];
R_temp[Cart::yz][Cart::xx][1]+=pma1*R_temp[Cart::z][Cart::xx][1]+wmp1*R_temp[Cart::z][Cart::xx][2];
R_temp[Cart::yz][Cart::xz][1]+=pma1*R_temp[Cart::z][Cart::xz][1]+wmp1*R_temp[Cart::z][Cart::xz][2];
R_temp[Cart::yz][Cart::zz][1]+=pma1*R_temp[Cart::z][Cart::zz][1]+wmp1*R_temp[Cart::z][Cart::zz][2];
R_temp[Cart::xx][Cart::yy][1]+=pma0*R_temp[Cart::x][Cart::yy][1]+wmp0*R_temp[Cart::x][Cart::yy][2];
R_temp[Cart::xx][Cart::xy][1]+=pma0*R_temp[Cart::x][Cart::xy][1]+wmp0*R_temp[Cart::x][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::x][Cart::y][2];
R_temp[Cart::xx][Cart::yz][1]+=pma0*R_temp[Cart::x][Cart::yz][1]+wmp0*R_temp[Cart::x][Cart::yz][2];
R_temp[Cart::xx][Cart::xx][1]+=pma0*R_temp[Cart::x][Cart::xx][1]+wmp0*R_temp[Cart::x][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::x][Cart::x][2];
R_temp[Cart::xx][Cart::xz][1]+=pma0*R_temp[Cart::x][Cart::xz][1]+wmp0*R_temp[Cart::x][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::x][Cart::z][2];
R_temp[Cart::xx][Cart::zz][1]+=pma0*R_temp[Cart::x][Cart::zz][1]+wmp0*R_temp[Cart::x][Cart::zz][2];
R_temp[Cart::xz][Cart::yy][1]+=pma0*R_temp[Cart::z][Cart::yy][1]+wmp0*R_temp[Cart::z][Cart::yy][2];
R_temp[Cart::xz][Cart::xy][1]+=pma0*R_temp[Cart::z][Cart::xy][1]+wmp0*R_temp[Cart::z][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::xz][Cart::yz][1]+=pma0*R_temp[Cart::z][Cart::yz][1]+wmp0*R_temp[Cart::z][Cart::yz][2];
R_temp[Cart::xz][Cart::xx][1]+=pma0*R_temp[Cart::z][Cart::xx][1]+wmp0*R_temp[Cart::z][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::z][Cart::x][2];
R_temp[Cart::xz][Cart::xz][1]+=pma0*R_temp[Cart::z][Cart::xz][1]+wmp0*R_temp[Cart::z][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::z][2];
R_temp[Cart::xz][Cart::zz][1]+=pma0*R_temp[Cart::z][Cart::zz][1]+wmp0*R_temp[Cart::z][Cart::zz][2];
R_temp[Cart::zz][Cart::yy][1]+=pma2*R_temp[Cart::z][Cart::yy][1]+wmp2*R_temp[Cart::z][Cart::yy][2];
R_temp[Cart::zz][Cart::xy][1]+=pma2*R_temp[Cart::z][Cart::xy][1]+wmp2*R_temp[Cart::z][Cart::xy][2];
R_temp[Cart::zz][Cart::yz][1]+=pma2*R_temp[Cart::z][Cart::yz][1]+wmp2*R_temp[Cart::z][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::y][2];
R_temp[Cart::zz][Cart::xx][1]+=pma2*R_temp[Cart::z][Cart::xx][1]+wmp2*R_temp[Cart::z][Cart::xx][2];
R_temp[Cart::zz][Cart::xz][1]+=pma2*R_temp[Cart::z][Cart::xz][1]+wmp2*R_temp[Cart::z][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::x][2];
R_temp[Cart::zz][Cart::zz][1]+=pma2*R_temp[Cart::z][Cart::zz][1]+wmp2*R_temp[Cart::z][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::z][Cart::z][2];
}
}
//------------------------------------------------------

//Integral d - s - d - m2
if (_mmax >2 ){
if (_lmax_alpha>1 && _lmax_gamma>1){
R_temp[Cart::yy][Cart::yy][2]+=pma1*R_temp[Cart::y][Cart::yy][2]+wmp1*R_temp[Cart::y][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::y][Cart::y][3];
R_temp[Cart::yy][Cart::xy][2]+=pma1*R_temp[Cart::y][Cart::xy][2]+wmp1*R_temp[Cart::y][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::y][Cart::x][3];
R_temp[Cart::yy][Cart::yz][2]+=pma1*R_temp[Cart::y][Cart::yz][2]+wmp1*R_temp[Cart::y][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::y][Cart::z][3];
R_temp[Cart::yy][Cart::xx][2]+=pma1*R_temp[Cart::y][Cart::xx][2]+wmp1*R_temp[Cart::y][Cart::xx][3];
R_temp[Cart::yy][Cart::xz][2]+=pma1*R_temp[Cart::y][Cart::xz][2]+wmp1*R_temp[Cart::y][Cart::xz][3];
R_temp[Cart::yy][Cart::zz][2]+=pma1*R_temp[Cart::y][Cart::zz][2]+wmp1*R_temp[Cart::y][Cart::zz][3];
R_temp[Cart::xy][Cart::yy][2]+=pma0*R_temp[Cart::y][Cart::yy][2]+wmp0*R_temp[Cart::y][Cart::yy][3];
R_temp[Cart::xy][Cart::xy][2]+=pma0*R_temp[Cart::y][Cart::xy][2]+wmp0*R_temp[Cart::y][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::y][Cart::y][3];
R_temp[Cart::xy][Cart::yz][2]+=pma0*R_temp[Cart::y][Cart::yz][2]+wmp0*R_temp[Cart::y][Cart::yz][3];
R_temp[Cart::xy][Cart::xx][2]+=pma0*R_temp[Cart::y][Cart::xx][2]+wmp0*R_temp[Cart::y][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::y][Cart::x][3];
R_temp[Cart::xy][Cart::xz][2]+=pma0*R_temp[Cart::y][Cart::xz][2]+wmp0*R_temp[Cart::y][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::y][Cart::z][3];
R_temp[Cart::xy][Cart::zz][2]+=pma0*R_temp[Cart::y][Cart::zz][2]+wmp0*R_temp[Cart::y][Cart::zz][3];
R_temp[Cart::yz][Cart::yy][2]+=pma1*R_temp[Cart::z][Cart::yy][2]+wmp1*R_temp[Cart::z][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::yz][Cart::xy][2]+=pma1*R_temp[Cart::z][Cart::xy][2]+wmp1*R_temp[Cart::z][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::z][Cart::x][3];
R_temp[Cart::yz][Cart::yz][2]+=pma1*R_temp[Cart::z][Cart::yz][2]+wmp1*R_temp[Cart::z][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::z][3];
R_temp[Cart::yz][Cart::xx][2]+=pma1*R_temp[Cart::z][Cart::xx][2]+wmp1*R_temp[Cart::z][Cart::xx][3];
R_temp[Cart::yz][Cart::xz][2]+=pma1*R_temp[Cart::z][Cart::xz][2]+wmp1*R_temp[Cart::z][Cart::xz][3];
R_temp[Cart::yz][Cart::zz][2]+=pma1*R_temp[Cart::z][Cart::zz][2]+wmp1*R_temp[Cart::z][Cart::zz][3];
R_temp[Cart::xx][Cart::yy][2]+=pma0*R_temp[Cart::x][Cart::yy][2]+wmp0*R_temp[Cart::x][Cart::yy][3];
R_temp[Cart::xx][Cart::xy][2]+=pma0*R_temp[Cart::x][Cart::xy][2]+wmp0*R_temp[Cart::x][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::x][Cart::y][3];
R_temp[Cart::xx][Cart::yz][2]+=pma0*R_temp[Cart::x][Cart::yz][2]+wmp0*R_temp[Cart::x][Cart::yz][3];
R_temp[Cart::xx][Cart::xx][2]+=pma0*R_temp[Cart::x][Cart::xx][2]+wmp0*R_temp[Cart::x][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::x][Cart::x][3];
R_temp[Cart::xx][Cart::xz][2]+=pma0*R_temp[Cart::x][Cart::xz][2]+wmp0*R_temp[Cart::x][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::x][Cart::z][3];
R_temp[Cart::xx][Cart::zz][2]+=pma0*R_temp[Cart::x][Cart::zz][2]+wmp0*R_temp[Cart::x][Cart::zz][3];
R_temp[Cart::xz][Cart::yy][2]+=pma0*R_temp[Cart::z][Cart::yy][2]+wmp0*R_temp[Cart::z][Cart::yy][3];
R_temp[Cart::xz][Cart::xy][2]+=pma0*R_temp[Cart::z][Cart::xy][2]+wmp0*R_temp[Cart::z][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::xz][Cart::yz][2]+=pma0*R_temp[Cart::z][Cart::yz][2]+wmp0*R_temp[Cart::z][Cart::yz][3];
R_temp[Cart::xz][Cart::xx][2]+=pma0*R_temp[Cart::z][Cart::xx][2]+wmp0*R_temp[Cart::z][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::z][Cart::x][3];
R_temp[Cart::xz][Cart::xz][2]+=pma0*R_temp[Cart::z][Cart::xz][2]+wmp0*R_temp[Cart::z][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::z][3];
R_temp[Cart::xz][Cart::zz][2]+=pma0*R_temp[Cart::z][Cart::zz][2]+wmp0*R_temp[Cart::z][Cart::zz][3];
R_temp[Cart::zz][Cart::yy][2]+=pma2*R_temp[Cart::z][Cart::yy][2]+wmp2*R_temp[Cart::z][Cart::yy][3];
R_temp[Cart::zz][Cart::xy][2]+=pma2*R_temp[Cart::z][Cart::xy][2]+wmp2*R_temp[Cart::z][Cart::xy][3];
R_temp[Cart::zz][Cart::yz][2]+=pma2*R_temp[Cart::z][Cart::yz][2]+wmp2*R_temp[Cart::z][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::y][3];
R_temp[Cart::zz][Cart::xx][2]+=pma2*R_temp[Cart::z][Cart::xx][2]+wmp2*R_temp[Cart::z][Cart::xx][3];
R_temp[Cart::zz][Cart::xz][2]+=pma2*R_temp[Cart::z][Cart::xz][2]+wmp2*R_temp[Cart::z][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::x][3];
R_temp[Cart::zz][Cart::zz][2]+=pma2*R_temp[Cart::z][Cart::zz][2]+wmp2*R_temp[Cart::z][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::z][Cart::z][3];
}
}
//------------------------------------------------------

//Integral d - s - d - m3
if (_mmax >3 ){
if (_lmax_alpha>1 && _lmax_gamma>1){
R_temp[Cart::yy][Cart::yy][3]+=pma1*R_temp[Cart::y][Cart::yy][3]+wmp1*R_temp[Cart::y][Cart::yy][4]+0.5/_decay*2*R_temp[Cart::y][Cart::y][4];
R_temp[Cart::yy][Cart::xy][3]+=pma1*R_temp[Cart::y][Cart::xy][3]+wmp1*R_temp[Cart::y][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::y][Cart::x][4];
R_temp[Cart::yy][Cart::yz][3]+=pma1*R_temp[Cart::y][Cart::yz][3]+wmp1*R_temp[Cart::y][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::y][Cart::z][4];
R_temp[Cart::yy][Cart::xx][3]+=pma1*R_temp[Cart::y][Cart::xx][3]+wmp1*R_temp[Cart::y][Cart::xx][4];
R_temp[Cart::yy][Cart::xz][3]+=pma1*R_temp[Cart::y][Cart::xz][3]+wmp1*R_temp[Cart::y][Cart::xz][4];
R_temp[Cart::yy][Cart::zz][3]+=pma1*R_temp[Cart::y][Cart::zz][3]+wmp1*R_temp[Cart::y][Cart::zz][4];
R_temp[Cart::xy][Cart::yy][3]+=pma0*R_temp[Cart::y][Cart::yy][3]+wmp0*R_temp[Cart::y][Cart::yy][4];
R_temp[Cart::xy][Cart::xy][3]+=pma0*R_temp[Cart::y][Cart::xy][3]+wmp0*R_temp[Cart::y][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::y][Cart::y][4];
R_temp[Cart::xy][Cart::yz][3]+=pma0*R_temp[Cart::y][Cart::yz][3]+wmp0*R_temp[Cart::y][Cart::yz][4];
R_temp[Cart::xy][Cart::xx][3]+=pma0*R_temp[Cart::y][Cart::xx][3]+wmp0*R_temp[Cart::y][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::y][Cart::x][4];
R_temp[Cart::xy][Cart::xz][3]+=pma0*R_temp[Cart::y][Cart::xz][3]+wmp0*R_temp[Cart::y][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::y][Cart::z][4];
R_temp[Cart::xy][Cart::zz][3]+=pma0*R_temp[Cart::y][Cart::zz][3]+wmp0*R_temp[Cart::y][Cart::zz][4];
R_temp[Cart::yz][Cart::yy][3]+=pma1*R_temp[Cart::z][Cart::yy][3]+wmp1*R_temp[Cart::z][Cart::yy][4]+0.5/_decay*2*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::yz][Cart::xy][3]+=pma1*R_temp[Cart::z][Cart::xy][3]+wmp1*R_temp[Cart::z][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::z][Cart::x][4];
R_temp[Cart::yz][Cart::yz][3]+=pma1*R_temp[Cart::z][Cart::yz][3]+wmp1*R_temp[Cart::z][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::z][4];
R_temp[Cart::yz][Cart::xx][3]+=pma1*R_temp[Cart::z][Cart::xx][3]+wmp1*R_temp[Cart::z][Cart::xx][4];
R_temp[Cart::yz][Cart::xz][3]+=pma1*R_temp[Cart::z][Cart::xz][3]+wmp1*R_temp[Cart::z][Cart::xz][4];
R_temp[Cart::yz][Cart::zz][3]+=pma1*R_temp[Cart::z][Cart::zz][3]+wmp1*R_temp[Cart::z][Cart::zz][4];
R_temp[Cart::xx][Cart::yy][3]+=pma0*R_temp[Cart::x][Cart::yy][3]+wmp0*R_temp[Cart::x][Cart::yy][4];
R_temp[Cart::xx][Cart::xy][3]+=pma0*R_temp[Cart::x][Cart::xy][3]+wmp0*R_temp[Cart::x][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::x][Cart::y][4];
R_temp[Cart::xx][Cart::yz][3]+=pma0*R_temp[Cart::x][Cart::yz][3]+wmp0*R_temp[Cart::x][Cart::yz][4];
R_temp[Cart::xx][Cart::xx][3]+=pma0*R_temp[Cart::x][Cart::xx][3]+wmp0*R_temp[Cart::x][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::x][Cart::x][4];
R_temp[Cart::xx][Cart::xz][3]+=pma0*R_temp[Cart::x][Cart::xz][3]+wmp0*R_temp[Cart::x][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::x][Cart::z][4];
R_temp[Cart::xx][Cart::zz][3]+=pma0*R_temp[Cart::x][Cart::zz][3]+wmp0*R_temp[Cart::x][Cart::zz][4];
R_temp[Cart::xz][Cart::yy][3]+=pma0*R_temp[Cart::z][Cart::yy][3]+wmp0*R_temp[Cart::z][Cart::yy][4];
R_temp[Cart::xz][Cart::xy][3]+=pma0*R_temp[Cart::z][Cart::xy][3]+wmp0*R_temp[Cart::z][Cart::xy][4]+0.5/_decay*1*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::xz][Cart::yz][3]+=pma0*R_temp[Cart::z][Cart::yz][3]+wmp0*R_temp[Cart::z][Cart::yz][4];
R_temp[Cart::xz][Cart::xx][3]+=pma0*R_temp[Cart::z][Cart::xx][3]+wmp0*R_temp[Cart::z][Cart::xx][4]+0.5/_decay*2*R_temp[Cart::z][Cart::x][4];
R_temp[Cart::xz][Cart::xz][3]+=pma0*R_temp[Cart::z][Cart::xz][3]+wmp0*R_temp[Cart::z][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::z][4];
R_temp[Cart::xz][Cart::zz][3]+=pma0*R_temp[Cart::z][Cart::zz][3]+wmp0*R_temp[Cart::z][Cart::zz][4];
R_temp[Cart::zz][Cart::yy][3]+=pma2*R_temp[Cart::z][Cart::yy][3]+wmp2*R_temp[Cart::z][Cart::yy][4];
R_temp[Cart::zz][Cart::xy][3]+=pma2*R_temp[Cart::z][Cart::xy][3]+wmp2*R_temp[Cart::z][Cart::xy][4];
R_temp[Cart::zz][Cart::yz][3]+=pma2*R_temp[Cart::z][Cart::yz][3]+wmp2*R_temp[Cart::z][Cart::yz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::y][4];
R_temp[Cart::zz][Cart::xx][3]+=pma2*R_temp[Cart::z][Cart::xx][3]+wmp2*R_temp[Cart::z][Cart::xx][4];
R_temp[Cart::zz][Cart::xz][3]+=pma2*R_temp[Cart::z][Cart::xz][3]+wmp2*R_temp[Cart::z][Cart::xz][4]+0.5/_decay*1*R_temp[Cart::z][Cart::x][4];
R_temp[Cart::zz][Cart::zz][3]+=pma2*R_temp[Cart::z][Cart::zz][3]+wmp2*R_temp[Cart::z][Cart::zz][4]+0.5/_decay*2*R_temp[Cart::z][Cart::z][4];
}
}
//------------------------------------------------------

//Integral d - s - f - m0
if (_lmax_alpha>1 && _lmax_gamma>2){
R_temp[Cart::yy][Cart::yyy][0]+=wmc1*R_temp[Cart::yy][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::yy][Cart::y][0]-cfak*R_temp[Cart::yy][Cart::y][1])+0.5/_decay*2*R_temp[Cart::y][Cart::yy][1];
R_temp[Cart::yy][Cart::xyy][0]+=wmc0*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::yy][Cart::yyz][0]+=wmc2*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::yy][Cart::xxy][0]+=wmc1*R_temp[Cart::yy][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::y][Cart::xx][1];
R_temp[Cart::yy][Cart::xyz][0]+=wmc0*R_temp[Cart::yy][Cart::yz][1];
R_temp[Cart::yy][Cart::yzz][0]+=wmc1*R_temp[Cart::yy][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::y][Cart::zz][1];
R_temp[Cart::yy][Cart::xxx][0]+=wmc0*R_temp[Cart::yy][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::yy][Cart::x][0]-cfak*R_temp[Cart::yy][Cart::x][1]);
R_temp[Cart::yy][Cart::xxz][0]+=wmc2*R_temp[Cart::yy][Cart::xx][1];
R_temp[Cart::yy][Cart::xzz][0]+=wmc0*R_temp[Cart::yy][Cart::zz][1];
R_temp[Cart::yy][Cart::zzz][0]+=wmc2*R_temp[Cart::yy][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::yy][Cart::z][0]-cfak*R_temp[Cart::yy][Cart::z][1]);
R_temp[Cart::xy][Cart::yyy][0]+=wmc1*R_temp[Cart::xy][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::xy][Cart::y][0]-cfak*R_temp[Cart::xy][Cart::y][1])+0.5/_decay*1*R_temp[Cart::x][Cart::yy][1];
R_temp[Cart::xy][Cart::xyy][0]+=wmc0*R_temp[Cart::xy][Cart::yy][1]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][1];
R_temp[Cart::xy][Cart::yyz][0]+=wmc2*R_temp[Cart::xy][Cart::yy][1];
R_temp[Cart::xy][Cart::xxy][0]+=wmc1*R_temp[Cart::xy][Cart::xx][1]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][1];
R_temp[Cart::xy][Cart::xyz][0]+=wmc0*R_temp[Cart::xy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::y][Cart::yz][1];
R_temp[Cart::xy][Cart::yzz][0]+=wmc1*R_temp[Cart::xy][Cart::zz][1]+0.5/_decay*1*R_temp[Cart::x][Cart::zz][1];
R_temp[Cart::xy][Cart::xxx][0]+=wmc0*R_temp[Cart::xy][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::xy][Cart::x][0]-cfak*R_temp[Cart::xy][Cart::x][1])+0.5/_decay*1*R_temp[Cart::y][Cart::xx][1];
R_temp[Cart::xy][Cart::xxz][0]+=wmc2*R_temp[Cart::xy][Cart::xx][1];
R_temp[Cart::xy][Cart::xzz][0]+=wmc0*R_temp[Cart::xy][Cart::zz][1]+0.5/_decay*1*R_temp[Cart::y][Cart::zz][1];
R_temp[Cart::xy][Cart::zzz][0]+=wmc2*R_temp[Cart::xy][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::xy][Cart::z][0]-cfak*R_temp[Cart::xy][Cart::z][1]);
R_temp[Cart::yz][Cart::yyy][0]+=wmc1*R_temp[Cart::yz][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::yz][Cart::y][0]-cfak*R_temp[Cart::yz][Cart::y][1])+0.5/_decay*1*R_temp[Cart::z][Cart::yy][1];
R_temp[Cart::yz][Cart::xyy][0]+=wmc0*R_temp[Cart::yz][Cart::yy][1];
R_temp[Cart::yz][Cart::yyz][0]+=wmc2*R_temp[Cart::yz][Cart::yy][1]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][1];
R_temp[Cart::yz][Cart::xxy][0]+=wmc1*R_temp[Cart::yz][Cart::xx][1]+0.5/_decay*1*R_temp[Cart::z][Cart::xx][1];
R_temp[Cart::yz][Cart::xyz][0]+=wmc0*R_temp[Cart::yz][Cart::yz][1];
R_temp[Cart::yz][Cart::yzz][0]+=wmc1*R_temp[Cart::yz][Cart::zz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][1];
R_temp[Cart::yz][Cart::xxx][0]+=wmc0*R_temp[Cart::yz][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::yz][Cart::x][0]-cfak*R_temp[Cart::yz][Cart::x][1]);
R_temp[Cart::yz][Cart::xxz][0]+=wmc2*R_temp[Cart::yz][Cart::xx][1]+0.5/_decay*1*R_temp[Cart::y][Cart::xx][1];
R_temp[Cart::yz][Cart::xzz][0]+=wmc0*R_temp[Cart::yz][Cart::zz][1];
R_temp[Cart::yz][Cart::zzz][0]+=wmc2*R_temp[Cart::yz][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::yz][Cart::z][0]-cfak*R_temp[Cart::yz][Cart::z][1])+0.5/_decay*1*R_temp[Cart::y][Cart::zz][1];
R_temp[Cart::xx][Cart::yyy][0]+=wmc1*R_temp[Cart::xx][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::xx][Cart::y][0]-cfak*R_temp[Cart::xx][Cart::y][1]);
R_temp[Cart::xx][Cart::xyy][0]+=wmc0*R_temp[Cart::xx][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::x][Cart::yy][1];
R_temp[Cart::xx][Cart::yyz][0]+=wmc2*R_temp[Cart::xx][Cart::yy][1];
R_temp[Cart::xx][Cart::xxy][0]+=wmc1*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xx][Cart::xyz][0]+=wmc0*R_temp[Cart::xx][Cart::yz][1]+0.5/_decay*2*R_temp[Cart::x][Cart::yz][1];
R_temp[Cart::xx][Cart::yzz][0]+=wmc1*R_temp[Cart::xx][Cart::zz][1];
R_temp[Cart::xx][Cart::xxx][0]+=wmc0*R_temp[Cart::xx][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::xx][Cart::x][0]-cfak*R_temp[Cart::xx][Cart::x][1])+0.5/_decay*2*R_temp[Cart::x][Cart::xx][1];
R_temp[Cart::xx][Cart::xxz][0]+=wmc2*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xx][Cart::xzz][0]+=wmc0*R_temp[Cart::xx][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::x][Cart::zz][1];
R_temp[Cart::xx][Cart::zzz][0]+=wmc2*R_temp[Cart::xx][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::xx][Cart::z][0]-cfak*R_temp[Cart::xx][Cart::z][1]);
R_temp[Cart::xz][Cart::yyy][0]+=wmc1*R_temp[Cart::xz][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::xz][Cart::y][0]-cfak*R_temp[Cart::xz][Cart::y][1]);
R_temp[Cart::xz][Cart::xyy][0]+=wmc0*R_temp[Cart::xz][Cart::yy][1]+0.5/_decay*1*R_temp[Cart::z][Cart::yy][1];
R_temp[Cart::xz][Cart::yyz][0]+=wmc2*R_temp[Cart::xz][Cart::yy][1]+0.5/_decay*1*R_temp[Cart::x][Cart::yy][1];
R_temp[Cart::xz][Cart::xxy][0]+=wmc1*R_temp[Cart::xz][Cart::xx][1];
R_temp[Cart::xz][Cart::xyz][0]+=wmc0*R_temp[Cart::xz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::yz][1];
R_temp[Cart::xz][Cart::yzz][0]+=wmc1*R_temp[Cart::xz][Cart::zz][1];
R_temp[Cart::xz][Cart::xxx][0]+=wmc0*R_temp[Cart::xz][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::xz][Cart::x][0]-cfak*R_temp[Cart::xz][Cart::x][1])+0.5/_decay*1*R_temp[Cart::z][Cart::xx][1];
R_temp[Cart::xz][Cart::xxz][0]+=wmc2*R_temp[Cart::xz][Cart::xx][1]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][1];
R_temp[Cart::xz][Cart::xzz][0]+=wmc0*R_temp[Cart::xz][Cart::zz][1]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][1];
R_temp[Cart::xz][Cart::zzz][0]+=wmc2*R_temp[Cart::xz][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::xz][Cart::z][0]-cfak*R_temp[Cart::xz][Cart::z][1])+0.5/_decay*1*R_temp[Cart::x][Cart::zz][1];
R_temp[Cart::zz][Cart::yyy][0]+=wmc1*R_temp[Cart::zz][Cart::yy][1]+1/_decay_gamma*(R_temp[Cart::zz][Cart::y][0]-cfak*R_temp[Cart::zz][Cart::y][1]);
R_temp[Cart::zz][Cart::xyy][0]+=wmc0*R_temp[Cart::zz][Cart::yy][1];
R_temp[Cart::zz][Cart::yyz][0]+=wmc2*R_temp[Cart::zz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::z][Cart::yy][1];
R_temp[Cart::zz][Cart::xxy][0]+=wmc1*R_temp[Cart::zz][Cart::xx][1];
R_temp[Cart::zz][Cart::xyz][0]+=wmc0*R_temp[Cart::zz][Cart::yz][1];
R_temp[Cart::zz][Cart::yzz][0]+=wmc1*R_temp[Cart::zz][Cart::zz][1];
R_temp[Cart::zz][Cart::xxx][0]+=wmc0*R_temp[Cart::zz][Cart::xx][1]+1/_decay_gamma*(R_temp[Cart::zz][Cart::x][0]-cfak*R_temp[Cart::zz][Cart::x][1]);
R_temp[Cart::zz][Cart::xxz][0]+=wmc2*R_temp[Cart::zz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::z][Cart::xx][1];
R_temp[Cart::zz][Cart::xzz][0]+=wmc0*R_temp[Cart::zz][Cart::zz][1];
R_temp[Cart::zz][Cart::zzz][0]+=wmc2*R_temp[Cart::zz][Cart::zz][1]+1/_decay_gamma*(R_temp[Cart::zz][Cart::z][0]-cfak*R_temp[Cart::zz][Cart::z][1])+0.5/_decay*2*R_temp[Cart::z][Cart::zz][1];
}
//------------------------------------------------------

//Integral d - s - f - m1
if (_mmax >1 ){
if (_lmax_alpha>1 && _lmax_gamma>2){
R_temp[Cart::yy][Cart::yyy][1]+=wmc1*R_temp[Cart::yy][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::yy][Cart::y][1]-cfak*R_temp[Cart::yy][Cart::y][2])+0.5/_decay*2*R_temp[Cart::y][Cart::yy][2];
R_temp[Cart::yy][Cart::xyy][1]+=wmc0*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::yy][Cart::yyz][1]+=wmc2*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::yy][Cart::xxy][1]+=wmc1*R_temp[Cart::yy][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::y][Cart::xx][2];
R_temp[Cart::yy][Cart::xyz][1]+=wmc0*R_temp[Cart::yy][Cart::yz][2];
R_temp[Cart::yy][Cart::yzz][1]+=wmc1*R_temp[Cart::yy][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::y][Cart::zz][2];
R_temp[Cart::yy][Cart::xxx][1]+=wmc0*R_temp[Cart::yy][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::yy][Cart::x][1]-cfak*R_temp[Cart::yy][Cart::x][2]);
R_temp[Cart::yy][Cart::xxz][1]+=wmc2*R_temp[Cart::yy][Cart::xx][2];
R_temp[Cart::yy][Cart::xzz][1]+=wmc0*R_temp[Cart::yy][Cart::zz][2];
R_temp[Cart::yy][Cart::zzz][1]+=wmc2*R_temp[Cart::yy][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::yy][Cart::z][1]-cfak*R_temp[Cart::yy][Cart::z][2]);
R_temp[Cart::xy][Cart::yyy][1]+=wmc1*R_temp[Cart::xy][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::xy][Cart::y][1]-cfak*R_temp[Cart::xy][Cart::y][2])+0.5/_decay*1*R_temp[Cart::x][Cart::yy][2];
R_temp[Cart::xy][Cart::xyy][1]+=wmc0*R_temp[Cart::xy][Cart::yy][2]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][2];
R_temp[Cart::xy][Cart::yyz][1]+=wmc2*R_temp[Cart::xy][Cart::yy][2];
R_temp[Cart::xy][Cart::xxy][1]+=wmc1*R_temp[Cart::xy][Cart::xx][2]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][2];
R_temp[Cart::xy][Cart::xyz][1]+=wmc0*R_temp[Cart::xy][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::y][Cart::yz][2];
R_temp[Cart::xy][Cart::yzz][1]+=wmc1*R_temp[Cart::xy][Cart::zz][2]+0.5/_decay*1*R_temp[Cart::x][Cart::zz][2];
R_temp[Cart::xy][Cart::xxx][1]+=wmc0*R_temp[Cart::xy][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::xy][Cart::x][1]-cfak*R_temp[Cart::xy][Cart::x][2])+0.5/_decay*1*R_temp[Cart::y][Cart::xx][2];
R_temp[Cart::xy][Cart::xxz][1]+=wmc2*R_temp[Cart::xy][Cart::xx][2];
R_temp[Cart::xy][Cart::xzz][1]+=wmc0*R_temp[Cart::xy][Cart::zz][2]+0.5/_decay*1*R_temp[Cart::y][Cart::zz][2];
R_temp[Cart::xy][Cart::zzz][1]+=wmc2*R_temp[Cart::xy][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::xy][Cart::z][1]-cfak*R_temp[Cart::xy][Cart::z][2]);
R_temp[Cart::yz][Cart::yyy][1]+=wmc1*R_temp[Cart::yz][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::yz][Cart::y][1]-cfak*R_temp[Cart::yz][Cart::y][2])+0.5/_decay*1*R_temp[Cart::z][Cart::yy][2];
R_temp[Cart::yz][Cart::xyy][1]+=wmc0*R_temp[Cart::yz][Cart::yy][2];
R_temp[Cart::yz][Cart::yyz][1]+=wmc2*R_temp[Cart::yz][Cart::yy][2]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][2];
R_temp[Cart::yz][Cart::xxy][1]+=wmc1*R_temp[Cart::yz][Cart::xx][2]+0.5/_decay*1*R_temp[Cart::z][Cart::xx][2];
R_temp[Cart::yz][Cart::xyz][1]+=wmc0*R_temp[Cart::yz][Cart::yz][2];
R_temp[Cart::yz][Cart::yzz][1]+=wmc1*R_temp[Cart::yz][Cart::zz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][2];
R_temp[Cart::yz][Cart::xxx][1]+=wmc0*R_temp[Cart::yz][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::yz][Cart::x][1]-cfak*R_temp[Cart::yz][Cart::x][2]);
R_temp[Cart::yz][Cart::xxz][1]+=wmc2*R_temp[Cart::yz][Cart::xx][2]+0.5/_decay*1*R_temp[Cart::y][Cart::xx][2];
R_temp[Cart::yz][Cart::xzz][1]+=wmc0*R_temp[Cart::yz][Cart::zz][2];
R_temp[Cart::yz][Cart::zzz][1]+=wmc2*R_temp[Cart::yz][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::yz][Cart::z][1]-cfak*R_temp[Cart::yz][Cart::z][2])+0.5/_decay*1*R_temp[Cart::y][Cart::zz][2];
R_temp[Cart::xx][Cart::yyy][1]+=wmc1*R_temp[Cart::xx][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::xx][Cart::y][1]-cfak*R_temp[Cart::xx][Cart::y][2]);
R_temp[Cart::xx][Cart::xyy][1]+=wmc0*R_temp[Cart::xx][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::x][Cart::yy][2];
R_temp[Cart::xx][Cart::yyz][1]+=wmc2*R_temp[Cart::xx][Cart::yy][2];
R_temp[Cart::xx][Cart::xxy][1]+=wmc1*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xx][Cart::xyz][1]+=wmc0*R_temp[Cart::xx][Cart::yz][2]+0.5/_decay*2*R_temp[Cart::x][Cart::yz][2];
R_temp[Cart::xx][Cart::yzz][1]+=wmc1*R_temp[Cart::xx][Cart::zz][2];
R_temp[Cart::xx][Cart::xxx][1]+=wmc0*R_temp[Cart::xx][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::xx][Cart::x][1]-cfak*R_temp[Cart::xx][Cart::x][2])+0.5/_decay*2*R_temp[Cart::x][Cart::xx][2];
R_temp[Cart::xx][Cart::xxz][1]+=wmc2*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xx][Cart::xzz][1]+=wmc0*R_temp[Cart::xx][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::x][Cart::zz][2];
R_temp[Cart::xx][Cart::zzz][1]+=wmc2*R_temp[Cart::xx][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::xx][Cart::z][1]-cfak*R_temp[Cart::xx][Cart::z][2]);
R_temp[Cart::xz][Cart::yyy][1]+=wmc1*R_temp[Cart::xz][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::xz][Cart::y][1]-cfak*R_temp[Cart::xz][Cart::y][2]);
R_temp[Cart::xz][Cart::xyy][1]+=wmc0*R_temp[Cart::xz][Cart::yy][2]+0.5/_decay*1*R_temp[Cart::z][Cart::yy][2];
R_temp[Cart::xz][Cart::yyz][1]+=wmc2*R_temp[Cart::xz][Cart::yy][2]+0.5/_decay*1*R_temp[Cart::x][Cart::yy][2];
R_temp[Cart::xz][Cart::xxy][1]+=wmc1*R_temp[Cart::xz][Cart::xx][2];
R_temp[Cart::xz][Cart::xyz][1]+=wmc0*R_temp[Cart::xz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::yz][2];
R_temp[Cart::xz][Cart::yzz][1]+=wmc1*R_temp[Cart::xz][Cart::zz][2];
R_temp[Cart::xz][Cart::xxx][1]+=wmc0*R_temp[Cart::xz][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::xz][Cart::x][1]-cfak*R_temp[Cart::xz][Cart::x][2])+0.5/_decay*1*R_temp[Cart::z][Cart::xx][2];
R_temp[Cart::xz][Cart::xxz][1]+=wmc2*R_temp[Cart::xz][Cart::xx][2]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][2];
R_temp[Cart::xz][Cart::xzz][1]+=wmc0*R_temp[Cart::xz][Cart::zz][2]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][2];
R_temp[Cart::xz][Cart::zzz][1]+=wmc2*R_temp[Cart::xz][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::xz][Cart::z][1]-cfak*R_temp[Cart::xz][Cart::z][2])+0.5/_decay*1*R_temp[Cart::x][Cart::zz][2];
R_temp[Cart::zz][Cart::yyy][1]+=wmc1*R_temp[Cart::zz][Cart::yy][2]+1/_decay_gamma*(R_temp[Cart::zz][Cart::y][1]-cfak*R_temp[Cart::zz][Cart::y][2]);
R_temp[Cart::zz][Cart::xyy][1]+=wmc0*R_temp[Cart::zz][Cart::yy][2];
R_temp[Cart::zz][Cart::yyz][1]+=wmc2*R_temp[Cart::zz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::z][Cart::yy][2];
R_temp[Cart::zz][Cart::xxy][1]+=wmc1*R_temp[Cart::zz][Cart::xx][2];
R_temp[Cart::zz][Cart::xyz][1]+=wmc0*R_temp[Cart::zz][Cart::yz][2];
R_temp[Cart::zz][Cart::yzz][1]+=wmc1*R_temp[Cart::zz][Cart::zz][2];
R_temp[Cart::zz][Cart::xxx][1]+=wmc0*R_temp[Cart::zz][Cart::xx][2]+1/_decay_gamma*(R_temp[Cart::zz][Cart::x][1]-cfak*R_temp[Cart::zz][Cart::x][2]);
R_temp[Cart::zz][Cart::xxz][1]+=wmc2*R_temp[Cart::zz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::z][Cart::xx][2];
R_temp[Cart::zz][Cart::xzz][1]+=wmc0*R_temp[Cart::zz][Cart::zz][2];
R_temp[Cart::zz][Cart::zzz][1]+=wmc2*R_temp[Cart::zz][Cart::zz][2]+1/_decay_gamma*(R_temp[Cart::zz][Cart::z][1]-cfak*R_temp[Cart::zz][Cart::z][2])+0.5/_decay*2*R_temp[Cart::z][Cart::zz][2];
}
}
//------------------------------------------------------

//Integral d - s - f - m2
if (_mmax >2 ){
if (_lmax_alpha>1 && _lmax_gamma>2){
R_temp[Cart::yy][Cart::yyy][2]+=wmc1*R_temp[Cart::yy][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::yy][Cart::y][2]-cfak*R_temp[Cart::yy][Cart::y][3])+0.5/_decay*2*R_temp[Cart::y][Cart::yy][3];
R_temp[Cart::yy][Cart::xyy][2]+=wmc0*R_temp[Cart::yy][Cart::yy][3];
R_temp[Cart::yy][Cart::yyz][2]+=wmc2*R_temp[Cart::yy][Cart::yy][3];
R_temp[Cart::yy][Cart::xxy][2]+=wmc1*R_temp[Cart::yy][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::y][Cart::xx][3];
R_temp[Cart::yy][Cart::xyz][2]+=wmc0*R_temp[Cart::yy][Cart::yz][3];
R_temp[Cart::yy][Cart::yzz][2]+=wmc1*R_temp[Cart::yy][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::y][Cart::zz][3];
R_temp[Cart::yy][Cart::xxx][2]+=wmc0*R_temp[Cart::yy][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::yy][Cart::x][2]-cfak*R_temp[Cart::yy][Cart::x][3]);
R_temp[Cart::yy][Cart::xxz][2]+=wmc2*R_temp[Cart::yy][Cart::xx][3];
R_temp[Cart::yy][Cart::xzz][2]+=wmc0*R_temp[Cart::yy][Cart::zz][3];
R_temp[Cart::yy][Cart::zzz][2]+=wmc2*R_temp[Cart::yy][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::yy][Cart::z][2]-cfak*R_temp[Cart::yy][Cart::z][3]);
R_temp[Cart::xy][Cart::yyy][2]+=wmc1*R_temp[Cart::xy][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::xy][Cart::y][2]-cfak*R_temp[Cart::xy][Cart::y][3])+0.5/_decay*1*R_temp[Cart::x][Cart::yy][3];
R_temp[Cart::xy][Cart::xyy][2]+=wmc0*R_temp[Cart::xy][Cart::yy][3]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][3];
R_temp[Cart::xy][Cart::yyz][2]+=wmc2*R_temp[Cart::xy][Cart::yy][3];
R_temp[Cart::xy][Cart::xxy][2]+=wmc1*R_temp[Cart::xy][Cart::xx][3]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][3];
R_temp[Cart::xy][Cart::xyz][2]+=wmc0*R_temp[Cart::xy][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::y][Cart::yz][3];
R_temp[Cart::xy][Cart::yzz][2]+=wmc1*R_temp[Cart::xy][Cart::zz][3]+0.5/_decay*1*R_temp[Cart::x][Cart::zz][3];
R_temp[Cart::xy][Cart::xxx][2]+=wmc0*R_temp[Cart::xy][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::xy][Cart::x][2]-cfak*R_temp[Cart::xy][Cart::x][3])+0.5/_decay*1*R_temp[Cart::y][Cart::xx][3];
R_temp[Cart::xy][Cart::xxz][2]+=wmc2*R_temp[Cart::xy][Cart::xx][3];
R_temp[Cart::xy][Cart::xzz][2]+=wmc0*R_temp[Cart::xy][Cart::zz][3]+0.5/_decay*1*R_temp[Cart::y][Cart::zz][3];
R_temp[Cart::xy][Cart::zzz][2]+=wmc2*R_temp[Cart::xy][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::xy][Cart::z][2]-cfak*R_temp[Cart::xy][Cart::z][3]);
R_temp[Cart::yz][Cart::yyy][2]+=wmc1*R_temp[Cart::yz][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::yz][Cart::y][2]-cfak*R_temp[Cart::yz][Cart::y][3])+0.5/_decay*1*R_temp[Cart::z][Cart::yy][3];
R_temp[Cart::yz][Cart::xyy][2]+=wmc0*R_temp[Cart::yz][Cart::yy][3];
R_temp[Cart::yz][Cart::yyz][2]+=wmc2*R_temp[Cart::yz][Cart::yy][3]+0.5/_decay*1*R_temp[Cart::y][Cart::yy][3];
R_temp[Cart::yz][Cart::xxy][2]+=wmc1*R_temp[Cart::yz][Cart::xx][3]+0.5/_decay*1*R_temp[Cart::z][Cart::xx][3];
R_temp[Cart::yz][Cart::xyz][2]+=wmc0*R_temp[Cart::yz][Cart::yz][3];
R_temp[Cart::yz][Cart::yzz][2]+=wmc1*R_temp[Cart::yz][Cart::zz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][3];
R_temp[Cart::yz][Cart::xxx][2]+=wmc0*R_temp[Cart::yz][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::yz][Cart::x][2]-cfak*R_temp[Cart::yz][Cart::x][3]);
R_temp[Cart::yz][Cart::xxz][2]+=wmc2*R_temp[Cart::yz][Cart::xx][3]+0.5/_decay*1*R_temp[Cart::y][Cart::xx][3];
R_temp[Cart::yz][Cart::xzz][2]+=wmc0*R_temp[Cart::yz][Cart::zz][3];
R_temp[Cart::yz][Cart::zzz][2]+=wmc2*R_temp[Cart::yz][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::yz][Cart::z][2]-cfak*R_temp[Cart::yz][Cart::z][3])+0.5/_decay*1*R_temp[Cart::y][Cart::zz][3];
R_temp[Cart::xx][Cart::yyy][2]+=wmc1*R_temp[Cart::xx][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::xx][Cart::y][2]-cfak*R_temp[Cart::xx][Cart::y][3]);
R_temp[Cart::xx][Cart::xyy][2]+=wmc0*R_temp[Cart::xx][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::x][Cart::yy][3];
R_temp[Cart::xx][Cart::yyz][2]+=wmc2*R_temp[Cart::xx][Cart::yy][3];
R_temp[Cart::xx][Cart::xxy][2]+=wmc1*R_temp[Cart::xx][Cart::xx][3];
R_temp[Cart::xx][Cart::xyz][2]+=wmc0*R_temp[Cart::xx][Cart::yz][3]+0.5/_decay*2*R_temp[Cart::x][Cart::yz][3];
R_temp[Cart::xx][Cart::yzz][2]+=wmc1*R_temp[Cart::xx][Cart::zz][3];
R_temp[Cart::xx][Cart::xxx][2]+=wmc0*R_temp[Cart::xx][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::xx][Cart::x][2]-cfak*R_temp[Cart::xx][Cart::x][3])+0.5/_decay*2*R_temp[Cart::x][Cart::xx][3];
R_temp[Cart::xx][Cart::xxz][2]+=wmc2*R_temp[Cart::xx][Cart::xx][3];
R_temp[Cart::xx][Cart::xzz][2]+=wmc0*R_temp[Cart::xx][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::x][Cart::zz][3];
R_temp[Cart::xx][Cart::zzz][2]+=wmc2*R_temp[Cart::xx][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::xx][Cart::z][2]-cfak*R_temp[Cart::xx][Cart::z][3]);
R_temp[Cart::xz][Cart::yyy][2]+=wmc1*R_temp[Cart::xz][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::xz][Cart::y][2]-cfak*R_temp[Cart::xz][Cart::y][3]);
R_temp[Cart::xz][Cart::xyy][2]+=wmc0*R_temp[Cart::xz][Cart::yy][3]+0.5/_decay*1*R_temp[Cart::z][Cart::yy][3];
R_temp[Cart::xz][Cart::yyz][2]+=wmc2*R_temp[Cart::xz][Cart::yy][3]+0.5/_decay*1*R_temp[Cart::x][Cart::yy][3];
R_temp[Cart::xz][Cart::xxy][2]+=wmc1*R_temp[Cart::xz][Cart::xx][3];
R_temp[Cart::xz][Cart::xyz][2]+=wmc0*R_temp[Cart::xz][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::yz][3];
R_temp[Cart::xz][Cart::yzz][2]+=wmc1*R_temp[Cart::xz][Cart::zz][3];
R_temp[Cart::xz][Cart::xxx][2]+=wmc0*R_temp[Cart::xz][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::xz][Cart::x][2]-cfak*R_temp[Cart::xz][Cart::x][3])+0.5/_decay*1*R_temp[Cart::z][Cart::xx][3];
R_temp[Cart::xz][Cart::xxz][2]+=wmc2*R_temp[Cart::xz][Cart::xx][3]+0.5/_decay*1*R_temp[Cart::x][Cart::xx][3];
R_temp[Cart::xz][Cart::xzz][2]+=wmc0*R_temp[Cart::xz][Cart::zz][3]+0.5/_decay*1*R_temp[Cart::z][Cart::zz][3];
R_temp[Cart::xz][Cart::zzz][2]+=wmc2*R_temp[Cart::xz][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::xz][Cart::z][2]-cfak*R_temp[Cart::xz][Cart::z][3])+0.5/_decay*1*R_temp[Cart::x][Cart::zz][3];
R_temp[Cart::zz][Cart::yyy][2]+=wmc1*R_temp[Cart::zz][Cart::yy][3]+1/_decay_gamma*(R_temp[Cart::zz][Cart::y][2]-cfak*R_temp[Cart::zz][Cart::y][3]);
R_temp[Cart::zz][Cart::xyy][2]+=wmc0*R_temp[Cart::zz][Cart::yy][3];
R_temp[Cart::zz][Cart::yyz][2]+=wmc2*R_temp[Cart::zz][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::z][Cart::yy][3];
R_temp[Cart::zz][Cart::xxy][2]+=wmc1*R_temp[Cart::zz][Cart::xx][3];
R_temp[Cart::zz][Cart::xyz][2]+=wmc0*R_temp[Cart::zz][Cart::yz][3];
R_temp[Cart::zz][Cart::yzz][2]+=wmc1*R_temp[Cart::zz][Cart::zz][3];
R_temp[Cart::zz][Cart::xxx][2]+=wmc0*R_temp[Cart::zz][Cart::xx][3]+1/_decay_gamma*(R_temp[Cart::zz][Cart::x][2]-cfak*R_temp[Cart::zz][Cart::x][3]);
R_temp[Cart::zz][Cart::xxz][2]+=wmc2*R_temp[Cart::zz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::z][Cart::xx][3];
R_temp[Cart::zz][Cart::xzz][2]+=wmc0*R_temp[Cart::zz][Cart::zz][3];
R_temp[Cart::zz][Cart::zzz][2]+=wmc2*R_temp[Cart::zz][Cart::zz][3]+1/_decay_gamma*(R_temp[Cart::zz][Cart::z][2]-cfak*R_temp[Cart::zz][Cart::z][3])+0.5/_decay*2*R_temp[Cart::z][Cart::zz][3];
}
}
//------------------------------------------------------

//Integral f - s - s - m0
if (_lmax_alpha>2){
R_temp[Cart::yyy][Cart::s][0]+=pma1*R_temp[Cart::yy][Cart::s][0]+wmp1*R_temp[Cart::yy][Cart::s][1]+1*rzeta*(R_temp[Cart::y][Cart::s][0]-gfak*R_temp[Cart::y][Cart::s][1]);
R_temp[Cart::xyy][Cart::s][0]+=pma0*R_temp[Cart::yy][Cart::s][0]+wmp0*R_temp[Cart::yy][Cart::s][1];
R_temp[Cart::yyz][Cart::s][0]+=pma2*R_temp[Cart::yy][Cart::s][0]+wmp2*R_temp[Cart::yy][Cart::s][1];
R_temp[Cart::xxy][Cart::s][0]+=pma1*R_temp[Cart::xx][Cart::s][0]+wmp1*R_temp[Cart::xx][Cart::s][1];
R_temp[Cart::xyz][Cart::s][0]+=pma0*R_temp[Cart::yz][Cart::s][0]+wmp0*R_temp[Cart::yz][Cart::s][1];
R_temp[Cart::yzz][Cart::s][0]+=pma1*R_temp[Cart::zz][Cart::s][0]+wmp1*R_temp[Cart::zz][Cart::s][1];
R_temp[Cart::xxx][Cart::s][0]+=pma0*R_temp[Cart::xx][Cart::s][0]+wmp0*R_temp[Cart::xx][Cart::s][1]+1*rzeta*(R_temp[Cart::x][Cart::s][0]-gfak*R_temp[Cart::x][Cart::s][1]);
R_temp[Cart::xxz][Cart::s][0]+=pma2*R_temp[Cart::xx][Cart::s][0]+wmp2*R_temp[Cart::xx][Cart::s][1];
R_temp[Cart::xzz][Cart::s][0]+=pma0*R_temp[Cart::zz][Cart::s][0]+wmp0*R_temp[Cart::zz][Cart::s][1];
R_temp[Cart::zzz][Cart::s][0]+=pma2*R_temp[Cart::zz][Cart::s][0]+wmp2*R_temp[Cart::zz][Cart::s][1]+1*rzeta*(R_temp[Cart::z][Cart::s][0]-gfak*R_temp[Cart::z][Cart::s][1]);
}
//------------------------------------------------------

//Integral f - s - s - m1
if (_mmax >1 ){
if (_lmax_alpha>2){
R_temp[Cart::yyy][Cart::s][1]+=pma1*R_temp[Cart::yy][Cart::s][1]+wmp1*R_temp[Cart::yy][Cart::s][2]+1*rzeta*(R_temp[Cart::y][Cart::s][1]-gfak*R_temp[Cart::y][Cart::s][2]);
R_temp[Cart::xyy][Cart::s][1]+=pma0*R_temp[Cart::yy][Cart::s][1]+wmp0*R_temp[Cart::yy][Cart::s][2];
R_temp[Cart::yyz][Cart::s][1]+=pma2*R_temp[Cart::yy][Cart::s][1]+wmp2*R_temp[Cart::yy][Cart::s][2];
R_temp[Cart::xxy][Cart::s][1]+=pma1*R_temp[Cart::xx][Cart::s][1]+wmp1*R_temp[Cart::xx][Cart::s][2];
R_temp[Cart::xyz][Cart::s][1]+=pma0*R_temp[Cart::yz][Cart::s][1]+wmp0*R_temp[Cart::yz][Cart::s][2];
R_temp[Cart::yzz][Cart::s][1]+=pma1*R_temp[Cart::zz][Cart::s][1]+wmp1*R_temp[Cart::zz][Cart::s][2];
R_temp[Cart::xxx][Cart::s][1]+=pma0*R_temp[Cart::xx][Cart::s][1]+wmp0*R_temp[Cart::xx][Cart::s][2]+1*rzeta*(R_temp[Cart::x][Cart::s][1]-gfak*R_temp[Cart::x][Cart::s][2]);
R_temp[Cart::xxz][Cart::s][1]+=pma2*R_temp[Cart::xx][Cart::s][1]+wmp2*R_temp[Cart::xx][Cart::s][2];
R_temp[Cart::xzz][Cart::s][1]+=pma0*R_temp[Cart::zz][Cart::s][1]+wmp0*R_temp[Cart::zz][Cart::s][2];
R_temp[Cart::zzz][Cart::s][1]+=pma2*R_temp[Cart::zz][Cart::s][1]+wmp2*R_temp[Cart::zz][Cart::s][2]+1*rzeta*(R_temp[Cart::z][Cart::s][1]-gfak*R_temp[Cart::z][Cart::s][2]);
}
}
//------------------------------------------------------

//Integral f - s - s - m2
if (_mmax >2 ){
if (_lmax_alpha>2){
R_temp[Cart::yyy][Cart::s][2]+=pma1*R_temp[Cart::yy][Cart::s][2]+wmp1*R_temp[Cart::yy][Cart::s][3]+1*rzeta*(R_temp[Cart::y][Cart::s][2]-gfak*R_temp[Cart::y][Cart::s][3]);
R_temp[Cart::xyy][Cart::s][2]+=pma0*R_temp[Cart::yy][Cart::s][2]+wmp0*R_temp[Cart::yy][Cart::s][3];
R_temp[Cart::yyz][Cart::s][2]+=pma2*R_temp[Cart::yy][Cart::s][2]+wmp2*R_temp[Cart::yy][Cart::s][3];
R_temp[Cart::xxy][Cart::s][2]+=pma1*R_temp[Cart::xx][Cart::s][2]+wmp1*R_temp[Cart::xx][Cart::s][3];
R_temp[Cart::xyz][Cart::s][2]+=pma0*R_temp[Cart::yz][Cart::s][2]+wmp0*R_temp[Cart::yz][Cart::s][3];
R_temp[Cart::yzz][Cart::s][2]+=pma1*R_temp[Cart::zz][Cart::s][2]+wmp1*R_temp[Cart::zz][Cart::s][3];
R_temp[Cart::xxx][Cart::s][2]+=pma0*R_temp[Cart::xx][Cart::s][2]+wmp0*R_temp[Cart::xx][Cart::s][3]+1*rzeta*(R_temp[Cart::x][Cart::s][2]-gfak*R_temp[Cart::x][Cart::s][3]);
R_temp[Cart::xxz][Cart::s][2]+=pma2*R_temp[Cart::xx][Cart::s][2]+wmp2*R_temp[Cart::xx][Cart::s][3];
R_temp[Cart::xzz][Cart::s][2]+=pma0*R_temp[Cart::zz][Cart::s][2]+wmp0*R_temp[Cart::zz][Cart::s][3];
R_temp[Cart::zzz][Cart::s][2]+=pma2*R_temp[Cart::zz][Cart::s][2]+wmp2*R_temp[Cart::zz][Cart::s][3]+1*rzeta*(R_temp[Cart::z][Cart::s][2]-gfak*R_temp[Cart::z][Cart::s][3]);
}
}
//------------------------------------------------------

//Integral f - s - s - m3
if (_mmax >3 ){
if (_lmax_alpha>2){
R_temp[Cart::yyy][Cart::s][3]+=pma1*R_temp[Cart::yy][Cart::s][3]+wmp1*R_temp[Cart::yy][Cart::s][4]+1*rzeta*(R_temp[Cart::y][Cart::s][3]-gfak*R_temp[Cart::y][Cart::s][4]);
R_temp[Cart::xyy][Cart::s][3]+=pma0*R_temp[Cart::yy][Cart::s][3]+wmp0*R_temp[Cart::yy][Cart::s][4];
R_temp[Cart::yyz][Cart::s][3]+=pma2*R_temp[Cart::yy][Cart::s][3]+wmp2*R_temp[Cart::yy][Cart::s][4];
R_temp[Cart::xxy][Cart::s][3]+=pma1*R_temp[Cart::xx][Cart::s][3]+wmp1*R_temp[Cart::xx][Cart::s][4];
R_temp[Cart::xyz][Cart::s][3]+=pma0*R_temp[Cart::yz][Cart::s][3]+wmp0*R_temp[Cart::yz][Cart::s][4];
R_temp[Cart::yzz][Cart::s][3]+=pma1*R_temp[Cart::zz][Cart::s][3]+wmp1*R_temp[Cart::zz][Cart::s][4];
R_temp[Cart::xxx][Cart::s][3]+=pma0*R_temp[Cart::xx][Cart::s][3]+wmp0*R_temp[Cart::xx][Cart::s][4]+1*rzeta*(R_temp[Cart::x][Cart::s][3]-gfak*R_temp[Cart::x][Cart::s][4]);
R_temp[Cart::xxz][Cart::s][3]+=pma2*R_temp[Cart::xx][Cart::s][3]+wmp2*R_temp[Cart::xx][Cart::s][4];
R_temp[Cart::xzz][Cart::s][3]+=pma0*R_temp[Cart::zz][Cart::s][3]+wmp0*R_temp[Cart::zz][Cart::s][4];
R_temp[Cart::zzz][Cart::s][3]+=pma2*R_temp[Cart::zz][Cart::s][3]+wmp2*R_temp[Cart::zz][Cart::s][4]+1*rzeta*(R_temp[Cart::z][Cart::s][3]-gfak*R_temp[Cart::z][Cart::s][4]);
}
}
//------------------------------------------------------

//Integral f - s - s - m4
if (_mmax >4 ){
if (_lmax_alpha>2){
R_temp[Cart::yyy][Cart::s][4]+=pma1*R_temp[Cart::yy][Cart::s][4]+wmp1*R_temp[Cart::yy][Cart::s][5]+1*rzeta*(R_temp[Cart::y][Cart::s][4]-gfak*R_temp[Cart::y][Cart::s][5]);
R_temp[Cart::xyy][Cart::s][4]+=pma0*R_temp[Cart::yy][Cart::s][4]+wmp0*R_temp[Cart::yy][Cart::s][5];
R_temp[Cart::yyz][Cart::s][4]+=pma2*R_temp[Cart::yy][Cart::s][4]+wmp2*R_temp[Cart::yy][Cart::s][5];
R_temp[Cart::xxy][Cart::s][4]+=pma1*R_temp[Cart::xx][Cart::s][4]+wmp1*R_temp[Cart::xx][Cart::s][5];
R_temp[Cart::xyz][Cart::s][4]+=pma0*R_temp[Cart::yz][Cart::s][4]+wmp0*R_temp[Cart::yz][Cart::s][5];
R_temp[Cart::yzz][Cart::s][4]+=pma1*R_temp[Cart::zz][Cart::s][4]+wmp1*R_temp[Cart::zz][Cart::s][5];
R_temp[Cart::xxx][Cart::s][4]+=pma0*R_temp[Cart::xx][Cart::s][4]+wmp0*R_temp[Cart::xx][Cart::s][5]+1*rzeta*(R_temp[Cart::x][Cart::s][4]-gfak*R_temp[Cart::x][Cart::s][5]);
R_temp[Cart::xxz][Cart::s][4]+=pma2*R_temp[Cart::xx][Cart::s][4]+wmp2*R_temp[Cart::xx][Cart::s][5];
R_temp[Cart::xzz][Cart::s][4]+=pma0*R_temp[Cart::zz][Cart::s][4]+wmp0*R_temp[Cart::zz][Cart::s][5];
R_temp[Cart::zzz][Cart::s][4]+=pma2*R_temp[Cart::zz][Cart::s][4]+wmp2*R_temp[Cart::zz][Cart::s][5]+1*rzeta*(R_temp[Cart::z][Cart::s][4]-gfak*R_temp[Cart::z][Cart::s][5]);
}
}
//------------------------------------------------------

//Integral f - s - p - m0
if (_lmax_alpha>2 && _lmax_gamma>0){
R_temp[Cart::yyy][Cart::y][0]+=pma1*R_temp[Cart::yy][Cart::y][0]+wmp1*R_temp[Cart::yy][Cart::y][1]+1*rzeta*(R_temp[Cart::y][Cart::y][0]-gfak*R_temp[Cart::y][Cart::y][1])+0.5/_decay*1*R_temp[Cart::yy][Cart::s][1];
R_temp[Cart::yyy][Cart::x][0]+=pma1*R_temp[Cart::yy][Cart::x][0]+wmp1*R_temp[Cart::yy][Cart::x][1]+1*rzeta*(R_temp[Cart::y][Cart::x][0]-gfak*R_temp[Cart::y][Cart::x][1]);
R_temp[Cart::yyy][Cart::z][0]+=pma1*R_temp[Cart::yy][Cart::z][0]+wmp1*R_temp[Cart::yy][Cart::z][1]+1*rzeta*(R_temp[Cart::y][Cart::z][0]-gfak*R_temp[Cart::y][Cart::z][1]);
R_temp[Cart::xyy][Cart::y][0]+=pma0*R_temp[Cart::yy][Cart::y][0]+wmp0*R_temp[Cart::yy][Cart::y][1];
R_temp[Cart::xyy][Cart::x][0]+=pma0*R_temp[Cart::yy][Cart::x][0]+wmp0*R_temp[Cart::yy][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][1];
R_temp[Cart::xyy][Cart::z][0]+=pma0*R_temp[Cart::yy][Cart::z][0]+wmp0*R_temp[Cart::yy][Cart::z][1];
R_temp[Cart::yyz][Cart::y][0]+=pma2*R_temp[Cart::yy][Cart::y][0]+wmp2*R_temp[Cart::yy][Cart::y][1];
R_temp[Cart::yyz][Cart::x][0]+=pma2*R_temp[Cart::yy][Cart::x][0]+wmp2*R_temp[Cart::yy][Cart::x][1];
R_temp[Cart::yyz][Cart::z][0]+=pma2*R_temp[Cart::yy][Cart::z][0]+wmp2*R_temp[Cart::yy][Cart::z][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][1];
R_temp[Cart::xxy][Cart::y][0]+=pma1*R_temp[Cart::xx][Cart::y][0]+wmp1*R_temp[Cart::xx][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][1];
R_temp[Cart::xxy][Cart::x][0]+=pma1*R_temp[Cart::xx][Cart::x][0]+wmp1*R_temp[Cart::xx][Cart::x][1];
R_temp[Cart::xxy][Cart::z][0]+=pma1*R_temp[Cart::xx][Cart::z][0]+wmp1*R_temp[Cart::xx][Cart::z][1];
R_temp[Cart::xyz][Cart::y][0]+=pma0*R_temp[Cart::yz][Cart::y][0]+wmp0*R_temp[Cart::yz][Cart::y][1];
R_temp[Cart::xyz][Cart::x][0]+=pma0*R_temp[Cart::yz][Cart::x][0]+wmp0*R_temp[Cart::yz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yz][Cart::s][1];
R_temp[Cart::xyz][Cart::z][0]+=pma0*R_temp[Cart::yz][Cart::z][0]+wmp0*R_temp[Cart::yz][Cart::z][1];
R_temp[Cart::yzz][Cart::y][0]+=pma1*R_temp[Cart::zz][Cart::y][0]+wmp1*R_temp[Cart::zz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][1];
R_temp[Cart::yzz][Cart::x][0]+=pma1*R_temp[Cart::zz][Cart::x][0]+wmp1*R_temp[Cart::zz][Cart::x][1];
R_temp[Cart::yzz][Cart::z][0]+=pma1*R_temp[Cart::zz][Cart::z][0]+wmp1*R_temp[Cart::zz][Cart::z][1];
R_temp[Cart::xxx][Cart::y][0]+=pma0*R_temp[Cart::xx][Cart::y][0]+wmp0*R_temp[Cart::xx][Cart::y][1]+1*rzeta*(R_temp[Cart::x][Cart::y][0]-gfak*R_temp[Cart::x][Cart::y][1]);
R_temp[Cart::xxx][Cart::x][0]+=pma0*R_temp[Cart::xx][Cart::x][0]+wmp0*R_temp[Cart::xx][Cart::x][1]+1*rzeta*(R_temp[Cart::x][Cart::x][0]-gfak*R_temp[Cart::x][Cart::x][1])+0.5/_decay*1*R_temp[Cart::xx][Cart::s][1];
R_temp[Cart::xxx][Cart::z][0]+=pma0*R_temp[Cart::xx][Cart::z][0]+wmp0*R_temp[Cart::xx][Cart::z][1]+1*rzeta*(R_temp[Cart::x][Cart::z][0]-gfak*R_temp[Cart::x][Cart::z][1]);
R_temp[Cart::xxz][Cart::y][0]+=pma2*R_temp[Cart::xx][Cart::y][0]+wmp2*R_temp[Cart::xx][Cart::y][1];
R_temp[Cart::xxz][Cart::x][0]+=pma2*R_temp[Cart::xx][Cart::x][0]+wmp2*R_temp[Cart::xx][Cart::x][1];
R_temp[Cart::xxz][Cart::z][0]+=pma2*R_temp[Cart::xx][Cart::z][0]+wmp2*R_temp[Cart::xx][Cart::z][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][1];
R_temp[Cart::xzz][Cart::y][0]+=pma0*R_temp[Cart::zz][Cart::y][0]+wmp0*R_temp[Cart::zz][Cart::y][1];
R_temp[Cart::xzz][Cart::x][0]+=pma0*R_temp[Cart::zz][Cart::x][0]+wmp0*R_temp[Cart::zz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][1];
R_temp[Cart::xzz][Cart::z][0]+=pma0*R_temp[Cart::zz][Cart::z][0]+wmp0*R_temp[Cart::zz][Cart::z][1];
R_temp[Cart::zzz][Cart::y][0]+=pma2*R_temp[Cart::zz][Cart::y][0]+wmp2*R_temp[Cart::zz][Cart::y][1]+1*rzeta*(R_temp[Cart::z][Cart::y][0]-gfak*R_temp[Cart::z][Cart::y][1]);
R_temp[Cart::zzz][Cart::x][0]+=pma2*R_temp[Cart::zz][Cart::x][0]+wmp2*R_temp[Cart::zz][Cart::x][1]+1*rzeta*(R_temp[Cart::z][Cart::x][0]-gfak*R_temp[Cart::z][Cart::x][1]);
R_temp[Cart::zzz][Cart::z][0]+=pma2*R_temp[Cart::zz][Cart::z][0]+wmp2*R_temp[Cart::zz][Cart::z][1]+1*rzeta*(R_temp[Cart::z][Cart::z][0]-gfak*R_temp[Cart::z][Cart::z][1])+0.5/_decay*1*R_temp[Cart::zz][Cart::s][1];
}
//------------------------------------------------------

//Integral f - s - p - m1
if (_mmax >1 ){
if (_lmax_alpha>2 && _lmax_gamma>0){
R_temp[Cart::yyy][Cart::y][1]+=pma1*R_temp[Cart::yy][Cart::y][1]+wmp1*R_temp[Cart::yy][Cart::y][2]+1*rzeta*(R_temp[Cart::y][Cart::y][1]-gfak*R_temp[Cart::y][Cart::y][2])+0.5/_decay*1*R_temp[Cart::yy][Cart::s][2];
R_temp[Cart::yyy][Cart::x][1]+=pma1*R_temp[Cart::yy][Cart::x][1]+wmp1*R_temp[Cart::yy][Cart::x][2]+1*rzeta*(R_temp[Cart::y][Cart::x][1]-gfak*R_temp[Cart::y][Cart::x][2]);
R_temp[Cart::yyy][Cart::z][1]+=pma1*R_temp[Cart::yy][Cart::z][1]+wmp1*R_temp[Cart::yy][Cart::z][2]+1*rzeta*(R_temp[Cart::y][Cart::z][1]-gfak*R_temp[Cart::y][Cart::z][2]);
R_temp[Cart::xyy][Cart::y][1]+=pma0*R_temp[Cart::yy][Cart::y][1]+wmp0*R_temp[Cart::yy][Cart::y][2];
R_temp[Cart::xyy][Cart::x][1]+=pma0*R_temp[Cart::yy][Cart::x][1]+wmp0*R_temp[Cart::yy][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][2];
R_temp[Cart::xyy][Cart::z][1]+=pma0*R_temp[Cart::yy][Cart::z][1]+wmp0*R_temp[Cart::yy][Cart::z][2];
R_temp[Cart::yyz][Cart::y][1]+=pma2*R_temp[Cart::yy][Cart::y][1]+wmp2*R_temp[Cart::yy][Cart::y][2];
R_temp[Cart::yyz][Cart::x][1]+=pma2*R_temp[Cart::yy][Cart::x][1]+wmp2*R_temp[Cart::yy][Cart::x][2];
R_temp[Cart::yyz][Cart::z][1]+=pma2*R_temp[Cart::yy][Cart::z][1]+wmp2*R_temp[Cart::yy][Cart::z][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][2];
R_temp[Cart::xxy][Cart::y][1]+=pma1*R_temp[Cart::xx][Cart::y][1]+wmp1*R_temp[Cart::xx][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][2];
R_temp[Cart::xxy][Cart::x][1]+=pma1*R_temp[Cart::xx][Cart::x][1]+wmp1*R_temp[Cart::xx][Cart::x][2];
R_temp[Cart::xxy][Cart::z][1]+=pma1*R_temp[Cart::xx][Cart::z][1]+wmp1*R_temp[Cart::xx][Cart::z][2];
R_temp[Cart::xyz][Cart::y][1]+=pma0*R_temp[Cart::yz][Cart::y][1]+wmp0*R_temp[Cart::yz][Cart::y][2];
R_temp[Cart::xyz][Cart::x][1]+=pma0*R_temp[Cart::yz][Cart::x][1]+wmp0*R_temp[Cart::yz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yz][Cart::s][2];
R_temp[Cart::xyz][Cart::z][1]+=pma0*R_temp[Cart::yz][Cart::z][1]+wmp0*R_temp[Cart::yz][Cart::z][2];
R_temp[Cart::yzz][Cart::y][1]+=pma1*R_temp[Cart::zz][Cart::y][1]+wmp1*R_temp[Cart::zz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][2];
R_temp[Cart::yzz][Cart::x][1]+=pma1*R_temp[Cart::zz][Cart::x][1]+wmp1*R_temp[Cart::zz][Cart::x][2];
R_temp[Cart::yzz][Cart::z][1]+=pma1*R_temp[Cart::zz][Cart::z][1]+wmp1*R_temp[Cart::zz][Cart::z][2];
R_temp[Cart::xxx][Cart::y][1]+=pma0*R_temp[Cart::xx][Cart::y][1]+wmp0*R_temp[Cart::xx][Cart::y][2]+1*rzeta*(R_temp[Cart::x][Cart::y][1]-gfak*R_temp[Cart::x][Cart::y][2]);
R_temp[Cart::xxx][Cart::x][1]+=pma0*R_temp[Cart::xx][Cart::x][1]+wmp0*R_temp[Cart::xx][Cart::x][2]+1*rzeta*(R_temp[Cart::x][Cart::x][1]-gfak*R_temp[Cart::x][Cart::x][2])+0.5/_decay*1*R_temp[Cart::xx][Cart::s][2];
R_temp[Cart::xxx][Cart::z][1]+=pma0*R_temp[Cart::xx][Cart::z][1]+wmp0*R_temp[Cart::xx][Cart::z][2]+1*rzeta*(R_temp[Cart::x][Cart::z][1]-gfak*R_temp[Cart::x][Cart::z][2]);
R_temp[Cart::xxz][Cart::y][1]+=pma2*R_temp[Cart::xx][Cart::y][1]+wmp2*R_temp[Cart::xx][Cart::y][2];
R_temp[Cart::xxz][Cart::x][1]+=pma2*R_temp[Cart::xx][Cart::x][1]+wmp2*R_temp[Cart::xx][Cart::x][2];
R_temp[Cart::xxz][Cart::z][1]+=pma2*R_temp[Cart::xx][Cart::z][1]+wmp2*R_temp[Cart::xx][Cart::z][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][2];
R_temp[Cart::xzz][Cart::y][1]+=pma0*R_temp[Cart::zz][Cart::y][1]+wmp0*R_temp[Cart::zz][Cart::y][2];
R_temp[Cart::xzz][Cart::x][1]+=pma0*R_temp[Cart::zz][Cart::x][1]+wmp0*R_temp[Cart::zz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][2];
R_temp[Cart::xzz][Cart::z][1]+=pma0*R_temp[Cart::zz][Cart::z][1]+wmp0*R_temp[Cart::zz][Cart::z][2];
R_temp[Cart::zzz][Cart::y][1]+=pma2*R_temp[Cart::zz][Cart::y][1]+wmp2*R_temp[Cart::zz][Cart::y][2]+1*rzeta*(R_temp[Cart::z][Cart::y][1]-gfak*R_temp[Cart::z][Cart::y][2]);
R_temp[Cart::zzz][Cart::x][1]+=pma2*R_temp[Cart::zz][Cart::x][1]+wmp2*R_temp[Cart::zz][Cart::x][2]+1*rzeta*(R_temp[Cart::z][Cart::x][1]-gfak*R_temp[Cart::z][Cart::x][2]);
R_temp[Cart::zzz][Cart::z][1]+=pma2*R_temp[Cart::zz][Cart::z][1]+wmp2*R_temp[Cart::zz][Cart::z][2]+1*rzeta*(R_temp[Cart::z][Cart::z][1]-gfak*R_temp[Cart::z][Cart::z][2])+0.5/_decay*1*R_temp[Cart::zz][Cart::s][2];
}
}
//------------------------------------------------------

//Integral f - s - p - m2
if (_mmax >2 ){
if (_lmax_alpha>2 && _lmax_gamma>0){
R_temp[Cart::yyy][Cart::y][2]+=pma1*R_temp[Cart::yy][Cart::y][2]+wmp1*R_temp[Cart::yy][Cart::y][3]+1*rzeta*(R_temp[Cart::y][Cart::y][2]-gfak*R_temp[Cart::y][Cart::y][3])+0.5/_decay*1*R_temp[Cart::yy][Cart::s][3];
R_temp[Cart::yyy][Cart::x][2]+=pma1*R_temp[Cart::yy][Cart::x][2]+wmp1*R_temp[Cart::yy][Cart::x][3]+1*rzeta*(R_temp[Cart::y][Cart::x][2]-gfak*R_temp[Cart::y][Cart::x][3]);
R_temp[Cart::yyy][Cart::z][2]+=pma1*R_temp[Cart::yy][Cart::z][2]+wmp1*R_temp[Cart::yy][Cart::z][3]+1*rzeta*(R_temp[Cart::y][Cart::z][2]-gfak*R_temp[Cart::y][Cart::z][3]);
R_temp[Cart::xyy][Cart::y][2]+=pma0*R_temp[Cart::yy][Cart::y][2]+wmp0*R_temp[Cart::yy][Cart::y][3];
R_temp[Cart::xyy][Cart::x][2]+=pma0*R_temp[Cart::yy][Cart::x][2]+wmp0*R_temp[Cart::yy][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][3];
R_temp[Cart::xyy][Cart::z][2]+=pma0*R_temp[Cart::yy][Cart::z][2]+wmp0*R_temp[Cart::yy][Cart::z][3];
R_temp[Cart::yyz][Cart::y][2]+=pma2*R_temp[Cart::yy][Cart::y][2]+wmp2*R_temp[Cart::yy][Cart::y][3];
R_temp[Cart::yyz][Cart::x][2]+=pma2*R_temp[Cart::yy][Cart::x][2]+wmp2*R_temp[Cart::yy][Cart::x][3];
R_temp[Cart::yyz][Cart::z][2]+=pma2*R_temp[Cart::yy][Cart::z][2]+wmp2*R_temp[Cart::yy][Cart::z][3]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][3];
R_temp[Cart::xxy][Cart::y][2]+=pma1*R_temp[Cart::xx][Cart::y][2]+wmp1*R_temp[Cart::xx][Cart::y][3]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][3];
R_temp[Cart::xxy][Cart::x][2]+=pma1*R_temp[Cart::xx][Cart::x][2]+wmp1*R_temp[Cart::xx][Cart::x][3];
R_temp[Cart::xxy][Cart::z][2]+=pma1*R_temp[Cart::xx][Cart::z][2]+wmp1*R_temp[Cart::xx][Cart::z][3];
R_temp[Cart::xyz][Cart::y][2]+=pma0*R_temp[Cart::yz][Cart::y][2]+wmp0*R_temp[Cart::yz][Cart::y][3];
R_temp[Cart::xyz][Cart::x][2]+=pma0*R_temp[Cart::yz][Cart::x][2]+wmp0*R_temp[Cart::yz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yz][Cart::s][3];
R_temp[Cart::xyz][Cart::z][2]+=pma0*R_temp[Cart::yz][Cart::z][2]+wmp0*R_temp[Cart::yz][Cart::z][3];
R_temp[Cart::yzz][Cart::y][2]+=pma1*R_temp[Cart::zz][Cart::y][2]+wmp1*R_temp[Cart::zz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][3];
R_temp[Cart::yzz][Cart::x][2]+=pma1*R_temp[Cart::zz][Cart::x][2]+wmp1*R_temp[Cart::zz][Cart::x][3];
R_temp[Cart::yzz][Cart::z][2]+=pma1*R_temp[Cart::zz][Cart::z][2]+wmp1*R_temp[Cart::zz][Cart::z][3];
R_temp[Cart::xxx][Cart::y][2]+=pma0*R_temp[Cart::xx][Cart::y][2]+wmp0*R_temp[Cart::xx][Cart::y][3]+1*rzeta*(R_temp[Cart::x][Cart::y][2]-gfak*R_temp[Cart::x][Cart::y][3]);
R_temp[Cart::xxx][Cart::x][2]+=pma0*R_temp[Cart::xx][Cart::x][2]+wmp0*R_temp[Cart::xx][Cart::x][3]+1*rzeta*(R_temp[Cart::x][Cart::x][2]-gfak*R_temp[Cart::x][Cart::x][3])+0.5/_decay*1*R_temp[Cart::xx][Cart::s][3];
R_temp[Cart::xxx][Cart::z][2]+=pma0*R_temp[Cart::xx][Cart::z][2]+wmp0*R_temp[Cart::xx][Cart::z][3]+1*rzeta*(R_temp[Cart::x][Cart::z][2]-gfak*R_temp[Cart::x][Cart::z][3]);
R_temp[Cart::xxz][Cart::y][2]+=pma2*R_temp[Cart::xx][Cart::y][2]+wmp2*R_temp[Cart::xx][Cart::y][3];
R_temp[Cart::xxz][Cart::x][2]+=pma2*R_temp[Cart::xx][Cart::x][2]+wmp2*R_temp[Cart::xx][Cart::x][3];
R_temp[Cart::xxz][Cart::z][2]+=pma2*R_temp[Cart::xx][Cart::z][2]+wmp2*R_temp[Cart::xx][Cart::z][3]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][3];
R_temp[Cart::xzz][Cart::y][2]+=pma0*R_temp[Cart::zz][Cart::y][2]+wmp0*R_temp[Cart::zz][Cart::y][3];
R_temp[Cart::xzz][Cart::x][2]+=pma0*R_temp[Cart::zz][Cart::x][2]+wmp0*R_temp[Cart::zz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][3];
R_temp[Cart::xzz][Cart::z][2]+=pma0*R_temp[Cart::zz][Cart::z][2]+wmp0*R_temp[Cart::zz][Cart::z][3];
R_temp[Cart::zzz][Cart::y][2]+=pma2*R_temp[Cart::zz][Cart::y][2]+wmp2*R_temp[Cart::zz][Cart::y][3]+1*rzeta*(R_temp[Cart::z][Cart::y][2]-gfak*R_temp[Cart::z][Cart::y][3]);
R_temp[Cart::zzz][Cart::x][2]+=pma2*R_temp[Cart::zz][Cart::x][2]+wmp2*R_temp[Cart::zz][Cart::x][3]+1*rzeta*(R_temp[Cart::z][Cart::x][2]-gfak*R_temp[Cart::z][Cart::x][3]);
R_temp[Cart::zzz][Cart::z][2]+=pma2*R_temp[Cart::zz][Cart::z][2]+wmp2*R_temp[Cart::zz][Cart::z][3]+1*rzeta*(R_temp[Cart::z][Cart::z][2]-gfak*R_temp[Cart::z][Cart::z][3])+0.5/_decay*1*R_temp[Cart::zz][Cart::s][3];
}
}
//------------------------------------------------------

//Integral f - s - p - m3
if (_mmax >3 ){
if (_lmax_alpha>2 && _lmax_gamma>0){
R_temp[Cart::yyy][Cart::y][3]+=pma1*R_temp[Cart::yy][Cart::y][3]+wmp1*R_temp[Cart::yy][Cart::y][4]+1*rzeta*(R_temp[Cart::y][Cart::y][3]-gfak*R_temp[Cart::y][Cart::y][4])+0.5/_decay*1*R_temp[Cart::yy][Cart::s][4];
R_temp[Cart::yyy][Cart::x][3]+=pma1*R_temp[Cart::yy][Cart::x][3]+wmp1*R_temp[Cart::yy][Cart::x][4]+1*rzeta*(R_temp[Cart::y][Cart::x][3]-gfak*R_temp[Cart::y][Cart::x][4]);
R_temp[Cart::yyy][Cart::z][3]+=pma1*R_temp[Cart::yy][Cart::z][3]+wmp1*R_temp[Cart::yy][Cart::z][4]+1*rzeta*(R_temp[Cart::y][Cart::z][3]-gfak*R_temp[Cart::y][Cart::z][4]);
R_temp[Cart::xyy][Cart::y][3]+=pma0*R_temp[Cart::yy][Cart::y][3]+wmp0*R_temp[Cart::yy][Cart::y][4];
R_temp[Cart::xyy][Cart::x][3]+=pma0*R_temp[Cart::yy][Cart::x][3]+wmp0*R_temp[Cart::yy][Cart::x][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][4];
R_temp[Cart::xyy][Cart::z][3]+=pma0*R_temp[Cart::yy][Cart::z][3]+wmp0*R_temp[Cart::yy][Cart::z][4];
R_temp[Cart::yyz][Cart::y][3]+=pma2*R_temp[Cart::yy][Cart::y][3]+wmp2*R_temp[Cart::yy][Cart::y][4];
R_temp[Cart::yyz][Cart::x][3]+=pma2*R_temp[Cart::yy][Cart::x][3]+wmp2*R_temp[Cart::yy][Cart::x][4];
R_temp[Cart::yyz][Cart::z][3]+=pma2*R_temp[Cart::yy][Cart::z][3]+wmp2*R_temp[Cart::yy][Cart::z][4]+0.5/_decay*1*R_temp[Cart::yy][Cart::s][4];
R_temp[Cart::xxy][Cart::y][3]+=pma1*R_temp[Cart::xx][Cart::y][3]+wmp1*R_temp[Cart::xx][Cart::y][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][4];
R_temp[Cart::xxy][Cart::x][3]+=pma1*R_temp[Cart::xx][Cart::x][3]+wmp1*R_temp[Cart::xx][Cart::x][4];
R_temp[Cart::xxy][Cart::z][3]+=pma1*R_temp[Cart::xx][Cart::z][3]+wmp1*R_temp[Cart::xx][Cart::z][4];
R_temp[Cart::xyz][Cart::y][3]+=pma0*R_temp[Cart::yz][Cart::y][3]+wmp0*R_temp[Cart::yz][Cart::y][4];
R_temp[Cart::xyz][Cart::x][3]+=pma0*R_temp[Cart::yz][Cart::x][3]+wmp0*R_temp[Cart::yz][Cart::x][4]+0.5/_decay*1*R_temp[Cart::yz][Cart::s][4];
R_temp[Cart::xyz][Cart::z][3]+=pma0*R_temp[Cart::yz][Cart::z][3]+wmp0*R_temp[Cart::yz][Cart::z][4];
R_temp[Cart::yzz][Cart::y][3]+=pma1*R_temp[Cart::zz][Cart::y][3]+wmp1*R_temp[Cart::zz][Cart::y][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][4];
R_temp[Cart::yzz][Cart::x][3]+=pma1*R_temp[Cart::zz][Cart::x][3]+wmp1*R_temp[Cart::zz][Cart::x][4];
R_temp[Cart::yzz][Cart::z][3]+=pma1*R_temp[Cart::zz][Cart::z][3]+wmp1*R_temp[Cart::zz][Cart::z][4];
R_temp[Cart::xxx][Cart::y][3]+=pma0*R_temp[Cart::xx][Cart::y][3]+wmp0*R_temp[Cart::xx][Cart::y][4]+1*rzeta*(R_temp[Cart::x][Cart::y][3]-gfak*R_temp[Cart::x][Cart::y][4]);
R_temp[Cart::xxx][Cart::x][3]+=pma0*R_temp[Cart::xx][Cart::x][3]+wmp0*R_temp[Cart::xx][Cart::x][4]+1*rzeta*(R_temp[Cart::x][Cart::x][3]-gfak*R_temp[Cart::x][Cart::x][4])+0.5/_decay*1*R_temp[Cart::xx][Cart::s][4];
R_temp[Cart::xxx][Cart::z][3]+=pma0*R_temp[Cart::xx][Cart::z][3]+wmp0*R_temp[Cart::xx][Cart::z][4]+1*rzeta*(R_temp[Cart::x][Cart::z][3]-gfak*R_temp[Cart::x][Cart::z][4]);
R_temp[Cart::xxz][Cart::y][3]+=pma2*R_temp[Cart::xx][Cart::y][3]+wmp2*R_temp[Cart::xx][Cart::y][4];
R_temp[Cart::xxz][Cart::x][3]+=pma2*R_temp[Cart::xx][Cart::x][3]+wmp2*R_temp[Cart::xx][Cart::x][4];
R_temp[Cart::xxz][Cart::z][3]+=pma2*R_temp[Cart::xx][Cart::z][3]+wmp2*R_temp[Cart::xx][Cart::z][4]+0.5/_decay*1*R_temp[Cart::xx][Cart::s][4];
R_temp[Cart::xzz][Cart::y][3]+=pma0*R_temp[Cart::zz][Cart::y][3]+wmp0*R_temp[Cart::zz][Cart::y][4];
R_temp[Cart::xzz][Cart::x][3]+=pma0*R_temp[Cart::zz][Cart::x][3]+wmp0*R_temp[Cart::zz][Cart::x][4]+0.5/_decay*1*R_temp[Cart::zz][Cart::s][4];
R_temp[Cart::xzz][Cart::z][3]+=pma0*R_temp[Cart::zz][Cart::z][3]+wmp0*R_temp[Cart::zz][Cart::z][4];
R_temp[Cart::zzz][Cart::y][3]+=pma2*R_temp[Cart::zz][Cart::y][3]+wmp2*R_temp[Cart::zz][Cart::y][4]+1*rzeta*(R_temp[Cart::z][Cart::y][3]-gfak*R_temp[Cart::z][Cart::y][4]);
R_temp[Cart::zzz][Cart::x][3]+=pma2*R_temp[Cart::zz][Cart::x][3]+wmp2*R_temp[Cart::zz][Cart::x][4]+1*rzeta*(R_temp[Cart::z][Cart::x][3]-gfak*R_temp[Cart::z][Cart::x][4]);
R_temp[Cart::zzz][Cart::z][3]+=pma2*R_temp[Cart::zz][Cart::z][3]+wmp2*R_temp[Cart::zz][Cart::z][4]+1*rzeta*(R_temp[Cart::z][Cart::z][3]-gfak*R_temp[Cart::z][Cart::z][4])+0.5/_decay*1*R_temp[Cart::zz][Cart::s][4];
}
}
//------------------------------------------------------

//Integral f - s - d - m0
if (_lmax_alpha>2 && _lmax_gamma>1){
R_temp[Cart::yyy][Cart::yy][0]+=pma1*R_temp[Cart::yy][Cart::yy][0]+wmp1*R_temp[Cart::yy][Cart::yy][1]+1*rzeta*(R_temp[Cart::y][Cart::yy][0]-gfak*R_temp[Cart::y][Cart::yy][1])+0.5/_decay*2*R_temp[Cart::yy][Cart::y][1];
R_temp[Cart::yyy][Cart::xy][0]+=pma1*R_temp[Cart::yy][Cart::xy][0]+wmp1*R_temp[Cart::yy][Cart::xy][1]+1*rzeta*(R_temp[Cart::y][Cart::xy][0]-gfak*R_temp[Cart::y][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::yy][Cart::x][1];
R_temp[Cart::yyy][Cart::yz][0]+=pma1*R_temp[Cart::yy][Cart::yz][0]+wmp1*R_temp[Cart::yy][Cart::yz][1]+1*rzeta*(R_temp[Cart::y][Cart::yz][0]-gfak*R_temp[Cart::y][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::yy][Cart::z][1];
R_temp[Cart::yyy][Cart::xx][0]+=pma1*R_temp[Cart::yy][Cart::xx][0]+wmp1*R_temp[Cart::yy][Cart::xx][1]+1*rzeta*(R_temp[Cart::y][Cart::xx][0]-gfak*R_temp[Cart::y][Cart::xx][1]);
R_temp[Cart::yyy][Cart::xz][0]+=pma1*R_temp[Cart::yy][Cart::xz][0]+wmp1*R_temp[Cart::yy][Cart::xz][1]+1*rzeta*(R_temp[Cart::y][Cart::xz][0]-gfak*R_temp[Cart::y][Cart::xz][1]);
R_temp[Cart::yyy][Cart::zz][0]+=pma1*R_temp[Cart::yy][Cart::zz][0]+wmp1*R_temp[Cart::yy][Cart::zz][1]+1*rzeta*(R_temp[Cart::y][Cart::zz][0]-gfak*R_temp[Cart::y][Cart::zz][1]);
R_temp[Cart::xyy][Cart::yy][0]+=pma0*R_temp[Cart::yy][Cart::yy][0]+wmp0*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::xyy][Cart::xy][0]+=pma0*R_temp[Cart::yy][Cart::xy][0]+wmp0*R_temp[Cart::yy][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][1];
R_temp[Cart::xyy][Cart::yz][0]+=pma0*R_temp[Cart::yy][Cart::yz][0]+wmp0*R_temp[Cart::yy][Cart::yz][1];
R_temp[Cart::xyy][Cart::xx][0]+=pma0*R_temp[Cart::yy][Cart::xx][0]+wmp0*R_temp[Cart::yy][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yy][Cart::x][1];
R_temp[Cart::xyy][Cart::xz][0]+=pma0*R_temp[Cart::yy][Cart::xz][0]+wmp0*R_temp[Cart::yy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::z][1];
R_temp[Cart::xyy][Cart::zz][0]+=pma0*R_temp[Cart::yy][Cart::zz][0]+wmp0*R_temp[Cart::yy][Cart::zz][1];
R_temp[Cart::yyz][Cart::yy][0]+=pma2*R_temp[Cart::yy][Cart::yy][0]+wmp2*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::yyz][Cart::xy][0]+=pma2*R_temp[Cart::yy][Cart::xy][0]+wmp2*R_temp[Cart::yy][Cart::xy][1];
R_temp[Cart::yyz][Cart::yz][0]+=pma2*R_temp[Cart::yy][Cart::yz][0]+wmp2*R_temp[Cart::yy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][1];
R_temp[Cart::yyz][Cart::xx][0]+=pma2*R_temp[Cart::yy][Cart::xx][0]+wmp2*R_temp[Cart::yy][Cart::xx][1];
R_temp[Cart::yyz][Cart::xz][0]+=pma2*R_temp[Cart::yy][Cart::xz][0]+wmp2*R_temp[Cart::yy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::x][1];
R_temp[Cart::yyz][Cart::zz][0]+=pma2*R_temp[Cart::yy][Cart::zz][0]+wmp2*R_temp[Cart::yy][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::yy][Cart::z][1];
R_temp[Cart::xxy][Cart::yy][0]+=pma1*R_temp[Cart::xx][Cart::yy][0]+wmp1*R_temp[Cart::xx][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xx][Cart::y][1];
R_temp[Cart::xxy][Cart::xy][0]+=pma1*R_temp[Cart::xx][Cart::xy][0]+wmp1*R_temp[Cart::xx][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][1];
R_temp[Cart::xxy][Cart::yz][0]+=pma1*R_temp[Cart::xx][Cart::yz][0]+wmp1*R_temp[Cart::xx][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::z][1];
R_temp[Cart::xxy][Cart::xx][0]+=pma1*R_temp[Cart::xx][Cart::xx][0]+wmp1*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xxy][Cart::xz][0]+=pma1*R_temp[Cart::xx][Cart::xz][0]+wmp1*R_temp[Cart::xx][Cart::xz][1];
R_temp[Cart::xxy][Cart::zz][0]+=pma1*R_temp[Cart::xx][Cart::zz][0]+wmp1*R_temp[Cart::xx][Cart::zz][1];
R_temp[Cart::xyz][Cart::yy][0]+=pma0*R_temp[Cart::yz][Cart::yy][0]+wmp0*R_temp[Cart::yz][Cart::yy][1];
R_temp[Cart::xyz][Cart::xy][0]+=pma0*R_temp[Cart::yz][Cart::xy][0]+wmp0*R_temp[Cart::yz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yz][Cart::y][1];
R_temp[Cart::xyz][Cart::yz][0]+=pma0*R_temp[Cart::yz][Cart::yz][0]+wmp0*R_temp[Cart::yz][Cart::yz][1];
R_temp[Cart::xyz][Cart::xx][0]+=pma0*R_temp[Cart::yz][Cart::xx][0]+wmp0*R_temp[Cart::yz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yz][Cart::x][1];
R_temp[Cart::xyz][Cart::xz][0]+=pma0*R_temp[Cart::yz][Cart::xz][0]+wmp0*R_temp[Cart::yz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yz][Cart::z][1];
R_temp[Cart::xyz][Cart::zz][0]+=pma0*R_temp[Cart::yz][Cart::zz][0]+wmp0*R_temp[Cart::yz][Cart::zz][1];
R_temp[Cart::yzz][Cart::yy][0]+=pma1*R_temp[Cart::zz][Cart::yy][0]+wmp1*R_temp[Cart::zz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::zz][Cart::y][1];
R_temp[Cart::yzz][Cart::xy][0]+=pma1*R_temp[Cart::zz][Cart::xy][0]+wmp1*R_temp[Cart::zz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::x][1];
R_temp[Cart::yzz][Cart::yz][0]+=pma1*R_temp[Cart::zz][Cart::yz][0]+wmp1*R_temp[Cart::zz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][1];
R_temp[Cart::yzz][Cart::xx][0]+=pma1*R_temp[Cart::zz][Cart::xx][0]+wmp1*R_temp[Cart::zz][Cart::xx][1];
R_temp[Cart::yzz][Cart::xz][0]+=pma1*R_temp[Cart::zz][Cart::xz][0]+wmp1*R_temp[Cart::zz][Cart::xz][1];
R_temp[Cart::yzz][Cart::zz][0]+=pma1*R_temp[Cart::zz][Cart::zz][0]+wmp1*R_temp[Cart::zz][Cart::zz][1];
R_temp[Cart::xxx][Cart::yy][0]+=pma0*R_temp[Cart::xx][Cart::yy][0]+wmp0*R_temp[Cart::xx][Cart::yy][1]+1*rzeta*(R_temp[Cart::x][Cart::yy][0]-gfak*R_temp[Cart::x][Cart::yy][1]);
R_temp[Cart::xxx][Cart::xy][0]+=pma0*R_temp[Cart::xx][Cart::xy][0]+wmp0*R_temp[Cart::xx][Cart::xy][1]+1*rzeta*(R_temp[Cart::x][Cart::xy][0]-gfak*R_temp[Cart::x][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::xx][Cart::y][1];
R_temp[Cart::xxx][Cart::yz][0]+=pma0*R_temp[Cart::xx][Cart::yz][0]+wmp0*R_temp[Cart::xx][Cart::yz][1]+1*rzeta*(R_temp[Cart::x][Cart::yz][0]-gfak*R_temp[Cart::x][Cart::yz][1]);
R_temp[Cart::xxx][Cart::xx][0]+=pma0*R_temp[Cart::xx][Cart::xx][0]+wmp0*R_temp[Cart::xx][Cart::xx][1]+1*rzeta*(R_temp[Cart::x][Cart::xx][0]-gfak*R_temp[Cart::x][Cart::xx][1])+0.5/_decay*2*R_temp[Cart::xx][Cart::x][1];
R_temp[Cart::xxx][Cart::xz][0]+=pma0*R_temp[Cart::xx][Cart::xz][0]+wmp0*R_temp[Cart::xx][Cart::xz][1]+1*rzeta*(R_temp[Cart::x][Cart::xz][0]-gfak*R_temp[Cart::x][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::xx][Cart::z][1];
R_temp[Cart::xxx][Cart::zz][0]+=pma0*R_temp[Cart::xx][Cart::zz][0]+wmp0*R_temp[Cart::xx][Cart::zz][1]+1*rzeta*(R_temp[Cart::x][Cart::zz][0]-gfak*R_temp[Cart::x][Cart::zz][1]);
R_temp[Cart::xxz][Cart::yy][0]+=pma2*R_temp[Cart::xx][Cart::yy][0]+wmp2*R_temp[Cart::xx][Cart::yy][1];
R_temp[Cart::xxz][Cart::xy][0]+=pma2*R_temp[Cart::xx][Cart::xy][0]+wmp2*R_temp[Cart::xx][Cart::xy][1];
R_temp[Cart::xxz][Cart::yz][0]+=pma2*R_temp[Cart::xx][Cart::yz][0]+wmp2*R_temp[Cart::xx][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::y][1];
R_temp[Cart::xxz][Cart::xx][0]+=pma2*R_temp[Cart::xx][Cart::xx][0]+wmp2*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xxz][Cart::xz][0]+=pma2*R_temp[Cart::xx][Cart::xz][0]+wmp2*R_temp[Cart::xx][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][1];
R_temp[Cart::xxz][Cart::zz][0]+=pma2*R_temp[Cart::xx][Cart::zz][0]+wmp2*R_temp[Cart::xx][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::xx][Cart::z][1];
R_temp[Cart::xzz][Cart::yy][0]+=pma0*R_temp[Cart::zz][Cart::yy][0]+wmp0*R_temp[Cart::zz][Cart::yy][1];
R_temp[Cart::xzz][Cart::xy][0]+=pma0*R_temp[Cart::zz][Cart::xy][0]+wmp0*R_temp[Cart::zz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::y][1];
R_temp[Cart::xzz][Cart::yz][0]+=pma0*R_temp[Cart::zz][Cart::yz][0]+wmp0*R_temp[Cart::zz][Cart::yz][1];
R_temp[Cart::xzz][Cart::xx][0]+=pma0*R_temp[Cart::zz][Cart::xx][0]+wmp0*R_temp[Cart::zz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::zz][Cart::x][1];
R_temp[Cart::xzz][Cart::xz][0]+=pma0*R_temp[Cart::zz][Cart::xz][0]+wmp0*R_temp[Cart::zz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][1];
R_temp[Cart::xzz][Cart::zz][0]+=pma0*R_temp[Cart::zz][Cart::zz][0]+wmp0*R_temp[Cart::zz][Cart::zz][1];
R_temp[Cart::zzz][Cart::yy][0]+=pma2*R_temp[Cart::zz][Cart::yy][0]+wmp2*R_temp[Cart::zz][Cart::yy][1]+1*rzeta*(R_temp[Cart::z][Cart::yy][0]-gfak*R_temp[Cart::z][Cart::yy][1]);
R_temp[Cart::zzz][Cart::xy][0]+=pma2*R_temp[Cart::zz][Cart::xy][0]+wmp2*R_temp[Cart::zz][Cart::xy][1]+1*rzeta*(R_temp[Cart::z][Cart::xy][0]-gfak*R_temp[Cart::z][Cart::xy][1]);
R_temp[Cart::zzz][Cart::yz][0]+=pma2*R_temp[Cart::zz][Cart::yz][0]+wmp2*R_temp[Cart::zz][Cart::yz][1]+1*rzeta*(R_temp[Cart::z][Cart::yz][0]-gfak*R_temp[Cart::z][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::zz][Cart::y][1];
R_temp[Cart::zzz][Cart::xx][0]+=pma2*R_temp[Cart::zz][Cart::xx][0]+wmp2*R_temp[Cart::zz][Cart::xx][1]+1*rzeta*(R_temp[Cart::z][Cart::xx][0]-gfak*R_temp[Cart::z][Cart::xx][1]);
R_temp[Cart::zzz][Cart::xz][0]+=pma2*R_temp[Cart::zz][Cart::xz][0]+wmp2*R_temp[Cart::zz][Cart::xz][1]+1*rzeta*(R_temp[Cart::z][Cart::xz][0]-gfak*R_temp[Cart::z][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::zz][Cart::x][1];
R_temp[Cart::zzz][Cart::zz][0]+=pma2*R_temp[Cart::zz][Cart::zz][0]+wmp2*R_temp[Cart::zz][Cart::zz][1]+1*rzeta*(R_temp[Cart::z][Cart::zz][0]-gfak*R_temp[Cart::z][Cart::zz][1])+0.5/_decay*2*R_temp[Cart::zz][Cart::z][1];
}
//------------------------------------------------------

//Integral f - s - d - m1
if (_mmax >1 ){
if (_lmax_alpha>2 && _lmax_gamma>1){
R_temp[Cart::yyy][Cart::yy][1]+=pma1*R_temp[Cart::yy][Cart::yy][1]+wmp1*R_temp[Cart::yy][Cart::yy][2]+1*rzeta*(R_temp[Cart::y][Cart::yy][1]-gfak*R_temp[Cart::y][Cart::yy][2])+0.5/_decay*2*R_temp[Cart::yy][Cart::y][2];
R_temp[Cart::yyy][Cart::xy][1]+=pma1*R_temp[Cart::yy][Cart::xy][1]+wmp1*R_temp[Cart::yy][Cart::xy][2]+1*rzeta*(R_temp[Cart::y][Cart::xy][1]-gfak*R_temp[Cart::y][Cart::xy][2])+0.5/_decay*1*R_temp[Cart::yy][Cart::x][2];
R_temp[Cart::yyy][Cart::yz][1]+=pma1*R_temp[Cart::yy][Cart::yz][1]+wmp1*R_temp[Cart::yy][Cart::yz][2]+1*rzeta*(R_temp[Cart::y][Cart::yz][1]-gfak*R_temp[Cart::y][Cart::yz][2])+0.5/_decay*1*R_temp[Cart::yy][Cart::z][2];
R_temp[Cart::yyy][Cart::xx][1]+=pma1*R_temp[Cart::yy][Cart::xx][1]+wmp1*R_temp[Cart::yy][Cart::xx][2]+1*rzeta*(R_temp[Cart::y][Cart::xx][1]-gfak*R_temp[Cart::y][Cart::xx][2]);
R_temp[Cart::yyy][Cart::xz][1]+=pma1*R_temp[Cart::yy][Cart::xz][1]+wmp1*R_temp[Cart::yy][Cart::xz][2]+1*rzeta*(R_temp[Cart::y][Cart::xz][1]-gfak*R_temp[Cart::y][Cart::xz][2]);
R_temp[Cart::yyy][Cart::zz][1]+=pma1*R_temp[Cart::yy][Cart::zz][1]+wmp1*R_temp[Cart::yy][Cart::zz][2]+1*rzeta*(R_temp[Cart::y][Cart::zz][1]-gfak*R_temp[Cart::y][Cart::zz][2]);
R_temp[Cart::xyy][Cart::yy][1]+=pma0*R_temp[Cart::yy][Cart::yy][1]+wmp0*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::xyy][Cart::xy][1]+=pma0*R_temp[Cart::yy][Cart::xy][1]+wmp0*R_temp[Cart::yy][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][2];
R_temp[Cart::xyy][Cart::yz][1]+=pma0*R_temp[Cart::yy][Cart::yz][1]+wmp0*R_temp[Cart::yy][Cart::yz][2];
R_temp[Cart::xyy][Cart::xx][1]+=pma0*R_temp[Cart::yy][Cart::xx][1]+wmp0*R_temp[Cart::yy][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yy][Cart::x][2];
R_temp[Cart::xyy][Cart::xz][1]+=pma0*R_temp[Cart::yy][Cart::xz][1]+wmp0*R_temp[Cart::yy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::z][2];
R_temp[Cart::xyy][Cart::zz][1]+=pma0*R_temp[Cart::yy][Cart::zz][1]+wmp0*R_temp[Cart::yy][Cart::zz][2];
R_temp[Cart::yyz][Cart::yy][1]+=pma2*R_temp[Cart::yy][Cart::yy][1]+wmp2*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::yyz][Cart::xy][1]+=pma2*R_temp[Cart::yy][Cart::xy][1]+wmp2*R_temp[Cart::yy][Cart::xy][2];
R_temp[Cart::yyz][Cart::yz][1]+=pma2*R_temp[Cart::yy][Cart::yz][1]+wmp2*R_temp[Cart::yy][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][2];
R_temp[Cart::yyz][Cart::xx][1]+=pma2*R_temp[Cart::yy][Cart::xx][1]+wmp2*R_temp[Cart::yy][Cart::xx][2];
R_temp[Cart::yyz][Cart::xz][1]+=pma2*R_temp[Cart::yy][Cart::xz][1]+wmp2*R_temp[Cart::yy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::x][2];
R_temp[Cart::yyz][Cart::zz][1]+=pma2*R_temp[Cart::yy][Cart::zz][1]+wmp2*R_temp[Cart::yy][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::yy][Cart::z][2];
R_temp[Cart::xxy][Cart::yy][1]+=pma1*R_temp[Cart::xx][Cart::yy][1]+wmp1*R_temp[Cart::xx][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xx][Cart::y][2];
R_temp[Cart::xxy][Cart::xy][1]+=pma1*R_temp[Cart::xx][Cart::xy][1]+wmp1*R_temp[Cart::xx][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][2];
R_temp[Cart::xxy][Cart::yz][1]+=pma1*R_temp[Cart::xx][Cart::yz][1]+wmp1*R_temp[Cart::xx][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::z][2];
R_temp[Cart::xxy][Cart::xx][1]+=pma1*R_temp[Cart::xx][Cart::xx][1]+wmp1*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xxy][Cart::xz][1]+=pma1*R_temp[Cart::xx][Cart::xz][1]+wmp1*R_temp[Cart::xx][Cart::xz][2];
R_temp[Cart::xxy][Cart::zz][1]+=pma1*R_temp[Cart::xx][Cart::zz][1]+wmp1*R_temp[Cart::xx][Cart::zz][2];
R_temp[Cart::xyz][Cart::yy][1]+=pma0*R_temp[Cart::yz][Cart::yy][1]+wmp0*R_temp[Cart::yz][Cart::yy][2];
R_temp[Cart::xyz][Cart::xy][1]+=pma0*R_temp[Cart::yz][Cart::xy][1]+wmp0*R_temp[Cart::yz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yz][Cart::y][2];
R_temp[Cart::xyz][Cart::yz][1]+=pma0*R_temp[Cart::yz][Cart::yz][1]+wmp0*R_temp[Cart::yz][Cart::yz][2];
R_temp[Cart::xyz][Cart::xx][1]+=pma0*R_temp[Cart::yz][Cart::xx][1]+wmp0*R_temp[Cart::yz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yz][Cart::x][2];
R_temp[Cart::xyz][Cart::xz][1]+=pma0*R_temp[Cart::yz][Cart::xz][1]+wmp0*R_temp[Cart::yz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yz][Cart::z][2];
R_temp[Cart::xyz][Cart::zz][1]+=pma0*R_temp[Cart::yz][Cart::zz][1]+wmp0*R_temp[Cart::yz][Cart::zz][2];
R_temp[Cart::yzz][Cart::yy][1]+=pma1*R_temp[Cart::zz][Cart::yy][1]+wmp1*R_temp[Cart::zz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::zz][Cart::y][2];
R_temp[Cart::yzz][Cart::xy][1]+=pma1*R_temp[Cart::zz][Cart::xy][1]+wmp1*R_temp[Cart::zz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::x][2];
R_temp[Cart::yzz][Cart::yz][1]+=pma1*R_temp[Cart::zz][Cart::yz][1]+wmp1*R_temp[Cart::zz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][2];
R_temp[Cart::yzz][Cart::xx][1]+=pma1*R_temp[Cart::zz][Cart::xx][1]+wmp1*R_temp[Cart::zz][Cart::xx][2];
R_temp[Cart::yzz][Cart::xz][1]+=pma1*R_temp[Cart::zz][Cart::xz][1]+wmp1*R_temp[Cart::zz][Cart::xz][2];
R_temp[Cart::yzz][Cart::zz][1]+=pma1*R_temp[Cart::zz][Cart::zz][1]+wmp1*R_temp[Cart::zz][Cart::zz][2];
R_temp[Cart::xxx][Cart::yy][1]+=pma0*R_temp[Cart::xx][Cart::yy][1]+wmp0*R_temp[Cart::xx][Cart::yy][2]+1*rzeta*(R_temp[Cart::x][Cart::yy][1]-gfak*R_temp[Cart::x][Cart::yy][2]);
R_temp[Cart::xxx][Cart::xy][1]+=pma0*R_temp[Cart::xx][Cart::xy][1]+wmp0*R_temp[Cart::xx][Cart::xy][2]+1*rzeta*(R_temp[Cart::x][Cart::xy][1]-gfak*R_temp[Cart::x][Cart::xy][2])+0.5/_decay*1*R_temp[Cart::xx][Cart::y][2];
R_temp[Cart::xxx][Cart::yz][1]+=pma0*R_temp[Cart::xx][Cart::yz][1]+wmp0*R_temp[Cart::xx][Cart::yz][2]+1*rzeta*(R_temp[Cart::x][Cart::yz][1]-gfak*R_temp[Cart::x][Cart::yz][2]);
R_temp[Cart::xxx][Cart::xx][1]+=pma0*R_temp[Cart::xx][Cart::xx][1]+wmp0*R_temp[Cart::xx][Cart::xx][2]+1*rzeta*(R_temp[Cart::x][Cart::xx][1]-gfak*R_temp[Cart::x][Cart::xx][2])+0.5/_decay*2*R_temp[Cart::xx][Cart::x][2];
R_temp[Cart::xxx][Cart::xz][1]+=pma0*R_temp[Cart::xx][Cart::xz][1]+wmp0*R_temp[Cart::xx][Cart::xz][2]+1*rzeta*(R_temp[Cart::x][Cart::xz][1]-gfak*R_temp[Cart::x][Cart::xz][2])+0.5/_decay*1*R_temp[Cart::xx][Cart::z][2];
R_temp[Cart::xxx][Cart::zz][1]+=pma0*R_temp[Cart::xx][Cart::zz][1]+wmp0*R_temp[Cart::xx][Cart::zz][2]+1*rzeta*(R_temp[Cart::x][Cart::zz][1]-gfak*R_temp[Cart::x][Cart::zz][2]);
R_temp[Cart::xxz][Cart::yy][1]+=pma2*R_temp[Cart::xx][Cart::yy][1]+wmp2*R_temp[Cart::xx][Cart::yy][2];
R_temp[Cart::xxz][Cart::xy][1]+=pma2*R_temp[Cart::xx][Cart::xy][1]+wmp2*R_temp[Cart::xx][Cart::xy][2];
R_temp[Cart::xxz][Cart::yz][1]+=pma2*R_temp[Cart::xx][Cart::yz][1]+wmp2*R_temp[Cart::xx][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::y][2];
R_temp[Cart::xxz][Cart::xx][1]+=pma2*R_temp[Cart::xx][Cart::xx][1]+wmp2*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xxz][Cart::xz][1]+=pma2*R_temp[Cart::xx][Cart::xz][1]+wmp2*R_temp[Cart::xx][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][2];
R_temp[Cart::xxz][Cart::zz][1]+=pma2*R_temp[Cart::xx][Cart::zz][1]+wmp2*R_temp[Cart::xx][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::xx][Cart::z][2];
R_temp[Cart::xzz][Cart::yy][1]+=pma0*R_temp[Cart::zz][Cart::yy][1]+wmp0*R_temp[Cart::zz][Cart::yy][2];
R_temp[Cart::xzz][Cart::xy][1]+=pma0*R_temp[Cart::zz][Cart::xy][1]+wmp0*R_temp[Cart::zz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::y][2];
R_temp[Cart::xzz][Cart::yz][1]+=pma0*R_temp[Cart::zz][Cart::yz][1]+wmp0*R_temp[Cart::zz][Cart::yz][2];
R_temp[Cart::xzz][Cart::xx][1]+=pma0*R_temp[Cart::zz][Cart::xx][1]+wmp0*R_temp[Cart::zz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::zz][Cart::x][2];
R_temp[Cart::xzz][Cart::xz][1]+=pma0*R_temp[Cart::zz][Cart::xz][1]+wmp0*R_temp[Cart::zz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][2];
R_temp[Cart::xzz][Cart::zz][1]+=pma0*R_temp[Cart::zz][Cart::zz][1]+wmp0*R_temp[Cart::zz][Cart::zz][2];
R_temp[Cart::zzz][Cart::yy][1]+=pma2*R_temp[Cart::zz][Cart::yy][1]+wmp2*R_temp[Cart::zz][Cart::yy][2]+1*rzeta*(R_temp[Cart::z][Cart::yy][1]-gfak*R_temp[Cart::z][Cart::yy][2]);
R_temp[Cart::zzz][Cart::xy][1]+=pma2*R_temp[Cart::zz][Cart::xy][1]+wmp2*R_temp[Cart::zz][Cart::xy][2]+1*rzeta*(R_temp[Cart::z][Cart::xy][1]-gfak*R_temp[Cart::z][Cart::xy][2]);
R_temp[Cart::zzz][Cart::yz][1]+=pma2*R_temp[Cart::zz][Cart::yz][1]+wmp2*R_temp[Cart::zz][Cart::yz][2]+1*rzeta*(R_temp[Cart::z][Cart::yz][1]-gfak*R_temp[Cart::z][Cart::yz][2])+0.5/_decay*1*R_temp[Cart::zz][Cart::y][2];
R_temp[Cart::zzz][Cart::xx][1]+=pma2*R_temp[Cart::zz][Cart::xx][1]+wmp2*R_temp[Cart::zz][Cart::xx][2]+1*rzeta*(R_temp[Cart::z][Cart::xx][1]-gfak*R_temp[Cart::z][Cart::xx][2]);
R_temp[Cart::zzz][Cart::xz][1]+=pma2*R_temp[Cart::zz][Cart::xz][1]+wmp2*R_temp[Cart::zz][Cart::xz][2]+1*rzeta*(R_temp[Cart::z][Cart::xz][1]-gfak*R_temp[Cart::z][Cart::xz][2])+0.5/_decay*1*R_temp[Cart::zz][Cart::x][2];
R_temp[Cart::zzz][Cart::zz][1]+=pma2*R_temp[Cart::zz][Cart::zz][1]+wmp2*R_temp[Cart::zz][Cart::zz][2]+1*rzeta*(R_temp[Cart::z][Cart::zz][1]-gfak*R_temp[Cart::z][Cart::zz][2])+0.5/_decay*2*R_temp[Cart::zz][Cart::z][2];
}
}
//------------------------------------------------------

//Integral f - s - d - m2
if (_mmax >2 ){
if (_lmax_alpha>2 && _lmax_gamma>1){
R_temp[Cart::yyy][Cart::yy][2]+=pma1*R_temp[Cart::yy][Cart::yy][2]+wmp1*R_temp[Cart::yy][Cart::yy][3]+1*rzeta*(R_temp[Cart::y][Cart::yy][2]-gfak*R_temp[Cart::y][Cart::yy][3])+0.5/_decay*2*R_temp[Cart::yy][Cart::y][3];
R_temp[Cart::yyy][Cart::xy][2]+=pma1*R_temp[Cart::yy][Cart::xy][2]+wmp1*R_temp[Cart::yy][Cart::xy][3]+1*rzeta*(R_temp[Cart::y][Cart::xy][2]-gfak*R_temp[Cart::y][Cart::xy][3])+0.5/_decay*1*R_temp[Cart::yy][Cart::x][3];
R_temp[Cart::yyy][Cart::yz][2]+=pma1*R_temp[Cart::yy][Cart::yz][2]+wmp1*R_temp[Cart::yy][Cart::yz][3]+1*rzeta*(R_temp[Cart::y][Cart::yz][2]-gfak*R_temp[Cart::y][Cart::yz][3])+0.5/_decay*1*R_temp[Cart::yy][Cart::z][3];
R_temp[Cart::yyy][Cart::xx][2]+=pma1*R_temp[Cart::yy][Cart::xx][2]+wmp1*R_temp[Cart::yy][Cart::xx][3]+1*rzeta*(R_temp[Cart::y][Cart::xx][2]-gfak*R_temp[Cart::y][Cart::xx][3]);
R_temp[Cart::yyy][Cart::xz][2]+=pma1*R_temp[Cart::yy][Cart::xz][2]+wmp1*R_temp[Cart::yy][Cart::xz][3]+1*rzeta*(R_temp[Cart::y][Cart::xz][2]-gfak*R_temp[Cart::y][Cart::xz][3]);
R_temp[Cart::yyy][Cart::zz][2]+=pma1*R_temp[Cart::yy][Cart::zz][2]+wmp1*R_temp[Cart::yy][Cart::zz][3]+1*rzeta*(R_temp[Cart::y][Cart::zz][2]-gfak*R_temp[Cart::y][Cart::zz][3]);
R_temp[Cart::xyy][Cart::yy][2]+=pma0*R_temp[Cart::yy][Cart::yy][2]+wmp0*R_temp[Cart::yy][Cart::yy][3];
R_temp[Cart::xyy][Cart::xy][2]+=pma0*R_temp[Cart::yy][Cart::xy][2]+wmp0*R_temp[Cart::yy][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][3];
R_temp[Cart::xyy][Cart::yz][2]+=pma0*R_temp[Cart::yy][Cart::yz][2]+wmp0*R_temp[Cart::yy][Cart::yz][3];
R_temp[Cart::xyy][Cart::xx][2]+=pma0*R_temp[Cart::yy][Cart::xx][2]+wmp0*R_temp[Cart::yy][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::yy][Cart::x][3];
R_temp[Cart::xyy][Cart::xz][2]+=pma0*R_temp[Cart::yy][Cart::xz][2]+wmp0*R_temp[Cart::yy][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yy][Cart::z][3];
R_temp[Cart::xyy][Cart::zz][2]+=pma0*R_temp[Cart::yy][Cart::zz][2]+wmp0*R_temp[Cart::yy][Cart::zz][3];
R_temp[Cart::yyz][Cart::yy][2]+=pma2*R_temp[Cart::yy][Cart::yy][2]+wmp2*R_temp[Cart::yy][Cart::yy][3];
R_temp[Cart::yyz][Cart::xy][2]+=pma2*R_temp[Cart::yy][Cart::xy][2]+wmp2*R_temp[Cart::yy][Cart::xy][3];
R_temp[Cart::yyz][Cart::yz][2]+=pma2*R_temp[Cart::yy][Cart::yz][2]+wmp2*R_temp[Cart::yy][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::yy][Cart::y][3];
R_temp[Cart::yyz][Cart::xx][2]+=pma2*R_temp[Cart::yy][Cart::xx][2]+wmp2*R_temp[Cart::yy][Cart::xx][3];
R_temp[Cart::yyz][Cart::xz][2]+=pma2*R_temp[Cart::yy][Cart::xz][2]+wmp2*R_temp[Cart::yy][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yy][Cart::x][3];
R_temp[Cart::yyz][Cart::zz][2]+=pma2*R_temp[Cart::yy][Cart::zz][2]+wmp2*R_temp[Cart::yy][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::yy][Cart::z][3];
R_temp[Cart::xxy][Cart::yy][2]+=pma1*R_temp[Cart::xx][Cart::yy][2]+wmp1*R_temp[Cart::xx][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::xx][Cart::y][3];
R_temp[Cart::xxy][Cart::xy][2]+=pma1*R_temp[Cart::xx][Cart::xy][2]+wmp1*R_temp[Cart::xx][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][3];
R_temp[Cart::xxy][Cart::yz][2]+=pma1*R_temp[Cart::xx][Cart::yz][2]+wmp1*R_temp[Cart::xx][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::xx][Cart::z][3];
R_temp[Cart::xxy][Cart::xx][2]+=pma1*R_temp[Cart::xx][Cart::xx][2]+wmp1*R_temp[Cart::xx][Cart::xx][3];
R_temp[Cart::xxy][Cart::xz][2]+=pma1*R_temp[Cart::xx][Cart::xz][2]+wmp1*R_temp[Cart::xx][Cart::xz][3];
R_temp[Cart::xxy][Cart::zz][2]+=pma1*R_temp[Cart::xx][Cart::zz][2]+wmp1*R_temp[Cart::xx][Cart::zz][3];
R_temp[Cart::xyz][Cart::yy][2]+=pma0*R_temp[Cart::yz][Cart::yy][2]+wmp0*R_temp[Cart::yz][Cart::yy][3];
R_temp[Cart::xyz][Cart::xy][2]+=pma0*R_temp[Cart::yz][Cart::xy][2]+wmp0*R_temp[Cart::yz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::yz][Cart::y][3];
R_temp[Cart::xyz][Cart::yz][2]+=pma0*R_temp[Cart::yz][Cart::yz][2]+wmp0*R_temp[Cart::yz][Cart::yz][3];
R_temp[Cart::xyz][Cart::xx][2]+=pma0*R_temp[Cart::yz][Cart::xx][2]+wmp0*R_temp[Cart::yz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::yz][Cart::x][3];
R_temp[Cart::xyz][Cart::xz][2]+=pma0*R_temp[Cart::yz][Cart::xz][2]+wmp0*R_temp[Cart::yz][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::yz][Cart::z][3];
R_temp[Cart::xyz][Cart::zz][2]+=pma0*R_temp[Cart::yz][Cart::zz][2]+wmp0*R_temp[Cart::yz][Cart::zz][3];
R_temp[Cart::yzz][Cart::yy][2]+=pma1*R_temp[Cart::zz][Cart::yy][2]+wmp1*R_temp[Cart::zz][Cart::yy][3]+0.5/_decay*2*R_temp[Cart::zz][Cart::y][3];
R_temp[Cart::yzz][Cart::xy][2]+=pma1*R_temp[Cart::zz][Cart::xy][2]+wmp1*R_temp[Cart::zz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::zz][Cart::x][3];
R_temp[Cart::yzz][Cart::yz][2]+=pma1*R_temp[Cart::zz][Cart::yz][2]+wmp1*R_temp[Cart::zz][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][3];
R_temp[Cart::yzz][Cart::xx][2]+=pma1*R_temp[Cart::zz][Cart::xx][2]+wmp1*R_temp[Cart::zz][Cart::xx][3];
R_temp[Cart::yzz][Cart::xz][2]+=pma1*R_temp[Cart::zz][Cart::xz][2]+wmp1*R_temp[Cart::zz][Cart::xz][3];
R_temp[Cart::yzz][Cart::zz][2]+=pma1*R_temp[Cart::zz][Cart::zz][2]+wmp1*R_temp[Cart::zz][Cart::zz][3];
R_temp[Cart::xxx][Cart::yy][2]+=pma0*R_temp[Cart::xx][Cart::yy][2]+wmp0*R_temp[Cart::xx][Cart::yy][3]+1*rzeta*(R_temp[Cart::x][Cart::yy][2]-gfak*R_temp[Cart::x][Cart::yy][3]);
R_temp[Cart::xxx][Cart::xy][2]+=pma0*R_temp[Cart::xx][Cart::xy][2]+wmp0*R_temp[Cart::xx][Cart::xy][3]+1*rzeta*(R_temp[Cart::x][Cart::xy][2]-gfak*R_temp[Cart::x][Cart::xy][3])+0.5/_decay*1*R_temp[Cart::xx][Cart::y][3];
R_temp[Cart::xxx][Cart::yz][2]+=pma0*R_temp[Cart::xx][Cart::yz][2]+wmp0*R_temp[Cart::xx][Cart::yz][3]+1*rzeta*(R_temp[Cart::x][Cart::yz][2]-gfak*R_temp[Cart::x][Cart::yz][3]);
R_temp[Cart::xxx][Cart::xx][2]+=pma0*R_temp[Cart::xx][Cart::xx][2]+wmp0*R_temp[Cart::xx][Cart::xx][3]+1*rzeta*(R_temp[Cart::x][Cart::xx][2]-gfak*R_temp[Cart::x][Cart::xx][3])+0.5/_decay*2*R_temp[Cart::xx][Cart::x][3];
R_temp[Cart::xxx][Cart::xz][2]+=pma0*R_temp[Cart::xx][Cart::xz][2]+wmp0*R_temp[Cart::xx][Cart::xz][3]+1*rzeta*(R_temp[Cart::x][Cart::xz][2]-gfak*R_temp[Cart::x][Cart::xz][3])+0.5/_decay*1*R_temp[Cart::xx][Cart::z][3];
R_temp[Cart::xxx][Cart::zz][2]+=pma0*R_temp[Cart::xx][Cart::zz][2]+wmp0*R_temp[Cart::xx][Cart::zz][3]+1*rzeta*(R_temp[Cart::x][Cart::zz][2]-gfak*R_temp[Cart::x][Cart::zz][3]);
R_temp[Cart::xxz][Cart::yy][2]+=pma2*R_temp[Cart::xx][Cart::yy][2]+wmp2*R_temp[Cart::xx][Cart::yy][3];
R_temp[Cart::xxz][Cart::xy][2]+=pma2*R_temp[Cart::xx][Cart::xy][2]+wmp2*R_temp[Cart::xx][Cart::xy][3];
R_temp[Cart::xxz][Cart::yz][2]+=pma2*R_temp[Cart::xx][Cart::yz][2]+wmp2*R_temp[Cart::xx][Cart::yz][3]+0.5/_decay*1*R_temp[Cart::xx][Cart::y][3];
R_temp[Cart::xxz][Cart::xx][2]+=pma2*R_temp[Cart::xx][Cart::xx][2]+wmp2*R_temp[Cart::xx][Cart::xx][3];
R_temp[Cart::xxz][Cart::xz][2]+=pma2*R_temp[Cart::xx][Cart::xz][2]+wmp2*R_temp[Cart::xx][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::xx][Cart::x][3];
R_temp[Cart::xxz][Cart::zz][2]+=pma2*R_temp[Cart::xx][Cart::zz][2]+wmp2*R_temp[Cart::xx][Cart::zz][3]+0.5/_decay*2*R_temp[Cart::xx][Cart::z][3];
R_temp[Cart::xzz][Cart::yy][2]+=pma0*R_temp[Cart::zz][Cart::yy][2]+wmp0*R_temp[Cart::zz][Cart::yy][3];
R_temp[Cart::xzz][Cart::xy][2]+=pma0*R_temp[Cart::zz][Cart::xy][2]+wmp0*R_temp[Cart::zz][Cart::xy][3]+0.5/_decay*1*R_temp[Cart::zz][Cart::y][3];
R_temp[Cart::xzz][Cart::yz][2]+=pma0*R_temp[Cart::zz][Cart::yz][2]+wmp0*R_temp[Cart::zz][Cart::yz][3];
R_temp[Cart::xzz][Cart::xx][2]+=pma0*R_temp[Cart::zz][Cart::xx][2]+wmp0*R_temp[Cart::zz][Cart::xx][3]+0.5/_decay*2*R_temp[Cart::zz][Cart::x][3];
R_temp[Cart::xzz][Cart::xz][2]+=pma0*R_temp[Cart::zz][Cart::xz][2]+wmp0*R_temp[Cart::zz][Cart::xz][3]+0.5/_decay*1*R_temp[Cart::zz][Cart::z][3];
R_temp[Cart::xzz][Cart::zz][2]+=pma0*R_temp[Cart::zz][Cart::zz][2]+wmp0*R_temp[Cart::zz][Cart::zz][3];
R_temp[Cart::zzz][Cart::yy][2]+=pma2*R_temp[Cart::zz][Cart::yy][2]+wmp2*R_temp[Cart::zz][Cart::yy][3]+1*rzeta*(R_temp[Cart::z][Cart::yy][2]-gfak*R_temp[Cart::z][Cart::yy][3]);
R_temp[Cart::zzz][Cart::xy][2]+=pma2*R_temp[Cart::zz][Cart::xy][2]+wmp2*R_temp[Cart::zz][Cart::xy][3]+1*rzeta*(R_temp[Cart::z][Cart::xy][2]-gfak*R_temp[Cart::z][Cart::xy][3]);
R_temp[Cart::zzz][Cart::yz][2]+=pma2*R_temp[Cart::zz][Cart::yz][2]+wmp2*R_temp[Cart::zz][Cart::yz][3]+1*rzeta*(R_temp[Cart::z][Cart::yz][2]-gfak*R_temp[Cart::z][Cart::yz][3])+0.5/_decay*1*R_temp[Cart::zz][Cart::y][3];
R_temp[Cart::zzz][Cart::xx][2]+=pma2*R_temp[Cart::zz][Cart::xx][2]+wmp2*R_temp[Cart::zz][Cart::xx][3]+1*rzeta*(R_temp[Cart::z][Cart::xx][2]-gfak*R_temp[Cart::z][Cart::xx][3]);
R_temp[Cart::zzz][Cart::xz][2]+=pma2*R_temp[Cart::zz][Cart::xz][2]+wmp2*R_temp[Cart::zz][Cart::xz][3]+1*rzeta*(R_temp[Cart::z][Cart::xz][2]-gfak*R_temp[Cart::z][Cart::xz][3])+0.5/_decay*1*R_temp[Cart::zz][Cart::x][3];
R_temp[Cart::zzz][Cart::zz][2]+=pma2*R_temp[Cart::zz][Cart::zz][2]+wmp2*R_temp[Cart::zz][Cart::zz][3]+1*rzeta*(R_temp[Cart::z][Cart::zz][2]-gfak*R_temp[Cart::z][Cart::zz][3])+0.5/_decay*2*R_temp[Cart::zz][Cart::z][3];
}
}
//------------------------------------------------------

//Integral f - s - f - m0
if (_lmax_alpha>2 && _lmax_gamma>2){
R_temp[Cart::yyy][Cart::yyy][0]+=pma1*R_temp[Cart::yy][Cart::yyy][0]+wmp1*R_temp[Cart::yy][Cart::yyy][1]+1*rzeta*(R_temp[Cart::y][Cart::yyy][0]-gfak*R_temp[Cart::y][Cart::yyy][1])+0.5/_decay*3*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::yyy][Cart::xyy][0]+=pma1*R_temp[Cart::yy][Cart::xyy][0]+wmp1*R_temp[Cart::yy][Cart::xyy][1]+1*rzeta*(R_temp[Cart::y][Cart::xyy][0]-gfak*R_temp[Cart::y][Cart::xyy][1])+0.5/_decay*2*R_temp[Cart::yy][Cart::xy][1];
R_temp[Cart::yyy][Cart::yyz][0]+=pma1*R_temp[Cart::yy][Cart::yyz][0]+wmp1*R_temp[Cart::yy][Cart::yyz][1]+1*rzeta*(R_temp[Cart::y][Cart::yyz][0]-gfak*R_temp[Cart::y][Cart::yyz][1])+0.5/_decay*2*R_temp[Cart::yy][Cart::yz][1];
R_temp[Cart::yyy][Cart::xxy][0]+=pma1*R_temp[Cart::yy][Cart::xxy][0]+wmp1*R_temp[Cart::yy][Cart::xxy][1]+1*rzeta*(R_temp[Cart::y][Cart::xxy][0]-gfak*R_temp[Cart::y][Cart::xxy][1])+0.5/_decay*1*R_temp[Cart::yy][Cart::xx][1];
R_temp[Cart::yyy][Cart::xyz][0]+=pma1*R_temp[Cart::yy][Cart::xyz][0]+wmp1*R_temp[Cart::yy][Cart::xyz][1]+1*rzeta*(R_temp[Cart::y][Cart::xyz][0]-gfak*R_temp[Cart::y][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::yy][Cart::xz][1];
R_temp[Cart::yyy][Cart::yzz][0]+=pma1*R_temp[Cart::yy][Cart::yzz][0]+wmp1*R_temp[Cart::yy][Cart::yzz][1]+1*rzeta*(R_temp[Cart::y][Cart::yzz][0]-gfak*R_temp[Cart::y][Cart::yzz][1])+0.5/_decay*1*R_temp[Cart::yy][Cart::zz][1];
R_temp[Cart::yyy][Cart::xxx][0]+=pma1*R_temp[Cart::yy][Cart::xxx][0]+wmp1*R_temp[Cart::yy][Cart::xxx][1]+1*rzeta*(R_temp[Cart::y][Cart::xxx][0]-gfak*R_temp[Cart::y][Cart::xxx][1]);
R_temp[Cart::yyy][Cart::xxz][0]+=pma1*R_temp[Cart::yy][Cart::xxz][0]+wmp1*R_temp[Cart::yy][Cart::xxz][1]+1*rzeta*(R_temp[Cart::y][Cart::xxz][0]-gfak*R_temp[Cart::y][Cart::xxz][1]);
R_temp[Cart::yyy][Cart::xzz][0]+=pma1*R_temp[Cart::yy][Cart::xzz][0]+wmp1*R_temp[Cart::yy][Cart::xzz][1]+1*rzeta*(R_temp[Cart::y][Cart::xzz][0]-gfak*R_temp[Cart::y][Cart::xzz][1]);
R_temp[Cart::yyy][Cart::zzz][0]+=pma1*R_temp[Cart::yy][Cart::zzz][0]+wmp1*R_temp[Cart::yy][Cart::zzz][1]+1*rzeta*(R_temp[Cart::y][Cart::zzz][0]-gfak*R_temp[Cart::y][Cart::zzz][1]);
R_temp[Cart::xyy][Cart::yyy][0]+=pma0*R_temp[Cart::yy][Cart::yyy][0]+wmp0*R_temp[Cart::yy][Cart::yyy][1];
R_temp[Cart::xyy][Cart::xyy][0]+=pma0*R_temp[Cart::yy][Cart::xyy][0]+wmp0*R_temp[Cart::yy][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::xyy][Cart::yyz][0]+=pma0*R_temp[Cart::yy][Cart::yyz][0]+wmp0*R_temp[Cart::yy][Cart::yyz][1];
R_temp[Cart::xyy][Cart::xxy][0]+=pma0*R_temp[Cart::yy][Cart::xxy][0]+wmp0*R_temp[Cart::yy][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yy][Cart::xy][1];
R_temp[Cart::xyy][Cart::xyz][0]+=pma0*R_temp[Cart::yy][Cart::xyz][0]+wmp0*R_temp[Cart::yy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::yz][1];
R_temp[Cart::xyy][Cart::yzz][0]+=pma0*R_temp[Cart::yy][Cart::yzz][0]+wmp0*R_temp[Cart::yy][Cart::yzz][1];
R_temp[Cart::xyy][Cart::xxx][0]+=pma0*R_temp[Cart::yy][Cart::xxx][0]+wmp0*R_temp[Cart::yy][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yy][Cart::xx][1];
R_temp[Cart::xyy][Cart::xxz][0]+=pma0*R_temp[Cart::yy][Cart::xxz][0]+wmp0*R_temp[Cart::yy][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yy][Cart::xz][1];
R_temp[Cart::xyy][Cart::xzz][0]+=pma0*R_temp[Cart::yy][Cart::xzz][0]+wmp0*R_temp[Cart::yy][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::zz][1];
R_temp[Cart::xyy][Cart::zzz][0]+=pma0*R_temp[Cart::yy][Cart::zzz][0]+wmp0*R_temp[Cart::yy][Cart::zzz][1];
R_temp[Cart::yyz][Cart::yyy][0]+=pma2*R_temp[Cart::yy][Cart::yyy][0]+wmp2*R_temp[Cart::yy][Cart::yyy][1];
R_temp[Cart::yyz][Cart::xyy][0]+=pma2*R_temp[Cart::yy][Cart::xyy][0]+wmp2*R_temp[Cart::yy][Cart::xyy][1];
R_temp[Cart::yyz][Cart::yyz][0]+=pma2*R_temp[Cart::yy][Cart::yyz][0]+wmp2*R_temp[Cart::yy][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::yy][1];
R_temp[Cart::yyz][Cart::xxy][0]+=pma2*R_temp[Cart::yy][Cart::xxy][0]+wmp2*R_temp[Cart::yy][Cart::xxy][1];
R_temp[Cart::yyz][Cart::xyz][0]+=pma2*R_temp[Cart::yy][Cart::xyz][0]+wmp2*R_temp[Cart::yy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::xy][1];
R_temp[Cart::yyz][Cart::yzz][0]+=pma2*R_temp[Cart::yy][Cart::yzz][0]+wmp2*R_temp[Cart::yy][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::yy][Cart::yz][1];
R_temp[Cart::yyz][Cart::xxx][0]+=pma2*R_temp[Cart::yy][Cart::xxx][0]+wmp2*R_temp[Cart::yy][Cart::xxx][1];
R_temp[Cart::yyz][Cart::xxz][0]+=pma2*R_temp[Cart::yy][Cart::xxz][0]+wmp2*R_temp[Cart::yy][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::yy][Cart::xx][1];
R_temp[Cart::yyz][Cart::xzz][0]+=pma2*R_temp[Cart::yy][Cart::xzz][0]+wmp2*R_temp[Cart::yy][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::yy][Cart::xz][1];
R_temp[Cart::yyz][Cart::zzz][0]+=pma2*R_temp[Cart::yy][Cart::zzz][0]+wmp2*R_temp[Cart::yy][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::yy][Cart::zz][1];
R_temp[Cart::xxy][Cart::yyy][0]+=pma1*R_temp[Cart::xx][Cart::yyy][0]+wmp1*R_temp[Cart::xx][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xx][Cart::yy][1];
R_temp[Cart::xxy][Cart::xyy][0]+=pma1*R_temp[Cart::xx][Cart::xyy][0]+wmp1*R_temp[Cart::xx][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xx][Cart::xy][1];
R_temp[Cart::xxy][Cart::yyz][0]+=pma1*R_temp[Cart::xx][Cart::yyz][0]+wmp1*R_temp[Cart::xx][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xx][Cart::yz][1];
R_temp[Cart::xxy][Cart::xxy][0]+=pma1*R_temp[Cart::xx][Cart::xxy][0]+wmp1*R_temp[Cart::xx][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xxy][Cart::xyz][0]+=pma1*R_temp[Cart::xx][Cart::xyz][0]+wmp1*R_temp[Cart::xx][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::xz][1];
R_temp[Cart::xxy][Cart::yzz][0]+=pma1*R_temp[Cart::xx][Cart::yzz][0]+wmp1*R_temp[Cart::xx][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::zz][1];
R_temp[Cart::xxy][Cart::xxx][0]+=pma1*R_temp[Cart::xx][Cart::xxx][0]+wmp1*R_temp[Cart::xx][Cart::xxx][1];
R_temp[Cart::xxy][Cart::xxz][0]+=pma1*R_temp[Cart::xx][Cart::xxz][0]+wmp1*R_temp[Cart::xx][Cart::xxz][1];
R_temp[Cart::xxy][Cart::xzz][0]+=pma1*R_temp[Cart::xx][Cart::xzz][0]+wmp1*R_temp[Cart::xx][Cart::xzz][1];
R_temp[Cart::xxy][Cart::zzz][0]+=pma1*R_temp[Cart::xx][Cart::zzz][0]+wmp1*R_temp[Cart::xx][Cart::zzz][1];
R_temp[Cart::xyz][Cart::yyy][0]+=pma0*R_temp[Cart::yz][Cart::yyy][0]+wmp0*R_temp[Cart::yz][Cart::yyy][1];
R_temp[Cart::xyz][Cart::xyy][0]+=pma0*R_temp[Cart::yz][Cart::xyy][0]+wmp0*R_temp[Cart::yz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yz][Cart::yy][1];
R_temp[Cart::xyz][Cart::yyz][0]+=pma0*R_temp[Cart::yz][Cart::yyz][0]+wmp0*R_temp[Cart::yz][Cart::yyz][1];
R_temp[Cart::xyz][Cart::xxy][0]+=pma0*R_temp[Cart::yz][Cart::xxy][0]+wmp0*R_temp[Cart::yz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yz][Cart::xy][1];
R_temp[Cart::xyz][Cart::xyz][0]+=pma0*R_temp[Cart::yz][Cart::xyz][0]+wmp0*R_temp[Cart::yz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yz][Cart::yz][1];
R_temp[Cart::xyz][Cart::yzz][0]+=pma0*R_temp[Cart::yz][Cart::yzz][0]+wmp0*R_temp[Cart::yz][Cart::yzz][1];
R_temp[Cart::xyz][Cart::xxx][0]+=pma0*R_temp[Cart::yz][Cart::xxx][0]+wmp0*R_temp[Cart::yz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yz][Cart::xx][1];
R_temp[Cart::xyz][Cart::xxz][0]+=pma0*R_temp[Cart::yz][Cart::xxz][0]+wmp0*R_temp[Cart::yz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yz][Cart::xz][1];
R_temp[Cart::xyz][Cart::xzz][0]+=pma0*R_temp[Cart::yz][Cart::xzz][0]+wmp0*R_temp[Cart::yz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yz][Cart::zz][1];
R_temp[Cart::xyz][Cart::zzz][0]+=pma0*R_temp[Cart::yz][Cart::zzz][0]+wmp0*R_temp[Cart::yz][Cart::zzz][1];
R_temp[Cart::yzz][Cart::yyy][0]+=pma1*R_temp[Cart::zz][Cart::yyy][0]+wmp1*R_temp[Cart::zz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::zz][Cart::yy][1];
R_temp[Cart::yzz][Cart::xyy][0]+=pma1*R_temp[Cart::zz][Cart::xyy][0]+wmp1*R_temp[Cart::zz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::zz][Cart::xy][1];
R_temp[Cart::yzz][Cart::yyz][0]+=pma1*R_temp[Cart::zz][Cart::yyz][0]+wmp1*R_temp[Cart::zz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::zz][Cart::yz][1];
R_temp[Cart::yzz][Cart::xxy][0]+=pma1*R_temp[Cart::zz][Cart::xxy][0]+wmp1*R_temp[Cart::zz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::xx][1];
R_temp[Cart::yzz][Cart::xyz][0]+=pma1*R_temp[Cart::zz][Cart::xyz][0]+wmp1*R_temp[Cart::zz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::xz][1];
R_temp[Cart::yzz][Cart::yzz][0]+=pma1*R_temp[Cart::zz][Cart::yzz][0]+wmp1*R_temp[Cart::zz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::zz][1];
R_temp[Cart::yzz][Cart::xxx][0]+=pma1*R_temp[Cart::zz][Cart::xxx][0]+wmp1*R_temp[Cart::zz][Cart::xxx][1];
R_temp[Cart::yzz][Cart::xxz][0]+=pma1*R_temp[Cart::zz][Cart::xxz][0]+wmp1*R_temp[Cart::zz][Cart::xxz][1];
R_temp[Cart::yzz][Cart::xzz][0]+=pma1*R_temp[Cart::zz][Cart::xzz][0]+wmp1*R_temp[Cart::zz][Cart::xzz][1];
R_temp[Cart::yzz][Cart::zzz][0]+=pma1*R_temp[Cart::zz][Cart::zzz][0]+wmp1*R_temp[Cart::zz][Cart::zzz][1];
R_temp[Cart::xxx][Cart::yyy][0]+=pma0*R_temp[Cart::xx][Cart::yyy][0]+wmp0*R_temp[Cart::xx][Cart::yyy][1]+1*rzeta*(R_temp[Cart::x][Cart::yyy][0]-gfak*R_temp[Cart::x][Cart::yyy][1]);
R_temp[Cart::xxx][Cart::xyy][0]+=pma0*R_temp[Cart::xx][Cart::xyy][0]+wmp0*R_temp[Cart::xx][Cart::xyy][1]+1*rzeta*(R_temp[Cart::x][Cart::xyy][0]-gfak*R_temp[Cart::x][Cart::xyy][1])+0.5/_decay*1*R_temp[Cart::xx][Cart::yy][1];
R_temp[Cart::xxx][Cart::yyz][0]+=pma0*R_temp[Cart::xx][Cart::yyz][0]+wmp0*R_temp[Cart::xx][Cart::yyz][1]+1*rzeta*(R_temp[Cart::x][Cart::yyz][0]-gfak*R_temp[Cart::x][Cart::yyz][1]);
R_temp[Cart::xxx][Cart::xxy][0]+=pma0*R_temp[Cart::xx][Cart::xxy][0]+wmp0*R_temp[Cart::xx][Cart::xxy][1]+1*rzeta*(R_temp[Cart::x][Cart::xxy][0]-gfak*R_temp[Cart::x][Cart::xxy][1])+0.5/_decay*2*R_temp[Cart::xx][Cart::xy][1];
R_temp[Cart::xxx][Cart::xyz][0]+=pma0*R_temp[Cart::xx][Cart::xyz][0]+wmp0*R_temp[Cart::xx][Cart::xyz][1]+1*rzeta*(R_temp[Cart::x][Cart::xyz][0]-gfak*R_temp[Cart::x][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::xx][Cart::yz][1];
R_temp[Cart::xxx][Cart::yzz][0]+=pma0*R_temp[Cart::xx][Cart::yzz][0]+wmp0*R_temp[Cart::xx][Cart::yzz][1]+1*rzeta*(R_temp[Cart::x][Cart::yzz][0]-gfak*R_temp[Cart::x][Cart::yzz][1]);
R_temp[Cart::xxx][Cart::xxx][0]+=pma0*R_temp[Cart::xx][Cart::xxx][0]+wmp0*R_temp[Cart::xx][Cart::xxx][1]+1*rzeta*(R_temp[Cart::x][Cart::xxx][0]-gfak*R_temp[Cart::x][Cart::xxx][1])+0.5/_decay*3*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xxx][Cart::xxz][0]+=pma0*R_temp[Cart::xx][Cart::xxz][0]+wmp0*R_temp[Cart::xx][Cart::xxz][1]+1*rzeta*(R_temp[Cart::x][Cart::xxz][0]-gfak*R_temp[Cart::x][Cart::xxz][1])+0.5/_decay*2*R_temp[Cart::xx][Cart::xz][1];
R_temp[Cart::xxx][Cart::xzz][0]+=pma0*R_temp[Cart::xx][Cart::xzz][0]+wmp0*R_temp[Cart::xx][Cart::xzz][1]+1*rzeta*(R_temp[Cart::x][Cart::xzz][0]-gfak*R_temp[Cart::x][Cart::xzz][1])+0.5/_decay*1*R_temp[Cart::xx][Cart::zz][1];
R_temp[Cart::xxx][Cart::zzz][0]+=pma0*R_temp[Cart::xx][Cart::zzz][0]+wmp0*R_temp[Cart::xx][Cart::zzz][1]+1*rzeta*(R_temp[Cart::x][Cart::zzz][0]-gfak*R_temp[Cart::x][Cart::zzz][1]);
R_temp[Cart::xxz][Cart::yyy][0]+=pma2*R_temp[Cart::xx][Cart::yyy][0]+wmp2*R_temp[Cart::xx][Cart::yyy][1];
R_temp[Cart::xxz][Cart::xyy][0]+=pma2*R_temp[Cart::xx][Cart::xyy][0]+wmp2*R_temp[Cart::xx][Cart::xyy][1];
R_temp[Cart::xxz][Cart::yyz][0]+=pma2*R_temp[Cart::xx][Cart::yyz][0]+wmp2*R_temp[Cart::xx][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::yy][1];
R_temp[Cart::xxz][Cart::xxy][0]+=pma2*R_temp[Cart::xx][Cart::xxy][0]+wmp2*R_temp[Cart::xx][Cart::xxy][1];
R_temp[Cart::xxz][Cart::xyz][0]+=pma2*R_temp[Cart::xx][Cart::xyz][0]+wmp2*R_temp[Cart::xx][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::xy][1];
R_temp[Cart::xxz][Cart::yzz][0]+=pma2*R_temp[Cart::xx][Cart::yzz][0]+wmp2*R_temp[Cart::xx][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::xx][Cart::yz][1];
R_temp[Cart::xxz][Cart::xxx][0]+=pma2*R_temp[Cart::xx][Cart::xxx][0]+wmp2*R_temp[Cart::xx][Cart::xxx][1];
R_temp[Cart::xxz][Cart::xxz][0]+=pma2*R_temp[Cart::xx][Cart::xxz][0]+wmp2*R_temp[Cart::xx][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::xx][Cart::xx][1];
R_temp[Cart::xxz][Cart::xzz][0]+=pma2*R_temp[Cart::xx][Cart::xzz][0]+wmp2*R_temp[Cart::xx][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::xx][Cart::xz][1];
R_temp[Cart::xxz][Cart::zzz][0]+=pma2*R_temp[Cart::xx][Cart::zzz][0]+wmp2*R_temp[Cart::xx][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::xx][Cart::zz][1];
R_temp[Cart::xzz][Cart::yyy][0]+=pma0*R_temp[Cart::zz][Cart::yyy][0]+wmp0*R_temp[Cart::zz][Cart::yyy][1];
R_temp[Cart::xzz][Cart::xyy][0]+=pma0*R_temp[Cart::zz][Cart::xyy][0]+wmp0*R_temp[Cart::zz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::yy][1];
R_temp[Cart::xzz][Cart::yyz][0]+=pma0*R_temp[Cart::zz][Cart::yyz][0]+wmp0*R_temp[Cart::zz][Cart::yyz][1];
R_temp[Cart::xzz][Cart::xxy][0]+=pma0*R_temp[Cart::zz][Cart::xxy][0]+wmp0*R_temp[Cart::zz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::zz][Cart::xy][1];
R_temp[Cart::xzz][Cart::xyz][0]+=pma0*R_temp[Cart::zz][Cart::xyz][0]+wmp0*R_temp[Cart::zz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::yz][1];
R_temp[Cart::xzz][Cart::yzz][0]+=pma0*R_temp[Cart::zz][Cart::yzz][0]+wmp0*R_temp[Cart::zz][Cart::yzz][1];
R_temp[Cart::xzz][Cart::xxx][0]+=pma0*R_temp[Cart::zz][Cart::xxx][0]+wmp0*R_temp[Cart::zz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::zz][Cart::xx][1];
R_temp[Cart::xzz][Cart::xxz][0]+=pma0*R_temp[Cart::zz][Cart::xxz][0]+wmp0*R_temp[Cart::zz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::zz][Cart::xz][1];
R_temp[Cart::xzz][Cart::xzz][0]+=pma0*R_temp[Cart::zz][Cart::xzz][0]+wmp0*R_temp[Cart::zz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::zz][Cart::zz][1];
R_temp[Cart::xzz][Cart::zzz][0]+=pma0*R_temp[Cart::zz][Cart::zzz][0]+wmp0*R_temp[Cart::zz][Cart::zzz][1];
R_temp[Cart::zzz][Cart::yyy][0]+=pma2*R_temp[Cart::zz][Cart::yyy][0]+wmp2*R_temp[Cart::zz][Cart::yyy][1]+1*rzeta*(R_temp[Cart::z][Cart::yyy][0]-gfak*R_temp[Cart::z][Cart::yyy][1]);
R_temp[Cart::zzz][Cart::xyy][0]+=pma2*R_temp[Cart::zz][Cart::xyy][0]+wmp2*R_temp[Cart::zz][Cart::xyy][1]+1*rzeta*(R_temp[Cart::z][Cart::xyy][0]-gfak*R_temp[Cart::z][Cart::xyy][1]);
R_temp[Cart::zzz][Cart::yyz][0]+=pma2*R_temp[Cart::zz][Cart::yyz][0]+wmp2*R_temp[Cart::zz][Cart::yyz][1]+1*rzeta*(R_temp[Cart::z][Cart::yyz][0]-gfak*R_temp[Cart::z][Cart::yyz][1])+0.5/_decay*1*R_temp[Cart::zz][Cart::yy][1];
R_temp[Cart::zzz][Cart::xxy][0]+=pma2*R_temp[Cart::zz][Cart::xxy][0]+wmp2*R_temp[Cart::zz][Cart::xxy][1]+1*rzeta*(R_temp[Cart::z][Cart::xxy][0]-gfak*R_temp[Cart::z][Cart::xxy][1]);
R_temp[Cart::zzz][Cart::xyz][0]+=pma2*R_temp[Cart::zz][Cart::xyz][0]+wmp2*R_temp[Cart::zz][Cart::xyz][1]+1*rzeta*(R_temp[Cart::z][Cart::xyz][0]-gfak*R_temp[Cart::z][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::zz][Cart::xy][1];
R_temp[Cart::zzz][Cart::yzz][0]+=pma2*R_temp[Cart::zz][Cart::yzz][0]+wmp2*R_temp[Cart::zz][Cart::yzz][1]+1*rzeta*(R_temp[Cart::z][Cart::yzz][0]-gfak*R_temp[Cart::z][Cart::yzz][1])+0.5/_decay*2*R_temp[Cart::zz][Cart::yz][1];
R_temp[Cart::zzz][Cart::xxx][0]+=pma2*R_temp[Cart::zz][Cart::xxx][0]+wmp2*R_temp[Cart::zz][Cart::xxx][1]+1*rzeta*(R_temp[Cart::z][Cart::xxx][0]-gfak*R_temp[Cart::z][Cart::xxx][1]);
R_temp[Cart::zzz][Cart::xxz][0]+=pma2*R_temp[Cart::zz][Cart::xxz][0]+wmp2*R_temp[Cart::zz][Cart::xxz][1]+1*rzeta*(R_temp[Cart::z][Cart::xxz][0]-gfak*R_temp[Cart::z][Cart::xxz][1])+0.5/_decay*1*R_temp[Cart::zz][Cart::xx][1];
R_temp[Cart::zzz][Cart::xzz][0]+=pma2*R_temp[Cart::zz][Cart::xzz][0]+wmp2*R_temp[Cart::zz][Cart::xzz][1]+1*rzeta*(R_temp[Cart::z][Cart::xzz][0]-gfak*R_temp[Cart::z][Cart::xzz][1])+0.5/_decay*2*R_temp[Cart::zz][Cart::xz][1];
R_temp[Cart::zzz][Cart::zzz][0]+=pma2*R_temp[Cart::zz][Cart::zzz][0]+wmp2*R_temp[Cart::zz][Cart::zzz][1]+1*rzeta*(R_temp[Cart::z][Cart::zzz][0]-gfak*R_temp[Cart::z][Cart::zzz][1])+0.5/_decay*3*R_temp[Cart::zz][Cart::zz][1];
}
//------------------------------------------------------

//Integral f - s - f - m1
if (_mmax >1 ){
if (_lmax_alpha>2 && _lmax_gamma>2){
R_temp[Cart::yyy][Cart::yyy][1]+=pma1*R_temp[Cart::yy][Cart::yyy][1]+wmp1*R_temp[Cart::yy][Cart::yyy][2]+1*rzeta*(R_temp[Cart::y][Cart::yyy][1]-gfak*R_temp[Cart::y][Cart::yyy][2])+0.5/_decay*3*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::yyy][Cart::xyy][1]+=pma1*R_temp[Cart::yy][Cart::xyy][1]+wmp1*R_temp[Cart::yy][Cart::xyy][2]+1*rzeta*(R_temp[Cart::y][Cart::xyy][1]-gfak*R_temp[Cart::y][Cart::xyy][2])+0.5/_decay*2*R_temp[Cart::yy][Cart::xy][2];
R_temp[Cart::yyy][Cart::yyz][1]+=pma1*R_temp[Cart::yy][Cart::yyz][1]+wmp1*R_temp[Cart::yy][Cart::yyz][2]+1*rzeta*(R_temp[Cart::y][Cart::yyz][1]-gfak*R_temp[Cart::y][Cart::yyz][2])+0.5/_decay*2*R_temp[Cart::yy][Cart::yz][2];
R_temp[Cart::yyy][Cart::xxy][1]+=pma1*R_temp[Cart::yy][Cart::xxy][1]+wmp1*R_temp[Cart::yy][Cart::xxy][2]+1*rzeta*(R_temp[Cart::y][Cart::xxy][1]-gfak*R_temp[Cart::y][Cart::xxy][2])+0.5/_decay*1*R_temp[Cart::yy][Cart::xx][2];
R_temp[Cart::yyy][Cart::xyz][1]+=pma1*R_temp[Cart::yy][Cart::xyz][1]+wmp1*R_temp[Cart::yy][Cart::xyz][2]+1*rzeta*(R_temp[Cart::y][Cart::xyz][1]-gfak*R_temp[Cart::y][Cart::xyz][2])+0.5/_decay*1*R_temp[Cart::yy][Cart::xz][2];
R_temp[Cart::yyy][Cart::yzz][1]+=pma1*R_temp[Cart::yy][Cart::yzz][1]+wmp1*R_temp[Cart::yy][Cart::yzz][2]+1*rzeta*(R_temp[Cart::y][Cart::yzz][1]-gfak*R_temp[Cart::y][Cart::yzz][2])+0.5/_decay*1*R_temp[Cart::yy][Cart::zz][2];
R_temp[Cart::yyy][Cart::xxx][1]+=pma1*R_temp[Cart::yy][Cart::xxx][1]+wmp1*R_temp[Cart::yy][Cart::xxx][2]+1*rzeta*(R_temp[Cart::y][Cart::xxx][1]-gfak*R_temp[Cart::y][Cart::xxx][2]);
R_temp[Cart::yyy][Cart::xxz][1]+=pma1*R_temp[Cart::yy][Cart::xxz][1]+wmp1*R_temp[Cart::yy][Cart::xxz][2]+1*rzeta*(R_temp[Cart::y][Cart::xxz][1]-gfak*R_temp[Cart::y][Cart::xxz][2]);
R_temp[Cart::yyy][Cart::xzz][1]+=pma1*R_temp[Cart::yy][Cart::xzz][1]+wmp1*R_temp[Cart::yy][Cart::xzz][2]+1*rzeta*(R_temp[Cart::y][Cart::xzz][1]-gfak*R_temp[Cart::y][Cart::xzz][2]);
R_temp[Cart::yyy][Cart::zzz][1]+=pma1*R_temp[Cart::yy][Cart::zzz][1]+wmp1*R_temp[Cart::yy][Cart::zzz][2]+1*rzeta*(R_temp[Cart::y][Cart::zzz][1]-gfak*R_temp[Cart::y][Cart::zzz][2]);
R_temp[Cart::xyy][Cart::yyy][1]+=pma0*R_temp[Cart::yy][Cart::yyy][1]+wmp0*R_temp[Cart::yy][Cart::yyy][2];
R_temp[Cart::xyy][Cart::xyy][1]+=pma0*R_temp[Cart::yy][Cart::xyy][1]+wmp0*R_temp[Cart::yy][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::xyy][Cart::yyz][1]+=pma0*R_temp[Cart::yy][Cart::yyz][1]+wmp0*R_temp[Cart::yy][Cart::yyz][2];
R_temp[Cart::xyy][Cart::xxy][1]+=pma0*R_temp[Cart::yy][Cart::xxy][1]+wmp0*R_temp[Cart::yy][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::yy][Cart::xy][2];
R_temp[Cart::xyy][Cart::xyz][1]+=pma0*R_temp[Cart::yy][Cart::xyz][1]+wmp0*R_temp[Cart::yy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::yz][2];
R_temp[Cart::xyy][Cart::yzz][1]+=pma0*R_temp[Cart::yy][Cart::yzz][1]+wmp0*R_temp[Cart::yy][Cart::yzz][2];
R_temp[Cart::xyy][Cart::xxx][1]+=pma0*R_temp[Cart::yy][Cart::xxx][1]+wmp0*R_temp[Cart::yy][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::yy][Cart::xx][2];
R_temp[Cart::xyy][Cart::xxz][1]+=pma0*R_temp[Cart::yy][Cart::xxz][1]+wmp0*R_temp[Cart::yy][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::yy][Cart::xz][2];
R_temp[Cart::xyy][Cart::xzz][1]+=pma0*R_temp[Cart::yy][Cart::xzz][1]+wmp0*R_temp[Cart::yy][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::zz][2];
R_temp[Cart::xyy][Cart::zzz][1]+=pma0*R_temp[Cart::yy][Cart::zzz][1]+wmp0*R_temp[Cart::yy][Cart::zzz][2];
R_temp[Cart::yyz][Cart::yyy][1]+=pma2*R_temp[Cart::yy][Cart::yyy][1]+wmp2*R_temp[Cart::yy][Cart::yyy][2];
R_temp[Cart::yyz][Cart::xyy][1]+=pma2*R_temp[Cart::yy][Cart::xyy][1]+wmp2*R_temp[Cart::yy][Cart::xyy][2];
R_temp[Cart::yyz][Cart::yyz][1]+=pma2*R_temp[Cart::yy][Cart::yyz][1]+wmp2*R_temp[Cart::yy][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::yy][2];
R_temp[Cart::yyz][Cart::xxy][1]+=pma2*R_temp[Cart::yy][Cart::xxy][1]+wmp2*R_temp[Cart::yy][Cart::xxy][2];
R_temp[Cart::yyz][Cart::xyz][1]+=pma2*R_temp[Cart::yy][Cart::xyz][1]+wmp2*R_temp[Cart::yy][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::xy][2];
R_temp[Cart::yyz][Cart::yzz][1]+=pma2*R_temp[Cart::yy][Cart::yzz][1]+wmp2*R_temp[Cart::yy][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::yy][Cart::yz][2];
R_temp[Cart::yyz][Cart::xxx][1]+=pma2*R_temp[Cart::yy][Cart::xxx][1]+wmp2*R_temp[Cart::yy][Cart::xxx][2];
R_temp[Cart::yyz][Cart::xxz][1]+=pma2*R_temp[Cart::yy][Cart::xxz][1]+wmp2*R_temp[Cart::yy][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::yy][Cart::xx][2];
R_temp[Cart::yyz][Cart::xzz][1]+=pma2*R_temp[Cart::yy][Cart::xzz][1]+wmp2*R_temp[Cart::yy][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::yy][Cart::xz][2];
R_temp[Cart::yyz][Cart::zzz][1]+=pma2*R_temp[Cart::yy][Cart::zzz][1]+wmp2*R_temp[Cart::yy][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::yy][Cart::zz][2];
R_temp[Cart::xxy][Cart::yyy][1]+=pma1*R_temp[Cart::xx][Cart::yyy][1]+wmp1*R_temp[Cart::xx][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::xx][Cart::yy][2];
R_temp[Cart::xxy][Cart::xyy][1]+=pma1*R_temp[Cart::xx][Cart::xyy][1]+wmp1*R_temp[Cart::xx][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::xx][Cart::xy][2];
R_temp[Cart::xxy][Cart::yyz][1]+=pma1*R_temp[Cart::xx][Cart::yyz][1]+wmp1*R_temp[Cart::xx][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::xx][Cart::yz][2];
R_temp[Cart::xxy][Cart::xxy][1]+=pma1*R_temp[Cart::xx][Cart::xxy][1]+wmp1*R_temp[Cart::xx][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xxy][Cart::xyz][1]+=pma1*R_temp[Cart::xx][Cart::xyz][1]+wmp1*R_temp[Cart::xx][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::xz][2];
R_temp[Cart::xxy][Cart::yzz][1]+=pma1*R_temp[Cart::xx][Cart::yzz][1]+wmp1*R_temp[Cart::xx][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::zz][2];
R_temp[Cart::xxy][Cart::xxx][1]+=pma1*R_temp[Cart::xx][Cart::xxx][1]+wmp1*R_temp[Cart::xx][Cart::xxx][2];
R_temp[Cart::xxy][Cart::xxz][1]+=pma1*R_temp[Cart::xx][Cart::xxz][1]+wmp1*R_temp[Cart::xx][Cart::xxz][2];
R_temp[Cart::xxy][Cart::xzz][1]+=pma1*R_temp[Cart::xx][Cart::xzz][1]+wmp1*R_temp[Cart::xx][Cart::xzz][2];
R_temp[Cart::xxy][Cart::zzz][1]+=pma1*R_temp[Cart::xx][Cart::zzz][1]+wmp1*R_temp[Cart::xx][Cart::zzz][2];
R_temp[Cart::xyz][Cart::yyy][1]+=pma0*R_temp[Cart::yz][Cart::yyy][1]+wmp0*R_temp[Cart::yz][Cart::yyy][2];
R_temp[Cart::xyz][Cart::xyy][1]+=pma0*R_temp[Cart::yz][Cart::xyy][1]+wmp0*R_temp[Cart::yz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::yz][Cart::yy][2];
R_temp[Cart::xyz][Cart::yyz][1]+=pma0*R_temp[Cart::yz][Cart::yyz][1]+wmp0*R_temp[Cart::yz][Cart::yyz][2];
R_temp[Cart::xyz][Cart::xxy][1]+=pma0*R_temp[Cart::yz][Cart::xxy][1]+wmp0*R_temp[Cart::yz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::yz][Cart::xy][2];
R_temp[Cart::xyz][Cart::xyz][1]+=pma0*R_temp[Cart::yz][Cart::xyz][1]+wmp0*R_temp[Cart::yz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::yz][Cart::yz][2];
R_temp[Cart::xyz][Cart::yzz][1]+=pma0*R_temp[Cart::yz][Cart::yzz][1]+wmp0*R_temp[Cart::yz][Cart::yzz][2];
R_temp[Cart::xyz][Cart::xxx][1]+=pma0*R_temp[Cart::yz][Cart::xxx][1]+wmp0*R_temp[Cart::yz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::yz][Cart::xx][2];
R_temp[Cart::xyz][Cart::xxz][1]+=pma0*R_temp[Cart::yz][Cart::xxz][1]+wmp0*R_temp[Cart::yz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::yz][Cart::xz][2];
R_temp[Cart::xyz][Cart::xzz][1]+=pma0*R_temp[Cart::yz][Cart::xzz][1]+wmp0*R_temp[Cart::yz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::yz][Cart::zz][2];
R_temp[Cart::xyz][Cart::zzz][1]+=pma0*R_temp[Cart::yz][Cart::zzz][1]+wmp0*R_temp[Cart::yz][Cart::zzz][2];
R_temp[Cart::yzz][Cart::yyy][1]+=pma1*R_temp[Cart::zz][Cart::yyy][1]+wmp1*R_temp[Cart::zz][Cart::yyy][2]+0.5/_decay*3*R_temp[Cart::zz][Cart::yy][2];
R_temp[Cart::yzz][Cart::xyy][1]+=pma1*R_temp[Cart::zz][Cart::xyy][1]+wmp1*R_temp[Cart::zz][Cart::xyy][2]+0.5/_decay*2*R_temp[Cart::zz][Cart::xy][2];
R_temp[Cart::yzz][Cart::yyz][1]+=pma1*R_temp[Cart::zz][Cart::yyz][1]+wmp1*R_temp[Cart::zz][Cart::yyz][2]+0.5/_decay*2*R_temp[Cart::zz][Cart::yz][2];
R_temp[Cart::yzz][Cart::xxy][1]+=pma1*R_temp[Cart::zz][Cart::xxy][1]+wmp1*R_temp[Cart::zz][Cart::xxy][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::xx][2];
R_temp[Cart::yzz][Cart::xyz][1]+=pma1*R_temp[Cart::zz][Cart::xyz][1]+wmp1*R_temp[Cart::zz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::xz][2];
R_temp[Cart::yzz][Cart::yzz][1]+=pma1*R_temp[Cart::zz][Cart::yzz][1]+wmp1*R_temp[Cart::zz][Cart::yzz][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::zz][2];
R_temp[Cart::yzz][Cart::xxx][1]+=pma1*R_temp[Cart::zz][Cart::xxx][1]+wmp1*R_temp[Cart::zz][Cart::xxx][2];
R_temp[Cart::yzz][Cart::xxz][1]+=pma1*R_temp[Cart::zz][Cart::xxz][1]+wmp1*R_temp[Cart::zz][Cart::xxz][2];
R_temp[Cart::yzz][Cart::xzz][1]+=pma1*R_temp[Cart::zz][Cart::xzz][1]+wmp1*R_temp[Cart::zz][Cart::xzz][2];
R_temp[Cart::yzz][Cart::zzz][1]+=pma1*R_temp[Cart::zz][Cart::zzz][1]+wmp1*R_temp[Cart::zz][Cart::zzz][2];
R_temp[Cart::xxx][Cart::yyy][1]+=pma0*R_temp[Cart::xx][Cart::yyy][1]+wmp0*R_temp[Cart::xx][Cart::yyy][2]+1*rzeta*(R_temp[Cart::x][Cart::yyy][1]-gfak*R_temp[Cart::x][Cart::yyy][2]);
R_temp[Cart::xxx][Cart::xyy][1]+=pma0*R_temp[Cart::xx][Cart::xyy][1]+wmp0*R_temp[Cart::xx][Cart::xyy][2]+1*rzeta*(R_temp[Cart::x][Cart::xyy][1]-gfak*R_temp[Cart::x][Cart::xyy][2])+0.5/_decay*1*R_temp[Cart::xx][Cart::yy][2];
R_temp[Cart::xxx][Cart::yyz][1]+=pma0*R_temp[Cart::xx][Cart::yyz][1]+wmp0*R_temp[Cart::xx][Cart::yyz][2]+1*rzeta*(R_temp[Cart::x][Cart::yyz][1]-gfak*R_temp[Cart::x][Cart::yyz][2]);
R_temp[Cart::xxx][Cart::xxy][1]+=pma0*R_temp[Cart::xx][Cart::xxy][1]+wmp0*R_temp[Cart::xx][Cart::xxy][2]+1*rzeta*(R_temp[Cart::x][Cart::xxy][1]-gfak*R_temp[Cart::x][Cart::xxy][2])+0.5/_decay*2*R_temp[Cart::xx][Cart::xy][2];
R_temp[Cart::xxx][Cart::xyz][1]+=pma0*R_temp[Cart::xx][Cart::xyz][1]+wmp0*R_temp[Cart::xx][Cart::xyz][2]+1*rzeta*(R_temp[Cart::x][Cart::xyz][1]-gfak*R_temp[Cart::x][Cart::xyz][2])+0.5/_decay*1*R_temp[Cart::xx][Cart::yz][2];
R_temp[Cart::xxx][Cart::yzz][1]+=pma0*R_temp[Cart::xx][Cart::yzz][1]+wmp0*R_temp[Cart::xx][Cart::yzz][2]+1*rzeta*(R_temp[Cart::x][Cart::yzz][1]-gfak*R_temp[Cart::x][Cart::yzz][2]);
R_temp[Cart::xxx][Cart::xxx][1]+=pma0*R_temp[Cart::xx][Cart::xxx][1]+wmp0*R_temp[Cart::xx][Cart::xxx][2]+1*rzeta*(R_temp[Cart::x][Cart::xxx][1]-gfak*R_temp[Cart::x][Cart::xxx][2])+0.5/_decay*3*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xxx][Cart::xxz][1]+=pma0*R_temp[Cart::xx][Cart::xxz][1]+wmp0*R_temp[Cart::xx][Cart::xxz][2]+1*rzeta*(R_temp[Cart::x][Cart::xxz][1]-gfak*R_temp[Cart::x][Cart::xxz][2])+0.5/_decay*2*R_temp[Cart::xx][Cart::xz][2];
R_temp[Cart::xxx][Cart::xzz][1]+=pma0*R_temp[Cart::xx][Cart::xzz][1]+wmp0*R_temp[Cart::xx][Cart::xzz][2]+1*rzeta*(R_temp[Cart::x][Cart::xzz][1]-gfak*R_temp[Cart::x][Cart::xzz][2])+0.5/_decay*1*R_temp[Cart::xx][Cart::zz][2];
R_temp[Cart::xxx][Cart::zzz][1]+=pma0*R_temp[Cart::xx][Cart::zzz][1]+wmp0*R_temp[Cart::xx][Cart::zzz][2]+1*rzeta*(R_temp[Cart::x][Cart::zzz][1]-gfak*R_temp[Cart::x][Cart::zzz][2]);
R_temp[Cart::xxz][Cart::yyy][1]+=pma2*R_temp[Cart::xx][Cart::yyy][1]+wmp2*R_temp[Cart::xx][Cart::yyy][2];
R_temp[Cart::xxz][Cart::xyy][1]+=pma2*R_temp[Cart::xx][Cart::xyy][1]+wmp2*R_temp[Cart::xx][Cart::xyy][2];
R_temp[Cart::xxz][Cart::yyz][1]+=pma2*R_temp[Cart::xx][Cart::yyz][1]+wmp2*R_temp[Cart::xx][Cart::yyz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::yy][2];
R_temp[Cart::xxz][Cart::xxy][1]+=pma2*R_temp[Cart::xx][Cart::xxy][1]+wmp2*R_temp[Cart::xx][Cart::xxy][2];
R_temp[Cart::xxz][Cart::xyz][1]+=pma2*R_temp[Cart::xx][Cart::xyz][1]+wmp2*R_temp[Cart::xx][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::xy][2];
R_temp[Cart::xxz][Cart::yzz][1]+=pma2*R_temp[Cart::xx][Cart::yzz][1]+wmp2*R_temp[Cart::xx][Cart::yzz][2]+0.5/_decay*2*R_temp[Cart::xx][Cart::yz][2];
R_temp[Cart::xxz][Cart::xxx][1]+=pma2*R_temp[Cart::xx][Cart::xxx][1]+wmp2*R_temp[Cart::xx][Cart::xxx][2];
R_temp[Cart::xxz][Cart::xxz][1]+=pma2*R_temp[Cart::xx][Cart::xxz][1]+wmp2*R_temp[Cart::xx][Cart::xxz][2]+0.5/_decay*1*R_temp[Cart::xx][Cart::xx][2];
R_temp[Cart::xxz][Cart::xzz][1]+=pma2*R_temp[Cart::xx][Cart::xzz][1]+wmp2*R_temp[Cart::xx][Cart::xzz][2]+0.5/_decay*2*R_temp[Cart::xx][Cart::xz][2];
R_temp[Cart::xxz][Cart::zzz][1]+=pma2*R_temp[Cart::xx][Cart::zzz][1]+wmp2*R_temp[Cart::xx][Cart::zzz][2]+0.5/_decay*3*R_temp[Cart::xx][Cart::zz][2];
R_temp[Cart::xzz][Cart::yyy][1]+=pma0*R_temp[Cart::zz][Cart::yyy][1]+wmp0*R_temp[Cart::zz][Cart::yyy][2];
R_temp[Cart::xzz][Cart::xyy][1]+=pma0*R_temp[Cart::zz][Cart::xyy][1]+wmp0*R_temp[Cart::zz][Cart::xyy][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::yy][2];
R_temp[Cart::xzz][Cart::yyz][1]+=pma0*R_temp[Cart::zz][Cart::yyz][1]+wmp0*R_temp[Cart::zz][Cart::yyz][2];
R_temp[Cart::xzz][Cart::xxy][1]+=pma0*R_temp[Cart::zz][Cart::xxy][1]+wmp0*R_temp[Cart::zz][Cart::xxy][2]+0.5/_decay*2*R_temp[Cart::zz][Cart::xy][2];
R_temp[Cart::xzz][Cart::xyz][1]+=pma0*R_temp[Cart::zz][Cart::xyz][1]+wmp0*R_temp[Cart::zz][Cart::xyz][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::yz][2];
R_temp[Cart::xzz][Cart::yzz][1]+=pma0*R_temp[Cart::zz][Cart::yzz][1]+wmp0*R_temp[Cart::zz][Cart::yzz][2];
R_temp[Cart::xzz][Cart::xxx][1]+=pma0*R_temp[Cart::zz][Cart::xxx][1]+wmp0*R_temp[Cart::zz][Cart::xxx][2]+0.5/_decay*3*R_temp[Cart::zz][Cart::xx][2];
R_temp[Cart::xzz][Cart::xxz][1]+=pma0*R_temp[Cart::zz][Cart::xxz][1]+wmp0*R_temp[Cart::zz][Cart::xxz][2]+0.5/_decay*2*R_temp[Cart::zz][Cart::xz][2];
R_temp[Cart::xzz][Cart::xzz][1]+=pma0*R_temp[Cart::zz][Cart::xzz][1]+wmp0*R_temp[Cart::zz][Cart::xzz][2]+0.5/_decay*1*R_temp[Cart::zz][Cart::zz][2];
R_temp[Cart::xzz][Cart::zzz][1]+=pma0*R_temp[Cart::zz][Cart::zzz][1]+wmp0*R_temp[Cart::zz][Cart::zzz][2];
R_temp[Cart::zzz][Cart::yyy][1]+=pma2*R_temp[Cart::zz][Cart::yyy][1]+wmp2*R_temp[Cart::zz][Cart::yyy][2]+1*rzeta*(R_temp[Cart::z][Cart::yyy][1]-gfak*R_temp[Cart::z][Cart::yyy][2]);
R_temp[Cart::zzz][Cart::xyy][1]+=pma2*R_temp[Cart::zz][Cart::xyy][1]+wmp2*R_temp[Cart::zz][Cart::xyy][2]+1*rzeta*(R_temp[Cart::z][Cart::xyy][1]-gfak*R_temp[Cart::z][Cart::xyy][2]);
R_temp[Cart::zzz][Cart::yyz][1]+=pma2*R_temp[Cart::zz][Cart::yyz][1]+wmp2*R_temp[Cart::zz][Cart::yyz][2]+1*rzeta*(R_temp[Cart::z][Cart::yyz][1]-gfak*R_temp[Cart::z][Cart::yyz][2])+0.5/_decay*1*R_temp[Cart::zz][Cart::yy][2];
R_temp[Cart::zzz][Cart::xxy][1]+=pma2*R_temp[Cart::zz][Cart::xxy][1]+wmp2*R_temp[Cart::zz][Cart::xxy][2]+1*rzeta*(R_temp[Cart::z][Cart::xxy][1]-gfak*R_temp[Cart::z][Cart::xxy][2]);
R_temp[Cart::zzz][Cart::xyz][1]+=pma2*R_temp[Cart::zz][Cart::xyz][1]+wmp2*R_temp[Cart::zz][Cart::xyz][2]+1*rzeta*(R_temp[Cart::z][Cart::xyz][1]-gfak*R_temp[Cart::z][Cart::xyz][2])+0.5/_decay*1*R_temp[Cart::zz][Cart::xy][2];
R_temp[Cart::zzz][Cart::yzz][1]+=pma2*R_temp[Cart::zz][Cart::yzz][1]+wmp2*R_temp[Cart::zz][Cart::yzz][2]+1*rzeta*(R_temp[Cart::z][Cart::yzz][1]-gfak*R_temp[Cart::z][Cart::yzz][2])+0.5/_decay*2*R_temp[Cart::zz][Cart::yz][2];
R_temp[Cart::zzz][Cart::xxx][1]+=pma2*R_temp[Cart::zz][Cart::xxx][1]+wmp2*R_temp[Cart::zz][Cart::xxx][2]+1*rzeta*(R_temp[Cart::z][Cart::xxx][1]-gfak*R_temp[Cart::z][Cart::xxx][2]);
R_temp[Cart::zzz][Cart::xxz][1]+=pma2*R_temp[Cart::zz][Cart::xxz][1]+wmp2*R_temp[Cart::zz][Cart::xxz][2]+1*rzeta*(R_temp[Cart::z][Cart::xxz][1]-gfak*R_temp[Cart::z][Cart::xxz][2])+0.5/_decay*1*R_temp[Cart::zz][Cart::xx][2];
R_temp[Cart::zzz][Cart::xzz][1]+=pma2*R_temp[Cart::zz][Cart::xzz][1]+wmp2*R_temp[Cart::zz][Cart::xzz][2]+1*rzeta*(R_temp[Cart::z][Cart::xzz][1]-gfak*R_temp[Cart::z][Cart::xzz][2])+0.5/_decay*2*R_temp[Cart::zz][Cart::xz][2];
R_temp[Cart::zzz][Cart::zzz][1]+=pma2*R_temp[Cart::zz][Cart::zzz][1]+wmp2*R_temp[Cart::zz][Cart::zzz][2]+1*rzeta*(R_temp[Cart::z][Cart::zzz][1]-gfak*R_temp[Cart::z][Cart::zzz][2])+0.5/_decay*3*R_temp[Cart::zz][Cart::zz][2];
}
}
//------------------------------------------------------

//Integral g - s - s - m0
if (_lmax_alpha>3){
R_temp[Cart::yyyy][Cart::s][0]+=pma1*R_temp[Cart::yyy][Cart::s][0]+wmp1*R_temp[Cart::yyy][Cart::s][1]+2*rzeta*(R_temp[Cart::yy][Cart::s][0]-gfak*R_temp[Cart::yy][Cart::s][1]);
R_temp[Cart::xyyy][Cart::s][0]+=pma0*R_temp[Cart::yyy][Cart::s][0]+wmp0*R_temp[Cart::yyy][Cart::s][1];
R_temp[Cart::yyyz][Cart::s][0]+=pma2*R_temp[Cart::yyy][Cart::s][0]+wmp2*R_temp[Cart::yyy][Cart::s][1];
R_temp[Cart::xxyy][Cart::s][0]+=pma0*R_temp[Cart::xyy][Cart::s][0]+wmp0*R_temp[Cart::xyy][Cart::s][1];
R_temp[Cart::xyyz][Cart::s][0]+=pma0*R_temp[Cart::yyz][Cart::s][0]+wmp0*R_temp[Cart::yyz][Cart::s][1];
R_temp[Cart::yyzz][Cart::s][0]+=pma1*R_temp[Cart::yzz][Cart::s][0]+wmp1*R_temp[Cart::yzz][Cart::s][1];
R_temp[Cart::xxxy][Cart::s][0]+=pma1*R_temp[Cart::xxx][Cart::s][0]+wmp1*R_temp[Cart::xxx][Cart::s][1];
R_temp[Cart::xxyz][Cart::s][0]+=pma1*R_temp[Cart::xxz][Cart::s][0]+wmp1*R_temp[Cart::xxz][Cart::s][1];
R_temp[Cart::xyzz][Cart::s][0]+=pma0*R_temp[Cart::yzz][Cart::s][0]+wmp0*R_temp[Cart::yzz][Cart::s][1];
R_temp[Cart::yzzz][Cart::s][0]+=pma1*R_temp[Cart::zzz][Cart::s][0]+wmp1*R_temp[Cart::zzz][Cart::s][1];
R_temp[Cart::xxxx][Cart::s][0]+=pma0*R_temp[Cart::xxx][Cart::s][0]+wmp0*R_temp[Cart::xxx][Cart::s][1]+2*rzeta*(R_temp[Cart::xx][Cart::s][0]-gfak*R_temp[Cart::xx][Cart::s][1]);
R_temp[Cart::xxxz][Cart::s][0]+=pma2*R_temp[Cart::xxx][Cart::s][0]+wmp2*R_temp[Cart::xxx][Cart::s][1];
R_temp[Cart::xxzz][Cart::s][0]+=pma0*R_temp[Cart::xzz][Cart::s][0]+wmp0*R_temp[Cart::xzz][Cart::s][1];
R_temp[Cart::xzzz][Cart::s][0]+=pma0*R_temp[Cart::zzz][Cart::s][0]+wmp0*R_temp[Cart::zzz][Cart::s][1];
R_temp[Cart::zzzz][Cart::s][0]+=pma2*R_temp[Cart::zzz][Cart::s][0]+wmp2*R_temp[Cart::zzz][Cart::s][1]+2*rzeta*(R_temp[Cart::zz][Cart::s][0]-gfak*R_temp[Cart::zz][Cart::s][1]);
}
//------------------------------------------------------

//Integral g - s - s - m1
if (_mmax >1 ){
if (_lmax_alpha>3){
R_temp[Cart::yyyy][Cart::s][1]+=pma1*R_temp[Cart::yyy][Cart::s][1]+wmp1*R_temp[Cart::yyy][Cart::s][2]+2*rzeta*(R_temp[Cart::yy][Cart::s][1]-gfak*R_temp[Cart::yy][Cart::s][2]);
R_temp[Cart::xyyy][Cart::s][1]+=pma0*R_temp[Cart::yyy][Cart::s][1]+wmp0*R_temp[Cart::yyy][Cart::s][2];
R_temp[Cart::yyyz][Cart::s][1]+=pma2*R_temp[Cart::yyy][Cart::s][1]+wmp2*R_temp[Cart::yyy][Cart::s][2];
R_temp[Cart::xxyy][Cart::s][1]+=pma0*R_temp[Cart::xyy][Cart::s][1]+wmp0*R_temp[Cart::xyy][Cart::s][2];
R_temp[Cart::xyyz][Cart::s][1]+=pma0*R_temp[Cart::yyz][Cart::s][1]+wmp0*R_temp[Cart::yyz][Cart::s][2];
R_temp[Cart::yyzz][Cart::s][1]+=pma1*R_temp[Cart::yzz][Cart::s][1]+wmp1*R_temp[Cart::yzz][Cart::s][2];
R_temp[Cart::xxxy][Cart::s][1]+=pma1*R_temp[Cart::xxx][Cart::s][1]+wmp1*R_temp[Cart::xxx][Cart::s][2];
R_temp[Cart::xxyz][Cart::s][1]+=pma1*R_temp[Cart::xxz][Cart::s][1]+wmp1*R_temp[Cart::xxz][Cart::s][2];
R_temp[Cart::xyzz][Cart::s][1]+=pma0*R_temp[Cart::yzz][Cart::s][1]+wmp0*R_temp[Cart::yzz][Cart::s][2];
R_temp[Cart::yzzz][Cart::s][1]+=pma1*R_temp[Cart::zzz][Cart::s][1]+wmp1*R_temp[Cart::zzz][Cart::s][2];
R_temp[Cart::xxxx][Cart::s][1]+=pma0*R_temp[Cart::xxx][Cart::s][1]+wmp0*R_temp[Cart::xxx][Cart::s][2]+2*rzeta*(R_temp[Cart::xx][Cart::s][1]-gfak*R_temp[Cart::xx][Cart::s][2]);
R_temp[Cart::xxxz][Cart::s][1]+=pma2*R_temp[Cart::xxx][Cart::s][1]+wmp2*R_temp[Cart::xxx][Cart::s][2];
R_temp[Cart::xxzz][Cart::s][1]+=pma0*R_temp[Cart::xzz][Cart::s][1]+wmp0*R_temp[Cart::xzz][Cart::s][2];
R_temp[Cart::xzzz][Cart::s][1]+=pma0*R_temp[Cart::zzz][Cart::s][1]+wmp0*R_temp[Cart::zzz][Cart::s][2];
R_temp[Cart::zzzz][Cart::s][1]+=pma2*R_temp[Cart::zzz][Cart::s][1]+wmp2*R_temp[Cart::zzz][Cart::s][2]+2*rzeta*(R_temp[Cart::zz][Cart::s][1]-gfak*R_temp[Cart::zz][Cart::s][2]);
}
}
//------------------------------------------------------

//Integral g - s - s - m2
if (_mmax >2 ){
if (_lmax_alpha>3){
R_temp[Cart::yyyy][Cart::s][2]+=pma1*R_temp[Cart::yyy][Cart::s][2]+wmp1*R_temp[Cart::yyy][Cart::s][3]+2*rzeta*(R_temp[Cart::yy][Cart::s][2]-gfak*R_temp[Cart::yy][Cart::s][3]);
R_temp[Cart::xyyy][Cart::s][2]+=pma0*R_temp[Cart::yyy][Cart::s][2]+wmp0*R_temp[Cart::yyy][Cart::s][3];
R_temp[Cart::yyyz][Cart::s][2]+=pma2*R_temp[Cart::yyy][Cart::s][2]+wmp2*R_temp[Cart::yyy][Cart::s][3];
R_temp[Cart::xxyy][Cart::s][2]+=pma0*R_temp[Cart::xyy][Cart::s][2]+wmp0*R_temp[Cart::xyy][Cart::s][3];
R_temp[Cart::xyyz][Cart::s][2]+=pma0*R_temp[Cart::yyz][Cart::s][2]+wmp0*R_temp[Cart::yyz][Cart::s][3];
R_temp[Cart::yyzz][Cart::s][2]+=pma1*R_temp[Cart::yzz][Cart::s][2]+wmp1*R_temp[Cart::yzz][Cart::s][3];
R_temp[Cart::xxxy][Cart::s][2]+=pma1*R_temp[Cart::xxx][Cart::s][2]+wmp1*R_temp[Cart::xxx][Cart::s][3];
R_temp[Cart::xxyz][Cart::s][2]+=pma1*R_temp[Cart::xxz][Cart::s][2]+wmp1*R_temp[Cart::xxz][Cart::s][3];
R_temp[Cart::xyzz][Cart::s][2]+=pma0*R_temp[Cart::yzz][Cart::s][2]+wmp0*R_temp[Cart::yzz][Cart::s][3];
R_temp[Cart::yzzz][Cart::s][2]+=pma1*R_temp[Cart::zzz][Cart::s][2]+wmp1*R_temp[Cart::zzz][Cart::s][3];
R_temp[Cart::xxxx][Cart::s][2]+=pma0*R_temp[Cart::xxx][Cart::s][2]+wmp0*R_temp[Cart::xxx][Cart::s][3]+2*rzeta*(R_temp[Cart::xx][Cart::s][2]-gfak*R_temp[Cart::xx][Cart::s][3]);
R_temp[Cart::xxxz][Cart::s][2]+=pma2*R_temp[Cart::xxx][Cart::s][2]+wmp2*R_temp[Cart::xxx][Cart::s][3];
R_temp[Cart::xxzz][Cart::s][2]+=pma0*R_temp[Cart::xzz][Cart::s][2]+wmp0*R_temp[Cart::xzz][Cart::s][3];
R_temp[Cart::xzzz][Cart::s][2]+=pma0*R_temp[Cart::zzz][Cart::s][2]+wmp0*R_temp[Cart::zzz][Cart::s][3];
R_temp[Cart::zzzz][Cart::s][2]+=pma2*R_temp[Cart::zzz][Cart::s][2]+wmp2*R_temp[Cart::zzz][Cart::s][3]+2*rzeta*(R_temp[Cart::zz][Cart::s][2]-gfak*R_temp[Cart::zz][Cart::s][3]);
}
}
//------------------------------------------------------

//Integral g - s - s - m3
if (_mmax >3 ){
if (_lmax_alpha>3){
R_temp[Cart::yyyy][Cart::s][3]+=pma1*R_temp[Cart::yyy][Cart::s][3]+wmp1*R_temp[Cart::yyy][Cart::s][4]+2*rzeta*(R_temp[Cart::yy][Cart::s][3]-gfak*R_temp[Cart::yy][Cart::s][4]);
R_temp[Cart::xyyy][Cart::s][3]+=pma0*R_temp[Cart::yyy][Cart::s][3]+wmp0*R_temp[Cart::yyy][Cart::s][4];
R_temp[Cart::yyyz][Cart::s][3]+=pma2*R_temp[Cart::yyy][Cart::s][3]+wmp2*R_temp[Cart::yyy][Cart::s][4];
R_temp[Cart::xxyy][Cart::s][3]+=pma0*R_temp[Cart::xyy][Cart::s][3]+wmp0*R_temp[Cart::xyy][Cart::s][4];
R_temp[Cart::xyyz][Cart::s][3]+=pma0*R_temp[Cart::yyz][Cart::s][3]+wmp0*R_temp[Cart::yyz][Cart::s][4];
R_temp[Cart::yyzz][Cart::s][3]+=pma1*R_temp[Cart::yzz][Cart::s][3]+wmp1*R_temp[Cart::yzz][Cart::s][4];
R_temp[Cart::xxxy][Cart::s][3]+=pma1*R_temp[Cart::xxx][Cart::s][3]+wmp1*R_temp[Cart::xxx][Cart::s][4];
R_temp[Cart::xxyz][Cart::s][3]+=pma1*R_temp[Cart::xxz][Cart::s][3]+wmp1*R_temp[Cart::xxz][Cart::s][4];
R_temp[Cart::xyzz][Cart::s][3]+=pma0*R_temp[Cart::yzz][Cart::s][3]+wmp0*R_temp[Cart::yzz][Cart::s][4];
R_temp[Cart::yzzz][Cart::s][3]+=pma1*R_temp[Cart::zzz][Cart::s][3]+wmp1*R_temp[Cart::zzz][Cart::s][4];
R_temp[Cart::xxxx][Cart::s][3]+=pma0*R_temp[Cart::xxx][Cart::s][3]+wmp0*R_temp[Cart::xxx][Cart::s][4]+2*rzeta*(R_temp[Cart::xx][Cart::s][3]-gfak*R_temp[Cart::xx][Cart::s][4]);
R_temp[Cart::xxxz][Cart::s][3]+=pma2*R_temp[Cart::xxx][Cart::s][3]+wmp2*R_temp[Cart::xxx][Cart::s][4];
R_temp[Cart::xxzz][Cart::s][3]+=pma0*R_temp[Cart::xzz][Cart::s][3]+wmp0*R_temp[Cart::xzz][Cart::s][4];
R_temp[Cart::xzzz][Cart::s][3]+=pma0*R_temp[Cart::zzz][Cart::s][3]+wmp0*R_temp[Cart::zzz][Cart::s][4];
R_temp[Cart::zzzz][Cart::s][3]+=pma2*R_temp[Cart::zzz][Cart::s][3]+wmp2*R_temp[Cart::zzz][Cart::s][4]+2*rzeta*(R_temp[Cart::zz][Cart::s][3]-gfak*R_temp[Cart::zz][Cart::s][4]);
}
}
//------------------------------------------------------

//Integral g - s - p - m0
if (_lmax_alpha>3 && _lmax_gamma>0){
R_temp[Cart::yyyy][Cart::y][0]+=pma1*R_temp[Cart::yyy][Cart::y][0]+wmp1*R_temp[Cart::yyy][Cart::y][1]+2*rzeta*(R_temp[Cart::yy][Cart::y][0]-gfak*R_temp[Cart::yy][Cart::y][1])+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][1];
R_temp[Cart::yyyy][Cart::x][0]+=pma1*R_temp[Cart::yyy][Cart::x][0]+wmp1*R_temp[Cart::yyy][Cart::x][1]+2*rzeta*(R_temp[Cart::yy][Cart::x][0]-gfak*R_temp[Cart::yy][Cart::x][1]);
R_temp[Cart::yyyy][Cart::z][0]+=pma1*R_temp[Cart::yyy][Cart::z][0]+wmp1*R_temp[Cart::yyy][Cart::z][1]+2*rzeta*(R_temp[Cart::yy][Cart::z][0]-gfak*R_temp[Cart::yy][Cart::z][1]);
R_temp[Cart::xyyy][Cart::y][0]+=pma0*R_temp[Cart::yyy][Cart::y][0]+wmp0*R_temp[Cart::yyy][Cart::y][1];
R_temp[Cart::xyyy][Cart::x][0]+=pma0*R_temp[Cart::yyy][Cart::x][0]+wmp0*R_temp[Cart::yyy][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][1];
R_temp[Cart::xyyy][Cart::z][0]+=pma0*R_temp[Cart::yyy][Cart::z][0]+wmp0*R_temp[Cart::yyy][Cart::z][1];
R_temp[Cart::yyyz][Cart::y][0]+=pma2*R_temp[Cart::yyy][Cart::y][0]+wmp2*R_temp[Cart::yyy][Cart::y][1];
R_temp[Cart::yyyz][Cart::x][0]+=pma2*R_temp[Cart::yyy][Cart::x][0]+wmp2*R_temp[Cart::yyy][Cart::x][1];
R_temp[Cart::yyyz][Cart::z][0]+=pma2*R_temp[Cart::yyy][Cart::z][0]+wmp2*R_temp[Cart::yyy][Cart::z][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][1];
R_temp[Cart::xxyy][Cart::y][0]+=pma0*R_temp[Cart::xyy][Cart::y][0]+wmp0*R_temp[Cart::xyy][Cart::y][1];
R_temp[Cart::xxyy][Cart::x][0]+=pma0*R_temp[Cart::xyy][Cart::x][0]+wmp0*R_temp[Cart::xyy][Cart::x][1]+0.5/_decay*1*R_temp[Cart::xyy][Cart::s][1];
R_temp[Cart::xxyy][Cart::z][0]+=pma0*R_temp[Cart::xyy][Cart::z][0]+wmp0*R_temp[Cart::xyy][Cart::z][1];
R_temp[Cart::xyyz][Cart::y][0]+=pma0*R_temp[Cart::yyz][Cart::y][0]+wmp0*R_temp[Cart::yyz][Cart::y][1];
R_temp[Cart::xyyz][Cart::x][0]+=pma0*R_temp[Cart::yyz][Cart::x][0]+wmp0*R_temp[Cart::yyz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yyz][Cart::s][1];
R_temp[Cart::xyyz][Cart::z][0]+=pma0*R_temp[Cart::yyz][Cart::z][0]+wmp0*R_temp[Cart::yyz][Cart::z][1];
R_temp[Cart::yyzz][Cart::y][0]+=pma1*R_temp[Cart::yzz][Cart::y][0]+wmp1*R_temp[Cart::yzz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][1];
R_temp[Cart::yyzz][Cart::x][0]+=pma1*R_temp[Cart::yzz][Cart::x][0]+wmp1*R_temp[Cart::yzz][Cart::x][1];
R_temp[Cart::yyzz][Cart::z][0]+=pma1*R_temp[Cart::yzz][Cart::z][0]+wmp1*R_temp[Cart::yzz][Cart::z][1];
R_temp[Cart::xxxy][Cart::y][0]+=pma1*R_temp[Cart::xxx][Cart::y][0]+wmp1*R_temp[Cart::xxx][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][1];
R_temp[Cart::xxxy][Cart::x][0]+=pma1*R_temp[Cart::xxx][Cart::x][0]+wmp1*R_temp[Cart::xxx][Cart::x][1];
R_temp[Cart::xxxy][Cart::z][0]+=pma1*R_temp[Cart::xxx][Cart::z][0]+wmp1*R_temp[Cart::xxx][Cart::z][1];
R_temp[Cart::xxyz][Cart::y][0]+=pma1*R_temp[Cart::xxz][Cart::y][0]+wmp1*R_temp[Cart::xxz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::xxz][Cart::s][1];
R_temp[Cart::xxyz][Cart::x][0]+=pma1*R_temp[Cart::xxz][Cart::x][0]+wmp1*R_temp[Cart::xxz][Cart::x][1];
R_temp[Cart::xxyz][Cart::z][0]+=pma1*R_temp[Cart::xxz][Cart::z][0]+wmp1*R_temp[Cart::xxz][Cart::z][1];
R_temp[Cart::xyzz][Cart::y][0]+=pma0*R_temp[Cart::yzz][Cart::y][0]+wmp0*R_temp[Cart::yzz][Cart::y][1];
R_temp[Cart::xyzz][Cart::x][0]+=pma0*R_temp[Cart::yzz][Cart::x][0]+wmp0*R_temp[Cart::yzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][1];
R_temp[Cart::xyzz][Cart::z][0]+=pma0*R_temp[Cart::yzz][Cart::z][0]+wmp0*R_temp[Cart::yzz][Cart::z][1];
R_temp[Cart::yzzz][Cart::y][0]+=pma1*R_temp[Cart::zzz][Cart::y][0]+wmp1*R_temp[Cart::zzz][Cart::y][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][1];
R_temp[Cart::yzzz][Cart::x][0]+=pma1*R_temp[Cart::zzz][Cart::x][0]+wmp1*R_temp[Cart::zzz][Cart::x][1];
R_temp[Cart::yzzz][Cart::z][0]+=pma1*R_temp[Cart::zzz][Cart::z][0]+wmp1*R_temp[Cart::zzz][Cart::z][1];
R_temp[Cart::xxxx][Cart::y][0]+=pma0*R_temp[Cart::xxx][Cart::y][0]+wmp0*R_temp[Cart::xxx][Cart::y][1]+2*rzeta*(R_temp[Cart::xx][Cart::y][0]-gfak*R_temp[Cart::xx][Cart::y][1]);
R_temp[Cart::xxxx][Cart::x][0]+=pma0*R_temp[Cart::xxx][Cart::x][0]+wmp0*R_temp[Cart::xxx][Cart::x][1]+2*rzeta*(R_temp[Cart::xx][Cart::x][0]-gfak*R_temp[Cart::xx][Cart::x][1])+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][1];
R_temp[Cart::xxxx][Cart::z][0]+=pma0*R_temp[Cart::xxx][Cart::z][0]+wmp0*R_temp[Cart::xxx][Cart::z][1]+2*rzeta*(R_temp[Cart::xx][Cart::z][0]-gfak*R_temp[Cart::xx][Cart::z][1]);
R_temp[Cart::xxxz][Cart::y][0]+=pma2*R_temp[Cart::xxx][Cart::y][0]+wmp2*R_temp[Cart::xxx][Cart::y][1];
R_temp[Cart::xxxz][Cart::x][0]+=pma2*R_temp[Cart::xxx][Cart::x][0]+wmp2*R_temp[Cart::xxx][Cart::x][1];
R_temp[Cart::xxxz][Cart::z][0]+=pma2*R_temp[Cart::xxx][Cart::z][0]+wmp2*R_temp[Cart::xxx][Cart::z][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][1];
R_temp[Cart::xxzz][Cart::y][0]+=pma0*R_temp[Cart::xzz][Cart::y][0]+wmp0*R_temp[Cart::xzz][Cart::y][1];
R_temp[Cart::xxzz][Cart::x][0]+=pma0*R_temp[Cart::xzz][Cart::x][0]+wmp0*R_temp[Cart::xzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::xzz][Cart::s][1];
R_temp[Cart::xxzz][Cart::z][0]+=pma0*R_temp[Cart::xzz][Cart::z][0]+wmp0*R_temp[Cart::xzz][Cart::z][1];
R_temp[Cart::xzzz][Cart::y][0]+=pma0*R_temp[Cart::zzz][Cart::y][0]+wmp0*R_temp[Cart::zzz][Cart::y][1];
R_temp[Cart::xzzz][Cart::x][0]+=pma0*R_temp[Cart::zzz][Cart::x][0]+wmp0*R_temp[Cart::zzz][Cart::x][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][1];
R_temp[Cart::xzzz][Cart::z][0]+=pma0*R_temp[Cart::zzz][Cart::z][0]+wmp0*R_temp[Cart::zzz][Cart::z][1];
R_temp[Cart::zzzz][Cart::y][0]+=pma2*R_temp[Cart::zzz][Cart::y][0]+wmp2*R_temp[Cart::zzz][Cart::y][1]+2*rzeta*(R_temp[Cart::zz][Cart::y][0]-gfak*R_temp[Cart::zz][Cart::y][1]);
R_temp[Cart::zzzz][Cart::x][0]+=pma2*R_temp[Cart::zzz][Cart::x][0]+wmp2*R_temp[Cart::zzz][Cart::x][1]+2*rzeta*(R_temp[Cart::zz][Cart::x][0]-gfak*R_temp[Cart::zz][Cart::x][1]);
R_temp[Cart::zzzz][Cart::z][0]+=pma2*R_temp[Cart::zzz][Cart::z][0]+wmp2*R_temp[Cart::zzz][Cart::z][1]+2*rzeta*(R_temp[Cart::zz][Cart::z][0]-gfak*R_temp[Cart::zz][Cart::z][1])+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][1];
}
//------------------------------------------------------

//Integral g - s - p - m1
if (_mmax >1 ){
if (_lmax_alpha>3 && _lmax_gamma>0){
R_temp[Cart::yyyy][Cart::y][1]+=pma1*R_temp[Cart::yyy][Cart::y][1]+wmp1*R_temp[Cart::yyy][Cart::y][2]+2*rzeta*(R_temp[Cart::yy][Cart::y][1]-gfak*R_temp[Cart::yy][Cart::y][2])+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][2];
R_temp[Cart::yyyy][Cart::x][1]+=pma1*R_temp[Cart::yyy][Cart::x][1]+wmp1*R_temp[Cart::yyy][Cart::x][2]+2*rzeta*(R_temp[Cart::yy][Cart::x][1]-gfak*R_temp[Cart::yy][Cart::x][2]);
R_temp[Cart::yyyy][Cart::z][1]+=pma1*R_temp[Cart::yyy][Cart::z][1]+wmp1*R_temp[Cart::yyy][Cart::z][2]+2*rzeta*(R_temp[Cart::yy][Cart::z][1]-gfak*R_temp[Cart::yy][Cart::z][2]);
R_temp[Cart::xyyy][Cart::y][1]+=pma0*R_temp[Cart::yyy][Cart::y][1]+wmp0*R_temp[Cart::yyy][Cart::y][2];
R_temp[Cart::xyyy][Cart::x][1]+=pma0*R_temp[Cart::yyy][Cart::x][1]+wmp0*R_temp[Cart::yyy][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][2];
R_temp[Cart::xyyy][Cart::z][1]+=pma0*R_temp[Cart::yyy][Cart::z][1]+wmp0*R_temp[Cart::yyy][Cart::z][2];
R_temp[Cart::yyyz][Cart::y][1]+=pma2*R_temp[Cart::yyy][Cart::y][1]+wmp2*R_temp[Cart::yyy][Cart::y][2];
R_temp[Cart::yyyz][Cart::x][1]+=pma2*R_temp[Cart::yyy][Cart::x][1]+wmp2*R_temp[Cart::yyy][Cart::x][2];
R_temp[Cart::yyyz][Cart::z][1]+=pma2*R_temp[Cart::yyy][Cart::z][1]+wmp2*R_temp[Cart::yyy][Cart::z][2]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][2];
R_temp[Cart::xxyy][Cart::y][1]+=pma0*R_temp[Cart::xyy][Cart::y][1]+wmp0*R_temp[Cart::xyy][Cart::y][2];
R_temp[Cart::xxyy][Cart::x][1]+=pma0*R_temp[Cart::xyy][Cart::x][1]+wmp0*R_temp[Cart::xyy][Cart::x][2]+0.5/_decay*1*R_temp[Cart::xyy][Cart::s][2];
R_temp[Cart::xxyy][Cart::z][1]+=pma0*R_temp[Cart::xyy][Cart::z][1]+wmp0*R_temp[Cart::xyy][Cart::z][2];
R_temp[Cart::xyyz][Cart::y][1]+=pma0*R_temp[Cart::yyz][Cart::y][1]+wmp0*R_temp[Cart::yyz][Cart::y][2];
R_temp[Cart::xyyz][Cart::x][1]+=pma0*R_temp[Cart::yyz][Cart::x][1]+wmp0*R_temp[Cart::yyz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yyz][Cart::s][2];
R_temp[Cart::xyyz][Cart::z][1]+=pma0*R_temp[Cart::yyz][Cart::z][1]+wmp0*R_temp[Cart::yyz][Cart::z][2];
R_temp[Cart::yyzz][Cart::y][1]+=pma1*R_temp[Cart::yzz][Cart::y][1]+wmp1*R_temp[Cart::yzz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][2];
R_temp[Cart::yyzz][Cart::x][1]+=pma1*R_temp[Cart::yzz][Cart::x][1]+wmp1*R_temp[Cart::yzz][Cart::x][2];
R_temp[Cart::yyzz][Cart::z][1]+=pma1*R_temp[Cart::yzz][Cart::z][1]+wmp1*R_temp[Cart::yzz][Cart::z][2];
R_temp[Cart::xxxy][Cart::y][1]+=pma1*R_temp[Cart::xxx][Cart::y][1]+wmp1*R_temp[Cart::xxx][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][2];
R_temp[Cart::xxxy][Cart::x][1]+=pma1*R_temp[Cart::xxx][Cart::x][1]+wmp1*R_temp[Cart::xxx][Cart::x][2];
R_temp[Cart::xxxy][Cart::z][1]+=pma1*R_temp[Cart::xxx][Cart::z][1]+wmp1*R_temp[Cart::xxx][Cart::z][2];
R_temp[Cart::xxyz][Cart::y][1]+=pma1*R_temp[Cart::xxz][Cart::y][1]+wmp1*R_temp[Cart::xxz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::xxz][Cart::s][2];
R_temp[Cart::xxyz][Cart::x][1]+=pma1*R_temp[Cart::xxz][Cart::x][1]+wmp1*R_temp[Cart::xxz][Cart::x][2];
R_temp[Cart::xxyz][Cart::z][1]+=pma1*R_temp[Cart::xxz][Cart::z][1]+wmp1*R_temp[Cart::xxz][Cart::z][2];
R_temp[Cart::xyzz][Cart::y][1]+=pma0*R_temp[Cart::yzz][Cart::y][1]+wmp0*R_temp[Cart::yzz][Cart::y][2];
R_temp[Cart::xyzz][Cart::x][1]+=pma0*R_temp[Cart::yzz][Cart::x][1]+wmp0*R_temp[Cart::yzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][2];
R_temp[Cart::xyzz][Cart::z][1]+=pma0*R_temp[Cart::yzz][Cart::z][1]+wmp0*R_temp[Cart::yzz][Cart::z][2];
R_temp[Cart::yzzz][Cart::y][1]+=pma1*R_temp[Cart::zzz][Cart::y][1]+wmp1*R_temp[Cart::zzz][Cart::y][2]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][2];
R_temp[Cart::yzzz][Cart::x][1]+=pma1*R_temp[Cart::zzz][Cart::x][1]+wmp1*R_temp[Cart::zzz][Cart::x][2];
R_temp[Cart::yzzz][Cart::z][1]+=pma1*R_temp[Cart::zzz][Cart::z][1]+wmp1*R_temp[Cart::zzz][Cart::z][2];
R_temp[Cart::xxxx][Cart::y][1]+=pma0*R_temp[Cart::xxx][Cart::y][1]+wmp0*R_temp[Cart::xxx][Cart::y][2]+2*rzeta*(R_temp[Cart::xx][Cart::y][1]-gfak*R_temp[Cart::xx][Cart::y][2]);
R_temp[Cart::xxxx][Cart::x][1]+=pma0*R_temp[Cart::xxx][Cart::x][1]+wmp0*R_temp[Cart::xxx][Cart::x][2]+2*rzeta*(R_temp[Cart::xx][Cart::x][1]-gfak*R_temp[Cart::xx][Cart::x][2])+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][2];
R_temp[Cart::xxxx][Cart::z][1]+=pma0*R_temp[Cart::xxx][Cart::z][1]+wmp0*R_temp[Cart::xxx][Cart::z][2]+2*rzeta*(R_temp[Cart::xx][Cart::z][1]-gfak*R_temp[Cart::xx][Cart::z][2]);
R_temp[Cart::xxxz][Cart::y][1]+=pma2*R_temp[Cart::xxx][Cart::y][1]+wmp2*R_temp[Cart::xxx][Cart::y][2];
R_temp[Cart::xxxz][Cart::x][1]+=pma2*R_temp[Cart::xxx][Cart::x][1]+wmp2*R_temp[Cart::xxx][Cart::x][2];
R_temp[Cart::xxxz][Cart::z][1]+=pma2*R_temp[Cart::xxx][Cart::z][1]+wmp2*R_temp[Cart::xxx][Cart::z][2]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][2];
R_temp[Cart::xxzz][Cart::y][1]+=pma0*R_temp[Cart::xzz][Cart::y][1]+wmp0*R_temp[Cart::xzz][Cart::y][2];
R_temp[Cart::xxzz][Cart::x][1]+=pma0*R_temp[Cart::xzz][Cart::x][1]+wmp0*R_temp[Cart::xzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::xzz][Cart::s][2];
R_temp[Cart::xxzz][Cart::z][1]+=pma0*R_temp[Cart::xzz][Cart::z][1]+wmp0*R_temp[Cart::xzz][Cart::z][2];
R_temp[Cart::xzzz][Cart::y][1]+=pma0*R_temp[Cart::zzz][Cart::y][1]+wmp0*R_temp[Cart::zzz][Cart::y][2];
R_temp[Cart::xzzz][Cart::x][1]+=pma0*R_temp[Cart::zzz][Cart::x][1]+wmp0*R_temp[Cart::zzz][Cart::x][2]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][2];
R_temp[Cart::xzzz][Cart::z][1]+=pma0*R_temp[Cart::zzz][Cart::z][1]+wmp0*R_temp[Cart::zzz][Cart::z][2];
R_temp[Cart::zzzz][Cart::y][1]+=pma2*R_temp[Cart::zzz][Cart::y][1]+wmp2*R_temp[Cart::zzz][Cart::y][2]+2*rzeta*(R_temp[Cart::zz][Cart::y][1]-gfak*R_temp[Cart::zz][Cart::y][2]);
R_temp[Cart::zzzz][Cart::x][1]+=pma2*R_temp[Cart::zzz][Cart::x][1]+wmp2*R_temp[Cart::zzz][Cart::x][2]+2*rzeta*(R_temp[Cart::zz][Cart::x][1]-gfak*R_temp[Cart::zz][Cart::x][2]);
R_temp[Cart::zzzz][Cart::z][1]+=pma2*R_temp[Cart::zzz][Cart::z][1]+wmp2*R_temp[Cart::zzz][Cart::z][2]+2*rzeta*(R_temp[Cart::zz][Cart::z][1]-gfak*R_temp[Cart::zz][Cart::z][2])+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][2];
}
}
//------------------------------------------------------

//Integral g - s - p - m2
if (_mmax >2 ){
if (_lmax_alpha>3 && _lmax_gamma>0){
R_temp[Cart::yyyy][Cart::y][2]+=pma1*R_temp[Cart::yyy][Cart::y][2]+wmp1*R_temp[Cart::yyy][Cart::y][3]+2*rzeta*(R_temp[Cart::yy][Cart::y][2]-gfak*R_temp[Cart::yy][Cart::y][3])+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][3];
R_temp[Cart::yyyy][Cart::x][2]+=pma1*R_temp[Cart::yyy][Cart::x][2]+wmp1*R_temp[Cart::yyy][Cart::x][3]+2*rzeta*(R_temp[Cart::yy][Cart::x][2]-gfak*R_temp[Cart::yy][Cart::x][3]);
R_temp[Cart::yyyy][Cart::z][2]+=pma1*R_temp[Cart::yyy][Cart::z][2]+wmp1*R_temp[Cart::yyy][Cart::z][3]+2*rzeta*(R_temp[Cart::yy][Cart::z][2]-gfak*R_temp[Cart::yy][Cart::z][3]);
R_temp[Cart::xyyy][Cart::y][2]+=pma0*R_temp[Cart::yyy][Cart::y][2]+wmp0*R_temp[Cart::yyy][Cart::y][3];
R_temp[Cart::xyyy][Cart::x][2]+=pma0*R_temp[Cart::yyy][Cart::x][2]+wmp0*R_temp[Cart::yyy][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][3];
R_temp[Cart::xyyy][Cart::z][2]+=pma0*R_temp[Cart::yyy][Cart::z][2]+wmp0*R_temp[Cart::yyy][Cart::z][3];
R_temp[Cart::yyyz][Cart::y][2]+=pma2*R_temp[Cart::yyy][Cart::y][2]+wmp2*R_temp[Cart::yyy][Cart::y][3];
R_temp[Cart::yyyz][Cart::x][2]+=pma2*R_temp[Cart::yyy][Cart::x][2]+wmp2*R_temp[Cart::yyy][Cart::x][3];
R_temp[Cart::yyyz][Cart::z][2]+=pma2*R_temp[Cart::yyy][Cart::z][2]+wmp2*R_temp[Cart::yyy][Cart::z][3]+0.5/_decay*1*R_temp[Cart::yyy][Cart::s][3];
R_temp[Cart::xxyy][Cart::y][2]+=pma0*R_temp[Cart::xyy][Cart::y][2]+wmp0*R_temp[Cart::xyy][Cart::y][3];
R_temp[Cart::xxyy][Cart::x][2]+=pma0*R_temp[Cart::xyy][Cart::x][2]+wmp0*R_temp[Cart::xyy][Cart::x][3]+0.5/_decay*1*R_temp[Cart::xyy][Cart::s][3];
R_temp[Cart::xxyy][Cart::z][2]+=pma0*R_temp[Cart::xyy][Cart::z][2]+wmp0*R_temp[Cart::xyy][Cart::z][3];
R_temp[Cart::xyyz][Cart::y][2]+=pma0*R_temp[Cart::yyz][Cart::y][2]+wmp0*R_temp[Cart::yyz][Cart::y][3];
R_temp[Cart::xyyz][Cart::x][2]+=pma0*R_temp[Cart::yyz][Cart::x][2]+wmp0*R_temp[Cart::yyz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yyz][Cart::s][3];
R_temp[Cart::xyyz][Cart::z][2]+=pma0*R_temp[Cart::yyz][Cart::z][2]+wmp0*R_temp[Cart::yyz][Cart::z][3];
R_temp[Cart::yyzz][Cart::y][2]+=pma1*R_temp[Cart::yzz][Cart::y][2]+wmp1*R_temp[Cart::yzz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][3];
R_temp[Cart::yyzz][Cart::x][2]+=pma1*R_temp[Cart::yzz][Cart::x][2]+wmp1*R_temp[Cart::yzz][Cart::x][3];
R_temp[Cart::yyzz][Cart::z][2]+=pma1*R_temp[Cart::yzz][Cart::z][2]+wmp1*R_temp[Cart::yzz][Cart::z][3];
R_temp[Cart::xxxy][Cart::y][2]+=pma1*R_temp[Cart::xxx][Cart::y][2]+wmp1*R_temp[Cart::xxx][Cart::y][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][3];
R_temp[Cart::xxxy][Cart::x][2]+=pma1*R_temp[Cart::xxx][Cart::x][2]+wmp1*R_temp[Cart::xxx][Cart::x][3];
R_temp[Cart::xxxy][Cart::z][2]+=pma1*R_temp[Cart::xxx][Cart::z][2]+wmp1*R_temp[Cart::xxx][Cart::z][3];
R_temp[Cart::xxyz][Cart::y][2]+=pma1*R_temp[Cart::xxz][Cart::y][2]+wmp1*R_temp[Cart::xxz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::xxz][Cart::s][3];
R_temp[Cart::xxyz][Cart::x][2]+=pma1*R_temp[Cart::xxz][Cart::x][2]+wmp1*R_temp[Cart::xxz][Cart::x][3];
R_temp[Cart::xxyz][Cart::z][2]+=pma1*R_temp[Cart::xxz][Cart::z][2]+wmp1*R_temp[Cart::xxz][Cart::z][3];
R_temp[Cart::xyzz][Cart::y][2]+=pma0*R_temp[Cart::yzz][Cart::y][2]+wmp0*R_temp[Cart::yzz][Cart::y][3];
R_temp[Cart::xyzz][Cart::x][2]+=pma0*R_temp[Cart::yzz][Cart::x][2]+wmp0*R_temp[Cart::yzz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::yzz][Cart::s][3];
R_temp[Cart::xyzz][Cart::z][2]+=pma0*R_temp[Cart::yzz][Cart::z][2]+wmp0*R_temp[Cart::yzz][Cart::z][3];
R_temp[Cart::yzzz][Cart::y][2]+=pma1*R_temp[Cart::zzz][Cart::y][2]+wmp1*R_temp[Cart::zzz][Cart::y][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][3];
R_temp[Cart::yzzz][Cart::x][2]+=pma1*R_temp[Cart::zzz][Cart::x][2]+wmp1*R_temp[Cart::zzz][Cart::x][3];
R_temp[Cart::yzzz][Cart::z][2]+=pma1*R_temp[Cart::zzz][Cart::z][2]+wmp1*R_temp[Cart::zzz][Cart::z][3];
R_temp[Cart::xxxx][Cart::y][2]+=pma0*R_temp[Cart::xxx][Cart::y][2]+wmp0*R_temp[Cart::xxx][Cart::y][3]+2*rzeta*(R_temp[Cart::xx][Cart::y][2]-gfak*R_temp[Cart::xx][Cart::y][3]);
R_temp[Cart::xxxx][Cart::x][2]+=pma0*R_temp[Cart::xxx][Cart::x][2]+wmp0*R_temp[Cart::xxx][Cart::x][3]+2*rzeta*(R_temp[Cart::xx][Cart::x][2]-gfak*R_temp[Cart::xx][Cart::x][3])+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][3];
R_temp[Cart::xxxx][Cart::z][2]+=pma0*R_temp[Cart::xxx][Cart::z][2]+wmp0*R_temp[Cart::xxx][Cart::z][3]+2*rzeta*(R_temp[Cart::xx][Cart::z][2]-gfak*R_temp[Cart::xx][Cart::z][3]);
R_temp[Cart::xxxz][Cart::y][2]+=pma2*R_temp[Cart::xxx][Cart::y][2]+wmp2*R_temp[Cart::xxx][Cart::y][3];
R_temp[Cart::xxxz][Cart::x][2]+=pma2*R_temp[Cart::xxx][Cart::x][2]+wmp2*R_temp[Cart::xxx][Cart::x][3];
R_temp[Cart::xxxz][Cart::z][2]+=pma2*R_temp[Cart::xxx][Cart::z][2]+wmp2*R_temp[Cart::xxx][Cart::z][3]+0.5/_decay*1*R_temp[Cart::xxx][Cart::s][3];
R_temp[Cart::xxzz][Cart::y][2]+=pma0*R_temp[Cart::xzz][Cart::y][2]+wmp0*R_temp[Cart::xzz][Cart::y][3];
R_temp[Cart::xxzz][Cart::x][2]+=pma0*R_temp[Cart::xzz][Cart::x][2]+wmp0*R_temp[Cart::xzz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::xzz][Cart::s][3];
R_temp[Cart::xxzz][Cart::z][2]+=pma0*R_temp[Cart::xzz][Cart::z][2]+wmp0*R_temp[Cart::xzz][Cart::z][3];
R_temp[Cart::xzzz][Cart::y][2]+=pma0*R_temp[Cart::zzz][Cart::y][2]+wmp0*R_temp[Cart::zzz][Cart::y][3];
R_temp[Cart::xzzz][Cart::x][2]+=pma0*R_temp[Cart::zzz][Cart::x][2]+wmp0*R_temp[Cart::zzz][Cart::x][3]+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][3];
R_temp[Cart::xzzz][Cart::z][2]+=pma0*R_temp[Cart::zzz][Cart::z][2]+wmp0*R_temp[Cart::zzz][Cart::z][3];
R_temp[Cart::zzzz][Cart::y][2]+=pma2*R_temp[Cart::zzz][Cart::y][2]+wmp2*R_temp[Cart::zzz][Cart::y][3]+2*rzeta*(R_temp[Cart::zz][Cart::y][2]-gfak*R_temp[Cart::zz][Cart::y][3]);
R_temp[Cart::zzzz][Cart::x][2]+=pma2*R_temp[Cart::zzz][Cart::x][2]+wmp2*R_temp[Cart::zzz][Cart::x][3]+2*rzeta*(R_temp[Cart::zz][Cart::x][2]-gfak*R_temp[Cart::zz][Cart::x][3]);
R_temp[Cart::zzzz][Cart::z][2]+=pma2*R_temp[Cart::zzz][Cart::z][2]+wmp2*R_temp[Cart::zzz][Cart::z][3]+2*rzeta*(R_temp[Cart::zz][Cart::z][2]-gfak*R_temp[Cart::zz][Cart::z][3])+0.5/_decay*1*R_temp[Cart::zzz][Cart::s][3];
}
}
//------------------------------------------------------

//Integral g - s - d - m0
if (_lmax_alpha>3 && _lmax_gamma>1){
R_temp[Cart::yyyy][Cart::yy][0]+=pma1*R_temp[Cart::yyy][Cart::yy][0]+wmp1*R_temp[Cart::yyy][Cart::yy][1]+2*rzeta*(R_temp[Cart::yy][Cart::yy][0]-gfak*R_temp[Cart::yy][Cart::yy][1])+0.5/_decay*2*R_temp[Cart::yyy][Cart::y][1];
R_temp[Cart::yyyy][Cart::xy][0]+=pma1*R_temp[Cart::yyy][Cart::xy][0]+wmp1*R_temp[Cart::yyy][Cart::xy][1]+2*rzeta*(R_temp[Cart::yy][Cart::xy][0]-gfak*R_temp[Cart::yy][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::yyy][Cart::x][1];
R_temp[Cart::yyyy][Cart::yz][0]+=pma1*R_temp[Cart::yyy][Cart::yz][0]+wmp1*R_temp[Cart::yyy][Cart::yz][1]+2*rzeta*(R_temp[Cart::yy][Cart::yz][0]-gfak*R_temp[Cart::yy][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::yyy][Cart::z][1];
R_temp[Cart::yyyy][Cart::xx][0]+=pma1*R_temp[Cart::yyy][Cart::xx][0]+wmp1*R_temp[Cart::yyy][Cart::xx][1]+2*rzeta*(R_temp[Cart::yy][Cart::xx][0]-gfak*R_temp[Cart::yy][Cart::xx][1]);
R_temp[Cart::yyyy][Cart::xz][0]+=pma1*R_temp[Cart::yyy][Cart::xz][0]+wmp1*R_temp[Cart::yyy][Cart::xz][1]+2*rzeta*(R_temp[Cart::yy][Cart::xz][0]-gfak*R_temp[Cart::yy][Cart::xz][1]);
R_temp[Cart::yyyy][Cart::zz][0]+=pma1*R_temp[Cart::yyy][Cart::zz][0]+wmp1*R_temp[Cart::yyy][Cart::zz][1]+2*rzeta*(R_temp[Cart::yy][Cart::zz][0]-gfak*R_temp[Cart::yy][Cart::zz][1]);
R_temp[Cart::xyyy][Cart::yy][0]+=pma0*R_temp[Cart::yyy][Cart::yy][0]+wmp0*R_temp[Cart::yyy][Cart::yy][1];
R_temp[Cart::xyyy][Cart::xy][0]+=pma0*R_temp[Cart::yyy][Cart::xy][0]+wmp0*R_temp[Cart::yyy][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::y][1];
R_temp[Cart::xyyy][Cart::yz][0]+=pma0*R_temp[Cart::yyy][Cart::yz][0]+wmp0*R_temp[Cart::yyy][Cart::yz][1];
R_temp[Cart::xyyy][Cart::xx][0]+=pma0*R_temp[Cart::yyy][Cart::xx][0]+wmp0*R_temp[Cart::yyy][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yyy][Cart::x][1];
R_temp[Cart::xyyy][Cart::xz][0]+=pma0*R_temp[Cart::yyy][Cart::xz][0]+wmp0*R_temp[Cart::yyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::z][1];
R_temp[Cart::xyyy][Cart::zz][0]+=pma0*R_temp[Cart::yyy][Cart::zz][0]+wmp0*R_temp[Cart::yyy][Cart::zz][1];
R_temp[Cart::yyyz][Cart::yy][0]+=pma2*R_temp[Cart::yyy][Cart::yy][0]+wmp2*R_temp[Cart::yyy][Cart::yy][1];
R_temp[Cart::yyyz][Cart::xy][0]+=pma2*R_temp[Cart::yyy][Cart::xy][0]+wmp2*R_temp[Cart::yyy][Cart::xy][1];
R_temp[Cart::yyyz][Cart::yz][0]+=pma2*R_temp[Cart::yyy][Cart::yz][0]+wmp2*R_temp[Cart::yyy][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::y][1];
R_temp[Cart::yyyz][Cart::xx][0]+=pma2*R_temp[Cart::yyy][Cart::xx][0]+wmp2*R_temp[Cart::yyy][Cart::xx][1];
R_temp[Cart::yyyz][Cart::xz][0]+=pma2*R_temp[Cart::yyy][Cart::xz][0]+wmp2*R_temp[Cart::yyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::x][1];
R_temp[Cart::yyyz][Cart::zz][0]+=pma2*R_temp[Cart::yyy][Cart::zz][0]+wmp2*R_temp[Cart::yyy][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::yyy][Cart::z][1];
R_temp[Cart::xxyy][Cart::yy][0]+=pma0*R_temp[Cart::xyy][Cart::yy][0]+wmp0*R_temp[Cart::xyy][Cart::yy][1];
R_temp[Cart::xxyy][Cart::xy][0]+=pma0*R_temp[Cart::xyy][Cart::xy][0]+wmp0*R_temp[Cart::xyy][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xyy][Cart::y][1];
R_temp[Cart::xxyy][Cart::yz][0]+=pma0*R_temp[Cart::xyy][Cart::yz][0]+wmp0*R_temp[Cart::xyy][Cart::yz][1];
R_temp[Cart::xxyy][Cart::xx][0]+=pma0*R_temp[Cart::xyy][Cart::xx][0]+wmp0*R_temp[Cart::xyy][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::xyy][Cart::x][1];
R_temp[Cart::xxyy][Cart::xz][0]+=pma0*R_temp[Cart::xyy][Cart::xz][0]+wmp0*R_temp[Cart::xyy][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xyy][Cart::z][1];
R_temp[Cart::xxyy][Cart::zz][0]+=pma0*R_temp[Cart::xyy][Cart::zz][0]+wmp0*R_temp[Cart::xyy][Cart::zz][1];
R_temp[Cart::xyyz][Cart::yy][0]+=pma0*R_temp[Cart::yyz][Cart::yy][0]+wmp0*R_temp[Cart::yyz][Cart::yy][1];
R_temp[Cart::xyyz][Cart::xy][0]+=pma0*R_temp[Cart::yyz][Cart::xy][0]+wmp0*R_temp[Cart::yyz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yyz][Cart::y][1];
R_temp[Cart::xyyz][Cart::yz][0]+=pma0*R_temp[Cart::yyz][Cart::yz][0]+wmp0*R_temp[Cart::yyz][Cart::yz][1];
R_temp[Cart::xyyz][Cart::xx][0]+=pma0*R_temp[Cart::yyz][Cart::xx][0]+wmp0*R_temp[Cart::yyz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yyz][Cart::x][1];
R_temp[Cart::xyyz][Cart::xz][0]+=pma0*R_temp[Cart::yyz][Cart::xz][0]+wmp0*R_temp[Cart::yyz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yyz][Cart::z][1];
R_temp[Cart::xyyz][Cart::zz][0]+=pma0*R_temp[Cart::yyz][Cart::zz][0]+wmp0*R_temp[Cart::yyz][Cart::zz][1];
R_temp[Cart::yyzz][Cart::yy][0]+=pma1*R_temp[Cart::yzz][Cart::yy][0]+wmp1*R_temp[Cart::yzz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::yzz][Cart::y][1];
R_temp[Cart::yyzz][Cart::xy][0]+=pma1*R_temp[Cart::yzz][Cart::xy][0]+wmp1*R_temp[Cart::yzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::x][1];
R_temp[Cart::yyzz][Cart::yz][0]+=pma1*R_temp[Cart::yzz][Cart::yz][0]+wmp1*R_temp[Cart::yzz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::z][1];
R_temp[Cart::yyzz][Cart::xx][0]+=pma1*R_temp[Cart::yzz][Cart::xx][0]+wmp1*R_temp[Cart::yzz][Cart::xx][1];
R_temp[Cart::yyzz][Cart::xz][0]+=pma1*R_temp[Cart::yzz][Cart::xz][0]+wmp1*R_temp[Cart::yzz][Cart::xz][1];
R_temp[Cart::yyzz][Cart::zz][0]+=pma1*R_temp[Cart::yzz][Cart::zz][0]+wmp1*R_temp[Cart::yzz][Cart::zz][1];
R_temp[Cart::xxxy][Cart::yy][0]+=pma1*R_temp[Cart::xxx][Cart::yy][0]+wmp1*R_temp[Cart::xxx][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxx][Cart::y][1];
R_temp[Cart::xxxy][Cart::xy][0]+=pma1*R_temp[Cart::xxx][Cart::xy][0]+wmp1*R_temp[Cart::xxx][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::x][1];
R_temp[Cart::xxxy][Cart::yz][0]+=pma1*R_temp[Cart::xxx][Cart::yz][0]+wmp1*R_temp[Cart::xxx][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::z][1];
R_temp[Cart::xxxy][Cart::xx][0]+=pma1*R_temp[Cart::xxx][Cart::xx][0]+wmp1*R_temp[Cart::xxx][Cart::xx][1];
R_temp[Cart::xxxy][Cart::xz][0]+=pma1*R_temp[Cart::xxx][Cart::xz][0]+wmp1*R_temp[Cart::xxx][Cart::xz][1];
R_temp[Cart::xxxy][Cart::zz][0]+=pma1*R_temp[Cart::xxx][Cart::zz][0]+wmp1*R_temp[Cart::xxx][Cart::zz][1];
R_temp[Cart::xxyz][Cart::yy][0]+=pma1*R_temp[Cart::xxz][Cart::yy][0]+wmp1*R_temp[Cart::xxz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::xxz][Cart::y][1];
R_temp[Cart::xxyz][Cart::xy][0]+=pma1*R_temp[Cart::xxz][Cart::xy][0]+wmp1*R_temp[Cart::xxz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xxz][Cart::x][1];
R_temp[Cart::xxyz][Cart::yz][0]+=pma1*R_temp[Cart::xxz][Cart::yz][0]+wmp1*R_temp[Cart::xxz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxz][Cart::z][1];
R_temp[Cart::xxyz][Cart::xx][0]+=pma1*R_temp[Cart::xxz][Cart::xx][0]+wmp1*R_temp[Cart::xxz][Cart::xx][1];
R_temp[Cart::xxyz][Cart::xz][0]+=pma1*R_temp[Cart::xxz][Cart::xz][0]+wmp1*R_temp[Cart::xxz][Cart::xz][1];
R_temp[Cart::xxyz][Cart::zz][0]+=pma1*R_temp[Cart::xxz][Cart::zz][0]+wmp1*R_temp[Cart::xxz][Cart::zz][1];
R_temp[Cart::xyzz][Cart::yy][0]+=pma0*R_temp[Cart::yzz][Cart::yy][0]+wmp0*R_temp[Cart::yzz][Cart::yy][1];
R_temp[Cart::xyzz][Cart::xy][0]+=pma0*R_temp[Cart::yzz][Cart::xy][0]+wmp0*R_temp[Cart::yzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::y][1];
R_temp[Cart::xyzz][Cart::yz][0]+=pma0*R_temp[Cart::yzz][Cart::yz][0]+wmp0*R_temp[Cart::yzz][Cart::yz][1];
R_temp[Cart::xyzz][Cart::xx][0]+=pma0*R_temp[Cart::yzz][Cart::xx][0]+wmp0*R_temp[Cart::yzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::yzz][Cart::x][1];
R_temp[Cart::xyzz][Cart::xz][0]+=pma0*R_temp[Cart::yzz][Cart::xz][0]+wmp0*R_temp[Cart::yzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::z][1];
R_temp[Cart::xyzz][Cart::zz][0]+=pma0*R_temp[Cart::yzz][Cart::zz][0]+wmp0*R_temp[Cart::yzz][Cart::zz][1];
R_temp[Cart::yzzz][Cart::yy][0]+=pma1*R_temp[Cart::zzz][Cart::yy][0]+wmp1*R_temp[Cart::zzz][Cart::yy][1]+0.5/_decay*2*R_temp[Cart::zzz][Cart::y][1];
R_temp[Cart::yzzz][Cart::xy][0]+=pma1*R_temp[Cart::zzz][Cart::xy][0]+wmp1*R_temp[Cart::zzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::x][1];
R_temp[Cart::yzzz][Cart::yz][0]+=pma1*R_temp[Cart::zzz][Cart::yz][0]+wmp1*R_temp[Cart::zzz][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::z][1];
R_temp[Cart::yzzz][Cart::xx][0]+=pma1*R_temp[Cart::zzz][Cart::xx][0]+wmp1*R_temp[Cart::zzz][Cart::xx][1];
R_temp[Cart::yzzz][Cart::xz][0]+=pma1*R_temp[Cart::zzz][Cart::xz][0]+wmp1*R_temp[Cart::zzz][Cart::xz][1];
R_temp[Cart::yzzz][Cart::zz][0]+=pma1*R_temp[Cart::zzz][Cart::zz][0]+wmp1*R_temp[Cart::zzz][Cart::zz][1];
R_temp[Cart::xxxx][Cart::yy][0]+=pma0*R_temp[Cart::xxx][Cart::yy][0]+wmp0*R_temp[Cart::xxx][Cart::yy][1]+2*rzeta*(R_temp[Cart::xx][Cart::yy][0]-gfak*R_temp[Cart::xx][Cart::yy][1]);
R_temp[Cart::xxxx][Cart::xy][0]+=pma0*R_temp[Cart::xxx][Cart::xy][0]+wmp0*R_temp[Cart::xxx][Cart::xy][1]+2*rzeta*(R_temp[Cart::xx][Cart::xy][0]-gfak*R_temp[Cart::xx][Cart::xy][1])+0.5/_decay*1*R_temp[Cart::xxx][Cart::y][1];
R_temp[Cart::xxxx][Cart::yz][0]+=pma0*R_temp[Cart::xxx][Cart::yz][0]+wmp0*R_temp[Cart::xxx][Cart::yz][1]+2*rzeta*(R_temp[Cart::xx][Cart::yz][0]-gfak*R_temp[Cart::xx][Cart::yz][1]);
R_temp[Cart::xxxx][Cart::xx][0]+=pma0*R_temp[Cart::xxx][Cart::xx][0]+wmp0*R_temp[Cart::xxx][Cart::xx][1]+2*rzeta*(R_temp[Cart::xx][Cart::xx][0]-gfak*R_temp[Cart::xx][Cart::xx][1])+0.5/_decay*2*R_temp[Cart::xxx][Cart::x][1];
R_temp[Cart::xxxx][Cart::xz][0]+=pma0*R_temp[Cart::xxx][Cart::xz][0]+wmp0*R_temp[Cart::xxx][Cart::xz][1]+2*rzeta*(R_temp[Cart::xx][Cart::xz][0]-gfak*R_temp[Cart::xx][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::xxx][Cart::z][1];
R_temp[Cart::xxxx][Cart::zz][0]+=pma0*R_temp[Cart::xxx][Cart::zz][0]+wmp0*R_temp[Cart::xxx][Cart::zz][1]+2*rzeta*(R_temp[Cart::xx][Cart::zz][0]-gfak*R_temp[Cart::xx][Cart::zz][1]);
R_temp[Cart::xxxz][Cart::yy][0]+=pma2*R_temp[Cart::xxx][Cart::yy][0]+wmp2*R_temp[Cart::xxx][Cart::yy][1];
R_temp[Cart::xxxz][Cart::xy][0]+=pma2*R_temp[Cart::xxx][Cart::xy][0]+wmp2*R_temp[Cart::xxx][Cart::xy][1];
R_temp[Cart::xxxz][Cart::yz][0]+=pma2*R_temp[Cart::xxx][Cart::yz][0]+wmp2*R_temp[Cart::xxx][Cart::yz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::y][1];
R_temp[Cart::xxxz][Cart::xx][0]+=pma2*R_temp[Cart::xxx][Cart::xx][0]+wmp2*R_temp[Cart::xxx][Cart::xx][1];
R_temp[Cart::xxxz][Cart::xz][0]+=pma2*R_temp[Cart::xxx][Cart::xz][0]+wmp2*R_temp[Cart::xxx][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::x][1];
R_temp[Cart::xxxz][Cart::zz][0]+=pma2*R_temp[Cart::xxx][Cart::zz][0]+wmp2*R_temp[Cart::xxx][Cart::zz][1]+0.5/_decay*2*R_temp[Cart::xxx][Cart::z][1];
R_temp[Cart::xxzz][Cart::yy][0]+=pma0*R_temp[Cart::xzz][Cart::yy][0]+wmp0*R_temp[Cart::xzz][Cart::yy][1];
R_temp[Cart::xxzz][Cart::xy][0]+=pma0*R_temp[Cart::xzz][Cart::xy][0]+wmp0*R_temp[Cart::xzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::xzz][Cart::y][1];
R_temp[Cart::xxzz][Cart::yz][0]+=pma0*R_temp[Cart::xzz][Cart::yz][0]+wmp0*R_temp[Cart::xzz][Cart::yz][1];
R_temp[Cart::xxzz][Cart::xx][0]+=pma0*R_temp[Cart::xzz][Cart::xx][0]+wmp0*R_temp[Cart::xzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::xzz][Cart::x][1];
R_temp[Cart::xxzz][Cart::xz][0]+=pma0*R_temp[Cart::xzz][Cart::xz][0]+wmp0*R_temp[Cart::xzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::xzz][Cart::z][1];
R_temp[Cart::xxzz][Cart::zz][0]+=pma0*R_temp[Cart::xzz][Cart::zz][0]+wmp0*R_temp[Cart::xzz][Cart::zz][1];
R_temp[Cart::xzzz][Cart::yy][0]+=pma0*R_temp[Cart::zzz][Cart::yy][0]+wmp0*R_temp[Cart::zzz][Cart::yy][1];
R_temp[Cart::xzzz][Cart::xy][0]+=pma0*R_temp[Cart::zzz][Cart::xy][0]+wmp0*R_temp[Cart::zzz][Cart::xy][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::y][1];
R_temp[Cart::xzzz][Cart::yz][0]+=pma0*R_temp[Cart::zzz][Cart::yz][0]+wmp0*R_temp[Cart::zzz][Cart::yz][1];
R_temp[Cart::xzzz][Cart::xx][0]+=pma0*R_temp[Cart::zzz][Cart::xx][0]+wmp0*R_temp[Cart::zzz][Cart::xx][1]+0.5/_decay*2*R_temp[Cart::zzz][Cart::x][1];
R_temp[Cart::xzzz][Cart::xz][0]+=pma0*R_temp[Cart::zzz][Cart::xz][0]+wmp0*R_temp[Cart::zzz][Cart::xz][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::z][1];
R_temp[Cart::xzzz][Cart::zz][0]+=pma0*R_temp[Cart::zzz][Cart::zz][0]+wmp0*R_temp[Cart::zzz][Cart::zz][1];
R_temp[Cart::zzzz][Cart::yy][0]+=pma2*R_temp[Cart::zzz][Cart::yy][0]+wmp2*R_temp[Cart::zzz][Cart::yy][1]+2*rzeta*(R_temp[Cart::zz][Cart::yy][0]-gfak*R_temp[Cart::zz][Cart::yy][1]);
R_temp[Cart::zzzz][Cart::xy][0]+=pma2*R_temp[Cart::zzz][Cart::xy][0]+wmp2*R_temp[Cart::zzz][Cart::xy][1]+2*rzeta*(R_temp[Cart::zz][Cart::xy][0]-gfak*R_temp[Cart::zz][Cart::xy][1]);
R_temp[Cart::zzzz][Cart::yz][0]+=pma2*R_temp[Cart::zzz][Cart::yz][0]+wmp2*R_temp[Cart::zzz][Cart::yz][1]+2*rzeta*(R_temp[Cart::zz][Cart::yz][0]-gfak*R_temp[Cart::zz][Cart::yz][1])+0.5/_decay*1*R_temp[Cart::zzz][Cart::y][1];
R_temp[Cart::zzzz][Cart::xx][0]+=pma2*R_temp[Cart::zzz][Cart::xx][0]+wmp2*R_temp[Cart::zzz][Cart::xx][1]+2*rzeta*(R_temp[Cart::zz][Cart::xx][0]-gfak*R_temp[Cart::zz][Cart::xx][1]);
R_temp[Cart::zzzz][Cart::xz][0]+=pma2*R_temp[Cart::zzz][Cart::xz][0]+wmp2*R_temp[Cart::zzz][Cart::xz][1]+2*rzeta*(R_temp[Cart::zz][Cart::xz][0]-gfak*R_temp[Cart::zz][Cart::xz][1])+0.5/_decay*1*R_temp[Cart::zzz][Cart::x][1];
R_temp[Cart::zzzz][Cart::zz][0]+=pma2*R_temp[Cart::zzz][Cart::zz][0]+wmp2*R_temp[Cart::zzz][Cart::zz][1]+2*rzeta*(R_temp[Cart::zz][Cart::zz][0]-gfak*R_temp[Cart::zz][Cart::zz][1])+0.5/_decay*2*R_temp[Cart::zzz][Cart::z][1];
}
//------------------------------------------------------

//Integral g - s - d - m1
if (_mmax >1 ){
if (_lmax_alpha>3 && _lmax_gamma>1){
R_temp[Cart::yyyy][Cart::yy][1]+=pma1*R_temp[Cart::yyy][Cart::yy][1]+wmp1*R_temp[Cart::yyy][Cart::yy][2]+2*rzeta*(R_temp[Cart::yy][Cart::yy][1]-gfak*R_temp[Cart::yy][Cart::yy][2])+0.5/_decay*2*R_temp[Cart::yyy][Cart::y][2];
R_temp[Cart::yyyy][Cart::xy][1]+=pma1*R_temp[Cart::yyy][Cart::xy][1]+wmp1*R_temp[Cart::yyy][Cart::xy][2]+2*rzeta*(R_temp[Cart::yy][Cart::xy][1]-gfak*R_temp[Cart::yy][Cart::xy][2])+0.5/_decay*1*R_temp[Cart::yyy][Cart::x][2];
R_temp[Cart::yyyy][Cart::yz][1]+=pma1*R_temp[Cart::yyy][Cart::yz][1]+wmp1*R_temp[Cart::yyy][Cart::yz][2]+2*rzeta*(R_temp[Cart::yy][Cart::yz][1]-gfak*R_temp[Cart::yy][Cart::yz][2])+0.5/_decay*1*R_temp[Cart::yyy][Cart::z][2];
R_temp[Cart::yyyy][Cart::xx][1]+=pma1*R_temp[Cart::yyy][Cart::xx][1]+wmp1*R_temp[Cart::yyy][Cart::xx][2]+2*rzeta*(R_temp[Cart::yy][Cart::xx][1]-gfak*R_temp[Cart::yy][Cart::xx][2]);
R_temp[Cart::yyyy][Cart::xz][1]+=pma1*R_temp[Cart::yyy][Cart::xz][1]+wmp1*R_temp[Cart::yyy][Cart::xz][2]+2*rzeta*(R_temp[Cart::yy][Cart::xz][1]-gfak*R_temp[Cart::yy][Cart::xz][2]);
R_temp[Cart::yyyy][Cart::zz][1]+=pma1*R_temp[Cart::yyy][Cart::zz][1]+wmp1*R_temp[Cart::yyy][Cart::zz][2]+2*rzeta*(R_temp[Cart::yy][Cart::zz][1]-gfak*R_temp[Cart::yy][Cart::zz][2]);
R_temp[Cart::xyyy][Cart::yy][1]+=pma0*R_temp[Cart::yyy][Cart::yy][1]+wmp0*R_temp[Cart::yyy][Cart::yy][2];
R_temp[Cart::xyyy][Cart::xy][1]+=pma0*R_temp[Cart::yyy][Cart::xy][1]+wmp0*R_temp[Cart::yyy][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yyy][Cart::y][2];
R_temp[Cart::xyyy][Cart::yz][1]+=pma0*R_temp[Cart::yyy][Cart::yz][1]+wmp0*R_temp[Cart::yyy][Cart::yz][2];
R_temp[Cart::xyyy][Cart::xx][1]+=pma0*R_temp[Cart::yyy][Cart::xx][1]+wmp0*R_temp[Cart::yyy][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yyy][Cart::x][2];
R_temp[Cart::xyyy][Cart::xz][1]+=pma0*R_temp[Cart::yyy][Cart::xz][1]+wmp0*R_temp[Cart::yyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyy][Cart::z][2];
R_temp[Cart::xyyy][Cart::zz][1]+=pma0*R_temp[Cart::yyy][Cart::zz][1]+wmp0*R_temp[Cart::yyy][Cart::zz][2];
R_temp[Cart::yyyz][Cart::yy][1]+=pma2*R_temp[Cart::yyy][Cart::yy][1]+wmp2*R_temp[Cart::yyy][Cart::yy][2];
R_temp[Cart::yyyz][Cart::xy][1]+=pma2*R_temp[Cart::yyy][Cart::xy][1]+wmp2*R_temp[Cart::yyy][Cart::xy][2];
R_temp[Cart::yyyz][Cart::yz][1]+=pma2*R_temp[Cart::yyy][Cart::yz][1]+wmp2*R_temp[Cart::yyy][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::yyy][Cart::y][2];
R_temp[Cart::yyyz][Cart::xx][1]+=pma2*R_temp[Cart::yyy][Cart::xx][1]+wmp2*R_temp[Cart::yyy][Cart::xx][2];
R_temp[Cart::yyyz][Cart::xz][1]+=pma2*R_temp[Cart::yyy][Cart::xz][1]+wmp2*R_temp[Cart::yyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyy][Cart::x][2];
R_temp[Cart::yyyz][Cart::zz][1]+=pma2*R_temp[Cart::yyy][Cart::zz][1]+wmp2*R_temp[Cart::yyy][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::yyy][Cart::z][2];
R_temp[Cart::xxyy][Cart::yy][1]+=pma0*R_temp[Cart::xyy][Cart::yy][1]+wmp0*R_temp[Cart::xyy][Cart::yy][2];
R_temp[Cart::xxyy][Cart::xy][1]+=pma0*R_temp[Cart::xyy][Cart::xy][1]+wmp0*R_temp[Cart::xyy][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xyy][Cart::y][2];
R_temp[Cart::xxyy][Cart::yz][1]+=pma0*R_temp[Cart::xyy][Cart::yz][1]+wmp0*R_temp[Cart::xyy][Cart::yz][2];
R_temp[Cart::xxyy][Cart::xx][1]+=pma0*R_temp[Cart::xyy][Cart::xx][1]+wmp0*R_temp[Cart::xyy][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::xyy][Cart::x][2];
R_temp[Cart::xxyy][Cart::xz][1]+=pma0*R_temp[Cart::xyy][Cart::xz][1]+wmp0*R_temp[Cart::xyy][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xyy][Cart::z][2];
R_temp[Cart::xxyy][Cart::zz][1]+=pma0*R_temp[Cart::xyy][Cart::zz][1]+wmp0*R_temp[Cart::xyy][Cart::zz][2];
R_temp[Cart::xyyz][Cart::yy][1]+=pma0*R_temp[Cart::yyz][Cart::yy][1]+wmp0*R_temp[Cart::yyz][Cart::yy][2];
R_temp[Cart::xyyz][Cart::xy][1]+=pma0*R_temp[Cart::yyz][Cart::xy][1]+wmp0*R_temp[Cart::yyz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yyz][Cart::y][2];
R_temp[Cart::xyyz][Cart::yz][1]+=pma0*R_temp[Cart::yyz][Cart::yz][1]+wmp0*R_temp[Cart::yyz][Cart::yz][2];
R_temp[Cart::xyyz][Cart::xx][1]+=pma0*R_temp[Cart::yyz][Cart::xx][1]+wmp0*R_temp[Cart::yyz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yyz][Cart::x][2];
R_temp[Cart::xyyz][Cart::xz][1]+=pma0*R_temp[Cart::yyz][Cart::xz][1]+wmp0*R_temp[Cart::yyz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yyz][Cart::z][2];
R_temp[Cart::xyyz][Cart::zz][1]+=pma0*R_temp[Cart::yyz][Cart::zz][1]+wmp0*R_temp[Cart::yyz][Cart::zz][2];
R_temp[Cart::yyzz][Cart::yy][1]+=pma1*R_temp[Cart::yzz][Cart::yy][1]+wmp1*R_temp[Cart::yzz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::yzz][Cart::y][2];
R_temp[Cart::yyzz][Cart::xy][1]+=pma1*R_temp[Cart::yzz][Cart::xy][1]+wmp1*R_temp[Cart::yzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yzz][Cart::x][2];
R_temp[Cart::yyzz][Cart::yz][1]+=pma1*R_temp[Cart::yzz][Cart::yz][1]+wmp1*R_temp[Cart::yzz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::yzz][Cart::z][2];
R_temp[Cart::yyzz][Cart::xx][1]+=pma1*R_temp[Cart::yzz][Cart::xx][1]+wmp1*R_temp[Cart::yzz][Cart::xx][2];
R_temp[Cart::yyzz][Cart::xz][1]+=pma1*R_temp[Cart::yzz][Cart::xz][1]+wmp1*R_temp[Cart::yzz][Cart::xz][2];
R_temp[Cart::yyzz][Cart::zz][1]+=pma1*R_temp[Cart::yzz][Cart::zz][1]+wmp1*R_temp[Cart::yzz][Cart::zz][2];
R_temp[Cart::xxxy][Cart::yy][1]+=pma1*R_temp[Cart::xxx][Cart::yy][1]+wmp1*R_temp[Cart::xxx][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xxx][Cart::y][2];
R_temp[Cart::xxxy][Cart::xy][1]+=pma1*R_temp[Cart::xxx][Cart::xy][1]+wmp1*R_temp[Cart::xxx][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xxx][Cart::x][2];
R_temp[Cart::xxxy][Cart::yz][1]+=pma1*R_temp[Cart::xxx][Cart::yz][1]+wmp1*R_temp[Cart::xxx][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxx][Cart::z][2];
R_temp[Cart::xxxy][Cart::xx][1]+=pma1*R_temp[Cart::xxx][Cart::xx][1]+wmp1*R_temp[Cart::xxx][Cart::xx][2];
R_temp[Cart::xxxy][Cart::xz][1]+=pma1*R_temp[Cart::xxx][Cart::xz][1]+wmp1*R_temp[Cart::xxx][Cart::xz][2];
R_temp[Cart::xxxy][Cart::zz][1]+=pma1*R_temp[Cart::xxx][Cart::zz][1]+wmp1*R_temp[Cart::xxx][Cart::zz][2];
R_temp[Cart::xxyz][Cart::yy][1]+=pma1*R_temp[Cart::xxz][Cart::yy][1]+wmp1*R_temp[Cart::xxz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::xxz][Cart::y][2];
R_temp[Cart::xxyz][Cart::xy][1]+=pma1*R_temp[Cart::xxz][Cart::xy][1]+wmp1*R_temp[Cart::xxz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xxz][Cart::x][2];
R_temp[Cart::xxyz][Cart::yz][1]+=pma1*R_temp[Cart::xxz][Cart::yz][1]+wmp1*R_temp[Cart::xxz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxz][Cart::z][2];
R_temp[Cart::xxyz][Cart::xx][1]+=pma1*R_temp[Cart::xxz][Cart::xx][1]+wmp1*R_temp[Cart::xxz][Cart::xx][2];
R_temp[Cart::xxyz][Cart::xz][1]+=pma1*R_temp[Cart::xxz][Cart::xz][1]+wmp1*R_temp[Cart::xxz][Cart::xz][2];
R_temp[Cart::xxyz][Cart::zz][1]+=pma1*R_temp[Cart::xxz][Cart::zz][1]+wmp1*R_temp[Cart::xxz][Cart::zz][2];
R_temp[Cart::xyzz][Cart::yy][1]+=pma0*R_temp[Cart::yzz][Cart::yy][1]+wmp0*R_temp[Cart::yzz][Cart::yy][2];
R_temp[Cart::xyzz][Cart::xy][1]+=pma0*R_temp[Cart::yzz][Cart::xy][1]+wmp0*R_temp[Cart::yzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::yzz][Cart::y][2];
R_temp[Cart::xyzz][Cart::yz][1]+=pma0*R_temp[Cart::yzz][Cart::yz][1]+wmp0*R_temp[Cart::yzz][Cart::yz][2];
R_temp[Cart::xyzz][Cart::xx][1]+=pma0*R_temp[Cart::yzz][Cart::xx][1]+wmp0*R_temp[Cart::yzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::yzz][Cart::x][2];
R_temp[Cart::xyzz][Cart::xz][1]+=pma0*R_temp[Cart::yzz][Cart::xz][1]+wmp0*R_temp[Cart::yzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::yzz][Cart::z][2];
R_temp[Cart::xyzz][Cart::zz][1]+=pma0*R_temp[Cart::yzz][Cart::zz][1]+wmp0*R_temp[Cart::yzz][Cart::zz][2];
R_temp[Cart::yzzz][Cart::yy][1]+=pma1*R_temp[Cart::zzz][Cart::yy][1]+wmp1*R_temp[Cart::zzz][Cart::yy][2]+0.5/_decay*2*R_temp[Cart::zzz][Cart::y][2];
R_temp[Cart::yzzz][Cart::xy][1]+=pma1*R_temp[Cart::zzz][Cart::xy][1]+wmp1*R_temp[Cart::zzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::zzz][Cart::x][2];
R_temp[Cart::yzzz][Cart::yz][1]+=pma1*R_temp[Cart::zzz][Cart::yz][1]+wmp1*R_temp[Cart::zzz][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::zzz][Cart::z][2];
R_temp[Cart::yzzz][Cart::xx][1]+=pma1*R_temp[Cart::zzz][Cart::xx][1]+wmp1*R_temp[Cart::zzz][Cart::xx][2];
R_temp[Cart::yzzz][Cart::xz][1]+=pma1*R_temp[Cart::zzz][Cart::xz][1]+wmp1*R_temp[Cart::zzz][Cart::xz][2];
R_temp[Cart::yzzz][Cart::zz][1]+=pma1*R_temp[Cart::zzz][Cart::zz][1]+wmp1*R_temp[Cart::zzz][Cart::zz][2];
R_temp[Cart::xxxx][Cart::yy][1]+=pma0*R_temp[Cart::xxx][Cart::yy][1]+wmp0*R_temp[Cart::xxx][Cart::yy][2]+2*rzeta*(R_temp[Cart::xx][Cart::yy][1]-gfak*R_temp[Cart::xx][Cart::yy][2]);
R_temp[Cart::xxxx][Cart::xy][1]+=pma0*R_temp[Cart::xxx][Cart::xy][1]+wmp0*R_temp[Cart::xxx][Cart::xy][2]+2*rzeta*(R_temp[Cart::xx][Cart::xy][1]-gfak*R_temp[Cart::xx][Cart::xy][2])+0.5/_decay*1*R_temp[Cart::xxx][Cart::y][2];
R_temp[Cart::xxxx][Cart::yz][1]+=pma0*R_temp[Cart::xxx][Cart::yz][1]+wmp0*R_temp[Cart::xxx][Cart::yz][2]+2*rzeta*(R_temp[Cart::xx][Cart::yz][1]-gfak*R_temp[Cart::xx][Cart::yz][2]);
R_temp[Cart::xxxx][Cart::xx][1]+=pma0*R_temp[Cart::xxx][Cart::xx][1]+wmp0*R_temp[Cart::xxx][Cart::xx][2]+2*rzeta*(R_temp[Cart::xx][Cart::xx][1]-gfak*R_temp[Cart::xx][Cart::xx][2])+0.5/_decay*2*R_temp[Cart::xxx][Cart::x][2];
R_temp[Cart::xxxx][Cart::xz][1]+=pma0*R_temp[Cart::xxx][Cart::xz][1]+wmp0*R_temp[Cart::xxx][Cart::xz][2]+2*rzeta*(R_temp[Cart::xx][Cart::xz][1]-gfak*R_temp[Cart::xx][Cart::xz][2])+0.5/_decay*1*R_temp[Cart::xxx][Cart::z][2];
R_temp[Cart::xxxx][Cart::zz][1]+=pma0*R_temp[Cart::xxx][Cart::zz][1]+wmp0*R_temp[Cart::xxx][Cart::zz][2]+2*rzeta*(R_temp[Cart::xx][Cart::zz][1]-gfak*R_temp[Cart::xx][Cart::zz][2]);
R_temp[Cart::xxxz][Cart::yy][1]+=pma2*R_temp[Cart::xxx][Cart::yy][1]+wmp2*R_temp[Cart::xxx][Cart::yy][2];
R_temp[Cart::xxxz][Cart::xy][1]+=pma2*R_temp[Cart::xxx][Cart::xy][1]+wmp2*R_temp[Cart::xxx][Cart::xy][2];
R_temp[Cart::xxxz][Cart::yz][1]+=pma2*R_temp[Cart::xxx][Cart::yz][1]+wmp2*R_temp[Cart::xxx][Cart::yz][2]+0.5/_decay*1*R_temp[Cart::xxx][Cart::y][2];
R_temp[Cart::xxxz][Cart::xx][1]+=pma2*R_temp[Cart::xxx][Cart::xx][1]+wmp2*R_temp[Cart::xxx][Cart::xx][2];
R_temp[Cart::xxxz][Cart::xz][1]+=pma2*R_temp[Cart::xxx][Cart::xz][1]+wmp2*R_temp[Cart::xxx][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xxx][Cart::x][2];
R_temp[Cart::xxxz][Cart::zz][1]+=pma2*R_temp[Cart::xxx][Cart::zz][1]+wmp2*R_temp[Cart::xxx][Cart::zz][2]+0.5/_decay*2*R_temp[Cart::xxx][Cart::z][2];
R_temp[Cart::xxzz][Cart::yy][1]+=pma0*R_temp[Cart::xzz][Cart::yy][1]+wmp0*R_temp[Cart::xzz][Cart::yy][2];
R_temp[Cart::xxzz][Cart::xy][1]+=pma0*R_temp[Cart::xzz][Cart::xy][1]+wmp0*R_temp[Cart::xzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::xzz][Cart::y][2];
R_temp[Cart::xxzz][Cart::yz][1]+=pma0*R_temp[Cart::xzz][Cart::yz][1]+wmp0*R_temp[Cart::xzz][Cart::yz][2];
R_temp[Cart::xxzz][Cart::xx][1]+=pma0*R_temp[Cart::xzz][Cart::xx][1]+wmp0*R_temp[Cart::xzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::xzz][Cart::x][2];
R_temp[Cart::xxzz][Cart::xz][1]+=pma0*R_temp[Cart::xzz][Cart::xz][1]+wmp0*R_temp[Cart::xzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::xzz][Cart::z][2];
R_temp[Cart::xxzz][Cart::zz][1]+=pma0*R_temp[Cart::xzz][Cart::zz][1]+wmp0*R_temp[Cart::xzz][Cart::zz][2];
R_temp[Cart::xzzz][Cart::yy][1]+=pma0*R_temp[Cart::zzz][Cart::yy][1]+wmp0*R_temp[Cart::zzz][Cart::yy][2];
R_temp[Cart::xzzz][Cart::xy][1]+=pma0*R_temp[Cart::zzz][Cart::xy][1]+wmp0*R_temp[Cart::zzz][Cart::xy][2]+0.5/_decay*1*R_temp[Cart::zzz][Cart::y][2];
R_temp[Cart::xzzz][Cart::yz][1]+=pma0*R_temp[Cart::zzz][Cart::yz][1]+wmp0*R_temp[Cart::zzz][Cart::yz][2];
R_temp[Cart::xzzz][Cart::xx][1]+=pma0*R_temp[Cart::zzz][Cart::xx][1]+wmp0*R_temp[Cart::zzz][Cart::xx][2]+0.5/_decay*2*R_temp[Cart::zzz][Cart::x][2];
R_temp[Cart::xzzz][Cart::xz][1]+=pma0*R_temp[Cart::zzz][Cart::xz][1]+wmp0*R_temp[Cart::zzz][Cart::xz][2]+0.5/_decay*1*R_temp[Cart::zzz][Cart::z][2];
R_temp[Cart::xzzz][Cart::zz][1]+=pma0*R_temp[Cart::zzz][Cart::zz][1]+wmp0*R_temp[Cart::zzz][Cart::zz][2];
R_temp[Cart::zzzz][Cart::yy][1]+=pma2*R_temp[Cart::zzz][Cart::yy][1]+wmp2*R_temp[Cart::zzz][Cart::yy][2]+2*rzeta*(R_temp[Cart::zz][Cart::yy][1]-gfak*R_temp[Cart::zz][Cart::yy][2]);
R_temp[Cart::zzzz][Cart::xy][1]+=pma2*R_temp[Cart::zzz][Cart::xy][1]+wmp2*R_temp[Cart::zzz][Cart::xy][2]+2*rzeta*(R_temp[Cart::zz][Cart::xy][1]-gfak*R_temp[Cart::zz][Cart::xy][2]);
R_temp[Cart::zzzz][Cart::yz][1]+=pma2*R_temp[Cart::zzz][Cart::yz][1]+wmp2*R_temp[Cart::zzz][Cart::yz][2]+2*rzeta*(R_temp[Cart::zz][Cart::yz][1]-gfak*R_temp[Cart::zz][Cart::yz][2])+0.5/_decay*1*R_temp[Cart::zzz][Cart::y][2];
R_temp[Cart::zzzz][Cart::xx][1]+=pma2*R_temp[Cart::zzz][Cart::xx][1]+wmp2*R_temp[Cart::zzz][Cart::xx][2]+2*rzeta*(R_temp[Cart::zz][Cart::xx][1]-gfak*R_temp[Cart::zz][Cart::xx][2]);
R_temp[Cart::zzzz][Cart::xz][1]+=pma2*R_temp[Cart::zzz][Cart::xz][1]+wmp2*R_temp[Cart::zzz][Cart::xz][2]+2*rzeta*(R_temp[Cart::zz][Cart::xz][1]-gfak*R_temp[Cart::zz][Cart::xz][2])+0.5/_decay*1*R_temp[Cart::zzz][Cart::x][2];
R_temp[Cart::zzzz][Cart::zz][1]+=pma2*R_temp[Cart::zzz][Cart::zz][1]+wmp2*R_temp[Cart::zzz][Cart::zz][2]+2*rzeta*(R_temp[Cart::zz][Cart::zz][1]-gfak*R_temp[Cart::zz][Cart::zz][2])+0.5/_decay*2*R_temp[Cart::zzz][Cart::z][2];
}
}
//------------------------------------------------------

//Integral g - s - f - m0
if (_lmax_alpha>3 && _lmax_gamma>2){
R_temp[Cart::yyyy][Cart::yyy][0]+=pma1*R_temp[Cart::yyy][Cart::yyy][0]+wmp1*R_temp[Cart::yyy][Cart::yyy][1]+2*rzeta*(R_temp[Cart::yy][Cart::yyy][0]-gfak*R_temp[Cart::yy][Cart::yyy][1])+0.5/_decay*3*R_temp[Cart::yyy][Cart::yy][1];
R_temp[Cart::yyyy][Cart::xyy][0]+=pma1*R_temp[Cart::yyy][Cart::xyy][0]+wmp1*R_temp[Cart::yyy][Cart::xyy][1]+2*rzeta*(R_temp[Cart::yy][Cart::xyy][0]-gfak*R_temp[Cart::yy][Cart::xyy][1])+0.5/_decay*2*R_temp[Cart::yyy][Cart::xy][1];
R_temp[Cart::yyyy][Cart::yyz][0]+=pma1*R_temp[Cart::yyy][Cart::yyz][0]+wmp1*R_temp[Cart::yyy][Cart::yyz][1]+2*rzeta*(R_temp[Cart::yy][Cart::yyz][0]-gfak*R_temp[Cart::yy][Cart::yyz][1])+0.5/_decay*2*R_temp[Cart::yyy][Cart::yz][1];
R_temp[Cart::yyyy][Cart::xxy][0]+=pma1*R_temp[Cart::yyy][Cart::xxy][0]+wmp1*R_temp[Cart::yyy][Cart::xxy][1]+2*rzeta*(R_temp[Cart::yy][Cart::xxy][0]-gfak*R_temp[Cart::yy][Cart::xxy][1])+0.5/_decay*1*R_temp[Cart::yyy][Cart::xx][1];
R_temp[Cart::yyyy][Cart::xyz][0]+=pma1*R_temp[Cart::yyy][Cart::xyz][0]+wmp1*R_temp[Cart::yyy][Cart::xyz][1]+2*rzeta*(R_temp[Cart::yy][Cart::xyz][0]-gfak*R_temp[Cart::yy][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::yyy][Cart::xz][1];
R_temp[Cart::yyyy][Cart::yzz][0]+=pma1*R_temp[Cart::yyy][Cart::yzz][0]+wmp1*R_temp[Cart::yyy][Cart::yzz][1]+2*rzeta*(R_temp[Cart::yy][Cart::yzz][0]-gfak*R_temp[Cart::yy][Cart::yzz][1])+0.5/_decay*1*R_temp[Cart::yyy][Cart::zz][1];
R_temp[Cart::yyyy][Cart::xxx][0]+=pma1*R_temp[Cart::yyy][Cart::xxx][0]+wmp1*R_temp[Cart::yyy][Cart::xxx][1]+2*rzeta*(R_temp[Cart::yy][Cart::xxx][0]-gfak*R_temp[Cart::yy][Cart::xxx][1]);
R_temp[Cart::yyyy][Cart::xxz][0]+=pma1*R_temp[Cart::yyy][Cart::xxz][0]+wmp1*R_temp[Cart::yyy][Cart::xxz][1]+2*rzeta*(R_temp[Cart::yy][Cart::xxz][0]-gfak*R_temp[Cart::yy][Cart::xxz][1]);
R_temp[Cart::yyyy][Cart::xzz][0]+=pma1*R_temp[Cart::yyy][Cart::xzz][0]+wmp1*R_temp[Cart::yyy][Cart::xzz][1]+2*rzeta*(R_temp[Cart::yy][Cart::xzz][0]-gfak*R_temp[Cart::yy][Cart::xzz][1]);
R_temp[Cart::yyyy][Cart::zzz][0]+=pma1*R_temp[Cart::yyy][Cart::zzz][0]+wmp1*R_temp[Cart::yyy][Cart::zzz][1]+2*rzeta*(R_temp[Cart::yy][Cart::zzz][0]-gfak*R_temp[Cart::yy][Cart::zzz][1]);
R_temp[Cart::xyyy][Cart::yyy][0]+=pma0*R_temp[Cart::yyy][Cart::yyy][0]+wmp0*R_temp[Cart::yyy][Cart::yyy][1];
R_temp[Cart::xyyy][Cart::xyy][0]+=pma0*R_temp[Cart::yyy][Cart::xyy][0]+wmp0*R_temp[Cart::yyy][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::yy][1];
R_temp[Cart::xyyy][Cart::yyz][0]+=pma0*R_temp[Cart::yyy][Cart::yyz][0]+wmp0*R_temp[Cart::yyy][Cart::yyz][1];
R_temp[Cart::xyyy][Cart::xxy][0]+=pma0*R_temp[Cart::yyy][Cart::xxy][0]+wmp0*R_temp[Cart::yyy][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yyy][Cart::xy][1];
R_temp[Cart::xyyy][Cart::xyz][0]+=pma0*R_temp[Cart::yyy][Cart::xyz][0]+wmp0*R_temp[Cart::yyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::yz][1];
R_temp[Cart::xyyy][Cart::yzz][0]+=pma0*R_temp[Cart::yyy][Cart::yzz][0]+wmp0*R_temp[Cart::yyy][Cart::yzz][1];
R_temp[Cart::xyyy][Cart::xxx][0]+=pma0*R_temp[Cart::yyy][Cart::xxx][0]+wmp0*R_temp[Cart::yyy][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yyy][Cart::xx][1];
R_temp[Cart::xyyy][Cart::xxz][0]+=pma0*R_temp[Cart::yyy][Cart::xxz][0]+wmp0*R_temp[Cart::yyy][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yyy][Cart::xz][1];
R_temp[Cart::xyyy][Cart::xzz][0]+=pma0*R_temp[Cart::yyy][Cart::xzz][0]+wmp0*R_temp[Cart::yyy][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::zz][1];
R_temp[Cart::xyyy][Cart::zzz][0]+=pma0*R_temp[Cart::yyy][Cart::zzz][0]+wmp0*R_temp[Cart::yyy][Cart::zzz][1];
R_temp[Cart::yyyz][Cart::yyy][0]+=pma2*R_temp[Cart::yyy][Cart::yyy][0]+wmp2*R_temp[Cart::yyy][Cart::yyy][1];
R_temp[Cart::yyyz][Cart::xyy][0]+=pma2*R_temp[Cart::yyy][Cart::xyy][0]+wmp2*R_temp[Cart::yyy][Cart::xyy][1];
R_temp[Cart::yyyz][Cart::yyz][0]+=pma2*R_temp[Cart::yyy][Cart::yyz][0]+wmp2*R_temp[Cart::yyy][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::yy][1];
R_temp[Cart::yyyz][Cart::xxy][0]+=pma2*R_temp[Cart::yyy][Cart::xxy][0]+wmp2*R_temp[Cart::yyy][Cart::xxy][1];
R_temp[Cart::yyyz][Cart::xyz][0]+=pma2*R_temp[Cart::yyy][Cart::xyz][0]+wmp2*R_temp[Cart::yyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::xy][1];
R_temp[Cart::yyyz][Cart::yzz][0]+=pma2*R_temp[Cart::yyy][Cart::yzz][0]+wmp2*R_temp[Cart::yyy][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::yyy][Cart::yz][1];
R_temp[Cart::yyyz][Cart::xxx][0]+=pma2*R_temp[Cart::yyy][Cart::xxx][0]+wmp2*R_temp[Cart::yyy][Cart::xxx][1];
R_temp[Cart::yyyz][Cart::xxz][0]+=pma2*R_temp[Cart::yyy][Cart::xxz][0]+wmp2*R_temp[Cart::yyy][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::yyy][Cart::xx][1];
R_temp[Cart::yyyz][Cart::xzz][0]+=pma2*R_temp[Cart::yyy][Cart::xzz][0]+wmp2*R_temp[Cart::yyy][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::yyy][Cart::xz][1];
R_temp[Cart::yyyz][Cart::zzz][0]+=pma2*R_temp[Cart::yyy][Cart::zzz][0]+wmp2*R_temp[Cart::yyy][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::yyy][Cart::zz][1];
R_temp[Cart::xxyy][Cart::yyy][0]+=pma0*R_temp[Cart::xyy][Cart::yyy][0]+wmp0*R_temp[Cart::xyy][Cart::yyy][1];
R_temp[Cart::xxyy][Cart::xyy][0]+=pma0*R_temp[Cart::xyy][Cart::xyy][0]+wmp0*R_temp[Cart::xyy][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::xyy][Cart::yy][1];
R_temp[Cart::xxyy][Cart::yyz][0]+=pma0*R_temp[Cart::xyy][Cart::yyz][0]+wmp0*R_temp[Cart::xyy][Cart::yyz][1];
R_temp[Cart::xxyy][Cart::xxy][0]+=pma0*R_temp[Cart::xyy][Cart::xxy][0]+wmp0*R_temp[Cart::xyy][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::xyy][Cart::xy][1];
R_temp[Cart::xxyy][Cart::xyz][0]+=pma0*R_temp[Cart::xyy][Cart::xyz][0]+wmp0*R_temp[Cart::xyy][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xyy][Cart::yz][1];
R_temp[Cart::xxyy][Cart::yzz][0]+=pma0*R_temp[Cart::xyy][Cart::yzz][0]+wmp0*R_temp[Cart::xyy][Cart::yzz][1];
R_temp[Cart::xxyy][Cart::xxx][0]+=pma0*R_temp[Cart::xyy][Cart::xxx][0]+wmp0*R_temp[Cart::xyy][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::xyy][Cart::xx][1];
R_temp[Cart::xxyy][Cart::xxz][0]+=pma0*R_temp[Cart::xyy][Cart::xxz][0]+wmp0*R_temp[Cart::xyy][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::xyy][Cart::xz][1];
R_temp[Cart::xxyy][Cart::xzz][0]+=pma0*R_temp[Cart::xyy][Cart::xzz][0]+wmp0*R_temp[Cart::xyy][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::xyy][Cart::zz][1];
R_temp[Cart::xxyy][Cart::zzz][0]+=pma0*R_temp[Cart::xyy][Cart::zzz][0]+wmp0*R_temp[Cart::xyy][Cart::zzz][1];
R_temp[Cart::xyyz][Cart::yyy][0]+=pma0*R_temp[Cart::yyz][Cart::yyy][0]+wmp0*R_temp[Cart::yyz][Cart::yyy][1];
R_temp[Cart::xyyz][Cart::xyy][0]+=pma0*R_temp[Cart::yyz][Cart::xyy][0]+wmp0*R_temp[Cart::yyz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yyz][Cart::yy][1];
R_temp[Cart::xyyz][Cart::yyz][0]+=pma0*R_temp[Cart::yyz][Cart::yyz][0]+wmp0*R_temp[Cart::yyz][Cart::yyz][1];
R_temp[Cart::xyyz][Cart::xxy][0]+=pma0*R_temp[Cart::yyz][Cart::xxy][0]+wmp0*R_temp[Cart::yyz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yyz][Cart::xy][1];
R_temp[Cart::xyyz][Cart::xyz][0]+=pma0*R_temp[Cart::yyz][Cart::xyz][0]+wmp0*R_temp[Cart::yyz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yyz][Cart::yz][1];
R_temp[Cart::xyyz][Cart::yzz][0]+=pma0*R_temp[Cart::yyz][Cart::yzz][0]+wmp0*R_temp[Cart::yyz][Cart::yzz][1];
R_temp[Cart::xyyz][Cart::xxx][0]+=pma0*R_temp[Cart::yyz][Cart::xxx][0]+wmp0*R_temp[Cart::yyz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yyz][Cart::xx][1];
R_temp[Cart::xyyz][Cart::xxz][0]+=pma0*R_temp[Cart::yyz][Cart::xxz][0]+wmp0*R_temp[Cart::yyz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yyz][Cart::xz][1];
R_temp[Cart::xyyz][Cart::xzz][0]+=pma0*R_temp[Cart::yyz][Cart::xzz][0]+wmp0*R_temp[Cart::yyz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yyz][Cart::zz][1];
R_temp[Cart::xyyz][Cart::zzz][0]+=pma0*R_temp[Cart::yyz][Cart::zzz][0]+wmp0*R_temp[Cart::yyz][Cart::zzz][1];
R_temp[Cart::yyzz][Cart::yyy][0]+=pma1*R_temp[Cart::yzz][Cart::yyy][0]+wmp1*R_temp[Cart::yzz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::yzz][Cart::yy][1];
R_temp[Cart::yyzz][Cart::xyy][0]+=pma1*R_temp[Cart::yzz][Cart::xyy][0]+wmp1*R_temp[Cart::yzz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::yzz][Cart::xy][1];
R_temp[Cart::yyzz][Cart::yyz][0]+=pma1*R_temp[Cart::yzz][Cart::yyz][0]+wmp1*R_temp[Cart::yzz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::yzz][Cart::yz][1];
R_temp[Cart::yyzz][Cart::xxy][0]+=pma1*R_temp[Cart::yzz][Cart::xxy][0]+wmp1*R_temp[Cart::yzz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::xx][1];
R_temp[Cart::yyzz][Cart::xyz][0]+=pma1*R_temp[Cart::yzz][Cart::xyz][0]+wmp1*R_temp[Cart::yzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::xz][1];
R_temp[Cart::yyzz][Cart::yzz][0]+=pma1*R_temp[Cart::yzz][Cart::yzz][0]+wmp1*R_temp[Cart::yzz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::zz][1];
R_temp[Cart::yyzz][Cart::xxx][0]+=pma1*R_temp[Cart::yzz][Cart::xxx][0]+wmp1*R_temp[Cart::yzz][Cart::xxx][1];
R_temp[Cart::yyzz][Cart::xxz][0]+=pma1*R_temp[Cart::yzz][Cart::xxz][0]+wmp1*R_temp[Cart::yzz][Cart::xxz][1];
R_temp[Cart::yyzz][Cart::xzz][0]+=pma1*R_temp[Cart::yzz][Cart::xzz][0]+wmp1*R_temp[Cart::yzz][Cart::xzz][1];
R_temp[Cart::yyzz][Cart::zzz][0]+=pma1*R_temp[Cart::yzz][Cart::zzz][0]+wmp1*R_temp[Cart::yzz][Cart::zzz][1];
R_temp[Cart::xxxy][Cart::yyy][0]+=pma1*R_temp[Cart::xxx][Cart::yyy][0]+wmp1*R_temp[Cart::xxx][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxx][Cart::yy][1];
R_temp[Cart::xxxy][Cart::xyy][0]+=pma1*R_temp[Cart::xxx][Cart::xyy][0]+wmp1*R_temp[Cart::xxx][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxx][Cart::xy][1];
R_temp[Cart::xxxy][Cart::yyz][0]+=pma1*R_temp[Cart::xxx][Cart::yyz][0]+wmp1*R_temp[Cart::xxx][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxx][Cart::yz][1];
R_temp[Cart::xxxy][Cart::xxy][0]+=pma1*R_temp[Cart::xxx][Cart::xxy][0]+wmp1*R_temp[Cart::xxx][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xx][1];
R_temp[Cart::xxxy][Cart::xyz][0]+=pma1*R_temp[Cart::xxx][Cart::xyz][0]+wmp1*R_temp[Cart::xxx][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xz][1];
R_temp[Cart::xxxy][Cart::yzz][0]+=pma1*R_temp[Cart::xxx][Cart::yzz][0]+wmp1*R_temp[Cart::xxx][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::zz][1];
R_temp[Cart::xxxy][Cart::xxx][0]+=pma1*R_temp[Cart::xxx][Cart::xxx][0]+wmp1*R_temp[Cart::xxx][Cart::xxx][1];
R_temp[Cart::xxxy][Cart::xxz][0]+=pma1*R_temp[Cart::xxx][Cart::xxz][0]+wmp1*R_temp[Cart::xxx][Cart::xxz][1];
R_temp[Cart::xxxy][Cart::xzz][0]+=pma1*R_temp[Cart::xxx][Cart::xzz][0]+wmp1*R_temp[Cart::xxx][Cart::xzz][1];
R_temp[Cart::xxxy][Cart::zzz][0]+=pma1*R_temp[Cart::xxx][Cart::zzz][0]+wmp1*R_temp[Cart::xxx][Cart::zzz][1];
R_temp[Cart::xxyz][Cart::yyy][0]+=pma1*R_temp[Cart::xxz][Cart::yyy][0]+wmp1*R_temp[Cart::xxz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::xxz][Cart::yy][1];
R_temp[Cart::xxyz][Cart::xyy][0]+=pma1*R_temp[Cart::xxz][Cart::xyy][0]+wmp1*R_temp[Cart::xxz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::xxz][Cart::xy][1];
R_temp[Cart::xxyz][Cart::yyz][0]+=pma1*R_temp[Cart::xxz][Cart::yyz][0]+wmp1*R_temp[Cart::xxz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::xxz][Cart::yz][1];
R_temp[Cart::xxyz][Cart::xxy][0]+=pma1*R_temp[Cart::xxz][Cart::xxy][0]+wmp1*R_temp[Cart::xxz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::xxz][Cart::xx][1];
R_temp[Cart::xxyz][Cart::xyz][0]+=pma1*R_temp[Cart::xxz][Cart::xyz][0]+wmp1*R_temp[Cart::xxz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxz][Cart::xz][1];
R_temp[Cart::xxyz][Cart::yzz][0]+=pma1*R_temp[Cart::xxz][Cart::yzz][0]+wmp1*R_temp[Cart::xxz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::xxz][Cart::zz][1];
R_temp[Cart::xxyz][Cart::xxx][0]+=pma1*R_temp[Cart::xxz][Cart::xxx][0]+wmp1*R_temp[Cart::xxz][Cart::xxx][1];
R_temp[Cart::xxyz][Cart::xxz][0]+=pma1*R_temp[Cart::xxz][Cart::xxz][0]+wmp1*R_temp[Cart::xxz][Cart::xxz][1];
R_temp[Cart::xxyz][Cart::xzz][0]+=pma1*R_temp[Cart::xxz][Cart::xzz][0]+wmp1*R_temp[Cart::xxz][Cart::xzz][1];
R_temp[Cart::xxyz][Cart::zzz][0]+=pma1*R_temp[Cart::xxz][Cart::zzz][0]+wmp1*R_temp[Cart::xxz][Cart::zzz][1];
R_temp[Cart::xyzz][Cart::yyy][0]+=pma0*R_temp[Cart::yzz][Cart::yyy][0]+wmp0*R_temp[Cart::yzz][Cart::yyy][1];
R_temp[Cart::xyzz][Cart::xyy][0]+=pma0*R_temp[Cart::yzz][Cart::xyy][0]+wmp0*R_temp[Cart::yzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::yy][1];
R_temp[Cart::xyzz][Cart::yyz][0]+=pma0*R_temp[Cart::yzz][Cart::yyz][0]+wmp0*R_temp[Cart::yzz][Cart::yyz][1];
R_temp[Cart::xyzz][Cart::xxy][0]+=pma0*R_temp[Cart::yzz][Cart::xxy][0]+wmp0*R_temp[Cart::yzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::yzz][Cart::xy][1];
R_temp[Cart::xyzz][Cart::xyz][0]+=pma0*R_temp[Cart::yzz][Cart::xyz][0]+wmp0*R_temp[Cart::yzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::yz][1];
R_temp[Cart::xyzz][Cart::yzz][0]+=pma0*R_temp[Cart::yzz][Cart::yzz][0]+wmp0*R_temp[Cart::yzz][Cart::yzz][1];
R_temp[Cart::xyzz][Cart::xxx][0]+=pma0*R_temp[Cart::yzz][Cart::xxx][0]+wmp0*R_temp[Cart::yzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::yzz][Cart::xx][1];
R_temp[Cart::xyzz][Cart::xxz][0]+=pma0*R_temp[Cart::yzz][Cart::xxz][0]+wmp0*R_temp[Cart::yzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::yzz][Cart::xz][1];
R_temp[Cart::xyzz][Cart::xzz][0]+=pma0*R_temp[Cart::yzz][Cart::xzz][0]+wmp0*R_temp[Cart::yzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::yzz][Cart::zz][1];
R_temp[Cart::xyzz][Cart::zzz][0]+=pma0*R_temp[Cart::yzz][Cart::zzz][0]+wmp0*R_temp[Cart::yzz][Cart::zzz][1];
R_temp[Cart::yzzz][Cart::yyy][0]+=pma1*R_temp[Cart::zzz][Cart::yyy][0]+wmp1*R_temp[Cart::zzz][Cart::yyy][1]+0.5/_decay*3*R_temp[Cart::zzz][Cart::yy][1];
R_temp[Cart::yzzz][Cart::xyy][0]+=pma1*R_temp[Cart::zzz][Cart::xyy][0]+wmp1*R_temp[Cart::zzz][Cart::xyy][1]+0.5/_decay*2*R_temp[Cart::zzz][Cart::xy][1];
R_temp[Cart::yzzz][Cart::yyz][0]+=pma1*R_temp[Cart::zzz][Cart::yyz][0]+wmp1*R_temp[Cart::zzz][Cart::yyz][1]+0.5/_decay*2*R_temp[Cart::zzz][Cart::yz][1];
R_temp[Cart::yzzz][Cart::xxy][0]+=pma1*R_temp[Cart::zzz][Cart::xxy][0]+wmp1*R_temp[Cart::zzz][Cart::xxy][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::xx][1];
R_temp[Cart::yzzz][Cart::xyz][0]+=pma1*R_temp[Cart::zzz][Cart::xyz][0]+wmp1*R_temp[Cart::zzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::xz][1];
R_temp[Cart::yzzz][Cart::yzz][0]+=pma1*R_temp[Cart::zzz][Cart::yzz][0]+wmp1*R_temp[Cart::zzz][Cart::yzz][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::zz][1];
R_temp[Cart::yzzz][Cart::xxx][0]+=pma1*R_temp[Cart::zzz][Cart::xxx][0]+wmp1*R_temp[Cart::zzz][Cart::xxx][1];
R_temp[Cart::yzzz][Cart::xxz][0]+=pma1*R_temp[Cart::zzz][Cart::xxz][0]+wmp1*R_temp[Cart::zzz][Cart::xxz][1];
R_temp[Cart::yzzz][Cart::xzz][0]+=pma1*R_temp[Cart::zzz][Cart::xzz][0]+wmp1*R_temp[Cart::zzz][Cart::xzz][1];
R_temp[Cart::yzzz][Cart::zzz][0]+=pma1*R_temp[Cart::zzz][Cart::zzz][0]+wmp1*R_temp[Cart::zzz][Cart::zzz][1];
R_temp[Cart::xxxx][Cart::yyy][0]+=pma0*R_temp[Cart::xxx][Cart::yyy][0]+wmp0*R_temp[Cart::xxx][Cart::yyy][1]+2*rzeta*(R_temp[Cart::xx][Cart::yyy][0]-gfak*R_temp[Cart::xx][Cart::yyy][1]);
R_temp[Cart::xxxx][Cart::xyy][0]+=pma0*R_temp[Cart::xxx][Cart::xyy][0]+wmp0*R_temp[Cart::xxx][Cart::xyy][1]+2*rzeta*(R_temp[Cart::xx][Cart::xyy][0]-gfak*R_temp[Cart::xx][Cart::xyy][1])+0.5/_decay*1*R_temp[Cart::xxx][Cart::yy][1];
R_temp[Cart::xxxx][Cart::yyz][0]+=pma0*R_temp[Cart::xxx][Cart::yyz][0]+wmp0*R_temp[Cart::xxx][Cart::yyz][1]+2*rzeta*(R_temp[Cart::xx][Cart::yyz][0]-gfak*R_temp[Cart::xx][Cart::yyz][1]);
R_temp[Cart::xxxx][Cart::xxy][0]+=pma0*R_temp[Cart::xxx][Cart::xxy][0]+wmp0*R_temp[Cart::xxx][Cart::xxy][1]+2*rzeta*(R_temp[Cart::xx][Cart::xxy][0]-gfak*R_temp[Cart::xx][Cart::xxy][1])+0.5/_decay*2*R_temp[Cart::xxx][Cart::xy][1];
R_temp[Cart::xxxx][Cart::xyz][0]+=pma0*R_temp[Cart::xxx][Cart::xyz][0]+wmp0*R_temp[Cart::xxx][Cart::xyz][1]+2*rzeta*(R_temp[Cart::xx][Cart::xyz][0]-gfak*R_temp[Cart::xx][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::xxx][Cart::yz][1];
R_temp[Cart::xxxx][Cart::yzz][0]+=pma0*R_temp[Cart::xxx][Cart::yzz][0]+wmp0*R_temp[Cart::xxx][Cart::yzz][1]+2*rzeta*(R_temp[Cart::xx][Cart::yzz][0]-gfak*R_temp[Cart::xx][Cart::yzz][1]);
R_temp[Cart::xxxx][Cart::xxx][0]+=pma0*R_temp[Cart::xxx][Cart::xxx][0]+wmp0*R_temp[Cart::xxx][Cart::xxx][1]+2*rzeta*(R_temp[Cart::xx][Cart::xxx][0]-gfak*R_temp[Cart::xx][Cart::xxx][1])+0.5/_decay*3*R_temp[Cart::xxx][Cart::xx][1];
R_temp[Cart::xxxx][Cart::xxz][0]+=pma0*R_temp[Cart::xxx][Cart::xxz][0]+wmp0*R_temp[Cart::xxx][Cart::xxz][1]+2*rzeta*(R_temp[Cart::xx][Cart::xxz][0]-gfak*R_temp[Cart::xx][Cart::xxz][1])+0.5/_decay*2*R_temp[Cart::xxx][Cart::xz][1];
R_temp[Cart::xxxx][Cart::xzz][0]+=pma0*R_temp[Cart::xxx][Cart::xzz][0]+wmp0*R_temp[Cart::xxx][Cart::xzz][1]+2*rzeta*(R_temp[Cart::xx][Cart::xzz][0]-gfak*R_temp[Cart::xx][Cart::xzz][1])+0.5/_decay*1*R_temp[Cart::xxx][Cart::zz][1];
R_temp[Cart::xxxx][Cart::zzz][0]+=pma0*R_temp[Cart::xxx][Cart::zzz][0]+wmp0*R_temp[Cart::xxx][Cart::zzz][1]+2*rzeta*(R_temp[Cart::xx][Cart::zzz][0]-gfak*R_temp[Cart::xx][Cart::zzz][1]);
R_temp[Cart::xxxz][Cart::yyy][0]+=pma2*R_temp[Cart::xxx][Cart::yyy][0]+wmp2*R_temp[Cart::xxx][Cart::yyy][1];
R_temp[Cart::xxxz][Cart::xyy][0]+=pma2*R_temp[Cart::xxx][Cart::xyy][0]+wmp2*R_temp[Cart::xxx][Cart::xyy][1];
R_temp[Cart::xxxz][Cart::yyz][0]+=pma2*R_temp[Cart::xxx][Cart::yyz][0]+wmp2*R_temp[Cart::xxx][Cart::yyz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::yy][1];
R_temp[Cart::xxxz][Cart::xxy][0]+=pma2*R_temp[Cart::xxx][Cart::xxy][0]+wmp2*R_temp[Cart::xxx][Cart::xxy][1];
R_temp[Cart::xxxz][Cart::xyz][0]+=pma2*R_temp[Cart::xxx][Cart::xyz][0]+wmp2*R_temp[Cart::xxx][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xy][1];
R_temp[Cart::xxxz][Cart::yzz][0]+=pma2*R_temp[Cart::xxx][Cart::yzz][0]+wmp2*R_temp[Cart::xxx][Cart::yzz][1]+0.5/_decay*2*R_temp[Cart::xxx][Cart::yz][1];
R_temp[Cart::xxxz][Cart::xxx][0]+=pma2*R_temp[Cart::xxx][Cart::xxx][0]+wmp2*R_temp[Cart::xxx][Cart::xxx][1];
R_temp[Cart::xxxz][Cart::xxz][0]+=pma2*R_temp[Cart::xxx][Cart::xxz][0]+wmp2*R_temp[Cart::xxx][Cart::xxz][1]+0.5/_decay*1*R_temp[Cart::xxx][Cart::xx][1];
R_temp[Cart::xxxz][Cart::xzz][0]+=pma2*R_temp[Cart::xxx][Cart::xzz][0]+wmp2*R_temp[Cart::xxx][Cart::xzz][1]+0.5/_decay*2*R_temp[Cart::xxx][Cart::xz][1];
R_temp[Cart::xxxz][Cart::zzz][0]+=pma2*R_temp[Cart::xxx][Cart::zzz][0]+wmp2*R_temp[Cart::xxx][Cart::zzz][1]+0.5/_decay*3*R_temp[Cart::xxx][Cart::zz][1];
R_temp[Cart::xxzz][Cart::yyy][0]+=pma0*R_temp[Cart::xzz][Cart::yyy][0]+wmp0*R_temp[Cart::xzz][Cart::yyy][1];
R_temp[Cart::xxzz][Cart::xyy][0]+=pma0*R_temp[Cart::xzz][Cart::xyy][0]+wmp0*R_temp[Cart::xzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::xzz][Cart::yy][1];
R_temp[Cart::xxzz][Cart::yyz][0]+=pma0*R_temp[Cart::xzz][Cart::yyz][0]+wmp0*R_temp[Cart::xzz][Cart::yyz][1];
R_temp[Cart::xxzz][Cart::xxy][0]+=pma0*R_temp[Cart::xzz][Cart::xxy][0]+wmp0*R_temp[Cart::xzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::xzz][Cart::xy][1];
R_temp[Cart::xxzz][Cart::xyz][0]+=pma0*R_temp[Cart::xzz][Cart::xyz][0]+wmp0*R_temp[Cart::xzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::xzz][Cart::yz][1];
R_temp[Cart::xxzz][Cart::yzz][0]+=pma0*R_temp[Cart::xzz][Cart::yzz][0]+wmp0*R_temp[Cart::xzz][Cart::yzz][1];
R_temp[Cart::xxzz][Cart::xxx][0]+=pma0*R_temp[Cart::xzz][Cart::xxx][0]+wmp0*R_temp[Cart::xzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::xzz][Cart::xx][1];
R_temp[Cart::xxzz][Cart::xxz][0]+=pma0*R_temp[Cart::xzz][Cart::xxz][0]+wmp0*R_temp[Cart::xzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::xzz][Cart::xz][1];
R_temp[Cart::xxzz][Cart::xzz][0]+=pma0*R_temp[Cart::xzz][Cart::xzz][0]+wmp0*R_temp[Cart::xzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::xzz][Cart::zz][1];
R_temp[Cart::xxzz][Cart::zzz][0]+=pma0*R_temp[Cart::xzz][Cart::zzz][0]+wmp0*R_temp[Cart::xzz][Cart::zzz][1];
R_temp[Cart::xzzz][Cart::yyy][0]+=pma0*R_temp[Cart::zzz][Cart::yyy][0]+wmp0*R_temp[Cart::zzz][Cart::yyy][1];
R_temp[Cart::xzzz][Cart::xyy][0]+=pma0*R_temp[Cart::zzz][Cart::xyy][0]+wmp0*R_temp[Cart::zzz][Cart::xyy][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::yy][1];
R_temp[Cart::xzzz][Cart::yyz][0]+=pma0*R_temp[Cart::zzz][Cart::yyz][0]+wmp0*R_temp[Cart::zzz][Cart::yyz][1];
R_temp[Cart::xzzz][Cart::xxy][0]+=pma0*R_temp[Cart::zzz][Cart::xxy][0]+wmp0*R_temp[Cart::zzz][Cart::xxy][1]+0.5/_decay*2*R_temp[Cart::zzz][Cart::xy][1];
R_temp[Cart::xzzz][Cart::xyz][0]+=pma0*R_temp[Cart::zzz][Cart::xyz][0]+wmp0*R_temp[Cart::zzz][Cart::xyz][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::yz][1];
R_temp[Cart::xzzz][Cart::yzz][0]+=pma0*R_temp[Cart::zzz][Cart::yzz][0]+wmp0*R_temp[Cart::zzz][Cart::yzz][1];
R_temp[Cart::xzzz][Cart::xxx][0]+=pma0*R_temp[Cart::zzz][Cart::xxx][0]+wmp0*R_temp[Cart::zzz][Cart::xxx][1]+0.5/_decay*3*R_temp[Cart::zzz][Cart::xx][1];
R_temp[Cart::xzzz][Cart::xxz][0]+=pma0*R_temp[Cart::zzz][Cart::xxz][0]+wmp0*R_temp[Cart::zzz][Cart::xxz][1]+0.5/_decay*2*R_temp[Cart::zzz][Cart::xz][1];
R_temp[Cart::xzzz][Cart::xzz][0]+=pma0*R_temp[Cart::zzz][Cart::xzz][0]+wmp0*R_temp[Cart::zzz][Cart::xzz][1]+0.5/_decay*1*R_temp[Cart::zzz][Cart::zz][1];
R_temp[Cart::xzzz][Cart::zzz][0]+=pma0*R_temp[Cart::zzz][Cart::zzz][0]+wmp0*R_temp[Cart::zzz][Cart::zzz][1];
R_temp[Cart::zzzz][Cart::yyy][0]+=pma2*R_temp[Cart::zzz][Cart::yyy][0]+wmp2*R_temp[Cart::zzz][Cart::yyy][1]+2*rzeta*(R_temp[Cart::zz][Cart::yyy][0]-gfak*R_temp[Cart::zz][Cart::yyy][1]);
R_temp[Cart::zzzz][Cart::xyy][0]+=pma2*R_temp[Cart::zzz][Cart::xyy][0]+wmp2*R_temp[Cart::zzz][Cart::xyy][1]+2*rzeta*(R_temp[Cart::zz][Cart::xyy][0]-gfak*R_temp[Cart::zz][Cart::xyy][1]);
R_temp[Cart::zzzz][Cart::yyz][0]+=pma2*R_temp[Cart::zzz][Cart::yyz][0]+wmp2*R_temp[Cart::zzz][Cart::yyz][1]+2*rzeta*(R_temp[Cart::zz][Cart::yyz][0]-gfak*R_temp[Cart::zz][Cart::yyz][1])+0.5/_decay*1*R_temp[Cart::zzz][Cart::yy][1];
R_temp[Cart::zzzz][Cart::xxy][0]+=pma2*R_temp[Cart::zzz][Cart::xxy][0]+wmp2*R_temp[Cart::zzz][Cart::xxy][1]+2*rzeta*(R_temp[Cart::zz][Cart::xxy][0]-gfak*R_temp[Cart::zz][Cart::xxy][1]);
R_temp[Cart::zzzz][Cart::xyz][0]+=pma2*R_temp[Cart::zzz][Cart::xyz][0]+wmp2*R_temp[Cart::zzz][Cart::xyz][1]+2*rzeta*(R_temp[Cart::zz][Cart::xyz][0]-gfak*R_temp[Cart::zz][Cart::xyz][1])+0.5/_decay*1*R_temp[Cart::zzz][Cart::xy][1];
R_temp[Cart::zzzz][Cart::yzz][0]+=pma2*R_temp[Cart::zzz][Cart::yzz][0]+wmp2*R_temp[Cart::zzz][Cart::yzz][1]+2*rzeta*(R_temp[Cart::zz][Cart::yzz][0]-gfak*R_temp[Cart::zz][Cart::yzz][1])+0.5/_decay*2*R_temp[Cart::zzz][Cart::yz][1];
R_temp[Cart::zzzz][Cart::xxx][0]+=pma2*R_temp[Cart::zzz][Cart::xxx][0]+wmp2*R_temp[Cart::zzz][Cart::xxx][1]+2*rzeta*(R_temp[Cart::zz][Cart::xxx][0]-gfak*R_temp[Cart::zz][Cart::xxx][1]);
R_temp[Cart::zzzz][Cart::xxz][0]+=pma2*R_temp[Cart::zzz][Cart::xxz][0]+wmp2*R_temp[Cart::zzz][Cart::xxz][1]+2*rzeta*(R_temp[Cart::zz][Cart::xxz][0]-gfak*R_temp[Cart::zz][Cart::xxz][1])+0.5/_decay*1*R_temp[Cart::zzz][Cart::xx][1];
R_temp[Cart::zzzz][Cart::xzz][0]+=pma2*R_temp[Cart::zzz][Cart::xzz][0]+wmp2*R_temp[Cart::zzz][Cart::xzz][1]+2*rzeta*(R_temp[Cart::zz][Cart::xzz][0]-gfak*R_temp[Cart::zz][Cart::xzz][1])+0.5/_decay*2*R_temp[Cart::zzz][Cart::xz][1];
R_temp[Cart::zzzz][Cart::zzz][0]+=pma2*R_temp[Cart::zzz][Cart::zzz][0]+wmp2*R_temp[Cart::zzz][Cart::zzz][1]+2*rzeta*(R_temp[Cart::zz][Cart::zzz][0]-gfak*R_temp[Cart::zz][Cart::zzz][1])+0.5/_decay*3*R_temp[Cart::zzz][Cart::zz][1];
}
//------------------------------------------------------


            cout << "halloende" << endl;


            
//copy into new array for 3D use.       
            
for (index i = 0; i != _ncombined; ++i) {

         for (index k = 0; k != _ngamma; ++k) {

                            R[i][0][k] = R_temp[i][k][0];
                        }

                }      




//Integral s - p - s
if (_lmax_beta>0){

R[Cart::s][Cart::y][Cart::s]=R[Cart::y][Cart::s][Cart::s]+amb1*R[Cart::s][Cart::s][Cart::s];
R[Cart::s][Cart::x][Cart::s]=R[Cart::x][Cart::s][Cart::s]+amb0*R[Cart::s][Cart::s][Cart::s];
R[Cart::s][Cart::z][Cart::s]=R[Cart::z][Cart::s][Cart::s]+amb2*R[Cart::s][Cart::s][Cart::s];
}
//------------------------------------------------------

//Integral p - p - s
if (_lmax_beta>0 && _lmax_alpha>0){

R[Cart::y][Cart::y][Cart::s]=R[Cart::yy][Cart::s][Cart::s]+amb1*R[Cart::y][Cart::s][Cart::s];
R[Cart::x][Cart::y][Cart::s]=R[Cart::xy][Cart::s][Cart::s]+amb1*R[Cart::x][Cart::s][Cart::s];
R[Cart::z][Cart::y][Cart::s]=R[Cart::yz][Cart::s][Cart::s]+amb1*R[Cart::z][Cart::s][Cart::s];
R[Cart::y][Cart::x][Cart::s]=R[Cart::xy][Cart::s][Cart::s]+amb0*R[Cart::y][Cart::s][Cart::s];
R[Cart::x][Cart::x][Cart::s]=R[Cart::xx][Cart::s][Cart::s]+amb0*R[Cart::x][Cart::s][Cart::s];
R[Cart::z][Cart::x][Cart::s]=R[Cart::xz][Cart::s][Cart::s]+amb0*R[Cart::z][Cart::s][Cart::s];
R[Cart::y][Cart::z][Cart::s]=R[Cart::yz][Cart::s][Cart::s]+amb2*R[Cart::y][Cart::s][Cart::s];
R[Cart::x][Cart::z][Cart::s]=R[Cart::xz][Cart::s][Cart::s]+amb2*R[Cart::x][Cart::s][Cart::s];
R[Cart::z][Cart::z][Cart::s]=R[Cart::zz][Cart::s][Cart::s]+amb2*R[Cart::z][Cart::s][Cart::s];
}
//------------------------------------------------------

//Integral d - p - s
if (_lmax_beta>0 && _lmax_alpha>1){

R[Cart::yy][Cart::y][Cart::s]=R[Cart::yyy][Cart::s][Cart::s]+amb1*R[Cart::yy][Cart::s][Cart::s];
R[Cart::xy][Cart::y][Cart::s]=R[Cart::xyy][Cart::s][Cart::s]+amb1*R[Cart::xy][Cart::s][Cart::s];
R[Cart::yz][Cart::y][Cart::s]=R[Cart::yyz][Cart::s][Cart::s]+amb1*R[Cart::yz][Cart::s][Cart::s];
R[Cart::xx][Cart::y][Cart::s]=R[Cart::xxy][Cart::s][Cart::s]+amb1*R[Cart::xx][Cart::s][Cart::s];
R[Cart::xz][Cart::y][Cart::s]=R[Cart::xyz][Cart::s][Cart::s]+amb1*R[Cart::xz][Cart::s][Cart::s];
R[Cart::zz][Cart::y][Cart::s]=R[Cart::yzz][Cart::s][Cart::s]+amb1*R[Cart::zz][Cart::s][Cart::s];
R[Cart::yy][Cart::x][Cart::s]=R[Cart::xyy][Cart::s][Cart::s]+amb0*R[Cart::yy][Cart::s][Cart::s];
R[Cart::xy][Cart::x][Cart::s]=R[Cart::xxy][Cart::s][Cart::s]+amb0*R[Cart::xy][Cart::s][Cart::s];
R[Cart::yz][Cart::x][Cart::s]=R[Cart::xyz][Cart::s][Cart::s]+amb0*R[Cart::yz][Cart::s][Cart::s];
R[Cart::xx][Cart::x][Cart::s]=R[Cart::xxx][Cart::s][Cart::s]+amb0*R[Cart::xx][Cart::s][Cart::s];
R[Cart::xz][Cart::x][Cart::s]=R[Cart::xxz][Cart::s][Cart::s]+amb0*R[Cart::xz][Cart::s][Cart::s];
R[Cart::zz][Cart::x][Cart::s]=R[Cart::xzz][Cart::s][Cart::s]+amb0*R[Cart::zz][Cart::s][Cart::s];
R[Cart::yy][Cart::z][Cart::s]=R[Cart::yyz][Cart::s][Cart::s]+amb2*R[Cart::yy][Cart::s][Cart::s];
R[Cart::xy][Cart::z][Cart::s]=R[Cart::xyz][Cart::s][Cart::s]+amb2*R[Cart::xy][Cart::s][Cart::s];
R[Cart::yz][Cart::z][Cart::s]=R[Cart::yzz][Cart::s][Cart::s]+amb2*R[Cart::yz][Cart::s][Cart::s];
R[Cart::xx][Cart::z][Cart::s]=R[Cart::xxz][Cart::s][Cart::s]+amb2*R[Cart::xx][Cart::s][Cart::s];
R[Cart::xz][Cart::z][Cart::s]=R[Cart::xzz][Cart::s][Cart::s]+amb2*R[Cart::xz][Cart::s][Cart::s];
R[Cart::zz][Cart::z][Cart::s]=R[Cart::zzz][Cart::s][Cart::s]+amb2*R[Cart::zz][Cart::s][Cart::s];
}
//------------------------------------------------------

//Integral f - p - s
if (_lmax_beta>0 && _lmax_alpha>2){

R[Cart::yyy][Cart::y][Cart::s]=R[Cart::yyyy][Cart::s][Cart::s]+amb1*R[Cart::yyy][Cart::s][Cart::s];
R[Cart::xyy][Cart::y][Cart::s]=R[Cart::xyyy][Cart::s][Cart::s]+amb1*R[Cart::xyy][Cart::s][Cart::s];
R[Cart::yyz][Cart::y][Cart::s]=R[Cart::yyyz][Cart::s][Cart::s]+amb1*R[Cart::yyz][Cart::s][Cart::s];
R[Cart::xxy][Cart::y][Cart::s]=R[Cart::xxyy][Cart::s][Cart::s]+amb1*R[Cart::xxy][Cart::s][Cart::s];
R[Cart::xyz][Cart::y][Cart::s]=R[Cart::xyyz][Cart::s][Cart::s]+amb1*R[Cart::xyz][Cart::s][Cart::s];
R[Cart::yzz][Cart::y][Cart::s]=R[Cart::yyzz][Cart::s][Cart::s]+amb1*R[Cart::yzz][Cart::s][Cart::s];
R[Cart::xxx][Cart::y][Cart::s]=R[Cart::xxxy][Cart::s][Cart::s]+amb1*R[Cart::xxx][Cart::s][Cart::s];
R[Cart::xxz][Cart::y][Cart::s]=R[Cart::xxyz][Cart::s][Cart::s]+amb1*R[Cart::xxz][Cart::s][Cart::s];
R[Cart::xzz][Cart::y][Cart::s]=R[Cart::xyzz][Cart::s][Cart::s]+amb1*R[Cart::xzz][Cart::s][Cart::s];
R[Cart::zzz][Cart::y][Cart::s]=R[Cart::yzzz][Cart::s][Cart::s]+amb1*R[Cart::zzz][Cart::s][Cart::s];
R[Cart::yyy][Cart::x][Cart::s]=R[Cart::xyyy][Cart::s][Cart::s]+amb0*R[Cart::yyy][Cart::s][Cart::s];
R[Cart::xyy][Cart::x][Cart::s]=R[Cart::xxyy][Cart::s][Cart::s]+amb0*R[Cart::xyy][Cart::s][Cart::s];
R[Cart::yyz][Cart::x][Cart::s]=R[Cart::xyyz][Cart::s][Cart::s]+amb0*R[Cart::yyz][Cart::s][Cart::s];
R[Cart::xxy][Cart::x][Cart::s]=R[Cart::xxxy][Cart::s][Cart::s]+amb0*R[Cart::xxy][Cart::s][Cart::s];
R[Cart::xyz][Cart::x][Cart::s]=R[Cart::xxyz][Cart::s][Cart::s]+amb0*R[Cart::xyz][Cart::s][Cart::s];
R[Cart::yzz][Cart::x][Cart::s]=R[Cart::xyzz][Cart::s][Cart::s]+amb0*R[Cart::yzz][Cart::s][Cart::s];
R[Cart::xxx][Cart::x][Cart::s]=R[Cart::xxxx][Cart::s][Cart::s]+amb0*R[Cart::xxx][Cart::s][Cart::s];
R[Cart::xxz][Cart::x][Cart::s]=R[Cart::xxxz][Cart::s][Cart::s]+amb0*R[Cart::xxz][Cart::s][Cart::s];
R[Cart::xzz][Cart::x][Cart::s]=R[Cart::xxzz][Cart::s][Cart::s]+amb0*R[Cart::xzz][Cart::s][Cart::s];
R[Cart::zzz][Cart::x][Cart::s]=R[Cart::xzzz][Cart::s][Cart::s]+amb0*R[Cart::zzz][Cart::s][Cart::s];
R[Cart::yyy][Cart::z][Cart::s]=R[Cart::yyyz][Cart::s][Cart::s]+amb2*R[Cart::yyy][Cart::s][Cart::s];
R[Cart::xyy][Cart::z][Cart::s]=R[Cart::xyyz][Cart::s][Cart::s]+amb2*R[Cart::xyy][Cart::s][Cart::s];
R[Cart::yyz][Cart::z][Cart::s]=R[Cart::yyzz][Cart::s][Cart::s]+amb2*R[Cart::yyz][Cart::s][Cart::s];
R[Cart::xxy][Cart::z][Cart::s]=R[Cart::xxyz][Cart::s][Cart::s]+amb2*R[Cart::xxy][Cart::s][Cart::s];
R[Cart::xyz][Cart::z][Cart::s]=R[Cart::xyzz][Cart::s][Cart::s]+amb2*R[Cart::xyz][Cart::s][Cart::s];
R[Cart::yzz][Cart::z][Cart::s]=R[Cart::yzzz][Cart::s][Cart::s]+amb2*R[Cart::yzz][Cart::s][Cart::s];
R[Cart::xxx][Cart::z][Cart::s]=R[Cart::xxxz][Cart::s][Cart::s]+amb2*R[Cart::xxx][Cart::s][Cart::s];
R[Cart::xxz][Cart::z][Cart::s]=R[Cart::xxzz][Cart::s][Cart::s]+amb2*R[Cart::xxz][Cart::s][Cart::s];
R[Cart::xzz][Cart::z][Cart::s]=R[Cart::xzzz][Cart::s][Cart::s]+amb2*R[Cart::xzz][Cart::s][Cart::s];
R[Cart::zzz][Cart::z][Cart::s]=R[Cart::zzzz][Cart::s][Cart::s]+amb2*R[Cart::zzz][Cart::s][Cart::s];
}
//------------------------------------------------------

//Integral s - p - p
if (_lmax_beta>0 && _lmax_gamma>0){

R[Cart::s][Cart::y][Cart::y]=R[Cart::y][Cart::s][Cart::y]+amb1*R[Cart::s][Cart::s][Cart::y];
R[Cart::s][Cart::x][Cart::y]=R[Cart::x][Cart::s][Cart::y]+amb0*R[Cart::s][Cart::s][Cart::y];
R[Cart::s][Cart::z][Cart::y]=R[Cart::z][Cart::s][Cart::y]+amb2*R[Cart::s][Cart::s][Cart::y];
R[Cart::s][Cart::y][Cart::x]=R[Cart::y][Cart::s][Cart::x]+amb1*R[Cart::s][Cart::s][Cart::x];
R[Cart::s][Cart::x][Cart::x]=R[Cart::x][Cart::s][Cart::x]+amb0*R[Cart::s][Cart::s][Cart::x];
R[Cart::s][Cart::z][Cart::x]=R[Cart::z][Cart::s][Cart::x]+amb2*R[Cart::s][Cart::s][Cart::x];
R[Cart::s][Cart::y][Cart::z]=R[Cart::y][Cart::s][Cart::z]+amb1*R[Cart::s][Cart::s][Cart::z];
R[Cart::s][Cart::x][Cart::z]=R[Cart::x][Cart::s][Cart::z]+amb0*R[Cart::s][Cart::s][Cart::z];
R[Cart::s][Cart::z][Cart::z]=R[Cart::z][Cart::s][Cart::z]+amb2*R[Cart::s][Cart::s][Cart::z];
}
//------------------------------------------------------

//Integral p - p - p
if (_lmax_beta>0 && _lmax_alpha>0 && _lmax_gamma>0){

R[Cart::y][Cart::y][Cart::y]=R[Cart::yy][Cart::s][Cart::y]+amb1*R[Cart::y][Cart::s][Cart::y];
R[Cart::x][Cart::y][Cart::y]=R[Cart::xy][Cart::s][Cart::y]+amb1*R[Cart::x][Cart::s][Cart::y];
R[Cart::z][Cart::y][Cart::y]=R[Cart::yz][Cart::s][Cart::y]+amb1*R[Cart::z][Cart::s][Cart::y];
R[Cart::y][Cart::x][Cart::y]=R[Cart::xy][Cart::s][Cart::y]+amb0*R[Cart::y][Cart::s][Cart::y];
R[Cart::x][Cart::x][Cart::y]=R[Cart::xx][Cart::s][Cart::y]+amb0*R[Cart::x][Cart::s][Cart::y];
R[Cart::z][Cart::x][Cart::y]=R[Cart::xz][Cart::s][Cart::y]+amb0*R[Cart::z][Cart::s][Cart::y];
R[Cart::y][Cart::z][Cart::y]=R[Cart::yz][Cart::s][Cart::y]+amb2*R[Cart::y][Cart::s][Cart::y];
R[Cart::x][Cart::z][Cart::y]=R[Cart::xz][Cart::s][Cart::y]+amb2*R[Cart::x][Cart::s][Cart::y];
R[Cart::z][Cart::z][Cart::y]=R[Cart::zz][Cart::s][Cart::y]+amb2*R[Cart::z][Cart::s][Cart::y];
R[Cart::y][Cart::y][Cart::x]=R[Cart::yy][Cart::s][Cart::x]+amb1*R[Cart::y][Cart::s][Cart::x];
R[Cart::x][Cart::y][Cart::x]=R[Cart::xy][Cart::s][Cart::x]+amb1*R[Cart::x][Cart::s][Cart::x];
R[Cart::z][Cart::y][Cart::x]=R[Cart::yz][Cart::s][Cart::x]+amb1*R[Cart::z][Cart::s][Cart::x];
R[Cart::y][Cart::x][Cart::x]=R[Cart::xy][Cart::s][Cart::x]+amb0*R[Cart::y][Cart::s][Cart::x];
R[Cart::x][Cart::x][Cart::x]=R[Cart::xx][Cart::s][Cart::x]+amb0*R[Cart::x][Cart::s][Cart::x];
R[Cart::z][Cart::x][Cart::x]=R[Cart::xz][Cart::s][Cart::x]+amb0*R[Cart::z][Cart::s][Cart::x];
R[Cart::y][Cart::z][Cart::x]=R[Cart::yz][Cart::s][Cart::x]+amb2*R[Cart::y][Cart::s][Cart::x];
R[Cart::x][Cart::z][Cart::x]=R[Cart::xz][Cart::s][Cart::x]+amb2*R[Cart::x][Cart::s][Cart::x];
R[Cart::z][Cart::z][Cart::x]=R[Cart::zz][Cart::s][Cart::x]+amb2*R[Cart::z][Cart::s][Cart::x];
R[Cart::y][Cart::y][Cart::z]=R[Cart::yy][Cart::s][Cart::z]+amb1*R[Cart::y][Cart::s][Cart::z];
R[Cart::x][Cart::y][Cart::z]=R[Cart::xy][Cart::s][Cart::z]+amb1*R[Cart::x][Cart::s][Cart::z];
R[Cart::z][Cart::y][Cart::z]=R[Cart::yz][Cart::s][Cart::z]+amb1*R[Cart::z][Cart::s][Cart::z];
R[Cart::y][Cart::x][Cart::z]=R[Cart::xy][Cart::s][Cart::z]+amb0*R[Cart::y][Cart::s][Cart::z];
R[Cart::x][Cart::x][Cart::z]=R[Cart::xx][Cart::s][Cart::z]+amb0*R[Cart::x][Cart::s][Cart::z];
R[Cart::z][Cart::x][Cart::z]=R[Cart::xz][Cart::s][Cart::z]+amb0*R[Cart::z][Cart::s][Cart::z];
R[Cart::y][Cart::z][Cart::z]=R[Cart::yz][Cart::s][Cart::z]+amb2*R[Cart::y][Cart::s][Cart::z];
R[Cart::x][Cart::z][Cart::z]=R[Cart::xz][Cart::s][Cart::z]+amb2*R[Cart::x][Cart::s][Cart::z];
R[Cart::z][Cart::z][Cart::z]=R[Cart::zz][Cart::s][Cart::z]+amb2*R[Cart::z][Cart::s][Cart::z];
}
//------------------------------------------------------

//Integral d - p - p
if (_lmax_beta>0 && _lmax_alpha>1 && _lmax_gamma>0){

R[Cart::yy][Cart::y][Cart::y]=R[Cart::yyy][Cart::s][Cart::y]+amb1*R[Cart::yy][Cart::s][Cart::y];
R[Cart::xy][Cart::y][Cart::y]=R[Cart::xyy][Cart::s][Cart::y]+amb1*R[Cart::xy][Cart::s][Cart::y];
R[Cart::yz][Cart::y][Cart::y]=R[Cart::yyz][Cart::s][Cart::y]+amb1*R[Cart::yz][Cart::s][Cart::y];
R[Cart::xx][Cart::y][Cart::y]=R[Cart::xxy][Cart::s][Cart::y]+amb1*R[Cart::xx][Cart::s][Cart::y];
R[Cart::xz][Cart::y][Cart::y]=R[Cart::xyz][Cart::s][Cart::y]+amb1*R[Cart::xz][Cart::s][Cart::y];
R[Cart::zz][Cart::y][Cart::y]=R[Cart::yzz][Cart::s][Cart::y]+amb1*R[Cart::zz][Cart::s][Cart::y];
R[Cart::yy][Cart::x][Cart::y]=R[Cart::xyy][Cart::s][Cart::y]+amb0*R[Cart::yy][Cart::s][Cart::y];
R[Cart::xy][Cart::x][Cart::y]=R[Cart::xxy][Cart::s][Cart::y]+amb0*R[Cart::xy][Cart::s][Cart::y];
R[Cart::yz][Cart::x][Cart::y]=R[Cart::xyz][Cart::s][Cart::y]+amb0*R[Cart::yz][Cart::s][Cart::y];
R[Cart::xx][Cart::x][Cart::y]=R[Cart::xxx][Cart::s][Cart::y]+amb0*R[Cart::xx][Cart::s][Cart::y];
R[Cart::xz][Cart::x][Cart::y]=R[Cart::xxz][Cart::s][Cart::y]+amb0*R[Cart::xz][Cart::s][Cart::y];
R[Cart::zz][Cart::x][Cart::y]=R[Cart::xzz][Cart::s][Cart::y]+amb0*R[Cart::zz][Cart::s][Cart::y];
R[Cart::yy][Cart::z][Cart::y]=R[Cart::yyz][Cart::s][Cart::y]+amb2*R[Cart::yy][Cart::s][Cart::y];
R[Cart::xy][Cart::z][Cart::y]=R[Cart::xyz][Cart::s][Cart::y]+amb2*R[Cart::xy][Cart::s][Cart::y];
R[Cart::yz][Cart::z][Cart::y]=R[Cart::yzz][Cart::s][Cart::y]+amb2*R[Cart::yz][Cart::s][Cart::y];
R[Cart::xx][Cart::z][Cart::y]=R[Cart::xxz][Cart::s][Cart::y]+amb2*R[Cart::xx][Cart::s][Cart::y];
R[Cart::xz][Cart::z][Cart::y]=R[Cart::xzz][Cart::s][Cart::y]+amb2*R[Cart::xz][Cart::s][Cart::y];
R[Cart::zz][Cart::z][Cart::y]=R[Cart::zzz][Cart::s][Cart::y]+amb2*R[Cart::zz][Cart::s][Cart::y];
R[Cart::yy][Cart::y][Cart::x]=R[Cart::yyy][Cart::s][Cart::x]+amb1*R[Cart::yy][Cart::s][Cart::x];
R[Cart::xy][Cart::y][Cart::x]=R[Cart::xyy][Cart::s][Cart::x]+amb1*R[Cart::xy][Cart::s][Cart::x];
R[Cart::yz][Cart::y][Cart::x]=R[Cart::yyz][Cart::s][Cart::x]+amb1*R[Cart::yz][Cart::s][Cart::x];
R[Cart::xx][Cart::y][Cart::x]=R[Cart::xxy][Cart::s][Cart::x]+amb1*R[Cart::xx][Cart::s][Cart::x];
R[Cart::xz][Cart::y][Cart::x]=R[Cart::xyz][Cart::s][Cart::x]+amb1*R[Cart::xz][Cart::s][Cart::x];
R[Cart::zz][Cart::y][Cart::x]=R[Cart::yzz][Cart::s][Cart::x]+amb1*R[Cart::zz][Cart::s][Cart::x];
R[Cart::yy][Cart::x][Cart::x]=R[Cart::xyy][Cart::s][Cart::x]+amb0*R[Cart::yy][Cart::s][Cart::x];
R[Cart::xy][Cart::x][Cart::x]=R[Cart::xxy][Cart::s][Cart::x]+amb0*R[Cart::xy][Cart::s][Cart::x];
R[Cart::yz][Cart::x][Cart::x]=R[Cart::xyz][Cart::s][Cart::x]+amb0*R[Cart::yz][Cart::s][Cart::x];
R[Cart::xx][Cart::x][Cart::x]=R[Cart::xxx][Cart::s][Cart::x]+amb0*R[Cart::xx][Cart::s][Cart::x];
R[Cart::xz][Cart::x][Cart::x]=R[Cart::xxz][Cart::s][Cart::x]+amb0*R[Cart::xz][Cart::s][Cart::x];
R[Cart::zz][Cart::x][Cart::x]=R[Cart::xzz][Cart::s][Cart::x]+amb0*R[Cart::zz][Cart::s][Cart::x];
R[Cart::yy][Cart::z][Cart::x]=R[Cart::yyz][Cart::s][Cart::x]+amb2*R[Cart::yy][Cart::s][Cart::x];
R[Cart::xy][Cart::z][Cart::x]=R[Cart::xyz][Cart::s][Cart::x]+amb2*R[Cart::xy][Cart::s][Cart::x];
R[Cart::yz][Cart::z][Cart::x]=R[Cart::yzz][Cart::s][Cart::x]+amb2*R[Cart::yz][Cart::s][Cart::x];
R[Cart::xx][Cart::z][Cart::x]=R[Cart::xxz][Cart::s][Cart::x]+amb2*R[Cart::xx][Cart::s][Cart::x];
R[Cart::xz][Cart::z][Cart::x]=R[Cart::xzz][Cart::s][Cart::x]+amb2*R[Cart::xz][Cart::s][Cart::x];
R[Cart::zz][Cart::z][Cart::x]=R[Cart::zzz][Cart::s][Cart::x]+amb2*R[Cart::zz][Cart::s][Cart::x];
R[Cart::yy][Cart::y][Cart::z]=R[Cart::yyy][Cart::s][Cart::z]+amb1*R[Cart::yy][Cart::s][Cart::z];
R[Cart::xy][Cart::y][Cart::z]=R[Cart::xyy][Cart::s][Cart::z]+amb1*R[Cart::xy][Cart::s][Cart::z];
R[Cart::yz][Cart::y][Cart::z]=R[Cart::yyz][Cart::s][Cart::z]+amb1*R[Cart::yz][Cart::s][Cart::z];
R[Cart::xx][Cart::y][Cart::z]=R[Cart::xxy][Cart::s][Cart::z]+amb1*R[Cart::xx][Cart::s][Cart::z];
R[Cart::xz][Cart::y][Cart::z]=R[Cart::xyz][Cart::s][Cart::z]+amb1*R[Cart::xz][Cart::s][Cart::z];
R[Cart::zz][Cart::y][Cart::z]=R[Cart::yzz][Cart::s][Cart::z]+amb1*R[Cart::zz][Cart::s][Cart::z];
R[Cart::yy][Cart::x][Cart::z]=R[Cart::xyy][Cart::s][Cart::z]+amb0*R[Cart::yy][Cart::s][Cart::z];
R[Cart::xy][Cart::x][Cart::z]=R[Cart::xxy][Cart::s][Cart::z]+amb0*R[Cart::xy][Cart::s][Cart::z];
R[Cart::yz][Cart::x][Cart::z]=R[Cart::xyz][Cart::s][Cart::z]+amb0*R[Cart::yz][Cart::s][Cart::z];
R[Cart::xx][Cart::x][Cart::z]=R[Cart::xxx][Cart::s][Cart::z]+amb0*R[Cart::xx][Cart::s][Cart::z];
R[Cart::xz][Cart::x][Cart::z]=R[Cart::xxz][Cart::s][Cart::z]+amb0*R[Cart::xz][Cart::s][Cart::z];
R[Cart::zz][Cart::x][Cart::z]=R[Cart::xzz][Cart::s][Cart::z]+amb0*R[Cart::zz][Cart::s][Cart::z];
R[Cart::yy][Cart::z][Cart::z]=R[Cart::yyz][Cart::s][Cart::z]+amb2*R[Cart::yy][Cart::s][Cart::z];
R[Cart::xy][Cart::z][Cart::z]=R[Cart::xyz][Cart::s][Cart::z]+amb2*R[Cart::xy][Cart::s][Cart::z];
R[Cart::yz][Cart::z][Cart::z]=R[Cart::yzz][Cart::s][Cart::z]+amb2*R[Cart::yz][Cart::s][Cart::z];
R[Cart::xx][Cart::z][Cart::z]=R[Cart::xxz][Cart::s][Cart::z]+amb2*R[Cart::xx][Cart::s][Cart::z];
R[Cart::xz][Cart::z][Cart::z]=R[Cart::xzz][Cart::s][Cart::z]+amb2*R[Cart::xz][Cart::s][Cart::z];
R[Cart::zz][Cart::z][Cart::z]=R[Cart::zzz][Cart::s][Cart::z]+amb2*R[Cart::zz][Cart::s][Cart::z];
}
//------------------------------------------------------

//Integral f - p - p
if (_lmax_beta>0 && _lmax_alpha>2 && _lmax_gamma>0){

R[Cart::yyy][Cart::y][Cart::y]=R[Cart::yyyy][Cart::s][Cart::y]+amb1*R[Cart::yyy][Cart::s][Cart::y];
R[Cart::xyy][Cart::y][Cart::y]=R[Cart::xyyy][Cart::s][Cart::y]+amb1*R[Cart::xyy][Cart::s][Cart::y];
R[Cart::yyz][Cart::y][Cart::y]=R[Cart::yyyz][Cart::s][Cart::y]+amb1*R[Cart::yyz][Cart::s][Cart::y];
R[Cart::xxy][Cart::y][Cart::y]=R[Cart::xxyy][Cart::s][Cart::y]+amb1*R[Cart::xxy][Cart::s][Cart::y];
R[Cart::xyz][Cart::y][Cart::y]=R[Cart::xyyz][Cart::s][Cart::y]+amb1*R[Cart::xyz][Cart::s][Cart::y];
R[Cart::yzz][Cart::y][Cart::y]=R[Cart::yyzz][Cart::s][Cart::y]+amb1*R[Cart::yzz][Cart::s][Cart::y];
R[Cart::xxx][Cart::y][Cart::y]=R[Cart::xxxy][Cart::s][Cart::y]+amb1*R[Cart::xxx][Cart::s][Cart::y];
R[Cart::xxz][Cart::y][Cart::y]=R[Cart::xxyz][Cart::s][Cart::y]+amb1*R[Cart::xxz][Cart::s][Cart::y];
R[Cart::xzz][Cart::y][Cart::y]=R[Cart::xyzz][Cart::s][Cart::y]+amb1*R[Cart::xzz][Cart::s][Cart::y];
R[Cart::zzz][Cart::y][Cart::y]=R[Cart::yzzz][Cart::s][Cart::y]+amb1*R[Cart::zzz][Cart::s][Cart::y];
R[Cart::yyy][Cart::x][Cart::y]=R[Cart::xyyy][Cart::s][Cart::y]+amb0*R[Cart::yyy][Cart::s][Cart::y];
R[Cart::xyy][Cart::x][Cart::y]=R[Cart::xxyy][Cart::s][Cart::y]+amb0*R[Cart::xyy][Cart::s][Cart::y];
R[Cart::yyz][Cart::x][Cart::y]=R[Cart::xyyz][Cart::s][Cart::y]+amb0*R[Cart::yyz][Cart::s][Cart::y];
R[Cart::xxy][Cart::x][Cart::y]=R[Cart::xxxy][Cart::s][Cart::y]+amb0*R[Cart::xxy][Cart::s][Cart::y];
R[Cart::xyz][Cart::x][Cart::y]=R[Cart::xxyz][Cart::s][Cart::y]+amb0*R[Cart::xyz][Cart::s][Cart::y];
R[Cart::yzz][Cart::x][Cart::y]=R[Cart::xyzz][Cart::s][Cart::y]+amb0*R[Cart::yzz][Cart::s][Cart::y];
R[Cart::xxx][Cart::x][Cart::y]=R[Cart::xxxx][Cart::s][Cart::y]+amb0*R[Cart::xxx][Cart::s][Cart::y];
R[Cart::xxz][Cart::x][Cart::y]=R[Cart::xxxz][Cart::s][Cart::y]+amb0*R[Cart::xxz][Cart::s][Cart::y];
R[Cart::xzz][Cart::x][Cart::y]=R[Cart::xxzz][Cart::s][Cart::y]+amb0*R[Cart::xzz][Cart::s][Cart::y];
R[Cart::zzz][Cart::x][Cart::y]=R[Cart::xzzz][Cart::s][Cart::y]+amb0*R[Cart::zzz][Cart::s][Cart::y];
R[Cart::yyy][Cart::z][Cart::y]=R[Cart::yyyz][Cart::s][Cart::y]+amb2*R[Cart::yyy][Cart::s][Cart::y];
R[Cart::xyy][Cart::z][Cart::y]=R[Cart::xyyz][Cart::s][Cart::y]+amb2*R[Cart::xyy][Cart::s][Cart::y];
R[Cart::yyz][Cart::z][Cart::y]=R[Cart::yyzz][Cart::s][Cart::y]+amb2*R[Cart::yyz][Cart::s][Cart::y];
R[Cart::xxy][Cart::z][Cart::y]=R[Cart::xxyz][Cart::s][Cart::y]+amb2*R[Cart::xxy][Cart::s][Cart::y];
R[Cart::xyz][Cart::z][Cart::y]=R[Cart::xyzz][Cart::s][Cart::y]+amb2*R[Cart::xyz][Cart::s][Cart::y];
R[Cart::yzz][Cart::z][Cart::y]=R[Cart::yzzz][Cart::s][Cart::y]+amb2*R[Cart::yzz][Cart::s][Cart::y];
R[Cart::xxx][Cart::z][Cart::y]=R[Cart::xxxz][Cart::s][Cart::y]+amb2*R[Cart::xxx][Cart::s][Cart::y];
R[Cart::xxz][Cart::z][Cart::y]=R[Cart::xxzz][Cart::s][Cart::y]+amb2*R[Cart::xxz][Cart::s][Cart::y];
R[Cart::xzz][Cart::z][Cart::y]=R[Cart::xzzz][Cart::s][Cart::y]+amb2*R[Cart::xzz][Cart::s][Cart::y];
R[Cart::zzz][Cart::z][Cart::y]=R[Cart::zzzz][Cart::s][Cart::y]+amb2*R[Cart::zzz][Cart::s][Cart::y];
R[Cart::yyy][Cart::y][Cart::x]=R[Cart::yyyy][Cart::s][Cart::x]+amb1*R[Cart::yyy][Cart::s][Cart::x];
R[Cart::xyy][Cart::y][Cart::x]=R[Cart::xyyy][Cart::s][Cart::x]+amb1*R[Cart::xyy][Cart::s][Cart::x];
R[Cart::yyz][Cart::y][Cart::x]=R[Cart::yyyz][Cart::s][Cart::x]+amb1*R[Cart::yyz][Cart::s][Cart::x];
R[Cart::xxy][Cart::y][Cart::x]=R[Cart::xxyy][Cart::s][Cart::x]+amb1*R[Cart::xxy][Cart::s][Cart::x];
R[Cart::xyz][Cart::y][Cart::x]=R[Cart::xyyz][Cart::s][Cart::x]+amb1*R[Cart::xyz][Cart::s][Cart::x];
R[Cart::yzz][Cart::y][Cart::x]=R[Cart::yyzz][Cart::s][Cart::x]+amb1*R[Cart::yzz][Cart::s][Cart::x];
R[Cart::xxx][Cart::y][Cart::x]=R[Cart::xxxy][Cart::s][Cart::x]+amb1*R[Cart::xxx][Cart::s][Cart::x];
R[Cart::xxz][Cart::y][Cart::x]=R[Cart::xxyz][Cart::s][Cart::x]+amb1*R[Cart::xxz][Cart::s][Cart::x];
R[Cart::xzz][Cart::y][Cart::x]=R[Cart::xyzz][Cart::s][Cart::x]+amb1*R[Cart::xzz][Cart::s][Cart::x];
R[Cart::zzz][Cart::y][Cart::x]=R[Cart::yzzz][Cart::s][Cart::x]+amb1*R[Cart::zzz][Cart::s][Cart::x];
R[Cart::yyy][Cart::x][Cart::x]=R[Cart::xyyy][Cart::s][Cart::x]+amb0*R[Cart::yyy][Cart::s][Cart::x];
R[Cart::xyy][Cart::x][Cart::x]=R[Cart::xxyy][Cart::s][Cart::x]+amb0*R[Cart::xyy][Cart::s][Cart::x];
R[Cart::yyz][Cart::x][Cart::x]=R[Cart::xyyz][Cart::s][Cart::x]+amb0*R[Cart::yyz][Cart::s][Cart::x];
R[Cart::xxy][Cart::x][Cart::x]=R[Cart::xxxy][Cart::s][Cart::x]+amb0*R[Cart::xxy][Cart::s][Cart::x];
R[Cart::xyz][Cart::x][Cart::x]=R[Cart::xxyz][Cart::s][Cart::x]+amb0*R[Cart::xyz][Cart::s][Cart::x];
R[Cart::yzz][Cart::x][Cart::x]=R[Cart::xyzz][Cart::s][Cart::x]+amb0*R[Cart::yzz][Cart::s][Cart::x];
R[Cart::xxx][Cart::x][Cart::x]=R[Cart::xxxx][Cart::s][Cart::x]+amb0*R[Cart::xxx][Cart::s][Cart::x];
R[Cart::xxz][Cart::x][Cart::x]=R[Cart::xxxz][Cart::s][Cart::x]+amb0*R[Cart::xxz][Cart::s][Cart::x];
R[Cart::xzz][Cart::x][Cart::x]=R[Cart::xxzz][Cart::s][Cart::x]+amb0*R[Cart::xzz][Cart::s][Cart::x];
R[Cart::zzz][Cart::x][Cart::x]=R[Cart::xzzz][Cart::s][Cart::x]+amb0*R[Cart::zzz][Cart::s][Cart::x];
R[Cart::yyy][Cart::z][Cart::x]=R[Cart::yyyz][Cart::s][Cart::x]+amb2*R[Cart::yyy][Cart::s][Cart::x];
R[Cart::xyy][Cart::z][Cart::x]=R[Cart::xyyz][Cart::s][Cart::x]+amb2*R[Cart::xyy][Cart::s][Cart::x];
R[Cart::yyz][Cart::z][Cart::x]=R[Cart::yyzz][Cart::s][Cart::x]+amb2*R[Cart::yyz][Cart::s][Cart::x];
R[Cart::xxy][Cart::z][Cart::x]=R[Cart::xxyz][Cart::s][Cart::x]+amb2*R[Cart::xxy][Cart::s][Cart::x];
R[Cart::xyz][Cart::z][Cart::x]=R[Cart::xyzz][Cart::s][Cart::x]+amb2*R[Cart::xyz][Cart::s][Cart::x];
R[Cart::yzz][Cart::z][Cart::x]=R[Cart::yzzz][Cart::s][Cart::x]+amb2*R[Cart::yzz][Cart::s][Cart::x];
R[Cart::xxx][Cart::z][Cart::x]=R[Cart::xxxz][Cart::s][Cart::x]+amb2*R[Cart::xxx][Cart::s][Cart::x];
R[Cart::xxz][Cart::z][Cart::x]=R[Cart::xxzz][Cart::s][Cart::x]+amb2*R[Cart::xxz][Cart::s][Cart::x];
R[Cart::xzz][Cart::z][Cart::x]=R[Cart::xzzz][Cart::s][Cart::x]+amb2*R[Cart::xzz][Cart::s][Cart::x];
R[Cart::zzz][Cart::z][Cart::x]=R[Cart::zzzz][Cart::s][Cart::x]+amb2*R[Cart::zzz][Cart::s][Cart::x];
R[Cart::yyy][Cart::y][Cart::z]=R[Cart::yyyy][Cart::s][Cart::z]+amb1*R[Cart::yyy][Cart::s][Cart::z];
R[Cart::xyy][Cart::y][Cart::z]=R[Cart::xyyy][Cart::s][Cart::z]+amb1*R[Cart::xyy][Cart::s][Cart::z];
R[Cart::yyz][Cart::y][Cart::z]=R[Cart::yyyz][Cart::s][Cart::z]+amb1*R[Cart::yyz][Cart::s][Cart::z];
R[Cart::xxy][Cart::y][Cart::z]=R[Cart::xxyy][Cart::s][Cart::z]+amb1*R[Cart::xxy][Cart::s][Cart::z];
R[Cart::xyz][Cart::y][Cart::z]=R[Cart::xyyz][Cart::s][Cart::z]+amb1*R[Cart::xyz][Cart::s][Cart::z];
R[Cart::yzz][Cart::y][Cart::z]=R[Cart::yyzz][Cart::s][Cart::z]+amb1*R[Cart::yzz][Cart::s][Cart::z];
R[Cart::xxx][Cart::y][Cart::z]=R[Cart::xxxy][Cart::s][Cart::z]+amb1*R[Cart::xxx][Cart::s][Cart::z];
R[Cart::xxz][Cart::y][Cart::z]=R[Cart::xxyz][Cart::s][Cart::z]+amb1*R[Cart::xxz][Cart::s][Cart::z];
R[Cart::xzz][Cart::y][Cart::z]=R[Cart::xyzz][Cart::s][Cart::z]+amb1*R[Cart::xzz][Cart::s][Cart::z];
R[Cart::zzz][Cart::y][Cart::z]=R[Cart::yzzz][Cart::s][Cart::z]+amb1*R[Cart::zzz][Cart::s][Cart::z];
R[Cart::yyy][Cart::x][Cart::z]=R[Cart::xyyy][Cart::s][Cart::z]+amb0*R[Cart::yyy][Cart::s][Cart::z];
R[Cart::xyy][Cart::x][Cart::z]=R[Cart::xxyy][Cart::s][Cart::z]+amb0*R[Cart::xyy][Cart::s][Cart::z];
R[Cart::yyz][Cart::x][Cart::z]=R[Cart::xyyz][Cart::s][Cart::z]+amb0*R[Cart::yyz][Cart::s][Cart::z];
R[Cart::xxy][Cart::x][Cart::z]=R[Cart::xxxy][Cart::s][Cart::z]+amb0*R[Cart::xxy][Cart::s][Cart::z];
R[Cart::xyz][Cart::x][Cart::z]=R[Cart::xxyz][Cart::s][Cart::z]+amb0*R[Cart::xyz][Cart::s][Cart::z];
R[Cart::yzz][Cart::x][Cart::z]=R[Cart::xyzz][Cart::s][Cart::z]+amb0*R[Cart::yzz][Cart::s][Cart::z];
R[Cart::xxx][Cart::x][Cart::z]=R[Cart::xxxx][Cart::s][Cart::z]+amb0*R[Cart::xxx][Cart::s][Cart::z];
R[Cart::xxz][Cart::x][Cart::z]=R[Cart::xxxz][Cart::s][Cart::z]+amb0*R[Cart::xxz][Cart::s][Cart::z];
R[Cart::xzz][Cart::x][Cart::z]=R[Cart::xxzz][Cart::s][Cart::z]+amb0*R[Cart::xzz][Cart::s][Cart::z];
R[Cart::zzz][Cart::x][Cart::z]=R[Cart::xzzz][Cart::s][Cart::z]+amb0*R[Cart::zzz][Cart::s][Cart::z];
R[Cart::yyy][Cart::z][Cart::z]=R[Cart::yyyz][Cart::s][Cart::z]+amb2*R[Cart::yyy][Cart::s][Cart::z];
R[Cart::xyy][Cart::z][Cart::z]=R[Cart::xyyz][Cart::s][Cart::z]+amb2*R[Cart::xyy][Cart::s][Cart::z];
R[Cart::yyz][Cart::z][Cart::z]=R[Cart::yyzz][Cart::s][Cart::z]+amb2*R[Cart::yyz][Cart::s][Cart::z];
R[Cart::xxy][Cart::z][Cart::z]=R[Cart::xxyz][Cart::s][Cart::z]+amb2*R[Cart::xxy][Cart::s][Cart::z];
R[Cart::xyz][Cart::z][Cart::z]=R[Cart::xyzz][Cart::s][Cart::z]+amb2*R[Cart::xyz][Cart::s][Cart::z];
R[Cart::yzz][Cart::z][Cart::z]=R[Cart::yzzz][Cart::s][Cart::z]+amb2*R[Cart::yzz][Cart::s][Cart::z];
R[Cart::xxx][Cart::z][Cart::z]=R[Cart::xxxz][Cart::s][Cart::z]+amb2*R[Cart::xxx][Cart::s][Cart::z];
R[Cart::xxz][Cart::z][Cart::z]=R[Cart::xxzz][Cart::s][Cart::z]+amb2*R[Cart::xxz][Cart::s][Cart::z];
R[Cart::xzz][Cart::z][Cart::z]=R[Cart::xzzz][Cart::s][Cart::z]+amb2*R[Cart::xzz][Cart::s][Cart::z];
R[Cart::zzz][Cart::z][Cart::z]=R[Cart::zzzz][Cart::s][Cart::z]+amb2*R[Cart::zzz][Cart::s][Cart::z];
}
//------------------------------------------------------

//Integral s - p - d
if (_lmax_beta>0 && _lmax_gamma>1){

R[Cart::s][Cart::y][Cart::yy]=R[Cart::y][Cart::s][Cart::yy]+amb1*R[Cart::s][Cart::s][Cart::yy];
R[Cart::s][Cart::x][Cart::yy]=R[Cart::x][Cart::s][Cart::yy]+amb0*R[Cart::s][Cart::s][Cart::yy];
R[Cart::s][Cart::z][Cart::yy]=R[Cart::z][Cart::s][Cart::yy]+amb2*R[Cart::s][Cart::s][Cart::yy];
R[Cart::s][Cart::y][Cart::xy]=R[Cart::y][Cart::s][Cart::xy]+amb1*R[Cart::s][Cart::s][Cart::xy];
R[Cart::s][Cart::x][Cart::xy]=R[Cart::x][Cart::s][Cart::xy]+amb0*R[Cart::s][Cart::s][Cart::xy];
R[Cart::s][Cart::z][Cart::xy]=R[Cart::z][Cart::s][Cart::xy]+amb2*R[Cart::s][Cart::s][Cart::xy];
R[Cart::s][Cart::y][Cart::yz]=R[Cart::y][Cart::s][Cart::yz]+amb1*R[Cart::s][Cart::s][Cart::yz];
R[Cart::s][Cart::x][Cart::yz]=R[Cart::x][Cart::s][Cart::yz]+amb0*R[Cart::s][Cart::s][Cart::yz];
R[Cart::s][Cart::z][Cart::yz]=R[Cart::z][Cart::s][Cart::yz]+amb2*R[Cart::s][Cart::s][Cart::yz];
R[Cart::s][Cart::y][Cart::xx]=R[Cart::y][Cart::s][Cart::xx]+amb1*R[Cart::s][Cart::s][Cart::xx];
R[Cart::s][Cart::x][Cart::xx]=R[Cart::x][Cart::s][Cart::xx]+amb0*R[Cart::s][Cart::s][Cart::xx];
R[Cart::s][Cart::z][Cart::xx]=R[Cart::z][Cart::s][Cart::xx]+amb2*R[Cart::s][Cart::s][Cart::xx];
R[Cart::s][Cart::y][Cart::xz]=R[Cart::y][Cart::s][Cart::xz]+amb1*R[Cart::s][Cart::s][Cart::xz];
R[Cart::s][Cart::x][Cart::xz]=R[Cart::x][Cart::s][Cart::xz]+amb0*R[Cart::s][Cart::s][Cart::xz];
R[Cart::s][Cart::z][Cart::xz]=R[Cart::z][Cart::s][Cart::xz]+amb2*R[Cart::s][Cart::s][Cart::xz];
R[Cart::s][Cart::y][Cart::zz]=R[Cart::y][Cart::s][Cart::zz]+amb1*R[Cart::s][Cart::s][Cart::zz];
R[Cart::s][Cart::x][Cart::zz]=R[Cart::x][Cart::s][Cart::zz]+amb0*R[Cart::s][Cart::s][Cart::zz];
R[Cart::s][Cart::z][Cart::zz]=R[Cart::z][Cart::s][Cart::zz]+amb2*R[Cart::s][Cart::s][Cart::zz];
}
//------------------------------------------------------

//Integral p - p - d
if (_lmax_beta>0 && _lmax_alpha>0 && _lmax_gamma>1){

R[Cart::y][Cart::y][Cart::yy]=R[Cart::yy][Cart::s][Cart::yy]+amb1*R[Cart::y][Cart::s][Cart::yy];
R[Cart::x][Cart::y][Cart::yy]=R[Cart::xy][Cart::s][Cart::yy]+amb1*R[Cart::x][Cart::s][Cart::yy];
R[Cart::z][Cart::y][Cart::yy]=R[Cart::yz][Cart::s][Cart::yy]+amb1*R[Cart::z][Cart::s][Cart::yy];
R[Cart::y][Cart::x][Cart::yy]=R[Cart::xy][Cart::s][Cart::yy]+amb0*R[Cart::y][Cart::s][Cart::yy];
R[Cart::x][Cart::x][Cart::yy]=R[Cart::xx][Cart::s][Cart::yy]+amb0*R[Cart::x][Cart::s][Cart::yy];
R[Cart::z][Cart::x][Cart::yy]=R[Cart::xz][Cart::s][Cart::yy]+amb0*R[Cart::z][Cart::s][Cart::yy];
R[Cart::y][Cart::z][Cart::yy]=R[Cart::yz][Cart::s][Cart::yy]+amb2*R[Cart::y][Cart::s][Cart::yy];
R[Cart::x][Cart::z][Cart::yy]=R[Cart::xz][Cart::s][Cart::yy]+amb2*R[Cart::x][Cart::s][Cart::yy];
R[Cart::z][Cart::z][Cart::yy]=R[Cart::zz][Cart::s][Cart::yy]+amb2*R[Cart::z][Cart::s][Cart::yy];
R[Cart::y][Cart::y][Cart::xy]=R[Cart::yy][Cart::s][Cart::xy]+amb1*R[Cart::y][Cart::s][Cart::xy];
R[Cart::x][Cart::y][Cart::xy]=R[Cart::xy][Cart::s][Cart::xy]+amb1*R[Cart::x][Cart::s][Cart::xy];
R[Cart::z][Cart::y][Cart::xy]=R[Cart::yz][Cart::s][Cart::xy]+amb1*R[Cart::z][Cart::s][Cart::xy];
R[Cart::y][Cart::x][Cart::xy]=R[Cart::xy][Cart::s][Cart::xy]+amb0*R[Cart::y][Cart::s][Cart::xy];
R[Cart::x][Cart::x][Cart::xy]=R[Cart::xx][Cart::s][Cart::xy]+amb0*R[Cart::x][Cart::s][Cart::xy];
R[Cart::z][Cart::x][Cart::xy]=R[Cart::xz][Cart::s][Cart::xy]+amb0*R[Cart::z][Cart::s][Cart::xy];
R[Cart::y][Cart::z][Cart::xy]=R[Cart::yz][Cart::s][Cart::xy]+amb2*R[Cart::y][Cart::s][Cart::xy];
R[Cart::x][Cart::z][Cart::xy]=R[Cart::xz][Cart::s][Cart::xy]+amb2*R[Cart::x][Cart::s][Cart::xy];
R[Cart::z][Cart::z][Cart::xy]=R[Cart::zz][Cart::s][Cart::xy]+amb2*R[Cart::z][Cart::s][Cart::xy];
R[Cart::y][Cart::y][Cart::yz]=R[Cart::yy][Cart::s][Cart::yz]+amb1*R[Cart::y][Cart::s][Cart::yz];
R[Cart::x][Cart::y][Cart::yz]=R[Cart::xy][Cart::s][Cart::yz]+amb1*R[Cart::x][Cart::s][Cart::yz];
R[Cart::z][Cart::y][Cart::yz]=R[Cart::yz][Cart::s][Cart::yz]+amb1*R[Cart::z][Cart::s][Cart::yz];
R[Cart::y][Cart::x][Cart::yz]=R[Cart::xy][Cart::s][Cart::yz]+amb0*R[Cart::y][Cart::s][Cart::yz];
R[Cart::x][Cart::x][Cart::yz]=R[Cart::xx][Cart::s][Cart::yz]+amb0*R[Cart::x][Cart::s][Cart::yz];
R[Cart::z][Cart::x][Cart::yz]=R[Cart::xz][Cart::s][Cart::yz]+amb0*R[Cart::z][Cart::s][Cart::yz];
R[Cart::y][Cart::z][Cart::yz]=R[Cart::yz][Cart::s][Cart::yz]+amb2*R[Cart::y][Cart::s][Cart::yz];
R[Cart::x][Cart::z][Cart::yz]=R[Cart::xz][Cart::s][Cart::yz]+amb2*R[Cart::x][Cart::s][Cart::yz];
R[Cart::z][Cart::z][Cart::yz]=R[Cart::zz][Cart::s][Cart::yz]+amb2*R[Cart::z][Cart::s][Cart::yz];
R[Cart::y][Cart::y][Cart::xx]=R[Cart::yy][Cart::s][Cart::xx]+amb1*R[Cart::y][Cart::s][Cart::xx];
R[Cart::x][Cart::y][Cart::xx]=R[Cart::xy][Cart::s][Cart::xx]+amb1*R[Cart::x][Cart::s][Cart::xx];
R[Cart::z][Cart::y][Cart::xx]=R[Cart::yz][Cart::s][Cart::xx]+amb1*R[Cart::z][Cart::s][Cart::xx];
R[Cart::y][Cart::x][Cart::xx]=R[Cart::xy][Cart::s][Cart::xx]+amb0*R[Cart::y][Cart::s][Cart::xx];
R[Cart::x][Cart::x][Cart::xx]=R[Cart::xx][Cart::s][Cart::xx]+amb0*R[Cart::x][Cart::s][Cart::xx];
R[Cart::z][Cart::x][Cart::xx]=R[Cart::xz][Cart::s][Cart::xx]+amb0*R[Cart::z][Cart::s][Cart::xx];
R[Cart::y][Cart::z][Cart::xx]=R[Cart::yz][Cart::s][Cart::xx]+amb2*R[Cart::y][Cart::s][Cart::xx];
R[Cart::x][Cart::z][Cart::xx]=R[Cart::xz][Cart::s][Cart::xx]+amb2*R[Cart::x][Cart::s][Cart::xx];
R[Cart::z][Cart::z][Cart::xx]=R[Cart::zz][Cart::s][Cart::xx]+amb2*R[Cart::z][Cart::s][Cart::xx];
R[Cart::y][Cart::y][Cart::xz]=R[Cart::yy][Cart::s][Cart::xz]+amb1*R[Cart::y][Cart::s][Cart::xz];
R[Cart::x][Cart::y][Cart::xz]=R[Cart::xy][Cart::s][Cart::xz]+amb1*R[Cart::x][Cart::s][Cart::xz];
R[Cart::z][Cart::y][Cart::xz]=R[Cart::yz][Cart::s][Cart::xz]+amb1*R[Cart::z][Cart::s][Cart::xz];
R[Cart::y][Cart::x][Cart::xz]=R[Cart::xy][Cart::s][Cart::xz]+amb0*R[Cart::y][Cart::s][Cart::xz];
R[Cart::x][Cart::x][Cart::xz]=R[Cart::xx][Cart::s][Cart::xz]+amb0*R[Cart::x][Cart::s][Cart::xz];
R[Cart::z][Cart::x][Cart::xz]=R[Cart::xz][Cart::s][Cart::xz]+amb0*R[Cart::z][Cart::s][Cart::xz];
R[Cart::y][Cart::z][Cart::xz]=R[Cart::yz][Cart::s][Cart::xz]+amb2*R[Cart::y][Cart::s][Cart::xz];
R[Cart::x][Cart::z][Cart::xz]=R[Cart::xz][Cart::s][Cart::xz]+amb2*R[Cart::x][Cart::s][Cart::xz];
R[Cart::z][Cart::z][Cart::xz]=R[Cart::zz][Cart::s][Cart::xz]+amb2*R[Cart::z][Cart::s][Cart::xz];
R[Cart::y][Cart::y][Cart::zz]=R[Cart::yy][Cart::s][Cart::zz]+amb1*R[Cart::y][Cart::s][Cart::zz];
R[Cart::x][Cart::y][Cart::zz]=R[Cart::xy][Cart::s][Cart::zz]+amb1*R[Cart::x][Cart::s][Cart::zz];
R[Cart::z][Cart::y][Cart::zz]=R[Cart::yz][Cart::s][Cart::zz]+amb1*R[Cart::z][Cart::s][Cart::zz];
R[Cart::y][Cart::x][Cart::zz]=R[Cart::xy][Cart::s][Cart::zz]+amb0*R[Cart::y][Cart::s][Cart::zz];
R[Cart::x][Cart::x][Cart::zz]=R[Cart::xx][Cart::s][Cart::zz]+amb0*R[Cart::x][Cart::s][Cart::zz];
R[Cart::z][Cart::x][Cart::zz]=R[Cart::xz][Cart::s][Cart::zz]+amb0*R[Cart::z][Cart::s][Cart::zz];
R[Cart::y][Cart::z][Cart::zz]=R[Cart::yz][Cart::s][Cart::zz]+amb2*R[Cart::y][Cart::s][Cart::zz];
R[Cart::x][Cart::z][Cart::zz]=R[Cart::xz][Cart::s][Cart::zz]+amb2*R[Cart::x][Cart::s][Cart::zz];
R[Cart::z][Cart::z][Cart::zz]=R[Cart::zz][Cart::s][Cart::zz]+amb2*R[Cart::z][Cart::s][Cart::zz];
}
//------------------------------------------------------

//Integral d - p - d
if (_lmax_beta>0 && _lmax_alpha>1 && _lmax_gamma>1){

R[Cart::yy][Cart::y][Cart::yy]=R[Cart::yyy][Cart::s][Cart::yy]+amb1*R[Cart::yy][Cart::s][Cart::yy];
R[Cart::xy][Cart::y][Cart::yy]=R[Cart::xyy][Cart::s][Cart::yy]+amb1*R[Cart::xy][Cart::s][Cart::yy];
R[Cart::yz][Cart::y][Cart::yy]=R[Cart::yyz][Cart::s][Cart::yy]+amb1*R[Cart::yz][Cart::s][Cart::yy];
R[Cart::xx][Cart::y][Cart::yy]=R[Cart::xxy][Cart::s][Cart::yy]+amb1*R[Cart::xx][Cart::s][Cart::yy];
R[Cart::xz][Cart::y][Cart::yy]=R[Cart::xyz][Cart::s][Cart::yy]+amb1*R[Cart::xz][Cart::s][Cart::yy];
R[Cart::zz][Cart::y][Cart::yy]=R[Cart::yzz][Cart::s][Cart::yy]+amb1*R[Cart::zz][Cart::s][Cart::yy];
R[Cart::yy][Cart::x][Cart::yy]=R[Cart::xyy][Cart::s][Cart::yy]+amb0*R[Cart::yy][Cart::s][Cart::yy];
R[Cart::xy][Cart::x][Cart::yy]=R[Cart::xxy][Cart::s][Cart::yy]+amb0*R[Cart::xy][Cart::s][Cart::yy];
R[Cart::yz][Cart::x][Cart::yy]=R[Cart::xyz][Cart::s][Cart::yy]+amb0*R[Cart::yz][Cart::s][Cart::yy];
R[Cart::xx][Cart::x][Cart::yy]=R[Cart::xxx][Cart::s][Cart::yy]+amb0*R[Cart::xx][Cart::s][Cart::yy];
R[Cart::xz][Cart::x][Cart::yy]=R[Cart::xxz][Cart::s][Cart::yy]+amb0*R[Cart::xz][Cart::s][Cart::yy];
R[Cart::zz][Cart::x][Cart::yy]=R[Cart::xzz][Cart::s][Cart::yy]+amb0*R[Cart::zz][Cart::s][Cart::yy];
R[Cart::yy][Cart::z][Cart::yy]=R[Cart::yyz][Cart::s][Cart::yy]+amb2*R[Cart::yy][Cart::s][Cart::yy];
R[Cart::xy][Cart::z][Cart::yy]=R[Cart::xyz][Cart::s][Cart::yy]+amb2*R[Cart::xy][Cart::s][Cart::yy];
R[Cart::yz][Cart::z][Cart::yy]=R[Cart::yzz][Cart::s][Cart::yy]+amb2*R[Cart::yz][Cart::s][Cart::yy];
R[Cart::xx][Cart::z][Cart::yy]=R[Cart::xxz][Cart::s][Cart::yy]+amb2*R[Cart::xx][Cart::s][Cart::yy];
R[Cart::xz][Cart::z][Cart::yy]=R[Cart::xzz][Cart::s][Cart::yy]+amb2*R[Cart::xz][Cart::s][Cart::yy];
R[Cart::zz][Cart::z][Cart::yy]=R[Cart::zzz][Cart::s][Cart::yy]+amb2*R[Cart::zz][Cart::s][Cart::yy];
R[Cart::yy][Cart::y][Cart::xy]=R[Cart::yyy][Cart::s][Cart::xy]+amb1*R[Cart::yy][Cart::s][Cart::xy];
R[Cart::xy][Cart::y][Cart::xy]=R[Cart::xyy][Cart::s][Cart::xy]+amb1*R[Cart::xy][Cart::s][Cart::xy];
R[Cart::yz][Cart::y][Cart::xy]=R[Cart::yyz][Cart::s][Cart::xy]+amb1*R[Cart::yz][Cart::s][Cart::xy];
R[Cart::xx][Cart::y][Cart::xy]=R[Cart::xxy][Cart::s][Cart::xy]+amb1*R[Cart::xx][Cart::s][Cart::xy];
R[Cart::xz][Cart::y][Cart::xy]=R[Cart::xyz][Cart::s][Cart::xy]+amb1*R[Cart::xz][Cart::s][Cart::xy];
R[Cart::zz][Cart::y][Cart::xy]=R[Cart::yzz][Cart::s][Cart::xy]+amb1*R[Cart::zz][Cart::s][Cart::xy];
R[Cart::yy][Cart::x][Cart::xy]=R[Cart::xyy][Cart::s][Cart::xy]+amb0*R[Cart::yy][Cart::s][Cart::xy];
R[Cart::xy][Cart::x][Cart::xy]=R[Cart::xxy][Cart::s][Cart::xy]+amb0*R[Cart::xy][Cart::s][Cart::xy];
R[Cart::yz][Cart::x][Cart::xy]=R[Cart::xyz][Cart::s][Cart::xy]+amb0*R[Cart::yz][Cart::s][Cart::xy];
R[Cart::xx][Cart::x][Cart::xy]=R[Cart::xxx][Cart::s][Cart::xy]+amb0*R[Cart::xx][Cart::s][Cart::xy];
R[Cart::xz][Cart::x][Cart::xy]=R[Cart::xxz][Cart::s][Cart::xy]+amb0*R[Cart::xz][Cart::s][Cart::xy];
R[Cart::zz][Cart::x][Cart::xy]=R[Cart::xzz][Cart::s][Cart::xy]+amb0*R[Cart::zz][Cart::s][Cart::xy];
R[Cart::yy][Cart::z][Cart::xy]=R[Cart::yyz][Cart::s][Cart::xy]+amb2*R[Cart::yy][Cart::s][Cart::xy];
R[Cart::xy][Cart::z][Cart::xy]=R[Cart::xyz][Cart::s][Cart::xy]+amb2*R[Cart::xy][Cart::s][Cart::xy];
R[Cart::yz][Cart::z][Cart::xy]=R[Cart::yzz][Cart::s][Cart::xy]+amb2*R[Cart::yz][Cart::s][Cart::xy];
R[Cart::xx][Cart::z][Cart::xy]=R[Cart::xxz][Cart::s][Cart::xy]+amb2*R[Cart::xx][Cart::s][Cart::xy];
R[Cart::xz][Cart::z][Cart::xy]=R[Cart::xzz][Cart::s][Cart::xy]+amb2*R[Cart::xz][Cart::s][Cart::xy];
R[Cart::zz][Cart::z][Cart::xy]=R[Cart::zzz][Cart::s][Cart::xy]+amb2*R[Cart::zz][Cart::s][Cart::xy];
R[Cart::yy][Cart::y][Cart::yz]=R[Cart::yyy][Cart::s][Cart::yz]+amb1*R[Cart::yy][Cart::s][Cart::yz];
R[Cart::xy][Cart::y][Cart::yz]=R[Cart::xyy][Cart::s][Cart::yz]+amb1*R[Cart::xy][Cart::s][Cart::yz];
R[Cart::yz][Cart::y][Cart::yz]=R[Cart::yyz][Cart::s][Cart::yz]+amb1*R[Cart::yz][Cart::s][Cart::yz];
R[Cart::xx][Cart::y][Cart::yz]=R[Cart::xxy][Cart::s][Cart::yz]+amb1*R[Cart::xx][Cart::s][Cart::yz];
R[Cart::xz][Cart::y][Cart::yz]=R[Cart::xyz][Cart::s][Cart::yz]+amb1*R[Cart::xz][Cart::s][Cart::yz];
R[Cart::zz][Cart::y][Cart::yz]=R[Cart::yzz][Cart::s][Cart::yz]+amb1*R[Cart::zz][Cart::s][Cart::yz];
R[Cart::yy][Cart::x][Cart::yz]=R[Cart::xyy][Cart::s][Cart::yz]+amb0*R[Cart::yy][Cart::s][Cart::yz];
R[Cart::xy][Cart::x][Cart::yz]=R[Cart::xxy][Cart::s][Cart::yz]+amb0*R[Cart::xy][Cart::s][Cart::yz];
R[Cart::yz][Cart::x][Cart::yz]=R[Cart::xyz][Cart::s][Cart::yz]+amb0*R[Cart::yz][Cart::s][Cart::yz];
R[Cart::xx][Cart::x][Cart::yz]=R[Cart::xxx][Cart::s][Cart::yz]+amb0*R[Cart::xx][Cart::s][Cart::yz];
R[Cart::xz][Cart::x][Cart::yz]=R[Cart::xxz][Cart::s][Cart::yz]+amb0*R[Cart::xz][Cart::s][Cart::yz];
R[Cart::zz][Cart::x][Cart::yz]=R[Cart::xzz][Cart::s][Cart::yz]+amb0*R[Cart::zz][Cart::s][Cart::yz];
R[Cart::yy][Cart::z][Cart::yz]=R[Cart::yyz][Cart::s][Cart::yz]+amb2*R[Cart::yy][Cart::s][Cart::yz];
R[Cart::xy][Cart::z][Cart::yz]=R[Cart::xyz][Cart::s][Cart::yz]+amb2*R[Cart::xy][Cart::s][Cart::yz];
R[Cart::yz][Cart::z][Cart::yz]=R[Cart::yzz][Cart::s][Cart::yz]+amb2*R[Cart::yz][Cart::s][Cart::yz];
R[Cart::xx][Cart::z][Cart::yz]=R[Cart::xxz][Cart::s][Cart::yz]+amb2*R[Cart::xx][Cart::s][Cart::yz];
R[Cart::xz][Cart::z][Cart::yz]=R[Cart::xzz][Cart::s][Cart::yz]+amb2*R[Cart::xz][Cart::s][Cart::yz];
R[Cart::zz][Cart::z][Cart::yz]=R[Cart::zzz][Cart::s][Cart::yz]+amb2*R[Cart::zz][Cart::s][Cart::yz];
R[Cart::yy][Cart::y][Cart::xx]=R[Cart::yyy][Cart::s][Cart::xx]+amb1*R[Cart::yy][Cart::s][Cart::xx];
R[Cart::xy][Cart::y][Cart::xx]=R[Cart::xyy][Cart::s][Cart::xx]+amb1*R[Cart::xy][Cart::s][Cart::xx];
R[Cart::yz][Cart::y][Cart::xx]=R[Cart::yyz][Cart::s][Cart::xx]+amb1*R[Cart::yz][Cart::s][Cart::xx];
R[Cart::xx][Cart::y][Cart::xx]=R[Cart::xxy][Cart::s][Cart::xx]+amb1*R[Cart::xx][Cart::s][Cart::xx];
R[Cart::xz][Cart::y][Cart::xx]=R[Cart::xyz][Cart::s][Cart::xx]+amb1*R[Cart::xz][Cart::s][Cart::xx];
R[Cart::zz][Cart::y][Cart::xx]=R[Cart::yzz][Cart::s][Cart::xx]+amb1*R[Cart::zz][Cart::s][Cart::xx];
R[Cart::yy][Cart::x][Cart::xx]=R[Cart::xyy][Cart::s][Cart::xx]+amb0*R[Cart::yy][Cart::s][Cart::xx];
R[Cart::xy][Cart::x][Cart::xx]=R[Cart::xxy][Cart::s][Cart::xx]+amb0*R[Cart::xy][Cart::s][Cart::xx];
R[Cart::yz][Cart::x][Cart::xx]=R[Cart::xyz][Cart::s][Cart::xx]+amb0*R[Cart::yz][Cart::s][Cart::xx];
R[Cart::xx][Cart::x][Cart::xx]=R[Cart::xxx][Cart::s][Cart::xx]+amb0*R[Cart::xx][Cart::s][Cart::xx];
R[Cart::xz][Cart::x][Cart::xx]=R[Cart::xxz][Cart::s][Cart::xx]+amb0*R[Cart::xz][Cart::s][Cart::xx];
R[Cart::zz][Cart::x][Cart::xx]=R[Cart::xzz][Cart::s][Cart::xx]+amb0*R[Cart::zz][Cart::s][Cart::xx];
R[Cart::yy][Cart::z][Cart::xx]=R[Cart::yyz][Cart::s][Cart::xx]+amb2*R[Cart::yy][Cart::s][Cart::xx];
R[Cart::xy][Cart::z][Cart::xx]=R[Cart::xyz][Cart::s][Cart::xx]+amb2*R[Cart::xy][Cart::s][Cart::xx];
R[Cart::yz][Cart::z][Cart::xx]=R[Cart::yzz][Cart::s][Cart::xx]+amb2*R[Cart::yz][Cart::s][Cart::xx];
R[Cart::xx][Cart::z][Cart::xx]=R[Cart::xxz][Cart::s][Cart::xx]+amb2*R[Cart::xx][Cart::s][Cart::xx];
R[Cart::xz][Cart::z][Cart::xx]=R[Cart::xzz][Cart::s][Cart::xx]+amb2*R[Cart::xz][Cart::s][Cart::xx];
R[Cart::zz][Cart::z][Cart::xx]=R[Cart::zzz][Cart::s][Cart::xx]+amb2*R[Cart::zz][Cart::s][Cart::xx];
R[Cart::yy][Cart::y][Cart::xz]=R[Cart::yyy][Cart::s][Cart::xz]+amb1*R[Cart::yy][Cart::s][Cart::xz];
R[Cart::xy][Cart::y][Cart::xz]=R[Cart::xyy][Cart::s][Cart::xz]+amb1*R[Cart::xy][Cart::s][Cart::xz];
R[Cart::yz][Cart::y][Cart::xz]=R[Cart::yyz][Cart::s][Cart::xz]+amb1*R[Cart::yz][Cart::s][Cart::xz];
R[Cart::xx][Cart::y][Cart::xz]=R[Cart::xxy][Cart::s][Cart::xz]+amb1*R[Cart::xx][Cart::s][Cart::xz];
R[Cart::xz][Cart::y][Cart::xz]=R[Cart::xyz][Cart::s][Cart::xz]+amb1*R[Cart::xz][Cart::s][Cart::xz];
R[Cart::zz][Cart::y][Cart::xz]=R[Cart::yzz][Cart::s][Cart::xz]+amb1*R[Cart::zz][Cart::s][Cart::xz];
R[Cart::yy][Cart::x][Cart::xz]=R[Cart::xyy][Cart::s][Cart::xz]+amb0*R[Cart::yy][Cart::s][Cart::xz];
R[Cart::xy][Cart::x][Cart::xz]=R[Cart::xxy][Cart::s][Cart::xz]+amb0*R[Cart::xy][Cart::s][Cart::xz];
R[Cart::yz][Cart::x][Cart::xz]=R[Cart::xyz][Cart::s][Cart::xz]+amb0*R[Cart::yz][Cart::s][Cart::xz];
R[Cart::xx][Cart::x][Cart::xz]=R[Cart::xxx][Cart::s][Cart::xz]+amb0*R[Cart::xx][Cart::s][Cart::xz];
R[Cart::xz][Cart::x][Cart::xz]=R[Cart::xxz][Cart::s][Cart::xz]+amb0*R[Cart::xz][Cart::s][Cart::xz];
R[Cart::zz][Cart::x][Cart::xz]=R[Cart::xzz][Cart::s][Cart::xz]+amb0*R[Cart::zz][Cart::s][Cart::xz];
R[Cart::yy][Cart::z][Cart::xz]=R[Cart::yyz][Cart::s][Cart::xz]+amb2*R[Cart::yy][Cart::s][Cart::xz];
R[Cart::xy][Cart::z][Cart::xz]=R[Cart::xyz][Cart::s][Cart::xz]+amb2*R[Cart::xy][Cart::s][Cart::xz];
R[Cart::yz][Cart::z][Cart::xz]=R[Cart::yzz][Cart::s][Cart::xz]+amb2*R[Cart::yz][Cart::s][Cart::xz];
R[Cart::xx][Cart::z][Cart::xz]=R[Cart::xxz][Cart::s][Cart::xz]+amb2*R[Cart::xx][Cart::s][Cart::xz];
R[Cart::xz][Cart::z][Cart::xz]=R[Cart::xzz][Cart::s][Cart::xz]+amb2*R[Cart::xz][Cart::s][Cart::xz];
R[Cart::zz][Cart::z][Cart::xz]=R[Cart::zzz][Cart::s][Cart::xz]+amb2*R[Cart::zz][Cart::s][Cart::xz];
R[Cart::yy][Cart::y][Cart::zz]=R[Cart::yyy][Cart::s][Cart::zz]+amb1*R[Cart::yy][Cart::s][Cart::zz];
R[Cart::xy][Cart::y][Cart::zz]=R[Cart::xyy][Cart::s][Cart::zz]+amb1*R[Cart::xy][Cart::s][Cart::zz];
R[Cart::yz][Cart::y][Cart::zz]=R[Cart::yyz][Cart::s][Cart::zz]+amb1*R[Cart::yz][Cart::s][Cart::zz];
R[Cart::xx][Cart::y][Cart::zz]=R[Cart::xxy][Cart::s][Cart::zz]+amb1*R[Cart::xx][Cart::s][Cart::zz];
R[Cart::xz][Cart::y][Cart::zz]=R[Cart::xyz][Cart::s][Cart::zz]+amb1*R[Cart::xz][Cart::s][Cart::zz];
R[Cart::zz][Cart::y][Cart::zz]=R[Cart::yzz][Cart::s][Cart::zz]+amb1*R[Cart::zz][Cart::s][Cart::zz];
R[Cart::yy][Cart::x][Cart::zz]=R[Cart::xyy][Cart::s][Cart::zz]+amb0*R[Cart::yy][Cart::s][Cart::zz];
R[Cart::xy][Cart::x][Cart::zz]=R[Cart::xxy][Cart::s][Cart::zz]+amb0*R[Cart::xy][Cart::s][Cart::zz];
R[Cart::yz][Cart::x][Cart::zz]=R[Cart::xyz][Cart::s][Cart::zz]+amb0*R[Cart::yz][Cart::s][Cart::zz];
R[Cart::xx][Cart::x][Cart::zz]=R[Cart::xxx][Cart::s][Cart::zz]+amb0*R[Cart::xx][Cart::s][Cart::zz];
R[Cart::xz][Cart::x][Cart::zz]=R[Cart::xxz][Cart::s][Cart::zz]+amb0*R[Cart::xz][Cart::s][Cart::zz];
R[Cart::zz][Cart::x][Cart::zz]=R[Cart::xzz][Cart::s][Cart::zz]+amb0*R[Cart::zz][Cart::s][Cart::zz];
R[Cart::yy][Cart::z][Cart::zz]=R[Cart::yyz][Cart::s][Cart::zz]+amb2*R[Cart::yy][Cart::s][Cart::zz];
R[Cart::xy][Cart::z][Cart::zz]=R[Cart::xyz][Cart::s][Cart::zz]+amb2*R[Cart::xy][Cart::s][Cart::zz];
R[Cart::yz][Cart::z][Cart::zz]=R[Cart::yzz][Cart::s][Cart::zz]+amb2*R[Cart::yz][Cart::s][Cart::zz];
R[Cart::xx][Cart::z][Cart::zz]=R[Cart::xxz][Cart::s][Cart::zz]+amb2*R[Cart::xx][Cart::s][Cart::zz];
R[Cart::xz][Cart::z][Cart::zz]=R[Cart::xzz][Cart::s][Cart::zz]+amb2*R[Cart::xz][Cart::s][Cart::zz];
R[Cart::zz][Cart::z][Cart::zz]=R[Cart::zzz][Cart::s][Cart::zz]+amb2*R[Cart::zz][Cart::s][Cart::zz];
}
//------------------------------------------------------

//Integral f - p - d
if (_lmax_beta>0 && _lmax_alpha>2 && _lmax_gamma>1){

R[Cart::yyy][Cart::y][Cart::yy]=R[Cart::yyyy][Cart::s][Cart::yy]+amb1*R[Cart::yyy][Cart::s][Cart::yy];
R[Cart::xyy][Cart::y][Cart::yy]=R[Cart::xyyy][Cart::s][Cart::yy]+amb1*R[Cart::xyy][Cart::s][Cart::yy];
R[Cart::yyz][Cart::y][Cart::yy]=R[Cart::yyyz][Cart::s][Cart::yy]+amb1*R[Cart::yyz][Cart::s][Cart::yy];
R[Cart::xxy][Cart::y][Cart::yy]=R[Cart::xxyy][Cart::s][Cart::yy]+amb1*R[Cart::xxy][Cart::s][Cart::yy];
R[Cart::xyz][Cart::y][Cart::yy]=R[Cart::xyyz][Cart::s][Cart::yy]+amb1*R[Cart::xyz][Cart::s][Cart::yy];
R[Cart::yzz][Cart::y][Cart::yy]=R[Cart::yyzz][Cart::s][Cart::yy]+amb1*R[Cart::yzz][Cart::s][Cart::yy];
R[Cart::xxx][Cart::y][Cart::yy]=R[Cart::xxxy][Cart::s][Cart::yy]+amb1*R[Cart::xxx][Cart::s][Cart::yy];
R[Cart::xxz][Cart::y][Cart::yy]=R[Cart::xxyz][Cart::s][Cart::yy]+amb1*R[Cart::xxz][Cart::s][Cart::yy];
R[Cart::xzz][Cart::y][Cart::yy]=R[Cart::xyzz][Cart::s][Cart::yy]+amb1*R[Cart::xzz][Cart::s][Cart::yy];
R[Cart::zzz][Cart::y][Cart::yy]=R[Cart::yzzz][Cart::s][Cart::yy]+amb1*R[Cart::zzz][Cart::s][Cart::yy];
R[Cart::yyy][Cart::x][Cart::yy]=R[Cart::xyyy][Cart::s][Cart::yy]+amb0*R[Cart::yyy][Cart::s][Cart::yy];
R[Cart::xyy][Cart::x][Cart::yy]=R[Cart::xxyy][Cart::s][Cart::yy]+amb0*R[Cart::xyy][Cart::s][Cart::yy];
R[Cart::yyz][Cart::x][Cart::yy]=R[Cart::xyyz][Cart::s][Cart::yy]+amb0*R[Cart::yyz][Cart::s][Cart::yy];
R[Cart::xxy][Cart::x][Cart::yy]=R[Cart::xxxy][Cart::s][Cart::yy]+amb0*R[Cart::xxy][Cart::s][Cart::yy];
R[Cart::xyz][Cart::x][Cart::yy]=R[Cart::xxyz][Cart::s][Cart::yy]+amb0*R[Cart::xyz][Cart::s][Cart::yy];
R[Cart::yzz][Cart::x][Cart::yy]=R[Cart::xyzz][Cart::s][Cart::yy]+amb0*R[Cart::yzz][Cart::s][Cart::yy];
R[Cart::xxx][Cart::x][Cart::yy]=R[Cart::xxxx][Cart::s][Cart::yy]+amb0*R[Cart::xxx][Cart::s][Cart::yy];
R[Cart::xxz][Cart::x][Cart::yy]=R[Cart::xxxz][Cart::s][Cart::yy]+amb0*R[Cart::xxz][Cart::s][Cart::yy];
R[Cart::xzz][Cart::x][Cart::yy]=R[Cart::xxzz][Cart::s][Cart::yy]+amb0*R[Cart::xzz][Cart::s][Cart::yy];
R[Cart::zzz][Cart::x][Cart::yy]=R[Cart::xzzz][Cart::s][Cart::yy]+amb0*R[Cart::zzz][Cart::s][Cart::yy];
R[Cart::yyy][Cart::z][Cart::yy]=R[Cart::yyyz][Cart::s][Cart::yy]+amb2*R[Cart::yyy][Cart::s][Cart::yy];
R[Cart::xyy][Cart::z][Cart::yy]=R[Cart::xyyz][Cart::s][Cart::yy]+amb2*R[Cart::xyy][Cart::s][Cart::yy];
R[Cart::yyz][Cart::z][Cart::yy]=R[Cart::yyzz][Cart::s][Cart::yy]+amb2*R[Cart::yyz][Cart::s][Cart::yy];
R[Cart::xxy][Cart::z][Cart::yy]=R[Cart::xxyz][Cart::s][Cart::yy]+amb2*R[Cart::xxy][Cart::s][Cart::yy];
R[Cart::xyz][Cart::z][Cart::yy]=R[Cart::xyzz][Cart::s][Cart::yy]+amb2*R[Cart::xyz][Cart::s][Cart::yy];
R[Cart::yzz][Cart::z][Cart::yy]=R[Cart::yzzz][Cart::s][Cart::yy]+amb2*R[Cart::yzz][Cart::s][Cart::yy];
R[Cart::xxx][Cart::z][Cart::yy]=R[Cart::xxxz][Cart::s][Cart::yy]+amb2*R[Cart::xxx][Cart::s][Cart::yy];
R[Cart::xxz][Cart::z][Cart::yy]=R[Cart::xxzz][Cart::s][Cart::yy]+amb2*R[Cart::xxz][Cart::s][Cart::yy];
R[Cart::xzz][Cart::z][Cart::yy]=R[Cart::xzzz][Cart::s][Cart::yy]+amb2*R[Cart::xzz][Cart::s][Cart::yy];
R[Cart::zzz][Cart::z][Cart::yy]=R[Cart::zzzz][Cart::s][Cart::yy]+amb2*R[Cart::zzz][Cart::s][Cart::yy];
R[Cart::yyy][Cart::y][Cart::xy]=R[Cart::yyyy][Cart::s][Cart::xy]+amb1*R[Cart::yyy][Cart::s][Cart::xy];
R[Cart::xyy][Cart::y][Cart::xy]=R[Cart::xyyy][Cart::s][Cart::xy]+amb1*R[Cart::xyy][Cart::s][Cart::xy];
R[Cart::yyz][Cart::y][Cart::xy]=R[Cart::yyyz][Cart::s][Cart::xy]+amb1*R[Cart::yyz][Cart::s][Cart::xy];
R[Cart::xxy][Cart::y][Cart::xy]=R[Cart::xxyy][Cart::s][Cart::xy]+amb1*R[Cart::xxy][Cart::s][Cart::xy];
R[Cart::xyz][Cart::y][Cart::xy]=R[Cart::xyyz][Cart::s][Cart::xy]+amb1*R[Cart::xyz][Cart::s][Cart::xy];
R[Cart::yzz][Cart::y][Cart::xy]=R[Cart::yyzz][Cart::s][Cart::xy]+amb1*R[Cart::yzz][Cart::s][Cart::xy];
R[Cart::xxx][Cart::y][Cart::xy]=R[Cart::xxxy][Cart::s][Cart::xy]+amb1*R[Cart::xxx][Cart::s][Cart::xy];
R[Cart::xxz][Cart::y][Cart::xy]=R[Cart::xxyz][Cart::s][Cart::xy]+amb1*R[Cart::xxz][Cart::s][Cart::xy];
R[Cart::xzz][Cart::y][Cart::xy]=R[Cart::xyzz][Cart::s][Cart::xy]+amb1*R[Cart::xzz][Cart::s][Cart::xy];
R[Cart::zzz][Cart::y][Cart::xy]=R[Cart::yzzz][Cart::s][Cart::xy]+amb1*R[Cart::zzz][Cart::s][Cart::xy];
R[Cart::yyy][Cart::x][Cart::xy]=R[Cart::xyyy][Cart::s][Cart::xy]+amb0*R[Cart::yyy][Cart::s][Cart::xy];
R[Cart::xyy][Cart::x][Cart::xy]=R[Cart::xxyy][Cart::s][Cart::xy]+amb0*R[Cart::xyy][Cart::s][Cart::xy];
R[Cart::yyz][Cart::x][Cart::xy]=R[Cart::xyyz][Cart::s][Cart::xy]+amb0*R[Cart::yyz][Cart::s][Cart::xy];
R[Cart::xxy][Cart::x][Cart::xy]=R[Cart::xxxy][Cart::s][Cart::xy]+amb0*R[Cart::xxy][Cart::s][Cart::xy];
R[Cart::xyz][Cart::x][Cart::xy]=R[Cart::xxyz][Cart::s][Cart::xy]+amb0*R[Cart::xyz][Cart::s][Cart::xy];
R[Cart::yzz][Cart::x][Cart::xy]=R[Cart::xyzz][Cart::s][Cart::xy]+amb0*R[Cart::yzz][Cart::s][Cart::xy];
R[Cart::xxx][Cart::x][Cart::xy]=R[Cart::xxxx][Cart::s][Cart::xy]+amb0*R[Cart::xxx][Cart::s][Cart::xy];
R[Cart::xxz][Cart::x][Cart::xy]=R[Cart::xxxz][Cart::s][Cart::xy]+amb0*R[Cart::xxz][Cart::s][Cart::xy];
R[Cart::xzz][Cart::x][Cart::xy]=R[Cart::xxzz][Cart::s][Cart::xy]+amb0*R[Cart::xzz][Cart::s][Cart::xy];
R[Cart::zzz][Cart::x][Cart::xy]=R[Cart::xzzz][Cart::s][Cart::xy]+amb0*R[Cart::zzz][Cart::s][Cart::xy];
R[Cart::yyy][Cart::z][Cart::xy]=R[Cart::yyyz][Cart::s][Cart::xy]+amb2*R[Cart::yyy][Cart::s][Cart::xy];
R[Cart::xyy][Cart::z][Cart::xy]=R[Cart::xyyz][Cart::s][Cart::xy]+amb2*R[Cart::xyy][Cart::s][Cart::xy];
R[Cart::yyz][Cart::z][Cart::xy]=R[Cart::yyzz][Cart::s][Cart::xy]+amb2*R[Cart::yyz][Cart::s][Cart::xy];
R[Cart::xxy][Cart::z][Cart::xy]=R[Cart::xxyz][Cart::s][Cart::xy]+amb2*R[Cart::xxy][Cart::s][Cart::xy];
R[Cart::xyz][Cart::z][Cart::xy]=R[Cart::xyzz][Cart::s][Cart::xy]+amb2*R[Cart::xyz][Cart::s][Cart::xy];
R[Cart::yzz][Cart::z][Cart::xy]=R[Cart::yzzz][Cart::s][Cart::xy]+amb2*R[Cart::yzz][Cart::s][Cart::xy];
R[Cart::xxx][Cart::z][Cart::xy]=R[Cart::xxxz][Cart::s][Cart::xy]+amb2*R[Cart::xxx][Cart::s][Cart::xy];
R[Cart::xxz][Cart::z][Cart::xy]=R[Cart::xxzz][Cart::s][Cart::xy]+amb2*R[Cart::xxz][Cart::s][Cart::xy];
R[Cart::xzz][Cart::z][Cart::xy]=R[Cart::xzzz][Cart::s][Cart::xy]+amb2*R[Cart::xzz][Cart::s][Cart::xy];
R[Cart::zzz][Cart::z][Cart::xy]=R[Cart::zzzz][Cart::s][Cart::xy]+amb2*R[Cart::zzz][Cart::s][Cart::xy];
R[Cart::yyy][Cart::y][Cart::yz]=R[Cart::yyyy][Cart::s][Cart::yz]+amb1*R[Cart::yyy][Cart::s][Cart::yz];
R[Cart::xyy][Cart::y][Cart::yz]=R[Cart::xyyy][Cart::s][Cart::yz]+amb1*R[Cart::xyy][Cart::s][Cart::yz];
R[Cart::yyz][Cart::y][Cart::yz]=R[Cart::yyyz][Cart::s][Cart::yz]+amb1*R[Cart::yyz][Cart::s][Cart::yz];
R[Cart::xxy][Cart::y][Cart::yz]=R[Cart::xxyy][Cart::s][Cart::yz]+amb1*R[Cart::xxy][Cart::s][Cart::yz];
R[Cart::xyz][Cart::y][Cart::yz]=R[Cart::xyyz][Cart::s][Cart::yz]+amb1*R[Cart::xyz][Cart::s][Cart::yz];
R[Cart::yzz][Cart::y][Cart::yz]=R[Cart::yyzz][Cart::s][Cart::yz]+amb1*R[Cart::yzz][Cart::s][Cart::yz];
R[Cart::xxx][Cart::y][Cart::yz]=R[Cart::xxxy][Cart::s][Cart::yz]+amb1*R[Cart::xxx][Cart::s][Cart::yz];
R[Cart::xxz][Cart::y][Cart::yz]=R[Cart::xxyz][Cart::s][Cart::yz]+amb1*R[Cart::xxz][Cart::s][Cart::yz];
R[Cart::xzz][Cart::y][Cart::yz]=R[Cart::xyzz][Cart::s][Cart::yz]+amb1*R[Cart::xzz][Cart::s][Cart::yz];
R[Cart::zzz][Cart::y][Cart::yz]=R[Cart::yzzz][Cart::s][Cart::yz]+amb1*R[Cart::zzz][Cart::s][Cart::yz];
R[Cart::yyy][Cart::x][Cart::yz]=R[Cart::xyyy][Cart::s][Cart::yz]+amb0*R[Cart::yyy][Cart::s][Cart::yz];
R[Cart::xyy][Cart::x][Cart::yz]=R[Cart::xxyy][Cart::s][Cart::yz]+amb0*R[Cart::xyy][Cart::s][Cart::yz];
R[Cart::yyz][Cart::x][Cart::yz]=R[Cart::xyyz][Cart::s][Cart::yz]+amb0*R[Cart::yyz][Cart::s][Cart::yz];
R[Cart::xxy][Cart::x][Cart::yz]=R[Cart::xxxy][Cart::s][Cart::yz]+amb0*R[Cart::xxy][Cart::s][Cart::yz];
R[Cart::xyz][Cart::x][Cart::yz]=R[Cart::xxyz][Cart::s][Cart::yz]+amb0*R[Cart::xyz][Cart::s][Cart::yz];
R[Cart::yzz][Cart::x][Cart::yz]=R[Cart::xyzz][Cart::s][Cart::yz]+amb0*R[Cart::yzz][Cart::s][Cart::yz];
R[Cart::xxx][Cart::x][Cart::yz]=R[Cart::xxxx][Cart::s][Cart::yz]+amb0*R[Cart::xxx][Cart::s][Cart::yz];
R[Cart::xxz][Cart::x][Cart::yz]=R[Cart::xxxz][Cart::s][Cart::yz]+amb0*R[Cart::xxz][Cart::s][Cart::yz];
R[Cart::xzz][Cart::x][Cart::yz]=R[Cart::xxzz][Cart::s][Cart::yz]+amb0*R[Cart::xzz][Cart::s][Cart::yz];
R[Cart::zzz][Cart::x][Cart::yz]=R[Cart::xzzz][Cart::s][Cart::yz]+amb0*R[Cart::zzz][Cart::s][Cart::yz];
R[Cart::yyy][Cart::z][Cart::yz]=R[Cart::yyyz][Cart::s][Cart::yz]+amb2*R[Cart::yyy][Cart::s][Cart::yz];
R[Cart::xyy][Cart::z][Cart::yz]=R[Cart::xyyz][Cart::s][Cart::yz]+amb2*R[Cart::xyy][Cart::s][Cart::yz];
R[Cart::yyz][Cart::z][Cart::yz]=R[Cart::yyzz][Cart::s][Cart::yz]+amb2*R[Cart::yyz][Cart::s][Cart::yz];
R[Cart::xxy][Cart::z][Cart::yz]=R[Cart::xxyz][Cart::s][Cart::yz]+amb2*R[Cart::xxy][Cart::s][Cart::yz];
R[Cart::xyz][Cart::z][Cart::yz]=R[Cart::xyzz][Cart::s][Cart::yz]+amb2*R[Cart::xyz][Cart::s][Cart::yz];
R[Cart::yzz][Cart::z][Cart::yz]=R[Cart::yzzz][Cart::s][Cart::yz]+amb2*R[Cart::yzz][Cart::s][Cart::yz];
R[Cart::xxx][Cart::z][Cart::yz]=R[Cart::xxxz][Cart::s][Cart::yz]+amb2*R[Cart::xxx][Cart::s][Cart::yz];
R[Cart::xxz][Cart::z][Cart::yz]=R[Cart::xxzz][Cart::s][Cart::yz]+amb2*R[Cart::xxz][Cart::s][Cart::yz];
R[Cart::xzz][Cart::z][Cart::yz]=R[Cart::xzzz][Cart::s][Cart::yz]+amb2*R[Cart::xzz][Cart::s][Cart::yz];
R[Cart::zzz][Cart::z][Cart::yz]=R[Cart::zzzz][Cart::s][Cart::yz]+amb2*R[Cart::zzz][Cart::s][Cart::yz];
R[Cart::yyy][Cart::y][Cart::xx]=R[Cart::yyyy][Cart::s][Cart::xx]+amb1*R[Cart::yyy][Cart::s][Cart::xx];
R[Cart::xyy][Cart::y][Cart::xx]=R[Cart::xyyy][Cart::s][Cart::xx]+amb1*R[Cart::xyy][Cart::s][Cart::xx];
R[Cart::yyz][Cart::y][Cart::xx]=R[Cart::yyyz][Cart::s][Cart::xx]+amb1*R[Cart::yyz][Cart::s][Cart::xx];
R[Cart::xxy][Cart::y][Cart::xx]=R[Cart::xxyy][Cart::s][Cart::xx]+amb1*R[Cart::xxy][Cart::s][Cart::xx];
R[Cart::xyz][Cart::y][Cart::xx]=R[Cart::xyyz][Cart::s][Cart::xx]+amb1*R[Cart::xyz][Cart::s][Cart::xx];
R[Cart::yzz][Cart::y][Cart::xx]=R[Cart::yyzz][Cart::s][Cart::xx]+amb1*R[Cart::yzz][Cart::s][Cart::xx];
R[Cart::xxx][Cart::y][Cart::xx]=R[Cart::xxxy][Cart::s][Cart::xx]+amb1*R[Cart::xxx][Cart::s][Cart::xx];
R[Cart::xxz][Cart::y][Cart::xx]=R[Cart::xxyz][Cart::s][Cart::xx]+amb1*R[Cart::xxz][Cart::s][Cart::xx];
R[Cart::xzz][Cart::y][Cart::xx]=R[Cart::xyzz][Cart::s][Cart::xx]+amb1*R[Cart::xzz][Cart::s][Cart::xx];
R[Cart::zzz][Cart::y][Cart::xx]=R[Cart::yzzz][Cart::s][Cart::xx]+amb1*R[Cart::zzz][Cart::s][Cart::xx];
R[Cart::yyy][Cart::x][Cart::xx]=R[Cart::xyyy][Cart::s][Cart::xx]+amb0*R[Cart::yyy][Cart::s][Cart::xx];
R[Cart::xyy][Cart::x][Cart::xx]=R[Cart::xxyy][Cart::s][Cart::xx]+amb0*R[Cart::xyy][Cart::s][Cart::xx];
R[Cart::yyz][Cart::x][Cart::xx]=R[Cart::xyyz][Cart::s][Cart::xx]+amb0*R[Cart::yyz][Cart::s][Cart::xx];
R[Cart::xxy][Cart::x][Cart::xx]=R[Cart::xxxy][Cart::s][Cart::xx]+amb0*R[Cart::xxy][Cart::s][Cart::xx];
R[Cart::xyz][Cart::x][Cart::xx]=R[Cart::xxyz][Cart::s][Cart::xx]+amb0*R[Cart::xyz][Cart::s][Cart::xx];
R[Cart::yzz][Cart::x][Cart::xx]=R[Cart::xyzz][Cart::s][Cart::xx]+amb0*R[Cart::yzz][Cart::s][Cart::xx];
R[Cart::xxx][Cart::x][Cart::xx]=R[Cart::xxxx][Cart::s][Cart::xx]+amb0*R[Cart::xxx][Cart::s][Cart::xx];
R[Cart::xxz][Cart::x][Cart::xx]=R[Cart::xxxz][Cart::s][Cart::xx]+amb0*R[Cart::xxz][Cart::s][Cart::xx];
R[Cart::xzz][Cart::x][Cart::xx]=R[Cart::xxzz][Cart::s][Cart::xx]+amb0*R[Cart::xzz][Cart::s][Cart::xx];
R[Cart::zzz][Cart::x][Cart::xx]=R[Cart::xzzz][Cart::s][Cart::xx]+amb0*R[Cart::zzz][Cart::s][Cart::xx];
R[Cart::yyy][Cart::z][Cart::xx]=R[Cart::yyyz][Cart::s][Cart::xx]+amb2*R[Cart::yyy][Cart::s][Cart::xx];
R[Cart::xyy][Cart::z][Cart::xx]=R[Cart::xyyz][Cart::s][Cart::xx]+amb2*R[Cart::xyy][Cart::s][Cart::xx];
R[Cart::yyz][Cart::z][Cart::xx]=R[Cart::yyzz][Cart::s][Cart::xx]+amb2*R[Cart::yyz][Cart::s][Cart::xx];
R[Cart::xxy][Cart::z][Cart::xx]=R[Cart::xxyz][Cart::s][Cart::xx]+amb2*R[Cart::xxy][Cart::s][Cart::xx];
R[Cart::xyz][Cart::z][Cart::xx]=R[Cart::xyzz][Cart::s][Cart::xx]+amb2*R[Cart::xyz][Cart::s][Cart::xx];
R[Cart::yzz][Cart::z][Cart::xx]=R[Cart::yzzz][Cart::s][Cart::xx]+amb2*R[Cart::yzz][Cart::s][Cart::xx];
R[Cart::xxx][Cart::z][Cart::xx]=R[Cart::xxxz][Cart::s][Cart::xx]+amb2*R[Cart::xxx][Cart::s][Cart::xx];
R[Cart::xxz][Cart::z][Cart::xx]=R[Cart::xxzz][Cart::s][Cart::xx]+amb2*R[Cart::xxz][Cart::s][Cart::xx];
R[Cart::xzz][Cart::z][Cart::xx]=R[Cart::xzzz][Cart::s][Cart::xx]+amb2*R[Cart::xzz][Cart::s][Cart::xx];
R[Cart::zzz][Cart::z][Cart::xx]=R[Cart::zzzz][Cart::s][Cart::xx]+amb2*R[Cart::zzz][Cart::s][Cart::xx];
R[Cart::yyy][Cart::y][Cart::xz]=R[Cart::yyyy][Cart::s][Cart::xz]+amb1*R[Cart::yyy][Cart::s][Cart::xz];
R[Cart::xyy][Cart::y][Cart::xz]=R[Cart::xyyy][Cart::s][Cart::xz]+amb1*R[Cart::xyy][Cart::s][Cart::xz];
R[Cart::yyz][Cart::y][Cart::xz]=R[Cart::yyyz][Cart::s][Cart::xz]+amb1*R[Cart::yyz][Cart::s][Cart::xz];
R[Cart::xxy][Cart::y][Cart::xz]=R[Cart::xxyy][Cart::s][Cart::xz]+amb1*R[Cart::xxy][Cart::s][Cart::xz];
R[Cart::xyz][Cart::y][Cart::xz]=R[Cart::xyyz][Cart::s][Cart::xz]+amb1*R[Cart::xyz][Cart::s][Cart::xz];
R[Cart::yzz][Cart::y][Cart::xz]=R[Cart::yyzz][Cart::s][Cart::xz]+amb1*R[Cart::yzz][Cart::s][Cart::xz];
R[Cart::xxx][Cart::y][Cart::xz]=R[Cart::xxxy][Cart::s][Cart::xz]+amb1*R[Cart::xxx][Cart::s][Cart::xz];
R[Cart::xxz][Cart::y][Cart::xz]=R[Cart::xxyz][Cart::s][Cart::xz]+amb1*R[Cart::xxz][Cart::s][Cart::xz];
R[Cart::xzz][Cart::y][Cart::xz]=R[Cart::xyzz][Cart::s][Cart::xz]+amb1*R[Cart::xzz][Cart::s][Cart::xz];
R[Cart::zzz][Cart::y][Cart::xz]=R[Cart::yzzz][Cart::s][Cart::xz]+amb1*R[Cart::zzz][Cart::s][Cart::xz];
R[Cart::yyy][Cart::x][Cart::xz]=R[Cart::xyyy][Cart::s][Cart::xz]+amb0*R[Cart::yyy][Cart::s][Cart::xz];
R[Cart::xyy][Cart::x][Cart::xz]=R[Cart::xxyy][Cart::s][Cart::xz]+amb0*R[Cart::xyy][Cart::s][Cart::xz];
R[Cart::yyz][Cart::x][Cart::xz]=R[Cart::xyyz][Cart::s][Cart::xz]+amb0*R[Cart::yyz][Cart::s][Cart::xz];
R[Cart::xxy][Cart::x][Cart::xz]=R[Cart::xxxy][Cart::s][Cart::xz]+amb0*R[Cart::xxy][Cart::s][Cart::xz];
R[Cart::xyz][Cart::x][Cart::xz]=R[Cart::xxyz][Cart::s][Cart::xz]+amb0*R[Cart::xyz][Cart::s][Cart::xz];
R[Cart::yzz][Cart::x][Cart::xz]=R[Cart::xyzz][Cart::s][Cart::xz]+amb0*R[Cart::yzz][Cart::s][Cart::xz];
R[Cart::xxx][Cart::x][Cart::xz]=R[Cart::xxxx][Cart::s][Cart::xz]+amb0*R[Cart::xxx][Cart::s][Cart::xz];
R[Cart::xxz][Cart::x][Cart::xz]=R[Cart::xxxz][Cart::s][Cart::xz]+amb0*R[Cart::xxz][Cart::s][Cart::xz];
R[Cart::xzz][Cart::x][Cart::xz]=R[Cart::xxzz][Cart::s][Cart::xz]+amb0*R[Cart::xzz][Cart::s][Cart::xz];
R[Cart::zzz][Cart::x][Cart::xz]=R[Cart::xzzz][Cart::s][Cart::xz]+amb0*R[Cart::zzz][Cart::s][Cart::xz];
R[Cart::yyy][Cart::z][Cart::xz]=R[Cart::yyyz][Cart::s][Cart::xz]+amb2*R[Cart::yyy][Cart::s][Cart::xz];
R[Cart::xyy][Cart::z][Cart::xz]=R[Cart::xyyz][Cart::s][Cart::xz]+amb2*R[Cart::xyy][Cart::s][Cart::xz];
R[Cart::yyz][Cart::z][Cart::xz]=R[Cart::yyzz][Cart::s][Cart::xz]+amb2*R[Cart::yyz][Cart::s][Cart::xz];
R[Cart::xxy][Cart::z][Cart::xz]=R[Cart::xxyz][Cart::s][Cart::xz]+amb2*R[Cart::xxy][Cart::s][Cart::xz];
R[Cart::xyz][Cart::z][Cart::xz]=R[Cart::xyzz][Cart::s][Cart::xz]+amb2*R[Cart::xyz][Cart::s][Cart::xz];
R[Cart::yzz][Cart::z][Cart::xz]=R[Cart::yzzz][Cart::s][Cart::xz]+amb2*R[Cart::yzz][Cart::s][Cart::xz];
R[Cart::xxx][Cart::z][Cart::xz]=R[Cart::xxxz][Cart::s][Cart::xz]+amb2*R[Cart::xxx][Cart::s][Cart::xz];
R[Cart::xxz][Cart::z][Cart::xz]=R[Cart::xxzz][Cart::s][Cart::xz]+amb2*R[Cart::xxz][Cart::s][Cart::xz];
R[Cart::xzz][Cart::z][Cart::xz]=R[Cart::xzzz][Cart::s][Cart::xz]+amb2*R[Cart::xzz][Cart::s][Cart::xz];
R[Cart::zzz][Cart::z][Cart::xz]=R[Cart::zzzz][Cart::s][Cart::xz]+amb2*R[Cart::zzz][Cart::s][Cart::xz];
R[Cart::yyy][Cart::y][Cart::zz]=R[Cart::yyyy][Cart::s][Cart::zz]+amb1*R[Cart::yyy][Cart::s][Cart::zz];
R[Cart::xyy][Cart::y][Cart::zz]=R[Cart::xyyy][Cart::s][Cart::zz]+amb1*R[Cart::xyy][Cart::s][Cart::zz];
R[Cart::yyz][Cart::y][Cart::zz]=R[Cart::yyyz][Cart::s][Cart::zz]+amb1*R[Cart::yyz][Cart::s][Cart::zz];
R[Cart::xxy][Cart::y][Cart::zz]=R[Cart::xxyy][Cart::s][Cart::zz]+amb1*R[Cart::xxy][Cart::s][Cart::zz];
R[Cart::xyz][Cart::y][Cart::zz]=R[Cart::xyyz][Cart::s][Cart::zz]+amb1*R[Cart::xyz][Cart::s][Cart::zz];
R[Cart::yzz][Cart::y][Cart::zz]=R[Cart::yyzz][Cart::s][Cart::zz]+amb1*R[Cart::yzz][Cart::s][Cart::zz];
R[Cart::xxx][Cart::y][Cart::zz]=R[Cart::xxxy][Cart::s][Cart::zz]+amb1*R[Cart::xxx][Cart::s][Cart::zz];
R[Cart::xxz][Cart::y][Cart::zz]=R[Cart::xxyz][Cart::s][Cart::zz]+amb1*R[Cart::xxz][Cart::s][Cart::zz];
R[Cart::xzz][Cart::y][Cart::zz]=R[Cart::xyzz][Cart::s][Cart::zz]+amb1*R[Cart::xzz][Cart::s][Cart::zz];
R[Cart::zzz][Cart::y][Cart::zz]=R[Cart::yzzz][Cart::s][Cart::zz]+amb1*R[Cart::zzz][Cart::s][Cart::zz];
R[Cart::yyy][Cart::x][Cart::zz]=R[Cart::xyyy][Cart::s][Cart::zz]+amb0*R[Cart::yyy][Cart::s][Cart::zz];
R[Cart::xyy][Cart::x][Cart::zz]=R[Cart::xxyy][Cart::s][Cart::zz]+amb0*R[Cart::xyy][Cart::s][Cart::zz];
R[Cart::yyz][Cart::x][Cart::zz]=R[Cart::xyyz][Cart::s][Cart::zz]+amb0*R[Cart::yyz][Cart::s][Cart::zz];
R[Cart::xxy][Cart::x][Cart::zz]=R[Cart::xxxy][Cart::s][Cart::zz]+amb0*R[Cart::xxy][Cart::s][Cart::zz];
R[Cart::xyz][Cart::x][Cart::zz]=R[Cart::xxyz][Cart::s][Cart::zz]+amb0*R[Cart::xyz][Cart::s][Cart::zz];
R[Cart::yzz][Cart::x][Cart::zz]=R[Cart::xyzz][Cart::s][Cart::zz]+amb0*R[Cart::yzz][Cart::s][Cart::zz];
R[Cart::xxx][Cart::x][Cart::zz]=R[Cart::xxxx][Cart::s][Cart::zz]+amb0*R[Cart::xxx][Cart::s][Cart::zz];
R[Cart::xxz][Cart::x][Cart::zz]=R[Cart::xxxz][Cart::s][Cart::zz]+amb0*R[Cart::xxz][Cart::s][Cart::zz];
R[Cart::xzz][Cart::x][Cart::zz]=R[Cart::xxzz][Cart::s][Cart::zz]+amb0*R[Cart::xzz][Cart::s][Cart::zz];
R[Cart::zzz][Cart::x][Cart::zz]=R[Cart::xzzz][Cart::s][Cart::zz]+amb0*R[Cart::zzz][Cart::s][Cart::zz];
R[Cart::yyy][Cart::z][Cart::zz]=R[Cart::yyyz][Cart::s][Cart::zz]+amb2*R[Cart::yyy][Cart::s][Cart::zz];
R[Cart::xyy][Cart::z][Cart::zz]=R[Cart::xyyz][Cart::s][Cart::zz]+amb2*R[Cart::xyy][Cart::s][Cart::zz];
R[Cart::yyz][Cart::z][Cart::zz]=R[Cart::yyzz][Cart::s][Cart::zz]+amb2*R[Cart::yyz][Cart::s][Cart::zz];
R[Cart::xxy][Cart::z][Cart::zz]=R[Cart::xxyz][Cart::s][Cart::zz]+amb2*R[Cart::xxy][Cart::s][Cart::zz];
R[Cart::xyz][Cart::z][Cart::zz]=R[Cart::xyzz][Cart::s][Cart::zz]+amb2*R[Cart::xyz][Cart::s][Cart::zz];
R[Cart::yzz][Cart::z][Cart::zz]=R[Cart::yzzz][Cart::s][Cart::zz]+amb2*R[Cart::yzz][Cart::s][Cart::zz];
R[Cart::xxx][Cart::z][Cart::zz]=R[Cart::xxxz][Cart::s][Cart::zz]+amb2*R[Cart::xxx][Cart::s][Cart::zz];
R[Cart::xxz][Cart::z][Cart::zz]=R[Cart::xxzz][Cart::s][Cart::zz]+amb2*R[Cart::xxz][Cart::s][Cart::zz];
R[Cart::xzz][Cart::z][Cart::zz]=R[Cart::xzzz][Cart::s][Cart::zz]+amb2*R[Cart::xzz][Cart::s][Cart::zz];
R[Cart::zzz][Cart::z][Cart::zz]=R[Cart::zzzz][Cart::s][Cart::zz]+amb2*R[Cart::zzz][Cart::s][Cart::zz];
}
//------------------------------------------------------

//Integral s - p - f
if (_lmax_beta>0 && _lmax_gamma>2){

R[Cart::s][Cart::y][Cart::yyy]=R[Cart::y][Cart::s][Cart::yyy]+amb1*R[Cart::s][Cart::s][Cart::yyy];
R[Cart::s][Cart::x][Cart::yyy]=R[Cart::x][Cart::s][Cart::yyy]+amb0*R[Cart::s][Cart::s][Cart::yyy];
R[Cart::s][Cart::z][Cart::yyy]=R[Cart::z][Cart::s][Cart::yyy]+amb2*R[Cart::s][Cart::s][Cart::yyy];
R[Cart::s][Cart::y][Cart::xyy]=R[Cart::y][Cart::s][Cart::xyy]+amb1*R[Cart::s][Cart::s][Cart::xyy];
R[Cart::s][Cart::x][Cart::xyy]=R[Cart::x][Cart::s][Cart::xyy]+amb0*R[Cart::s][Cart::s][Cart::xyy];
R[Cart::s][Cart::z][Cart::xyy]=R[Cart::z][Cart::s][Cart::xyy]+amb2*R[Cart::s][Cart::s][Cart::xyy];
R[Cart::s][Cart::y][Cart::yyz]=R[Cart::y][Cart::s][Cart::yyz]+amb1*R[Cart::s][Cart::s][Cart::yyz];
R[Cart::s][Cart::x][Cart::yyz]=R[Cart::x][Cart::s][Cart::yyz]+amb0*R[Cart::s][Cart::s][Cart::yyz];
R[Cart::s][Cart::z][Cart::yyz]=R[Cart::z][Cart::s][Cart::yyz]+amb2*R[Cart::s][Cart::s][Cart::yyz];
R[Cart::s][Cart::y][Cart::xxy]=R[Cart::y][Cart::s][Cart::xxy]+amb1*R[Cart::s][Cart::s][Cart::xxy];
R[Cart::s][Cart::x][Cart::xxy]=R[Cart::x][Cart::s][Cart::xxy]+amb0*R[Cart::s][Cart::s][Cart::xxy];
R[Cart::s][Cart::z][Cart::xxy]=R[Cart::z][Cart::s][Cart::xxy]+amb2*R[Cart::s][Cart::s][Cart::xxy];
R[Cart::s][Cart::y][Cart::xyz]=R[Cart::y][Cart::s][Cart::xyz]+amb1*R[Cart::s][Cart::s][Cart::xyz];
R[Cart::s][Cart::x][Cart::xyz]=R[Cart::x][Cart::s][Cart::xyz]+amb0*R[Cart::s][Cart::s][Cart::xyz];
R[Cart::s][Cart::z][Cart::xyz]=R[Cart::z][Cart::s][Cart::xyz]+amb2*R[Cart::s][Cart::s][Cart::xyz];
R[Cart::s][Cart::y][Cart::yzz]=R[Cart::y][Cart::s][Cart::yzz]+amb1*R[Cart::s][Cart::s][Cart::yzz];
R[Cart::s][Cart::x][Cart::yzz]=R[Cart::x][Cart::s][Cart::yzz]+amb0*R[Cart::s][Cart::s][Cart::yzz];
R[Cart::s][Cart::z][Cart::yzz]=R[Cart::z][Cart::s][Cart::yzz]+amb2*R[Cart::s][Cart::s][Cart::yzz];
R[Cart::s][Cart::y][Cart::xxx]=R[Cart::y][Cart::s][Cart::xxx]+amb1*R[Cart::s][Cart::s][Cart::xxx];
R[Cart::s][Cart::x][Cart::xxx]=R[Cart::x][Cart::s][Cart::xxx]+amb0*R[Cart::s][Cart::s][Cart::xxx];
R[Cart::s][Cart::z][Cart::xxx]=R[Cart::z][Cart::s][Cart::xxx]+amb2*R[Cart::s][Cart::s][Cart::xxx];
R[Cart::s][Cart::y][Cart::xxz]=R[Cart::y][Cart::s][Cart::xxz]+amb1*R[Cart::s][Cart::s][Cart::xxz];
R[Cart::s][Cart::x][Cart::xxz]=R[Cart::x][Cart::s][Cart::xxz]+amb0*R[Cart::s][Cart::s][Cart::xxz];
R[Cart::s][Cart::z][Cart::xxz]=R[Cart::z][Cart::s][Cart::xxz]+amb2*R[Cart::s][Cart::s][Cart::xxz];
R[Cart::s][Cart::y][Cart::xzz]=R[Cart::y][Cart::s][Cart::xzz]+amb1*R[Cart::s][Cart::s][Cart::xzz];
R[Cart::s][Cart::x][Cart::xzz]=R[Cart::x][Cart::s][Cart::xzz]+amb0*R[Cart::s][Cart::s][Cart::xzz];
R[Cart::s][Cart::z][Cart::xzz]=R[Cart::z][Cart::s][Cart::xzz]+amb2*R[Cart::s][Cart::s][Cart::xzz];
R[Cart::s][Cart::y][Cart::zzz]=R[Cart::y][Cart::s][Cart::zzz]+amb1*R[Cart::s][Cart::s][Cart::zzz];
R[Cart::s][Cart::x][Cart::zzz]=R[Cart::x][Cart::s][Cart::zzz]+amb0*R[Cart::s][Cart::s][Cart::zzz];
R[Cart::s][Cart::z][Cart::zzz]=R[Cart::z][Cart::s][Cart::zzz]+amb2*R[Cart::s][Cart::s][Cart::zzz];
}
//------------------------------------------------------

//Integral p - p - f
if (_lmax_beta>0 && _lmax_alpha>0 && _lmax_gamma>2){

R[Cart::y][Cart::y][Cart::yyy]=R[Cart::yy][Cart::s][Cart::yyy]+amb1*R[Cart::y][Cart::s][Cart::yyy];
R[Cart::x][Cart::y][Cart::yyy]=R[Cart::xy][Cart::s][Cart::yyy]+amb1*R[Cart::x][Cart::s][Cart::yyy];
R[Cart::z][Cart::y][Cart::yyy]=R[Cart::yz][Cart::s][Cart::yyy]+amb1*R[Cart::z][Cart::s][Cart::yyy];
R[Cart::y][Cart::x][Cart::yyy]=R[Cart::xy][Cart::s][Cart::yyy]+amb0*R[Cart::y][Cart::s][Cart::yyy];
R[Cart::x][Cart::x][Cart::yyy]=R[Cart::xx][Cart::s][Cart::yyy]+amb0*R[Cart::x][Cart::s][Cart::yyy];
R[Cart::z][Cart::x][Cart::yyy]=R[Cart::xz][Cart::s][Cart::yyy]+amb0*R[Cart::z][Cart::s][Cart::yyy];
R[Cart::y][Cart::z][Cart::yyy]=R[Cart::yz][Cart::s][Cart::yyy]+amb2*R[Cart::y][Cart::s][Cart::yyy];
R[Cart::x][Cart::z][Cart::yyy]=R[Cart::xz][Cart::s][Cart::yyy]+amb2*R[Cart::x][Cart::s][Cart::yyy];
R[Cart::z][Cart::z][Cart::yyy]=R[Cart::zz][Cart::s][Cart::yyy]+amb2*R[Cart::z][Cart::s][Cart::yyy];
R[Cart::y][Cart::y][Cart::xyy]=R[Cart::yy][Cart::s][Cart::xyy]+amb1*R[Cart::y][Cart::s][Cart::xyy];
R[Cart::x][Cart::y][Cart::xyy]=R[Cart::xy][Cart::s][Cart::xyy]+amb1*R[Cart::x][Cart::s][Cart::xyy];
R[Cart::z][Cart::y][Cart::xyy]=R[Cart::yz][Cart::s][Cart::xyy]+amb1*R[Cart::z][Cart::s][Cart::xyy];
R[Cart::y][Cart::x][Cart::xyy]=R[Cart::xy][Cart::s][Cart::xyy]+amb0*R[Cart::y][Cart::s][Cart::xyy];
R[Cart::x][Cart::x][Cart::xyy]=R[Cart::xx][Cart::s][Cart::xyy]+amb0*R[Cart::x][Cart::s][Cart::xyy];
R[Cart::z][Cart::x][Cart::xyy]=R[Cart::xz][Cart::s][Cart::xyy]+amb0*R[Cart::z][Cart::s][Cart::xyy];
R[Cart::y][Cart::z][Cart::xyy]=R[Cart::yz][Cart::s][Cart::xyy]+amb2*R[Cart::y][Cart::s][Cart::xyy];
R[Cart::x][Cart::z][Cart::xyy]=R[Cart::xz][Cart::s][Cart::xyy]+amb2*R[Cart::x][Cart::s][Cart::xyy];
R[Cart::z][Cart::z][Cart::xyy]=R[Cart::zz][Cart::s][Cart::xyy]+amb2*R[Cart::z][Cart::s][Cart::xyy];
R[Cart::y][Cart::y][Cart::yyz]=R[Cart::yy][Cart::s][Cart::yyz]+amb1*R[Cart::y][Cart::s][Cart::yyz];
R[Cart::x][Cart::y][Cart::yyz]=R[Cart::xy][Cart::s][Cart::yyz]+amb1*R[Cart::x][Cart::s][Cart::yyz];
R[Cart::z][Cart::y][Cart::yyz]=R[Cart::yz][Cart::s][Cart::yyz]+amb1*R[Cart::z][Cart::s][Cart::yyz];
R[Cart::y][Cart::x][Cart::yyz]=R[Cart::xy][Cart::s][Cart::yyz]+amb0*R[Cart::y][Cart::s][Cart::yyz];
R[Cart::x][Cart::x][Cart::yyz]=R[Cart::xx][Cart::s][Cart::yyz]+amb0*R[Cart::x][Cart::s][Cart::yyz];
R[Cart::z][Cart::x][Cart::yyz]=R[Cart::xz][Cart::s][Cart::yyz]+amb0*R[Cart::z][Cart::s][Cart::yyz];
R[Cart::y][Cart::z][Cart::yyz]=R[Cart::yz][Cart::s][Cart::yyz]+amb2*R[Cart::y][Cart::s][Cart::yyz];
R[Cart::x][Cart::z][Cart::yyz]=R[Cart::xz][Cart::s][Cart::yyz]+amb2*R[Cart::x][Cart::s][Cart::yyz];
R[Cart::z][Cart::z][Cart::yyz]=R[Cart::zz][Cart::s][Cart::yyz]+amb2*R[Cart::z][Cart::s][Cart::yyz];
R[Cart::y][Cart::y][Cart::xxy]=R[Cart::yy][Cart::s][Cart::xxy]+amb1*R[Cart::y][Cart::s][Cart::xxy];
R[Cart::x][Cart::y][Cart::xxy]=R[Cart::xy][Cart::s][Cart::xxy]+amb1*R[Cart::x][Cart::s][Cart::xxy];
R[Cart::z][Cart::y][Cart::xxy]=R[Cart::yz][Cart::s][Cart::xxy]+amb1*R[Cart::z][Cart::s][Cart::xxy];
R[Cart::y][Cart::x][Cart::xxy]=R[Cart::xy][Cart::s][Cart::xxy]+amb0*R[Cart::y][Cart::s][Cart::xxy];
R[Cart::x][Cart::x][Cart::xxy]=R[Cart::xx][Cart::s][Cart::xxy]+amb0*R[Cart::x][Cart::s][Cart::xxy];
R[Cart::z][Cart::x][Cart::xxy]=R[Cart::xz][Cart::s][Cart::xxy]+amb0*R[Cart::z][Cart::s][Cart::xxy];
R[Cart::y][Cart::z][Cart::xxy]=R[Cart::yz][Cart::s][Cart::xxy]+amb2*R[Cart::y][Cart::s][Cart::xxy];
R[Cart::x][Cart::z][Cart::xxy]=R[Cart::xz][Cart::s][Cart::xxy]+amb2*R[Cart::x][Cart::s][Cart::xxy];
R[Cart::z][Cart::z][Cart::xxy]=R[Cart::zz][Cart::s][Cart::xxy]+amb2*R[Cart::z][Cart::s][Cart::xxy];
R[Cart::y][Cart::y][Cart::xyz]=R[Cart::yy][Cart::s][Cart::xyz]+amb1*R[Cart::y][Cart::s][Cart::xyz];
R[Cart::x][Cart::y][Cart::xyz]=R[Cart::xy][Cart::s][Cart::xyz]+amb1*R[Cart::x][Cart::s][Cart::xyz];
R[Cart::z][Cart::y][Cart::xyz]=R[Cart::yz][Cart::s][Cart::xyz]+amb1*R[Cart::z][Cart::s][Cart::xyz];
R[Cart::y][Cart::x][Cart::xyz]=R[Cart::xy][Cart::s][Cart::xyz]+amb0*R[Cart::y][Cart::s][Cart::xyz];
R[Cart::x][Cart::x][Cart::xyz]=R[Cart::xx][Cart::s][Cart::xyz]+amb0*R[Cart::x][Cart::s][Cart::xyz];
R[Cart::z][Cart::x][Cart::xyz]=R[Cart::xz][Cart::s][Cart::xyz]+amb0*R[Cart::z][Cart::s][Cart::xyz];
R[Cart::y][Cart::z][Cart::xyz]=R[Cart::yz][Cart::s][Cart::xyz]+amb2*R[Cart::y][Cart::s][Cart::xyz];
R[Cart::x][Cart::z][Cart::xyz]=R[Cart::xz][Cart::s][Cart::xyz]+amb2*R[Cart::x][Cart::s][Cart::xyz];
R[Cart::z][Cart::z][Cart::xyz]=R[Cart::zz][Cart::s][Cart::xyz]+amb2*R[Cart::z][Cart::s][Cart::xyz];
R[Cart::y][Cart::y][Cart::yzz]=R[Cart::yy][Cart::s][Cart::yzz]+amb1*R[Cart::y][Cart::s][Cart::yzz];
R[Cart::x][Cart::y][Cart::yzz]=R[Cart::xy][Cart::s][Cart::yzz]+amb1*R[Cart::x][Cart::s][Cart::yzz];
R[Cart::z][Cart::y][Cart::yzz]=R[Cart::yz][Cart::s][Cart::yzz]+amb1*R[Cart::z][Cart::s][Cart::yzz];
R[Cart::y][Cart::x][Cart::yzz]=R[Cart::xy][Cart::s][Cart::yzz]+amb0*R[Cart::y][Cart::s][Cart::yzz];
R[Cart::x][Cart::x][Cart::yzz]=R[Cart::xx][Cart::s][Cart::yzz]+amb0*R[Cart::x][Cart::s][Cart::yzz];
R[Cart::z][Cart::x][Cart::yzz]=R[Cart::xz][Cart::s][Cart::yzz]+amb0*R[Cart::z][Cart::s][Cart::yzz];
R[Cart::y][Cart::z][Cart::yzz]=R[Cart::yz][Cart::s][Cart::yzz]+amb2*R[Cart::y][Cart::s][Cart::yzz];
R[Cart::x][Cart::z][Cart::yzz]=R[Cart::xz][Cart::s][Cart::yzz]+amb2*R[Cart::x][Cart::s][Cart::yzz];
R[Cart::z][Cart::z][Cart::yzz]=R[Cart::zz][Cart::s][Cart::yzz]+amb2*R[Cart::z][Cart::s][Cart::yzz];
R[Cart::y][Cart::y][Cart::xxx]=R[Cart::yy][Cart::s][Cart::xxx]+amb1*R[Cart::y][Cart::s][Cart::xxx];
R[Cart::x][Cart::y][Cart::xxx]=R[Cart::xy][Cart::s][Cart::xxx]+amb1*R[Cart::x][Cart::s][Cart::xxx];
R[Cart::z][Cart::y][Cart::xxx]=R[Cart::yz][Cart::s][Cart::xxx]+amb1*R[Cart::z][Cart::s][Cart::xxx];
R[Cart::y][Cart::x][Cart::xxx]=R[Cart::xy][Cart::s][Cart::xxx]+amb0*R[Cart::y][Cart::s][Cart::xxx];
R[Cart::x][Cart::x][Cart::xxx]=R[Cart::xx][Cart::s][Cart::xxx]+amb0*R[Cart::x][Cart::s][Cart::xxx];
R[Cart::z][Cart::x][Cart::xxx]=R[Cart::xz][Cart::s][Cart::xxx]+amb0*R[Cart::z][Cart::s][Cart::xxx];
R[Cart::y][Cart::z][Cart::xxx]=R[Cart::yz][Cart::s][Cart::xxx]+amb2*R[Cart::y][Cart::s][Cart::xxx];
R[Cart::x][Cart::z][Cart::xxx]=R[Cart::xz][Cart::s][Cart::xxx]+amb2*R[Cart::x][Cart::s][Cart::xxx];
R[Cart::z][Cart::z][Cart::xxx]=R[Cart::zz][Cart::s][Cart::xxx]+amb2*R[Cart::z][Cart::s][Cart::xxx];
R[Cart::y][Cart::y][Cart::xxz]=R[Cart::yy][Cart::s][Cart::xxz]+amb1*R[Cart::y][Cart::s][Cart::xxz];
R[Cart::x][Cart::y][Cart::xxz]=R[Cart::xy][Cart::s][Cart::xxz]+amb1*R[Cart::x][Cart::s][Cart::xxz];
R[Cart::z][Cart::y][Cart::xxz]=R[Cart::yz][Cart::s][Cart::xxz]+amb1*R[Cart::z][Cart::s][Cart::xxz];
R[Cart::y][Cart::x][Cart::xxz]=R[Cart::xy][Cart::s][Cart::xxz]+amb0*R[Cart::y][Cart::s][Cart::xxz];
R[Cart::x][Cart::x][Cart::xxz]=R[Cart::xx][Cart::s][Cart::xxz]+amb0*R[Cart::x][Cart::s][Cart::xxz];
R[Cart::z][Cart::x][Cart::xxz]=R[Cart::xz][Cart::s][Cart::xxz]+amb0*R[Cart::z][Cart::s][Cart::xxz];
R[Cart::y][Cart::z][Cart::xxz]=R[Cart::yz][Cart::s][Cart::xxz]+amb2*R[Cart::y][Cart::s][Cart::xxz];
R[Cart::x][Cart::z][Cart::xxz]=R[Cart::xz][Cart::s][Cart::xxz]+amb2*R[Cart::x][Cart::s][Cart::xxz];
R[Cart::z][Cart::z][Cart::xxz]=R[Cart::zz][Cart::s][Cart::xxz]+amb2*R[Cart::z][Cart::s][Cart::xxz];
R[Cart::y][Cart::y][Cart::xzz]=R[Cart::yy][Cart::s][Cart::xzz]+amb1*R[Cart::y][Cart::s][Cart::xzz];
R[Cart::x][Cart::y][Cart::xzz]=R[Cart::xy][Cart::s][Cart::xzz]+amb1*R[Cart::x][Cart::s][Cart::xzz];
R[Cart::z][Cart::y][Cart::xzz]=R[Cart::yz][Cart::s][Cart::xzz]+amb1*R[Cart::z][Cart::s][Cart::xzz];
R[Cart::y][Cart::x][Cart::xzz]=R[Cart::xy][Cart::s][Cart::xzz]+amb0*R[Cart::y][Cart::s][Cart::xzz];
R[Cart::x][Cart::x][Cart::xzz]=R[Cart::xx][Cart::s][Cart::xzz]+amb0*R[Cart::x][Cart::s][Cart::xzz];
R[Cart::z][Cart::x][Cart::xzz]=R[Cart::xz][Cart::s][Cart::xzz]+amb0*R[Cart::z][Cart::s][Cart::xzz];
R[Cart::y][Cart::z][Cart::xzz]=R[Cart::yz][Cart::s][Cart::xzz]+amb2*R[Cart::y][Cart::s][Cart::xzz];
R[Cart::x][Cart::z][Cart::xzz]=R[Cart::xz][Cart::s][Cart::xzz]+amb2*R[Cart::x][Cart::s][Cart::xzz];
R[Cart::z][Cart::z][Cart::xzz]=R[Cart::zz][Cart::s][Cart::xzz]+amb2*R[Cart::z][Cart::s][Cart::xzz];
R[Cart::y][Cart::y][Cart::zzz]=R[Cart::yy][Cart::s][Cart::zzz]+amb1*R[Cart::y][Cart::s][Cart::zzz];
R[Cart::x][Cart::y][Cart::zzz]=R[Cart::xy][Cart::s][Cart::zzz]+amb1*R[Cart::x][Cart::s][Cart::zzz];
R[Cart::z][Cart::y][Cart::zzz]=R[Cart::yz][Cart::s][Cart::zzz]+amb1*R[Cart::z][Cart::s][Cart::zzz];
R[Cart::y][Cart::x][Cart::zzz]=R[Cart::xy][Cart::s][Cart::zzz]+amb0*R[Cart::y][Cart::s][Cart::zzz];
R[Cart::x][Cart::x][Cart::zzz]=R[Cart::xx][Cart::s][Cart::zzz]+amb0*R[Cart::x][Cart::s][Cart::zzz];
R[Cart::z][Cart::x][Cart::zzz]=R[Cart::xz][Cart::s][Cart::zzz]+amb0*R[Cart::z][Cart::s][Cart::zzz];
R[Cart::y][Cart::z][Cart::zzz]=R[Cart::yz][Cart::s][Cart::zzz]+amb2*R[Cart::y][Cart::s][Cart::zzz];
R[Cart::x][Cart::z][Cart::zzz]=R[Cart::xz][Cart::s][Cart::zzz]+amb2*R[Cart::x][Cart::s][Cart::zzz];
R[Cart::z][Cart::z][Cart::zzz]=R[Cart::zz][Cart::s][Cart::zzz]+amb2*R[Cart::z][Cart::s][Cart::zzz];
}
//------------------------------------------------------

//Integral d - p - f
if (_lmax_beta>0 && _lmax_alpha>1 && _lmax_gamma>2){

R[Cart::yy][Cart::y][Cart::yyy]=R[Cart::yyy][Cart::s][Cart::yyy]+amb1*R[Cart::yy][Cart::s][Cart::yyy];
R[Cart::xy][Cart::y][Cart::yyy]=R[Cart::xyy][Cart::s][Cart::yyy]+amb1*R[Cart::xy][Cart::s][Cart::yyy];
R[Cart::yz][Cart::y][Cart::yyy]=R[Cart::yyz][Cart::s][Cart::yyy]+amb1*R[Cart::yz][Cart::s][Cart::yyy];
R[Cart::xx][Cart::y][Cart::yyy]=R[Cart::xxy][Cart::s][Cart::yyy]+amb1*R[Cart::xx][Cart::s][Cart::yyy];
R[Cart::xz][Cart::y][Cart::yyy]=R[Cart::xyz][Cart::s][Cart::yyy]+amb1*R[Cart::xz][Cart::s][Cart::yyy];
R[Cart::zz][Cart::y][Cart::yyy]=R[Cart::yzz][Cart::s][Cart::yyy]+amb1*R[Cart::zz][Cart::s][Cart::yyy];
R[Cart::yy][Cart::x][Cart::yyy]=R[Cart::xyy][Cart::s][Cart::yyy]+amb0*R[Cart::yy][Cart::s][Cart::yyy];
R[Cart::xy][Cart::x][Cart::yyy]=R[Cart::xxy][Cart::s][Cart::yyy]+amb0*R[Cart::xy][Cart::s][Cart::yyy];
R[Cart::yz][Cart::x][Cart::yyy]=R[Cart::xyz][Cart::s][Cart::yyy]+amb0*R[Cart::yz][Cart::s][Cart::yyy];
R[Cart::xx][Cart::x][Cart::yyy]=R[Cart::xxx][Cart::s][Cart::yyy]+amb0*R[Cart::xx][Cart::s][Cart::yyy];
R[Cart::xz][Cart::x][Cart::yyy]=R[Cart::xxz][Cart::s][Cart::yyy]+amb0*R[Cart::xz][Cart::s][Cart::yyy];
R[Cart::zz][Cart::x][Cart::yyy]=R[Cart::xzz][Cart::s][Cart::yyy]+amb0*R[Cart::zz][Cart::s][Cart::yyy];
R[Cart::yy][Cart::z][Cart::yyy]=R[Cart::yyz][Cart::s][Cart::yyy]+amb2*R[Cart::yy][Cart::s][Cart::yyy];
R[Cart::xy][Cart::z][Cart::yyy]=R[Cart::xyz][Cart::s][Cart::yyy]+amb2*R[Cart::xy][Cart::s][Cart::yyy];
R[Cart::yz][Cart::z][Cart::yyy]=R[Cart::yzz][Cart::s][Cart::yyy]+amb2*R[Cart::yz][Cart::s][Cart::yyy];
R[Cart::xx][Cart::z][Cart::yyy]=R[Cart::xxz][Cart::s][Cart::yyy]+amb2*R[Cart::xx][Cart::s][Cart::yyy];
R[Cart::xz][Cart::z][Cart::yyy]=R[Cart::xzz][Cart::s][Cart::yyy]+amb2*R[Cart::xz][Cart::s][Cart::yyy];
R[Cart::zz][Cart::z][Cart::yyy]=R[Cart::zzz][Cart::s][Cart::yyy]+amb2*R[Cart::zz][Cart::s][Cart::yyy];
R[Cart::yy][Cart::y][Cart::xyy]=R[Cart::yyy][Cart::s][Cart::xyy]+amb1*R[Cart::yy][Cart::s][Cart::xyy];
R[Cart::xy][Cart::y][Cart::xyy]=R[Cart::xyy][Cart::s][Cart::xyy]+amb1*R[Cart::xy][Cart::s][Cart::xyy];
R[Cart::yz][Cart::y][Cart::xyy]=R[Cart::yyz][Cart::s][Cart::xyy]+amb1*R[Cart::yz][Cart::s][Cart::xyy];
R[Cart::xx][Cart::y][Cart::xyy]=R[Cart::xxy][Cart::s][Cart::xyy]+amb1*R[Cart::xx][Cart::s][Cart::xyy];
R[Cart::xz][Cart::y][Cart::xyy]=R[Cart::xyz][Cart::s][Cart::xyy]+amb1*R[Cart::xz][Cart::s][Cart::xyy];
R[Cart::zz][Cart::y][Cart::xyy]=R[Cart::yzz][Cart::s][Cart::xyy]+amb1*R[Cart::zz][Cart::s][Cart::xyy];
R[Cart::yy][Cart::x][Cart::xyy]=R[Cart::xyy][Cart::s][Cart::xyy]+amb0*R[Cart::yy][Cart::s][Cart::xyy];
R[Cart::xy][Cart::x][Cart::xyy]=R[Cart::xxy][Cart::s][Cart::xyy]+amb0*R[Cart::xy][Cart::s][Cart::xyy];
R[Cart::yz][Cart::x][Cart::xyy]=R[Cart::xyz][Cart::s][Cart::xyy]+amb0*R[Cart::yz][Cart::s][Cart::xyy];
R[Cart::xx][Cart::x][Cart::xyy]=R[Cart::xxx][Cart::s][Cart::xyy]+amb0*R[Cart::xx][Cart::s][Cart::xyy];
R[Cart::xz][Cart::x][Cart::xyy]=R[Cart::xxz][Cart::s][Cart::xyy]+amb0*R[Cart::xz][Cart::s][Cart::xyy];
R[Cart::zz][Cart::x][Cart::xyy]=R[Cart::xzz][Cart::s][Cart::xyy]+amb0*R[Cart::zz][Cart::s][Cart::xyy];
R[Cart::yy][Cart::z][Cart::xyy]=R[Cart::yyz][Cart::s][Cart::xyy]+amb2*R[Cart::yy][Cart::s][Cart::xyy];
R[Cart::xy][Cart::z][Cart::xyy]=R[Cart::xyz][Cart::s][Cart::xyy]+amb2*R[Cart::xy][Cart::s][Cart::xyy];
R[Cart::yz][Cart::z][Cart::xyy]=R[Cart::yzz][Cart::s][Cart::xyy]+amb2*R[Cart::yz][Cart::s][Cart::xyy];
R[Cart::xx][Cart::z][Cart::xyy]=R[Cart::xxz][Cart::s][Cart::xyy]+amb2*R[Cart::xx][Cart::s][Cart::xyy];
R[Cart::xz][Cart::z][Cart::xyy]=R[Cart::xzz][Cart::s][Cart::xyy]+amb2*R[Cart::xz][Cart::s][Cart::xyy];
R[Cart::zz][Cart::z][Cart::xyy]=R[Cart::zzz][Cart::s][Cart::xyy]+amb2*R[Cart::zz][Cart::s][Cart::xyy];
R[Cart::yy][Cart::y][Cart::yyz]=R[Cart::yyy][Cart::s][Cart::yyz]+amb1*R[Cart::yy][Cart::s][Cart::yyz];
R[Cart::xy][Cart::y][Cart::yyz]=R[Cart::xyy][Cart::s][Cart::yyz]+amb1*R[Cart::xy][Cart::s][Cart::yyz];
R[Cart::yz][Cart::y][Cart::yyz]=R[Cart::yyz][Cart::s][Cart::yyz]+amb1*R[Cart::yz][Cart::s][Cart::yyz];
R[Cart::xx][Cart::y][Cart::yyz]=R[Cart::xxy][Cart::s][Cart::yyz]+amb1*R[Cart::xx][Cart::s][Cart::yyz];
R[Cart::xz][Cart::y][Cart::yyz]=R[Cart::xyz][Cart::s][Cart::yyz]+amb1*R[Cart::xz][Cart::s][Cart::yyz];
R[Cart::zz][Cart::y][Cart::yyz]=R[Cart::yzz][Cart::s][Cart::yyz]+amb1*R[Cart::zz][Cart::s][Cart::yyz];
R[Cart::yy][Cart::x][Cart::yyz]=R[Cart::xyy][Cart::s][Cart::yyz]+amb0*R[Cart::yy][Cart::s][Cart::yyz];
R[Cart::xy][Cart::x][Cart::yyz]=R[Cart::xxy][Cart::s][Cart::yyz]+amb0*R[Cart::xy][Cart::s][Cart::yyz];
R[Cart::yz][Cart::x][Cart::yyz]=R[Cart::xyz][Cart::s][Cart::yyz]+amb0*R[Cart::yz][Cart::s][Cart::yyz];
R[Cart::xx][Cart::x][Cart::yyz]=R[Cart::xxx][Cart::s][Cart::yyz]+amb0*R[Cart::xx][Cart::s][Cart::yyz];
R[Cart::xz][Cart::x][Cart::yyz]=R[Cart::xxz][Cart::s][Cart::yyz]+amb0*R[Cart::xz][Cart::s][Cart::yyz];
R[Cart::zz][Cart::x][Cart::yyz]=R[Cart::xzz][Cart::s][Cart::yyz]+amb0*R[Cart::zz][Cart::s][Cart::yyz];
R[Cart::yy][Cart::z][Cart::yyz]=R[Cart::yyz][Cart::s][Cart::yyz]+amb2*R[Cart::yy][Cart::s][Cart::yyz];
R[Cart::xy][Cart::z][Cart::yyz]=R[Cart::xyz][Cart::s][Cart::yyz]+amb2*R[Cart::xy][Cart::s][Cart::yyz];
R[Cart::yz][Cart::z][Cart::yyz]=R[Cart::yzz][Cart::s][Cart::yyz]+amb2*R[Cart::yz][Cart::s][Cart::yyz];
R[Cart::xx][Cart::z][Cart::yyz]=R[Cart::xxz][Cart::s][Cart::yyz]+amb2*R[Cart::xx][Cart::s][Cart::yyz];
R[Cart::xz][Cart::z][Cart::yyz]=R[Cart::xzz][Cart::s][Cart::yyz]+amb2*R[Cart::xz][Cart::s][Cart::yyz];
R[Cart::zz][Cart::z][Cart::yyz]=R[Cart::zzz][Cart::s][Cart::yyz]+amb2*R[Cart::zz][Cart::s][Cart::yyz];
R[Cart::yy][Cart::y][Cart::xxy]=R[Cart::yyy][Cart::s][Cart::xxy]+amb1*R[Cart::yy][Cart::s][Cart::xxy];
R[Cart::xy][Cart::y][Cart::xxy]=R[Cart::xyy][Cart::s][Cart::xxy]+amb1*R[Cart::xy][Cart::s][Cart::xxy];
R[Cart::yz][Cart::y][Cart::xxy]=R[Cart::yyz][Cart::s][Cart::xxy]+amb1*R[Cart::yz][Cart::s][Cart::xxy];
R[Cart::xx][Cart::y][Cart::xxy]=R[Cart::xxy][Cart::s][Cart::xxy]+amb1*R[Cart::xx][Cart::s][Cart::xxy];
R[Cart::xz][Cart::y][Cart::xxy]=R[Cart::xyz][Cart::s][Cart::xxy]+amb1*R[Cart::xz][Cart::s][Cart::xxy];
R[Cart::zz][Cart::y][Cart::xxy]=R[Cart::yzz][Cart::s][Cart::xxy]+amb1*R[Cart::zz][Cart::s][Cart::xxy];
R[Cart::yy][Cart::x][Cart::xxy]=R[Cart::xyy][Cart::s][Cart::xxy]+amb0*R[Cart::yy][Cart::s][Cart::xxy];
R[Cart::xy][Cart::x][Cart::xxy]=R[Cart::xxy][Cart::s][Cart::xxy]+amb0*R[Cart::xy][Cart::s][Cart::xxy];
R[Cart::yz][Cart::x][Cart::xxy]=R[Cart::xyz][Cart::s][Cart::xxy]+amb0*R[Cart::yz][Cart::s][Cart::xxy];
R[Cart::xx][Cart::x][Cart::xxy]=R[Cart::xxx][Cart::s][Cart::xxy]+amb0*R[Cart::xx][Cart::s][Cart::xxy];
R[Cart::xz][Cart::x][Cart::xxy]=R[Cart::xxz][Cart::s][Cart::xxy]+amb0*R[Cart::xz][Cart::s][Cart::xxy];
R[Cart::zz][Cart::x][Cart::xxy]=R[Cart::xzz][Cart::s][Cart::xxy]+amb0*R[Cart::zz][Cart::s][Cart::xxy];
R[Cart::yy][Cart::z][Cart::xxy]=R[Cart::yyz][Cart::s][Cart::xxy]+amb2*R[Cart::yy][Cart::s][Cart::xxy];
R[Cart::xy][Cart::z][Cart::xxy]=R[Cart::xyz][Cart::s][Cart::xxy]+amb2*R[Cart::xy][Cart::s][Cart::xxy];
R[Cart::yz][Cart::z][Cart::xxy]=R[Cart::yzz][Cart::s][Cart::xxy]+amb2*R[Cart::yz][Cart::s][Cart::xxy];
R[Cart::xx][Cart::z][Cart::xxy]=R[Cart::xxz][Cart::s][Cart::xxy]+amb2*R[Cart::xx][Cart::s][Cart::xxy];
R[Cart::xz][Cart::z][Cart::xxy]=R[Cart::xzz][Cart::s][Cart::xxy]+amb2*R[Cart::xz][Cart::s][Cart::xxy];
R[Cart::zz][Cart::z][Cart::xxy]=R[Cart::zzz][Cart::s][Cart::xxy]+amb2*R[Cart::zz][Cart::s][Cart::xxy];
R[Cart::yy][Cart::y][Cart::xyz]=R[Cart::yyy][Cart::s][Cart::xyz]+amb1*R[Cart::yy][Cart::s][Cart::xyz];
R[Cart::xy][Cart::y][Cart::xyz]=R[Cart::xyy][Cart::s][Cart::xyz]+amb1*R[Cart::xy][Cart::s][Cart::xyz];
R[Cart::yz][Cart::y][Cart::xyz]=R[Cart::yyz][Cart::s][Cart::xyz]+amb1*R[Cart::yz][Cart::s][Cart::xyz];
R[Cart::xx][Cart::y][Cart::xyz]=R[Cart::xxy][Cart::s][Cart::xyz]+amb1*R[Cart::xx][Cart::s][Cart::xyz];
R[Cart::xz][Cart::y][Cart::xyz]=R[Cart::xyz][Cart::s][Cart::xyz]+amb1*R[Cart::xz][Cart::s][Cart::xyz];
R[Cart::zz][Cart::y][Cart::xyz]=R[Cart::yzz][Cart::s][Cart::xyz]+amb1*R[Cart::zz][Cart::s][Cart::xyz];
R[Cart::yy][Cart::x][Cart::xyz]=R[Cart::xyy][Cart::s][Cart::xyz]+amb0*R[Cart::yy][Cart::s][Cart::xyz];
R[Cart::xy][Cart::x][Cart::xyz]=R[Cart::xxy][Cart::s][Cart::xyz]+amb0*R[Cart::xy][Cart::s][Cart::xyz];
R[Cart::yz][Cart::x][Cart::xyz]=R[Cart::xyz][Cart::s][Cart::xyz]+amb0*R[Cart::yz][Cart::s][Cart::xyz];
R[Cart::xx][Cart::x][Cart::xyz]=R[Cart::xxx][Cart::s][Cart::xyz]+amb0*R[Cart::xx][Cart::s][Cart::xyz];
R[Cart::xz][Cart::x][Cart::xyz]=R[Cart::xxz][Cart::s][Cart::xyz]+amb0*R[Cart::xz][Cart::s][Cart::xyz];
R[Cart::zz][Cart::x][Cart::xyz]=R[Cart::xzz][Cart::s][Cart::xyz]+amb0*R[Cart::zz][Cart::s][Cart::xyz];
R[Cart::yy][Cart::z][Cart::xyz]=R[Cart::yyz][Cart::s][Cart::xyz]+amb2*R[Cart::yy][Cart::s][Cart::xyz];
R[Cart::xy][Cart::z][Cart::xyz]=R[Cart::xyz][Cart::s][Cart::xyz]+amb2*R[Cart::xy][Cart::s][Cart::xyz];
R[Cart::yz][Cart::z][Cart::xyz]=R[Cart::yzz][Cart::s][Cart::xyz]+amb2*R[Cart::yz][Cart::s][Cart::xyz];
R[Cart::xx][Cart::z][Cart::xyz]=R[Cart::xxz][Cart::s][Cart::xyz]+amb2*R[Cart::xx][Cart::s][Cart::xyz];
R[Cart::xz][Cart::z][Cart::xyz]=R[Cart::xzz][Cart::s][Cart::xyz]+amb2*R[Cart::xz][Cart::s][Cart::xyz];
R[Cart::zz][Cart::z][Cart::xyz]=R[Cart::zzz][Cart::s][Cart::xyz]+amb2*R[Cart::zz][Cart::s][Cart::xyz];
R[Cart::yy][Cart::y][Cart::yzz]=R[Cart::yyy][Cart::s][Cart::yzz]+amb1*R[Cart::yy][Cart::s][Cart::yzz];
R[Cart::xy][Cart::y][Cart::yzz]=R[Cart::xyy][Cart::s][Cart::yzz]+amb1*R[Cart::xy][Cart::s][Cart::yzz];
R[Cart::yz][Cart::y][Cart::yzz]=R[Cart::yyz][Cart::s][Cart::yzz]+amb1*R[Cart::yz][Cart::s][Cart::yzz];
R[Cart::xx][Cart::y][Cart::yzz]=R[Cart::xxy][Cart::s][Cart::yzz]+amb1*R[Cart::xx][Cart::s][Cart::yzz];
R[Cart::xz][Cart::y][Cart::yzz]=R[Cart::xyz][Cart::s][Cart::yzz]+amb1*R[Cart::xz][Cart::s][Cart::yzz];
R[Cart::zz][Cart::y][Cart::yzz]=R[Cart::yzz][Cart::s][Cart::yzz]+amb1*R[Cart::zz][Cart::s][Cart::yzz];
R[Cart::yy][Cart::x][Cart::yzz]=R[Cart::xyy][Cart::s][Cart::yzz]+amb0*R[Cart::yy][Cart::s][Cart::yzz];
R[Cart::xy][Cart::x][Cart::yzz]=R[Cart::xxy][Cart::s][Cart::yzz]+amb0*R[Cart::xy][Cart::s][Cart::yzz];
R[Cart::yz][Cart::x][Cart::yzz]=R[Cart::xyz][Cart::s][Cart::yzz]+amb0*R[Cart::yz][Cart::s][Cart::yzz];
R[Cart::xx][Cart::x][Cart::yzz]=R[Cart::xxx][Cart::s][Cart::yzz]+amb0*R[Cart::xx][Cart::s][Cart::yzz];
R[Cart::xz][Cart::x][Cart::yzz]=R[Cart::xxz][Cart::s][Cart::yzz]+amb0*R[Cart::xz][Cart::s][Cart::yzz];
R[Cart::zz][Cart::x][Cart::yzz]=R[Cart::xzz][Cart::s][Cart::yzz]+amb0*R[Cart::zz][Cart::s][Cart::yzz];
R[Cart::yy][Cart::z][Cart::yzz]=R[Cart::yyz][Cart::s][Cart::yzz]+amb2*R[Cart::yy][Cart::s][Cart::yzz];
R[Cart::xy][Cart::z][Cart::yzz]=R[Cart::xyz][Cart::s][Cart::yzz]+amb2*R[Cart::xy][Cart::s][Cart::yzz];
R[Cart::yz][Cart::z][Cart::yzz]=R[Cart::yzz][Cart::s][Cart::yzz]+amb2*R[Cart::yz][Cart::s][Cart::yzz];
R[Cart::xx][Cart::z][Cart::yzz]=R[Cart::xxz][Cart::s][Cart::yzz]+amb2*R[Cart::xx][Cart::s][Cart::yzz];
R[Cart::xz][Cart::z][Cart::yzz]=R[Cart::xzz][Cart::s][Cart::yzz]+amb2*R[Cart::xz][Cart::s][Cart::yzz];
R[Cart::zz][Cart::z][Cart::yzz]=R[Cart::zzz][Cart::s][Cart::yzz]+amb2*R[Cart::zz][Cart::s][Cart::yzz];
R[Cart::yy][Cart::y][Cart::xxx]=R[Cart::yyy][Cart::s][Cart::xxx]+amb1*R[Cart::yy][Cart::s][Cart::xxx];
R[Cart::xy][Cart::y][Cart::xxx]=R[Cart::xyy][Cart::s][Cart::xxx]+amb1*R[Cart::xy][Cart::s][Cart::xxx];
R[Cart::yz][Cart::y][Cart::xxx]=R[Cart::yyz][Cart::s][Cart::xxx]+amb1*R[Cart::yz][Cart::s][Cart::xxx];
R[Cart::xx][Cart::y][Cart::xxx]=R[Cart::xxy][Cart::s][Cart::xxx]+amb1*R[Cart::xx][Cart::s][Cart::xxx];
R[Cart::xz][Cart::y][Cart::xxx]=R[Cart::xyz][Cart::s][Cart::xxx]+amb1*R[Cart::xz][Cart::s][Cart::xxx];
R[Cart::zz][Cart::y][Cart::xxx]=R[Cart::yzz][Cart::s][Cart::xxx]+amb1*R[Cart::zz][Cart::s][Cart::xxx];
R[Cart::yy][Cart::x][Cart::xxx]=R[Cart::xyy][Cart::s][Cart::xxx]+amb0*R[Cart::yy][Cart::s][Cart::xxx];
R[Cart::xy][Cart::x][Cart::xxx]=R[Cart::xxy][Cart::s][Cart::xxx]+amb0*R[Cart::xy][Cart::s][Cart::xxx];
R[Cart::yz][Cart::x][Cart::xxx]=R[Cart::xyz][Cart::s][Cart::xxx]+amb0*R[Cart::yz][Cart::s][Cart::xxx];
R[Cart::xx][Cart::x][Cart::xxx]=R[Cart::xxx][Cart::s][Cart::xxx]+amb0*R[Cart::xx][Cart::s][Cart::xxx];
R[Cart::xz][Cart::x][Cart::xxx]=R[Cart::xxz][Cart::s][Cart::xxx]+amb0*R[Cart::xz][Cart::s][Cart::xxx];
R[Cart::zz][Cart::x][Cart::xxx]=R[Cart::xzz][Cart::s][Cart::xxx]+amb0*R[Cart::zz][Cart::s][Cart::xxx];
R[Cart::yy][Cart::z][Cart::xxx]=R[Cart::yyz][Cart::s][Cart::xxx]+amb2*R[Cart::yy][Cart::s][Cart::xxx];
R[Cart::xy][Cart::z][Cart::xxx]=R[Cart::xyz][Cart::s][Cart::xxx]+amb2*R[Cart::xy][Cart::s][Cart::xxx];
R[Cart::yz][Cart::z][Cart::xxx]=R[Cart::yzz][Cart::s][Cart::xxx]+amb2*R[Cart::yz][Cart::s][Cart::xxx];
R[Cart::xx][Cart::z][Cart::xxx]=R[Cart::xxz][Cart::s][Cart::xxx]+amb2*R[Cart::xx][Cart::s][Cart::xxx];
R[Cart::xz][Cart::z][Cart::xxx]=R[Cart::xzz][Cart::s][Cart::xxx]+amb2*R[Cart::xz][Cart::s][Cart::xxx];
R[Cart::zz][Cart::z][Cart::xxx]=R[Cart::zzz][Cart::s][Cart::xxx]+amb2*R[Cart::zz][Cart::s][Cart::xxx];
R[Cart::yy][Cart::y][Cart::xxz]=R[Cart::yyy][Cart::s][Cart::xxz]+amb1*R[Cart::yy][Cart::s][Cart::xxz];
R[Cart::xy][Cart::y][Cart::xxz]=R[Cart::xyy][Cart::s][Cart::xxz]+amb1*R[Cart::xy][Cart::s][Cart::xxz];
R[Cart::yz][Cart::y][Cart::xxz]=R[Cart::yyz][Cart::s][Cart::xxz]+amb1*R[Cart::yz][Cart::s][Cart::xxz];
R[Cart::xx][Cart::y][Cart::xxz]=R[Cart::xxy][Cart::s][Cart::xxz]+amb1*R[Cart::xx][Cart::s][Cart::xxz];
R[Cart::xz][Cart::y][Cart::xxz]=R[Cart::xyz][Cart::s][Cart::xxz]+amb1*R[Cart::xz][Cart::s][Cart::xxz];
R[Cart::zz][Cart::y][Cart::xxz]=R[Cart::yzz][Cart::s][Cart::xxz]+amb1*R[Cart::zz][Cart::s][Cart::xxz];
R[Cart::yy][Cart::x][Cart::xxz]=R[Cart::xyy][Cart::s][Cart::xxz]+amb0*R[Cart::yy][Cart::s][Cart::xxz];
R[Cart::xy][Cart::x][Cart::xxz]=R[Cart::xxy][Cart::s][Cart::xxz]+amb0*R[Cart::xy][Cart::s][Cart::xxz];
R[Cart::yz][Cart::x][Cart::xxz]=R[Cart::xyz][Cart::s][Cart::xxz]+amb0*R[Cart::yz][Cart::s][Cart::xxz];
R[Cart::xx][Cart::x][Cart::xxz]=R[Cart::xxx][Cart::s][Cart::xxz]+amb0*R[Cart::xx][Cart::s][Cart::xxz];
R[Cart::xz][Cart::x][Cart::xxz]=R[Cart::xxz][Cart::s][Cart::xxz]+amb0*R[Cart::xz][Cart::s][Cart::xxz];
R[Cart::zz][Cart::x][Cart::xxz]=R[Cart::xzz][Cart::s][Cart::xxz]+amb0*R[Cart::zz][Cart::s][Cart::xxz];
R[Cart::yy][Cart::z][Cart::xxz]=R[Cart::yyz][Cart::s][Cart::xxz]+amb2*R[Cart::yy][Cart::s][Cart::xxz];
R[Cart::xy][Cart::z][Cart::xxz]=R[Cart::xyz][Cart::s][Cart::xxz]+amb2*R[Cart::xy][Cart::s][Cart::xxz];
R[Cart::yz][Cart::z][Cart::xxz]=R[Cart::yzz][Cart::s][Cart::xxz]+amb2*R[Cart::yz][Cart::s][Cart::xxz];
R[Cart::xx][Cart::z][Cart::xxz]=R[Cart::xxz][Cart::s][Cart::xxz]+amb2*R[Cart::xx][Cart::s][Cart::xxz];
R[Cart::xz][Cart::z][Cart::xxz]=R[Cart::xzz][Cart::s][Cart::xxz]+amb2*R[Cart::xz][Cart::s][Cart::xxz];
R[Cart::zz][Cart::z][Cart::xxz]=R[Cart::zzz][Cart::s][Cart::xxz]+amb2*R[Cart::zz][Cart::s][Cart::xxz];
R[Cart::yy][Cart::y][Cart::xzz]=R[Cart::yyy][Cart::s][Cart::xzz]+amb1*R[Cart::yy][Cart::s][Cart::xzz];
R[Cart::xy][Cart::y][Cart::xzz]=R[Cart::xyy][Cart::s][Cart::xzz]+amb1*R[Cart::xy][Cart::s][Cart::xzz];
R[Cart::yz][Cart::y][Cart::xzz]=R[Cart::yyz][Cart::s][Cart::xzz]+amb1*R[Cart::yz][Cart::s][Cart::xzz];
R[Cart::xx][Cart::y][Cart::xzz]=R[Cart::xxy][Cart::s][Cart::xzz]+amb1*R[Cart::xx][Cart::s][Cart::xzz];
R[Cart::xz][Cart::y][Cart::xzz]=R[Cart::xyz][Cart::s][Cart::xzz]+amb1*R[Cart::xz][Cart::s][Cart::xzz];
R[Cart::zz][Cart::y][Cart::xzz]=R[Cart::yzz][Cart::s][Cart::xzz]+amb1*R[Cart::zz][Cart::s][Cart::xzz];
R[Cart::yy][Cart::x][Cart::xzz]=R[Cart::xyy][Cart::s][Cart::xzz]+amb0*R[Cart::yy][Cart::s][Cart::xzz];
R[Cart::xy][Cart::x][Cart::xzz]=R[Cart::xxy][Cart::s][Cart::xzz]+amb0*R[Cart::xy][Cart::s][Cart::xzz];
R[Cart::yz][Cart::x][Cart::xzz]=R[Cart::xyz][Cart::s][Cart::xzz]+amb0*R[Cart::yz][Cart::s][Cart::xzz];
R[Cart::xx][Cart::x][Cart::xzz]=R[Cart::xxx][Cart::s][Cart::xzz]+amb0*R[Cart::xx][Cart::s][Cart::xzz];
R[Cart::xz][Cart::x][Cart::xzz]=R[Cart::xxz][Cart::s][Cart::xzz]+amb0*R[Cart::xz][Cart::s][Cart::xzz];
R[Cart::zz][Cart::x][Cart::xzz]=R[Cart::xzz][Cart::s][Cart::xzz]+amb0*R[Cart::zz][Cart::s][Cart::xzz];
R[Cart::yy][Cart::z][Cart::xzz]=R[Cart::yyz][Cart::s][Cart::xzz]+amb2*R[Cart::yy][Cart::s][Cart::xzz];
R[Cart::xy][Cart::z][Cart::xzz]=R[Cart::xyz][Cart::s][Cart::xzz]+amb2*R[Cart::xy][Cart::s][Cart::xzz];
R[Cart::yz][Cart::z][Cart::xzz]=R[Cart::yzz][Cart::s][Cart::xzz]+amb2*R[Cart::yz][Cart::s][Cart::xzz];
R[Cart::xx][Cart::z][Cart::xzz]=R[Cart::xxz][Cart::s][Cart::xzz]+amb2*R[Cart::xx][Cart::s][Cart::xzz];
R[Cart::xz][Cart::z][Cart::xzz]=R[Cart::xzz][Cart::s][Cart::xzz]+amb2*R[Cart::xz][Cart::s][Cart::xzz];
R[Cart::zz][Cart::z][Cart::xzz]=R[Cart::zzz][Cart::s][Cart::xzz]+amb2*R[Cart::zz][Cart::s][Cart::xzz];
R[Cart::yy][Cart::y][Cart::zzz]=R[Cart::yyy][Cart::s][Cart::zzz]+amb1*R[Cart::yy][Cart::s][Cart::zzz];
R[Cart::xy][Cart::y][Cart::zzz]=R[Cart::xyy][Cart::s][Cart::zzz]+amb1*R[Cart::xy][Cart::s][Cart::zzz];
R[Cart::yz][Cart::y][Cart::zzz]=R[Cart::yyz][Cart::s][Cart::zzz]+amb1*R[Cart::yz][Cart::s][Cart::zzz];
R[Cart::xx][Cart::y][Cart::zzz]=R[Cart::xxy][Cart::s][Cart::zzz]+amb1*R[Cart::xx][Cart::s][Cart::zzz];
R[Cart::xz][Cart::y][Cart::zzz]=R[Cart::xyz][Cart::s][Cart::zzz]+amb1*R[Cart::xz][Cart::s][Cart::zzz];
R[Cart::zz][Cart::y][Cart::zzz]=R[Cart::yzz][Cart::s][Cart::zzz]+amb1*R[Cart::zz][Cart::s][Cart::zzz];
R[Cart::yy][Cart::x][Cart::zzz]=R[Cart::xyy][Cart::s][Cart::zzz]+amb0*R[Cart::yy][Cart::s][Cart::zzz];
R[Cart::xy][Cart::x][Cart::zzz]=R[Cart::xxy][Cart::s][Cart::zzz]+amb0*R[Cart::xy][Cart::s][Cart::zzz];
R[Cart::yz][Cart::x][Cart::zzz]=R[Cart::xyz][Cart::s][Cart::zzz]+amb0*R[Cart::yz][Cart::s][Cart::zzz];
R[Cart::xx][Cart::x][Cart::zzz]=R[Cart::xxx][Cart::s][Cart::zzz]+amb0*R[Cart::xx][Cart::s][Cart::zzz];
R[Cart::xz][Cart::x][Cart::zzz]=R[Cart::xxz][Cart::s][Cart::zzz]+amb0*R[Cart::xz][Cart::s][Cart::zzz];
R[Cart::zz][Cart::x][Cart::zzz]=R[Cart::xzz][Cart::s][Cart::zzz]+amb0*R[Cart::zz][Cart::s][Cart::zzz];
R[Cart::yy][Cart::z][Cart::zzz]=R[Cart::yyz][Cart::s][Cart::zzz]+amb2*R[Cart::yy][Cart::s][Cart::zzz];
R[Cart::xy][Cart::z][Cart::zzz]=R[Cart::xyz][Cart::s][Cart::zzz]+amb2*R[Cart::xy][Cart::s][Cart::zzz];
R[Cart::yz][Cart::z][Cart::zzz]=R[Cart::yzz][Cart::s][Cart::zzz]+amb2*R[Cart::yz][Cart::s][Cart::zzz];
R[Cart::xx][Cart::z][Cart::zzz]=R[Cart::xxz][Cart::s][Cart::zzz]+amb2*R[Cart::xx][Cart::s][Cart::zzz];
R[Cart::xz][Cart::z][Cart::zzz]=R[Cart::xzz][Cart::s][Cart::zzz]+amb2*R[Cart::xz][Cart::s][Cart::zzz];
R[Cart::zz][Cart::z][Cart::zzz]=R[Cart::zzz][Cart::s][Cart::zzz]+amb2*R[Cart::zz][Cart::s][Cart::zzz];
}
//------------------------------------------------------

//Integral f - p - f
if (_lmax_beta>0 && _lmax_alpha>2 && _lmax_gamma>2){

R[Cart::yyy][Cart::y][Cart::yyy]=R[Cart::yyyy][Cart::s][Cart::yyy]+amb1*R[Cart::yyy][Cart::s][Cart::yyy];
R[Cart::xyy][Cart::y][Cart::yyy]=R[Cart::xyyy][Cart::s][Cart::yyy]+amb1*R[Cart::xyy][Cart::s][Cart::yyy];
R[Cart::yyz][Cart::y][Cart::yyy]=R[Cart::yyyz][Cart::s][Cart::yyy]+amb1*R[Cart::yyz][Cart::s][Cart::yyy];
R[Cart::xxy][Cart::y][Cart::yyy]=R[Cart::xxyy][Cart::s][Cart::yyy]+amb1*R[Cart::xxy][Cart::s][Cart::yyy];
R[Cart::xyz][Cart::y][Cart::yyy]=R[Cart::xyyz][Cart::s][Cart::yyy]+amb1*R[Cart::xyz][Cart::s][Cart::yyy];
R[Cart::yzz][Cart::y][Cart::yyy]=R[Cart::yyzz][Cart::s][Cart::yyy]+amb1*R[Cart::yzz][Cart::s][Cart::yyy];
R[Cart::xxx][Cart::y][Cart::yyy]=R[Cart::xxxy][Cart::s][Cart::yyy]+amb1*R[Cart::xxx][Cart::s][Cart::yyy];
R[Cart::xxz][Cart::y][Cart::yyy]=R[Cart::xxyz][Cart::s][Cart::yyy]+amb1*R[Cart::xxz][Cart::s][Cart::yyy];
R[Cart::xzz][Cart::y][Cart::yyy]=R[Cart::xyzz][Cart::s][Cart::yyy]+amb1*R[Cart::xzz][Cart::s][Cart::yyy];
R[Cart::zzz][Cart::y][Cart::yyy]=R[Cart::yzzz][Cart::s][Cart::yyy]+amb1*R[Cart::zzz][Cart::s][Cart::yyy];
R[Cart::yyy][Cart::x][Cart::yyy]=R[Cart::xyyy][Cart::s][Cart::yyy]+amb0*R[Cart::yyy][Cart::s][Cart::yyy];
R[Cart::xyy][Cart::x][Cart::yyy]=R[Cart::xxyy][Cart::s][Cart::yyy]+amb0*R[Cart::xyy][Cart::s][Cart::yyy];
R[Cart::yyz][Cart::x][Cart::yyy]=R[Cart::xyyz][Cart::s][Cart::yyy]+amb0*R[Cart::yyz][Cart::s][Cart::yyy];
R[Cart::xxy][Cart::x][Cart::yyy]=R[Cart::xxxy][Cart::s][Cart::yyy]+amb0*R[Cart::xxy][Cart::s][Cart::yyy];
R[Cart::xyz][Cart::x][Cart::yyy]=R[Cart::xxyz][Cart::s][Cart::yyy]+amb0*R[Cart::xyz][Cart::s][Cart::yyy];
R[Cart::yzz][Cart::x][Cart::yyy]=R[Cart::xyzz][Cart::s][Cart::yyy]+amb0*R[Cart::yzz][Cart::s][Cart::yyy];
R[Cart::xxx][Cart::x][Cart::yyy]=R[Cart::xxxx][Cart::s][Cart::yyy]+amb0*R[Cart::xxx][Cart::s][Cart::yyy];
R[Cart::xxz][Cart::x][Cart::yyy]=R[Cart::xxxz][Cart::s][Cart::yyy]+amb0*R[Cart::xxz][Cart::s][Cart::yyy];
R[Cart::xzz][Cart::x][Cart::yyy]=R[Cart::xxzz][Cart::s][Cart::yyy]+amb0*R[Cart::xzz][Cart::s][Cart::yyy];
R[Cart::zzz][Cart::x][Cart::yyy]=R[Cart::xzzz][Cart::s][Cart::yyy]+amb0*R[Cart::zzz][Cart::s][Cart::yyy];
R[Cart::yyy][Cart::z][Cart::yyy]=R[Cart::yyyz][Cart::s][Cart::yyy]+amb2*R[Cart::yyy][Cart::s][Cart::yyy];
R[Cart::xyy][Cart::z][Cart::yyy]=R[Cart::xyyz][Cart::s][Cart::yyy]+amb2*R[Cart::xyy][Cart::s][Cart::yyy];
R[Cart::yyz][Cart::z][Cart::yyy]=R[Cart::yyzz][Cart::s][Cart::yyy]+amb2*R[Cart::yyz][Cart::s][Cart::yyy];
R[Cart::xxy][Cart::z][Cart::yyy]=R[Cart::xxyz][Cart::s][Cart::yyy]+amb2*R[Cart::xxy][Cart::s][Cart::yyy];
R[Cart::xyz][Cart::z][Cart::yyy]=R[Cart::xyzz][Cart::s][Cart::yyy]+amb2*R[Cart::xyz][Cart::s][Cart::yyy];
R[Cart::yzz][Cart::z][Cart::yyy]=R[Cart::yzzz][Cart::s][Cart::yyy]+amb2*R[Cart::yzz][Cart::s][Cart::yyy];
R[Cart::xxx][Cart::z][Cart::yyy]=R[Cart::xxxz][Cart::s][Cart::yyy]+amb2*R[Cart::xxx][Cart::s][Cart::yyy];
R[Cart::xxz][Cart::z][Cart::yyy]=R[Cart::xxzz][Cart::s][Cart::yyy]+amb2*R[Cart::xxz][Cart::s][Cart::yyy];
R[Cart::xzz][Cart::z][Cart::yyy]=R[Cart::xzzz][Cart::s][Cart::yyy]+amb2*R[Cart::xzz][Cart::s][Cart::yyy];
R[Cart::zzz][Cart::z][Cart::yyy]=R[Cart::zzzz][Cart::s][Cart::yyy]+amb2*R[Cart::zzz][Cart::s][Cart::yyy];
R[Cart::yyy][Cart::y][Cart::xyy]=R[Cart::yyyy][Cart::s][Cart::xyy]+amb1*R[Cart::yyy][Cart::s][Cart::xyy];
R[Cart::xyy][Cart::y][Cart::xyy]=R[Cart::xyyy][Cart::s][Cart::xyy]+amb1*R[Cart::xyy][Cart::s][Cart::xyy];
R[Cart::yyz][Cart::y][Cart::xyy]=R[Cart::yyyz][Cart::s][Cart::xyy]+amb1*R[Cart::yyz][Cart::s][Cart::xyy];
R[Cart::xxy][Cart::y][Cart::xyy]=R[Cart::xxyy][Cart::s][Cart::xyy]+amb1*R[Cart::xxy][Cart::s][Cart::xyy];
R[Cart::xyz][Cart::y][Cart::xyy]=R[Cart::xyyz][Cart::s][Cart::xyy]+amb1*R[Cart::xyz][Cart::s][Cart::xyy];
R[Cart::yzz][Cart::y][Cart::xyy]=R[Cart::yyzz][Cart::s][Cart::xyy]+amb1*R[Cart::yzz][Cart::s][Cart::xyy];
R[Cart::xxx][Cart::y][Cart::xyy]=R[Cart::xxxy][Cart::s][Cart::xyy]+amb1*R[Cart::xxx][Cart::s][Cart::xyy];
R[Cart::xxz][Cart::y][Cart::xyy]=R[Cart::xxyz][Cart::s][Cart::xyy]+amb1*R[Cart::xxz][Cart::s][Cart::xyy];
R[Cart::xzz][Cart::y][Cart::xyy]=R[Cart::xyzz][Cart::s][Cart::xyy]+amb1*R[Cart::xzz][Cart::s][Cart::xyy];
R[Cart::zzz][Cart::y][Cart::xyy]=R[Cart::yzzz][Cart::s][Cart::xyy]+amb1*R[Cart::zzz][Cart::s][Cart::xyy];
R[Cart::yyy][Cart::x][Cart::xyy]=R[Cart::xyyy][Cart::s][Cart::xyy]+amb0*R[Cart::yyy][Cart::s][Cart::xyy];
R[Cart::xyy][Cart::x][Cart::xyy]=R[Cart::xxyy][Cart::s][Cart::xyy]+amb0*R[Cart::xyy][Cart::s][Cart::xyy];
R[Cart::yyz][Cart::x][Cart::xyy]=R[Cart::xyyz][Cart::s][Cart::xyy]+amb0*R[Cart::yyz][Cart::s][Cart::xyy];
R[Cart::xxy][Cart::x][Cart::xyy]=R[Cart::xxxy][Cart::s][Cart::xyy]+amb0*R[Cart::xxy][Cart::s][Cart::xyy];
R[Cart::xyz][Cart::x][Cart::xyy]=R[Cart::xxyz][Cart::s][Cart::xyy]+amb0*R[Cart::xyz][Cart::s][Cart::xyy];
R[Cart::yzz][Cart::x][Cart::xyy]=R[Cart::xyzz][Cart::s][Cart::xyy]+amb0*R[Cart::yzz][Cart::s][Cart::xyy];
R[Cart::xxx][Cart::x][Cart::xyy]=R[Cart::xxxx][Cart::s][Cart::xyy]+amb0*R[Cart::xxx][Cart::s][Cart::xyy];
R[Cart::xxz][Cart::x][Cart::xyy]=R[Cart::xxxz][Cart::s][Cart::xyy]+amb0*R[Cart::xxz][Cart::s][Cart::xyy];
R[Cart::xzz][Cart::x][Cart::xyy]=R[Cart::xxzz][Cart::s][Cart::xyy]+amb0*R[Cart::xzz][Cart::s][Cart::xyy];
R[Cart::zzz][Cart::x][Cart::xyy]=R[Cart::xzzz][Cart::s][Cart::xyy]+amb0*R[Cart::zzz][Cart::s][Cart::xyy];
R[Cart::yyy][Cart::z][Cart::xyy]=R[Cart::yyyz][Cart::s][Cart::xyy]+amb2*R[Cart::yyy][Cart::s][Cart::xyy];
R[Cart::xyy][Cart::z][Cart::xyy]=R[Cart::xyyz][Cart::s][Cart::xyy]+amb2*R[Cart::xyy][Cart::s][Cart::xyy];
R[Cart::yyz][Cart::z][Cart::xyy]=R[Cart::yyzz][Cart::s][Cart::xyy]+amb2*R[Cart::yyz][Cart::s][Cart::xyy];
R[Cart::xxy][Cart::z][Cart::xyy]=R[Cart::xxyz][Cart::s][Cart::xyy]+amb2*R[Cart::xxy][Cart::s][Cart::xyy];
R[Cart::xyz][Cart::z][Cart::xyy]=R[Cart::xyzz][Cart::s][Cart::xyy]+amb2*R[Cart::xyz][Cart::s][Cart::xyy];
R[Cart::yzz][Cart::z][Cart::xyy]=R[Cart::yzzz][Cart::s][Cart::xyy]+amb2*R[Cart::yzz][Cart::s][Cart::xyy];
R[Cart::xxx][Cart::z][Cart::xyy]=R[Cart::xxxz][Cart::s][Cart::xyy]+amb2*R[Cart::xxx][Cart::s][Cart::xyy];
R[Cart::xxz][Cart::z][Cart::xyy]=R[Cart::xxzz][Cart::s][Cart::xyy]+amb2*R[Cart::xxz][Cart::s][Cart::xyy];
R[Cart::xzz][Cart::z][Cart::xyy]=R[Cart::xzzz][Cart::s][Cart::xyy]+amb2*R[Cart::xzz][Cart::s][Cart::xyy];
R[Cart::zzz][Cart::z][Cart::xyy]=R[Cart::zzzz][Cart::s][Cart::xyy]+amb2*R[Cart::zzz][Cart::s][Cart::xyy];
R[Cart::yyy][Cart::y][Cart::yyz]=R[Cart::yyyy][Cart::s][Cart::yyz]+amb1*R[Cart::yyy][Cart::s][Cart::yyz];
R[Cart::xyy][Cart::y][Cart::yyz]=R[Cart::xyyy][Cart::s][Cart::yyz]+amb1*R[Cart::xyy][Cart::s][Cart::yyz];
R[Cart::yyz][Cart::y][Cart::yyz]=R[Cart::yyyz][Cart::s][Cart::yyz]+amb1*R[Cart::yyz][Cart::s][Cart::yyz];
R[Cart::xxy][Cart::y][Cart::yyz]=R[Cart::xxyy][Cart::s][Cart::yyz]+amb1*R[Cart::xxy][Cart::s][Cart::yyz];
R[Cart::xyz][Cart::y][Cart::yyz]=R[Cart::xyyz][Cart::s][Cart::yyz]+amb1*R[Cart::xyz][Cart::s][Cart::yyz];
R[Cart::yzz][Cart::y][Cart::yyz]=R[Cart::yyzz][Cart::s][Cart::yyz]+amb1*R[Cart::yzz][Cart::s][Cart::yyz];
R[Cart::xxx][Cart::y][Cart::yyz]=R[Cart::xxxy][Cart::s][Cart::yyz]+amb1*R[Cart::xxx][Cart::s][Cart::yyz];
R[Cart::xxz][Cart::y][Cart::yyz]=R[Cart::xxyz][Cart::s][Cart::yyz]+amb1*R[Cart::xxz][Cart::s][Cart::yyz];
R[Cart::xzz][Cart::y][Cart::yyz]=R[Cart::xyzz][Cart::s][Cart::yyz]+amb1*R[Cart::xzz][Cart::s][Cart::yyz];
R[Cart::zzz][Cart::y][Cart::yyz]=R[Cart::yzzz][Cart::s][Cart::yyz]+amb1*R[Cart::zzz][Cart::s][Cart::yyz];
R[Cart::yyy][Cart::x][Cart::yyz]=R[Cart::xyyy][Cart::s][Cart::yyz]+amb0*R[Cart::yyy][Cart::s][Cart::yyz];
R[Cart::xyy][Cart::x][Cart::yyz]=R[Cart::xxyy][Cart::s][Cart::yyz]+amb0*R[Cart::xyy][Cart::s][Cart::yyz];
R[Cart::yyz][Cart::x][Cart::yyz]=R[Cart::xyyz][Cart::s][Cart::yyz]+amb0*R[Cart::yyz][Cart::s][Cart::yyz];
R[Cart::xxy][Cart::x][Cart::yyz]=R[Cart::xxxy][Cart::s][Cart::yyz]+amb0*R[Cart::xxy][Cart::s][Cart::yyz];
R[Cart::xyz][Cart::x][Cart::yyz]=R[Cart::xxyz][Cart::s][Cart::yyz]+amb0*R[Cart::xyz][Cart::s][Cart::yyz];
R[Cart::yzz][Cart::x][Cart::yyz]=R[Cart::xyzz][Cart::s][Cart::yyz]+amb0*R[Cart::yzz][Cart::s][Cart::yyz];
R[Cart::xxx][Cart::x][Cart::yyz]=R[Cart::xxxx][Cart::s][Cart::yyz]+amb0*R[Cart::xxx][Cart::s][Cart::yyz];
R[Cart::xxz][Cart::x][Cart::yyz]=R[Cart::xxxz][Cart::s][Cart::yyz]+amb0*R[Cart::xxz][Cart::s][Cart::yyz];
R[Cart::xzz][Cart::x][Cart::yyz]=R[Cart::xxzz][Cart::s][Cart::yyz]+amb0*R[Cart::xzz][Cart::s][Cart::yyz];
R[Cart::zzz][Cart::x][Cart::yyz]=R[Cart::xzzz][Cart::s][Cart::yyz]+amb0*R[Cart::zzz][Cart::s][Cart::yyz];
R[Cart::yyy][Cart::z][Cart::yyz]=R[Cart::yyyz][Cart::s][Cart::yyz]+amb2*R[Cart::yyy][Cart::s][Cart::yyz];
R[Cart::xyy][Cart::z][Cart::yyz]=R[Cart::xyyz][Cart::s][Cart::yyz]+amb2*R[Cart::xyy][Cart::s][Cart::yyz];
R[Cart::yyz][Cart::z][Cart::yyz]=R[Cart::yyzz][Cart::s][Cart::yyz]+amb2*R[Cart::yyz][Cart::s][Cart::yyz];
R[Cart::xxy][Cart::z][Cart::yyz]=R[Cart::xxyz][Cart::s][Cart::yyz]+amb2*R[Cart::xxy][Cart::s][Cart::yyz];
R[Cart::xyz][Cart::z][Cart::yyz]=R[Cart::xyzz][Cart::s][Cart::yyz]+amb2*R[Cart::xyz][Cart::s][Cart::yyz];
R[Cart::yzz][Cart::z][Cart::yyz]=R[Cart::yzzz][Cart::s][Cart::yyz]+amb2*R[Cart::yzz][Cart::s][Cart::yyz];
R[Cart::xxx][Cart::z][Cart::yyz]=R[Cart::xxxz][Cart::s][Cart::yyz]+amb2*R[Cart::xxx][Cart::s][Cart::yyz];
R[Cart::xxz][Cart::z][Cart::yyz]=R[Cart::xxzz][Cart::s][Cart::yyz]+amb2*R[Cart::xxz][Cart::s][Cart::yyz];
R[Cart::xzz][Cart::z][Cart::yyz]=R[Cart::xzzz][Cart::s][Cart::yyz]+amb2*R[Cart::xzz][Cart::s][Cart::yyz];
R[Cart::zzz][Cart::z][Cart::yyz]=R[Cart::zzzz][Cart::s][Cart::yyz]+amb2*R[Cart::zzz][Cart::s][Cart::yyz];
R[Cart::yyy][Cart::y][Cart::xxy]=R[Cart::yyyy][Cart::s][Cart::xxy]+amb1*R[Cart::yyy][Cart::s][Cart::xxy];
R[Cart::xyy][Cart::y][Cart::xxy]=R[Cart::xyyy][Cart::s][Cart::xxy]+amb1*R[Cart::xyy][Cart::s][Cart::xxy];
R[Cart::yyz][Cart::y][Cart::xxy]=R[Cart::yyyz][Cart::s][Cart::xxy]+amb1*R[Cart::yyz][Cart::s][Cart::xxy];
R[Cart::xxy][Cart::y][Cart::xxy]=R[Cart::xxyy][Cart::s][Cart::xxy]+amb1*R[Cart::xxy][Cart::s][Cart::xxy];
R[Cart::xyz][Cart::y][Cart::xxy]=R[Cart::xyyz][Cart::s][Cart::xxy]+amb1*R[Cart::xyz][Cart::s][Cart::xxy];
R[Cart::yzz][Cart::y][Cart::xxy]=R[Cart::yyzz][Cart::s][Cart::xxy]+amb1*R[Cart::yzz][Cart::s][Cart::xxy];
R[Cart::xxx][Cart::y][Cart::xxy]=R[Cart::xxxy][Cart::s][Cart::xxy]+amb1*R[Cart::xxx][Cart::s][Cart::xxy];
R[Cart::xxz][Cart::y][Cart::xxy]=R[Cart::xxyz][Cart::s][Cart::xxy]+amb1*R[Cart::xxz][Cart::s][Cart::xxy];
R[Cart::xzz][Cart::y][Cart::xxy]=R[Cart::xyzz][Cart::s][Cart::xxy]+amb1*R[Cart::xzz][Cart::s][Cart::xxy];
R[Cart::zzz][Cart::y][Cart::xxy]=R[Cart::yzzz][Cart::s][Cart::xxy]+amb1*R[Cart::zzz][Cart::s][Cart::xxy];
R[Cart::yyy][Cart::x][Cart::xxy]=R[Cart::xyyy][Cart::s][Cart::xxy]+amb0*R[Cart::yyy][Cart::s][Cart::xxy];
R[Cart::xyy][Cart::x][Cart::xxy]=R[Cart::xxyy][Cart::s][Cart::xxy]+amb0*R[Cart::xyy][Cart::s][Cart::xxy];
R[Cart::yyz][Cart::x][Cart::xxy]=R[Cart::xyyz][Cart::s][Cart::xxy]+amb0*R[Cart::yyz][Cart::s][Cart::xxy];
R[Cart::xxy][Cart::x][Cart::xxy]=R[Cart::xxxy][Cart::s][Cart::xxy]+amb0*R[Cart::xxy][Cart::s][Cart::xxy];
R[Cart::xyz][Cart::x][Cart::xxy]=R[Cart::xxyz][Cart::s][Cart::xxy]+amb0*R[Cart::xyz][Cart::s][Cart::xxy];
R[Cart::yzz][Cart::x][Cart::xxy]=R[Cart::xyzz][Cart::s][Cart::xxy]+amb0*R[Cart::yzz][Cart::s][Cart::xxy];
R[Cart::xxx][Cart::x][Cart::xxy]=R[Cart::xxxx][Cart::s][Cart::xxy]+amb0*R[Cart::xxx][Cart::s][Cart::xxy];
R[Cart::xxz][Cart::x][Cart::xxy]=R[Cart::xxxz][Cart::s][Cart::xxy]+amb0*R[Cart::xxz][Cart::s][Cart::xxy];
R[Cart::xzz][Cart::x][Cart::xxy]=R[Cart::xxzz][Cart::s][Cart::xxy]+amb0*R[Cart::xzz][Cart::s][Cart::xxy];
R[Cart::zzz][Cart::x][Cart::xxy]=R[Cart::xzzz][Cart::s][Cart::xxy]+amb0*R[Cart::zzz][Cart::s][Cart::xxy];
R[Cart::yyy][Cart::z][Cart::xxy]=R[Cart::yyyz][Cart::s][Cart::xxy]+amb2*R[Cart::yyy][Cart::s][Cart::xxy];
R[Cart::xyy][Cart::z][Cart::xxy]=R[Cart::xyyz][Cart::s][Cart::xxy]+amb2*R[Cart::xyy][Cart::s][Cart::xxy];
R[Cart::yyz][Cart::z][Cart::xxy]=R[Cart::yyzz][Cart::s][Cart::xxy]+amb2*R[Cart::yyz][Cart::s][Cart::xxy];
R[Cart::xxy][Cart::z][Cart::xxy]=R[Cart::xxyz][Cart::s][Cart::xxy]+amb2*R[Cart::xxy][Cart::s][Cart::xxy];
R[Cart::xyz][Cart::z][Cart::xxy]=R[Cart::xyzz][Cart::s][Cart::xxy]+amb2*R[Cart::xyz][Cart::s][Cart::xxy];
R[Cart::yzz][Cart::z][Cart::xxy]=R[Cart::yzzz][Cart::s][Cart::xxy]+amb2*R[Cart::yzz][Cart::s][Cart::xxy];
R[Cart::xxx][Cart::z][Cart::xxy]=R[Cart::xxxz][Cart::s][Cart::xxy]+amb2*R[Cart::xxx][Cart::s][Cart::xxy];
R[Cart::xxz][Cart::z][Cart::xxy]=R[Cart::xxzz][Cart::s][Cart::xxy]+amb2*R[Cart::xxz][Cart::s][Cart::xxy];
R[Cart::xzz][Cart::z][Cart::xxy]=R[Cart::xzzz][Cart::s][Cart::xxy]+amb2*R[Cart::xzz][Cart::s][Cart::xxy];
R[Cart::zzz][Cart::z][Cart::xxy]=R[Cart::zzzz][Cart::s][Cart::xxy]+amb2*R[Cart::zzz][Cart::s][Cart::xxy];
R[Cart::yyy][Cart::y][Cart::xyz]=R[Cart::yyyy][Cart::s][Cart::xyz]+amb1*R[Cart::yyy][Cart::s][Cart::xyz];
R[Cart::xyy][Cart::y][Cart::xyz]=R[Cart::xyyy][Cart::s][Cart::xyz]+amb1*R[Cart::xyy][Cart::s][Cart::xyz];
R[Cart::yyz][Cart::y][Cart::xyz]=R[Cart::yyyz][Cart::s][Cart::xyz]+amb1*R[Cart::yyz][Cart::s][Cart::xyz];
R[Cart::xxy][Cart::y][Cart::xyz]=R[Cart::xxyy][Cart::s][Cart::xyz]+amb1*R[Cart::xxy][Cart::s][Cart::xyz];
R[Cart::xyz][Cart::y][Cart::xyz]=R[Cart::xyyz][Cart::s][Cart::xyz]+amb1*R[Cart::xyz][Cart::s][Cart::xyz];
R[Cart::yzz][Cart::y][Cart::xyz]=R[Cart::yyzz][Cart::s][Cart::xyz]+amb1*R[Cart::yzz][Cart::s][Cart::xyz];
R[Cart::xxx][Cart::y][Cart::xyz]=R[Cart::xxxy][Cart::s][Cart::xyz]+amb1*R[Cart::xxx][Cart::s][Cart::xyz];
R[Cart::xxz][Cart::y][Cart::xyz]=R[Cart::xxyz][Cart::s][Cart::xyz]+amb1*R[Cart::xxz][Cart::s][Cart::xyz];
R[Cart::xzz][Cart::y][Cart::xyz]=R[Cart::xyzz][Cart::s][Cart::xyz]+amb1*R[Cart::xzz][Cart::s][Cart::xyz];
R[Cart::zzz][Cart::y][Cart::xyz]=R[Cart::yzzz][Cart::s][Cart::xyz]+amb1*R[Cart::zzz][Cart::s][Cart::xyz];
R[Cart::yyy][Cart::x][Cart::xyz]=R[Cart::xyyy][Cart::s][Cart::xyz]+amb0*R[Cart::yyy][Cart::s][Cart::xyz];
R[Cart::xyy][Cart::x][Cart::xyz]=R[Cart::xxyy][Cart::s][Cart::xyz]+amb0*R[Cart::xyy][Cart::s][Cart::xyz];
R[Cart::yyz][Cart::x][Cart::xyz]=R[Cart::xyyz][Cart::s][Cart::xyz]+amb0*R[Cart::yyz][Cart::s][Cart::xyz];
R[Cart::xxy][Cart::x][Cart::xyz]=R[Cart::xxxy][Cart::s][Cart::xyz]+amb0*R[Cart::xxy][Cart::s][Cart::xyz];
R[Cart::xyz][Cart::x][Cart::xyz]=R[Cart::xxyz][Cart::s][Cart::xyz]+amb0*R[Cart::xyz][Cart::s][Cart::xyz];
R[Cart::yzz][Cart::x][Cart::xyz]=R[Cart::xyzz][Cart::s][Cart::xyz]+amb0*R[Cart::yzz][Cart::s][Cart::xyz];
R[Cart::xxx][Cart::x][Cart::xyz]=R[Cart::xxxx][Cart::s][Cart::xyz]+amb0*R[Cart::xxx][Cart::s][Cart::xyz];
R[Cart::xxz][Cart::x][Cart::xyz]=R[Cart::xxxz][Cart::s][Cart::xyz]+amb0*R[Cart::xxz][Cart::s][Cart::xyz];
R[Cart::xzz][Cart::x][Cart::xyz]=R[Cart::xxzz][Cart::s][Cart::xyz]+amb0*R[Cart::xzz][Cart::s][Cart::xyz];
R[Cart::zzz][Cart::x][Cart::xyz]=R[Cart::xzzz][Cart::s][Cart::xyz]+amb0*R[Cart::zzz][Cart::s][Cart::xyz];
R[Cart::yyy][Cart::z][Cart::xyz]=R[Cart::yyyz][Cart::s][Cart::xyz]+amb2*R[Cart::yyy][Cart::s][Cart::xyz];
R[Cart::xyy][Cart::z][Cart::xyz]=R[Cart::xyyz][Cart::s][Cart::xyz]+amb2*R[Cart::xyy][Cart::s][Cart::xyz];
R[Cart::yyz][Cart::z][Cart::xyz]=R[Cart::yyzz][Cart::s][Cart::xyz]+amb2*R[Cart::yyz][Cart::s][Cart::xyz];
R[Cart::xxy][Cart::z][Cart::xyz]=R[Cart::xxyz][Cart::s][Cart::xyz]+amb2*R[Cart::xxy][Cart::s][Cart::xyz];
R[Cart::xyz][Cart::z][Cart::xyz]=R[Cart::xyzz][Cart::s][Cart::xyz]+amb2*R[Cart::xyz][Cart::s][Cart::xyz];
R[Cart::yzz][Cart::z][Cart::xyz]=R[Cart::yzzz][Cart::s][Cart::xyz]+amb2*R[Cart::yzz][Cart::s][Cart::xyz];
R[Cart::xxx][Cart::z][Cart::xyz]=R[Cart::xxxz][Cart::s][Cart::xyz]+amb2*R[Cart::xxx][Cart::s][Cart::xyz];
R[Cart::xxz][Cart::z][Cart::xyz]=R[Cart::xxzz][Cart::s][Cart::xyz]+amb2*R[Cart::xxz][Cart::s][Cart::xyz];
R[Cart::xzz][Cart::z][Cart::xyz]=R[Cart::xzzz][Cart::s][Cart::xyz]+amb2*R[Cart::xzz][Cart::s][Cart::xyz];
R[Cart::zzz][Cart::z][Cart::xyz]=R[Cart::zzzz][Cart::s][Cart::xyz]+amb2*R[Cart::zzz][Cart::s][Cart::xyz];
R[Cart::yyy][Cart::y][Cart::yzz]=R[Cart::yyyy][Cart::s][Cart::yzz]+amb1*R[Cart::yyy][Cart::s][Cart::yzz];
R[Cart::xyy][Cart::y][Cart::yzz]=R[Cart::xyyy][Cart::s][Cart::yzz]+amb1*R[Cart::xyy][Cart::s][Cart::yzz];
R[Cart::yyz][Cart::y][Cart::yzz]=R[Cart::yyyz][Cart::s][Cart::yzz]+amb1*R[Cart::yyz][Cart::s][Cart::yzz];
R[Cart::xxy][Cart::y][Cart::yzz]=R[Cart::xxyy][Cart::s][Cart::yzz]+amb1*R[Cart::xxy][Cart::s][Cart::yzz];
R[Cart::xyz][Cart::y][Cart::yzz]=R[Cart::xyyz][Cart::s][Cart::yzz]+amb1*R[Cart::xyz][Cart::s][Cart::yzz];
R[Cart::yzz][Cart::y][Cart::yzz]=R[Cart::yyzz][Cart::s][Cart::yzz]+amb1*R[Cart::yzz][Cart::s][Cart::yzz];
R[Cart::xxx][Cart::y][Cart::yzz]=R[Cart::xxxy][Cart::s][Cart::yzz]+amb1*R[Cart::xxx][Cart::s][Cart::yzz];
R[Cart::xxz][Cart::y][Cart::yzz]=R[Cart::xxyz][Cart::s][Cart::yzz]+amb1*R[Cart::xxz][Cart::s][Cart::yzz];
R[Cart::xzz][Cart::y][Cart::yzz]=R[Cart::xyzz][Cart::s][Cart::yzz]+amb1*R[Cart::xzz][Cart::s][Cart::yzz];
R[Cart::zzz][Cart::y][Cart::yzz]=R[Cart::yzzz][Cart::s][Cart::yzz]+amb1*R[Cart::zzz][Cart::s][Cart::yzz];
R[Cart::yyy][Cart::x][Cart::yzz]=R[Cart::xyyy][Cart::s][Cart::yzz]+amb0*R[Cart::yyy][Cart::s][Cart::yzz];
R[Cart::xyy][Cart::x][Cart::yzz]=R[Cart::xxyy][Cart::s][Cart::yzz]+amb0*R[Cart::xyy][Cart::s][Cart::yzz];
R[Cart::yyz][Cart::x][Cart::yzz]=R[Cart::xyyz][Cart::s][Cart::yzz]+amb0*R[Cart::yyz][Cart::s][Cart::yzz];
R[Cart::xxy][Cart::x][Cart::yzz]=R[Cart::xxxy][Cart::s][Cart::yzz]+amb0*R[Cart::xxy][Cart::s][Cart::yzz];
R[Cart::xyz][Cart::x][Cart::yzz]=R[Cart::xxyz][Cart::s][Cart::yzz]+amb0*R[Cart::xyz][Cart::s][Cart::yzz];
R[Cart::yzz][Cart::x][Cart::yzz]=R[Cart::xyzz][Cart::s][Cart::yzz]+amb0*R[Cart::yzz][Cart::s][Cart::yzz];
R[Cart::xxx][Cart::x][Cart::yzz]=R[Cart::xxxx][Cart::s][Cart::yzz]+amb0*R[Cart::xxx][Cart::s][Cart::yzz];
R[Cart::xxz][Cart::x][Cart::yzz]=R[Cart::xxxz][Cart::s][Cart::yzz]+amb0*R[Cart::xxz][Cart::s][Cart::yzz];
R[Cart::xzz][Cart::x][Cart::yzz]=R[Cart::xxzz][Cart::s][Cart::yzz]+amb0*R[Cart::xzz][Cart::s][Cart::yzz];
R[Cart::zzz][Cart::x][Cart::yzz]=R[Cart::xzzz][Cart::s][Cart::yzz]+amb0*R[Cart::zzz][Cart::s][Cart::yzz];
R[Cart::yyy][Cart::z][Cart::yzz]=R[Cart::yyyz][Cart::s][Cart::yzz]+amb2*R[Cart::yyy][Cart::s][Cart::yzz];
R[Cart::xyy][Cart::z][Cart::yzz]=R[Cart::xyyz][Cart::s][Cart::yzz]+amb2*R[Cart::xyy][Cart::s][Cart::yzz];
R[Cart::yyz][Cart::z][Cart::yzz]=R[Cart::yyzz][Cart::s][Cart::yzz]+amb2*R[Cart::yyz][Cart::s][Cart::yzz];
R[Cart::xxy][Cart::z][Cart::yzz]=R[Cart::xxyz][Cart::s][Cart::yzz]+amb2*R[Cart::xxy][Cart::s][Cart::yzz];
R[Cart::xyz][Cart::z][Cart::yzz]=R[Cart::xyzz][Cart::s][Cart::yzz]+amb2*R[Cart::xyz][Cart::s][Cart::yzz];
R[Cart::yzz][Cart::z][Cart::yzz]=R[Cart::yzzz][Cart::s][Cart::yzz]+amb2*R[Cart::yzz][Cart::s][Cart::yzz];
R[Cart::xxx][Cart::z][Cart::yzz]=R[Cart::xxxz][Cart::s][Cart::yzz]+amb2*R[Cart::xxx][Cart::s][Cart::yzz];
R[Cart::xxz][Cart::z][Cart::yzz]=R[Cart::xxzz][Cart::s][Cart::yzz]+amb2*R[Cart::xxz][Cart::s][Cart::yzz];
R[Cart::xzz][Cart::z][Cart::yzz]=R[Cart::xzzz][Cart::s][Cart::yzz]+amb2*R[Cart::xzz][Cart::s][Cart::yzz];
R[Cart::zzz][Cart::z][Cart::yzz]=R[Cart::zzzz][Cart::s][Cart::yzz]+amb2*R[Cart::zzz][Cart::s][Cart::yzz];
R[Cart::yyy][Cart::y][Cart::xxx]=R[Cart::yyyy][Cart::s][Cart::xxx]+amb1*R[Cart::yyy][Cart::s][Cart::xxx];
R[Cart::xyy][Cart::y][Cart::xxx]=R[Cart::xyyy][Cart::s][Cart::xxx]+amb1*R[Cart::xyy][Cart::s][Cart::xxx];
R[Cart::yyz][Cart::y][Cart::xxx]=R[Cart::yyyz][Cart::s][Cart::xxx]+amb1*R[Cart::yyz][Cart::s][Cart::xxx];
R[Cart::xxy][Cart::y][Cart::xxx]=R[Cart::xxyy][Cart::s][Cart::xxx]+amb1*R[Cart::xxy][Cart::s][Cart::xxx];
R[Cart::xyz][Cart::y][Cart::xxx]=R[Cart::xyyz][Cart::s][Cart::xxx]+amb1*R[Cart::xyz][Cart::s][Cart::xxx];
R[Cart::yzz][Cart::y][Cart::xxx]=R[Cart::yyzz][Cart::s][Cart::xxx]+amb1*R[Cart::yzz][Cart::s][Cart::xxx];
R[Cart::xxx][Cart::y][Cart::xxx]=R[Cart::xxxy][Cart::s][Cart::xxx]+amb1*R[Cart::xxx][Cart::s][Cart::xxx];
R[Cart::xxz][Cart::y][Cart::xxx]=R[Cart::xxyz][Cart::s][Cart::xxx]+amb1*R[Cart::xxz][Cart::s][Cart::xxx];
R[Cart::xzz][Cart::y][Cart::xxx]=R[Cart::xyzz][Cart::s][Cart::xxx]+amb1*R[Cart::xzz][Cart::s][Cart::xxx];
R[Cart::zzz][Cart::y][Cart::xxx]=R[Cart::yzzz][Cart::s][Cart::xxx]+amb1*R[Cart::zzz][Cart::s][Cart::xxx];
R[Cart::yyy][Cart::x][Cart::xxx]=R[Cart::xyyy][Cart::s][Cart::xxx]+amb0*R[Cart::yyy][Cart::s][Cart::xxx];
R[Cart::xyy][Cart::x][Cart::xxx]=R[Cart::xxyy][Cart::s][Cart::xxx]+amb0*R[Cart::xyy][Cart::s][Cart::xxx];
R[Cart::yyz][Cart::x][Cart::xxx]=R[Cart::xyyz][Cart::s][Cart::xxx]+amb0*R[Cart::yyz][Cart::s][Cart::xxx];
R[Cart::xxy][Cart::x][Cart::xxx]=R[Cart::xxxy][Cart::s][Cart::xxx]+amb0*R[Cart::xxy][Cart::s][Cart::xxx];
R[Cart::xyz][Cart::x][Cart::xxx]=R[Cart::xxyz][Cart::s][Cart::xxx]+amb0*R[Cart::xyz][Cart::s][Cart::xxx];
R[Cart::yzz][Cart::x][Cart::xxx]=R[Cart::xyzz][Cart::s][Cart::xxx]+amb0*R[Cart::yzz][Cart::s][Cart::xxx];
R[Cart::xxx][Cart::x][Cart::xxx]=R[Cart::xxxx][Cart::s][Cart::xxx]+amb0*R[Cart::xxx][Cart::s][Cart::xxx];
R[Cart::xxz][Cart::x][Cart::xxx]=R[Cart::xxxz][Cart::s][Cart::xxx]+amb0*R[Cart::xxz][Cart::s][Cart::xxx];
R[Cart::xzz][Cart::x][Cart::xxx]=R[Cart::xxzz][Cart::s][Cart::xxx]+amb0*R[Cart::xzz][Cart::s][Cart::xxx];
R[Cart::zzz][Cart::x][Cart::xxx]=R[Cart::xzzz][Cart::s][Cart::xxx]+amb0*R[Cart::zzz][Cart::s][Cart::xxx];
R[Cart::yyy][Cart::z][Cart::xxx]=R[Cart::yyyz][Cart::s][Cart::xxx]+amb2*R[Cart::yyy][Cart::s][Cart::xxx];
R[Cart::xyy][Cart::z][Cart::xxx]=R[Cart::xyyz][Cart::s][Cart::xxx]+amb2*R[Cart::xyy][Cart::s][Cart::xxx];
R[Cart::yyz][Cart::z][Cart::xxx]=R[Cart::yyzz][Cart::s][Cart::xxx]+amb2*R[Cart::yyz][Cart::s][Cart::xxx];
R[Cart::xxy][Cart::z][Cart::xxx]=R[Cart::xxyz][Cart::s][Cart::xxx]+amb2*R[Cart::xxy][Cart::s][Cart::xxx];
R[Cart::xyz][Cart::z][Cart::xxx]=R[Cart::xyzz][Cart::s][Cart::xxx]+amb2*R[Cart::xyz][Cart::s][Cart::xxx];
R[Cart::yzz][Cart::z][Cart::xxx]=R[Cart::yzzz][Cart::s][Cart::xxx]+amb2*R[Cart::yzz][Cart::s][Cart::xxx];
R[Cart::xxx][Cart::z][Cart::xxx]=R[Cart::xxxz][Cart::s][Cart::xxx]+amb2*R[Cart::xxx][Cart::s][Cart::xxx];
R[Cart::xxz][Cart::z][Cart::xxx]=R[Cart::xxzz][Cart::s][Cart::xxx]+amb2*R[Cart::xxz][Cart::s][Cart::xxx];
R[Cart::xzz][Cart::z][Cart::xxx]=R[Cart::xzzz][Cart::s][Cart::xxx]+amb2*R[Cart::xzz][Cart::s][Cart::xxx];
R[Cart::zzz][Cart::z][Cart::xxx]=R[Cart::zzzz][Cart::s][Cart::xxx]+amb2*R[Cart::zzz][Cart::s][Cart::xxx];
R[Cart::yyy][Cart::y][Cart::xxz]=R[Cart::yyyy][Cart::s][Cart::xxz]+amb1*R[Cart::yyy][Cart::s][Cart::xxz];
R[Cart::xyy][Cart::y][Cart::xxz]=R[Cart::xyyy][Cart::s][Cart::xxz]+amb1*R[Cart::xyy][Cart::s][Cart::xxz];
R[Cart::yyz][Cart::y][Cart::xxz]=R[Cart::yyyz][Cart::s][Cart::xxz]+amb1*R[Cart::yyz][Cart::s][Cart::xxz];
R[Cart::xxy][Cart::y][Cart::xxz]=R[Cart::xxyy][Cart::s][Cart::xxz]+amb1*R[Cart::xxy][Cart::s][Cart::xxz];
R[Cart::xyz][Cart::y][Cart::xxz]=R[Cart::xyyz][Cart::s][Cart::xxz]+amb1*R[Cart::xyz][Cart::s][Cart::xxz];
R[Cart::yzz][Cart::y][Cart::xxz]=R[Cart::yyzz][Cart::s][Cart::xxz]+amb1*R[Cart::yzz][Cart::s][Cart::xxz];
R[Cart::xxx][Cart::y][Cart::xxz]=R[Cart::xxxy][Cart::s][Cart::xxz]+amb1*R[Cart::xxx][Cart::s][Cart::xxz];
R[Cart::xxz][Cart::y][Cart::xxz]=R[Cart::xxyz][Cart::s][Cart::xxz]+amb1*R[Cart::xxz][Cart::s][Cart::xxz];
R[Cart::xzz][Cart::y][Cart::xxz]=R[Cart::xyzz][Cart::s][Cart::xxz]+amb1*R[Cart::xzz][Cart::s][Cart::xxz];
R[Cart::zzz][Cart::y][Cart::xxz]=R[Cart::yzzz][Cart::s][Cart::xxz]+amb1*R[Cart::zzz][Cart::s][Cart::xxz];
R[Cart::yyy][Cart::x][Cart::xxz]=R[Cart::xyyy][Cart::s][Cart::xxz]+amb0*R[Cart::yyy][Cart::s][Cart::xxz];
R[Cart::xyy][Cart::x][Cart::xxz]=R[Cart::xxyy][Cart::s][Cart::xxz]+amb0*R[Cart::xyy][Cart::s][Cart::xxz];
R[Cart::yyz][Cart::x][Cart::xxz]=R[Cart::xyyz][Cart::s][Cart::xxz]+amb0*R[Cart::yyz][Cart::s][Cart::xxz];
R[Cart::xxy][Cart::x][Cart::xxz]=R[Cart::xxxy][Cart::s][Cart::xxz]+amb0*R[Cart::xxy][Cart::s][Cart::xxz];
R[Cart::xyz][Cart::x][Cart::xxz]=R[Cart::xxyz][Cart::s][Cart::xxz]+amb0*R[Cart::xyz][Cart::s][Cart::xxz];
R[Cart::yzz][Cart::x][Cart::xxz]=R[Cart::xyzz][Cart::s][Cart::xxz]+amb0*R[Cart::yzz][Cart::s][Cart::xxz];
R[Cart::xxx][Cart::x][Cart::xxz]=R[Cart::xxxx][Cart::s][Cart::xxz]+amb0*R[Cart::xxx][Cart::s][Cart::xxz];
R[Cart::xxz][Cart::x][Cart::xxz]=R[Cart::xxxz][Cart::s][Cart::xxz]+amb0*R[Cart::xxz][Cart::s][Cart::xxz];
R[Cart::xzz][Cart::x][Cart::xxz]=R[Cart::xxzz][Cart::s][Cart::xxz]+amb0*R[Cart::xzz][Cart::s][Cart::xxz];
R[Cart::zzz][Cart::x][Cart::xxz]=R[Cart::xzzz][Cart::s][Cart::xxz]+amb0*R[Cart::zzz][Cart::s][Cart::xxz];
R[Cart::yyy][Cart::z][Cart::xxz]=R[Cart::yyyz][Cart::s][Cart::xxz]+amb2*R[Cart::yyy][Cart::s][Cart::xxz];
R[Cart::xyy][Cart::z][Cart::xxz]=R[Cart::xyyz][Cart::s][Cart::xxz]+amb2*R[Cart::xyy][Cart::s][Cart::xxz];
R[Cart::yyz][Cart::z][Cart::xxz]=R[Cart::yyzz][Cart::s][Cart::xxz]+amb2*R[Cart::yyz][Cart::s][Cart::xxz];
R[Cart::xxy][Cart::z][Cart::xxz]=R[Cart::xxyz][Cart::s][Cart::xxz]+amb2*R[Cart::xxy][Cart::s][Cart::xxz];
R[Cart::xyz][Cart::z][Cart::xxz]=R[Cart::xyzz][Cart::s][Cart::xxz]+amb2*R[Cart::xyz][Cart::s][Cart::xxz];
R[Cart::yzz][Cart::z][Cart::xxz]=R[Cart::yzzz][Cart::s][Cart::xxz]+amb2*R[Cart::yzz][Cart::s][Cart::xxz];
R[Cart::xxx][Cart::z][Cart::xxz]=R[Cart::xxxz][Cart::s][Cart::xxz]+amb2*R[Cart::xxx][Cart::s][Cart::xxz];
R[Cart::xxz][Cart::z][Cart::xxz]=R[Cart::xxzz][Cart::s][Cart::xxz]+amb2*R[Cart::xxz][Cart::s][Cart::xxz];
R[Cart::xzz][Cart::z][Cart::xxz]=R[Cart::xzzz][Cart::s][Cart::xxz]+amb2*R[Cart::xzz][Cart::s][Cart::xxz];
R[Cart::zzz][Cart::z][Cart::xxz]=R[Cart::zzzz][Cart::s][Cart::xxz]+amb2*R[Cart::zzz][Cart::s][Cart::xxz];
R[Cart::yyy][Cart::y][Cart::xzz]=R[Cart::yyyy][Cart::s][Cart::xzz]+amb1*R[Cart::yyy][Cart::s][Cart::xzz];
R[Cart::xyy][Cart::y][Cart::xzz]=R[Cart::xyyy][Cart::s][Cart::xzz]+amb1*R[Cart::xyy][Cart::s][Cart::xzz];
R[Cart::yyz][Cart::y][Cart::xzz]=R[Cart::yyyz][Cart::s][Cart::xzz]+amb1*R[Cart::yyz][Cart::s][Cart::xzz];
R[Cart::xxy][Cart::y][Cart::xzz]=R[Cart::xxyy][Cart::s][Cart::xzz]+amb1*R[Cart::xxy][Cart::s][Cart::xzz];
R[Cart::xyz][Cart::y][Cart::xzz]=R[Cart::xyyz][Cart::s][Cart::xzz]+amb1*R[Cart::xyz][Cart::s][Cart::xzz];
R[Cart::yzz][Cart::y][Cart::xzz]=R[Cart::yyzz][Cart::s][Cart::xzz]+amb1*R[Cart::yzz][Cart::s][Cart::xzz];
R[Cart::xxx][Cart::y][Cart::xzz]=R[Cart::xxxy][Cart::s][Cart::xzz]+amb1*R[Cart::xxx][Cart::s][Cart::xzz];
R[Cart::xxz][Cart::y][Cart::xzz]=R[Cart::xxyz][Cart::s][Cart::xzz]+amb1*R[Cart::xxz][Cart::s][Cart::xzz];
R[Cart::xzz][Cart::y][Cart::xzz]=R[Cart::xyzz][Cart::s][Cart::xzz]+amb1*R[Cart::xzz][Cart::s][Cart::xzz];
R[Cart::zzz][Cart::y][Cart::xzz]=R[Cart::yzzz][Cart::s][Cart::xzz]+amb1*R[Cart::zzz][Cart::s][Cart::xzz];
R[Cart::yyy][Cart::x][Cart::xzz]=R[Cart::xyyy][Cart::s][Cart::xzz]+amb0*R[Cart::yyy][Cart::s][Cart::xzz];
R[Cart::xyy][Cart::x][Cart::xzz]=R[Cart::xxyy][Cart::s][Cart::xzz]+amb0*R[Cart::xyy][Cart::s][Cart::xzz];
R[Cart::yyz][Cart::x][Cart::xzz]=R[Cart::xyyz][Cart::s][Cart::xzz]+amb0*R[Cart::yyz][Cart::s][Cart::xzz];
R[Cart::xxy][Cart::x][Cart::xzz]=R[Cart::xxxy][Cart::s][Cart::xzz]+amb0*R[Cart::xxy][Cart::s][Cart::xzz];
R[Cart::xyz][Cart::x][Cart::xzz]=R[Cart::xxyz][Cart::s][Cart::xzz]+amb0*R[Cart::xyz][Cart::s][Cart::xzz];
R[Cart::yzz][Cart::x][Cart::xzz]=R[Cart::xyzz][Cart::s][Cart::xzz]+amb0*R[Cart::yzz][Cart::s][Cart::xzz];
R[Cart::xxx][Cart::x][Cart::xzz]=R[Cart::xxxx][Cart::s][Cart::xzz]+amb0*R[Cart::xxx][Cart::s][Cart::xzz];
R[Cart::xxz][Cart::x][Cart::xzz]=R[Cart::xxxz][Cart::s][Cart::xzz]+amb0*R[Cart::xxz][Cart::s][Cart::xzz];
R[Cart::xzz][Cart::x][Cart::xzz]=R[Cart::xxzz][Cart::s][Cart::xzz]+amb0*R[Cart::xzz][Cart::s][Cart::xzz];
R[Cart::zzz][Cart::x][Cart::xzz]=R[Cart::xzzz][Cart::s][Cart::xzz]+amb0*R[Cart::zzz][Cart::s][Cart::xzz];
R[Cart::yyy][Cart::z][Cart::xzz]=R[Cart::yyyz][Cart::s][Cart::xzz]+amb2*R[Cart::yyy][Cart::s][Cart::xzz];
R[Cart::xyy][Cart::z][Cart::xzz]=R[Cart::xyyz][Cart::s][Cart::xzz]+amb2*R[Cart::xyy][Cart::s][Cart::xzz];
R[Cart::yyz][Cart::z][Cart::xzz]=R[Cart::yyzz][Cart::s][Cart::xzz]+amb2*R[Cart::yyz][Cart::s][Cart::xzz];
R[Cart::xxy][Cart::z][Cart::xzz]=R[Cart::xxyz][Cart::s][Cart::xzz]+amb2*R[Cart::xxy][Cart::s][Cart::xzz];
R[Cart::xyz][Cart::z][Cart::xzz]=R[Cart::xyzz][Cart::s][Cart::xzz]+amb2*R[Cart::xyz][Cart::s][Cart::xzz];
R[Cart::yzz][Cart::z][Cart::xzz]=R[Cart::yzzz][Cart::s][Cart::xzz]+amb2*R[Cart::yzz][Cart::s][Cart::xzz];
R[Cart::xxx][Cart::z][Cart::xzz]=R[Cart::xxxz][Cart::s][Cart::xzz]+amb2*R[Cart::xxx][Cart::s][Cart::xzz];
R[Cart::xxz][Cart::z][Cart::xzz]=R[Cart::xxzz][Cart::s][Cart::xzz]+amb2*R[Cart::xxz][Cart::s][Cart::xzz];
R[Cart::xzz][Cart::z][Cart::xzz]=R[Cart::xzzz][Cart::s][Cart::xzz]+amb2*R[Cart::xzz][Cart::s][Cart::xzz];
R[Cart::zzz][Cart::z][Cart::xzz]=R[Cart::zzzz][Cart::s][Cart::xzz]+amb2*R[Cart::zzz][Cart::s][Cart::xzz];
R[Cart::yyy][Cart::y][Cart::zzz]=R[Cart::yyyy][Cart::s][Cart::zzz]+amb1*R[Cart::yyy][Cart::s][Cart::zzz];
R[Cart::xyy][Cart::y][Cart::zzz]=R[Cart::xyyy][Cart::s][Cart::zzz]+amb1*R[Cart::xyy][Cart::s][Cart::zzz];
R[Cart::yyz][Cart::y][Cart::zzz]=R[Cart::yyyz][Cart::s][Cart::zzz]+amb1*R[Cart::yyz][Cart::s][Cart::zzz];
R[Cart::xxy][Cart::y][Cart::zzz]=R[Cart::xxyy][Cart::s][Cart::zzz]+amb1*R[Cart::xxy][Cart::s][Cart::zzz];
R[Cart::xyz][Cart::y][Cart::zzz]=R[Cart::xyyz][Cart::s][Cart::zzz]+amb1*R[Cart::xyz][Cart::s][Cart::zzz];
R[Cart::yzz][Cart::y][Cart::zzz]=R[Cart::yyzz][Cart::s][Cart::zzz]+amb1*R[Cart::yzz][Cart::s][Cart::zzz];
R[Cart::xxx][Cart::y][Cart::zzz]=R[Cart::xxxy][Cart::s][Cart::zzz]+amb1*R[Cart::xxx][Cart::s][Cart::zzz];
R[Cart::xxz][Cart::y][Cart::zzz]=R[Cart::xxyz][Cart::s][Cart::zzz]+amb1*R[Cart::xxz][Cart::s][Cart::zzz];
R[Cart::xzz][Cart::y][Cart::zzz]=R[Cart::xyzz][Cart::s][Cart::zzz]+amb1*R[Cart::xzz][Cart::s][Cart::zzz];
R[Cart::zzz][Cart::y][Cart::zzz]=R[Cart::yzzz][Cart::s][Cart::zzz]+amb1*R[Cart::zzz][Cart::s][Cart::zzz];
R[Cart::yyy][Cart::x][Cart::zzz]=R[Cart::xyyy][Cart::s][Cart::zzz]+amb0*R[Cart::yyy][Cart::s][Cart::zzz];
R[Cart::xyy][Cart::x][Cart::zzz]=R[Cart::xxyy][Cart::s][Cart::zzz]+amb0*R[Cart::xyy][Cart::s][Cart::zzz];
R[Cart::yyz][Cart::x][Cart::zzz]=R[Cart::xyyz][Cart::s][Cart::zzz]+amb0*R[Cart::yyz][Cart::s][Cart::zzz];
R[Cart::xxy][Cart::x][Cart::zzz]=R[Cart::xxxy][Cart::s][Cart::zzz]+amb0*R[Cart::xxy][Cart::s][Cart::zzz];
R[Cart::xyz][Cart::x][Cart::zzz]=R[Cart::xxyz][Cart::s][Cart::zzz]+amb0*R[Cart::xyz][Cart::s][Cart::zzz];
R[Cart::yzz][Cart::x][Cart::zzz]=R[Cart::xyzz][Cart::s][Cart::zzz]+amb0*R[Cart::yzz][Cart::s][Cart::zzz];
R[Cart::xxx][Cart::x][Cart::zzz]=R[Cart::xxxx][Cart::s][Cart::zzz]+amb0*R[Cart::xxx][Cart::s][Cart::zzz];
R[Cart::xxz][Cart::x][Cart::zzz]=R[Cart::xxxz][Cart::s][Cart::zzz]+amb0*R[Cart::xxz][Cart::s][Cart::zzz];
R[Cart::xzz][Cart::x][Cart::zzz]=R[Cart::xxzz][Cart::s][Cart::zzz]+amb0*R[Cart::xzz][Cart::s][Cart::zzz];
R[Cart::zzz][Cart::x][Cart::zzz]=R[Cart::xzzz][Cart::s][Cart::zzz]+amb0*R[Cart::zzz][Cart::s][Cart::zzz];
R[Cart::yyy][Cart::z][Cart::zzz]=R[Cart::yyyz][Cart::s][Cart::zzz]+amb2*R[Cart::yyy][Cart::s][Cart::zzz];
R[Cart::xyy][Cart::z][Cart::zzz]=R[Cart::xyyz][Cart::s][Cart::zzz]+amb2*R[Cart::xyy][Cart::s][Cart::zzz];
R[Cart::yyz][Cart::z][Cart::zzz]=R[Cart::yyzz][Cart::s][Cart::zzz]+amb2*R[Cart::yyz][Cart::s][Cart::zzz];
R[Cart::xxy][Cart::z][Cart::zzz]=R[Cart::xxyz][Cart::s][Cart::zzz]+amb2*R[Cart::xxy][Cart::s][Cart::zzz];
R[Cart::xyz][Cart::z][Cart::zzz]=R[Cart::xyzz][Cart::s][Cart::zzz]+amb2*R[Cart::xyz][Cart::s][Cart::zzz];
R[Cart::yzz][Cart::z][Cart::zzz]=R[Cart::yzzz][Cart::s][Cart::zzz]+amb2*R[Cart::yzz][Cart::s][Cart::zzz];
R[Cart::xxx][Cart::z][Cart::zzz]=R[Cart::xxxz][Cart::s][Cart::zzz]+amb2*R[Cart::xxx][Cart::s][Cart::zzz];
R[Cart::xxz][Cart::z][Cart::zzz]=R[Cart::xxzz][Cart::s][Cart::zzz]+amb2*R[Cart::xxz][Cart::s][Cart::zzz];
R[Cart::xzz][Cart::z][Cart::zzz]=R[Cart::xzzz][Cart::s][Cart::zzz]+amb2*R[Cart::xzz][Cart::s][Cart::zzz];
R[Cart::zzz][Cart::z][Cart::zzz]=R[Cart::zzzz][Cart::s][Cart::zzz]+amb2*R[Cart::zzz][Cart::s][Cart::zzz];
}
//------------------------------------------------------

//Integral s - d - s
if (_lmax_beta>1){

R[Cart::s][Cart::yy][Cart::s]=R[Cart::y][Cart::y][Cart::s]+amb1*R[Cart::s][Cart::y][Cart::s];
R[Cart::s][Cart::xy][Cart::s]=R[Cart::x][Cart::y][Cart::s]+amb0*R[Cart::s][Cart::y][Cart::s];
R[Cart::s][Cart::yz][Cart::s]=R[Cart::y][Cart::z][Cart::s]+amb1*R[Cart::s][Cart::z][Cart::s];
R[Cart::s][Cart::xx][Cart::s]=R[Cart::x][Cart::x][Cart::s]+amb0*R[Cart::s][Cart::x][Cart::s];
R[Cart::s][Cart::xz][Cart::s]=R[Cart::x][Cart::z][Cart::s]+amb0*R[Cart::s][Cart::z][Cart::s];
R[Cart::s][Cart::zz][Cart::s]=R[Cart::z][Cart::z][Cart::s]+amb2*R[Cart::s][Cart::z][Cart::s];
}
//------------------------------------------------------

//Integral p - d - s
if (_lmax_beta>1 && _lmax_alpha>0){

R[Cart::y][Cart::yy][Cart::s]=R[Cart::yy][Cart::y][Cart::s]+amb1*R[Cart::y][Cart::y][Cart::s];
R[Cart::x][Cart::yy][Cart::s]=R[Cart::xy][Cart::y][Cart::s]+amb1*R[Cart::x][Cart::y][Cart::s];
R[Cart::z][Cart::yy][Cart::s]=R[Cart::yz][Cart::y][Cart::s]+amb1*R[Cart::z][Cart::y][Cart::s];
R[Cart::y][Cart::xy][Cart::s]=R[Cart::xy][Cart::y][Cart::s]+amb0*R[Cart::y][Cart::y][Cart::s];
R[Cart::x][Cart::xy][Cart::s]=R[Cart::xx][Cart::y][Cart::s]+amb0*R[Cart::x][Cart::y][Cart::s];
R[Cart::z][Cart::xy][Cart::s]=R[Cart::xz][Cart::y][Cart::s]+amb0*R[Cart::z][Cart::y][Cart::s];
R[Cart::y][Cart::yz][Cart::s]=R[Cart::yy][Cart::z][Cart::s]+amb1*R[Cart::y][Cart::z][Cart::s];
R[Cart::x][Cart::yz][Cart::s]=R[Cart::xy][Cart::z][Cart::s]+amb1*R[Cart::x][Cart::z][Cart::s];
R[Cart::z][Cart::yz][Cart::s]=R[Cart::yz][Cart::z][Cart::s]+amb1*R[Cart::z][Cart::z][Cart::s];
R[Cart::y][Cart::xx][Cart::s]=R[Cart::xy][Cart::x][Cart::s]+amb0*R[Cart::y][Cart::x][Cart::s];
R[Cart::x][Cart::xx][Cart::s]=R[Cart::xx][Cart::x][Cart::s]+amb0*R[Cart::x][Cart::x][Cart::s];
R[Cart::z][Cart::xx][Cart::s]=R[Cart::xz][Cart::x][Cart::s]+amb0*R[Cart::z][Cart::x][Cart::s];
R[Cart::y][Cart::xz][Cart::s]=R[Cart::xy][Cart::z][Cart::s]+amb0*R[Cart::y][Cart::z][Cart::s];
R[Cart::x][Cart::xz][Cart::s]=R[Cart::xx][Cart::z][Cart::s]+amb0*R[Cart::x][Cart::z][Cart::s];
R[Cart::z][Cart::xz][Cart::s]=R[Cart::xz][Cart::z][Cart::s]+amb0*R[Cart::z][Cart::z][Cart::s];
R[Cart::y][Cart::zz][Cart::s]=R[Cart::yz][Cart::z][Cart::s]+amb2*R[Cart::y][Cart::z][Cart::s];
R[Cart::x][Cart::zz][Cart::s]=R[Cart::xz][Cart::z][Cart::s]+amb2*R[Cart::x][Cart::z][Cart::s];
R[Cart::z][Cart::zz][Cart::s]=R[Cart::zz][Cart::z][Cart::s]+amb2*R[Cart::z][Cart::z][Cart::s];
}
//------------------------------------------------------

//Integral d - d - s
if (_lmax_beta>1 && _lmax_alpha>1){

R[Cart::yy][Cart::yy][Cart::s]=R[Cart::yyy][Cart::y][Cart::s]+amb1*R[Cart::yy][Cart::y][Cart::s];
R[Cart::xy][Cart::yy][Cart::s]=R[Cart::xyy][Cart::y][Cart::s]+amb1*R[Cart::xy][Cart::y][Cart::s];
R[Cart::yz][Cart::yy][Cart::s]=R[Cart::yyz][Cart::y][Cart::s]+amb1*R[Cart::yz][Cart::y][Cart::s];
R[Cart::xx][Cart::yy][Cart::s]=R[Cart::xxy][Cart::y][Cart::s]+amb1*R[Cart::xx][Cart::y][Cart::s];
R[Cart::xz][Cart::yy][Cart::s]=R[Cart::xyz][Cart::y][Cart::s]+amb1*R[Cart::xz][Cart::y][Cart::s];
R[Cart::zz][Cart::yy][Cart::s]=R[Cart::yzz][Cart::y][Cart::s]+amb1*R[Cart::zz][Cart::y][Cart::s];
R[Cart::yy][Cart::xy][Cart::s]=R[Cart::xyy][Cart::y][Cart::s]+amb0*R[Cart::yy][Cart::y][Cart::s];
R[Cart::xy][Cart::xy][Cart::s]=R[Cart::xxy][Cart::y][Cart::s]+amb0*R[Cart::xy][Cart::y][Cart::s];
R[Cart::yz][Cart::xy][Cart::s]=R[Cart::xyz][Cart::y][Cart::s]+amb0*R[Cart::yz][Cart::y][Cart::s];
R[Cart::xx][Cart::xy][Cart::s]=R[Cart::xxx][Cart::y][Cart::s]+amb0*R[Cart::xx][Cart::y][Cart::s];
R[Cart::xz][Cart::xy][Cart::s]=R[Cart::xxz][Cart::y][Cart::s]+amb0*R[Cart::xz][Cart::y][Cart::s];
R[Cart::zz][Cart::xy][Cart::s]=R[Cart::xzz][Cart::y][Cart::s]+amb0*R[Cart::zz][Cart::y][Cart::s];
R[Cart::yy][Cart::yz][Cart::s]=R[Cart::yyy][Cart::z][Cart::s]+amb1*R[Cart::yy][Cart::z][Cart::s];
R[Cart::xy][Cart::yz][Cart::s]=R[Cart::xyy][Cart::z][Cart::s]+amb1*R[Cart::xy][Cart::z][Cart::s];
R[Cart::yz][Cart::yz][Cart::s]=R[Cart::yyz][Cart::z][Cart::s]+amb1*R[Cart::yz][Cart::z][Cart::s];
R[Cart::xx][Cart::yz][Cart::s]=R[Cart::xxy][Cart::z][Cart::s]+amb1*R[Cart::xx][Cart::z][Cart::s];
R[Cart::xz][Cart::yz][Cart::s]=R[Cart::xyz][Cart::z][Cart::s]+amb1*R[Cart::xz][Cart::z][Cart::s];
R[Cart::zz][Cart::yz][Cart::s]=R[Cart::yzz][Cart::z][Cart::s]+amb1*R[Cart::zz][Cart::z][Cart::s];
R[Cart::yy][Cart::xx][Cart::s]=R[Cart::xyy][Cart::x][Cart::s]+amb0*R[Cart::yy][Cart::x][Cart::s];
R[Cart::xy][Cart::xx][Cart::s]=R[Cart::xxy][Cart::x][Cart::s]+amb0*R[Cart::xy][Cart::x][Cart::s];
R[Cart::yz][Cart::xx][Cart::s]=R[Cart::xyz][Cart::x][Cart::s]+amb0*R[Cart::yz][Cart::x][Cart::s];
R[Cart::xx][Cart::xx][Cart::s]=R[Cart::xxx][Cart::x][Cart::s]+amb0*R[Cart::xx][Cart::x][Cart::s];
R[Cart::xz][Cart::xx][Cart::s]=R[Cart::xxz][Cart::x][Cart::s]+amb0*R[Cart::xz][Cart::x][Cart::s];
R[Cart::zz][Cart::xx][Cart::s]=R[Cart::xzz][Cart::x][Cart::s]+amb0*R[Cart::zz][Cart::x][Cart::s];
R[Cart::yy][Cart::xz][Cart::s]=R[Cart::xyy][Cart::z][Cart::s]+amb0*R[Cart::yy][Cart::z][Cart::s];
R[Cart::xy][Cart::xz][Cart::s]=R[Cart::xxy][Cart::z][Cart::s]+amb0*R[Cart::xy][Cart::z][Cart::s];
R[Cart::yz][Cart::xz][Cart::s]=R[Cart::xyz][Cart::z][Cart::s]+amb0*R[Cart::yz][Cart::z][Cart::s];
R[Cart::xx][Cart::xz][Cart::s]=R[Cart::xxx][Cart::z][Cart::s]+amb0*R[Cart::xx][Cart::z][Cart::s];
R[Cart::xz][Cart::xz][Cart::s]=R[Cart::xxz][Cart::z][Cart::s]+amb0*R[Cart::xz][Cart::z][Cart::s];
R[Cart::zz][Cart::xz][Cart::s]=R[Cart::xzz][Cart::z][Cart::s]+amb0*R[Cart::zz][Cart::z][Cart::s];
R[Cart::yy][Cart::zz][Cart::s]=R[Cart::yyz][Cart::z][Cart::s]+amb2*R[Cart::yy][Cart::z][Cart::s];
R[Cart::xy][Cart::zz][Cart::s]=R[Cart::xyz][Cart::z][Cart::s]+amb2*R[Cart::xy][Cart::z][Cart::s];
R[Cart::yz][Cart::zz][Cart::s]=R[Cart::yzz][Cart::z][Cart::s]+amb2*R[Cart::yz][Cart::z][Cart::s];
R[Cart::xx][Cart::zz][Cart::s]=R[Cart::xxz][Cart::z][Cart::s]+amb2*R[Cart::xx][Cart::z][Cart::s];
R[Cart::xz][Cart::zz][Cart::s]=R[Cart::xzz][Cart::z][Cart::s]+amb2*R[Cart::xz][Cart::z][Cart::s];
R[Cart::zz][Cart::zz][Cart::s]=R[Cart::zzz][Cart::z][Cart::s]+amb2*R[Cart::zz][Cart::z][Cart::s];
}
//------------------------------------------------------

//Integral s - d - p
if (_lmax_beta>1 && _lmax_gamma>0){

R[Cart::s][Cart::yy][Cart::y]=R[Cart::y][Cart::y][Cart::y]+amb1*R[Cart::s][Cart::y][Cart::y];
R[Cart::s][Cart::xy][Cart::y]=R[Cart::x][Cart::y][Cart::y]+amb0*R[Cart::s][Cart::y][Cart::y];
R[Cart::s][Cart::yz][Cart::y]=R[Cart::y][Cart::z][Cart::y]+amb1*R[Cart::s][Cart::z][Cart::y];
R[Cart::s][Cart::xx][Cart::y]=R[Cart::x][Cart::x][Cart::y]+amb0*R[Cart::s][Cart::x][Cart::y];
R[Cart::s][Cart::xz][Cart::y]=R[Cart::x][Cart::z][Cart::y]+amb0*R[Cart::s][Cart::z][Cart::y];
R[Cart::s][Cart::zz][Cart::y]=R[Cart::z][Cart::z][Cart::y]+amb2*R[Cart::s][Cart::z][Cart::y];
R[Cart::s][Cart::yy][Cart::x]=R[Cart::y][Cart::y][Cart::x]+amb1*R[Cart::s][Cart::y][Cart::x];
R[Cart::s][Cart::xy][Cart::x]=R[Cart::x][Cart::y][Cart::x]+amb0*R[Cart::s][Cart::y][Cart::x];
R[Cart::s][Cart::yz][Cart::x]=R[Cart::y][Cart::z][Cart::x]+amb1*R[Cart::s][Cart::z][Cart::x];
R[Cart::s][Cart::xx][Cart::x]=R[Cart::x][Cart::x][Cart::x]+amb0*R[Cart::s][Cart::x][Cart::x];
R[Cart::s][Cart::xz][Cart::x]=R[Cart::x][Cart::z][Cart::x]+amb0*R[Cart::s][Cart::z][Cart::x];
R[Cart::s][Cart::zz][Cart::x]=R[Cart::z][Cart::z][Cart::x]+amb2*R[Cart::s][Cart::z][Cart::x];
R[Cart::s][Cart::yy][Cart::z]=R[Cart::y][Cart::y][Cart::z]+amb1*R[Cart::s][Cart::y][Cart::z];
R[Cart::s][Cart::xy][Cart::z]=R[Cart::x][Cart::y][Cart::z]+amb0*R[Cart::s][Cart::y][Cart::z];
R[Cart::s][Cart::yz][Cart::z]=R[Cart::y][Cart::z][Cart::z]+amb1*R[Cart::s][Cart::z][Cart::z];
R[Cart::s][Cart::xx][Cart::z]=R[Cart::x][Cart::x][Cart::z]+amb0*R[Cart::s][Cart::x][Cart::z];
R[Cart::s][Cart::xz][Cart::z]=R[Cart::x][Cart::z][Cart::z]+amb0*R[Cart::s][Cart::z][Cart::z];
R[Cart::s][Cart::zz][Cart::z]=R[Cart::z][Cart::z][Cart::z]+amb2*R[Cart::s][Cart::z][Cart::z];
}
//------------------------------------------------------

//Integral p - d - p
if (_lmax_beta>1 && _lmax_alpha>0 && _lmax_gamma>0){

R[Cart::y][Cart::yy][Cart::y]=R[Cart::yy][Cart::y][Cart::y]+amb1*R[Cart::y][Cart::y][Cart::y];
R[Cart::x][Cart::yy][Cart::y]=R[Cart::xy][Cart::y][Cart::y]+amb1*R[Cart::x][Cart::y][Cart::y];
R[Cart::z][Cart::yy][Cart::y]=R[Cart::yz][Cart::y][Cart::y]+amb1*R[Cart::z][Cart::y][Cart::y];
R[Cart::y][Cart::xy][Cart::y]=R[Cart::xy][Cart::y][Cart::y]+amb0*R[Cart::y][Cart::y][Cart::y];
R[Cart::x][Cart::xy][Cart::y]=R[Cart::xx][Cart::y][Cart::y]+amb0*R[Cart::x][Cart::y][Cart::y];
R[Cart::z][Cart::xy][Cart::y]=R[Cart::xz][Cart::y][Cart::y]+amb0*R[Cart::z][Cart::y][Cart::y];
R[Cart::y][Cart::yz][Cart::y]=R[Cart::yy][Cart::z][Cart::y]+amb1*R[Cart::y][Cart::z][Cart::y];
R[Cart::x][Cart::yz][Cart::y]=R[Cart::xy][Cart::z][Cart::y]+amb1*R[Cart::x][Cart::z][Cart::y];
R[Cart::z][Cart::yz][Cart::y]=R[Cart::yz][Cart::z][Cart::y]+amb1*R[Cart::z][Cart::z][Cart::y];
R[Cart::y][Cart::xx][Cart::y]=R[Cart::xy][Cart::x][Cart::y]+amb0*R[Cart::y][Cart::x][Cart::y];
R[Cart::x][Cart::xx][Cart::y]=R[Cart::xx][Cart::x][Cart::y]+amb0*R[Cart::x][Cart::x][Cart::y];
R[Cart::z][Cart::xx][Cart::y]=R[Cart::xz][Cart::x][Cart::y]+amb0*R[Cart::z][Cart::x][Cart::y];
R[Cart::y][Cart::xz][Cart::y]=R[Cart::xy][Cart::z][Cart::y]+amb0*R[Cart::y][Cart::z][Cart::y];
R[Cart::x][Cart::xz][Cart::y]=R[Cart::xx][Cart::z][Cart::y]+amb0*R[Cart::x][Cart::z][Cart::y];
R[Cart::z][Cart::xz][Cart::y]=R[Cart::xz][Cart::z][Cart::y]+amb0*R[Cart::z][Cart::z][Cart::y];
R[Cart::y][Cart::zz][Cart::y]=R[Cart::yz][Cart::z][Cart::y]+amb2*R[Cart::y][Cart::z][Cart::y];
R[Cart::x][Cart::zz][Cart::y]=R[Cart::xz][Cart::z][Cart::y]+amb2*R[Cart::x][Cart::z][Cart::y];
R[Cart::z][Cart::zz][Cart::y]=R[Cart::zz][Cart::z][Cart::y]+amb2*R[Cart::z][Cart::z][Cart::y];
R[Cart::y][Cart::yy][Cart::x]=R[Cart::yy][Cart::y][Cart::x]+amb1*R[Cart::y][Cart::y][Cart::x];
R[Cart::x][Cart::yy][Cart::x]=R[Cart::xy][Cart::y][Cart::x]+amb1*R[Cart::x][Cart::y][Cart::x];
R[Cart::z][Cart::yy][Cart::x]=R[Cart::yz][Cart::y][Cart::x]+amb1*R[Cart::z][Cart::y][Cart::x];
R[Cart::y][Cart::xy][Cart::x]=R[Cart::xy][Cart::y][Cart::x]+amb0*R[Cart::y][Cart::y][Cart::x];
R[Cart::x][Cart::xy][Cart::x]=R[Cart::xx][Cart::y][Cart::x]+amb0*R[Cart::x][Cart::y][Cart::x];
R[Cart::z][Cart::xy][Cart::x]=R[Cart::xz][Cart::y][Cart::x]+amb0*R[Cart::z][Cart::y][Cart::x];
R[Cart::y][Cart::yz][Cart::x]=R[Cart::yy][Cart::z][Cart::x]+amb1*R[Cart::y][Cart::z][Cart::x];
R[Cart::x][Cart::yz][Cart::x]=R[Cart::xy][Cart::z][Cart::x]+amb1*R[Cart::x][Cart::z][Cart::x];
R[Cart::z][Cart::yz][Cart::x]=R[Cart::yz][Cart::z][Cart::x]+amb1*R[Cart::z][Cart::z][Cart::x];
R[Cart::y][Cart::xx][Cart::x]=R[Cart::xy][Cart::x][Cart::x]+amb0*R[Cart::y][Cart::x][Cart::x];
R[Cart::x][Cart::xx][Cart::x]=R[Cart::xx][Cart::x][Cart::x]+amb0*R[Cart::x][Cart::x][Cart::x];
R[Cart::z][Cart::xx][Cart::x]=R[Cart::xz][Cart::x][Cart::x]+amb0*R[Cart::z][Cart::x][Cart::x];
R[Cart::y][Cart::xz][Cart::x]=R[Cart::xy][Cart::z][Cart::x]+amb0*R[Cart::y][Cart::z][Cart::x];
R[Cart::x][Cart::xz][Cart::x]=R[Cart::xx][Cart::z][Cart::x]+amb0*R[Cart::x][Cart::z][Cart::x];
R[Cart::z][Cart::xz][Cart::x]=R[Cart::xz][Cart::z][Cart::x]+amb0*R[Cart::z][Cart::z][Cart::x];
R[Cart::y][Cart::zz][Cart::x]=R[Cart::yz][Cart::z][Cart::x]+amb2*R[Cart::y][Cart::z][Cart::x];
R[Cart::x][Cart::zz][Cart::x]=R[Cart::xz][Cart::z][Cart::x]+amb2*R[Cart::x][Cart::z][Cart::x];
R[Cart::z][Cart::zz][Cart::x]=R[Cart::zz][Cart::z][Cart::x]+amb2*R[Cart::z][Cart::z][Cart::x];
R[Cart::y][Cart::yy][Cart::z]=R[Cart::yy][Cart::y][Cart::z]+amb1*R[Cart::y][Cart::y][Cart::z];
R[Cart::x][Cart::yy][Cart::z]=R[Cart::xy][Cart::y][Cart::z]+amb1*R[Cart::x][Cart::y][Cart::z];
R[Cart::z][Cart::yy][Cart::z]=R[Cart::yz][Cart::y][Cart::z]+amb1*R[Cart::z][Cart::y][Cart::z];
R[Cart::y][Cart::xy][Cart::z]=R[Cart::xy][Cart::y][Cart::z]+amb0*R[Cart::y][Cart::y][Cart::z];
R[Cart::x][Cart::xy][Cart::z]=R[Cart::xx][Cart::y][Cart::z]+amb0*R[Cart::x][Cart::y][Cart::z];
R[Cart::z][Cart::xy][Cart::z]=R[Cart::xz][Cart::y][Cart::z]+amb0*R[Cart::z][Cart::y][Cart::z];
R[Cart::y][Cart::yz][Cart::z]=R[Cart::yy][Cart::z][Cart::z]+amb1*R[Cart::y][Cart::z][Cart::z];
R[Cart::x][Cart::yz][Cart::z]=R[Cart::xy][Cart::z][Cart::z]+amb1*R[Cart::x][Cart::z][Cart::z];
R[Cart::z][Cart::yz][Cart::z]=R[Cart::yz][Cart::z][Cart::z]+amb1*R[Cart::z][Cart::z][Cart::z];
R[Cart::y][Cart::xx][Cart::z]=R[Cart::xy][Cart::x][Cart::z]+amb0*R[Cart::y][Cart::x][Cart::z];
R[Cart::x][Cart::xx][Cart::z]=R[Cart::xx][Cart::x][Cart::z]+amb0*R[Cart::x][Cart::x][Cart::z];
R[Cart::z][Cart::xx][Cart::z]=R[Cart::xz][Cart::x][Cart::z]+amb0*R[Cart::z][Cart::x][Cart::z];
R[Cart::y][Cart::xz][Cart::z]=R[Cart::xy][Cart::z][Cart::z]+amb0*R[Cart::y][Cart::z][Cart::z];
R[Cart::x][Cart::xz][Cart::z]=R[Cart::xx][Cart::z][Cart::z]+amb0*R[Cart::x][Cart::z][Cart::z];
R[Cart::z][Cart::xz][Cart::z]=R[Cart::xz][Cart::z][Cart::z]+amb0*R[Cart::z][Cart::z][Cart::z];
R[Cart::y][Cart::zz][Cart::z]=R[Cart::yz][Cart::z][Cart::z]+amb2*R[Cart::y][Cart::z][Cart::z];
R[Cart::x][Cart::zz][Cart::z]=R[Cart::xz][Cart::z][Cart::z]+amb2*R[Cart::x][Cart::z][Cart::z];
R[Cart::z][Cart::zz][Cart::z]=R[Cart::zz][Cart::z][Cart::z]+amb2*R[Cart::z][Cart::z][Cart::z];
}
//------------------------------------------------------

//Integral d - d - p
if (_lmax_beta>1 && _lmax_alpha>1 && _lmax_gamma>0){

R[Cart::yy][Cart::yy][Cart::y]=R[Cart::yyy][Cart::y][Cart::y]+amb1*R[Cart::yy][Cart::y][Cart::y];
R[Cart::xy][Cart::yy][Cart::y]=R[Cart::xyy][Cart::y][Cart::y]+amb1*R[Cart::xy][Cart::y][Cart::y];
R[Cart::yz][Cart::yy][Cart::y]=R[Cart::yyz][Cart::y][Cart::y]+amb1*R[Cart::yz][Cart::y][Cart::y];
R[Cart::xx][Cart::yy][Cart::y]=R[Cart::xxy][Cart::y][Cart::y]+amb1*R[Cart::xx][Cart::y][Cart::y];
R[Cart::xz][Cart::yy][Cart::y]=R[Cart::xyz][Cart::y][Cart::y]+amb1*R[Cart::xz][Cart::y][Cart::y];
R[Cart::zz][Cart::yy][Cart::y]=R[Cart::yzz][Cart::y][Cart::y]+amb1*R[Cart::zz][Cart::y][Cart::y];
R[Cart::yy][Cart::xy][Cart::y]=R[Cart::xyy][Cart::y][Cart::y]+amb0*R[Cart::yy][Cart::y][Cart::y];
R[Cart::xy][Cart::xy][Cart::y]=R[Cart::xxy][Cart::y][Cart::y]+amb0*R[Cart::xy][Cart::y][Cart::y];
R[Cart::yz][Cart::xy][Cart::y]=R[Cart::xyz][Cart::y][Cart::y]+amb0*R[Cart::yz][Cart::y][Cart::y];
R[Cart::xx][Cart::xy][Cart::y]=R[Cart::xxx][Cart::y][Cart::y]+amb0*R[Cart::xx][Cart::y][Cart::y];
R[Cart::xz][Cart::xy][Cart::y]=R[Cart::xxz][Cart::y][Cart::y]+amb0*R[Cart::xz][Cart::y][Cart::y];
R[Cart::zz][Cart::xy][Cart::y]=R[Cart::xzz][Cart::y][Cart::y]+amb0*R[Cart::zz][Cart::y][Cart::y];
R[Cart::yy][Cart::yz][Cart::y]=R[Cart::yyy][Cart::z][Cart::y]+amb1*R[Cart::yy][Cart::z][Cart::y];
R[Cart::xy][Cart::yz][Cart::y]=R[Cart::xyy][Cart::z][Cart::y]+amb1*R[Cart::xy][Cart::z][Cart::y];
R[Cart::yz][Cart::yz][Cart::y]=R[Cart::yyz][Cart::z][Cart::y]+amb1*R[Cart::yz][Cart::z][Cart::y];
R[Cart::xx][Cart::yz][Cart::y]=R[Cart::xxy][Cart::z][Cart::y]+amb1*R[Cart::xx][Cart::z][Cart::y];
R[Cart::xz][Cart::yz][Cart::y]=R[Cart::xyz][Cart::z][Cart::y]+amb1*R[Cart::xz][Cart::z][Cart::y];
R[Cart::zz][Cart::yz][Cart::y]=R[Cart::yzz][Cart::z][Cart::y]+amb1*R[Cart::zz][Cart::z][Cart::y];
R[Cart::yy][Cart::xx][Cart::y]=R[Cart::xyy][Cart::x][Cart::y]+amb0*R[Cart::yy][Cart::x][Cart::y];
R[Cart::xy][Cart::xx][Cart::y]=R[Cart::xxy][Cart::x][Cart::y]+amb0*R[Cart::xy][Cart::x][Cart::y];
R[Cart::yz][Cart::xx][Cart::y]=R[Cart::xyz][Cart::x][Cart::y]+amb0*R[Cart::yz][Cart::x][Cart::y];
R[Cart::xx][Cart::xx][Cart::y]=R[Cart::xxx][Cart::x][Cart::y]+amb0*R[Cart::xx][Cart::x][Cart::y];
R[Cart::xz][Cart::xx][Cart::y]=R[Cart::xxz][Cart::x][Cart::y]+amb0*R[Cart::xz][Cart::x][Cart::y];
R[Cart::zz][Cart::xx][Cart::y]=R[Cart::xzz][Cart::x][Cart::y]+amb0*R[Cart::zz][Cart::x][Cart::y];
R[Cart::yy][Cart::xz][Cart::y]=R[Cart::xyy][Cart::z][Cart::y]+amb0*R[Cart::yy][Cart::z][Cart::y];
R[Cart::xy][Cart::xz][Cart::y]=R[Cart::xxy][Cart::z][Cart::y]+amb0*R[Cart::xy][Cart::z][Cart::y];
R[Cart::yz][Cart::xz][Cart::y]=R[Cart::xyz][Cart::z][Cart::y]+amb0*R[Cart::yz][Cart::z][Cart::y];
R[Cart::xx][Cart::xz][Cart::y]=R[Cart::xxx][Cart::z][Cart::y]+amb0*R[Cart::xx][Cart::z][Cart::y];
R[Cart::xz][Cart::xz][Cart::y]=R[Cart::xxz][Cart::z][Cart::y]+amb0*R[Cart::xz][Cart::z][Cart::y];
R[Cart::zz][Cart::xz][Cart::y]=R[Cart::xzz][Cart::z][Cart::y]+amb0*R[Cart::zz][Cart::z][Cart::y];
R[Cart::yy][Cart::zz][Cart::y]=R[Cart::yyz][Cart::z][Cart::y]+amb2*R[Cart::yy][Cart::z][Cart::y];
R[Cart::xy][Cart::zz][Cart::y]=R[Cart::xyz][Cart::z][Cart::y]+amb2*R[Cart::xy][Cart::z][Cart::y];
R[Cart::yz][Cart::zz][Cart::y]=R[Cart::yzz][Cart::z][Cart::y]+amb2*R[Cart::yz][Cart::z][Cart::y];
R[Cart::xx][Cart::zz][Cart::y]=R[Cart::xxz][Cart::z][Cart::y]+amb2*R[Cart::xx][Cart::z][Cart::y];
R[Cart::xz][Cart::zz][Cart::y]=R[Cart::xzz][Cart::z][Cart::y]+amb2*R[Cart::xz][Cart::z][Cart::y];
R[Cart::zz][Cart::zz][Cart::y]=R[Cart::zzz][Cart::z][Cart::y]+amb2*R[Cart::zz][Cart::z][Cart::y];
R[Cart::yy][Cart::yy][Cart::x]=R[Cart::yyy][Cart::y][Cart::x]+amb1*R[Cart::yy][Cart::y][Cart::x];
R[Cart::xy][Cart::yy][Cart::x]=R[Cart::xyy][Cart::y][Cart::x]+amb1*R[Cart::xy][Cart::y][Cart::x];
R[Cart::yz][Cart::yy][Cart::x]=R[Cart::yyz][Cart::y][Cart::x]+amb1*R[Cart::yz][Cart::y][Cart::x];
R[Cart::xx][Cart::yy][Cart::x]=R[Cart::xxy][Cart::y][Cart::x]+amb1*R[Cart::xx][Cart::y][Cart::x];
R[Cart::xz][Cart::yy][Cart::x]=R[Cart::xyz][Cart::y][Cart::x]+amb1*R[Cart::xz][Cart::y][Cart::x];
R[Cart::zz][Cart::yy][Cart::x]=R[Cart::yzz][Cart::y][Cart::x]+amb1*R[Cart::zz][Cart::y][Cart::x];
R[Cart::yy][Cart::xy][Cart::x]=R[Cart::xyy][Cart::y][Cart::x]+amb0*R[Cart::yy][Cart::y][Cart::x];
R[Cart::xy][Cart::xy][Cart::x]=R[Cart::xxy][Cart::y][Cart::x]+amb0*R[Cart::xy][Cart::y][Cart::x];
R[Cart::yz][Cart::xy][Cart::x]=R[Cart::xyz][Cart::y][Cart::x]+amb0*R[Cart::yz][Cart::y][Cart::x];
R[Cart::xx][Cart::xy][Cart::x]=R[Cart::xxx][Cart::y][Cart::x]+amb0*R[Cart::xx][Cart::y][Cart::x];
R[Cart::xz][Cart::xy][Cart::x]=R[Cart::xxz][Cart::y][Cart::x]+amb0*R[Cart::xz][Cart::y][Cart::x];
R[Cart::zz][Cart::xy][Cart::x]=R[Cart::xzz][Cart::y][Cart::x]+amb0*R[Cart::zz][Cart::y][Cart::x];
R[Cart::yy][Cart::yz][Cart::x]=R[Cart::yyy][Cart::z][Cart::x]+amb1*R[Cart::yy][Cart::z][Cart::x];
R[Cart::xy][Cart::yz][Cart::x]=R[Cart::xyy][Cart::z][Cart::x]+amb1*R[Cart::xy][Cart::z][Cart::x];
R[Cart::yz][Cart::yz][Cart::x]=R[Cart::yyz][Cart::z][Cart::x]+amb1*R[Cart::yz][Cart::z][Cart::x];
R[Cart::xx][Cart::yz][Cart::x]=R[Cart::xxy][Cart::z][Cart::x]+amb1*R[Cart::xx][Cart::z][Cart::x];
R[Cart::xz][Cart::yz][Cart::x]=R[Cart::xyz][Cart::z][Cart::x]+amb1*R[Cart::xz][Cart::z][Cart::x];
R[Cart::zz][Cart::yz][Cart::x]=R[Cart::yzz][Cart::z][Cart::x]+amb1*R[Cart::zz][Cart::z][Cart::x];
R[Cart::yy][Cart::xx][Cart::x]=R[Cart::xyy][Cart::x][Cart::x]+amb0*R[Cart::yy][Cart::x][Cart::x];
R[Cart::xy][Cart::xx][Cart::x]=R[Cart::xxy][Cart::x][Cart::x]+amb0*R[Cart::xy][Cart::x][Cart::x];
R[Cart::yz][Cart::xx][Cart::x]=R[Cart::xyz][Cart::x][Cart::x]+amb0*R[Cart::yz][Cart::x][Cart::x];
R[Cart::xx][Cart::xx][Cart::x]=R[Cart::xxx][Cart::x][Cart::x]+amb0*R[Cart::xx][Cart::x][Cart::x];
R[Cart::xz][Cart::xx][Cart::x]=R[Cart::xxz][Cart::x][Cart::x]+amb0*R[Cart::xz][Cart::x][Cart::x];
R[Cart::zz][Cart::xx][Cart::x]=R[Cart::xzz][Cart::x][Cart::x]+amb0*R[Cart::zz][Cart::x][Cart::x];
R[Cart::yy][Cart::xz][Cart::x]=R[Cart::xyy][Cart::z][Cart::x]+amb0*R[Cart::yy][Cart::z][Cart::x];
R[Cart::xy][Cart::xz][Cart::x]=R[Cart::xxy][Cart::z][Cart::x]+amb0*R[Cart::xy][Cart::z][Cart::x];
R[Cart::yz][Cart::xz][Cart::x]=R[Cart::xyz][Cart::z][Cart::x]+amb0*R[Cart::yz][Cart::z][Cart::x];
R[Cart::xx][Cart::xz][Cart::x]=R[Cart::xxx][Cart::z][Cart::x]+amb0*R[Cart::xx][Cart::z][Cart::x];
R[Cart::xz][Cart::xz][Cart::x]=R[Cart::xxz][Cart::z][Cart::x]+amb0*R[Cart::xz][Cart::z][Cart::x];
R[Cart::zz][Cart::xz][Cart::x]=R[Cart::xzz][Cart::z][Cart::x]+amb0*R[Cart::zz][Cart::z][Cart::x];
R[Cart::yy][Cart::zz][Cart::x]=R[Cart::yyz][Cart::z][Cart::x]+amb2*R[Cart::yy][Cart::z][Cart::x];
R[Cart::xy][Cart::zz][Cart::x]=R[Cart::xyz][Cart::z][Cart::x]+amb2*R[Cart::xy][Cart::z][Cart::x];
R[Cart::yz][Cart::zz][Cart::x]=R[Cart::yzz][Cart::z][Cart::x]+amb2*R[Cart::yz][Cart::z][Cart::x];
R[Cart::xx][Cart::zz][Cart::x]=R[Cart::xxz][Cart::z][Cart::x]+amb2*R[Cart::xx][Cart::z][Cart::x];
R[Cart::xz][Cart::zz][Cart::x]=R[Cart::xzz][Cart::z][Cart::x]+amb2*R[Cart::xz][Cart::z][Cart::x];
R[Cart::zz][Cart::zz][Cart::x]=R[Cart::zzz][Cart::z][Cart::x]+amb2*R[Cart::zz][Cart::z][Cart::x];
R[Cart::yy][Cart::yy][Cart::z]=R[Cart::yyy][Cart::y][Cart::z]+amb1*R[Cart::yy][Cart::y][Cart::z];
R[Cart::xy][Cart::yy][Cart::z]=R[Cart::xyy][Cart::y][Cart::z]+amb1*R[Cart::xy][Cart::y][Cart::z];
R[Cart::yz][Cart::yy][Cart::z]=R[Cart::yyz][Cart::y][Cart::z]+amb1*R[Cart::yz][Cart::y][Cart::z];
R[Cart::xx][Cart::yy][Cart::z]=R[Cart::xxy][Cart::y][Cart::z]+amb1*R[Cart::xx][Cart::y][Cart::z];
R[Cart::xz][Cart::yy][Cart::z]=R[Cart::xyz][Cart::y][Cart::z]+amb1*R[Cart::xz][Cart::y][Cart::z];
R[Cart::zz][Cart::yy][Cart::z]=R[Cart::yzz][Cart::y][Cart::z]+amb1*R[Cart::zz][Cart::y][Cart::z];
R[Cart::yy][Cart::xy][Cart::z]=R[Cart::xyy][Cart::y][Cart::z]+amb0*R[Cart::yy][Cart::y][Cart::z];
R[Cart::xy][Cart::xy][Cart::z]=R[Cart::xxy][Cart::y][Cart::z]+amb0*R[Cart::xy][Cart::y][Cart::z];
R[Cart::yz][Cart::xy][Cart::z]=R[Cart::xyz][Cart::y][Cart::z]+amb0*R[Cart::yz][Cart::y][Cart::z];
R[Cart::xx][Cart::xy][Cart::z]=R[Cart::xxx][Cart::y][Cart::z]+amb0*R[Cart::xx][Cart::y][Cart::z];
R[Cart::xz][Cart::xy][Cart::z]=R[Cart::xxz][Cart::y][Cart::z]+amb0*R[Cart::xz][Cart::y][Cart::z];
R[Cart::zz][Cart::xy][Cart::z]=R[Cart::xzz][Cart::y][Cart::z]+amb0*R[Cart::zz][Cart::y][Cart::z];
R[Cart::yy][Cart::yz][Cart::z]=R[Cart::yyy][Cart::z][Cart::z]+amb1*R[Cart::yy][Cart::z][Cart::z];
R[Cart::xy][Cart::yz][Cart::z]=R[Cart::xyy][Cart::z][Cart::z]+amb1*R[Cart::xy][Cart::z][Cart::z];
R[Cart::yz][Cart::yz][Cart::z]=R[Cart::yyz][Cart::z][Cart::z]+amb1*R[Cart::yz][Cart::z][Cart::z];
R[Cart::xx][Cart::yz][Cart::z]=R[Cart::xxy][Cart::z][Cart::z]+amb1*R[Cart::xx][Cart::z][Cart::z];
R[Cart::xz][Cart::yz][Cart::z]=R[Cart::xyz][Cart::z][Cart::z]+amb1*R[Cart::xz][Cart::z][Cart::z];
R[Cart::zz][Cart::yz][Cart::z]=R[Cart::yzz][Cart::z][Cart::z]+amb1*R[Cart::zz][Cart::z][Cart::z];
R[Cart::yy][Cart::xx][Cart::z]=R[Cart::xyy][Cart::x][Cart::z]+amb0*R[Cart::yy][Cart::x][Cart::z];
R[Cart::xy][Cart::xx][Cart::z]=R[Cart::xxy][Cart::x][Cart::z]+amb0*R[Cart::xy][Cart::x][Cart::z];
R[Cart::yz][Cart::xx][Cart::z]=R[Cart::xyz][Cart::x][Cart::z]+amb0*R[Cart::yz][Cart::x][Cart::z];
R[Cart::xx][Cart::xx][Cart::z]=R[Cart::xxx][Cart::x][Cart::z]+amb0*R[Cart::xx][Cart::x][Cart::z];
R[Cart::xz][Cart::xx][Cart::z]=R[Cart::xxz][Cart::x][Cart::z]+amb0*R[Cart::xz][Cart::x][Cart::z];
R[Cart::zz][Cart::xx][Cart::z]=R[Cart::xzz][Cart::x][Cart::z]+amb0*R[Cart::zz][Cart::x][Cart::z];
R[Cart::yy][Cart::xz][Cart::z]=R[Cart::xyy][Cart::z][Cart::z]+amb0*R[Cart::yy][Cart::z][Cart::z];
R[Cart::xy][Cart::xz][Cart::z]=R[Cart::xxy][Cart::z][Cart::z]+amb0*R[Cart::xy][Cart::z][Cart::z];
R[Cart::yz][Cart::xz][Cart::z]=R[Cart::xyz][Cart::z][Cart::z]+amb0*R[Cart::yz][Cart::z][Cart::z];
R[Cart::xx][Cart::xz][Cart::z]=R[Cart::xxx][Cart::z][Cart::z]+amb0*R[Cart::xx][Cart::z][Cart::z];
R[Cart::xz][Cart::xz][Cart::z]=R[Cart::xxz][Cart::z][Cart::z]+amb0*R[Cart::xz][Cart::z][Cart::z];
R[Cart::zz][Cart::xz][Cart::z]=R[Cart::xzz][Cart::z][Cart::z]+amb0*R[Cart::zz][Cart::z][Cart::z];
R[Cart::yy][Cart::zz][Cart::z]=R[Cart::yyz][Cart::z][Cart::z]+amb2*R[Cart::yy][Cart::z][Cart::z];
R[Cart::xy][Cart::zz][Cart::z]=R[Cart::xyz][Cart::z][Cart::z]+amb2*R[Cart::xy][Cart::z][Cart::z];
R[Cart::yz][Cart::zz][Cart::z]=R[Cart::yzz][Cart::z][Cart::z]+amb2*R[Cart::yz][Cart::z][Cart::z];
R[Cart::xx][Cart::zz][Cart::z]=R[Cart::xxz][Cart::z][Cart::z]+amb2*R[Cart::xx][Cart::z][Cart::z];
R[Cart::xz][Cart::zz][Cart::z]=R[Cart::xzz][Cart::z][Cart::z]+amb2*R[Cart::xz][Cart::z][Cart::z];
R[Cart::zz][Cart::zz][Cart::z]=R[Cart::zzz][Cart::z][Cart::z]+amb2*R[Cart::zz][Cart::z][Cart::z];
}
//------------------------------------------------------

//Integral s - d - d
if (_lmax_beta>1 && _lmax_gamma>1){

R[Cart::s][Cart::yy][Cart::yy]=R[Cart::y][Cart::y][Cart::yy]+amb1*R[Cart::s][Cart::y][Cart::yy];
R[Cart::s][Cart::xy][Cart::yy]=R[Cart::x][Cart::y][Cart::yy]+amb0*R[Cart::s][Cart::y][Cart::yy];
R[Cart::s][Cart::yz][Cart::yy]=R[Cart::y][Cart::z][Cart::yy]+amb1*R[Cart::s][Cart::z][Cart::yy];
R[Cart::s][Cart::xx][Cart::yy]=R[Cart::x][Cart::x][Cart::yy]+amb0*R[Cart::s][Cart::x][Cart::yy];
R[Cart::s][Cart::xz][Cart::yy]=R[Cart::x][Cart::z][Cart::yy]+amb0*R[Cart::s][Cart::z][Cart::yy];
R[Cart::s][Cart::zz][Cart::yy]=R[Cart::z][Cart::z][Cart::yy]+amb2*R[Cart::s][Cart::z][Cart::yy];
R[Cart::s][Cart::yy][Cart::xy]=R[Cart::y][Cart::y][Cart::xy]+amb1*R[Cart::s][Cart::y][Cart::xy];
R[Cart::s][Cart::xy][Cart::xy]=R[Cart::x][Cart::y][Cart::xy]+amb0*R[Cart::s][Cart::y][Cart::xy];
R[Cart::s][Cart::yz][Cart::xy]=R[Cart::y][Cart::z][Cart::xy]+amb1*R[Cart::s][Cart::z][Cart::xy];
R[Cart::s][Cart::xx][Cart::xy]=R[Cart::x][Cart::x][Cart::xy]+amb0*R[Cart::s][Cart::x][Cart::xy];
R[Cart::s][Cart::xz][Cart::xy]=R[Cart::x][Cart::z][Cart::xy]+amb0*R[Cart::s][Cart::z][Cart::xy];
R[Cart::s][Cart::zz][Cart::xy]=R[Cart::z][Cart::z][Cart::xy]+amb2*R[Cart::s][Cart::z][Cart::xy];
R[Cart::s][Cart::yy][Cart::yz]=R[Cart::y][Cart::y][Cart::yz]+amb1*R[Cart::s][Cart::y][Cart::yz];
R[Cart::s][Cart::xy][Cart::yz]=R[Cart::x][Cart::y][Cart::yz]+amb0*R[Cart::s][Cart::y][Cart::yz];
R[Cart::s][Cart::yz][Cart::yz]=R[Cart::y][Cart::z][Cart::yz]+amb1*R[Cart::s][Cart::z][Cart::yz];
R[Cart::s][Cart::xx][Cart::yz]=R[Cart::x][Cart::x][Cart::yz]+amb0*R[Cart::s][Cart::x][Cart::yz];
R[Cart::s][Cart::xz][Cart::yz]=R[Cart::x][Cart::z][Cart::yz]+amb0*R[Cart::s][Cart::z][Cart::yz];
R[Cart::s][Cart::zz][Cart::yz]=R[Cart::z][Cart::z][Cart::yz]+amb2*R[Cart::s][Cart::z][Cart::yz];
R[Cart::s][Cart::yy][Cart::xx]=R[Cart::y][Cart::y][Cart::xx]+amb1*R[Cart::s][Cart::y][Cart::xx];
R[Cart::s][Cart::xy][Cart::xx]=R[Cart::x][Cart::y][Cart::xx]+amb0*R[Cart::s][Cart::y][Cart::xx];
R[Cart::s][Cart::yz][Cart::xx]=R[Cart::y][Cart::z][Cart::xx]+amb1*R[Cart::s][Cart::z][Cart::xx];
R[Cart::s][Cart::xx][Cart::xx]=R[Cart::x][Cart::x][Cart::xx]+amb0*R[Cart::s][Cart::x][Cart::xx];
R[Cart::s][Cart::xz][Cart::xx]=R[Cart::x][Cart::z][Cart::xx]+amb0*R[Cart::s][Cart::z][Cart::xx];
R[Cart::s][Cart::zz][Cart::xx]=R[Cart::z][Cart::z][Cart::xx]+amb2*R[Cart::s][Cart::z][Cart::xx];
R[Cart::s][Cart::yy][Cart::xz]=R[Cart::y][Cart::y][Cart::xz]+amb1*R[Cart::s][Cart::y][Cart::xz];
R[Cart::s][Cart::xy][Cart::xz]=R[Cart::x][Cart::y][Cart::xz]+amb0*R[Cart::s][Cart::y][Cart::xz];
R[Cart::s][Cart::yz][Cart::xz]=R[Cart::y][Cart::z][Cart::xz]+amb1*R[Cart::s][Cart::z][Cart::xz];
R[Cart::s][Cart::xx][Cart::xz]=R[Cart::x][Cart::x][Cart::xz]+amb0*R[Cart::s][Cart::x][Cart::xz];
R[Cart::s][Cart::xz][Cart::xz]=R[Cart::x][Cart::z][Cart::xz]+amb0*R[Cart::s][Cart::z][Cart::xz];
R[Cart::s][Cart::zz][Cart::xz]=R[Cart::z][Cart::z][Cart::xz]+amb2*R[Cart::s][Cart::z][Cart::xz];
R[Cart::s][Cart::yy][Cart::zz]=R[Cart::y][Cart::y][Cart::zz]+amb1*R[Cart::s][Cart::y][Cart::zz];
R[Cart::s][Cart::xy][Cart::zz]=R[Cart::x][Cart::y][Cart::zz]+amb0*R[Cart::s][Cart::y][Cart::zz];
R[Cart::s][Cart::yz][Cart::zz]=R[Cart::y][Cart::z][Cart::zz]+amb1*R[Cart::s][Cart::z][Cart::zz];
R[Cart::s][Cart::xx][Cart::zz]=R[Cart::x][Cart::x][Cart::zz]+amb0*R[Cart::s][Cart::x][Cart::zz];
R[Cart::s][Cart::xz][Cart::zz]=R[Cart::x][Cart::z][Cart::zz]+amb0*R[Cart::s][Cart::z][Cart::zz];
R[Cart::s][Cart::zz][Cart::zz]=R[Cart::z][Cart::z][Cart::zz]+amb2*R[Cart::s][Cart::z][Cart::zz];
}
//------------------------------------------------------

//Integral p - d - d
if (_lmax_beta>1 && _lmax_alpha>0 && _lmax_gamma>1){

R[Cart::y][Cart::yy][Cart::yy]=R[Cart::yy][Cart::y][Cart::yy]+amb1*R[Cart::y][Cart::y][Cart::yy];
R[Cart::x][Cart::yy][Cart::yy]=R[Cart::xy][Cart::y][Cart::yy]+amb1*R[Cart::x][Cart::y][Cart::yy];
R[Cart::z][Cart::yy][Cart::yy]=R[Cart::yz][Cart::y][Cart::yy]+amb1*R[Cart::z][Cart::y][Cart::yy];
R[Cart::y][Cart::xy][Cart::yy]=R[Cart::xy][Cart::y][Cart::yy]+amb0*R[Cart::y][Cart::y][Cart::yy];
R[Cart::x][Cart::xy][Cart::yy]=R[Cart::xx][Cart::y][Cart::yy]+amb0*R[Cart::x][Cart::y][Cart::yy];
R[Cart::z][Cart::xy][Cart::yy]=R[Cart::xz][Cart::y][Cart::yy]+amb0*R[Cart::z][Cart::y][Cart::yy];
R[Cart::y][Cart::yz][Cart::yy]=R[Cart::yy][Cart::z][Cart::yy]+amb1*R[Cart::y][Cart::z][Cart::yy];
R[Cart::x][Cart::yz][Cart::yy]=R[Cart::xy][Cart::z][Cart::yy]+amb1*R[Cart::x][Cart::z][Cart::yy];
R[Cart::z][Cart::yz][Cart::yy]=R[Cart::yz][Cart::z][Cart::yy]+amb1*R[Cart::z][Cart::z][Cart::yy];
R[Cart::y][Cart::xx][Cart::yy]=R[Cart::xy][Cart::x][Cart::yy]+amb0*R[Cart::y][Cart::x][Cart::yy];
R[Cart::x][Cart::xx][Cart::yy]=R[Cart::xx][Cart::x][Cart::yy]+amb0*R[Cart::x][Cart::x][Cart::yy];
R[Cart::z][Cart::xx][Cart::yy]=R[Cart::xz][Cart::x][Cart::yy]+amb0*R[Cart::z][Cart::x][Cart::yy];
R[Cart::y][Cart::xz][Cart::yy]=R[Cart::xy][Cart::z][Cart::yy]+amb0*R[Cart::y][Cart::z][Cart::yy];
R[Cart::x][Cart::xz][Cart::yy]=R[Cart::xx][Cart::z][Cart::yy]+amb0*R[Cart::x][Cart::z][Cart::yy];
R[Cart::z][Cart::xz][Cart::yy]=R[Cart::xz][Cart::z][Cart::yy]+amb0*R[Cart::z][Cart::z][Cart::yy];
R[Cart::y][Cart::zz][Cart::yy]=R[Cart::yz][Cart::z][Cart::yy]+amb2*R[Cart::y][Cart::z][Cart::yy];
R[Cart::x][Cart::zz][Cart::yy]=R[Cart::xz][Cart::z][Cart::yy]+amb2*R[Cart::x][Cart::z][Cart::yy];
R[Cart::z][Cart::zz][Cart::yy]=R[Cart::zz][Cart::z][Cart::yy]+amb2*R[Cart::z][Cart::z][Cart::yy];
R[Cart::y][Cart::yy][Cart::xy]=R[Cart::yy][Cart::y][Cart::xy]+amb1*R[Cart::y][Cart::y][Cart::xy];
R[Cart::x][Cart::yy][Cart::xy]=R[Cart::xy][Cart::y][Cart::xy]+amb1*R[Cart::x][Cart::y][Cart::xy];
R[Cart::z][Cart::yy][Cart::xy]=R[Cart::yz][Cart::y][Cart::xy]+amb1*R[Cart::z][Cart::y][Cart::xy];
R[Cart::y][Cart::xy][Cart::xy]=R[Cart::xy][Cart::y][Cart::xy]+amb0*R[Cart::y][Cart::y][Cart::xy];
R[Cart::x][Cart::xy][Cart::xy]=R[Cart::xx][Cart::y][Cart::xy]+amb0*R[Cart::x][Cart::y][Cart::xy];
R[Cart::z][Cart::xy][Cart::xy]=R[Cart::xz][Cart::y][Cart::xy]+amb0*R[Cart::z][Cart::y][Cart::xy];
R[Cart::y][Cart::yz][Cart::xy]=R[Cart::yy][Cart::z][Cart::xy]+amb1*R[Cart::y][Cart::z][Cart::xy];
R[Cart::x][Cart::yz][Cart::xy]=R[Cart::xy][Cart::z][Cart::xy]+amb1*R[Cart::x][Cart::z][Cart::xy];
R[Cart::z][Cart::yz][Cart::xy]=R[Cart::yz][Cart::z][Cart::xy]+amb1*R[Cart::z][Cart::z][Cart::xy];
R[Cart::y][Cart::xx][Cart::xy]=R[Cart::xy][Cart::x][Cart::xy]+amb0*R[Cart::y][Cart::x][Cart::xy];
R[Cart::x][Cart::xx][Cart::xy]=R[Cart::xx][Cart::x][Cart::xy]+amb0*R[Cart::x][Cart::x][Cart::xy];
R[Cart::z][Cart::xx][Cart::xy]=R[Cart::xz][Cart::x][Cart::xy]+amb0*R[Cart::z][Cart::x][Cart::xy];
R[Cart::y][Cart::xz][Cart::xy]=R[Cart::xy][Cart::z][Cart::xy]+amb0*R[Cart::y][Cart::z][Cart::xy];
R[Cart::x][Cart::xz][Cart::xy]=R[Cart::xx][Cart::z][Cart::xy]+amb0*R[Cart::x][Cart::z][Cart::xy];
R[Cart::z][Cart::xz][Cart::xy]=R[Cart::xz][Cart::z][Cart::xy]+amb0*R[Cart::z][Cart::z][Cart::xy];
R[Cart::y][Cart::zz][Cart::xy]=R[Cart::yz][Cart::z][Cart::xy]+amb2*R[Cart::y][Cart::z][Cart::xy];
R[Cart::x][Cart::zz][Cart::xy]=R[Cart::xz][Cart::z][Cart::xy]+amb2*R[Cart::x][Cart::z][Cart::xy];
R[Cart::z][Cart::zz][Cart::xy]=R[Cart::zz][Cart::z][Cart::xy]+amb2*R[Cart::z][Cart::z][Cart::xy];
R[Cart::y][Cart::yy][Cart::yz]=R[Cart::yy][Cart::y][Cart::yz]+amb1*R[Cart::y][Cart::y][Cart::yz];
R[Cart::x][Cart::yy][Cart::yz]=R[Cart::xy][Cart::y][Cart::yz]+amb1*R[Cart::x][Cart::y][Cart::yz];
R[Cart::z][Cart::yy][Cart::yz]=R[Cart::yz][Cart::y][Cart::yz]+amb1*R[Cart::z][Cart::y][Cart::yz];
R[Cart::y][Cart::xy][Cart::yz]=R[Cart::xy][Cart::y][Cart::yz]+amb0*R[Cart::y][Cart::y][Cart::yz];
R[Cart::x][Cart::xy][Cart::yz]=R[Cart::xx][Cart::y][Cart::yz]+amb0*R[Cart::x][Cart::y][Cart::yz];
R[Cart::z][Cart::xy][Cart::yz]=R[Cart::xz][Cart::y][Cart::yz]+amb0*R[Cart::z][Cart::y][Cart::yz];
R[Cart::y][Cart::yz][Cart::yz]=R[Cart::yy][Cart::z][Cart::yz]+amb1*R[Cart::y][Cart::z][Cart::yz];
R[Cart::x][Cart::yz][Cart::yz]=R[Cart::xy][Cart::z][Cart::yz]+amb1*R[Cart::x][Cart::z][Cart::yz];
R[Cart::z][Cart::yz][Cart::yz]=R[Cart::yz][Cart::z][Cart::yz]+amb1*R[Cart::z][Cart::z][Cart::yz];
R[Cart::y][Cart::xx][Cart::yz]=R[Cart::xy][Cart::x][Cart::yz]+amb0*R[Cart::y][Cart::x][Cart::yz];
R[Cart::x][Cart::xx][Cart::yz]=R[Cart::xx][Cart::x][Cart::yz]+amb0*R[Cart::x][Cart::x][Cart::yz];
R[Cart::z][Cart::xx][Cart::yz]=R[Cart::xz][Cart::x][Cart::yz]+amb0*R[Cart::z][Cart::x][Cart::yz];
R[Cart::y][Cart::xz][Cart::yz]=R[Cart::xy][Cart::z][Cart::yz]+amb0*R[Cart::y][Cart::z][Cart::yz];
R[Cart::x][Cart::xz][Cart::yz]=R[Cart::xx][Cart::z][Cart::yz]+amb0*R[Cart::x][Cart::z][Cart::yz];
R[Cart::z][Cart::xz][Cart::yz]=R[Cart::xz][Cart::z][Cart::yz]+amb0*R[Cart::z][Cart::z][Cart::yz];
R[Cart::y][Cart::zz][Cart::yz]=R[Cart::yz][Cart::z][Cart::yz]+amb2*R[Cart::y][Cart::z][Cart::yz];
R[Cart::x][Cart::zz][Cart::yz]=R[Cart::xz][Cart::z][Cart::yz]+amb2*R[Cart::x][Cart::z][Cart::yz];
R[Cart::z][Cart::zz][Cart::yz]=R[Cart::zz][Cart::z][Cart::yz]+amb2*R[Cart::z][Cart::z][Cart::yz];
R[Cart::y][Cart::yy][Cart::xx]=R[Cart::yy][Cart::y][Cart::xx]+amb1*R[Cart::y][Cart::y][Cart::xx];
R[Cart::x][Cart::yy][Cart::xx]=R[Cart::xy][Cart::y][Cart::xx]+amb1*R[Cart::x][Cart::y][Cart::xx];
R[Cart::z][Cart::yy][Cart::xx]=R[Cart::yz][Cart::y][Cart::xx]+amb1*R[Cart::z][Cart::y][Cart::xx];
R[Cart::y][Cart::xy][Cart::xx]=R[Cart::xy][Cart::y][Cart::xx]+amb0*R[Cart::y][Cart::y][Cart::xx];
R[Cart::x][Cart::xy][Cart::xx]=R[Cart::xx][Cart::y][Cart::xx]+amb0*R[Cart::x][Cart::y][Cart::xx];
R[Cart::z][Cart::xy][Cart::xx]=R[Cart::xz][Cart::y][Cart::xx]+amb0*R[Cart::z][Cart::y][Cart::xx];
R[Cart::y][Cart::yz][Cart::xx]=R[Cart::yy][Cart::z][Cart::xx]+amb1*R[Cart::y][Cart::z][Cart::xx];
R[Cart::x][Cart::yz][Cart::xx]=R[Cart::xy][Cart::z][Cart::xx]+amb1*R[Cart::x][Cart::z][Cart::xx];
R[Cart::z][Cart::yz][Cart::xx]=R[Cart::yz][Cart::z][Cart::xx]+amb1*R[Cart::z][Cart::z][Cart::xx];
R[Cart::y][Cart::xx][Cart::xx]=R[Cart::xy][Cart::x][Cart::xx]+amb0*R[Cart::y][Cart::x][Cart::xx];
R[Cart::x][Cart::xx][Cart::xx]=R[Cart::xx][Cart::x][Cart::xx]+amb0*R[Cart::x][Cart::x][Cart::xx];
R[Cart::z][Cart::xx][Cart::xx]=R[Cart::xz][Cart::x][Cart::xx]+amb0*R[Cart::z][Cart::x][Cart::xx];
R[Cart::y][Cart::xz][Cart::xx]=R[Cart::xy][Cart::z][Cart::xx]+amb0*R[Cart::y][Cart::z][Cart::xx];
R[Cart::x][Cart::xz][Cart::xx]=R[Cart::xx][Cart::z][Cart::xx]+amb0*R[Cart::x][Cart::z][Cart::xx];
R[Cart::z][Cart::xz][Cart::xx]=R[Cart::xz][Cart::z][Cart::xx]+amb0*R[Cart::z][Cart::z][Cart::xx];
R[Cart::y][Cart::zz][Cart::xx]=R[Cart::yz][Cart::z][Cart::xx]+amb2*R[Cart::y][Cart::z][Cart::xx];
R[Cart::x][Cart::zz][Cart::xx]=R[Cart::xz][Cart::z][Cart::xx]+amb2*R[Cart::x][Cart::z][Cart::xx];
R[Cart::z][Cart::zz][Cart::xx]=R[Cart::zz][Cart::z][Cart::xx]+amb2*R[Cart::z][Cart::z][Cart::xx];
R[Cart::y][Cart::yy][Cart::xz]=R[Cart::yy][Cart::y][Cart::xz]+amb1*R[Cart::y][Cart::y][Cart::xz];
R[Cart::x][Cart::yy][Cart::xz]=R[Cart::xy][Cart::y][Cart::xz]+amb1*R[Cart::x][Cart::y][Cart::xz];
R[Cart::z][Cart::yy][Cart::xz]=R[Cart::yz][Cart::y][Cart::xz]+amb1*R[Cart::z][Cart::y][Cart::xz];
R[Cart::y][Cart::xy][Cart::xz]=R[Cart::xy][Cart::y][Cart::xz]+amb0*R[Cart::y][Cart::y][Cart::xz];
R[Cart::x][Cart::xy][Cart::xz]=R[Cart::xx][Cart::y][Cart::xz]+amb0*R[Cart::x][Cart::y][Cart::xz];
R[Cart::z][Cart::xy][Cart::xz]=R[Cart::xz][Cart::y][Cart::xz]+amb0*R[Cart::z][Cart::y][Cart::xz];
R[Cart::y][Cart::yz][Cart::xz]=R[Cart::yy][Cart::z][Cart::xz]+amb1*R[Cart::y][Cart::z][Cart::xz];
R[Cart::x][Cart::yz][Cart::xz]=R[Cart::xy][Cart::z][Cart::xz]+amb1*R[Cart::x][Cart::z][Cart::xz];
R[Cart::z][Cart::yz][Cart::xz]=R[Cart::yz][Cart::z][Cart::xz]+amb1*R[Cart::z][Cart::z][Cart::xz];
R[Cart::y][Cart::xx][Cart::xz]=R[Cart::xy][Cart::x][Cart::xz]+amb0*R[Cart::y][Cart::x][Cart::xz];
R[Cart::x][Cart::xx][Cart::xz]=R[Cart::xx][Cart::x][Cart::xz]+amb0*R[Cart::x][Cart::x][Cart::xz];
R[Cart::z][Cart::xx][Cart::xz]=R[Cart::xz][Cart::x][Cart::xz]+amb0*R[Cart::z][Cart::x][Cart::xz];
R[Cart::y][Cart::xz][Cart::xz]=R[Cart::xy][Cart::z][Cart::xz]+amb0*R[Cart::y][Cart::z][Cart::xz];
R[Cart::x][Cart::xz][Cart::xz]=R[Cart::xx][Cart::z][Cart::xz]+amb0*R[Cart::x][Cart::z][Cart::xz];
R[Cart::z][Cart::xz][Cart::xz]=R[Cart::xz][Cart::z][Cart::xz]+amb0*R[Cart::z][Cart::z][Cart::xz];
R[Cart::y][Cart::zz][Cart::xz]=R[Cart::yz][Cart::z][Cart::xz]+amb2*R[Cart::y][Cart::z][Cart::xz];
R[Cart::x][Cart::zz][Cart::xz]=R[Cart::xz][Cart::z][Cart::xz]+amb2*R[Cart::x][Cart::z][Cart::xz];
R[Cart::z][Cart::zz][Cart::xz]=R[Cart::zz][Cart::z][Cart::xz]+amb2*R[Cart::z][Cart::z][Cart::xz];
R[Cart::y][Cart::yy][Cart::zz]=R[Cart::yy][Cart::y][Cart::zz]+amb1*R[Cart::y][Cart::y][Cart::zz];
R[Cart::x][Cart::yy][Cart::zz]=R[Cart::xy][Cart::y][Cart::zz]+amb1*R[Cart::x][Cart::y][Cart::zz];
R[Cart::z][Cart::yy][Cart::zz]=R[Cart::yz][Cart::y][Cart::zz]+amb1*R[Cart::z][Cart::y][Cart::zz];
R[Cart::y][Cart::xy][Cart::zz]=R[Cart::xy][Cart::y][Cart::zz]+amb0*R[Cart::y][Cart::y][Cart::zz];
R[Cart::x][Cart::xy][Cart::zz]=R[Cart::xx][Cart::y][Cart::zz]+amb0*R[Cart::x][Cart::y][Cart::zz];
R[Cart::z][Cart::xy][Cart::zz]=R[Cart::xz][Cart::y][Cart::zz]+amb0*R[Cart::z][Cart::y][Cart::zz];
R[Cart::y][Cart::yz][Cart::zz]=R[Cart::yy][Cart::z][Cart::zz]+amb1*R[Cart::y][Cart::z][Cart::zz];
R[Cart::x][Cart::yz][Cart::zz]=R[Cart::xy][Cart::z][Cart::zz]+amb1*R[Cart::x][Cart::z][Cart::zz];
R[Cart::z][Cart::yz][Cart::zz]=R[Cart::yz][Cart::z][Cart::zz]+amb1*R[Cart::z][Cart::z][Cart::zz];
R[Cart::y][Cart::xx][Cart::zz]=R[Cart::xy][Cart::x][Cart::zz]+amb0*R[Cart::y][Cart::x][Cart::zz];
R[Cart::x][Cart::xx][Cart::zz]=R[Cart::xx][Cart::x][Cart::zz]+amb0*R[Cart::x][Cart::x][Cart::zz];
R[Cart::z][Cart::xx][Cart::zz]=R[Cart::xz][Cart::x][Cart::zz]+amb0*R[Cart::z][Cart::x][Cart::zz];
R[Cart::y][Cart::xz][Cart::zz]=R[Cart::xy][Cart::z][Cart::zz]+amb0*R[Cart::y][Cart::z][Cart::zz];
R[Cart::x][Cart::xz][Cart::zz]=R[Cart::xx][Cart::z][Cart::zz]+amb0*R[Cart::x][Cart::z][Cart::zz];
R[Cart::z][Cart::xz][Cart::zz]=R[Cart::xz][Cart::z][Cart::zz]+amb0*R[Cart::z][Cart::z][Cart::zz];
R[Cart::y][Cart::zz][Cart::zz]=R[Cart::yz][Cart::z][Cart::zz]+amb2*R[Cart::y][Cart::z][Cart::zz];
R[Cart::x][Cart::zz][Cart::zz]=R[Cart::xz][Cart::z][Cart::zz]+amb2*R[Cart::x][Cart::z][Cart::zz];
R[Cart::z][Cart::zz][Cart::zz]=R[Cart::zz][Cart::z][Cart::zz]+amb2*R[Cart::z][Cart::z][Cart::zz];
}
//------------------------------------------------------

//Integral d - d - d
if (_lmax_beta>1 && _lmax_alpha>1 && _lmax_gamma>1){

R[Cart::yy][Cart::yy][Cart::yy]=R[Cart::yyy][Cart::y][Cart::yy]+amb1*R[Cart::yy][Cart::y][Cart::yy];
R[Cart::xy][Cart::yy][Cart::yy]=R[Cart::xyy][Cart::y][Cart::yy]+amb1*R[Cart::xy][Cart::y][Cart::yy];
R[Cart::yz][Cart::yy][Cart::yy]=R[Cart::yyz][Cart::y][Cart::yy]+amb1*R[Cart::yz][Cart::y][Cart::yy];
R[Cart::xx][Cart::yy][Cart::yy]=R[Cart::xxy][Cart::y][Cart::yy]+amb1*R[Cart::xx][Cart::y][Cart::yy];
R[Cart::xz][Cart::yy][Cart::yy]=R[Cart::xyz][Cart::y][Cart::yy]+amb1*R[Cart::xz][Cart::y][Cart::yy];
R[Cart::zz][Cart::yy][Cart::yy]=R[Cart::yzz][Cart::y][Cart::yy]+amb1*R[Cart::zz][Cart::y][Cart::yy];
R[Cart::yy][Cart::xy][Cart::yy]=R[Cart::xyy][Cart::y][Cart::yy]+amb0*R[Cart::yy][Cart::y][Cart::yy];
R[Cart::xy][Cart::xy][Cart::yy]=R[Cart::xxy][Cart::y][Cart::yy]+amb0*R[Cart::xy][Cart::y][Cart::yy];
R[Cart::yz][Cart::xy][Cart::yy]=R[Cart::xyz][Cart::y][Cart::yy]+amb0*R[Cart::yz][Cart::y][Cart::yy];
R[Cart::xx][Cart::xy][Cart::yy]=R[Cart::xxx][Cart::y][Cart::yy]+amb0*R[Cart::xx][Cart::y][Cart::yy];
R[Cart::xz][Cart::xy][Cart::yy]=R[Cart::xxz][Cart::y][Cart::yy]+amb0*R[Cart::xz][Cart::y][Cart::yy];
R[Cart::zz][Cart::xy][Cart::yy]=R[Cart::xzz][Cart::y][Cart::yy]+amb0*R[Cart::zz][Cart::y][Cart::yy];
R[Cart::yy][Cart::yz][Cart::yy]=R[Cart::yyy][Cart::z][Cart::yy]+amb1*R[Cart::yy][Cart::z][Cart::yy];
R[Cart::xy][Cart::yz][Cart::yy]=R[Cart::xyy][Cart::z][Cart::yy]+amb1*R[Cart::xy][Cart::z][Cart::yy];
R[Cart::yz][Cart::yz][Cart::yy]=R[Cart::yyz][Cart::z][Cart::yy]+amb1*R[Cart::yz][Cart::z][Cart::yy];
R[Cart::xx][Cart::yz][Cart::yy]=R[Cart::xxy][Cart::z][Cart::yy]+amb1*R[Cart::xx][Cart::z][Cart::yy];
R[Cart::xz][Cart::yz][Cart::yy]=R[Cart::xyz][Cart::z][Cart::yy]+amb1*R[Cart::xz][Cart::z][Cart::yy];
R[Cart::zz][Cart::yz][Cart::yy]=R[Cart::yzz][Cart::z][Cart::yy]+amb1*R[Cart::zz][Cart::z][Cart::yy];
R[Cart::yy][Cart::xx][Cart::yy]=R[Cart::xyy][Cart::x][Cart::yy]+amb0*R[Cart::yy][Cart::x][Cart::yy];
R[Cart::xy][Cart::xx][Cart::yy]=R[Cart::xxy][Cart::x][Cart::yy]+amb0*R[Cart::xy][Cart::x][Cart::yy];
R[Cart::yz][Cart::xx][Cart::yy]=R[Cart::xyz][Cart::x][Cart::yy]+amb0*R[Cart::yz][Cart::x][Cart::yy];
R[Cart::xx][Cart::xx][Cart::yy]=R[Cart::xxx][Cart::x][Cart::yy]+amb0*R[Cart::xx][Cart::x][Cart::yy];
R[Cart::xz][Cart::xx][Cart::yy]=R[Cart::xxz][Cart::x][Cart::yy]+amb0*R[Cart::xz][Cart::x][Cart::yy];
R[Cart::zz][Cart::xx][Cart::yy]=R[Cart::xzz][Cart::x][Cart::yy]+amb0*R[Cart::zz][Cart::x][Cart::yy];
R[Cart::yy][Cart::xz][Cart::yy]=R[Cart::xyy][Cart::z][Cart::yy]+amb0*R[Cart::yy][Cart::z][Cart::yy];
R[Cart::xy][Cart::xz][Cart::yy]=R[Cart::xxy][Cart::z][Cart::yy]+amb0*R[Cart::xy][Cart::z][Cart::yy];
R[Cart::yz][Cart::xz][Cart::yy]=R[Cart::xyz][Cart::z][Cart::yy]+amb0*R[Cart::yz][Cart::z][Cart::yy];
R[Cart::xx][Cart::xz][Cart::yy]=R[Cart::xxx][Cart::z][Cart::yy]+amb0*R[Cart::xx][Cart::z][Cart::yy];
R[Cart::xz][Cart::xz][Cart::yy]=R[Cart::xxz][Cart::z][Cart::yy]+amb0*R[Cart::xz][Cart::z][Cart::yy];
R[Cart::zz][Cart::xz][Cart::yy]=R[Cart::xzz][Cart::z][Cart::yy]+amb0*R[Cart::zz][Cart::z][Cart::yy];
R[Cart::yy][Cart::zz][Cart::yy]=R[Cart::yyz][Cart::z][Cart::yy]+amb2*R[Cart::yy][Cart::z][Cart::yy];
R[Cart::xy][Cart::zz][Cart::yy]=R[Cart::xyz][Cart::z][Cart::yy]+amb2*R[Cart::xy][Cart::z][Cart::yy];
R[Cart::yz][Cart::zz][Cart::yy]=R[Cart::yzz][Cart::z][Cart::yy]+amb2*R[Cart::yz][Cart::z][Cart::yy];
R[Cart::xx][Cart::zz][Cart::yy]=R[Cart::xxz][Cart::z][Cart::yy]+amb2*R[Cart::xx][Cart::z][Cart::yy];
R[Cart::xz][Cart::zz][Cart::yy]=R[Cart::xzz][Cart::z][Cart::yy]+amb2*R[Cart::xz][Cart::z][Cart::yy];
R[Cart::zz][Cart::zz][Cart::yy]=R[Cart::zzz][Cart::z][Cart::yy]+amb2*R[Cart::zz][Cart::z][Cart::yy];
R[Cart::yy][Cart::yy][Cart::xy]=R[Cart::yyy][Cart::y][Cart::xy]+amb1*R[Cart::yy][Cart::y][Cart::xy];
R[Cart::xy][Cart::yy][Cart::xy]=R[Cart::xyy][Cart::y][Cart::xy]+amb1*R[Cart::xy][Cart::y][Cart::xy];
R[Cart::yz][Cart::yy][Cart::xy]=R[Cart::yyz][Cart::y][Cart::xy]+amb1*R[Cart::yz][Cart::y][Cart::xy];
R[Cart::xx][Cart::yy][Cart::xy]=R[Cart::xxy][Cart::y][Cart::xy]+amb1*R[Cart::xx][Cart::y][Cart::xy];
R[Cart::xz][Cart::yy][Cart::xy]=R[Cart::xyz][Cart::y][Cart::xy]+amb1*R[Cart::xz][Cart::y][Cart::xy];
R[Cart::zz][Cart::yy][Cart::xy]=R[Cart::yzz][Cart::y][Cart::xy]+amb1*R[Cart::zz][Cart::y][Cart::xy];
R[Cart::yy][Cart::xy][Cart::xy]=R[Cart::xyy][Cart::y][Cart::xy]+amb0*R[Cart::yy][Cart::y][Cart::xy];
R[Cart::xy][Cart::xy][Cart::xy]=R[Cart::xxy][Cart::y][Cart::xy]+amb0*R[Cart::xy][Cart::y][Cart::xy];
R[Cart::yz][Cart::xy][Cart::xy]=R[Cart::xyz][Cart::y][Cart::xy]+amb0*R[Cart::yz][Cart::y][Cart::xy];
R[Cart::xx][Cart::xy][Cart::xy]=R[Cart::xxx][Cart::y][Cart::xy]+amb0*R[Cart::xx][Cart::y][Cart::xy];
R[Cart::xz][Cart::xy][Cart::xy]=R[Cart::xxz][Cart::y][Cart::xy]+amb0*R[Cart::xz][Cart::y][Cart::xy];
R[Cart::zz][Cart::xy][Cart::xy]=R[Cart::xzz][Cart::y][Cart::xy]+amb0*R[Cart::zz][Cart::y][Cart::xy];
R[Cart::yy][Cart::yz][Cart::xy]=R[Cart::yyy][Cart::z][Cart::xy]+amb1*R[Cart::yy][Cart::z][Cart::xy];
R[Cart::xy][Cart::yz][Cart::xy]=R[Cart::xyy][Cart::z][Cart::xy]+amb1*R[Cart::xy][Cart::z][Cart::xy];
R[Cart::yz][Cart::yz][Cart::xy]=R[Cart::yyz][Cart::z][Cart::xy]+amb1*R[Cart::yz][Cart::z][Cart::xy];
R[Cart::xx][Cart::yz][Cart::xy]=R[Cart::xxy][Cart::z][Cart::xy]+amb1*R[Cart::xx][Cart::z][Cart::xy];
R[Cart::xz][Cart::yz][Cart::xy]=R[Cart::xyz][Cart::z][Cart::xy]+amb1*R[Cart::xz][Cart::z][Cart::xy];
R[Cart::zz][Cart::yz][Cart::xy]=R[Cart::yzz][Cart::z][Cart::xy]+amb1*R[Cart::zz][Cart::z][Cart::xy];
R[Cart::yy][Cart::xx][Cart::xy]=R[Cart::xyy][Cart::x][Cart::xy]+amb0*R[Cart::yy][Cart::x][Cart::xy];
R[Cart::xy][Cart::xx][Cart::xy]=R[Cart::xxy][Cart::x][Cart::xy]+amb0*R[Cart::xy][Cart::x][Cart::xy];
R[Cart::yz][Cart::xx][Cart::xy]=R[Cart::xyz][Cart::x][Cart::xy]+amb0*R[Cart::yz][Cart::x][Cart::xy];
R[Cart::xx][Cart::xx][Cart::xy]=R[Cart::xxx][Cart::x][Cart::xy]+amb0*R[Cart::xx][Cart::x][Cart::xy];
R[Cart::xz][Cart::xx][Cart::xy]=R[Cart::xxz][Cart::x][Cart::xy]+amb0*R[Cart::xz][Cart::x][Cart::xy];
R[Cart::zz][Cart::xx][Cart::xy]=R[Cart::xzz][Cart::x][Cart::xy]+amb0*R[Cart::zz][Cart::x][Cart::xy];
R[Cart::yy][Cart::xz][Cart::xy]=R[Cart::xyy][Cart::z][Cart::xy]+amb0*R[Cart::yy][Cart::z][Cart::xy];
R[Cart::xy][Cart::xz][Cart::xy]=R[Cart::xxy][Cart::z][Cart::xy]+amb0*R[Cart::xy][Cart::z][Cart::xy];
R[Cart::yz][Cart::xz][Cart::xy]=R[Cart::xyz][Cart::z][Cart::xy]+amb0*R[Cart::yz][Cart::z][Cart::xy];
R[Cart::xx][Cart::xz][Cart::xy]=R[Cart::xxx][Cart::z][Cart::xy]+amb0*R[Cart::xx][Cart::z][Cart::xy];
R[Cart::xz][Cart::xz][Cart::xy]=R[Cart::xxz][Cart::z][Cart::xy]+amb0*R[Cart::xz][Cart::z][Cart::xy];
R[Cart::zz][Cart::xz][Cart::xy]=R[Cart::xzz][Cart::z][Cart::xy]+amb0*R[Cart::zz][Cart::z][Cart::xy];
R[Cart::yy][Cart::zz][Cart::xy]=R[Cart::yyz][Cart::z][Cart::xy]+amb2*R[Cart::yy][Cart::z][Cart::xy];
R[Cart::xy][Cart::zz][Cart::xy]=R[Cart::xyz][Cart::z][Cart::xy]+amb2*R[Cart::xy][Cart::z][Cart::xy];
R[Cart::yz][Cart::zz][Cart::xy]=R[Cart::yzz][Cart::z][Cart::xy]+amb2*R[Cart::yz][Cart::z][Cart::xy];
R[Cart::xx][Cart::zz][Cart::xy]=R[Cart::xxz][Cart::z][Cart::xy]+amb2*R[Cart::xx][Cart::z][Cart::xy];
R[Cart::xz][Cart::zz][Cart::xy]=R[Cart::xzz][Cart::z][Cart::xy]+amb2*R[Cart::xz][Cart::z][Cart::xy];
R[Cart::zz][Cart::zz][Cart::xy]=R[Cart::zzz][Cart::z][Cart::xy]+amb2*R[Cart::zz][Cart::z][Cart::xy];
R[Cart::yy][Cart::yy][Cart::yz]=R[Cart::yyy][Cart::y][Cart::yz]+amb1*R[Cart::yy][Cart::y][Cart::yz];
R[Cart::xy][Cart::yy][Cart::yz]=R[Cart::xyy][Cart::y][Cart::yz]+amb1*R[Cart::xy][Cart::y][Cart::yz];
R[Cart::yz][Cart::yy][Cart::yz]=R[Cart::yyz][Cart::y][Cart::yz]+amb1*R[Cart::yz][Cart::y][Cart::yz];
R[Cart::xx][Cart::yy][Cart::yz]=R[Cart::xxy][Cart::y][Cart::yz]+amb1*R[Cart::xx][Cart::y][Cart::yz];
R[Cart::xz][Cart::yy][Cart::yz]=R[Cart::xyz][Cart::y][Cart::yz]+amb1*R[Cart::xz][Cart::y][Cart::yz];
R[Cart::zz][Cart::yy][Cart::yz]=R[Cart::yzz][Cart::y][Cart::yz]+amb1*R[Cart::zz][Cart::y][Cart::yz];
R[Cart::yy][Cart::xy][Cart::yz]=R[Cart::xyy][Cart::y][Cart::yz]+amb0*R[Cart::yy][Cart::y][Cart::yz];
R[Cart::xy][Cart::xy][Cart::yz]=R[Cart::xxy][Cart::y][Cart::yz]+amb0*R[Cart::xy][Cart::y][Cart::yz];
R[Cart::yz][Cart::xy][Cart::yz]=R[Cart::xyz][Cart::y][Cart::yz]+amb0*R[Cart::yz][Cart::y][Cart::yz];
R[Cart::xx][Cart::xy][Cart::yz]=R[Cart::xxx][Cart::y][Cart::yz]+amb0*R[Cart::xx][Cart::y][Cart::yz];
R[Cart::xz][Cart::xy][Cart::yz]=R[Cart::xxz][Cart::y][Cart::yz]+amb0*R[Cart::xz][Cart::y][Cart::yz];
R[Cart::zz][Cart::xy][Cart::yz]=R[Cart::xzz][Cart::y][Cart::yz]+amb0*R[Cart::zz][Cart::y][Cart::yz];
R[Cart::yy][Cart::yz][Cart::yz]=R[Cart::yyy][Cart::z][Cart::yz]+amb1*R[Cart::yy][Cart::z][Cart::yz];
R[Cart::xy][Cart::yz][Cart::yz]=R[Cart::xyy][Cart::z][Cart::yz]+amb1*R[Cart::xy][Cart::z][Cart::yz];
R[Cart::yz][Cart::yz][Cart::yz]=R[Cart::yyz][Cart::z][Cart::yz]+amb1*R[Cart::yz][Cart::z][Cart::yz];
R[Cart::xx][Cart::yz][Cart::yz]=R[Cart::xxy][Cart::z][Cart::yz]+amb1*R[Cart::xx][Cart::z][Cart::yz];
R[Cart::xz][Cart::yz][Cart::yz]=R[Cart::xyz][Cart::z][Cart::yz]+amb1*R[Cart::xz][Cart::z][Cart::yz];
R[Cart::zz][Cart::yz][Cart::yz]=R[Cart::yzz][Cart::z][Cart::yz]+amb1*R[Cart::zz][Cart::z][Cart::yz];
R[Cart::yy][Cart::xx][Cart::yz]=R[Cart::xyy][Cart::x][Cart::yz]+amb0*R[Cart::yy][Cart::x][Cart::yz];
R[Cart::xy][Cart::xx][Cart::yz]=R[Cart::xxy][Cart::x][Cart::yz]+amb0*R[Cart::xy][Cart::x][Cart::yz];
R[Cart::yz][Cart::xx][Cart::yz]=R[Cart::xyz][Cart::x][Cart::yz]+amb0*R[Cart::yz][Cart::x][Cart::yz];
R[Cart::xx][Cart::xx][Cart::yz]=R[Cart::xxx][Cart::x][Cart::yz]+amb0*R[Cart::xx][Cart::x][Cart::yz];
R[Cart::xz][Cart::xx][Cart::yz]=R[Cart::xxz][Cart::x][Cart::yz]+amb0*R[Cart::xz][Cart::x][Cart::yz];
R[Cart::zz][Cart::xx][Cart::yz]=R[Cart::xzz][Cart::x][Cart::yz]+amb0*R[Cart::zz][Cart::x][Cart::yz];
R[Cart::yy][Cart::xz][Cart::yz]=R[Cart::xyy][Cart::z][Cart::yz]+amb0*R[Cart::yy][Cart::z][Cart::yz];
R[Cart::xy][Cart::xz][Cart::yz]=R[Cart::xxy][Cart::z][Cart::yz]+amb0*R[Cart::xy][Cart::z][Cart::yz];
R[Cart::yz][Cart::xz][Cart::yz]=R[Cart::xyz][Cart::z][Cart::yz]+amb0*R[Cart::yz][Cart::z][Cart::yz];
R[Cart::xx][Cart::xz][Cart::yz]=R[Cart::xxx][Cart::z][Cart::yz]+amb0*R[Cart::xx][Cart::z][Cart::yz];
R[Cart::xz][Cart::xz][Cart::yz]=R[Cart::xxz][Cart::z][Cart::yz]+amb0*R[Cart::xz][Cart::z][Cart::yz];
R[Cart::zz][Cart::xz][Cart::yz]=R[Cart::xzz][Cart::z][Cart::yz]+amb0*R[Cart::zz][Cart::z][Cart::yz];
R[Cart::yy][Cart::zz][Cart::yz]=R[Cart::yyz][Cart::z][Cart::yz]+amb2*R[Cart::yy][Cart::z][Cart::yz];
R[Cart::xy][Cart::zz][Cart::yz]=R[Cart::xyz][Cart::z][Cart::yz]+amb2*R[Cart::xy][Cart::z][Cart::yz];
R[Cart::yz][Cart::zz][Cart::yz]=R[Cart::yzz][Cart::z][Cart::yz]+amb2*R[Cart::yz][Cart::z][Cart::yz];
R[Cart::xx][Cart::zz][Cart::yz]=R[Cart::xxz][Cart::z][Cart::yz]+amb2*R[Cart::xx][Cart::z][Cart::yz];
R[Cart::xz][Cart::zz][Cart::yz]=R[Cart::xzz][Cart::z][Cart::yz]+amb2*R[Cart::xz][Cart::z][Cart::yz];
R[Cart::zz][Cart::zz][Cart::yz]=R[Cart::zzz][Cart::z][Cart::yz]+amb2*R[Cart::zz][Cart::z][Cart::yz];
R[Cart::yy][Cart::yy][Cart::xx]=R[Cart::yyy][Cart::y][Cart::xx]+amb1*R[Cart::yy][Cart::y][Cart::xx];
R[Cart::xy][Cart::yy][Cart::xx]=R[Cart::xyy][Cart::y][Cart::xx]+amb1*R[Cart::xy][Cart::y][Cart::xx];
R[Cart::yz][Cart::yy][Cart::xx]=R[Cart::yyz][Cart::y][Cart::xx]+amb1*R[Cart::yz][Cart::y][Cart::xx];
R[Cart::xx][Cart::yy][Cart::xx]=R[Cart::xxy][Cart::y][Cart::xx]+amb1*R[Cart::xx][Cart::y][Cart::xx];
R[Cart::xz][Cart::yy][Cart::xx]=R[Cart::xyz][Cart::y][Cart::xx]+amb1*R[Cart::xz][Cart::y][Cart::xx];
R[Cart::zz][Cart::yy][Cart::xx]=R[Cart::yzz][Cart::y][Cart::xx]+amb1*R[Cart::zz][Cart::y][Cart::xx];
R[Cart::yy][Cart::xy][Cart::xx]=R[Cart::xyy][Cart::y][Cart::xx]+amb0*R[Cart::yy][Cart::y][Cart::xx];
R[Cart::xy][Cart::xy][Cart::xx]=R[Cart::xxy][Cart::y][Cart::xx]+amb0*R[Cart::xy][Cart::y][Cart::xx];
R[Cart::yz][Cart::xy][Cart::xx]=R[Cart::xyz][Cart::y][Cart::xx]+amb0*R[Cart::yz][Cart::y][Cart::xx];
R[Cart::xx][Cart::xy][Cart::xx]=R[Cart::xxx][Cart::y][Cart::xx]+amb0*R[Cart::xx][Cart::y][Cart::xx];
R[Cart::xz][Cart::xy][Cart::xx]=R[Cart::xxz][Cart::y][Cart::xx]+amb0*R[Cart::xz][Cart::y][Cart::xx];
R[Cart::zz][Cart::xy][Cart::xx]=R[Cart::xzz][Cart::y][Cart::xx]+amb0*R[Cart::zz][Cart::y][Cart::xx];
R[Cart::yy][Cart::yz][Cart::xx]=R[Cart::yyy][Cart::z][Cart::xx]+amb1*R[Cart::yy][Cart::z][Cart::xx];
R[Cart::xy][Cart::yz][Cart::xx]=R[Cart::xyy][Cart::z][Cart::xx]+amb1*R[Cart::xy][Cart::z][Cart::xx];
R[Cart::yz][Cart::yz][Cart::xx]=R[Cart::yyz][Cart::z][Cart::xx]+amb1*R[Cart::yz][Cart::z][Cart::xx];
R[Cart::xx][Cart::yz][Cart::xx]=R[Cart::xxy][Cart::z][Cart::xx]+amb1*R[Cart::xx][Cart::z][Cart::xx];
R[Cart::xz][Cart::yz][Cart::xx]=R[Cart::xyz][Cart::z][Cart::xx]+amb1*R[Cart::xz][Cart::z][Cart::xx];
R[Cart::zz][Cart::yz][Cart::xx]=R[Cart::yzz][Cart::z][Cart::xx]+amb1*R[Cart::zz][Cart::z][Cart::xx];
R[Cart::yy][Cart::xx][Cart::xx]=R[Cart::xyy][Cart::x][Cart::xx]+amb0*R[Cart::yy][Cart::x][Cart::xx];
R[Cart::xy][Cart::xx][Cart::xx]=R[Cart::xxy][Cart::x][Cart::xx]+amb0*R[Cart::xy][Cart::x][Cart::xx];
R[Cart::yz][Cart::xx][Cart::xx]=R[Cart::xyz][Cart::x][Cart::xx]+amb0*R[Cart::yz][Cart::x][Cart::xx];
R[Cart::xx][Cart::xx][Cart::xx]=R[Cart::xxx][Cart::x][Cart::xx]+amb0*R[Cart::xx][Cart::x][Cart::xx];
R[Cart::xz][Cart::xx][Cart::xx]=R[Cart::xxz][Cart::x][Cart::xx]+amb0*R[Cart::xz][Cart::x][Cart::xx];
R[Cart::zz][Cart::xx][Cart::xx]=R[Cart::xzz][Cart::x][Cart::xx]+amb0*R[Cart::zz][Cart::x][Cart::xx];
R[Cart::yy][Cart::xz][Cart::xx]=R[Cart::xyy][Cart::z][Cart::xx]+amb0*R[Cart::yy][Cart::z][Cart::xx];
R[Cart::xy][Cart::xz][Cart::xx]=R[Cart::xxy][Cart::z][Cart::xx]+amb0*R[Cart::xy][Cart::z][Cart::xx];
R[Cart::yz][Cart::xz][Cart::xx]=R[Cart::xyz][Cart::z][Cart::xx]+amb0*R[Cart::yz][Cart::z][Cart::xx];
R[Cart::xx][Cart::xz][Cart::xx]=R[Cart::xxx][Cart::z][Cart::xx]+amb0*R[Cart::xx][Cart::z][Cart::xx];
R[Cart::xz][Cart::xz][Cart::xx]=R[Cart::xxz][Cart::z][Cart::xx]+amb0*R[Cart::xz][Cart::z][Cart::xx];
R[Cart::zz][Cart::xz][Cart::xx]=R[Cart::xzz][Cart::z][Cart::xx]+amb0*R[Cart::zz][Cart::z][Cart::xx];
R[Cart::yy][Cart::zz][Cart::xx]=R[Cart::yyz][Cart::z][Cart::xx]+amb2*R[Cart::yy][Cart::z][Cart::xx];
R[Cart::xy][Cart::zz][Cart::xx]=R[Cart::xyz][Cart::z][Cart::xx]+amb2*R[Cart::xy][Cart::z][Cart::xx];
R[Cart::yz][Cart::zz][Cart::xx]=R[Cart::yzz][Cart::z][Cart::xx]+amb2*R[Cart::yz][Cart::z][Cart::xx];
R[Cart::xx][Cart::zz][Cart::xx]=R[Cart::xxz][Cart::z][Cart::xx]+amb2*R[Cart::xx][Cart::z][Cart::xx];
R[Cart::xz][Cart::zz][Cart::xx]=R[Cart::xzz][Cart::z][Cart::xx]+amb2*R[Cart::xz][Cart::z][Cart::xx];
R[Cart::zz][Cart::zz][Cart::xx]=R[Cart::zzz][Cart::z][Cart::xx]+amb2*R[Cart::zz][Cart::z][Cart::xx];
R[Cart::yy][Cart::yy][Cart::xz]=R[Cart::yyy][Cart::y][Cart::xz]+amb1*R[Cart::yy][Cart::y][Cart::xz];
R[Cart::xy][Cart::yy][Cart::xz]=R[Cart::xyy][Cart::y][Cart::xz]+amb1*R[Cart::xy][Cart::y][Cart::xz];
R[Cart::yz][Cart::yy][Cart::xz]=R[Cart::yyz][Cart::y][Cart::xz]+amb1*R[Cart::yz][Cart::y][Cart::xz];
R[Cart::xx][Cart::yy][Cart::xz]=R[Cart::xxy][Cart::y][Cart::xz]+amb1*R[Cart::xx][Cart::y][Cart::xz];
R[Cart::xz][Cart::yy][Cart::xz]=R[Cart::xyz][Cart::y][Cart::xz]+amb1*R[Cart::xz][Cart::y][Cart::xz];
R[Cart::zz][Cart::yy][Cart::xz]=R[Cart::yzz][Cart::y][Cart::xz]+amb1*R[Cart::zz][Cart::y][Cart::xz];
R[Cart::yy][Cart::xy][Cart::xz]=R[Cart::xyy][Cart::y][Cart::xz]+amb0*R[Cart::yy][Cart::y][Cart::xz];
R[Cart::xy][Cart::xy][Cart::xz]=R[Cart::xxy][Cart::y][Cart::xz]+amb0*R[Cart::xy][Cart::y][Cart::xz];
R[Cart::yz][Cart::xy][Cart::xz]=R[Cart::xyz][Cart::y][Cart::xz]+amb0*R[Cart::yz][Cart::y][Cart::xz];
R[Cart::xx][Cart::xy][Cart::xz]=R[Cart::xxx][Cart::y][Cart::xz]+amb0*R[Cart::xx][Cart::y][Cart::xz];
R[Cart::xz][Cart::xy][Cart::xz]=R[Cart::xxz][Cart::y][Cart::xz]+amb0*R[Cart::xz][Cart::y][Cart::xz];
R[Cart::zz][Cart::xy][Cart::xz]=R[Cart::xzz][Cart::y][Cart::xz]+amb0*R[Cart::zz][Cart::y][Cart::xz];
R[Cart::yy][Cart::yz][Cart::xz]=R[Cart::yyy][Cart::z][Cart::xz]+amb1*R[Cart::yy][Cart::z][Cart::xz];
R[Cart::xy][Cart::yz][Cart::xz]=R[Cart::xyy][Cart::z][Cart::xz]+amb1*R[Cart::xy][Cart::z][Cart::xz];
R[Cart::yz][Cart::yz][Cart::xz]=R[Cart::yyz][Cart::z][Cart::xz]+amb1*R[Cart::yz][Cart::z][Cart::xz];
R[Cart::xx][Cart::yz][Cart::xz]=R[Cart::xxy][Cart::z][Cart::xz]+amb1*R[Cart::xx][Cart::z][Cart::xz];
R[Cart::xz][Cart::yz][Cart::xz]=R[Cart::xyz][Cart::z][Cart::xz]+amb1*R[Cart::xz][Cart::z][Cart::xz];
R[Cart::zz][Cart::yz][Cart::xz]=R[Cart::yzz][Cart::z][Cart::xz]+amb1*R[Cart::zz][Cart::z][Cart::xz];
R[Cart::yy][Cart::xx][Cart::xz]=R[Cart::xyy][Cart::x][Cart::xz]+amb0*R[Cart::yy][Cart::x][Cart::xz];
R[Cart::xy][Cart::xx][Cart::xz]=R[Cart::xxy][Cart::x][Cart::xz]+amb0*R[Cart::xy][Cart::x][Cart::xz];
R[Cart::yz][Cart::xx][Cart::xz]=R[Cart::xyz][Cart::x][Cart::xz]+amb0*R[Cart::yz][Cart::x][Cart::xz];
R[Cart::xx][Cart::xx][Cart::xz]=R[Cart::xxx][Cart::x][Cart::xz]+amb0*R[Cart::xx][Cart::x][Cart::xz];
R[Cart::xz][Cart::xx][Cart::xz]=R[Cart::xxz][Cart::x][Cart::xz]+amb0*R[Cart::xz][Cart::x][Cart::xz];
R[Cart::zz][Cart::xx][Cart::xz]=R[Cart::xzz][Cart::x][Cart::xz]+amb0*R[Cart::zz][Cart::x][Cart::xz];
R[Cart::yy][Cart::xz][Cart::xz]=R[Cart::xyy][Cart::z][Cart::xz]+amb0*R[Cart::yy][Cart::z][Cart::xz];
R[Cart::xy][Cart::xz][Cart::xz]=R[Cart::xxy][Cart::z][Cart::xz]+amb0*R[Cart::xy][Cart::z][Cart::xz];
R[Cart::yz][Cart::xz][Cart::xz]=R[Cart::xyz][Cart::z][Cart::xz]+amb0*R[Cart::yz][Cart::z][Cart::xz];
R[Cart::xx][Cart::xz][Cart::xz]=R[Cart::xxx][Cart::z][Cart::xz]+amb0*R[Cart::xx][Cart::z][Cart::xz];
R[Cart::xz][Cart::xz][Cart::xz]=R[Cart::xxz][Cart::z][Cart::xz]+amb0*R[Cart::xz][Cart::z][Cart::xz];
R[Cart::zz][Cart::xz][Cart::xz]=R[Cart::xzz][Cart::z][Cart::xz]+amb0*R[Cart::zz][Cart::z][Cart::xz];
R[Cart::yy][Cart::zz][Cart::xz]=R[Cart::yyz][Cart::z][Cart::xz]+amb2*R[Cart::yy][Cart::z][Cart::xz];
R[Cart::xy][Cart::zz][Cart::xz]=R[Cart::xyz][Cart::z][Cart::xz]+amb2*R[Cart::xy][Cart::z][Cart::xz];
R[Cart::yz][Cart::zz][Cart::xz]=R[Cart::yzz][Cart::z][Cart::xz]+amb2*R[Cart::yz][Cart::z][Cart::xz];
R[Cart::xx][Cart::zz][Cart::xz]=R[Cart::xxz][Cart::z][Cart::xz]+amb2*R[Cart::xx][Cart::z][Cart::xz];
R[Cart::xz][Cart::zz][Cart::xz]=R[Cart::xzz][Cart::z][Cart::xz]+amb2*R[Cart::xz][Cart::z][Cart::xz];
R[Cart::zz][Cart::zz][Cart::xz]=R[Cart::zzz][Cart::z][Cart::xz]+amb2*R[Cart::zz][Cart::z][Cart::xz];
R[Cart::yy][Cart::yy][Cart::zz]=R[Cart::yyy][Cart::y][Cart::zz]+amb1*R[Cart::yy][Cart::y][Cart::zz];
R[Cart::xy][Cart::yy][Cart::zz]=R[Cart::xyy][Cart::y][Cart::zz]+amb1*R[Cart::xy][Cart::y][Cart::zz];
R[Cart::yz][Cart::yy][Cart::zz]=R[Cart::yyz][Cart::y][Cart::zz]+amb1*R[Cart::yz][Cart::y][Cart::zz];
R[Cart::xx][Cart::yy][Cart::zz]=R[Cart::xxy][Cart::y][Cart::zz]+amb1*R[Cart::xx][Cart::y][Cart::zz];
R[Cart::xz][Cart::yy][Cart::zz]=R[Cart::xyz][Cart::y][Cart::zz]+amb1*R[Cart::xz][Cart::y][Cart::zz];
R[Cart::zz][Cart::yy][Cart::zz]=R[Cart::yzz][Cart::y][Cart::zz]+amb1*R[Cart::zz][Cart::y][Cart::zz];
R[Cart::yy][Cart::xy][Cart::zz]=R[Cart::xyy][Cart::y][Cart::zz]+amb0*R[Cart::yy][Cart::y][Cart::zz];
R[Cart::xy][Cart::xy][Cart::zz]=R[Cart::xxy][Cart::y][Cart::zz]+amb0*R[Cart::xy][Cart::y][Cart::zz];
R[Cart::yz][Cart::xy][Cart::zz]=R[Cart::xyz][Cart::y][Cart::zz]+amb0*R[Cart::yz][Cart::y][Cart::zz];
R[Cart::xx][Cart::xy][Cart::zz]=R[Cart::xxx][Cart::y][Cart::zz]+amb0*R[Cart::xx][Cart::y][Cart::zz];
R[Cart::xz][Cart::xy][Cart::zz]=R[Cart::xxz][Cart::y][Cart::zz]+amb0*R[Cart::xz][Cart::y][Cart::zz];
R[Cart::zz][Cart::xy][Cart::zz]=R[Cart::xzz][Cart::y][Cart::zz]+amb0*R[Cart::zz][Cart::y][Cart::zz];
R[Cart::yy][Cart::yz][Cart::zz]=R[Cart::yyy][Cart::z][Cart::zz]+amb1*R[Cart::yy][Cart::z][Cart::zz];
R[Cart::xy][Cart::yz][Cart::zz]=R[Cart::xyy][Cart::z][Cart::zz]+amb1*R[Cart::xy][Cart::z][Cart::zz];
R[Cart::yz][Cart::yz][Cart::zz]=R[Cart::yyz][Cart::z][Cart::zz]+amb1*R[Cart::yz][Cart::z][Cart::zz];
R[Cart::xx][Cart::yz][Cart::zz]=R[Cart::xxy][Cart::z][Cart::zz]+amb1*R[Cart::xx][Cart::z][Cart::zz];
R[Cart::xz][Cart::yz][Cart::zz]=R[Cart::xyz][Cart::z][Cart::zz]+amb1*R[Cart::xz][Cart::z][Cart::zz];
R[Cart::zz][Cart::yz][Cart::zz]=R[Cart::yzz][Cart::z][Cart::zz]+amb1*R[Cart::zz][Cart::z][Cart::zz];
R[Cart::yy][Cart::xx][Cart::zz]=R[Cart::xyy][Cart::x][Cart::zz]+amb0*R[Cart::yy][Cart::x][Cart::zz];
R[Cart::xy][Cart::xx][Cart::zz]=R[Cart::xxy][Cart::x][Cart::zz]+amb0*R[Cart::xy][Cart::x][Cart::zz];
R[Cart::yz][Cart::xx][Cart::zz]=R[Cart::xyz][Cart::x][Cart::zz]+amb0*R[Cart::yz][Cart::x][Cart::zz];
R[Cart::xx][Cart::xx][Cart::zz]=R[Cart::xxx][Cart::x][Cart::zz]+amb0*R[Cart::xx][Cart::x][Cart::zz];
R[Cart::xz][Cart::xx][Cart::zz]=R[Cart::xxz][Cart::x][Cart::zz]+amb0*R[Cart::xz][Cart::x][Cart::zz];
R[Cart::zz][Cart::xx][Cart::zz]=R[Cart::xzz][Cart::x][Cart::zz]+amb0*R[Cart::zz][Cart::x][Cart::zz];
R[Cart::yy][Cart::xz][Cart::zz]=R[Cart::xyy][Cart::z][Cart::zz]+amb0*R[Cart::yy][Cart::z][Cart::zz];
R[Cart::xy][Cart::xz][Cart::zz]=R[Cart::xxy][Cart::z][Cart::zz]+amb0*R[Cart::xy][Cart::z][Cart::zz];
R[Cart::yz][Cart::xz][Cart::zz]=R[Cart::xyz][Cart::z][Cart::zz]+amb0*R[Cart::yz][Cart::z][Cart::zz];
R[Cart::xx][Cart::xz][Cart::zz]=R[Cart::xxx][Cart::z][Cart::zz]+amb0*R[Cart::xx][Cart::z][Cart::zz];
R[Cart::xz][Cart::xz][Cart::zz]=R[Cart::xxz][Cart::z][Cart::zz]+amb0*R[Cart::xz][Cart::z][Cart::zz];
R[Cart::zz][Cart::xz][Cart::zz]=R[Cart::xzz][Cart::z][Cart::zz]+amb0*R[Cart::zz][Cart::z][Cart::zz];
R[Cart::yy][Cart::zz][Cart::zz]=R[Cart::yyz][Cart::z][Cart::zz]+amb2*R[Cart::yy][Cart::z][Cart::zz];
R[Cart::xy][Cart::zz][Cart::zz]=R[Cart::xyz][Cart::z][Cart::zz]+amb2*R[Cart::xy][Cart::z][Cart::zz];
R[Cart::yz][Cart::zz][Cart::zz]=R[Cart::yzz][Cart::z][Cart::zz]+amb2*R[Cart::yz][Cart::z][Cart::zz];
R[Cart::xx][Cart::zz][Cart::zz]=R[Cart::xxz][Cart::z][Cart::zz]+amb2*R[Cart::xx][Cart::z][Cart::zz];
R[Cart::xz][Cart::zz][Cart::zz]=R[Cart::xzz][Cart::z][Cart::zz]+amb2*R[Cart::xz][Cart::z][Cart::zz];
R[Cart::zz][Cart::zz][Cart::zz]=R[Cart::zzz][Cart::z][Cart::zz]+amb2*R[Cart::zz][Cart::z][Cart::zz];
}
//------------------------------------------------------

//Integral s - d - f
if (_lmax_beta>1 && _lmax_gamma>2){

R[Cart::s][Cart::yy][Cart::yyy]=R[Cart::y][Cart::y][Cart::yyy]+amb1*R[Cart::s][Cart::y][Cart::yyy];
R[Cart::s][Cart::xy][Cart::yyy]=R[Cart::x][Cart::y][Cart::yyy]+amb0*R[Cart::s][Cart::y][Cart::yyy];
R[Cart::s][Cart::yz][Cart::yyy]=R[Cart::y][Cart::z][Cart::yyy]+amb1*R[Cart::s][Cart::z][Cart::yyy];
R[Cart::s][Cart::xx][Cart::yyy]=R[Cart::x][Cart::x][Cart::yyy]+amb0*R[Cart::s][Cart::x][Cart::yyy];
R[Cart::s][Cart::xz][Cart::yyy]=R[Cart::x][Cart::z][Cart::yyy]+amb0*R[Cart::s][Cart::z][Cart::yyy];
R[Cart::s][Cart::zz][Cart::yyy]=R[Cart::z][Cart::z][Cart::yyy]+amb2*R[Cart::s][Cart::z][Cart::yyy];
R[Cart::s][Cart::yy][Cart::xyy]=R[Cart::y][Cart::y][Cart::xyy]+amb1*R[Cart::s][Cart::y][Cart::xyy];
R[Cart::s][Cart::xy][Cart::xyy]=R[Cart::x][Cart::y][Cart::xyy]+amb0*R[Cart::s][Cart::y][Cart::xyy];
R[Cart::s][Cart::yz][Cart::xyy]=R[Cart::y][Cart::z][Cart::xyy]+amb1*R[Cart::s][Cart::z][Cart::xyy];
R[Cart::s][Cart::xx][Cart::xyy]=R[Cart::x][Cart::x][Cart::xyy]+amb0*R[Cart::s][Cart::x][Cart::xyy];
R[Cart::s][Cart::xz][Cart::xyy]=R[Cart::x][Cart::z][Cart::xyy]+amb0*R[Cart::s][Cart::z][Cart::xyy];
R[Cart::s][Cart::zz][Cart::xyy]=R[Cart::z][Cart::z][Cart::xyy]+amb2*R[Cart::s][Cart::z][Cart::xyy];
R[Cart::s][Cart::yy][Cart::yyz]=R[Cart::y][Cart::y][Cart::yyz]+amb1*R[Cart::s][Cart::y][Cart::yyz];
R[Cart::s][Cart::xy][Cart::yyz]=R[Cart::x][Cart::y][Cart::yyz]+amb0*R[Cart::s][Cart::y][Cart::yyz];
R[Cart::s][Cart::yz][Cart::yyz]=R[Cart::y][Cart::z][Cart::yyz]+amb1*R[Cart::s][Cart::z][Cart::yyz];
R[Cart::s][Cart::xx][Cart::yyz]=R[Cart::x][Cart::x][Cart::yyz]+amb0*R[Cart::s][Cart::x][Cart::yyz];
R[Cart::s][Cart::xz][Cart::yyz]=R[Cart::x][Cart::z][Cart::yyz]+amb0*R[Cart::s][Cart::z][Cart::yyz];
R[Cart::s][Cart::zz][Cart::yyz]=R[Cart::z][Cart::z][Cart::yyz]+amb2*R[Cart::s][Cart::z][Cart::yyz];
R[Cart::s][Cart::yy][Cart::xxy]=R[Cart::y][Cart::y][Cart::xxy]+amb1*R[Cart::s][Cart::y][Cart::xxy];
R[Cart::s][Cart::xy][Cart::xxy]=R[Cart::x][Cart::y][Cart::xxy]+amb0*R[Cart::s][Cart::y][Cart::xxy];
R[Cart::s][Cart::yz][Cart::xxy]=R[Cart::y][Cart::z][Cart::xxy]+amb1*R[Cart::s][Cart::z][Cart::xxy];
R[Cart::s][Cart::xx][Cart::xxy]=R[Cart::x][Cart::x][Cart::xxy]+amb0*R[Cart::s][Cart::x][Cart::xxy];
R[Cart::s][Cart::xz][Cart::xxy]=R[Cart::x][Cart::z][Cart::xxy]+amb0*R[Cart::s][Cart::z][Cart::xxy];
R[Cart::s][Cart::zz][Cart::xxy]=R[Cart::z][Cart::z][Cart::xxy]+amb2*R[Cart::s][Cart::z][Cart::xxy];
R[Cart::s][Cart::yy][Cart::xyz]=R[Cart::y][Cart::y][Cart::xyz]+amb1*R[Cart::s][Cart::y][Cart::xyz];
R[Cart::s][Cart::xy][Cart::xyz]=R[Cart::x][Cart::y][Cart::xyz]+amb0*R[Cart::s][Cart::y][Cart::xyz];
R[Cart::s][Cart::yz][Cart::xyz]=R[Cart::y][Cart::z][Cart::xyz]+amb1*R[Cart::s][Cart::z][Cart::xyz];
R[Cart::s][Cart::xx][Cart::xyz]=R[Cart::x][Cart::x][Cart::xyz]+amb0*R[Cart::s][Cart::x][Cart::xyz];
R[Cart::s][Cart::xz][Cart::xyz]=R[Cart::x][Cart::z][Cart::xyz]+amb0*R[Cart::s][Cart::z][Cart::xyz];
R[Cart::s][Cart::zz][Cart::xyz]=R[Cart::z][Cart::z][Cart::xyz]+amb2*R[Cart::s][Cart::z][Cart::xyz];
R[Cart::s][Cart::yy][Cart::yzz]=R[Cart::y][Cart::y][Cart::yzz]+amb1*R[Cart::s][Cart::y][Cart::yzz];
R[Cart::s][Cart::xy][Cart::yzz]=R[Cart::x][Cart::y][Cart::yzz]+amb0*R[Cart::s][Cart::y][Cart::yzz];
R[Cart::s][Cart::yz][Cart::yzz]=R[Cart::y][Cart::z][Cart::yzz]+amb1*R[Cart::s][Cart::z][Cart::yzz];
R[Cart::s][Cart::xx][Cart::yzz]=R[Cart::x][Cart::x][Cart::yzz]+amb0*R[Cart::s][Cart::x][Cart::yzz];
R[Cart::s][Cart::xz][Cart::yzz]=R[Cart::x][Cart::z][Cart::yzz]+amb0*R[Cart::s][Cart::z][Cart::yzz];
R[Cart::s][Cart::zz][Cart::yzz]=R[Cart::z][Cart::z][Cart::yzz]+amb2*R[Cart::s][Cart::z][Cart::yzz];
R[Cart::s][Cart::yy][Cart::xxx]=R[Cart::y][Cart::y][Cart::xxx]+amb1*R[Cart::s][Cart::y][Cart::xxx];
R[Cart::s][Cart::xy][Cart::xxx]=R[Cart::x][Cart::y][Cart::xxx]+amb0*R[Cart::s][Cart::y][Cart::xxx];
R[Cart::s][Cart::yz][Cart::xxx]=R[Cart::y][Cart::z][Cart::xxx]+amb1*R[Cart::s][Cart::z][Cart::xxx];
R[Cart::s][Cart::xx][Cart::xxx]=R[Cart::x][Cart::x][Cart::xxx]+amb0*R[Cart::s][Cart::x][Cart::xxx];
R[Cart::s][Cart::xz][Cart::xxx]=R[Cart::x][Cart::z][Cart::xxx]+amb0*R[Cart::s][Cart::z][Cart::xxx];
R[Cart::s][Cart::zz][Cart::xxx]=R[Cart::z][Cart::z][Cart::xxx]+amb2*R[Cart::s][Cart::z][Cart::xxx];
R[Cart::s][Cart::yy][Cart::xxz]=R[Cart::y][Cart::y][Cart::xxz]+amb1*R[Cart::s][Cart::y][Cart::xxz];
R[Cart::s][Cart::xy][Cart::xxz]=R[Cart::x][Cart::y][Cart::xxz]+amb0*R[Cart::s][Cart::y][Cart::xxz];
R[Cart::s][Cart::yz][Cart::xxz]=R[Cart::y][Cart::z][Cart::xxz]+amb1*R[Cart::s][Cart::z][Cart::xxz];
R[Cart::s][Cart::xx][Cart::xxz]=R[Cart::x][Cart::x][Cart::xxz]+amb0*R[Cart::s][Cart::x][Cart::xxz];
R[Cart::s][Cart::xz][Cart::xxz]=R[Cart::x][Cart::z][Cart::xxz]+amb0*R[Cart::s][Cart::z][Cart::xxz];
R[Cart::s][Cart::zz][Cart::xxz]=R[Cart::z][Cart::z][Cart::xxz]+amb2*R[Cart::s][Cart::z][Cart::xxz];
R[Cart::s][Cart::yy][Cart::xzz]=R[Cart::y][Cart::y][Cart::xzz]+amb1*R[Cart::s][Cart::y][Cart::xzz];
R[Cart::s][Cart::xy][Cart::xzz]=R[Cart::x][Cart::y][Cart::xzz]+amb0*R[Cart::s][Cart::y][Cart::xzz];
R[Cart::s][Cart::yz][Cart::xzz]=R[Cart::y][Cart::z][Cart::xzz]+amb1*R[Cart::s][Cart::z][Cart::xzz];
R[Cart::s][Cart::xx][Cart::xzz]=R[Cart::x][Cart::x][Cart::xzz]+amb0*R[Cart::s][Cart::x][Cart::xzz];
R[Cart::s][Cart::xz][Cart::xzz]=R[Cart::x][Cart::z][Cart::xzz]+amb0*R[Cart::s][Cart::z][Cart::xzz];
R[Cart::s][Cart::zz][Cart::xzz]=R[Cart::z][Cart::z][Cart::xzz]+amb2*R[Cart::s][Cart::z][Cart::xzz];
R[Cart::s][Cart::yy][Cart::zzz]=R[Cart::y][Cart::y][Cart::zzz]+amb1*R[Cart::s][Cart::y][Cart::zzz];
R[Cart::s][Cart::xy][Cart::zzz]=R[Cart::x][Cart::y][Cart::zzz]+amb0*R[Cart::s][Cart::y][Cart::zzz];
R[Cart::s][Cart::yz][Cart::zzz]=R[Cart::y][Cart::z][Cart::zzz]+amb1*R[Cart::s][Cart::z][Cart::zzz];
R[Cart::s][Cart::xx][Cart::zzz]=R[Cart::x][Cart::x][Cart::zzz]+amb0*R[Cart::s][Cart::x][Cart::zzz];
R[Cart::s][Cart::xz][Cart::zzz]=R[Cart::x][Cart::z][Cart::zzz]+amb0*R[Cart::s][Cart::z][Cart::zzz];
R[Cart::s][Cart::zz][Cart::zzz]=R[Cart::z][Cart::z][Cart::zzz]+amb2*R[Cart::s][Cart::z][Cart::zzz];
}
//------------------------------------------------------

//Integral p - d - f
if (_lmax_beta>1 && _lmax_alpha>0 && _lmax_gamma>2){

R[Cart::y][Cart::yy][Cart::yyy]=R[Cart::yy][Cart::y][Cart::yyy]+amb1*R[Cart::y][Cart::y][Cart::yyy];
R[Cart::x][Cart::yy][Cart::yyy]=R[Cart::xy][Cart::y][Cart::yyy]+amb1*R[Cart::x][Cart::y][Cart::yyy];
R[Cart::z][Cart::yy][Cart::yyy]=R[Cart::yz][Cart::y][Cart::yyy]+amb1*R[Cart::z][Cart::y][Cart::yyy];
R[Cart::y][Cart::xy][Cart::yyy]=R[Cart::xy][Cart::y][Cart::yyy]+amb0*R[Cart::y][Cart::y][Cart::yyy];
R[Cart::x][Cart::xy][Cart::yyy]=R[Cart::xx][Cart::y][Cart::yyy]+amb0*R[Cart::x][Cart::y][Cart::yyy];
R[Cart::z][Cart::xy][Cart::yyy]=R[Cart::xz][Cart::y][Cart::yyy]+amb0*R[Cart::z][Cart::y][Cart::yyy];
R[Cart::y][Cart::yz][Cart::yyy]=R[Cart::yy][Cart::z][Cart::yyy]+amb1*R[Cart::y][Cart::z][Cart::yyy];
R[Cart::x][Cart::yz][Cart::yyy]=R[Cart::xy][Cart::z][Cart::yyy]+amb1*R[Cart::x][Cart::z][Cart::yyy];
R[Cart::z][Cart::yz][Cart::yyy]=R[Cart::yz][Cart::z][Cart::yyy]+amb1*R[Cart::z][Cart::z][Cart::yyy];
R[Cart::y][Cart::xx][Cart::yyy]=R[Cart::xy][Cart::x][Cart::yyy]+amb0*R[Cart::y][Cart::x][Cart::yyy];
R[Cart::x][Cart::xx][Cart::yyy]=R[Cart::xx][Cart::x][Cart::yyy]+amb0*R[Cart::x][Cart::x][Cart::yyy];
R[Cart::z][Cart::xx][Cart::yyy]=R[Cart::xz][Cart::x][Cart::yyy]+amb0*R[Cart::z][Cart::x][Cart::yyy];
R[Cart::y][Cart::xz][Cart::yyy]=R[Cart::xy][Cart::z][Cart::yyy]+amb0*R[Cart::y][Cart::z][Cart::yyy];
R[Cart::x][Cart::xz][Cart::yyy]=R[Cart::xx][Cart::z][Cart::yyy]+amb0*R[Cart::x][Cart::z][Cart::yyy];
R[Cart::z][Cart::xz][Cart::yyy]=R[Cart::xz][Cart::z][Cart::yyy]+amb0*R[Cart::z][Cart::z][Cart::yyy];
R[Cart::y][Cart::zz][Cart::yyy]=R[Cart::yz][Cart::z][Cart::yyy]+amb2*R[Cart::y][Cart::z][Cart::yyy];
R[Cart::x][Cart::zz][Cart::yyy]=R[Cart::xz][Cart::z][Cart::yyy]+amb2*R[Cart::x][Cart::z][Cart::yyy];
R[Cart::z][Cart::zz][Cart::yyy]=R[Cart::zz][Cart::z][Cart::yyy]+amb2*R[Cart::z][Cart::z][Cart::yyy];
R[Cart::y][Cart::yy][Cart::xyy]=R[Cart::yy][Cart::y][Cart::xyy]+amb1*R[Cart::y][Cart::y][Cart::xyy];
R[Cart::x][Cart::yy][Cart::xyy]=R[Cart::xy][Cart::y][Cart::xyy]+amb1*R[Cart::x][Cart::y][Cart::xyy];
R[Cart::z][Cart::yy][Cart::xyy]=R[Cart::yz][Cart::y][Cart::xyy]+amb1*R[Cart::z][Cart::y][Cart::xyy];
R[Cart::y][Cart::xy][Cart::xyy]=R[Cart::xy][Cart::y][Cart::xyy]+amb0*R[Cart::y][Cart::y][Cart::xyy];
R[Cart::x][Cart::xy][Cart::xyy]=R[Cart::xx][Cart::y][Cart::xyy]+amb0*R[Cart::x][Cart::y][Cart::xyy];
R[Cart::z][Cart::xy][Cart::xyy]=R[Cart::xz][Cart::y][Cart::xyy]+amb0*R[Cart::z][Cart::y][Cart::xyy];
R[Cart::y][Cart::yz][Cart::xyy]=R[Cart::yy][Cart::z][Cart::xyy]+amb1*R[Cart::y][Cart::z][Cart::xyy];
R[Cart::x][Cart::yz][Cart::xyy]=R[Cart::xy][Cart::z][Cart::xyy]+amb1*R[Cart::x][Cart::z][Cart::xyy];
R[Cart::z][Cart::yz][Cart::xyy]=R[Cart::yz][Cart::z][Cart::xyy]+amb1*R[Cart::z][Cart::z][Cart::xyy];
R[Cart::y][Cart::xx][Cart::xyy]=R[Cart::xy][Cart::x][Cart::xyy]+amb0*R[Cart::y][Cart::x][Cart::xyy];
R[Cart::x][Cart::xx][Cart::xyy]=R[Cart::xx][Cart::x][Cart::xyy]+amb0*R[Cart::x][Cart::x][Cart::xyy];
R[Cart::z][Cart::xx][Cart::xyy]=R[Cart::xz][Cart::x][Cart::xyy]+amb0*R[Cart::z][Cart::x][Cart::xyy];
R[Cart::y][Cart::xz][Cart::xyy]=R[Cart::xy][Cart::z][Cart::xyy]+amb0*R[Cart::y][Cart::z][Cart::xyy];
R[Cart::x][Cart::xz][Cart::xyy]=R[Cart::xx][Cart::z][Cart::xyy]+amb0*R[Cart::x][Cart::z][Cart::xyy];
R[Cart::z][Cart::xz][Cart::xyy]=R[Cart::xz][Cart::z][Cart::xyy]+amb0*R[Cart::z][Cart::z][Cart::xyy];
R[Cart::y][Cart::zz][Cart::xyy]=R[Cart::yz][Cart::z][Cart::xyy]+amb2*R[Cart::y][Cart::z][Cart::xyy];
R[Cart::x][Cart::zz][Cart::xyy]=R[Cart::xz][Cart::z][Cart::xyy]+amb2*R[Cart::x][Cart::z][Cart::xyy];
R[Cart::z][Cart::zz][Cart::xyy]=R[Cart::zz][Cart::z][Cart::xyy]+amb2*R[Cart::z][Cart::z][Cart::xyy];
R[Cart::y][Cart::yy][Cart::yyz]=R[Cart::yy][Cart::y][Cart::yyz]+amb1*R[Cart::y][Cart::y][Cart::yyz];
R[Cart::x][Cart::yy][Cart::yyz]=R[Cart::xy][Cart::y][Cart::yyz]+amb1*R[Cart::x][Cart::y][Cart::yyz];
R[Cart::z][Cart::yy][Cart::yyz]=R[Cart::yz][Cart::y][Cart::yyz]+amb1*R[Cart::z][Cart::y][Cart::yyz];
R[Cart::y][Cart::xy][Cart::yyz]=R[Cart::xy][Cart::y][Cart::yyz]+amb0*R[Cart::y][Cart::y][Cart::yyz];
R[Cart::x][Cart::xy][Cart::yyz]=R[Cart::xx][Cart::y][Cart::yyz]+amb0*R[Cart::x][Cart::y][Cart::yyz];
R[Cart::z][Cart::xy][Cart::yyz]=R[Cart::xz][Cart::y][Cart::yyz]+amb0*R[Cart::z][Cart::y][Cart::yyz];
R[Cart::y][Cart::yz][Cart::yyz]=R[Cart::yy][Cart::z][Cart::yyz]+amb1*R[Cart::y][Cart::z][Cart::yyz];
R[Cart::x][Cart::yz][Cart::yyz]=R[Cart::xy][Cart::z][Cart::yyz]+amb1*R[Cart::x][Cart::z][Cart::yyz];
R[Cart::z][Cart::yz][Cart::yyz]=R[Cart::yz][Cart::z][Cart::yyz]+amb1*R[Cart::z][Cart::z][Cart::yyz];
R[Cart::y][Cart::xx][Cart::yyz]=R[Cart::xy][Cart::x][Cart::yyz]+amb0*R[Cart::y][Cart::x][Cart::yyz];
R[Cart::x][Cart::xx][Cart::yyz]=R[Cart::xx][Cart::x][Cart::yyz]+amb0*R[Cart::x][Cart::x][Cart::yyz];
R[Cart::z][Cart::xx][Cart::yyz]=R[Cart::xz][Cart::x][Cart::yyz]+amb0*R[Cart::z][Cart::x][Cart::yyz];
R[Cart::y][Cart::xz][Cart::yyz]=R[Cart::xy][Cart::z][Cart::yyz]+amb0*R[Cart::y][Cart::z][Cart::yyz];
R[Cart::x][Cart::xz][Cart::yyz]=R[Cart::xx][Cart::z][Cart::yyz]+amb0*R[Cart::x][Cart::z][Cart::yyz];
R[Cart::z][Cart::xz][Cart::yyz]=R[Cart::xz][Cart::z][Cart::yyz]+amb0*R[Cart::z][Cart::z][Cart::yyz];
R[Cart::y][Cart::zz][Cart::yyz]=R[Cart::yz][Cart::z][Cart::yyz]+amb2*R[Cart::y][Cart::z][Cart::yyz];
R[Cart::x][Cart::zz][Cart::yyz]=R[Cart::xz][Cart::z][Cart::yyz]+amb2*R[Cart::x][Cart::z][Cart::yyz];
R[Cart::z][Cart::zz][Cart::yyz]=R[Cart::zz][Cart::z][Cart::yyz]+amb2*R[Cart::z][Cart::z][Cart::yyz];
R[Cart::y][Cart::yy][Cart::xxy]=R[Cart::yy][Cart::y][Cart::xxy]+amb1*R[Cart::y][Cart::y][Cart::xxy];
R[Cart::x][Cart::yy][Cart::xxy]=R[Cart::xy][Cart::y][Cart::xxy]+amb1*R[Cart::x][Cart::y][Cart::xxy];
R[Cart::z][Cart::yy][Cart::xxy]=R[Cart::yz][Cart::y][Cart::xxy]+amb1*R[Cart::z][Cart::y][Cart::xxy];
R[Cart::y][Cart::xy][Cart::xxy]=R[Cart::xy][Cart::y][Cart::xxy]+amb0*R[Cart::y][Cart::y][Cart::xxy];
R[Cart::x][Cart::xy][Cart::xxy]=R[Cart::xx][Cart::y][Cart::xxy]+amb0*R[Cart::x][Cart::y][Cart::xxy];
R[Cart::z][Cart::xy][Cart::xxy]=R[Cart::xz][Cart::y][Cart::xxy]+amb0*R[Cart::z][Cart::y][Cart::xxy];
R[Cart::y][Cart::yz][Cart::xxy]=R[Cart::yy][Cart::z][Cart::xxy]+amb1*R[Cart::y][Cart::z][Cart::xxy];
R[Cart::x][Cart::yz][Cart::xxy]=R[Cart::xy][Cart::z][Cart::xxy]+amb1*R[Cart::x][Cart::z][Cart::xxy];
R[Cart::z][Cart::yz][Cart::xxy]=R[Cart::yz][Cart::z][Cart::xxy]+amb1*R[Cart::z][Cart::z][Cart::xxy];
R[Cart::y][Cart::xx][Cart::xxy]=R[Cart::xy][Cart::x][Cart::xxy]+amb0*R[Cart::y][Cart::x][Cart::xxy];
R[Cart::x][Cart::xx][Cart::xxy]=R[Cart::xx][Cart::x][Cart::xxy]+amb0*R[Cart::x][Cart::x][Cart::xxy];
R[Cart::z][Cart::xx][Cart::xxy]=R[Cart::xz][Cart::x][Cart::xxy]+amb0*R[Cart::z][Cart::x][Cart::xxy];
R[Cart::y][Cart::xz][Cart::xxy]=R[Cart::xy][Cart::z][Cart::xxy]+amb0*R[Cart::y][Cart::z][Cart::xxy];
R[Cart::x][Cart::xz][Cart::xxy]=R[Cart::xx][Cart::z][Cart::xxy]+amb0*R[Cart::x][Cart::z][Cart::xxy];
R[Cart::z][Cart::xz][Cart::xxy]=R[Cart::xz][Cart::z][Cart::xxy]+amb0*R[Cart::z][Cart::z][Cart::xxy];
R[Cart::y][Cart::zz][Cart::xxy]=R[Cart::yz][Cart::z][Cart::xxy]+amb2*R[Cart::y][Cart::z][Cart::xxy];
R[Cart::x][Cart::zz][Cart::xxy]=R[Cart::xz][Cart::z][Cart::xxy]+amb2*R[Cart::x][Cart::z][Cart::xxy];
R[Cart::z][Cart::zz][Cart::xxy]=R[Cart::zz][Cart::z][Cart::xxy]+amb2*R[Cart::z][Cart::z][Cart::xxy];
R[Cart::y][Cart::yy][Cart::xyz]=R[Cart::yy][Cart::y][Cart::xyz]+amb1*R[Cart::y][Cart::y][Cart::xyz];
R[Cart::x][Cart::yy][Cart::xyz]=R[Cart::xy][Cart::y][Cart::xyz]+amb1*R[Cart::x][Cart::y][Cart::xyz];
R[Cart::z][Cart::yy][Cart::xyz]=R[Cart::yz][Cart::y][Cart::xyz]+amb1*R[Cart::z][Cart::y][Cart::xyz];
R[Cart::y][Cart::xy][Cart::xyz]=R[Cart::xy][Cart::y][Cart::xyz]+amb0*R[Cart::y][Cart::y][Cart::xyz];
R[Cart::x][Cart::xy][Cart::xyz]=R[Cart::xx][Cart::y][Cart::xyz]+amb0*R[Cart::x][Cart::y][Cart::xyz];
R[Cart::z][Cart::xy][Cart::xyz]=R[Cart::xz][Cart::y][Cart::xyz]+amb0*R[Cart::z][Cart::y][Cart::xyz];
R[Cart::y][Cart::yz][Cart::xyz]=R[Cart::yy][Cart::z][Cart::xyz]+amb1*R[Cart::y][Cart::z][Cart::xyz];
R[Cart::x][Cart::yz][Cart::xyz]=R[Cart::xy][Cart::z][Cart::xyz]+amb1*R[Cart::x][Cart::z][Cart::xyz];
R[Cart::z][Cart::yz][Cart::xyz]=R[Cart::yz][Cart::z][Cart::xyz]+amb1*R[Cart::z][Cart::z][Cart::xyz];
R[Cart::y][Cart::xx][Cart::xyz]=R[Cart::xy][Cart::x][Cart::xyz]+amb0*R[Cart::y][Cart::x][Cart::xyz];
R[Cart::x][Cart::xx][Cart::xyz]=R[Cart::xx][Cart::x][Cart::xyz]+amb0*R[Cart::x][Cart::x][Cart::xyz];
R[Cart::z][Cart::xx][Cart::xyz]=R[Cart::xz][Cart::x][Cart::xyz]+amb0*R[Cart::z][Cart::x][Cart::xyz];
R[Cart::y][Cart::xz][Cart::xyz]=R[Cart::xy][Cart::z][Cart::xyz]+amb0*R[Cart::y][Cart::z][Cart::xyz];
R[Cart::x][Cart::xz][Cart::xyz]=R[Cart::xx][Cart::z][Cart::xyz]+amb0*R[Cart::x][Cart::z][Cart::xyz];
R[Cart::z][Cart::xz][Cart::xyz]=R[Cart::xz][Cart::z][Cart::xyz]+amb0*R[Cart::z][Cart::z][Cart::xyz];
R[Cart::y][Cart::zz][Cart::xyz]=R[Cart::yz][Cart::z][Cart::xyz]+amb2*R[Cart::y][Cart::z][Cart::xyz];
R[Cart::x][Cart::zz][Cart::xyz]=R[Cart::xz][Cart::z][Cart::xyz]+amb2*R[Cart::x][Cart::z][Cart::xyz];
R[Cart::z][Cart::zz][Cart::xyz]=R[Cart::zz][Cart::z][Cart::xyz]+amb2*R[Cart::z][Cart::z][Cart::xyz];
R[Cart::y][Cart::yy][Cart::yzz]=R[Cart::yy][Cart::y][Cart::yzz]+amb1*R[Cart::y][Cart::y][Cart::yzz];
R[Cart::x][Cart::yy][Cart::yzz]=R[Cart::xy][Cart::y][Cart::yzz]+amb1*R[Cart::x][Cart::y][Cart::yzz];
R[Cart::z][Cart::yy][Cart::yzz]=R[Cart::yz][Cart::y][Cart::yzz]+amb1*R[Cart::z][Cart::y][Cart::yzz];
R[Cart::y][Cart::xy][Cart::yzz]=R[Cart::xy][Cart::y][Cart::yzz]+amb0*R[Cart::y][Cart::y][Cart::yzz];
R[Cart::x][Cart::xy][Cart::yzz]=R[Cart::xx][Cart::y][Cart::yzz]+amb0*R[Cart::x][Cart::y][Cart::yzz];
R[Cart::z][Cart::xy][Cart::yzz]=R[Cart::xz][Cart::y][Cart::yzz]+amb0*R[Cart::z][Cart::y][Cart::yzz];
R[Cart::y][Cart::yz][Cart::yzz]=R[Cart::yy][Cart::z][Cart::yzz]+amb1*R[Cart::y][Cart::z][Cart::yzz];
R[Cart::x][Cart::yz][Cart::yzz]=R[Cart::xy][Cart::z][Cart::yzz]+amb1*R[Cart::x][Cart::z][Cart::yzz];
R[Cart::z][Cart::yz][Cart::yzz]=R[Cart::yz][Cart::z][Cart::yzz]+amb1*R[Cart::z][Cart::z][Cart::yzz];
R[Cart::y][Cart::xx][Cart::yzz]=R[Cart::xy][Cart::x][Cart::yzz]+amb0*R[Cart::y][Cart::x][Cart::yzz];
R[Cart::x][Cart::xx][Cart::yzz]=R[Cart::xx][Cart::x][Cart::yzz]+amb0*R[Cart::x][Cart::x][Cart::yzz];
R[Cart::z][Cart::xx][Cart::yzz]=R[Cart::xz][Cart::x][Cart::yzz]+amb0*R[Cart::z][Cart::x][Cart::yzz];
R[Cart::y][Cart::xz][Cart::yzz]=R[Cart::xy][Cart::z][Cart::yzz]+amb0*R[Cart::y][Cart::z][Cart::yzz];
R[Cart::x][Cart::xz][Cart::yzz]=R[Cart::xx][Cart::z][Cart::yzz]+amb0*R[Cart::x][Cart::z][Cart::yzz];
R[Cart::z][Cart::xz][Cart::yzz]=R[Cart::xz][Cart::z][Cart::yzz]+amb0*R[Cart::z][Cart::z][Cart::yzz];
R[Cart::y][Cart::zz][Cart::yzz]=R[Cart::yz][Cart::z][Cart::yzz]+amb2*R[Cart::y][Cart::z][Cart::yzz];
R[Cart::x][Cart::zz][Cart::yzz]=R[Cart::xz][Cart::z][Cart::yzz]+amb2*R[Cart::x][Cart::z][Cart::yzz];
R[Cart::z][Cart::zz][Cart::yzz]=R[Cart::zz][Cart::z][Cart::yzz]+amb2*R[Cart::z][Cart::z][Cart::yzz];
R[Cart::y][Cart::yy][Cart::xxx]=R[Cart::yy][Cart::y][Cart::xxx]+amb1*R[Cart::y][Cart::y][Cart::xxx];
R[Cart::x][Cart::yy][Cart::xxx]=R[Cart::xy][Cart::y][Cart::xxx]+amb1*R[Cart::x][Cart::y][Cart::xxx];
R[Cart::z][Cart::yy][Cart::xxx]=R[Cart::yz][Cart::y][Cart::xxx]+amb1*R[Cart::z][Cart::y][Cart::xxx];
R[Cart::y][Cart::xy][Cart::xxx]=R[Cart::xy][Cart::y][Cart::xxx]+amb0*R[Cart::y][Cart::y][Cart::xxx];
R[Cart::x][Cart::xy][Cart::xxx]=R[Cart::xx][Cart::y][Cart::xxx]+amb0*R[Cart::x][Cart::y][Cart::xxx];
R[Cart::z][Cart::xy][Cart::xxx]=R[Cart::xz][Cart::y][Cart::xxx]+amb0*R[Cart::z][Cart::y][Cart::xxx];
R[Cart::y][Cart::yz][Cart::xxx]=R[Cart::yy][Cart::z][Cart::xxx]+amb1*R[Cart::y][Cart::z][Cart::xxx];
R[Cart::x][Cart::yz][Cart::xxx]=R[Cart::xy][Cart::z][Cart::xxx]+amb1*R[Cart::x][Cart::z][Cart::xxx];
R[Cart::z][Cart::yz][Cart::xxx]=R[Cart::yz][Cart::z][Cart::xxx]+amb1*R[Cart::z][Cart::z][Cart::xxx];
R[Cart::y][Cart::xx][Cart::xxx]=R[Cart::xy][Cart::x][Cart::xxx]+amb0*R[Cart::y][Cart::x][Cart::xxx];
R[Cart::x][Cart::xx][Cart::xxx]=R[Cart::xx][Cart::x][Cart::xxx]+amb0*R[Cart::x][Cart::x][Cart::xxx];
R[Cart::z][Cart::xx][Cart::xxx]=R[Cart::xz][Cart::x][Cart::xxx]+amb0*R[Cart::z][Cart::x][Cart::xxx];
R[Cart::y][Cart::xz][Cart::xxx]=R[Cart::xy][Cart::z][Cart::xxx]+amb0*R[Cart::y][Cart::z][Cart::xxx];
R[Cart::x][Cart::xz][Cart::xxx]=R[Cart::xx][Cart::z][Cart::xxx]+amb0*R[Cart::x][Cart::z][Cart::xxx];
R[Cart::z][Cart::xz][Cart::xxx]=R[Cart::xz][Cart::z][Cart::xxx]+amb0*R[Cart::z][Cart::z][Cart::xxx];
R[Cart::y][Cart::zz][Cart::xxx]=R[Cart::yz][Cart::z][Cart::xxx]+amb2*R[Cart::y][Cart::z][Cart::xxx];
R[Cart::x][Cart::zz][Cart::xxx]=R[Cart::xz][Cart::z][Cart::xxx]+amb2*R[Cart::x][Cart::z][Cart::xxx];
R[Cart::z][Cart::zz][Cart::xxx]=R[Cart::zz][Cart::z][Cart::xxx]+amb2*R[Cart::z][Cart::z][Cart::xxx];
R[Cart::y][Cart::yy][Cart::xxz]=R[Cart::yy][Cart::y][Cart::xxz]+amb1*R[Cart::y][Cart::y][Cart::xxz];
R[Cart::x][Cart::yy][Cart::xxz]=R[Cart::xy][Cart::y][Cart::xxz]+amb1*R[Cart::x][Cart::y][Cart::xxz];
R[Cart::z][Cart::yy][Cart::xxz]=R[Cart::yz][Cart::y][Cart::xxz]+amb1*R[Cart::z][Cart::y][Cart::xxz];
R[Cart::y][Cart::xy][Cart::xxz]=R[Cart::xy][Cart::y][Cart::xxz]+amb0*R[Cart::y][Cart::y][Cart::xxz];
R[Cart::x][Cart::xy][Cart::xxz]=R[Cart::xx][Cart::y][Cart::xxz]+amb0*R[Cart::x][Cart::y][Cart::xxz];
R[Cart::z][Cart::xy][Cart::xxz]=R[Cart::xz][Cart::y][Cart::xxz]+amb0*R[Cart::z][Cart::y][Cart::xxz];
R[Cart::y][Cart::yz][Cart::xxz]=R[Cart::yy][Cart::z][Cart::xxz]+amb1*R[Cart::y][Cart::z][Cart::xxz];
R[Cart::x][Cart::yz][Cart::xxz]=R[Cart::xy][Cart::z][Cart::xxz]+amb1*R[Cart::x][Cart::z][Cart::xxz];
R[Cart::z][Cart::yz][Cart::xxz]=R[Cart::yz][Cart::z][Cart::xxz]+amb1*R[Cart::z][Cart::z][Cart::xxz];
R[Cart::y][Cart::xx][Cart::xxz]=R[Cart::xy][Cart::x][Cart::xxz]+amb0*R[Cart::y][Cart::x][Cart::xxz];
R[Cart::x][Cart::xx][Cart::xxz]=R[Cart::xx][Cart::x][Cart::xxz]+amb0*R[Cart::x][Cart::x][Cart::xxz];
R[Cart::z][Cart::xx][Cart::xxz]=R[Cart::xz][Cart::x][Cart::xxz]+amb0*R[Cart::z][Cart::x][Cart::xxz];
R[Cart::y][Cart::xz][Cart::xxz]=R[Cart::xy][Cart::z][Cart::xxz]+amb0*R[Cart::y][Cart::z][Cart::xxz];
R[Cart::x][Cart::xz][Cart::xxz]=R[Cart::xx][Cart::z][Cart::xxz]+amb0*R[Cart::x][Cart::z][Cart::xxz];
R[Cart::z][Cart::xz][Cart::xxz]=R[Cart::xz][Cart::z][Cart::xxz]+amb0*R[Cart::z][Cart::z][Cart::xxz];
R[Cart::y][Cart::zz][Cart::xxz]=R[Cart::yz][Cart::z][Cart::xxz]+amb2*R[Cart::y][Cart::z][Cart::xxz];
R[Cart::x][Cart::zz][Cart::xxz]=R[Cart::xz][Cart::z][Cart::xxz]+amb2*R[Cart::x][Cart::z][Cart::xxz];
R[Cart::z][Cart::zz][Cart::xxz]=R[Cart::zz][Cart::z][Cart::xxz]+amb2*R[Cart::z][Cart::z][Cart::xxz];
R[Cart::y][Cart::yy][Cart::xzz]=R[Cart::yy][Cart::y][Cart::xzz]+amb1*R[Cart::y][Cart::y][Cart::xzz];
R[Cart::x][Cart::yy][Cart::xzz]=R[Cart::xy][Cart::y][Cart::xzz]+amb1*R[Cart::x][Cart::y][Cart::xzz];
R[Cart::z][Cart::yy][Cart::xzz]=R[Cart::yz][Cart::y][Cart::xzz]+amb1*R[Cart::z][Cart::y][Cart::xzz];
R[Cart::y][Cart::xy][Cart::xzz]=R[Cart::xy][Cart::y][Cart::xzz]+amb0*R[Cart::y][Cart::y][Cart::xzz];
R[Cart::x][Cart::xy][Cart::xzz]=R[Cart::xx][Cart::y][Cart::xzz]+amb0*R[Cart::x][Cart::y][Cart::xzz];
R[Cart::z][Cart::xy][Cart::xzz]=R[Cart::xz][Cart::y][Cart::xzz]+amb0*R[Cart::z][Cart::y][Cart::xzz];
R[Cart::y][Cart::yz][Cart::xzz]=R[Cart::yy][Cart::z][Cart::xzz]+amb1*R[Cart::y][Cart::z][Cart::xzz];
R[Cart::x][Cart::yz][Cart::xzz]=R[Cart::xy][Cart::z][Cart::xzz]+amb1*R[Cart::x][Cart::z][Cart::xzz];
R[Cart::z][Cart::yz][Cart::xzz]=R[Cart::yz][Cart::z][Cart::xzz]+amb1*R[Cart::z][Cart::z][Cart::xzz];
R[Cart::y][Cart::xx][Cart::xzz]=R[Cart::xy][Cart::x][Cart::xzz]+amb0*R[Cart::y][Cart::x][Cart::xzz];
R[Cart::x][Cart::xx][Cart::xzz]=R[Cart::xx][Cart::x][Cart::xzz]+amb0*R[Cart::x][Cart::x][Cart::xzz];
R[Cart::z][Cart::xx][Cart::xzz]=R[Cart::xz][Cart::x][Cart::xzz]+amb0*R[Cart::z][Cart::x][Cart::xzz];
R[Cart::y][Cart::xz][Cart::xzz]=R[Cart::xy][Cart::z][Cart::xzz]+amb0*R[Cart::y][Cart::z][Cart::xzz];
R[Cart::x][Cart::xz][Cart::xzz]=R[Cart::xx][Cart::z][Cart::xzz]+amb0*R[Cart::x][Cart::z][Cart::xzz];
R[Cart::z][Cart::xz][Cart::xzz]=R[Cart::xz][Cart::z][Cart::xzz]+amb0*R[Cart::z][Cart::z][Cart::xzz];
R[Cart::y][Cart::zz][Cart::xzz]=R[Cart::yz][Cart::z][Cart::xzz]+amb2*R[Cart::y][Cart::z][Cart::xzz];
R[Cart::x][Cart::zz][Cart::xzz]=R[Cart::xz][Cart::z][Cart::xzz]+amb2*R[Cart::x][Cart::z][Cart::xzz];
R[Cart::z][Cart::zz][Cart::xzz]=R[Cart::zz][Cart::z][Cart::xzz]+amb2*R[Cart::z][Cart::z][Cart::xzz];
R[Cart::y][Cart::yy][Cart::zzz]=R[Cart::yy][Cart::y][Cart::zzz]+amb1*R[Cart::y][Cart::y][Cart::zzz];
R[Cart::x][Cart::yy][Cart::zzz]=R[Cart::xy][Cart::y][Cart::zzz]+amb1*R[Cart::x][Cart::y][Cart::zzz];
R[Cart::z][Cart::yy][Cart::zzz]=R[Cart::yz][Cart::y][Cart::zzz]+amb1*R[Cart::z][Cart::y][Cart::zzz];
R[Cart::y][Cart::xy][Cart::zzz]=R[Cart::xy][Cart::y][Cart::zzz]+amb0*R[Cart::y][Cart::y][Cart::zzz];
R[Cart::x][Cart::xy][Cart::zzz]=R[Cart::xx][Cart::y][Cart::zzz]+amb0*R[Cart::x][Cart::y][Cart::zzz];
R[Cart::z][Cart::xy][Cart::zzz]=R[Cart::xz][Cart::y][Cart::zzz]+amb0*R[Cart::z][Cart::y][Cart::zzz];
R[Cart::y][Cart::yz][Cart::zzz]=R[Cart::yy][Cart::z][Cart::zzz]+amb1*R[Cart::y][Cart::z][Cart::zzz];
R[Cart::x][Cart::yz][Cart::zzz]=R[Cart::xy][Cart::z][Cart::zzz]+amb1*R[Cart::x][Cart::z][Cart::zzz];
R[Cart::z][Cart::yz][Cart::zzz]=R[Cart::yz][Cart::z][Cart::zzz]+amb1*R[Cart::z][Cart::z][Cart::zzz];
R[Cart::y][Cart::xx][Cart::zzz]=R[Cart::xy][Cart::x][Cart::zzz]+amb0*R[Cart::y][Cart::x][Cart::zzz];
R[Cart::x][Cart::xx][Cart::zzz]=R[Cart::xx][Cart::x][Cart::zzz]+amb0*R[Cart::x][Cart::x][Cart::zzz];
R[Cart::z][Cart::xx][Cart::zzz]=R[Cart::xz][Cart::x][Cart::zzz]+amb0*R[Cart::z][Cart::x][Cart::zzz];
R[Cart::y][Cart::xz][Cart::zzz]=R[Cart::xy][Cart::z][Cart::zzz]+amb0*R[Cart::y][Cart::z][Cart::zzz];
R[Cart::x][Cart::xz][Cart::zzz]=R[Cart::xx][Cart::z][Cart::zzz]+amb0*R[Cart::x][Cart::z][Cart::zzz];
R[Cart::z][Cart::xz][Cart::zzz]=R[Cart::xz][Cart::z][Cart::zzz]+amb0*R[Cart::z][Cart::z][Cart::zzz];
R[Cart::y][Cart::zz][Cart::zzz]=R[Cart::yz][Cart::z][Cart::zzz]+amb2*R[Cart::y][Cart::z][Cart::zzz];
R[Cart::x][Cart::zz][Cart::zzz]=R[Cart::xz][Cart::z][Cart::zzz]+amb2*R[Cart::x][Cart::z][Cart::zzz];
R[Cart::z][Cart::zz][Cart::zzz]=R[Cart::zz][Cart::z][Cart::zzz]+amb2*R[Cart::z][Cart::z][Cart::zzz];
}
//------------------------------------------------------

//Integral d - d - f
if (_lmax_beta>1 && _lmax_alpha>1 && _lmax_gamma>2){

R[Cart::yy][Cart::yy][Cart::yyy]=R[Cart::yyy][Cart::y][Cart::yyy]+amb1*R[Cart::yy][Cart::y][Cart::yyy];
R[Cart::xy][Cart::yy][Cart::yyy]=R[Cart::xyy][Cart::y][Cart::yyy]+amb1*R[Cart::xy][Cart::y][Cart::yyy];
R[Cart::yz][Cart::yy][Cart::yyy]=R[Cart::yyz][Cart::y][Cart::yyy]+amb1*R[Cart::yz][Cart::y][Cart::yyy];
R[Cart::xx][Cart::yy][Cart::yyy]=R[Cart::xxy][Cart::y][Cart::yyy]+amb1*R[Cart::xx][Cart::y][Cart::yyy];
R[Cart::xz][Cart::yy][Cart::yyy]=R[Cart::xyz][Cart::y][Cart::yyy]+amb1*R[Cart::xz][Cart::y][Cart::yyy];
R[Cart::zz][Cart::yy][Cart::yyy]=R[Cart::yzz][Cart::y][Cart::yyy]+amb1*R[Cart::zz][Cart::y][Cart::yyy];
R[Cart::yy][Cart::xy][Cart::yyy]=R[Cart::xyy][Cart::y][Cart::yyy]+amb0*R[Cart::yy][Cart::y][Cart::yyy];
R[Cart::xy][Cart::xy][Cart::yyy]=R[Cart::xxy][Cart::y][Cart::yyy]+amb0*R[Cart::xy][Cart::y][Cart::yyy];
R[Cart::yz][Cart::xy][Cart::yyy]=R[Cart::xyz][Cart::y][Cart::yyy]+amb0*R[Cart::yz][Cart::y][Cart::yyy];
R[Cart::xx][Cart::xy][Cart::yyy]=R[Cart::xxx][Cart::y][Cart::yyy]+amb0*R[Cart::xx][Cart::y][Cart::yyy];
R[Cart::xz][Cart::xy][Cart::yyy]=R[Cart::xxz][Cart::y][Cart::yyy]+amb0*R[Cart::xz][Cart::y][Cart::yyy];
R[Cart::zz][Cart::xy][Cart::yyy]=R[Cart::xzz][Cart::y][Cart::yyy]+amb0*R[Cart::zz][Cart::y][Cart::yyy];
R[Cart::yy][Cart::yz][Cart::yyy]=R[Cart::yyy][Cart::z][Cart::yyy]+amb1*R[Cart::yy][Cart::z][Cart::yyy];
R[Cart::xy][Cart::yz][Cart::yyy]=R[Cart::xyy][Cart::z][Cart::yyy]+amb1*R[Cart::xy][Cart::z][Cart::yyy];
R[Cart::yz][Cart::yz][Cart::yyy]=R[Cart::yyz][Cart::z][Cart::yyy]+amb1*R[Cart::yz][Cart::z][Cart::yyy];
R[Cart::xx][Cart::yz][Cart::yyy]=R[Cart::xxy][Cart::z][Cart::yyy]+amb1*R[Cart::xx][Cart::z][Cart::yyy];
R[Cart::xz][Cart::yz][Cart::yyy]=R[Cart::xyz][Cart::z][Cart::yyy]+amb1*R[Cart::xz][Cart::z][Cart::yyy];
R[Cart::zz][Cart::yz][Cart::yyy]=R[Cart::yzz][Cart::z][Cart::yyy]+amb1*R[Cart::zz][Cart::z][Cart::yyy];
R[Cart::yy][Cart::xx][Cart::yyy]=R[Cart::xyy][Cart::x][Cart::yyy]+amb0*R[Cart::yy][Cart::x][Cart::yyy];
R[Cart::xy][Cart::xx][Cart::yyy]=R[Cart::xxy][Cart::x][Cart::yyy]+amb0*R[Cart::xy][Cart::x][Cart::yyy];
R[Cart::yz][Cart::xx][Cart::yyy]=R[Cart::xyz][Cart::x][Cart::yyy]+amb0*R[Cart::yz][Cart::x][Cart::yyy];
R[Cart::xx][Cart::xx][Cart::yyy]=R[Cart::xxx][Cart::x][Cart::yyy]+amb0*R[Cart::xx][Cart::x][Cart::yyy];
R[Cart::xz][Cart::xx][Cart::yyy]=R[Cart::xxz][Cart::x][Cart::yyy]+amb0*R[Cart::xz][Cart::x][Cart::yyy];
R[Cart::zz][Cart::xx][Cart::yyy]=R[Cart::xzz][Cart::x][Cart::yyy]+amb0*R[Cart::zz][Cart::x][Cart::yyy];
R[Cart::yy][Cart::xz][Cart::yyy]=R[Cart::xyy][Cart::z][Cart::yyy]+amb0*R[Cart::yy][Cart::z][Cart::yyy];
R[Cart::xy][Cart::xz][Cart::yyy]=R[Cart::xxy][Cart::z][Cart::yyy]+amb0*R[Cart::xy][Cart::z][Cart::yyy];
R[Cart::yz][Cart::xz][Cart::yyy]=R[Cart::xyz][Cart::z][Cart::yyy]+amb0*R[Cart::yz][Cart::z][Cart::yyy];
R[Cart::xx][Cart::xz][Cart::yyy]=R[Cart::xxx][Cart::z][Cart::yyy]+amb0*R[Cart::xx][Cart::z][Cart::yyy];
R[Cart::xz][Cart::xz][Cart::yyy]=R[Cart::xxz][Cart::z][Cart::yyy]+amb0*R[Cart::xz][Cart::z][Cart::yyy];
R[Cart::zz][Cart::xz][Cart::yyy]=R[Cart::xzz][Cart::z][Cart::yyy]+amb0*R[Cart::zz][Cart::z][Cart::yyy];
R[Cart::yy][Cart::zz][Cart::yyy]=R[Cart::yyz][Cart::z][Cart::yyy]+amb2*R[Cart::yy][Cart::z][Cart::yyy];
R[Cart::xy][Cart::zz][Cart::yyy]=R[Cart::xyz][Cart::z][Cart::yyy]+amb2*R[Cart::xy][Cart::z][Cart::yyy];
R[Cart::yz][Cart::zz][Cart::yyy]=R[Cart::yzz][Cart::z][Cart::yyy]+amb2*R[Cart::yz][Cart::z][Cart::yyy];
R[Cart::xx][Cart::zz][Cart::yyy]=R[Cart::xxz][Cart::z][Cart::yyy]+amb2*R[Cart::xx][Cart::z][Cart::yyy];
R[Cart::xz][Cart::zz][Cart::yyy]=R[Cart::xzz][Cart::z][Cart::yyy]+amb2*R[Cart::xz][Cart::z][Cart::yyy];
R[Cart::zz][Cart::zz][Cart::yyy]=R[Cart::zzz][Cart::z][Cart::yyy]+amb2*R[Cart::zz][Cart::z][Cart::yyy];
R[Cart::yy][Cart::yy][Cart::xyy]=R[Cart::yyy][Cart::y][Cart::xyy]+amb1*R[Cart::yy][Cart::y][Cart::xyy];
R[Cart::xy][Cart::yy][Cart::xyy]=R[Cart::xyy][Cart::y][Cart::xyy]+amb1*R[Cart::xy][Cart::y][Cart::xyy];
R[Cart::yz][Cart::yy][Cart::xyy]=R[Cart::yyz][Cart::y][Cart::xyy]+amb1*R[Cart::yz][Cart::y][Cart::xyy];
R[Cart::xx][Cart::yy][Cart::xyy]=R[Cart::xxy][Cart::y][Cart::xyy]+amb1*R[Cart::xx][Cart::y][Cart::xyy];
R[Cart::xz][Cart::yy][Cart::xyy]=R[Cart::xyz][Cart::y][Cart::xyy]+amb1*R[Cart::xz][Cart::y][Cart::xyy];
R[Cart::zz][Cart::yy][Cart::xyy]=R[Cart::yzz][Cart::y][Cart::xyy]+amb1*R[Cart::zz][Cart::y][Cart::xyy];
R[Cart::yy][Cart::xy][Cart::xyy]=R[Cart::xyy][Cart::y][Cart::xyy]+amb0*R[Cart::yy][Cart::y][Cart::xyy];
R[Cart::xy][Cart::xy][Cart::xyy]=R[Cart::xxy][Cart::y][Cart::xyy]+amb0*R[Cart::xy][Cart::y][Cart::xyy];
R[Cart::yz][Cart::xy][Cart::xyy]=R[Cart::xyz][Cart::y][Cart::xyy]+amb0*R[Cart::yz][Cart::y][Cart::xyy];
R[Cart::xx][Cart::xy][Cart::xyy]=R[Cart::xxx][Cart::y][Cart::xyy]+amb0*R[Cart::xx][Cart::y][Cart::xyy];
R[Cart::xz][Cart::xy][Cart::xyy]=R[Cart::xxz][Cart::y][Cart::xyy]+amb0*R[Cart::xz][Cart::y][Cart::xyy];
R[Cart::zz][Cart::xy][Cart::xyy]=R[Cart::xzz][Cart::y][Cart::xyy]+amb0*R[Cart::zz][Cart::y][Cart::xyy];
R[Cart::yy][Cart::yz][Cart::xyy]=R[Cart::yyy][Cart::z][Cart::xyy]+amb1*R[Cart::yy][Cart::z][Cart::xyy];
R[Cart::xy][Cart::yz][Cart::xyy]=R[Cart::xyy][Cart::z][Cart::xyy]+amb1*R[Cart::xy][Cart::z][Cart::xyy];
R[Cart::yz][Cart::yz][Cart::xyy]=R[Cart::yyz][Cart::z][Cart::xyy]+amb1*R[Cart::yz][Cart::z][Cart::xyy];
R[Cart::xx][Cart::yz][Cart::xyy]=R[Cart::xxy][Cart::z][Cart::xyy]+amb1*R[Cart::xx][Cart::z][Cart::xyy];
R[Cart::xz][Cart::yz][Cart::xyy]=R[Cart::xyz][Cart::z][Cart::xyy]+amb1*R[Cart::xz][Cart::z][Cart::xyy];
R[Cart::zz][Cart::yz][Cart::xyy]=R[Cart::yzz][Cart::z][Cart::xyy]+amb1*R[Cart::zz][Cart::z][Cart::xyy];
R[Cart::yy][Cart::xx][Cart::xyy]=R[Cart::xyy][Cart::x][Cart::xyy]+amb0*R[Cart::yy][Cart::x][Cart::xyy];
R[Cart::xy][Cart::xx][Cart::xyy]=R[Cart::xxy][Cart::x][Cart::xyy]+amb0*R[Cart::xy][Cart::x][Cart::xyy];
R[Cart::yz][Cart::xx][Cart::xyy]=R[Cart::xyz][Cart::x][Cart::xyy]+amb0*R[Cart::yz][Cart::x][Cart::xyy];
R[Cart::xx][Cart::xx][Cart::xyy]=R[Cart::xxx][Cart::x][Cart::xyy]+amb0*R[Cart::xx][Cart::x][Cart::xyy];
R[Cart::xz][Cart::xx][Cart::xyy]=R[Cart::xxz][Cart::x][Cart::xyy]+amb0*R[Cart::xz][Cart::x][Cart::xyy];
R[Cart::zz][Cart::xx][Cart::xyy]=R[Cart::xzz][Cart::x][Cart::xyy]+amb0*R[Cart::zz][Cart::x][Cart::xyy];
R[Cart::yy][Cart::xz][Cart::xyy]=R[Cart::xyy][Cart::z][Cart::xyy]+amb0*R[Cart::yy][Cart::z][Cart::xyy];
R[Cart::xy][Cart::xz][Cart::xyy]=R[Cart::xxy][Cart::z][Cart::xyy]+amb0*R[Cart::xy][Cart::z][Cart::xyy];
R[Cart::yz][Cart::xz][Cart::xyy]=R[Cart::xyz][Cart::z][Cart::xyy]+amb0*R[Cart::yz][Cart::z][Cart::xyy];
R[Cart::xx][Cart::xz][Cart::xyy]=R[Cart::xxx][Cart::z][Cart::xyy]+amb0*R[Cart::xx][Cart::z][Cart::xyy];
R[Cart::xz][Cart::xz][Cart::xyy]=R[Cart::xxz][Cart::z][Cart::xyy]+amb0*R[Cart::xz][Cart::z][Cart::xyy];
R[Cart::zz][Cart::xz][Cart::xyy]=R[Cart::xzz][Cart::z][Cart::xyy]+amb0*R[Cart::zz][Cart::z][Cart::xyy];
R[Cart::yy][Cart::zz][Cart::xyy]=R[Cart::yyz][Cart::z][Cart::xyy]+amb2*R[Cart::yy][Cart::z][Cart::xyy];
R[Cart::xy][Cart::zz][Cart::xyy]=R[Cart::xyz][Cart::z][Cart::xyy]+amb2*R[Cart::xy][Cart::z][Cart::xyy];
R[Cart::yz][Cart::zz][Cart::xyy]=R[Cart::yzz][Cart::z][Cart::xyy]+amb2*R[Cart::yz][Cart::z][Cart::xyy];
R[Cart::xx][Cart::zz][Cart::xyy]=R[Cart::xxz][Cart::z][Cart::xyy]+amb2*R[Cart::xx][Cart::z][Cart::xyy];
R[Cart::xz][Cart::zz][Cart::xyy]=R[Cart::xzz][Cart::z][Cart::xyy]+amb2*R[Cart::xz][Cart::z][Cart::xyy];
R[Cart::zz][Cart::zz][Cart::xyy]=R[Cart::zzz][Cart::z][Cart::xyy]+amb2*R[Cart::zz][Cart::z][Cart::xyy];
R[Cart::yy][Cart::yy][Cart::yyz]=R[Cart::yyy][Cart::y][Cart::yyz]+amb1*R[Cart::yy][Cart::y][Cart::yyz];
R[Cart::xy][Cart::yy][Cart::yyz]=R[Cart::xyy][Cart::y][Cart::yyz]+amb1*R[Cart::xy][Cart::y][Cart::yyz];
R[Cart::yz][Cart::yy][Cart::yyz]=R[Cart::yyz][Cart::y][Cart::yyz]+amb1*R[Cart::yz][Cart::y][Cart::yyz];
R[Cart::xx][Cart::yy][Cart::yyz]=R[Cart::xxy][Cart::y][Cart::yyz]+amb1*R[Cart::xx][Cart::y][Cart::yyz];
R[Cart::xz][Cart::yy][Cart::yyz]=R[Cart::xyz][Cart::y][Cart::yyz]+amb1*R[Cart::xz][Cart::y][Cart::yyz];
R[Cart::zz][Cart::yy][Cart::yyz]=R[Cart::yzz][Cart::y][Cart::yyz]+amb1*R[Cart::zz][Cart::y][Cart::yyz];
R[Cart::yy][Cart::xy][Cart::yyz]=R[Cart::xyy][Cart::y][Cart::yyz]+amb0*R[Cart::yy][Cart::y][Cart::yyz];
R[Cart::xy][Cart::xy][Cart::yyz]=R[Cart::xxy][Cart::y][Cart::yyz]+amb0*R[Cart::xy][Cart::y][Cart::yyz];
R[Cart::yz][Cart::xy][Cart::yyz]=R[Cart::xyz][Cart::y][Cart::yyz]+amb0*R[Cart::yz][Cart::y][Cart::yyz];
R[Cart::xx][Cart::xy][Cart::yyz]=R[Cart::xxx][Cart::y][Cart::yyz]+amb0*R[Cart::xx][Cart::y][Cart::yyz];
R[Cart::xz][Cart::xy][Cart::yyz]=R[Cart::xxz][Cart::y][Cart::yyz]+amb0*R[Cart::xz][Cart::y][Cart::yyz];
R[Cart::zz][Cart::xy][Cart::yyz]=R[Cart::xzz][Cart::y][Cart::yyz]+amb0*R[Cart::zz][Cart::y][Cart::yyz];
R[Cart::yy][Cart::yz][Cart::yyz]=R[Cart::yyy][Cart::z][Cart::yyz]+amb1*R[Cart::yy][Cart::z][Cart::yyz];
R[Cart::xy][Cart::yz][Cart::yyz]=R[Cart::xyy][Cart::z][Cart::yyz]+amb1*R[Cart::xy][Cart::z][Cart::yyz];
R[Cart::yz][Cart::yz][Cart::yyz]=R[Cart::yyz][Cart::z][Cart::yyz]+amb1*R[Cart::yz][Cart::z][Cart::yyz];
R[Cart::xx][Cart::yz][Cart::yyz]=R[Cart::xxy][Cart::z][Cart::yyz]+amb1*R[Cart::xx][Cart::z][Cart::yyz];
R[Cart::xz][Cart::yz][Cart::yyz]=R[Cart::xyz][Cart::z][Cart::yyz]+amb1*R[Cart::xz][Cart::z][Cart::yyz];
R[Cart::zz][Cart::yz][Cart::yyz]=R[Cart::yzz][Cart::z][Cart::yyz]+amb1*R[Cart::zz][Cart::z][Cart::yyz];
R[Cart::yy][Cart::xx][Cart::yyz]=R[Cart::xyy][Cart::x][Cart::yyz]+amb0*R[Cart::yy][Cart::x][Cart::yyz];
R[Cart::xy][Cart::xx][Cart::yyz]=R[Cart::xxy][Cart::x][Cart::yyz]+amb0*R[Cart::xy][Cart::x][Cart::yyz];
R[Cart::yz][Cart::xx][Cart::yyz]=R[Cart::xyz][Cart::x][Cart::yyz]+amb0*R[Cart::yz][Cart::x][Cart::yyz];
R[Cart::xx][Cart::xx][Cart::yyz]=R[Cart::xxx][Cart::x][Cart::yyz]+amb0*R[Cart::xx][Cart::x][Cart::yyz];
R[Cart::xz][Cart::xx][Cart::yyz]=R[Cart::xxz][Cart::x][Cart::yyz]+amb0*R[Cart::xz][Cart::x][Cart::yyz];
R[Cart::zz][Cart::xx][Cart::yyz]=R[Cart::xzz][Cart::x][Cart::yyz]+amb0*R[Cart::zz][Cart::x][Cart::yyz];
R[Cart::yy][Cart::xz][Cart::yyz]=R[Cart::xyy][Cart::z][Cart::yyz]+amb0*R[Cart::yy][Cart::z][Cart::yyz];
R[Cart::xy][Cart::xz][Cart::yyz]=R[Cart::xxy][Cart::z][Cart::yyz]+amb0*R[Cart::xy][Cart::z][Cart::yyz];
R[Cart::yz][Cart::xz][Cart::yyz]=R[Cart::xyz][Cart::z][Cart::yyz]+amb0*R[Cart::yz][Cart::z][Cart::yyz];
R[Cart::xx][Cart::xz][Cart::yyz]=R[Cart::xxx][Cart::z][Cart::yyz]+amb0*R[Cart::xx][Cart::z][Cart::yyz];
R[Cart::xz][Cart::xz][Cart::yyz]=R[Cart::xxz][Cart::z][Cart::yyz]+amb0*R[Cart::xz][Cart::z][Cart::yyz];
R[Cart::zz][Cart::xz][Cart::yyz]=R[Cart::xzz][Cart::z][Cart::yyz]+amb0*R[Cart::zz][Cart::z][Cart::yyz];
R[Cart::yy][Cart::zz][Cart::yyz]=R[Cart::yyz][Cart::z][Cart::yyz]+amb2*R[Cart::yy][Cart::z][Cart::yyz];
R[Cart::xy][Cart::zz][Cart::yyz]=R[Cart::xyz][Cart::z][Cart::yyz]+amb2*R[Cart::xy][Cart::z][Cart::yyz];
R[Cart::yz][Cart::zz][Cart::yyz]=R[Cart::yzz][Cart::z][Cart::yyz]+amb2*R[Cart::yz][Cart::z][Cart::yyz];
R[Cart::xx][Cart::zz][Cart::yyz]=R[Cart::xxz][Cart::z][Cart::yyz]+amb2*R[Cart::xx][Cart::z][Cart::yyz];
R[Cart::xz][Cart::zz][Cart::yyz]=R[Cart::xzz][Cart::z][Cart::yyz]+amb2*R[Cart::xz][Cart::z][Cart::yyz];
R[Cart::zz][Cart::zz][Cart::yyz]=R[Cart::zzz][Cart::z][Cart::yyz]+amb2*R[Cart::zz][Cart::z][Cart::yyz];
R[Cart::yy][Cart::yy][Cart::xxy]=R[Cart::yyy][Cart::y][Cart::xxy]+amb1*R[Cart::yy][Cart::y][Cart::xxy];
R[Cart::xy][Cart::yy][Cart::xxy]=R[Cart::xyy][Cart::y][Cart::xxy]+amb1*R[Cart::xy][Cart::y][Cart::xxy];
R[Cart::yz][Cart::yy][Cart::xxy]=R[Cart::yyz][Cart::y][Cart::xxy]+amb1*R[Cart::yz][Cart::y][Cart::xxy];
R[Cart::xx][Cart::yy][Cart::xxy]=R[Cart::xxy][Cart::y][Cart::xxy]+amb1*R[Cart::xx][Cart::y][Cart::xxy];
R[Cart::xz][Cart::yy][Cart::xxy]=R[Cart::xyz][Cart::y][Cart::xxy]+amb1*R[Cart::xz][Cart::y][Cart::xxy];
R[Cart::zz][Cart::yy][Cart::xxy]=R[Cart::yzz][Cart::y][Cart::xxy]+amb1*R[Cart::zz][Cart::y][Cart::xxy];
R[Cart::yy][Cart::xy][Cart::xxy]=R[Cart::xyy][Cart::y][Cart::xxy]+amb0*R[Cart::yy][Cart::y][Cart::xxy];
R[Cart::xy][Cart::xy][Cart::xxy]=R[Cart::xxy][Cart::y][Cart::xxy]+amb0*R[Cart::xy][Cart::y][Cart::xxy];
R[Cart::yz][Cart::xy][Cart::xxy]=R[Cart::xyz][Cart::y][Cart::xxy]+amb0*R[Cart::yz][Cart::y][Cart::xxy];
R[Cart::xx][Cart::xy][Cart::xxy]=R[Cart::xxx][Cart::y][Cart::xxy]+amb0*R[Cart::xx][Cart::y][Cart::xxy];
R[Cart::xz][Cart::xy][Cart::xxy]=R[Cart::xxz][Cart::y][Cart::xxy]+amb0*R[Cart::xz][Cart::y][Cart::xxy];
R[Cart::zz][Cart::xy][Cart::xxy]=R[Cart::xzz][Cart::y][Cart::xxy]+amb0*R[Cart::zz][Cart::y][Cart::xxy];
R[Cart::yy][Cart::yz][Cart::xxy]=R[Cart::yyy][Cart::z][Cart::xxy]+amb1*R[Cart::yy][Cart::z][Cart::xxy];
R[Cart::xy][Cart::yz][Cart::xxy]=R[Cart::xyy][Cart::z][Cart::xxy]+amb1*R[Cart::xy][Cart::z][Cart::xxy];
R[Cart::yz][Cart::yz][Cart::xxy]=R[Cart::yyz][Cart::z][Cart::xxy]+amb1*R[Cart::yz][Cart::z][Cart::xxy];
R[Cart::xx][Cart::yz][Cart::xxy]=R[Cart::xxy][Cart::z][Cart::xxy]+amb1*R[Cart::xx][Cart::z][Cart::xxy];
R[Cart::xz][Cart::yz][Cart::xxy]=R[Cart::xyz][Cart::z][Cart::xxy]+amb1*R[Cart::xz][Cart::z][Cart::xxy];
R[Cart::zz][Cart::yz][Cart::xxy]=R[Cart::yzz][Cart::z][Cart::xxy]+amb1*R[Cart::zz][Cart::z][Cart::xxy];
R[Cart::yy][Cart::xx][Cart::xxy]=R[Cart::xyy][Cart::x][Cart::xxy]+amb0*R[Cart::yy][Cart::x][Cart::xxy];
R[Cart::xy][Cart::xx][Cart::xxy]=R[Cart::xxy][Cart::x][Cart::xxy]+amb0*R[Cart::xy][Cart::x][Cart::xxy];
R[Cart::yz][Cart::xx][Cart::xxy]=R[Cart::xyz][Cart::x][Cart::xxy]+amb0*R[Cart::yz][Cart::x][Cart::xxy];
R[Cart::xx][Cart::xx][Cart::xxy]=R[Cart::xxx][Cart::x][Cart::xxy]+amb0*R[Cart::xx][Cart::x][Cart::xxy];
R[Cart::xz][Cart::xx][Cart::xxy]=R[Cart::xxz][Cart::x][Cart::xxy]+amb0*R[Cart::xz][Cart::x][Cart::xxy];
R[Cart::zz][Cart::xx][Cart::xxy]=R[Cart::xzz][Cart::x][Cart::xxy]+amb0*R[Cart::zz][Cart::x][Cart::xxy];
R[Cart::yy][Cart::xz][Cart::xxy]=R[Cart::xyy][Cart::z][Cart::xxy]+amb0*R[Cart::yy][Cart::z][Cart::xxy];
R[Cart::xy][Cart::xz][Cart::xxy]=R[Cart::xxy][Cart::z][Cart::xxy]+amb0*R[Cart::xy][Cart::z][Cart::xxy];
R[Cart::yz][Cart::xz][Cart::xxy]=R[Cart::xyz][Cart::z][Cart::xxy]+amb0*R[Cart::yz][Cart::z][Cart::xxy];
R[Cart::xx][Cart::xz][Cart::xxy]=R[Cart::xxx][Cart::z][Cart::xxy]+amb0*R[Cart::xx][Cart::z][Cart::xxy];
R[Cart::xz][Cart::xz][Cart::xxy]=R[Cart::xxz][Cart::z][Cart::xxy]+amb0*R[Cart::xz][Cart::z][Cart::xxy];
R[Cart::zz][Cart::xz][Cart::xxy]=R[Cart::xzz][Cart::z][Cart::xxy]+amb0*R[Cart::zz][Cart::z][Cart::xxy];
R[Cart::yy][Cart::zz][Cart::xxy]=R[Cart::yyz][Cart::z][Cart::xxy]+amb2*R[Cart::yy][Cart::z][Cart::xxy];
R[Cart::xy][Cart::zz][Cart::xxy]=R[Cart::xyz][Cart::z][Cart::xxy]+amb2*R[Cart::xy][Cart::z][Cart::xxy];
R[Cart::yz][Cart::zz][Cart::xxy]=R[Cart::yzz][Cart::z][Cart::xxy]+amb2*R[Cart::yz][Cart::z][Cart::xxy];
R[Cart::xx][Cart::zz][Cart::xxy]=R[Cart::xxz][Cart::z][Cart::xxy]+amb2*R[Cart::xx][Cart::z][Cart::xxy];
R[Cart::xz][Cart::zz][Cart::xxy]=R[Cart::xzz][Cart::z][Cart::xxy]+amb2*R[Cart::xz][Cart::z][Cart::xxy];
R[Cart::zz][Cart::zz][Cart::xxy]=R[Cart::zzz][Cart::z][Cart::xxy]+amb2*R[Cart::zz][Cart::z][Cart::xxy];
R[Cart::yy][Cart::yy][Cart::xyz]=R[Cart::yyy][Cart::y][Cart::xyz]+amb1*R[Cart::yy][Cart::y][Cart::xyz];
R[Cart::xy][Cart::yy][Cart::xyz]=R[Cart::xyy][Cart::y][Cart::xyz]+amb1*R[Cart::xy][Cart::y][Cart::xyz];
R[Cart::yz][Cart::yy][Cart::xyz]=R[Cart::yyz][Cart::y][Cart::xyz]+amb1*R[Cart::yz][Cart::y][Cart::xyz];
R[Cart::xx][Cart::yy][Cart::xyz]=R[Cart::xxy][Cart::y][Cart::xyz]+amb1*R[Cart::xx][Cart::y][Cart::xyz];
R[Cart::xz][Cart::yy][Cart::xyz]=R[Cart::xyz][Cart::y][Cart::xyz]+amb1*R[Cart::xz][Cart::y][Cart::xyz];
R[Cart::zz][Cart::yy][Cart::xyz]=R[Cart::yzz][Cart::y][Cart::xyz]+amb1*R[Cart::zz][Cart::y][Cart::xyz];
R[Cart::yy][Cart::xy][Cart::xyz]=R[Cart::xyy][Cart::y][Cart::xyz]+amb0*R[Cart::yy][Cart::y][Cart::xyz];
R[Cart::xy][Cart::xy][Cart::xyz]=R[Cart::xxy][Cart::y][Cart::xyz]+amb0*R[Cart::xy][Cart::y][Cart::xyz];
R[Cart::yz][Cart::xy][Cart::xyz]=R[Cart::xyz][Cart::y][Cart::xyz]+amb0*R[Cart::yz][Cart::y][Cart::xyz];
R[Cart::xx][Cart::xy][Cart::xyz]=R[Cart::xxx][Cart::y][Cart::xyz]+amb0*R[Cart::xx][Cart::y][Cart::xyz];
R[Cart::xz][Cart::xy][Cart::xyz]=R[Cart::xxz][Cart::y][Cart::xyz]+amb0*R[Cart::xz][Cart::y][Cart::xyz];
R[Cart::zz][Cart::xy][Cart::xyz]=R[Cart::xzz][Cart::y][Cart::xyz]+amb0*R[Cart::zz][Cart::y][Cart::xyz];
R[Cart::yy][Cart::yz][Cart::xyz]=R[Cart::yyy][Cart::z][Cart::xyz]+amb1*R[Cart::yy][Cart::z][Cart::xyz];
R[Cart::xy][Cart::yz][Cart::xyz]=R[Cart::xyy][Cart::z][Cart::xyz]+amb1*R[Cart::xy][Cart::z][Cart::xyz];
R[Cart::yz][Cart::yz][Cart::xyz]=R[Cart::yyz][Cart::z][Cart::xyz]+amb1*R[Cart::yz][Cart::z][Cart::xyz];
R[Cart::xx][Cart::yz][Cart::xyz]=R[Cart::xxy][Cart::z][Cart::xyz]+amb1*R[Cart::xx][Cart::z][Cart::xyz];
R[Cart::xz][Cart::yz][Cart::xyz]=R[Cart::xyz][Cart::z][Cart::xyz]+amb1*R[Cart::xz][Cart::z][Cart::xyz];
R[Cart::zz][Cart::yz][Cart::xyz]=R[Cart::yzz][Cart::z][Cart::xyz]+amb1*R[Cart::zz][Cart::z][Cart::xyz];
R[Cart::yy][Cart::xx][Cart::xyz]=R[Cart::xyy][Cart::x][Cart::xyz]+amb0*R[Cart::yy][Cart::x][Cart::xyz];
R[Cart::xy][Cart::xx][Cart::xyz]=R[Cart::xxy][Cart::x][Cart::xyz]+amb0*R[Cart::xy][Cart::x][Cart::xyz];
R[Cart::yz][Cart::xx][Cart::xyz]=R[Cart::xyz][Cart::x][Cart::xyz]+amb0*R[Cart::yz][Cart::x][Cart::xyz];
R[Cart::xx][Cart::xx][Cart::xyz]=R[Cart::xxx][Cart::x][Cart::xyz]+amb0*R[Cart::xx][Cart::x][Cart::xyz];
R[Cart::xz][Cart::xx][Cart::xyz]=R[Cart::xxz][Cart::x][Cart::xyz]+amb0*R[Cart::xz][Cart::x][Cart::xyz];
R[Cart::zz][Cart::xx][Cart::xyz]=R[Cart::xzz][Cart::x][Cart::xyz]+amb0*R[Cart::zz][Cart::x][Cart::xyz];
R[Cart::yy][Cart::xz][Cart::xyz]=R[Cart::xyy][Cart::z][Cart::xyz]+amb0*R[Cart::yy][Cart::z][Cart::xyz];
R[Cart::xy][Cart::xz][Cart::xyz]=R[Cart::xxy][Cart::z][Cart::xyz]+amb0*R[Cart::xy][Cart::z][Cart::xyz];
R[Cart::yz][Cart::xz][Cart::xyz]=R[Cart::xyz][Cart::z][Cart::xyz]+amb0*R[Cart::yz][Cart::z][Cart::xyz];
R[Cart::xx][Cart::xz][Cart::xyz]=R[Cart::xxx][Cart::z][Cart::xyz]+amb0*R[Cart::xx][Cart::z][Cart::xyz];
R[Cart::xz][Cart::xz][Cart::xyz]=R[Cart::xxz][Cart::z][Cart::xyz]+amb0*R[Cart::xz][Cart::z][Cart::xyz];
R[Cart::zz][Cart::xz][Cart::xyz]=R[Cart::xzz][Cart::z][Cart::xyz]+amb0*R[Cart::zz][Cart::z][Cart::xyz];
R[Cart::yy][Cart::zz][Cart::xyz]=R[Cart::yyz][Cart::z][Cart::xyz]+amb2*R[Cart::yy][Cart::z][Cart::xyz];
R[Cart::xy][Cart::zz][Cart::xyz]=R[Cart::xyz][Cart::z][Cart::xyz]+amb2*R[Cart::xy][Cart::z][Cart::xyz];
R[Cart::yz][Cart::zz][Cart::xyz]=R[Cart::yzz][Cart::z][Cart::xyz]+amb2*R[Cart::yz][Cart::z][Cart::xyz];
R[Cart::xx][Cart::zz][Cart::xyz]=R[Cart::xxz][Cart::z][Cart::xyz]+amb2*R[Cart::xx][Cart::z][Cart::xyz];
R[Cart::xz][Cart::zz][Cart::xyz]=R[Cart::xzz][Cart::z][Cart::xyz]+amb2*R[Cart::xz][Cart::z][Cart::xyz];
R[Cart::zz][Cart::zz][Cart::xyz]=R[Cart::zzz][Cart::z][Cart::xyz]+amb2*R[Cart::zz][Cart::z][Cart::xyz];
R[Cart::yy][Cart::yy][Cart::yzz]=R[Cart::yyy][Cart::y][Cart::yzz]+amb1*R[Cart::yy][Cart::y][Cart::yzz];
R[Cart::xy][Cart::yy][Cart::yzz]=R[Cart::xyy][Cart::y][Cart::yzz]+amb1*R[Cart::xy][Cart::y][Cart::yzz];
R[Cart::yz][Cart::yy][Cart::yzz]=R[Cart::yyz][Cart::y][Cart::yzz]+amb1*R[Cart::yz][Cart::y][Cart::yzz];
R[Cart::xx][Cart::yy][Cart::yzz]=R[Cart::xxy][Cart::y][Cart::yzz]+amb1*R[Cart::xx][Cart::y][Cart::yzz];
R[Cart::xz][Cart::yy][Cart::yzz]=R[Cart::xyz][Cart::y][Cart::yzz]+amb1*R[Cart::xz][Cart::y][Cart::yzz];
R[Cart::zz][Cart::yy][Cart::yzz]=R[Cart::yzz][Cart::y][Cart::yzz]+amb1*R[Cart::zz][Cart::y][Cart::yzz];
R[Cart::yy][Cart::xy][Cart::yzz]=R[Cart::xyy][Cart::y][Cart::yzz]+amb0*R[Cart::yy][Cart::y][Cart::yzz];
R[Cart::xy][Cart::xy][Cart::yzz]=R[Cart::xxy][Cart::y][Cart::yzz]+amb0*R[Cart::xy][Cart::y][Cart::yzz];
R[Cart::yz][Cart::xy][Cart::yzz]=R[Cart::xyz][Cart::y][Cart::yzz]+amb0*R[Cart::yz][Cart::y][Cart::yzz];
R[Cart::xx][Cart::xy][Cart::yzz]=R[Cart::xxx][Cart::y][Cart::yzz]+amb0*R[Cart::xx][Cart::y][Cart::yzz];
R[Cart::xz][Cart::xy][Cart::yzz]=R[Cart::xxz][Cart::y][Cart::yzz]+amb0*R[Cart::xz][Cart::y][Cart::yzz];
R[Cart::zz][Cart::xy][Cart::yzz]=R[Cart::xzz][Cart::y][Cart::yzz]+amb0*R[Cart::zz][Cart::y][Cart::yzz];
R[Cart::yy][Cart::yz][Cart::yzz]=R[Cart::yyy][Cart::z][Cart::yzz]+amb1*R[Cart::yy][Cart::z][Cart::yzz];
R[Cart::xy][Cart::yz][Cart::yzz]=R[Cart::xyy][Cart::z][Cart::yzz]+amb1*R[Cart::xy][Cart::z][Cart::yzz];
R[Cart::yz][Cart::yz][Cart::yzz]=R[Cart::yyz][Cart::z][Cart::yzz]+amb1*R[Cart::yz][Cart::z][Cart::yzz];
R[Cart::xx][Cart::yz][Cart::yzz]=R[Cart::xxy][Cart::z][Cart::yzz]+amb1*R[Cart::xx][Cart::z][Cart::yzz];
R[Cart::xz][Cart::yz][Cart::yzz]=R[Cart::xyz][Cart::z][Cart::yzz]+amb1*R[Cart::xz][Cart::z][Cart::yzz];
R[Cart::zz][Cart::yz][Cart::yzz]=R[Cart::yzz][Cart::z][Cart::yzz]+amb1*R[Cart::zz][Cart::z][Cart::yzz];
R[Cart::yy][Cart::xx][Cart::yzz]=R[Cart::xyy][Cart::x][Cart::yzz]+amb0*R[Cart::yy][Cart::x][Cart::yzz];
R[Cart::xy][Cart::xx][Cart::yzz]=R[Cart::xxy][Cart::x][Cart::yzz]+amb0*R[Cart::xy][Cart::x][Cart::yzz];
R[Cart::yz][Cart::xx][Cart::yzz]=R[Cart::xyz][Cart::x][Cart::yzz]+amb0*R[Cart::yz][Cart::x][Cart::yzz];
R[Cart::xx][Cart::xx][Cart::yzz]=R[Cart::xxx][Cart::x][Cart::yzz]+amb0*R[Cart::xx][Cart::x][Cart::yzz];
R[Cart::xz][Cart::xx][Cart::yzz]=R[Cart::xxz][Cart::x][Cart::yzz]+amb0*R[Cart::xz][Cart::x][Cart::yzz];
R[Cart::zz][Cart::xx][Cart::yzz]=R[Cart::xzz][Cart::x][Cart::yzz]+amb0*R[Cart::zz][Cart::x][Cart::yzz];
R[Cart::yy][Cart::xz][Cart::yzz]=R[Cart::xyy][Cart::z][Cart::yzz]+amb0*R[Cart::yy][Cart::z][Cart::yzz];
R[Cart::xy][Cart::xz][Cart::yzz]=R[Cart::xxy][Cart::z][Cart::yzz]+amb0*R[Cart::xy][Cart::z][Cart::yzz];
R[Cart::yz][Cart::xz][Cart::yzz]=R[Cart::xyz][Cart::z][Cart::yzz]+amb0*R[Cart::yz][Cart::z][Cart::yzz];
R[Cart::xx][Cart::xz][Cart::yzz]=R[Cart::xxx][Cart::z][Cart::yzz]+amb0*R[Cart::xx][Cart::z][Cart::yzz];
R[Cart::xz][Cart::xz][Cart::yzz]=R[Cart::xxz][Cart::z][Cart::yzz]+amb0*R[Cart::xz][Cart::z][Cart::yzz];
R[Cart::zz][Cart::xz][Cart::yzz]=R[Cart::xzz][Cart::z][Cart::yzz]+amb0*R[Cart::zz][Cart::z][Cart::yzz];
R[Cart::yy][Cart::zz][Cart::yzz]=R[Cart::yyz][Cart::z][Cart::yzz]+amb2*R[Cart::yy][Cart::z][Cart::yzz];
R[Cart::xy][Cart::zz][Cart::yzz]=R[Cart::xyz][Cart::z][Cart::yzz]+amb2*R[Cart::xy][Cart::z][Cart::yzz];
R[Cart::yz][Cart::zz][Cart::yzz]=R[Cart::yzz][Cart::z][Cart::yzz]+amb2*R[Cart::yz][Cart::z][Cart::yzz];
R[Cart::xx][Cart::zz][Cart::yzz]=R[Cart::xxz][Cart::z][Cart::yzz]+amb2*R[Cart::xx][Cart::z][Cart::yzz];
R[Cart::xz][Cart::zz][Cart::yzz]=R[Cart::xzz][Cart::z][Cart::yzz]+amb2*R[Cart::xz][Cart::z][Cart::yzz];
R[Cart::zz][Cart::zz][Cart::yzz]=R[Cart::zzz][Cart::z][Cart::yzz]+amb2*R[Cart::zz][Cart::z][Cart::yzz];
R[Cart::yy][Cart::yy][Cart::xxx]=R[Cart::yyy][Cart::y][Cart::xxx]+amb1*R[Cart::yy][Cart::y][Cart::xxx];
R[Cart::xy][Cart::yy][Cart::xxx]=R[Cart::xyy][Cart::y][Cart::xxx]+amb1*R[Cart::xy][Cart::y][Cart::xxx];
R[Cart::yz][Cart::yy][Cart::xxx]=R[Cart::yyz][Cart::y][Cart::xxx]+amb1*R[Cart::yz][Cart::y][Cart::xxx];
R[Cart::xx][Cart::yy][Cart::xxx]=R[Cart::xxy][Cart::y][Cart::xxx]+amb1*R[Cart::xx][Cart::y][Cart::xxx];
R[Cart::xz][Cart::yy][Cart::xxx]=R[Cart::xyz][Cart::y][Cart::xxx]+amb1*R[Cart::xz][Cart::y][Cart::xxx];
R[Cart::zz][Cart::yy][Cart::xxx]=R[Cart::yzz][Cart::y][Cart::xxx]+amb1*R[Cart::zz][Cart::y][Cart::xxx];
R[Cart::yy][Cart::xy][Cart::xxx]=R[Cart::xyy][Cart::y][Cart::xxx]+amb0*R[Cart::yy][Cart::y][Cart::xxx];
R[Cart::xy][Cart::xy][Cart::xxx]=R[Cart::xxy][Cart::y][Cart::xxx]+amb0*R[Cart::xy][Cart::y][Cart::xxx];
R[Cart::yz][Cart::xy][Cart::xxx]=R[Cart::xyz][Cart::y][Cart::xxx]+amb0*R[Cart::yz][Cart::y][Cart::xxx];
R[Cart::xx][Cart::xy][Cart::xxx]=R[Cart::xxx][Cart::y][Cart::xxx]+amb0*R[Cart::xx][Cart::y][Cart::xxx];
R[Cart::xz][Cart::xy][Cart::xxx]=R[Cart::xxz][Cart::y][Cart::xxx]+amb0*R[Cart::xz][Cart::y][Cart::xxx];
R[Cart::zz][Cart::xy][Cart::xxx]=R[Cart::xzz][Cart::y][Cart::xxx]+amb0*R[Cart::zz][Cart::y][Cart::xxx];
R[Cart::yy][Cart::yz][Cart::xxx]=R[Cart::yyy][Cart::z][Cart::xxx]+amb1*R[Cart::yy][Cart::z][Cart::xxx];
R[Cart::xy][Cart::yz][Cart::xxx]=R[Cart::xyy][Cart::z][Cart::xxx]+amb1*R[Cart::xy][Cart::z][Cart::xxx];
R[Cart::yz][Cart::yz][Cart::xxx]=R[Cart::yyz][Cart::z][Cart::xxx]+amb1*R[Cart::yz][Cart::z][Cart::xxx];
R[Cart::xx][Cart::yz][Cart::xxx]=R[Cart::xxy][Cart::z][Cart::xxx]+amb1*R[Cart::xx][Cart::z][Cart::xxx];
R[Cart::xz][Cart::yz][Cart::xxx]=R[Cart::xyz][Cart::z][Cart::xxx]+amb1*R[Cart::xz][Cart::z][Cart::xxx];
R[Cart::zz][Cart::yz][Cart::xxx]=R[Cart::yzz][Cart::z][Cart::xxx]+amb1*R[Cart::zz][Cart::z][Cart::xxx];
R[Cart::yy][Cart::xx][Cart::xxx]=R[Cart::xyy][Cart::x][Cart::xxx]+amb0*R[Cart::yy][Cart::x][Cart::xxx];
R[Cart::xy][Cart::xx][Cart::xxx]=R[Cart::xxy][Cart::x][Cart::xxx]+amb0*R[Cart::xy][Cart::x][Cart::xxx];
R[Cart::yz][Cart::xx][Cart::xxx]=R[Cart::xyz][Cart::x][Cart::xxx]+amb0*R[Cart::yz][Cart::x][Cart::xxx];
R[Cart::xx][Cart::xx][Cart::xxx]=R[Cart::xxx][Cart::x][Cart::xxx]+amb0*R[Cart::xx][Cart::x][Cart::xxx];
R[Cart::xz][Cart::xx][Cart::xxx]=R[Cart::xxz][Cart::x][Cart::xxx]+amb0*R[Cart::xz][Cart::x][Cart::xxx];
R[Cart::zz][Cart::xx][Cart::xxx]=R[Cart::xzz][Cart::x][Cart::xxx]+amb0*R[Cart::zz][Cart::x][Cart::xxx];
R[Cart::yy][Cart::xz][Cart::xxx]=R[Cart::xyy][Cart::z][Cart::xxx]+amb0*R[Cart::yy][Cart::z][Cart::xxx];
R[Cart::xy][Cart::xz][Cart::xxx]=R[Cart::xxy][Cart::z][Cart::xxx]+amb0*R[Cart::xy][Cart::z][Cart::xxx];
R[Cart::yz][Cart::xz][Cart::xxx]=R[Cart::xyz][Cart::z][Cart::xxx]+amb0*R[Cart::yz][Cart::z][Cart::xxx];
R[Cart::xx][Cart::xz][Cart::xxx]=R[Cart::xxx][Cart::z][Cart::xxx]+amb0*R[Cart::xx][Cart::z][Cart::xxx];
R[Cart::xz][Cart::xz][Cart::xxx]=R[Cart::xxz][Cart::z][Cart::xxx]+amb0*R[Cart::xz][Cart::z][Cart::xxx];
R[Cart::zz][Cart::xz][Cart::xxx]=R[Cart::xzz][Cart::z][Cart::xxx]+amb0*R[Cart::zz][Cart::z][Cart::xxx];
R[Cart::yy][Cart::zz][Cart::xxx]=R[Cart::yyz][Cart::z][Cart::xxx]+amb2*R[Cart::yy][Cart::z][Cart::xxx];
R[Cart::xy][Cart::zz][Cart::xxx]=R[Cart::xyz][Cart::z][Cart::xxx]+amb2*R[Cart::xy][Cart::z][Cart::xxx];
R[Cart::yz][Cart::zz][Cart::xxx]=R[Cart::yzz][Cart::z][Cart::xxx]+amb2*R[Cart::yz][Cart::z][Cart::xxx];
R[Cart::xx][Cart::zz][Cart::xxx]=R[Cart::xxz][Cart::z][Cart::xxx]+amb2*R[Cart::xx][Cart::z][Cart::xxx];
R[Cart::xz][Cart::zz][Cart::xxx]=R[Cart::xzz][Cart::z][Cart::xxx]+amb2*R[Cart::xz][Cart::z][Cart::xxx];
R[Cart::zz][Cart::zz][Cart::xxx]=R[Cart::zzz][Cart::z][Cart::xxx]+amb2*R[Cart::zz][Cart::z][Cart::xxx];
R[Cart::yy][Cart::yy][Cart::xxz]=R[Cart::yyy][Cart::y][Cart::xxz]+amb1*R[Cart::yy][Cart::y][Cart::xxz];
R[Cart::xy][Cart::yy][Cart::xxz]=R[Cart::xyy][Cart::y][Cart::xxz]+amb1*R[Cart::xy][Cart::y][Cart::xxz];
R[Cart::yz][Cart::yy][Cart::xxz]=R[Cart::yyz][Cart::y][Cart::xxz]+amb1*R[Cart::yz][Cart::y][Cart::xxz];
R[Cart::xx][Cart::yy][Cart::xxz]=R[Cart::xxy][Cart::y][Cart::xxz]+amb1*R[Cart::xx][Cart::y][Cart::xxz];
R[Cart::xz][Cart::yy][Cart::xxz]=R[Cart::xyz][Cart::y][Cart::xxz]+amb1*R[Cart::xz][Cart::y][Cart::xxz];
R[Cart::zz][Cart::yy][Cart::xxz]=R[Cart::yzz][Cart::y][Cart::xxz]+amb1*R[Cart::zz][Cart::y][Cart::xxz];
R[Cart::yy][Cart::xy][Cart::xxz]=R[Cart::xyy][Cart::y][Cart::xxz]+amb0*R[Cart::yy][Cart::y][Cart::xxz];
R[Cart::xy][Cart::xy][Cart::xxz]=R[Cart::xxy][Cart::y][Cart::xxz]+amb0*R[Cart::xy][Cart::y][Cart::xxz];
R[Cart::yz][Cart::xy][Cart::xxz]=R[Cart::xyz][Cart::y][Cart::xxz]+amb0*R[Cart::yz][Cart::y][Cart::xxz];
R[Cart::xx][Cart::xy][Cart::xxz]=R[Cart::xxx][Cart::y][Cart::xxz]+amb0*R[Cart::xx][Cart::y][Cart::xxz];
R[Cart::xz][Cart::xy][Cart::xxz]=R[Cart::xxz][Cart::y][Cart::xxz]+amb0*R[Cart::xz][Cart::y][Cart::xxz];
R[Cart::zz][Cart::xy][Cart::xxz]=R[Cart::xzz][Cart::y][Cart::xxz]+amb0*R[Cart::zz][Cart::y][Cart::xxz];
R[Cart::yy][Cart::yz][Cart::xxz]=R[Cart::yyy][Cart::z][Cart::xxz]+amb1*R[Cart::yy][Cart::z][Cart::xxz];
R[Cart::xy][Cart::yz][Cart::xxz]=R[Cart::xyy][Cart::z][Cart::xxz]+amb1*R[Cart::xy][Cart::z][Cart::xxz];
R[Cart::yz][Cart::yz][Cart::xxz]=R[Cart::yyz][Cart::z][Cart::xxz]+amb1*R[Cart::yz][Cart::z][Cart::xxz];
R[Cart::xx][Cart::yz][Cart::xxz]=R[Cart::xxy][Cart::z][Cart::xxz]+amb1*R[Cart::xx][Cart::z][Cart::xxz];
R[Cart::xz][Cart::yz][Cart::xxz]=R[Cart::xyz][Cart::z][Cart::xxz]+amb1*R[Cart::xz][Cart::z][Cart::xxz];
R[Cart::zz][Cart::yz][Cart::xxz]=R[Cart::yzz][Cart::z][Cart::xxz]+amb1*R[Cart::zz][Cart::z][Cart::xxz];
R[Cart::yy][Cart::xx][Cart::xxz]=R[Cart::xyy][Cart::x][Cart::xxz]+amb0*R[Cart::yy][Cart::x][Cart::xxz];
R[Cart::xy][Cart::xx][Cart::xxz]=R[Cart::xxy][Cart::x][Cart::xxz]+amb0*R[Cart::xy][Cart::x][Cart::xxz];
R[Cart::yz][Cart::xx][Cart::xxz]=R[Cart::xyz][Cart::x][Cart::xxz]+amb0*R[Cart::yz][Cart::x][Cart::xxz];
R[Cart::xx][Cart::xx][Cart::xxz]=R[Cart::xxx][Cart::x][Cart::xxz]+amb0*R[Cart::xx][Cart::x][Cart::xxz];
R[Cart::xz][Cart::xx][Cart::xxz]=R[Cart::xxz][Cart::x][Cart::xxz]+amb0*R[Cart::xz][Cart::x][Cart::xxz];
R[Cart::zz][Cart::xx][Cart::xxz]=R[Cart::xzz][Cart::x][Cart::xxz]+amb0*R[Cart::zz][Cart::x][Cart::xxz];
R[Cart::yy][Cart::xz][Cart::xxz]=R[Cart::xyy][Cart::z][Cart::xxz]+amb0*R[Cart::yy][Cart::z][Cart::xxz];
R[Cart::xy][Cart::xz][Cart::xxz]=R[Cart::xxy][Cart::z][Cart::xxz]+amb0*R[Cart::xy][Cart::z][Cart::xxz];
R[Cart::yz][Cart::xz][Cart::xxz]=R[Cart::xyz][Cart::z][Cart::xxz]+amb0*R[Cart::yz][Cart::z][Cart::xxz];
R[Cart::xx][Cart::xz][Cart::xxz]=R[Cart::xxx][Cart::z][Cart::xxz]+amb0*R[Cart::xx][Cart::z][Cart::xxz];
R[Cart::xz][Cart::xz][Cart::xxz]=R[Cart::xxz][Cart::z][Cart::xxz]+amb0*R[Cart::xz][Cart::z][Cart::xxz];
R[Cart::zz][Cart::xz][Cart::xxz]=R[Cart::xzz][Cart::z][Cart::xxz]+amb0*R[Cart::zz][Cart::z][Cart::xxz];
R[Cart::yy][Cart::zz][Cart::xxz]=R[Cart::yyz][Cart::z][Cart::xxz]+amb2*R[Cart::yy][Cart::z][Cart::xxz];
R[Cart::xy][Cart::zz][Cart::xxz]=R[Cart::xyz][Cart::z][Cart::xxz]+amb2*R[Cart::xy][Cart::z][Cart::xxz];
R[Cart::yz][Cart::zz][Cart::xxz]=R[Cart::yzz][Cart::z][Cart::xxz]+amb2*R[Cart::yz][Cart::z][Cart::xxz];
R[Cart::xx][Cart::zz][Cart::xxz]=R[Cart::xxz][Cart::z][Cart::xxz]+amb2*R[Cart::xx][Cart::z][Cart::xxz];
R[Cart::xz][Cart::zz][Cart::xxz]=R[Cart::xzz][Cart::z][Cart::xxz]+amb2*R[Cart::xz][Cart::z][Cart::xxz];
R[Cart::zz][Cart::zz][Cart::xxz]=R[Cart::zzz][Cart::z][Cart::xxz]+amb2*R[Cart::zz][Cart::z][Cart::xxz];
R[Cart::yy][Cart::yy][Cart::xzz]=R[Cart::yyy][Cart::y][Cart::xzz]+amb1*R[Cart::yy][Cart::y][Cart::xzz];
R[Cart::xy][Cart::yy][Cart::xzz]=R[Cart::xyy][Cart::y][Cart::xzz]+amb1*R[Cart::xy][Cart::y][Cart::xzz];
R[Cart::yz][Cart::yy][Cart::xzz]=R[Cart::yyz][Cart::y][Cart::xzz]+amb1*R[Cart::yz][Cart::y][Cart::xzz];
R[Cart::xx][Cart::yy][Cart::xzz]=R[Cart::xxy][Cart::y][Cart::xzz]+amb1*R[Cart::xx][Cart::y][Cart::xzz];
R[Cart::xz][Cart::yy][Cart::xzz]=R[Cart::xyz][Cart::y][Cart::xzz]+amb1*R[Cart::xz][Cart::y][Cart::xzz];
R[Cart::zz][Cart::yy][Cart::xzz]=R[Cart::yzz][Cart::y][Cart::xzz]+amb1*R[Cart::zz][Cart::y][Cart::xzz];
R[Cart::yy][Cart::xy][Cart::xzz]=R[Cart::xyy][Cart::y][Cart::xzz]+amb0*R[Cart::yy][Cart::y][Cart::xzz];
R[Cart::xy][Cart::xy][Cart::xzz]=R[Cart::xxy][Cart::y][Cart::xzz]+amb0*R[Cart::xy][Cart::y][Cart::xzz];
R[Cart::yz][Cart::xy][Cart::xzz]=R[Cart::xyz][Cart::y][Cart::xzz]+amb0*R[Cart::yz][Cart::y][Cart::xzz];
R[Cart::xx][Cart::xy][Cart::xzz]=R[Cart::xxx][Cart::y][Cart::xzz]+amb0*R[Cart::xx][Cart::y][Cart::xzz];
R[Cart::xz][Cart::xy][Cart::xzz]=R[Cart::xxz][Cart::y][Cart::xzz]+amb0*R[Cart::xz][Cart::y][Cart::xzz];
R[Cart::zz][Cart::xy][Cart::xzz]=R[Cart::xzz][Cart::y][Cart::xzz]+amb0*R[Cart::zz][Cart::y][Cart::xzz];
R[Cart::yy][Cart::yz][Cart::xzz]=R[Cart::yyy][Cart::z][Cart::xzz]+amb1*R[Cart::yy][Cart::z][Cart::xzz];
R[Cart::xy][Cart::yz][Cart::xzz]=R[Cart::xyy][Cart::z][Cart::xzz]+amb1*R[Cart::xy][Cart::z][Cart::xzz];
R[Cart::yz][Cart::yz][Cart::xzz]=R[Cart::yyz][Cart::z][Cart::xzz]+amb1*R[Cart::yz][Cart::z][Cart::xzz];
R[Cart::xx][Cart::yz][Cart::xzz]=R[Cart::xxy][Cart::z][Cart::xzz]+amb1*R[Cart::xx][Cart::z][Cart::xzz];
R[Cart::xz][Cart::yz][Cart::xzz]=R[Cart::xyz][Cart::z][Cart::xzz]+amb1*R[Cart::xz][Cart::z][Cart::xzz];
R[Cart::zz][Cart::yz][Cart::xzz]=R[Cart::yzz][Cart::z][Cart::xzz]+amb1*R[Cart::zz][Cart::z][Cart::xzz];
R[Cart::yy][Cart::xx][Cart::xzz]=R[Cart::xyy][Cart::x][Cart::xzz]+amb0*R[Cart::yy][Cart::x][Cart::xzz];
R[Cart::xy][Cart::xx][Cart::xzz]=R[Cart::xxy][Cart::x][Cart::xzz]+amb0*R[Cart::xy][Cart::x][Cart::xzz];
R[Cart::yz][Cart::xx][Cart::xzz]=R[Cart::xyz][Cart::x][Cart::xzz]+amb0*R[Cart::yz][Cart::x][Cart::xzz];
R[Cart::xx][Cart::xx][Cart::xzz]=R[Cart::xxx][Cart::x][Cart::xzz]+amb0*R[Cart::xx][Cart::x][Cart::xzz];
R[Cart::xz][Cart::xx][Cart::xzz]=R[Cart::xxz][Cart::x][Cart::xzz]+amb0*R[Cart::xz][Cart::x][Cart::xzz];
R[Cart::zz][Cart::xx][Cart::xzz]=R[Cart::xzz][Cart::x][Cart::xzz]+amb0*R[Cart::zz][Cart::x][Cart::xzz];
R[Cart::yy][Cart::xz][Cart::xzz]=R[Cart::xyy][Cart::z][Cart::xzz]+amb0*R[Cart::yy][Cart::z][Cart::xzz];
R[Cart::xy][Cart::xz][Cart::xzz]=R[Cart::xxy][Cart::z][Cart::xzz]+amb0*R[Cart::xy][Cart::z][Cart::xzz];
R[Cart::yz][Cart::xz][Cart::xzz]=R[Cart::xyz][Cart::z][Cart::xzz]+amb0*R[Cart::yz][Cart::z][Cart::xzz];
R[Cart::xx][Cart::xz][Cart::xzz]=R[Cart::xxx][Cart::z][Cart::xzz]+amb0*R[Cart::xx][Cart::z][Cart::xzz];
R[Cart::xz][Cart::xz][Cart::xzz]=R[Cart::xxz][Cart::z][Cart::xzz]+amb0*R[Cart::xz][Cart::z][Cart::xzz];
R[Cart::zz][Cart::xz][Cart::xzz]=R[Cart::xzz][Cart::z][Cart::xzz]+amb0*R[Cart::zz][Cart::z][Cart::xzz];
R[Cart::yy][Cart::zz][Cart::xzz]=R[Cart::yyz][Cart::z][Cart::xzz]+amb2*R[Cart::yy][Cart::z][Cart::xzz];
R[Cart::xy][Cart::zz][Cart::xzz]=R[Cart::xyz][Cart::z][Cart::xzz]+amb2*R[Cart::xy][Cart::z][Cart::xzz];
R[Cart::yz][Cart::zz][Cart::xzz]=R[Cart::yzz][Cart::z][Cart::xzz]+amb2*R[Cart::yz][Cart::z][Cart::xzz];
R[Cart::xx][Cart::zz][Cart::xzz]=R[Cart::xxz][Cart::z][Cart::xzz]+amb2*R[Cart::xx][Cart::z][Cart::xzz];
R[Cart::xz][Cart::zz][Cart::xzz]=R[Cart::xzz][Cart::z][Cart::xzz]+amb2*R[Cart::xz][Cart::z][Cart::xzz];
R[Cart::zz][Cart::zz][Cart::xzz]=R[Cart::zzz][Cart::z][Cart::xzz]+amb2*R[Cart::zz][Cart::z][Cart::xzz];
R[Cart::yy][Cart::yy][Cart::zzz]=R[Cart::yyy][Cart::y][Cart::zzz]+amb1*R[Cart::yy][Cart::y][Cart::zzz];
R[Cart::xy][Cart::yy][Cart::zzz]=R[Cart::xyy][Cart::y][Cart::zzz]+amb1*R[Cart::xy][Cart::y][Cart::zzz];
R[Cart::yz][Cart::yy][Cart::zzz]=R[Cart::yyz][Cart::y][Cart::zzz]+amb1*R[Cart::yz][Cart::y][Cart::zzz];
R[Cart::xx][Cart::yy][Cart::zzz]=R[Cart::xxy][Cart::y][Cart::zzz]+amb1*R[Cart::xx][Cart::y][Cart::zzz];
R[Cart::xz][Cart::yy][Cart::zzz]=R[Cart::xyz][Cart::y][Cart::zzz]+amb1*R[Cart::xz][Cart::y][Cart::zzz];
R[Cart::zz][Cart::yy][Cart::zzz]=R[Cart::yzz][Cart::y][Cart::zzz]+amb1*R[Cart::zz][Cart::y][Cart::zzz];
R[Cart::yy][Cart::xy][Cart::zzz]=R[Cart::xyy][Cart::y][Cart::zzz]+amb0*R[Cart::yy][Cart::y][Cart::zzz];
R[Cart::xy][Cart::xy][Cart::zzz]=R[Cart::xxy][Cart::y][Cart::zzz]+amb0*R[Cart::xy][Cart::y][Cart::zzz];
R[Cart::yz][Cart::xy][Cart::zzz]=R[Cart::xyz][Cart::y][Cart::zzz]+amb0*R[Cart::yz][Cart::y][Cart::zzz];
R[Cart::xx][Cart::xy][Cart::zzz]=R[Cart::xxx][Cart::y][Cart::zzz]+amb0*R[Cart::xx][Cart::y][Cart::zzz];
R[Cart::xz][Cart::xy][Cart::zzz]=R[Cart::xxz][Cart::y][Cart::zzz]+amb0*R[Cart::xz][Cart::y][Cart::zzz];
R[Cart::zz][Cart::xy][Cart::zzz]=R[Cart::xzz][Cart::y][Cart::zzz]+amb0*R[Cart::zz][Cart::y][Cart::zzz];
R[Cart::yy][Cart::yz][Cart::zzz]=R[Cart::yyy][Cart::z][Cart::zzz]+amb1*R[Cart::yy][Cart::z][Cart::zzz];
R[Cart::xy][Cart::yz][Cart::zzz]=R[Cart::xyy][Cart::z][Cart::zzz]+amb1*R[Cart::xy][Cart::z][Cart::zzz];
R[Cart::yz][Cart::yz][Cart::zzz]=R[Cart::yyz][Cart::z][Cart::zzz]+amb1*R[Cart::yz][Cart::z][Cart::zzz];
R[Cart::xx][Cart::yz][Cart::zzz]=R[Cart::xxy][Cart::z][Cart::zzz]+amb1*R[Cart::xx][Cart::z][Cart::zzz];
R[Cart::xz][Cart::yz][Cart::zzz]=R[Cart::xyz][Cart::z][Cart::zzz]+amb1*R[Cart::xz][Cart::z][Cart::zzz];
R[Cart::zz][Cart::yz][Cart::zzz]=R[Cart::yzz][Cart::z][Cart::zzz]+amb1*R[Cart::zz][Cart::z][Cart::zzz];
R[Cart::yy][Cart::xx][Cart::zzz]=R[Cart::xyy][Cart::x][Cart::zzz]+amb0*R[Cart::yy][Cart::x][Cart::zzz];
R[Cart::xy][Cart::xx][Cart::zzz]=R[Cart::xxy][Cart::x][Cart::zzz]+amb0*R[Cart::xy][Cart::x][Cart::zzz];
R[Cart::yz][Cart::xx][Cart::zzz]=R[Cart::xyz][Cart::x][Cart::zzz]+amb0*R[Cart::yz][Cart::x][Cart::zzz];
R[Cart::xx][Cart::xx][Cart::zzz]=R[Cart::xxx][Cart::x][Cart::zzz]+amb0*R[Cart::xx][Cart::x][Cart::zzz];
R[Cart::xz][Cart::xx][Cart::zzz]=R[Cart::xxz][Cart::x][Cart::zzz]+amb0*R[Cart::xz][Cart::x][Cart::zzz];
R[Cart::zz][Cart::xx][Cart::zzz]=R[Cart::xzz][Cart::x][Cart::zzz]+amb0*R[Cart::zz][Cart::x][Cart::zzz];
R[Cart::yy][Cart::xz][Cart::zzz]=R[Cart::xyy][Cart::z][Cart::zzz]+amb0*R[Cart::yy][Cart::z][Cart::zzz];
R[Cart::xy][Cart::xz][Cart::zzz]=R[Cart::xxy][Cart::z][Cart::zzz]+amb0*R[Cart::xy][Cart::z][Cart::zzz];
R[Cart::yz][Cart::xz][Cart::zzz]=R[Cart::xyz][Cart::z][Cart::zzz]+amb0*R[Cart::yz][Cart::z][Cart::zzz];
R[Cart::xx][Cart::xz][Cart::zzz]=R[Cart::xxx][Cart::z][Cart::zzz]+amb0*R[Cart::xx][Cart::z][Cart::zzz];
R[Cart::xz][Cart::xz][Cart::zzz]=R[Cart::xxz][Cart::z][Cart::zzz]+amb0*R[Cart::xz][Cart::z][Cart::zzz];
R[Cart::zz][Cart::xz][Cart::zzz]=R[Cart::xzz][Cart::z][Cart::zzz]+amb0*R[Cart::zz][Cart::z][Cart::zzz];
R[Cart::yy][Cart::zz][Cart::zzz]=R[Cart::yyz][Cart::z][Cart::zzz]+amb2*R[Cart::yy][Cart::z][Cart::zzz];
R[Cart::xy][Cart::zz][Cart::zzz]=R[Cart::xyz][Cart::z][Cart::zzz]+amb2*R[Cart::xy][Cart::z][Cart::zzz];
R[Cart::yz][Cart::zz][Cart::zzz]=R[Cart::yzz][Cart::z][Cart::zzz]+amb2*R[Cart::yz][Cart::z][Cart::zzz];
R[Cart::xx][Cart::zz][Cart::zzz]=R[Cart::xxz][Cart::z][Cart::zzz]+amb2*R[Cart::xx][Cart::z][Cart::zzz];
R[Cart::xz][Cart::zz][Cart::zzz]=R[Cart::xzz][Cart::z][Cart::zzz]+amb2*R[Cart::xz][Cart::z][Cart::zzz];
R[Cart::zz][Cart::zz][Cart::zzz]=R[Cart::zzz][Cart::z][Cart::zzz]+amb2*R[Cart::zz][Cart::z][Cart::zzz];
}
//------------------------------------------------------





   // data is now stored in unnormalized cartesian Gaussians in the multiarray
            // Now, weird-looking construction since multiarray is not accessible for ub::prod
            //              s  px  py  pz dxz dyz dxy d3z2-r2 dx2-y2  f1  f2  f3  f4  f5  f6  f7
            int istart[] = {0, 1, 2, 3, 5, 6, 4, 7, 7, 12, 10, 11, 11, 10, 19, 15};
            int istop[] = {0, 1, 2, 3, 5, 6, 4, 9, 8, 17, 16, 18, 13, 14, 19, 17};

            // ub::vector<ub::matrix<double> >& _subvector
            // which ones do we want to store
            int _offset_beta = _shell_beta->getOffset();
            int _offset_alpha = _shell_alpha->getOffset();
            int _offset_gamma = _shell_gamma->getOffset();

            // prepare transformation matrices
            int _ntrafo_beta = _shell_beta->getNumFunc() + _offset_beta;
            int _ntrafo_alpha = _shell_alpha->getNumFunc() + _offset_alpha;
            int _ntrafo_gamma = _shell_gamma->getNumFunc() + _offset_gamma;

            ub::matrix<double> _trafo_beta = ub::zero_matrix<double>(_ntrafo_beta, _nbeta);
            ub::matrix<double> _trafo_alpha = ub::zero_matrix<double>(_ntrafo_alpha, _nalpha);
            ub::matrix<double> _trafo_gamma = ub::zero_matrix<double>(_ntrafo_gamma, _ngamma);

            
            std::vector<double> _contractions_alpha = (*italpha)->contraction;
            std::vector<double> _contractions_gamma = (*itgamma)->contraction;
            std::vector<double> _contractions_beta    = (*itbeta)->contraction;
            
            // get transformation matrices
            this->getTrafo(_trafo_beta, _lmax_beta, _decay_beta, _contractions_beta);
            this->getTrafo(_trafo_alpha, _lmax_alpha, _decay_alpha, _contractions_alpha);
            this->getTrafo(_trafo_gamma, _lmax_gamma, _decay_gamma, _contractions_gamma);

            // transform from unnormalized cartesians to normalized sphericals
            // container with indices starting at zero
            ma_type R_sph;
            R_sph.resize(extents[ _ntrafo_alpha ][ _ntrafo_beta ][ _ntrafo_gamma ]);

            for (int _i_beta = 0; _i_beta < _ntrafo_beta; _i_beta++) {
                for (int _i_alpha = 0; _i_alpha < _ntrafo_alpha; _i_alpha++) {
                    for (int _i_gamma = 0; _i_gamma < _ntrafo_gamma; _i_gamma++) {

                        R_sph[ _i_alpha ][ _i_beta ][ _i_gamma ] = 0.0;

                        for (int _i_beta_t = istart[ _i_beta ]; _i_beta_t <= istop[ _i_beta ]; _i_beta_t++) {
                            for (int _i_alpha_t = istart[ _i_alpha ]; _i_alpha_t <= istop[ _i_alpha ]; _i_alpha_t++) {
                                for (int _i_gamma_t = istart[ _i_gamma ]; _i_gamma_t <= istop[ _i_gamma ]; _i_gamma_t++) {

                                    R_sph[ _i_alpha ][ _i_beta ][ _i_gamma ] += R[ _i_alpha_t ][ _i_beta_t][ _i_gamma_t]
                                            * _trafo_alpha(_i_alpha, _i_alpha_t) * _trafo_beta(_i_beta, _i_beta_t) * _trafo_gamma(_i_gamma, _i_gamma_t);


                                }
                            }
                        }
                    }
                }
            }

            
            if(alphabetaswitch==true){
            cout << "switched back" << endl;    
            // only store the parts, we need
           
               for (int _i_gamma = 0; _i_gamma < _shell_gamma->getNumFunc(); _i_gamma++) {
                for (int _i_alpha = 0; _i_alpha < _shell_alpha->getNumFunc(); _i_alpha++) {
                    for (int _i_beta = 0; _i_beta < _shell_beta->getNumFunc(); _i_beta++) {


                        int _i_index = _shell_alpha->getNumFunc() * _i_gamma + _i_alpha;

                        _subvector(_i_beta, _i_index) += R_sph[ _offset_alpha + _i_alpha ][ _offset_beta + _i_beta ][ _offset_gamma + _i_gamma ];

                    }
                }
            } 
            }
            
            else{

            
               for (int _i_gamma = 0; _i_gamma < _shell_gamma->getNumFunc(); _i_gamma++) {
                for (int _i_alpha = 0; _i_alpha < _shell_alpha->getNumFunc(); _i_alpha++) {
                    for (int _i_beta = 0; _i_beta < _shell_beta->getNumFunc(); _i_beta++) {


                        int _i_index = _shell_beta->getNumFunc() * _i_gamma + _i_beta;

                        _subvector(_i_alpha, _i_index) += R_sph[ _offset_alpha + _i_alpha ][ _offset_beta + _i_beta ][ _offset_gamma + _i_gamma ];

                    }
                }
            } 
            }
            
            
            
            
            


                }
            }
        }

                        
 
                    cout << "ende" << endl;
    
       return _does_contribute;     
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

            for (unsigned m = 1; m < _FmT.size(); m++ ){
                _FmT[m] = (2*m-1) * _FmT[m-1]/(2.0*_T) - exp(-_T)/(2.0*_T) ;
            }
        }

        if ( _T < 1e-10 ){
           for ( unsigned m=0; m < _FmT.size(); m++){
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
