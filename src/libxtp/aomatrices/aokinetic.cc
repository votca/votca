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

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>

#include <vector>



namespace votca { namespace xtp {


    
    void AOKinetic::FillBlock( Eigen::Block<Eigen::MatrixXd> & _matrix,const AOShell* _shell_row,const AOShell* _shell_col) {
       
       
        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();
        
        
        if (_lmax_col >4 || _lmax_row >4){
            std::cerr << "Orbitals higher than g are not yet implemented. This should not have happened!" << std::flush;
             exit(1);
        }

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
        
             // get shell positions
        const tools::vec& _pos_row = _shell_row->getPos();
        const tools::vec& _pos_col = _shell_col->getPos();
        const tools::vec  _diff    = _pos_row - _pos_col;
       
          
        double _distsq = (_diff*_diff);   
        
     
        
        int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};
        
        int nx[] = { 0,
              1, 0, 0,
              2, 1, 1, 0, 0, 0,
              3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
              4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0 };

 int ny[] = { 0,
              0, 1, 0,
              0, 1, 0, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
              0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0 };

 int nz[] = { 0,
              0, 0, 1,
              0, 0, 1, 0, 1, 2,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
              0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4 };


 int i_less_x[] = {  0,
                     0,  0,  0,
                     1,  2,  3,  0,  0,  0,
                     4,  5,  6,  7,  8,  9,  0,  0,  0,  0,
                    10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  0,  0,  0,  0,  0 };

 int i_less_y[] = {  0,
                     0,  0,  0,
                     0,  1,  0,  2,  3,  0,
                     0,  4,  0,  5,  6,  0,  7,  8,  9,  0,
                     0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,  0 };

 int i_less_z[] = {  0,
                     0,  0,  0,
                     0,  0,  1,  0,  2,  3,
                     0,  0,  4,  0,  5,  6,  0,  7,  8,  9,
                     0,  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19 };

       // cout << "row shell is " << _shell_row->getSize() << " -fold contracted!" << endl;
        //cout << "col shell is " << _shell_col->getSize() << " -fold contracted!" << endl;
        
       
        // iterate over Gaussians in this _shell_row
        for ( AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
            const double _decay_row = itr->getDecay();
            
            for ( AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
           
            
           

            // get decay constants 
            const double _decay_col = itc->getDecay();
            
            // some helpers
            const double _fak  = 0.5/(_decay_row + _decay_col);
            const double rzeta = 2.0 * _fak;
            
            const double xi=rzeta* _decay_row * _decay_col;
            
            
            // check if distance between postions is big, then skip step   
            double _exparg = xi *_distsq;
	    if ( _exparg > 30.0 ) { continue; }
            
    
            
            const double PmA0 = rzeta*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
            const double PmA1 = rzeta*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
            const double PmA2 = rzeta*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

            const double PmB0 = rzeta*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
            const double PmB1 = rzeta*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
            const double PmB2 = rzeta*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
             const double xi2 = 2.*xi; //////////////
             const double _fak_a = rzeta * _decay_col; /////////////
             const double _fak_b = rzeta * _decay_row; ////////////
            
            // matrix for kinetic energies
            Eigen::MatrixXd kin = Eigen::MatrixXd::Zero(_nrows,_ncols);
            //matrix for unnormalized overlap integrals
            Eigen::MatrixXd ol = Eigen::MatrixXd::Zero(_nrows,_ncols);
        
            // s-s overlap integral
            ol(Cart::s,Cart::s) = pow(rzeta,1.5)*pow(4.0*_decay_row*_decay_col,0.75) * exp(-_exparg);
            // s-s- kinetic energy integral
            kin(Cart::s,Cart::s)= ol(Cart::s,Cart::s)*xi*(3-2*xi*_distsq);
            
//Integrals     p - s
if (_lmax_row > 0) {
  ol(Cart::x,0) = PmA0*ol(0,0);
  ol(Cart::y,0) = PmA1*ol(0,0);
  ol(Cart::z,0) = PmA2*ol(0,0);
}
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 1) {
  double term = _fak*ol(0,0);
  ol(Cart::xx,0) = PmA0*ol(Cart::x,0) + term;
  ol(Cart::xy,0) = PmA0*ol(Cart::y,0);
  ol(Cart::xz,0) = PmA0*ol(Cart::z,0);
  ol(Cart::yy,0) = PmA1*ol(Cart::y,0) + term;
  ol(Cart::yz,0) = PmA1*ol(Cart::z,0);
  ol(Cart::zz,0) = PmA2*ol(Cart::z,0) + term;
}
//------------------------------------------------------

//Integrals     f - s
if (_lmax_row > 2) {
  ol(Cart::xxx,0) = PmA0*ol(Cart::xx,0) + 2*_fak*ol(Cart::x,0);
  ol(Cart::xxy,0) = PmA1*ol(Cart::xx,0);
  ol(Cart::xxz,0) = PmA2*ol(Cart::xx,0);
  ol(Cart::xyy,0) = PmA0*ol(Cart::yy,0);
  ol(Cart::xyz,0) = PmA0*ol(Cart::yz,0);
  ol(Cart::xzz,0) = PmA0*ol(Cart::zz,0);
  ol(Cart::yyy,0) = PmA1*ol(Cart::yy,0) + 2*_fak*ol(Cart::y,0);
  ol(Cart::yyz,0) = PmA2*ol(Cart::yy,0);
  ol(Cart::yzz,0) = PmA1*ol(Cart::zz,0);
  ol(Cart::zzz,0) = PmA2*ol(Cart::zz,0) + 2*_fak*ol(Cart::z,0);
}
//------------------------------------------------------

//Integrals     g - s
if (_lmax_row > 3) {
  double term_xx = _fak*ol(Cart::xx,0);
  double term_yy = _fak*ol(Cart::yy,0);
  double term_zz = _fak*ol(Cart::zz,0);
  ol(Cart::xxxx,0) = PmA0*ol(Cart::xxx,0) + 3*term_xx;
  ol(Cart::xxxy,0) = PmA1*ol(Cart::xxx,0);
  ol(Cart::xxxz,0) = PmA2*ol(Cart::xxx,0);
  ol(Cart::xxyy,0) = PmA0*ol(Cart::xyy,0) + term_yy;
  ol(Cart::xxyz,0) = PmA1*ol(Cart::xxz,0);
  ol(Cart::xxzz,0) = PmA0*ol(Cart::xzz,0) + term_zz;
  ol(Cart::xyyy,0) = PmA0*ol(Cart::yyy,0);
  ol(Cart::xyyz,0) = PmA0*ol(Cart::yyz,0);
  ol(Cart::xyzz,0) = PmA0*ol(Cart::yzz,0);
  ol(Cart::xzzz,0) = PmA0*ol(Cart::zzz,0);
  ol(Cart::yyyy,0) = PmA1*ol(Cart::yyy,0) + 3*term_yy;
  ol(Cart::yyyz,0) = PmA2*ol(Cart::yyy,0);
  ol(Cart::yyzz,0) = PmA1*ol(Cart::yzz,0) + term_zz;
  ol(Cart::yzzz,0) = PmA1*ol(Cart::zzz,0);
  ol(Cart::zzzz,0) = PmA2*ol(Cart::zzz,0) + 3*term_zz;
}
//------------------------------------------------------



if (_lmax_col > 0) {

  //Integrals     s - p
  ol(0,Cart::x) = PmB0*ol(0,0);
  ol(0,Cart::y) = PmB1*ol(0,0);
  ol(0,Cart::z) = PmB2*ol(0,0);
  //------------------------------------------------------

  //Integrals     p - p     d - p     f - p     g - p
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    ol(_i,Cart::x) = PmB0*ol(_i,0) + nx[_i]*_fak*ol(i_less_x[_i],0);
    ol(_i,Cart::y) = PmB1*ol(_i,0) + ny[_i]*_fak*ol(i_less_y[_i],0);
    ol(_i,Cart::z) = PmB2*ol(_i,0) + nz[_i]*_fak*ol(i_less_z[_i],0);
  }
  //------------------------------------------------------

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

  //Integrals     s - d
  double term = _fak*ol(0,0);
  ol(0,Cart::xx) = PmB0*ol(0,Cart::x) + term;
  ol(0,Cart::xy) = PmB0*ol(0,Cart::y);
  ol(0,Cart::xz) = PmB0*ol(0,Cart::z);
  ol(0,Cart::yy) = PmB1*ol(0,Cart::y) + term;
  ol(0,Cart::yz) = PmB1*ol(0,Cart::z);
  ol(0,Cart::zz) = PmB2*ol(0,Cart::z) + term;
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    double term = _fak*ol(_i,0);
    ol(_i,Cart::xx) = PmB0*ol(_i,Cart::x) + nx[_i]*_fak*ol(i_less_x[_i],Cart::x) + term;
    ol(_i,Cart::xy) = PmB0*ol(_i,Cart::y) + nx[_i]*_fak*ol(i_less_x[_i],Cart::y);
    ol(_i,Cart::xz) = PmB0*ol(_i,Cart::z) + nx[_i]*_fak*ol(i_less_x[_i],Cart::z);
    ol(_i,Cart::yy) = PmB1*ol(_i,Cart::y) + ny[_i]*_fak*ol(i_less_y[_i],Cart::y) + term;
    ol(_i,Cart::yz) = PmB1*ol(_i,Cart::z) + ny[_i]*_fak*ol(i_less_y[_i],Cart::z);
    ol(_i,Cart::zz) = PmB2*ol(_i,Cart::z) + nz[_i]*_fak*ol(i_less_z[_i],Cart::z) + term;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

  //Integrals     s - f
  ol(0,Cart::xxx) = PmB0*ol(0,Cart::xx) + 2*_fak*ol(0,Cart::x);
  ol(0,Cart::xxy) = PmB1*ol(0,Cart::xx);
  ol(0,Cart::xxz) = PmB2*ol(0,Cart::xx);
  ol(0,Cart::xyy) = PmB0*ol(0,Cart::yy);
  ol(0,Cart::xyz) = PmB0*ol(0,Cart::yz);
  ol(0,Cart::xzz) = PmB0*ol(0,Cart::zz);
  ol(0,Cart::yyy) = PmB1*ol(0,Cart::yy) + 2*_fak*ol(0,Cart::y);
  ol(0,Cart::yyz) = PmB2*ol(0,Cart::yy);
  ol(0,Cart::yzz) = PmB1*ol(0,Cart::zz);
  ol(0,Cart::zzz) = PmB2*ol(0,Cart::zz) + 2*_fak*ol(0,Cart::z);
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    double term_x = 2*_fak*ol(_i,Cart::x);
    double term_y = 2*_fak*ol(_i,Cart::y);
    double term_z = 2*_fak*ol(_i,Cart::z);
    ol(_i,Cart::xxx) = PmB0*ol(_i,Cart::xx) + nx[_i]*_fak*ol(i_less_x[_i],Cart::xx) + term_x;
    ol(_i,Cart::xxy) = PmB1*ol(_i,Cart::xx) + ny[_i]*_fak*ol(i_less_y[_i],Cart::xx);
    ol(_i,Cart::xxz) = PmB2*ol(_i,Cart::xx) + nz[_i]*_fak*ol(i_less_z[_i],Cart::xx);
    ol(_i,Cart::xyy) = PmB0*ol(_i,Cart::yy) + nx[_i]*_fak*ol(i_less_x[_i],Cart::yy);
    ol(_i,Cart::xyz) = PmB0*ol(_i,Cart::yz) + nx[_i]*_fak*ol(i_less_x[_i],Cart::yz);
    ol(_i,Cart::xzz) = PmB0*ol(_i,Cart::zz) + nx[_i]*_fak*ol(i_less_x[_i],Cart::zz);
    ol(_i,Cart::yyy) = PmB1*ol(_i,Cart::yy) + ny[_i]*_fak*ol(i_less_y[_i],Cart::yy) + term_y;
    ol(_i,Cart::yyz) = PmB2*ol(_i,Cart::yy) + nz[_i]*_fak*ol(i_less_z[_i],Cart::yy);
    ol(_i,Cart::yzz) = PmB1*ol(_i,Cart::zz) + ny[_i]*_fak*ol(i_less_y[_i],Cart::zz);
    ol(_i,Cart::zzz) = PmB2*ol(_i,Cart::zz) + nz[_i]*_fak*ol(i_less_z[_i],Cart::zz) + term_z;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - g
  double term_xx = _fak*ol(0,Cart::xx);
  double term_yy = _fak*ol(0,Cart::yy);
  double term_zz = _fak*ol(0,Cart::zz);
  ol(0,Cart::xxxx) = PmB0*ol(0,Cart::xxx) + 3*term_xx;
  ol(0,Cart::xxxy) = PmB1*ol(0,Cart::xxx);
  ol(0,Cart::xxxz) = PmB2*ol(0,Cart::xxx);
  ol(0,Cart::xxyy) = PmB0*ol(0,Cart::xyy) + term_yy;
  ol(0,Cart::xxyz) = PmB1*ol(0,Cart::xxz);
  ol(0,Cart::xxzz) = PmB0*ol(0,Cart::xzz) + term_zz;
  ol(0,Cart::xyyy) = PmB0*ol(0,Cart::yyy);
  ol(0,Cart::xyyz) = PmB0*ol(0,Cart::yyz);
  ol(0,Cart::xyzz) = PmB0*ol(0,Cart::yzz);
  ol(0,Cart::xzzz) = PmB0*ol(0,Cart::zzz);
  ol(0,Cart::yyyy) = PmB1*ol(0,Cart::yyy) + 3*term_yy;
  ol(0,Cart::yyyz) = PmB2*ol(0,Cart::yyy);
  ol(0,Cart::yyzz) = PmB1*ol(0,Cart::yzz) + term_zz;
  ol(0,Cart::yzzz) = PmB1*ol(0,Cart::zzz);
  ol(0,Cart::zzzz) = PmB2*ol(0,Cart::zzz) + 3*term_zz;
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    double term_xx = _fak*ol(_i,Cart::xx);
    double term_yy = _fak*ol(_i,Cart::yy);
    double term_zz = _fak*ol(_i,Cart::zz);
    ol(_i,Cart::xxxx) = PmB0*ol(_i,Cart::xxx) + nx[_i]*_fak*ol(i_less_x[_i],Cart::xxx) + 3*term_xx;
    ol(_i,Cart::xxxy) = PmB1*ol(_i,Cart::xxx) + ny[_i]*_fak*ol(i_less_y[_i],Cart::xxx);
    ol(_i,Cart::xxxz) = PmB2*ol(_i,Cart::xxx) + nz[_i]*_fak*ol(i_less_z[_i],Cart::xxx);
    ol(_i,Cart::xxyy) = PmB0*ol(_i,Cart::xyy) + nx[_i]*_fak*ol(i_less_x[_i],Cart::xyy) + term_yy;
    ol(_i,Cart::xxyz) = PmB1*ol(_i,Cart::xxz) + ny[_i]*_fak*ol(i_less_y[_i],Cart::xxz);
    ol(_i,Cart::xxzz) = PmB0*ol(_i,Cart::xzz) + nx[_i]*_fak*ol(i_less_x[_i],Cart::xzz) + term_zz;
    ol(_i,Cart::xyyy) = PmB0*ol(_i,Cart::yyy) + nx[_i]*_fak*ol(i_less_x[_i],Cart::yyy);
    ol(_i,Cart::xyyz) = PmB0*ol(_i,Cart::yyz) + nx[_i]*_fak*ol(i_less_x[_i],Cart::yyz);
    ol(_i,Cart::xyzz) = PmB0*ol(_i,Cart::yzz) + nx[_i]*_fak*ol(i_less_x[_i],Cart::yzz);
    ol(_i,Cart::xzzz) = PmB0*ol(_i,Cart::zzz) + nx[_i]*_fak*ol(i_less_x[_i],Cart::zzz);
    ol(_i,Cart::yyyy) = PmB1*ol(_i,Cart::yyy) + ny[_i]*_fak*ol(i_less_y[_i],Cart::yyy) + 3*term_yy;
    ol(_i,Cart::yyyz) = PmB2*ol(_i,Cart::yyy) + nz[_i]*_fak*ol(i_less_z[_i],Cart::yyy);
    ol(_i,Cart::yyzz) = PmB1*ol(_i,Cart::yzz) + ny[_i]*_fak*ol(i_less_y[_i],Cart::yzz) + term_zz;
    ol(_i,Cart::yzzz) = PmB1*ol(_i,Cart::zzz) + ny[_i]*_fak*ol(i_less_y[_i],Cart::zzz);
    ol(_i,Cart::zzzz) = PmB2*ol(_i,Cart::zzz) + nz[_i]*_fak*ol(i_less_z[_i],Cart::zzz) + 3*term_zz;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)




//Integrals     p - s
if (_lmax_row > 0) {
  kin(Cart::x,0) = xi2*ol(Cart::x,0) + PmA0*kin(0,0);
  kin(Cart::y,0) = xi2*ol(Cart::y,0) + PmA1*kin(0,0);
  kin(Cart::z,0) = xi2*ol(Cart::z,0) + PmA2*kin(0,0);
}
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 1) {
  double term = _fak*kin(0,0)-_fak_a*ol(0,0);
  kin(Cart::xx,0) = xi2*ol(Cart::xx,0) + PmA0*kin(Cart::x,0) + term;
  kin(Cart::xy,0) = xi2*ol(Cart::xy,0) + PmA0*kin(Cart::y,0);
  kin(Cart::xz,0) = xi2*ol(Cart::xz,0) + PmA0*kin(Cart::z,0);
  kin(Cart::yy,0) = xi2*ol(Cart::yy,0) + PmA1*kin(Cart::y,0) + term;
  kin(Cart::yz,0) = xi2*ol(Cart::yz,0) + PmA1*kin(Cart::z,0);
  kin(Cart::zz,0) = xi2*ol(Cart::zz,0) + PmA2*kin(Cart::z,0) + term;
}
//------------------------------------------------------

//Integrals     f - s
if (_lmax_row > 2) {
  kin(Cart::xxx,0) = xi2*ol(Cart::xxx,0) + PmA0*kin(Cart::xx,0) + 2*(_fak*kin(Cart::x,0)-_fak_a*ol(Cart::x,0));
  kin(Cart::xxy,0) = xi2*ol(Cart::xxy,0) + PmA1*kin(Cart::xx,0);
  kin(Cart::xxz,0) = xi2*ol(Cart::xxz,0) + PmA2*kin(Cart::xx,0);
  kin(Cart::xyy,0) = xi2*ol(Cart::xyy,0) + PmA0*kin(Cart::yy,0);
  kin(Cart::xyz,0) = xi2*ol(Cart::xyz,0) + PmA0*kin(Cart::yz,0);
  kin(Cart::xzz,0) = xi2*ol(Cart::xzz,0) + PmA0*kin(Cart::zz,0);
  kin(Cart::yyy,0) = xi2*ol(Cart::yyy,0) + PmA1*kin(Cart::yy,0) + 2*(_fak*kin(Cart::y,0)-_fak_a*ol(Cart::y,0));
  kin(Cart::yyz,0) = xi2*ol(Cart::yyz,0) + PmA2*kin(Cart::yy,0);
  kin(Cart::yzz,0) = xi2*ol(Cart::yzz,0) + PmA1*kin(Cart::zz,0);
  kin(Cart::zzz,0) = xi2*ol(Cart::zzz,0) + PmA2*kin(Cart::zz,0) + 2*(_fak*kin(Cart::z,0)-_fak_a*ol(Cart::z,0));
}
//------------------------------------------------------

//Integrals     g - s
if (_lmax_row > 3) {
  double term_xx = _fak*kin(Cart::xx,0)-_fak_a*ol(Cart::xx,0);
  double term_yy = _fak*kin(Cart::yy,0)-_fak_a*ol(Cart::yy,0);
  double term_zz = _fak*kin(Cart::zz,0)-_fak_a*ol(Cart::zz,0);
  kin(Cart::xxxx,0) = xi2*ol(Cart::xxxx,0) + PmA0*kin(Cart::xxx,0) + 3*term_xx;
  kin(Cart::xxxy,0) = xi2*ol(Cart::xxxy,0) + PmA1*kin(Cart::xxx,0);
  kin(Cart::xxxz,0) = xi2*ol(Cart::xxxz,0) + PmA2*kin(Cart::xxx,0);
  kin(Cart::xxyy,0) = xi2*ol(Cart::xxyy,0) + PmA0*kin(Cart::xyy,0) + term_yy;
  kin(Cart::xxyz,0) = xi2*ol(Cart::xxyz,0) + PmA1*kin(Cart::xxz,0);
  kin(Cart::xxzz,0) = xi2*ol(Cart::xxzz,0) + PmA0*kin(Cart::xzz,0) + term_zz;
  kin(Cart::xyyy,0) = xi2*ol(Cart::xyyy,0) + PmA0*kin(Cart::yyy,0);
  kin(Cart::xyyz,0) = xi2*ol(Cart::xyyz,0) + PmA0*kin(Cart::yyz,0);
  kin(Cart::xyzz,0) = xi2*ol(Cart::xyzz,0) + PmA0*kin(Cart::yzz,0);
  kin(Cart::xzzz,0) = xi2*ol(Cart::xzzz,0) + PmA0*kin(Cart::zzz,0);
  kin(Cart::yyyy,0) = xi2*ol(Cart::yyyy,0) + PmA1*kin(Cart::yyy,0) + 3*term_yy;
  kin(Cart::yyyz,0) = xi2*ol(Cart::yyyz,0) + PmA2*kin(Cart::yyy,0);
  kin(Cart::yyzz,0) = xi2*ol(Cart::yyzz,0) + PmA1*kin(Cart::yzz,0) + term_zz;
  kin(Cart::yzzz,0) = xi2*ol(Cart::yzzz,0) + PmA1*kin(Cart::zzz,0);
  kin(Cart::zzzz,0) = xi2*ol(Cart::zzzz,0) + PmA2*kin(Cart::zzz,0) + 3*term_zz;
}
//------------------------------------------------------



if (_lmax_col > 0) {

  //Integrals     s - p
  kin(0,Cart::x) = xi2*ol(0,Cart::x) + PmB0*kin(0,0);
  kin(0,Cart::y) = xi2*ol(0,Cart::y) + PmB1*kin(0,0);
  kin(0,Cart::z) = xi2*ol(0,Cart::z) + PmB2*kin(0,0);
  //------------------------------------------------------

  //Integrals     p - p     d - p     f - p     g - p
  for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
    kin(_i,Cart::x) = xi2*ol(_i,Cart::x) + PmB0*kin(_i,0) + nx[_i]*_fak*kin(i_less_x[_i],0);
    kin(_i,Cart::y) = xi2*ol(_i,Cart::y) + PmB1*kin(_i,0) + ny[_i]*_fak*kin(i_less_y[_i],0);
    kin(_i,Cart::z) = xi2*ol(_i,Cart::z) + PmB2*kin(_i,0) + nz[_i]*_fak*kin(i_less_z[_i],0);
  }
  //------------------------------------------------------

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

  //Integrals     s - d
  double term = _fak*kin(0,0)-_fak_b*ol(0,0);
  kin(0,Cart::xx) = xi2*ol(0,Cart::xx) + PmB0*kin(0,Cart::x) + term;
  kin(0,Cart::xy) = xi2*ol(0,Cart::xy) + PmB0*kin(0,Cart::y);
  kin(0,Cart::xz) = xi2*ol(0,Cart::xz) + PmB0*kin(0,Cart::z);
  kin(0,Cart::yy) = xi2*ol(0,Cart::yy) + PmB1*kin(0,Cart::y) + term;
  kin(0,Cart::yz) = xi2*ol(0,Cart::yz) + PmB1*kin(0,Cart::z);
  kin(0,Cart::zz) = xi2*ol(0,Cart::zz) + PmB2*kin(0,Cart::z) + term;
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d
  for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
    double term = _fak*kin(_i,0)-_fak_b*ol(_i,0);
    kin(_i,Cart::xx) = xi2*ol(_i,Cart::xx) + PmB0*kin(_i,Cart::x) + nx[_i]*_fak*kin(i_less_x[_i],Cart::x) + term;
    kin(_i,Cart::xy) = xi2*ol(_i,Cart::xy) + PmB0*kin(_i,Cart::y) + nx[_i]*_fak*kin(i_less_x[_i],Cart::y);
    kin(_i,Cart::xz) = xi2*ol(_i,Cart::xz) + PmB0*kin(_i,Cart::z) + nx[_i]*_fak*kin(i_less_x[_i],Cart::z);
    kin(_i,Cart::yy) = xi2*ol(_i,Cart::yy) + PmB1*kin(_i,Cart::y) + ny[_i]*_fak*kin(i_less_y[_i],Cart::y) + term;
    kin(_i,Cart::yz) = xi2*ol(_i,Cart::yz) + PmB1*kin(_i,Cart::z) + ny[_i]*_fak*kin(i_less_y[_i],Cart::z);
    kin(_i,Cart::zz) = xi2*ol(_i,Cart::zz) + PmB2*kin(_i,Cart::z) + nz[_i]*_fak*kin(i_less_z[_i],Cart::z) + term;
//    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

  //Integrals     s - f
  kin(0,Cart::xxx) = xi2*ol(0,Cart::xxx) + PmB0*kin(0,Cart::xx) + 2*(_fak*kin(0,Cart::x)-_fak_b*ol(0,Cart::x));
  kin(0,Cart::xxy) = xi2*ol(0,Cart::xxy) + PmB1*kin(0,Cart::xx);
  kin(0,Cart::xxz) = xi2*ol(0,Cart::xxz) + PmB2*kin(0,Cart::xx);
  kin(0,Cart::xyy) = xi2*ol(0,Cart::xyy) + PmB0*kin(0,Cart::yy);
  kin(0,Cart::xyz) = xi2*ol(0,Cart::xyz) + PmB0*kin(0,Cart::yz);
  kin(0,Cart::xzz) = xi2*ol(0,Cart::xzz) + PmB0*kin(0,Cart::zz);
  kin(0,Cart::yyy) = xi2*ol(0,Cart::yyy) + PmB1*kin(0,Cart::yy) + 2*(_fak*kin(0,Cart::y)-_fak_b*ol(0,Cart::y));
  kin(0,Cart::yyz) = xi2*ol(0,Cart::yyz) + PmB2*kin(0,Cart::yy);
  kin(0,Cart::yzz) = xi2*ol(0,Cart::yzz) + PmB1*kin(0,Cart::zz);
  kin(0,Cart::zzz) = xi2*ol(0,Cart::zzz) + PmB2*kin(0,Cart::zz) + 2*(_fak*kin(0,Cart::z)-_fak_b*ol(0,Cart::z));
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f
  for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
    double term_x = 2*(_fak*kin(_i,Cart::x)-_fak_b*ol(_i,Cart::x));
    double term_y = 2*(_fak*kin(_i,Cart::y)-_fak_b*ol(_i,Cart::y));
    double term_z = 2*(_fak*kin(_i,Cart::z)-_fak_b*ol(_i,Cart::z));
    kin(_i,Cart::xxx) = xi2*ol(_i,Cart::xxx) + PmB0*kin(_i,Cart::xx) + nx[_i]*_fak*kin(i_less_x[_i],Cart::xx) + term_x;
    kin(_i,Cart::xxy) = xi2*ol(_i,Cart::xxy) + PmB1*kin(_i,Cart::xx) + ny[_i]*_fak*kin(i_less_y[_i],Cart::xx);
    kin(_i,Cart::xxz) = xi2*ol(_i,Cart::xxz) + PmB2*kin(_i,Cart::xx) + nz[_i]*_fak*kin(i_less_z[_i],Cart::xx);
    kin(_i,Cart::xyy) = xi2*ol(_i,Cart::xyy) + PmB0*kin(_i,Cart::yy) + nx[_i]*_fak*kin(i_less_x[_i],Cart::yy);
    kin(_i,Cart::xyz) = xi2*ol(_i,Cart::xyz) + PmB0*kin(_i,Cart::yz) + nx[_i]*_fak*kin(i_less_x[_i],Cart::yz);
    kin(_i,Cart::xzz) = xi2*ol(_i,Cart::xzz) + PmB0*kin(_i,Cart::zz) + nx[_i]*_fak*kin(i_less_x[_i],Cart::zz);
    kin(_i,Cart::yyy) = xi2*ol(_i,Cart::yyy) + PmB1*kin(_i,Cart::yy) + ny[_i]*_fak*kin(i_less_y[_i],Cart::yy) + term_y;
    kin(_i,Cart::yyz) = xi2*ol(_i,Cart::yyz) + PmB2*kin(_i,Cart::yy) + nz[_i]*_fak*kin(i_less_z[_i],Cart::yy);
    kin(_i,Cart::yzz) = xi2*ol(_i,Cart::yzz) + PmB1*kin(_i,Cart::zz) + ny[_i]*_fak*kin(i_less_y[_i],Cart::zz);
    kin(_i,Cart::zzz) = xi2*ol(_i,Cart::zzz) + PmB2*kin(_i,Cart::zz) + nz[_i]*_fak*kin(i_less_z[_i],Cart::zz) + term_z;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - g
  double term_xx = _fak*kin(0,Cart::xx)-_fak_b*ol(0,Cart::xx);
  double term_yy = _fak*kin(0,Cart::yy)-_fak_b*ol(0,Cart::yy);
  double term_zz = _fak*kin(0,Cart::zz)-_fak_b*ol(0,Cart::zz);
  kin(0,Cart::xxxx) = xi2*ol(0,Cart::xxxx) + PmB0*kin(0,Cart::xxx) + 3*term_xx;
  kin(0,Cart::xxxy) = xi2*ol(0,Cart::xxxy) + PmB1*kin(0,Cart::xxx);
  kin(0,Cart::xxxz) = xi2*ol(0,Cart::xxxz) + PmB2*kin(0,Cart::xxx);
  kin(0,Cart::xxyy) = xi2*ol(0,Cart::xxyy) + PmB0*kin(0,Cart::xyy) + term_yy;
  kin(0,Cart::xxyz) = xi2*ol(0,Cart::xxyz) + PmB1*kin(0,Cart::xxz);
  kin(0,Cart::xxzz) = xi2*ol(0,Cart::xxzz) + PmB0*kin(0,Cart::xzz) + term_zz;
  kin(0,Cart::xyyy) = xi2*ol(0,Cart::xyyy) + PmB0*kin(0,Cart::yyy);
  kin(0,Cart::xyyz) = xi2*ol(0,Cart::xyyz) + PmB0*kin(0,Cart::yyz);
  kin(0,Cart::xyzz) = xi2*ol(0,Cart::xyzz) + PmB0*kin(0,Cart::yzz);
  kin(0,Cart::xzzz) = xi2*ol(0,Cart::xzzz) + PmB0*kin(0,Cart::zzz);
  kin(0,Cart::yyyy) = xi2*ol(0,Cart::yyyy) + PmB1*kin(0,Cart::yyy) + 3*term_yy;
  kin(0,Cart::yyyz) = xi2*ol(0,Cart::yyyz) + PmB2*kin(0,Cart::yyy);
  kin(0,Cart::yyzz) = xi2*ol(0,Cart::yyzz) + PmB1*kin(0,Cart::yzz) + term_zz;
  kin(0,Cart::yzzz) = xi2*ol(0,Cart::yzzz) + PmB1*kin(0,Cart::zzz);
  kin(0,Cart::zzzz) = xi2*ol(0,Cart::zzzz) + PmB2*kin(0,Cart::zzz) + 3*term_zz;
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g
  for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
    double term_xx = _fak*kin(_i,Cart::xx)-_fak_b*ol(_i,Cart::xx);
    double term_yy = _fak*kin(_i,Cart::yy)-_fak_b*ol(_i,Cart::yy);
    double term_zz = _fak*kin(_i,Cart::zz)-_fak_b*ol(_i,Cart::zz);
    kin(_i,Cart::xxxx) = xi2*ol(_i,Cart::xxxx) + PmB0*kin(_i,Cart::xxx) + nx[_i]*_fak*kin(i_less_x[_i],Cart::xxx) + 3*term_xx;
    kin(_i,Cart::xxxy) = xi2*ol(_i,Cart::xxxy) + PmB1*kin(_i,Cart::xxx) + ny[_i]*_fak*kin(i_less_y[_i],Cart::xxx);
    kin(_i,Cart::xxxz) = xi2*ol(_i,Cart::xxxz) + PmB2*kin(_i,Cart::xxx) + nz[_i]*_fak*kin(i_less_z[_i],Cart::xxx);
    kin(_i,Cart::xxyy) = xi2*ol(_i,Cart::xxyy) + PmB0*kin(_i,Cart::xyy) + nx[_i]*_fak*kin(i_less_x[_i],Cart::xyy) + term_yy;
    kin(_i,Cart::xxyz) = xi2*ol(_i,Cart::xxyz) + PmB1*kin(_i,Cart::xxz) + ny[_i]*_fak*kin(i_less_y[_i],Cart::xxz);
    kin(_i,Cart::xxzz) = xi2*ol(_i,Cart::xxzz) + PmB0*kin(_i,Cart::xzz) + nx[_i]*_fak*kin(i_less_x[_i],Cart::xzz) + term_zz;
    kin(_i,Cart::xyyy) = xi2*ol(_i,Cart::xyyy) + PmB0*kin(_i,Cart::yyy) + nx[_i]*_fak*kin(i_less_x[_i],Cart::yyy);
    kin(_i,Cart::xyyz) = xi2*ol(_i,Cart::xyyz) + PmB0*kin(_i,Cart::yyz) + nx[_i]*_fak*kin(i_less_x[_i],Cart::yyz);
    kin(_i,Cart::xyzz) = xi2*ol(_i,Cart::xyzz) + PmB0*kin(_i,Cart::yzz) + nx[_i]*_fak*kin(i_less_x[_i],Cart::yzz);
    kin(_i,Cart::xzzz) = xi2*ol(_i,Cart::xzzz) + PmB0*kin(_i,Cart::zzz) + nx[_i]*_fak*kin(i_less_x[_i],Cart::zzz);
    kin(_i,Cart::yyyy) = xi2*ol(_i,Cart::yyyy) + PmB1*kin(_i,Cart::yyy) + ny[_i]*_fak*kin(i_less_y[_i],Cart::yyy) + 3*term_yy;
    kin(_i,Cart::yyyz) = xi2*ol(_i,Cart::yyyz) + PmB2*kin(_i,Cart::yyy) + nz[_i]*_fak*kin(i_less_z[_i],Cart::yyy);
    kin(_i,Cart::yyzz) = xi2*ol(_i,Cart::yyzz) + PmB1*kin(_i,Cart::yzz) + ny[_i]*_fak*kin(i_less_y[_i],Cart::yzz) + term_zz;
    kin(_i,Cart::yzzz) = xi2*ol(_i,Cart::yzzz) + PmB1*kin(_i,Cart::zzz) + ny[_i]*_fak*kin(i_less_y[_i],Cart::zzz);
    kin(_i,Cart::zzzz) = xi2*ol(_i,Cart::zzzz) + PmB2*kin(_i,Cart::zzz) + nz[_i]*_fak*kin(i_less_z[_i],Cart::zzz) + 3*term_zz;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)

                // normalization and cartesian -> spherical factors
             Eigen::MatrixXd _kin_sph = getTrafo(*itr)*kin*getTrafo(*itc).transpose();
        // save to _matrix
        
        for ( unsigned i = 0; i< _matrix.rows(); i++ ) {
            for (unsigned j = 0; j < _matrix.cols(); j++) {
                _matrix(i,j) += _kin_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
            }
        }
        
        
        
        
        
                }//col
            }//row
 return;
        }
    }
}
    
