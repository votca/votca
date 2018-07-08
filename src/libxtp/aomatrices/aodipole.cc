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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>
#include <string>
#include <vector>



namespace votca { namespace xtp {
    

    
    void AODipole::FillBlock( std::vector< Eigen::Block<Eigen::MatrixXd> >& _matrix,const AOShell* _shell_row,const AOShell* _shell_col) {

        
        /* Calculating the AO matrix of the gradient operator requires 
         * the raw overlap matrix (i.e. in unnormalized cartesians) 
        
         */

        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();
        
        if ( _lmax_col > 4 ) {
            std::cerr << "Momentum transition dipoles only implemented for S,P,D,F,G functions in DFT basis!" << std::flush;
            exit(1);
        }

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
    
        // initialize local matrix block for unnormalized cartesians
        std::vector< Eigen::MatrixXd > _dip;
        for (int _i_comp = 0; _i_comp < 3; _i_comp++){
            _dip.push_back(Eigen::MatrixXd::Zero(_nrows,_ncols));
        }
        
        // initialize local matrix block for unnormalized cartesians of overlap
        // int _ncols_ol = this->getBlockSize( _lmax_col +1 ); 
        
        Eigen::MatrixXd _ol =Eigen::MatrixXd::Zero(_nrows,_ncols);
        
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



        // some helpers
       
        // definition of a center around which the moment should be calculated
        tools::vec _center(0.0); // here: origin, can be changed later
        tools::vec  _pmc(0.0);
        
        
        // iterate over Gaussians in this _shell_row
        for ( AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
            // get decay constant
            const double _decay_row = itr->getDecay();
            
            for ( AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
                //get decay constant
                const double _decay_col = itc->getDecay();
        
       
                const double _fak  = 0.5/(_decay_row + _decay_col);
                const double _fak2 = 2.0 * _fak;

                double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
                // check if distance between postions is big, then skip step   
       
                if ( _exparg > 30.0 ) { continue; }
        


        const double PmA0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
        const double PmA1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
        const double PmA2 = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

        const double PmB0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
        const double PmB1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
        const double PmB2 = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();        
        
        
        _pmc= _fak2*(_decay_row * _pos_row + _decay_col * _pos_col)-_center;
   
        // calculate s-s- overlap matrix element
        _ol(0,0) = pow(4.0*_decay_row*_decay_col,0.75) * pow(_fak2,1.5)*exp(-_fak2 * _decay_row * _decay_col *_distsq); // s-s element

        // s-s dipole moment integrals
        for ( int _i_comp = 0 ; _i_comp < 3; _i_comp++ ){
            _dip[_i_comp](0,0) = _pmc[_i_comp]*_ol(0,0);
        }

//Integrals     p - s
if (_lmax_row > 0) {
  _ol(Cart::x,0) = PmA0*_ol(0,0);
  _ol(Cart::y,0) = PmA1*_ol(0,0);
  _ol(Cart::z,0) = PmA2*_ol(0,0);
}
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 1) {
  double term = _fak*_ol(0,0);
  _ol(Cart::xx,0) = PmA0*_ol(Cart::x,0) + term;
  _ol(Cart::xy,0) = PmA0*_ol(Cart::y,0);
  _ol(Cart::xz,0) = PmA0*_ol(Cart::z,0);
  _ol(Cart::yy,0) = PmA1*_ol(Cart::y,0) + term;
  _ol(Cart::yz,0) = PmA1*_ol(Cart::z,0);
  _ol(Cart::zz,0) = PmA2*_ol(Cart::z,0) + term;
}
//------------------------------------------------------

//Integrals     f - s
if (_lmax_row > 2) {
  _ol(Cart::xxx,0) = PmA0*_ol(Cart::xx,0) + 2*_fak*_ol(Cart::x,0);
  _ol(Cart::xxy,0) = PmA1*_ol(Cart::xx,0);
  _ol(Cart::xxz,0) = PmA2*_ol(Cart::xx,0);
  _ol(Cart::xyy,0) = PmA0*_ol(Cart::yy,0);
  _ol(Cart::xyz,0) = PmA0*_ol(Cart::yz,0);
  _ol(Cart::xzz,0) = PmA0*_ol(Cart::zz,0);
  _ol(Cart::yyy,0) = PmA1*_ol(Cart::yy,0) + 2*_fak*_ol(Cart::y,0);
  _ol(Cart::yyz,0) = PmA2*_ol(Cart::yy,0);
  _ol(Cart::yzz,0) = PmA1*_ol(Cart::zz,0);
  _ol(Cart::zzz,0) = PmA2*_ol(Cart::zz,0) + 2*_fak*_ol(Cart::z,0);
}
//------------------------------------------------------

//Integrals     g - s
if (_lmax_row > 3) {
  double term_xx = _fak*_ol(Cart::xx,0);
  double term_yy = _fak*_ol(Cart::yy,0);
  double term_zz = _fak*_ol(Cart::zz,0);
  _ol(Cart::xxxx,0) = PmA0*_ol(Cart::xxx,0) + 3*term_xx;
  _ol(Cart::xxxy,0) = PmA1*_ol(Cart::xxx,0);
  _ol(Cart::xxxz,0) = PmA2*_ol(Cart::xxx,0);
  _ol(Cart::xxyy,0) = PmA0*_ol(Cart::xyy,0) + term_yy;
  _ol(Cart::xxyz,0) = PmA1*_ol(Cart::xxz,0);
  _ol(Cart::xxzz,0) = PmA0*_ol(Cart::xzz,0) + term_zz;
  _ol(Cart::xyyy,0) = PmA0*_ol(Cart::yyy,0);
  _ol(Cart::xyyz,0) = PmA0*_ol(Cart::yyz,0);
  _ol(Cart::xyzz,0) = PmA0*_ol(Cart::yzz,0);
  _ol(Cart::xzzz,0) = PmA0*_ol(Cart::zzz,0);
  _ol(Cart::yyyy,0) = PmA1*_ol(Cart::yyy,0) + 3*term_yy;
  _ol(Cart::yyyz,0) = PmA2*_ol(Cart::yyy,0);
  _ol(Cart::yyzz,0) = PmA1*_ol(Cart::yzz,0) + term_zz;
  _ol(Cart::yzzz,0) = PmA1*_ol(Cart::zzz,0);
  _ol(Cart::zzzz,0) = PmA2*_ol(Cart::zzz,0) + 3*term_zz;
}
//------------------------------------------------------



if (_lmax_col > 0) {

  //Integrals     s - p
  _ol(0,Cart::x) = PmB0*_ol(0,0);
  _ol(0,Cart::y) = PmB1*_ol(0,0);
  _ol(0,Cart::z) = PmB2*_ol(0,0);
  //------------------------------------------------------

  //Integrals     p - p     d - p     f - p     g - p
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    _ol(_i,Cart::x) = PmB0*_ol(_i,0) + nx[_i]*_fak*_ol(i_less_x[_i],0);
    _ol(_i,Cart::y) = PmB1*_ol(_i,0) + ny[_i]*_fak*_ol(i_less_y[_i],0);
    _ol(_i,Cart::z) = PmB2*_ol(_i,0) + nz[_i]*_fak*_ol(i_less_z[_i],0);
  }
  //------------------------------------------------------

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

  //Integrals     s - d
  double term = _fak*_ol(0,0);
  _ol(0,Cart::xx) = PmB0*_ol(0,Cart::x) + term;
  _ol(0,Cart::xy) = PmB0*_ol(0,Cart::y);
  _ol(0,Cart::xz) = PmB0*_ol(0,Cart::z);
  _ol(0,Cart::yy) = PmB1*_ol(0,Cart::y) + term;
  _ol(0,Cart::yz) = PmB1*_ol(0,Cart::z);
  _ol(0,Cart::zz) = PmB2*_ol(0,Cart::z) + term;
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    double term = _fak*_ol(_i,0);
    _ol(_i,Cart::xx) = PmB0*_ol(_i,Cart::x) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::x) + term;
    _ol(_i,Cart::xy) = PmB0*_ol(_i,Cart::y) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::y);
    _ol(_i,Cart::xz) = PmB0*_ol(_i,Cart::z) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::z);
    _ol(_i,Cart::yy) = PmB1*_ol(_i,Cart::y) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::y) + term;
    _ol(_i,Cart::yz) = PmB1*_ol(_i,Cart::z) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::z);
    _ol(_i,Cart::zz) = PmB2*_ol(_i,Cart::z) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::z) + term;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

  //Integrals     s - f
  _ol(0,Cart::xxx) = PmB0*_ol(0,Cart::xx) + 2*_fak*_ol(0,Cart::x);
  _ol(0,Cart::xxy) = PmB1*_ol(0,Cart::xx);
  _ol(0,Cart::xxz) = PmB2*_ol(0,Cart::xx);
  _ol(0,Cart::xyy) = PmB0*_ol(0,Cart::yy);
  _ol(0,Cart::xyz) = PmB0*_ol(0,Cart::yz);
  _ol(0,Cart::xzz) = PmB0*_ol(0,Cart::zz);
  _ol(0,Cart::yyy) = PmB1*_ol(0,Cart::yy) + 2*_fak*_ol(0,Cart::y);
  _ol(0,Cart::yyz) = PmB2*_ol(0,Cart::yy);
  _ol(0,Cart::yzz) = PmB1*_ol(0,Cart::zz);
  _ol(0,Cart::zzz) = PmB2*_ol(0,Cart::zz) + 2*_fak*_ol(0,Cart::z);
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    double term_x = 2*_fak*_ol(_i,Cart::x);
    double term_y = 2*_fak*_ol(_i,Cart::y);
    double term_z = 2*_fak*_ol(_i,Cart::z);
    _ol(_i,Cart::xxx) = PmB0*_ol(_i,Cart::xx) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::xx) + term_x;
    _ol(_i,Cart::xxy) = PmB1*_ol(_i,Cart::xx) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::xx);
    _ol(_i,Cart::xxz) = PmB2*_ol(_i,Cart::xx) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::xx);
    _ol(_i,Cart::xyy) = PmB0*_ol(_i,Cart::yy) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::yy);
    _ol(_i,Cart::xyz) = PmB0*_ol(_i,Cart::yz) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::yz);
    _ol(_i,Cart::xzz) = PmB0*_ol(_i,Cart::zz) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::zz);
    _ol(_i,Cart::yyy) = PmB1*_ol(_i,Cart::yy) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::yy) + term_y;
    _ol(_i,Cart::yyz) = PmB2*_ol(_i,Cart::yy) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::yy);
    _ol(_i,Cart::yzz) = PmB1*_ol(_i,Cart::zz) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::zz);
    _ol(_i,Cart::zzz) = PmB2*_ol(_i,Cart::zz) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::zz) + term_z;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - g
  double term_xx = _fak*_ol(0,Cart::xx);
  double term_yy = _fak*_ol(0,Cart::yy);
  double term_zz = _fak*_ol(0,Cart::zz);
  _ol(0,Cart::xxxx) = PmB0*_ol(0,Cart::xxx) + 3*term_xx;
  _ol(0,Cart::xxxy) = PmB1*_ol(0,Cart::xxx);
  _ol(0,Cart::xxxz) = PmB2*_ol(0,Cart::xxx);
  _ol(0,Cart::xxyy) = PmB0*_ol(0,Cart::xyy) + term_yy;
  _ol(0,Cart::xxyz) = PmB1*_ol(0,Cart::xxz);
  _ol(0,Cart::xxzz) = PmB0*_ol(0,Cart::xzz) + term_zz;
  _ol(0,Cart::xyyy) = PmB0*_ol(0,Cart::yyy);
  _ol(0,Cart::xyyz) = PmB0*_ol(0,Cart::yyz);
  _ol(0,Cart::xyzz) = PmB0*_ol(0,Cart::yzz);
  _ol(0,Cart::xzzz) = PmB0*_ol(0,Cart::zzz);
  _ol(0,Cart::yyyy) = PmB1*_ol(0,Cart::yyy) + 3*term_yy;
  _ol(0,Cart::yyyz) = PmB2*_ol(0,Cart::yyy);
  _ol(0,Cart::yyzz) = PmB1*_ol(0,Cart::yzz) + term_zz;
  _ol(0,Cart::yzzz) = PmB1*_ol(0,Cart::zzz);
  _ol(0,Cart::zzzz) = PmB2*_ol(0,Cart::zzz) + 3*term_zz;
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    double term_xx = _fak*_ol(_i,Cart::xx);
    double term_yy = _fak*_ol(_i,Cart::yy);
    double term_zz = _fak*_ol(_i,Cart::zz);
    _ol(_i,Cart::xxxx) = PmB0*_ol(_i,Cart::xxx) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::xxx) + 3*term_xx;
    _ol(_i,Cart::xxxy) = PmB1*_ol(_i,Cart::xxx) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::xxx);
    _ol(_i,Cart::xxxz) = PmB2*_ol(_i,Cart::xxx) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::xxx);
    _ol(_i,Cart::xxyy) = PmB0*_ol(_i,Cart::xyy) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::xyy) + term_yy;
    _ol(_i,Cart::xxyz) = PmB1*_ol(_i,Cart::xxz) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::xxz);
    _ol(_i,Cart::xxzz) = PmB0*_ol(_i,Cart::xzz) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::xzz) + term_zz;
    _ol(_i,Cart::xyyy) = PmB0*_ol(_i,Cart::yyy) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::yyy);
    _ol(_i,Cart::xyyz) = PmB0*_ol(_i,Cart::yyz) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::yyz);
    _ol(_i,Cart::xyzz) = PmB0*_ol(_i,Cart::yzz) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::yzz);
    _ol(_i,Cart::xzzz) = PmB0*_ol(_i,Cart::zzz) + nx[_i]*_fak*_ol(i_less_x[_i],Cart::zzz);
    _ol(_i,Cart::yyyy) = PmB1*_ol(_i,Cart::yyy) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::yyy) + 3*term_yy;
    _ol(_i,Cart::yyyz) = PmB2*_ol(_i,Cart::yyy) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::yyy);
    _ol(_i,Cart::yyzz) = PmB1*_ol(_i,Cart::yzz) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::yzz) + term_zz;
    _ol(_i,Cart::yzzz) = PmB1*_ol(_i,Cart::zzz) + ny[_i]*_fak*_ol(i_less_y[_i],Cart::zzz);
    _ol(_i,Cart::zzzz) = PmB2*_ol(_i,Cart::zzz) + nz[_i]*_fak*_ol(i_less_z[_i],Cart::zzz) + 3*term_zz;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)




//Integrals     p - s
if (_lmax_row > 0) {
  for (int _k =  0; _k < 3; _k++) {
    _dip[_k](Cart::x,0) = PmA0*_dip[_k](0,0) + (_k==0)*_fak*_ol(0,0);
    _dip[_k](Cart::y,0) = PmA1*_dip[_k](0,0) + (_k==1)*_fak*_ol(0,0);
    _dip[_k](Cart::z,0) = PmA2*_dip[_k](0,0) + (_k==2)*_fak*_ol(0,0);
  }
}
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 1) {
  for (int _k =  0; _k < 3; _k++) {
    double term = _fak*_dip[_k](0,0);
    _dip[_k](Cart::xx,0) = PmA0*_dip[_k](Cart::x,0) + (_k==0)*_fak*_ol(Cart::x,0) + term;
    _dip[_k](Cart::xy,0) = PmA0*_dip[_k](Cart::y,0) + (_k==0)*_fak*_ol(Cart::y,0);
    _dip[_k](Cart::xz,0) = PmA0*_dip[_k](Cart::z,0) + (_k==0)*_fak*_ol(Cart::z,0);
    _dip[_k](Cart::yy,0) = PmA1*_dip[_k](Cart::y,0) + (_k==1)*_fak*_ol(Cart::y,0) + term;
    _dip[_k](Cart::yz,0) = PmA1*_dip[_k](Cart::z,0) + (_k==1)*_fak*_ol(Cart::z,0);
    _dip[_k](Cart::zz,0) = PmA2*_dip[_k](Cart::z,0) + (_k==2)*_fak*_ol(Cart::z,0) + term;
  }
}
//------------------------------------------------------

//Integrals     f - s
if (_lmax_row > 2) {
  for (int _k =  0; _k < 3; _k++) {
    _dip[_k](Cart::xxx,0) = PmA0*_dip[_k](Cart::xx,0) + (_k==0)*_fak*_ol(Cart::xx,0) + 2*_fak*_dip[_k](Cart::x,0);
    _dip[_k](Cart::xxy,0) = PmA1*_dip[_k](Cart::xx,0) + (_k==1)*_fak*_ol(Cart::xx,0);
    _dip[_k](Cart::xxz,0) = PmA2*_dip[_k](Cart::xx,0) + (_k==2)*_fak*_ol(Cart::xx,0);
    _dip[_k](Cart::xyy,0) = PmA0*_dip[_k](Cart::yy,0) + (_k==0)*_fak*_ol(Cart::yy,0);
    _dip[_k](Cart::xyz,0) = PmA0*_dip[_k](Cart::yz,0) + (_k==0)*_fak*_ol(Cart::yz,0);
    _dip[_k](Cart::xzz,0) = PmA0*_dip[_k](Cart::zz,0) + (_k==0)*_fak*_ol(Cart::zz,0);
    _dip[_k](Cart::yyy,0) = PmA1*_dip[_k](Cart::yy,0) + (_k==1)*_fak*_ol(Cart::yy,0) + 2*_fak*_dip[_k](Cart::y,0);
    _dip[_k](Cart::yyz,0) = PmA2*_dip[_k](Cart::yy,0) + (_k==2)*_fak*_ol(Cart::yy,0);
    _dip[_k](Cart::yzz,0) = PmA1*_dip[_k](Cart::zz,0) + (_k==1)*_fak*_ol(Cart::zz,0);
    _dip[_k](Cart::zzz,0) = PmA2*_dip[_k](Cart::zz,0) + (_k==2)*_fak*_ol(Cart::zz,0) + 2*_fak*_dip[_k](Cart::z,0);
  }
}
//------------------------------------------------------

//Integrals     g - s
if (_lmax_row > 3) {
  for (int _k =  0; _k < 3; _k++) {
    double term_xx = _fak*_dip[_k](Cart::xx,0);
    double term_yy = _fak*_dip[_k](Cart::yy,0);
    double term_zz = _fak*_dip[_k](Cart::zz,0);
    _dip[_k](Cart::xxxx,0) = PmA0*_dip[_k](Cart::xxx,0) + (_k==0)*_fak*_ol(Cart::xxx,0) + 3*term_xx;
    _dip[_k](Cart::xxxy,0) = PmA1*_dip[_k](Cart::xxx,0) + (_k==1)*_fak*_ol(Cart::xxx,0);
    _dip[_k](Cart::xxxz,0) = PmA2*_dip[_k](Cart::xxx,0) + (_k==2)*_fak*_ol(Cart::xxx,0);
    _dip[_k](Cart::xxyy,0) = PmA0*_dip[_k](Cart::xyy,0) + (_k==0)*_fak*_ol(Cart::xyy,0) + term_yy;
    _dip[_k](Cart::xxyz,0) = PmA1*_dip[_k](Cart::xxz,0) + (_k==1)*_fak*_ol(Cart::xxz,0);
    _dip[_k](Cart::xxzz,0) = PmA0*_dip[_k](Cart::xzz,0) + (_k==0)*_fak*_ol(Cart::xzz,0) + term_zz;
    _dip[_k](Cart::xyyy,0) = PmA0*_dip[_k](Cart::yyy,0) + (_k==0)*_fak*_ol(Cart::yyy,0);
    _dip[_k](Cart::xyyz,0) = PmA0*_dip[_k](Cart::yyz,0) + (_k==0)*_fak*_ol(Cart::yyz,0);
    _dip[_k](Cart::xyzz,0) = PmA0*_dip[_k](Cart::yzz,0) + (_k==0)*_fak*_ol(Cart::yzz,0);
    _dip[_k](Cart::xzzz,0) = PmA0*_dip[_k](Cart::zzz,0) + (_k==0)*_fak*_ol(Cart::zzz,0);
    _dip[_k](Cart::yyyy,0) = PmA1*_dip[_k](Cart::yyy,0) + (_k==1)*_fak*_ol(Cart::yyy,0) + 3*term_yy;
    _dip[_k](Cart::yyyz,0) = PmA2*_dip[_k](Cart::yyy,0) + (_k==2)*_fak*_ol(Cart::yyy,0);
    _dip[_k](Cart::yyzz,0) = PmA1*_dip[_k](Cart::yzz,0) + (_k==1)*_fak*_ol(Cart::yzz,0) + term_zz;
    _dip[_k](Cart::yzzz,0) = PmA1*_dip[_k](Cart::zzz,0) + (_k==1)*_fak*_ol(Cart::zzz,0);
    _dip[_k](Cart::zzzz,0) = PmA2*_dip[_k](Cart::zzz,0) + (_k==2)*_fak*_ol(Cart::zzz,0) + 3*term_zz;
  }
}
//------------------------------------------------------



if (_lmax_col > 0) {

  //Integrals     s - p
  for (int _k =  0; _k < 3; _k++) {
    _dip[_k](0,Cart::x) = PmB0*_dip[_k](0,0) + (_k==0)*_fak*_ol(0,0);
    _dip[_k](0,Cart::y) = PmB1*_dip[_k](0,0) + (_k==1)*_fak*_ol(0,0);
    _dip[_k](0,Cart::z) = PmB2*_dip[_k](0,0) + (_k==2)*_fak*_ol(0,0);
  }
  //------------------------------------------------------

  //Integrals     p - p     d - p     f - p     g - p
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    for (int _k =  0; _k < 3; _k++) {
      _dip[_k](_i,Cart::x) = PmB0*_dip[_k](_i,0) + (_k==0)*_fak*_ol(_i,0) + nx[_i]*_fak*_dip[_k](i_less_x[_i],0);
      _dip[_k](_i,Cart::y) = PmB1*_dip[_k](_i,0) + (_k==1)*_fak*_ol(_i,0) + ny[_i]*_fak*_dip[_k](i_less_y[_i],0);
      _dip[_k](_i,Cart::z) = PmB2*_dip[_k](_i,0) + (_k==2)*_fak*_ol(_i,0) + nz[_i]*_fak*_dip[_k](i_less_z[_i],0);
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

  //Integrals     s - d
  for (int _k =  0; _k < 3; _k++) {
    double term = _fak*_dip[_k](0,0);
    _dip[_k](0,Cart::xx) = PmB0*_dip[_k](0,Cart::x) + (_k==0)*_fak*_ol(0,Cart::x) + term;
    _dip[_k](0,Cart::xy) = PmB0*_dip[_k](0,Cart::y) + (_k==0)*_fak*_ol(0,Cart::y);
    _dip[_k](0,Cart::xz) = PmB0*_dip[_k](0,Cart::z) + (_k==0)*_fak*_ol(0,Cart::z);
    _dip[_k](0,Cart::yy) = PmB1*_dip[_k](0,Cart::y) + (_k==1)*_fak*_ol(0,Cart::y) + term;
    _dip[_k](0,Cart::yz) = PmB1*_dip[_k](0,Cart::z) + (_k==1)*_fak*_ol(0,Cart::z);
    _dip[_k](0,Cart::zz) = PmB2*_dip[_k](0,Cart::z) + (_k==2)*_fak*_ol(0,Cart::z) + term;
  }
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    for (int _k =  0; _k < 3; _k++) {
      double term = _fak*_dip[_k](_i,0);
      _dip[_k](_i,Cart::xx) = PmB0*_dip[_k](_i,Cart::x) + (_k==0)*_fak*_ol(_i,Cart::x) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::x) + term;
      _dip[_k](_i,Cart::xy) = PmB0*_dip[_k](_i,Cart::y) + (_k==0)*_fak*_ol(_i,Cart::y) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::y);
      _dip[_k](_i,Cart::xz) = PmB0*_dip[_k](_i,Cart::z) + (_k==0)*_fak*_ol(_i,Cart::z) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::z);
      _dip[_k](_i,Cart::yy) = PmB1*_dip[_k](_i,Cart::y) + (_k==1)*_fak*_ol(_i,Cart::y) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::y) + term;
      _dip[_k](_i,Cart::yz) = PmB1*_dip[_k](_i,Cart::z) + (_k==1)*_fak*_ol(_i,Cart::z) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::z);
      _dip[_k](_i,Cart::zz) = PmB2*_dip[_k](_i,Cart::z) + (_k==2)*_fak*_ol(_i,Cart::z) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::z) + term;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

  //Integrals     s - f
  for (int _k =  0; _k < 3; _k++) {
    _dip[_k](0,Cart::xxx) = PmB0*_dip[_k](0,Cart::xx) + (_k==0)*_fak*_ol(0,Cart::xx) + 2*_fak*_dip[_k](0,Cart::x);
    _dip[_k](0,Cart::xxy) = PmB1*_dip[_k](0,Cart::xx) + (_k==1)*_fak*_ol(0,Cart::xx);
    _dip[_k](0,Cart::xxz) = PmB2*_dip[_k](0,Cart::xx) + (_k==2)*_fak*_ol(0,Cart::xx);
    _dip[_k](0,Cart::xyy) = PmB0*_dip[_k](0,Cart::yy) + (_k==0)*_fak*_ol(0,Cart::yy);
    _dip[_k](0,Cart::xyz) = PmB0*_dip[_k](0,Cart::yz) + (_k==0)*_fak*_ol(0,Cart::yz);
    _dip[_k](0,Cart::xzz) = PmB0*_dip[_k](0,Cart::zz) + (_k==0)*_fak*_ol(0,Cart::zz);
    _dip[_k](0,Cart::yyy) = PmB1*_dip[_k](0,Cart::yy) + (_k==1)*_fak*_ol(0,Cart::yy) + 2*_fak*_dip[_k](0,Cart::y);
    _dip[_k](0,Cart::yyz) = PmB2*_dip[_k](0,Cart::yy) + (_k==2)*_fak*_ol(0,Cart::yy);
    _dip[_k](0,Cart::yzz) = PmB1*_dip[_k](0,Cart::zz) + (_k==1)*_fak*_ol(0,Cart::zz);
    _dip[_k](0,Cart::zzz) = PmB2*_dip[_k](0,Cart::zz) + (_k==2)*_fak*_ol(0,Cart::zz) + 2*_fak*_dip[_k](0,Cart::z);
  }
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    for (int _k =  0; _k < 3; _k++) {
      double term_x = 2*_fak*_dip[_k](_i,Cart::x);
      double term_y = 2*_fak*_dip[_k](_i,Cart::y);
      double term_z = 2*_fak*_dip[_k](_i,Cart::z);
      _dip[_k](_i,Cart::xxx) = PmB0*_dip[_k](_i,Cart::xx) + (_k==0)*_fak*_ol(_i,Cart::xx) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::xx) + term_x;
      _dip[_k](_i,Cart::xxy) = PmB1*_dip[_k](_i,Cart::xx) + (_k==1)*_fak*_ol(_i,Cart::xx) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::xx);
      _dip[_k](_i,Cart::xxz) = PmB2*_dip[_k](_i,Cart::xx) + (_k==2)*_fak*_ol(_i,Cart::xx) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::xx);
      _dip[_k](_i,Cart::xyy) = PmB0*_dip[_k](_i,Cart::yy) + (_k==0)*_fak*_ol(_i,Cart::yy) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::yy);
      _dip[_k](_i,Cart::xyz) = PmB0*_dip[_k](_i,Cart::yz) + (_k==0)*_fak*_ol(_i,Cart::yz) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::yz);
      _dip[_k](_i,Cart::xzz) = PmB0*_dip[_k](_i,Cart::zz) + (_k==0)*_fak*_ol(_i,Cart::zz) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::zz);
      _dip[_k](_i,Cart::yyy) = PmB1*_dip[_k](_i,Cart::yy) + (_k==1)*_fak*_ol(_i,Cart::yy) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::yy) + term_y;
      _dip[_k](_i,Cart::yyz) = PmB2*_dip[_k](_i,Cart::yy) + (_k==2)*_fak*_ol(_i,Cart::yy) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::yy);
      _dip[_k](_i,Cart::yzz) = PmB1*_dip[_k](_i,Cart::zz) + (_k==1)*_fak*_ol(_i,Cart::zz) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::zz);
      _dip[_k](_i,Cart::zzz) = PmB2*_dip[_k](_i,Cart::zz) + (_k==2)*_fak*_ol(_i,Cart::zz) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::zz) + term_z;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - g
  for (int _k =  0; _k < 3; _k++) {
    double term_xx = _fak*_dip[_k](0,Cart::xx);
    double term_yy = _fak*_dip[_k](0,Cart::yy);
    double term_zz = _fak*_dip[_k](0,Cart::zz);
    _dip[_k](0,Cart::xxxx) = PmB0*_dip[_k](0,Cart::xxx) + (_k==0)*_fak*_ol(0,Cart::xxx) + 3*term_xx;
    _dip[_k](0,Cart::xxxy) = PmB1*_dip[_k](0,Cart::xxx) + (_k==1)*_fak*_ol(0,Cart::xxx);
    _dip[_k](0,Cart::xxxz) = PmB2*_dip[_k](0,Cart::xxx) + (_k==2)*_fak*_ol(0,Cart::xxx);
    _dip[_k](0,Cart::xxyy) = PmB0*_dip[_k](0,Cart::xyy) + (_k==0)*_fak*_ol(0,Cart::xyy) + term_yy;
    _dip[_k](0,Cart::xxyz) = PmB1*_dip[_k](0,Cart::xxz) + (_k==1)*_fak*_ol(0,Cart::xxz);
    _dip[_k](0,Cart::xxzz) = PmB0*_dip[_k](0,Cart::xzz) + (_k==0)*_fak*_ol(0,Cart::xzz) + term_zz;
    _dip[_k](0,Cart::xyyy) = PmB0*_dip[_k](0,Cart::yyy) + (_k==0)*_fak*_ol(0,Cart::yyy);
    _dip[_k](0,Cart::xyyz) = PmB0*_dip[_k](0,Cart::yyz) + (_k==0)*_fak*_ol(0,Cart::yyz);
    _dip[_k](0,Cart::xyzz) = PmB0*_dip[_k](0,Cart::yzz) + (_k==0)*_fak*_ol(0,Cart::yzz);
    _dip[_k](0,Cart::xzzz) = PmB0*_dip[_k](0,Cart::zzz) + (_k==0)*_fak*_ol(0,Cart::zzz);
    _dip[_k](0,Cart::yyyy) = PmB1*_dip[_k](0,Cart::yyy) + (_k==1)*_fak*_ol(0,Cart::yyy) + 3*term_yy;
    _dip[_k](0,Cart::yyyz) = PmB2*_dip[_k](0,Cart::yyy) + (_k==2)*_fak*_ol(0,Cart::yyy);
    _dip[_k](0,Cart::yyzz) = PmB1*_dip[_k](0,Cart::yzz) + (_k==1)*_fak*_ol(0,Cart::yzz) + term_zz;
    _dip[_k](0,Cart::yzzz) = PmB1*_dip[_k](0,Cart::zzz) + (_k==1)*_fak*_ol(0,Cart::zzz);
    _dip[_k](0,Cart::zzzz) = PmB2*_dip[_k](0,Cart::zzz) + (_k==2)*_fak*_ol(0,Cart::zzz) + 3*term_zz;
  }
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    for (int _k =  0; _k < 3; _k++) {
      double term_xx = _fak*_dip[_k](_i,Cart::xx);
      double term_yy = _fak*_dip[_k](_i,Cart::yy);
      double term_zz = _fak*_dip[_k](_i,Cart::zz);
      _dip[_k](_i,Cart::xxxx) = PmB0*_dip[_k](_i,Cart::xxx) + (_k==0)*_fak*_ol(_i,Cart::xxx) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::xxx) + 3*term_xx;
      _dip[_k](_i,Cart::xxxy) = PmB1*_dip[_k](_i,Cart::xxx) + (_k==1)*_fak*_ol(_i,Cart::xxx) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::xxx);
      _dip[_k](_i,Cart::xxxz) = PmB2*_dip[_k](_i,Cart::xxx) + (_k==2)*_fak*_ol(_i,Cart::xxx) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::xxx);
      _dip[_k](_i,Cart::xxyy) = PmB0*_dip[_k](_i,Cart::xyy) + (_k==0)*_fak*_ol(_i,Cart::xyy) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::xyy) + term_yy;
      _dip[_k](_i,Cart::xxyz) = PmB1*_dip[_k](_i,Cart::xxz) + (_k==1)*_fak*_ol(_i,Cart::xxz) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::xxz);
      _dip[_k](_i,Cart::xxzz) = PmB0*_dip[_k](_i,Cart::xzz) + (_k==0)*_fak*_ol(_i,Cart::xzz) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::xzz) + term_zz;
      _dip[_k](_i,Cart::xyyy) = PmB0*_dip[_k](_i,Cart::yyy) + (_k==0)*_fak*_ol(_i,Cart::yyy) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::yyy);
      _dip[_k](_i,Cart::xyyz) = PmB0*_dip[_k](_i,Cart::yyz) + (_k==0)*_fak*_ol(_i,Cart::yyz) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::yyz);
      _dip[_k](_i,Cart::xyzz) = PmB0*_dip[_k](_i,Cart::yzz) + (_k==0)*_fak*_ol(_i,Cart::yzz) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::yzz);
      _dip[_k](_i,Cart::xzzz) = PmB0*_dip[_k](_i,Cart::zzz) + (_k==0)*_fak*_ol(_i,Cart::zzz) + nx[_i]*_fak*_dip[_k](i_less_x[_i],Cart::zzz);
      _dip[_k](_i,Cart::yyyy) = PmB1*_dip[_k](_i,Cart::yyy) + (_k==1)*_fak*_ol(_i,Cart::yyy) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::yyy) + 3*term_yy;
      _dip[_k](_i,Cart::yyyz) = PmB2*_dip[_k](_i,Cart::yyy) + (_k==2)*_fak*_ol(_i,Cart::yyy) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::yyy);
      _dip[_k](_i,Cart::yyzz) = PmB1*_dip[_k](_i,Cart::yzz) + (_k==1)*_fak*_ol(_i,Cart::yzz) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::yzz) + term_zz;
      _dip[_k](_i,Cart::yzzz) = PmB1*_dip[_k](_i,Cart::zzz) + (_k==1)*_fak*_ol(_i,Cart::zzz) + ny[_i]*_fak*_dip[_k](i_less_y[_i],Cart::zzz);
      _dip[_k](_i,Cart::zzzz) = PmB2*_dip[_k](_i,Cart::zzz) + (_k==2)*_fak*_ol(_i,Cart::zzz) + nz[_i]*_fak*_dip[_k](i_less_z[_i],Cart::zzz) + 3*term_zz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)


       
        
        Eigen::MatrixXd _trafo_row = getTrafo(*itr);
        Eigen::MatrixXd _trafo_col= getTrafo(*itc);       
       
        // cartesian -> spherical
       
        for ( int _i_comp = 0; _i_comp < 3; _i_comp++){

            Eigen::MatrixXd _dip_sph = _trafo_row*_dip[ _i_comp ] * _trafo_col.transpose() ;
            
            // save to _matrix
            for ( unsigned i = 0; i< _matrix[0].rows(); i++ ) {
                for (unsigned j = 0; j < _matrix[0].cols(); j++){
                    _matrix[ _i_comp ](i,j) += _dip_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
                }
            }
        }
        
       
            }// _shell_col Gaussians
        }// _shell_row Gaussians
    }
    
  
        
    
    
}}

