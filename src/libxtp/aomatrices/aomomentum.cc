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

#include <vector>



namespace votca { namespace xtp {


    
    void AOMomentum::FillBlock( std::vector< Eigen::Block<Eigen::MatrixXd> >& _matrix, const AOShell* _shell_row,const AOShell* _shell_col) {

        
        /* Calculating the AO matrix of the gradient operator requires 
         * the raw overlap matrix (i.e. in unnormalized cartesians) 
         * with lmax of the column shell increased by one:
         * 
         *         phi_(ijk) = x^i y^j z^k exp(-beta r^2)
         * => d/dx phi_(ijk) = (i*x^(i-1) - 2beta x^(i+1)) y^j z^k exp(-beta r^2)
         *    d/dy phi_(ijk) = x^i (j*y^(j-1) - 2beta y^(j+1)) z^k exp(-beta r^2)
         *    d/dz phi_(ijk) = x^i y^j (k*z^(k-1) - 2beta z^(k+1)) exp(-beta r^2)
         * 
         * i.e.:   d/dx phi_s  = d/dx phi_(000) = -2beta phi_(100) = -2beta phi_px
         *         d/dy phi_px = d/dy phi_(100) = -2beta phi_(110) = -2beta phi_dxy
         *         d/dz phi_pz = d/dz phi_(001) = phi_(000) - 2beta phi_(002) 
         *                                      = phi_s     - 2beta phi_dxx 
         * 
         * and with that
         *         <s|d/dx|s>  = -2beta <s|px>
         *         <s|d/dy|px> = -2beta <s|dxy>
         *         <s|d/dz|pz> = <s|s> - 2beta <s|dxx>
         *         ...
         * 
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
        std::vector< Eigen::MatrixXd > _mom;
        for (int _i_comp = 0; _i_comp < 3; _i_comp++){
            _mom.push_back(Eigen::MatrixXd ::Zero(_nrows,_ncols));
        }

        std::vector< Eigen::MatrixXd > _2nd_mom; //////////////
        for (int _i_comp = 0; _i_comp < 6; _i_comp++){ //////////////
            _2nd_mom.push_back(Eigen::MatrixXd ::Zero(_nrows,_ncols)); //////////////
        } //////////////
        
        // initialize local matrix block for unnormalized cartesians of overlap
        int _nrows_ol = this->getBlockSize( _lmax_row +1 ); //////////////
        int _ncols_ol = this->getBlockSize( _lmax_col +1 );
        // make copy of shell_col and change type, lmax
        //AOShell _shell_col_local = (*_shell_col);


        Eigen::MatrixXd _ol = Eigen::MatrixXd::Zero(_nrows_ol,_ncols_ol); //////////////
        
        // get shell positions
        const tools::vec& _pos_row = _shell_row->getPos();
        const tools::vec& _pos_col = _shell_col->getPos();
        const tools::vec  _diff    = _pos_row - _pos_col;

        double _distsq = (_diff*_diff);





 int n_orbitals[] = {1, 4, 10, 20, 35, 56, 84};


 int nx[] = {
 0,
 1, 0, 0,
 2, 1, 1, 0, 0, 0,
 3, 2, 2, 1, 1, 1, 0, 0, 0, 0,
 4, 3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0,
 5, 4, 4, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0
 };

 int ny[] = {
 0,
 0, 1, 0,
 0, 1, 0, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0,
 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0
 };

 int nz[] = {
 0,
 0, 0, 1,
 0, 0, 1, 0, 1, 2,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4,
 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5
 };


 int i_less_x[] = {
  0,
  0,  0,  0,
  1,  2,  3,  0,  0,  0,
  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,
 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  0,  0,  0,  0,  0,
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,  0,  0,  0,  0,  0,  0
 };

 int i_less_y[] = {
  0,
  0,  0,  0,
  0,  1,  0,  2,  3,  0,
  0,  4,  0,  5,  6,  0,  7,  8,  9,  0,
  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,  0,
  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34,  0
 };

 int i_less_z[] = {
  0,
  0,  0,  0,
  0,  0,  1,  0,  2,  3,
  0,  0,  4,  0,  5,  6,  0,  7,  8,  9,
  0,  0, 10,  0, 11, 12,  0, 13, 14, 15,  0, 16, 17, 18, 19,
  0,  0, 20,  0, 21, 22,  0, 23, 24, 25,  0, 26, 27, 28, 29,  0, 30, 31, 32, 33, 34
 };


 int i_more_x[] = {
  1,
  4,  5,  6,
 10, 11, 12, 13, 14, 15,
 20, 21, 22, 23, 24, 25, 26, 27, 28, 29,
 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49
 };

 int i_more_y[] = {
  2,
  5,  7,  8,
 11, 13, 14, 16, 17, 18,
 21, 23, 24, 26, 27, 28, 30, 31, 32, 33,
 36, 38, 39, 41, 42, 43, 45, 46, 47, 48, 50, 51, 52, 53, 54
 };

 int i_more_z[] = {
  3,
  6,  8,  9,
 12, 14, 15, 17, 18, 19,
 22, 24, 25, 27, 28, 29, 31, 32, 33, 34,
 37, 39, 40, 42, 43, 44, 46, 47, 48, 49, 51, 52, 53, 54, 55
 };


     
       // iterate over Gaussians in this _shell_row   
        for ( AOShell::GaussianIterator itr = _shell_row->begin(); itr != _shell_row->end(); ++itr){
            // iterate over Gaussians in this _shell_col
            // get decay constant
            const double _decay_row = itr->getDecay();
            
            for ( AOShell::GaussianIterator itc = _shell_col->begin(); itc != _shell_col->end(); ++itc){
                //get decay constant
                const double _decay_col = itc->getDecay();
        
        // some helpers
        
        
        const double _fak  = 0.5/(_decay_row + _decay_col);
        const double _fak2 = 2.0 * _fak;
        double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
           
       /// check if distance between postions is big, then skip step   
       
        if ( _exparg > 30.0 ) { continue; }        
      
        const double PmA0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
        const double PmA1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
        const double PmA2 = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

        const double PmB0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
        const double PmB1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
        const double PmB2 = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();



        // calculate s-s- overlap matrix element
        _ol(0,0) = pow(4.0*_decay_row*_decay_col,0.75) * pow(_fak2,1.5)*exp(-_fak2 * _decay_row * _decay_col *_distsq); // s-s element




//Integrals     p - s
_ol(Cart::x,0) = PmA0*_ol(0,0);
_ol(Cart::y,0) = PmA1*_ol(0,0);
_ol(Cart::z,0) = PmA2*_ol(0,0);
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 0) {
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
if (_lmax_row > 1) {
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
if (_lmax_row > 2) {
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

//Integrals     h - s
if (_lmax_row > 3) {
  double term_xxx = _fak*_ol(Cart::xxx,0);
  double term_yyy = _fak*_ol(Cart::yyy,0);
  double term_zzz = _fak*_ol(Cart::zzz,0);
  _ol(Cart::xxxxx,0) = PmA0*_ol(Cart::xxxx,0) + 4*term_xxx;
  _ol(Cart::xxxxy,0) = PmA1*_ol(Cart::xxxx,0);
  _ol(Cart::xxxxz,0) = PmA2*_ol(Cart::xxxx,0);
  _ol(Cart::xxxyy,0) = PmA1*_ol(Cart::xxxy,0) + term_xxx;
  _ol(Cart::xxxyz,0) = PmA1*_ol(Cart::xxxz,0);
  _ol(Cart::xxxzz,0) = PmA2*_ol(Cart::xxxz,0) + term_xxx;
  _ol(Cart::xxyyy,0) = PmA0*_ol(Cart::xyyy,0) + term_yyy;
  _ol(Cart::xxyyz,0) = PmA2*_ol(Cart::xxyy,0);
  _ol(Cart::xxyzz,0) = PmA1*_ol(Cart::xxzz,0);
  _ol(Cart::xxzzz,0) = PmA0*_ol(Cart::xzzz,0) + term_zzz;
  _ol(Cart::xyyyy,0) = PmA0*_ol(Cart::yyyy,0);
  _ol(Cart::xyyyz,0) = PmA0*_ol(Cart::yyyz,0);
  _ol(Cart::xyyzz,0) = PmA0*_ol(Cart::yyzz,0);
  _ol(Cart::xyzzz,0) = PmA0*_ol(Cart::yzzz,0);
  _ol(Cart::xzzzz,0) = PmA0*_ol(Cart::zzzz,0);
  _ol(Cart::yyyyy,0) = PmA1*_ol(Cart::yyyy,0) + 4*term_yyy;
  _ol(Cart::yyyyz,0) = PmA2*_ol(Cart::yyyy,0);
  _ol(Cart::yyyzz,0) = PmA2*_ol(Cart::yyyz,0) + term_yyy;
  _ol(Cart::yyzzz,0) = PmA1*_ol(Cart::yzzz,0) + term_zzz;
  _ol(Cart::yzzzz,0) = PmA1*_ol(Cart::zzzz,0);
  _ol(Cart::zzzzz,0) = PmA2*_ol(Cart::zzzz,0) + 4*term_zzz;
}
//------------------------------------------------------






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


if (_lmax_col > 0) {

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

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

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
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_x = 2*_fak*_ol(_i,Cart::x);
    double term_y = 2*_fak*_ol(_i,Cart::y);
    double term_z = 2*_fak*_ol(_i,Cart::z);
    _ol(_i,Cart::xxx) = PmB0*_ol(_i,Cart::xx) + nx_i*_fak*_ol(ilx_i,Cart::xx) + term_x;
    _ol(_i,Cart::xxy) = PmB1*_ol(_i,Cart::xx) + ny_i*_fak*_ol(ily_i,Cart::xx);
    _ol(_i,Cart::xxz) = PmB2*_ol(_i,Cart::xx) + nz_i*_fak*_ol(ilz_i,Cart::xx);
    _ol(_i,Cart::xyy) = PmB0*_ol(_i,Cart::yy) + nx_i*_fak*_ol(ilx_i,Cart::yy);
    _ol(_i,Cart::xyz) = PmB0*_ol(_i,Cart::yz) + nx_i*_fak*_ol(ilx_i,Cart::yz);
    _ol(_i,Cart::xzz) = PmB0*_ol(_i,Cart::zz) + nx_i*_fak*_ol(ilx_i,Cart::zz);
    _ol(_i,Cart::yyy) = PmB1*_ol(_i,Cart::yy) + ny_i*_fak*_ol(ily_i,Cart::yy) + term_y;
    _ol(_i,Cart::yyz) = PmB2*_ol(_i,Cart::yy) + nz_i*_fak*_ol(ilz_i,Cart::yy);
    _ol(_i,Cart::yzz) = PmB1*_ol(_i,Cart::zz) + ny_i*_fak*_ol(ily_i,Cart::zz);
    _ol(_i,Cart::zzz) = PmB2*_ol(_i,Cart::zz) + nz_i*_fak*_ol(ilz_i,Cart::zz) + term_z;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

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
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_xx = _fak*_ol(_i,Cart::xx);
    double term_yy = _fak*_ol(_i,Cart::yy);
    double term_zz = _fak*_ol(_i,Cart::zz);
    _ol(_i,Cart::xxxx) = PmB0*_ol(_i,Cart::xxx) + nx_i*_fak*_ol(ilx_i,Cart::xxx) + 3*term_xx;
    _ol(_i,Cart::xxxy) = PmB1*_ol(_i,Cart::xxx) + ny_i*_fak*_ol(ily_i,Cart::xxx);
    _ol(_i,Cart::xxxz) = PmB2*_ol(_i,Cart::xxx) + nz_i*_fak*_ol(ilz_i,Cart::xxx);
    _ol(_i,Cart::xxyy) = PmB0*_ol(_i,Cart::xyy) + nx_i*_fak*_ol(ilx_i,Cart::xyy) + term_yy;
    _ol(_i,Cart::xxyz) = PmB1*_ol(_i,Cart::xxz) + ny_i*_fak*_ol(ily_i,Cart::xxz);
    _ol(_i,Cart::xxzz) = PmB0*_ol(_i,Cart::xzz) + nx_i*_fak*_ol(ilx_i,Cart::xzz) + term_zz;
    _ol(_i,Cart::xyyy) = PmB0*_ol(_i,Cart::yyy) + nx_i*_fak*_ol(ilx_i,Cart::yyy);
    _ol(_i,Cart::xyyz) = PmB0*_ol(_i,Cart::yyz) + nx_i*_fak*_ol(ilx_i,Cart::yyz);
    _ol(_i,Cart::xyzz) = PmB0*_ol(_i,Cart::yzz) + nx_i*_fak*_ol(ilx_i,Cart::yzz);
    _ol(_i,Cart::xzzz) = PmB0*_ol(_i,Cart::zzz) + nx_i*_fak*_ol(ilx_i,Cart::zzz);
    _ol(_i,Cart::yyyy) = PmB1*_ol(_i,Cart::yyy) + ny_i*_fak*_ol(ily_i,Cart::yyy) + 3*term_yy;
    _ol(_i,Cart::yyyz) = PmB2*_ol(_i,Cart::yyy) + nz_i*_fak*_ol(ilz_i,Cart::yyy);
    _ol(_i,Cart::yyzz) = PmB1*_ol(_i,Cart::yzz) + ny_i*_fak*_ol(ily_i,Cart::yzz) + term_zz;
    _ol(_i,Cart::yzzz) = PmB1*_ol(_i,Cart::zzz) + ny_i*_fak*_ol(ily_i,Cart::zzz);
    _ol(_i,Cart::zzzz) = PmB2*_ol(_i,Cart::zzz) + nz_i*_fak*_ol(ilz_i,Cart::zzz) + 3*term_zz;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - h
  double term_xxx = _fak*_ol(0,Cart::xxx);
  double term_yyy = _fak*_ol(0,Cart::yyy);
  double term_zzz = _fak*_ol(0,Cart::zzz);
  _ol(0,Cart::xxxxx) = PmB0*_ol(0,Cart::xxxx) + 4*term_xxx;
  _ol(0,Cart::xxxxy) = PmB1*_ol(0,Cart::xxxx);
  _ol(0,Cart::xxxxz) = PmB2*_ol(0,Cart::xxxx);
  _ol(0,Cart::xxxyy) = PmB1*_ol(0,Cart::xxxy) + term_xxx;
  _ol(0,Cart::xxxyz) = PmB1*_ol(0,Cart::xxxz);
  _ol(0,Cart::xxxzz) = PmB2*_ol(0,Cart::xxxz) + term_xxx;
  _ol(0,Cart::xxyyy) = PmB0*_ol(0,Cart::xyyy) + term_yyy;
  _ol(0,Cart::xxyyz) = PmB2*_ol(0,Cart::xxyy);
  _ol(0,Cart::xxyzz) = PmB1*_ol(0,Cart::xxzz);
  _ol(0,Cart::xxzzz) = PmB0*_ol(0,Cart::xzzz) + term_zzz;
  _ol(0,Cart::xyyyy) = PmB0*_ol(0,Cart::yyyy);
  _ol(0,Cart::xyyyz) = PmB0*_ol(0,Cart::yyyz);
  _ol(0,Cart::xyyzz) = PmB0*_ol(0,Cart::yyzz);
  _ol(0,Cart::xyzzz) = PmB0*_ol(0,Cart::yzzz);
  _ol(0,Cart::xzzzz) = PmB0*_ol(0,Cart::zzzz);
  _ol(0,Cart::yyyyy) = PmB1*_ol(0,Cart::yyyy) + 4*term_yyy;
  _ol(0,Cart::yyyyz) = PmB2*_ol(0,Cart::yyyy);
  _ol(0,Cart::yyyzz) = PmB2*_ol(0,Cart::yyyz) + term_yyy;
  _ol(0,Cart::yyzzz) = PmB1*_ol(0,Cart::yzzz) + term_zzz;
  _ol(0,Cart::yzzzz) = PmB1*_ol(0,Cart::zzzz);
  _ol(0,Cart::zzzzz) = PmB2*_ol(0,Cart::zzzz) + 4*term_zzz;
  //------------------------------------------------------

  //Integrals     p - h     d - h     f - h     g - h
  for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
    int nx_i = nx[_i];
    int ny_i = ny[_i];
    int nz_i = nz[_i];
    int ilx_i = i_less_x[_i];
    int ily_i = i_less_y[_i];
    int ilz_i = i_less_z[_i];
    double term_xxx = _fak*_ol(_i,Cart::xxx);
    double term_yyy = _fak*_ol(_i,Cart::yyy);
    double term_zzz = _fak*_ol(_i,Cart::zzz);
    _ol(_i,Cart::xxxxx) = PmB0*_ol(_i,Cart::xxxx) + nx_i*_fak*_ol(ilx_i,Cart::xxxx) + 4*term_xxx;
    _ol(_i,Cart::xxxxy) = PmB1*_ol(_i,Cart::xxxx) + ny_i*_fak*_ol(ily_i,Cart::xxxx);
    _ol(_i,Cart::xxxxz) = PmB2*_ol(_i,Cart::xxxx) + nz_i*_fak*_ol(ilz_i,Cart::xxxx);
    _ol(_i,Cart::xxxyy) = PmB1*_ol(_i,Cart::xxxy) + ny_i*_fak*_ol(ily_i,Cart::xxxy) + term_xxx;
    _ol(_i,Cart::xxxyz) = PmB1*_ol(_i,Cart::xxxz) + ny_i*_fak*_ol(ily_i,Cart::xxxz);
    _ol(_i,Cart::xxxzz) = PmB2*_ol(_i,Cart::xxxz) + nz_i*_fak*_ol(ilz_i,Cart::xxxz) + term_xxx;
    _ol(_i,Cart::xxyyy) = PmB0*_ol(_i,Cart::xyyy) + nx_i*_fak*_ol(ilx_i,Cart::xyyy) + term_yyy;
    _ol(_i,Cart::xxyyz) = PmB2*_ol(_i,Cart::xxyy) + nz_i*_fak*_ol(ilz_i,Cart::xxyy);
    _ol(_i,Cart::xxyzz) = PmB1*_ol(_i,Cart::xxzz) + ny_i*_fak*_ol(ily_i,Cart::xxzz);
    _ol(_i,Cart::xxzzz) = PmB0*_ol(_i,Cart::xzzz) + nx_i*_fak*_ol(ilx_i,Cart::xzzz) + term_zzz;
    _ol(_i,Cart::xyyyy) = PmB0*_ol(_i,Cart::yyyy) + nx_i*_fak*_ol(ilx_i,Cart::yyyy);
    _ol(_i,Cart::xyyyz) = PmB0*_ol(_i,Cart::yyyz) + nx_i*_fak*_ol(ilx_i,Cart::yyyz);
    _ol(_i,Cart::xyyzz) = PmB0*_ol(_i,Cart::yyzz) + nx_i*_fak*_ol(ilx_i,Cart::yyzz);
    _ol(_i,Cart::xyzzz) = PmB0*_ol(_i,Cart::yzzz) + nx_i*_fak*_ol(ilx_i,Cart::yzzz);
    _ol(_i,Cart::xzzzz) = PmB0*_ol(_i,Cart::zzzz) + nx_i*_fak*_ol(ilx_i,Cart::zzzz);
    _ol(_i,Cart::yyyyy) = PmB1*_ol(_i,Cart::yyyy) + ny_i*_fak*_ol(ily_i,Cart::yyyy) + 4*term_yyy;
    _ol(_i,Cart::yyyyz) = PmB2*_ol(_i,Cart::yyyy) + nz_i*_fak*_ol(ilz_i,Cart::yyyy);
    _ol(_i,Cart::yyyzz) = PmB2*_ol(_i,Cart::yyyz) + nz_i*_fak*_ol(ilz_i,Cart::yyyz) + term_yyy;
    _ol(_i,Cart::yyzzz) = PmB1*_ol(_i,Cart::yzzz) + ny_i*_fak*_ol(ily_i,Cart::yzzz) + term_zzz;
    _ol(_i,Cart::yzzzz) = PmB1*_ol(_i,Cart::zzzz) + ny_i*_fak*_ol(ily_i,Cart::zzzz);
    _ol(_i,Cart::zzzzz) = PmB2*_ol(_i,Cart::zzzz) + nz_i*_fak*_ol(ilz_i,Cart::zzzz) + 4*term_zzz;
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)




double _2_alpha = 2.0 * _decay_row;
double _2_beta = 2.0 * _decay_col;
for (int _i = 0; _i < _ncols; _i++) {

  int nx_i = nx[_i];
  int ny_i = ny[_i];
  int nz_i = nz[_i];
  int ilx_i = i_less_x[_i];
  int ily_i = i_less_y[_i];
  int ilz_i = i_less_z[_i];
  int imx_i = i_more_x[_i];
  int imy_i = i_more_y[_i];
  int imz_i = i_more_z[_i];

  for (int _j = 0; _j < _nrows; _j++) {

    _mom[0](_j,_i) = nx_i*_ol(_j,ilx_i) - _2_beta*_ol(_j,imx_i);
    _mom[1](_j,_i) = ny_i*_ol(_j,ily_i) - _2_beta*_ol(_j,imy_i);
    _mom[2](_j,_i) = nz_i*_ol(_j,ilz_i) - _2_beta*_ol(_j,imz_i);

    int nx_j = nx[_j];
    int ny_j = ny[_j];
    int nz_j = nz[_j];
    int ilx_j = i_less_x[_j];
    int ily_j = i_less_y[_j];
    int ilz_j = i_less_z[_j];
    int imx_j = i_more_x[_j];
    int imy_j = i_more_y[_j];
    int imz_j = i_more_z[_j];

    _2nd_mom[0](_j,_i) = nx_j * (_2_beta*_ol(ilx_j,imx_i) - nx_i*_ol(ilx_j,ilx_i)) - _2_alpha * (_2_beta*_ol(imx_j,imx_i) - nx_i*_ol(imx_j,ilx_i)); // d2/(dxdx)
    _2nd_mom[1](_j,_i) = nx_j * (_2_beta*_ol(ilx_j,imy_i) - ny_i*_ol(ilx_j,ily_i)) - _2_alpha * (_2_beta*_ol(imx_j,imy_i) - ny_i*_ol(imx_j,ily_i)); // d2/(dxdy)
    _2nd_mom[2](_j,_i) = nx_j * (_2_beta*_ol(ilx_j,imz_i) - nz_i*_ol(ilx_j,ilz_i)) - _2_alpha * (_2_beta*_ol(imx_j,imz_i) - nz_i*_ol(imx_j,ilz_i)); // d2/(dxdz)

    _2nd_mom[3](_j,_i) = ny_j * (_2_beta*_ol(ily_j,imy_i) - ny_i*_ol(ily_j,ily_i)) - _2_alpha * (_2_beta*_ol(imy_j,imy_i) - ny_i*_ol(imy_j,ily_i)); // d2/(dydy)
    _2nd_mom[4](_j,_i) = ny_j * (_2_beta*_ol(ily_j,imz_i) - nz_i*_ol(ily_j,ilz_i)) - _2_alpha * (_2_beta*_ol(imy_j,imz_i) - nz_i*_ol(imy_j,ilz_i)); // d2/(dydz)

    _2nd_mom[5](_j,_i) = nz_j * (_2_beta*_ol(ilz_j,imz_i) - nz_i*_ol(ilz_j,ilz_i)) - _2_alpha * (_2_beta*_ol(imz_j,imz_i) - nz_i*_ol(imz_j,ilz_i)); // d2/(dzdz)

  }
}

      Eigen::MatrixXd _trafo_row = getTrafo(*itr);
      Eigen::MatrixXd _trafo_col = getTrafo(*itc);
          // cartesian -> spherical
      for (int _i_comp = 0; _i_comp < 3; _i_comp++) {
        Eigen::MatrixXd _mom_sph = _trafo_row * _mom[ _i_comp ] * _trafo_col.transpose();
        // save to _matrix
        for (unsigned i = 0; i < _matrix[0].rows(); i++) {
          for (unsigned j = 0; j < _matrix[0].cols(); j++) {
            _matrix[ _i_comp ](i, j) += _mom_sph(i + _shell_row->getOffset(), j + _shell_col->getOffset());
          }
        }
      }
        
        
       
            }// _shell_col Gaussians
        }// _shell_row Gaussians
    }
    
  
        
    
    
}}

