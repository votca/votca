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
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <votca/tools/elements.h>
#include <votca/tools/constants.h>

#include "votca/xtp/qmatom.h"







namespace votca { namespace xtp {

    

    
    void AOESP::FillBlock( Eigen::Block<Eigen::MatrixXd>& _matrix,const AOShell* _shell_row,const AOShell* _shell_col) {
        /*cout << "\nAO block: "<< endl;
        cout << "\t row: " << _shell_row->getType() << " at " << _shell_row->getPos() << endl;
        cout << "\t col: " << _shell_col->getType() << " at " << _shell_col->getPos() << endl;*/
        const double pi = boost::math::constants::pi<double>();
       
        
        
        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();
        int _lsum = _lmax_row + _lmax_col;
        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
    
        // initialize local matrix block for unnormalized cartesians
        Eigen::MatrixXd nuc   = Eigen::MatrixXd::Zero(_nrows,_ncols);
        

        //cout << nuc.size1() << ":" << nuc.size2() << endl;
        
        /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
         * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION 
         * COEFFICIENTS, AND ADD TO matrix(i,j)
         */
        
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
      
        
        // get shell positions
        const tools::vec& _pos_row = _shell_row->getPos();
        const tools::vec& _pos_col = _shell_col->getPos();
        const tools::vec  _diff    = _pos_row - _pos_col;
        // initialize some helper
      
        double _distsq = (_diff*_diff); 
        
         
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
        
        // some helpers
       


        const double zeta = _decay_row + _decay_col;

        double PmA0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
        double PmA1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
        double PmA2 = _fak2*(_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

        double PmB0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
        double PmB1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
        double PmB2 = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
        double PmC0 = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _r.getX();
        double PmC1 = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _r.getY();
        double PmC2 = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _r.getZ();
        
        
        const double _U = zeta*(PmC0*PmC0+PmC1*PmC1+PmC2*PmC2);
        
       
        const std::vector<double>_FmU=XIntegrate(_lsum+1,_U );
        //cout << endl;
        
        
        // (s-s element normiert )
        double _prefactor = 2*sqrt(1.0/pi)*pow(4.0*_decay_row*_decay_col,0.75) * _fak2 * exp(-_exparg);
        nuc(Cart::s,Cart::s)   = _prefactor * _FmU[0];
        
        typedef boost::multi_array<double, 3> ma_type;
                         ma_type nuc3(boost::extents[_nrows][_ncols][_lsum+1]);
                         typedef ma_type::index index;

                           for (index i = 0; i < _nrows; ++i) {
                               for (index j = 0; j < _ncols; ++j) {
                                   for (index k = 0; k < _lsum+1; ++k) {
                                       nuc3[i][j][k] = 0.;
                                   }
                               }
                           }


            for (int _i = 0; _i < _lsum+1; _i++){ //////////////////////
                nuc3[0][0][_i] = _prefactor*_FmU[_i];
            }
        
     
//Integrals     p - s
if (_lmax_row > 0) {
  for (int m = 0; m < _lsum; m++) {
    nuc3[Cart::x][0][m] = PmA0*nuc3[0][0][m] - PmC0*nuc3[0][0][m+1];
    nuc3[Cart::y][0][m] = PmA1*nuc3[0][0][m] - PmC1*nuc3[0][0][m+1];
    nuc3[Cart::z][0][m] = PmA2*nuc3[0][0][m] - PmC2*nuc3[0][0][m+1];
  }
}
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 1) {
  for (int m = 0; m < _lsum-1; m++) {
    double term = _fak*(nuc3[0][0][m]-nuc3[0][0][m+1]);
    nuc3[Cart::xx][0][m] = PmA0*nuc3[Cart::x][0][m] - PmC0*nuc3[Cart::x][0][m+1] + term;
    nuc3[Cart::xy][0][m] = PmA0*nuc3[Cart::y][0][m] - PmC0*nuc3[Cart::y][0][m+1];
    nuc3[Cart::xz][0][m] = PmA0*nuc3[Cart::z][0][m] - PmC0*nuc3[Cart::z][0][m+1];
    nuc3[Cart::yy][0][m] = PmA1*nuc3[Cart::y][0][m] - PmC1*nuc3[Cart::y][0][m+1] + term;
    nuc3[Cart::yz][0][m] = PmA1*nuc3[Cart::z][0][m] - PmC1*nuc3[Cart::z][0][m+1];
    nuc3[Cart::zz][0][m] = PmA2*nuc3[Cart::z][0][m] - PmC2*nuc3[Cart::z][0][m+1] + term;
  }
}
//------------------------------------------------------

//Integrals     f - s
if (_lmax_row > 2) {
  for (int m = 0; m < _lsum-2; m++) {
    nuc3[Cart::xxx][0][m] = PmA0*nuc3[Cart::xx][0][m] - PmC0*nuc3[Cart::xx][0][m+1] + 2*_fak*(nuc3[Cart::x][0][m]-nuc3[Cart::x][0][m+1]);
    nuc3[Cart::xxy][0][m] = PmA1*nuc3[Cart::xx][0][m] - PmC1*nuc3[Cart::xx][0][m+1];
    nuc3[Cart::xxz][0][m] = PmA2*nuc3[Cart::xx][0][m] - PmC2*nuc3[Cart::xx][0][m+1];
    nuc3[Cart::xyy][0][m] = PmA0*nuc3[Cart::yy][0][m] - PmC0*nuc3[Cart::yy][0][m+1];
    nuc3[Cart::xyz][0][m] = PmA0*nuc3[Cart::yz][0][m] - PmC0*nuc3[Cart::yz][0][m+1];
    nuc3[Cart::xzz][0][m] = PmA0*nuc3[Cart::zz][0][m] - PmC0*nuc3[Cart::zz][0][m+1];
    nuc3[Cart::yyy][0][m] = PmA1*nuc3[Cart::yy][0][m] - PmC1*nuc3[Cart::yy][0][m+1] + 2*_fak*(nuc3[Cart::y][0][m]-nuc3[Cart::y][0][m+1]);
    nuc3[Cart::yyz][0][m] = PmA2*nuc3[Cart::yy][0][m] - PmC2*nuc3[Cart::yy][0][m+1];
    nuc3[Cart::yzz][0][m] = PmA1*nuc3[Cart::zz][0][m] - PmC1*nuc3[Cart::zz][0][m+1];
    nuc3[Cart::zzz][0][m] = PmA2*nuc3[Cart::zz][0][m] - PmC2*nuc3[Cart::zz][0][m+1] + 2*_fak*(nuc3[Cart::z][0][m]-nuc3[Cart::z][0][m+1]);
  }
}
//------------------------------------------------------

//Integrals     g - s
if (_lmax_row > 3) {
  for (int m = 0; m < _lsum-3; m++) {
    double term_xx = _fak*(nuc3[Cart::xx][0][m]-nuc3[Cart::xx][0][m+1]);
    double term_yy = _fak*(nuc3[Cart::yy][0][m]-nuc3[Cart::yy][0][m+1]);
    double term_zz = _fak*(nuc3[Cart::zz][0][m]-nuc3[Cart::zz][0][m+1]);
    nuc3[Cart::xxxx][0][m] = PmA0*nuc3[Cart::xxx][0][m] - PmC0*nuc3[Cart::xxx][0][m+1] + 3*term_xx;
    nuc3[Cart::xxxy][0][m] = PmA1*nuc3[Cart::xxx][0][m] - PmC1*nuc3[Cart::xxx][0][m+1];
    nuc3[Cart::xxxz][0][m] = PmA2*nuc3[Cart::xxx][0][m] - PmC2*nuc3[Cart::xxx][0][m+1];
    nuc3[Cart::xxyy][0][m] = PmA0*nuc3[Cart::xyy][0][m] - PmC0*nuc3[Cart::xyy][0][m+1] + term_yy;
    nuc3[Cart::xxyz][0][m] = PmA1*nuc3[Cart::xxz][0][m] - PmC1*nuc3[Cart::xxz][0][m+1];
    nuc3[Cart::xxzz][0][m] = PmA0*nuc3[Cart::xzz][0][m] - PmC0*nuc3[Cart::xzz][0][m+1] + term_zz;
    nuc3[Cart::xyyy][0][m] = PmA0*nuc3[Cart::yyy][0][m] - PmC0*nuc3[Cart::yyy][0][m+1];
    nuc3[Cart::xyyz][0][m] = PmA0*nuc3[Cart::yyz][0][m] - PmC0*nuc3[Cart::yyz][0][m+1];
    nuc3[Cart::xyzz][0][m] = PmA0*nuc3[Cart::yzz][0][m] - PmC0*nuc3[Cart::yzz][0][m+1];
    nuc3[Cart::xzzz][0][m] = PmA0*nuc3[Cart::zzz][0][m] - PmC0*nuc3[Cart::zzz][0][m+1];
    nuc3[Cart::yyyy][0][m] = PmA1*nuc3[Cart::yyy][0][m] - PmC1*nuc3[Cart::yyy][0][m+1] + 3*term_yy;
    nuc3[Cart::yyyz][0][m] = PmA2*nuc3[Cart::yyy][0][m] - PmC2*nuc3[Cart::yyy][0][m+1];
    nuc3[Cart::yyzz][0][m] = PmA1*nuc3[Cart::yzz][0][m] - PmC1*nuc3[Cart::yzz][0][m+1] + term_zz;
    nuc3[Cart::yzzz][0][m] = PmA1*nuc3[Cart::zzz][0][m] - PmC1*nuc3[Cart::zzz][0][m+1];
    nuc3[Cart::zzzz][0][m] = PmA2*nuc3[Cart::zzz][0][m] - PmC2*nuc3[Cart::zzz][0][m+1] + 3*term_zz;
  }
}
//------------------------------------------------------




if (_lmax_col > 0) {

  //Integrals     s - p
  for (int m = 0; m < _lmax_col; m++) {
    nuc3[0][Cart::x][m] = PmB0*nuc3[0][0][m] - PmC0*nuc3[0][0][m+1];
    nuc3[0][Cart::y][m] = PmB1*nuc3[0][0][m] - PmC1*nuc3[0][0][m+1];
    nuc3[0][Cart::z][m] = PmB2*nuc3[0][0][m] - PmC2*nuc3[0][0][m+1];
  }
  //------------------------------------------------------

  //Integrals     p - p
  if (_lmax_row > 0) {
    for (int m = 0; m < _lmax_col; m++) {
      double term = _fak*(nuc3[0][0][m]-nuc3[0][0][m+1]);
      for (int _i =  1; _i < 4; _i++) {
        nuc3[_i][Cart::x][m] = PmB0*nuc3[_i][0][m] - PmC0*nuc3[_i][0][m+1] + nx[_i]*term;
        nuc3[_i][Cart::y][m] = PmB1*nuc3[_i][0][m] - PmC1*nuc3[_i][0][m+1] + ny[_i]*term;
        nuc3[_i][Cart::z][m] = PmB2*nuc3[_i][0][m] - PmC2*nuc3[_i][0][m+1] + nz[_i]*term;
      }
    }
  }
  //------------------------------------------------------

  //Integrals     d - p     f - p     g - p
  for (int m = 0; m < _lmax_col; m++) {
    for (int _i = 4; _i < n_orbitals[_lmax_row]; _i++) {
      nuc3[_i][Cart::x][m] = PmB0*nuc3[_i][0][m] - PmC0*nuc3[_i][0][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][0][m] - nuc3[i_less_x[_i]][0][m+1]);
      nuc3[_i][Cart::y][m] = PmB1*nuc3[_i][0][m] - PmC1*nuc3[_i][0][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][0][m] - nuc3[i_less_y[_i]][0][m+1]);
      nuc3[_i][Cart::z][m] = PmB2*nuc3[_i][0][m] - PmC2*nuc3[_i][0][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][0][m] - nuc3[i_less_z[_i]][0][m+1]);
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

  //Integrals     s - d
  for (int m = 0; m < _lmax_col-1; m++) {
    double term = _fak*(nuc3[0][0][m]-nuc3[0][0][m+1]);
    nuc3[0][Cart::xx][m] = PmB0*nuc3[0][Cart::x][m] - PmC0*nuc3[0][Cart::x][m+1] + term;
    nuc3[0][Cart::xy][m] = PmB0*nuc3[0][Cart::y][m] - PmC0*nuc3[0][Cart::y][m+1];
    nuc3[0][Cart::xz][m] = PmB0*nuc3[0][Cart::z][m] - PmC0*nuc3[0][Cart::z][m+1];
    nuc3[0][Cart::yy][m] = PmB1*nuc3[0][Cart::y][m] - PmC1*nuc3[0][Cart::y][m+1] + term;
    nuc3[0][Cart::yz][m] = PmB1*nuc3[0][Cart::z][m] - PmC1*nuc3[0][Cart::z][m+1];
    nuc3[0][Cart::zz][m] = PmB2*nuc3[0][Cart::z][m] - PmC2*nuc3[0][Cart::z][m+1] + term;
  }
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d
  for (int m = 0; m < _lmax_col-1; m++) {
    for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
      double term = _fak*(nuc3[_i][0][m]-nuc3[_i][0][m+1]);
      nuc3[_i][Cart::xx][m] = PmB0*nuc3[_i][Cart::x][m] - PmC0*nuc3[_i][Cart::x][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::x][m] - nuc3[i_less_x[_i]][Cart::x][m+1]) + term;
      nuc3[_i][Cart::xy][m] = PmB0*nuc3[_i][Cart::y][m] - PmC0*nuc3[_i][Cart::y][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::y][m] - nuc3[i_less_x[_i]][Cart::y][m+1]);
      nuc3[_i][Cart::xz][m] = PmB0*nuc3[_i][Cart::z][m] - PmC0*nuc3[_i][Cart::z][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::z][m] - nuc3[i_less_x[_i]][Cart::z][m+1]);
      nuc3[_i][Cart::yy][m] = PmB1*nuc3[_i][Cart::y][m] - PmC1*nuc3[_i][Cart::y][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::y][m] - nuc3[i_less_y[_i]][Cart::y][m+1]) + term;
      nuc3[_i][Cart::yz][m] = PmB1*nuc3[_i][Cart::z][m] - PmC1*nuc3[_i][Cart::z][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::z][m] - nuc3[i_less_y[_i]][Cart::z][m+1]);
      nuc3[_i][Cart::zz][m] = PmB2*nuc3[_i][Cart::z][m] - PmC2*nuc3[_i][Cart::z][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::z][m] - nuc3[i_less_z[_i]][Cart::z][m+1]) + term;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

  //Integrals     s - f
  for (int m = 0; m < _lmax_col-2; m++) {
    nuc3[0][Cart::xxx][m] = PmB0*nuc3[0][Cart::xx][m] - PmC0*nuc3[0][Cart::xx][m+1] + 2*_fak*(nuc3[0][Cart::x][m]-nuc3[0][Cart::x][m+1]);
    nuc3[0][Cart::xxy][m] = PmB1*nuc3[0][Cart::xx][m] - PmC1*nuc3[0][Cart::xx][m+1];
    nuc3[0][Cart::xxz][m] = PmB2*nuc3[0][Cart::xx][m] - PmC2*nuc3[0][Cart::xx][m+1];
    nuc3[0][Cart::xyy][m] = PmB0*nuc3[0][Cart::yy][m] - PmC0*nuc3[0][Cart::yy][m+1];
    nuc3[0][Cart::xyz][m] = PmB0*nuc3[0][Cart::yz][m] - PmC0*nuc3[0][Cart::yz][m+1];
    nuc3[0][Cart::xzz][m] = PmB0*nuc3[0][Cart::zz][m] - PmC0*nuc3[0][Cart::zz][m+1];
    nuc3[0][Cart::yyy][m] = PmB1*nuc3[0][Cart::yy][m] - PmC1*nuc3[0][Cart::yy][m+1] + 2*_fak*(nuc3[0][Cart::y][m]-nuc3[0][Cart::y][m+1]);
    nuc3[0][Cart::yyz][m] = PmB2*nuc3[0][Cart::yy][m] - PmC2*nuc3[0][Cart::yy][m+1];
    nuc3[0][Cart::yzz][m] = PmB1*nuc3[0][Cart::zz][m] - PmC1*nuc3[0][Cart::zz][m+1];
    nuc3[0][Cart::zzz][m] = PmB2*nuc3[0][Cart::zz][m] - PmC2*nuc3[0][Cart::zz][m+1] + 2*_fak*(nuc3[0][Cart::z][m]-nuc3[0][Cart::z][m+1]);
  }
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f
  for (int m = 0; m < _lmax_col-2; m++) {
    for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
      double term_x = 2*_fak*(nuc3[_i][Cart::x][m]-nuc3[_i][Cart::x][m+1]);
      double term_y = 2*_fak*(nuc3[_i][Cart::y][m]-nuc3[_i][Cart::y][m+1]);
      double term_z = 2*_fak*(nuc3[_i][Cart::z][m]-nuc3[_i][Cart::z][m+1]);
      nuc3[_i][Cart::xxx][m] = PmB0*nuc3[_i][Cart::xx][m] - PmC0*nuc3[_i][Cart::xx][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::xx][m] - nuc3[i_less_x[_i]][Cart::xx][m+1]) + term_x;
      nuc3[_i][Cart::xxy][m] = PmB1*nuc3[_i][Cart::xx][m] - PmC1*nuc3[_i][Cart::xx][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::xx][m] - nuc3[i_less_y[_i]][Cart::xx][m+1]);
      nuc3[_i][Cart::xxz][m] = PmB2*nuc3[_i][Cart::xx][m] - PmC2*nuc3[_i][Cart::xx][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::xx][m] - nuc3[i_less_z[_i]][Cart::xx][m+1]);
      nuc3[_i][Cart::xyy][m] = PmB0*nuc3[_i][Cart::yy][m] - PmC0*nuc3[_i][Cart::yy][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::yy][m] - nuc3[i_less_x[_i]][Cart::yy][m+1]);
      nuc3[_i][Cart::xyz][m] = PmB0*nuc3[_i][Cart::yz][m] - PmC0*nuc3[_i][Cart::yz][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::yz][m] - nuc3[i_less_x[_i]][Cart::yz][m+1]);
      nuc3[_i][Cart::xzz][m] = PmB0*nuc3[_i][Cart::zz][m] - PmC0*nuc3[_i][Cart::zz][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::zz][m] - nuc3[i_less_x[_i]][Cart::zz][m+1]);
      nuc3[_i][Cart::yyy][m] = PmB1*nuc3[_i][Cart::yy][m] - PmC1*nuc3[_i][Cart::yy][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::yy][m] - nuc3[i_less_y[_i]][Cart::yy][m+1]) + term_y;
      nuc3[_i][Cart::yyz][m] = PmB2*nuc3[_i][Cart::yy][m] - PmC2*nuc3[_i][Cart::yy][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::yy][m] - nuc3[i_less_z[_i]][Cart::yy][m+1]);
      nuc3[_i][Cart::yzz][m] = PmB1*nuc3[_i][Cart::zz][m] - PmC1*nuc3[_i][Cart::zz][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::zz][m] - nuc3[i_less_y[_i]][Cart::zz][m+1]);
      nuc3[_i][Cart::zzz][m] = PmB2*nuc3[_i][Cart::zz][m] - PmC2*nuc3[_i][Cart::zz][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::zz][m] - nuc3[i_less_z[_i]][Cart::zz][m+1]) + term_z;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - g
  for (int m = 0; m < _lmax_col-3; m++) {
    double term_xx = _fak*(nuc3[0][Cart::xx][m]-nuc3[0][Cart::xx][m+1]);
    double term_yy = _fak*(nuc3[0][Cart::yy][m]-nuc3[0][Cart::yy][m+1]);
    double term_zz = _fak*(nuc3[0][Cart::zz][m]-nuc3[0][Cart::zz][m+1]);
    nuc3[0][Cart::xxxx][m] = PmB0*nuc3[0][Cart::xxx][m] - PmC0*nuc3[0][Cart::xxx][m+1] + 3*term_xx;
    nuc3[0][Cart::xxxy][m] = PmB1*nuc3[0][Cart::xxx][m] - PmC1*nuc3[0][Cart::xxx][m+1];
    nuc3[0][Cart::xxxz][m] = PmB2*nuc3[0][Cart::xxx][m] - PmC2*nuc3[0][Cart::xxx][m+1];
    nuc3[0][Cart::xxyy][m] = PmB0*nuc3[0][Cart::xyy][m] - PmC0*nuc3[0][Cart::xyy][m+1] + term_yy;
    nuc3[0][Cart::xxyz][m] = PmB1*nuc3[0][Cart::xxz][m] - PmC1*nuc3[0][Cart::xxz][m+1];
    nuc3[0][Cart::xxzz][m] = PmB0*nuc3[0][Cart::xzz][m] - PmC0*nuc3[0][Cart::xzz][m+1] + term_zz;
    nuc3[0][Cart::xyyy][m] = PmB0*nuc3[0][Cart::yyy][m] - PmC0*nuc3[0][Cart::yyy][m+1];
    nuc3[0][Cart::xyyz][m] = PmB0*nuc3[0][Cart::yyz][m] - PmC0*nuc3[0][Cart::yyz][m+1];
    nuc3[0][Cart::xyzz][m] = PmB0*nuc3[0][Cart::yzz][m] - PmC0*nuc3[0][Cart::yzz][m+1];
    nuc3[0][Cart::xzzz][m] = PmB0*nuc3[0][Cart::zzz][m] - PmC0*nuc3[0][Cart::zzz][m+1];
    nuc3[0][Cart::yyyy][m] = PmB1*nuc3[0][Cart::yyy][m] - PmC1*nuc3[0][Cart::yyy][m+1] + 3*term_yy;
    nuc3[0][Cart::yyyz][m] = PmB2*nuc3[0][Cart::yyy][m] - PmC2*nuc3[0][Cart::yyy][m+1];
    nuc3[0][Cart::yyzz][m] = PmB1*nuc3[0][Cart::yzz][m] - PmC1*nuc3[0][Cart::yzz][m+1] + term_zz;
    nuc3[0][Cart::yzzz][m] = PmB1*nuc3[0][Cart::zzz][m] - PmC1*nuc3[0][Cart::zzz][m+1];
    nuc3[0][Cart::zzzz][m] = PmB2*nuc3[0][Cart::zzz][m] - PmC2*nuc3[0][Cart::zzz][m+1] + 3*term_zz;
  }
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g
  for (int m = 0; m < _lmax_col-3; m++) {
    for (int _i = 1; _i < n_orbitals[_lmax_row]; _i++) {
      double term_xx = _fak*(nuc3[_i][Cart::xx][m]-nuc3[_i][Cart::xx][m+1]);
      double term_yy = _fak*(nuc3[_i][Cart::yy][m]-nuc3[_i][Cart::yy][m+1]);
      double term_zz = _fak*(nuc3[_i][Cart::zz][m]-nuc3[_i][Cart::zz][m+1]);
      nuc3[_i][Cart::xxxx][m] = PmB0*nuc3[_i][Cart::xxx][m] - PmC0*nuc3[_i][Cart::xxx][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::xxx][m] - nuc3[i_less_x[_i]][Cart::xxx][m+1]) + 3*term_xx;
      nuc3[_i][Cart::xxxy][m] = PmB1*nuc3[_i][Cart::xxx][m] - PmC1*nuc3[_i][Cart::xxx][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::xxx][m] - nuc3[i_less_y[_i]][Cart::xxx][m+1]);
      nuc3[_i][Cart::xxxz][m] = PmB2*nuc3[_i][Cart::xxx][m] - PmC2*nuc3[_i][Cart::xxx][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::xxx][m] - nuc3[i_less_z[_i]][Cart::xxx][m+1]);
      nuc3[_i][Cart::xxyy][m] = PmB0*nuc3[_i][Cart::xyy][m] - PmC0*nuc3[_i][Cart::xyy][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::xyy][m] - nuc3[i_less_x[_i]][Cart::xyy][m+1]) + term_yy;
      nuc3[_i][Cart::xxyz][m] = PmB1*nuc3[_i][Cart::xxz][m] - PmC1*nuc3[_i][Cart::xxz][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::xxz][m] - nuc3[i_less_y[_i]][Cart::xxz][m+1]);
      nuc3[_i][Cart::xxzz][m] = PmB0*nuc3[_i][Cart::xzz][m] - PmC0*nuc3[_i][Cart::xzz][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::xzz][m] - nuc3[i_less_x[_i]][Cart::xzz][m+1]) + term_zz;
      nuc3[_i][Cart::xyyy][m] = PmB0*nuc3[_i][Cart::yyy][m] - PmC0*nuc3[_i][Cart::yyy][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::yyy][m] - nuc3[i_less_x[_i]][Cart::yyy][m+1]);
      nuc3[_i][Cart::xyyz][m] = PmB0*nuc3[_i][Cart::yyz][m] - PmC0*nuc3[_i][Cart::yyz][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::yyz][m] - nuc3[i_less_x[_i]][Cart::yyz][m+1]);
      nuc3[_i][Cart::xyzz][m] = PmB0*nuc3[_i][Cart::yzz][m] - PmC0*nuc3[_i][Cart::yzz][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::yzz][m] - nuc3[i_less_x[_i]][Cart::yzz][m+1]);
      nuc3[_i][Cart::xzzz][m] = PmB0*nuc3[_i][Cart::zzz][m] - PmC0*nuc3[_i][Cart::zzz][m+1] + nx[_i]*_fak*(nuc3[i_less_x[_i]][Cart::zzz][m] - nuc3[i_less_x[_i]][Cart::zzz][m+1]);
      nuc3[_i][Cart::yyyy][m] = PmB1*nuc3[_i][Cart::yyy][m] - PmC1*nuc3[_i][Cart::yyy][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::yyy][m] - nuc3[i_less_y[_i]][Cart::yyy][m+1]) + 3*term_yy;
      nuc3[_i][Cart::yyyz][m] = PmB2*nuc3[_i][Cart::yyy][m] - PmC2*nuc3[_i][Cart::yyy][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::yyy][m] - nuc3[i_less_z[_i]][Cart::yyy][m+1]);
      nuc3[_i][Cart::yyzz][m] = PmB1*nuc3[_i][Cart::yzz][m] - PmC1*nuc3[_i][Cart::yzz][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::yzz][m] - nuc3[i_less_y[_i]][Cart::yzz][m+1]) + term_zz;
      nuc3[_i][Cart::yzzz][m] = PmB1*nuc3[_i][Cart::zzz][m] - PmC1*nuc3[_i][Cart::zzz][m+1] + ny[_i]*_fak*(nuc3[i_less_y[_i]][Cart::zzz][m] - nuc3[i_less_y[_i]][Cart::zzz][m+1]);
      nuc3[_i][Cart::zzzz][m] = PmB2*nuc3[_i][Cart::zzz][m] - PmC2*nuc3[_i][Cart::zzz][m+1] + nz[_i]*_fak*(nuc3[i_less_z[_i]][Cart::zzz][m] - nuc3[i_less_z[_i]][Cart::zzz][m+1]) + 3*term_zz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)
                         
                         
         for (int i = 0; i < _nrows; i++) {
            for (int j = 0; j < _ncols; j++) {
               nuc(i,j) = nuc3[i][j][0];
            }
         }                

        
       
        
         Eigen::MatrixXd _nuc_sph = getTrafo(*itr)*nuc*getTrafo(*itc).transpose();
        // save to _matrix
        
        for ( unsigned i = 0; i< _matrix.rows(); i++ ) {
            for (unsigned j = 0; j < _matrix.cols(); j++) {
                _matrix(i,j) += _nuc_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
            }
        }
        
      
            }// _shell_col Gaussians
        }// _shell_row Gaussians
         return;
    }
    
  // Calculates the electrostatic potential matrix element between two basis functions, for an array of atomcores.

    void AOESP::Fillnucpotential(const AOBasis& aobasis, std::vector<QMAtom*>& _atoms) {
            tools::Elements _elements;
            _nuclearpotential = Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

            for (unsigned j = 0; j < _atoms.size(); j++) {
                double Znuc = _atoms[j]->getNuccharge();
                _aomatrix = Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
                setPosition(_atoms[j]->getPos());
                Fill(aobasis);
                _nuclearpotential -= (Znuc) * _aomatrix;        
            }
            return;
        }

        void AOESP::Fillextpotential(const AOBasis& aobasis,const std::vector<ctp::PolarSeg*> & _sites) {
            
            _externalpotential = Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());

            for (unsigned int i = 0; i < _sites.size(); i++) {
                for (ctp::PolarSeg::const_iterator it = _sites[i]->begin(); it < _sites[i]->end(); ++it) {
                    tools::vec positionofsite = (*it)->getPos() * tools::conv::nm2bohr;
                    _aomatrix = Eigen::MatrixXd::Zero(aobasis.AOBasisSize(), aobasis.AOBasisSize());
                    setPosition(positionofsite);
                    Fill(aobasis);
                    _externalpotential -= (*it)->getQ00() * _aomatrix;
                }
            }
            return;
        }    
    
}}

