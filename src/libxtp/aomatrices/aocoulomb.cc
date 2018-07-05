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


 

    void AOCoulomb::FillBlock(Eigen::Block<Eigen::MatrixXd>& _matrix,const  AOShell* _shell_row,const AOShell* _shell_col) {
      
            // shell info, only lmax tells how far to go
            const int _lmax_row = _shell_row->getLmax();
            const int _lmax_col = _shell_col->getLmax();

            // set size of internal block for recursion
            int _nrows = this->getBlockSize(_lmax_row);
            int _ncols = this->getBlockSize(_lmax_col);
            const int _mmax = _lmax_row + _lmax_col; 
            const int _nextra = _mmax +1;

            // get shell positions
            const tools::vec& _pos_row = _shell_row->getPos();
            const tools::vec& _pos_col = _shell_col->getPos();
            const tools::vec _diff = _pos_row - _pos_col;
            double _distsq = (_diff.getX() * _diff.getX()) + (_diff.getY() * _diff.getY()) + (_diff.getZ() * _diff.getZ());
            
            const double pi = boost::math::constants::pi<double>();
             // some helpers
            std::vector<double> _wmp=std::vector<double>(3);
            std::vector<double> _wmq=std::vector<double>(3);

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
         

           
            
          
        // iterate over Gaussians in this _shell_row
            for ( AOShell::GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
                const double _decay_row = itr->getDecay();
                const double r_decay_row = 0.5/_decay_row;
                const double powfactor_row=itr->getPowfactor();
                for ( AOShell::GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
                    
                     // get decay constants 
                        const double _decay_col = itc->getDecay();
                        const double r_decay_col = 0.5/_decay_col; 
                       const double powfactor_col=itc->getPowfactor();
                      
                         
                         tensor3d _cou(boost::extents[_nrows][_ncols][_nextra]);
                         
                         
                                  
                           for (index3d i = 0; i != _nrows; ++i) {
                               for (index3d j = 0; j != _ncols; ++j) {
                                   for (index3d k = 0; k != _nextra; ++k) {
                                       _cou[i][j][k] = 0.0;
                                   }
                               }
                           }

                         
                         
            const double _decay = _decay_row + _decay_col; 
            const double r_decay = 0.5/_decay; 
            const double r_decay_2 = 2.*r_decay; 
            const double fac_a_ac = _decay_row/_decay; 
            const double fac_c_ac = _decay_col/_decay; 
 
            
            const double wmp0 = r_decay_2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_row.getX();
            const double wmp1 = r_decay_2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_row.getY();
            const double wmp2 = r_decay_2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_row.getZ();

            const double wmq0 = r_decay_2 * (_decay_row * _pos_row.getX() + _decay_col * _pos_col.getX()) - _pos_col.getX();
            const double wmq1 = r_decay_2 * (_decay_row * _pos_row.getY() + _decay_col * _pos_col.getY()) - _pos_col.getY();
            const double wmq2 = r_decay_2 * (_decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ()) - _pos_col.getZ();

            
            const double _T = fac_a_ac * _decay_col * _distsq;


            double _fak = 2.0 * pow(pi, 2.5) / (_decay_row * _decay_col * sqrt(_decay_row + _decay_col));
            _fak = _fak *  powfactor_col*powfactor_row;

         
            const std::vector<double> _FmT=XIntegrate(_nextra, _T);

            // get initial data from _FmT -> s-s element
            for (index3d i = 0; i != _nextra; ++i) {
                _cou[0][0][i] = _fak * _FmT[i];
            }


         
//Integrals     p - s
if (_lmax_row > 0) {
  for (int m = 0; m < _mmax; m++) {
    _cou[Cart::x][0][m] = wmp0*_cou[0][0][m+1];
    _cou[Cart::y][0][m] = wmp1*_cou[0][0][m+1];
    _cou[Cart::z][0][m] = wmp2*_cou[0][0][m+1];
  }
}
//------------------------------------------------------

//Integrals     d - s
if (_lmax_row > 1) {
  for (int m = 0; m < _mmax-1; m++) {
    double term = r_decay_row*(_cou[0][0][m]-fac_c_ac*_cou[0][0][m+1]);
    _cou[Cart::xx][0][m] = wmp0*_cou[Cart::x][0][m+1] + term;
    _cou[Cart::xy][0][m] = wmp0*_cou[Cart::y][0][m+1];
    _cou[Cart::xz][0][m] = wmp0*_cou[Cart::z][0][m+1];
    _cou[Cart::yy][0][m] = wmp1*_cou[Cart::y][0][m+1] + term;
    _cou[Cart::yz][0][m] = wmp1*_cou[Cart::z][0][m+1];
    _cou[Cart::zz][0][m] = wmp2*_cou[Cart::z][0][m+1] + term;
  }
}
//------------------------------------------------------

//Integrals     f - s
if (_lmax_row > 2) {
  for (int m = 0; m < _mmax-2; m++) {
    _cou[Cart::xxx][0][m] = wmp0*_cou[Cart::xx][0][m+1] + 2*r_decay_row*(_cou[Cart::x][0][m]-fac_c_ac*_cou[Cart::x][0][m+1]);
    _cou[Cart::xxy][0][m] = wmp1*_cou[Cart::xx][0][m+1];
    _cou[Cart::xxz][0][m] = wmp2*_cou[Cart::xx][0][m+1];
    _cou[Cart::xyy][0][m] = wmp0*_cou[Cart::yy][0][m+1];
    _cou[Cart::xyz][0][m] = wmp0*_cou[Cart::yz][0][m+1];
    _cou[Cart::xzz][0][m] = wmp0*_cou[Cart::zz][0][m+1];
    _cou[Cart::yyy][0][m] = wmp1*_cou[Cart::yy][0][m+1] + 2*r_decay_row*(_cou[Cart::y][0][m]-fac_c_ac*_cou[Cart::y][0][m+1]);
    _cou[Cart::yyz][0][m] = wmp2*_cou[Cart::yy][0][m+1];
    _cou[Cart::yzz][0][m] = wmp1*_cou[Cart::zz][0][m+1];
    _cou[Cart::zzz][0][m] = wmp2*_cou[Cart::zz][0][m+1] + 2*r_decay_row*(_cou[Cart::z][0][m]-fac_c_ac*_cou[Cart::z][0][m+1]);
  }
}
//------------------------------------------------------

//Integrals     g - s
if (_lmax_row > 3) {
  for (int m = 0; m < _mmax-3; m++) {
    double term_xx = r_decay_row*(_cou[Cart::xx][0][m]-fac_c_ac*_cou[Cart::xx][0][m+1]);
    double term_yy = r_decay_row*(_cou[Cart::yy][0][m]-fac_c_ac*_cou[Cart::yy][0][m+1]);
    double term_zz = r_decay_row*(_cou[Cart::zz][0][m]-fac_c_ac*_cou[Cart::zz][0][m+1]);
    _cou[Cart::xxxx][0][m] = wmp0*_cou[Cart::xxx][0][m+1] + 3*term_xx;
    _cou[Cart::xxxy][0][m] = wmp1*_cou[Cart::xxx][0][m+1];
    _cou[Cart::xxxz][0][m] = wmp2*_cou[Cart::xxx][0][m+1];
    _cou[Cart::xxyy][0][m] = wmp0*_cou[Cart::xyy][0][m+1] + term_yy;
    _cou[Cart::xxyz][0][m] = wmp1*_cou[Cart::xxz][0][m+1];
    _cou[Cart::xxzz][0][m] = wmp0*_cou[Cart::xzz][0][m+1] + term_zz;
    _cou[Cart::xyyy][0][m] = wmp0*_cou[Cart::yyy][0][m+1];
    _cou[Cart::xyyz][0][m] = wmp0*_cou[Cart::yyz][0][m+1];
    _cou[Cart::xyzz][0][m] = wmp0*_cou[Cart::yzz][0][m+1];
    _cou[Cart::xzzz][0][m] = wmp0*_cou[Cart::zzz][0][m+1];
    _cou[Cart::yyyy][0][m] = wmp1*_cou[Cart::yyy][0][m+1] + 3*term_yy;
    _cou[Cart::yyyz][0][m] = wmp2*_cou[Cart::yyy][0][m+1];
    _cou[Cart::yyzz][0][m] = wmp1*_cou[Cart::yzz][0][m+1] + term_zz;
    _cou[Cart::yzzz][0][m] = wmp1*_cou[Cart::zzz][0][m+1];
    _cou[Cart::zzzz][0][m] = wmp2*_cou[Cart::zzz][0][m+1] + 3*term_zz;
  }
}
//------------------------------------------------------

//Integrals     h - s
if (_lmax_row > 4) {
  for (int m = 0; m < _mmax-4; m++) {
    double term_xxx = r_decay_row*(_cou[Cart::xxx][0][m]-fac_c_ac*_cou[Cart::xxx][0][m+1]);
    double term_yyy = r_decay_row*(_cou[Cart::yyy][0][m]-fac_c_ac*_cou[Cart::yyy][0][m+1]);
    double term_zzz = r_decay_row*(_cou[Cart::zzz][0][m]-fac_c_ac*_cou[Cart::zzz][0][m+1]);
    _cou[Cart::xxxxx][0][m] = wmp0*_cou[Cart::xxxx][0][m+1] + 4*term_xxx;
    _cou[Cart::xxxxy][0][m] = wmp1*_cou[Cart::xxxx][0][m+1];
    _cou[Cart::xxxxz][0][m] = wmp2*_cou[Cart::xxxx][0][m+1];
    _cou[Cart::xxxyy][0][m] = wmp1*_cou[Cart::xxxy][0][m+1] + term_xxx;
    _cou[Cart::xxxyz][0][m] = wmp1*_cou[Cart::xxxz][0][m+1];
    _cou[Cart::xxxzz][0][m] = wmp2*_cou[Cart::xxxz][0][m+1] + term_xxx;
    _cou[Cart::xxyyy][0][m] = wmp0*_cou[Cart::xyyy][0][m+1] + term_yyy;
    _cou[Cart::xxyyz][0][m] = wmp2*_cou[Cart::xxyy][0][m+1];
    _cou[Cart::xxyzz][0][m] = wmp1*_cou[Cart::xxzz][0][m+1];
    _cou[Cart::xxzzz][0][m] = wmp0*_cou[Cart::xzzz][0][m+1] + term_zzz;
    _cou[Cart::xyyyy][0][m] = wmp0*_cou[Cart::yyyy][0][m+1];
    _cou[Cart::xyyyz][0][m] = wmp0*_cou[Cart::yyyz][0][m+1];
    _cou[Cart::xyyzz][0][m] = wmp0*_cou[Cart::yyzz][0][m+1];
    _cou[Cart::xyzzz][0][m] = wmp0*_cou[Cart::yzzz][0][m+1];
    _cou[Cart::xzzzz][0][m] = wmp0*_cou[Cart::zzzz][0][m+1];
    _cou[Cart::yyyyy][0][m] = wmp1*_cou[Cart::yyyy][0][m+1] + 4*term_yyy;
    _cou[Cart::yyyyz][0][m] = wmp2*_cou[Cart::yyyy][0][m+1];
    _cou[Cart::yyyzz][0][m] = wmp2*_cou[Cart::yyyz][0][m+1] + term_yyy;
    _cou[Cart::yyzzz][0][m] = wmp1*_cou[Cart::yzzz][0][m+1] + term_zzz;
    _cou[Cart::yzzzz][0][m] = wmp1*_cou[Cart::zzzz][0][m+1];
    _cou[Cart::zzzzz][0][m] = wmp2*_cou[Cart::zzzz][0][m+1] + 4*term_zzz;
  }
}
//------------------------------------------------------

//Integrals     i - s
if (_lmax_row > 5) {
  for (int m = 0; m < _mmax-5; m++) {
    double term_xxxx = r_decay_row*(_cou[Cart::xxxx][0][m]-fac_c_ac*_cou[Cart::xxxx][0][m+1]);
    double term_xyyy = r_decay_row*(_cou[Cart::xyyy][0][m]-fac_c_ac*_cou[Cart::xyyy][0][m+1]);
    double term_xzzz = r_decay_row*(_cou[Cart::xzzz][0][m]-fac_c_ac*_cou[Cart::xzzz][0][m+1]);
    double term_yyyy = r_decay_row*(_cou[Cart::yyyy][0][m]-fac_c_ac*_cou[Cart::yyyy][0][m+1]);
    double term_yyzz = r_decay_row*(_cou[Cart::yyzz][0][m]-fac_c_ac*_cou[Cart::yyzz][0][m+1]);
    double term_yzzz = r_decay_row*(_cou[Cart::yzzz][0][m]-fac_c_ac*_cou[Cart::yzzz][0][m+1]);
    double term_zzzz = r_decay_row*(_cou[Cart::zzzz][0][m]-fac_c_ac*_cou[Cart::zzzz][0][m+1]);
    _cou[Cart::xxxxxx][0][m] = wmp0*_cou[Cart::xxxxx][0][m+1] + 5*term_xxxx;
    _cou[Cart::xxxxxy][0][m] = wmp1*_cou[Cart::xxxxx][0][m+1];
    _cou[Cart::xxxxxz][0][m] = wmp2*_cou[Cart::xxxxx][0][m+1];
    _cou[Cart::xxxxyy][0][m] = wmp1*_cou[Cart::xxxxy][0][m+1] + term_xxxx;
    _cou[Cart::xxxxyz][0][m] = wmp1*_cou[Cart::xxxxz][0][m+1];
    _cou[Cart::xxxxzz][0][m] = wmp2*_cou[Cart::xxxxz][0][m+1] + term_xxxx;
    _cou[Cart::xxxyyy][0][m] = wmp0*_cou[Cart::xxyyy][0][m+1] + 2*term_xyyy;
    _cou[Cart::xxxyyz][0][m] = wmp2*_cou[Cart::xxxyy][0][m+1];
    _cou[Cart::xxxyzz][0][m] = wmp1*_cou[Cart::xxxzz][0][m+1];
    _cou[Cart::xxxzzz][0][m] = wmp0*_cou[Cart::xxzzz][0][m+1] + 2*term_xzzz;
    _cou[Cart::xxyyyy][0][m] = wmp0*_cou[Cart::xyyyy][0][m+1] + term_yyyy;
    _cou[Cart::xxyyyz][0][m] = wmp2*_cou[Cart::xxyyy][0][m+1];
    _cou[Cart::xxyyzz][0][m] = wmp0*_cou[Cart::xyyzz][0][m+1] + term_yyzz;
    _cou[Cart::xxyzzz][0][m] = wmp1*_cou[Cart::xxzzz][0][m+1];
    _cou[Cart::xxzzzz][0][m] = wmp0*_cou[Cart::xzzzz][0][m+1] + term_zzzz;
    _cou[Cart::xyyyyy][0][m] = wmp0*_cou[Cart::yyyyy][0][m+1];
    _cou[Cart::xyyyyz][0][m] = wmp0*_cou[Cart::yyyyz][0][m+1];
    _cou[Cart::xyyyzz][0][m] = wmp0*_cou[Cart::yyyzz][0][m+1];
    _cou[Cart::xyyzzz][0][m] = wmp0*_cou[Cart::yyzzz][0][m+1];
    _cou[Cart::xyzzzz][0][m] = wmp0*_cou[Cart::yzzzz][0][m+1];
    _cou[Cart::xzzzzz][0][m] = wmp0*_cou[Cart::zzzzz][0][m+1];
    _cou[Cart::yyyyyy][0][m] = wmp1*_cou[Cart::yyyyy][0][m+1] + 5*term_yyyy;
    _cou[Cart::yyyyyz][0][m] = wmp2*_cou[Cart::yyyyy][0][m+1];
    _cou[Cart::yyyyzz][0][m] = wmp2*_cou[Cart::yyyyz][0][m+1] + term_yyyy;
    _cou[Cart::yyyzzz][0][m] = wmp1*_cou[Cart::yyzzz][0][m+1] + 2*term_yzzz;
    _cou[Cart::yyzzzz][0][m] = wmp1*_cou[Cart::yzzzz][0][m+1] + term_zzzz;
    _cou[Cart::yzzzzz][0][m] = wmp1*_cou[Cart::zzzzz][0][m+1];
    _cou[Cart::zzzzzz][0][m] = wmp2*_cou[Cart::zzzzz][0][m+1] + 5*term_zzzz;
  }
}
//------------------------------------------------------


if (_lmax_col > 0) {

  //Integrals     s - p
  for (int m = 0; m < _lmax_col; m++) {
    _cou[0][Cart::x][m] = wmq0*_cou[0][0][m+1];
    _cou[0][Cart::y][m] = wmq1*_cou[0][0][m+1];
    _cou[0][Cart::z][m] = wmq2*_cou[0][0][m+1];
  }
  //------------------------------------------------------

  //Integrals     p - p
  if (_lmax_row > 0) {
    for (int m = 0; m < _lmax_col; m++) {
      double term = r_decay*_cou[0][0][m+1];
      for (int _i =  1; _i < 4; _i++) {
        _cou[_i][Cart::x][m] = wmq0*_cou[_i][0][m+1] + nx[_i]*term;
        _cou[_i][Cart::y][m] = wmq1*_cou[_i][0][m+1] + ny[_i]*term;
        _cou[_i][Cart::z][m] = wmq2*_cou[_i][0][m+1] + nz[_i]*term;
      }
    }
  }
  //------------------------------------------------------

  //Integrals     d - p     f - p     g - p     h - p     i - p
  for (int _i_row = 2; _i_row < _lmax_row+1; _i_row++) {
    for (int m = 0; m < _lmax_col; m++) {
      for (int _i =  4; _i < n_orbitals[_lmax_row]; _i++) {
        _cou[_i][Cart::x][m] = wmq0*_cou[_i][0][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][0][m+1];
        _cou[_i][Cart::y][m] = wmq1*_cou[_i][0][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][0][m+1];
        _cou[_i][Cart::z][m] = wmq2*_cou[_i][0][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][0][m+1];
      }
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 0)


if (_lmax_col > 1) {

  //Integrals     s - d
  for (int m = 0; m < _lmax_col-1; m++) {
    double term = r_decay_col*(_cou[0][0][m]-fac_a_ac*_cou[0][0][m+1]);
    _cou[0][Cart::xx][m] = wmq0*_cou[0][Cart::x][m+1] + term;
    _cou[0][Cart::xy][m] = wmq0*_cou[0][Cart::y][m+1];
    _cou[0][Cart::xz][m] = wmq0*_cou[0][Cart::z][m+1];
    _cou[0][Cart::yy][m] = wmq1*_cou[0][Cart::y][m+1] + term;
    _cou[0][Cart::yz][m] = wmq1*_cou[0][Cart::z][m+1];
    _cou[0][Cart::zz][m] = wmq2*_cou[0][Cart::z][m+1] + term;
  }
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d     h - d     i - d
  for (int m = 0; m < _lmax_col-1; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
      double term = r_decay_col*(_cou[_i][0][m]-fac_a_ac*_cou[_i][0][m+1]);
      _cou[_i][Cart::xx][m] = wmq0*_cou[_i][Cart::x][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::x][m+1] + term;
      _cou[_i][Cart::xy][m] = wmq0*_cou[_i][Cart::y][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::y][m+1];
      _cou[_i][Cart::xz][m] = wmq0*_cou[_i][Cart::z][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::z][m+1];
      _cou[_i][Cart::yy][m] = wmq1*_cou[_i][Cart::y][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::y][m+1] + term;
      _cou[_i][Cart::yz][m] = wmq1*_cou[_i][Cart::z][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::z][m+1];
      _cou[_i][Cart::zz][m] = wmq2*_cou[_i][Cart::z][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::z][m+1] + term;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 1)


if (_lmax_col > 2) {

  //Integrals     s - f
  for (int m = 0; m < _lmax_col-2; m++) {
    _cou[0][Cart::xxx][m] = wmq0*_cou[0][Cart::xx][m+1] + 2*r_decay_col*(_cou[0][Cart::x][m]-fac_a_ac*_cou[0][Cart::x][m+1]);
    _cou[0][Cart::xxy][m] = wmq1*_cou[0][Cart::xx][m+1];
    _cou[0][Cart::xxz][m] = wmq2*_cou[0][Cart::xx][m+1];
    _cou[0][Cart::xyy][m] = wmq0*_cou[0][Cart::yy][m+1];
    _cou[0][Cart::xyz][m] = wmq0*_cou[0][Cart::yz][m+1];
    _cou[0][Cart::xzz][m] = wmq0*_cou[0][Cart::zz][m+1];
    _cou[0][Cart::yyy][m] = wmq1*_cou[0][Cart::yy][m+1] + 2*r_decay_col*(_cou[0][Cart::y][m]-fac_a_ac*_cou[0][Cart::y][m+1]);
    _cou[0][Cart::yyz][m] = wmq2*_cou[0][Cart::yy][m+1];
    _cou[0][Cart::yzz][m] = wmq1*_cou[0][Cart::zz][m+1];
    _cou[0][Cart::zzz][m] = wmq2*_cou[0][Cart::zz][m+1] + 2*r_decay_col*(_cou[0][Cart::z][m]-fac_a_ac*_cou[0][Cart::z][m+1]);
  }
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f     h - f     i - f
  for (int m = 0; m < _lmax_col-2; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
      double term_x = 2*r_decay_col*(_cou[_i][Cart::x][m]-fac_a_ac*_cou[_i][Cart::x][m+1]);
      double term_y = 2*r_decay_col*(_cou[_i][Cart::y][m]-fac_a_ac*_cou[_i][Cart::y][m+1]);
      double term_z = 2*r_decay_col*(_cou[_i][Cart::z][m]-fac_a_ac*_cou[_i][Cart::z][m+1]);
      _cou[_i][Cart::xxx][m] = wmq0*_cou[_i][Cart::xx][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xx][m+1] + term_x;
      _cou[_i][Cart::xxy][m] = wmq1*_cou[_i][Cart::xx][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xx][m+1];
      _cou[_i][Cart::xxz][m] = wmq2*_cou[_i][Cart::xx][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xx][m+1];
      _cou[_i][Cart::xyy][m] = wmq0*_cou[_i][Cart::yy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yy][m+1];
      _cou[_i][Cart::xyz][m] = wmq0*_cou[_i][Cart::yz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yz][m+1];
      _cou[_i][Cart::xzz][m] = wmq0*_cou[_i][Cart::zz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::zz][m+1];
      _cou[_i][Cart::yyy][m] = wmq1*_cou[_i][Cart::yy][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yy][m+1] + term_y;
      _cou[_i][Cart::yyz][m] = wmq2*_cou[_i][Cart::yy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::yy][m+1];
      _cou[_i][Cart::yzz][m] = wmq1*_cou[_i][Cart::zz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::zz][m+1];
      _cou[_i][Cart::zzz][m] = wmq2*_cou[_i][Cart::zz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::zz][m+1] + term_z;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 2)


if (_lmax_col > 3) {

  //Integrals     s - g
  for (int m = 0; m < _lmax_col-3; m++) {
    double term_xx = r_decay_col*(_cou[0][Cart::xx][m]-fac_a_ac*_cou[0][Cart::xx][m+1]);
    double term_yy = r_decay_col*(_cou[0][Cart::yy][m]-fac_a_ac*_cou[0][Cart::yy][m+1]);
    double term_zz = r_decay_col*(_cou[0][Cart::zz][m]-fac_a_ac*_cou[0][Cart::zz][m+1]);
    _cou[0][Cart::xxxx][m] = wmq0*_cou[0][Cart::xxx][m+1] + 3*term_xx;
    _cou[0][Cart::xxxy][m] = wmq1*_cou[0][Cart::xxx][m+1];
    _cou[0][Cart::xxxz][m] = wmq2*_cou[0][Cart::xxx][m+1];
    _cou[0][Cart::xxyy][m] = wmq0*_cou[0][Cart::xyy][m+1] + term_yy;
    _cou[0][Cart::xxyz][m] = wmq1*_cou[0][Cart::xxz][m+1];
    _cou[0][Cart::xxzz][m] = wmq0*_cou[0][Cart::xzz][m+1] + term_zz;
    _cou[0][Cart::xyyy][m] = wmq0*_cou[0][Cart::yyy][m+1];
    _cou[0][Cart::xyyz][m] = wmq0*_cou[0][Cart::yyz][m+1];
    _cou[0][Cart::xyzz][m] = wmq0*_cou[0][Cart::yzz][m+1];
    _cou[0][Cart::xzzz][m] = wmq0*_cou[0][Cart::zzz][m+1];
    _cou[0][Cart::yyyy][m] = wmq1*_cou[0][Cart::yyy][m+1] + 3*term_yy;
    _cou[0][Cart::yyyz][m] = wmq2*_cou[0][Cart::yyy][m+1];
    _cou[0][Cart::yyzz][m] = wmq1*_cou[0][Cart::yzz][m+1] + term_zz;
    _cou[0][Cart::yzzz][m] = wmq1*_cou[0][Cart::zzz][m+1];
    _cou[0][Cart::zzzz][m] = wmq2*_cou[0][Cart::zzz][m+1] + 3*term_zz;
  }
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g     h - g     i - g
  for (int m = 0; m < _lmax_col-3; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
      double term_xx = r_decay_col*(_cou[_i][Cart::xx][m]-fac_a_ac*_cou[_i][Cart::xx][m+1]);
      double term_yy = r_decay_col*(_cou[_i][Cart::yy][m]-fac_a_ac*_cou[_i][Cart::yy][m+1]);
      double term_zz = r_decay_col*(_cou[_i][Cart::zz][m]-fac_a_ac*_cou[_i][Cart::zz][m+1]);
      _cou[_i][Cart::xxxx][m] = wmq0*_cou[_i][Cart::xxx][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xxx][m+1] + 3*term_xx;
      _cou[_i][Cart::xxxy][m] = wmq1*_cou[_i][Cart::xxx][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxx][m+1];
      _cou[_i][Cart::xxxz][m] = wmq2*_cou[_i][Cart::xxx][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxx][m+1];
      _cou[_i][Cart::xxyy][m] = wmq0*_cou[_i][Cart::xyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xyy][m+1] + term_yy;
      _cou[_i][Cart::xxyz][m] = wmq1*_cou[_i][Cart::xxz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxz][m+1];
      _cou[_i][Cart::xxzz][m] = wmq0*_cou[_i][Cart::xzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xzz][m+1] + term_zz;
      _cou[_i][Cart::xyyy][m] = wmq0*_cou[_i][Cart::yyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyy][m+1];
      _cou[_i][Cart::xyyz][m] = wmq0*_cou[_i][Cart::yyz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyz][m+1];
      _cou[_i][Cart::xyzz][m] = wmq0*_cou[_i][Cart::yzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yzz][m+1];
      _cou[_i][Cart::xzzz][m] = wmq0*_cou[_i][Cart::zzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::zzz][m+1];
      _cou[_i][Cart::yyyy][m] = wmq1*_cou[_i][Cart::yyy][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yyy][m+1] + 3*term_yy;
      _cou[_i][Cart::yyyz][m] = wmq2*_cou[_i][Cart::yyy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::yyy][m+1];
      _cou[_i][Cart::yyzz][m] = wmq1*_cou[_i][Cart::yzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yzz][m+1] + term_zz;
      _cou[_i][Cart::yzzz][m] = wmq1*_cou[_i][Cart::zzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::zzz][m+1];
      _cou[_i][Cart::zzzz][m] = wmq2*_cou[_i][Cart::zzz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::zzz][m+1] + 3*term_zz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 3)


if (_lmax_col > 4) {

  //Integrals     s - h
  for (int m = 0; m < _lmax_col-4; m++) {
    double term_xxx = r_decay_col*(_cou[0][Cart::xxx][m]-fac_a_ac*_cou[0][Cart::xxx][m+1]);
    double term_yyy = r_decay_col*(_cou[0][Cart::yyy][m]-fac_a_ac*_cou[0][Cart::yyy][m+1]);
    double term_zzz = r_decay_col*(_cou[0][Cart::zzz][m]-fac_a_ac*_cou[0][Cart::zzz][m+1]);
    _cou[0][Cart::xxxxx][m] = wmq0*_cou[0][Cart::xxxx][m+1] + 4*term_xxx;
    _cou[0][Cart::xxxxy][m] = wmq1*_cou[0][Cart::xxxx][m+1];
    _cou[0][Cart::xxxxz][m] = wmq2*_cou[0][Cart::xxxx][m+1];
    _cou[0][Cart::xxxyy][m] = wmq1*_cou[0][Cart::xxxy][m+1] + term_xxx;
    _cou[0][Cart::xxxyz][m] = wmq1*_cou[0][Cart::xxxz][m+1];
    _cou[0][Cart::xxxzz][m] = wmq2*_cou[0][Cart::xxxz][m+1] + term_xxx;
    _cou[0][Cart::xxyyy][m] = wmq0*_cou[0][Cart::xyyy][m+1] + term_yyy;
    _cou[0][Cart::xxyyz][m] = wmq2*_cou[0][Cart::xxyy][m+1];
    _cou[0][Cart::xxyzz][m] = wmq1*_cou[0][Cart::xxzz][m+1];
    _cou[0][Cart::xxzzz][m] = wmq0*_cou[0][Cart::xzzz][m+1] + term_zzz;
    _cou[0][Cart::xyyyy][m] = wmq0*_cou[0][Cart::yyyy][m+1];
    _cou[0][Cart::xyyyz][m] = wmq0*_cou[0][Cart::yyyz][m+1];
    _cou[0][Cart::xyyzz][m] = wmq0*_cou[0][Cart::yyzz][m+1];
    _cou[0][Cart::xyzzz][m] = wmq0*_cou[0][Cart::yzzz][m+1];
    _cou[0][Cart::xzzzz][m] = wmq0*_cou[0][Cart::zzzz][m+1];
    _cou[0][Cart::yyyyy][m] = wmq1*_cou[0][Cart::yyyy][m+1] + 4*term_yyy;
    _cou[0][Cart::yyyyz][m] = wmq2*_cou[0][Cart::yyyy][m+1];
    _cou[0][Cart::yyyzz][m] = wmq2*_cou[0][Cart::yyyz][m+1] + term_yyy;
    _cou[0][Cart::yyzzz][m] = wmq1*_cou[0][Cart::yzzz][m+1] + term_zzz;
    _cou[0][Cart::yzzzz][m] = wmq1*_cou[0][Cart::zzzz][m+1];
    _cou[0][Cart::zzzzz][m] = wmq2*_cou[0][Cart::zzzz][m+1] + 4*term_zzz;
  }
  //------------------------------------------------------

  //Integrals     p - h     d - h     f - h     g - h     h - h     i - h
  for (int m = 0; m < _lmax_col-4; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
      double term_xxx = r_decay_col*(_cou[_i][Cart::xxx][m]-fac_a_ac*_cou[_i][Cart::xxx][m+1]);
      double term_yyy = r_decay_col*(_cou[_i][Cart::yyy][m]-fac_a_ac*_cou[_i][Cart::yyy][m+1]);
      double term_zzz = r_decay_col*(_cou[_i][Cart::zzz][m]-fac_a_ac*_cou[_i][Cart::zzz][m+1]);
      _cou[_i][Cart::xxxxx][m] = wmq0*_cou[_i][Cart::xxxx][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xxxx][m+1] + 4*term_xxx;
      _cou[_i][Cart::xxxxy][m] = wmq1*_cou[_i][Cart::xxxx][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxx][m+1];
      _cou[_i][Cart::xxxxz][m] = wmq2*_cou[_i][Cart::xxxx][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxxx][m+1];
      _cou[_i][Cart::xxxyy][m] = wmq1*_cou[_i][Cart::xxxy][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxy][m+1] + term_xxx;
      _cou[_i][Cart::xxxyz][m] = wmq1*_cou[_i][Cart::xxxz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxz][m+1];
      _cou[_i][Cart::xxxzz][m] = wmq2*_cou[_i][Cart::xxxz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxxz][m+1] + term_xxx;
      _cou[_i][Cart::xxyyy][m] = wmq0*_cou[_i][Cart::xyyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xyyy][m+1] + term_yyy;
      _cou[_i][Cart::xxyyz][m] = wmq2*_cou[_i][Cart::xxyy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxyy][m+1];
      _cou[_i][Cart::xxyzz][m] = wmq1*_cou[_i][Cart::xxzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxzz][m+1];
      _cou[_i][Cart::xxzzz][m] = wmq0*_cou[_i][Cart::xzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xzzz][m+1] + term_zzz;
      _cou[_i][Cart::xyyyy][m] = wmq0*_cou[_i][Cart::yyyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyyy][m+1];
      _cou[_i][Cart::xyyyz][m] = wmq0*_cou[_i][Cart::yyyz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyyz][m+1];
      _cou[_i][Cart::xyyzz][m] = wmq0*_cou[_i][Cart::yyzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyzz][m+1];
      _cou[_i][Cart::xyzzz][m] = wmq0*_cou[_i][Cart::yzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yzzz][m+1];
      _cou[_i][Cart::xzzzz][m] = wmq0*_cou[_i][Cart::zzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::zzzz][m+1];
      _cou[_i][Cart::yyyyy][m] = wmq1*_cou[_i][Cart::yyyy][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yyyy][m+1] + 4*term_yyy;
      _cou[_i][Cart::yyyyz][m] = wmq2*_cou[_i][Cart::yyyy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::yyyy][m+1];
      _cou[_i][Cart::yyyzz][m] = wmq2*_cou[_i][Cart::yyyz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::yyyz][m+1] + term_yyy;
      _cou[_i][Cart::yyzzz][m] = wmq1*_cou[_i][Cart::yzzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yzzz][m+1] + term_zzz;
      _cou[_i][Cart::yzzzz][m] = wmq1*_cou[_i][Cart::zzzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::zzzz][m+1];
      _cou[_i][Cart::zzzzz][m] = wmq2*_cou[_i][Cart::zzzz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::zzzz][m+1] + 4*term_zzz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 4)


if (_lmax_col > 5) {

  //Integrals     s - i
  for (int m = 0; m < _lmax_col-5; m++) {
    double term_xxxx = r_decay_col*(_cou[0][Cart::xxxx][m]-fac_a_ac*_cou[0][Cart::xxxx][m+1]);
    double term_xyyy = r_decay_col*(_cou[0][Cart::xyyy][m]-fac_a_ac*_cou[0][Cart::xyyy][m+1]);
    double term_xzzz = r_decay_col*(_cou[0][Cart::xzzz][m]-fac_a_ac*_cou[0][Cart::xzzz][m+1]);
    double term_yyyy = r_decay_col*(_cou[0][Cart::yyyy][m]-fac_a_ac*_cou[0][Cart::yyyy][m+1]);
    double term_yyzz = r_decay_col*(_cou[0][Cart::yyzz][m]-fac_a_ac*_cou[0][Cart::yyzz][m+1]);
    double term_yzzz = r_decay_col*(_cou[0][Cart::yzzz][m]-fac_a_ac*_cou[0][Cart::yzzz][m+1]);
    double term_zzzz = r_decay_col*(_cou[0][Cart::zzzz][m]-fac_a_ac*_cou[0][Cart::zzzz][m+1]);
    _cou[0][Cart::xxxxxx][m] = wmq0*_cou[0][Cart::xxxxx][m+1] + 5*term_xxxx;
    _cou[0][Cart::xxxxxy][m] = wmq1*_cou[0][Cart::xxxxx][m+1];
    _cou[0][Cart::xxxxxz][m] = wmq2*_cou[0][Cart::xxxxx][m+1];
    _cou[0][Cart::xxxxyy][m] = wmq1*_cou[0][Cart::xxxxy][m+1] + term_xxxx;
    _cou[0][Cart::xxxxyz][m] = wmq1*_cou[0][Cart::xxxxz][m+1];
    _cou[0][Cart::xxxxzz][m] = wmq2*_cou[0][Cart::xxxxz][m+1] + term_xxxx;
    _cou[0][Cart::xxxyyy][m] = wmq0*_cou[0][Cart::xxyyy][m+1] + 2*term_xyyy;
    _cou[0][Cart::xxxyyz][m] = wmq2*_cou[0][Cart::xxxyy][m+1];
    _cou[0][Cart::xxxyzz][m] = wmq1*_cou[0][Cart::xxxzz][m+1];
    _cou[0][Cart::xxxzzz][m] = wmq0*_cou[0][Cart::xxzzz][m+1] + 2*term_xzzz;
    _cou[0][Cart::xxyyyy][m] = wmq0*_cou[0][Cart::xyyyy][m+1] + term_yyyy;
    _cou[0][Cart::xxyyyz][m] = wmq2*_cou[0][Cart::xxyyy][m+1];
    _cou[0][Cart::xxyyzz][m] = wmq0*_cou[0][Cart::xyyzz][m+1] + term_yyzz;
    _cou[0][Cart::xxyzzz][m] = wmq1*_cou[0][Cart::xxzzz][m+1];
    _cou[0][Cart::xxzzzz][m] = wmq0*_cou[0][Cart::xzzzz][m+1] + term_zzzz;
    _cou[0][Cart::xyyyyy][m] = wmq0*_cou[0][Cart::yyyyy][m+1];
    _cou[0][Cart::xyyyyz][m] = wmq0*_cou[0][Cart::yyyyz][m+1];
    _cou[0][Cart::xyyyzz][m] = wmq0*_cou[0][Cart::yyyzz][m+1];
    _cou[0][Cart::xyyzzz][m] = wmq0*_cou[0][Cart::yyzzz][m+1];
    _cou[0][Cart::xyzzzz][m] = wmq0*_cou[0][Cart::yzzzz][m+1];
    _cou[0][Cart::xzzzzz][m] = wmq0*_cou[0][Cart::zzzzz][m+1];
    _cou[0][Cart::yyyyyy][m] = wmq1*_cou[0][Cart::yyyyy][m+1] + 5*term_yyyy;
    _cou[0][Cart::yyyyyz][m] = wmq2*_cou[0][Cart::yyyyy][m+1];
    _cou[0][Cart::yyyyzz][m] = wmq2*_cou[0][Cart::yyyyz][m+1] + term_yyyy;
    _cou[0][Cart::yyyzzz][m] = wmq1*_cou[0][Cart::yyzzz][m+1] + 2*term_yzzz;
    _cou[0][Cart::yyzzzz][m] = wmq1*_cou[0][Cart::yzzzz][m+1] + term_zzzz;
    _cou[0][Cart::yzzzzz][m] = wmq1*_cou[0][Cart::zzzzz][m+1];
    _cou[0][Cart::zzzzzz][m] = wmq2*_cou[0][Cart::zzzzz][m+1] + 5*term_zzzz;
  }
  //------------------------------------------------------

  //Integrals     p - i     d - i     f - i     g - i     h - i     i - i
  for (int m = 0; m < _lmax_col-5; m++) {
    for (int _i =  1; _i < n_orbitals[_lmax_row]; _i++) {
      double term_xxxx = r_decay_col*(_cou[_i][Cart::xxxx][m]-fac_a_ac*_cou[_i][Cart::xxxx][m+1]);
      double term_xyyy = r_decay_col*(_cou[_i][Cart::xyyy][m]-fac_a_ac*_cou[_i][Cart::xyyy][m+1]);
      double term_xzzz = r_decay_col*(_cou[_i][Cart::xzzz][m]-fac_a_ac*_cou[_i][Cart::xzzz][m+1]);
      double term_yyyy = r_decay_col*(_cou[_i][Cart::yyyy][m]-fac_a_ac*_cou[_i][Cart::yyyy][m+1]);
      double term_yyzz = r_decay_col*(_cou[_i][Cart::yyzz][m]-fac_a_ac*_cou[_i][Cart::yyzz][m+1]);
      double term_yzzz = r_decay_col*(_cou[_i][Cart::yzzz][m]-fac_a_ac*_cou[_i][Cart::yzzz][m+1]);
      double term_zzzz = r_decay_col*(_cou[_i][Cart::zzzz][m]-fac_a_ac*_cou[_i][Cart::zzzz][m+1]);
      _cou[_i][Cart::xxxxxx][m] = wmq0*_cou[_i][Cart::xxxxx][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xxxxx][m+1] + 5*term_xxxx;
      _cou[_i][Cart::xxxxxy][m] = wmq1*_cou[_i][Cart::xxxxx][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxxx][m+1];
      _cou[_i][Cart::xxxxxz][m] = wmq2*_cou[_i][Cart::xxxxx][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxxxx][m+1];
      _cou[_i][Cart::xxxxyy][m] = wmq1*_cou[_i][Cart::xxxxy][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxxy][m+1] + term_xxxx;
      _cou[_i][Cart::xxxxyz][m] = wmq1*_cou[_i][Cart::xxxxz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxxz][m+1];
      _cou[_i][Cart::xxxxzz][m] = wmq2*_cou[_i][Cart::xxxxz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxxxz][m+1] + term_xxxx;
      _cou[_i][Cart::xxxyyy][m] = wmq0*_cou[_i][Cart::xxyyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xxyyy][m+1] + 2*term_xyyy;
      _cou[_i][Cart::xxxyyz][m] = wmq2*_cou[_i][Cart::xxxyy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxxyy][m+1];
      _cou[_i][Cart::xxxyzz][m] = wmq1*_cou[_i][Cart::xxxzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxxzz][m+1];
      _cou[_i][Cart::xxxzzz][m] = wmq0*_cou[_i][Cart::xxzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xxzzz][m+1] + 2*term_xzzz;
      _cou[_i][Cart::xxyyyy][m] = wmq0*_cou[_i][Cart::xyyyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xyyyy][m+1] + term_yyyy;
      _cou[_i][Cart::xxyyyz][m] = wmq2*_cou[_i][Cart::xxyyy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::xxyyy][m+1];
      _cou[_i][Cart::xxyyzz][m] = wmq0*_cou[_i][Cart::xyyzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xyyzz][m+1] + term_yyzz;
      _cou[_i][Cart::xxyzzz][m] = wmq1*_cou[_i][Cart::xxzzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::xxzzz][m+1];
      _cou[_i][Cart::xxzzzz][m] = wmq0*_cou[_i][Cart::xzzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::xzzzz][m+1] + term_zzzz;
      _cou[_i][Cart::xyyyyy][m] = wmq0*_cou[_i][Cart::yyyyy][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyyyy][m+1];
      _cou[_i][Cart::xyyyyz][m] = wmq0*_cou[_i][Cart::yyyyz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyyyz][m+1];
      _cou[_i][Cart::xyyyzz][m] = wmq0*_cou[_i][Cart::yyyzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyyzz][m+1];
      _cou[_i][Cart::xyyzzz][m] = wmq0*_cou[_i][Cart::yyzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yyzzz][m+1];
      _cou[_i][Cart::xyzzzz][m] = wmq0*_cou[_i][Cart::yzzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::yzzzz][m+1];
      _cou[_i][Cart::xzzzzz][m] = wmq0*_cou[_i][Cart::zzzzz][m+1] + nx[_i]*r_decay*_cou[i_less_x[_i]][Cart::zzzzz][m+1];
      _cou[_i][Cart::yyyyyy][m] = wmq1*_cou[_i][Cart::yyyyy][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yyyyy][m+1] + 5*term_yyyy;
      _cou[_i][Cart::yyyyyz][m] = wmq2*_cou[_i][Cart::yyyyy][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::yyyyy][m+1];
      _cou[_i][Cart::yyyyzz][m] = wmq2*_cou[_i][Cart::yyyyz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::yyyyz][m+1] + term_yyyy;
      _cou[_i][Cart::yyyzzz][m] = wmq1*_cou[_i][Cart::yyzzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yyzzz][m+1] + 2*term_yzzz;
      _cou[_i][Cart::yyzzzz][m] = wmq1*_cou[_i][Cart::yzzzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::yzzzz][m+1] + term_zzzz;
      _cou[_i][Cart::yzzzzz][m] = wmq1*_cou[_i][Cart::zzzzz][m+1] + ny[_i]*r_decay*_cou[i_less_y[_i]][Cart::zzzzz][m+1];
      _cou[_i][Cart::zzzzzz][m] = wmq2*_cou[_i][Cart::zzzzz][m+1] + nz[_i]*r_decay*_cou[i_less_z[_i]][Cart::zzzzz][m+1] + 5*term_zzzz;
    }
  }
  //------------------------------------------------------

} // end if (_lmax_col > 5)
 
 
         
            // normalization and cartesian -> spherical factors
            int _ntrafo_row = _shell_row->getNumFunc() + _shell_row->getOffset();
            int _ntrafo_col = _shell_col->getNumFunc() + _shell_col->getOffset();

            

            // get transformation matrices including contraction coefficients
          const std::vector<double>& _contractions_row = itr->getContraction();
          const std::vector<double>& _contractions_col = itc->getContraction();

          

            // put _cou[i][j][0] into eigen matrix
            Eigen::MatrixXd _coumat = Eigen::MatrixXd::Zero(_nrows, _ncols);
            for (unsigned i = 0; i < _coumat.rows(); i++) {
                for (unsigned j = 0; j < _coumat.cols(); j++) {
                    _coumat(i, j) = _cou[i][j][0];
                }
            }

            Eigen::MatrixXd _cou_tmp = Eigen::MatrixXd::Zero(_ntrafo_row, _ncols);
            
            
              // s-functions
            double factor = _contractions_row[0];
            for (int _i =  0; _i < _ncols; _i++) {
              _cou_tmp(0,_i) = factor * _coumat(0,_i); /// Y 0,0
            }

            if (_lmax_row > 0) {
              // p-functions
              factor = 2.0*sqrt(_decay_row)*_contractions_row[1];
              for (int _i =  0; _i < _ncols; _i++) {
                _cou_tmp(1,_i) = factor * _coumat(3,_i); /// Y 1,0
                _cou_tmp(2,_i) = factor * _coumat(2,_i); /// Y 1,-1
                _cou_tmp(3,_i) = factor * _coumat(1,_i); /// Y 1,1
              }
            }

            if (_lmax_row > 1) {
              // d-functions
              factor = 2.0*_decay_row*_contractions_row[2];
              double factor_1 =  factor/sqrt(3.0);
              for (int _i =  0; _i < _ncols; _i++) {
                _cou_tmp(4,_i) = factor_1 * ( 2.0*_coumat(Cart::zz,_i) - _coumat(Cart::xx,_i) - _coumat(Cart::yy,_i) );  /// d3z2-r2  Y 2,0
                _cou_tmp(5,_i) = 2.*factor * _coumat(Cart::yz,_i);  /// dyz  Y 2,-1
                _cou_tmp(6,_i) = 2.*factor * _coumat(Cart::xz,_i);  /// dxz  Y 2,1
                _cou_tmp(7,_i) = 2.*factor * _coumat(Cart::xy,_i);  /// dxy  Y 2,-2
                _cou_tmp(8,_i) = factor * ( _coumat(Cart::xx,_i) - _coumat(Cart::yy,_i) );  /// dx2-y2  Y 2,2
              }
            }

            if (_lmax_row > 2) {
              // f-functions
              factor = 2.0*pow(_decay_row,1.5)*_contractions_row[3];
              double factor_1 = factor*2./sqrt(15.);
              double factor_2 = factor*sqrt(2.)/sqrt(5.);
              double factor_3 = factor*sqrt(2.)/sqrt(3.);
              for (int _i =  0; _i < _ncols; _i++) {
                _cou_tmp(9,_i) = factor_1 * ( 2.*_coumat(Cart::zzz,_i) - 3.*_coumat(Cart::xxz,_i) - 3.* _coumat(Cart::yyz,_i) ); /// Y 3,0
                _cou_tmp(10,_i) = factor_2 * ( 4.*_coumat(Cart::yzz,_i) - _coumat(Cart::xxy,_i) - _coumat(Cart::yyy,_i) ); /// Y 3,-1
                _cou_tmp(11,_i) = factor_2 * ( 4.*_coumat(Cart::xzz,_i) - _coumat(Cart::xxx,_i) - _coumat(Cart::xyy,_i) ); /// Y 3,1
                _cou_tmp(12,_i) = 4.*factor * _coumat(Cart::xyz,_i); /// Y 3,-2
                _cou_tmp(13,_i) = 2.*factor * ( _coumat(Cart::xxz,_i) - _coumat(Cart::yyz,_i) ); /// Y 3,2
                _cou_tmp(14,_i) = factor_3 * ( 3.*_coumat(Cart::xxy,_i) - _coumat(Cart::yyy,_i) ); /// Y 3,-3
                _cou_tmp(15,_i) = factor_3 * ( _coumat(Cart::xxx,_i) - 3.*_coumat(Cart::xyy,_i) ); /// Y 3,3
              }
            }

            if (_lmax_row > 3) {
              // g-functions
              factor = 2./sqrt(3.)*_decay_row*_decay_row*_contractions_row[4];
              double factor_1 = factor/sqrt(35.);
              double factor_2 = factor*4./sqrt(14.);
              double factor_3 = factor*2./sqrt(7.);
              double factor_4 = factor*2.*sqrt(2.);
              for (int _i =  0; _i < _ncols; _i++) {
                _cou_tmp(16,_i) = factor_1 * (    3.*(_coumat(Cart::xxxx,_i) + _coumat(Cart::yyyy,_i))
                                                 + 6.*_coumat(Cart::xxyy,_i)
                                               - 24.*(_coumat(Cart::xxzz,_i) + _coumat(Cart::yyzz,_i))
                                                 + 8.*_coumat(Cart::zzzz,_i) );                               /// Y 4,0
                _cou_tmp(17,_i) = factor_2 * ( -3.*(_coumat(Cart::xxyz,_i) + _coumat(Cart::yyyz,_i))
                                               + 4.*_coumat(Cart::yzzz,_i) );                                 /// Y 4,-1
                _cou_tmp(18,_i) = factor_2 * ( -3.*(_coumat(Cart::xxxz,_i) + _coumat(Cart::xyyz,_i))
                                               + 4.*_coumat(Cart::xzzz,_i) );                                 /// Y 4,1
                _cou_tmp(19,_i) = 2.*factor_3 * (    -_coumat(Cart::xxxy,_i)
                                                     - _coumat(Cart::xyyy,_i)
                                                  + 6.*_coumat(Cart::xyzz,_i) );                              /// Y 4,-2
                _cou_tmp(20,_i) = factor_3 * (      -_coumat(Cart::xxxx,_i)
                                               + 6.*(_coumat(Cart::xxzz,_i) - _coumat(Cart::yyzz,_i))
                                                  + _coumat(Cart::yyyy,_i) );                                 /// Y 4,2
                _cou_tmp(21,_i) = factor_4 * ( 3.*_coumat(Cart::xxyz,_i) 
                                                - _coumat(Cart::yyyz,_i) );                                   /// Y 4,-3
                _cou_tmp(22,_i) = factor_4 * (      _coumat(Cart::xxxz,_i) 
                                               - 3.*_coumat(Cart::xyyz,_i) );                                 /// Y 4,3
                _cou_tmp(23,_i) = 4.*factor * (   _coumat(Cart::xxxy,_i)
                                                - _coumat(Cart::xyyy,_i) );                                   /// Y 4,-4
                _cou_tmp(24,_i) = factor * (      _coumat(Cart::xxxx,_i) 
                                             - 6.*_coumat(Cart::xxyy,_i)
                                                + _coumat(Cart::yyyy,_i) );                                   /// Y 4,4
              }
            }

            if (_lmax_row > 4) {
              // h-functions
              factor = (2./3.)*pow(_decay_row,2.5)*_contractions_row[5];
              double factor_1 = factor*2./sqrt(105.);
              double factor_2 = factor*2./sqrt(7.);
              double factor_3 = factor*sqrt(6.)/3.;
              double factor_4 = factor*2.*sqrt(3.);
              double factor_5 = factor*.2*sqrt(30.);
              for (int _i =  0; _i < _ncols; _i++) {
                _cou_tmp(25,_i) = factor_1 * (   15.*(_coumat(Cart::xxxxz,_i) + _coumat(Cart::yyyyz,_i))
                                                + 30.*_coumat(Cart::xxyyz,_i)
                                               - 40.*(_coumat(Cart::xxzzz,_i) + _coumat(Cart::yyzzz,_i))
                                                 + 8.*_coumat(Cart::zzzzz,_i) );                              /// Y 5,0

                _cou_tmp(26,_i) = factor_2 * (        _coumat(Cart::xxxxy,_i)
                                                 + 2.*_coumat(Cart::xxyyy,_i)
                                               - 12.*(_coumat(Cart::xxyzz,_i) + _coumat(Cart::yyyzz,_i))
                                                    + _coumat(Cart::yyyyy,_i)
                                                 + 8.*_coumat(Cart::yzzzz,_i) );                              /// Y 5,-1

                _cou_tmp(27,_i) = factor_2 * (        _coumat(Cart::xxxxx,_i)
                                                 + 2.*_coumat(Cart::xxxyy,_i)
                                               - 12.*(_coumat(Cart::xxxzz,_i) + _coumat(Cart::xyyzz,_i))
                                                    + _coumat(Cart::xyyyy,_i)
                                                 + 8.*_coumat(Cart::xzzzz,_i) );                              /// Y 5,1

                _cou_tmp(28,_i) = 8.*factor * (     -_coumat(Cart::xxxyz,_i)
                                                   - _coumat(Cart::xyyyz,_i)
                                                + 2.*_coumat(Cart::xyzzz,_i) );                               /// Y 5,-2

                _cou_tmp(29,_i) = 4.*factor * (      -_coumat(Cart::xxxxz,_i)
                                                + 2.*(_coumat(Cart::xxzzz,_i) - _coumat(Cart::yyzzz,_i))
                                                    + _coumat(Cart::yyyyz,_i) );                              /// Y 5,2

                _cou_tmp(30,_i) = factor_3 * (   -3.*_coumat(Cart::xxxxy,_i)
                                                - 2.*_coumat(Cart::xxyyy,_i)
                                               + 24.*_coumat(Cart::xxyzz,_i)
                                                   + _coumat(Cart::yyyyy,_i)
                                                - 8.*_coumat(Cart::yyyzz,_i) );                               /// Y 5,-3

                _cou_tmp(31,_i) = factor_3 * (      -_coumat(Cart::xxxxx,_i)
                                                + 2.*_coumat(Cart::xxxyy,_i)
                                                + 8.*_coumat(Cart::xxxzz,_i)
                                                + 3.*_coumat(Cart::xyyyy,_i)
                                               - 24.*_coumat(Cart::xyyzz,_i) );                               /// Y 5,3

                _cou_tmp(32,_i) = 4.*factor_4 * (   _coumat(Cart::xxxyz,_i)
                                                  - _coumat(Cart::xyyyz,_i) );                                /// Y 5,-4

                _cou_tmp(33,_i) = factor_4 * (      _coumat(Cart::xxxxz,_i)
                                               - 6.*_coumat(Cart::xxyyz,_i)
                                                  + _coumat(Cart::yyyyz,_i) );                                /// Y 5,4

                _cou_tmp(34,_i) = factor_5 * (    5.*_coumat(Cart::xxxxy,_i)
                                               - 10.*_coumat(Cart::xxyyy,_i)
                                                   + _coumat(Cart::yyyyy,_i) );                               /// Y 5,-5

                _cou_tmp(35,_i) = factor_5 * (       _coumat(Cart::xxxxx,_i)
                                               - 10.*_coumat(Cart::xxxyy,_i)
                                                + 5.*_coumat(Cart::xyyyy,_i) );                               /// Y 5,5
              }
            }


            if (_lmax_row > 5) {
              // i-functions
              factor = (2./3.)*_decay_row*_decay_row*_decay_row*_contractions_row[6];
              double factor_1 = factor*2./sqrt(1155.);
              double factor_2 = factor*4./sqrt(55.);
              double factor_3 = factor*sqrt(22.)/11.;
              double factor_4 = factor*2.*sqrt(165.)/55.;
              double factor_5 = factor*.4*sqrt(30.);
              double factor_6 = factor*.2*sqrt(10.);
              for (int _i =  0; _i < _ncols; _i++) {
                _cou_tmp(36,_i) = factor_1 * (    -5.*(_coumat(Cart::xxxxxx,_i) + _coumat(Cart::yyyyyy,_i))
                                                - 15.*(_coumat(Cart::xxxxyy,_i) + _coumat(Cart::xxyyyy,_i))
                                                + 90.*(_coumat(Cart::xxxxzz,_i) + _coumat(Cart::yyyyzz,_i))
                                                + 180.*_coumat(Cart::xxyyzz,_i)
                                               - 120.*(_coumat(Cart::xxzzzz,_i) + _coumat(Cart::yyzzzz,_i))
                                                 + 16.*_coumat(Cart::zzzzzz,_i) );                                /// Y 6,0

                _cou_tmp(37,_i) = factor_2 * (    5.*(_coumat(Cart::xxxxyz,_i) + _coumat(Cart::yyyyyz,_i))
                                                + 10.*_coumat(Cart::xxyyyz,_i)
                                               - 20.*(_coumat(Cart::xxyzzz,_i) + _coumat(Cart::yyyzzz,_i))
                                                 + 8.*_coumat(Cart::yzzzzz,_i) );                                 /// Y 6,-1

                _cou_tmp(38,_i) = factor_2 * (    5.*(_coumat(Cart::xxxxxz,_i) + _coumat(Cart::xyyyyz,_i))
                                                + 10.*_coumat(Cart::xxxyyz,_i)
                                               - 20.*(_coumat(Cart::xxxzzz,_i) + _coumat(Cart::xyyzzz,_i))
                                                 + 8.*_coumat(Cart::xzzzzz,_i) );                                 /// Y 6,1

                _cou_tmp(39,_i) = 2.*factor_3 * (        _coumat(Cart::xxxxxy,_i)
                                                    + 2.*_coumat(Cart::xxxyyy,_i)
                                                  - 16.*(_coumat(Cart::xxxyzz,_i) + _coumat(Cart::xyyyzz,_i) - _coumat(Cart::xyzzzz,_i))
                                                       + _coumat(Cart::xyyyyy,_i) );                              /// Y 6,-2

                _cou_tmp(40,_i) = factor_3 * (        _coumat(Cart::xxxxxy,_i)
                                                    + _coumat(Cart::xxxxyy,_i)
                                               - 16.*(_coumat(Cart::xxxxzz,_i) - _coumat(Cart::xxzzzz,_i)
                                                                               - _coumat(Cart::yyyyzz,_i) + _coumat(Cart::yyzzzz,_i))
                                                    - _coumat(Cart::xxyyyy,_i)
                                                    - _coumat(Cart::yyyyyy,_i) );                                 /// Y 6,2

                _cou_tmp(41,_i) = 2.*factor_3 * (   -9.*_coumat(Cart::xxxxyz,_i)
                                                   - 6.*_coumat(Cart::xxyyyz,_i)
                                                  + 24.*_coumat(Cart::xxyzzz,_i)
                                                   + 3.*_coumat(Cart::yyyyyz,_i)
                                                   - 8.*_coumat(Cart::yyyzzz,_i) );                               /// Y 6,-3

                _cou_tmp(42,_i) = 2.*factor_3 * (   -3.*_coumat(Cart::xxxxxz,_i)
                                                   + 6.*_coumat(Cart::xxxyyz,_i)
                                                   + 8.*_coumat(Cart::xxxzzz,_i)
                                                   + 9.*_coumat(Cart::xyyyyz,_i)
                                                  - 24.*_coumat(Cart::xyyzzz,_i) );                               /// Y 6,3

                _cou_tmp(43,_i) = 4.*factor_4 * (       -_coumat(Cart::xxxxxy,_i)
                                                  + 10.*(_coumat(Cart::xxxyzz,_i) - _coumat(Cart::xyyyzz,_i))
                                                       + _coumat(Cart::xyyyyy,_i) );                              /// Y 6,-4

                _cou_tmp(44,_i) = factor_4 * (       -_coumat(Cart::xxxxxx,_i)
                                                + 5.*(_coumat(Cart::xxxxyy,_i) + _coumat(Cart::xxyyyy,_i))
                                               + 10.*(_coumat(Cart::xxxxzz,_i) + _coumat(Cart::yyyyzz,_i))
                                                - 60.*_coumat(Cart::xxyyzz,_i)
                                                   -  _coumat(Cart::yyyyyy,_i) );                                 /// Y 6,4

                _cou_tmp(45,_i) = factor_5 * (    5.*_coumat(Cart::xxxxyz,_i)
                                               - 10.*_coumat(Cart::xxyyyz,_i)
                                                   + _coumat(Cart::yyyyyz,_i) );                                  /// Y 6,-5

                _cou_tmp(46,_i) = factor_5 * (       _coumat(Cart::xxxxxz,_i)
                                               - 10.*_coumat(Cart::xxxyyz,_i)
                                                + 5.*_coumat(Cart::xyyyyz,_i) );                                  /// Y 6,5

                _cou_tmp(47,_i) = 2.*factor_6 * (    3.*_coumat(Cart::xxxxxy,_i)
                                                  - 10.*_coumat(Cart::xxxyyy,_i)
                                                   + 3.*_coumat(Cart::xyyyyy,_i) );                               /// Y 6,-6

                _cou_tmp(48,_i) = factor_6 * (        _coumat(Cart::xxxxxx,_i)
                                               - 15.*(_coumat(Cart::xxxxyy,_i) - _coumat(Cart::xxyyyy,_i))
                                                    - _coumat(Cart::yyyyyy,_i) );                                 /// Y 6,6

              }
            }

                
                
            

    Eigen::MatrixXd _cou_sph = Eigen::MatrixXd::Zero(_ntrafo_row, _ntrafo_col);  ////////////////////////////////////

        
              // s-functions
            factor = _contractions_col[0];
            for (int _i =  0; _i < _ntrafo_row; _i++) {
              _cou_sph(_i,0) = factor * _cou_tmp(_i,0); /// Y 0,0
            }

            if (_lmax_col > 0) {
              // p-functions
              factor = 2.0*sqrt(_decay_col)*_contractions_col[1];
              for (int _i =  0; _i < _ntrafo_row; _i++) {
                _cou_sph(_i,1) = factor * _cou_tmp(_i,3); /// Y 1,0
                _cou_sph(_i,2) = factor * _cou_tmp(_i,2); /// Y 1,-1
                _cou_sph(_i,3) = factor * _cou_tmp(_i,1); /// Y 1,1
              }
            }

            if (_lmax_col > 1) {
              // d-functions
              factor = 2.0*_decay_col*_contractions_col[2];
              double factor_1 =  factor/sqrt(3.0);
              for (int _i =  0; _i < _ntrafo_row; _i++) {
                _cou_sph(_i,4) = factor_1 * ( 2.0*_cou_tmp(_i,Cart::zz) - _cou_tmp(_i,Cart::xx) - _cou_tmp(_i,Cart::yy) );  /// d3z2-r2  Y 2,0
                _cou_sph(_i,5) = 2.*factor * _cou_tmp(_i,Cart::yz);  /// dyz  Y 2,-1
                _cou_sph(_i,6) = 2.*factor * _cou_tmp(_i,Cart::xz);  /// dxz  Y 2,1
                _cou_sph(_i,7) = 2.*factor * _cou_tmp(_i,Cart::xy);  /// dxy  Y 2,-2
                _cou_sph(_i,8) = factor * ( _cou_tmp(_i,Cart::xx) - _cou_tmp(_i,Cart::yy) );  /// dx2-y2  Y 2,2
              }
            }

            if (_lmax_col > 2) {
              // f-functions
              factor = 2.0*pow(_decay_col,1.5)*_contractions_col[3];
              double factor_1 = factor*2./sqrt(15.);
              double factor_2 = factor*sqrt(2.)/sqrt(5.);
              double factor_3 = factor*sqrt(2.)/sqrt(3.);
              for (int _i =  0; _i < _ntrafo_row; _i++) {
                _cou_sph(_i,9) = factor_1 * ( 2.*_cou_tmp(_i,Cart::zzz) - 3.*_cou_tmp(_i,Cart::xxz) - 3.* _cou_tmp(_i,Cart::yyz) ); /// Y 3,0
                _cou_sph(_i,10) = factor_2 * ( 4.*_cou_tmp(_i,Cart::yzz) - _cou_tmp(_i,Cart::xxy) - _cou_tmp(_i,Cart::yyy) ); /// Y 3,-1
                _cou_sph(_i,11) = factor_2 * ( 4.*_cou_tmp(_i,Cart::xzz) - _cou_tmp(_i,Cart::xxx) - _cou_tmp(_i,Cart::xyy) ); /// Y 3,1
                _cou_sph(_i,12) = 4.*factor * _cou_tmp(_i,Cart::xyz); /// Y 3,-2
                _cou_sph(_i,13) = 2.*factor * ( _cou_tmp(_i,Cart::xxz) - _cou_tmp(_i,Cart::yyz) ); /// Y 3,2
                _cou_sph(_i,14) = factor_3 * ( 3.*_cou_tmp(_i,Cart::xxy) - _cou_tmp(_i,Cart::yyy) ); /// Y 3,-3
                _cou_sph(_i,15) = factor_3 * ( _cou_tmp(_i,Cart::xxx) - 3.*_cou_tmp(_i,Cart::xyy) ); /// Y 3,3
              }
            }

            if (_lmax_col > 3) {
              // g-functions
              factor = 2./sqrt(3.)*_decay_col*_decay_col*_contractions_col[4];
              double factor_1 = factor/sqrt(35.);
              double factor_2 = factor*4./sqrt(14.);
              double factor_3 = factor*2./sqrt(7.);
              double factor_4 = factor*2.*sqrt(2.);
              for (int _i =  0; _i < _ntrafo_row; _i++) {
                _cou_sph(_i,16) = factor_1 * (    3.*(_cou_tmp(_i,Cart::xxxx) + _cou_tmp(_i,Cart::yyyy))
                                                 + 6.*_cou_tmp(_i,Cart::xxyy)
                                               - 24.*(_cou_tmp(_i,Cart::xxzz) + _cou_tmp(_i,Cart::yyzz))
                                                 + 8.*_cou_tmp(_i,Cart::zzzz) );                               /// Y 4,0
                _cou_sph(_i,17) = factor_2 * ( -3.*(_cou_tmp(_i,Cart::xxyz) + _cou_tmp(_i,Cart::yyyz))
                                               + 4.*_cou_tmp(_i,Cart::yzzz) );                                 /// Y 4,-1
                _cou_sph(_i,18) = factor_2 * ( -3.*(_cou_tmp(_i,Cart::xxxz) + _cou_tmp(_i,Cart::xyyz))
                                               + 4.*_cou_tmp(_i,Cart::xzzz) );                                 /// Y 4,1
                _cou_sph(_i,19) = 2.*factor_3 * (    -_cou_tmp(_i,Cart::xxxy)
                                                     - _cou_tmp(_i,Cart::xyyy)
                                                  + 6.*_cou_tmp(_i,Cart::xyzz) );                              /// Y 4,-2
                _cou_sph(_i,20) = factor_3 * (      -_cou_tmp(_i,Cart::xxxx)
                                               + 6.*(_cou_tmp(_i,Cart::xxzz) - _cou_tmp(_i,Cart::yyzz))
                                                  + _cou_tmp(_i,Cart::yyyy) );                                 /// Y 4,2
                _cou_sph(_i,21) = factor_4 * ( 3.*_cou_tmp(_i,Cart::xxyz) 
                                                - _cou_tmp(_i,Cart::yyyz) );                                   /// Y 4,-3
                _cou_sph(_i,22) = factor_4 * (      _cou_tmp(_i,Cart::xxxz) 
                                               - 3.*_cou_tmp(_i,Cart::xyyz) );                                 /// Y 4,3
                _cou_sph(_i,23) = 4.*factor * (   _cou_tmp(_i,Cart::xxxy)
                                                - _cou_tmp(_i,Cart::xyyy) );                                   /// Y 4,-4
                _cou_sph(_i,24) = factor * (      _cou_tmp(_i,Cart::xxxx) 
                                             - 6.*_cou_tmp(_i,Cart::xxyy)
                                                + _cou_tmp(_i,Cart::yyyy) );                                   /// Y 4,4
              }
            }

            if (_lmax_col > 4) {
              // h-functions
              factor = (2./3.)*pow(_decay_col,2.5)*_contractions_col[5];
              double factor_1 = factor*2./sqrt(105.);
              double factor_2 = factor*2./sqrt(7.);
              double factor_3 = factor*sqrt(6.)/3.;
              double factor_4 = factor*2.*sqrt(3.);
              double factor_5 = factor*.2*sqrt(30.);
              for (int _i =  0; _i < _ntrafo_row; _i++) {
                _cou_sph(_i,25) = factor_1 * (   15.*(_cou_tmp(_i,Cart::xxxxz) + _cou_tmp(_i,Cart::yyyyz))
                                                + 30.*_cou_tmp(_i,Cart::xxyyz)
                                               - 40.*(_cou_tmp(_i,Cart::xxzzz) + _cou_tmp(_i,Cart::yyzzz))
                                                 + 8.*_cou_tmp(_i,Cart::zzzzz) );                              /// Y 5,0

                _cou_sph(_i,26) = factor_2 * (        _cou_tmp(_i,Cart::xxxxy)
                                                 + 2.*_cou_tmp(_i,Cart::xxyyy)
                                               - 12.*(_cou_tmp(_i,Cart::xxyzz) + _cou_tmp(_i,Cart::yyyzz))
                                                    + _cou_tmp(_i,Cart::yyyyy)
                                                 + 8.*_cou_tmp(_i,Cart::yzzzz) );                              /// Y 5,-1

                _cou_sph(_i,27) = factor_2 * (        _cou_tmp(_i,Cart::xxxxx)
                                                 + 2.*_cou_tmp(_i,Cart::xxxyy)
                                               - 12.*(_cou_tmp(_i,Cart::xxxzz) + _cou_tmp(_i,Cart::xyyzz))
                                                    + _cou_tmp(_i,Cart::xyyyy)
                                                 + 8.*_cou_tmp(_i,Cart::xzzzz) );                              /// Y 5,1

                _cou_sph(_i,28) = 8.*factor * (     -_cou_tmp(_i,Cart::xxxyz)
                                                   - _cou_tmp(_i,Cart::xyyyz)
                                                + 2.*_cou_tmp(_i,Cart::xyzzz) );                               /// Y 5,-2

                _cou_sph(_i,29) = 4.*factor * (      -_cou_tmp(_i,Cart::xxxxz)
                                                + 2.*(_cou_tmp(_i,Cart::xxzzz) - _cou_tmp(_i,Cart::yyzzz))
                                                    + _cou_tmp(_i,Cart::yyyyz) );                              /// Y 5,2

                _cou_sph(_i,30) = factor_3 * (   -3.*_cou_tmp(_i,Cart::xxxxy)
                                                - 2.*_cou_tmp(_i,Cart::xxyyy)
                                               + 24.*_cou_tmp(_i,Cart::xxyzz)
                                                   + _cou_tmp(_i,Cart::yyyyy)
                                                - 8.*_cou_tmp(_i,Cart::yyyzz) );                               /// Y 5,-3

                _cou_sph(_i,31) = factor_3 * (      -_cou_tmp(_i,Cart::xxxxx)
                                                + 2.*_cou_tmp(_i,Cart::xxxyy)
                                                + 8.*_cou_tmp(_i,Cart::xxxzz)
                                                + 3.*_cou_tmp(_i,Cart::xyyyy)
                                               - 24.*_cou_tmp(_i,Cart::xyyzz) );                               /// Y 5,3

                _cou_sph(_i,32) = 4.*factor_4 * (   _cou_tmp(_i,Cart::xxxyz)
                                                  - _cou_tmp(_i,Cart::xyyyz) );                                /// Y 5,-4

                _cou_sph(_i,33) = factor_4 * (      _cou_tmp(_i,Cart::xxxxz)
                                               - 6.*_cou_tmp(_i,Cart::xxyyz)
                                                  + _cou_tmp(_i,Cart::yyyyz) );                                /// Y 5,4

                _cou_sph(_i,34) = factor_5 * (    5.*_cou_tmp(_i,Cart::xxxxy)
                                               - 10.*_cou_tmp(_i,Cart::xxyyy)
                                                   + _cou_tmp(_i,Cart::yyyyy) );                               /// Y 5,-5

                _cou_sph(_i,35) = factor_5 * (       _cou_tmp(_i,Cart::xxxxx)
                                               - 10.*_cou_tmp(_i,Cart::xxxyy)
                                                + 5.*_cou_tmp(_i,Cart::xyyyy) );                               /// Y 5,5
              }
            }


            if (_lmax_col > 5) {
              // i-functions
              factor = (2./3.)*_decay_col*_decay_col*_decay_col*_contractions_col[6];
              double factor_1 = factor*2./sqrt(1155.);
              double factor_2 = factor*4./sqrt(55.);
              double factor_3 = factor*sqrt(22.)/11.;
              double factor_4 = factor*2.*sqrt(165.)/55.;
              double factor_5 = factor*.4*sqrt(30.);
              double factor_6 = factor*.2*sqrt(10.);
              for (int _i =  0; _i < _ntrafo_row; _i++) {
                _cou_sph(_i,36) = factor_1 * (    -5.*(_cou_tmp(_i,Cart::xxxxxx) + _cou_tmp(_i,Cart::yyyyyy))
                                                - 15.*(_cou_tmp(_i,Cart::xxxxyy) + _cou_tmp(_i,Cart::xxyyyy))
                                                + 90.*(_cou_tmp(_i,Cart::xxxxzz) + _cou_tmp(_i,Cart::yyyyzz))
                                                + 180.*_cou_tmp(_i,Cart::xxyyzz)
                                               - 120.*(_cou_tmp(_i,Cart::xxzzzz) + _cou_tmp(_i,Cart::yyzzzz))
                                                 + 16.*_cou_tmp(_i,Cart::zzzzzz) );                                /// Y 6,0

                _cou_sph(_i,37) = factor_2 * (    5.*(_cou_tmp(_i,Cart::xxxxyz) + _cou_tmp(_i,Cart::yyyyyz))
                                                + 10.*_cou_tmp(_i,Cart::xxyyyz)
                                               - 20.*(_cou_tmp(_i,Cart::xxyzzz) + _cou_tmp(_i,Cart::yyyzzz))
                                                 + 8.*_cou_tmp(_i,Cart::yzzzzz) );                                 /// Y 6,-1

                _cou_sph(_i,38) = factor_2 * (    5.*(_cou_tmp(_i,Cart::xxxxxz) + _cou_tmp(_i,Cart::xyyyyz))
                                                + 10.*_cou_tmp(_i,Cart::xxxyyz)
                                               - 20.*(_cou_tmp(_i,Cart::xxxzzz) + _cou_tmp(_i,Cart::xyyzzz))
                                                 + 8.*_cou_tmp(_i,Cart::xzzzzz) );                                 /// Y 6,1

                _cou_sph(_i,39) = 2.*factor_3 * (        _cou_tmp(_i,Cart::xxxxxy)
                                                    + 2.*_cou_tmp(_i,Cart::xxxyyy)
                                                  - 16.*(_cou_tmp(_i,Cart::xxxyzz) + _cou_tmp(_i,Cart::xyyyzz) - _cou_tmp(_i,Cart::xyzzzz))
                                                       + _cou_tmp(_i,Cart::xyyyyy) );                              /// Y 6,-2

                _cou_sph(_i,40) = factor_3 * (        _cou_tmp(_i,Cart::xxxxxy)
                                                    + _cou_tmp(_i,Cart::xxxxyy)
                                               - 16.*(_cou_tmp(_i,Cart::xxxxzz) - _cou_tmp(_i,Cart::xxzzzz)
                                                                                - _cou_tmp(_i,Cart::yyyyzz) + _cou_tmp(_i,Cart::yyzzzz))
                                                    - _cou_tmp(_i,Cart::xxyyyy)
                                                    - _cou_tmp(_i,Cart::yyyyyy) );                                 /// Y 6,2

                _cou_sph(_i,41) = 2.*factor_3 * (   -9.*_cou_tmp(_i,Cart::xxxxyz)
                                                   - 6.*_cou_tmp(_i,Cart::xxyyyz)
                                                  + 24.*_cou_tmp(_i,Cart::xxyzzz)
                                                   + 3.*_cou_tmp(_i,Cart::yyyyyz)
                                                   - 8.*_cou_tmp(_i,Cart::yyyzzz) );                               /// Y 6,-3

                _cou_sph(_i,42) = 2.*factor_3 * (   -3.*_cou_tmp(_i,Cart::xxxxxz)
                                                   + 6.*_cou_tmp(_i,Cart::xxxyyz)
                                                   + 8.*_cou_tmp(_i,Cart::xxxzzz)
                                                   + 9.*_cou_tmp(_i,Cart::xyyyyz)
                                                  - 24.*_cou_tmp(_i,Cart::xyyzzz) );                               /// Y 6,3

                _cou_sph(_i,43) = 4.*factor_4 * (       -_cou_tmp(_i,Cart::xxxxxy)
                                                  + 10.*(_cou_tmp(_i,Cart::xxxyzz) - _cou_tmp(_i,Cart::xyyyzz))
                                                       + _cou_tmp(_i,Cart::xyyyyy) );                              /// Y 6,-4

                _cou_sph(_i,44) = factor_4 * (       -_cou_tmp(_i,Cart::xxxxxx)
                                                + 5.*(_cou_tmp(_i,Cart::xxxxyy) + _cou_tmp(_i,Cart::xxyyyy))
                                               + 10.*(_cou_tmp(_i,Cart::xxxxzz) + _cou_tmp(_i,Cart::yyyyzz))
                                                - 60.*_cou_tmp(_i,Cart::xxyyzz)
                                                   -  _cou_tmp(_i,Cart::yyyyyy) );                                 /// Y 6,4

                _cou_sph(_i,45) = factor_5 * (    5.*_cou_tmp(_i,Cart::xxxxyz)
                                               - 10.*_cou_tmp(_i,Cart::xxyyyz)
                                                   + _cou_tmp(_i,Cart::yyyyyz) );                                  /// Y 6,-5

                _cou_sph(_i,46) = factor_5 * (       _cou_tmp(_i,Cart::xxxxxz)
                                               - 10.*_cou_tmp(_i,Cart::xxxyyz)
                                                + 5.*_cou_tmp(_i,Cart::xyyyyz) );                                  /// Y 6,5

                _cou_sph(_i,47) = 2.*factor_6 * (    3.*_cou_tmp(_i,Cart::xxxxxy)
                                                  - 10.*_cou_tmp(_i,Cart::xxxyyy)
                                                   + 3.*_cou_tmp(_i,Cart::xyyyyy) );                               /// Y 6,-6

                _cou_sph(_i,48) = factor_6 * (        _cou_tmp(_i,Cart::xxxxxx)
                                               - 15.*(_cou_tmp(_i,Cart::xxxxyy) - _cou_tmp(_i,Cart::xxyyyy))
                                                    - _cou_tmp(_i,Cart::yyyyyy) );                                 /// Y 6,6

              }
            }



            // save to _matrix
            for (unsigned i = 0; i < _matrix.rows(); i++) {
                for (unsigned j = 0; j < _matrix.cols(); j++) {
                    _matrix(i, j) += _cou_sph(i + _shell_row->getOffset(), j + _shell_col->getOffset());
                }
            }

                } // _shell_col Gaussians
            } // _shell_row Gaussians
           return; 
            }    
    

    
    Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt_GWBSE(const AOOverlap& _auxoverlap, double etol){
        
        
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eo(_auxoverlap.Matrix());
    Eigen::MatrixXd Ssqrt=eo.operatorSqrt();
    //This converts V into (S1/2 V S1/2)-1/2 S1/2, which is needed to construct 4c integrals,
       
      Eigen::MatrixXd ortho=Ssqrt*_aomatrix*Ssqrt;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ortho); 
      Eigen::VectorXd diagonal=Eigen::VectorXd::Zero(es.eigenvalues().size());
      removedfunctions=0;
      for (unsigned i=0;i<diagonal.size();++i){
          if(es.eigenvalues()(i)<etol){
              removedfunctions++;
          }else{
              diagonal(i)=1.0/std::sqrt(es.eigenvalues()(i));
          }
      }
      
      Eigen::MatrixXd Vm1=es.eigenvectors() * diagonal.asDiagonal() * es.eigenvectors().transpose();
      Eigen::MatrixXd result=(Vm1*Ssqrt).transpose();
    return result;
    }
    
     Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt(double etol){
       Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(_aomatrix);
       Eigen::VectorXd diagonal=Eigen::VectorXd::Zero(es.eigenvalues().size());
      removedfunctions=0;
      for (unsigned i=0;i<diagonal.size();++i){
          if(es.eigenvalues()(i)<etol){
              removedfunctions++;
          }else{
              diagonal(i)=1.0/std::sqrt(es.eigenvalues()(i));
          }
      }
           
     return es.eigenvectors() * diagonal.asDiagonal() * es.eigenvectors().transpose();
    }
    
    
    

    
}}

