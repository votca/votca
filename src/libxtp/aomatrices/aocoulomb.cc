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
 * olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */


#include <votca/xtp/aomatrix.h>
#include <votca/xtp/aobasis.h>
#include <vector>



namespace votca { namespace xtp {


 

    void AOCoulomb::FillBlock(Eigen::Block<Eigen::MatrixXd>& matrix,const  AOShell* shell_row,const AOShell* shell_col) {
      
            // shell info, only lmax tells how far to go
            const int lmax_row = shell_row->getLmax();
            const int lmax_col = shell_col->getLmax();

            // set size of internal block for recursion
            int nrows = this->getBlockSize(lmax_row);
            int ncols = this->getBlockSize(lmax_col);
            const int mmax = lmax_row + lmax_col; 
            const int nextra = mmax +1;

            // get shell positions
            const tools::vec& pos_row = shell_row->getPos();
            const tools::vec& pos_col = shell_col->getPos();
            const tools::vec diff = pos_row - pos_col;
            double distsq = (diff.getX() * diff.getX()) + (diff.getY() * diff.getY()) + (diff.getZ() * diff.getZ());
            
            const double pi = boost::math::constants::pi<double>();
             // some helpers
            std::vector<double> wmp=std::vector<double>(3);
            std::vector<double> wmq=std::vector<double>(3);

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
         

           
            
          
        // iterate over Gaussians in this shell_row
            for ( AOShell::GaussianIterator itr = shell_row->begin(); itr != shell_row->end(); ++itr){
            // iterate over Gaussians in this shell_col
                const double decay_row = itr->getDecay();
                const double rdecay_row = 0.5/decay_row;
                const double powfactor_row=itr->getPowfactor();
                for ( AOShell::GaussianIterator itc = shell_col->begin(); itc != shell_col->end(); ++itc){
                    
                     // get decay constants 
                        const double decay_col = itc->getDecay();
                        const double rdecay_col = 0.5/decay_col; 
                       const double powfactor_col=itc->getPowfactor();
                      
                         
                         tensor3d cou(boost::extents[nrows][ncols][nextra]);
                                  
                           for (index3d i = 0; i != nrows; ++i) {
                               for (index3d j = 0; j != ncols; ++j) {
                                   for (index3d k = 0; k != nextra; ++k) {
                                       cou[i][j][k] = 0.0;
                                   }
                               }
                           }

                         
                         
            const double decay = decay_row + decay_col; 
            const double r_decay = 0.5/decay; 
            const double r_decay_2 = 2.*r_decay; 
            const double fac_a_ac = decay_row/decay; 
            const double fac_c_ac = decay_col/decay; 
 
            
            const double wmp0 = r_decay_2 * (decay_row * pos_row.getX() + decay_col * pos_col.getX()) - pos_row.getX();
            const double wmp1 = r_decay_2 * (decay_row * pos_row.getY() + decay_col * pos_col.getY()) - pos_row.getY();
            const double wmp2 = r_decay_2 * (decay_row * pos_row.getZ() + decay_col * pos_col.getZ()) - pos_row.getZ();

            const double wmq0 = r_decay_2 * (decay_row * pos_row.getX() + decay_col * pos_col.getX()) - pos_col.getX();
            const double wmq1 = r_decay_2 * (decay_row * pos_row.getY() + decay_col * pos_col.getY()) - pos_col.getY();
            const double wmq2 = r_decay_2 * (decay_row * pos_row.getZ() + decay_col * pos_col.getZ()) - pos_col.getZ();

            
            const double T = fac_a_ac * decay_col * distsq;


            double fak = 2.0 * pow(pi, 2.5) / (decay_row * decay_col * sqrt(decay_row + decay_col));
            fak = fak *  powfactor_col*powfactor_row;

         
            const std::vector<double> FmT=XIntegrate(nextra, T);

            // get initial data from FmT -> s-s element
            for (index3d i = 0; i != nextra; ++i) {
                cou[0][0][i] = fak * FmT[i];
            }


         
//Integrals     p - s
if (lmax_row > 0) {
  for (int m = 0; m < mmax; m++) {
    cou[Cart::x][0][m] = wmp0*cou[0][0][m+1];
    cou[Cart::y][0][m] = wmp1*cou[0][0][m+1];
    cou[Cart::z][0][m] = wmp2*cou[0][0][m+1];
  }
}
//------------------------------------------------------

//Integrals     d - s
if (lmax_row > 1) {
  for (int m = 0; m < mmax-1; m++) {
    double term = rdecay_row*(cou[0][0][m]-fac_c_ac*cou[0][0][m+1]);
    cou[Cart::xx][0][m] = wmp0*cou[Cart::x][0][m+1] + term;
    cou[Cart::xy][0][m] = wmp0*cou[Cart::y][0][m+1];
    cou[Cart::xz][0][m] = wmp0*cou[Cart::z][0][m+1];
    cou[Cart::yy][0][m] = wmp1*cou[Cart::y][0][m+1] + term;
    cou[Cart::yz][0][m] = wmp1*cou[Cart::z][0][m+1];
    cou[Cart::zz][0][m] = wmp2*cou[Cart::z][0][m+1] + term;
  }
}
//------------------------------------------------------

//Integrals     f - s
if (lmax_row > 2) {
  for (int m = 0; m < mmax-2; m++) {
    cou[Cart::xxx][0][m] = wmp0*cou[Cart::xx][0][m+1] + 2*rdecay_row*(cou[Cart::x][0][m]-fac_c_ac*cou[Cart::x][0][m+1]);
    cou[Cart::xxy][0][m] = wmp1*cou[Cart::xx][0][m+1];
    cou[Cart::xxz][0][m] = wmp2*cou[Cart::xx][0][m+1];
    cou[Cart::xyy][0][m] = wmp0*cou[Cart::yy][0][m+1];
    cou[Cart::xyz][0][m] = wmp0*cou[Cart::yz][0][m+1];
    cou[Cart::xzz][0][m] = wmp0*cou[Cart::zz][0][m+1];
    cou[Cart::yyy][0][m] = wmp1*cou[Cart::yy][0][m+1] + 2*rdecay_row*(cou[Cart::y][0][m]-fac_c_ac*cou[Cart::y][0][m+1]);
    cou[Cart::yyz][0][m] = wmp2*cou[Cart::yy][0][m+1];
    cou[Cart::yzz][0][m] = wmp1*cou[Cart::zz][0][m+1];
    cou[Cart::zzz][0][m] = wmp2*cou[Cart::zz][0][m+1] + 2*rdecay_row*(cou[Cart::z][0][m]-fac_c_ac*cou[Cart::z][0][m+1]);
  }
}
//------------------------------------------------------

//Integrals     g - s
if (lmax_row > 3) {
  for (int m = 0; m < mmax-3; m++) {
    double term_xx = rdecay_row*(cou[Cart::xx][0][m]-fac_c_ac*cou[Cart::xx][0][m+1]);
    double term_yy = rdecay_row*(cou[Cart::yy][0][m]-fac_c_ac*cou[Cart::yy][0][m+1]);
    double term_zz = rdecay_row*(cou[Cart::zz][0][m]-fac_c_ac*cou[Cart::zz][0][m+1]);
    cou[Cart::xxxx][0][m] = wmp0*cou[Cart::xxx][0][m+1] + 3*term_xx;
    cou[Cart::xxxy][0][m] = wmp1*cou[Cart::xxx][0][m+1];
    cou[Cart::xxxz][0][m] = wmp2*cou[Cart::xxx][0][m+1];
    cou[Cart::xxyy][0][m] = wmp0*cou[Cart::xyy][0][m+1] + term_yy;
    cou[Cart::xxyz][0][m] = wmp1*cou[Cart::xxz][0][m+1];
    cou[Cart::xxzz][0][m] = wmp0*cou[Cart::xzz][0][m+1] + term_zz;
    cou[Cart::xyyy][0][m] = wmp0*cou[Cart::yyy][0][m+1];
    cou[Cart::xyyz][0][m] = wmp0*cou[Cart::yyz][0][m+1];
    cou[Cart::xyzz][0][m] = wmp0*cou[Cart::yzz][0][m+1];
    cou[Cart::xzzz][0][m] = wmp0*cou[Cart::zzz][0][m+1];
    cou[Cart::yyyy][0][m] = wmp1*cou[Cart::yyy][0][m+1] + 3*term_yy;
    cou[Cart::yyyz][0][m] = wmp2*cou[Cart::yyy][0][m+1];
    cou[Cart::yyzz][0][m] = wmp1*cou[Cart::yzz][0][m+1] + term_zz;
    cou[Cart::yzzz][0][m] = wmp1*cou[Cart::zzz][0][m+1];
    cou[Cart::zzzz][0][m] = wmp2*cou[Cart::zzz][0][m+1] + 3*term_zz;
  }
}
//------------------------------------------------------

//Integrals     h - s
if (lmax_row > 4) {
  for (int m = 0; m < mmax-4; m++) {
    double term_xxx = rdecay_row*(cou[Cart::xxx][0][m]-fac_c_ac*cou[Cart::xxx][0][m+1]);
    double term_yyy = rdecay_row*(cou[Cart::yyy][0][m]-fac_c_ac*cou[Cart::yyy][0][m+1]);
    double term_zzz = rdecay_row*(cou[Cart::zzz][0][m]-fac_c_ac*cou[Cart::zzz][0][m+1]);
    cou[Cart::xxxxx][0][m] = wmp0*cou[Cart::xxxx][0][m+1] + 4*term_xxx;
    cou[Cart::xxxxy][0][m] = wmp1*cou[Cart::xxxx][0][m+1];
    cou[Cart::xxxxz][0][m] = wmp2*cou[Cart::xxxx][0][m+1];
    cou[Cart::xxxyy][0][m] = wmp1*cou[Cart::xxxy][0][m+1] + term_xxx;
    cou[Cart::xxxyz][0][m] = wmp1*cou[Cart::xxxz][0][m+1];
    cou[Cart::xxxzz][0][m] = wmp2*cou[Cart::xxxz][0][m+1] + term_xxx;
    cou[Cart::xxyyy][0][m] = wmp0*cou[Cart::xyyy][0][m+1] + term_yyy;
    cou[Cart::xxyyz][0][m] = wmp2*cou[Cart::xxyy][0][m+1];
    cou[Cart::xxyzz][0][m] = wmp1*cou[Cart::xxzz][0][m+1];
    cou[Cart::xxzzz][0][m] = wmp0*cou[Cart::xzzz][0][m+1] + term_zzz;
    cou[Cart::xyyyy][0][m] = wmp0*cou[Cart::yyyy][0][m+1];
    cou[Cart::xyyyz][0][m] = wmp0*cou[Cart::yyyz][0][m+1];
    cou[Cart::xyyzz][0][m] = wmp0*cou[Cart::yyzz][0][m+1];
    cou[Cart::xyzzz][0][m] = wmp0*cou[Cart::yzzz][0][m+1];
    cou[Cart::xzzzz][0][m] = wmp0*cou[Cart::zzzz][0][m+1];
    cou[Cart::yyyyy][0][m] = wmp1*cou[Cart::yyyy][0][m+1] + 4*term_yyy;
    cou[Cart::yyyyz][0][m] = wmp2*cou[Cart::yyyy][0][m+1];
    cou[Cart::yyyzz][0][m] = wmp2*cou[Cart::yyyz][0][m+1] + term_yyy;
    cou[Cart::yyzzz][0][m] = wmp1*cou[Cart::yzzz][0][m+1] + term_zzz;
    cou[Cart::yzzzz][0][m] = wmp1*cou[Cart::zzzz][0][m+1];
    cou[Cart::zzzzz][0][m] = wmp2*cou[Cart::zzzz][0][m+1] + 4*term_zzz;
  }
}
//------------------------------------------------------

//Integrals     i - s
if (lmax_row > 5) {
  for (int m = 0; m < mmax-5; m++) {
    double term_xxxx = rdecay_row*(cou[Cart::xxxx][0][m]-fac_c_ac*cou[Cart::xxxx][0][m+1]);
    double term_xyyy = rdecay_row*(cou[Cart::xyyy][0][m]-fac_c_ac*cou[Cart::xyyy][0][m+1]);
    double term_xzzz = rdecay_row*(cou[Cart::xzzz][0][m]-fac_c_ac*cou[Cart::xzzz][0][m+1]);
    double term_yyyy = rdecay_row*(cou[Cart::yyyy][0][m]-fac_c_ac*cou[Cart::yyyy][0][m+1]);
    double term_yyzz = rdecay_row*(cou[Cart::yyzz][0][m]-fac_c_ac*cou[Cart::yyzz][0][m+1]);
    double term_yzzz = rdecay_row*(cou[Cart::yzzz][0][m]-fac_c_ac*cou[Cart::yzzz][0][m+1]);
    double term_zzzz = rdecay_row*(cou[Cart::zzzz][0][m]-fac_c_ac*cou[Cart::zzzz][0][m+1]);
    cou[Cart::xxxxxx][0][m] = wmp0*cou[Cart::xxxxx][0][m+1] + 5*term_xxxx;
    cou[Cart::xxxxxy][0][m] = wmp1*cou[Cart::xxxxx][0][m+1];
    cou[Cart::xxxxxz][0][m] = wmp2*cou[Cart::xxxxx][0][m+1];
    cou[Cart::xxxxyy][0][m] = wmp1*cou[Cart::xxxxy][0][m+1] + term_xxxx;
    cou[Cart::xxxxyz][0][m] = wmp1*cou[Cart::xxxxz][0][m+1];
    cou[Cart::xxxxzz][0][m] = wmp2*cou[Cart::xxxxz][0][m+1] + term_xxxx;
    cou[Cart::xxxyyy][0][m] = wmp0*cou[Cart::xxyyy][0][m+1] + 2*term_xyyy;
    cou[Cart::xxxyyz][0][m] = wmp2*cou[Cart::xxxyy][0][m+1];
    cou[Cart::xxxyzz][0][m] = wmp1*cou[Cart::xxxzz][0][m+1];
    cou[Cart::xxxzzz][0][m] = wmp0*cou[Cart::xxzzz][0][m+1] + 2*term_xzzz;
    cou[Cart::xxyyyy][0][m] = wmp0*cou[Cart::xyyyy][0][m+1] + term_yyyy;
    cou[Cart::xxyyyz][0][m] = wmp2*cou[Cart::xxyyy][0][m+1];
    cou[Cart::xxyyzz][0][m] = wmp0*cou[Cart::xyyzz][0][m+1] + term_yyzz;
    cou[Cart::xxyzzz][0][m] = wmp1*cou[Cart::xxzzz][0][m+1];
    cou[Cart::xxzzzz][0][m] = wmp0*cou[Cart::xzzzz][0][m+1] + term_zzzz;
    cou[Cart::xyyyyy][0][m] = wmp0*cou[Cart::yyyyy][0][m+1];
    cou[Cart::xyyyyz][0][m] = wmp0*cou[Cart::yyyyz][0][m+1];
    cou[Cart::xyyyzz][0][m] = wmp0*cou[Cart::yyyzz][0][m+1];
    cou[Cart::xyyzzz][0][m] = wmp0*cou[Cart::yyzzz][0][m+1];
    cou[Cart::xyzzzz][0][m] = wmp0*cou[Cart::yzzzz][0][m+1];
    cou[Cart::xzzzzz][0][m] = wmp0*cou[Cart::zzzzz][0][m+1];
    cou[Cart::yyyyyy][0][m] = wmp1*cou[Cart::yyyyy][0][m+1] + 5*term_yyyy;
    cou[Cart::yyyyyz][0][m] = wmp2*cou[Cart::yyyyy][0][m+1];
    cou[Cart::yyyyzz][0][m] = wmp2*cou[Cart::yyyyz][0][m+1] + term_yyyy;
    cou[Cart::yyyzzz][0][m] = wmp1*cou[Cart::yyzzz][0][m+1] + 2*term_yzzz;
    cou[Cart::yyzzzz][0][m] = wmp1*cou[Cart::yzzzz][0][m+1] + term_zzzz;
    cou[Cart::yzzzzz][0][m] = wmp1*cou[Cart::zzzzz][0][m+1];
    cou[Cart::zzzzzz][0][m] = wmp2*cou[Cart::zzzzz][0][m+1] + 5*term_zzzz;
  }
}
//------------------------------------------------------


if (lmax_col > 0) {

  //Integrals     s - p
  for (int m = 0; m < lmax_col; m++) {
    cou[0][Cart::x][m] = wmq0*cou[0][0][m+1];
    cou[0][Cart::y][m] = wmq1*cou[0][0][m+1];
    cou[0][Cart::z][m] = wmq2*cou[0][0][m+1];
  }
  //------------------------------------------------------

  //Integrals     p - p
  if (lmax_row > 0) {
    for (int m = 0; m < lmax_col; m++) {
      double term = r_decay*cou[0][0][m+1];
      for (int i =  1; i < 4; i++) {
        cou[i][Cart::x][m] = wmq0*cou[i][0][m+1] + nx[i]*term;
        cou[i][Cart::y][m] = wmq1*cou[i][0][m+1] + ny[i]*term;
        cou[i][Cart::z][m] = wmq2*cou[i][0][m+1] + nz[i]*term;
      }
    }
  }
  //------------------------------------------------------

  //Integrals     d - p     f - p     g - p     h - p     i - p
  for (int i_row = 2; i_row < lmax_row+1; i_row++) {
    for (int m = 0; m < lmax_col; m++) {
      for (int i =  4; i < n_orbitals[lmax_row]; i++) {
        cou[i][Cart::x][m] = wmq0*cou[i][0][m+1] + nx[i]*r_decay*cou[i_less_x[i]][0][m+1];
        cou[i][Cart::y][m] = wmq1*cou[i][0][m+1] + ny[i]*r_decay*cou[i_less_y[i]][0][m+1];
        cou[i][Cart::z][m] = wmq2*cou[i][0][m+1] + nz[i]*r_decay*cou[i_less_z[i]][0][m+1];
      }
    }
  }
  //------------------------------------------------------

} // end if (lmax_col > 0)


if (lmax_col > 1) {

  //Integrals     s - d
  for (int m = 0; m < lmax_col-1; m++) {
    double term = rdecay_col*(cou[0][0][m]-fac_a_ac*cou[0][0][m+1]);
    cou[0][Cart::xx][m] = wmq0*cou[0][Cart::x][m+1] + term;
    cou[0][Cart::xy][m] = wmq0*cou[0][Cart::y][m+1];
    cou[0][Cart::xz][m] = wmq0*cou[0][Cart::z][m+1];
    cou[0][Cart::yy][m] = wmq1*cou[0][Cart::y][m+1] + term;
    cou[0][Cart::yz][m] = wmq1*cou[0][Cart::z][m+1];
    cou[0][Cart::zz][m] = wmq2*cou[0][Cart::z][m+1] + term;
  }
  //------------------------------------------------------

  //Integrals     p - d     d - d     f - d     g - d     h - d     i - d
  for (int m = 0; m < lmax_col-1; m++) {
    for (int i =  1; i < n_orbitals[lmax_row]; i++) {
      double term = rdecay_col*(cou[i][0][m]-fac_a_ac*cou[i][0][m+1]);
      cou[i][Cart::xx][m] = wmq0*cou[i][Cart::x][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::x][m+1] + term;
      cou[i][Cart::xy][m] = wmq0*cou[i][Cart::y][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::y][m+1];
      cou[i][Cart::xz][m] = wmq0*cou[i][Cart::z][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::z][m+1];
      cou[i][Cart::yy][m] = wmq1*cou[i][Cart::y][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::y][m+1] + term;
      cou[i][Cart::yz][m] = wmq1*cou[i][Cart::z][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::z][m+1];
      cou[i][Cart::zz][m] = wmq2*cou[i][Cart::z][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::z][m+1] + term;
    }
  }
  //------------------------------------------------------

} // end if (lmax_col > 1)


if (lmax_col > 2) {

  //Integrals     s - f
  for (int m = 0; m < lmax_col-2; m++) {
    cou[0][Cart::xxx][m] = wmq0*cou[0][Cart::xx][m+1] + 2*rdecay_col*(cou[0][Cart::x][m]-fac_a_ac*cou[0][Cart::x][m+1]);
    cou[0][Cart::xxy][m] = wmq1*cou[0][Cart::xx][m+1];
    cou[0][Cart::xxz][m] = wmq2*cou[0][Cart::xx][m+1];
    cou[0][Cart::xyy][m] = wmq0*cou[0][Cart::yy][m+1];
    cou[0][Cart::xyz][m] = wmq0*cou[0][Cart::yz][m+1];
    cou[0][Cart::xzz][m] = wmq0*cou[0][Cart::zz][m+1];
    cou[0][Cart::yyy][m] = wmq1*cou[0][Cart::yy][m+1] + 2*rdecay_col*(cou[0][Cart::y][m]-fac_a_ac*cou[0][Cart::y][m+1]);
    cou[0][Cart::yyz][m] = wmq2*cou[0][Cart::yy][m+1];
    cou[0][Cart::yzz][m] = wmq1*cou[0][Cart::zz][m+1];
    cou[0][Cart::zzz][m] = wmq2*cou[0][Cart::zz][m+1] + 2*rdecay_col*(cou[0][Cart::z][m]-fac_a_ac*cou[0][Cart::z][m+1]);
  }
  //------------------------------------------------------

  //Integrals     p - f     d - f     f - f     g - f     h - f     i - f
  for (int m = 0; m < lmax_col-2; m++) {
    for (int i =  1; i < n_orbitals[lmax_row]; i++) {
      double term_x = 2*rdecay_col*(cou[i][Cart::x][m]-fac_a_ac*cou[i][Cart::x][m+1]);
      double term_y = 2*rdecay_col*(cou[i][Cart::y][m]-fac_a_ac*cou[i][Cart::y][m+1]);
      double term_z = 2*rdecay_col*(cou[i][Cart::z][m]-fac_a_ac*cou[i][Cart::z][m+1]);
      cou[i][Cart::xxx][m] = wmq0*cou[i][Cart::xx][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xx][m+1] + term_x;
      cou[i][Cart::xxy][m] = wmq1*cou[i][Cart::xx][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xx][m+1];
      cou[i][Cart::xxz][m] = wmq2*cou[i][Cart::xx][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xx][m+1];
      cou[i][Cart::xyy][m] = wmq0*cou[i][Cart::yy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yy][m+1];
      cou[i][Cart::xyz][m] = wmq0*cou[i][Cart::yz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yz][m+1];
      cou[i][Cart::xzz][m] = wmq0*cou[i][Cart::zz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::zz][m+1];
      cou[i][Cart::yyy][m] = wmq1*cou[i][Cart::yy][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yy][m+1] + term_y;
      cou[i][Cart::yyz][m] = wmq2*cou[i][Cart::yy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::yy][m+1];
      cou[i][Cart::yzz][m] = wmq1*cou[i][Cart::zz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::zz][m+1];
      cou[i][Cart::zzz][m] = wmq2*cou[i][Cart::zz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::zz][m+1] + term_z;
    }
  }
  //------------------------------------------------------

} // end if (lmax_col > 2)


if (lmax_col > 3) {

  //Integrals     s - g
  for (int m = 0; m < lmax_col-3; m++) {
    double term_xx = rdecay_col*(cou[0][Cart::xx][m]-fac_a_ac*cou[0][Cart::xx][m+1]);
    double term_yy = rdecay_col*(cou[0][Cart::yy][m]-fac_a_ac*cou[0][Cart::yy][m+1]);
    double term_zz = rdecay_col*(cou[0][Cart::zz][m]-fac_a_ac*cou[0][Cart::zz][m+1]);
    cou[0][Cart::xxxx][m] = wmq0*cou[0][Cart::xxx][m+1] + 3*term_xx;
    cou[0][Cart::xxxy][m] = wmq1*cou[0][Cart::xxx][m+1];
    cou[0][Cart::xxxz][m] = wmq2*cou[0][Cart::xxx][m+1];
    cou[0][Cart::xxyy][m] = wmq0*cou[0][Cart::xyy][m+1] + term_yy;
    cou[0][Cart::xxyz][m] = wmq1*cou[0][Cart::xxz][m+1];
    cou[0][Cart::xxzz][m] = wmq0*cou[0][Cart::xzz][m+1] + term_zz;
    cou[0][Cart::xyyy][m] = wmq0*cou[0][Cart::yyy][m+1];
    cou[0][Cart::xyyz][m] = wmq0*cou[0][Cart::yyz][m+1];
    cou[0][Cart::xyzz][m] = wmq0*cou[0][Cart::yzz][m+1];
    cou[0][Cart::xzzz][m] = wmq0*cou[0][Cart::zzz][m+1];
    cou[0][Cart::yyyy][m] = wmq1*cou[0][Cart::yyy][m+1] + 3*term_yy;
    cou[0][Cart::yyyz][m] = wmq2*cou[0][Cart::yyy][m+1];
    cou[0][Cart::yyzz][m] = wmq1*cou[0][Cart::yzz][m+1] + term_zz;
    cou[0][Cart::yzzz][m] = wmq1*cou[0][Cart::zzz][m+1];
    cou[0][Cart::zzzz][m] = wmq2*cou[0][Cart::zzz][m+1] + 3*term_zz;
  }
  //------------------------------------------------------

  //Integrals     p - g     d - g     f - g     g - g     h - g     i - g
  for (int m = 0; m < lmax_col-3; m++) {
    for (int i =  1; i < n_orbitals[lmax_row]; i++) {
      double term_xx = rdecay_col*(cou[i][Cart::xx][m]-fac_a_ac*cou[i][Cart::xx][m+1]);
      double term_yy = rdecay_col*(cou[i][Cart::yy][m]-fac_a_ac*cou[i][Cart::yy][m+1]);
      double term_zz = rdecay_col*(cou[i][Cart::zz][m]-fac_a_ac*cou[i][Cart::zz][m+1]);
      cou[i][Cart::xxxx][m] = wmq0*cou[i][Cart::xxx][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xxx][m+1] + 3*term_xx;
      cou[i][Cart::xxxy][m] = wmq1*cou[i][Cart::xxx][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxx][m+1];
      cou[i][Cart::xxxz][m] = wmq2*cou[i][Cart::xxx][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxx][m+1];
      cou[i][Cart::xxyy][m] = wmq0*cou[i][Cart::xyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xyy][m+1] + term_yy;
      cou[i][Cart::xxyz][m] = wmq1*cou[i][Cart::xxz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxz][m+1];
      cou[i][Cart::xxzz][m] = wmq0*cou[i][Cart::xzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xzz][m+1] + term_zz;
      cou[i][Cart::xyyy][m] = wmq0*cou[i][Cart::yyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyy][m+1];
      cou[i][Cart::xyyz][m] = wmq0*cou[i][Cart::yyz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyz][m+1];
      cou[i][Cart::xyzz][m] = wmq0*cou[i][Cart::yzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yzz][m+1];
      cou[i][Cart::xzzz][m] = wmq0*cou[i][Cart::zzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::zzz][m+1];
      cou[i][Cart::yyyy][m] = wmq1*cou[i][Cart::yyy][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yyy][m+1] + 3*term_yy;
      cou[i][Cart::yyyz][m] = wmq2*cou[i][Cart::yyy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::yyy][m+1];
      cou[i][Cart::yyzz][m] = wmq1*cou[i][Cart::yzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yzz][m+1] + term_zz;
      cou[i][Cart::yzzz][m] = wmq1*cou[i][Cart::zzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::zzz][m+1];
      cou[i][Cart::zzzz][m] = wmq2*cou[i][Cart::zzz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::zzz][m+1] + 3*term_zz;
    }
  }
  //------------------------------------------------------

} // end if (lmax_col > 3)


if (lmax_col > 4) {

  //Integrals     s - h
  for (int m = 0; m < lmax_col-4; m++) {
    double term_xxx = rdecay_col*(cou[0][Cart::xxx][m]-fac_a_ac*cou[0][Cart::xxx][m+1]);
    double term_yyy = rdecay_col*(cou[0][Cart::yyy][m]-fac_a_ac*cou[0][Cart::yyy][m+1]);
    double term_zzz = rdecay_col*(cou[0][Cart::zzz][m]-fac_a_ac*cou[0][Cart::zzz][m+1]);
    cou[0][Cart::xxxxx][m] = wmq0*cou[0][Cart::xxxx][m+1] + 4*term_xxx;
    cou[0][Cart::xxxxy][m] = wmq1*cou[0][Cart::xxxx][m+1];
    cou[0][Cart::xxxxz][m] = wmq2*cou[0][Cart::xxxx][m+1];
    cou[0][Cart::xxxyy][m] = wmq1*cou[0][Cart::xxxy][m+1] + term_xxx;
    cou[0][Cart::xxxyz][m] = wmq1*cou[0][Cart::xxxz][m+1];
    cou[0][Cart::xxxzz][m] = wmq2*cou[0][Cart::xxxz][m+1] + term_xxx;
    cou[0][Cart::xxyyy][m] = wmq0*cou[0][Cart::xyyy][m+1] + term_yyy;
    cou[0][Cart::xxyyz][m] = wmq2*cou[0][Cart::xxyy][m+1];
    cou[0][Cart::xxyzz][m] = wmq1*cou[0][Cart::xxzz][m+1];
    cou[0][Cart::xxzzz][m] = wmq0*cou[0][Cart::xzzz][m+1] + term_zzz;
    cou[0][Cart::xyyyy][m] = wmq0*cou[0][Cart::yyyy][m+1];
    cou[0][Cart::xyyyz][m] = wmq0*cou[0][Cart::yyyz][m+1];
    cou[0][Cart::xyyzz][m] = wmq0*cou[0][Cart::yyzz][m+1];
    cou[0][Cart::xyzzz][m] = wmq0*cou[0][Cart::yzzz][m+1];
    cou[0][Cart::xzzzz][m] = wmq0*cou[0][Cart::zzzz][m+1];
    cou[0][Cart::yyyyy][m] = wmq1*cou[0][Cart::yyyy][m+1] + 4*term_yyy;
    cou[0][Cart::yyyyz][m] = wmq2*cou[0][Cart::yyyy][m+1];
    cou[0][Cart::yyyzz][m] = wmq2*cou[0][Cart::yyyz][m+1] + term_yyy;
    cou[0][Cart::yyzzz][m] = wmq1*cou[0][Cart::yzzz][m+1] + term_zzz;
    cou[0][Cart::yzzzz][m] = wmq1*cou[0][Cart::zzzz][m+1];
    cou[0][Cart::zzzzz][m] = wmq2*cou[0][Cart::zzzz][m+1] + 4*term_zzz;
  }
  //------------------------------------------------------

  //Integrals     p - h     d - h     f - h     g - h     h - h     i - h
  for (int m = 0; m < lmax_col-4; m++) {
    for (int i =  1; i < n_orbitals[lmax_row]; i++) {
      double term_xxx = rdecay_col*(cou[i][Cart::xxx][m]-fac_a_ac*cou[i][Cart::xxx][m+1]);
      double term_yyy = rdecay_col*(cou[i][Cart::yyy][m]-fac_a_ac*cou[i][Cart::yyy][m+1]);
      double term_zzz = rdecay_col*(cou[i][Cart::zzz][m]-fac_a_ac*cou[i][Cart::zzz][m+1]);
      cou[i][Cart::xxxxx][m] = wmq0*cou[i][Cart::xxxx][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xxxx][m+1] + 4*term_xxx;
      cou[i][Cart::xxxxy][m] = wmq1*cou[i][Cart::xxxx][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxx][m+1];
      cou[i][Cart::xxxxz][m] = wmq2*cou[i][Cart::xxxx][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxxx][m+1];
      cou[i][Cart::xxxyy][m] = wmq1*cou[i][Cart::xxxy][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxy][m+1] + term_xxx;
      cou[i][Cart::xxxyz][m] = wmq1*cou[i][Cart::xxxz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxz][m+1];
      cou[i][Cart::xxxzz][m] = wmq2*cou[i][Cart::xxxz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxxz][m+1] + term_xxx;
      cou[i][Cart::xxyyy][m] = wmq0*cou[i][Cart::xyyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xyyy][m+1] + term_yyy;
      cou[i][Cart::xxyyz][m] = wmq2*cou[i][Cart::xxyy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxyy][m+1];
      cou[i][Cart::xxyzz][m] = wmq1*cou[i][Cart::xxzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxzz][m+1];
      cou[i][Cart::xxzzz][m] = wmq0*cou[i][Cart::xzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xzzz][m+1] + term_zzz;
      cou[i][Cart::xyyyy][m] = wmq0*cou[i][Cart::yyyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyyy][m+1];
      cou[i][Cart::xyyyz][m] = wmq0*cou[i][Cart::yyyz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyyz][m+1];
      cou[i][Cart::xyyzz][m] = wmq0*cou[i][Cart::yyzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyzz][m+1];
      cou[i][Cart::xyzzz][m] = wmq0*cou[i][Cart::yzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yzzz][m+1];
      cou[i][Cart::xzzzz][m] = wmq0*cou[i][Cart::zzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::zzzz][m+1];
      cou[i][Cart::yyyyy][m] = wmq1*cou[i][Cart::yyyy][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yyyy][m+1] + 4*term_yyy;
      cou[i][Cart::yyyyz][m] = wmq2*cou[i][Cart::yyyy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::yyyy][m+1];
      cou[i][Cart::yyyzz][m] = wmq2*cou[i][Cart::yyyz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::yyyz][m+1] + term_yyy;
      cou[i][Cart::yyzzz][m] = wmq1*cou[i][Cart::yzzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yzzz][m+1] + term_zzz;
      cou[i][Cart::yzzzz][m] = wmq1*cou[i][Cart::zzzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::zzzz][m+1];
      cou[i][Cart::zzzzz][m] = wmq2*cou[i][Cart::zzzz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::zzzz][m+1] + 4*term_zzz;
    }
  }
  //------------------------------------------------------

} // end if (lmax_col > 4)


if (lmax_col > 5) {

  //Integrals     s - i
  for (int m = 0; m < lmax_col-5; m++) {
    double term_xxxx = rdecay_col*(cou[0][Cart::xxxx][m]-fac_a_ac*cou[0][Cart::xxxx][m+1]);
    double term_xyyy = rdecay_col*(cou[0][Cart::xyyy][m]-fac_a_ac*cou[0][Cart::xyyy][m+1]);
    double term_xzzz = rdecay_col*(cou[0][Cart::xzzz][m]-fac_a_ac*cou[0][Cart::xzzz][m+1]);
    double term_yyyy = rdecay_col*(cou[0][Cart::yyyy][m]-fac_a_ac*cou[0][Cart::yyyy][m+1]);
    double term_yyzz = rdecay_col*(cou[0][Cart::yyzz][m]-fac_a_ac*cou[0][Cart::yyzz][m+1]);
    double term_yzzz = rdecay_col*(cou[0][Cart::yzzz][m]-fac_a_ac*cou[0][Cart::yzzz][m+1]);
    double term_zzzz = rdecay_col*(cou[0][Cart::zzzz][m]-fac_a_ac*cou[0][Cart::zzzz][m+1]);
    cou[0][Cart::xxxxxx][m] = wmq0*cou[0][Cart::xxxxx][m+1] + 5*term_xxxx;
    cou[0][Cart::xxxxxy][m] = wmq1*cou[0][Cart::xxxxx][m+1];
    cou[0][Cart::xxxxxz][m] = wmq2*cou[0][Cart::xxxxx][m+1];
    cou[0][Cart::xxxxyy][m] = wmq1*cou[0][Cart::xxxxy][m+1] + term_xxxx;
    cou[0][Cart::xxxxyz][m] = wmq1*cou[0][Cart::xxxxz][m+1];
    cou[0][Cart::xxxxzz][m] = wmq2*cou[0][Cart::xxxxz][m+1] + term_xxxx;
    cou[0][Cart::xxxyyy][m] = wmq0*cou[0][Cart::xxyyy][m+1] + 2*term_xyyy;
    cou[0][Cart::xxxyyz][m] = wmq2*cou[0][Cart::xxxyy][m+1];
    cou[0][Cart::xxxyzz][m] = wmq1*cou[0][Cart::xxxzz][m+1];
    cou[0][Cart::xxxzzz][m] = wmq0*cou[0][Cart::xxzzz][m+1] + 2*term_xzzz;
    cou[0][Cart::xxyyyy][m] = wmq0*cou[0][Cart::xyyyy][m+1] + term_yyyy;
    cou[0][Cart::xxyyyz][m] = wmq2*cou[0][Cart::xxyyy][m+1];
    cou[0][Cart::xxyyzz][m] = wmq0*cou[0][Cart::xyyzz][m+1] + term_yyzz;
    cou[0][Cart::xxyzzz][m] = wmq1*cou[0][Cart::xxzzz][m+1];
    cou[0][Cart::xxzzzz][m] = wmq0*cou[0][Cart::xzzzz][m+1] + term_zzzz;
    cou[0][Cart::xyyyyy][m] = wmq0*cou[0][Cart::yyyyy][m+1];
    cou[0][Cart::xyyyyz][m] = wmq0*cou[0][Cart::yyyyz][m+1];
    cou[0][Cart::xyyyzz][m] = wmq0*cou[0][Cart::yyyzz][m+1];
    cou[0][Cart::xyyzzz][m] = wmq0*cou[0][Cart::yyzzz][m+1];
    cou[0][Cart::xyzzzz][m] = wmq0*cou[0][Cart::yzzzz][m+1];
    cou[0][Cart::xzzzzz][m] = wmq0*cou[0][Cart::zzzzz][m+1];
    cou[0][Cart::yyyyyy][m] = wmq1*cou[0][Cart::yyyyy][m+1] + 5*term_yyyy;
    cou[0][Cart::yyyyyz][m] = wmq2*cou[0][Cart::yyyyy][m+1];
    cou[0][Cart::yyyyzz][m] = wmq2*cou[0][Cart::yyyyz][m+1] + term_yyyy;
    cou[0][Cart::yyyzzz][m] = wmq1*cou[0][Cart::yyzzz][m+1] + 2*term_yzzz;
    cou[0][Cart::yyzzzz][m] = wmq1*cou[0][Cart::yzzzz][m+1] + term_zzzz;
    cou[0][Cart::yzzzzz][m] = wmq1*cou[0][Cart::zzzzz][m+1];
    cou[0][Cart::zzzzzz][m] = wmq2*cou[0][Cart::zzzzz][m+1] + 5*term_zzzz;
  }
  //------------------------------------------------------

  //Integrals     p - i     d - i     f - i     g - i     h - i     i - i
  for (int m = 0; m < lmax_col-5; m++) {
    for (int i =  1; i < n_orbitals[lmax_row]; i++) {
      double term_xxxx = rdecay_col*(cou[i][Cart::xxxx][m]-fac_a_ac*cou[i][Cart::xxxx][m+1]);
      double term_xyyy = rdecay_col*(cou[i][Cart::xyyy][m]-fac_a_ac*cou[i][Cart::xyyy][m+1]);
      double term_xzzz = rdecay_col*(cou[i][Cart::xzzz][m]-fac_a_ac*cou[i][Cart::xzzz][m+1]);
      double term_yyyy = rdecay_col*(cou[i][Cart::yyyy][m]-fac_a_ac*cou[i][Cart::yyyy][m+1]);
      double term_yyzz = rdecay_col*(cou[i][Cart::yyzz][m]-fac_a_ac*cou[i][Cart::yyzz][m+1]);
      double term_yzzz = rdecay_col*(cou[i][Cart::yzzz][m]-fac_a_ac*cou[i][Cart::yzzz][m+1]);
      double term_zzzz = rdecay_col*(cou[i][Cart::zzzz][m]-fac_a_ac*cou[i][Cart::zzzz][m+1]);
      cou[i][Cart::xxxxxx][m] = wmq0*cou[i][Cart::xxxxx][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xxxxx][m+1] + 5*term_xxxx;
      cou[i][Cart::xxxxxy][m] = wmq1*cou[i][Cart::xxxxx][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxxx][m+1];
      cou[i][Cart::xxxxxz][m] = wmq2*cou[i][Cart::xxxxx][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxxxx][m+1];
      cou[i][Cart::xxxxyy][m] = wmq1*cou[i][Cart::xxxxy][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxxy][m+1] + term_xxxx;
      cou[i][Cart::xxxxyz][m] = wmq1*cou[i][Cart::xxxxz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxxz][m+1];
      cou[i][Cart::xxxxzz][m] = wmq2*cou[i][Cart::xxxxz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxxxz][m+1] + term_xxxx;
      cou[i][Cart::xxxyyy][m] = wmq0*cou[i][Cart::xxyyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xxyyy][m+1] + 2*term_xyyy;
      cou[i][Cart::xxxyyz][m] = wmq2*cou[i][Cart::xxxyy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxxyy][m+1];
      cou[i][Cart::xxxyzz][m] = wmq1*cou[i][Cart::xxxzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxxzz][m+1];
      cou[i][Cart::xxxzzz][m] = wmq0*cou[i][Cart::xxzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xxzzz][m+1] + 2*term_xzzz;
      cou[i][Cart::xxyyyy][m] = wmq0*cou[i][Cart::xyyyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xyyyy][m+1] + term_yyyy;
      cou[i][Cart::xxyyyz][m] = wmq2*cou[i][Cart::xxyyy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::xxyyy][m+1];
      cou[i][Cart::xxyyzz][m] = wmq0*cou[i][Cart::xyyzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xyyzz][m+1] + term_yyzz;
      cou[i][Cart::xxyzzz][m] = wmq1*cou[i][Cart::xxzzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::xxzzz][m+1];
      cou[i][Cart::xxzzzz][m] = wmq0*cou[i][Cart::xzzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::xzzzz][m+1] + term_zzzz;
      cou[i][Cart::xyyyyy][m] = wmq0*cou[i][Cart::yyyyy][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyyyy][m+1];
      cou[i][Cart::xyyyyz][m] = wmq0*cou[i][Cart::yyyyz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyyyz][m+1];
      cou[i][Cart::xyyyzz][m] = wmq0*cou[i][Cart::yyyzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyyzz][m+1];
      cou[i][Cart::xyyzzz][m] = wmq0*cou[i][Cart::yyzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yyzzz][m+1];
      cou[i][Cart::xyzzzz][m] = wmq0*cou[i][Cart::yzzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::yzzzz][m+1];
      cou[i][Cart::xzzzzz][m] = wmq0*cou[i][Cart::zzzzz][m+1] + nx[i]*r_decay*cou[i_less_x[i]][Cart::zzzzz][m+1];
      cou[i][Cart::yyyyyy][m] = wmq1*cou[i][Cart::yyyyy][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yyyyy][m+1] + 5*term_yyyy;
      cou[i][Cart::yyyyyz][m] = wmq2*cou[i][Cart::yyyyy][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::yyyyy][m+1];
      cou[i][Cart::yyyyzz][m] = wmq2*cou[i][Cart::yyyyz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::yyyyz][m+1] + term_yyyy;
      cou[i][Cart::yyyzzz][m] = wmq1*cou[i][Cart::yyzzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yyzzz][m+1] + 2*term_yzzz;
      cou[i][Cart::yyzzzz][m] = wmq1*cou[i][Cart::yzzzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::yzzzz][m+1] + term_zzzz;
      cou[i][Cart::yzzzzz][m] = wmq1*cou[i][Cart::zzzzz][m+1] + ny[i]*r_decay*cou[i_less_y[i]][Cart::zzzzz][m+1];
      cou[i][Cart::zzzzzz][m] = wmq2*cou[i][Cart::zzzzz][m+1] + nz[i]*r_decay*cou[i_less_z[i]][Cart::zzzzz][m+1] + 5*term_zzzz;
    }
  }
  //------------------------------------------------------

} // end if (lmax_col > 5)
 
 
         
            // normalization and cartesian -> spherical factors
            int ntrafo_row = shell_row->getNumFunc() + shell_row->getOffset();
            int ntrafo_col = shell_col->getNumFunc() + shell_col->getOffset();

            

            // get transformation matrices including contraction coefficients
          const std::vector<double>& contractions_row = itr->getContraction();
          const std::vector<double>& contractions_col = itc->getContraction();

          

            // put cou[i][j][0] into eigen matrix
            Eigen::MatrixXd coumat = Eigen::MatrixXd::Zero(nrows, ncols);
            for (unsigned i = 0; i < coumat.rows(); i++) {
                for (unsigned j = 0; j < coumat.cols(); j++) {
                    coumat(i, j) = cou[i][j][0];
                }
            }

            Eigen::MatrixXd cou_tmp = Eigen::MatrixXd::Zero(ntrafo_row, ncols);
            
            
              // s-functions
            double factor = contractions_row[0];
            for (int i =  0; i < ncols; i++) {
              cou_tmp(0,i) = factor * coumat(0,i); /// Y 0,0
            }

            if (lmax_row > 0) {
              // p-functions
              factor = 2.0*sqrt(decay_row)*contractions_row[1];
              for (int i =  0; i < ncols; i++) {
                cou_tmp(1,i) = factor * coumat(3,i); /// Y 1,0
                cou_tmp(2,i) = factor * coumat(2,i); /// Y 1,-1
                cou_tmp(3,i) = factor * coumat(1,i); /// Y 1,1
              }
            }

            if (lmax_row > 1) {
              // d-functions
              factor = 2.0*decay_row*contractions_row[2];
              double factor_1 =  factor/sqrt(3.0);
              for (int i =  0; i < ncols; i++) {
                cou_tmp(4,i) = factor_1 * ( 2.0*coumat(Cart::zz,i) - coumat(Cart::xx,i) - coumat(Cart::yy,i) );  /// d3z2-r2  Y 2,0
                cou_tmp(5,i) = 2.*factor * coumat(Cart::yz,i);  /// dyz  Y 2,-1
                cou_tmp(6,i) = 2.*factor * coumat(Cart::xz,i);  /// dxz  Y 2,1
                cou_tmp(7,i) = 2.*factor * coumat(Cart::xy,i);  /// dxy  Y 2,-2
                cou_tmp(8,i) = factor * ( coumat(Cart::xx,i) - coumat(Cart::yy,i) );  /// dx2-y2  Y 2,2
              }
            }

            if (lmax_row > 2) {
              // f-functions
              factor = 2.0*pow(decay_row,1.5)*contractions_row[3];
              double factor_1 = factor*2./sqrt(15.);
              double factor_2 = factor*sqrt(2.)/sqrt(5.);
              double factor_3 = factor*sqrt(2.)/sqrt(3.);
              for (int i =  0; i < ncols; i++) {
                cou_tmp(9,i) = factor_1 * ( 2.*coumat(Cart::zzz,i) - 3.*coumat(Cart::xxz,i) - 3.* coumat(Cart::yyz,i) ); /// Y 3,0
                cou_tmp(10,i) = factor_2 * ( 4.*coumat(Cart::yzz,i) - coumat(Cart::xxy,i) - coumat(Cart::yyy,i) ); /// Y 3,-1
                cou_tmp(11,i) = factor_2 * ( 4.*coumat(Cart::xzz,i) - coumat(Cart::xxx,i) - coumat(Cart::xyy,i) ); /// Y 3,1
                cou_tmp(12,i) = 4.*factor * coumat(Cart::xyz,i); /// Y 3,-2
                cou_tmp(13,i) = 2.*factor * ( coumat(Cart::xxz,i) - coumat(Cart::yyz,i) ); /// Y 3,2
                cou_tmp(14,i) = factor_3 * ( 3.*coumat(Cart::xxy,i) - coumat(Cart::yyy,i) ); /// Y 3,-3
                cou_tmp(15,i) = factor_3 * ( coumat(Cart::xxx,i) - 3.*coumat(Cart::xyy,i) ); /// Y 3,3
              }
            }

            if (lmax_row > 3) {
              // g-functions
              factor = 2./sqrt(3.)*decay_row*decay_row*contractions_row[4];
              double factor_1 = factor/sqrt(35.);
              double factor_2 = factor*4./sqrt(14.);
              double factor_3 = factor*2./sqrt(7.);
              double factor_4 = factor*2.*sqrt(2.);
              for (int i =  0; i < ncols; i++) {
                cou_tmp(16,i) = factor_1 * (    3.*(coumat(Cart::xxxx,i) + coumat(Cart::yyyy,i))
                                                 + 6.*coumat(Cart::xxyy,i)
                                               - 24.*(coumat(Cart::xxzz,i) + coumat(Cart::yyzz,i))
                                                 + 8.*coumat(Cart::zzzz,i) );                               /// Y 4,0
                cou_tmp(17,i) = factor_2 * ( -3.*(coumat(Cart::xxyz,i) + coumat(Cart::yyyz,i))
                                               + 4.*coumat(Cart::yzzz,i) );                                 /// Y 4,-1
                cou_tmp(18,i) = factor_2 * ( -3.*(coumat(Cart::xxxz,i) + coumat(Cart::xyyz,i))
                                               + 4.*coumat(Cart::xzzz,i) );                                 /// Y 4,1
                cou_tmp(19,i) = 2.*factor_3 * (    -coumat(Cart::xxxy,i)
                                                     - coumat(Cart::xyyy,i)
                                                  + 6.*coumat(Cart::xyzz,i) );                              /// Y 4,-2
                cou_tmp(20,i) = factor_3 * (      -coumat(Cart::xxxx,i)
                                               + 6.*(coumat(Cart::xxzz,i) - coumat(Cart::yyzz,i))
                                                  + coumat(Cart::yyyy,i) );                                 /// Y 4,2
                cou_tmp(21,i) = factor_4 * ( 3.*coumat(Cart::xxyz,i) 
                                                - coumat(Cart::yyyz,i) );                                   /// Y 4,-3
                cou_tmp(22,i) = factor_4 * (      coumat(Cart::xxxz,i) 
                                               - 3.*coumat(Cart::xyyz,i) );                                 /// Y 4,3
                cou_tmp(23,i) = 4.*factor * (   coumat(Cart::xxxy,i)
                                                - coumat(Cart::xyyy,i) );                                   /// Y 4,-4
                cou_tmp(24,i) = factor * (      coumat(Cart::xxxx,i) 
                                             - 6.*coumat(Cart::xxyy,i)
                                                + coumat(Cart::yyyy,i) );                                   /// Y 4,4
              }
            }

            if (lmax_row > 4) {
              // h-functions
              factor = (2./3.)*pow(decay_row,2.5)*contractions_row[5];
              double factor_1 = factor*2./sqrt(105.);
              double factor_2 = factor*2./sqrt(7.);
              double factor_3 = factor*sqrt(6.)/3.;
              double factor_4 = factor*2.*sqrt(3.);
              double factor_5 = factor*.2*sqrt(30.);
              for (int i =  0; i < ncols; i++) {
                cou_tmp(25,i) = factor_1 * (   15.*(coumat(Cart::xxxxz,i) + coumat(Cart::yyyyz,i))
                                                + 30.*coumat(Cart::xxyyz,i)
                                               - 40.*(coumat(Cart::xxzzz,i) + coumat(Cart::yyzzz,i))
                                                 + 8.*coumat(Cart::zzzzz,i) );                              /// Y 5,0

                cou_tmp(26,i) = factor_2 * (        coumat(Cart::xxxxy,i)
                                                 + 2.*coumat(Cart::xxyyy,i)
                                               - 12.*(coumat(Cart::xxyzz,i) + coumat(Cart::yyyzz,i))
                                                    + coumat(Cart::yyyyy,i)
                                                 + 8.*coumat(Cart::yzzzz,i) );                              /// Y 5,-1

                cou_tmp(27,i) = factor_2 * (        coumat(Cart::xxxxx,i)
                                                 + 2.*coumat(Cart::xxxyy,i)
                                               - 12.*(coumat(Cart::xxxzz,i) + coumat(Cart::xyyzz,i))
                                                    + coumat(Cart::xyyyy,i)
                                                 + 8.*coumat(Cart::xzzzz,i) );                              /// Y 5,1

                cou_tmp(28,i) = 8.*factor * (     -coumat(Cart::xxxyz,i)
                                                   - coumat(Cart::xyyyz,i)
                                                + 2.*coumat(Cart::xyzzz,i) );                               /// Y 5,-2

                cou_tmp(29,i) = 4.*factor * (      -coumat(Cart::xxxxz,i)
                                                + 2.*(coumat(Cart::xxzzz,i) - coumat(Cart::yyzzz,i))
                                                    + coumat(Cart::yyyyz,i) );                              /// Y 5,2

                cou_tmp(30,i) = factor_3 * (   -3.*coumat(Cart::xxxxy,i)
                                                - 2.*coumat(Cart::xxyyy,i)
                                               + 24.*coumat(Cart::xxyzz,i)
                                                   + coumat(Cart::yyyyy,i)
                                                - 8.*coumat(Cart::yyyzz,i) );                               /// Y 5,-3

                cou_tmp(31,i) = factor_3 * (      -coumat(Cart::xxxxx,i)
                                                + 2.*coumat(Cart::xxxyy,i)
                                                + 8.*coumat(Cart::xxxzz,i)
                                                + 3.*coumat(Cart::xyyyy,i)
                                               - 24.*coumat(Cart::xyyzz,i) );                               /// Y 5,3

                cou_tmp(32,i) = 4.*factor_4 * (   coumat(Cart::xxxyz,i)
                                                  - coumat(Cart::xyyyz,i) );                                /// Y 5,-4

                cou_tmp(33,i) = factor_4 * (      coumat(Cart::xxxxz,i)
                                               - 6.*coumat(Cart::xxyyz,i)
                                                  + coumat(Cart::yyyyz,i) );                                /// Y 5,4

                cou_tmp(34,i) = factor_5 * (    5.*coumat(Cart::xxxxy,i)
                                               - 10.*coumat(Cart::xxyyy,i)
                                                   + coumat(Cart::yyyyy,i) );                               /// Y 5,-5

                cou_tmp(35,i) = factor_5 * (       coumat(Cart::xxxxx,i)
                                               - 10.*coumat(Cart::xxxyy,i)
                                                + 5.*coumat(Cart::xyyyy,i) );                               /// Y 5,5
              }
            }


            if (lmax_row > 5) {
              // i-functions
              factor = (2./3.)*decay_row*decay_row*decay_row*contractions_row[6];
              double factor_1 = factor*2./sqrt(1155.);
              double factor_2 = factor*4./sqrt(55.);
              double factor_3 = factor*sqrt(22.)/11.;
              double factor_4 = factor*2.*sqrt(165.)/55.;
              double factor_5 = factor*.4*sqrt(30.);
              double factor_6 = factor*.2*sqrt(10.);
              for (int i =  0; i < ncols; i++) {
                cou_tmp(36,i) = factor_1 * (    -5.*(coumat(Cart::xxxxxx,i) + coumat(Cart::yyyyyy,i))
                                                - 15.*(coumat(Cart::xxxxyy,i) + coumat(Cart::xxyyyy,i))
                                                + 90.*(coumat(Cart::xxxxzz,i) + coumat(Cart::yyyyzz,i))
                                                + 180.*coumat(Cart::xxyyzz,i)
                                               - 120.*(coumat(Cart::xxzzzz,i) + coumat(Cart::yyzzzz,i))
                                                 + 16.*coumat(Cart::zzzzzz,i) );                                /// Y 6,0

                cou_tmp(37,i) = factor_2 * (    5.*(coumat(Cart::xxxxyz,i) + coumat(Cart::yyyyyz,i))
                                                + 10.*coumat(Cart::xxyyyz,i)
                                               - 20.*(coumat(Cart::xxyzzz,i) + coumat(Cart::yyyzzz,i))
                                                 + 8.*coumat(Cart::yzzzzz,i) );                                 /// Y 6,-1

                cou_tmp(38,i) = factor_2 * (    5.*(coumat(Cart::xxxxxz,i) + coumat(Cart::xyyyyz,i))
                                                + 10.*coumat(Cart::xxxyyz,i)
                                               - 20.*(coumat(Cart::xxxzzz,i) + coumat(Cart::xyyzzz,i))
                                                 + 8.*coumat(Cart::xzzzzz,i) );                                 /// Y 6,1

                cou_tmp(39,i) = 2.*factor_3 * (        coumat(Cart::xxxxxy,i)
                                                    + 2.*coumat(Cart::xxxyyy,i)
                                                  - 16.*(coumat(Cart::xxxyzz,i) + coumat(Cart::xyyyzz,i) - coumat(Cart::xyzzzz,i))
                                                       + coumat(Cart::xyyyyy,i) );                              /// Y 6,-2

                cou_tmp(40,i) = factor_3 * (        coumat(Cart::xxxxxy,i)
                                                    + coumat(Cart::xxxxyy,i)
                                               - 16.*(coumat(Cart::xxxxzz,i) - coumat(Cart::xxzzzz,i)
                                                                               - coumat(Cart::yyyyzz,i) + coumat(Cart::yyzzzz,i))
                                                    - coumat(Cart::xxyyyy,i)
                                                    - coumat(Cart::yyyyyy,i) );                                 /// Y 6,2

                cou_tmp(41,i) = 2.*factor_3 * (   -9.*coumat(Cart::xxxxyz,i)
                                                   - 6.*coumat(Cart::xxyyyz,i)
                                                  + 24.*coumat(Cart::xxyzzz,i)
                                                   + 3.*coumat(Cart::yyyyyz,i)
                                                   - 8.*coumat(Cart::yyyzzz,i) );                               /// Y 6,-3

                cou_tmp(42,i) = 2.*factor_3 * (   -3.*coumat(Cart::xxxxxz,i)
                                                   + 6.*coumat(Cart::xxxyyz,i)
                                                   + 8.*coumat(Cart::xxxzzz,i)
                                                   + 9.*coumat(Cart::xyyyyz,i)
                                                  - 24.*coumat(Cart::xyyzzz,i) );                               /// Y 6,3

                cou_tmp(43,i) = 4.*factor_4 * (       -coumat(Cart::xxxxxy,i)
                                                  + 10.*(coumat(Cart::xxxyzz,i) - coumat(Cart::xyyyzz,i))
                                                       + coumat(Cart::xyyyyy,i) );                              /// Y 6,-4

                cou_tmp(44,i) = factor_4 * (       -coumat(Cart::xxxxxx,i)
                                                + 5.*(coumat(Cart::xxxxyy,i) + coumat(Cart::xxyyyy,i))
                                               + 10.*(coumat(Cart::xxxxzz,i) + coumat(Cart::yyyyzz,i))
                                                - 60.*coumat(Cart::xxyyzz,i)
                                                   -  coumat(Cart::yyyyyy,i) );                                 /// Y 6,4

                cou_tmp(45,i) = factor_5 * (    5.*coumat(Cart::xxxxyz,i)
                                               - 10.*coumat(Cart::xxyyyz,i)
                                                   + coumat(Cart::yyyyyz,i) );                                  /// Y 6,-5

                cou_tmp(46,i) = factor_5 * (       coumat(Cart::xxxxxz,i)
                                               - 10.*coumat(Cart::xxxyyz,i)
                                                + 5.*coumat(Cart::xyyyyz,i) );                                  /// Y 6,5

                cou_tmp(47,i) = 2.*factor_6 * (    3.*coumat(Cart::xxxxxy,i)
                                                  - 10.*coumat(Cart::xxxyyy,i)
                                                   + 3.*coumat(Cart::xyyyyy,i) );                               /// Y 6,-6

                cou_tmp(48,i) = factor_6 * (        coumat(Cart::xxxxxx,i)
                                               - 15.*(coumat(Cart::xxxxyy,i) - coumat(Cart::xxyyyy,i))
                                                    - coumat(Cart::yyyyyy,i) );                                 /// Y 6,6

              }
            }

                
                
            

    Eigen::MatrixXd cou_sph = Eigen::MatrixXd::Zero(ntrafo_row, ntrafo_col);  ////////////////////////////////////

        
              // s-functions
            factor = contractions_col[0];
            for (int i =  0; i < ntrafo_row; i++) {
              cou_sph(i,0) = factor * cou_tmp(i,0); /// Y 0,0
            }

            if (lmax_col > 0) {
              // p-functions
              factor = 2.0*sqrt(decay_col)*contractions_col[1];
              for (int i =  0; i < ntrafo_row; i++) {
                cou_sph(i,1) = factor * cou_tmp(i,3); /// Y 1,0
                cou_sph(i,2) = factor * cou_tmp(i,2); /// Y 1,-1
                cou_sph(i,3) = factor * cou_tmp(i,1); /// Y 1,1
              }
            }

            if (lmax_col > 1) {
              // d-functions
              factor = 2.0*decay_col*contractions_col[2];
              double factor_1 =  factor/sqrt(3.0);
              for (int i =  0; i < ntrafo_row; i++) {
                cou_sph(i,4) = factor_1 * ( 2.0*cou_tmp(i,Cart::zz) - cou_tmp(i,Cart::xx) - cou_tmp(i,Cart::yy) );  /// d3z2-r2  Y 2,0
                cou_sph(i,5) = 2.*factor * cou_tmp(i,Cart::yz);  /// dyz  Y 2,-1
                cou_sph(i,6) = 2.*factor * cou_tmp(i,Cart::xz);  /// dxz  Y 2,1
                cou_sph(i,7) = 2.*factor * cou_tmp(i,Cart::xy);  /// dxy  Y 2,-2
                cou_sph(i,8) = factor * ( cou_tmp(i,Cart::xx) - cou_tmp(i,Cart::yy) );  /// dx2-y2  Y 2,2
              }
            }

            if (lmax_col > 2) {
              // f-functions
              factor = 2.0*pow(decay_col,1.5)*contractions_col[3];
              double factor_1 = factor*2./sqrt(15.);
              double factor_2 = factor*sqrt(2.)/sqrt(5.);
              double factor_3 = factor*sqrt(2.)/sqrt(3.);
              for (int i =  0; i < ntrafo_row; i++) {
                cou_sph(i,9) = factor_1 * ( 2.*cou_tmp(i,Cart::zzz) - 3.*cou_tmp(i,Cart::xxz) - 3.* cou_tmp(i,Cart::yyz) ); /// Y 3,0
                cou_sph(i,10) = factor_2 * ( 4.*cou_tmp(i,Cart::yzz) - cou_tmp(i,Cart::xxy) - cou_tmp(i,Cart::yyy) ); /// Y 3,-1
                cou_sph(i,11) = factor_2 * ( 4.*cou_tmp(i,Cart::xzz) - cou_tmp(i,Cart::xxx) - cou_tmp(i,Cart::xyy) ); /// Y 3,1
                cou_sph(i,12) = 4.*factor * cou_tmp(i,Cart::xyz); /// Y 3,-2
                cou_sph(i,13) = 2.*factor * ( cou_tmp(i,Cart::xxz) - cou_tmp(i,Cart::yyz) ); /// Y 3,2
                cou_sph(i,14) = factor_3 * ( 3.*cou_tmp(i,Cart::xxy) - cou_tmp(i,Cart::yyy) ); /// Y 3,-3
                cou_sph(i,15) = factor_3 * ( cou_tmp(i,Cart::xxx) - 3.*cou_tmp(i,Cart::xyy) ); /// Y 3,3
              }
            }

            if (lmax_col > 3) {
              // g-functions
              factor = 2./sqrt(3.)*decay_col*decay_col*contractions_col[4];
              double factor_1 = factor/sqrt(35.);
              double factor_2 = factor*4./sqrt(14.);
              double factor_3 = factor*2./sqrt(7.);
              double factor_4 = factor*2.*sqrt(2.);
              for (int i =  0; i < ntrafo_row; i++) {
                cou_sph(i,16) = factor_1 * (    3.*(cou_tmp(i,Cart::xxxx) + cou_tmp(i,Cart::yyyy))
                                                 + 6.*cou_tmp(i,Cart::xxyy)
                                               - 24.*(cou_tmp(i,Cart::xxzz) + cou_tmp(i,Cart::yyzz))
                                                 + 8.*cou_tmp(i,Cart::zzzz) );                               /// Y 4,0
                cou_sph(i,17) = factor_2 * ( -3.*(cou_tmp(i,Cart::xxyz) + cou_tmp(i,Cart::yyyz))
                                               + 4.*cou_tmp(i,Cart::yzzz) );                                 /// Y 4,-1
                cou_sph(i,18) = factor_2 * ( -3.*(cou_tmp(i,Cart::xxxz) + cou_tmp(i,Cart::xyyz))
                                               + 4.*cou_tmp(i,Cart::xzzz) );                                 /// Y 4,1
                cou_sph(i,19) = 2.*factor_3 * (    -cou_tmp(i,Cart::xxxy)
                                                     - cou_tmp(i,Cart::xyyy)
                                                  + 6.*cou_tmp(i,Cart::xyzz) );                              /// Y 4,-2
                cou_sph(i,20) = factor_3 * (      -cou_tmp(i,Cart::xxxx)
                                               + 6.*(cou_tmp(i,Cart::xxzz) - cou_tmp(i,Cart::yyzz))
                                                  + cou_tmp(i,Cart::yyyy) );                                 /// Y 4,2
                cou_sph(i,21) = factor_4 * ( 3.*cou_tmp(i,Cart::xxyz) 
                                                - cou_tmp(i,Cart::yyyz) );                                   /// Y 4,-3
                cou_sph(i,22) = factor_4 * (      cou_tmp(i,Cart::xxxz) 
                                               - 3.*cou_tmp(i,Cart::xyyz) );                                 /// Y 4,3
                cou_sph(i,23) = 4.*factor * (   cou_tmp(i,Cart::xxxy)
                                                - cou_tmp(i,Cart::xyyy) );                                   /// Y 4,-4
                cou_sph(i,24) = factor * (      cou_tmp(i,Cart::xxxx) 
                                             - 6.*cou_tmp(i,Cart::xxyy)
                                                + cou_tmp(i,Cart::yyyy) );                                   /// Y 4,4
              }
            }

            if (lmax_col > 4) {
              // h-functions
              factor = (2./3.)*pow(decay_col,2.5)*contractions_col[5];
              double factor_1 = factor*2./sqrt(105.);
              double factor_2 = factor*2./sqrt(7.);
              double factor_3 = factor*sqrt(6.)/3.;
              double factor_4 = factor*2.*sqrt(3.);
              double factor_5 = factor*.2*sqrt(30.);
              for (int i =  0; i < ntrafo_row; i++) {
                cou_sph(i,25) = factor_1 * (   15.*(cou_tmp(i,Cart::xxxxz) + cou_tmp(i,Cart::yyyyz))
                                                + 30.*cou_tmp(i,Cart::xxyyz)
                                               - 40.*(cou_tmp(i,Cart::xxzzz) + cou_tmp(i,Cart::yyzzz))
                                                 + 8.*cou_tmp(i,Cart::zzzzz) );                              /// Y 5,0

                cou_sph(i,26) = factor_2 * (        cou_tmp(i,Cart::xxxxy)
                                                 + 2.*cou_tmp(i,Cart::xxyyy)
                                               - 12.*(cou_tmp(i,Cart::xxyzz) + cou_tmp(i,Cart::yyyzz))
                                                    + cou_tmp(i,Cart::yyyyy)
                                                 + 8.*cou_tmp(i,Cart::yzzzz) );                              /// Y 5,-1

                cou_sph(i,27) = factor_2 * (        cou_tmp(i,Cart::xxxxx)
                                                 + 2.*cou_tmp(i,Cart::xxxyy)
                                               - 12.*(cou_tmp(i,Cart::xxxzz) + cou_tmp(i,Cart::xyyzz))
                                                    + cou_tmp(i,Cart::xyyyy)
                                                 + 8.*cou_tmp(i,Cart::xzzzz) );                              /// Y 5,1

                cou_sph(i,28) = 8.*factor * (     -cou_tmp(i,Cart::xxxyz)
                                                   - cou_tmp(i,Cart::xyyyz)
                                                + 2.*cou_tmp(i,Cart::xyzzz) );                               /// Y 5,-2

                cou_sph(i,29) = 4.*factor * (      -cou_tmp(i,Cart::xxxxz)
                                                + 2.*(cou_tmp(i,Cart::xxzzz) - cou_tmp(i,Cart::yyzzz))
                                                    + cou_tmp(i,Cart::yyyyz) );                              /// Y 5,2

                cou_sph(i,30) = factor_3 * (   -3.*cou_tmp(i,Cart::xxxxy)
                                                - 2.*cou_tmp(i,Cart::xxyyy)
                                               + 24.*cou_tmp(i,Cart::xxyzz)
                                                   + cou_tmp(i,Cart::yyyyy)
                                                - 8.*cou_tmp(i,Cart::yyyzz) );                               /// Y 5,-3

                cou_sph(i,31) = factor_3 * (      -cou_tmp(i,Cart::xxxxx)
                                                + 2.*cou_tmp(i,Cart::xxxyy)
                                                + 8.*cou_tmp(i,Cart::xxxzz)
                                                + 3.*cou_tmp(i,Cart::xyyyy)
                                               - 24.*cou_tmp(i,Cart::xyyzz) );                               /// Y 5,3

                cou_sph(i,32) = 4.*factor_4 * (   cou_tmp(i,Cart::xxxyz)
                                                  - cou_tmp(i,Cart::xyyyz) );                                /// Y 5,-4

                cou_sph(i,33) = factor_4 * (      cou_tmp(i,Cart::xxxxz)
                                               - 6.*cou_tmp(i,Cart::xxyyz)
                                                  + cou_tmp(i,Cart::yyyyz) );                                /// Y 5,4

                cou_sph(i,34) = factor_5 * (    5.*cou_tmp(i,Cart::xxxxy)
                                               - 10.*cou_tmp(i,Cart::xxyyy)
                                                   + cou_tmp(i,Cart::yyyyy) );                               /// Y 5,-5

                cou_sph(i,35) = factor_5 * (       cou_tmp(i,Cart::xxxxx)
                                               - 10.*cou_tmp(i,Cart::xxxyy)
                                                + 5.*cou_tmp(i,Cart::xyyyy) );                               /// Y 5,5
              }
            }


            if (lmax_col > 5) {
              // i-functions
              factor = (2./3.)*decay_col*decay_col*decay_col*contractions_col[6];
              double factor_1 = factor*2./sqrt(1155.);
              double factor_2 = factor*4./sqrt(55.);
              double factor_3 = factor*sqrt(22.)/11.;
              double factor_4 = factor*2.*sqrt(165.)/55.;
              double factor_5 = factor*.4*sqrt(30.);
              double factor_6 = factor*.2*sqrt(10.);
              for (int i =  0; i < ntrafo_row; i++) {
                cou_sph(i,36) = factor_1 * (    -5.*(cou_tmp(i,Cart::xxxxxx) + cou_tmp(i,Cart::yyyyyy))
                                                - 15.*(cou_tmp(i,Cart::xxxxyy) + cou_tmp(i,Cart::xxyyyy))
                                                + 90.*(cou_tmp(i,Cart::xxxxzz) + cou_tmp(i,Cart::yyyyzz))
                                                + 180.*cou_tmp(i,Cart::xxyyzz)
                                               - 120.*(cou_tmp(i,Cart::xxzzzz) + cou_tmp(i,Cart::yyzzzz))
                                                 + 16.*cou_tmp(i,Cart::zzzzzz) );                                /// Y 6,0

                cou_sph(i,37) = factor_2 * (    5.*(cou_tmp(i,Cart::xxxxyz) + cou_tmp(i,Cart::yyyyyz))
                                                + 10.*cou_tmp(i,Cart::xxyyyz)
                                               - 20.*(cou_tmp(i,Cart::xxyzzz) + cou_tmp(i,Cart::yyyzzz))
                                                 + 8.*cou_tmp(i,Cart::yzzzzz) );                                 /// Y 6,-1

                cou_sph(i,38) = factor_2 * (    5.*(cou_tmp(i,Cart::xxxxxz) + cou_tmp(i,Cart::xyyyyz))
                                                + 10.*cou_tmp(i,Cart::xxxyyz)
                                               - 20.*(cou_tmp(i,Cart::xxxzzz) + cou_tmp(i,Cart::xyyzzz))
                                                 + 8.*cou_tmp(i,Cart::xzzzzz) );                                 /// Y 6,1

                cou_sph(i,39) = 2.*factor_3 * (        cou_tmp(i,Cart::xxxxxy)
                                                    + 2.*cou_tmp(i,Cart::xxxyyy)
                                                  - 16.*(cou_tmp(i,Cart::xxxyzz) + cou_tmp(i,Cart::xyyyzz) - cou_tmp(i,Cart::xyzzzz))
                                                       + cou_tmp(i,Cart::xyyyyy) );                              /// Y 6,-2

                cou_sph(i,40) = factor_3 * (        cou_tmp(i,Cart::xxxxxy)
                                                    + cou_tmp(i,Cart::xxxxyy)
                                               - 16.*(cou_tmp(i,Cart::xxxxzz) - cou_tmp(i,Cart::xxzzzz)
                                                                                - cou_tmp(i,Cart::yyyyzz) + cou_tmp(i,Cart::yyzzzz))
                                                    - cou_tmp(i,Cart::xxyyyy)
                                                    - cou_tmp(i,Cart::yyyyyy) );                                 /// Y 6,2

                cou_sph(i,41) = 2.*factor_3 * (   -9.*cou_tmp(i,Cart::xxxxyz)
                                                   - 6.*cou_tmp(i,Cart::xxyyyz)
                                                  + 24.*cou_tmp(i,Cart::xxyzzz)
                                                   + 3.*cou_tmp(i,Cart::yyyyyz)
                                                   - 8.*cou_tmp(i,Cart::yyyzzz) );                               /// Y 6,-3

                cou_sph(i,42) = 2.*factor_3 * (   -3.*cou_tmp(i,Cart::xxxxxz)
                                                   + 6.*cou_tmp(i,Cart::xxxyyz)
                                                   + 8.*cou_tmp(i,Cart::xxxzzz)
                                                   + 9.*cou_tmp(i,Cart::xyyyyz)
                                                  - 24.*cou_tmp(i,Cart::xyyzzz) );                               /// Y 6,3

                cou_sph(i,43) = 4.*factor_4 * (       -cou_tmp(i,Cart::xxxxxy)
                                                  + 10.*(cou_tmp(i,Cart::xxxyzz) - cou_tmp(i,Cart::xyyyzz))
                                                       + cou_tmp(i,Cart::xyyyyy) );                              /// Y 6,-4

                cou_sph(i,44) = factor_4 * (       -cou_tmp(i,Cart::xxxxxx)
                                                + 5.*(cou_tmp(i,Cart::xxxxyy) + cou_tmp(i,Cart::xxyyyy))
                                               + 10.*(cou_tmp(i,Cart::xxxxzz) + cou_tmp(i,Cart::yyyyzz))
                                                - 60.*cou_tmp(i,Cart::xxyyzz)
                                                   -  cou_tmp(i,Cart::yyyyyy) );                                 /// Y 6,4

                cou_sph(i,45) = factor_5 * (    5.*cou_tmp(i,Cart::xxxxyz)
                                               - 10.*cou_tmp(i,Cart::xxyyyz)
                                                   + cou_tmp(i,Cart::yyyyyz) );                                  /// Y 6,-5

                cou_sph(i,46) = factor_5 * (       cou_tmp(i,Cart::xxxxxz)
                                               - 10.*cou_tmp(i,Cart::xxxyyz)
                                                + 5.*cou_tmp(i,Cart::xyyyyz) );                                  /// Y 6,5

                cou_sph(i,47) = 2.*factor_6 * (    3.*cou_tmp(i,Cart::xxxxxy)
                                                  - 10.*cou_tmp(i,Cart::xxxyyy)
                                                   + 3.*cou_tmp(i,Cart::xyyyyy) );                               /// Y 6,-6

                cou_sph(i,48) = factor_6 * (        cou_tmp(i,Cart::xxxxxx)
                                               - 15.*(cou_tmp(i,Cart::xxxxyy) - cou_tmp(i,Cart::xxyyyy))
                                                    - cou_tmp(i,Cart::yyyyyy) );                                 /// Y 6,6

              }
            }



            // save to matrix
            for (unsigned i = 0; i < matrix.rows(); i++) {
                for (unsigned j = 0; j < matrix.cols(); j++) {
                    matrix(i, j) += cou_sph(i + shell_row->getOffset(), j + shell_col->getOffset());
                }
            }

                } // shell_col Gaussians
            } // shell_row Gaussians
           return; 
            }    
    

    //This converts V into ((S-1/2 V S-1/2)-1/2 S-1/2)T, which is needed to construct 4c integrals,
    Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap, double etol){
        
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eo(auxoverlap.Matrix());
      removedfunctions=0;
      Eigen::VectorXd diagonal_overlap=Eigen::VectorXd::Zero(eo.eigenvalues().size());
     for (unsigned i=0;i<diagonal_overlap.size();++i){
          if(eo.eigenvalues()(i)<etol){
              removedfunctions++;
          }else{
              diagonal_overlap(i)=1.0/std::sqrt(eo.eigenvalues()(i));
          }
      }
      Eigen::MatrixXd Ssqrt=eo.eigenvectors() * diagonal_overlap.asDiagonal() * eo.eigenvectors().transpose();

      Eigen::MatrixXd ortho=Ssqrt*_aomatrix*Ssqrt;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(ortho); 
      Eigen::VectorXd diagonal=Eigen::VectorXd::Zero(es.eigenvalues().size());
    
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

