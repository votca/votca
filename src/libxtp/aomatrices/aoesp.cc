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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/xtp/votca_xtp_config.h>

#include <votca/xtp/aomatrix.h>

#include <votca/xtp/aobasis.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>
#include <votca/xtp/logger.h>
#include <votca/tools/linalg.h>
#include <votca/xtp/elements.h>
//#include <boost/timer/timer.hpp>

using namespace std;
using namespace votca::tools;



namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
    
    

    
    void AOESP::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col , AOBasis* ecp) {
        /*cout << "\nAO block: "<< endl;
        cout << "\t row: " << _shell_row->getType() << " at " << _shell_row->getPos() << endl;
        cout << "\t col: " << _shell_col->getType() << " at " << _shell_col->getPos() << endl;*/
        const double pi = boost::math::constants::pi<double>();
       
        
        // cout << _gridpoint << endl;
        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
    
        // initialize local matrix block for unnormalized cartesians
        ub::matrix<double> nuc   = ub::zero_matrix<double>(_nrows,_ncols);
        ub::matrix<double> nucm1 = ub::zero_matrix<double>(_nrows,_ncols);
        ub::matrix<double> nucm2 = ub::zero_matrix<double>(_nrows,_ncols);
        ub::matrix<double> nucm3 = ub::zero_matrix<double>(_nrows,_ncols);
        ub::matrix<double> nucm4 = ub::zero_matrix<double>(_nrows,_ncols);

        //cout << nuc.size1() << ":" << nuc.size2() << endl;
        
        /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
         * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION 
         * COEFFICIENTS, AND ADD TO matrix(i,j)
         */
        
      
       
      
        
        // get shell positions
        const vec& _pos_row = _shell_row->getPos();
        const vec& _pos_col = _shell_col->getPos();
        const vec  _diff    = _pos_row - _pos_col;
        // initialize some helper
        vector<double> PmA (3,0.0);
        vector<double> PmB (3,0.0);
        vector<double> PmC (3,0.0);
        double _distsq = (_diff.getX()*_diff.getX()) + (_diff.getY()*_diff.getY()) + (_diff.getZ()*_diff.getZ()); 
        
         typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
        // iterate over Gaussians in this _shell_row
        for ( GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
            // get decay constant
            const double& _decay_row = (*itr)->decay;
            
            for ( GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
                //get decay constant
                const double& _decay_col = (*itc)->decay;
        
                const double _fak  = 0.5/(_decay_row + _decay_col);
                const double _fak2 = 2.0 * _fak;
                
                
                double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
        
       // check if distance between postions is big, then skip step   
       
                if ( _exparg > 30.0 ) { continue; }
        
        // some helpers
       


        const double zeta = _decay_row + _decay_col;

        PmA[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
        PmA[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
        PmA[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

        PmB[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
        PmB[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
        PmB[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
        PmC[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _gridpoint[0];
        PmC[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _gridpoint[1];
        PmC[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _gridpoint[2];
        
        
        const double _U = zeta*(PmC[0]*PmC[0]+PmC[1]*PmC[1]+PmC[2]*PmC[2]);
        
        vector<double> _FmU(5, 0.0); // that size needs to be checked!

        XIntegrate(_FmU,_U );
        //cout << endl;
        
        
        // (s-s element normiert )
        double _prefactor = 2*sqrt(1.0/pi)*pow(4.0*_decay_row*_decay_col,0.75) * _fak2 * exp(-_exparg);
        nuc(Cart::s,Cart::s)   = _prefactor * _FmU[0];
        nucm1(Cart::s,Cart::s) = _prefactor * _FmU[1];
        nucm2(Cart::s,Cart::s) = _prefactor * _FmU[2];
        nucm3(Cart::s,Cart::s) = _prefactor * _FmU[3];
        nucm4(Cart::s,Cart::s) = _prefactor * _FmU[4];
        double _lsum=_lmax_row + _lmax_col;
        
        if (_lmax_col >2 || _lmax_row >2){
            cerr << "Orbitals higher than d are not yet implemented. This should not have happened!" << flush;
             exit(1);
        }
        

        // s-p-0
        if ( _lmax_col > 0 ) {
                nuc(Cart::s,Cart::y) =PmB[1]*nuc(Cart::s,Cart::s)-PmC[1]*nucm1(Cart::s,Cart::s);
                nuc(Cart::s,Cart::x) =PmB[0]*nuc(Cart::s,Cart::s)-PmC[0]*nucm1(Cart::s,Cart::s);
                nuc(Cart::s,Cart::z) =PmB[2]*nuc(Cart::s,Cart::s)-PmC[2]*nucm1(Cart::s,Cart::s);


          // cout << "\t setting s-p" << flush;
  
        }
        
        // p-s-0
        if ( _lmax_row > 0 ) {
           //cout << "\t setting p-s" << flush;
                nuc(Cart::y,Cart::s) =PmA[1]*nuc(Cart::s,Cart::s)-PmC[1]*nucm1(Cart::s,Cart::s);
                nuc(Cart::x,Cart::s) =PmA[0]*nuc(Cart::s,Cart::s)-PmC[0]*nucm1(Cart::s,Cart::s);
                nuc(Cart::z,Cart::s) =PmA[2]*nuc(Cart::s,Cart::s)-PmC[2]*nucm1(Cart::s,Cart::s);

        }
        

        if ( _lsum > 1 ) {
          // cout << "\t setting p-p" << endl; 
            
            //m=1
            //s-p-1
           if (_lmax_col>0){
            
                nucm1(Cart::s,Cart::y) =PmB[1]*nucm1(Cart::s,Cart::s)-PmC[1]*nucm2(Cart::s,Cart::s);
                nucm1(Cart::s,Cart::x) =PmB[0]*nucm1(Cart::s,Cart::s)-PmC[0]*nucm2(Cart::s,Cart::s);
                nucm1(Cart::s,Cart::z) =PmB[2]*nucm1(Cart::s,Cart::s)-PmC[2]*nucm2(Cart::s,Cart::s);
           }
              
            // p-s-1
             if (_lmax_row>0){   
                nucm1(Cart::y,Cart::s) =PmA[1]*nucm1(Cart::s,Cart::s)-PmC[1]*nucm2(Cart::s,Cart::s);
                nucm1(Cart::x,Cart::s) =PmA[0]*nucm1(Cart::s,Cart::s)-PmC[0]*nucm2(Cart::s,Cart::s);
                nucm1(Cart::z,Cart::s) =PmA[2]*nucm1(Cart::s,Cart::s)-PmC[2]*nucm2(Cart::s,Cart::s);
            }      
        }
                
        if ( _lmax_row > 0 && _lmax_col > 0 ) {
           //cout << "\t setting p-p" << endl; 
            
            // p-p-0 
                nuc(Cart::y,Cart::y) =PmA[1]*nuc(Cart::s,Cart::y)-PmC[1]*nucm1(Cart::s,Cart::y)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));
                nuc(Cart::y,Cart::x) =PmA[1]*nuc(Cart::s,Cart::x)-PmC[1]*nucm1(Cart::s,Cart::x);
                nuc(Cart::y,Cart::z) =PmA[1]*nuc(Cart::s,Cart::z)-PmC[1]*nucm1(Cart::s,Cart::z);
                nuc(Cart::x,Cart::y) =PmA[0]*nuc(Cart::s,Cart::y)-PmC[0]*nucm1(Cart::s,Cart::y);
                nuc(Cart::x,Cart::x) =PmA[0]*nuc(Cart::s,Cart::x)-PmC[0]*nucm1(Cart::s,Cart::x)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));
                nuc(Cart::x,Cart::z) =PmA[0]*nuc(Cart::s,Cart::z)-PmC[0]*nucm1(Cart::s,Cart::z);
                nuc(Cart::z,Cart::y) =PmA[2]*nuc(Cart::s,Cart::y)-PmC[2]*nucm1(Cart::s,Cart::y);
                nuc(Cart::z,Cart::x) =PmA[2]*nuc(Cart::s,Cart::x)-PmC[2]*nucm1(Cart::s,Cart::x);
                nuc(Cart::z,Cart::z) =PmA[2]*nuc(Cart::s,Cart::z)-PmC[2]*nucm1(Cart::s,Cart::z)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));

        } 
      
        // s-d
       if ( _lmax_col > 1){
            //cout << "\t setting s-d" << endl;
          // s-d-0
                nuc(Cart::s,Cart::yy) =PmB[1]*nuc(Cart::s,Cart::y)-PmC[1]*nucm1(Cart::s,Cart::y)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));
                nuc(Cart::s,Cart::xy) =PmB[0]*nuc(Cart::s,Cart::y)-PmC[0]*nucm1(Cart::s,Cart::y);
                nuc(Cart::s,Cart::yz) =PmB[1]*nuc(Cart::s,Cart::z)-PmC[1]*nucm1(Cart::s,Cart::z);
                nuc(Cart::s,Cart::xx) =PmB[0]*nuc(Cart::s,Cart::x)-PmC[0]*nucm1(Cart::s,Cart::x)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));
                nuc(Cart::s,Cart::xz) =PmB[0]*nuc(Cart::s,Cart::z)-PmC[0]*nucm1(Cart::s,Cart::z);
                nuc(Cart::s,Cart::zz) =PmB[2]*nuc(Cart::s,Cart::z)-PmC[2]*nucm1(Cart::s,Cart::z)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));

            
        }
        
        
         // d-s
        if ( _lmax_row > 1){
           //cout << "\t setting d-s" << endl;
                nuc(Cart::yy,Cart::s) =PmA[1]*nuc(Cart::y,Cart::s)-PmC[1]*nucm1(Cart::y,Cart::s)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));
                nuc(Cart::xy,Cart::s) =PmA[0]*nuc(Cart::y,Cart::s)-PmC[0]*nucm1(Cart::y,Cart::s);
                nuc(Cart::yz,Cart::s) =PmA[1]*nuc(Cart::z,Cart::s)-PmC[1]*nucm1(Cart::z,Cart::s);
                nuc(Cart::xx,Cart::s) =PmA[0]*nuc(Cart::x,Cart::s)-PmC[0]*nucm1(Cart::x,Cart::s)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));
                nuc(Cart::xz,Cart::s) =PmA[0]*nuc(Cart::z,Cart::s)-PmC[0]*nucm1(Cart::z,Cart::s);
                nuc(Cart::zz,Cart::s) =PmA[2]*nuc(Cart::z,Cart::s)-PmC[2]*nucm1(Cart::z,Cart::s)+_fak*(nuc(Cart::s,Cart::s)-nucm1(Cart::s,Cart::s));

         
        }
        
        
        if ( _lsum > 2 ){
            //cout << "\t setting p-d" << endl;
            if ( _lmax_col > 0){
            //s-p-2
                nucm2(Cart::s,Cart::y) =PmB[1]*nucm2(Cart::s,Cart::s)-PmC[1]*nucm3(Cart::s,Cart::s);
                nucm2(Cart::s,Cart::x) =PmB[0]*nucm2(Cart::s,Cart::s)-PmC[0]*nucm3(Cart::s,Cart::s);
                nucm2(Cart::s,Cart::z) =PmB[2]*nucm2(Cart::s,Cart::s)-PmC[2]*nucm3(Cart::s,Cart::s);
            }
            if ( _lmax_row > 0){    
            //p-s-2
                nucm2(Cart::y,Cart::s) =PmA[1]*nucm2(Cart::s,Cart::s)-PmC[1]*nucm3(Cart::s,Cart::s);
                nucm2(Cart::x,Cart::s) =PmA[0]*nucm2(Cart::s,Cart::s)-PmC[0]*nucm3(Cart::s,Cart::s);
                nucm2(Cart::z,Cart::s) =PmA[2]*nucm2(Cart::s,Cart::s)-PmC[2]*nucm3(Cart::s,Cart::s);
            }

            if ( _lmax_row > 0 && _lmax_col > 0 ) {
            //p-p-1
                nucm1(Cart::y,Cart::y) =PmA[1]*nucm1(Cart::s,Cart::y)-PmC[1]*nucm2(Cart::s,Cart::y)+_fak*(nucm1(Cart::s,Cart::s)-nucm2(Cart::s,Cart::s));
                nucm1(Cart::y,Cart::x) =PmA[1]*nucm1(Cart::s,Cart::x)-PmC[1]*nucm2(Cart::s,Cart::x);
                nucm1(Cart::y,Cart::z) =PmA[1]*nucm1(Cart::s,Cart::z)-PmC[1]*nucm2(Cart::s,Cart::z);
                nucm1(Cart::x,Cart::y) =PmA[0]*nucm1(Cart::s,Cart::y)-PmC[0]*nucm2(Cart::s,Cart::y);
                nucm1(Cart::x,Cart::x) =PmA[0]*nucm1(Cart::s,Cart::x)-PmC[0]*nucm2(Cart::s,Cart::x)+_fak*(nucm1(Cart::s,Cart::s)-nucm2(Cart::s,Cart::s));
                nucm1(Cart::x,Cart::z) =PmA[0]*nucm1(Cart::s,Cart::z)-PmC[0]*nucm2(Cart::s,Cart::z);
                nucm1(Cart::z,Cart::y) =PmA[2]*nucm1(Cart::s,Cart::y)-PmC[2]*nucm2(Cart::s,Cart::y);
                nucm1(Cart::z,Cart::x) =PmA[2]*nucm1(Cart::s,Cart::x)-PmC[2]*nucm2(Cart::s,Cart::x);
                nucm1(Cart::z,Cart::z) =PmA[2]*nucm1(Cart::s,Cart::z)-PmC[2]*nucm2(Cart::s,Cart::z)+_fak*(nucm1(Cart::s,Cart::s)-nucm2(Cart::s,Cart::s));
        }
        }
        
        
        
        
         //p-d
        if ( _lmax_row > 0 && _lmax_col > 1){
            //cout << "\t setting p-d" << endl;
            
        
            // p-d-0
                nuc(Cart::y,Cart::yy) =PmB[1]*nuc(Cart::y,Cart::y)-PmC[1]*nucm1(Cart::y,Cart::y)+_fak*(nuc(Cart::y,Cart::s)-nucm1(Cart::y,Cart::s))+_fak*(nuc(Cart::s,Cart::y)-nucm1(Cart::s,Cart::y));
                nuc(Cart::y,Cart::xy) =PmB[0]*nuc(Cart::y,Cart::y)-PmC[0]*nucm1(Cart::y,Cart::y);
                nuc(Cart::y,Cart::yz) =PmB[1]*nuc(Cart::y,Cart::z)-PmC[1]*nucm1(Cart::y,Cart::z)+_fak*(nuc(Cart::s,Cart::z)-nucm1(Cart::s,Cart::z));
                nuc(Cart::y,Cart::xx) =PmB[0]*nuc(Cart::y,Cart::x)-PmC[0]*nucm1(Cart::y,Cart::x)+_fak*(nuc(Cart::y,Cart::s)-nucm1(Cart::y,Cart::s));
                nuc(Cart::y,Cart::xz) =PmB[0]*nuc(Cart::y,Cart::z)-PmC[0]*nucm1(Cart::y,Cart::z);
                nuc(Cart::y,Cart::zz) =PmB[2]*nuc(Cart::y,Cart::z)-PmC[2]*nucm1(Cart::y,Cart::z)+_fak*(nuc(Cart::y,Cart::s)-nucm1(Cart::y,Cart::s));
                nuc(Cart::x,Cart::yy) =PmB[1]*nuc(Cart::x,Cart::y)-PmC[1]*nucm1(Cart::x,Cart::y)+_fak*(nuc(Cart::x,Cart::s)-nucm1(Cart::x,Cart::s));
                nuc(Cart::x,Cart::xy) =PmB[0]*nuc(Cart::x,Cart::y)-PmC[0]*nucm1(Cart::x,Cart::y)+_fak*(nuc(Cart::s,Cart::y)-nucm1(Cart::s,Cart::y));
                nuc(Cart::x,Cart::yz) =PmB[1]*nuc(Cart::x,Cart::z)-PmC[1]*nucm1(Cart::x,Cart::z);
                nuc(Cart::x,Cart::xx) =PmB[0]*nuc(Cart::x,Cart::x)-PmC[0]*nucm1(Cart::x,Cart::x)+_fak*(nuc(Cart::x,Cart::s)-nucm1(Cart::x,Cart::s))+_fak*(nuc(Cart::s,Cart::x)-nucm1(Cart::s,Cart::x));
                nuc(Cart::x,Cart::xz) =PmB[0]*nuc(Cart::x,Cart::z)-PmC[0]*nucm1(Cart::x,Cart::z)+_fak*(nuc(Cart::s,Cart::z)-nucm1(Cart::s,Cart::z));
                nuc(Cart::x,Cart::zz) =PmB[2]*nuc(Cart::x,Cart::z)-PmC[2]*nucm1(Cart::x,Cart::z)+_fak*(nuc(Cart::x,Cart::s)-nucm1(Cart::x,Cart::s));
                nuc(Cart::z,Cart::yy) =PmB[1]*nuc(Cart::z,Cart::y)-PmC[1]*nucm1(Cart::z,Cart::y)+_fak*(nuc(Cart::z,Cart::s)-nucm1(Cart::z,Cart::s));
                nuc(Cart::z,Cart::xy) =PmB[0]*nuc(Cart::z,Cart::y)-PmC[0]*nucm1(Cart::z,Cart::y);
                nuc(Cart::z,Cart::yz) =PmB[1]*nuc(Cart::z,Cart::z)-PmC[1]*nucm1(Cart::z,Cart::z);
                nuc(Cart::z,Cart::xx) =PmB[0]*nuc(Cart::z,Cart::x)-PmC[0]*nucm1(Cart::z,Cart::x)+_fak*(nuc(Cart::z,Cart::s)-nucm1(Cart::z,Cart::s));
                nuc(Cart::z,Cart::xz) =PmB[0]*nuc(Cart::z,Cart::z)-PmC[0]*nucm1(Cart::z,Cart::z);
                nuc(Cart::z,Cart::zz) =PmB[2]*nuc(Cart::z,Cart::z)-PmC[2]*nucm1(Cart::z,Cart::z)+_fak*(nuc(Cart::z,Cart::s)-nucm1(Cart::z,Cart::s))+_fak*(nuc(Cart::s,Cart::z)-nucm1(Cart::s,Cart::z));

         
        }

       
        
        
        // d-p
        if ( _lmax_row >1 && _lmax_col > 0){
           //cout << "\t setting d-p" << endl;
            
                nuc(Cart::yy,Cart::y) =PmA[1]*nuc(Cart::y,Cart::y)-PmC[1]*nucm1(Cart::y,Cart::y)+_fak*(nuc(Cart::s,Cart::y)-nucm1(Cart::s,Cart::y))+_fak*(nuc(Cart::y,Cart::s)-nucm1(Cart::y,Cart::s));
                nuc(Cart::yy,Cart::x) =PmA[1]*nuc(Cart::y,Cart::x)-PmC[1]*nucm1(Cart::y,Cart::x)+_fak*(nuc(Cart::s,Cart::x)-nucm1(Cart::s,Cart::x));
                nuc(Cart::yy,Cart::z) =PmA[1]*nuc(Cart::y,Cart::z)-PmC[1]*nucm1(Cart::y,Cart::z)+_fak*(nuc(Cart::s,Cart::z)-nucm1(Cart::s,Cart::z));
                nuc(Cart::xy,Cart::y) =PmA[0]*nuc(Cart::y,Cart::y)-PmC[0]*nucm1(Cart::y,Cart::y);
                nuc(Cart::xy,Cart::x) =PmA[0]*nuc(Cart::y,Cart::x)-PmC[0]*nucm1(Cart::y,Cart::x)+_fak*(nuc(Cart::y,Cart::s)-nucm1(Cart::y,Cart::s));
                nuc(Cart::xy,Cart::z) =PmA[0]*nuc(Cart::y,Cart::z)-PmC[0]*nucm1(Cart::y,Cart::z);
                nuc(Cart::yz,Cart::y) =PmA[1]*nuc(Cart::z,Cart::y)-PmC[1]*nucm1(Cart::z,Cart::y)+_fak*(nuc(Cart::z,Cart::s)-nucm1(Cart::z,Cart::s));
                nuc(Cart::yz,Cart::x) =PmA[1]*nuc(Cart::z,Cart::x)-PmC[1]*nucm1(Cart::z,Cart::x);
                nuc(Cart::yz,Cart::z) =PmA[1]*nuc(Cart::z,Cart::z)-PmC[1]*nucm1(Cart::z,Cart::z);
                nuc(Cart::xx,Cart::y) =PmA[0]*nuc(Cart::x,Cart::y)-PmC[0]*nucm1(Cart::x,Cart::y)+_fak*(nuc(Cart::s,Cart::y)-nucm1(Cart::s,Cart::y));
                nuc(Cart::xx,Cart::x) =PmA[0]*nuc(Cart::x,Cart::x)-PmC[0]*nucm1(Cart::x,Cart::x)+_fak*(nuc(Cart::s,Cart::x)-nucm1(Cart::s,Cart::x))+_fak*(nuc(Cart::x,Cart::s)-nucm1(Cart::x,Cart::s));
                nuc(Cart::xx,Cart::z) =PmA[0]*nuc(Cart::x,Cart::z)-PmC[0]*nucm1(Cart::x,Cart::z)+_fak*(nuc(Cart::s,Cart::z)-nucm1(Cart::s,Cart::z));
                nuc(Cart::xz,Cart::y) =PmA[0]*nuc(Cart::z,Cart::y)-PmC[0]*nucm1(Cart::z,Cart::y);
                nuc(Cart::xz,Cart::x) =PmA[0]*nuc(Cart::z,Cart::x)-PmC[0]*nucm1(Cart::z,Cart::x)+_fak*(nuc(Cart::z,Cart::s)-nucm1(Cart::z,Cart::s));
                nuc(Cart::xz,Cart::z) =PmA[0]*nuc(Cart::z,Cart::z)-PmC[0]*nucm1(Cart::z,Cart::z);
                nuc(Cart::zz,Cart::y) =PmA[2]*nuc(Cart::z,Cart::y)-PmC[2]*nucm1(Cart::z,Cart::y)+_fak*(nuc(Cart::s,Cart::y)-nucm1(Cart::s,Cart::y));
                nuc(Cart::zz,Cart::x) =PmA[2]*nuc(Cart::z,Cart::x)-PmC[2]*nucm1(Cart::z,Cart::x)+_fak*(nuc(Cart::s,Cart::x)-nucm1(Cart::s,Cart::x));
                nuc(Cart::zz,Cart::z) =PmA[2]*nuc(Cart::z,Cart::z)-PmC[2]*nucm1(Cart::z,Cart::z)+_fak*(nuc(Cart::s,Cart::z)-nucm1(Cart::s,Cart::z))+_fak*(nuc(Cart::z,Cart::s)-nucm1(Cart::z,Cart::s));


        }
        
         
        if ( _lsum > 2 && _lmax_col>1 && _lmax_row>1){
            
            
           if (  _lmax_col > 0){
            //s-p-3
            
            nucm3(Cart::s, Cart::y) = PmB[1] * nucm3(Cart::s, Cart::s) - PmC[1] * nucm4(Cart::s, Cart::s);
            nucm3(Cart::s, Cart::x) = PmB[0] * nucm3(Cart::s, Cart::s) - PmC[0] * nucm4(Cart::s, Cart::s);
            nucm3(Cart::s, Cart::z) = PmB[2] * nucm3(Cart::s, Cart::s) - PmC[2] * nucm4(Cart::s, Cart::s);
            }
            
             if ( _lmax_row > 0 && _lmax_col > 0 ) {
            //p-p-2
            
            nucm2(Cart::y,Cart::y) =PmA[1]*nucm2(Cart::s,Cart::y)-PmC[1]*nucm3(Cart::s,Cart::y)+_fak*(nucm2(Cart::s,Cart::s)-nucm3(Cart::s,Cart::s));
            nucm2(Cart::y,Cart::x) =PmA[1]*nucm2(Cart::s,Cart::x)-PmC[1]*nucm3(Cart::s,Cart::x);
            nucm2(Cart::y,Cart::z) =PmA[1]*nucm2(Cart::s,Cart::z)-PmC[1]*nucm3(Cart::s,Cart::z);
            nucm2(Cart::x,Cart::y) =PmA[0]*nucm2(Cart::s,Cart::y)-PmC[0]*nucm3(Cart::s,Cart::y);
            nucm2(Cart::x,Cart::x) =PmA[0]*nucm2(Cart::s,Cart::x)-PmC[0]*nucm3(Cart::s,Cart::x)+_fak*(nucm2(Cart::s,Cart::s)-nucm3(Cart::s,Cart::s));
            nucm2(Cart::x,Cart::z) =PmA[0]*nucm2(Cart::s,Cart::z)-PmC[0]*nucm3(Cart::s,Cart::z);
            nucm2(Cart::z,Cart::y) =PmA[2]*nucm2(Cart::s,Cart::y)-PmC[2]*nucm3(Cart::s,Cart::y);
            nucm2(Cart::z,Cart::x) =PmA[2]*nucm2(Cart::s,Cart::x)-PmC[2]*nucm3(Cart::s,Cart::x);
            nucm2(Cart::z,Cart::z) =PmA[2]*nucm2(Cart::s,Cart::z)-PmC[2]*nucm3(Cart::s,Cart::z)+_fak*(nucm2(Cart::s,Cart::s)-nucm3(Cart::s,Cart::s));
             }
            
            
              if ( _lmax_row > 0 && _lmax_col > 1){
                  //s-d-1
             
            nucm1(Cart::s,Cart::yy) =PmB[1]*nucm1(Cart::s,Cart::y)-PmC[1]*nucm2(Cart::s,Cart::y)+_fak*(nucm1(Cart::s,Cart::s)-nucm2(Cart::s,Cart::s));
            nucm1(Cart::s,Cart::xy) =PmB[0]*nucm1(Cart::s,Cart::y)-PmC[0]*nucm2(Cart::s,Cart::y);
            nucm1(Cart::s,Cart::yz) =PmB[1]*nucm1(Cart::s,Cart::z)-PmC[1]*nucm2(Cart::s,Cart::z);
            nucm1(Cart::s,Cart::xx) =PmB[0]*nucm1(Cart::s,Cart::x)-PmC[0]*nucm2(Cart::s,Cart::x)+_fak*(nucm1(Cart::s,Cart::s)-nucm2(Cart::s,Cart::s));
            nucm1(Cart::s,Cart::xz) =PmB[0]*nucm1(Cart::s,Cart::z)-PmC[0]*nucm2(Cart::s,Cart::z);
            nucm1(Cart::s,Cart::zz) =PmB[2]*nucm1(Cart::s,Cart::z)-PmC[2]*nucm2(Cart::s,Cart::z)+_fak*(nucm1(Cart::s,Cart::s)-nucm2(Cart::s,Cart::s));
        }
            
             if ( _lmax_row > 0 && _lmax_col > 1){
            //p-d-1
            nucm1(Cart::y,Cart::yy) =PmB[1]*nucm1(Cart::y,Cart::y)-PmC[1]*nucm2(Cart::y,Cart::y)+_fak*(nucm1(Cart::y,Cart::s)-nucm2(Cart::y,Cart::s))+_fak*(nucm1(Cart::s,Cart::y)-nucm2(Cart::s,Cart::y));
            nucm1(Cart::y,Cart::xy) =PmB[0]*nucm1(Cart::y,Cart::y)-PmC[0]*nucm2(Cart::y,Cart::y);
            nucm1(Cart::y,Cart::yz) =PmB[1]*nucm1(Cart::y,Cart::z)-PmC[1]*nucm2(Cart::y,Cart::z)+_fak*(nucm1(Cart::s,Cart::z)-nucm2(Cart::s,Cart::z));
            nucm1(Cart::y,Cart::xx) =PmB[0]*nucm1(Cart::y,Cart::x)-PmC[0]*nucm2(Cart::y,Cart::x)+_fak*(nucm1(Cart::y,Cart::s)-nucm2(Cart::y,Cart::s));
            nucm1(Cart::y,Cart::xz) =PmB[0]*nucm1(Cart::y,Cart::z)-PmC[0]*nucm2(Cart::y,Cart::z);
            nucm1(Cart::y,Cart::zz) =PmB[2]*nucm1(Cart::y,Cart::z)-PmC[2]*nucm2(Cart::y,Cart::z)+_fak*(nucm1(Cart::y,Cart::s)-nucm2(Cart::y,Cart::s));
            nucm1(Cart::x,Cart::yy) =PmB[1]*nucm1(Cart::x,Cart::y)-PmC[1]*nucm2(Cart::x,Cart::y)+_fak*(nucm1(Cart::x,Cart::s)-nucm2(Cart::x,Cart::s));
            nucm1(Cart::x,Cart::xy) =PmB[0]*nucm1(Cart::x,Cart::y)-PmC[0]*nucm2(Cart::x,Cart::y)+_fak*(nucm1(Cart::s,Cart::y)-nucm2(Cart::s,Cart::y));
            nucm1(Cart::x,Cart::yz) =PmB[1]*nucm1(Cart::x,Cart::z)-PmC[1]*nucm2(Cart::x,Cart::z);
            nucm1(Cart::x,Cart::xx) =PmB[0]*nucm1(Cart::x,Cart::x)-PmC[0]*nucm2(Cart::x,Cart::x)+_fak*(nucm1(Cart::x,Cart::s)-nucm2(Cart::x,Cart::s))+_fak*(nucm1(Cart::s,Cart::x)-nucm2(Cart::s,Cart::x));
            nucm1(Cart::x,Cart::xz) =PmB[0]*nucm1(Cart::x,Cart::z)-PmC[0]*nucm2(Cart::x,Cart::z)+_fak*(nucm1(Cart::s,Cart::z)-nucm2(Cart::s,Cart::z));
            nucm1(Cart::x,Cart::zz) =PmB[2]*nucm1(Cart::x,Cart::z)-PmC[2]*nucm2(Cart::x,Cart::z)+_fak*(nucm1(Cart::x,Cart::s)-nucm2(Cart::x,Cart::s));
            nucm1(Cart::z,Cart::yy) =PmB[1]*nucm1(Cart::z,Cart::y)-PmC[1]*nucm2(Cart::z,Cart::y)+_fak*(nucm1(Cart::z,Cart::s)-nucm2(Cart::z,Cart::s));
            nucm1(Cart::z,Cart::xy) =PmB[0]*nucm1(Cart::z,Cart::y)-PmC[0]*nucm2(Cart::z,Cart::y);
            nucm1(Cart::z,Cart::yz) =PmB[1]*nucm1(Cart::z,Cart::z)-PmC[1]*nucm2(Cart::z,Cart::z);
            nucm1(Cart::z,Cart::xx) =PmB[0]*nucm1(Cart::z,Cart::x)-PmC[0]*nucm2(Cart::z,Cart::x)+_fak*(nucm1(Cart::z,Cart::s)-nucm2(Cart::z,Cart::s));
            nucm1(Cart::z,Cart::xz) =PmB[0]*nucm1(Cart::z,Cart::z)-PmC[0]*nucm2(Cart::z,Cart::z);
            nucm1(Cart::z,Cart::zz) =PmB[2]*nucm1(Cart::z,Cart::z)-PmC[2]*nucm2(Cart::z,Cart::z)+_fak*(nucm1(Cart::z,Cart::s)-nucm2(Cart::z,Cart::s))+_fak*(nucm1(Cart::s,Cart::z)-nucm2(Cart::s,Cart::z));
             }
            
        }
     
            
            
            
            
            
            
        
        
        // d-d
        if ( _lmax_row > 1 && _lmax_col > 1 ){
             // cout << "\t setting d-d" << endl;
 
            
            //d-d-0
            nuc(Cart::yy,Cart::yy) =PmA[1]*nuc(Cart::y,Cart::yy)-PmC[1]*nucm1(Cart::y,Cart::yy)+_fak*(nuc(Cart::s,Cart::yy)-nucm1(Cart::s,Cart::yy))+_fak2*(nuc(Cart::y,Cart::y)-nucm1(Cart::y,Cart::y));
            nuc(Cart::yy,Cart::xy) =PmA[1]*nuc(Cart::y,Cart::xy)-PmC[1]*nucm1(Cart::y,Cart::xy)+_fak*(nuc(Cart::s,Cart::xy)-nucm1(Cart::s,Cart::xy))+_fak*(nuc(Cart::y,Cart::x)-nucm1(Cart::y,Cart::x));
            nuc(Cart::yy,Cart::yz) =PmA[1]*nuc(Cart::y,Cart::yz)-PmC[1]*nucm1(Cart::y,Cart::yz)+_fak*(nuc(Cart::s,Cart::yz)-nucm1(Cart::s,Cart::yz))+_fak*(nuc(Cart::y,Cart::z)-nucm1(Cart::y,Cart::z));
            nuc(Cart::yy,Cart::xx) =PmA[1]*nuc(Cart::y,Cart::xx)-PmC[1]*nucm1(Cart::y,Cart::xx)+_fak*(nuc(Cart::s,Cart::xx)-nucm1(Cart::s,Cart::xx));
            nuc(Cart::yy,Cart::xz) =PmA[1]*nuc(Cart::y,Cart::xz)-PmC[1]*nucm1(Cart::y,Cart::xz)+_fak*(nuc(Cart::s,Cart::xz)-nucm1(Cart::s,Cart::xz));
            nuc(Cart::yy,Cart::zz) =PmA[1]*nuc(Cart::y,Cart::zz)-PmC[1]*nucm1(Cart::y,Cart::zz)+_fak*(nuc(Cart::s,Cart::zz)-nucm1(Cart::s,Cart::zz));
            nuc(Cart::xy,Cart::yy) =PmA[0]*nuc(Cart::y,Cart::yy)-PmC[0]*nucm1(Cart::y,Cart::yy);
            nuc(Cart::xy,Cart::xy) =PmA[0]*nuc(Cart::y,Cart::xy)-PmC[0]*nucm1(Cart::y,Cart::xy)+_fak*(nuc(Cart::y,Cart::y)-nucm1(Cart::y,Cart::y));
            nuc(Cart::xy,Cart::yz) =PmA[0]*nuc(Cart::y,Cart::yz)-PmC[0]*nucm1(Cart::y,Cart::yz);
            nuc(Cart::xy,Cart::xx) =PmA[0]*nuc(Cart::y,Cart::xx)-PmC[0]*nucm1(Cart::y,Cart::xx)+_fak2*(nuc(Cart::y,Cart::x)-nucm1(Cart::y,Cart::x));
            nuc(Cart::xy,Cart::xz) =PmA[0]*nuc(Cart::y,Cart::xz)-PmC[0]*nucm1(Cart::y,Cart::xz)+_fak*(nuc(Cart::y,Cart::z)-nucm1(Cart::y,Cart::z));
            nuc(Cart::xy,Cart::zz) =PmA[0]*nuc(Cart::y,Cart::zz)-PmC[0]*nucm1(Cart::y,Cart::zz);
            nuc(Cart::yz,Cart::yy) =PmA[1]*nuc(Cart::z,Cart::yy)-PmC[1]*nucm1(Cart::z,Cart::yy)+_fak2*(nuc(Cart::z,Cart::y)-nucm1(Cart::z,Cart::y));
            nuc(Cart::yz,Cart::xy) =PmA[1]*nuc(Cart::z,Cart::xy)-PmC[1]*nucm1(Cart::z,Cart::xy)+_fak*(nuc(Cart::z,Cart::x)-nucm1(Cart::z,Cart::x));
            nuc(Cart::yz,Cart::yz) =PmA[1]*nuc(Cart::z,Cart::yz)-PmC[1]*nucm1(Cart::z,Cart::yz)+_fak*(nuc(Cart::z,Cart::z)-nucm1(Cart::z,Cart::z));
            nuc(Cart::yz,Cart::xx) =PmA[1]*nuc(Cart::z,Cart::xx)-PmC[1]*nucm1(Cart::z,Cart::xx);
            nuc(Cart::yz,Cart::xz) =PmA[1]*nuc(Cart::z,Cart::xz)-PmC[1]*nucm1(Cart::z,Cart::xz);
            nuc(Cart::yz,Cart::zz) =PmA[1]*nuc(Cart::z,Cart::zz)-PmC[1]*nucm1(Cart::z,Cart::zz);
            nuc(Cart::xx,Cart::yy) =PmA[0]*nuc(Cart::x,Cart::yy)-PmC[0]*nucm1(Cart::x,Cart::yy)+_fak*(nuc(Cart::s,Cart::yy)-nucm1(Cart::s,Cart::yy));
            nuc(Cart::xx,Cart::xy) =PmA[0]*nuc(Cart::x,Cart::xy)-PmC[0]*nucm1(Cart::x,Cart::xy)+_fak*(nuc(Cart::s,Cart::xy)-nucm1(Cart::s,Cart::xy))+_fak*(nuc(Cart::x,Cart::y)-nucm1(Cart::x,Cart::y));
            nuc(Cart::xx,Cart::yz) =PmA[0]*nuc(Cart::x,Cart::yz)-PmC[0]*nucm1(Cart::x,Cart::yz)+_fak*(nuc(Cart::s,Cart::yz)-nucm1(Cart::s,Cart::yz));
            nuc(Cart::xx,Cart::xx) =PmA[0]*nuc(Cart::x,Cart::xx)-PmC[0]*nucm1(Cart::x,Cart::xx)+_fak*(nuc(Cart::s,Cart::xx)-nucm1(Cart::s,Cart::xx))+_fak2*(nuc(Cart::x,Cart::x)-nucm1(Cart::x,Cart::x));
            nuc(Cart::xx,Cart::xz) =PmA[0]*nuc(Cart::x,Cart::xz)-PmC[0]*nucm1(Cart::x,Cart::xz)+_fak*(nuc(Cart::s,Cart::xz)-nucm1(Cart::s,Cart::xz))+_fak*(nuc(Cart::x,Cart::z)-nucm1(Cart::x,Cart::z));
            nuc(Cart::xx,Cart::zz) =PmA[0]*nuc(Cart::x,Cart::zz)-PmC[0]*nucm1(Cart::x,Cart::zz)+_fak*(nuc(Cart::s,Cart::zz)-nucm1(Cart::s,Cart::zz));
            nuc(Cart::xz,Cart::yy) =PmA[0]*nuc(Cart::z,Cart::yy)-PmC[0]*nucm1(Cart::z,Cart::yy);
            nuc(Cart::xz,Cart::xy) =PmA[0]*nuc(Cart::z,Cart::xy)-PmC[0]*nucm1(Cart::z,Cart::xy)+_fak*(nuc(Cart::z,Cart::y)-nucm1(Cart::z,Cart::y));
            nuc(Cart::xz,Cart::yz) =PmA[0]*nuc(Cart::z,Cart::yz)-PmC[0]*nucm1(Cart::z,Cart::yz);
            nuc(Cart::xz,Cart::xx) =PmA[0]*nuc(Cart::z,Cart::xx)-PmC[0]*nucm1(Cart::z,Cart::xx)+_fak2*(nuc(Cart::z,Cart::x)-nucm1(Cart::z,Cart::x));
            nuc(Cart::xz,Cart::xz) =PmA[0]*nuc(Cart::z,Cart::xz)-PmC[0]*nucm1(Cart::z,Cart::xz)+_fak*(nuc(Cart::z,Cart::z)-nucm1(Cart::z,Cart::z));
            nuc(Cart::xz,Cart::zz) =PmA[0]*nuc(Cart::z,Cart::zz)-PmC[0]*nucm1(Cart::z,Cart::zz);
            nuc(Cart::zz,Cart::yy) =PmA[2]*nuc(Cart::z,Cart::yy)-PmC[2]*nucm1(Cart::z,Cart::yy)+_fak*(nuc(Cart::s,Cart::yy)-nucm1(Cart::s,Cart::yy));
            nuc(Cart::zz,Cart::xy) =PmA[2]*nuc(Cart::z,Cart::xy)-PmC[2]*nucm1(Cart::z,Cart::xy)+_fak*(nuc(Cart::s,Cart::xy)-nucm1(Cart::s,Cart::xy));
            nuc(Cart::zz,Cart::yz) =PmA[2]*nuc(Cart::z,Cart::yz)-PmC[2]*nucm1(Cart::z,Cart::yz)+_fak*(nuc(Cart::s,Cart::yz)-nucm1(Cart::s,Cart::yz))+_fak*(nuc(Cart::z,Cart::y)-nucm1(Cart::z,Cart::y));
            nuc(Cart::zz,Cart::xx) =PmA[2]*nuc(Cart::z,Cart::xx)-PmC[2]*nucm1(Cart::z,Cart::xx)+_fak*(nuc(Cart::s,Cart::xx)-nucm1(Cart::s,Cart::xx));
            nuc(Cart::zz,Cart::xz) =PmA[2]*nuc(Cart::z,Cart::xz)-PmC[2]*nucm1(Cart::z,Cart::xz)+_fak*(nuc(Cart::s,Cart::xz)-nucm1(Cart::s,Cart::xz))+_fak*(nuc(Cart::z,Cart::x)-nucm1(Cart::z,Cart::x));
            nuc(Cart::zz,Cart::zz) =PmA[2]*nuc(Cart::z,Cart::zz)-PmC[2]*nucm1(Cart::z,Cart::zz)+_fak*(nuc(Cart::s,Cart::zz)-nucm1(Cart::s,Cart::zz))+_fak2*(nuc(Cart::z,Cart::z)-nucm1(Cart::z,Cart::z));


            
        }
        
       // boost::timer::cpu_times t11 = cpu_t.elapsed();
        
        //cout << "Done with unnormalized matrix " << endl;
        
        // normalization and cartesian -> spherical factors
        int _ntrafo_row = _shell_row->getNumFunc() + _shell_row->getOffset();
        int _ntrafo_col = _shell_col->getNumFunc() + _shell_col->getOffset();
        
        //cout << " _ntrafo_row " << _ntrafo_row << ":" << _shell_row->getType() << endl;
        //cout << " _ntrafo_col " << _ntrafo_col << ":" << _shell_col->getType() << endl;
        ub::matrix<double> _trafo_row = ub::zero_matrix<double>(_ntrafo_row,_nrows);
        ub::matrix<double> _trafo_col = ub::zero_matrix<double>(_ntrafo_col,_ncols);

        // get transformation matrices including contraction coefficients
        std::vector<double> _contractions_row = (*itr)->contraction;
        std::vector<double> _contractions_col = (*itc)->contraction;
        this->getTrafo( _trafo_row, _lmax_row, _decay_row, _contractions_row);
        this->getTrafo( _trafo_col, _lmax_col, _decay_col, _contractions_col);
        

        // cartesian -> spherical
             
        ub::matrix<double> _nuc_tmp = ub::prod( _trafo_row, nuc );
        ub::matrix<double> _trafo_col_tposed = ub::trans( _trafo_col );
        ub::matrix<double> _nuc_sph = ub::prod( _nuc_tmp, _trafo_col_tposed );
        // save to _matrix
        // if (_lmax_row > 1 || _lmax_col > 1 ){
        //_matrix = ub::project(_nuc_sph, ub::range(_shell_row->getOffset(), _matrix.size1() + 1), ub::range(_shell_col->getOffset(), _matrix.size2()));
        //}
        //else {
        for ( unsigned i = 0; i< _matrix.size1(); i++ ) {
            for (unsigned j = 0; j < _matrix.size2(); j++){
                _matrix(i,j) += _nuc_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
            }
        }
        
        //}
        //nuc.clear();
            }// _shell_col Gaussians
        }// _shell_row Gaussians
    }
    
  // Calculates the electrostatic potential matrix element between two basis functions, for an array of atomcores.
    void AOESP::Fillnucpotential( AOBasis* aobasis, std::vector<QMAtom*>& _atoms){
    Elements _elements;
    _nuclearpotential=ub::zero_matrix<double>(aobasis->AOBasisSize(),aobasis->AOBasisSize());
    ub::vector<double> positionofatom=ub::zero_vector<double>(3);
   for ( unsigned j = 0; j < _atoms.size(); j++){

            
            positionofatom(0) = _atoms[j]->x*1.8897259886;
            positionofatom(1) = _atoms[j]->y*1.8897259886;
            positionofatom(2) = _atoms[j]->z*1.8897259886;
             //cout << "NUC POS" << positionofatom(0) << " " << positionofatom(1) << " " << positionofatom(2) << " " << endl;
	    double Znuc = _elements.getNucCrg(_atoms[j]->type);
            //cout << "NUCLEAR CHARGE" << Znuc << endl;
            _aomatrix = ub::zero_matrix<double>( aobasis->AOBasisSize(),aobasis->AOBasisSize() );
            Fill(aobasis,positionofatom);
            //Print("TMAT");
            double Zecp =2 ;
            cout << " Warning, Zecp set to " << Zecp << endl;
            _nuclearpotential+=(-1)*(Znuc-Zecp)*_aomatrix;
           // cout << "nucpotential(0,0) " << _nuclearpotential(0,0)<< endl;
    
    }
    
    }    
    
}}

