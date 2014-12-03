/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICEN_olE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 * _olee the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include <votca/ctp/aomatrix.h>

#include <votca/ctp/aobasis.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>
#include <votca/ctp/logger.h>
#include <votca/tools/linalg.h>
#include <votca/ctp/elements.h>
//#include <boost/timer/timer.hpp>

using namespace std;
using namespace votca::tools;



namespace votca { namespace ctp {
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
        nuc(Cartesian::s,Cartesian::s)   = _prefactor * _FmU[0];
        nucm1(Cartesian::s,Cartesian::s) = _prefactor * _FmU[1];
        nucm2(Cartesian::s,Cartesian::s) = _prefactor * _FmU[2];
        nucm3(Cartesian::s,Cartesian::s) = _prefactor * _FmU[3];
        nucm4(Cartesian::s,Cartesian::s) = _prefactor * _FmU[4];
        double _lsum=_lmax_row + _lmax_col;
        
        if (_lmax_col >2 || _lmax_row >2){
            cerr << "Orbitals higher than d are not yet implemented. This should not have happened!" << flush;
             exit(1);
        }
        

        // s-p-0
        if ( _lmax_col > 0 ) {
                nuc(Cartesian::s,Cartesian::y) =PmB[1]*nuc(Cartesian::s,Cartesian::s)-PmC[1]*nucm1(Cartesian::s,Cartesian::s);
                nuc(Cartesian::s,Cartesian::x) =PmB[0]*nuc(Cartesian::s,Cartesian::s)-PmC[0]*nucm1(Cartesian::s,Cartesian::s);
                nuc(Cartesian::s,Cartesian::z) =PmB[2]*nuc(Cartesian::s,Cartesian::s)-PmC[2]*nucm1(Cartesian::s,Cartesian::s);


          // cout << "\t setting s-p" << flush;
  
        }
        
        // p-s-0
        if ( _lmax_row > 0 ) {
           //cout << "\t setting p-s" << flush;
                nuc(Cartesian::y,Cartesian::s) =PmA[1]*nuc(Cartesian::s,Cartesian::s)-PmC[1]*nucm1(Cartesian::s,Cartesian::s);
                nuc(Cartesian::x,Cartesian::s) =PmA[0]*nuc(Cartesian::s,Cartesian::s)-PmC[0]*nucm1(Cartesian::s,Cartesian::s);
                nuc(Cartesian::z,Cartesian::s) =PmA[2]*nuc(Cartesian::s,Cartesian::s)-PmC[2]*nucm1(Cartesian::s,Cartesian::s);

        }
        

        if ( _lsum > 1 ) {
          // cout << "\t setting p-p" << endl; 
            
            //m=1
            //s-p-1
           if (_lmax_col>0){
            
                nucm1(Cartesian::s,Cartesian::y) =PmB[1]*nucm1(Cartesian::s,Cartesian::s)-PmC[1]*nucm2(Cartesian::s,Cartesian::s);
                nucm1(Cartesian::s,Cartesian::x) =PmB[0]*nucm1(Cartesian::s,Cartesian::s)-PmC[0]*nucm2(Cartesian::s,Cartesian::s);
                nucm1(Cartesian::s,Cartesian::z) =PmB[2]*nucm1(Cartesian::s,Cartesian::s)-PmC[2]*nucm2(Cartesian::s,Cartesian::s);
           }
              
            // p-s-1
             if (_lmax_row>0){   
                nucm1(Cartesian::y,Cartesian::s) =PmA[1]*nucm1(Cartesian::s,Cartesian::s)-PmC[1]*nucm2(Cartesian::s,Cartesian::s);
                nucm1(Cartesian::x,Cartesian::s) =PmA[0]*nucm1(Cartesian::s,Cartesian::s)-PmC[0]*nucm2(Cartesian::s,Cartesian::s);
                nucm1(Cartesian::z,Cartesian::s) =PmA[2]*nucm1(Cartesian::s,Cartesian::s)-PmC[2]*nucm2(Cartesian::s,Cartesian::s);
            }      
        }
                
        if ( _lmax_row > 0 && _lmax_col > 0 ) {
           //cout << "\t setting p-p" << endl; 
            
            // p-p-0 
                nuc(Cartesian::y,Cartesian::y) =PmA[1]*nuc(Cartesian::s,Cartesian::y)-PmC[1]*nucm1(Cartesian::s,Cartesian::y)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));
                nuc(Cartesian::y,Cartesian::x) =PmA[1]*nuc(Cartesian::s,Cartesian::x)-PmC[1]*nucm1(Cartesian::s,Cartesian::x);
                nuc(Cartesian::y,Cartesian::z) =PmA[1]*nuc(Cartesian::s,Cartesian::z)-PmC[1]*nucm1(Cartesian::s,Cartesian::z);
                nuc(Cartesian::x,Cartesian::y) =PmA[0]*nuc(Cartesian::s,Cartesian::y)-PmC[0]*nucm1(Cartesian::s,Cartesian::y);
                nuc(Cartesian::x,Cartesian::x) =PmA[0]*nuc(Cartesian::s,Cartesian::x)-PmC[0]*nucm1(Cartesian::s,Cartesian::x)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));
                nuc(Cartesian::x,Cartesian::z) =PmA[0]*nuc(Cartesian::s,Cartesian::z)-PmC[0]*nucm1(Cartesian::s,Cartesian::z);
                nuc(Cartesian::z,Cartesian::y) =PmA[2]*nuc(Cartesian::s,Cartesian::y)-PmC[2]*nucm1(Cartesian::s,Cartesian::y);
                nuc(Cartesian::z,Cartesian::x) =PmA[2]*nuc(Cartesian::s,Cartesian::x)-PmC[2]*nucm1(Cartesian::s,Cartesian::x);
                nuc(Cartesian::z,Cartesian::z) =PmA[2]*nuc(Cartesian::s,Cartesian::z)-PmC[2]*nucm1(Cartesian::s,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));

        } 
      
        // s-d
       if ( _lmax_col > 1){
            //cout << "\t setting s-d" << endl;
          // s-d-0
                nuc(Cartesian::s,Cartesian::yy) =PmB[1]*nuc(Cartesian::s,Cartesian::y)-PmC[1]*nucm1(Cartesian::s,Cartesian::y)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));
                nuc(Cartesian::s,Cartesian::xy) =PmB[0]*nuc(Cartesian::s,Cartesian::y)-PmC[0]*nucm1(Cartesian::s,Cartesian::y);
                nuc(Cartesian::s,Cartesian::yz) =PmB[1]*nuc(Cartesian::s,Cartesian::z)-PmC[1]*nucm1(Cartesian::s,Cartesian::z);
                nuc(Cartesian::s,Cartesian::xx) =PmB[0]*nuc(Cartesian::s,Cartesian::x)-PmC[0]*nucm1(Cartesian::s,Cartesian::x)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));
                nuc(Cartesian::s,Cartesian::xz) =PmB[0]*nuc(Cartesian::s,Cartesian::z)-PmC[0]*nucm1(Cartesian::s,Cartesian::z);
                nuc(Cartesian::s,Cartesian::zz) =PmB[2]*nuc(Cartesian::s,Cartesian::z)-PmC[2]*nucm1(Cartesian::s,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));

            
        }
        
        
         // d-s
        if ( _lmax_row > 1){
           //cout << "\t setting d-s" << endl;
                nuc(Cartesian::yy,Cartesian::s) =PmA[1]*nuc(Cartesian::y,Cartesian::s)-PmC[1]*nucm1(Cartesian::y,Cartesian::s)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));
                nuc(Cartesian::xy,Cartesian::s) =PmA[0]*nuc(Cartesian::y,Cartesian::s)-PmC[0]*nucm1(Cartesian::y,Cartesian::s);
                nuc(Cartesian::yz,Cartesian::s) =PmA[1]*nuc(Cartesian::z,Cartesian::s)-PmC[1]*nucm1(Cartesian::z,Cartesian::s);
                nuc(Cartesian::xx,Cartesian::s) =PmA[0]*nuc(Cartesian::x,Cartesian::s)-PmC[0]*nucm1(Cartesian::x,Cartesian::s)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));
                nuc(Cartesian::xz,Cartesian::s) =PmA[0]*nuc(Cartesian::z,Cartesian::s)-PmC[0]*nucm1(Cartesian::z,Cartesian::s);
                nuc(Cartesian::zz,Cartesian::s) =PmA[2]*nuc(Cartesian::z,Cartesian::s)-PmC[2]*nucm1(Cartesian::z,Cartesian::s)+_fak*(nuc(Cartesian::s,Cartesian::s)-nucm1(Cartesian::s,Cartesian::s));

         
        }
        
        
        if ( _lsum > 2 ){
            //cout << "\t setting p-d" << endl;
            if ( _lmax_col > 0){
            //s-p-2
                nucm2(Cartesian::s,Cartesian::y) =PmB[1]*nucm2(Cartesian::s,Cartesian::s)-PmC[1]*nucm3(Cartesian::s,Cartesian::s);
                nucm2(Cartesian::s,Cartesian::x) =PmB[0]*nucm2(Cartesian::s,Cartesian::s)-PmC[0]*nucm3(Cartesian::s,Cartesian::s);
                nucm2(Cartesian::s,Cartesian::z) =PmB[2]*nucm2(Cartesian::s,Cartesian::s)-PmC[2]*nucm3(Cartesian::s,Cartesian::s);
            }
            if ( _lmax_row > 0){    
            //p-s-2
                nucm2(Cartesian::y,Cartesian::s) =PmA[1]*nucm2(Cartesian::s,Cartesian::s)-PmC[1]*nucm3(Cartesian::s,Cartesian::s);
                nucm2(Cartesian::x,Cartesian::s) =PmA[0]*nucm2(Cartesian::s,Cartesian::s)-PmC[0]*nucm3(Cartesian::s,Cartesian::s);
                nucm2(Cartesian::z,Cartesian::s) =PmA[2]*nucm2(Cartesian::s,Cartesian::s)-PmC[2]*nucm3(Cartesian::s,Cartesian::s);
            }

            if ( _lmax_row > 0 && _lmax_col > 0 ) {
            //p-p-1
                nucm1(Cartesian::y,Cartesian::y) =PmA[1]*nucm1(Cartesian::s,Cartesian::y)-PmC[1]*nucm2(Cartesian::s,Cartesian::y)+_fak*(nucm1(Cartesian::s,Cartesian::s)-nucm2(Cartesian::s,Cartesian::s));
                nucm1(Cartesian::y,Cartesian::x) =PmA[1]*nucm1(Cartesian::s,Cartesian::x)-PmC[1]*nucm2(Cartesian::s,Cartesian::x);
                nucm1(Cartesian::y,Cartesian::z) =PmA[1]*nucm1(Cartesian::s,Cartesian::z)-PmC[1]*nucm2(Cartesian::s,Cartesian::z);
                nucm1(Cartesian::x,Cartesian::y) =PmA[0]*nucm1(Cartesian::s,Cartesian::y)-PmC[0]*nucm2(Cartesian::s,Cartesian::y);
                nucm1(Cartesian::x,Cartesian::x) =PmA[0]*nucm1(Cartesian::s,Cartesian::x)-PmC[0]*nucm2(Cartesian::s,Cartesian::x)+_fak*(nucm1(Cartesian::s,Cartesian::s)-nucm2(Cartesian::s,Cartesian::s));
                nucm1(Cartesian::x,Cartesian::z) =PmA[0]*nucm1(Cartesian::s,Cartesian::z)-PmC[0]*nucm2(Cartesian::s,Cartesian::z);
                nucm1(Cartesian::z,Cartesian::y) =PmA[2]*nucm1(Cartesian::s,Cartesian::y)-PmC[2]*nucm2(Cartesian::s,Cartesian::y);
                nucm1(Cartesian::z,Cartesian::x) =PmA[2]*nucm1(Cartesian::s,Cartesian::x)-PmC[2]*nucm2(Cartesian::s,Cartesian::x);
                nucm1(Cartesian::z,Cartesian::z) =PmA[2]*nucm1(Cartesian::s,Cartesian::z)-PmC[2]*nucm2(Cartesian::s,Cartesian::z)+_fak*(nucm1(Cartesian::s,Cartesian::s)-nucm2(Cartesian::s,Cartesian::s));
        }
        }
        
        
        
        
         //p-d
        if ( _lmax_row > 0 && _lmax_col > 1){
            //cout << "\t setting p-d" << endl;
            
        
            // p-d-0
                nuc(Cartesian::y,Cartesian::yy) =PmB[1]*nuc(Cartesian::y,Cartesian::y)-PmC[1]*nucm1(Cartesian::y,Cartesian::y)+_fak*(nuc(Cartesian::y,Cartesian::s)-nucm1(Cartesian::y,Cartesian::s))+_fak*(nuc(Cartesian::s,Cartesian::y)-nucm1(Cartesian::s,Cartesian::y));
                nuc(Cartesian::y,Cartesian::xy) =PmB[0]*nuc(Cartesian::y,Cartesian::y)-PmC[0]*nucm1(Cartesian::y,Cartesian::y);
                nuc(Cartesian::y,Cartesian::yz) =PmB[1]*nuc(Cartesian::y,Cartesian::z)-PmC[1]*nucm1(Cartesian::y,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::z)-nucm1(Cartesian::s,Cartesian::z));
                nuc(Cartesian::y,Cartesian::xx) =PmB[0]*nuc(Cartesian::y,Cartesian::x)-PmC[0]*nucm1(Cartesian::y,Cartesian::x)+_fak*(nuc(Cartesian::y,Cartesian::s)-nucm1(Cartesian::y,Cartesian::s));
                nuc(Cartesian::y,Cartesian::xz) =PmB[0]*nuc(Cartesian::y,Cartesian::z)-PmC[0]*nucm1(Cartesian::y,Cartesian::z);
                nuc(Cartesian::y,Cartesian::zz) =PmB[2]*nuc(Cartesian::y,Cartesian::z)-PmC[2]*nucm1(Cartesian::y,Cartesian::z)+_fak*(nuc(Cartesian::y,Cartesian::s)-nucm1(Cartesian::y,Cartesian::s));
                nuc(Cartesian::x,Cartesian::yy) =PmB[1]*nuc(Cartesian::x,Cartesian::y)-PmC[1]*nucm1(Cartesian::x,Cartesian::y)+_fak*(nuc(Cartesian::x,Cartesian::s)-nucm1(Cartesian::x,Cartesian::s));
                nuc(Cartesian::x,Cartesian::xy) =PmB[0]*nuc(Cartesian::x,Cartesian::y)-PmC[0]*nucm1(Cartesian::x,Cartesian::y)+_fak*(nuc(Cartesian::s,Cartesian::y)-nucm1(Cartesian::s,Cartesian::y));
                nuc(Cartesian::x,Cartesian::yz) =PmB[1]*nuc(Cartesian::x,Cartesian::z)-PmC[1]*nucm1(Cartesian::x,Cartesian::z);
                nuc(Cartesian::x,Cartesian::xx) =PmB[0]*nuc(Cartesian::x,Cartesian::x)-PmC[0]*nucm1(Cartesian::x,Cartesian::x)+_fak*(nuc(Cartesian::x,Cartesian::s)-nucm1(Cartesian::x,Cartesian::s))+_fak*(nuc(Cartesian::s,Cartesian::x)-nucm1(Cartesian::s,Cartesian::x));
                nuc(Cartesian::x,Cartesian::xz) =PmB[0]*nuc(Cartesian::x,Cartesian::z)-PmC[0]*nucm1(Cartesian::x,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::z)-nucm1(Cartesian::s,Cartesian::z));
                nuc(Cartesian::x,Cartesian::zz) =PmB[2]*nuc(Cartesian::x,Cartesian::z)-PmC[2]*nucm1(Cartesian::x,Cartesian::z)+_fak*(nuc(Cartesian::x,Cartesian::s)-nucm1(Cartesian::x,Cartesian::s));
                nuc(Cartesian::z,Cartesian::yy) =PmB[1]*nuc(Cartesian::z,Cartesian::y)-PmC[1]*nucm1(Cartesian::z,Cartesian::y)+_fak*(nuc(Cartesian::z,Cartesian::s)-nucm1(Cartesian::z,Cartesian::s));
                nuc(Cartesian::z,Cartesian::xy) =PmB[0]*nuc(Cartesian::z,Cartesian::y)-PmC[0]*nucm1(Cartesian::z,Cartesian::y);
                nuc(Cartesian::z,Cartesian::yz) =PmB[1]*nuc(Cartesian::z,Cartesian::z)-PmC[1]*nucm1(Cartesian::z,Cartesian::z);
                nuc(Cartesian::z,Cartesian::xx) =PmB[0]*nuc(Cartesian::z,Cartesian::x)-PmC[0]*nucm1(Cartesian::z,Cartesian::x)+_fak*(nuc(Cartesian::z,Cartesian::s)-nucm1(Cartesian::z,Cartesian::s));
                nuc(Cartesian::z,Cartesian::xz) =PmB[0]*nuc(Cartesian::z,Cartesian::z)-PmC[0]*nucm1(Cartesian::z,Cartesian::z);
                nuc(Cartesian::z,Cartesian::zz) =PmB[2]*nuc(Cartesian::z,Cartesian::z)-PmC[2]*nucm1(Cartesian::z,Cartesian::z)+_fak*(nuc(Cartesian::z,Cartesian::s)-nucm1(Cartesian::z,Cartesian::s))+_fak*(nuc(Cartesian::s,Cartesian::z)-nucm1(Cartesian::s,Cartesian::z));

         
        }

       
        
        
        // d-p
        if ( _lmax_row >1 && _lmax_col > 0){
           //cout << "\t setting d-p" << endl;
            
                nuc(Cartesian::yy,Cartesian::y) =PmA[1]*nuc(Cartesian::y,Cartesian::y)-PmC[1]*nucm1(Cartesian::y,Cartesian::y)+_fak*(nuc(Cartesian::s,Cartesian::y)-nucm1(Cartesian::s,Cartesian::y))+_fak*(nuc(Cartesian::y,Cartesian::s)-nucm1(Cartesian::y,Cartesian::s));
                nuc(Cartesian::yy,Cartesian::x) =PmA[1]*nuc(Cartesian::y,Cartesian::x)-PmC[1]*nucm1(Cartesian::y,Cartesian::x)+_fak*(nuc(Cartesian::s,Cartesian::x)-nucm1(Cartesian::s,Cartesian::x));
                nuc(Cartesian::yy,Cartesian::z) =PmA[1]*nuc(Cartesian::y,Cartesian::z)-PmC[1]*nucm1(Cartesian::y,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::z)-nucm1(Cartesian::s,Cartesian::z));
                nuc(Cartesian::xy,Cartesian::y) =PmA[0]*nuc(Cartesian::y,Cartesian::y)-PmC[0]*nucm1(Cartesian::y,Cartesian::y);
                nuc(Cartesian::xy,Cartesian::x) =PmA[0]*nuc(Cartesian::y,Cartesian::x)-PmC[0]*nucm1(Cartesian::y,Cartesian::x)+_fak*(nuc(Cartesian::y,Cartesian::s)-nucm1(Cartesian::y,Cartesian::s));
                nuc(Cartesian::xy,Cartesian::z) =PmA[0]*nuc(Cartesian::y,Cartesian::z)-PmC[0]*nucm1(Cartesian::y,Cartesian::z);
                nuc(Cartesian::yz,Cartesian::y) =PmA[1]*nuc(Cartesian::z,Cartesian::y)-PmC[1]*nucm1(Cartesian::z,Cartesian::y)+_fak*(nuc(Cartesian::z,Cartesian::s)-nucm1(Cartesian::z,Cartesian::s));
                nuc(Cartesian::yz,Cartesian::x) =PmA[1]*nuc(Cartesian::z,Cartesian::x)-PmC[1]*nucm1(Cartesian::z,Cartesian::x);
                nuc(Cartesian::yz,Cartesian::z) =PmA[1]*nuc(Cartesian::z,Cartesian::z)-PmC[1]*nucm1(Cartesian::z,Cartesian::z);
                nuc(Cartesian::xx,Cartesian::y) =PmA[0]*nuc(Cartesian::x,Cartesian::y)-PmC[0]*nucm1(Cartesian::x,Cartesian::y)+_fak*(nuc(Cartesian::s,Cartesian::y)-nucm1(Cartesian::s,Cartesian::y));
                nuc(Cartesian::xx,Cartesian::x) =PmA[0]*nuc(Cartesian::x,Cartesian::x)-PmC[0]*nucm1(Cartesian::x,Cartesian::x)+_fak*(nuc(Cartesian::s,Cartesian::x)-nucm1(Cartesian::s,Cartesian::x))+_fak*(nuc(Cartesian::x,Cartesian::s)-nucm1(Cartesian::x,Cartesian::s));
                nuc(Cartesian::xx,Cartesian::z) =PmA[0]*nuc(Cartesian::x,Cartesian::z)-PmC[0]*nucm1(Cartesian::x,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::z)-nucm1(Cartesian::s,Cartesian::z));
                nuc(Cartesian::xz,Cartesian::y) =PmA[0]*nuc(Cartesian::z,Cartesian::y)-PmC[0]*nucm1(Cartesian::z,Cartesian::y);
                nuc(Cartesian::xz,Cartesian::x) =PmA[0]*nuc(Cartesian::z,Cartesian::x)-PmC[0]*nucm1(Cartesian::z,Cartesian::x)+_fak*(nuc(Cartesian::z,Cartesian::s)-nucm1(Cartesian::z,Cartesian::s));
                nuc(Cartesian::xz,Cartesian::z) =PmA[0]*nuc(Cartesian::z,Cartesian::z)-PmC[0]*nucm1(Cartesian::z,Cartesian::z);
                nuc(Cartesian::zz,Cartesian::y) =PmA[2]*nuc(Cartesian::z,Cartesian::y)-PmC[2]*nucm1(Cartesian::z,Cartesian::y)+_fak*(nuc(Cartesian::s,Cartesian::y)-nucm1(Cartesian::s,Cartesian::y));
                nuc(Cartesian::zz,Cartesian::x) =PmA[2]*nuc(Cartesian::z,Cartesian::x)-PmC[2]*nucm1(Cartesian::z,Cartesian::x)+_fak*(nuc(Cartesian::s,Cartesian::x)-nucm1(Cartesian::s,Cartesian::x));
                nuc(Cartesian::zz,Cartesian::z) =PmA[2]*nuc(Cartesian::z,Cartesian::z)-PmC[2]*nucm1(Cartesian::z,Cartesian::z)+_fak*(nuc(Cartesian::s,Cartesian::z)-nucm1(Cartesian::s,Cartesian::z))+_fak*(nuc(Cartesian::z,Cartesian::s)-nucm1(Cartesian::z,Cartesian::s));


        }
        
         
        if ( _lsum > 2 && _lmax_col>1 && _lmax_row>1){
            
            
           if (  _lmax_col > 0){
            //s-p-3
            
            nucm3(Cartesian::s, Cartesian::y) = PmB[1] * nucm3(Cartesian::s, Cartesian::s) - PmC[1] * nucm4(Cartesian::s, Cartesian::s);
            nucm3(Cartesian::s, Cartesian::x) = PmB[0] * nucm3(Cartesian::s, Cartesian::s) - PmC[0] * nucm4(Cartesian::s, Cartesian::s);
            nucm3(Cartesian::s, Cartesian::z) = PmB[2] * nucm3(Cartesian::s, Cartesian::s) - PmC[2] * nucm4(Cartesian::s, Cartesian::s);
            }
            
             if ( _lmax_row > 0 && _lmax_col > 0 ) {
            //p-p-2
            
            nucm2(Cartesian::y,Cartesian::y) =PmA[1]*nucm2(Cartesian::s,Cartesian::y)-PmC[1]*nucm3(Cartesian::s,Cartesian::y)+_fak*(nucm2(Cartesian::s,Cartesian::s)-nucm3(Cartesian::s,Cartesian::s));
            nucm2(Cartesian::y,Cartesian::x) =PmA[1]*nucm2(Cartesian::s,Cartesian::x)-PmC[1]*nucm3(Cartesian::s,Cartesian::x);
            nucm2(Cartesian::y,Cartesian::z) =PmA[1]*nucm2(Cartesian::s,Cartesian::z)-PmC[1]*nucm3(Cartesian::s,Cartesian::z);
            nucm2(Cartesian::x,Cartesian::y) =PmA[0]*nucm2(Cartesian::s,Cartesian::y)-PmC[0]*nucm3(Cartesian::s,Cartesian::y);
            nucm2(Cartesian::x,Cartesian::x) =PmA[0]*nucm2(Cartesian::s,Cartesian::x)-PmC[0]*nucm3(Cartesian::s,Cartesian::x)+_fak*(nucm2(Cartesian::s,Cartesian::s)-nucm3(Cartesian::s,Cartesian::s));
            nucm2(Cartesian::x,Cartesian::z) =PmA[0]*nucm2(Cartesian::s,Cartesian::z)-PmC[0]*nucm3(Cartesian::s,Cartesian::z);
            nucm2(Cartesian::z,Cartesian::y) =PmA[2]*nucm2(Cartesian::s,Cartesian::y)-PmC[2]*nucm3(Cartesian::s,Cartesian::y);
            nucm2(Cartesian::z,Cartesian::x) =PmA[2]*nucm2(Cartesian::s,Cartesian::x)-PmC[2]*nucm3(Cartesian::s,Cartesian::x);
            nucm2(Cartesian::z,Cartesian::z) =PmA[2]*nucm2(Cartesian::s,Cartesian::z)-PmC[2]*nucm3(Cartesian::s,Cartesian::z)+_fak*(nucm2(Cartesian::s,Cartesian::s)-nucm3(Cartesian::s,Cartesian::s));
             }
            
            
              if ( _lmax_row > 0 && _lmax_col > 1){
                  //s-d-1
             
            nucm1(Cartesian::s,Cartesian::yy) =PmB[1]*nucm1(Cartesian::s,Cartesian::y)-PmC[1]*nucm2(Cartesian::s,Cartesian::y)+_fak*(nucm1(Cartesian::s,Cartesian::s)-nucm2(Cartesian::s,Cartesian::s));
            nucm1(Cartesian::s,Cartesian::xy) =PmB[0]*nucm1(Cartesian::s,Cartesian::y)-PmC[0]*nucm2(Cartesian::s,Cartesian::y);
            nucm1(Cartesian::s,Cartesian::yz) =PmB[1]*nucm1(Cartesian::s,Cartesian::z)-PmC[1]*nucm2(Cartesian::s,Cartesian::z);
            nucm1(Cartesian::s,Cartesian::xx) =PmB[0]*nucm1(Cartesian::s,Cartesian::x)-PmC[0]*nucm2(Cartesian::s,Cartesian::x)+_fak*(nucm1(Cartesian::s,Cartesian::s)-nucm2(Cartesian::s,Cartesian::s));
            nucm1(Cartesian::s,Cartesian::xz) =PmB[0]*nucm1(Cartesian::s,Cartesian::z)-PmC[0]*nucm2(Cartesian::s,Cartesian::z);
            nucm1(Cartesian::s,Cartesian::zz) =PmB[2]*nucm1(Cartesian::s,Cartesian::z)-PmC[2]*nucm2(Cartesian::s,Cartesian::z)+_fak*(nucm1(Cartesian::s,Cartesian::s)-nucm2(Cartesian::s,Cartesian::s));
        }
            
             if ( _lmax_row > 0 && _lmax_col > 1){
            //p-d-1
            nucm1(Cartesian::y,Cartesian::yy) =PmB[1]*nucm1(Cartesian::y,Cartesian::y)-PmC[1]*nucm2(Cartesian::y,Cartesian::y)+_fak*(nucm1(Cartesian::y,Cartesian::s)-nucm2(Cartesian::y,Cartesian::s))+_fak*(nucm1(Cartesian::s,Cartesian::y)-nucm2(Cartesian::s,Cartesian::y));
            nucm1(Cartesian::y,Cartesian::xy) =PmB[0]*nucm1(Cartesian::y,Cartesian::y)-PmC[0]*nucm2(Cartesian::y,Cartesian::y);
            nucm1(Cartesian::y,Cartesian::yz) =PmB[1]*nucm1(Cartesian::y,Cartesian::z)-PmC[1]*nucm2(Cartesian::y,Cartesian::z)+_fak*(nucm1(Cartesian::s,Cartesian::z)-nucm2(Cartesian::s,Cartesian::z));
            nucm1(Cartesian::y,Cartesian::xx) =PmB[0]*nucm1(Cartesian::y,Cartesian::x)-PmC[0]*nucm2(Cartesian::y,Cartesian::x)+_fak*(nucm1(Cartesian::y,Cartesian::s)-nucm2(Cartesian::y,Cartesian::s));
            nucm1(Cartesian::y,Cartesian::xz) =PmB[0]*nucm1(Cartesian::y,Cartesian::z)-PmC[0]*nucm2(Cartesian::y,Cartesian::z);
            nucm1(Cartesian::y,Cartesian::zz) =PmB[2]*nucm1(Cartesian::y,Cartesian::z)-PmC[2]*nucm2(Cartesian::y,Cartesian::z)+_fak*(nucm1(Cartesian::y,Cartesian::s)-nucm2(Cartesian::y,Cartesian::s));
            nucm1(Cartesian::x,Cartesian::yy) =PmB[1]*nucm1(Cartesian::x,Cartesian::y)-PmC[1]*nucm2(Cartesian::x,Cartesian::y)+_fak*(nucm1(Cartesian::x,Cartesian::s)-nucm2(Cartesian::x,Cartesian::s));
            nucm1(Cartesian::x,Cartesian::xy) =PmB[0]*nucm1(Cartesian::x,Cartesian::y)-PmC[0]*nucm2(Cartesian::x,Cartesian::y)+_fak*(nucm1(Cartesian::s,Cartesian::y)-nucm2(Cartesian::s,Cartesian::y));
            nucm1(Cartesian::x,Cartesian::yz) =PmB[1]*nucm1(Cartesian::x,Cartesian::z)-PmC[1]*nucm2(Cartesian::x,Cartesian::z);
            nucm1(Cartesian::x,Cartesian::xx) =PmB[0]*nucm1(Cartesian::x,Cartesian::x)-PmC[0]*nucm2(Cartesian::x,Cartesian::x)+_fak*(nucm1(Cartesian::x,Cartesian::s)-nucm2(Cartesian::x,Cartesian::s))+_fak*(nucm1(Cartesian::s,Cartesian::x)-nucm2(Cartesian::s,Cartesian::x));
            nucm1(Cartesian::x,Cartesian::xz) =PmB[0]*nucm1(Cartesian::x,Cartesian::z)-PmC[0]*nucm2(Cartesian::x,Cartesian::z)+_fak*(nucm1(Cartesian::s,Cartesian::z)-nucm2(Cartesian::s,Cartesian::z));
            nucm1(Cartesian::x,Cartesian::zz) =PmB[2]*nucm1(Cartesian::x,Cartesian::z)-PmC[2]*nucm2(Cartesian::x,Cartesian::z)+_fak*(nucm1(Cartesian::x,Cartesian::s)-nucm2(Cartesian::x,Cartesian::s));
            nucm1(Cartesian::z,Cartesian::yy) =PmB[1]*nucm1(Cartesian::z,Cartesian::y)-PmC[1]*nucm2(Cartesian::z,Cartesian::y)+_fak*(nucm1(Cartesian::z,Cartesian::s)-nucm2(Cartesian::z,Cartesian::s));
            nucm1(Cartesian::z,Cartesian::xy) =PmB[0]*nucm1(Cartesian::z,Cartesian::y)-PmC[0]*nucm2(Cartesian::z,Cartesian::y);
            nucm1(Cartesian::z,Cartesian::yz) =PmB[1]*nucm1(Cartesian::z,Cartesian::z)-PmC[1]*nucm2(Cartesian::z,Cartesian::z);
            nucm1(Cartesian::z,Cartesian::xx) =PmB[0]*nucm1(Cartesian::z,Cartesian::x)-PmC[0]*nucm2(Cartesian::z,Cartesian::x)+_fak*(nucm1(Cartesian::z,Cartesian::s)-nucm2(Cartesian::z,Cartesian::s));
            nucm1(Cartesian::z,Cartesian::xz) =PmB[0]*nucm1(Cartesian::z,Cartesian::z)-PmC[0]*nucm2(Cartesian::z,Cartesian::z);
            nucm1(Cartesian::z,Cartesian::zz) =PmB[2]*nucm1(Cartesian::z,Cartesian::z)-PmC[2]*nucm2(Cartesian::z,Cartesian::z)+_fak*(nucm1(Cartesian::z,Cartesian::s)-nucm2(Cartesian::z,Cartesian::s))+_fak*(nucm1(Cartesian::s,Cartesian::z)-nucm2(Cartesian::s,Cartesian::z));
             }
            
        }
     
            
            
            
            
            
            
        
        
        // d-d
        if ( _lmax_row > 1 && _lmax_col > 1 ){
             // cout << "\t setting d-d" << endl;
 
            
            //d-d-0
            nuc(Cartesian::yy,Cartesian::yy) =PmA[1]*nuc(Cartesian::y,Cartesian::yy)-PmC[1]*nucm1(Cartesian::y,Cartesian::yy)+_fak*(nuc(Cartesian::s,Cartesian::yy)-nucm1(Cartesian::s,Cartesian::yy))+_fak2*(nuc(Cartesian::y,Cartesian::y)-nucm1(Cartesian::y,Cartesian::y));
            nuc(Cartesian::yy,Cartesian::xy) =PmA[1]*nuc(Cartesian::y,Cartesian::xy)-PmC[1]*nucm1(Cartesian::y,Cartesian::xy)+_fak*(nuc(Cartesian::s,Cartesian::xy)-nucm1(Cartesian::s,Cartesian::xy))+_fak*(nuc(Cartesian::y,Cartesian::x)-nucm1(Cartesian::y,Cartesian::x));
            nuc(Cartesian::yy,Cartesian::yz) =PmA[1]*nuc(Cartesian::y,Cartesian::yz)-PmC[1]*nucm1(Cartesian::y,Cartesian::yz)+_fak*(nuc(Cartesian::s,Cartesian::yz)-nucm1(Cartesian::s,Cartesian::yz))+_fak*(nuc(Cartesian::y,Cartesian::z)-nucm1(Cartesian::y,Cartesian::z));
            nuc(Cartesian::yy,Cartesian::xx) =PmA[1]*nuc(Cartesian::y,Cartesian::xx)-PmC[1]*nucm1(Cartesian::y,Cartesian::xx)+_fak*(nuc(Cartesian::s,Cartesian::xx)-nucm1(Cartesian::s,Cartesian::xx));
            nuc(Cartesian::yy,Cartesian::xz) =PmA[1]*nuc(Cartesian::y,Cartesian::xz)-PmC[1]*nucm1(Cartesian::y,Cartesian::xz)+_fak*(nuc(Cartesian::s,Cartesian::xz)-nucm1(Cartesian::s,Cartesian::xz));
            nuc(Cartesian::yy,Cartesian::zz) =PmA[1]*nuc(Cartesian::y,Cartesian::zz)-PmC[1]*nucm1(Cartesian::y,Cartesian::zz)+_fak*(nuc(Cartesian::s,Cartesian::zz)-nucm1(Cartesian::s,Cartesian::zz));
            nuc(Cartesian::xy,Cartesian::yy) =PmA[0]*nuc(Cartesian::y,Cartesian::yy)-PmC[0]*nucm1(Cartesian::y,Cartesian::yy);
            nuc(Cartesian::xy,Cartesian::xy) =PmA[0]*nuc(Cartesian::y,Cartesian::xy)-PmC[0]*nucm1(Cartesian::y,Cartesian::xy)+_fak*(nuc(Cartesian::y,Cartesian::y)-nucm1(Cartesian::y,Cartesian::y));
            nuc(Cartesian::xy,Cartesian::yz) =PmA[0]*nuc(Cartesian::y,Cartesian::yz)-PmC[0]*nucm1(Cartesian::y,Cartesian::yz);
            nuc(Cartesian::xy,Cartesian::xx) =PmA[0]*nuc(Cartesian::y,Cartesian::xx)-PmC[0]*nucm1(Cartesian::y,Cartesian::xx)+_fak2*(nuc(Cartesian::y,Cartesian::x)-nucm1(Cartesian::y,Cartesian::x));
            nuc(Cartesian::xy,Cartesian::xz) =PmA[0]*nuc(Cartesian::y,Cartesian::xz)-PmC[0]*nucm1(Cartesian::y,Cartesian::xz)+_fak*(nuc(Cartesian::y,Cartesian::z)-nucm1(Cartesian::y,Cartesian::z));
            nuc(Cartesian::xy,Cartesian::zz) =PmA[0]*nuc(Cartesian::y,Cartesian::zz)-PmC[0]*nucm1(Cartesian::y,Cartesian::zz);
            nuc(Cartesian::yz,Cartesian::yy) =PmA[1]*nuc(Cartesian::z,Cartesian::yy)-PmC[1]*nucm1(Cartesian::z,Cartesian::yy)+_fak2*(nuc(Cartesian::z,Cartesian::y)-nucm1(Cartesian::z,Cartesian::y));
            nuc(Cartesian::yz,Cartesian::xy) =PmA[1]*nuc(Cartesian::z,Cartesian::xy)-PmC[1]*nucm1(Cartesian::z,Cartesian::xy)+_fak*(nuc(Cartesian::z,Cartesian::x)-nucm1(Cartesian::z,Cartesian::x));
            nuc(Cartesian::yz,Cartesian::yz) =PmA[1]*nuc(Cartesian::z,Cartesian::yz)-PmC[1]*nucm1(Cartesian::z,Cartesian::yz)+_fak*(nuc(Cartesian::z,Cartesian::z)-nucm1(Cartesian::z,Cartesian::z));
            nuc(Cartesian::yz,Cartesian::xx) =PmA[1]*nuc(Cartesian::z,Cartesian::xx)-PmC[1]*nucm1(Cartesian::z,Cartesian::xx);
            nuc(Cartesian::yz,Cartesian::xz) =PmA[1]*nuc(Cartesian::z,Cartesian::xz)-PmC[1]*nucm1(Cartesian::z,Cartesian::xz);
            nuc(Cartesian::yz,Cartesian::zz) =PmA[1]*nuc(Cartesian::z,Cartesian::zz)-PmC[1]*nucm1(Cartesian::z,Cartesian::zz);
            nuc(Cartesian::xx,Cartesian::yy) =PmA[0]*nuc(Cartesian::x,Cartesian::yy)-PmC[0]*nucm1(Cartesian::x,Cartesian::yy)+_fak*(nuc(Cartesian::s,Cartesian::yy)-nucm1(Cartesian::s,Cartesian::yy));
            nuc(Cartesian::xx,Cartesian::xy) =PmA[0]*nuc(Cartesian::x,Cartesian::xy)-PmC[0]*nucm1(Cartesian::x,Cartesian::xy)+_fak*(nuc(Cartesian::s,Cartesian::xy)-nucm1(Cartesian::s,Cartesian::xy))+_fak*(nuc(Cartesian::x,Cartesian::y)-nucm1(Cartesian::x,Cartesian::y));
            nuc(Cartesian::xx,Cartesian::yz) =PmA[0]*nuc(Cartesian::x,Cartesian::yz)-PmC[0]*nucm1(Cartesian::x,Cartesian::yz)+_fak*(nuc(Cartesian::s,Cartesian::yz)-nucm1(Cartesian::s,Cartesian::yz));
            nuc(Cartesian::xx,Cartesian::xx) =PmA[0]*nuc(Cartesian::x,Cartesian::xx)-PmC[0]*nucm1(Cartesian::x,Cartesian::xx)+_fak*(nuc(Cartesian::s,Cartesian::xx)-nucm1(Cartesian::s,Cartesian::xx))+_fak2*(nuc(Cartesian::x,Cartesian::x)-nucm1(Cartesian::x,Cartesian::x));
            nuc(Cartesian::xx,Cartesian::xz) =PmA[0]*nuc(Cartesian::x,Cartesian::xz)-PmC[0]*nucm1(Cartesian::x,Cartesian::xz)+_fak*(nuc(Cartesian::s,Cartesian::xz)-nucm1(Cartesian::s,Cartesian::xz))+_fak*(nuc(Cartesian::x,Cartesian::z)-nucm1(Cartesian::x,Cartesian::z));
            nuc(Cartesian::xx,Cartesian::zz) =PmA[0]*nuc(Cartesian::x,Cartesian::zz)-PmC[0]*nucm1(Cartesian::x,Cartesian::zz)+_fak*(nuc(Cartesian::s,Cartesian::zz)-nucm1(Cartesian::s,Cartesian::zz));
            nuc(Cartesian::xz,Cartesian::yy) =PmA[0]*nuc(Cartesian::z,Cartesian::yy)-PmC[0]*nucm1(Cartesian::z,Cartesian::yy);
            nuc(Cartesian::xz,Cartesian::xy) =PmA[0]*nuc(Cartesian::z,Cartesian::xy)-PmC[0]*nucm1(Cartesian::z,Cartesian::xy)+_fak*(nuc(Cartesian::z,Cartesian::y)-nucm1(Cartesian::z,Cartesian::y));
            nuc(Cartesian::xz,Cartesian::yz) =PmA[0]*nuc(Cartesian::z,Cartesian::yz)-PmC[0]*nucm1(Cartesian::z,Cartesian::yz);
            nuc(Cartesian::xz,Cartesian::xx) =PmA[0]*nuc(Cartesian::z,Cartesian::xx)-PmC[0]*nucm1(Cartesian::z,Cartesian::xx)+_fak2*(nuc(Cartesian::z,Cartesian::x)-nucm1(Cartesian::z,Cartesian::x));
            nuc(Cartesian::xz,Cartesian::xz) =PmA[0]*nuc(Cartesian::z,Cartesian::xz)-PmC[0]*nucm1(Cartesian::z,Cartesian::xz)+_fak*(nuc(Cartesian::z,Cartesian::z)-nucm1(Cartesian::z,Cartesian::z));
            nuc(Cartesian::xz,Cartesian::zz) =PmA[0]*nuc(Cartesian::z,Cartesian::zz)-PmC[0]*nucm1(Cartesian::z,Cartesian::zz);
            nuc(Cartesian::zz,Cartesian::yy) =PmA[2]*nuc(Cartesian::z,Cartesian::yy)-PmC[2]*nucm1(Cartesian::z,Cartesian::yy)+_fak*(nuc(Cartesian::s,Cartesian::yy)-nucm1(Cartesian::s,Cartesian::yy));
            nuc(Cartesian::zz,Cartesian::xy) =PmA[2]*nuc(Cartesian::z,Cartesian::xy)-PmC[2]*nucm1(Cartesian::z,Cartesian::xy)+_fak*(nuc(Cartesian::s,Cartesian::xy)-nucm1(Cartesian::s,Cartesian::xy));
            nuc(Cartesian::zz,Cartesian::yz) =PmA[2]*nuc(Cartesian::z,Cartesian::yz)-PmC[2]*nucm1(Cartesian::z,Cartesian::yz)+_fak*(nuc(Cartesian::s,Cartesian::yz)-nucm1(Cartesian::s,Cartesian::yz))+_fak*(nuc(Cartesian::z,Cartesian::y)-nucm1(Cartesian::z,Cartesian::y));
            nuc(Cartesian::zz,Cartesian::xx) =PmA[2]*nuc(Cartesian::z,Cartesian::xx)-PmC[2]*nucm1(Cartesian::z,Cartesian::xx)+_fak*(nuc(Cartesian::s,Cartesian::xx)-nucm1(Cartesian::s,Cartesian::xx));
            nuc(Cartesian::zz,Cartesian::xz) =PmA[2]*nuc(Cartesian::z,Cartesian::xz)-PmC[2]*nucm1(Cartesian::z,Cartesian::xz)+_fak*(nuc(Cartesian::s,Cartesian::xz)-nucm1(Cartesian::s,Cartesian::xz))+_fak*(nuc(Cartesian::z,Cartesian::x)-nucm1(Cartesian::z,Cartesian::x));
            nuc(Cartesian::zz,Cartesian::zz) =PmA[2]*nuc(Cartesian::z,Cartesian::zz)-PmC[2]*nucm1(Cartesian::z,Cartesian::zz)+_fak*(nuc(Cartesian::s,Cartesian::zz)-nucm1(Cartesian::s,Cartesian::zz))+_fak2*(nuc(Cartesian::z,Cartesian::z)-nucm1(Cartesian::z,Cartesian::z));


            
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
        for ( int i = 0; i< _matrix.size1(); i++ ) {
            for (int j = 0; j < _matrix.size2(); j++){
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
   for ( int j = 0; j < _atoms.size(); j++){

            
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

