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
#include <votca/tools/linalg.h>

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
#include <votca/ctp/logger.h>
#include <votca/tools/linalg.h>


using namespace votca::tools;

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;

    
    void AOOverlap::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, AOBasis* ecp ) {
        /*cout << "\nAO block: "<< endl;
        cout << "\t row: " << _shell_row->getType() << " at " << _shell_row->getPos() << endl;
        cout << "\t col: " << _shell_col->getType() << " at " << _shell_col->getPos() << endl;*/
       
        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
    
        
        /* FOR CONTRACTED FUNCTIONS, ADD LOOP OVER ALL DECAYS IN CONTRACTION
         * MULTIPLY THE TRANSFORMATION MATRICES BY APPROPRIATE CONTRACTION 
         * COEFFICIENTS, AND ADD TO matrix(i,j)
         */

        // get shell positions
        const vec& _pos_row = _shell_row->getPos();
        const vec& _pos_col = _shell_col->getPos();
        const vec  _diff    = _pos_row - _pos_col;
        std::vector<double> _pma (3,0.0);
        std::vector<double> _pmb (3,0.0);
          
        double _distsq = (_diff.getX()*_diff.getX()) + (_diff.getY()*_diff.getY()) + (_diff.getZ()*_diff.getZ());   


       // cout << "row shell is " << _shell_row->getSize() << " -fold contracted!" << endl;
        //cout << "col shell is " << _shell_col->getSize() << " -fold contracted!" << endl;
        
        typedef std::vector< AOGaussianPrimitive* >::iterator GaussianIterator;
        // iterate over Gaussians in this _shell_row
        for ( GaussianIterator itr = _shell_row->firstGaussian(); itr != _shell_row->lastGaussian(); ++itr){
            // iterate over Gaussians in this _shell_col
            const double& _decay_row = (*itr)->decay;
            
            for ( GaussianIterator itc = _shell_col->firstGaussian(); itc != _shell_col->lastGaussian(); ++itc){
           
            
           

            // get decay constants 
            const double& _decay_col = (*itc)->decay;
            
            // some helpers
            const double _fak  = 0.5/(_decay_row + _decay_col);
            const double _fak2 = 2.0 * _fak;
            
            // check if distance between postions is big, then skip step   
            double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
            if ( _exparg > 30.0 ) { continue; }
             // initialize local matrix block for unnormalized cartesians
            ub::matrix<double> _ol = ub::zero_matrix<double>(_nrows,_ncols);
        
            // some helpers
            const double _fak3 = 3.0 * _fak;
            const double _fak4 = 4.0 * _fak;

            
            if ( _distsq > 0.0001  ){
            
            _pma[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
            _pma[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
            _pma[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

            _pmb[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
            _pmb[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
            _pmb[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
            }else{
                _pma[0] = 0.0;
                _pma[1] = 0.0;
                _pma[2] = 0.0;
                
                _pmb[0] = 0.0;
                _pmb[1] = 0.0;
                _pmb[2] = 0.0;
                
                
            }
            

        
        // calculate matrix elements
        _ol(0,0) = pow(4.0*_decay_row*_decay_col,0.75) * pow(_fak2,1.5)*exp(-_exparg); // s-s element
        //cout << "\t setting s-s: " << _ol(0,0) << endl;
        
        // s-p
        if ( _lmax_col > 0 ) {
          // cout << "\t setting s-p" << flush;
           _ol(0,1) = _pmb[0]*_ol(0,0); // s-px
           _ol(0,2) = _pmb[1]*_ol(0,0); // s-py
           _ol(0,3) = _pmb[2]*_ol(0,0); // s-pz
        }
        
        // p-s
        if ( _lmax_row > 0 ) {
           //cout << "\t setting p-s" << flush;

           _ol(1,0) = _pma[0]*_ol(0,0); // px-s
           _ol(2,0) = _pma[1]*_ol(0,0); // py-s
           _ol(3,0) = _pma[2]*_ol(0,0); // pz-s
        }
        
        // p-p
        if ( _lmax_row > 0 && _lmax_col > 0 ) {
           //cout << "\t setting p-p" << endl;            
           _ol(1,1) = _pma[0]*_ol(0,1) + _fak * _ol(0,0); // px-px
           _ol(1,2) = _pma[0]*_ol(0,2); // px-py
           _ol(1,3) = _pma[0]*_ol(0,3); // px-pz
           _ol(2,1) = _pma[1]*_ol(0,1); // py-px
           _ol(2,2) = _pma[1]*_ol(0,2) + _fak * _ol(0,0); // py-py
           _ol(2,3) = _pma[1]*_ol(0,3); // py-pz
           _ol(3,1) = _pma[2]*_ol(0,1); // pz-px
           _ol(3,2) = _pma[2]*_ol(0,2); // pz-py
           _ol(3,3) = _pma[2]*_ol(0,3) + _fak * _ol(0,0); // pz-pz
        }
        /* 
         * d-function elements are six cartesians, in order 
         * dxy, dxz, dyz, dxx, dyy, dzz
         * we would like to have five spherical Gaussians in order
         * dxz, dyz, dxy, d3z2-r2, dx2-y2
         */
        // s-d
        if ( _lmax_col > 1){
            //cout << "\t setting s-d" << endl;
            _ol(0,4) = _pmb[1]*_ol(0,1); // s-dxy
            _ol(0,5) = _pmb[2]*_ol(0,1); // s-dxz
            _ol(0,6) = _pmb[2]*_ol(0,2); // s-dyz
            _ol(0,7) = _pmb[0]*_ol(0,1) + _fak*_ol(0,0); // s-dxx
            _ol(0,8) = _pmb[1]*_ol(0,2) + _fak*_ol(0,0); // s-dyy
            _ol(0,9) = _pmb[2]*_ol(0,3) + _fak*_ol(0,0); // s-dzz
        }
        
        // p-d
        if ( _lmax_row > 0 && _lmax_col > 1){
            //cout << "\t setting p-d" << endl;
            
            _ol(1,4) = _pma[0]*_ol(0,4) + _fak * _ol(0,2);
            _ol(1,5) = _pma[0]*_ol(0,5) + _fak * _ol(0,3);
            _ol(1,6) = _pma[0]*_ol(0,6);
            _ol(1,7) = _pma[0]*_ol(0,7) + _fak2 * _ol(0,1);
            _ol(1,8) = _pma[0]*_ol(0,8);
            _ol(1,9) = _pma[0]*_ol(0,9);
 
            _ol(2,4) = _pma[1]*_ol(0,4) + _fak * _ol(0,1);
            _ol(2,5) = _pma[1]*_ol(0,5);
            _ol(2,6) = _pma[1]*_ol(0,6) + _fak * _ol(0,3);
            _ol(2,7) = _pma[1]*_ol(0,7);
            _ol(2,8) = _pma[1]*_ol(0,8) + _fak2 * _ol(0,2);
            _ol(2,9) = _pma[1]*_ol(0,9);

            _ol(3,4) = _pma[2]*_ol(0,4);
            _ol(3,5) = _pma[2]*_ol(0,5) + _fak * _ol(0,1);
            _ol(3,6) = _pma[2]*_ol(0,6) + _fak * _ol(0,2);
            _ol(3,7) = _pma[2]*_ol(0,7);
            _ol(3,8) = _pma[2]*_ol(0,8);
            _ol(3,9) = _pma[2]*_ol(0,9) + _fak2 * _ol(0,3);
        }

        // d-s
        if ( _lmax_row > 1){
           //cout << "\t setting d-s" << endl;
            _ol(4,0) = _pma[1]*_ol(1,0); // dxy-s
            _ol(5,0) = _pma[2]*_ol(1,0); // dxz-s
            _ol(6,0) = _pma[2]*_ol(2,0); // dyz-s
            _ol(7,0) = _pma[0]*_ol(1,0) + _fak * _ol(0,0); // dxx-s
            _ol(8,0) = _pma[1]*_ol(2,0) + _fak * _ol(0,0); // dyy-s
            _ol(9,0) = _pma[2]*_ol(3,0) + _fak * _ol(0,0); // dzz-s
        }
        
        
        // d-p
        if ( _lmax_row >1 && _lmax_col > 0){
           //cout << "\t setting d-p" << endl;

             _ol(4,1) = _pma[1]*_ol(1,1);
             _ol(5,1) = _pma[2]*_ol(1,1);
             _ol(6,1) = _pma[2]*_ol(2,1);
             
             _ol(7,1) = _pma[0]*_ol(1,1) + _fak  * ( _ol(0,1) + _ol(1,0) );
             _ol(8,1) = _pma[1]*_ol(2,1) + _fak  * _ol(0,1);
             _ol(9,1) = _pma[2]*_ol(3,1) + _fak  * _ol(0,1);

             _ol(4,2) = _pma[1]*_ol(1,2) + _fak  * _ol(1,0);
             _ol(5,2) = _pma[2]*_ol(1,2);
             _ol(6,2) = _pma[2]*_ol(2,2);
             
             _ol(7,2) = _pma[0]*_ol(1,2) + _fak  * _ol(0,2);
             _ol(8,2) = _pma[1]*_ol(2,2) + _fak  * ( _ol(0,2) + _ol (2,0) );
             _ol(9,2) = _pma[2]*_ol(3,2) + _fak  * _ol(0,2);

             _ol(4,3) = _pma[1]*_ol(1,3);
             _ol(5,3) = _pma[2]*_ol(1,3) + _fak  * _ol(1,0);
             _ol(6,3) = _pma[2]*_ol(2,3) + _fak  * _ol(2,0);
             _ol(7,3) = _pma[0]*_ol(1,3) + _fak  * _ol(0,3);
             _ol(8,3) = _pma[1]*_ol(2,3) + _fak  * _ol(0,3);
             _ol(9,3) = _pma[2]*_ol(3,3) + _fak  * ( _ol(0,3) + _ol(3,0) );
           
        }
        
        // d-d
        if ( _lmax_row > 1 && _lmax_col > 1 ){
             // cout << "\t setting d-d" << endl;
            
             _ol(4,4) = _pma[1]*_ol(1,4) + _fak * _ol(1,1);
             _ol(5,4) = _pma[2]*_ol(1,4);
             _ol(6,4) = _pma[2]*_ol(2,4);
             _ol(7,4) = _pma[0]*_ol(1,4) + _fak * (_ol(0,4) + _ol(1,2) );
             _ol(8,4) = _pma[1]*_ol(2,4) + _fak * (_ol(0,4) + _ol(2,1) );
             _ol(9,4) = _pma[2]*_ol(3,4) + _fak * _ol(0,4);

             _ol(4,5) = _pma[1]*_ol(1,5);
             _ol(5,5) = _pma[2]*_ol(1,5) + _fak * _ol(1,1);
             _ol(6,5) = _pma[2]*_ol(2,5) + _fak * _ol(2,1);
             _ol(7,5) = _pma[0]*_ol(1,5) + _fak * (_ol(0,5) + _ol(1,3) );
             _ol(8,5) = _pma[1]*_ol(2,5) + _fak * _ol(0,5);
             _ol(9,5) = _pma[2]*_ol(3,5) + _fak * (_ol(0,5) + _ol(3,1) );

             _ol(4,6) = _pma[1]*_ol(1,6) + _fak * _ol(1,3);
             _ol(5,6) = _pma[2]*_ol(1,6) + _fak * _ol(1,2);
             _ol(6,6) = _pma[2]*_ol(2,6) + _fak * _ol(2,2);
             _ol(7,6) = _pma[0]*_ol(1,6) + _fak * _ol(0,6);
             _ol(8,6) = _pma[1]*_ol(2,6) + _fak * (_ol(0,6) + _ol(2,3) );
             _ol(9,6) = _pma[2]*_ol(3,6) + _fak * (_ol(0,6) + _ol(3,2) );

             _ol(4,7) = _pma[1]*_ol(1,7);
             _ol(5,7) = _pma[2]*_ol(1,7);
             _ol(6,7) = _pma[2]*_ol(2,7);
             _ol(7,7) = _pma[0]*_ol(1,7) + _fak * (_ol(0,7) + 2.0*_ol(1,1) );
             _ol(8,7) = _pma[1]*_ol(2,7) + _fak * _ol(0,7);
             _ol(9,7) = _pma[2]*_ol(3,7) + _fak * _ol(0,7);

             _ol(4,8) = _pma[1]*_ol(1,8) + _fak2 * _ol(1,2);
             _ol(5,8) = _pma[2]*_ol(1,8);
             _ol(6,8) = _pma[2]*_ol(2,8);
             _ol(7,8) = _pma[0]*_ol(1,8) + _fak * _ol(0,8);
             _ol(8,8) = _pma[1]*_ol(2,8) + _fak * (_ol(0,8) + 2.0*_ol(2,2) );
             _ol(9,8) = _pma[2]*_ol(3,8) + _fak * _ol(0,8);

             _ol(4,9) = _pma[1]*_ol(1,9);
             _ol(5,9) = _pma[2]*_ol(1,9) + _fak2 * _ol(1,3);
             _ol(6,9) = _pma[2]*_ol(2,9) + _fak2 * _ol(2,3);
             _ol(7,9) = _pma[0]*_ol(1,9) + _fak * _ol( 0,9);
             _ol(8,9) = _pma[1]*_ol(2,9) + _fak * _ol(0,9);
             _ol(9,9) = _pma[2]*_ol(3,9) + _fak * (_ol(0,9) + 2.0*_ol(3,3) );
            
            
        }

        // s-f 
        if ( _lmax_col > 2 ){
             _ol(0,10) = _pmb[0]*_ol(0,7) + _fak2* _ol(0,1);
             _ol(0,11) = _pmb[1]*_ol(0,8) + _fak2* _ol(0,2);
             _ol(0,12) = _pmb[2]*_ol(0,9) + _fak2* _ol(0,3);
             _ol(0,13) = _pmb[0]*_ol(0,4) + _fak * _ol(0,2);
             _ol(0,14) = _pmb[1]*_ol(0,4) + _fak * _ol(0,1);
             _ol(0,15) = _pmb[0]*_ol(0,5) + _fak * _ol(0,3);
             _ol(0,16) = _pmb[2]*_ol(0,5) + _fak * _ol(0,1);
             _ol(0,17) = _pmb[1]*_ol(0,6) + _fak * _ol(0,3);
             _ol(0,18) = _pmb[2]*_ol(0,6) + _fak * _ol(0,2);
             _ol(0,19) = _pmb[2]*_ol(0,4);
        }

        // f-s
        if ( _lmax_row > 2){
             _ol(10,0) = _pma[0]*_ol(7,0) + _fak2* _ol( 1,0);
             _ol(11,0) = _pma[1]*_ol(8,0) + _fak2* _ol( 2,0);
             _ol(12,0) = _pma[2]*_ol(9,0) + _fak2* _ol( 3,0);
             _ol(13,0) = _pma[0]*_ol(4,0) + _fak * _ol( 2,0);
             _ol(14,0) = _pma[1]*_ol(4,0) + _fak * _ol( 1,0);
             _ol(15,0) = _pma[0]*_ol(5,0) + _fak * _ol( 3,0);
             _ol(16,0) = _pma[2]*_ol(5,0) + _fak * _ol( 1,0);
             _ol(17,0) = _pma[1]*_ol(6,0) + _fak * _ol( 3,0);
             _ol(18,0) = _pma[2]*_ol(6,0) + _fak * _ol( 2,0);
             _ol(19,0) = _pma[2]*_ol(4,0);
        }
        
        // p-f
        if ( _lmax_row > 0 && _lmax_col > 2 ){
             _ol( 1,10) = _pma[0]*_ol( 0,10) + _fak3* _ol( 0, 7);
             _ol( 2,10) = _pma[1]*_ol( 0,10);
             _ol( 3,10) = _pma[2]*_ol( 0,10);
             _ol( 1,11) = _pma[0]*_ol( 0,11);
             _ol( 2,11) = _pma[1]*_ol( 0,11) + _fak3* _ol( 0, 8);
             _ol( 3,11) = _pma[2]*_ol( 0,11);
             _ol( 1,12) = _pma[0]*_ol( 0,12);
             _ol( 2,12) = _pma[1]*_ol( 0,12);
             _ol( 3,12) = _pma[2]*_ol( 0,12) + _fak3* _ol( 0,9);
             _ol( 1,13) = _pma[0]*_ol( 0,13) + _fak2* _ol( 0, 4);
             _ol( 2,13) = _pma[1]*_ol( 0,13) + _fak * _ol( 0, 7);
             _ol( 3,13) = _pma[2]*_ol( 0,13);
             _ol( 1,14) = _pma[0]*_ol( 0,14) + _fak * _ol( 0, 8);
             _ol( 2,14) = _pma[1]*_ol( 0,14) + _fak2* _ol( 0, 4);
             _ol( 3,14) = _pma[2]*_ol( 0,14);
             _ol( 1,15) = _pma[0]*_ol( 0,15) + _fak2* _ol( 0, 5);
             _ol( 2,15) = _pma[1]*_ol( 0,15);
             _ol( 3,15) = _pma[2]*_ol( 0,15) + _fak * _ol( 0, 7);
             _ol( 1,16) = _pma[0]*_ol( 0,16) + _fak * _ol( 0,9);
             _ol( 2,16) = _pma[1]*_ol( 0,16);
             _ol( 3,16) = _pma[2]*_ol( 0,16) + _fak2* _ol( 0, 5);
             _ol( 1,17) = _pma[0]*_ol( 0,17);
             _ol( 2,17) = _pma[1]*_ol( 0,17) + _fak2* _ol( 0, 6);
             _ol( 3,17) = _pma[2]*_ol( 0,17) + _fak * _ol( 0, 8);
             _ol( 1,18) = _pma[0]*_ol( 0,18);
             _ol( 2,18) = _pma[1]*_ol( 0,18) + _fak * _ol( 0,9);
             _ol( 3,18) = _pma[2]*_ol( 0,18) + _fak2* _ol( 0, 6);
             _ol( 1,19) = _pma[0]*_ol( 0,19) + _fak * _ol( 0, 6);
             _ol( 2,19) = _pma[1]*_ol( 0,19) + _fak * _ol( 0, 5);
             _ol( 3,19) = _pma[2]*_ol( 0,19) + _fak * _ol( 0, 4);            
        }
   
        // f-p
        if (_lmax_row > 2 && _lmax_col > 0 ){
             _ol(13, 1) = _pma[0]*_ol( 4, 1) + _fak * (_ol( 2, 1) + _ol( 4, 0) );
             _ol(14, 1) = _pma[1]*_ol( 4, 1) + _fak * _ol( 1, 1);
             _ol(19, 1) = _pma[2]*_ol( 4, 1);
             _ol(13, 2) = _pma[0]*_ol( 4, 2) + _fak * _ol( 2, 2);
             _ol(14, 2) = _pma[1]*_ol( 4, 2) + _fak * (_ol( 1, 2) + _ol( 4, 0) );
             _ol(19, 2) = _pma[2]*_ol( 4, 2);
             _ol(13, 3) = _pma[0]*_ol( 4, 3) + _fak * _ol( 2, 3);
             _ol(14, 3) = _pma[1]*_ol( 4, 3) + _fak * _ol( 1, 3);
             _ol(19, 3) = _pma[2]*_ol( 4, 3) + _fak * _ol( 4, 0);
             _ol(15, 1) = _pma[0]*_ol( 5, 1) + _fak * (_ol( 3, 1) + _ol( 5, 0) );
             _ol(16, 1) = _pma[2]*_ol( 5, 1) + _fak * _ol( 1, 1);
             _ol(15, 2) = _pma[0]*_ol( 5, 2) + _fak * _ol( 3, 2);
             _ol(16, 2) = _pma[2]*_ol( 5, 2) + _fak * _ol( 1, 2);
             _ol(15, 3) = _pma[0]*_ol( 5, 3) + _fak * _ol( 3, 3);
             _ol(16, 3) = _pma[2]*_ol( 5, 3) + _fak * (_ol( 1, 3) + _ol( 5, 0) );
             _ol(17, 1) = _pma[1]*_ol( 6, 1) + _fak * _ol( 3, 1);
             _ol(18, 1) = _pma[2]*_ol( 6, 1) + _fak * _ol( 2, 1);
             _ol(17, 2) = _pma[1]*_ol( 6, 2) + _fak * (_ol( 3, 2) + _ol( 6, 0) );
             _ol(18, 2) = _pma[2]*_ol( 6, 2) + _fak * _ol( 2, 2);
             _ol(17, 3) = _pma[1]*_ol( 6, 3) + _fak * _ol( 3, 3);
             _ol(18, 3) = _pma[2]*_ol( 6, 3) + _fak * (_ol( 2, 3) + _ol( 6, 0) );
             _ol(10, 1) = _pma[0]*_ol( 7, 1) + _fak * (2.0*_ol( 1, 1) + _ol( 7, 0) );
             _ol(10, 2) = _pma[0]*_ol( 7, 2) + _fak2* _ol( 1, 2);
             _ol(10, 3) = _pma[0]*_ol( 7, 3) + _fak2* _ol( 1, 3);
             _ol(11, 1) = _pma[1]*_ol( 8, 1) + _fak2* _ol( 2, 1);
             _ol(11, 2) = _pma[1]*_ol( 8, 2) + _fak * (2.0*_ol( 2, 2) + _ol( 8, 0) );
             _ol(11, 3) = _pma[1]*_ol( 8, 3) + _fak2* _ol( 2, 3);
             _ol(12, 1) = _pma[2]*_ol( 9, 1) + _fak2* _ol( 3, 1);
             _ol(12, 2) = _pma[2]*_ol( 9, 2) + _fak2* _ol( 3, 2);
             _ol(12, 3) = _pma[2]*_ol( 9, 3) + _fak * (2.0*_ol( 3, 3) + _ol(9, 0) );            
        }
        
        // d-f
        if ( _lmax_row > 1 && _lmax_col >2 ){
             _ol( 7,10) = _pma[0]*_ol( 1,10) + _fak * (_ol( 0,10) + 3.0*_ol( 1, 7) );
             _ol( 4,10) = _pma[1]*_ol( 1,10);
             _ol( 5,10) = _pma[2]*_ol( 1,10);
             _ol( 7,11) = _pma[0]*_ol( 1,11) + _fak * _ol( 0,1);
             _ol( 4,11) = _pma[1]*_ol( 1,11) + _fak3* _ol( 1, 8);
             _ol( 5,11) = _pma[2]*_ol( 1,11);
             _ol( 7,12) = _pma[0]*_ol( 1,12) + _fak * _ol( 0,12);
             _ol( 4,12) = _pma[1]*_ol( 1,12);
             _ol( 5,12) = _pma[2]*_ol( 1,12) + _fak3* _ol( 1,9);
             _ol( 7,13) = _pma[0]*_ol( 1,13) + _fak * (_ol( 0,13) + 2.0*_ol( 1, 4) );
             _ol( 4,13) = _pma[1]*_ol( 1,13) + _fak * _ol( 1, 7);
             _ol( 5,13) = _pma[2]*_ol( 1,13);
             _ol( 7,14) = _pma[0]*_ol( 1,14) + _fak * (_ol( 0,14) + _ol( 1, 8) );
             _ol( 4,14) = _pma[1]*_ol( 1,14) + _fak2* _ol( 1, 4);
             _ol( 5,14) = _pma[2]*_ol( 1,14);
             _ol( 7,15) = _pma[0]*_ol( 1,15) + _fak * (_ol( 0,15) + 2.0*_ol( 1, 5) );
             _ol( 4,15) = _pma[1]*_ol( 1,15);
             _ol( 5,15) = _pma[2]*_ol( 1,15) + _fak * _ol( 1, 7);
             _ol( 7,16) = _pma[0]*_ol( 1,16) + _fak * (_ol( 0,16) + _ol( 1,9) );
             _ol( 4,16) = _pma[1]*_ol( 1,16);
             _ol( 5,16) = _pma[2]*_ol( 1,16) + _fak2* _ol( 1, 5);
             _ol( 7,17) = _pma[0]*_ol( 1,17) + _fak * _ol( 0,17);
             _ol( 4,17) = _pma[1]*_ol( 1,17) + _fak2* _ol( 1, 6);
             _ol( 5,17) = _pma[2]*_ol( 1,17) + _fak * _ol( 1, 8);
             _ol( 7,18) = _pma[0]*_ol( 1,18) + _fak * _ol( 0,18);
             _ol( 4,18) = _pma[1]*_ol( 1,18) + _fak * _ol( 1,9);
             _ol( 5,18) = _pma[2]*_ol( 1,18) + _fak2* _ol( 1, 6);
             _ol( 7,19) = _pma[0]*_ol( 1,19) + _fak * (_ol( 0,19) + _ol( 1, 6) );
             _ol( 4,19) = _pma[1]*_ol( 1,19) + _fak * _ol( 1, 5);
             _ol( 5,19) = _pma[2]*_ol( 1,19) + _fak * _ol( 1, 4);
             _ol( 8,10) = _pma[1]*_ol( 2,10) + _fak * _ol( 0,10);
             _ol( 6,10) = _pma[2]*_ol( 2,10);
             _ol( 8,11) = _pma[1]*_ol( 2,11) + _fak * (_ol( 0,11) + 3.0*_ol( 2, 8) );
             _ol( 6,11) = _pma[2]*_ol( 2,11);
             _ol( 8,12) = _pma[1]*_ol( 2,12) + _fak * _ol( 0,12);
             _ol( 6,12) = _pma[2]*_ol( 2,12) + _fak3* _ol( 2,9);
             _ol( 8,13) = _pma[1]*_ol( 2,13) + _fak * (_ol( 0,13) + _ol( 2, 7) );
             _ol( 6,13) = _pma[2]*_ol( 2,13);
             _ol( 8,14) = _pma[1]*_ol( 2,14) + _fak * (_ol( 0,14) + 2.0*_ol( 2, 4) );
             _ol( 6,14) = _pma[2]*_ol( 2,14);
             _ol( 8,15) = _pma[1]*_ol( 2,15) + _fak * _ol( 0,15);
             _ol( 6,15) = _pma[2]*_ol( 2,15) + _fak * _ol( 2, 7);
             _ol( 8,16) = _pma[1]*_ol( 2,16) + _fak * _ol( 0,16);
             _ol( 6,16) = _pma[2]*_ol( 2,16) + _fak2* _ol( 2, 5);
             _ol( 8,17) = _pma[1]*_ol( 2,17) + _fak * (_ol( 0,17) + 2.0*_ol( 2, 6) );
             _ol( 6,17) = _pma[2]*_ol( 2,17) + _fak * _ol( 2, 8);
             _ol( 8,18) = _pma[1]*_ol( 2,18) + _fak * (_ol( 0,18) + _ol( 2,9) );
             _ol( 6,18) = _pma[2]*_ol( 2,18) + _fak2* _ol( 2, 6);
             _ol( 8,19) = _pma[1]*_ol( 2,19) + _fak * (_ol( 0,19) + _ol( 2, 5) );
             _ol( 6,19) = _pma[2]*_ol( 2,19) + _fak * _ol( 2, 4);
             _ol(9,10) = _pma[2]*_ol( 3,10) + _fak * _ol( 0,10);
             _ol(9,11) = _pma[2]*_ol( 3,11) + _fak * _ol( 0,11);
             _ol(9,12) = _pma[2]*_ol( 3,12) + _fak * (_ol( 0,12) + 3.0*_ol( 3,9) );
             _ol(9,13) = _pma[2]*_ol( 3,13) + _fak * _ol( 0,13);
             _ol(9,14) = _pma[2]*_ol( 3,14) + _fak * _ol( 0,14);
             _ol(9,15) = _pma[2]*_ol( 3,15) + _fak * (_ol( 0,15) + _ol( 3, 7) );
             _ol(9,16) = _pma[2]*_ol( 3,16) + _fak * (_ol( 0,16) + 2.0*_ol( 3, 5) );
             _ol(9,17) = _pma[2]*_ol( 3,17) + _fak * (_ol( 0,17) + _ol( 3, 8) );
             _ol(9,18) = _pma[2]*_ol( 3,18) + _fak * (_ol( 0,18) + 2.0*_ol( 3, 5) );
             _ol(9,19) = _pma[2]*_ol( 3,19) + _fak * (_ol( 0,19) + _ol( 3, 4) );
        }
        // f-d
        if ( _lmax_row > 2 && _lmax_col > 1 ){
             _ol(13, 4) = _pma[0]*_ol( 4, 4) + _fak * (_ol( 2, 4) + _ol( 4, 2) );
             _ol(14, 4) = _pma[1]*_ol( 4, 4) + _fak * (_ol( 1, 4) + _ol( 4, 1) );
             _ol(19, 4) = _pma[2]*_ol( 4, 4);
             _ol(13, 5) = _pma[0]*_ol( 4, 5) + _fak * (_ol( 2, 5) + _ol( 4, 3) );
             _ol(14, 5) = _pma[1]*_ol( 4, 5) + _fak * _ol( 1, 5);
             _ol(19, 5) = _pma[2]*_ol( 4, 5) + _fak * _ol( 4, 1);
             _ol(13, 6) = _pma[0]*_ol( 4, 6) + _fak * _ol( 2, 6);
             _ol(14, 6) = _pma[1]*_ol( 4, 6) + _fak * (_ol( 1, 6) + _ol( 4, 3) );
             _ol(19, 6) = _pma[2]*_ol( 4, 6) + _fak * _ol( 4, 2);
             _ol(13, 7) = _pma[0]*_ol( 4, 7) + _fak * (_ol( 2, 7) + 2.0*_ol( 4, 1) );
             _ol(14, 7) = _pma[1]*_ol( 4, 7) + _fak * _ol( 1, 7);
             _ol(19, 7) = _pma[2]*_ol( 4, 7);
             _ol(13, 8) = _pma[0]*_ol( 4, 8) + _fak * _ol( 2, 8);
             _ol(14, 8) = _pma[1]*_ol( 4, 8) + _fak * (_ol( 1, 8) + 2.0*_ol( 4, 2) );
             _ol(19, 8) = _pma[2]*_ol( 4, 8);
             _ol(13,9) = _pma[0]*_ol( 4,9) + _fak * _ol( 2,9);
             _ol(14,9) = _pma[1]*_ol( 4,9) + _fak * _ol( 1,9);
             _ol(19,9) = _pma[2]*_ol( 4,9) + _fak2* _ol( 4, 3);
             _ol(15, 4) = _pma[0]*_ol( 5, 4) + _fak * (_ol( 3, 4) + _ol( 5, 2) );
             _ol(16, 4) = _pma[2]*_ol( 5, 4) + _fak * _ol( 1, 4);
             _ol(15, 5) = _pma[0]*_ol( 5, 5) + _fak * (_ol( 3, 5) + _ol( 5, 3) );
             _ol(16, 5) = _pma[2]*_ol( 5, 5) + _fak * (_ol( 1, 5) + _ol( 5, 1) );
             _ol(15, 6) = _pma[0]*_ol( 5, 6) + _fak * _ol( 3, 6);
             _ol(16, 6) = _pma[2]*_ol( 5, 6) + _fak * (_ol( 1, 6) + _ol( 5, 2) );
             _ol(15, 7) = _pma[0]*_ol( 5, 7) + _fak * (_ol( 3, 7) + 2.0*_ol( 5, 1) );
             _ol(16, 7) = _pma[2]*_ol( 5, 7) + _fak * _ol( 1, 7);
             _ol(15, 8) = _pma[0]*_ol( 5, 8) + _fak * _ol( 3, 8);
             _ol(16, 8) = _pma[2]*_ol( 5, 8) + _fak * _ol( 1, 8);
             _ol(15,9) = _pma[0]*_ol( 5,9) + _fak * _ol( 3,9);
             _ol(16,9) = _pma[2]*_ol( 5,9) + _fak * (_ol( 1,9) + 2.0*_ol( 5, 3) );
             _ol(17, 4) = _pma[1]*_ol( 6, 4) + _fak * (_ol( 3, 4) + _ol( 6, 1) );
             _ol(18, 4) = _pma[2]*_ol( 6, 4) + _fak * _ol( 2, 4);
             _ol(17, 5) = _pma[1]*_ol( 6, 5) + _fak * _ol( 3, 5);
             _ol(18, 5) = _pma[2]*_ol( 6, 5) + _fak * (_ol( 2, 5) + _ol( 6, 1) );
             _ol(17, 6) = _pma[1]*_ol( 6, 6) + _fak * (_ol( 3, 6) + _ol( 6, 3) );
             _ol(18, 6) = _pma[2]*_ol( 6, 6) + _fak * (_ol( 2, 5) + _ol( 5, 2) );
             _ol(17, 7) = _pma[1]*_ol( 6, 7) + _fak * _ol( 3, 7);
             _ol(18, 7) = _pma[2]*_ol( 6, 7) + _fak * _ol( 2, 7);
             _ol(17, 8) = _pma[1]*_ol( 6, 8) + _fak * (_ol( 3, 8) + 2.0*_ol( 6, 2) );
             _ol(18, 8) = _pma[2]*_ol( 6, 8) + _fak * _ol( 2, 8);
             _ol(17,9) = _pma[1]*_ol( 6,9) + _fak * _ol( 3,9);
             _ol(18,9) = _pma[2]*_ol( 6,9) + _fak * (_ol( 2,9) + 2.0*_ol( 6, 3) );
             _ol(10, 4) = _pma[0]*_ol( 7, 4) + _fak * (2.0*_ol( 1, 4) + _ol( 7, 2) );
             _ol(10, 5) = _pma[0]*_ol( 7, 5) + _fak * (2.0*_ol( 1, 5) + _ol( 7, 3) );
             _ol(10, 6) = _pma[0]*_ol( 7, 6) + _fak2* _ol( 1, 6);
             _ol(10, 7) = _pma[0]*_ol( 7, 7) + _fak * (2.0*_ol( 1, 7) + 2.0*_ol( 7, 1));
             _ol(10, 8) = _pma[0]*_ol( 7, 8) + _fak2* _ol( 1, 8);
             _ol(10,9) = _pma[0]*_ol( 7,9) + _fak2* _ol( 1,9);
             _ol(11, 4) = _pma[1]*_ol( 8, 4) + _fak * (2.0*_ol( 2, 4) + _ol( 8, 1) );
             _ol(11, 5) = _pma[1]*_ol( 8, 5) + _fak2* _ol( 2, 5);
             _ol(11, 6) = _pma[1]*_ol( 8, 6) + _fak * (2.0*_ol( 2, 6) + _ol( 8, 3) );
             _ol(11, 7) = _pma[1]*_ol( 8, 7) + _fak2* _ol( 2, 7);
             _ol(11, 8) = _pma[1]*_ol( 8, 8) + _fak * (2.0*_ol( 2, 8) + 2.0*_ol( 8, 2));
             _ol(11,9) = _pma[1]*_ol( 8,9) + _fak2* _ol( 2,9);
             _ol(12, 4) = _pma[2]*_ol(9, 4) + _fak2* _ol( 3, 4);
             _ol(12, 5) = _pma[2]*_ol(9, 5) + _fak * (2.0*_ol( 3, 5) + _ol(9, 1) );
             _ol(12, 6) = _pma[2]*_ol(9, 6) + _fak * (2.0*_ol( 3, 6) + _ol(9, 2) );
             _ol(12, 7) = _pma[2]*_ol(9, 7) + _fak2* _ol( 3, 7);
             _ol(12, 8) = _pma[2]*_ol(9, 8) + _fak2* _ol( 3, 8);
             _ol(12,9) = _pma[2]*_ol(9,9) + _fak * (2.0*_ol( 3,9) + 2.0*_ol(9, 3));
        }
        // f-f
        if ( _lmax_row > 2 && _lmax_col > 2 ){
             _ol(13,10) = _pma[0]*_ol( 4,10) + _fak * (_ol( 2,10) + 3.0*_ol( 4, 7) );
             _ol(14,10) = _pma[1]*_ol( 4,10) + _fak * _ol( 1,10);
             _ol(19,10) = _pma[2]*_ol( 4,10);
             _ol(13,11) = _pma[0]*_ol( 4,11) + _fak * _ol( 2,11);
             _ol(14,11) = _pma[1]*_ol( 4,11) + _fak * (_ol( 1,11) + 3.0*_ol( 4, 8) );
             _ol(19,11) = _pma[2]*_ol( 4,11);
             _ol(13,12) = _pma[0]*_ol( 4,12) + _fak * _ol( 2,12);
             _ol(14,12) = _pma[1]*_ol( 4,12) + _fak * _ol( 1,12);
             _ol(19,12) = _pma[2]*_ol( 4,12) + _fak3* _ol( 4,9);
             _ol(13,13) = _pma[0]*_ol( 4,13) + _fak * (_ol( 2,13) + 2.0*_ol( 4,4) );
             _ol(14,13) = _pma[1]*_ol( 4,13) + _fak * (_ol( 1,13) + _ol( 4, 7) );
             _ol(19,13) = _pma[2]*_ol( 4,13);
             _ol(13,14) = _pma[0]*_ol( 4,14) + _fak * (_ol( 2,14) + _ol( 4, 8) );
             _ol(14,14) = _pma[1]*_ol( 4,14) + _fak * (_ol( 1,14) + 2.0*_ol( 4, 4) );
             _ol(19,14) = _pma[2]*_ol( 4,14);
             _ol(13,15) = _pma[0]*_ol( 4,15) + _fak * (_ol( 2,15) + 2.0*_ol( 4, 5) );
             _ol(14,15) = _pma[1]*_ol( 4,15) + _fak * _ol( 1,15);
             _ol(19,15) = _pma[2]*_ol( 4,15) + _fak * _ol( 4, 7);
             _ol(13,16) = _pma[0]*_ol( 4,16) + _fak * (_ol( 2,16) + _ol( 4,9) );
             _ol(14,16) = _pma[1]*_ol( 4,16) + _fak * _ol( 1,16);
             _ol(19,16) = _pma[2]*_ol( 4,16) + _fak2* _ol( 4, 5);
             _ol(13,17) = _pma[0]*_ol( 4,17) + _fak * _ol( 2,17);
             _ol(14,17) = _pma[1]*_ol( 4,17) + _fak * (_ol( 1,17) + 2.0*_ol( 4, 6) );
             _ol(19,17) = _pma[2]*_ol( 4,17) + _fak * _ol( 4, 8);
             _ol(13,18) = _pma[0]*_ol( 4,18) + _fak * _ol( 2,18);
             _ol(14,18) = _pma[1]*_ol( 4,18) + _fak * (_ol( 1,18) + _ol( 4,9) );
             _ol(19,18) = _pma[2]*_ol( 4,18) + _fak2* _ol( 4, 6);
             _ol(13,19) = _pma[0]*_ol( 4,19) + _fak * (_ol( 2,19) + _ol( 4, 6) );
             _ol(14,19) = _pma[1]*_ol( 4,19) + _fak * (_ol( 1,19) + _ol( 4, 5) );
             _ol(19,19) = _pma[2]*_ol( 4,19) + _fak * _ol( 4, 4);
             _ol(15,10) = _pma[0]*_ol( 5,10) + _fak * (_ol( 3,10) + 3.0*_ol( 5, 7) );
             _ol(16,10) = _pma[2]*_ol( 5,10) + _fak * _ol( 1,10);
             _ol(15,11) = _pma[0]*_ol( 5,11) + _fak * _ol( 3,11);
             _ol(16,11) = _pma[2]*_ol( 5,11) + _fak * _ol( 1,11);
             _ol(15,12) = _pma[0]*_ol( 5,12) + _fak * _ol( 3,12);
             _ol(16,12) = _pma[2]*_ol( 5,12) + _fak * (_ol( 1,12) + 3.0*_ol( 5,9) );
             _ol(15,13) = _pma[0]*_ol( 5,13) + _fak * (_ol( 3,13) + 2.0*_ol( 5, 4) );
             _ol(16,13) = _pma[2]*_ol( 5,13) + _fak * _ol( 1,13);
             _ol(15,14) = _pma[0]*_ol( 5,14) + _fak * (_ol( 3,14) + _ol( 5, 8) );
             _ol(16,14) = _pma[2]*_ol( 5,14) + _fak * _ol( 1,14);
             _ol(15,15) = _pma[0]*_ol( 5,15) + _fak * (_ol( 3,15) + 2.0*_ol( 5, 5) );
             _ol(16,15) = _pma[2]*_ol( 5,15) + _fak * (_ol( 1,15) + _ol( 5, 7) );
             _ol(15,16) = _pma[0]*_ol( 5,16) + _fak * (_ol( 3,16) + _ol( 5,9) );
             _ol(16,16) = _pma[2]*_ol( 5,16) + _fak * (_ol( 1,16) + 2.0*_ol( 5, 5) );
             _ol(15,17) = _pma[0]*_ol( 5,17) + _fak * _ol( 3,17);
             _ol(16,17) = _pma[2]*_ol( 5,17) + _fak * (_ol( 1,17) + _ol( 5, 8) );
             _ol(15,18) = _pma[0]*_ol( 5,18) + _fak * _ol( 3,18);
             _ol(16,18) = _pma[2]*_ol( 5,18) + _fak * (_ol( 1,18) + 2.0*_ol( 5, 6) );
             _ol(15,19) = _pma[0]*_ol( 5,19) + _fak * (_ol( 3,19) + _ol( 5, 6) );
             _ol(16,19) = _pma[2]*_ol( 5,19) + _fak * (_ol( 1,19) + _ol( 5, 4) );
             _ol(17,10) = _pma[1]*_ol( 6,10) + _fak * _ol( 3,10);
             _ol(18,10) = _pma[2]*_ol( 6,10) + _fak * _ol( 2,10);
             _ol(17,11) = _pma[1]*_ol( 6,11) + _fak * (_ol( 3,11) + 3.0*_ol( 6, 8) );
             _ol(18,11) = _pma[2]*_ol( 6,11) + _fak * _ol( 2,11);
             _ol(17,12) = _pma[1]*_ol( 6,12) + _fak * _ol( 3,12);
             _ol(18,12) = _pma[2]*_ol( 6,12) + _fak * (_ol( 2,12) + 3.0*_ol( 6,9) );
             _ol(17,13) = _pma[1]*_ol( 6,13) + _fak * (_ol( 3,13) + _ol( 6, 7) );
             _ol(18,13) = _pma[2]*_ol( 6,13) + _fak * _ol( 2,13);
             _ol(17,14) = _pma[1]*_ol( 6,14) + _fak * (_ol( 3,14) + 2.0*_ol( 6, 4) );
             _ol(18,14) = _pma[2]*_ol( 6,14) + _fak * _ol( 2,14);
             _ol(17,15) = _pma[1]*_ol( 6,15) + _fak * _ol( 3,15);
             _ol(18,15) = _pma[2]*_ol( 6,15) + _fak * (_ol( 2,15) + _ol( 6, 7) );
             _ol(17,16) = _pma[1]*_ol( 6,16) + _fak * _ol( 3,16);
             _ol(18,16) = _pma[2]*_ol( 6,16) + _fak * (_ol( 2,16) + 2.0*_ol( 6, 5) );
             _ol(17,17) = _pma[1]*_ol( 6,17) + _fak * (_ol( 3,17) + 2.0*_ol( 6, 6) );
             _ol(18,17) = _pma[2]*_ol( 6,17) + _fak * (_ol( 2,17) + _ol( 6, 8) );
             _ol(17,18) = _pma[1]*_ol( 6,18) + _fak * (_ol( 3,18) + _ol( 6,9) );
             _ol(18,18) = _pma[2]*_ol( 6,18) + _fak * (_ol( 2,18) + 2.0*_ol( 6, 6) );
             _ol(17,19) = _pma[1]*_ol( 6,19) + _fak * (_ol( 3,19) + _ol( 6, 5) );
             _ol(18,19) = _pma[2]*_ol( 6,19) + _fak * (_ol( 2,19) + _ol( 6, 4) );
             _ol(10,10) = _pma[0]*_ol( 7,10) + _fak * (2.0*_ol( 1,10) + 3.0*_ol( 7, 7));
             _ol(10,11) = _pma[0]*_ol( 7,11) + _fak2* _ol( 1,11);
             _ol(10,12) = _pma[0]*_ol( 7,12) + _fak2* _ol( 1,12);
             _ol(10,13) = _pma[0]*_ol( 7,13) + _fak * (2.0*_ol( 1,13) + 2.0*_ol( 7, 4));
             _ol(10,14) = _pma[0]*_ol( 7,14) + _fak * (2.0*_ol( 1,14) + _ol( 7, 8) );
             _ol(10,15) = _pma[0]*_ol( 7,15) + _fak * (2.0*_ol( 1,15) + 2.0*_ol( 7, 5));
             _ol(10,16) = _pma[0]*_ol( 7,16) + _fak * (2.0*_ol( 1,16) + _ol( 7,9) );
             _ol(10,17) = _pma[0]*_ol( 7,17) + _fak2* _ol( 1,17);
             _ol(10,18) = _pma[0]*_ol( 7,18) + _fak2* _ol( 1,18);
             _ol(10,19) = _pma[0]*_ol( 7,19) + _fak * (2.0*_ol( 1,19) + _ol( 7, 6) );
             _ol(11,10) = _pma[1]*_ol( 8,10) + _fak2* _ol( 2,10);
             _ol(11,11) = _pma[1]*_ol( 8,11) + _fak * (2.0*_ol( 2,11) + 3.0*_ol( 8, 8));
             _ol(11,12) = _pma[1]*_ol( 8,12) + _fak2* _ol( 2,12);
             _ol(11,13) = _pma[1]*_ol( 8,13) + _fak * (2.0*_ol( 2,13) + _ol( 8, 7) );
             _ol(11,14) = _pma[1]*_ol( 8,14) + _fak * (2.0*_ol( 2,14) + 2.0*_ol( 8, 4));
             _ol(11,15) = _pma[1]*_ol( 8,15) + _fak2* _ol( 2,15);
             _ol(11,16) = _pma[1]*_ol( 8,16) + _fak2* _ol( 2,16);
             _ol(11,17) = _pma[1]*_ol( 8,17) + _fak * (2.0*_ol( 2,17) + 2.0*_ol( 8, 6));
             _ol(11,18) = _pma[1]*_ol( 8,18) + _fak * (2.0*_ol( 2,18) + _ol( 8,9) );
             _ol(11,19) = _pma[1]*_ol( 8,19) + _fak * (2.0*_ol( 2,19) + _ol( 8, 5) );
             _ol(12,10) = _pma[2]*_ol(9,10) + _fak2* _ol( 3,10);
             _ol(12,11) = _pma[2]*_ol(9,11) + _fak2* _ol( 3,11);
             _ol(12,12) = _pma[2]*_ol(9,12) + _fak * (2.0*_ol( 3,12) + 3.0*_ol(9,9));
             _ol(12,13) = _pma[2]*_ol(9,13) + _fak2* _ol( 3,13);
             _ol(12,14) = _pma[2]*_ol(9,14) + _fak2* _ol( 3,14);
             _ol(12,15) = _pma[2]*_ol(9,15) + _fak * (2.0*_ol( 3,15) + _ol(9, 7) );
             _ol(12,16) = _pma[2]*_ol(9,16) + _fak * (2.0*_ol( 3,16) + 2.0*_ol(9, 5));
             _ol(12,17) = _pma[2]*_ol(9,17) + _fak * (2.0*_ol( 3,17) + _ol(9, 8) );
             _ol(12,18) = _pma[2]*_ol(9,18) + _fak * (2.0*_ol( 3,18) + 2.0*_ol(9, 6));
             _ol(12,19) = _pma[2]*_ol(9,19) + _fak * (2.0*_ol( 3,19) + _ol(9, 4) );
        }
        // s-g
        
     if( (_lmax_col)>3) {
     _ol(0,20) = _pmb[0]*_ol(0,10) + _fak3* _ol(0,7);
     _ol(0,23) = _pmb[1]*_ol(0,10);
     _ol(0,25) = _pmb[2]*_ol(0,10);
     _ol(0,24) = _pmb[0]*_ol(0,11);
     _ol(0,21) = _pmb[1]*_ol(0,11) + _fak3* _ol(0,8);
     _ol(0,27) = _pmb[2]*_ol(0,11);
     _ol(0,26) = _pmb[0]*_ol(0,12);
     _ol(0,28) = _pmb[1]*_ol(0,12);
     _ol(0,22) = _pmb[2]*_ol(0,12) + _fak3* _ol(0,9);
     _ol(0,31) = _pmb[1]*_ol(0,13) + _fak * _ol(0,7);
     _ol(0,32) = _pmb[2]*_ol(0,13);
     _ol(0,33) = _pmb[2]*_ol(0,14);
     _ol(0,29) = _pmb[2]*_ol(0,15) + _fak * _ol(0,7);
     _ol(0,34) = _pmb[1]*_ol(0,16);
     _ol(0,30) = _pmb[2]*_ol(0,17) + _fak * _ol(0,8);
    }
        // g-s
        
         if((_lmax_row)>3 ) {
     _ol(20,0) = _pma[0]*_ol(10,0) + _fak3* _ol(7,0);
     _ol(23,0) = _pma[1]*_ol(10,0);
     _ol(25,0) = _pma[2]*_ol(10,0);
     _ol(24,0) = _pma[0]*_ol(11,0);
     _ol(21,0) = _pma[1]*_ol(11,0) + _fak3* _ol(8,0);
     _ol(27,0) = _pma[2]*_ol(11,0);
     _ol(26,0) = _pma[0]*_ol(12,0);
     _ol(28,0) = _pma[1]*_ol(12,0);
     _ol(22,0) = _pma[2]*_ol(12,0) + _fak3* _ol(9,0);
     _ol(31,0) = _pma[1]*_ol(13,0) + _fak * _ol(7,0);
     _ol(32,0) = _pma[2]*_ol(13,0);
     _ol(33,0) = _pma[2]*_ol(14,0);
     _ol(29,0) = _pma[2]*_ol(15,0) + _fak * _ol(7,0);
     _ol(34,0) = _pma[1]*_ol(16,0);
     _ol(30,0) = _pma[2]*_ol(17,0) + _fak * _ol(8,0);
  }
        // p-g
        
    if((_lmax_row)>0 && (_lmax_col)>3) {
     _ol(1,20) = _pma[0]*_ol(0,20) + _fak4* _ol(0,10);
     _ol(2,20) = _pma[1]*_ol(0,20);
     _ol(3,20) = _pma[2]*_ol(0,20);
     _ol(1,21) = _pma[0]*_ol(0,21);
     _ol(2,21) = _pma[1]*_ol(0,21) + _fak4* _ol(0,11);
     _ol(3,21) = _pma[2]*_ol(0,21);
     _ol(1,22) = _pma[0]*_ol(0,22);
     _ol(2,22) = _pma[1]*_ol(0,22);
     _ol(3,22) = _pma[2]*_ol(0,22) + _fak4* _ol(0,12);
     _ol(1,23) = _pma[0]*_ol(0,23) + _fak3* _ol(0,13);
     _ol(2,23) = _pma[1]*_ol(0,23) + _fak * _ol(0,10);
     _ol(3,23) = _pma[2]*_ol(0,23);
     _ol(1,24) = _pma[0]*_ol(0,24) + _fak * _ol(0,11);
     _ol(2,24) = _pma[1]*_ol(0,24) + _fak3* _ol(0,14);
     _ol(3,24) = _pma[2]*_ol(0,24);
     _ol(1,25) = _pma[0]*_ol(0,25) + _fak3* _ol(0,15);
     _ol(2,25) = _pma[1]*_ol(0,25);
     _ol(3,25) = _pma[2]*_ol(0,25) + _fak * _ol(0,10);
     _ol(1,26) = _pma[0]*_ol(0,26) + _fak * _ol(0,12);
     _ol(2,26) = _pma[1]*_ol(0,26);
     _ol(3,26) = _pma[2]*_ol(0,26) + _fak3* _ol(0,16);
     _ol(1,27) = _pma[0]*_ol(0,27);
     _ol(2,27) = _pma[1]*_ol(0,27) + _fak3* _ol(0,17);
     _ol(3,27) = _pma[2]*_ol(0,27) + _fak * _ol(0,11);
     _ol(1,28) = _pma[0]*_ol(0,28);
     _ol(2,28) = _pma[1]*_ol(0,28) + _fak * _ol(0,12);
     _ol(3,28) = _pma[2]*_ol(0,28) + _fak3* _ol(0,18);
     _ol(1,29) = _pma[0]*_ol(0,29) + _fak2* _ol(0,16);
     _ol(2,29) = _pma[1]*_ol(0,29);
     _ol(3,29) = _pma[2]*_ol(0,29) + _fak2* _ol(0,15);
     _ol(1,30) = _pma[0]*_ol(0,30);
     _ol(2,30) = _pma[1]*_ol(0,30) + _fak2* _ol(0,18);
     _ol(3,30) = _pma[2]*_ol(0,30) + _fak2* _ol(0,17);
     _ol(1,31) = _pma[0]*_ol(0,31) + _fak2* _ol(0,14);
     _ol(2,31) = _pma[1]*_ol(0,31) + _fak2* _ol(0,13);
     _ol(3,31) = _pma[2]*_ol(0,31);
     _ol(1,32) = _pma[0]*_ol(0,32) + _fak2* _ol(0,19);
     _ol(2,32) = _pma[1]*_ol(0,32) + _fak * _ol(0,15);
     _ol(3,32) = _pma[2]*_ol(0,32) + _fak * _ol(0,13);
     _ol(1,33) = _pma[0]*_ol(0,33) + _fak * _ol(0,17);
     _ol(2,33) = _pma[1]*_ol(0,33) + _fak2* _ol(0,19);
     _ol(3,33) = _pma[2]*_ol(0,33) + _fak * _ol(0,14);
     _ol(1,34) = _pma[0]*_ol(0,34) + _fak * _ol(0,18);
     _ol(2,34) = _pma[1]*_ol(0,34) + _fak * _ol(0,16);
     _ol(3,34) = _pma[2]*_ol(0,34) + _fak2* _ol(0,19);
  }    
        
        // g-p
        
        
         
  if((_lmax_row)>3 && (_lmax_col)>0) {
     _ol(20,1) = _pma[0]*_ol(10,1) + _fak * (3.0*_ol(7,1) + _ol(10,0) );
     _ol(23,1) = _pma[1]*_ol(10,1);
     _ol(25,1) = _pma[2]*_ol(10,1);
     _ol(20,2) = _pma[0]*_ol(10,2) + _fak3* _ol(7,2);
     _ol(23,2) = _pma[1]*_ol(10,2) + _fak * _ol(10,0);
     _ol(25,2) = _pma[2]*_ol(10,2);
     _ol(20,3) = _pma[0]*_ol(10,3) + _fak3* _ol(7,3);
     _ol(23,3) = _pma[1]*_ol(10,3);
     _ol(25,3) = _pma[2]*_ol(10,3) + _fak * _ol(10,0);
     _ol(24,1) = _pma[0]*_ol(11,1) + _fak * _ol(11,0);
     _ol(21,1) = _pma[1]*_ol(11,1) + _fak3* _ol(8,1);
     _ol(27,1) = _pma[2]*_ol(11,1);
     _ol(24,2) = _pma[0]*_ol(11,2);
     _ol(21,2) = _pma[1]*_ol(11,2) + _fak * (3.0*_ol(8,2) + _ol(11,0) );
     _ol(27,2) = _pma[2]*_ol(11,2);
     _ol(24,3) = _pma[0]*_ol(11,3);
     _ol(21,3) = _pma[1]*_ol(11,3) + _fak3* _ol(8,3);
     _ol(27,3) = _pma[2]*_ol(11,3) + _fak * _ol(11,0);
     _ol(26,1) = _pma[0]*_ol(12,1) + _fak * _ol(12,0);
     _ol(28,1) = _pma[1]*_ol(12,1);
     _ol(22,1) = _pma[2]*_ol(12,1) + _fak3* _ol(9,1);
     _ol(26,2) = _pma[0]*_ol(12,2);
     _ol(28,2) = _pma[1]*_ol(12,2) + _fak * _ol(12,0);
     _ol(22,2) = _pma[2]*_ol(12,2) + _fak3* _ol(9,2);
     _ol(26,3) = _pma[0]*_ol(12,3);
     _ol(28,3) = _pma[1]*_ol(12,3);
     _ol(22,3) = _pma[2]*_ol(12,3) + _fak * (3.0*_ol(9,3) + _ol(12,0) );
     _ol(31,1) = _pma[1]*_ol(13,1) + _fak * _ol(7,1);
     _ol(32,1) = _pma[2]*_ol(13,1);
     _ol(31,2) = _pma[1]*_ol(13,2) + _fak * (_ol(7,2) + _ol(13,0) );
     _ol(32,2) = _pma[2]*_ol(13,2);
     _ol(31,3) = _pma[1]*_ol(13,3) + _fak * _ol(7,3);
     _ol(32,3) = _pma[2]*_ol(13,3) + _fak * _ol(13,0);
     _ol(33,1) = _pma[2]*_ol(14,1);
     _ol(33,2) = _pma[2]*_ol(14,2);
     _ol(33,3) = _pma[2]*_ol(14,3) + _fak * _ol(14,0);
     _ol(29,1) = _pma[2]*_ol(15,1) + _fak * _ol(7,1);
     _ol(29,2) = _pma[2]*_ol(15,2) + _fak * _ol(7,2);
     _ol(29,3) = _pma[2]*_ol(15,3) + _fak * (_ol(7,3) + _ol(15,0) );
     _ol(34,1) = _pma[1]*_ol(16,1);
     _ol(34,2) = _pma[1]*_ol(16,2) + _fak * _ol(16,0);
     _ol(34,3) = _pma[1]*_ol(16,3);
     _ol(30,1) = _pma[2]*_ol(17,1) + _fak * _ol(8,1);
     _ol(30,2) = _pma[2]*_ol(17,2) + _fak * _ol(8,2);
     _ol(30,3) = _pma[2]*_ol(17,3) + _fak * (_ol(8,3) + _ol(17,0) );
  }
        // d-g
        
        if((_lmax_row)>1 && (_lmax_col)>3) {
     _ol(7,20) = _pma[0]*_ol(1,20) + _fak * (_ol(0,20) + 4.0*_ol(1,10) );
     _ol(4,20) = _pma[1]*_ol(1,20);
     _ol(5,20) = _pma[2]*_ol(1,20);
     _ol(7,21) = _pma[0]*_ol(1,21) + _fak * _ol(0,21);
     _ol(4,21) = _pma[1]*_ol(1,21) + _fak4* _ol(1,11);
     _ol(5,21) = _pma[2]*_ol(1,21);
     _ol(7,22) = _pma[0]*_ol(1,22) + _fak * _ol(0,22);
     _ol(4,22) = _pma[1]*_ol(1,22);
     _ol(5,22) = _pma[2]*_ol(1,22) + _fak4* _ol(1,12);
     _ol(7,23) = _pma[0]*_ol(1,23) + _fak * (_ol(0,23) + 3.0*_ol(1,13) );
     _ol(4,23) = _pma[1]*_ol(1,23) + _fak * _ol(1,10);
     _ol(5,23) = _pma[2]*_ol(1,23);
     _ol(7,24) = _pma[0]*_ol(1,24) + _fak * (_ol(0,24) + _ol(1,11) );
     _ol(4,24) = _pma[1]*_ol(1,24) + _fak3* _ol(1,14);
     _ol(5,24) = _pma[2]*_ol(1,24);
     _ol(7,25) = _pma[0]*_ol(1,25) + _fak * (_ol(0,25) + 3.0*_ol(1,15) );
     _ol(4,25) = _pma[1]*_ol(1,25);
     _ol(5,25) = _pma[2]*_ol(1,25) + _fak * _ol(1,10);
     _ol(7,26) = _pma[0]*_ol(1,26) + _fak * (_ol(0,26) + _ol(1,12) );
     _ol(4,26) = _pma[1]*_ol(1,26);
     _ol(5,26) = _pma[2]*_ol(1,26) + _fak3* _ol(1,16);
     _ol(7,27) = _pma[0]*_ol(1,27) + _fak * _ol(0,27);
     _ol(4,27) = _pma[1]*_ol(1,27) + _fak3* _ol(1,17);
     _ol(5,27) = _pma[2]*_ol(1,27) + _fak * _ol(1,11);
     _ol(7,28) = _pma[0]*_ol(1,28) + _fak * _ol(0,28);
     _ol(4,28) = _pma[1]*_ol(1,28) + _fak * _ol(1,12);
     _ol(5,28) = _pma[2]*_ol(1,28) + _fak3* _ol(1,18);
     _ol(7,29) = _pma[0]*_ol(1,29) + _fak * (_ol(0,29) + 2.0*_ol(1,16) );
     _ol(4,29) = _pma[1]*_ol(1,29);
     _ol(5,29) = _pma[2]*_ol(1,29) + _fak2* _ol(1,15);
     _ol(7,30) = _pma[0]*_ol(1,30) + _fak * _ol(0,30);
     _ol(4,30) = _pma[1]*_ol(1,30) + _fak2* _ol(1,18);
     _ol(5,30) = _pma[2]*_ol(1,30) + _fak2* _ol(1,17);
     _ol(7,31) = _pma[0]*_ol(1,31) + _fak * (_ol(0,31) + 2.0*_ol(1,14) );
     _ol(4,31) = _pma[1]*_ol(1,31) + _fak2* _ol(1,13);
     _ol(5,31) = _pma[2]*_ol(1,31);
     _ol(7,32) = _pma[0]*_ol(1,32) + _fak * (_ol(0,32) + 2.0*_ol(1,19) );
     _ol(4,32) = _pma[1]*_ol(1,32) + _fak * _ol(1,15);
     _ol(5,32) = _pma[2]*_ol(1,32) + _fak * _ol(1,13);
     _ol(7,33) = _pma[0]*_ol(1,33) + _fak * (_ol(0,33) + _ol(1,17) );
     _ol(4,33) = _pma[1]*_ol(1,33) + _fak2* _ol(1,19);
     _ol(5,33) = _pma[2]*_ol(1,33) + _fak * _ol(1,14);
     _ol(7,34) = _pma[0]*_ol(1,34) + _fak * (_ol(0,34) + _ol(1,18) );
     _ol(4,34) = _pma[1]*_ol(1,34) + _fak * _ol(1,16);
     _ol(5,34) = _pma[2]*_ol(1,34) + _fak2* _ol(1,19);
     _ol(8,20) = _pma[1]*_ol(2,20) + _fak * _ol(0,20);
     _ol(6,20) = _pma[2]*_ol(2,20);
     _ol(8,21) = _pma[1]*_ol(2,21) + _fak * (_ol(0,21) + 4.0*_ol(2,11) );
     _ol(6,21) = _pma[2]*_ol(2,21);
     _ol(8,22) = _pma[1]*_ol(2,22) + _fak * _ol(0,22);
     _ol(6,22) = _pma[2]*_ol(2,22) + _fak4* _ol(2,12);
     _ol(8,23) = _pma[1]*_ol(2,23) + _fak * (_ol(0,23) + _ol(2,10) );
     _ol(6,23) = _pma[2]*_ol(2,23);
     _ol(8,24) = _pma[1]*_ol(2,24) + _fak * (_ol(0,24) + 3.0*_ol(2,14) );
     _ol(6,24) = _pma[2]*_ol(2,24);
     _ol(8,25) = _pma[1]*_ol(2,25) + _fak * _ol(0,25);
     _ol(6,25) = _pma[2]*_ol(2,25) + _fak * _ol(2,10);
     _ol(8,26) = _pma[1]*_ol(2,26) + _fak * _ol(0,26);
     _ol(6,26) = _pma[2]*_ol(2,26) + _fak3* _ol(2,16);
     _ol(8,27) = _pma[1]*_ol(2,27) + _fak * (_ol(0,27) + 3.0*_ol(2,17) );
     _ol(6,27) = _pma[2]*_ol(2,27) + _fak * _ol(2,11);
     _ol(8,28) = _pma[1]*_ol(2,28) + _fak * (_ol(0,28) + _ol(2,12) );
     _ol(6,28) = _pma[2]*_ol(2,28) + _fak3* _ol(2,18);
     _ol(8,29) = _pma[1]*_ol(2,29) + _fak * _ol(0,29);
     _ol(6,29) = _pma[2]*_ol(2,29) + _fak2* _ol(2,15);
     _ol(8,30) = _pma[1]*_ol(2,30) + _fak * (_ol(0,30) + 2.0*_ol(2,18) );
     _ol(6,30) = _pma[2]*_ol(2,30) + _fak2* _ol(2,17);
     _ol(8,31) = _pma[1]*_ol(2,31) + _fak * (_ol(0,31) + 2.0*_ol(2,13) );
     _ol(6,31) = _pma[2]*_ol(2,31);
     _ol(8,32) = _pma[1]*_ol(2,32) + _fak * (_ol(0,32) + _ol(2,15) );
     _ol(6,32) = _pma[2]*_ol(2,32) + _fak * _ol(2,13);
     _ol(8,33) = _pma[1]*_ol(2,33) + _fak * (_ol(0,33) + 2.0*_ol(2,19) );
     _ol(6,33) = _pma[2]*_ol(2,33) + _fak * _ol(2,14);
     _ol(8,34) = _pma[1]*_ol(2,34) + _fak * (_ol(0,34) + _ol(2,16) );
     _ol(6,34) = _pma[2]*_ol(2,34) + _fak2* _ol(2,19);
     _ol(9,20) = _pma[2]*_ol(3,20) + _fak * _ol(0,20);
     _ol(9,21) = _pma[2]*_ol(3,21) + _fak * _ol(0,21);
     _ol(9,22) = _pma[2]*_ol(3,22) + _fak * (_ol(0,22) + 4.0*_ol(3,12) );
     _ol(9,23) = _pma[2]*_ol(3,23) + _fak * _ol(0,23);
     _ol(9,24) = _pma[2]*_ol(3,24) + _fak * _ol(0,24);
     _ol(9,25) = _pma[2]*_ol(3,25) + _fak * (_ol(0,25) + _ol(3,10) );
     _ol(9,26) = _pma[2]*_ol(3,26) + _fak * (_ol(0,26) + 3.0*_ol(3,16) );
     _ol(9,27) = _pma[2]*_ol(3,27) + _fak * (_ol(0,27) + _ol(3,11) );
     _ol(9,28) = _pma[2]*_ol(3,28) + _fak * (_ol(0,28) + 3.0*_ol(3,18) );
     _ol(9,29) = _pma[2]*_ol(3,29) + _fak * (_ol(0,29) + 2.0*_ol(3,15) );
     _ol(9,30) = _pma[2]*_ol(3,30) + _fak * (_ol(0,30) + 2.0*_ol(3,17) );
     _ol(9,31) = _pma[2]*_ol(3,31) + _fak * _ol(0,31);
     _ol(9,32) = _pma[2]*_ol(3,32) + _fak * (_ol(0,32) + _ol(3,13) );
     _ol(9,33) = _pma[2]*_ol(3,33) + _fak * (_ol(0,33) + _ol(3,14) );
     _ol(9,34) = _pma[2]*_ol(3,34) + _fak * (_ol(0,34) + 2.0*_ol(3,19) );
  }
        // g-d
        
         if((_lmax_row)>3 && (_lmax_col)>1) {
     _ol(20,4) = _pma[0]*_ol(10,4) + _fak * (3.0*_ol(7,4) + _ol(10,2) );
     _ol(23,4) = _pma[1]*_ol(10,4) + _fak * _ol(10,1);
     _ol(25,4) = _pma[2]*_ol(10,4);
     _ol(20,5) = _pma[0]*_ol(10,5) + _fak * (3.0*_ol(7,5) + _ol(10,3) );
     _ol(23,5) = _pma[1]*_ol(10,5);
     _ol(25,5) = _pma[2]*_ol(10,5) + _fak * _ol(10,1);
     _ol(20,6) = _pma[0]*_ol(10,6) + _fak3* _ol(7,6);
     _ol(23,6) = _pma[1]*_ol(10,6) + _fak * _ol(10,3);
     _ol(25,6) = _pma[2]*_ol(10,6) + _fak * _ol(10,2);
     _ol(20,7) = _pma[0]*_ol(10,7) + _fak * (3.0*_ol(7,7) + 2.0*_ol(10,1));
     _ol(23,7) = _pma[1]*_ol(10,7);
     _ol(25,7) = _pma[2]*_ol(10,7);
     _ol(20,8) = _pma[0]*_ol(10,8) + _fak3* _ol(7,8);
     _ol(23,8) = _pma[1]*_ol(10,8) + _fak2* _ol(10,2);
     _ol(25,8) = _pma[2]*_ol(10,8);
     _ol(20,9) = _pma[0]*_ol(10,9) + _fak3* _ol(7,9);
     _ol(23,9) = _pma[1]*_ol(10,9);
     _ol(25,9) = _pma[2]*_ol(10,9) + _fak2* _ol(10,3);
     _ol(24,4) = _pma[0]*_ol(11,4) + _fak * _ol(11,2);
     _ol(21,4) = _pma[1]*_ol(11,4) + _fak * (3.0*_ol(8,4) + _ol(11,1) );
     _ol(27,4) = _pma[2]*_ol(11,4);
     _ol(24,5) = _pma[0]*_ol(11,5) + _fak * _ol(11,3);
     _ol(21,5) = _pma[1]*_ol(11,5) + _fak3* _ol(8,5);
     _ol(27,5) = _pma[2]*_ol(11,5) + _fak * _ol(11,1);
     _ol(24,6) = _pma[0]*_ol(11,6);
     _ol(21,6) = _pma[1]*_ol(11,6) + _fak * (3.0*_ol(8,6) + _ol(11,3) );
     _ol(27,6) = _pma[2]*_ol(11,6) + _fak * _ol(11,2);
     _ol(24,7) = _pma[0]*_ol(11,7) + _fak2* _ol(11,1);
     _ol(21,7) = _pma[1]*_ol(11,7) + _fak3* _ol(8,7);
     _ol(27,7) = _pma[2]*_ol(11,7);
     _ol(24,8) = _pma[0]*_ol(11,8);
     _ol(21,8) = _pma[1]*_ol(11,8) + _fak * (3.0*_ol(8,8) + 2.0*_ol(11,2));
     _ol(27,8) = _pma[2]*_ol(11,8);
     _ol(24,9) = _pma[0]*_ol(11,9);
     _ol(21,9) = _pma[1]*_ol(11,9) + _fak3* _ol(8,9);
     _ol(27,9) = _pma[2]*_ol(11,9) + _fak2* _ol(11,3);
     _ol(26,4) = _pma[0]*_ol(12,4) + _fak * _ol(12,2);
     _ol(28,4) = _pma[1]*_ol(12,4) + _fak * _ol(12,1);
     _ol(22,4) = _pma[2]*_ol(12,4) + _fak3* _ol(9,4);
     _ol(26,5) = _pma[0]*_ol(12,5) + _fak * _ol(12,3);
     _ol(28,5) = _pma[1]*_ol(12,5);
     _ol(22,5) = _pma[2]*_ol(12,5) + _fak * (3.0*_ol(9,5) + _ol(12,1) );
     _ol(26,6) = _pma[0]*_ol(12,6);
     _ol(28,6) = _pma[1]*_ol(12,6) + _fak * _ol(12,3);
     _ol(22,6) = _pma[2]*_ol(12,6) + _fak * (3.0*_ol(9,6) + _ol(12,2) );
     _ol(26,7) = _pma[0]*_ol(12,7) + _fak2* _ol(12,1);
     _ol(28,7) = _pma[1]*_ol(12,7);
     _ol(22,7) = _pma[2]*_ol(12,7) + _fak3* _ol(9,7);
     _ol(26,8) = _pma[0]*_ol(12,8);
     _ol(28,8) = _pma[1]*_ol(12,8) + _fak2* _ol(12,2);
     _ol(22,8) = _pma[2]*_ol(12,8) + _fak3* _ol(9,8);
     _ol(26,9) = _pma[0]*_ol(12,9);
     _ol(28,9) = _pma[1]*_ol(12,9);
     _ol(22,9) = _pma[2]*_ol(12,9) + _fak * (3.0*_ol(9,9) + 2.0*_ol(12,3));
     _ol(31,4) = _pma[1]*_ol(13,4) + _fak * (_ol(7,4) + _ol(13,1) );
     _ol(32,4) = _pma[2]*_ol(13,4);
     _ol(31,5) = _pma[1]*_ol(13,5) + _fak * _ol(7,5);
     _ol(32,5) = _pma[2]*_ol(13,5) + _fak * _ol(13,1);
     _ol(31,6) = _pma[1]*_ol(13,6) + _fak * (_ol(7,6) + _ol(13,3) );
     _ol(32,6) = _pma[2]*_ol(13,6) + _fak * _ol(13,2);
     _ol(31,7) = _pma[1]*_ol(13,7) + _fak * _ol(7,7);
     _ol(32,7) = _pma[2]*_ol(13,7);
     _ol(31,8) = _pma[1]*_ol(13,8) + _fak * (_ol(7,8) + 2.0*_ol(13,2) );
     _ol(32,8) = _pma[2]*_ol(13,8);
     _ol(31,9) = _pma[1]*_ol(13,9) + _fak * _ol(7,9);
     _ol(32,9) = _pma[2]*_ol(13,9) + _fak2* _ol(13,3);
     _ol(33,4) = _pma[2]*_ol(14,4);
     _ol(33,5) = _pma[2]*_ol(14,5) + _fak * _ol(14,1);
     _ol(33,6) = _pma[2]*_ol(14,6) + _fak * _ol(14,2);
     _ol(33,7) = _pma[2]*_ol(14,7);
     _ol(33,8) = _pma[2]*_ol(14,8);
     _ol(33,9) = _pma[2]*_ol(14,9) + _fak2* _ol(14,3);
     _ol(29,4) = _pma[2]*_ol(15,4) + _fak * _ol(7,4);
     _ol(29,5) = _pma[2]*_ol(15,5) + _fak * (_ol(7,5) + _ol(15,1) );
     _ol(29,6) = _pma[2]*_ol(15,6) + _fak * (_ol(7,6) + _ol(15,2) );
     _ol(29,7) = _pma[2]*_ol(15,7) + _fak * _ol(7,7);
     _ol(29,8) = _pma[2]*_ol(15,8) + _fak * _ol(7,8);
     _ol(29,9) = _pma[2]*_ol(15,9) + _fak * (_ol(7,9) + 2.0*_ol(15,3) );
     _ol(34,4) = _pma[1]*_ol(16,4) + _fak * _ol(16,1);
     _ol(34,5) = _pma[1]*_ol(16,5);
     _ol(34,6) = _pma[1]*_ol(16,6) + _fak * _ol(16,3);
     _ol(34,7) = _pma[1]*_ol(16,7);
     _ol(34,8) = _pma[1]*_ol(16,8) + _fak2* _ol(16,2);
     _ol(34,9) = _pma[1]*_ol(16,9);
     _ol(30,4) = _pma[2]*_ol(17,4) + _fak * _ol(8,4);
     _ol(30,5) = _pma[2]*_ol(17,5) + _fak * (_ol(8,5) + _ol(17,1) );
     _ol(30,6) = _pma[2]*_ol(17,6) + _fak * (_ol(8,6) + _ol(17,2) );
     _ol(30,7) = _pma[2]*_ol(17,7) + _fak * _ol(8,7);
     _ol(30,8) = _pma[2]*_ol(17,8) + _fak * _ol(8,8);
     _ol(30,9) = _pma[2]*_ol(17,9) + _fak * (_ol(8,9) + 2.0*_ol(17,3) );
  }
        // f-g
        
        if((_lmax_row)>2 && (_lmax_col)>3) {
     _ol(13,20) = _pma[0]*_ol(4,20) + _fak * (_ol(2,20) + 4.0*_ol(4,10) );
     _ol(14,20) = _pma[1]*_ol(4,20) + _fak * _ol(1,20);
     _ol(19,20) = _pma[2]*_ol(4,20);
     _ol(13,21) = _pma[0]*_ol(4,21) + _fak * _ol(2,21);
     _ol(14,21) = _pma[1]*_ol(4,21) + _fak * (_ol(1,21) + 4.0*_ol(4,11) );
     _ol(19,21) = _pma[2]*_ol(4,21);
     _ol(13,22) = _pma[0]*_ol(4,22) + _fak * _ol(2,22);
     _ol(14,22) = _pma[1]*_ol(4,22) + _fak * _ol(1,22);
     _ol(19,22) = _pma[2]*_ol(4,22) + _fak4* _ol(4,12);
     _ol(13,23) = _pma[0]*_ol(4,23) + _fak * (_ol(2,23) + 3.0*_ol(4,13) );
     _ol(14,23) = _pma[1]*_ol(4,23) + _fak * (_ol(1,23) + _ol(4,10) );
     _ol(19,23) = _pma[2]*_ol(4,23);
     _ol(13,24) = _pma[0]*_ol(4,24) + _fak * (_ol(2,24) + _ol(4,11) );
     _ol(14,24) = _pma[1]*_ol(4,24) + _fak * (_ol(1,24) + 3.0*_ol(4,14) );
     _ol(19,24) = _pma[2]*_ol(4,24);
     _ol(13,25) = _pma[0]*_ol(4,25) + _fak * (_ol(2,25) + 3.0*_ol(4,15) );
     _ol(14,25) = _pma[1]*_ol(4,25) + _fak * _ol(1,25);
     _ol(19,25) = _pma[2]*_ol(4,25) + _fak * _ol(4,10);
     _ol(13,26) = _pma[0]*_ol(4,26) + _fak * (_ol(2,26) + _ol(4,12) );
     _ol(14,26) = _pma[1]*_ol(4,26) + _fak * _ol(1,26);
     _ol(19,26) = _pma[2]*_ol(4,26) + _fak3* _ol(4,16);
     _ol(13,27) = _pma[0]*_ol(4,27) + _fak * _ol(2,27);
     _ol(14,27) = _pma[1]*_ol(4,27) + _fak * (_ol(1,27) + 3.0*_ol(4,17) );
     _ol(19,27) = _pma[2]*_ol(4,27) + _fak * _ol(4,11);
     _ol(13,28) = _pma[0]*_ol(4,28) + _fak * _ol(2,28);
     _ol(14,28) = _pma[1]*_ol(4,28) + _fak * (_ol(1,28) + _ol(4,12) );
     _ol(19,28) = _pma[2]*_ol(4,28) + _fak3* _ol(4,18);
     _ol(13,29) = _pma[0]*_ol(4,29) + _fak * (_ol(2,29) + 2.0*_ol(4,16) );
     _ol(14,29) = _pma[1]*_ol(4,29) + _fak * _ol(1,29);
     _ol(19,29) = _pma[2]*_ol(4,29) + _fak2* _ol(4,15);
     _ol(13,30) = _pma[0]*_ol(4,30) + _fak * _ol(2,30);
     _ol(14,30) = _pma[1]*_ol(4,30) + _fak * (_ol(1,30) + 2.0*_ol(4,18) );
     _ol(19,30) = _pma[2]*_ol(4,30) + _fak2* _ol(4,17);
     _ol(13,31) = _pma[0]*_ol(4,31) + _fak * (_ol(2,31) + 2.0*_ol(4,14) );
     _ol(14,31) = _pma[1]*_ol(4,31) + _fak * (_ol(1,31) + 2.0*_ol(4,13) );
     _ol(19,31) = _pma[2]*_ol(4,31);
     _ol(13,32) = _pma[0]*_ol(4,32) + _fak * (_ol(2,32) + 2.0*_ol(4,19) );
     _ol(14,32) = _pma[1]*_ol(4,32) + _fak * (_ol(1,32) + _ol(4,15) );
     _ol(19,32) = _pma[2]*_ol(4,32) + _fak * _ol(4,13);
     _ol(13,33) = _pma[0]*_ol(4,33) + _fak * (_ol(2,33) + _ol(4,17) );
     _ol(14,33) = _pma[1]*_ol(4,33) + _fak * (_ol(1,33) + 2.0*_ol(4,19) );
     _ol(19,33) = _pma[2]*_ol(4,33) + _fak * _ol(4,14);
     _ol(13,34) = _pma[0]*_ol(4,34) + _fak * (_ol(2,34) + _ol(4,18) );
     _ol(14,34) = _pma[1]*_ol(4,34) + _fak * (_ol(1,34) + _ol(4,16) );
     _ol(19,34) = _pma[2]*_ol(4,34) + _fak2* _ol(4,19);
     _ol(15,20) = _pma[0]*_ol(5,20) + _fak * (_ol(3,20) + 4.0*_ol(5,10) );
     _ol(16,20) = _pma[2]*_ol(5,20) + _fak * _ol(1,20);
     _ol(15,21) = _pma[0]*_ol(5,21) + _fak * _ol(3,21);
     _ol(16,21) = _pma[2]*_ol(5,21) + _fak * _ol(1,21);
     _ol(15,22) = _pma[0]*_ol(5,22) + _fak * _ol(3,22);
     _ol(16,22) = _pma[2]*_ol(5,22) + _fak * (_ol(1,22) + 4.0*_ol(5,12) );
     _ol(15,23) = _pma[0]*_ol(5,23) + _fak * (_ol(3,23) + 3.0*_ol(5,13) );
     _ol(16,23) = _pma[2]*_ol(5,23) + _fak * _ol(1,23);
     _ol(15,24) = _pma[0]*_ol(5,24) + _fak * (_ol(3,24) + _ol(5,11) );
     _ol(16,24) = _pma[2]*_ol(5,24) + _fak * _ol(1,24);
     _ol(15,25) = _pma[0]*_ol(5,25) + _fak * (_ol(3,25) + 3.0*_ol(5,15) );
     _ol(16,25) = _pma[2]*_ol(5,25) + _fak * (_ol(1,25) + _ol(5,10) );
     _ol(15,26) = _pma[0]*_ol(5,26) + _fak * (_ol(3,26) + _ol(5,12) );
     _ol(16,26) = _pma[2]*_ol(5,26) + _fak * (_ol(1,26) + 3.0*_ol(5,16) );
     _ol(15,27) = _pma[0]*_ol(5,27) + _fak * _ol(3,27);
     _ol(16,27) = _pma[2]*_ol(5,27) + _fak * (_ol(1,27) + _ol(5,11) );
     _ol(15,28) = _pma[0]*_ol(5,28) + _fak * _ol(3,28);
     _ol(16,28) = _pma[2]*_ol(5,28) + _fak * (_ol(1,28) + 3.0*_ol(5,18) );
     _ol(15,29) = _pma[0]*_ol(5,29) + _fak * (_ol(3,29) + 2.0*_ol(5,16) );
     _ol(16,29) = _pma[2]*_ol(5,29) + _fak * (_ol(1,29) + 2.0*_ol(5,15) );
     _ol(15,30) = _pma[0]*_ol(5,30) + _fak * _ol(3,30);
     _ol(16,30) = _pma[2]*_ol(5,30) + _fak * (_ol(1,30) + 2.0*_ol(5,17) );
     _ol(15,31) = _pma[0]*_ol(5,31) + _fak * (_ol(3,31) + 2.0*_ol(5,14) );
     _ol(16,31) = _pma[2]*_ol(5,31) + _fak * _ol(1,31);
     _ol(15,32) = _pma[0]*_ol(5,32) + _fak * (_ol(3,32) + 2.0*_ol(5,19) );
     _ol(16,32) = _pma[2]*_ol(5,32) + _fak * (_ol(1,32) + _ol(5,13) );
     _ol(15,33) = _pma[0]*_ol(5,33) + _fak * (_ol(3,33) + _ol(5,17) );
     _ol(16,33) = _pma[2]*_ol(5,33) + _fak * (_ol(1,33) + _ol(5,14) );
     _ol(15,34) = _pma[0]*_ol(5,34) + _fak * (_ol(3,34) + _ol(5,18) );
     _ol(16,34) = _pma[2]*_ol(5,34) + _fak * (_ol(1,34) + 2.0*_ol(5,19) );
     _ol(17,20) = _pma[1]*_ol(6,20) + _fak * _ol(3,20);
     _ol(18,20) = _pma[2]*_ol(6,20) + _fak * _ol(2,20);
     _ol(17,21) = _pma[1]*_ol(6,21) + _fak * (_ol(3,21) + 4.0*_ol(6,11) );
     _ol(18,21) = _pma[2]*_ol(6,21) + _fak * _ol(2,21);
     _ol(17,22) = _pma[1]*_ol(6,22) + _fak * _ol(3,22);
     _ol(18,22) = _pma[2]*_ol(6,22) + _fak * (_ol(2,22) + 4.0*_ol(6,12) );
     _ol(17,23) = _pma[1]*_ol(6,23) + _fak * (_ol(3,23) + _ol(6,10) );
     _ol(18,23) = _pma[2]*_ol(6,23) + _fak * _ol(2,23);
     _ol(17,24) = _pma[1]*_ol(6,24) + _fak * (_ol(3,24) + 3.0*_ol(6,14) );
     _ol(18,24) = _pma[2]*_ol(6,24) + _fak * _ol(2,24);
     _ol(17,25) = _pma[1]*_ol(6,25) + _fak * _ol(3,25);
     _ol(18,25) = _pma[2]*_ol(6,25) + _fak * (_ol(2,25) + _ol(6,10) );
     _ol(17,26) = _pma[1]*_ol(6,26) + _fak * _ol(3,26);
     _ol(18,26) = _pma[2]*_ol(6,26) + _fak * (_ol(2,26) + 3.0*_ol(6,16) );
     _ol(17,27) = _pma[1]*_ol(6,27) + _fak * (_ol(3,27) + 3.0*_ol(6,17) );
     _ol(18,27) = _pma[2]*_ol(6,27) + _fak * (_ol(2,27) + _ol(6,11) );
     _ol(17,28) = _pma[1]*_ol(6,28) + _fak * (_ol(3,28) + _ol(6,12) );
     _ol(18,28) = _pma[2]*_ol(6,28) + _fak * (_ol(2,28) + 3.0*_ol(6,18) );
     _ol(17,29) = _pma[1]*_ol(6,29) + _fak * _ol(3,29);
     _ol(18,29) = _pma[2]*_ol(6,29) + _fak * (_ol(2,29) + 2.0*_ol(6,15) );
     _ol(17,30) = _pma[1]*_ol(6,30) + _fak * (_ol(3,30) + 2.0*_ol(6,18) );
     _ol(18,30) = _pma[2]*_ol(6,30) + _fak * (_ol(2,30) + 2.0*_ol(6,17) );
     _ol(17,31) = _pma[1]*_ol(6,31) + _fak * (_ol(3,31) + 2.0*_ol(6,13) );
     _ol(18,31) = _pma[2]*_ol(6,31) + _fak * _ol(2,31);
     _ol(17,32) = _pma[1]*_ol(6,32) + _fak * (_ol(3,32) + _ol(6,15) );
     _ol(18,32) = _pma[2]*_ol(6,32) + _fak * (_ol(2,32) + _ol(6,13) );
     _ol(17,33) = _pma[1]*_ol(6,33) + _fak * (_ol(3,33) + 2.0*_ol(6,19) );
     _ol(18,33) = _pma[2]*_ol(6,33) + _fak * (_ol(2,33) + _ol(6,14) );
     _ol(17,34) = _pma[1]*_ol(6,34) + _fak * (_ol(3,34) + _ol(6,16) );
     _ol(18,34) = _pma[2]*_ol(6,34) + _fak * (_ol(2,34) + 2.0*_ol(6,19) );
     _ol(10,20) = _pma[0]*_ol(7,20) + _fak * (2.0*_ol(1,20) + 4.0*_ol(7,10));
     _ol(10,21) = _pma[0]*_ol(7,21) + _fak2* _ol(1,21);
     _ol(10,22) = _pma[0]*_ol(7,22) + _fak2* _ol(1,22);
     _ol(10,23) = _pma[0]*_ol(7,23) + _fak * (2.0*_ol(1,23) + 3.0*_ol(7,13));
     _ol(10,24) = _pma[0]*_ol(7,24) + _fak * (2.0*_ol(1,24) + _ol(7,11) );
     _ol(10,25) = _pma[0]*_ol(7,25) + _fak * (2.0*_ol(1,25) + 3.0*_ol(7,15));
     _ol(10,26) = _pma[0]*_ol(7,26) + _fak * (2.0*_ol(1,26) + _ol(7,12) );
     _ol(10,27) = _pma[0]*_ol(7,27) + _fak2* _ol(1,27);
     _ol(10,28) = _pma[0]*_ol(7,28) + _fak2* _ol(1,28);
     _ol(10,29) = _pma[0]*_ol(7,29) + _fak * (2.0*_ol(1,29) + 2.0*_ol(7,16));
     _ol(10,30) = _pma[0]*_ol(7,30) + _fak2* _ol(1,30);
     _ol(10,31) = _pma[0]*_ol(7,31) + _fak * (2.0*_ol(1,31) + 2.0*_ol(7,14));
     _ol(10,32) = _pma[0]*_ol(7,32) + _fak * (2.0*_ol(1,32) + 2.0*_ol(7,19));
     _ol(10,33) = _pma[0]*_ol(7,33) + _fak * (2.0*_ol(1,33) + _ol(7,17) );
     _ol(10,34) = _pma[0]*_ol(7,34) + _fak * (2.0*_ol(1,34) + _ol(7,18) );
     _ol(11,20) = _pma[1]*_ol(8,20) + _fak2* _ol(2,20);
     _ol(11,21) = _pma[1]*_ol(8,21) + _fak * (2.0*_ol(2,21) + 4.0*_ol(8,11));
     _ol(11,22) = _pma[1]*_ol(8,22) + _fak2* _ol(2,22);
     _ol(11,23) = _pma[1]*_ol(8,23) + _fak * (2.0*_ol(2,23) + _ol(8,10) );
     _ol(11,24) = _pma[1]*_ol(8,24) + _fak * (2.0*_ol(2,24) + 3.0*_ol(8,14));
     _ol(11,25) = _pma[1]*_ol(8,25) + _fak2* _ol(2,25);
     _ol(11,26) = _pma[1]*_ol(8,26) + _fak2* _ol(2,26);
     _ol(11,27) = _pma[1]*_ol(8,27) + _fak * (2.0*_ol(2,27) + 3.0*_ol(8,17));
     _ol(11,28) = _pma[1]*_ol(8,28) + _fak * (2.0*_ol(2,28) + _ol(8,12) );
     _ol(11,29) = _pma[1]*_ol(8,29) + _fak2* _ol(2,29);
     _ol(11,30) = _pma[1]*_ol(8,30) + _fak * (2.0*_ol(2,30) + 2.0*_ol(8,18));
     _ol(11,31) = _pma[1]*_ol(8,31) + _fak * (2.0*_ol(2,31) + 2.0*_ol(8,13));
     _ol(11,32) = _pma[1]*_ol(8,32) + _fak * (2.0*_ol(2,32) + _ol(8,15) );
     _ol(11,33) = _pma[1]*_ol(8,33) + _fak * (2.0*_ol(2,33) + 2.0*_ol(8,19));
     _ol(11,34) = _pma[1]*_ol(8,34) + _fak * (2.0*_ol(2,34) + _ol(8,16) );
     _ol(12,20) = _pma[2]*_ol(9,20) + _fak2* _ol(3,20);
     _ol(12,21) = _pma[2]*_ol(9,21) + _fak2* _ol(3,21);
     _ol(12,22) = _pma[2]*_ol(9,22) + _fak * (2.0*_ol(3,22) + 4.0*_ol(9,12));
     _ol(12,23) = _pma[2]*_ol(9,23) + _fak2* _ol(3,23);
     _ol(12,24) = _pma[2]*_ol(9,24) + _fak2* _ol(3,24);
     _ol(12,25) = _pma[2]*_ol(9,25) + _fak * (2.0*_ol(3,25) + _ol(9,10) );
     _ol(12,26) = _pma[2]*_ol(9,26) + _fak * (2.0*_ol(3,26) + 3.0*_ol(9,16));
     _ol(12,27) = _pma[2]*_ol(9,27) + _fak * (2.0*_ol(3,27) + _ol(9,11) );
     _ol(12,28) = _pma[2]*_ol(9,28) + _fak * (2.0*_ol(3,28) + 3.0*_ol(9,18));
     _ol(12,29) = _pma[2]*_ol(9,29) + _fak * (2.0*_ol(3,29) + 2.0*_ol(9,15));
     _ol(12,30) = _pma[2]*_ol(9,30) + _fak * (2.0*_ol(3,30) + 2.0*_ol(9,17));
     _ol(12,31) = _pma[2]*_ol(9,31) + _fak2* _ol(3,31);
     _ol(12,32) = _pma[2]*_ol(9,32) + _fak * (2.0*_ol(3,32) + _ol(9,13) );
     _ol(12,33) = _pma[2]*_ol(9,33) + _fak * (2.0*_ol(3,33) + _ol(9,14) );
     _ol(12,34) = _pma[2]*_ol(9,34) + _fak * (2.0*_ol(3,34) + 2.0*_ol(9,19));
  }
        // g-f
         if((_lmax_row)>3 && (_lmax_col)>2) {
     _ol(20,10) = _pma[0]*_ol(10,10) + _fak * (3.0*_ol(7,10) + 3.0*_ol(10,7));
     _ol(23,10) = _pma[1]*_ol(10,10);
     _ol(25,10) = _pma[2]*_ol(10,10);
     _ol(20,11) = _pma[0]*_ol(10,11) + _fak3* _ol(7,11);
     _ol(23,11) = _pma[1]*_ol(10,11) + _fak3* _ol(10,8);
     _ol(25,11) = _pma[2]*_ol(10,11);
     _ol(20,12) = _pma[0]*_ol(10,12) + _fak3* _ol(7,12);
     _ol(23,12) = _pma[1]*_ol(10,12);
     _ol(25,12) = _pma[2]*_ol(10,12) + _fak3* _ol(10,9);
     _ol(20,13) = _pma[0]*_ol(10,13) + _fak * (3.0*_ol(7,13) + 2.0*_ol(10,4));
     _ol(23,13) = _pma[1]*_ol(10,13) + _fak * _ol(10,7);
     _ol(25,13) = _pma[2]*_ol(10,13);
     _ol(20,14) = _pma[0]*_ol(10,14) + _fak * (3.0*_ol(7,14) + _ol(10,8) );
     _ol(23,14) = _pma[1]*_ol(10,14) + _fak2* _ol(10,4);
     _ol(25,14) = _pma[2]*_ol(10,14);
     _ol(20,15) = _pma[0]*_ol(10,15) + _fak * (3.0*_ol(7,15) + 2.0*_ol(10,5));
     _ol(23,15) = _pma[1]*_ol(10,15);
     _ol(25,15) = _pma[2]*_ol(10,15) + _fak * _ol(10,7);
     _ol(20,16) = _pma[0]*_ol(10,16) + _fak * (3.0*_ol(7,16) + _ol(10,9) );
     _ol(23,16) = _pma[1]*_ol(10,16);
     _ol(25,16) = _pma[2]*_ol(10,16) + _fak2* _ol(10,5);
     _ol(20,17) = _pma[0]*_ol(10,17) + _fak3* _ol(7,17);
     _ol(23,17) = _pma[1]*_ol(10,17) + _fak2* _ol(10,6);
     _ol(25,17) = _pma[2]*_ol(10,17) + _fak * _ol(10,8);
     _ol(20,18) = _pma[0]*_ol(10,18) + _fak3* _ol(7,18);
     _ol(23,18) = _pma[1]*_ol(10,18) + _fak * _ol(10,9);
     _ol(25,18) = _pma[2]*_ol(10,18) + _fak2* _ol(10,6);
     _ol(20,19) = _pma[0]*_ol(10,19) + _fak * (3.0*_ol(7,19) + _ol(10,6) );
     _ol(23,19) = _pma[1]*_ol(10,19) + _fak * _ol(10,5);
     _ol(25,19) = _pma[2]*_ol(10,19) + _fak * _ol(10,4);
     _ol(24,10) = _pma[0]*_ol(11,10) + _fak3* _ol(11,7);
     _ol(21,10) = _pma[1]*_ol(11,10) + _fak3* _ol(8,10);
     _ol(27,10) = _pma[2]*_ol(11,10);
     _ol(24,11) = _pma[0]*_ol(11,11);
     _ol(21,11) = _pma[1]*_ol(11,11) + _fak * (3.0*_ol(8,11) + 3.0*_ol(11,8));
     _ol(27,11) = _pma[2]*_ol(11,11);
     _ol(24,12) = _pma[0]*_ol(11,12);
     _ol(21,12) = _pma[1]*_ol(11,12) + _fak3* _ol(8,12);
     _ol(27,12) = _pma[2]*_ol(11,12) + _fak3* _ol(11,9);
     _ol(24,13) = _pma[0]*_ol(11,13) + _fak2* _ol(11,4);
     _ol(21,13) = _pma[1]*_ol(11,13) + _fak * (3.0*_ol(8,13) + _ol(11,7) );
     _ol(27,13) = _pma[2]*_ol(11,13);
     _ol(24,14) = _pma[0]*_ol(11,14) + _fak * _ol(11,8);
     _ol(21,14) = _pma[1]*_ol(11,14) + _fak * (3.0*_ol(8,14) + 2.0*_ol(11,4));
     _ol(27,14) = _pma[2]*_ol(11,14);
     _ol(24,15) = _pma[0]*_ol(11,15) + _fak2* _ol(11,5);
     _ol(21,15) = _pma[1]*_ol(11,15) + _fak3* _ol(8,15);
     _ol(27,15) = _pma[2]*_ol(11,15) + _fak * _ol(11,7);
     _ol(24,16) = _pma[0]*_ol(11,16) + _fak * _ol(11,9);
     _ol(21,16) = _pma[1]*_ol(11,16) + _fak3* _ol(8,16);
     _ol(27,16) = _pma[2]*_ol(11,16) + _fak2* _ol(11,5);
     _ol(24,17) = _pma[0]*_ol(11,17);
     _ol(21,17) = _pma[1]*_ol(11,17) + _fak * (3.0*_ol(8,17) + 2.0*_ol(11,6));
     _ol(27,17) = _pma[2]*_ol(11,17) + _fak * _ol(11,8);
     _ol(24,18) = _pma[0]*_ol(11,18);
     _ol(21,18) = _pma[1]*_ol(11,18) + _fak * (3.0*_ol(8,18) + _ol(11,9) );
     _ol(27,18) = _pma[2]*_ol(11,18) + _fak2* _ol(11,6);
     _ol(24,19) = _pma[0]*_ol(11,19) + _fak * _ol(11,6);
     _ol(21,19) = _pma[1]*_ol(11,19) + _fak * (3.0*_ol(8,19) + _ol(11,5) );
     _ol(27,19) = _pma[2]*_ol(11,19) + _fak * _ol(11,4);
     _ol(26,10) = _pma[0]*_ol(12,10) + _fak3* _ol(12,7);
     _ol(28,10) = _pma[1]*_ol(12,10);
     _ol(22,10) = _pma[2]*_ol(12,10) + _fak3* _ol(9,10);
     _ol(26,11) = _pma[0]*_ol(12,11);
     _ol(28,11) = _pma[1]*_ol(12,11) + _fak3* _ol(12,8);
     _ol(22,11) = _pma[2]*_ol(12,11) + _fak3* _ol(9,11);
     _ol(26,12) = _pma[0]*_ol(12,12);
     _ol(28,12) = _pma[1]*_ol(12,12);
     _ol(22,12) = _pma[2]*_ol(12,12) + _fak * (3.0*_ol(9,12) + 3.0*_ol(12,9));
     _ol(26,13) = _pma[0]*_ol(12,13) + _fak2* _ol(12,4);
     _ol(28,13) = _pma[1]*_ol(12,13) + _fak * _ol(12,7);
     _ol(22,13) = _pma[2]*_ol(12,13) + _fak3* _ol(9,13);
     _ol(26,14) = _pma[0]*_ol(12,14) + _fak * _ol(12,8);
     _ol(28,14) = _pma[1]*_ol(12,14) + _fak2* _ol(12,4);
     _ol(22,14) = _pma[2]*_ol(12,14) + _fak3* _ol(9,14);
     _ol(26,15) = _pma[0]*_ol(12,15) + _fak2* _ol(12,5);
     _ol(28,15) = _pma[1]*_ol(12,15);
     _ol(22,15) = _pma[2]*_ol(12,15) + _fak * (3.0*_ol(9,15) + _ol(12,7) );
     _ol(26,16) = _pma[0]*_ol(12,16) + _fak * _ol(12,9);
     _ol(28,16) = _pma[1]*_ol(12,16);
     _ol(22,16) = _pma[2]*_ol(12,16) + _fak * (3.0*_ol(9,16) + 2.0*_ol(12,5));
     _ol(26,17) = _pma[0]*_ol(12,17);
     _ol(28,17) = _pma[1]*_ol(12,17) + _fak2* _ol(12,6);
     _ol(22,17) = _pma[2]*_ol(12,17) + _fak * (3.0*_ol(9,17) + _ol(12,8) );
     _ol(26,18) = _pma[0]*_ol(12,18);
     _ol(28,18) = _pma[1]*_ol(12,18) + _fak * _ol(12,9);
     _ol(22,18) = _pma[2]*_ol(12,18) + _fak * (3.0*_ol(9,18) + 2.0*_ol(12,6));
     _ol(26,19) = _pma[0]*_ol(12,19) + _fak * _ol(12,6);
     _ol(28,19) = _pma[1]*_ol(12,19) + _fak * _ol(12,5);
     _ol(22,19) = _pma[2]*_ol(12,19) + _fak * (3.0*_ol(9,19) + _ol(12,4) );
     _ol(31,10) = _pma[1]*_ol(13,10) + _fak * _ol(7,10);
     _ol(32,10) = _pma[2]*_ol(13,10);
     _ol(31,11) = _pma[1]*_ol(13,11) + _fak * (_ol(7,11) + 3.0*_ol(13,8) );
     _ol(32,11) = _pma[2]*_ol(13,11);
     _ol(31,12) = _pma[1]*_ol(13,12) + _fak * _ol(7,12);
     _ol(32,12) = _pma[2]*_ol(13,12) + _fak3* _ol(13,9);
     _ol(31,13) = _pma[1]*_ol(13,13) + _fak * (_ol(7,13) + _ol(13,7) );
     _ol(32,13) = _pma[2]*_ol(13,13);
     _ol(31,14) = _pma[1]*_ol(13,14) + _fak * (_ol(7,14) + 2.0*_ol(13,4) );
     _ol(32,14) = _pma[2]*_ol(13,14);
     _ol(31,15) = _pma[1]*_ol(13,15) + _fak * _ol(7,15);
     _ol(32,15) = _pma[2]*_ol(13,15) + _fak * _ol(13,7);
     _ol(31,16) = _pma[1]*_ol(13,16) + _fak * _ol(7,16);
     _ol(32,16) = _pma[2]*_ol(13,16) + _fak2* _ol(13,5);
     _ol(31,17) = _pma[1]*_ol(13,17) + _fak * (_ol(7,17) + 2.0*_ol(13,6) );
     _ol(32,17) = _pma[2]*_ol(13,17) + _fak * _ol(13,8);
     _ol(31,18) = _pma[1]*_ol(13,18) + _fak * (_ol(7,18) + _ol(13,9) );
     _ol(32,18) = _pma[2]*_ol(13,18) + _fak2* _ol(13,6);
     _ol(31,19) = _pma[1]*_ol(13,19) + _fak * (_ol(7,19) + _ol(13,5) );
     _ol(32,19) = _pma[2]*_ol(13,19) + _fak * _ol(13,4);
     _ol(33,10) = _pma[2]*_ol(14,10);
     _ol(33,11) = _pma[2]*_ol(14,11);
     _ol(33,12) = _pma[2]*_ol(14,12) + _fak3* _ol(14,9);
     _ol(33,13) = _pma[2]*_ol(14,13);
     _ol(33,14) = _pma[2]*_ol(14,14);
     _ol(33,15) = _pma[2]*_ol(14,15) + _fak * _ol(14,7);
     _ol(33,16) = _pma[2]*_ol(14,16) + _fak2* _ol(14,5);
     _ol(33,17) = _pma[2]*_ol(14,17) + _fak * _ol(14,8);
     _ol(33,18) = _pma[2]*_ol(14,18) + _fak2* _ol(14,6);
     _ol(33,19) = _pma[2]*_ol(14,19) + _fak * _ol(14,4);
     _ol(29,10) = _pma[2]*_ol(15,10) + _fak * _ol(7,10);
     _ol(29,11) = _pma[2]*_ol(15,11) + _fak * _ol(7,11);
     _ol(29,12) = _pma[2]*_ol(15,12) + _fak * (_ol(7,12) + 3.0*_ol(15,9) );
     _ol(29,13) = _pma[2]*_ol(15,13) + _fak * _ol(7,13);
     _ol(29,14) = _pma[2]*_ol(15,14) + _fak * _ol(7,14);
     _ol(29,15) = _pma[2]*_ol(15,15) + _fak * (_ol(7,15) + _ol(15,7) );
     _ol(29,16) = _pma[2]*_ol(15,16) + _fak * (_ol(7,16) + 2.0*_ol(15,5) );
     _ol(29,17) = _pma[2]*_ol(15,17) + _fak * (_ol(7,17) + _ol(15,8) );
     _ol(29,18) = _pma[2]*_ol(15,18) + _fak * (_ol(7,18) + 2.0*_ol(15,6) );
     _ol(29,19) = _pma[2]*_ol(15,19) + _fak * (_ol(7,19) + _ol(15,4) );
     _ol(34,10) = _pma[1]*_ol(16,10);
     _ol(34,11) = _pma[1]*_ol(16,11) + _fak3* _ol(16,8);
     _ol(34,12) = _pma[1]*_ol(16,12);
     _ol(34,13) = _pma[1]*_ol(16,13) + _fak * _ol(16,7);
     _ol(34,14) = _pma[1]*_ol(16,14) + _fak2* _ol(16,4);
     _ol(34,15) = _pma[1]*_ol(16,15);
     _ol(34,16) = _pma[1]*_ol(16,16);
     _ol(34,17) = _pma[1]*_ol(16,17) + _fak2* _ol(16,6);
     _ol(34,18) = _pma[1]*_ol(16,18) + _fak * _ol(16,9);
     _ol(34,19) = _pma[1]*_ol(16,19) + _fak * _ol(16,5);
     _ol(30,10) = _pma[2]*_ol(17,10) + _fak * _ol(8,10);
     _ol(30,11) = _pma[2]*_ol(17,11) + _fak * _ol(8,11);
     _ol(30,12) = _pma[2]*_ol(17,12) + _fak * (_ol(8,12) + 3.0*_ol(17,9) );
     _ol(30,13) = _pma[2]*_ol(17,13) + _fak * _ol(8,13);
     _ol(30,14) = _pma[2]*_ol(17,14) + _fak * _ol(8,14);
     _ol(30,15) = _pma[2]*_ol(17,15) + _fak * (_ol(8,15) + _ol(17,7) );
     _ol(30,16) = _pma[2]*_ol(17,16) + _fak * (_ol(8,16) + 2.0*_ol(17,5) );
     _ol(30,17) = _pma[2]*_ol(17,17) + _fak * (_ol(8,17) + _ol(17,8) );
     _ol(30,18) = _pma[2]*_ol(17,18) + _fak * (_ol(8,18) + 2.0*_ol(17,6) );
     _ol(30,19) = _pma[2]*_ol(17,19) + _fak * (_ol(8,19) + _ol(17,4) );
  }
        // g-g
          if((_lmax_row)>3 && (_lmax_col)>3) {
     _ol(20,20) = _pma[0]*_ol(10,20) + _fak * (3.0*_ol(7,20) + 4.0*_ol(10,10));
     _ol(23,20) = _pma[1]*_ol(10,20);
     _ol(25,20) = _pma[2]*_ol(10,20);
     _ol(20,21) = _pma[0]*_ol(10,21) + _fak3* _ol(7,21);
     _ol(23,21) = _pma[1]*_ol(10,21) + _fak4* _ol(10,11);
     _ol(25,21) = _pma[2]*_ol(10,21);
     _ol(20,22) = _pma[0]*_ol(10,22) + _fak3* _ol(7,22);
     _ol(23,22) = _pma[1]*_ol(10,22);
     _ol(25,22) = _pma[2]*_ol(10,22) + _fak4* _ol(10,12);
     _ol(20,23) = _pma[0]*_ol(10,23) + _fak * (3.0*_ol(7,23) + 3.0*_ol(10,13));
     _ol(23,23) = _pma[1]*_ol(10,23) + _fak * _ol(10,10);
     _ol(25,23) = _pma[2]*_ol(10,23);
     _ol(20,24) = _pma[0]*_ol(10,24) + _fak * (3.0*_ol(7,24) + _ol(10,11) );
     _ol(23,24) = _pma[1]*_ol(10,24) + _fak3* _ol(10,14);
     _ol(25,24) = _pma[2]*_ol(10,24);
     _ol(20,25) = _pma[0]*_ol(10,25) + _fak * (3.0*_ol(7,25) + 3.0*_ol(10,15));
     _ol(23,25) = _pma[1]*_ol(10,25);
     _ol(25,25) = _pma[2]*_ol(10,25) + _fak * _ol(10,10);
     _ol(20,26) = _pma[0]*_ol(10,26) + _fak * (3.0*_ol(7,26) + _ol(10,12) );
     _ol(23,26) = _pma[1]*_ol(10,26);
     _ol(25,26) = _pma[2]*_ol(10,26) + _fak3* _ol(10,16);
     _ol(20,27) = _pma[0]*_ol(10,27) + _fak3* _ol(7,27);
     _ol(23,27) = _pma[1]*_ol(10,27) + _fak3* _ol(10,17);
     _ol(25,27) = _pma[2]*_ol(10,27) + _fak * _ol(10,11);
     _ol(20,28) = _pma[0]*_ol(10,28) + _fak3* _ol(7,28);
     _ol(23,28) = _pma[1]*_ol(10,28) + _fak * _ol(10,12);
     _ol(25,28) = _pma[2]*_ol(10,28) + _fak3* _ol(10,18);
     _ol(20,29) = _pma[0]*_ol(10,29) + _fak * (3.0*_ol(7,29) + 2.0*_ol(10,16));
     _ol(23,29) = _pma[1]*_ol(10,29);
     _ol(25,29) = _pma[2]*_ol(10,29) + _fak2* _ol(10,15);
     _ol(20,30) = _pma[0]*_ol(10,30) + _fak3* _ol(7,30);
     _ol(23,30) = _pma[1]*_ol(10,30) + _fak2* _ol(10,18);
     _ol(25,30) = _pma[2]*_ol(10,30) + _fak2* _ol(10,17);
     _ol(20,31) = _pma[0]*_ol(10,31) + _fak * (3.0*_ol(7,31) + 2.0*_ol(10,14));
     _ol(23,31) = _pma[1]*_ol(10,31) + _fak2* _ol(10,13);
     _ol(25,31) = _pma[2]*_ol(10,31);
     _ol(20,32) = _pma[0]*_ol(10,32) + _fak * (3.0*_ol(7,32) + 2.0*_ol(10,19));
     _ol(23,32) = _pma[1]*_ol(10,32) + _fak * _ol(10,15);
     _ol(25,32) = _pma[2]*_ol(10,32) + _fak * _ol(10,13);
     _ol(20,33) = _pma[0]*_ol(10,33) + _fak * (3.0*_ol(7,33) + _ol(10,17) );
     _ol(23,33) = _pma[1]*_ol(10,33) + _fak2* _ol(10,19);
     _ol(25,33) = _pma[2]*_ol(10,33) + _fak * _ol(10,14);
     _ol(20,34) = _pma[0]*_ol(10,34) + _fak * (3.0*_ol(7,34) + _ol(10,18) );
     _ol(23,34) = _pma[1]*_ol(10,34) + _fak * _ol(10,16);
     _ol(25,34) = _pma[2]*_ol(10,34) + _fak2* _ol(10,19);
     _ol(24,20) = _pma[0]*_ol(11,20) + _fak4* _ol(11,10);
     _ol(21,20) = _pma[1]*_ol(11,20) + _fak3* _ol(8,20);
     _ol(27,20) = _pma[2]*_ol(11,20);
     _ol(24,21) = _pma[0]*_ol(11,21);
     _ol(21,21) = _pma[1]*_ol(11,21) + _fak * (3.0*_ol(8,21) + 4.0*_ol(11,11));
     _ol(27,21) = _pma[2]*_ol(11,21);
     _ol(24,22) = _pma[0]*_ol(11,22);
     _ol(21,22) = _pma[1]*_ol(11,22) + _fak3* _ol(8,22);
     _ol(27,22) = _pma[2]*_ol(11,22) + _fak4* _ol(11,12);
     _ol(24,23) = _pma[0]*_ol(11,23) + _fak3* _ol(11,13);
     _ol(21,23) = _pma[1]*_ol(11,23) + _fak * (3.0*_ol(8,23) + _ol(11,10) );
     _ol(27,23) = _pma[2]*_ol(11,23);
     _ol(24,24) = _pma[0]*_ol(11,24) + _fak * _ol(11,11);
     _ol(21,24) = _pma[1]*_ol(11,24) + _fak * (3.0*_ol(8,24) + 3.0*_ol(11,14));
     _ol(27,24) = _pma[2]*_ol(11,24);
     _ol(24,25) = _pma[0]*_ol(11,25) + _fak3* _ol(11,15);
     _ol(21,25) = _pma[1]*_ol(11,25) + _fak3* _ol(8,25);
     _ol(27,25) = _pma[2]*_ol(11,25) + _fak * _ol(11,10);
     _ol(24,26) = _pma[0]*_ol(11,26) + _fak * _ol(11,12);
     _ol(21,26) = _pma[1]*_ol(11,26) + _fak3* _ol(8,26);
     _ol(27,26) = _pma[2]*_ol(11,26) + _fak3* _ol(11,16);
     _ol(24,27) = _pma[0]*_ol(11,27);
     _ol(21,27) = _pma[1]*_ol(11,27) + _fak * (3.0*_ol(8,27) + 3.0*_ol(11,17));
     _ol(27,27) = _pma[2]*_ol(11,27) + _fak * _ol(11,11);
     _ol(24,28) = _pma[0]*_ol(11,28);
     _ol(21,28) = _pma[1]*_ol(11,28) + _fak * (3.0*_ol(8,28) + _ol(11,12) );
     _ol(27,28) = _pma[2]*_ol(11,28) + _fak3* _ol(11,18);
     _ol(24,29) = _pma[0]*_ol(11,29) + _fak2* _ol(11,16);
     _ol(21,29) = _pma[1]*_ol(11,29) + _fak3* _ol(8,29);
     _ol(27,29) = _pma[2]*_ol(11,29) + _fak2* _ol(11,15);
     _ol(24,30) = _pma[0]*_ol(11,30);
     _ol(21,30) = _pma[1]*_ol(11,30) + _fak * (3.0*_ol(8,30) + 2.0*_ol(11,18));
     _ol(27,30) = _pma[2]*_ol(11,30) + _fak2* _ol(11,17);
     _ol(24,31) = _pma[0]*_ol(11,31) + _fak2* _ol(11,14);
     _ol(21,31) = _pma[1]*_ol(11,31) + _fak * (3.0*_ol(8,31) + 2.0*_ol(11,13));
     _ol(27,31) = _pma[2]*_ol(11,31);
     _ol(24,32) = _pma[0]*_ol(11,32) + _fak2* _ol(11,19);
     _ol(21,32) = _pma[1]*_ol(11,32) + _fak * (3.0*_ol(8,32) + _ol(11,15) );
     _ol(27,32) = _pma[2]*_ol(11,32) + _fak * _ol(11,13);
     _ol(24,33) = _pma[0]*_ol(11,33) + _fak * _ol(11,17);
     _ol(21,33) = _pma[1]*_ol(11,33) + _fak * (3.0*_ol(8,33) + 2.0*_ol(11,19));
     _ol(27,33) = _pma[2]*_ol(11,33) + _fak * _ol(11,14);
     _ol(24,34) = _pma[0]*_ol(11,34) + _fak * _ol(11,18);
     _ol(21,34) = _pma[1]*_ol(11,34) + _fak * (3.0*_ol(8,34) + _ol(11,16) );
     _ol(27,34) = _pma[2]*_ol(11,34) + _fak2* _ol(11,19);
     _ol(26,20) = _pma[0]*_ol(12,20) + _fak4* _ol(12,10);
     _ol(28,20) = _pma[1]*_ol(12,20);
     _ol(22,20) = _pma[2]*_ol(12,20) + _fak3* _ol(9,20);
     _ol(26,21) = _pma[0]*_ol(12,21);
     _ol(28,21) = _pma[1]*_ol(12,21) + _fak4* _ol(12,11);
     _ol(22,21) = _pma[2]*_ol(12,21) + _fak3* _ol(9,21);
     _ol(26,22) = _pma[0]*_ol(12,22);
     _ol(28,22) = _pma[1]*_ol(12,22);
     _ol(22,22) = _pma[2]*_ol(12,22) + _fak * (3.0*_ol(9,22) + 4.0*_ol(12,12));
     _ol(26,23) = _pma[0]*_ol(12,23) + _fak3* _ol(12,13);
     _ol(28,23) = _pma[1]*_ol(12,23) + _fak * _ol(12,10);
     _ol(22,23) = _pma[2]*_ol(12,23) + _fak3* _ol(9,23);
     _ol(26,24) = _pma[0]*_ol(12,24) + _fak * _ol(12,11);
     _ol(28,24) = _pma[1]*_ol(12,24) + _fak3* _ol(12,14);
     _ol(22,24) = _pma[2]*_ol(12,24) + _fak3* _ol(9,24);
     _ol(26,25) = _pma[0]*_ol(12,25) + _fak3* _ol(12,15);
     _ol(28,25) = _pma[1]*_ol(12,25);
     _ol(22,25) = _pma[2]*_ol(12,25) + _fak * (3.0*_ol(9,25) + _ol(12,10) );
     _ol(26,26) = _pma[0]*_ol(12,26) + _fak * _ol(12,12);
     _ol(28,26) = _pma[1]*_ol(12,26);
     _ol(22,26) = _pma[2]*_ol(12,26) + _fak * (3.0*_ol(9,26) + 3.0*_ol(12,16));
     _ol(26,27) = _pma[0]*_ol(12,27);
     _ol(28,27) = _pma[1]*_ol(12,27) + _fak3* _ol(12,17);
     _ol(22,27) = _pma[2]*_ol(12,27) + _fak * (3.0*_ol(9,27) + _ol(12,11) );
     _ol(26,28) = _pma[0]*_ol(12,28);
     _ol(28,28) = _pma[1]*_ol(12,28) + _fak * _ol(12,12);
     _ol(22,28) = _pma[2]*_ol(12,28) + _fak * (3.0*_ol(9,28) + 3.0*_ol(12,18));
     _ol(26,29) = _pma[0]*_ol(12,29) + _fak2* _ol(12,16);
     _ol(28,29) = _pma[1]*_ol(12,29);
     _ol(22,29) = _pma[2]*_ol(12,29) + _fak * (3.0*_ol(9,29) + 2.0*_ol(12,15));
     _ol(26,30) = _pma[0]*_ol(12,30);
     _ol(28,30) = _pma[1]*_ol(12,30) + _fak2* _ol(12,18);
     _ol(22,30) = _pma[2]*_ol(12,30) + _fak * (3.0*_ol(9,30) + 2.0*_ol(12,17));
     _ol(26,31) = _pma[0]*_ol(12,31) + _fak2* _ol(12,14);
     _ol(28,31) = _pma[1]*_ol(12,31) + _fak2* _ol(12,13);
     _ol(22,31) = _pma[2]*_ol(12,31) + _fak3* _ol(9,31);
     _ol(26,32) = _pma[0]*_ol(12,32) + _fak2* _ol(12,19);
     _ol(28,32) = _pma[1]*_ol(12,32) + _fak * _ol(12,15);
     _ol(22,32) = _pma[2]*_ol(12,32) + _fak * (3.0*_ol(9,32) + _ol(12,13) );
     _ol(26,33) = _pma[0]*_ol(12,33) + _fak * _ol(12,17);
     _ol(28,33) = _pma[1]*_ol(12,33) + _fak2* _ol(12,19);
     _ol(22,33) = _pma[2]*_ol(12,33) + _fak * (3.0*_ol(9,33) + _ol(12,14) );
     _ol(26,34) = _pma[0]*_ol(12,34) + _fak * _ol(12,18);
     _ol(28,34) = _pma[1]*_ol(12,34) + _fak * _ol(12,16);
     _ol(22,34) = _pma[2]*_ol(12,34) + _fak * (3.0*_ol(9,34) + 2.0*_ol(12,19));
     _ol(31,20) = _pma[1]*_ol(13,20) + _fak * _ol(7,20);
     _ol(32,20) = _pma[2]*_ol(13,20);
     _ol(31,21) = _pma[1]*_ol(13,21) + _fak * (_ol(7,21) + 4.0*_ol(13,11) );
     _ol(32,21) = _pma[2]*_ol(13,21);
     _ol(31,22) = _pma[1]*_ol(13,22) + _fak * _ol(7,22);
     _ol(32,22) = _pma[2]*_ol(13,22) + _fak4* _ol(13,12);
     _ol(31,23) = _pma[1]*_ol(13,23) + _fak * (_ol(7,23) + _ol(13,10) );
     _ol(32,23) = _pma[2]*_ol(13,23);
     _ol(31,24) = _pma[1]*_ol(13,24) + _fak * (_ol(7,24) + 3.0*_ol(13,14) );
     _ol(32,24) = _pma[2]*_ol(13,24);
     _ol(31,25) = _pma[1]*_ol(13,25) + _fak * _ol(7,25);
     _ol(32,25) = _pma[2]*_ol(13,25) + _fak * _ol(13,10);
     _ol(31,26) = _pma[1]*_ol(13,26) + _fak * _ol(7,26);
     _ol(32,26) = _pma[2]*_ol(13,26) + _fak3* _ol(13,16);
     _ol(31,27) = _pma[1]*_ol(13,27) + _fak * (_ol(7,27) + 3.0*_ol(13,17) );
     _ol(32,27) = _pma[2]*_ol(13,27) + _fak * _ol(13,11);
     _ol(31,28) = _pma[1]*_ol(13,28) + _fak * (_ol(7,28) + _ol(13,12) );
     _ol(32,28) = _pma[2]*_ol(13,28) + _fak3* _ol(13,18);
     _ol(31,29) = _pma[1]*_ol(13,29) + _fak * _ol(7,29);
     _ol(32,29) = _pma[2]*_ol(13,29) + _fak2* _ol(13,15);
     _ol(31,30) = _pma[1]*_ol(13,30) + _fak * (_ol(7,30) + 2.0*_ol(13,18) );
     _ol(32,30) = _pma[2]*_ol(13,30) + _fak2* _ol(13,17);
     _ol(31,31) = _pma[1]*_ol(13,31) + _fak * (_ol(7,31) + 2.0*_ol(13,13) );
     _ol(32,31) = _pma[2]*_ol(13,31);
     _ol(31,32) = _pma[1]*_ol(13,32) + _fak * (_ol(7,32) + _ol(13,15) );
     _ol(32,32) = _pma[2]*_ol(13,32) + _fak * _ol(13,13);
     _ol(31,33) = _pma[1]*_ol(13,33) + _fak * (_ol(7,33) + 2.0*_ol(13,19) );
     _ol(32,33) = _pma[2]*_ol(13,33) + _fak * _ol(13,14);
     _ol(31,34) = _pma[1]*_ol(13,34) + _fak * (_ol(7,34) + _ol(13,16) );
     _ol(32,34) = _pma[2]*_ol(13,34) + _fak2* _ol(13,19);
     _ol(33,20) = _pma[2]*_ol(14,20);
     _ol(33,21) = _pma[2]*_ol(14,21);
     _ol(33,22) = _pma[2]*_ol(14,22) + _fak4* _ol(14,12);
     _ol(33,23) = _pma[2]*_ol(14,23);
     _ol(33,24) = _pma[2]*_ol(14,24);
     _ol(33,25) = _pma[2]*_ol(14,25) + _fak * _ol(14,10);
     _ol(33,26) = _pma[2]*_ol(14,26) + _fak3* _ol(14,16);
     _ol(33,27) = _pma[2]*_ol(14,27) + _fak * _ol(14,11);
     _ol(33,28) = _pma[2]*_ol(14,28) + _fak3* _ol(14,18);
     _ol(33,29) = _pma[2]*_ol(14,29) + _fak2* _ol(14,15);
     _ol(33,30) = _pma[2]*_ol(14,30) + _fak2* _ol(14,17);
     _ol(33,31) = _pma[2]*_ol(14,31);
     _ol(33,32) = _pma[2]*_ol(14,32) + _fak * _ol(14,13);
     _ol(33,33) = _pma[2]*_ol(14,33) + _fak * _ol(14,14);
     _ol(33,34) = _pma[2]*_ol(14,34) + _fak2* _ol(14,19);
     _ol(29,20) = _pma[2]*_ol(15,20) + _fak * _ol(7,20);
     _ol(29,21) = _pma[2]*_ol(15,21) + _fak * _ol(7,21);
     _ol(29,22) = _pma[2]*_ol(15,22) + _fak * (_ol(7,22) + 4.0*_ol(15,12) );
     _ol(29,23) = _pma[2]*_ol(15,23) + _fak * _ol(7,23);
     _ol(29,24) = _pma[2]*_ol(15,24) + _fak * _ol(7,24);
     _ol(29,25) = _pma[2]*_ol(15,25) + _fak * (_ol(7,25) + _ol(15,10) );
     _ol(29,26) = _pma[2]*_ol(15,26) + _fak * (_ol(7,26) + 3.0*_ol(15,16) );
     _ol(29,27) = _pma[2]*_ol(15,27) + _fak * (_ol(7,27) + _ol(15,11) );
     _ol(29,28) = _pma[2]*_ol(15,28) + _fak * (_ol(7,28) + 3.0*_ol(15,18) );
     _ol(29,29) = _pma[2]*_ol(15,29) + _fak * (_ol(7,29) + 2.0*_ol(15,15) );
     _ol(29,30) = _pma[2]*_ol(15,30) + _fak * (_ol(7,30) + 2.0*_ol(15,17) );
     _ol(29,31) = _pma[2]*_ol(15,31) + _fak * _ol(7,31);
     _ol(29,32) = _pma[2]*_ol(15,32) + _fak * (_ol(7,32) + _ol(15,13) );
     _ol(29,33) = _pma[2]*_ol(15,33) + _fak * (_ol(7,33) + _ol(15,14) );
     _ol(29,34) = _pma[2]*_ol(15,34) + _fak * (_ol(7,34) + 2.0*_ol(15,19) );
     _ol(34,20) = _pma[1]*_ol(16,20);
     _ol(34,21) = _pma[1]*_ol(16,21) + _fak4* _ol(16,11);
     _ol(34,22) = _pma[1]*_ol(16,22);
     _ol(34,23) = _pma[1]*_ol(16,23) + _fak * _ol(16,10);
     _ol(34,24) = _pma[1]*_ol(16,24) + _fak3* _ol(16,14);
     _ol(34,25) = _pma[1]*_ol(16,25);
     _ol(34,26) = _pma[1]*_ol(16,26);
     _ol(34,27) = _pma[1]*_ol(16,27) + _fak3* _ol(16,17);
     _ol(34,28) = _pma[1]*_ol(16,28) + _fak * _ol(16,12);
     _ol(34,29) = _pma[1]*_ol(16,29);
     _ol(34,30) = _pma[1]*_ol(16,30) + _fak2* _ol(16,18);
     _ol(34,31) = _pma[1]*_ol(16,31) + _fak2* _ol(16,13);
     _ol(34,32) = _pma[1]*_ol(16,32) + _fak * _ol(16,15);
     _ol(34,33) = _pma[1]*_ol(16,33) + _fak2* _ol(16,19);
     _ol(34,34) = _pma[1]*_ol(16,34) + _fak * _ol(16,16);
     _ol(30,20) = _pma[2]*_ol(17,20) + _fak * _ol(8,20);
     _ol(30,21) = _pma[2]*_ol(17,21) + _fak * _ol(8,21);
     _ol(30,22) = _pma[2]*_ol(17,22) + _fak * (_ol(8,22) + 4.0*_ol(17,12) );
     _ol(30,23) = _pma[2]*_ol(17,23) + _fak * _ol(8,23);
     _ol(30,24) = _pma[2]*_ol(17,24) + _fak * _ol(8,24);
     _ol(30,25) = _pma[2]*_ol(17,25) + _fak * (_ol(8,25) + _ol(17,10) );
     _ol(30,26) = _pma[2]*_ol(17,26) + _fak * (_ol(8,26) + 3.0*_ol(17,16) );
     _ol(30,27) = _pma[2]*_ol(17,27) + _fak * (_ol(8,27) + _ol(17,11) );
     _ol(30,28) = _pma[2]*_ol(17,28) + _fak * (_ol(8,28) + 3.0*_ol(17,18) );
     _ol(30,29) = _pma[2]*_ol(17,29) + _fak * (_ol(8,29) + 2.0*_ol(17,15) );
     _ol(30,30) = _pma[2]*_ol(17,30) + _fak * (_ol(8,30) + 2.0*_ol(17,17) );
     _ol(30,31) = _pma[2]*_ol(17,31) + _fak * _ol(8,31);
     _ol(30,32) = _pma[2]*_ol(17,32) + _fak * (_ol(8,32) + _ol(17,13) );
     _ol(30,33) = _pma[2]*_ol(17,33) + _fak * (_ol(8,33) + _ol(17,14) );
     _ol(30,34) = _pma[2]*_ol(17,34) + _fak * (_ol(8,34) + 2.0*_ol(17,19) );
  }
        
        
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
             
        ub::matrix<double> _ol_tmp = ub::prod( _trafo_row, _ol );
        ub::matrix<double> _trafo_col_tposed = ub::trans( _trafo_col );
        ub::matrix<double> _ol_sph = ub::prod( _ol_tmp, _trafo_col_tposed );
        // save to _matrix
        for ( unsigned i = 0; i< _matrix.size1(); i++ ) {
            for (unsigned j = 0; j < _matrix.size2(); j++){
                _matrix(i,j) += _ol_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
            }
        }
        
        
        _ol.clear();
            } // _shell_col Gaussians
        } // _shell_row Gaussians
    }
    
  
        
    
    
}}

