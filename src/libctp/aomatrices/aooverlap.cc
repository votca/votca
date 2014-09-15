/* 
 *            Copyright 2009-2012 The VOTCA Development Team
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

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;

    
    void AOOverlap::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col , bool _raw) {
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
        vector<double> _pma (3,0.0);
        vector<double> _pmb (3,0.0);
          
        double _distsq = (_diff.getX()*_diff.getX()) + (_diff.getY()*_diff.getY()) + (_diff.getZ()*_diff.getZ());    

       // cout << "row shell is " << _shell_row->getSize() << " -fold contracted!" << endl;
        //cout << "col shell is " << _shell_col->getSize() << " -fold contracted!" << endl;
        
        typedef vector< AOGaussianPrimitive* >::iterator GaussianIterator;
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

            _pma[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
            _pma[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
            _pma[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

            _pmb[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
            _pmb[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
            _pmb[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
            

        
        // calculate matrix elements
        _ol(Cartesian::s,Cartesian::s) = pow(4.0*_decay_row*_decay_col,0.75) * pow(_fak2,1.5)*exp(-_fak2 * _decay_row * _decay_col *_distsq); // s-s element
        //cout << "\t setting s-s: " << _ol(Cartesian::s,Cartesian::s) << endl;
        
        // s-p
        if ( _lmax_col > 0 ) {
          // cout << "\t setting s-p" << flush;
           _ol(Cartesian::s,Cartesian::x) = _pmb[0]*_ol(Cartesian::s,Cartesian::s); // s-px
           _ol(Cartesian::s,Cartesian::y) = _pmb[1]*_ol(Cartesian::s,Cartesian::s); // s-py
           _ol(Cartesian::s,Cartesian::z) = _pmb[2]*_ol(Cartesian::s,Cartesian::s); // s-pz
        }
        
        // p-s
        if ( _lmax_row > 0 ) {
           //cout << "\t setting p-s" << flush;

           _ol(Cartesian::p,Cartesian::s) = _pma[0]*_ol(Cartesian::s,Cartesian::s); // px-s
           _ol(Cartesian::y,Cartesian::s) = _pma[1]*_ol(Cartesian::s,Cartesian::s); // py-s
           _ol(Cartesian::z,Cartesian::s) = _pma[2]*_ol(Cartesian::s,Cartesian::s); // pz-s
        }
        
        // p-p
        if ( _lmax_row > 0 && _lmax_col > 0 ) {
           //cout << "\t setting p-p" << endl;            
           _ol(Cartesian::x,Cartesian::x) = _pma[0]*_ol(Cartesian::s,Cartesian::x) + _fak * _ol(Cartesian::s,Cartesian::s); // px-px
           _ol(Cartesian::x,Cartesian::y) = _pma[0]*_ol(Cartesian::s,Cartesian::y); // px-py
           _ol(Cartesian::x,Cartesian::z) = _pma[0]*_ol(Cartesian::s,Cartesian::z); // px-pz
           _ol(Cartesian::y,Cartesian::x) = _pma[1]*_ol(Cartesian::s,Cartesian::x); // py-px
           _ol(Cartesian::y,Cartesian::y) = _pma[1]*_ol(Cartesian::s,Cartesian::y) + _fak * _ol(Cartesian::s,Cartesian::s); // py-py
           _ol(Cartesian::y,Cartesian::z) = _pma[1]*_ol(Cartesian::s,Cartesian::z); // py-pz
           _ol(Cartesian::z,Cartesian::x) = _pma[2]*_ol(Cartesian::s,Cartesian::x); // pz-px
           _ol(Cartesian::z,Cartesian::y) = _pma[2]*_ol(Cartesian::s,Cartesian::y); // pz-py
           _ol(Cartesian::z,Cartesian::z) = _pma[2]*_ol(Cartesian::s,Cartesian::z) + _fak * _ol(Cartesian::s,Cartesian::s); // pz-pz
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
            _ol(Cartesian::s,Cartesian::xy) = _pmb[1]*_ol(Cartesian::s,Cartesian::x); // s-dxy
            _ol(Cartesian::s,Cartesian::xz) = _pmb[2]*_ol(Cartesian::s,Cartesian::x); // s-dxz
            _ol(Cartesian::s,Cartesian::yz) = _pmb[2]*_ol(Cartesian::s,Cartesian::y); // s-dyz
            _ol(Cartesian::s,Cartesian::xx) = _pmb[0]*_ol(Cartesian::s,Cartesian::x) + _fak*_ol(Cartesian::s,Cartesian::s); // s-dxx
            _ol(Cartesian::s,Cartesian::yy) = _pmb[1]*_ol(Cartesian::s,Cartesian::y) + _fak*_ol(Cartesian::s,Cartesian::s); // s-dyy
            _ol(Cartesian::s,Cartesian::zz) = _pmb[2]*_ol(Cartesian::s,Cartesian::z) + _fak*_ol(Cartesian::s,Cartesian::s); // s-dzz
        }
        
        // p-d
        if ( _lmax_row > 0 && _lmax_col > 1){
            //cout << "\t setting p-d" << endl;
            
            _ol(Cartesian::x,Cartesian::xy) = _pma[0]*_ol(Cartesian::s,Cartesian::xy) + _fak * _ol(Cartesian::s,Cartesian::y);
            _ol(Cartesian::x,Cartesian::xz) = _pma[0]*_ol(Cartesian::s,Cartesian::xz) + _fak * _ol(Cartesian::s,Cartesian::z);
            _ol(Cartesian::x,Cartesian::yz) = _pma[0]*_ol(Cartesian::s,Cartesian::yz);
            _ol(Cartesian::x,Cartesian::xx) = _pma[0]*_ol(Cartesian::s,Cartesian::xx) + _fak2 * _ol(Cartesian::s,Cartesian::x);
            _ol(Cartesian::x,Cartesian::yy) = _pma[0]*_ol(Cartesian::s,Cartesian::yy);
            _ol(Cartesian::x,Cartesian::zz) = _pma[0]*_ol(Cartesian::s,Cartesian::zz);
 
            _ol(Cartesian::y,Cartesian::xy) = _pma[1]*_ol(Cartesian::s,Cartesian::xy) + _fak * _ol(Cartesian::s,Cartesian::x);
            _ol(Cartesian::y,Cartesian::xz) = _pma[1]*_ol(Cartesian::s,Cartesian::xz);
            _ol(Cartesian::y,Cartesian::yz) = _pma[1]*_ol(Cartesian::s,Cartesian::yz) + _fak * _ol(Cartesian::s,Cartesian::z);
            _ol(Cartesian::y,Cartesian::xx) = _pma[1]*_ol(Cartesian::s,Cartesian::xx);
            _ol(Cartesian::y,Cartesian::yy) = _pma[1]*_ol(Cartesian::s,Cartesian::yy) + _fak2 * _ol(Cartesian::s,Cartesian::y);
            _ol(Cartesian::y,Cartesian::zz) = _pma[1]*_ol(Cartesian::s,Cartesian::zz);

            _ol(Cartesian::z,Cartesian::xy) = _pma[2]*_ol(Cartesian::s,Cartesian::xy);
            _ol(Cartesian::z,Cartesian::xz) = _pma[2]*_ol(Cartesian::s,Cartesian::xz) + _fak * _ol(Cartesian::s,Cartesian::x);
            _ol(Cartesian::z,Cartesian::yz) = _pma[2]*_ol(Cartesian::s,Cartesian::yz) + _fak * _ol(Cartesian::s,Cartesian::y);
            _ol(Cartesian::z,Cartesian::xx) = _pma[2]*_ol(Cartesian::s,Cartesian::xx);
            _ol(Cartesian::z,Cartesian::yy) = _pma[2]*_ol(Cartesian::s,Cartesian::yy);
            _ol(Cartesian::z,Cartesian::zz) = _pma[2]*_ol(Cartesian::s,Cartesian::zz) + _fak2 * _ol(Cartesian::s,Cartesian::z);
        }

        // d-s
        if ( _lmax_row > 1){
           //cout << "\t setting d-s" << endl;
            _ol(Cartesian::xy,Cartesian::s) = _pma[1]*_ol(Cartesian::p,Cartesian::s); // dxy-s
            _ol(Cartesian::xy,Cartesian::s) = _pma[2]*_ol(Cartesian::p,Cartesian::s); // dxz-s
            _ol(Cartesian::yz,Cartesian::s) = _pma[2]*_ol(Cartesian::y,Cartesian::s); // dyz-s
            _ol(Cartesian::xx,Cartesian::s) = _pma[0]*_ol(Cartesian::p,Cartesian::s) + _fak * _ol(Cartesian::s,Cartesian::s); // dxx-s
            _ol(Cartesian::yy,Cartesian::s) = _pma[1]*_ol(Cartesian::y,Cartesian::s) + _fak * _ol(Cartesian::s,Cartesian::s); // dyy-s
            _ol(Cartesian::zz,Cartesian::s) = _pma[2]*_ol(Cartesian::z,Cartesian::s) + _fak * _ol(Cartesian::s,Cartesian::s); // dzz-s
        }
        
        
        // d-p
        if ( _lmax_row >1 && _lmax_col > 0){
           //cout << "\t setting d-p" << endl;

             _ol(Cartesian::xy,Cartesian::x) = _pma[1]*_ol(Cartesian::x,Cartesian::x);
             _ol(Cartesian::xy,Cartesian::x) = _pma[2]*_ol(Cartesian::x,Cartesian::x);
             _ol(Cartesian::yz,Cartesian::x) = _pma[2]*_ol(Cartesian::y,Cartesian::x);
             
             _ol(Cartesian::xx,Cartesian::x) = _pma[0]*_ol(Cartesian::x,Cartesian::x) + _fak  * ( _ol(Cartesian::s,Cartesian::x) + _ol(Cartesian::p,Cartesian::s) );
             _ol(Cartesian::yy,Cartesian::x) = _pma[1]*_ol(Cartesian::y,Cartesian::x) + _fak  * _ol(Cartesian::s,Cartesian::x);
             _ol(Cartesian::zz,Cartesian::x) = _pma[2]*_ol(Cartesian::z,Cartesian::x) + _fak  * _ol(Cartesian::s,Cartesian::x);

             _ol(Cartesian::xy,Cartesian::y) = _pma[1]*_ol(Cartesian::x,Cartesian::y) + _fak  * _ol(Cartesian::p,Cartesian::s);
             _ol(Cartesian::xy,Cartesian::y) = _pma[2]*_ol(Cartesian::x,Cartesian::y);
             _ol(Cartesian::yz,Cartesian::y) = _pma[2]*_ol(Cartesian::y,Cartesian::y);
             
             _ol(Cartesian::xx,Cartesian::y) = _pma[0]*_ol(Cartesian::x,Cartesian::y) + _fak  * _ol(Cartesian::s,Cartesian::y);
             _ol(Cartesian::yy,Cartesian::y) = _pma[1]*_ol(Cartesian::y,Cartesian::y) + _fak  * ( _ol(Cartesian::s,Cartesian::y) + _ol (2,Cartesian::s) );
             _ol(Cartesian::zz,Cartesian::y) = _pma[2]*_ol(Cartesian::z,Cartesian::y) + _fak  * _ol(Cartesian::s,Cartesian::y);

             _ol(Cartesian::xy,Cartesian::z) = _pma[1]*_ol(Cartesian::x,Cartesian::z);
             _ol(Cartesian::xy,Cartesian::z) = _pma[2]*_ol(Cartesian::x,Cartesian::z) + _fak  * _ol(Cartesian::p,Cartesian::s);
             _ol(Cartesian::yz,Cartesian::z) = _pma[2]*_ol(Cartesian::y,Cartesian::z) + _fak  * _ol(Cartesian::y,Cartesian::s);
             _ol(Cartesian::xx,Cartesian::z) = _pma[0]*_ol(Cartesian::x,Cartesian::z) + _fak  * _ol(Cartesian::s,Cartesian::z);
             _ol(Cartesian::yy,Cartesian::z) = _pma[1]*_ol(Cartesian::y,Cartesian::z) + _fak  * _ol(Cartesian::s,Cartesian::z);
             _ol(Cartesian::zz,Cartesian::z) = _pma[2]*_ol(Cartesian::z,Cartesian::z) + _fak  * ( _ol(Cartesian::s,Cartesian::z) + _ol(Cartesian::z,Cartesian::s) );
           
        }
        
        // d-d
        if ( _lmax_row > 1 && _lmax_col > 1 ){
             // cout << "\t setting d-d" << endl;
            
             _ol(Cartesian::xy,Cartesian::xy) = _pma[1]*_ol(Cartesian::x,Cartesian::xy) + _fak * _ol(Cartesian::x,Cartesian::x);
             _ol(Cartesian::xy,Cartesian::xy) = _pma[2]*_ol(Cartesian::x,Cartesian::xy);
             _ol(Cartesian::yz,Cartesian::xy) = _pma[2]*_ol(Cartesian::y,Cartesian::xy);
             _ol(Cartesian::xx,Cartesian::xy) = _pma[0]*_ol(Cartesian::x,Cartesian::xy) + _fak * (_ol(Cartesian::s,Cartesian::xy) + _ol(Cartesian::x,Cartesian::y) );
             _ol(Cartesian::yy,Cartesian::xy) = _pma[1]*_ol(Cartesian::y,Cartesian::xy) + _fak * (_ol(Cartesian::s,Cartesian::xy) + _ol(Cartesian::y,Cartesian::x) );
             _ol(Cartesian::zz,Cartesian::xy) = _pma[2]*_ol(Cartesian::z,Cartesian::xy) + _fak * _ol(Cartesian::s,Cartesian::xy);

             _ol(Cartesian::xy,Cartesian::xz) = _pma[1]*_ol(Cartesian::x,Cartesian::xz);
             _ol(Cartesian::xy,Cartesian::xz) = _pma[2]*_ol(Cartesian::x,Cartesian::xz) + _fak * _ol(Cartesian::x,Cartesian::x);
             _ol(Cartesian::yz,Cartesian::xz) = _pma[2]*_ol(Cartesian::y,Cartesian::xz) + _fak * _ol(Cartesian::y,Cartesian::x);
             _ol(Cartesian::xx,Cartesian::xz) = _pma[0]*_ol(Cartesian::x,Cartesian::xz) + _fak * (_ol(Cartesian::s,Cartesian::xz) + _ol(Cartesian::x,Cartesian::z) );
             _ol(Cartesian::yy,Cartesian::xz) = _pma[1]*_ol(Cartesian::y,Cartesian::xz) + _fak * _ol(Cartesian::s,Cartesian::xz);
             _ol(Cartesian::zz,Cartesian::xz) = _pma[2]*_ol(Cartesian::z,Cartesian::xz) + _fak * (_ol(Cartesian::s,Cartesian::xz) + _ol(Cartesian::z,Cartesian::x) );

             _ol(Cartesian::xy,Cartesian::yz) = _pma[1]*_ol(Cartesian::x,Cartesian::yz) + _fak * _ol(Cartesian::x,Cartesian::z);
             _ol(Cartesian::xy,Cartesian::yz) = _pma[2]*_ol(Cartesian::x,Cartesian::yz) + _fak * _ol(Cartesian::x,Cartesian::y);
             _ol(Cartesian::yz,Cartesian::yz) = _pma[2]*_ol(Cartesian::y,Cartesian::yz) + _fak * _ol(Cartesian::y,Cartesian::y);
             _ol(Cartesian::xx,Cartesian::yz) = _pma[0]*_ol(Cartesian::x,Cartesian::yz) + _fak * _ol(Cartesian::s,Cartesian::yz);
             _ol(Cartesian::yy,Cartesian::yz) = _pma[1]*_ol(Cartesian::y,Cartesian::yz) + _fak * (_ol(Cartesian::s,Cartesian::yz) + _ol(Cartesian::y,Cartesian::z) );
             _ol(Cartesian::zz,Cartesian::yz) = _pma[2]*_ol(Cartesian::z,Cartesian::yz) + _fak * (_ol(Cartesian::s,Cartesian::yz) + _ol(Cartesian::z,Cartesian::y) );

             _ol(Cartesian::xy,Cartesian::xx) = _pma[1]*_ol(Cartesian::x,Cartesian::xx);
             _ol(Cartesian::xy,Cartesian::xx) = _pma[2]*_ol(Cartesian::x,Cartesian::xx);
             _ol(Cartesian::yz,Cartesian::xx) = _pma[2]*_ol(Cartesian::y,Cartesian::xx);
             _ol(Cartesian::xx,Cartesian::xx) = _pma[0]*_ol(Cartesian::x,Cartesian::xx) + _fak * (_ol(Cartesian::s,Cartesian::xx) + 2.0*_ol(Cartesian::x,Cartesian::x) );
             _ol(Cartesian::yy,Cartesian::xx) = _pma[1]*_ol(Cartesian::y,Cartesian::xx) + _fak * _ol(Cartesian::s,Cartesian::xx);
             _ol(Cartesian::zz,Cartesian::xx) = _pma[2]*_ol(Cartesian::z,Cartesian::xx) + _fak * _ol(Cartesian::s,Cartesian::xx);

             _ol(Cartesian::xy,Cartesian::yy) = _pma[1]*_ol(Cartesian::x,Cartesian::yy) + _fak2 * _ol(Cartesian::x,Cartesian::y);
             _ol(Cartesian::xy,Cartesian::yy) = _pma[2]*_ol(Cartesian::x,Cartesian::yy);
             _ol(Cartesian::yz,Cartesian::yy) = _pma[2]*_ol(Cartesian::y,Cartesian::yy);
             _ol(Cartesian::xx,Cartesian::yy) = _pma[0]*_ol(Cartesian::x,Cartesian::yy) + _fak * _ol(Cartesian::s,Cartesian::yy);
             _ol(Cartesian::yy,Cartesian::yy) = _pma[1]*_ol(Cartesian::y,Cartesian::yy) + _fak * (_ol(Cartesian::s,Cartesian::yy) + 2.0*_ol(Cartesian::y,Cartesian::y) );
             _ol(Cartesian::zz,Cartesian::yy) = _pma[2]*_ol(Cartesian::z,Cartesian::yy) + _fak * _ol(Cartesian::s,Cartesian::yy);

             _ol(Cartesian::xy,Cartesian::zz) = _pma[1]*_ol(Cartesian::x,Cartesian::zz);
             _ol(Cartesian::xy,Cartesian::zz) = _pma[2]*_ol(Cartesian::x,Cartesian::zz) + _fak2 * _ol(Cartesian::x,Cartesian::z);
             _ol(Cartesian::yz,Cartesian::zz) = _pma[2]*_ol(Cartesian::y,Cartesian::zz) + _fak2 * _ol(Cartesian::y,Cartesian::z);
             _ol(Cartesian::xx,Cartesian::zz) = _pma[0]*_ol(Cartesian::x,Cartesian::zz) + _fak * _ol( 0,Cartesian::zz);
             _ol(Cartesian::yy,Cartesian::zz) = _pma[1]*_ol(Cartesian::y,Cartesian::zz) + _fak * _ol(Cartesian::s,Cartesian::zz);
             _ol(Cartesian::zz,Cartesian::zz) = _pma[2]*_ol(Cartesian::z,Cartesian::zz) + _fak * (_ol(Cartesian::s,Cartesian::zz) + 2.0*_ol(Cartesian::z,Cartesian::z) );
            
            
        }

        // s-f 
        if ( _lmax_col > 2 ){
             _ol(Cartesian::s,10) = _pmb[0]*_ol(Cartesian::s,Cartesian::xx) + _fak2* _ol(Cartesian::s,Cartesian::x);
             _ol(Cartesian::s,11) = _pmb[1]*_ol(Cartesian::s,Cartesian::yy) + _fak2* _ol(Cartesian::s,Cartesian::y);
             _ol(Cartesian::s,12) = _pmb[2]*_ol(Cartesian::s,Cartesian::zz) + _fak2* _ol(Cartesian::s,Cartesian::z);
             _ol(Cartesian::s,13) = _pmb[0]*_ol(Cartesian::s,Cartesian::xy) + _fak * _ol(Cartesian::s,Cartesian::y);
             _ol(Cartesian::s,14) = _pmb[1]*_ol(Cartesian::s,Cartesian::xy) + _fak * _ol(Cartesian::s,Cartesian::x);
             _ol(Cartesian::s,15) = _pmb[0]*_ol(Cartesian::s,Cartesian::xz) + _fak * _ol(Cartesian::s,Cartesian::z);
             _ol(Cartesian::s,16) = _pmb[2]*_ol(Cartesian::s,Cartesian::xz) + _fak * _ol(Cartesian::s,Cartesian::x);
             _ol(Cartesian::s,17) = _pmb[1]*_ol(Cartesian::s,Cartesian::yz) + _fak * _ol(Cartesian::s,Cartesian::z);
             _ol(Cartesian::s,18) = _pmb[2]*_ol(Cartesian::s,Cartesian::yz) + _fak * _ol(Cartesian::s,Cartesian::y);
             _ol(Cartesian::s,19) = _pmb[2]*_ol(Cartesian::s,Cartesian::xy);
        }

        // f-s
        if ( _lmax_row > 2){
             _ol(10,Cartesian::s) = _pma[0]*_ol(Cartesian::xx,Cartesian::s) + _fak2* _ol( 1,Cartesian::s);
             _ol(11,Cartesian::s) = _pma[1]*_ol(Cartesian::yy,Cartesian::s) + _fak2* _ol( 2,Cartesian::s);
             _ol(12,Cartesian::s) = _pma[2]*_ol(Cartesian::zz,Cartesian::s) + _fak2* _ol( 3,Cartesian::s);
             _ol(13,Cartesian::s) = _pma[0]*_ol(Cartesian::xy,Cartesian::s) + _fak * _ol( 2,Cartesian::s);
             _ol(14,Cartesian::s) = _pma[1]*_ol(Cartesian::xy,Cartesian::s) + _fak * _ol( 1,Cartesian::s);
             _ol(15,Cartesian::s) = _pma[0]*_ol(Cartesian::xy,Cartesian::s) + _fak * _ol( 3,Cartesian::s);
             _ol(16,Cartesian::s) = _pma[2]*_ol(Cartesian::xy,Cartesian::s) + _fak * _ol( 1,Cartesian::s);
             _ol(17,Cartesian::s) = _pma[1]*_ol(Cartesian::yz,Cartesian::s) + _fak * _ol( 3,Cartesian::s);
             _ol(18,Cartesian::s) = _pma[2]*_ol(Cartesian::yz,Cartesian::s) + _fak * _ol( 2,Cartesian::s);
             _ol(19,Cartesian::s) = _pma[2]*_ol(Cartesian::xy,Cartesian::s);
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
             _ol( 3,12) = _pma[2]*_ol( 0,12) + _fak3* _ol( 0,Cartesian::zz);
             _ol( 1,13) = _pma[0]*_ol( 0,13) + _fak2* _ol( 0, 4);
             _ol( 2,13) = _pma[1]*_ol( 0,13) + _fak * _ol( 0, 7);
             _ol( 3,13) = _pma[2]*_ol( 0,13);
             _ol( 1,14) = _pma[0]*_ol( 0,14) + _fak * _ol( 0, 8);
             _ol( 2,14) = _pma[1]*_ol( 0,14) + _fak2* _ol( 0, 4);
             _ol( 3,14) = _pma[2]*_ol( 0,14);
             _ol( 1,15) = _pma[0]*_ol( 0,15) + _fak2* _ol( 0, 5);
             _ol( 2,15) = _pma[1]*_ol( 0,15);
             _ol( 3,15) = _pma[2]*_ol( 0,15) + _fak * _ol( 0, 7);
             _ol( 1,16) = _pma[0]*_ol( 0,16) + _fak * _ol( 0,Cartesian::zz);
             _ol( 2,16) = _pma[1]*_ol( 0,16);
             _ol( 3,16) = _pma[2]*_ol( 0,16) + _fak2* _ol( 0, 5);
             _ol( 1,17) = _pma[0]*_ol( 0,17);
             _ol( 2,17) = _pma[1]*_ol( 0,17) + _fak2* _ol( 0, 6);
             _ol( 3,17) = _pma[2]*_ol( 0,17) + _fak * _ol( 0, 8);
             _ol( 1,18) = _pma[0]*_ol( 0,18);
             _ol( 2,18) = _pma[1]*_ol( 0,18) + _fak * _ol( 0,Cartesian::zz);
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
             _ol(12, 3) = _pma[2]*_ol( 9, 3) + _fak * (2.0*_ol( 3, 3) + _ol(Cartesian::zz, 0) );            
        }
        
        // d-f
        if ( _lmax_row > 1 && _lmax_col >2 ){
             _ol( 7,10) = _pma[0]*_ol( 1,10) + _fak * (_ol( 0,10) + 3.0*_ol( 1, 7) );
             _ol( 4,10) = _pma[1]*_ol( 1,10);
             _ol( 5,10) = _pma[2]*_ol( 1,10);
             _ol( 7,11) = _pma[0]*_ol( 1,11) + _fak * _ol( 0,Cartesian::x);
             _ol( 4,11) = _pma[1]*_ol( 1,11) + _fak3* _ol( 1, 8);
             _ol( 5,11) = _pma[2]*_ol( 1,11);
             _ol( 7,12) = _pma[0]*_ol( 1,12) + _fak * _ol( 0,12);
             _ol( 4,12) = _pma[1]*_ol( 1,12);
             _ol( 5,12) = _pma[2]*_ol( 1,12) + _fak3* _ol( 1,Cartesian::zz);
             _ol( 7,13) = _pma[0]*_ol( 1,13) + _fak * (_ol( 0,13) + 2.0*_ol( 1, 4) );
             _ol( 4,13) = _pma[1]*_ol( 1,13) + _fak * _ol( 1, 7);
             _ol( 5,13) = _pma[2]*_ol( 1,13);
             _ol( 7,14) = _pma[0]*_ol( 1,14) + _fak * (_ol( 0,14) + _ol( 1, 8) );
             _ol( 4,14) = _pma[1]*_ol( 1,14) + _fak2* _ol( 1, 4);
             _ol( 5,14) = _pma[2]*_ol( 1,14);
             _ol( 7,15) = _pma[0]*_ol( 1,15) + _fak * (_ol( 0,15) + 2.0*_ol( 1, 5) );
             _ol( 4,15) = _pma[1]*_ol( 1,15);
             _ol( 5,15) = _pma[2]*_ol( 1,15) + _fak * _ol( 1, 7);
             _ol( 7,16) = _pma[0]*_ol( 1,16) + _fak * (_ol( 0,16) + _ol( 1,Cartesian::zz) );
             _ol( 4,16) = _pma[1]*_ol( 1,16);
             _ol( 5,16) = _pma[2]*_ol( 1,16) + _fak2* _ol( 1, 5);
             _ol( 7,17) = _pma[0]*_ol( 1,17) + _fak * _ol( 0,17);
             _ol( 4,17) = _pma[1]*_ol( 1,17) + _fak2* _ol( 1, 6);
             _ol( 5,17) = _pma[2]*_ol( 1,17) + _fak * _ol( 1, 8);
             _ol( 7,18) = _pma[0]*_ol( 1,18) + _fak * _ol( 0,18);
             _ol( 4,18) = _pma[1]*_ol( 1,18) + _fak * _ol( 1,Cartesian::zz);
             _ol( 5,18) = _pma[2]*_ol( 1,18) + _fak2* _ol( 1, 6);
             _ol( 7,19) = _pma[0]*_ol( 1,19) + _fak * (_ol( 0,19) + _ol( 1, 6) );
             _ol( 4,19) = _pma[1]*_ol( 1,19) + _fak * _ol( 1, 5);
             _ol( 5,19) = _pma[2]*_ol( 1,19) + _fak * _ol( 1, 4);
             _ol( 8,10) = _pma[1]*_ol( 2,10) + _fak * _ol( 0,10);
             _ol( 6,10) = _pma[2]*_ol( 2,10);
             _ol( 8,11) = _pma[1]*_ol( 2,11) + _fak * (_ol( 0,11) + 3.0*_ol( 2, 8) );
             _ol( 6,11) = _pma[2]*_ol( 2,11);
             _ol( 8,12) = _pma[1]*_ol( 2,12) + _fak * _ol( 0,12);
             _ol( 6,12) = _pma[2]*_ol( 2,12) + _fak3* _ol( 2,Cartesian::zz);
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
             _ol( 8,18) = _pma[1]*_ol( 2,18) + _fak * (_ol( 0,18) + _ol( 2,Cartesian::zz) );
             _ol( 6,18) = _pma[2]*_ol( 2,18) + _fak2* _ol( 2, 6);
             _ol( 8,19) = _pma[1]*_ol( 2,19) + _fak * (_ol( 0,19) + _ol( 2, 5) );
             _ol( 6,19) = _pma[2]*_ol( 2,19) + _fak * _ol( 2, 4);
             _ol(Cartesian::zz,10) = _pma[2]*_ol( 3,10) + _fak * _ol( 0,10);
             _ol(Cartesian::zz,11) = _pma[2]*_ol( 3,11) + _fak * _ol( 0,11);
             _ol(Cartesian::zz,12) = _pma[2]*_ol( 3,12) + _fak * (_ol( 0,12) + 3.0*_ol( 3,Cartesian::zz) );
             _ol(Cartesian::zz,13) = _pma[2]*_ol( 3,13) + _fak * _ol( 0,13);
             _ol(Cartesian::zz,14) = _pma[2]*_ol( 3,14) + _fak * _ol( 0,14);
             _ol(Cartesian::zz,15) = _pma[2]*_ol( 3,15) + _fak * (_ol( 0,15) + _ol( 3, 7) );
             _ol(Cartesian::zz,16) = _pma[2]*_ol( 3,16) + _fak * (_ol( 0,16) + 2.0*_ol( 3, 5) );
             _ol(Cartesian::zz,17) = _pma[2]*_ol( 3,17) + _fak * (_ol( 0,17) + _ol( 3, 8) );
             _ol(Cartesian::zz,18) = _pma[2]*_ol( 3,18) + _fak * (_ol( 0,18) + 2.0*_ol( 3, 5) );
             _ol(Cartesian::zz,19) = _pma[2]*_ol( 3,19) + _fak * (_ol( 0,19) + _ol( 3, 4) );
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
             _ol(13,Cartesian::zz) = _pma[0]*_ol( 4,Cartesian::zz) + _fak * _ol( 2,Cartesian::zz);
             _ol(14,Cartesian::zz) = _pma[1]*_ol( 4,Cartesian::zz) + _fak * _ol( 1,Cartesian::zz);
             _ol(19,Cartesian::zz) = _pma[2]*_ol( 4,Cartesian::zz) + _fak2* _ol( 4, 3);
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
             _ol(15,Cartesian::zz) = _pma[0]*_ol( 5,Cartesian::zz) + _fak * _ol( 3,Cartesian::zz);
             _ol(16,Cartesian::zz) = _pma[2]*_ol( 5,Cartesian::zz) + _fak * (_ol( 1,Cartesian::zz) + 2.0*_ol( 5, 3) );
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
             _ol(17,Cartesian::zz) = _pma[1]*_ol( 6,Cartesian::zz) + _fak * _ol( 3,Cartesian::zz);
             _ol(18,Cartesian::zz) = _pma[2]*_ol( 6,Cartesian::zz) + _fak * (_ol( 2,Cartesian::zz) + 2.0*_ol( 6, 3) );
             _ol(10, 4) = _pma[0]*_ol( 7, 4) + _fak * (2.0*_ol( 1, 4) + _ol( 7, 2) );
             _ol(10, 5) = _pma[0]*_ol( 7, 5) + _fak * (2.0*_ol( 1, 5) + _ol( 7, 3) );
             _ol(10, 6) = _pma[0]*_ol( 7, 6) + _fak2* _ol( 1, 6);
             _ol(10, 7) = _pma[0]*_ol( 7, 7) + _fak * (2.0*_ol( 1, 7) + 2.0*_ol( 7, 1));
             _ol(10, 8) = _pma[0]*_ol( 7, 8) + _fak2* _ol( 1, 8);
             _ol(10,Cartesian::zz) = _pma[0]*_ol( 7,Cartesian::zz) + _fak2* _ol( 1,Cartesian::zz);
             _ol(11, 4) = _pma[1]*_ol( 8, 4) + _fak * (2.0*_ol( 2, 4) + _ol( 8, 1) );
             _ol(11, 5) = _pma[1]*_ol( 8, 5) + _fak2* _ol( 2, 5);
             _ol(11, 6) = _pma[1]*_ol( 8, 6) + _fak * (2.0*_ol( 2, 6) + _ol( 8, 3) );
             _ol(11, 7) = _pma[1]*_ol( 8, 7) + _fak2* _ol( 2, 7);
             _ol(11, 8) = _pma[1]*_ol( 8, 8) + _fak * (2.0*_ol( 2, 8) + 2.0*_ol( 8, 2));
             _ol(11,Cartesian::zz) = _pma[1]*_ol( 8,Cartesian::zz) + _fak2* _ol( 2,Cartesian::zz);
             _ol(12, 4) = _pma[2]*_ol(Cartesian::zz, 4) + _fak2* _ol( 3, 4);
             _ol(12, 5) = _pma[2]*_ol(Cartesian::zz, 5) + _fak * (2.0*_ol( 3, 5) + _ol(Cartesian::zz, 1) );
             _ol(12, 6) = _pma[2]*_ol(Cartesian::zz, 6) + _fak * (2.0*_ol( 3, 6) + _ol(Cartesian::zz, 2) );
             _ol(12, 7) = _pma[2]*_ol(Cartesian::zz, 7) + _fak2* _ol( 3, 7);
             _ol(12, 8) = _pma[2]*_ol(Cartesian::zz, 8) + _fak2* _ol( 3, 8);
             _ol(12,Cartesian::zz) = _pma[2]*_ol(Cartesian::zz,Cartesian::zz) + _fak * (2.0*_ol( 3,Cartesian::zz) + 2.0*_ol(Cartesian::zz, 3));
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
             _ol(19,12) = _pma[2]*_ol( 4,12) + _fak3* _ol( 4,Cartesian::zz);
             _ol(13,13) = _pma[0]*_ol( 4,13) + _fak * (_ol( 2,13) + 2.0*_ol( 4,Cartesian::xy) );
             _ol(14,13) = _pma[1]*_ol( 4,13) + _fak * (_ol( 1,13) + _ol( 4, 7) );
             _ol(19,13) = _pma[2]*_ol( 4,13);
             _ol(13,14) = _pma[0]*_ol( 4,14) + _fak * (_ol( 2,14) + _ol( 4, 8) );
             _ol(14,14) = _pma[1]*_ol( 4,14) + _fak * (_ol( 1,14) + 2.0*_ol( 4, 4) );
             _ol(19,14) = _pma[2]*_ol( 4,14);
             _ol(13,15) = _pma[0]*_ol( 4,15) + _fak * (_ol( 2,15) + 2.0*_ol( 4, 5) );
             _ol(14,15) = _pma[1]*_ol( 4,15) + _fak * _ol( 1,15);
             _ol(19,15) = _pma[2]*_ol( 4,15) + _fak * _ol( 4, 7);
             _ol(13,16) = _pma[0]*_ol( 4,16) + _fak * (_ol( 2,16) + _ol( 4,Cartesian::zz) );
             _ol(14,16) = _pma[1]*_ol( 4,16) + _fak * _ol( 1,16);
             _ol(19,16) = _pma[2]*_ol( 4,16) + _fak2* _ol( 4, 5);
             _ol(13,17) = _pma[0]*_ol( 4,17) + _fak * _ol( 2,17);
             _ol(14,17) = _pma[1]*_ol( 4,17) + _fak * (_ol( 1,17) + 2.0*_ol( 4, 6) );
             _ol(19,17) = _pma[2]*_ol( 4,17) + _fak * _ol( 4, 8);
             _ol(13,18) = _pma[0]*_ol( 4,18) + _fak * _ol( 2,18);
             _ol(14,18) = _pma[1]*_ol( 4,18) + _fak * (_ol( 1,18) + _ol( 4,Cartesian::zz) );
             _ol(19,18) = _pma[2]*_ol( 4,18) + _fak2* _ol( 4, 6);
             _ol(13,19) = _pma[0]*_ol( 4,19) + _fak * (_ol( 2,19) + _ol( 4, 6) );
             _ol(14,19) = _pma[1]*_ol( 4,19) + _fak * (_ol( 1,19) + _ol( 4, 5) );
             _ol(19,19) = _pma[2]*_ol( 4,19) + _fak * _ol( 4, 4);
             _ol(15,10) = _pma[0]*_ol( 5,10) + _fak * (_ol( 3,10) + 3.0*_ol( 5, 7) );
             _ol(16,10) = _pma[2]*_ol( 5,10) + _fak * _ol( 1,10);
             _ol(15,11) = _pma[0]*_ol( 5,11) + _fak * _ol( 3,11);
             _ol(16,11) = _pma[2]*_ol( 5,11) + _fak * _ol( 1,11);
             _ol(15,12) = _pma[0]*_ol( 5,12) + _fak * _ol( 3,12);
             _ol(16,12) = _pma[2]*_ol( 5,12) + _fak * (_ol( 1,12) + 3.0*_ol( 5,Cartesian::zz) );
             _ol(15,13) = _pma[0]*_ol( 5,13) + _fak * (_ol( 3,13) + 2.0*_ol( 5, 4) );
             _ol(16,13) = _pma[2]*_ol( 5,13) + _fak * _ol( 1,13);
             _ol(15,14) = _pma[0]*_ol( 5,14) + _fak * (_ol( 3,14) + _ol( 5, 8) );
             _ol(16,14) = _pma[2]*_ol( 5,14) + _fak * _ol( 1,14);
             _ol(15,15) = _pma[0]*_ol( 5,15) + _fak * (_ol( 3,15) + 2.0*_ol( 5, 5) );
             _ol(16,15) = _pma[2]*_ol( 5,15) + _fak * (_ol( 1,15) + _ol( 5, 7) );
             _ol(15,16) = _pma[0]*_ol( 5,16) + _fak * (_ol( 3,16) + _ol( 5,Cartesian::zz) );
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
             _ol(18,12) = _pma[2]*_ol( 6,12) + _fak * (_ol( 2,12) + 3.0*_ol( 6,Cartesian::zz) );
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
             _ol(17,18) = _pma[1]*_ol( 6,18) + _fak * (_ol( 3,18) + _ol( 6,Cartesian::zz) );
             _ol(18,18) = _pma[2]*_ol( 6,18) + _fak * (_ol( 2,18) + 2.0*_ol( 6, 6) );
             _ol(17,19) = _pma[1]*_ol( 6,19) + _fak * (_ol( 3,19) + _ol( 6, 5) );
             _ol(18,19) = _pma[2]*_ol( 6,19) + _fak * (_ol( 2,19) + _ol( 6, 4) );
             _ol(10,10) = _pma[0]*_ol( 7,10) + _fak * (2.0*_ol( 1,10) + 3.0*_ol( 7, 7));
             _ol(10,11) = _pma[0]*_ol( 7,11) + _fak2* _ol( 1,11);
             _ol(10,12) = _pma[0]*_ol( 7,12) + _fak2* _ol( 1,12);
             _ol(10,13) = _pma[0]*_ol( 7,13) + _fak * (2.0*_ol( 1,13) + 2.0*_ol( 7, 4));
             _ol(10,14) = _pma[0]*_ol( 7,14) + _fak * (2.0*_ol( 1,14) + _ol( 7, 8) );
             _ol(10,15) = _pma[0]*_ol( 7,15) + _fak * (2.0*_ol( 1,15) + 2.0*_ol( 7, 5));
             _ol(10,16) = _pma[0]*_ol( 7,16) + _fak * (2.0*_ol( 1,16) + _ol( 7,Cartesian::zz) );
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
             _ol(11,18) = _pma[1]*_ol( 8,18) + _fak * (2.0*_ol( 2,18) + _ol( 8,Cartesian::zz) );
             _ol(11,19) = _pma[1]*_ol( 8,19) + _fak * (2.0*_ol( 2,19) + _ol( 8, 5) );
             _ol(12,10) = _pma[2]*_ol(Cartesian::zz,10) + _fak2* _ol( 3,10);
             _ol(12,11) = _pma[2]*_ol(Cartesian::zz,11) + _fak2* _ol( 3,11);
             _ol(12,12) = _pma[2]*_ol(Cartesian::zz,12) + _fak * (2.0*_ol( 3,12) + 3.0*_ol(Cartesian::zz,Cartesian::zz));
             _ol(12,13) = _pma[2]*_ol(Cartesian::zz,13) + _fak2* _ol( 3,13);
             _ol(12,14) = _pma[2]*_ol(Cartesian::zz,14) + _fak2* _ol( 3,14);
             _ol(12,15) = _pma[2]*_ol(Cartesian::zz,15) + _fak * (2.0*_ol( 3,15) + _ol(Cartesian::zz, 7) );
             _ol(12,16) = _pma[2]*_ol(Cartesian::zz,16) + _fak * (2.0*_ol( 3,16) + 2.0*_ol(Cartesian::zz, 5));
             _ol(12,17) = _pma[2]*_ol(Cartesian::zz,17) + _fak * (2.0*_ol( 3,17) + _ol(Cartesian::zz, 8) );
             _ol(12,18) = _pma[2]*_ol(Cartesian::zz,18) + _fak * (2.0*_ol( 3,18) + 2.0*_ol(Cartesian::zz, 6));
             _ol(12,19) = _pma[2]*_ol(Cartesian::zz,19) + _fak * (2.0*_ol( 3,19) + _ol(Cartesian::zz, 4) );
        }
        // s-g
        // g-s
        // p-g
        // g-p
        // d-g
        // g-d
        // f-g
        // g-f
        // g-g
        
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
        for ( int i = 0; i< _matrix.size1(); i++ ) {
            for (int j = 0; j < _matrix.size2(); j++){
                _matrix(i,j) += _ol_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
            }
        }
        
        
        _ol.clear();
            } // _shell_col Gaussians
        } // _shell_row Gaussians
    }
    
  
        
    
    
}}

