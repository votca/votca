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

using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;

    
    void AOMomentum::FillBlock( std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix, AOShell* _shell_row, AOShell* _shell_col , AOBasis* ecp) {

        
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
        
        if ( _lmax_col > 2 ) {
            cerr << "Momentum transition dipoles only implemented for S,P,D functions in DFT basis!" << flush;
            exit(1);
        }

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
    
        // initialize local matrix block for unnormalized cartesians
        std::vector< ub::matrix<double> > _mom;
        for (int _i_comp = 0; _i_comp < 3; _i_comp++){
            _mom.push_back(ub::zero_matrix<double>(_nrows,_ncols));
        }
        
        // initialize local matrix block for unnormalized cartesians of overlap
        int _ncols_ol = this->getBlockSize( _lmax_col +1 ); 
        // make copy of shell_col and change type, lmax
        //AOShell _shell_col_local = (*_shell_col);
        
        ub::matrix<double> _ol = ub::zero_matrix<double>(_nrows,_ncols_ol);
        
        
        // get shell positions
        const vec& _pos_row = _shell_row->getPos();
        const vec& _pos_col = _shell_col->getPos();
        const vec  _diff    = _pos_row - _pos_col;
        vector<double> _pma (3,0.0);
        vector<double> _pmb (3,0.0);
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
        
        // some helpers
        
        
        const double _fak  = 0.5/(_decay_row + _decay_col);
        const double _fak2 = 2.0 * _fak;
        double _exparg = _fak2 * _decay_row * _decay_col *_distsq;
           
       /// check if distance between postions is big, then skip step   
       
        if ( _exparg > 30.0 ) { continue; }

        _pma[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
        _pma[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
        _pma[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

        _pmb[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
        _pmb[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
        _pmb[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
      

        
        
        const double _fak3 = 3.0 * _fak;
        //const double _fak4 = 4.0 * _fak;
        
        // calculate s-s- overlap matrix element
        _ol(0,0) = pow(4.0*_decay_row*_decay_col,0.75) * pow(_fak2,1.5)*exp(-_fak2 * _decay_row * _decay_col *_distsq); // s-s element

        

        // s-p momentum integrals
        if ( _lmax_col +1 > 0 ) {
     
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
        if ( _lmax_row > 0 && _lmax_col + 1 > 0 ) {
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
        if ( _lmax_col +1 > 1){
            //cout << "\t setting s-d" << endl;
            _ol(0,4) = _pmb[1]*_ol(0,1); // s-dxy
            _ol(0,5) = _pmb[2]*_ol(0,1); // s-dxz
            _ol(0,6) = _pmb[2]*_ol(0,2); // s-dyz
            _ol(0,7) = _pmb[0]*_ol(0,1) + _fak*_ol(0,0); // s-dxx
            _ol(0,8) = _pmb[1]*_ol(0,2) + _fak*_ol(0,0); // s-dyy
            _ol(0,9) = _pmb[2]*_ol(0,3) + _fak*_ol(0,0); // s-dzz
        }
        
        // p-d
        if ( _lmax_row > 0 && _lmax_col +1 > 1){
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
        if ( _lmax_row >1 && _lmax_col +1 > 0){
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
        if ( _lmax_row > 1 && _lmax_col +1 > 1 ){
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
        if ( _lmax_col +1 > 2 ){
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
        if ( _lmax_row > 0 && _lmax_col +1 > 2 ){
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
        if (_lmax_row > 2 && _lmax_col +1 > 0 ){
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
        if ( _lmax_row > 1 && _lmax_col +1 >2 ){
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
        if ( _lmax_row > 2 && _lmax_col +1 > 1 ){
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
        if ( _lmax_row > 2 && _lmax_col +1 > 2 ){
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
        // g-s
        // p-g
        // g-p
        // d-g
        // g-d
        // f-g
        // g-f
        // g-g
        
        // after overlap matrix is prepared, we can get the gradient/momentum matrix
       
        /* matrix storage _matrix[m](i,j)
         * with:
         *    m: cartesian component of momentum 0-x; 1-y; 2-z
         *    i: row basis function
         *    j: col basis function
         */

            double _two_beta = 2.0 * _decay_col;
            for (int _i_row = 0; _i_row < _nrows; _i_row++) {

                // x-component
                _mom[0](_i_row, 0) = -_two_beta * _ol(_i_row, 1);
                if (_lmax_col > 0) {
                    _mom[0](_i_row, 1) = _ol(_i_row, 0) - _two_beta * _ol(_i_row, 7);
                    _mom[0](_i_row, 2) = -_two_beta * _ol(_i_row, 4);
                    _mom[0](_i_row, 3) = -_two_beta * _ol(_i_row, 5);
                }
                if (_lmax_col > 1) {
                    _mom[0](_i_row, 4) = _ol(_i_row, 2) - _two_beta * _ol(_i_row, 13);
                    _mom[0](_i_row, 5) = _ol(_i_row, 3) - _two_beta * _ol(_i_row, 15);
                    _mom[0](_i_row, 6) = -_two_beta * _ol(_i_row, 19);
                    _mom[0](_i_row, 7) = 2.0 * _ol(_i_row, 1) - _two_beta * _ol(_i_row, 10);
                    _mom[0](_i_row, 8) = -_two_beta * _ol(_i_row, 14);
                    _mom[0](_i_row, 9) = -_two_beta * _ol(_i_row, 16);
                }

                // y-component
                _mom[1](_i_row, 0) = -_two_beta * _ol(_i_row, 2);
                if (_lmax_col > 0) {
                    _mom[1](_i_row, 1) = -_two_beta * _ol(_i_row, 4);
                    _mom[1](_i_row, 2) = _ol(_i_row, 0) - _two_beta * _ol(_i_row, 8);
                    _mom[1](_i_row, 3) = -_two_beta * _ol(_i_row, 6);
                }
                if (_lmax_col > 1) {
                    _mom[1](_i_row, 4) = _ol(_i_row, 1) - _two_beta * _ol(_i_row, 14);
                    _mom[1](_i_row, 5) = -_two_beta * _ol(_i_row, 19);
                    _mom[1](_i_row, 6) = _ol(_i_row, 3) - _two_beta * _ol(_i_row, 17);
                    _mom[1](_i_row, 7) = -_two_beta * _ol(_i_row, 13);
                    _mom[1](_i_row, 8) = 2.0 * _ol(_i_row, 2) - _two_beta * _ol(_i_row, 11);
                    _mom[1](_i_row, 9) = -_two_beta * _ol(_i_row, 18);
                }

                
                // z-component
                _mom[2](_i_row, 0) = -_two_beta * _ol(_i_row, 3);
                if (_lmax_col > 0) {
                    _mom[2](_i_row, 1) = -_two_beta * _ol(_i_row, 5);
                    _mom[2](_i_row, 2) = -_two_beta * _ol(_i_row, 6);
                    _mom[2](_i_row, 3) = _ol(_i_row, 0) - _two_beta * _ol(_i_row, 9);
                }
                if (_lmax_col > 1) {
                    _mom[2](_i_row, 4) = -_two_beta * _ol(_i_row, 19);
                    _mom[2](_i_row, 5) = _ol(_i_row, 1) - _two_beta * _ol(_i_row, 16);
                    _mom[2](_i_row, 6) = _ol(_i_row, 2) - _two_beta * _ol(_i_row, 18);
                    _mom[2](_i_row, 7) = -_two_beta * _ol(_i_row, 15);
                    _mom[2](_i_row, 8) = -_two_beta * _ol(_i_row, 17);
                    _mom[2](_i_row, 9) = 2.0 * _ol(_i_row, 3) - _two_beta * _ol(_i_row, 12);
                }

            }
        
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
        ub::matrix<double> _trafo_col_tposed = ub::trans( _trafo_col );

        // cartesian -> spherical
       
        for ( int _i_comp = 0; _i_comp < 3; _i_comp++){

            ub::matrix<double> _mom_tmp = ub::prod( _trafo_row, _mom[ _i_comp ] );

            ub::matrix<double> _mom_sph = ub::prod( _mom_tmp, _trafo_col_tposed );
            
            // save to _matrix
            for ( unsigned i = 0; i< _matrix[0].size1(); i++ ) {
                for (unsigned j = 0; j < _matrix[0].size2(); j++){
                    _matrix[ _i_comp ](i,j) += _mom_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
                }
            }
        }
        
        _ol.clear();
            }// _shell_col Gaussians
        }// _shell_row Gaussians
    }
    
  
        
    
    
}}

