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

    
    void AODipole::FillBlock( std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix, AOShell* _shell_row, AOShell* _shell_col , AOBasis* ecp) {

        
        /* Calculating the AO matrix of the gradient operator requires 
         * the raw overlap matrix (i.e. in unnormalized cartesians) 
        
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
        std::vector< ub::matrix<double> > _dip;
        for (int _i_comp = 0; _i_comp < 3; _i_comp++){
            _dip.push_back(ub::zero_matrix<double>(_nrows,_ncols));
        }
        
        // initialize local matrix block for unnormalized cartesians of overlap
        // int _ncols_ol = this->getBlockSize( _lmax_col +1 ); 
        
        ub::matrix<double> _ol = ub::zero_matrix<double>(_nrows,_ncols);
        
         // get shell positions
        const vec& _pos_row = _shell_row->getPos();
        const vec& _pos_col = _shell_col->getPos();
        const vec  _diff    = _pos_row - _pos_col;
        double _distsq = (_diff.getX()*_diff.getX()) + (_diff.getY()*_diff.getY()) + (_diff.getZ()*_diff.getZ()); 
         
        
        
        // some helpers
       
        vector<double> _pma (3,0.0);
        vector<double> _pmb (3,0.0);
        // definition of a center around which the moment should be calculated
        vector<double> _center(3,0.0); // here: origin, can be changed later
        vector<double> _pmc(3,0.0);
        
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
        
        //const double _fak3 = 3.0 * _fak;
        //const double _fak4 = 4.0 * _fak;

        _pma[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
        _pma[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
        _pma[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

        _pmb[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
        _pmb[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
        _pmb[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
        
        
        _pmc[0] = _fak2*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _center[0];
        _pmc[1] = _fak2*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _center[1];
        _pmc[2] = _fak2*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _center[2];
        
        
        

        
        
        // calculate s-s- overlap matrix element
        _ol(0,0) = pow(4.0*_decay_row*_decay_col,0.75) * pow(_fak2,1.5)*exp(-_fak2 * _decay_row * _decay_col *_distsq); // s-s element

        // s-s dipole moment integrals
        for ( int _i_comp = 0 ; _i_comp < 3; _i_comp++ ){
            _dip[_i_comp](0,0) = _pmc[_i_comp]*_ol(0,0);
        }
        
        // s-p dipole integrals
        if ( _lmax_col > 0 ) {
     
            // overlap
            _ol(0,1) = _pmb[0]*_ol(0,0); // s-px
            _ol(0,2) = _pmb[1]*_ol(0,0); // s-py
            _ol(0,3) = _pmb[2]*_ol(0,0); // s-pz
            
            // dipole 
            // x-components
            _dip[0](0,1) = _pmb[0]*_dip[0](0,0) + _fak * _ol(0,0); // s-x-px
            _dip[0](0,2) = _pmb[1]*_dip[0](0,0) ; // s-x-py
            _dip[0](0,3) = _pmb[2]*_dip[0](0,0) ; // s-x-pz

            // y-components
            _dip[1](0,1) = _pmb[0]*_dip[1](0,0); // s-y-px
            _dip[1](0,2) = _pmb[1]*_dip[1](0,0) + _fak * _ol(0,0); // s-y-py
            _dip[1](0,3) = _pmb[2]*_dip[1](0,0); // s-y-pz

            // z-components
            _dip[2](0,1) = _pmb[0]*_dip[2](0,0); // s-z-px
            _dip[2](0,2) = _pmb[1]*_dip[2](0,0); // s-z-py
            _dip[2](0,3) = _pmb[2]*_dip[2](0,0) + _fak * _ol(0,0); // s-z-pz
            
            
        }
        
        // p-s
        if ( _lmax_row > 0 ) {
           //cout << "\t setting p-s" << flush;

           _ol(1,0) = _pma[0]*_ol(0,0); // px-s
           _ol(2,0) = _pma[1]*_ol(0,0); // py-s
           _ol(3,0) = _pma[2]*_ol(0,0); // pz-s


            // dipole 
            // x-components
            _dip[0](1,0) = _pma[0]*_dip[0](0,0) + _fak * _ol(0,0); // px-x-s
            _dip[0](2,0) = _pma[1]*_dip[0](0,0) ; // py-x-s
            _dip[0](3,0) = _pma[2]*_dip[0](0,0) ; // pz-x-s

            // y-components
            _dip[1](1,0) = _pma[0]*_dip[1](0,0); // px-y-s
            _dip[1](2,0) = _pma[1]*_dip[1](0,0) + _fak * _ol(0,0); // py-y-s
            _dip[1](3,0) = _pma[2]*_dip[1](0,0); // pz-y-s

            // z-components
            _dip[2](1,0) = _pma[0]*_dip[2](0,0); // px-z-s
            _dip[2](2,0) = _pma[1]*_dip[2](0,0); // py-z-s
            _dip[2](3,0) = _pma[2]*_dip[2](0,0) + _fak * _ol(0,0); // pz-z-s

        }
        
        // p-p
        if ( _lmax_row > 0 && _lmax_col  > 0 ) {
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


           // dipole
           // x-components
           _dip[0](1,1) = _pmb[0]*_dip[0](1,0) + _fak * (_dip[0](0,0) + _ol(1,0)); // px-x-px, ok
           _dip[0](1,2) = _pmb[1]*_dip[0](1,0); // px-x-py, ok
           _dip[0](1,3) = _pmb[2]*_dip[0](1,0); //px-x-pz, ok
           
           _dip[0](2,1) = _pmb[0]*_dip[0](2,0) + _fak * _ol(2,0); // py-x-px, ok
           _dip[0](2,2) = _pmb[1]*_dip[0](2,0) + _fak * _dip[0](0,0); //py-x-py, ok
           _dip[0](2,3) = _pmb[2]*_dip[0](2,0); // py-x-pz, ok
           
           _dip[0](3,1) = _pmb[0]*_dip[0](3,0) + _fak * _ol(3,0); // pz-x-px, ok
           _dip[0](3,2) = _pmb[1]*_dip[0](3,0); // pz-x-py, ok
           _dip[0](3,3) = _pmb[2]*_dip[0](3,0) + _fak * _dip[0](0,0); // pz-x-pz, ok
           
           // y-components
           _dip[1](1,1) = _pmb[0]*_dip[1](1,0) + _fak * _dip[1](0,0); // px-y-px, ok
           _dip[1](1,2) = _pmb[1]*_dip[1](1,0) + _fak * _ol(1,0) ; // px-y-py, ok
           _dip[1](1,3) = _pmb[2]*_dip[1](1,0); //px-y-pz, ok
           
           _dip[1](2,1) = _pmb[0]*_dip[1](2,0); // py-y-px, ok
           _dip[1](2,2) = _pmb[1]*_dip[1](2,0) + _fak * ( _dip[1](0,0) + _ol(2,0)); //py-y-py, ok
           _dip[1](2,3) = _pmb[2]*_dip[1](2,0); // py-y-pz, ok
           
           _dip[1](3,1) = _pmb[0]*_dip[1](3,0); // pz-y-px, ok
           _dip[1](3,2) = _pmb[1]*_dip[1](3,0) + _fak * _ol(3,0); // pz-y-py,ok
           _dip[1](3,3) = _pmb[2]*_dip[1](3,0) + _fak * _dip[1](0,0); // pz-y-pz,ok
           
                      
            // z-components
           _dip[2](1,1) = _pmb[0]*_dip[2](1,0) + _fak * _dip[2](0,0) ; // px-z-px, ok
           _dip[2](1,2) = _pmb[1]*_dip[2](1,0); // px-z-py, ok
           _dip[2](1,3) = _pmb[2]*_dip[2](1,0) + _fak * _ol(1,0); //px-z-pz, ok
           
           _dip[2](2,1) = _pmb[0]*_dip[2](2,0); // py-z-px, ok
           _dip[2](2,2) = _pmb[1]*_dip[2](2,0) + _fak * _dip[2](0,0); //py-z-py, ok
           _dip[2](2,3) = _pmb[2]*_dip[2](2,0) + _fak * _ol(2,0); // py-z-pz,ok
           
           _dip[2](3,1) = _pmb[0]*_dip[2](3,0); // pz-z-px, ok
           _dip[2](3,2) = _pmb[1]*_dip[2](3,0); // pz-z-py, ok
           _dip[2](3,3) = _pmb[2]*_dip[2](3,0) + _fak * ( _dip[2](0,0) + _ol(3,0)); // pz-z-pz,ok

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
            
            
            // x-components
            _dip[0](0,4) = _pmb[1]*_dip[0](0,1);// s-x-dxy, ok
            _dip[0](0,5) = _pmb[2]*_dip[0](0,1);// s-x-dxz, ok
            _dip[0](0,6) = _pmb[2]*_dip[0](0,2);// s-x-dyz, ok
            _dip[0](0,7) = _pmb[0]*_dip[0](0,1) + _fak * ( _dip[0](0,0) + _ol(0,1)  );// s-x-dxx, ok
            _dip[0](0,8) = _pmb[1]*_dip[0](0,2) + _fak * _dip[0](0,0);// s-x-dyy, ok
            _dip[0](0,9) = _pmb[2]*_dip[0](0,3) + _fak * _dip[0](0,0);// s-x-dzz, ok
            
            // y-components
            _dip[1](0,4) = _pmb[1]*_dip[1](0,1) + _fak * _ol(0,1);// s-y-dxy, ok
            _dip[1](0,5) = _pmb[2]*_dip[1](0,1);// s-y-dxz, ok
            _dip[1](0,6) = _pmb[2]*_dip[1](0,2);// s-y-dyz, ok
            _dip[1](0,7) = _pmb[0]*_dip[1](0,1) + _fak * _dip[1](0,0) ;// s-y-dxx, ok
            _dip[1](0,8) = _pmb[1]*_dip[1](0,2) + _fak * (_dip[1](0,0) + _ol(0,2)) ;// s-y-dyy, ok
            _dip[1](0,9) = _pmb[2]*_dip[1](0,3) + _fak * _dip[1](0,0);// s-y-dzz, ok           
            
            // z-components
            _dip[2](0,4) = _pmb[1]*_dip[2](0,1);// s-z-dxy, ok
            _dip[2](0,5) = _pmb[2]*_dip[2](0,1) + _fak * _ol(0,1);// s-z-dxz, ok
            _dip[2](0,6) = _pmb[2]*_dip[2](0,2) + _fak * _ol(0,2);// s-z-dyz, ok
            _dip[2](0,7) = _pmb[0]*_dip[2](0,1) + _fak * _dip[2](0,0) ;// s-z-dxx, ok
            _dip[2](0,8) = _pmb[1]*_dip[2](0,2) + _fak * _dip[2](0,0) ;// s-z-dyy, ok
            _dip[2](0,9) = _pmb[2]*_dip[2](0,3) + _fak * (_dip[2](0,0) + _ol(0,3) );// s-z-dzz, ok
            
            
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
            
            // x-components
            _dip[0](1,4) = _pma[0] * _dip[0](0,4) + _fak * ( _dip[0](0,2) + _ol(0,4) ); // px-x-dxy, ok
            _dip[0](1,5) = _pma[0] * _dip[0](0,5) + _fak * ( _dip[0](0,3) + _ol(0,5) ); // px-x-dxz, ok
            _dip[0](1,6) = _pma[0] * _dip[0](0,6) + _fak * _ol(0,6) ; // px-x-dyz, ok
            _dip[0](1,7) = _pma[0] * _dip[0](0,7) + _fak * ( _dip[0](0,1) + _dip[0](0,1) + _ol(0,7)); // px-x-dxx, ok
            _dip[0](1,8) = _pma[0] * _dip[0](0,8) + _fak * _ol(0,8); // px-x-dyy, ok
            _dip[0](1,9) = _pma[0] * _dip[0](0,9) + _fak * _ol(0,9); // px-x-dzz, ok
                    
            _dip[0](2,4) = _pma[1] * _dip[0](0,4) + _fak * _dip[0](0,1); // py-x-dxy, ok
            _dip[0](2,5) = _pma[1] * _dip[0](0,5); // py-x-dxz, ok
            _dip[0](2,6) = _pma[1] * _dip[0](0,6) + _fak * _dip[0](0,3) ; // py-x-dyz,ok
            _dip[0](2,7) = _pma[1] * _dip[0](0,7); // py-x-dxx, ok
            _dip[0](2,8) = _pma[1] * _dip[0](0,8) + _fak * ( _dip[0](0,2) *2.0 ); // py-x-dyy, ok
            _dip[0](2,9) = _pma[1] * _dip[0](0,9); // py-x-dzz, ok
            
            _dip[0](3,4) = _pma[2] * _dip[0](0,4); // pz-x-dxy, ok
            _dip[0](3,5) = _pma[2] * _dip[0](0,5) + _fak * _dip[0](0,1); // pz-x-dxz, ok
            _dip[0](3,6) = _pma[2] * _dip[0](0,6) + _fak * _dip[0](0,2) ; // pz-x-dyz, ok
            _dip[0](3,7) = _pma[2] * _dip[0](0,7); // pz-x-dxx, ok
            _dip[0](3,8) = _pma[2] * _dip[0](0,8); // pz-x-dyy, ok
            _dip[0](3,9) = _pma[2] * _dip[0](0,9) + _fak * ( _dip[0](0,3) * 2.0 ); // pz-x-dzz, ok
            
            // y-components
            _dip[1](1,4) = _pma[0] * _dip[1](0,4) + _fak * _dip[1](0,2); // px-y-dxy, ok
            _dip[1](1,5) = _pma[0] * _dip[1](0,5) + _fak * _dip[1](0,3); // px-y-dxz, ok
            _dip[1](1,6) = _pma[0] * _dip[1](0,6); // px-y-dyz, ok
            _dip[1](1,7) = _pma[0] * _dip[1](0,7) + _fak * ( _dip[1](0,1) + _dip[1](0,1) ); // px-y-dxx, ok
            _dip[1](1,8) = _pma[0] * _dip[1](0,8); // px-y-dyy, ok
            _dip[1](1,9) = _pma[0] * _dip[1](0,9); // px-y-dzz, ok
                    
            _dip[1](2,4) = _pma[1] * _dip[1](0,4) + _fak * ( _dip[1](0,1) + _ol(0,4) ); // py-y-dxy, ok
            _dip[1](2,5) = _pma[1] * _dip[1](0,5) + _fak * _ol(0,5); // py-y-dxz, ok
            _dip[1](2,6) = _pma[1] * _dip[1](0,6) + _fak * ( _dip[1](0,3) + _ol(0,6)) ; // py-y-dyz, ok
            _dip[1](2,7) = _pma[1] * _dip[1](0,7) + _fak * _ol(0,7); // py-y-dxx, ok
            _dip[1](2,8) = _pma[1] * _dip[1](0,8) + _fak * ( _dip[1](0,2) *2.0 + _ol(0,8) ); // py-y-dyy, ok
            _dip[1](2,9) = _pma[1] * _dip[1](0,9) + _fak * _ol(0,9); // py-y-dzz, ok
            
            _dip[1](3,4) = _pma[2] * _dip[1](0,4); // pz-y-dxy, ok
            _dip[1](3,5) = _pma[2] * _dip[1](0,5) + _fak * _dip[1](0,1); // pz-y-dxz, ok
            _dip[1](3,6) = _pma[2] * _dip[1](0,6) + _fak * _dip[1](0,2) ; // pz-y-dyz, ok
            _dip[1](3,7) = _pma[2] * _dip[1](0,7); // pz-y-dxx, ok
            _dip[1](3,8) = _pma[2] * _dip[1](0,8); // pz-y-dyy, ok 
            _dip[1](3,9) = _pma[2] * _dip[1](0,9) + _fak * ( _dip[1](0,3) * 2.0 ); // pz-y-dzz, ok
            


            // z-components
            _dip[2](1,4) = _pma[0] * _dip[2](0,4) + _fak * _dip[2](0,2); // px-z-dxy, ok
            _dip[2](1,5) = _pma[0] * _dip[2](0,5) + _fak * _dip[2](0,3); // px-z-dxz, ok
            _dip[2](1,6) = _pma[0] * _dip[2](0,6); // px-z-dyz, ok
            _dip[2](1,7) = _pma[0] * _dip[2](0,7) + _fak * ( _dip[2](0,1) + _dip[2](0,1) ); // px-z-dxx, ok
            _dip[2](1,8) = _pma[0] * _dip[2](0,8); // px-z-dyy, ok
            _dip[2](1,9) = _pma[0] * _dip[2](0,9); // px-z-dzz, ok
            
            _dip[2](2,4) = _pma[1] * _dip[2](0,4) + _fak *  _dip[2](0,1); // py-z-dxy, ok
            _dip[2](2,5) = _pma[1] * _dip[2](0,5); // py-z-dxz, ok
            _dip[2](2,6) = _pma[1] * _dip[2](0,6) + _fak *  _dip[2](0,3); // py-z-dyz, ok
            _dip[2](2,7) = _pma[1] * _dip[2](0,7); // py-z-dxx, ok
            _dip[2](2,8) = _pma[1] * _dip[2](0,8) + _fak * ( _dip[2](0,2) *2.0  ); // py-z-dyy, ok
            _dip[2](2,9) = _pma[1] * _dip[2](0,9); // py-z-dzz, ok
            
            _dip[2](3,4) = _pma[2] * _dip[2](0,4) + _fak * _ol(0,4); // pz-z-dxy, ok
            _dip[2](3,5) = _pma[2] * _dip[2](0,5) + _fak * (_dip[2](0,1) + _ol(0,5)); // pz-z-dxz, ok
            _dip[2](3,6) = _pma[2] * _dip[2](0,6) + _fak * (_dip[2](0,2) + _ol(0,6)); // pz-z-dyz, ok
            _dip[2](3,7) = _pma[2] * _dip[2](0,7) + _fak * _ol(0,7); // pz-z-dxx, ok
            _dip[2](3,8) = _pma[2] * _dip[2](0,8) + _fak * _ol(0,8); // pz-z-dyy, ok
            _dip[2](3,9) = _pma[2] * _dip[2](0,9) + _fak * ( _dip[2](0,3) * 2.0 + _ol(0,9)); // pz-z-dzz, ok
            
            
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
            
            // x-components
            _dip[0](4,0) = _pma[1] * _dip[0](1,0) ; // dxy-x-s, ok
            _dip[0](5,0) = _pma[2] * _dip[0](1,0) ; // dxz-x-s, ok
            _dip[0](6,0) = _pma[2] * _dip[0](2,0) ; // dyz-x-s, ok
            _dip[0](7,0) = _pma[0] * _dip[0](1,0) + _fak * ( _dip[0](0,0) + _ol(1,0) ); // dxx-x-s, ok
            _dip[0](8,0) = _pma[1] * _dip[0](2,0) + _fak * _dip[0](0,0); // dyy-x-s, ok
            _dip[0](9,0) = _pma[2] * _dip[0](3,0) + _fak * _dip[0](0,0); // dzz-x-s, ok
            
            
            // y-components
            _dip[1](4,0) = _pma[1] * _dip[1](1,0) + _fak * _ol(1,0); // dxy-y-s, ok
            _dip[1](5,0) = _pma[2] * _dip[1](1,0) ; // dxz-y-s, ok
            _dip[1](6,0) = _pma[2] * _dip[1](2,0) ; // dyz-y-s, ok
            _dip[1](7,0) = _pma[0] * _dip[1](1,0) + _fak * _dip[1](0,0); // dxx-y-s, ok
            _dip[1](8,0) = _pma[1] * _dip[1](2,0) + _fak * ( _dip[1](0,0) + _ol(2,0) ); // dyy-y-s, ok
            _dip[1](9,0) = _pma[2] * _dip[1](3,0) + _fak * _dip[1](0,0); // dzz-y-s, ok
            
             // z-components
            _dip[2](4,0) = _pma[1] * _dip[2](1,0); // dxy-z-s, ok
            _dip[2](5,0) = _pma[2] * _dip[2](1,0) + _fak * _ol(1,0) ; // dxz-z-s, ok
            _dip[2](6,0) = _pma[2] * _dip[2](2,0) + _fak * _ol(2,0); // dyz-z-s, ok
            _dip[2](7,0) = _pma[0] * _dip[2](1,0) + _fak * _dip[2](0,0); // dxx-z-s, ok
            _dip[2](8,0) = _pma[1] * _dip[2](2,0) + _fak * _dip[2](0,0); // dyy-z-s, ok
            _dip[2](9,0) = _pma[2] * _dip[2](3,0) + _fak * (_dip[2](0,0) + _ol(3,0)); // dzz-z-s, ok
            
            
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

             
             // x-components
             _dip[0](4,1) = _pmb[0] * _dip[0](4,0) + _fak * ( _dip[0](2,0) + _ol(4,0) ); // dxy-x-px, ok
             _dip[0](5,1) = _pmb[0] * _dip[0](5,0) + _fak * ( _dip[0](3,0) + _ol(5,0) ); // dxz-x-px, ok
             _dip[0](6,1) = _pmb[0] * _dip[0](6,0) + _fak * _ol(6,0); // dyz-x-px, ok
             _dip[0](7,1) = _pmb[0] * _dip[0](7,0) + _fak * ( 2.0*_dip[0](1,0) + _ol(7,0) ); // dxx-x-px, ok
             _dip[0](8,1) = _pmb[0] * _dip[0](8,0) + _fak * _ol(8,0); // dyy-x-px, ok
             _dip[0](9,1) = _pmb[0] * _dip[0](9,0) + _fak * _ol(9,0); // dzz-x-px, ok

             _dip[0](4,2) = _pmb[1] * _dip[0](4,0) + _fak * _dip[0](1,0) ; // dxy-x-py, ok
             _dip[0](5,2) = _pmb[1] * _dip[0](5,0); // dxz-x-py, ok
             _dip[0](6,2) = _pmb[1] * _dip[0](6,0) + _fak * _dip[0](3,0); // dyz-x-py, ok
             _dip[0](7,2) = _pmb[1] * _dip[0](7,0); // dxx-x-py, ok
             _dip[0](8,2) = _pmb[1] * _dip[0](8,0) + _fak *2.0* _dip[0](2,0); // dyy-x-py, ok
             _dip[0](9,2) = _pmb[1] * _dip[0](9,0); // dzz-x-py, ok
             
             _dip[0](4,3) = _pmb[2] * _dip[0](4,0); // dxy-x-pz, ok
             _dip[0](5,3) = _pmb[2] * _dip[0](5,0) + _fak * _dip[0](1,0) ; // dxz-x-pz, ok
             _dip[0](6,3) = _pmb[2] * _dip[0](6,0) + _fak * _dip[0](2,0) ; // dyz-x-pz, ok
             _dip[0](7,3) = _pmb[2] * _dip[0](7,0); // dxx-x-pz, ok
             _dip[0](8,3) = _pmb[2] * _dip[0](8,0); // dyy-x-pz, ok
             _dip[0](9,3) = _pmb[2] * _dip[0](9,0) + _fak *  2.0*_dip[0](3,0) ; // dzz-x-pz, ok

             
             // y-components
             _dip[1](4,1) = _pmb[0] * _dip[1](4,0) + _fak * _dip[1](2,0); // dxy-y-px. ok
             _dip[1](5,1) = _pmb[0] * _dip[1](5,0) + _fak * _dip[1](3,0); // dxz-y-px, ok
             _dip[1](6,1) = _pmb[0] * _dip[1](6,0); // dyz-y-px, ok
             _dip[1](7,1) = _pmb[0] * _dip[1](7,0) + _fak * 2.0 * _dip[1](1,0); // dxx-y-px, ok
             _dip[1](8,1) = _pmb[0] * _dip[1](8,0); // dyy-y-px, ok
             _dip[1](9,1) = _pmb[0] * _dip[1](9,0); // dzz-y-px, ok

             _dip[1](4,2) = _pmb[1] * _dip[1](4,0) + _fak * ( _dip[1](1,0) + _ol(4,0) ) ; // dxy-y-py, ok
             _dip[1](5,2) = _pmb[1] * _dip[1](5,0) + _fak * _ol(5,0); // dxz-y-py, ok
             _dip[1](6,2) = _pmb[1] * _dip[1](6,0) + _fak * ( _dip[1](3,0) + _ol(6,0) ) ; // dyz-y-py, ok
             _dip[1](7,2) = _pmb[1] * _dip[1](7,0) + _fak * _ol(7,0); // dxx-y-py, ok
             _dip[1](8,2) = _pmb[1] * _dip[1](8,0) + _fak * (2.0* _dip[1](2,0) + _ol(8,0)) ; // dyy-y-py, ok
             _dip[1](9,2) = _pmb[1] * _dip[1](9,0) + _fak * _ol(9,0); // dzz-y-py, ok
             
             
             _dip[1](4,3) = _pmb[2] * _dip[1](4,0); // dxy-y-pz, ok
             _dip[1](5,3) = _pmb[2] * _dip[1](5,0) + _fak * _dip[1](1,0) ; // dxz-y-pz, ok
             _dip[1](6,3) = _pmb[2] * _dip[1](6,0) + _fak * _dip[1](2,0) ; // dyz-y-pz, ok
             _dip[1](7,3) = _pmb[2] * _dip[1](7,0); // dxx-y-pz, ok
             _dip[1](8,3) = _pmb[2] * _dip[1](8,0); // dyy-y-pz, ok
             _dip[1](9,3) = _pmb[2] * _dip[1](9,0) + _fak *  2.0*_dip[1](3,0) ; // dzz-y-pz, ok



             // z-components
             _dip[2](4,1) = _pmb[0] * _dip[2](4,0) + _fak * _dip[2](2,0); // dxy-z-px, ok
             _dip[2](5,1) = _pmb[0] * _dip[2](5,0) + _fak * _dip[2](3,0); // dxz-z-px, ok
             _dip[2](6,1) = _pmb[0] * _dip[2](6,0); // dyz-z-px, ok
             _dip[2](7,1) = _pmb[0] * _dip[2](7,0) + _fak * 2.0 * _dip[2](1,0); // dxx-z-px, ok
             _dip[2](8,1) = _pmb[0] * _dip[2](8,0); // dyy-z-px, ok
             _dip[2](9,1) = _pmb[0] * _dip[2](9,0); // dzz-z-px, ok

             _dip[2](4,2) = _pmb[1] * _dip[2](4,0) + _fak * _dip[2](1,0) ; // dxy-z-py, ok
             _dip[2](5,2) = _pmb[1] * _dip[2](5,0); // dxz-z-py, ok
             _dip[2](6,2) = _pmb[1] * _dip[2](6,0) + _fak * _dip[2](3,0) ; // dyz-z-py, ok
             _dip[2](7,2) = _pmb[1] * _dip[2](7,0); // dxx-z-py, ok
             _dip[2](8,2) = _pmb[1] * _dip[2](8,0) + _fak * 2.0* _dip[2](2,0) ; // dyy-z-py, ok
             _dip[2](9,2) = _pmb[1] * _dip[2](9,0); // dzz-z-py, ok
             
             _dip[2](4,3) = _pmb[2] * _dip[2](4,0) + _fak * _ol(4,0); // dxy-z-pz, ok
             _dip[2](5,3) = _pmb[2] * _dip[2](5,0) + _fak * ( _dip[2](1,0) + _ol(5,0) ) ; // dxz-z-pz, ok
             _dip[2](6,3) = _pmb[2] * _dip[2](6,0) + _fak * ( _dip[2](2,0) + _ol(6,0) ); // dyz-z-pz, ok
             _dip[2](7,3) = _pmb[2] * _dip[2](7,0) + _fak * _ol(7,0) ; // dxx-z-pz, ok
             _dip[2](8,3) = _pmb[2] * _dip[2](8,0) + _fak * _ol(8,0) ; // dyy-z-pz, ok
             _dip[2](9,3) = _pmb[2] * _dip[2](9,0) + _fak * ( 2.0*_dip[2](3,0) + _ol(9,0) ) ; // dzz-z-pz, ok

             
        }
        
        // d-d
        if ( _lmax_row > 1 && _lmax_col  > 1 ){
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
            
             
             // x-components
             _dip[0](4,4) = _pmb[1] * _dip[0](4,1) + _fak *  _dip[0](1,1)  ; // dxy-x-dxy, ok
             _dip[0](4,5) = _pmb[2] * _dip[0](4,1); // dxy-x-dxz, ok
             _dip[0](4,6) = _pmb[2] * _dip[0](4,2); // dxy-x-dyz, ok
             _dip[0](4,7) = _pmb[0] * _dip[0](4,1) + _fak * ( _dip[0](2,1) + _dip[0](4,0) + _ol(4,1) ); // dxy-x-dxx, ok 
             _dip[0](4,8) = _pmb[1] * _dip[0](4,2) + _fak * ( _dip[0](1,2) + _dip[0](4,0) ); // dxy-x-dyy, ok
             _dip[0](4,9) = _pmb[2] * _dip[0](4,3) + _fak * _dip[0](4,0) ; // dxy-x-dzz, ok

             _dip[0](5,4) = _pmb[1] * _dip[0](5,1); // dxz-x-dxy, ok
             _dip[0](5,5) = _pmb[2] * _dip[0](5,1) + _fak *  _dip[0](1,1) ; // dxz-x-dxz, ok
             _dip[0](5,6) = _pmb[2] * _dip[0](5,2) + _fak *  _dip[0](1,2); // dxz-x-dyz, ok
             _dip[0](5,7) = _pmb[0] * _dip[0](5,1) + _fak * ( _dip[0](3,1) + _dip[0](5,0) + _ol(5,1) ); // dxz-x-dxx, ok
             _dip[0](5,8) = _pmb[1] * _dip[0](5,2) + _fak * _dip[0](5,0); // dxz-x-dyy, ok
             _dip[0](5,9) = _pmb[2] * _dip[0](5,3) + _fak * ( _dip[0](1,3) + _dip[0](5,0) ) ; // dxz-x-dzz, ok
             
             _dip[0](6,4) = _pmb[1] * _dip[0](6,1) + _fak * _dip[0](3,1); // dyz-x-dxy, ok
             _dip[0](6,5) = _pmb[2] * _dip[0](6,1) + _fak * _dip[0](2,1) ; // dyz-x-dxz, ok
             _dip[0](6,6) = _pmb[2] * _dip[0](6,2) + _fak * _dip[0](2,2); // dyz-x-dyz, ok
             _dip[0](6,7) = _pmb[0] * _dip[0](6,1) + _fak * ( _dip[0](6,0) + _ol(6,1) ); // dyz-x-dxx, ok
             _dip[0](6,8) = _pmb[1] * _dip[0](6,2) + _fak * ( _dip[0](3,2) + _dip[0](6,0) ); // dyz-x-dyy, ok
             _dip[0](6,9) = _pmb[2] * _dip[0](6,3) + _fak * ( _dip[0](2,3) + _dip[0](6,0) ) ; // dyz-x-dzz, ok
             
             _dip[0](7,4) = _pmb[1] * _dip[0](7,1); // dxx-x-dxy, ok
             _dip[0](7,5) = _pmb[2] * _dip[0](7,1); // dxx-x-dxz, ok
             _dip[0](7,6) = _pmb[2] * _dip[0](7,2); // dxx-x-dyz, ok
             _dip[0](7,7) = _pmb[0] * _dip[0](7,1) + _fak * ( 2.0*_dip[0](1,1) + _dip[0](7,0) +  _ol(7,1) ) ; // dxx-x-dxx, ok
             _dip[0](7,8) = _pmb[1] * _dip[0](7,2) + _fak *  _dip[0](7,0) ; // dxx-x-dyy, ok
             _dip[0](7,9) = _pmb[2] * _dip[0](7,3) + _fak *  _dip[0](7,0) ; // dxx-x-dzz, ok

             _dip[0](8,4) = _pmb[1] * _dip[0](8,1) + _fak * 2.0 * _dip[0](2,1) ; // dyy-x-dxy, ok
             _dip[0](8,5) = _pmb[2] * _dip[0](8,1); // dyy-x-dxz, ok
             _dip[0](8,6) = _pmb[2] * _dip[0](8,2); // dyy-x-dyz, ok
             _dip[0](8,7) = _pmb[0] * _dip[0](8,1) + _fak *  ( _dip[0](8,0) + _ol(8,1) ) ; // dyy-x-dxx, ok
             _dip[0](8,8) = _pmb[1] * _dip[0](8,2) + _fak *  (2.0 * _dip[0](2,2) + _dip[0](8,0) ) ; // dyy-x-dyy, ok
             _dip[0](8,9) = _pmb[2] * _dip[0](8,3) + _fak * _dip[0](8,0) ; // dyy-x-dzz, ok

             _dip[0](9,4) = _pmb[1] * _dip[0](9,1); // dzz-x-dxy, ok
             _dip[0](9,5) = _pmb[2] * _dip[0](9,1) + _fak * 2.0 * _dip[0](3,1); // dzz-x-dxz, ok
             _dip[0](9,6) = _pmb[2] * _dip[0](9,2) + _fak * 2.0 * _dip[0](3,2); // dzz-x-dyz, ok
             _dip[0](9,7) = _pmb[0] * _dip[0](9,1) + _fak * ( _dip[0](9,0) + _ol(9,1) ) ; // dzz-x-dxx, ok
             _dip[0](9,8) = _pmb[1] * _dip[0](9,2) + _fak * _dip[0](9,0) ; // dzz-x-dyy, ok
             _dip[0](9,9) = _pmb[2] * _dip[0](9,3) + _fak * (2.0 * _dip[0](3,3) + _dip[0](9,0)); // dzz-x-dzz

             
             // y-components
             _dip[1](4,4) = _pmb[1] * _dip[1](4,1) + _fak * ( _dip[1](1,1) + _ol(4,1) ) ; // dxy-y-dxy, ok
             _dip[1](4,5) = _pmb[2] * _dip[1](4,1); // dxy-y-dxz, ok
             _dip[1](4,6) = _pmb[2] * _dip[1](4,2); // dxy-y-dyz, ok
             _dip[1](4,7) = _pmb[0] * _dip[1](4,1) + _fak * ( _dip[1](2,1) + _dip[1](4,0) ); // dxy-y-dxx, ok 
             _dip[1](4,8) = _pmb[1] * _dip[1](4,2) + _fak * ( _dip[1](1,2) + _dip[1](4,0) + _ol(4,2) ); // dxy-y-dyy, ok
             _dip[1](4,9) = _pmb[2] * _dip[1](4,3) + _fak * _dip[1](4,0) ; // dxy-y-dzz, ok             
          
             
             _dip[1](5,4) = _pmb[1] * _dip[1](5,1) + _fak * _ol(5,1); // dxz-y-dxy, ok
             _dip[1](5,5) = _pmb[2] * _dip[1](5,1) + _fak * _dip[1](1,1) ; // dxz-y-dxz, ok
             _dip[1](5,6) = _pmb[2] * _dip[1](5,2) + _fak * _dip[1](1,2); // dxz-y-dyz, ok
             _dip[1](5,7) = _pmb[0] * _dip[1](5,1) + _fak * ( _dip[1](3,1) + _dip[1](5,0) ); // dxz-y-dxx, ok
             _dip[1](5,8) = _pmb[1] * _dip[1](5,2) + _fak * ( _dip[1](5,0) + _ol(5,2) ); // dxz-y-dyy, ok
             _dip[1](5,9) = _pmb[2] * _dip[1](5,3) + _fak * ( _dip[1](1,3) + _dip[1](5,0) ) ; // dxz-y-dzz, ok             
             
             
             _dip[1](6,4) = _pmb[1] * _dip[1](6,1) + _fak * ( _dip[1](3,1) + _ol(6,1) ); // dyz-y-dxy, ok
             _dip[1](6,5) = _pmb[2] * _dip[1](6,1) + _fak * _dip[1](2,1) ; // dyz-y-dxz, ok
             _dip[1](6,6) = _pmb[2] * _dip[1](6,2) + _fak * _dip[1](2,2); // dyz-y-dyz, ok
             _dip[1](6,7) = _pmb[0] * _dip[1](6,1) + _fak * _dip[1](6,0) ; // dyz-y-dxx, ok
             _dip[1](6,8) = _pmb[1] * _dip[1](6,2) + _fak * ( _dip[1](3,2) + _dip[1](6,0) + _ol(6,2) ); // dyz-y-dyy, ok
             _dip[1](6,9) = _pmb[2] * _dip[1](6,3) + _fak * ( _dip[1](2,3) + _dip[1](6,0) ) ; // dyz-y-dzz, ok           
             
             
             _dip[1](7,4) = _pmb[1] * _dip[1](7,1) + _fak * _ol(7,1); // dxx-y-dxy, ok
             _dip[1](7,5) = _pmb[2] * _dip[1](7,1); // dxx-y-dxz, ok
             _dip[1](7,6) = _pmb[2] * _dip[1](7,2); // dxx-y-dyz, ok
             _dip[1](7,7) = _pmb[0] * _dip[1](7,1) + _fak * ( 2.0*_dip[1](1,1) + _dip[1](7,0) ) ; // dxx-y-dxx, ok
             _dip[1](7,8) = _pmb[1] * _dip[1](7,2) + _fak * ( _dip[1](7,0) + _ol(7,2) ); // dxx-y-dyy, ok
             _dip[1](7,9) = _pmb[2] * _dip[1](7,3) + _fak *  _dip[1](7,0) ; // dxx-y-dzz, ok 
             
             
             _dip[1](8,4) = _pmb[1] * _dip[1](8,1) + _fak * ( 2.0 * _dip[1](2,1) + _ol(8,1) ) ; // dyy-y-dxy, ok
             _dip[1](8,5) = _pmb[2] * _dip[1](8,1); // dyy-y-dxz, ok
             _dip[1](8,6) = _pmb[2] * _dip[1](8,2); // dyy-y-dyz, ok
             _dip[1](8,7) = _pmb[0] * _dip[1](8,1) + _fak * _dip[1](8,0) ; // dyy-y-dxx, ok
             _dip[1](8,8) = _pmb[1] * _dip[1](8,2) + _fak *  (2.0 * _dip[1](2,2) + _dip[1](8,0) + _ol(8,2) ) ; // dyy-y-dyy, ok
             _dip[1](8,9) = _pmb[2] * _dip[1](8,3) + _fak * _dip[1](8,0) ; // dyy-y-dzz, ok

             _dip[1](9,4) = _pmb[1] * _dip[1](9,1) + _fak * _ol(9,1); // dzz-y-dxy, ok
             _dip[1](9,5) = _pmb[2] * _dip[1](9,1) + _fak * 2.0 * _dip[1](3,1); // dzz-y-dxz, ok
             _dip[1](9,6) = _pmb[2] * _dip[1](9,2) + _fak * 2.0 * _dip[1](3,2); // dzz-y-dyz, ok
             _dip[1](9,7) = _pmb[0] * _dip[1](9,1) + _fak * _dip[1](9,0) ; // dzz-y-dxx, ok
             _dip[1](9,8) = _pmb[1] * _dip[1](9,2) + _fak * ( _dip[1](9,0) + _ol(9,2) ) ; // dzz-y-dyy, ok
             _dip[1](9,9) = _pmb[2] * _dip[1](9,3) + _fak * (2.0 * _dip[1](3,3) + _dip[1](9,0)); // dzz-y-dzz             

             
             // z-components
             _dip[2](4,4) = _pmb[1] * _dip[2](4,1) + _fak * _dip[2](1,1) ; // dxy-z-dxy, ok
             _dip[2](4,5) = _pmb[2] * _dip[2](4,1) + _fak * _ol(4,1); // dxy-z-dxz, ok
             _dip[2](4,6) = _pmb[2] * _dip[2](4,2) + _fak * _ol(4,2); // dxy-z-dyz, ok
             _dip[2](4,7) = _pmb[0] * _dip[2](4,1) + _fak * ( _dip[2](2,1) + _dip[2](4,0) ); // dxy-z-dxx, ok 
             _dip[2](4,8) = _pmb[1] * _dip[2](4,2) + _fak * ( _dip[2](1,2) + _dip[2](4,0) ); // dxy-z-dyy, ok
             _dip[2](4,9) = _pmb[2] * _dip[2](4,3) + _fak * ( _dip[2](4,0) + _ol(4,3) ); // dxy-z-dzz, ok   

             _dip[2](5,4) = _pmb[1] * _dip[2](5,1); // dxz-z-dxy, ok
             _dip[2](5,5) = _pmb[2] * _dip[2](5,1) + _fak * ( _dip[2](1,1) + _ol(5,1) ) ; // dxz-z-dxz, ok
             _dip[2](5,6) = _pmb[2] * _dip[2](5,2) + _fak * ( _dip[2](1,2) + _ol(5,2) ) ; // dxz-z-dyz, ok
             _dip[2](5,7) = _pmb[0] * _dip[2](5,1) + _fak * ( _dip[2](3,1) + _dip[2](5,0) ); // dxz-z-dxx, ok
             _dip[2](5,8) = _pmb[1] * _dip[2](5,2) + _fak * _dip[2](5,0) ; // dxz-z-dyy, ok
             _dip[2](5,9) = _pmb[2] * _dip[2](5,3) + _fak * ( _dip[2](1,3) + _dip[2](5,0) + _ol(5,3) ) ; // dxz-z-dzz, ok                


             _dip[2](6,4) = _pmb[1] * _dip[2](6,1) + _fak * _dip[2](3,1) ; // dyz-z-dxy, ok
             _dip[2](6,5) = _pmb[2] * _dip[2](6,1) + _fak * ( _dip[2](2,1) + _ol(6,1) ); // dyz-z-dxz, ok
             _dip[2](6,6) = _pmb[2] * _dip[2](6,2) + _fak * ( _dip[2](2,2) + _ol(6,2) ); // dyz-z-dyz, ok
             _dip[2](6,7) = _pmb[0] * _dip[2](6,1) + _fak * _dip[2](6,0) ; // dyz-z-dxx, ok
             _dip[2](6,8) = _pmb[1] * _dip[2](6,2) + _fak * ( _dip[2](3,2) + _dip[2](6,0) ); // dyz-z-dyy, ok
             _dip[2](6,9) = _pmb[2] * _dip[2](6,3) + _fak * ( _dip[2](2,3) + _dip[2](6,0) + _ol(6,3) ) ; // dyz-z-dzz, ok                   


                          
             _dip[2](7,4) = _pmb[1] * _dip[2](7,1); // dxx-z-dxy, ok
             _dip[2](7,5) = _pmb[2] * _dip[2](7,1) + _fak * _ol(7,1); // dxx-z-dxz, ok
             _dip[2](7,6) = _pmb[2] * _dip[2](7,2) + _fak * _ol(7,2); // dxx-z-dyz, ok
             _dip[2](7,7) = _pmb[0] * _dip[2](7,1) + _fak * ( 2.0*_dip[2](1,1) + _dip[2](7,0) ) ; // dxx-z-dxx, ok
             _dip[2](7,8) = _pmb[1] * _dip[2](7,2) + _fak * _dip[2](7,0); // dxx-z-dyy, ok
             _dip[2](7,9) = _pmb[2] * _dip[2](7,3) + _fak * ( _dip[2](7,0) + _ol(7,3) ) ; // dxx-z-dzz, ok 
             
             _dip[2](8,4) = _pmb[1] * _dip[2](8,1) + _fak * 2.0 * _dip[2](2,1); // dyy-z-dxy, ok
             _dip[2](8,5) = _pmb[2] * _dip[2](8,1) + _fak * _ol(8,1); // dyy-z-dxz, ok
             _dip[2](8,6) = _pmb[2] * _dip[2](8,2) + _fak * _ol(8,2); // dyy-z-dyz, ok
             _dip[2](8,7) = _pmb[0] * _dip[2](8,1) + _fak * _dip[2](8,0) ; // dyy-z-dxx, ok
             _dip[2](8,8) = _pmb[1] * _dip[2](8,2) + _fak *  (2.0 * _dip[2](2,2) + _dip[2](8,0) ) ; // dyy-z-dyy, ok
             _dip[2](8,9) = _pmb[2] * _dip[2](8,3) + _fak * ( _dip[2](8,0) + _ol(8,3) ) ; // dyy-z-dzz, ok
             
             _dip[2](9,4) = _pmb[1] * _dip[2](9,1); // dzz-z-dxy, ok
             _dip[2](9,5) = _pmb[2] * _dip[2](9,1) + _fak * ( 2.0 * _dip[2](3,1) + _ol(9,1) ); // dzz-z-dxz, ok
             _dip[2](9,6) = _pmb[2] * _dip[2](9,2) + _fak * ( 2.0 * _dip[2](3,2) + _ol(9,2) ); // dzz-z-dyz, ok
             _dip[2](9,7) = _pmb[0] * _dip[2](9,1) + _fak * _dip[2](9,0) ; // dzz-z-dxx, ok
             _dip[2](9,8) = _pmb[1] * _dip[2](9,2) + _fak * _dip[2](9,0)  ; // dzz-z-dyy, ok
             _dip[2](9,9) = _pmb[2] * _dip[2](9,3) + _fak * (2.0 * _dip[2](3,3) + _dip[2](9,0) + _ol(9,3) ); // dzz-y-dzz               

             
        }
/*
        // s-f 
        if ( _lmax_col  > 2 ){
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
        if ( _lmax_row > 0 && _lmax_col  > 2 ){
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
        if (_lmax_row > 2 && _lmax_col  > 0 ){
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
        if ( _lmax_row > 1 && _lmax_col  >2 ){
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
        if ( _lmax_row > 2 && _lmax_col  > 1 ){
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
        if ( _lmax_row > 2 && _lmax_col  > 2 ){
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
        */
        
        
        
   
        
        // normalization and cartesian -> spherical factors
        int _ntrafo_row = _shell_row->getNumFunc() + _shell_row->getOffset();
        int _ntrafo_col = _shell_col->getNumFunc() + _shell_col->getOffset();
        
        //cout << " _ntrafo_row " << _ntrafo_row << ":" << _shell_row->getType() << endl;
        //cout << " _ntrafo_col " << _ntrafo_col << ":" << _shell_col->getType() << endl;
        ub::matrix<double> _trafo_row = ub::zero_matrix<double>(_ntrafo_row,_nrows);
        ub::matrix<double> _trafo_col = ub::zero_matrix<double>(_ntrafo_col,_ncols);

        // get transformation matrices
           // get transformation matrices including contraction coefficients
        std::vector<double> _contractions_row = (*itr)->contraction;
        std::vector<double> _contractions_col = (*itc)->contraction;
        
        this->getTrafo( _trafo_row, _lmax_row, _decay_row, _contractions_row);
        this->getTrafo( _trafo_col, _lmax_col, _decay_col, _contractions_col);
        ub::matrix<double> _trafo_col_tposed = ub::trans( _trafo_col );

        // cartesian -> spherical
       
        for ( int _i_comp = 0; _i_comp < 3; _i_comp++){

            ub::matrix<double> _dip_tmp = ub::prod( _trafo_row, _dip[ _i_comp ] );

            ub::matrix<double> _dip_sph = ub::prod( _dip_tmp, _trafo_col_tposed );
            
            // save to _matrix
            for ( unsigned i = 0; i< _matrix[0].size1(); i++ ) {
                for (unsigned j = 0; j < _matrix[0].size2(); j++){
                    _matrix[ _i_comp ](i,j) += _dip_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
                }
            }
        }
        
        _ol.clear();
            }// _shell_col Gaussians
        }// _shell_row Gaussians
    }
    
  
        
    
    
}}

