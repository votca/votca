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

    
    void AOKinetic::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, AOBasis* ecp) {
        //const double pi = boost::math::constants::pi<double>();
            /*cout << "\nAO block: "<< endl;
        cout << "\t row: " << _shell_row->getType() << " at " << _shell_row->getPos() << endl;
        cout << "\t col: " << _shell_col->getType() << " at " << _shell_col->getPos() << endl;*/
       
        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();
        
        
        if (_lmax_col >2 || _lmax_row >2){
            cerr << "Orbitals higher than d are not yet implemented. This should not have happened!" << flush;
             exit(1);
        }

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
        
             // get shell positions
        const vec& _pos_row = _shell_row->getPos();
        const vec& _pos_col = _shell_col->getPos();
        const vec  _diff    = _pos_row - _pos_col;
       
          
        double _distsq = (_diff.getX()*_diff.getX()) + (_diff.getY()*_diff.getY()) + (_diff.getZ()*_diff.getZ());   
        
        vector<double> PmA (3,0.0);
        vector<double> PmB (3,0.0);

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
            const double rzeta = 2.0 * _fak;
            const double rzetaA = 1/_decay_row;
            const double rzetaB = 1/_decay_col;
            const double xi=rzeta* _decay_row * _decay_col;
            
            
            // check if distance between postions is big, then skip step   
            double _exparg = xi *_distsq;
	    if ( _exparg > 30.0 ) { continue; }
            
    
            
            PmA[0] = rzeta*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_row.getX();
            PmA[1] = rzeta*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_row.getY();
            PmA[2] = rzeta*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_row.getZ();

            PmB[0] = rzeta*( _decay_row * _pos_row.getX() + _decay_col * _pos_col.getX() ) - _pos_col.getX();
            PmB[1] = rzeta*( _decay_row * _pos_row.getY() + _decay_col * _pos_col.getY() ) - _pos_col.getY();
            PmB[2] = rzeta*( _decay_row * _pos_row.getZ() + _decay_col * _pos_col.getZ() ) - _pos_col.getZ();
        
            // matrix for kinetic energies
            ub::matrix<double> kin = ub::zero_matrix<double>(_nrows,_ncols);
            //matrix for unnormalized overlap integrals
            ub::matrix<double> ol = ub::zero_matrix<double>(_nrows,_ncols);
        
            // s-s overlap integral
            ol(Cart::s,Cart::s) = pow(rzeta,1.5)*pow(4.0*_decay_row*_decay_col,0.75) * exp(-_exparg);
            // s-s- kinetic energy integral
            kin(Cart::s,Cart::s)= ol(Cart::s,Cart::s)*xi*(3-2*xi*_distsq);
            
          //s-p
            if ( _lmax_col>0 ){
            //std::cout << "\t setting s-p|" << std::flush;      
            ol(Cart::s,Cart::y) = PmB[1]*ol(Cart::s,Cart::s);
            ol(Cart::s,Cart::x) = PmB[0]*ol(Cart::s,Cart::s);
            ol(Cart::s,Cart::z) = PmB[2]*ol(Cart::s,Cart::s);
            
            kin(Cart::s,Cart::y) = PmB[1]*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::s,Cart::y));
            kin(Cart::s,Cart::x) = PmB[0]*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::s,Cart::x));
            kin(Cart::s,Cart::z) = PmB[2]*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::s,Cart::z));
            }
            
            //p-s
            if ( _lmax_row>0 ){
            //std::cout << "\t setting p-s|" << std::flush;
            ol(Cart::y,Cart::s) = PmA[1]*ol(Cart::s,Cart::s);
            ol(Cart::x,Cart::s) = PmA[0]*ol(Cart::s,Cart::s);
            ol(Cart::z,Cart::s) = PmA[2]*ol(Cart::s,Cart::s);
       
            kin(Cart::y,Cart::s) = PmA[1]*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::y,Cart::s));
            kin(Cart::x,Cart::s) = PmA[0]*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::x,Cart::s));
            kin(Cart::z,Cart::s) = PmA[2]*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::z,Cart::s));
            }
            
            //p-p
            if ( _lmax_row>0 && _lmax_col>0 ){
            ol(Cart::y,Cart::y) = PmA[1]*ol(Cart::s,Cart::y)+1*_fak*ol(Cart::s,Cart::s);
            ol(Cart::y,Cart::x) = PmA[1]*ol(Cart::s,Cart::x);
            ol(Cart::y,Cart::z) = PmA[1]*ol(Cart::s,Cart::z);
            ol(Cart::x,Cart::y) = PmA[0]*ol(Cart::s,Cart::y);
            ol(Cart::x,Cart::x) = PmA[0]*ol(Cart::s,Cart::x)+1*_fak*ol(Cart::s,Cart::s);
            ol(Cart::x,Cart::z) = PmA[0]*ol(Cart::s,Cart::z);
            ol(Cart::z,Cart::y) = PmA[2]*ol(Cart::s,Cart::y);
            ol(Cart::z,Cart::x) = PmA[2]*ol(Cart::s,Cart::x);
            ol(Cart::z,Cart::z) = PmA[2]*ol(Cart::s,Cart::z)+1*_fak*ol(Cart::s,Cart::s);
     
            //std::cout << "\t setting p-p|" << std::flush;
            kin(Cart::y,Cart::y) = PmA[1]*kin(Cart::s,Cart::y)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::y,Cart::y));
            kin(Cart::y,Cart::x) = PmA[1]*kin(Cart::s,Cart::x)+2*xi*(ol(Cart::y,Cart::x));
            kin(Cart::y,Cart::z) = PmA[1]*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::y,Cart::z));
            kin(Cart::x,Cart::y) = PmA[0]*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::x,Cart::y));
            kin(Cart::x,Cart::x) = PmA[0]*kin(Cart::s,Cart::x)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::x,Cart::x));
            kin(Cart::x,Cart::z) = PmA[0]*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::x,Cart::z));
            kin(Cart::z,Cart::y) = PmA[2]*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::z,Cart::y));
            kin(Cart::z,Cart::x) = PmA[2]*kin(Cart::s,Cart::x)+2*xi*(ol(Cart::z,Cart::x));
            kin(Cart::z,Cart::z) = PmA[2]*kin(Cart::s,Cart::z)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::z,Cart::z));
            }
            
            //s-d
            if ( _lmax_col>1 ){
            ol(Cart::s,Cart::yy) = PmB[1]*ol(Cart::s,Cart::y)+1*_fak*ol(Cart::s,Cart::s);
            ol(Cart::s,Cart::xy) = PmB[0]*ol(Cart::s,Cart::y);
            ol(Cart::s,Cart::yz) = PmB[1]*ol(Cart::s,Cart::z);
            ol(Cart::s,Cart::xx) = PmB[0]*ol(Cart::s,Cart::x)+1*_fak*ol(Cart::s,Cart::s);
            ol(Cart::s,Cart::xz) = PmB[0]*ol(Cart::s,Cart::z);
            ol(Cart::s,Cart::zz) = PmB[2]*ol(Cart::s,Cart::z)+1*_fak*ol(Cart::s,Cart::s);    
                
            //std::cout << "\t setting s-d|" << std::flush;
            kin(Cart::s,Cart::yy) = PmB[1]*kin(Cart::s,Cart::y)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::s,Cart::yy)-0.5*rzetaB*ol(Cart::s,Cart::s));
            kin(Cart::s,Cart::xy) = PmB[0]*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::s,Cart::xy));
            kin(Cart::s,Cart::yz) = PmB[1]*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::s,Cart::yz));
            kin(Cart::s,Cart::xx) = PmB[0]*kin(Cart::s,Cart::x)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::s,Cart::xx)-0.5*rzetaB*ol(Cart::s,Cart::s));
            kin(Cart::s,Cart::xz) = PmB[0]*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::s,Cart::xz));
            kin(Cart::s,Cart::zz) = PmB[2]*kin(Cart::s,Cart::z)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::s,Cart::zz)-0.5*rzetaB*ol(Cart::s,Cart::s));            
            }

            //d-s
            if ( _lmax_row>1 ){
            //std::cout << "\t setting d-s|" << std::flush;
            ol(Cart::yy,Cart::s) = PmA[1]*ol(Cart::y,Cart::s)+1*_fak*ol(Cart::s,Cart::s);
            ol(Cart::xy,Cart::s) = PmA[0]*ol(Cart::y,Cart::s);
            ol(Cart::yz,Cart::s) = PmA[1]*ol(Cart::z,Cart::s);
            ol(Cart::xx,Cart::s) = PmA[0]*ol(Cart::x,Cart::s)+1*_fak*ol(Cart::s,Cart::s);
            ol(Cart::xz,Cart::s) = PmA[0]*ol(Cart::z,Cart::s);
            ol(Cart::zz,Cart::s) = PmA[2]*ol(Cart::z,Cart::s)+1*_fak*ol(Cart::s,Cart::s);
            
            kin(Cart::yy,Cart::s) = PmA[1]*kin(Cart::y,Cart::s)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::yy,Cart::s)-0.5*rzetaA*ol(Cart::s,Cart::s));
            kin(Cart::xy,Cart::s) = PmA[0]*kin(Cart::y,Cart::s)+2*xi*(ol(Cart::xy,Cart::s));
            kin(Cart::yz,Cart::s) = PmA[1]*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::yz,Cart::s));
            kin(Cart::xx,Cart::s) = PmA[0]*kin(Cart::x,Cart::s)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::xx,Cart::s)-0.5*rzetaA*ol(Cart::s,Cart::s));
            kin(Cart::xz,Cart::s) = PmA[0]*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::xz,Cart::s));
            kin(Cart::zz,Cart::s) = PmA[2]*kin(Cart::z,Cart::s)+1*_fak*kin(Cart::s,Cart::s)+2*xi*(ol(Cart::zz,Cart::s)-0.5*rzetaA*ol(Cart::s,Cart::s));
            }
            
            //p-d
            if ( _lmax_row>0 && _lmax_col>1 ){
            //std::cout << "\t setting p-d|" << std::flush;
            ol(Cart::y,Cart::yy) = PmB[1]*ol(Cart::y,Cart::y)+1*_fak*ol(Cart::y,Cart::s)+1*_fak*ol(Cart::s,Cart::y);
            ol(Cart::y,Cart::xy) = PmB[0]*ol(Cart::y,Cart::y);
            ol(Cart::y,Cart::yz) = PmB[1]*ol(Cart::y,Cart::z)+1*_fak*ol(Cart::s,Cart::z);
            ol(Cart::y,Cart::xx) = PmB[0]*ol(Cart::y,Cart::x)+1*_fak*ol(Cart::y,Cart::s);
            ol(Cart::y,Cart::xz) = PmB[0]*ol(Cart::y,Cart::z);
            ol(Cart::y,Cart::zz) = PmB[2]*ol(Cart::y,Cart::z)+1*_fak*ol(Cart::y,Cart::s);
            ol(Cart::x,Cart::yy) = PmB[1]*ol(Cart::x,Cart::y)+1*_fak*ol(Cart::x,Cart::s);
            ol(Cart::x,Cart::xy) = PmB[0]*ol(Cart::x,Cart::y)+1*_fak*ol(Cart::s,Cart::y);
            ol(Cart::x,Cart::yz) = PmB[1]*ol(Cart::x,Cart::z);
            ol(Cart::x,Cart::xx) = PmB[0]*ol(Cart::x,Cart::x)+1*_fak*ol(Cart::x,Cart::s)+1*_fak*ol(Cart::s,Cart::x);
            ol(Cart::x,Cart::xz) = PmB[0]*ol(Cart::x,Cart::z)+1*_fak*ol(Cart::s,Cart::z);
            ol(Cart::x,Cart::zz) = PmB[2]*ol(Cart::x,Cart::z)+1*_fak*ol(Cart::x,Cart::s);
            ol(Cart::z,Cart::yy) = PmB[1]*ol(Cart::z,Cart::y)+1*_fak*ol(Cart::z,Cart::s);
            ol(Cart::z,Cart::xy) = PmB[0]*ol(Cart::z,Cart::y);
            ol(Cart::z,Cart::yz) = PmB[1]*ol(Cart::z,Cart::z);
            ol(Cart::z,Cart::xx) = PmB[0]*ol(Cart::z,Cart::x)+1*_fak*ol(Cart::z,Cart::s);
            ol(Cart::z,Cart::xz) = PmB[0]*ol(Cart::z,Cart::z);
            ol(Cart::z,Cart::zz) = PmB[2]*ol(Cart::z,Cart::z)+1*_fak*ol(Cart::z,Cart::s)+1*_fak*ol(Cart::s,Cart::z);
            
            kin(Cart::y,Cart::yy) = PmB[1]*kin(Cart::y,Cart::y)+1*_fak*kin(Cart::y,Cart::s)+1*_fak*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::y,Cart::yy)-0.5*rzetaB*ol(Cart::y,Cart::s));
            kin(Cart::y,Cart::xy) = PmB[0]*kin(Cart::y,Cart::y)+2*xi*(ol(Cart::y,Cart::xy));
            kin(Cart::y,Cart::yz) = PmB[1]*kin(Cart::y,Cart::z)+1*_fak*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::y,Cart::yz));
            kin(Cart::y,Cart::xx) = PmB[0]*kin(Cart::y,Cart::x)+1*_fak*kin(Cart::y,Cart::s)+2*xi*(ol(Cart::y,Cart::xx)-0.5*rzetaB*ol(Cart::y,Cart::s));
            kin(Cart::y,Cart::xz) = PmB[0]*kin(Cart::y,Cart::z)+2*xi*(ol(Cart::y,Cart::xz));
            kin(Cart::y,Cart::zz) = PmB[2]*kin(Cart::y,Cart::z)+1*_fak*kin(Cart::y,Cart::s)+2*xi*(ol(Cart::y,Cart::zz)-0.5*rzetaB*ol(Cart::y,Cart::s));
            kin(Cart::x,Cart::yy) = PmB[1]*kin(Cart::x,Cart::y)+1*_fak*kin(Cart::x,Cart::s)+2*xi*(ol(Cart::x,Cart::yy)-0.5*rzetaB*ol(Cart::x,Cart::s));
            kin(Cart::x,Cart::xy) = PmB[0]*kin(Cart::x,Cart::y)+1*_fak*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::x,Cart::xy));
            kin(Cart::x,Cart::yz) = PmB[1]*kin(Cart::x,Cart::z)+2*xi*(ol(Cart::x,Cart::yz));
            kin(Cart::x,Cart::xx) = PmB[0]*kin(Cart::x,Cart::x)+1*_fak*kin(Cart::x,Cart::s)+1*_fak*kin(Cart::s,Cart::x)+2*xi*(ol(Cart::x,Cart::xx)-0.5*rzetaB*ol(Cart::x,Cart::s));
            kin(Cart::x,Cart::xz) = PmB[0]*kin(Cart::x,Cart::z)+1*_fak*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::x,Cart::xz));
            kin(Cart::x,Cart::zz) = PmB[2]*kin(Cart::x,Cart::z)+1*_fak*kin(Cart::x,Cart::s)+2*xi*(ol(Cart::x,Cart::zz)-0.5*rzetaB*ol(Cart::x,Cart::s));
            kin(Cart::z,Cart::yy) = PmB[1]*kin(Cart::z,Cart::y)+1*_fak*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::z,Cart::yy)-0.5*rzetaB*ol(Cart::z,Cart::s));
            kin(Cart::z,Cart::xy) = PmB[0]*kin(Cart::z,Cart::y)+2*xi*(ol(Cart::z,Cart::xy));
            kin(Cart::z,Cart::yz) = PmB[1]*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::z,Cart::yz));
            kin(Cart::z,Cart::xx) = PmB[0]*kin(Cart::z,Cart::x)+1*_fak*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::z,Cart::xx)-0.5*rzetaB*ol(Cart::z,Cart::s));
            kin(Cart::z,Cart::xz) = PmB[0]*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::z,Cart::xz));
            kin(Cart::z,Cart::zz) = PmB[2]*kin(Cart::z,Cart::z)+1*_fak*kin(Cart::z,Cart::s)+1*_fak*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::z,Cart::zz)-0.5*rzetaB*ol(Cart::z,Cart::s));
            }
            
            //d-p
            if ( _lmax_row>1 && _lmax_col>0 ){
            //std::cout << "\t setting d-p|" << std::flush;
            ol(Cart::yy,Cart::y) = PmA[1]*ol(Cart::y,Cart::y)+1*_fak*ol(Cart::s,Cart::y)+1*_fak*ol(Cart::y,Cart::s);
            ol(Cart::yy,Cart::x) = PmA[1]*ol(Cart::y,Cart::x)+1*_fak*ol(Cart::s,Cart::x);
            ol(Cart::yy,Cart::z) = PmA[1]*ol(Cart::y,Cart::z)+1*_fak*ol(Cart::s,Cart::z);
            ol(Cart::xy,Cart::y) = PmA[0]*ol(Cart::y,Cart::y);
            ol(Cart::xy,Cart::x) = PmA[0]*ol(Cart::y,Cart::x)+1*_fak*ol(Cart::y,Cart::s);
            ol(Cart::xy,Cart::z) = PmA[0]*ol(Cart::y,Cart::z);
            ol(Cart::yz,Cart::y) = PmA[1]*ol(Cart::z,Cart::y)+1*_fak*ol(Cart::z,Cart::s);
            ol(Cart::yz,Cart::x) = PmA[1]*ol(Cart::z,Cart::x);
            ol(Cart::yz,Cart::z) = PmA[1]*ol(Cart::z,Cart::z);
            ol(Cart::xx,Cart::y) = PmA[0]*ol(Cart::x,Cart::y)+1*_fak*ol(Cart::s,Cart::y);
            ol(Cart::xx,Cart::x) = PmA[0]*ol(Cart::x,Cart::x)+1*_fak*ol(Cart::s,Cart::x)+1*_fak*ol(Cart::x,Cart::s);
            ol(Cart::xx,Cart::z) = PmA[0]*ol(Cart::x,Cart::z)+1*_fak*ol(Cart::s,Cart::z);
            ol(Cart::xz,Cart::y) = PmA[0]*ol(Cart::z,Cart::y);
            ol(Cart::xz,Cart::x) = PmA[0]*ol(Cart::z,Cart::x)+1*_fak*ol(Cart::z,Cart::s);
            ol(Cart::xz,Cart::z) = PmA[0]*ol(Cart::z,Cart::z);
            ol(Cart::zz,Cart::y) = PmA[2]*ol(Cart::z,Cart::y)+1*_fak*ol(Cart::s,Cart::y);
            ol(Cart::zz,Cart::x) = PmA[2]*ol(Cart::z,Cart::x)+1*_fak*ol(Cart::s,Cart::x);
            ol(Cart::zz,Cart::z) = PmA[2]*ol(Cart::z,Cart::z)+1*_fak*ol(Cart::s,Cart::z)+1*_fak*ol(Cart::z,Cart::s);
            
            kin(Cart::yy,Cart::y) = PmA[1]*kin(Cart::y,Cart::y)+1*_fak*kin(Cart::s,Cart::y)+1*_fak*kin(Cart::y,Cart::s)+2*xi*(ol(Cart::yy,Cart::y)-0.5*rzetaA*ol(Cart::s,Cart::y));
            kin(Cart::yy,Cart::x) = PmA[1]*kin(Cart::y,Cart::x)+1*_fak*kin(Cart::s,Cart::x)+2*xi*(ol(Cart::yy,Cart::x)-0.5*rzetaA*ol(Cart::s,Cart::x));
            kin(Cart::yy,Cart::z) = PmA[1]*kin(Cart::y,Cart::z)+1*_fak*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::yy,Cart::z)-0.5*rzetaA*ol(Cart::s,Cart::z));
            kin(Cart::xy,Cart::y) = PmA[0]*kin(Cart::y,Cart::y)+2*xi*(ol(Cart::xy,Cart::y));
            kin(Cart::xy,Cart::x) = PmA[0]*kin(Cart::y,Cart::x)+1*_fak*kin(Cart::y,Cart::s)+2*xi*(ol(Cart::xy,Cart::x));
            kin(Cart::xy,Cart::z) = PmA[0]*kin(Cart::y,Cart::z)+2*xi*(ol(Cart::xy,Cart::z));
            kin(Cart::yz,Cart::y) = PmA[1]*kin(Cart::z,Cart::y)+1*_fak*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::yz,Cart::y));
            kin(Cart::yz,Cart::x) = PmA[1]*kin(Cart::z,Cart::x)+2*xi*(ol(Cart::yz,Cart::x));
            kin(Cart::yz,Cart::z) = PmA[1]*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::yz,Cart::z));
            kin(Cart::xx,Cart::y) = PmA[0]*kin(Cart::x,Cart::y)+1*_fak*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::xx,Cart::y)-0.5*rzetaA*ol(Cart::s,Cart::y));
            kin(Cart::xx,Cart::x) = PmA[0]*kin(Cart::x,Cart::x)+1*_fak*kin(Cart::s,Cart::x)+1*_fak*kin(Cart::x,Cart::s)+2*xi*(ol(Cart::xx,Cart::x)-0.5*rzetaA*ol(Cart::s,Cart::x));
            kin(Cart::xx,Cart::z) = PmA[0]*kin(Cart::x,Cart::z)+1*_fak*kin(Cart::s,Cart::z)+2*xi*(ol(Cart::xx,Cart::z)-0.5*rzetaA*ol(Cart::s,Cart::z));
            kin(Cart::xz,Cart::y) = PmA[0]*kin(Cart::z,Cart::y)+2*xi*(ol(Cart::xz,Cart::y));
            kin(Cart::xz,Cart::x) = PmA[0]*kin(Cart::z,Cart::x)+1*_fak*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::xz,Cart::x));
            kin(Cart::xz,Cart::z) = PmA[0]*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::xz,Cart::z));
            kin(Cart::zz,Cart::y) = PmA[2]*kin(Cart::z,Cart::y)+1*_fak*kin(Cart::s,Cart::y)+2*xi*(ol(Cart::zz,Cart::y)-0.5*rzetaA*ol(Cart::s,Cart::y));
            kin(Cart::zz,Cart::x) = PmA[2]*kin(Cart::z,Cart::x)+1*_fak*kin(Cart::s,Cart::x)+2*xi*(ol(Cart::zz,Cart::x)-0.5*rzetaA*ol(Cart::s,Cart::x));
            kin(Cart::zz,Cart::z) = PmA[2]*kin(Cart::z,Cart::z)+1*_fak*kin(Cart::s,Cart::z)+1*_fak*kin(Cart::z,Cart::s)+2*xi*(ol(Cart::zz,Cart::z)-0.5*rzetaA*ol(Cart::s,Cart::z));
            }

            //d-d
            if ( _lmax_row>1 && _lmax_col>1 ){
            //std::cout << "\t setting d-d|" << std::flush;
            ol(Cart::yy,Cart::yy) = PmA[1]*ol(Cart::y,Cart::yy)+1*_fak*ol(Cart::s,Cart::yy)+2*_fak*ol(Cart::y,Cart::y);
            ol(Cart::yy,Cart::xy) = PmA[1]*ol(Cart::y,Cart::xy)+1*_fak*ol(Cart::s,Cart::xy)+1*_fak*ol(Cart::y,Cart::x);
            ol(Cart::yy,Cart::yz) = PmA[1]*ol(Cart::y,Cart::yz)+1*_fak*ol(Cart::s,Cart::yz)+1*_fak*ol(Cart::y,Cart::z);
            ol(Cart::yy,Cart::xx) = PmA[1]*ol(Cart::y,Cart::xx)+1*_fak*ol(Cart::s,Cart::xx);
            ol(Cart::yy,Cart::xz) = PmA[1]*ol(Cart::y,Cart::xz)+1*_fak*ol(Cart::s,Cart::xz);
            ol(Cart::yy,Cart::zz) = PmA[1]*ol(Cart::y,Cart::zz)+1*_fak*ol(Cart::s,Cart::zz);
            ol(Cart::xy,Cart::yy) = PmA[0]*ol(Cart::y,Cart::yy);
            ol(Cart::xy,Cart::xy) = PmA[0]*ol(Cart::y,Cart::xy)+1*_fak*ol(Cart::y,Cart::y);
            ol(Cart::xy,Cart::yz) = PmA[0]*ol(Cart::y,Cart::yz);
            ol(Cart::xy,Cart::xx) = PmA[0]*ol(Cart::y,Cart::xx)+2*_fak*ol(Cart::y,Cart::x);
            ol(Cart::xy,Cart::xz) = PmA[0]*ol(Cart::y,Cart::xz)+1*_fak*ol(Cart::y,Cart::z);
            ol(Cart::xy,Cart::zz) = PmA[0]*ol(Cart::y,Cart::zz);
            ol(Cart::yz,Cart::yy) = PmA[1]*ol(Cart::z,Cart::yy)+2*_fak*ol(Cart::z,Cart::y);
            ol(Cart::yz,Cart::xy) = PmA[1]*ol(Cart::z,Cart::xy)+1*_fak*ol(Cart::z,Cart::x);
            ol(Cart::yz,Cart::yz) = PmA[1]*ol(Cart::z,Cart::yz)+1*_fak*ol(Cart::z,Cart::z);
            ol(Cart::yz,Cart::xx) = PmA[1]*ol(Cart::z,Cart::xx);
            ol(Cart::yz,Cart::xz) = PmA[1]*ol(Cart::z,Cart::xz);
            ol(Cart::yz,Cart::zz) = PmA[1]*ol(Cart::z,Cart::zz);
            ol(Cart::xx,Cart::yy) = PmA[0]*ol(Cart::x,Cart::yy)+1*_fak*ol(Cart::s,Cart::yy);
            ol(Cart::xx,Cart::xy) = PmA[0]*ol(Cart::x,Cart::xy)+1*_fak*ol(Cart::s,Cart::xy)+1*_fak*ol(Cart::x,Cart::y);
            ol(Cart::xx,Cart::yz) = PmA[0]*ol(Cart::x,Cart::yz)+1*_fak*ol(Cart::s,Cart::yz);
            ol(Cart::xx,Cart::xx) = PmA[0]*ol(Cart::x,Cart::xx)+1*_fak*ol(Cart::s,Cart::xx)+2*_fak*ol(Cart::x,Cart::x);
            ol(Cart::xx,Cart::xz) = PmA[0]*ol(Cart::x,Cart::xz)+1*_fak*ol(Cart::s,Cart::xz)+1*_fak*ol(Cart::x,Cart::z);
            ol(Cart::xx,Cart::zz) = PmA[0]*ol(Cart::x,Cart::zz)+1*_fak*ol(Cart::s,Cart::zz);
            ol(Cart::xz,Cart::yy) = PmA[0]*ol(Cart::z,Cart::yy);
            ol(Cart::xz,Cart::xy) = PmA[0]*ol(Cart::z,Cart::xy)+1*_fak*ol(Cart::z,Cart::y);
            ol(Cart::xz,Cart::yz) = PmA[0]*ol(Cart::z,Cart::yz);
            ol(Cart::xz,Cart::xx) = PmA[0]*ol(Cart::z,Cart::xx)+2*_fak*ol(Cart::z,Cart::x);
            ol(Cart::xz,Cart::xz) = PmA[0]*ol(Cart::z,Cart::xz)+1*_fak*ol(Cart::z,Cart::z);
            ol(Cart::xz,Cart::zz) = PmA[0]*ol(Cart::z,Cart::zz);
            ol(Cart::zz,Cart::yy) = PmA[2]*ol(Cart::z,Cart::yy)+1*_fak*ol(Cart::s,Cart::yy);
            ol(Cart::zz,Cart::xy) = PmA[2]*ol(Cart::z,Cart::xy)+1*_fak*ol(Cart::s,Cart::xy);
            ol(Cart::zz,Cart::yz) = PmA[2]*ol(Cart::z,Cart::yz)+1*_fak*ol(Cart::s,Cart::yz)+1*_fak*ol(Cart::z,Cart::y);
            ol(Cart::zz,Cart::xx) = PmA[2]*ol(Cart::z,Cart::xx)+1*_fak*ol(Cart::s,Cart::xx);
            ol(Cart::zz,Cart::xz) = PmA[2]*ol(Cart::z,Cart::xz)+1*_fak*ol(Cart::s,Cart::xz)+1*_fak*ol(Cart::z,Cart::x);
            ol(Cart::zz,Cart::zz) = PmA[2]*ol(Cart::z,Cart::zz)+1*_fak*ol(Cart::s,Cart::zz)+2*_fak*ol(Cart::z,Cart::z);

            kin(Cart::yy,Cart::yy) = PmA[1]*kin(Cart::y,Cart::yy)+1*_fak*kin(Cart::s,Cart::yy)+2*_fak*kin(Cart::y,Cart::y)+2*xi*(ol(Cart::yy,Cart::yy)-0.5*rzetaA*ol(Cart::s,Cart::yy));
            kin(Cart::yy,Cart::xy) = PmA[1]*kin(Cart::y,Cart::xy)+1*_fak*kin(Cart::s,Cart::xy)+1*_fak*kin(Cart::y,Cart::x)+2*xi*(ol(Cart::yy,Cart::xy)-0.5*rzetaA*ol(Cart::s,Cart::xy));
            kin(Cart::yy,Cart::yz) = PmA[1]*kin(Cart::y,Cart::yz)+1*_fak*kin(Cart::s,Cart::yz)+1*_fak*kin(Cart::y,Cart::z)+2*xi*(ol(Cart::yy,Cart::yz)-0.5*rzetaA*ol(Cart::s,Cart::yz));
            kin(Cart::yy,Cart::xx) = PmA[1]*kin(Cart::y,Cart::xx)+1*_fak*kin(Cart::s,Cart::xx)+2*xi*(ol(Cart::yy,Cart::xx)-0.5*rzetaA*ol(Cart::s,Cart::xx));
            kin(Cart::yy,Cart::xz) = PmA[1]*kin(Cart::y,Cart::xz)+1*_fak*kin(Cart::s,Cart::xz)+2*xi*(ol(Cart::yy,Cart::xz)-0.5*rzetaA*ol(Cart::s,Cart::xz));
            kin(Cart::yy,Cart::zz) = PmA[1]*kin(Cart::y,Cart::zz)+1*_fak*kin(Cart::s,Cart::zz)+2*xi*(ol(Cart::yy,Cart::zz)-0.5*rzetaA*ol(Cart::s,Cart::zz));
            kin(Cart::xy,Cart::yy) = PmA[0]*kin(Cart::y,Cart::yy)+2*xi*(ol(Cart::xy,Cart::yy));
            kin(Cart::xy,Cart::xy) = PmA[0]*kin(Cart::y,Cart::xy)+1*_fak*kin(Cart::y,Cart::y)+2*xi*(ol(Cart::xy,Cart::xy));
            kin(Cart::xy,Cart::yz) = PmA[0]*kin(Cart::y,Cart::yz)+2*xi*(ol(Cart::xy,Cart::yz));
            kin(Cart::xy,Cart::xx) = PmA[0]*kin(Cart::y,Cart::xx)+2*_fak*kin(Cart::y,Cart::x)+2*xi*(ol(Cart::xy,Cart::xx));
            kin(Cart::xy,Cart::xz) = PmA[0]*kin(Cart::y,Cart::xz)+1*_fak*kin(Cart::y,Cart::z)+2*xi*(ol(Cart::xy,Cart::xz));
            kin(Cart::xy,Cart::zz) = PmA[0]*kin(Cart::y,Cart::zz)+2*xi*(ol(Cart::xy,Cart::zz));
            kin(Cart::yz,Cart::yy) = PmA[1]*kin(Cart::z,Cart::yy)+2*_fak*kin(Cart::z,Cart::y)+2*xi*(ol(Cart::yz,Cart::yy));
            kin(Cart::yz,Cart::xy) = PmA[1]*kin(Cart::z,Cart::xy)+1*_fak*kin(Cart::z,Cart::x)+2*xi*(ol(Cart::yz,Cart::xy));
            kin(Cart::yz,Cart::yz) = PmA[1]*kin(Cart::z,Cart::yz)+1*_fak*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::yz,Cart::yz));
            kin(Cart::yz,Cart::xx) = PmA[1]*kin(Cart::z,Cart::xx)+2*xi*(ol(Cart::yz,Cart::xx));
            kin(Cart::yz,Cart::xz) = PmA[1]*kin(Cart::z,Cart::xz)+2*xi*(ol(Cart::yz,Cart::xz));
            kin(Cart::yz,Cart::zz) = PmA[1]*kin(Cart::z,Cart::zz)+2*xi*(ol(Cart::yz,Cart::zz));
            kin(Cart::xx,Cart::yy) = PmA[0]*kin(Cart::x,Cart::yy)+1*_fak*kin(Cart::s,Cart::yy)+2*xi*(ol(Cart::xx,Cart::yy)-0.5*rzetaA*ol(Cart::s,Cart::yy));
            kin(Cart::xx,Cart::xy) = PmA[0]*kin(Cart::x,Cart::xy)+1*_fak*kin(Cart::s,Cart::xy)+1*_fak*kin(Cart::x,Cart::y)+2*xi*(ol(Cart::xx,Cart::xy)-0.5*rzetaA*ol(Cart::s,Cart::xy));
            kin(Cart::xx,Cart::yz) = PmA[0]*kin(Cart::x,Cart::yz)+1*_fak*kin(Cart::s,Cart::yz)+2*xi*(ol(Cart::xx,Cart::yz)-0.5*rzetaA*ol(Cart::s,Cart::yz));
            kin(Cart::xx,Cart::xx) = PmA[0]*kin(Cart::x,Cart::xx)+1*_fak*kin(Cart::s,Cart::xx)+2*_fak*kin(Cart::x,Cart::x)+2*xi*(ol(Cart::xx,Cart::xx)-0.5*rzetaA*ol(Cart::s,Cart::xx));
            kin(Cart::xx,Cart::xz) = PmA[0]*kin(Cart::x,Cart::xz)+1*_fak*kin(Cart::s,Cart::xz)+1*_fak*kin(Cart::x,Cart::z)+2*xi*(ol(Cart::xx,Cart::xz)-0.5*rzetaA*ol(Cart::s,Cart::xz));
            kin(Cart::xx,Cart::zz) = PmA[0]*kin(Cart::x,Cart::zz)+1*_fak*kin(Cart::s,Cart::zz)+2*xi*(ol(Cart::xx,Cart::zz)-0.5*rzetaA*ol(Cart::s,Cart::zz));
            kin(Cart::xz,Cart::yy) = PmA[0]*kin(Cart::z,Cart::yy)+2*xi*(ol(Cart::xz,Cart::yy));
            kin(Cart::xz,Cart::xy) = PmA[0]*kin(Cart::z,Cart::xy)+1*_fak*kin(Cart::z,Cart::y)+2*xi*(ol(Cart::xz,Cart::xy));
            kin(Cart::xz,Cart::yz) = PmA[0]*kin(Cart::z,Cart::yz)+2*xi*(ol(Cart::xz,Cart::yz));
            kin(Cart::xz,Cart::xx) = PmA[0]*kin(Cart::z,Cart::xx)+2*_fak*kin(Cart::z,Cart::x)+2*xi*(ol(Cart::xz,Cart::xx));
            kin(Cart::xz,Cart::xz) = PmA[0]*kin(Cart::z,Cart::xz)+1*_fak*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::xz,Cart::xz));
            kin(Cart::xz,Cart::zz) = PmA[0]*kin(Cart::z,Cart::zz)+2*xi*(ol(Cart::xz,Cart::zz));
            kin(Cart::zz,Cart::yy) = PmA[2]*kin(Cart::z,Cart::yy)+1*_fak*kin(Cart::s,Cart::yy)+2*xi*(ol(Cart::zz,Cart::yy)-0.5*rzetaA*ol(Cart::s,Cart::yy));
            kin(Cart::zz,Cart::xy) = PmA[2]*kin(Cart::z,Cart::xy)+1*_fak*kin(Cart::s,Cart::xy)+2*xi*(ol(Cart::zz,Cart::xy)-0.5*rzetaA*ol(Cart::s,Cart::xy));
            kin(Cart::zz,Cart::yz) = PmA[2]*kin(Cart::z,Cart::yz)+1*_fak*kin(Cart::s,Cart::yz)+1*_fak*kin(Cart::z,Cart::y)+2*xi*(ol(Cart::zz,Cart::yz)-0.5*rzetaA*ol(Cart::s,Cart::yz));
            kin(Cart::zz,Cart::xx) = PmA[2]*kin(Cart::z,Cart::xx)+1*_fak*kin(Cart::s,Cart::xx)+2*xi*(ol(Cart::zz,Cart::xx)-0.5*rzetaA*ol(Cart::s,Cart::xx));
            kin(Cart::zz,Cart::xz) = PmA[2]*kin(Cart::z,Cart::xz)+1*_fak*kin(Cart::s,Cart::xz)+1*_fak*kin(Cart::z,Cart::x)+2*xi*(ol(Cart::zz,Cart::xz)-0.5*rzetaA*ol(Cart::s,Cart::xz));
            kin(Cart::zz,Cart::zz) = PmA[2]*kin(Cart::z,Cart::zz)+1*_fak*kin(Cart::s,Cart::zz)+2*_fak*kin(Cart::z,Cart::z)+2*xi*(ol(Cart::zz,Cart::zz)-0.5*rzetaA*ol(Cart::s,Cart::zz));
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


             // cartesian -> spherical

             ub::matrix<double> kin_tmp = ub::prod( _trafo_row, kin );
             ub::matrix<double> _trafo_col_tposed = ub::trans( _trafo_col );
             ub::matrix<double> kin_sph = ub::prod( kin_tmp, _trafo_col_tposed );
             // save to _matrix
             for ( unsigned i = 0; i< _matrix.size1(); i++ ) {
                 for (unsigned j = 0; j < _matrix.size2(); j++){
                     _matrix(i,j) += kin_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
                    }
            }
        
        
        
        
        
                }//col
            }//row
        }
    }
}
    
