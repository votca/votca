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

    
    void AOKinetic::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, AOBasis* ecp) {
        const double pi = boost::math::constants::pi<double>();
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
            ol(Cartesian::s,Cartesian::s) = pow(rzeta,1.5)*pow(4.0*_decay_row*_decay_col,0.75) * exp(-_exparg);
            // s-s- kinetic energy integral
            kin(Cartesian::s,Cartesian::s)= ol(Cartesian::s,Cartesian::s)*xi*(3-2*xi*_distsq);
            
          //s-p
            if ( _lmax_col>0 ){
            //std::cout << "\t setting s-p|" << std::flush;      
            ol(Cartesian::s,Cartesian::y) = PmB[1]*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::s,Cartesian::x) = PmB[0]*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::s,Cartesian::z) = PmB[2]*ol(Cartesian::s,Cartesian::s);
            
            kin(Cartesian::s,Cartesian::y) = PmB[1]*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::s,Cartesian::y));
            kin(Cartesian::s,Cartesian::x) = PmB[0]*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::s,Cartesian::x));
            kin(Cartesian::s,Cartesian::z) = PmB[2]*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::s,Cartesian::z));
            }
            
            //p-s
            if ( _lmax_row>0 ){
            //std::cout << "\t setting p-s|" << std::flush;
            ol(Cartesian::y,Cartesian::s) = PmA[1]*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::x,Cartesian::s) = PmA[0]*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::z,Cartesian::s) = PmA[2]*ol(Cartesian::s,Cartesian::s);
       
            kin(Cartesian::y,Cartesian::s) = PmA[1]*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::y,Cartesian::s));
            kin(Cartesian::x,Cartesian::s) = PmA[0]*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::x,Cartesian::s));
            kin(Cartesian::z,Cartesian::s) = PmA[2]*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::z,Cartesian::s));
            }
            
            //p-p
            if ( _lmax_row>0 && _lmax_col>0 ){
            ol(Cartesian::y,Cartesian::y) = PmA[1]*ol(Cartesian::s,Cartesian::y)+1*_fak*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::y,Cartesian::x) = PmA[1]*ol(Cartesian::s,Cartesian::x);
            ol(Cartesian::y,Cartesian::z) = PmA[1]*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::x,Cartesian::y) = PmA[0]*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::x,Cartesian::x) = PmA[0]*ol(Cartesian::s,Cartesian::x)+1*_fak*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::x,Cartesian::z) = PmA[0]*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::z,Cartesian::y) = PmA[2]*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::z,Cartesian::x) = PmA[2]*ol(Cartesian::s,Cartesian::x);
            ol(Cartesian::z,Cartesian::z) = PmA[2]*ol(Cartesian::s,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::s);
     
            //std::cout << "\t setting p-p|" << std::flush;
            kin(Cartesian::y,Cartesian::y) = PmA[1]*kin(Cartesian::s,Cartesian::y)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::y,Cartesian::y));
            kin(Cartesian::y,Cartesian::x) = PmA[1]*kin(Cartesian::s,Cartesian::x)+2*xi*(ol(Cartesian::y,Cartesian::x));
            kin(Cartesian::y,Cartesian::z) = PmA[1]*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::y,Cartesian::z));
            kin(Cartesian::x,Cartesian::y) = PmA[0]*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::x,Cartesian::y));
            kin(Cartesian::x,Cartesian::x) = PmA[0]*kin(Cartesian::s,Cartesian::x)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::x,Cartesian::x));
            kin(Cartesian::x,Cartesian::z) = PmA[0]*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::x,Cartesian::z));
            kin(Cartesian::z,Cartesian::y) = PmA[2]*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::z,Cartesian::y));
            kin(Cartesian::z,Cartesian::x) = PmA[2]*kin(Cartesian::s,Cartesian::x)+2*xi*(ol(Cartesian::z,Cartesian::x));
            kin(Cartesian::z,Cartesian::z) = PmA[2]*kin(Cartesian::s,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::z,Cartesian::z));
            }
            
            //s-d
            if ( _lmax_col>1 ){
            ol(Cartesian::s,Cartesian::yy) = PmB[1]*ol(Cartesian::s,Cartesian::y)+1*_fak*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::s,Cartesian::xy) = PmB[0]*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::s,Cartesian::yz) = PmB[1]*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::s,Cartesian::xx) = PmB[0]*ol(Cartesian::s,Cartesian::x)+1*_fak*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::s,Cartesian::xz) = PmB[0]*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::s,Cartesian::zz) = PmB[2]*ol(Cartesian::s,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::s);    
                
            //std::cout << "\t setting s-d|" << std::flush;
            kin(Cartesian::s,Cartesian::yy) = PmB[1]*kin(Cartesian::s,Cartesian::y)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::s,Cartesian::yy)-0.5*rzetaB*ol(Cartesian::s,Cartesian::s));
            kin(Cartesian::s,Cartesian::xy) = PmB[0]*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::s,Cartesian::xy));
            kin(Cartesian::s,Cartesian::yz) = PmB[1]*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::s,Cartesian::yz));
            kin(Cartesian::s,Cartesian::xx) = PmB[0]*kin(Cartesian::s,Cartesian::x)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::s,Cartesian::xx)-0.5*rzetaB*ol(Cartesian::s,Cartesian::s));
            kin(Cartesian::s,Cartesian::xz) = PmB[0]*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::s,Cartesian::xz));
            kin(Cartesian::s,Cartesian::zz) = PmB[2]*kin(Cartesian::s,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::s,Cartesian::zz)-0.5*rzetaB*ol(Cartesian::s,Cartesian::s));            
            }

            //d-s
            if ( _lmax_row>1 ){
            //std::cout << "\t setting d-s|" << std::flush;
            ol(Cartesian::yy,Cartesian::s) = PmA[1]*ol(Cartesian::y,Cartesian::s)+1*_fak*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::xy,Cartesian::s) = PmA[0]*ol(Cartesian::y,Cartesian::s);
            ol(Cartesian::yz,Cartesian::s) = PmA[1]*ol(Cartesian::z,Cartesian::s);
            ol(Cartesian::xx,Cartesian::s) = PmA[0]*ol(Cartesian::x,Cartesian::s)+1*_fak*ol(Cartesian::s,Cartesian::s);
            ol(Cartesian::xz,Cartesian::s) = PmA[0]*ol(Cartesian::z,Cartesian::s);
            ol(Cartesian::zz,Cartesian::s) = PmA[2]*ol(Cartesian::z,Cartesian::s)+1*_fak*ol(Cartesian::s,Cartesian::s);
            
            kin(Cartesian::yy,Cartesian::s) = PmA[1]*kin(Cartesian::y,Cartesian::s)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::yy,Cartesian::s)-0.5*rzetaA*ol(Cartesian::s,Cartesian::s));
            kin(Cartesian::xy,Cartesian::s) = PmA[0]*kin(Cartesian::y,Cartesian::s)+2*xi*(ol(Cartesian::xy,Cartesian::s));
            kin(Cartesian::yz,Cartesian::s) = PmA[1]*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::yz,Cartesian::s));
            kin(Cartesian::xx,Cartesian::s) = PmA[0]*kin(Cartesian::x,Cartesian::s)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::xx,Cartesian::s)-0.5*rzetaA*ol(Cartesian::s,Cartesian::s));
            kin(Cartesian::xz,Cartesian::s) = PmA[0]*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::xz,Cartesian::s));
            kin(Cartesian::zz,Cartesian::s) = PmA[2]*kin(Cartesian::z,Cartesian::s)+1*_fak*kin(Cartesian::s,Cartesian::s)+2*xi*(ol(Cartesian::zz,Cartesian::s)-0.5*rzetaA*ol(Cartesian::s,Cartesian::s));
            }
            
            //p-d
            if ( _lmax_row>0 && _lmax_col>1 ){
            //std::cout << "\t setting p-d|" << std::flush;
            ol(Cartesian::y,Cartesian::yy) = PmB[1]*ol(Cartesian::y,Cartesian::y)+1*_fak*ol(Cartesian::y,Cartesian::s)+1*_fak*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::y,Cartesian::xy) = PmB[0]*ol(Cartesian::y,Cartesian::y);
            ol(Cartesian::y,Cartesian::yz) = PmB[1]*ol(Cartesian::y,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::y,Cartesian::xx) = PmB[0]*ol(Cartesian::y,Cartesian::x)+1*_fak*ol(Cartesian::y,Cartesian::s);
            ol(Cartesian::y,Cartesian::xz) = PmB[0]*ol(Cartesian::y,Cartesian::z);
            ol(Cartesian::y,Cartesian::zz) = PmB[2]*ol(Cartesian::y,Cartesian::z)+1*_fak*ol(Cartesian::y,Cartesian::s);
            ol(Cartesian::x,Cartesian::yy) = PmB[1]*ol(Cartesian::x,Cartesian::y)+1*_fak*ol(Cartesian::x,Cartesian::s);
            ol(Cartesian::x,Cartesian::xy) = PmB[0]*ol(Cartesian::x,Cartesian::y)+1*_fak*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::x,Cartesian::yz) = PmB[1]*ol(Cartesian::x,Cartesian::z);
            ol(Cartesian::x,Cartesian::xx) = PmB[0]*ol(Cartesian::x,Cartesian::x)+1*_fak*ol(Cartesian::x,Cartesian::s)+1*_fak*ol(Cartesian::s,Cartesian::x);
            ol(Cartesian::x,Cartesian::xz) = PmB[0]*ol(Cartesian::x,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::x,Cartesian::zz) = PmB[2]*ol(Cartesian::x,Cartesian::z)+1*_fak*ol(Cartesian::x,Cartesian::s);
            ol(Cartesian::z,Cartesian::yy) = PmB[1]*ol(Cartesian::z,Cartesian::y)+1*_fak*ol(Cartesian::z,Cartesian::s);
            ol(Cartesian::z,Cartesian::xy) = PmB[0]*ol(Cartesian::z,Cartesian::y);
            ol(Cartesian::z,Cartesian::yz) = PmB[1]*ol(Cartesian::z,Cartesian::z);
            ol(Cartesian::z,Cartesian::xx) = PmB[0]*ol(Cartesian::z,Cartesian::x)+1*_fak*ol(Cartesian::z,Cartesian::s);
            ol(Cartesian::z,Cartesian::xz) = PmB[0]*ol(Cartesian::z,Cartesian::z);
            ol(Cartesian::z,Cartesian::zz) = PmB[2]*ol(Cartesian::z,Cartesian::z)+1*_fak*ol(Cartesian::z,Cartesian::s)+1*_fak*ol(Cartesian::s,Cartesian::z);
            
            kin(Cartesian::y,Cartesian::yy) = PmB[1]*kin(Cartesian::y,Cartesian::y)+1*_fak*kin(Cartesian::y,Cartesian::s)+1*_fak*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::y,Cartesian::yy)-0.5*rzetaB*ol(Cartesian::y,Cartesian::s));
            kin(Cartesian::y,Cartesian::xy) = PmB[0]*kin(Cartesian::y,Cartesian::y)+2*xi*(ol(Cartesian::y,Cartesian::xy));
            kin(Cartesian::y,Cartesian::yz) = PmB[1]*kin(Cartesian::y,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::y,Cartesian::yz));
            kin(Cartesian::y,Cartesian::xx) = PmB[0]*kin(Cartesian::y,Cartesian::x)+1*_fak*kin(Cartesian::y,Cartesian::s)+2*xi*(ol(Cartesian::y,Cartesian::xx)-0.5*rzetaB*ol(Cartesian::y,Cartesian::s));
            kin(Cartesian::y,Cartesian::xz) = PmB[0]*kin(Cartesian::y,Cartesian::z)+2*xi*(ol(Cartesian::y,Cartesian::xz));
            kin(Cartesian::y,Cartesian::zz) = PmB[2]*kin(Cartesian::y,Cartesian::z)+1*_fak*kin(Cartesian::y,Cartesian::s)+2*xi*(ol(Cartesian::y,Cartesian::zz)-0.5*rzetaB*ol(Cartesian::y,Cartesian::s));
            kin(Cartesian::x,Cartesian::yy) = PmB[1]*kin(Cartesian::x,Cartesian::y)+1*_fak*kin(Cartesian::x,Cartesian::s)+2*xi*(ol(Cartesian::x,Cartesian::yy)-0.5*rzetaB*ol(Cartesian::x,Cartesian::s));
            kin(Cartesian::x,Cartesian::xy) = PmB[0]*kin(Cartesian::x,Cartesian::y)+1*_fak*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::x,Cartesian::xy));
            kin(Cartesian::x,Cartesian::yz) = PmB[1]*kin(Cartesian::x,Cartesian::z)+2*xi*(ol(Cartesian::x,Cartesian::yz));
            kin(Cartesian::x,Cartesian::xx) = PmB[0]*kin(Cartesian::x,Cartesian::x)+1*_fak*kin(Cartesian::x,Cartesian::s)+1*_fak*kin(Cartesian::s,Cartesian::x)+2*xi*(ol(Cartesian::x,Cartesian::xx)-0.5*rzetaB*ol(Cartesian::x,Cartesian::s));
            kin(Cartesian::x,Cartesian::xz) = PmB[0]*kin(Cartesian::x,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::x,Cartesian::xz));
            kin(Cartesian::x,Cartesian::zz) = PmB[2]*kin(Cartesian::x,Cartesian::z)+1*_fak*kin(Cartesian::x,Cartesian::s)+2*xi*(ol(Cartesian::x,Cartesian::zz)-0.5*rzetaB*ol(Cartesian::x,Cartesian::s));
            kin(Cartesian::z,Cartesian::yy) = PmB[1]*kin(Cartesian::z,Cartesian::y)+1*_fak*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::z,Cartesian::yy)-0.5*rzetaB*ol(Cartesian::z,Cartesian::s));
            kin(Cartesian::z,Cartesian::xy) = PmB[0]*kin(Cartesian::z,Cartesian::y)+2*xi*(ol(Cartesian::z,Cartesian::xy));
            kin(Cartesian::z,Cartesian::yz) = PmB[1]*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::z,Cartesian::yz));
            kin(Cartesian::z,Cartesian::xx) = PmB[0]*kin(Cartesian::z,Cartesian::x)+1*_fak*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::z,Cartesian::xx)-0.5*rzetaB*ol(Cartesian::z,Cartesian::s));
            kin(Cartesian::z,Cartesian::xz) = PmB[0]*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::z,Cartesian::xz));
            kin(Cartesian::z,Cartesian::zz) = PmB[2]*kin(Cartesian::z,Cartesian::z)+1*_fak*kin(Cartesian::z,Cartesian::s)+1*_fak*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::z,Cartesian::zz)-0.5*rzetaB*ol(Cartesian::z,Cartesian::s));
            }
            
            //d-p
            if ( _lmax_row>1 && _lmax_col>0 ){
            //std::cout << "\t setting d-p|" << std::flush;
            ol(Cartesian::yy,Cartesian::y) = PmA[1]*ol(Cartesian::y,Cartesian::y)+1*_fak*ol(Cartesian::s,Cartesian::y)+1*_fak*ol(Cartesian::y,Cartesian::s);
            ol(Cartesian::yy,Cartesian::x) = PmA[1]*ol(Cartesian::y,Cartesian::x)+1*_fak*ol(Cartesian::s,Cartesian::x);
            ol(Cartesian::yy,Cartesian::z) = PmA[1]*ol(Cartesian::y,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::xy,Cartesian::y) = PmA[0]*ol(Cartesian::y,Cartesian::y);
            ol(Cartesian::xy,Cartesian::x) = PmA[0]*ol(Cartesian::y,Cartesian::x)+1*_fak*ol(Cartesian::y,Cartesian::s);
            ol(Cartesian::xy,Cartesian::z) = PmA[0]*ol(Cartesian::y,Cartesian::z);
            ol(Cartesian::yz,Cartesian::y) = PmA[1]*ol(Cartesian::z,Cartesian::y)+1*_fak*ol(Cartesian::z,Cartesian::s);
            ol(Cartesian::yz,Cartesian::x) = PmA[1]*ol(Cartesian::z,Cartesian::x);
            ol(Cartesian::yz,Cartesian::z) = PmA[1]*ol(Cartesian::z,Cartesian::z);
            ol(Cartesian::xx,Cartesian::y) = PmA[0]*ol(Cartesian::x,Cartesian::y)+1*_fak*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::xx,Cartesian::x) = PmA[0]*ol(Cartesian::x,Cartesian::x)+1*_fak*ol(Cartesian::s,Cartesian::x)+1*_fak*ol(Cartesian::x,Cartesian::s);
            ol(Cartesian::xx,Cartesian::z) = PmA[0]*ol(Cartesian::x,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::z);
            ol(Cartesian::xz,Cartesian::y) = PmA[0]*ol(Cartesian::z,Cartesian::y);
            ol(Cartesian::xz,Cartesian::x) = PmA[0]*ol(Cartesian::z,Cartesian::x)+1*_fak*ol(Cartesian::z,Cartesian::s);
            ol(Cartesian::xz,Cartesian::z) = PmA[0]*ol(Cartesian::z,Cartesian::z);
            ol(Cartesian::zz,Cartesian::y) = PmA[2]*ol(Cartesian::z,Cartesian::y)+1*_fak*ol(Cartesian::s,Cartesian::y);
            ol(Cartesian::zz,Cartesian::x) = PmA[2]*ol(Cartesian::z,Cartesian::x)+1*_fak*ol(Cartesian::s,Cartesian::x);
            ol(Cartesian::zz,Cartesian::z) = PmA[2]*ol(Cartesian::z,Cartesian::z)+1*_fak*ol(Cartesian::s,Cartesian::z)+1*_fak*ol(Cartesian::z,Cartesian::s);
            
            kin(Cartesian::yy,Cartesian::y) = PmA[1]*kin(Cartesian::y,Cartesian::y)+1*_fak*kin(Cartesian::s,Cartesian::y)+1*_fak*kin(Cartesian::y,Cartesian::s)+2*xi*(ol(Cartesian::yy,Cartesian::y)-0.5*rzetaA*ol(Cartesian::s,Cartesian::y));
            kin(Cartesian::yy,Cartesian::x) = PmA[1]*kin(Cartesian::y,Cartesian::x)+1*_fak*kin(Cartesian::s,Cartesian::x)+2*xi*(ol(Cartesian::yy,Cartesian::x)-0.5*rzetaA*ol(Cartesian::s,Cartesian::x));
            kin(Cartesian::yy,Cartesian::z) = PmA[1]*kin(Cartesian::y,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::yy,Cartesian::z)-0.5*rzetaA*ol(Cartesian::s,Cartesian::z));
            kin(Cartesian::xy,Cartesian::y) = PmA[0]*kin(Cartesian::y,Cartesian::y)+2*xi*(ol(Cartesian::xy,Cartesian::y));
            kin(Cartesian::xy,Cartesian::x) = PmA[0]*kin(Cartesian::y,Cartesian::x)+1*_fak*kin(Cartesian::y,Cartesian::s)+2*xi*(ol(Cartesian::xy,Cartesian::x));
            kin(Cartesian::xy,Cartesian::z) = PmA[0]*kin(Cartesian::y,Cartesian::z)+2*xi*(ol(Cartesian::xy,Cartesian::z));
            kin(Cartesian::yz,Cartesian::y) = PmA[1]*kin(Cartesian::z,Cartesian::y)+1*_fak*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::yz,Cartesian::y));
            kin(Cartesian::yz,Cartesian::x) = PmA[1]*kin(Cartesian::z,Cartesian::x)+2*xi*(ol(Cartesian::yz,Cartesian::x));
            kin(Cartesian::yz,Cartesian::z) = PmA[1]*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::yz,Cartesian::z));
            kin(Cartesian::xx,Cartesian::y) = PmA[0]*kin(Cartesian::x,Cartesian::y)+1*_fak*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::xx,Cartesian::y)-0.5*rzetaA*ol(Cartesian::s,Cartesian::y));
            kin(Cartesian::xx,Cartesian::x) = PmA[0]*kin(Cartesian::x,Cartesian::x)+1*_fak*kin(Cartesian::s,Cartesian::x)+1*_fak*kin(Cartesian::x,Cartesian::s)+2*xi*(ol(Cartesian::xx,Cartesian::x)-0.5*rzetaA*ol(Cartesian::s,Cartesian::x));
            kin(Cartesian::xx,Cartesian::z) = PmA[0]*kin(Cartesian::x,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::z)+2*xi*(ol(Cartesian::xx,Cartesian::z)-0.5*rzetaA*ol(Cartesian::s,Cartesian::z));
            kin(Cartesian::xz,Cartesian::y) = PmA[0]*kin(Cartesian::z,Cartesian::y)+2*xi*(ol(Cartesian::xz,Cartesian::y));
            kin(Cartesian::xz,Cartesian::x) = PmA[0]*kin(Cartesian::z,Cartesian::x)+1*_fak*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::xz,Cartesian::x));
            kin(Cartesian::xz,Cartesian::z) = PmA[0]*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::xz,Cartesian::z));
            kin(Cartesian::zz,Cartesian::y) = PmA[2]*kin(Cartesian::z,Cartesian::y)+1*_fak*kin(Cartesian::s,Cartesian::y)+2*xi*(ol(Cartesian::zz,Cartesian::y)-0.5*rzetaA*ol(Cartesian::s,Cartesian::y));
            kin(Cartesian::zz,Cartesian::x) = PmA[2]*kin(Cartesian::z,Cartesian::x)+1*_fak*kin(Cartesian::s,Cartesian::x)+2*xi*(ol(Cartesian::zz,Cartesian::x)-0.5*rzetaA*ol(Cartesian::s,Cartesian::x));
            kin(Cartesian::zz,Cartesian::z) = PmA[2]*kin(Cartesian::z,Cartesian::z)+1*_fak*kin(Cartesian::s,Cartesian::z)+1*_fak*kin(Cartesian::z,Cartesian::s)+2*xi*(ol(Cartesian::zz,Cartesian::z)-0.5*rzetaA*ol(Cartesian::s,Cartesian::z));
            }

            //d-d
            if ( _lmax_row>1 && _lmax_col>1 ){
            //std::cout << "\t setting d-d|" << std::flush;
            ol(Cartesian::yy,Cartesian::yy) = PmA[1]*ol(Cartesian::y,Cartesian::yy)+1*_fak*ol(Cartesian::s,Cartesian::yy)+2*_fak*ol(Cartesian::y,Cartesian::y);
            ol(Cartesian::yy,Cartesian::xy) = PmA[1]*ol(Cartesian::y,Cartesian::xy)+1*_fak*ol(Cartesian::s,Cartesian::xy)+1*_fak*ol(Cartesian::y,Cartesian::x);
            ol(Cartesian::yy,Cartesian::yz) = PmA[1]*ol(Cartesian::y,Cartesian::yz)+1*_fak*ol(Cartesian::s,Cartesian::yz)+1*_fak*ol(Cartesian::y,Cartesian::z);
            ol(Cartesian::yy,Cartesian::xx) = PmA[1]*ol(Cartesian::y,Cartesian::xx)+1*_fak*ol(Cartesian::s,Cartesian::xx);
            ol(Cartesian::yy,Cartesian::xz) = PmA[1]*ol(Cartesian::y,Cartesian::xz)+1*_fak*ol(Cartesian::s,Cartesian::xz);
            ol(Cartesian::yy,Cartesian::zz) = PmA[1]*ol(Cartesian::y,Cartesian::zz)+1*_fak*ol(Cartesian::s,Cartesian::zz);
            ol(Cartesian::xy,Cartesian::yy) = PmA[0]*ol(Cartesian::y,Cartesian::yy);
            ol(Cartesian::xy,Cartesian::xy) = PmA[0]*ol(Cartesian::y,Cartesian::xy)+1*_fak*ol(Cartesian::y,Cartesian::y);
            ol(Cartesian::xy,Cartesian::yz) = PmA[0]*ol(Cartesian::y,Cartesian::yz);
            ol(Cartesian::xy,Cartesian::xx) = PmA[0]*ol(Cartesian::y,Cartesian::xx)+2*_fak*ol(Cartesian::y,Cartesian::x);
            ol(Cartesian::xy,Cartesian::xz) = PmA[0]*ol(Cartesian::y,Cartesian::xz)+1*_fak*ol(Cartesian::y,Cartesian::z);
            ol(Cartesian::xy,Cartesian::zz) = PmA[0]*ol(Cartesian::y,Cartesian::zz);
            ol(Cartesian::yz,Cartesian::yy) = PmA[1]*ol(Cartesian::z,Cartesian::yy)+2*_fak*ol(Cartesian::z,Cartesian::y);
            ol(Cartesian::yz,Cartesian::xy) = PmA[1]*ol(Cartesian::z,Cartesian::xy)+1*_fak*ol(Cartesian::z,Cartesian::x);
            ol(Cartesian::yz,Cartesian::yz) = PmA[1]*ol(Cartesian::z,Cartesian::yz)+1*_fak*ol(Cartesian::z,Cartesian::z);
            ol(Cartesian::yz,Cartesian::xx) = PmA[1]*ol(Cartesian::z,Cartesian::xx);
            ol(Cartesian::yz,Cartesian::xz) = PmA[1]*ol(Cartesian::z,Cartesian::xz);
            ol(Cartesian::yz,Cartesian::zz) = PmA[1]*ol(Cartesian::z,Cartesian::zz);
            ol(Cartesian::xx,Cartesian::yy) = PmA[0]*ol(Cartesian::x,Cartesian::yy)+1*_fak*ol(Cartesian::s,Cartesian::yy);
            ol(Cartesian::xx,Cartesian::xy) = PmA[0]*ol(Cartesian::x,Cartesian::xy)+1*_fak*ol(Cartesian::s,Cartesian::xy)+1*_fak*ol(Cartesian::x,Cartesian::y);
            ol(Cartesian::xx,Cartesian::yz) = PmA[0]*ol(Cartesian::x,Cartesian::yz)+1*_fak*ol(Cartesian::s,Cartesian::yz);
            ol(Cartesian::xx,Cartesian::xx) = PmA[0]*ol(Cartesian::x,Cartesian::xx)+1*_fak*ol(Cartesian::s,Cartesian::xx)+2*_fak*ol(Cartesian::x,Cartesian::x);
            ol(Cartesian::xx,Cartesian::xz) = PmA[0]*ol(Cartesian::x,Cartesian::xz)+1*_fak*ol(Cartesian::s,Cartesian::xz)+1*_fak*ol(Cartesian::x,Cartesian::z);
            ol(Cartesian::xx,Cartesian::zz) = PmA[0]*ol(Cartesian::x,Cartesian::zz)+1*_fak*ol(Cartesian::s,Cartesian::zz);
            ol(Cartesian::xz,Cartesian::yy) = PmA[0]*ol(Cartesian::z,Cartesian::yy);
            ol(Cartesian::xz,Cartesian::xy) = PmA[0]*ol(Cartesian::z,Cartesian::xy)+1*_fak*ol(Cartesian::z,Cartesian::y);
            ol(Cartesian::xz,Cartesian::yz) = PmA[0]*ol(Cartesian::z,Cartesian::yz);
            ol(Cartesian::xz,Cartesian::xx) = PmA[0]*ol(Cartesian::z,Cartesian::xx)+2*_fak*ol(Cartesian::z,Cartesian::x);
            ol(Cartesian::xz,Cartesian::xz) = PmA[0]*ol(Cartesian::z,Cartesian::xz)+1*_fak*ol(Cartesian::z,Cartesian::z);
            ol(Cartesian::xz,Cartesian::zz) = PmA[0]*ol(Cartesian::z,Cartesian::zz);
            ol(Cartesian::zz,Cartesian::yy) = PmA[2]*ol(Cartesian::z,Cartesian::yy)+1*_fak*ol(Cartesian::s,Cartesian::yy);
            ol(Cartesian::zz,Cartesian::xy) = PmA[2]*ol(Cartesian::z,Cartesian::xy)+1*_fak*ol(Cartesian::s,Cartesian::xy);
            ol(Cartesian::zz,Cartesian::yz) = PmA[2]*ol(Cartesian::z,Cartesian::yz)+1*_fak*ol(Cartesian::s,Cartesian::yz)+1*_fak*ol(Cartesian::z,Cartesian::y);
            ol(Cartesian::zz,Cartesian::xx) = PmA[2]*ol(Cartesian::z,Cartesian::xx)+1*_fak*ol(Cartesian::s,Cartesian::xx);
            ol(Cartesian::zz,Cartesian::xz) = PmA[2]*ol(Cartesian::z,Cartesian::xz)+1*_fak*ol(Cartesian::s,Cartesian::xz)+1*_fak*ol(Cartesian::z,Cartesian::x);
            ol(Cartesian::zz,Cartesian::zz) = PmA[2]*ol(Cartesian::z,Cartesian::zz)+1*_fak*ol(Cartesian::s,Cartesian::zz)+2*_fak*ol(Cartesian::z,Cartesian::z);

            kin(Cartesian::yy,Cartesian::yy) = PmA[1]*kin(Cartesian::y,Cartesian::yy)+1*_fak*kin(Cartesian::s,Cartesian::yy)+2*_fak*kin(Cartesian::y,Cartesian::y)+2*xi*(ol(Cartesian::yy,Cartesian::yy)-0.5*rzetaA*ol(Cartesian::s,Cartesian::yy));
            kin(Cartesian::yy,Cartesian::xy) = PmA[1]*kin(Cartesian::y,Cartesian::xy)+1*_fak*kin(Cartesian::s,Cartesian::xy)+1*_fak*kin(Cartesian::y,Cartesian::x)+2*xi*(ol(Cartesian::yy,Cartesian::xy)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xy));
            kin(Cartesian::yy,Cartesian::yz) = PmA[1]*kin(Cartesian::y,Cartesian::yz)+1*_fak*kin(Cartesian::s,Cartesian::yz)+1*_fak*kin(Cartesian::y,Cartesian::z)+2*xi*(ol(Cartesian::yy,Cartesian::yz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::yz));
            kin(Cartesian::yy,Cartesian::xx) = PmA[1]*kin(Cartesian::y,Cartesian::xx)+1*_fak*kin(Cartesian::s,Cartesian::xx)+2*xi*(ol(Cartesian::yy,Cartesian::xx)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xx));
            kin(Cartesian::yy,Cartesian::xz) = PmA[1]*kin(Cartesian::y,Cartesian::xz)+1*_fak*kin(Cartesian::s,Cartesian::xz)+2*xi*(ol(Cartesian::yy,Cartesian::xz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xz));
            kin(Cartesian::yy,Cartesian::zz) = PmA[1]*kin(Cartesian::y,Cartesian::zz)+1*_fak*kin(Cartesian::s,Cartesian::zz)+2*xi*(ol(Cartesian::yy,Cartesian::zz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::zz));
            kin(Cartesian::xy,Cartesian::yy) = PmA[0]*kin(Cartesian::y,Cartesian::yy)+2*xi*(ol(Cartesian::xy,Cartesian::yy));
            kin(Cartesian::xy,Cartesian::xy) = PmA[0]*kin(Cartesian::y,Cartesian::xy)+1*_fak*kin(Cartesian::y,Cartesian::y)+2*xi*(ol(Cartesian::xy,Cartesian::xy));
            kin(Cartesian::xy,Cartesian::yz) = PmA[0]*kin(Cartesian::y,Cartesian::yz)+2*xi*(ol(Cartesian::xy,Cartesian::yz));
            kin(Cartesian::xy,Cartesian::xx) = PmA[0]*kin(Cartesian::y,Cartesian::xx)+2*_fak*kin(Cartesian::y,Cartesian::x)+2*xi*(ol(Cartesian::xy,Cartesian::xx));
            kin(Cartesian::xy,Cartesian::xz) = PmA[0]*kin(Cartesian::y,Cartesian::xz)+1*_fak*kin(Cartesian::y,Cartesian::z)+2*xi*(ol(Cartesian::xy,Cartesian::xz));
            kin(Cartesian::xy,Cartesian::zz) = PmA[0]*kin(Cartesian::y,Cartesian::zz)+2*xi*(ol(Cartesian::xy,Cartesian::zz));
            kin(Cartesian::yz,Cartesian::yy) = PmA[1]*kin(Cartesian::z,Cartesian::yy)+2*_fak*kin(Cartesian::z,Cartesian::y)+2*xi*(ol(Cartesian::yz,Cartesian::yy));
            kin(Cartesian::yz,Cartesian::xy) = PmA[1]*kin(Cartesian::z,Cartesian::xy)+1*_fak*kin(Cartesian::z,Cartesian::x)+2*xi*(ol(Cartesian::yz,Cartesian::xy));
            kin(Cartesian::yz,Cartesian::yz) = PmA[1]*kin(Cartesian::z,Cartesian::yz)+1*_fak*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::yz,Cartesian::yz));
            kin(Cartesian::yz,Cartesian::xx) = PmA[1]*kin(Cartesian::z,Cartesian::xx)+2*xi*(ol(Cartesian::yz,Cartesian::xx));
            kin(Cartesian::yz,Cartesian::xz) = PmA[1]*kin(Cartesian::z,Cartesian::xz)+2*xi*(ol(Cartesian::yz,Cartesian::xz));
            kin(Cartesian::yz,Cartesian::zz) = PmA[1]*kin(Cartesian::z,Cartesian::zz)+2*xi*(ol(Cartesian::yz,Cartesian::zz));
            kin(Cartesian::xx,Cartesian::yy) = PmA[0]*kin(Cartesian::x,Cartesian::yy)+1*_fak*kin(Cartesian::s,Cartesian::yy)+2*xi*(ol(Cartesian::xx,Cartesian::yy)-0.5*rzetaA*ol(Cartesian::s,Cartesian::yy));
            kin(Cartesian::xx,Cartesian::xy) = PmA[0]*kin(Cartesian::x,Cartesian::xy)+1*_fak*kin(Cartesian::s,Cartesian::xy)+1*_fak*kin(Cartesian::x,Cartesian::y)+2*xi*(ol(Cartesian::xx,Cartesian::xy)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xy));
            kin(Cartesian::xx,Cartesian::yz) = PmA[0]*kin(Cartesian::x,Cartesian::yz)+1*_fak*kin(Cartesian::s,Cartesian::yz)+2*xi*(ol(Cartesian::xx,Cartesian::yz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::yz));
            kin(Cartesian::xx,Cartesian::xx) = PmA[0]*kin(Cartesian::x,Cartesian::xx)+1*_fak*kin(Cartesian::s,Cartesian::xx)+2*_fak*kin(Cartesian::x,Cartesian::x)+2*xi*(ol(Cartesian::xx,Cartesian::xx)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xx));
            kin(Cartesian::xx,Cartesian::xz) = PmA[0]*kin(Cartesian::x,Cartesian::xz)+1*_fak*kin(Cartesian::s,Cartesian::xz)+1*_fak*kin(Cartesian::x,Cartesian::z)+2*xi*(ol(Cartesian::xx,Cartesian::xz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xz));
            kin(Cartesian::xx,Cartesian::zz) = PmA[0]*kin(Cartesian::x,Cartesian::zz)+1*_fak*kin(Cartesian::s,Cartesian::zz)+2*xi*(ol(Cartesian::xx,Cartesian::zz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::zz));
            kin(Cartesian::xz,Cartesian::yy) = PmA[0]*kin(Cartesian::z,Cartesian::yy)+2*xi*(ol(Cartesian::xz,Cartesian::yy));
            kin(Cartesian::xz,Cartesian::xy) = PmA[0]*kin(Cartesian::z,Cartesian::xy)+1*_fak*kin(Cartesian::z,Cartesian::y)+2*xi*(ol(Cartesian::xz,Cartesian::xy));
            kin(Cartesian::xz,Cartesian::yz) = PmA[0]*kin(Cartesian::z,Cartesian::yz)+2*xi*(ol(Cartesian::xz,Cartesian::yz));
            kin(Cartesian::xz,Cartesian::xx) = PmA[0]*kin(Cartesian::z,Cartesian::xx)+2*_fak*kin(Cartesian::z,Cartesian::x)+2*xi*(ol(Cartesian::xz,Cartesian::xx));
            kin(Cartesian::xz,Cartesian::xz) = PmA[0]*kin(Cartesian::z,Cartesian::xz)+1*_fak*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::xz,Cartesian::xz));
            kin(Cartesian::xz,Cartesian::zz) = PmA[0]*kin(Cartesian::z,Cartesian::zz)+2*xi*(ol(Cartesian::xz,Cartesian::zz));
            kin(Cartesian::zz,Cartesian::yy) = PmA[2]*kin(Cartesian::z,Cartesian::yy)+1*_fak*kin(Cartesian::s,Cartesian::yy)+2*xi*(ol(Cartesian::zz,Cartesian::yy)-0.5*rzetaA*ol(Cartesian::s,Cartesian::yy));
            kin(Cartesian::zz,Cartesian::xy) = PmA[2]*kin(Cartesian::z,Cartesian::xy)+1*_fak*kin(Cartesian::s,Cartesian::xy)+2*xi*(ol(Cartesian::zz,Cartesian::xy)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xy));
            kin(Cartesian::zz,Cartesian::yz) = PmA[2]*kin(Cartesian::z,Cartesian::yz)+1*_fak*kin(Cartesian::s,Cartesian::yz)+1*_fak*kin(Cartesian::z,Cartesian::y)+2*xi*(ol(Cartesian::zz,Cartesian::yz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::yz));
            kin(Cartesian::zz,Cartesian::xx) = PmA[2]*kin(Cartesian::z,Cartesian::xx)+1*_fak*kin(Cartesian::s,Cartesian::xx)+2*xi*(ol(Cartesian::zz,Cartesian::xx)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xx));
            kin(Cartesian::zz,Cartesian::xz) = PmA[2]*kin(Cartesian::z,Cartesian::xz)+1*_fak*kin(Cartesian::s,Cartesian::xz)+1*_fak*kin(Cartesian::z,Cartesian::x)+2*xi*(ol(Cartesian::zz,Cartesian::xz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::xz));
            kin(Cartesian::zz,Cartesian::zz) = PmA[2]*kin(Cartesian::z,Cartesian::zz)+1*_fak*kin(Cartesian::s,Cartesian::zz)+2*_fak*kin(Cartesian::z,Cartesian::z)+2*xi*(ol(Cartesian::zz,Cartesian::zz)-0.5*rzetaA*ol(Cartesian::s,Cartesian::zz));
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
             for ( int i = 0; i< _matrix.size1(); i++ ) {
                 for (int j = 0; j < _matrix.size2(); j++){
                     _matrix(i,j) += kin_sph(i+_shell_row->getOffset(),j+_shell_col->getOffset());
                    }
            }
        
        
        
        
        
                }//col
            }//row
        }
    }
}
    
