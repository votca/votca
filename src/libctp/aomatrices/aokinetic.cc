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

    
    void AOKinetic::FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col , bool _raw) {
        
            /*cout << "\nAO block: "<< endl;
        cout << "\t row: " << _shell_row->getType() << " at " << _shell_row->getPos() << endl;
        cout << "\t col: " << _shell_col->getType() << " at " << _shell_col->getPos() << endl;*/
       
        // shell info, only lmax tells how far to go
        int _lmax_row = _shell_row->getLmax();
        int _lmax_col = _shell_col->getLmax();

        // set size of internal block for recursion
        int _nrows = this->getBlockSize( _lmax_row ); 
        int _ncols = this->getBlockSize( _lmax_col ); 
        
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
        
            // matrix for kinetic energies
            ub::matrix<double> kin = ub::zero_matrix<double>(_ntrafo_row,_nrows);
            //matrix for unnormalized overlap integrals
            ub::matrix<double> ol = ub::zero_matrix<double>(_ntrafo_row,_nrows);
        
            
            
            
            
            
            
            
            
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
        
        
        
        
        
                }
            }
        }
    }
}
    