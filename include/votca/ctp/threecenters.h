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

#ifndef __CTP_THREECENTERS__H
#define	__CTP_THREECENTERS__H

#include <votca/ctp/aobasis.h>
#include <votca/ctp/segment.h>
#include <string>
#include <map>
#include <vector>
#include <votca/tools/property.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>
#include "basisset.h"
//#include "linalg_tools.h"

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    class TCMatrix {
    public:
        ub::matrix<double> _aomatrix; 

        typedef boost::multi_array<double, 3> ma_type;
        typedef boost::multi_array_types::extent_range range;
        typedef ma_type::index index;
        ma_type::extent_gen extents;
        ma_type _array;
        
        void Initialize( int _basissize, int _mmin, int _mmax, int _nmin, int _nmax ) {
            
            this->_array.resize( extents[ _basissize ][ range( _mmin-1, _mmax )  ][ range(_nmin-1, _nmax )]);
    
            // --- initialize it to zero
            for( index i = this->_array.index_bases()[0] ; i < this->_array.shape()[0]  ; ++i) {
               for( index j =  this->_array.index_bases()[1] ; j < this->_array.shape()[1] ; ++j){
                  for( index k = this->_array.index_bases()[2] ; k < this->_array.shape()[2] ; ++k){
                      this->_array[i][j][k] = 0.0;
                  }
               }
            } 
        }
        
        void Fill( AOBasis& gwbasis, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals );
        
        int getBlockSize( int size );
        
        void getTrafo( ub::matrix<double>& _trafo, int _lmax, const double& _decay );
        
        void PrintIndexToFunction( AOBasis* aobasis);
        
        // matrix print 
        void Print( string _ident);
        
        // block fill prototype
        virtual void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col) {} ;
        
    };
    

}}

#endif	/* AOMATRIX_H */

