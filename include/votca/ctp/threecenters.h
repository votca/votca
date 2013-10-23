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

              
        // likely better: go via vector of matrices
        ub::vector< ub::matrix<double> > _matrix;
        
        // band summation indices
        int mmin;
        int mmax;
        int nmin;
        int nmax;
        int ntotal;
        int mtotal;
        
        
        void Initialize ( int _basissize, int mmin, int mmax, int nmin, int nmax){

            // here as stoarage indices starting from zero
            this->mmin   = mmin -1 ;
            this->mmax   = mmax -1 ;
            this->nmin   = nmin -1 ;
            this->nmax   = nmax -1 ;
            this->mtotal = mmax - mmin +1;
            this->ntotal = nmax - nmin +1;
            
            // vector has _basissize elements
            this->_matrix.resize( _basissize );
            
            // each element is a m-by-n matrix, initialize to zero
            for ( int i = 0; i < _basissize; i++){
                this->_matrix(i) = ub::zero_matrix<double>(mtotal,ntotal);
            }
            
        }
        
        /* void Initialize( int _basissize, int _mmin, int _mmax, int _nmin, int _nmax ) {
            
            this->_array.resize( extents[ _basissize ][ range( _mmin-1, _mmax )  ][ range(_nmin-1, _nmax )]);
    
            // --- initialize it to zero
            for( index i = this->_array.index_bases()[0] ; i < this->_array.shape()[0]  ; ++i) {
               for( index j =  this->_array.index_bases()[1] ; j < this->_array.shape()[1] ; ++j){
                  for( index k = this->_array.index_bases()[2] ; k < this->_array.shape()[2] ; ++k){
                      this->_array[i][j][k] = 0.0;
                  }
               }
            } 
        }*/
        
        void Fill( AOBasis& gwbasis, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals );
        
        int getBlockSize( int size );
        
        void getTrafo( ub::matrix<double>& _trafo, int _lmax, const double& _decay );
        
        void PrintIndexToFunction( AOBasis* aobasis);
        
        // matrix print 
        void Print( string _ident);
        
        // block fill prototype
        //void FillBlock(ub::vector_range< ub::vector< ub::matrix<double> > >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        void FillBlock(ub::vector_range< ub::vector< ub::matrix<double> > >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        
        //bool FillThreeCenterOLBlock( ub::vector< ub::matrix<double> >& _subvector, AOShell* _shell, AOShell* _shell_row, AOShell* _shell_col);
        bool FillThreeCenterOLBlock(  ub::matrix<double> & _subvector, AOShell* _shell, AOShell* _shell_row, AOShell* _shell_col);
    };
    

}}

#endif	/* AOMATRIX_H */

