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
#include "aomatrix.h"
//#include "linalg_tools.h"

using namespace std;
using namespace votca::tools;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    class TCMatrix {
    public:
    
        ub::vector< ub::matrix<double> >& matrix() { return this->_matrix ; }

        int get_mmin() { return this->mmin ;}
        int get_mmax() { return this->mmax ;}
        int get_nmin() { return this->nmin ;}
        int get_nmax() { return this->nmax ;}
        int get_mtot() { return this->mtotal ;}
        int get_ntot() { return this->ntotal ;}

        void set_mmin( int i ) { this->mmin = i ;}
        void set_mmax( int i ) { this->mmax = i ;}
        void set_nmin( int i ) { this->nmin = i ;}
        void set_nmax( int i ) { this->nmax = i ;}
        void set_mtot( int i ) { this->mtotal = i ;}
        void set_ntot( int i ) { this->ntotal = i ;}
        
        
        void Initialize ( int _basissize, int mmin, int mmax, int nmin, int nmax){

            // here as storage indices starting from zero
            set_mmin( mmin - 1 );
            set_mmax( mmax - 1 );
            set_nmin( nmin - 1 );
            set_nmax( nmax - 1 );
            set_mtot( mmax - mmin +1 );
            set_ntot( nmax - nmin +1 );

            /* let's try a different storage that is more convenient for
             later access 
            
            // vector has _basissize elements
            this->_matrix.resize( _basissize );
            
            // each element is a m-by-n matrix, initialize to zero
            for ( int i = 0; i < _basissize; i++){
                this->_matrix(i) = ub::zero_matrix<double>(mtotal,ntotal);
            }

             */

            
            // vector has mtotal elements
            this->_matrix.resize( this->get_mtot() );
            
            // each element is a gwabasis-by-n matrix, initialize to zero
            for ( int i = 0; i < this->get_mtot() ; i++){
                this->_matrix(i) = ub::zero_matrix<double>(_basissize,ntotal);
            }


            
        }
        
        void Fill( AOBasis& gwbasis, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals );

        // matrix print 
        void Print( string _ident);

        void Symmetrize( ub::matrix<double>& coulomb  );

        ub::matrix<double> matrixProd( int m, ub::matrix<double>& matrix);
        
    private:
        
        // store vector of matrices
        ub::vector< ub::matrix<double> > _matrix;
        
        // band summation indices
        int mmin;
        int mmax;
        int nmin;
        int nmax;
        int ntotal;
        int mtotal;
        
        
        int getBlockSize( int size );
        
        void getTrafo( ub::matrix<double>& _trafo, int _lmax, const double& _decay );
        
        void FillBlock(ub::vector_range< ub::vector< ub::matrix<double> > >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        
        void FillBlock(ub::vector< ub::matrix<double> >& _matrix,  AOShell* _shell, AOBasis& dftbasis, ub::matrix<double>& _dft_orbitals ) ;
        
        bool FillThreeCenterOLBlock(  ub::matrix<double> & _subvector, AOShell* _shell, AOShell* _shell_row, AOShell* _shell_col);


        
        
    };
    

}}

#endif	/* AOMATRIX_H */

