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

#ifndef __CTP_AOMATRIX__H
#define	__CTP_AOMATRIX__H

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

namespace Cartesian {
        enum Index {
                s, x, y, z,  xy, xz, yz, xx, yy, zz
                };
}

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
    
    
    
    /* "superclass" AOSuperMatrix contains all common functionality for
     * atomic orbital matrix types
     */
        class AOSuperMatrix{
    public:
        
        int getBlockSize( int size );
        
        void getTrafo( ub::matrix<double>& _trafo, int _lmax, const double& _decay );
        void getTrafo( ub::matrix<double>& _trafo, int _lmax, const double& _decay , std::vector<double> contractions);
        
        void PrintIndexToFunction( AOBasis* aobasis);
        
        
    };
    
    
    // base class for 1D atomic orbital matrix types (overlap, Coulomb, ESP)
    class AOMatrix : public AOSuperMatrix {
    public:
        ub::matrix<double> _aomatrix; 
        ub::vector<double> _gridpoint;
        
        void Initialize( int size ) {
            this->_aomatrix = ub::zero_matrix<double>(size,size);
        }
        
        void Fill( AOBasis* aobasis, ub::vector<double> r = ub::zero_vector<double>(3) );
        
        // matrix print 
        void Print( string _ident);
        // integrate F
        void XIntegrate( vector<double>& _FmT, const double& _T );
        // block fill prototype
        virtual void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, bool _raw = false) {} ;

        // ~AOMatrix(){};

    };
    
    
    
    /* base class class for 3D atomic orbital matrix types 
     * (in principle, we could make nD and 1D and 3D are just special types)
     */
    class AOMatrix3D : public AOSuperMatrix {
    public:
        std::vector<ub::matrix<double> > _aomatrix; 
        
        void Initialize( int size ) {
            _aomatrix.resize(3);
            for (int i = 0; i < 3 ; i++){
              _aomatrix[ i ] = ub::zero_matrix<double>(size,size);
            }
        }



        // matrix print 
        void Print( string _ident);

        
        void Fill( AOBasis* aobasis );

        // block fill prototype
        virtual void FillBlock(std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix, AOShell* _shell_row, AOShell* _shell_col, bool _raw = false) {} ;

        
        void Cleanup();
        
      //  ~AOMatrix3D();
        
    };
    
    
    
    /* derived class for atomic orbital gradient matrices, required for
     * momentum transition dipoles
     */
    class AOMomentum : public AOMatrix3D { 
        
        //block fill for gradient/momentum operator, implementation in aomomentum.cc
        void FillBlock( std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix, AOShell* _shell_row, AOShell* _shell_col , bool _raw = false);
        
        
    };
    
    
    
    /* derived class for atomic orbital electrical dipole matrices, required for
     * electical transition dipoles
     */
    class AODipole : public AOMatrix3D { 
        
        //block fill for gradient/momentum operator, implementation in aomomentum.cc
        void FillBlock( std::vector< ub::matrix_range< ub::matrix<double> > >& _matrix, AOShell* _shell_row, AOShell* _shell_col , bool _raw = false);
        
        
    };
    
    
    // derived class for atomic orbital nuclear potential
    class AOESP : public AOMatrix{
    public:
        //block fill for overlap, implementation in aoesp.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, bool _raw = false );
        //void Print();
        
        // ~AOESP();
        
        
    };
    
    
    // derived class for atomic orbital overlap
    class AOOverlap : public AOMatrix{
    public:
        //block fill for overlap, implementation in aooverlap.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, bool _raw = false );
        //void Print();
        
	//        ~AOOverlap();
        
    };
    
    // inline AOOverlap::~AOOverlap(){
        
    //_aomatrix.clear();
    //_aomatrix.resize(0,0);
        
    //}
    
    //derived class for atomic orbital Coulomb interaction
    class AOCoulomb : public AOMatrix{
    public:
        int getExtraBlockSize( int lmax_row, int lmax_col  );
        void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col, bool _raw = false);
        void Symmetrize( AOOverlap& _overlap , AOBasis& _basis, AOOverlap& _overlap_inverse , AOOverlap& _gwoverlap_cholesky_inverse );
        
  
    };
}}

#endif	/* AOMATRIX_H */

