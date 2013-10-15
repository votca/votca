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

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    class AOMatrix {
    public:
        ub::matrix<double> _aomatrix; 
        
        void Initialize( int size ) {
            this->_aomatrix = ub::zero_matrix<double>(size,size);
        }
        
        void Fill( AOBasis* aobasis );
        
        int getBlockSize( int size );
        
        void getTrafo( ub::matrix<double>& _trafo, int _lmax, const double& _decay );

        /* Matrix inversion routine.
           Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
        //template<class T>
        bool InvertMatrix(const ub::matrix<double>& input, ub::matrix<double>& inverse);

        
        void PrintIndexToFunction( AOBasis* aobasis);
        
        // matrix print 
        void Print( string _ident);
        
        // block fill prototype
        virtual void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col) {} ;
        
    };
    
    class AOOverlap : public AOMatrix{
    public:
        //block fill for overlap, implementation in aomatrix.cc
        void FillBlock( ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col );
        //void Print();
    };
    
    
    class AOCoulomb : public AOMatrix{
    public:
        int getExtraBlockSize( int lmax_row, int lmax_col  );
        void FillBlock(ub::matrix_range< ub::matrix<double> >& _matrix, AOShell* _shell_row, AOShell* _shell_col);
        void XIntegrate( vector<double>& _FmT, const double& _T );
        void Symmetrize( AOOverlap& _overlap , AOBasis& _basis);
        
       // some helpers
       vector<double> _wmp ;
       vector<double> _wmq ;
       double _distsq;
  
       double _faka  ;
       double _faka2 ;
       double _faka3 ;
       double _fakaca ;
       double _fakaac ;
       double _fakac  ;
       double _fakac2 ;
       double _fakac3 ;
       double _fakac4 ;
       double _fakc   ;
       double _fakc2  ;
       double _fakc3  ;
       double _fakca  ;
       double _fakca2 ;
       double _fakca3 ;
       double _fakca4 ;
        
        
        bool has_s_p;
        bool has_p_s;
        bool has_s_p_1;
        bool has_p_s_1;
        bool has_p_p;
        bool has_d_s;
        bool has_s_d;
        bool has_s_p_2;
        bool has_p_p_1;
        bool has_d_p;
        bool has_p_d;
        bool has_s_d_1;
        bool has_s_p_3;
        bool has_s_d_2;
        bool has_p_d_1;
        bool has_d_d;
        bool has_p_s_2;
        bool has_d_s_1;
        bool has_f_s;
        bool has_s_f;
        bool has_p_p_2;
        bool has_d_p_1;
        bool has_f_p;
        bool has_s_f_1;
        bool has_p_f;
        bool has_d_f;
        bool has_s_p_4;
        bool has_s_d_3;
        bool has_s_f_2;
        bool has_p_f_1;
        bool has_f_d;
        bool has_p_d_2;
        bool has_d_d_1;
        bool has_s_p_5;
        bool has_s_d_4;
        bool has_s_f_3;
        bool has_p_f_2;
        bool has_d_f_1;
        bool has_f_f;
        
        
        bool calc_s_p( boost::multi_array<double, 3>& _cou );
        bool calc_p_s( boost::multi_array<double, 3>& _cou );
        bool calc_s_p_1( boost::multi_array<double, 3>& _cou );
        bool calc_p_s_1( boost::multi_array<double, 3>& _cou );
        bool calc_p_p( boost::multi_array<double, 3>& _cou );
        bool calc_d_s( boost::multi_array<double, 3>& _cou );
        bool calc_s_d( boost::multi_array<double, 3>& _cou );
        bool calc_s_p_2( boost::multi_array<double, 3>& _cou );  
        bool calc_p_p_1( boost::multi_array<double, 3>& _cou );
        bool calc_d_p( boost::multi_array<double, 3>& _cou );  
        bool calc_p_d( boost::multi_array<double, 3>& _cou );
        bool calc_s_d_1( boost::multi_array<double, 3>& _cou ); 
        bool calc_s_p_3( boost::multi_array<double, 3>& _cou ); 
        bool calc_s_d_2( boost::multi_array<double, 3>& _cou );
        bool calc_p_d_1( boost::multi_array<double, 3>& _cou ); 
        bool calc_d_d( boost::multi_array<double, 3>& _cou );
        bool calc_p_s_2( boost::multi_array<double, 3>& _cou );
        bool calc_d_s_1( boost::multi_array<double, 3>& _cou );
        bool calc_f_s( boost::multi_array<double, 3>& _cou );
        bool calc_s_f( boost::multi_array<double, 3>& _cou );
        bool calc_d_p_1( boost::multi_array<double, 3>& _cou );
        bool calc_p_p_2( boost::multi_array<double, 3>& _cou );
        bool calc_f_p( boost::multi_array<double, 3>& _cou );
        bool calc_s_f_1( boost::multi_array<double, 3>& _cou );
        bool calc_p_f( boost::multi_array<double, 3>& _cou );
        bool calc_d_f( boost::multi_array<double, 3>& _cou );
        bool calc_s_p_4( boost::multi_array<double, 3>& _cou );
        bool calc_s_d_3( boost::multi_array<double, 3>& _cou );
        bool calc_s_f_2( boost::multi_array<double, 3>& _cou );
        bool calc_p_f_1( boost::multi_array<double, 3>& _cou );
        bool calc_p_d_2( boost::multi_array<double, 3>& _cou );
        bool calc_d_d_1( boost::multi_array<double, 3>& _cou );
        bool calc_f_d( boost::multi_array<double, 3>& _cou );
        bool calc_s_p_5( boost::multi_array<double, 3>& _cou );
        bool calc_s_d_4( boost::multi_array<double, 3>& _cou );
        bool calc_s_f_3( boost::multi_array<double, 3>& _cou );
        bool calc_p_f_2( boost::multi_array<double, 3>& _cou );
        bool calc_d_f_1( boost::multi_array<double, 3>& _cou );
        bool calc_f_f( boost::multi_array<double, 3>& _cou );
    };
}}

#endif	/* AOMATRIX_H */

