/* 
 *            Copyright 2009-2012 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICEN_olE-2.0
 *
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/ctp/aomatrix.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/threecenters.h>
#include <votca/ctp/ERIs.h>
#include <string>
#include <vector>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <votca/tools/linalg.h>
#include <omp.h>

using namespace std;
using namespace votca::tools;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;
        
        
        
        void ERIs::Initialize (AOBasis &_dftbasis, AOBasis &_auxbasis,  AOOverlap &_auxAOverlap, AOCoulomb &_auxAOcoulomb){

           
           
            _threecenter.Fill( _auxbasis, _dftbasis );
            ub::matrix<double> _inverse=ub::zero_matrix<double>( _auxAOverlap.Dimension(), _auxAOverlap.Dimension());
            linalg_invert( _auxAOverlap.Matrix() , _inverse);
            ub::matrix<double> _temp=ub::prod(_auxAOcoulomb.Matrix(),_inverse);
            _Vcoulomb=ub::prod(_inverse,_temp);
            
        
        
        
        
        
        
        
        }
        
        
        void ERIs::CalculateERIs (ub::matrix<double> &DMAT){
           for ( int i=0,i < _three) 
            
            
        }
        
        
        
    }}