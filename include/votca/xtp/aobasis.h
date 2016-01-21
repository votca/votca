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

#ifndef __XTP_AOBASIS__H
#define	__XTP_AOBASIS__H

#include <votca/tools/property.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/qmatom.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/basisset.h>


using namespace std;
using namespace votca::tools;

namespace votca { namespace xtp {
namespace ub = boost::numeric::ublas;

class AOShell;



/*
 * A collection of shells associated with a specific element  
 */
class AOBasis 
{   
public:
        AOBasis( ) { ; }
        ~AOBasis(); 
       // template< class T >
       //static void ReorderMOs(ub::matrix<T> &v, vector<int> const &order )  {
       //static void ReorderMOs(ub::matrix<T> &v, string start, string target )  { 
       void ReorderMOs(ub::matrix<double> &v, string start, string target ); 
       
      //template< class T >   
      void MultiplyMOs(ub::matrix<double> &v, vector<int> const &multiplier );

    // void AOBasisFill( BasisSet* bs , vector<Segment* > segments);
    void AOBasisFill( BasisSet* bs , vector<QMAtom* > segments, int fragbreak = -1);
    void ECPFill( BasisSet* bs , vector<QMAtom* > segments); 
    
    int NumFuncShell( string shell );
    int NumFuncShell_cartesian( string shell );
    int OffsetFuncShell( string shell );
    int OffsetFuncShell_cartesian( string shell );
    int AOBasisSize() {return _AOBasisSize; }
    
    typedef vector< AOShell* >::iterator AOShellIterator;
    AOShellIterator firstShell() { return _aoshells.begin(); }
    AOShellIterator lastShell(){ return _aoshells.end(); }

    // string getType() { return _type; }
    // int    getNumFunc() { return _numFunc ;}
        
    AOShell* getShell( AOShellIterator it ) { return (*it); }
    
    AOShell* getShell( int idx ){ return _aoshells[idx] ;}
    
    AOShell* addShell( string shellType, double shellScale, int shellFunc, int startIndex, int offset, vec pos, string name, int index ); 
  
    
    vector<AOShell*> getShells() { return _aoshells; }
    
 
   
   int _AOBasisSize;
   int _AOBasisFragA;
   int _AOBasisFragB;
   
   bool _is_stable;
   
    vector<AOShell*> _aoshells;

    // void getReorderVector( string& package, vector<int>& neworder );
    void getReorderVector( string& start, string& target, vector<int>& neworder );
    //void addReorderShell( string& package, string& shell, vector<int>& neworder );
    void addReorderShell( string& start, string& target, string& shell, vector<int>& neworder );
    //void getMultiplierVector( string& package, vector<int>& multiplier );
    void getMultiplierVector( string& start, string& target, vector<int>& multiplier );
    //void addMultiplierShell( string& package, string& shell, vector<int>& multiplier );  
    void addMultiplierShell( string& start, string& target, string& shell, vector<int>& multiplier );  
    
    
    void getTransformationCartToSpherical( string& package, ub::matrix<double>& _trafomatrix );
    void addTrafoCartShell(  AOShell* shell , ub::matrix_range< ub::matrix<double> >& _submatrix );
    
    int getMaxFunctions ( );
};


 
}}

#endif	/* AOBASIS_H */

