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

#include <votca/ctp/qmatom.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/math/constants/constants.hpp>
#include <votca/xtp/basisset.h>



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
      
       void ReorderMOs(ub::matrix<double> &v, std::string start, std::string target ); 
     
      void MultiplyMOs(ub::matrix<double> &v, std::vector<int> const &multiplier );

    // void AOBasisFill( BasisSet* bs , std::vector<Segment* > segments);
    void AOBasisFill( BasisSet* bs , std::vector<ctp::QMAtom* > segments, int fragbreak = -1);
    void ECPFill( BasisSet* bs , std::vector<ctp::QMAtom* > segments); 
    
    
    int AOBasisSize() {return _AOBasisSize; }
    
    typedef std::vector< AOShell* >::iterator AOShellIterator;
    AOShellIterator firstShell() { return _aoshells.begin(); }
    AOShellIterator lastShell(){ return _aoshells.end(); }

    
    
        
    AOShell* getShell( AOShellIterator it ) { return (*it); }
    
    AOShell* getShell( int idx ){ return _aoshells[idx] ;}
    
    AOShell* addShell( std::string shellType,int Lmax,int Lmin, double shellScale, int shellFunc, int startIndex, int offset, vec pos, std::string name, int index ); 
  
    
    std::vector<AOShell*> getShells() { return _aoshells; }
    
 
   
   int _AOBasisSize;
   int _AOBasisFragA;
   int _AOBasisFragB;
   
   bool _is_stable;
   
    std::vector<AOShell*> _aoshells;

  
    void getReorderVector( std::string& start, std::string& target, std::vector<int>& neworder );
   
    void addReorderShell( std::string& start, std::string& target, std::string& shell, std::vector<int>& neworder );
  
    void getMultiplierVector( std::string& start, std::string& target, std::vector<int>& multiplier );
    
    void addMultiplierShell( std::string& start, std::string& target, std::string& shell, std::vector<int>& multiplier );  
    
    
    void getTransformationCartToSpherical( std::string& package, ub::matrix<double>& _trafomatrix );
    void addTrafoCartShell(  AOShell* shell , ub::matrix_range< ub::matrix<double> >& _submatrix );
    
    int getMaxFunctions ( );
};


 
}}

#endif	/* AOBASIS_H */

