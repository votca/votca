/* 
 *            Copyright 2009-2017 The VOTCA Development Team
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
#include <boost/numeric/ublas/symmetric.hpp>




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
      
       void ReorderMOs(ub::matrix<double> &v,const std::string& start, const std::string& target ); 
       
       void ReorderMatrix(ub::symmetric_matrix<double> &v,const string& start,const string& target );
     
      

    void AOBasisFill( BasisSet* bs , std::vector<ctp::QMAtom* > segments, int fragbreak = -1);
    void ECPFill( BasisSet* bs , std::vector<ctp::QMAtom* > segments); 
    
    
    unsigned int AOBasisSize() const {return _AOBasisSize; }
    
    typedef std::vector< AOShell* >::const_iterator AOShellIterator;
    AOShellIterator firstShell() const{ return _aoshells.begin(); }
    AOShellIterator lastShell() const{ return _aoshells.end(); }

    ub::matrix<double> getTransformationCartToSpherical(const std::string& package);
    
        
    const AOShell* getShell( AOShellIterator it ) const{ return (*it); }
    
    const AOShell* getShell( int idx )const{ return _aoshells[idx] ;}
    
    AOShell* addShell( std::string shellType,int Lmax,int Lmin, double shellScale, int shellFunc, int startIndex, int offset, vec pos, std::string name, int index ); 
  
    
    const std::vector<AOShell*>& getShells() const{ return _aoshells; }
    
    unsigned getNumofShells() const{return _aoshells.size();}
   

   int _AOBasisFragA;
   int _AOBasisFragB;
   private:
       
       
  void MultiplyMOs(ub::matrix<double> &v, std::vector<int> const &multiplier );
   
    std::vector<AOShell*> _aoshells;

    vector<int> invertOrder(const vector<int>& order );
    
    std::vector<int> getReorderVector(const std::string& start,const std::string& target );
   
    void addReorderShell(const std::string& start,const std::string& target,const std::string& shell, std::vector<int>& neworder );
  
    std::vector<int> getMultiplierVector(const std::string& start,const std::string& target );
    
    void addMultiplierShell(const std::string& start,const std::string& target,const std::string& shell, std::vector<int>& multiplier );  
    
    
    
    void addTrafoCartShell( const AOShell* shell , ub::matrix_range< ub::matrix<double> >& _submatrix );
    
    int getMaxFunctions ( );
    
private:
    unsigned int _AOBasisSize;
    
};


 
}}

#endif	/* AOBASIS_H */

