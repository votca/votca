/* 
 *            Copyright 2009-2018 The VOTCA Development Team
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


#include <boost/math/constants/constants.hpp>
#include <votca/xtp/basisset.h>
#include <votca/tools/vec.h>
#include <votca/xtp/eigen.h>



namespace votca { namespace xtp {

class AOShell;
class QMAtom;


/**
 * \brief Container to hold Basisfunctions for all atoms 
 * 
 * It is constructed from a vector of QMAtoms and a BasisSet.
 */
class AOBasis 
{   
public:
    
    ~AOBasis();//has to be declared, deletes std::vector<*Shell>_aoshells
    void ReorderMOs(Eigen::MatrixXd &v,const std::string& start, const std::string& target ); 
       
    void ReorderMatrix(Eigen::MatrixXd &v,const std::string& start,const std::string& target );

    void AOBasisFill( const BasisSet& bs , std::vector<QMAtom* > segments, int fragbreak = -1);
    void ECPFill( const BasisSet& bs , std::vector<QMAtom* > segments); 
    
    unsigned AOBasisSize() const {return _AOBasisSize; }
    
    typedef std::vector< AOShell* >::const_iterator AOShellIterator;
    AOShellIterator firstShell() const{ return _aoshells.begin(); }
    AOShellIterator lastShell() const{ return _aoshells.end(); }

    Eigen::MatrixXd getTransformationCartToSpherical(const std::string& package);
    
    const AOShell* getShell( int idx )const{ return _aoshells[idx] ;}

    const std::vector<AOShell*>& getShells() const{ return _aoshells; }
    
    const std::vector<const AOShell*> getShellsperAtom(int AtomId)const;
    
    int getFuncperAtom(int AtomId) const;
    
    unsigned getNumofShells() const{return _aoshells.size();}
   
    int getAOBasisFragA() const{return _AOBasisFragA;}
    
   int getAOBasisFragB() const{return _AOBasisFragB;}
  

private:
    
  AOShell* addShell( const Shell& shell, const QMAtom& atom, int startIndex );  
  
  AOShell* addECPShell( const Shell& shell, const QMAtom& atom, int startIndex,bool nonlocal);  
    
       
  void MultiplyMOs(Eigen::MatrixXd &v, std::vector<int> const &multiplier );
   
    std::vector<AOShell*> _aoshells;

    std::vector<int> invertOrder(const std::vector<int>& order );
    
    std::vector<int> getReorderVector(const std::string& start,const std::string& target );
   
    void addReorderShell(const std::string& start,const std::string& target,const std::string& shell, std::vector<int>& neworder );
  
    std::vector<int> getMultiplierVector(const std::string& start,const std::string& target );
    
    void addMultiplierShell(const std::string& start,const std::string& target,const std::string& shell, std::vector<int>& multiplier );  
  
    void addTrafoCartShell( const AOShell* shell , Eigen::Block<Eigen::MatrixXd>& _submatrix );
    
    
   int _AOBasisFragA;
   int _AOBasisFragB;
    unsigned int _AOBasisSize;
    
};


 
}}

#endif	/* AOBASIS_H */

