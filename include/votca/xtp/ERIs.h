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

#ifndef _VOTCA_XTP_ERIS_H
#define	_VOTCA_XTP_ERIS_H



#include <votca/xtp/threecenters.h>
#include <votca/tools/linalg.h>




namespace votca { namespace xtp {
    namespace ub = boost::numeric::ublas;
/**
* \brief Takes a density matrix and and an auxillary basis set and calculates the electron repulsion integrals. 
*
* 
* 
*/    
    class ERIs{    
        
    public:
      
        
        
        void Initialize(AOBasis &_dftbasis, AOBasis &_auxbasis, const ub::matrix<double> &inverse_Coulomb);
        void Initialize_4c_small_molecule(AOBasis &_dftbasis); ///////////
      
        const ub::matrix<double>& getEXX() const{return _EXXs;}
        
        const ub::matrix<double>& getERIs() const{return _ERIs;}
        double& getERIsenergy(){return _ERIsenergy;}
        
        void CalculateERIs(const ub::matrix<double> &DMAT);
        void CalculateERIs_4c_small_molecule(const ub::matrix<double> &DMAT); ///////////////////////////////////////
        void CalculateEXX_4c_small_molecule(const ub::matrix<double> &DMAT);
        
        int getSize1(){return _ERIs.size1();}
        int getSize2(){return _ERIs.size2();}
        
        void printERIs();
        
    private:
        ub::matrix<double> _inverse_Coulomb;
        
        TCMatrix_dft _threecenter;
        FCMatrix_dft _fourcenter; ////////////////////////
       
        ub::matrix<double> _ERIs;
        ub::matrix<double> _EXXs;
        double _ERIsenergy;
        double _EXXenergy;
        void CalculateEnergy(const ub::vector<double> &dmatasarray);
        void CalculateEXXEnergy(const ub::vector<double> &dmatasarray);
    };
    
    
    
    

}}

#endif	/* ERIS_H */

