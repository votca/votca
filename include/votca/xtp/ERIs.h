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


#include <votca/xtp/threecenter.h>
#include <votca/xtp/fourcenter.h>

namespace votca { namespace xtp {
   
/**
* \brief Takes a density matrix and and an auxillary basis set and calculates the electron repulsion integrals. 
*
* 
* 
*/    
    class ERIs{    
        
    public:
      
        
        
        void Initialize(AOBasis &_dftbasis, AOBasis &_auxbasis, const Eigen::MatrixXd &inverse_Coulomb);
        void Initialize_4c_small_molecule(AOBasis &_dftbasis); 
      
        const Eigen::MatrixXd& getEXX() const{return _EXXs;}
        
        const Eigen::MatrixXd& getERIs() const{return _ERIs;}
        double& getERIsenergy(){return _ERIsenergy;}
        double& getEXXsenergy(){return _EXXsenergy;}
        
        void CalculateERIs(const Eigen::MatrixXd &DMAT);
        void CalculateEXX(const Eigen::MatrixXd &DMAT);
        void CalculateEXX(const Eigen::Block<Eigen::MatrixXd>& occMos,const Eigen::MatrixXd &DMAT);
        void CalculateERIs_4c_small_molecule(const Eigen::MatrixXd &DMAT); 
        void CalculateEXX_4c_small_molecule(const Eigen::MatrixXd &DMAT);
        
        int getSize1(){return _ERIs.rows();}
        int getSize2(){return _ERIs.cols();}
        
        void printERIs();
        
    private:
        
        
        TCMatrix_dft _threecenter;
        FCMatrix _fourcenter; 
       
        Eigen::MatrixXd _ERIs;
        Eigen::MatrixXd _EXXs;
        double _ERIsenergy;
        double _EXXsenergy;
        void CalculateEnergy(const Eigen::MatrixXd &DMAT);
        void CalculateEXXEnergy(const Eigen::MatrixXd &DMAT);
    };
    
    
    
    

}}

#endif	/* ERIS_H */

