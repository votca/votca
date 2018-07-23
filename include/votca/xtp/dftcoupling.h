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

#ifndef _VOTCA_XTP_DFTCOUPLING_H
#define	_VOTCA_XTP_DFTCOUPLING_H

#include <votca/xtp/couplingbase.h>

namespace votca { namespace xtp {


/**
* \brief Evaluates electronic coupling elements
*
* B. Baumeier, J. Kirkpatrick, D. Andrienko, 
* Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
*/

class DFTcoupling :CouplingBase
{
public:

    std::string  Identify() { return "dftcoupling"; }
    DFTcoupling():
    _degeneracy(0.0),_numberofstates(0){;}

    void CalculateCouplings(const Orbitals& orbitalsA, 
                               const Orbitals& orbitalsB, 
                               const Orbitals& orbitalsAB);
    
    void Initialize(tools::Property&);
    
    void Addoutput(tools::Property & type_summary,const Orbitals& orbitalsA, 
                               const Orbitals& orbitalsB);
    

private:
    
    
    void WriteToProperty(tools::Property& type_summary, const Orbitals& orbitalsA,
            const Orbitals& orbitalsB, int a, int b);
    double getCouplingElement( int levelA, int levelB,  
                               const Orbitals& orbitalsA,  
                               const Orbitals& orbitalsB
                                )const;
    
    
    std::pair<int,int> DetermineRangeOfStates(const Orbitals& orbital)const;
    
    Eigen::MatrixXd JAB;
    
    double _degeneracy;
    double _numberofstates;
    
    
    std::pair<int,int> Range_orbA;
    std::pair<int,int> Range_orbB;
    
  


};

}}

#endif	/* _VOTCA_XTP_DFTCOUPLING_H */


