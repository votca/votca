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

#include <votca/xtp/orbitals.h>
#include <votca/ctp/logger.h>



namespace votca { namespace xtp {


/**
* \brief Evaluates electronic coupling elements
*
* B. Baumeier, J. Kirkpatrick, D. Andrienko, 
* Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
*/

class DFTcoupling 
{
public:

    DFTcoupling() {};
   ~DFTcoupling() {};



    Eigen::MatrixXd CalculateIntegrals(Orbitals& _orbitalsA, 
                               Orbitals& _orbitalsB, 
                               Orbitals& _orbitalsAB);
    
    double getCouplingElement( int levelA, int levelB,  
                               Orbitals& _orbitalsA,  
                               Orbitals& _orbitalsB, 
                               Eigen::MatrixXd* _JAB,
                               double _energy_difference = 0
                                );
    
    void setLogger( ctp::Logger* pLog ) { _pLog = pLog; }
    
private:
    
    ctp::Logger *_pLog;
    
  


};

}}

#endif	/* _VOTCA_XTP_DFTCOUPLING_H */


