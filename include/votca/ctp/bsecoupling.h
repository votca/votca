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

#include <votca/ctp/orbitals.h>
#include <votca/ctp/logger.h>

#ifndef _VOTCA_CTP_BSECOUPLING_H
#define	_VOTCA_CTP_BSECOUPLING_H

namespace votca { namespace ctp {


/**
* \brief Evaluates electronic coupling elements
*
* B. Baumeier, J. Kirkpatrick, D. Andrienko, 
* Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
*/

class BSECoupling 
{
public:

    BSECoupling() {};
   ~BSECoupling() {};


    bool CalculateCouplings(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB, 
                               boost::numeric::ublas::matrix<float>* _JAB_singlet, boost::numeric::ublas::matrix<float>* _JAB_triplet, string _type);  
    
    
    bool ProjectExcitons(boost::numeric::ublas::matrix<float>& _kap,
                         boost::numeric::ublas::matrix<float>& _kbp, 
                         boost::numeric::ublas::matrix<float>& _bseA, 
                         boost::numeric::ublas::matrix<float>& _bseB, 
                         boost::numeric::ublas::matrix<float>& _H, 
                         boost::numeric::ublas::matrix<float>& _J );

     
    float getCouplingElement( int levelA, int levelB,  
                               Orbitals* _orbitalsA,  
                               Orbitals* _orbitalsB, 
                               boost::numeric::ublas::matrix<float>* _JAB,
                               double _energy_difference = 0
                                );
    
    void setLogger( Logger* pLog ) { _pLog = pLog; }
    
private:
    
    Logger *_pLog;
    
    void SQRTOverlap(boost::numeric::ublas::symmetric_matrix<double> &S, 
                     boost::numeric::ublas::matrix<double> &Sm2);


};

}}

#endif	/* _VOTCA_CTP_BSECOUPLING_H */


