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

#include <votca/xtp/orbitals.h>
#include <votca/xtp/logger.h>

#ifndef _VOTCA_XTP_OVERLAP_H
#define	_VOTCA_XTP_OVERLAP_H

namespace votca { namespace xtp {


/**
* \brief Evaluates electronic coupling elements
*
* B. Baumeier, J. Kirkpatrick, D. Andrienko, 
* Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
*/

class Overlap 
{
public:

    Overlap() {};
   ~Overlap() {};


    bool CalculateIntegralsOptimized(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB, 
                               boost::numeric::ublas::matrix<double>* _JAB);  

    bool CalculateIntegrals(   Orbitals* _orbitalsA, 
                               Orbitals* _orbitalsB, 
                               Orbitals* _orbitalsAB, 
                               boost::numeric::ublas::matrix<double>* _JAB);
    
    double getCouplingElement( int levelA, int levelB,  
                               Orbitals* _orbitalsA,  
                               Orbitals* _orbitalsB, 
                               boost::numeric::ublas::matrix<double>* _JAB,
                               double _energy_difference = 0
                                );
    
    void setLogger( Logger* pLog ) { _pLog = pLog; }
    
private:
    
    Logger *_pLog;
    
    void SQRTOverlap(boost::numeric::ublas::symmetric_matrix<double> &S, 
                     boost::numeric::ublas::matrix<double> &Sm2);


};

}}

#endif	/* _VOTCA_XTP_OVERLAP_H */


