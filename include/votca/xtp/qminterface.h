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

#ifndef __QMMINTERFACE__H
#define	__QMMINTERFACE__H

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/tools/linalg.h>

#include <votca/ctp/xjob.h>
#include <votca/ctp/xinductor.h>

// add gwbse header for excited state support
#include <votca/xtp/gwbse.h>
#include <votca/xtp/qmpackagefactory.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/gdma.h>


namespace votca { namespace xtp {

    
// ========================================================================== //
// QM-MM INTERFACE CLASS - CONVERTS BETWEEN QMATOMS <> POLAR OBJECTS          //
// ========================================================================== //
    
class QMMInterface
{
public:
    
    QMMInterface() { _polar_table = ctp::POLAR_TABLE(); };
   ~QMMInterface() {};
    
    // CONVERSION QM -> MM
    ctp::APolarSite *Convert(ctp::QMAtom *atm, int id = -1);
    
    ctp::PolarSeg *Convert(std::vector<ctp::QMAtom*> &atms);
  
    
private:
    
    // Allocates polarizabilities in A**3 to element types
    std::map<std::string,double> _polar_table;
    
};




}}

#endif
