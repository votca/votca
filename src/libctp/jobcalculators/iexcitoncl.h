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

#ifndef _CALC_COUPLING_EXCL_H
#define	_CALC_COUPLING_EXCL_H

#include <votca/tools/property.h>

#include <votca/ctp/parallelxjobcalc.h>
#include <boost/numeric/ublas/io.hpp>
#include <sys/stat.h>
#include <boost/filesystem.hpp>

namespace votca { namespace ctp {
    
/**
* \brief Evaluates Transition Charge distributions classically
*
* Evaluates the electrostatic classical coupling between molecules in 
* their excited states.
* and a dimer in GAUSSIAN, NWChem, or TURBOMOLE format.
* 
*  B. Baumeier, J. Kirkpatrick, D. Andrienko, 
*  Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
* Callname: idft
*/

class IEXCITON : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    IEXCITON() {};
   ~IEXCITON() {};
   
    void    Initialize(tools::Property *options );
    
    string  Identify() { return "iexcitoncl"; }
    
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *Thread);

    void WriteJobFile(Topology *top);


    

private:

    
    
    // what to do
    bool _do_polar;
    
    // what to write in the storage

        

        
    void ParseOptionsXML( tools::Property *opt);    
    
    
    /** 
     * \brief Guess for a dimer based on monomer orbitals
     * 
     * Given two monomer orbitals (A and B) constructs a guess for dimer
     * orbitals: | A 0 | and energies: [EA, EB]
     *           | 0 B |
     */

    
 
        
};

}}
#endif	/* _CALC_INTEGRALS_DFT_H */
