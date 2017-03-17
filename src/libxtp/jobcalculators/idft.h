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

#ifndef _CALC_INTEGRALS_DFT_H
#define	_CALC_INTEGRALS_DFT_H

#include <votca/tools/property.h>

#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/overlap.h>

#include <boost/numeric/ublas/io.hpp>
#include <sys/stat.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    
/**
* \brief Density-functional-based electronic coupling elements
*
* Evaluates DFT-based electronic coupling elements for all conjugated
* segments from the neighbor list. Requires molecular orbitals of two monomers
* and a dimer in GAUSSIAN, NWChem, or TURBOMOLE format.
* 
*  B. Baumeier, J. Kirkpatrick, D. Andrienko, 
*  Phys. Chem. Chem. Phys., 12, 11103-11113, 2010
* 
* Callname: idft
*/

class IDFT : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    IDFT() {};
   ~IDFT() {};
   
    void    Initialize(tools::Property *options );
    
    string  Identify() { return "idft"; }
    
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *Thread);

    void WriteJobFile(Topology *top);

    void ReadJobFile( Topology *top );
    
/*  
    void    EvalPair(Topology *top, QMPair *pair, int slot);
    void    CleanUp();
*/

private:

    int                 _max_occupied_levels;
    int                 _max_unoccupied_levels;     
    int                 _trim_factor;
    
    string              _package;
    Property            _package_options; 
    
    // what to do
    bool                _do_input;
    bool                _do_run;
    bool                _do_parse;
    bool                _do_project;
    bool                _do_trim;
    bool                _do_extract;
    
    // what to write in the storage
    bool                _store_orbitals;
    bool                _store_overlap;
    bool                _store_integrals;
    
    double              _energy_difference;    
        
    string              _outParent;
        
    void ParseOptionsXML( tools::Property *opt);    
    
    
    /** 
     * \brief Guess for a dimer based on monomer orbitals
     * 
     * Given two monomer orbitals (A and B) constructs a guess for dimer
     * orbitals: | A 0 | and energies: [EA, EB]
     *           | 0 B |
     */
    void PrepareGuess(         Orbitals *_orbitalsA, 
                               Orbitals *_orbitalsB, 
                               Orbitals *_orbitalsAB, 
                               Logger *log = NULL );
    
    void LoadOrbitals(string file_name, Orbitals* orbitals, Logger *log = NULL );
        
};

}}
#endif	/* _CALC_INTEGRALS_DFT_H */
