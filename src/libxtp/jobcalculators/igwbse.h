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

#ifndef _CALC_INTEGRALS_BSE_H
#define	_CALC_INTEGRALS_BSE_H

#include <votca/tools/property.h>

#include <votca/xtp/parallelxjobcalc.h>
#include <votca/xtp/orbitals.h>
//#include <votca/xtp/overlap.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/bsecoupling.h>
#include <boost/numeric/ublas/io.hpp>
#include <sys/stat.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    
/**
* \brief GWBSE-based exciton coupling elements
*
* Evaluates GWBSE-based exciton coupling elements for all conjugated
* segments from the neighbor list. Requires molecular orbitals of two monomers
* and a dimer in GAUSSIAN, NWChem, or TURBOMOLE (sometime) format.
* 
* Callname: igwbse
*/

class IGWBSE : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    IGWBSE() {};
   ~IGWBSE() {};
   
    void    Initialize(tools::Property *options );
    
    string  Identify() { return "igwbse"; }
    
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *Thread);

    void WriteJobFile(Topology *top);

    void ReadJobFile( Topology *top );
    
/*  
    void    EvalPair(Topology *top, QMPair *pair, int slot);
    void    CleanUp();
*/

private:

    int                 _number_excitons;
    //int                 _max_unoccupied_levels;     
    //int                 _trim_factor;
    
    string              _package;
    Property            _package_options; 
    Property            _gwbse_options; 
    
    string              _spintype;
    bool                _do_singlets;
    bool                _do_triplets;
    //GWBSE               _gwbse;
    //BSECoupling         _bsecoupling; 

    // what to do
    bool                _do_dft_input;
    bool                _do_dft_run;
    bool                _do_dft_parse;
    bool                _do_gwbse;
    bool                _do_coupling;
    bool                _do_trim;
    
    // what to write in the storage
    bool                _store_orbitals;
    bool                _store_overlap;
    bool                _store_integrals;
    bool                _store_ehint;
    
    double              _energy_difference;    
        
    string              _outParent;
    
    
    // parsing options
    std::map<std::string, int> _singlet_levels;
    std::map<std::string, int> _triplet_levels;

    void ParseOptionsXML( tools::Property *opt);    
    
    std::map<std::string, int> FillParseMaps(string Mapstring);
    
    
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
#endif	/* _CALC_INTEGRALS_BSE_H */
