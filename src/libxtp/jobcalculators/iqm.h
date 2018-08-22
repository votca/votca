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

#ifndef _CALC_XTP_IQM_H
#define	_CALC_XTP_IQM_H

#include <votca/tools/property.h>

#include <votca/ctp/parallelxjobcalc.h>
#include <votca/xtp/orbitals.h>
#include <votca/xtp/dftcoupling.h>
#include <votca/xtp/gwbse.h>
#include <votca/xtp/bsecoupling.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>

namespace votca { namespace xtp {
    
/**
* \brief DFT & GWBSE-based coupling elements
*
* Evaluates DFT & GWBSE-based coupling elements for all conjugated
* segments from the neighbor list. Requires molecular orbitals of two monomers
* and a dimer in GAUSSIAN, NWChem, or ORCAformat.
* 
* Callname: iqm
*/

class IQM : public ctp::ParallelXJobCalc< vector<ctp::Job*>, ctp::Job*, ctp::Job::JobResult >
{
public:
   
    void    Initialize(tools::Property *options ); 
    string  Identify() { return "iqm"; }   
    ctp::Job::JobResult EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *Thread);  
    void WriteJobFile(ctp::Topology *top);
    void ReadJobFile( ctp::Topology *top );

private:
    
    double GetBSECouplingFromProp(tools::Property& bseprop, int stateA, int stateB);
    double GetDFTCouplingFromProp(tools::Property& dftprop, int stateA, int stateB);
    void SetJobToFailed(ctp::Job::JobResult& jres, ctp::Logger* pLog, const string& errormessage);
    void addLinkers(std::vector< ctp::Segment* > &segments, ctp::Topology *top);
    bool isLinker(const std::string& name);
    void WriteCoordinatesToOrbitalsPBC(ctp::QMPair& pair, Orbitals& orbitals);
    void ParseOptionsXML( tools::Property &opt);    
    std::map<std::string, int> FillParseMaps(const string& Mapstring);
    
    string              _package;
    Property            _dftpackage_options; 
    Property            _gwbse_options; 
    Property            _bsecoupling_options; 
    Property            _dftcoupling_options; 

    // what to do
    bool                _do_dft_input;
    bool                _do_dft_run;
    bool                _do_dft_parse;
    bool                _do_dftcoupling;
    bool                _do_gwbse;
    bool                _do_bsecoupling;
    
    std::vector< std::string > _linker_names;
    
    // what to write in the storage
    bool                _store_dft;
    bool                _store_singlets;
    bool                _store_triplets;
    bool                _store_ehint;
    bool                _write_orbfile;
    
    //double              _energy_difference;    
        
    string              _outParent;
      
    // parsing options
    std::map<std::string, QMState> _singlet_levels;
    std::map<std::string, QMState> _triplet_levels;
    
    std::map<std::string, QMState> _hole_levels;
    std::map<std::string, QMState> _electron_levels;

        
};

}}
#endif	
