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

#include <votca/xtp/parallelxjobcalc.h>
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

class IQM : public ParallelXJobCalc< std::vector<Job*>, Job*, Job::JobResult >
{
public:
   
    void    Initialize(tools::Property *options ); 
    std::string  Identify() { return "iqm"; }   
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *Thread);  
    void WriteJobFile(Topology *top);
    void ReadJobFile( Topology *top );

private:
    
    double GetBSECouplingFromProp(tools::Property& bseprop,const QMState& stateA,const QMState& stateB);
    double GetDFTCouplingFromProp(tools::Property& dftprop, int stateA, int stateB);
    void SetJobToFailed(Job::JobResult& jres, Logger* pLog, const std::string& errormessage);
    void WriteLoggerToFile(const std::string& logfile, Logger& logger);
    void addLinkers(std::vector< Segment* > &segments, Topology *top);
    bool isLinker(const std::string& name);
    void WriteCoordinatesToOrbitalsPBC(QMPair& pair, Orbitals& orbitals);
    void ParseOptionsXML(tools::Property &opt);    
    std::map<std::string, QMState> FillParseMaps(const std::string& Mapstring);
    
    QMState GetElementFromMap(const std::map<std::string, QMState>& elementmap,const std::string& elementname )const;
    
    std::string              _package;
    tools::Property            _dftpackage_options; 
    tools::Property            _gwbse_options; 
    tools::Property            _bsecoupling_options; 
    tools::Property            _dftcoupling_options; 

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
      
    // parsing options
    std::map<std::string, QMState> _singlet_levels;
    std::map<std::string, QMState> _triplet_levels;
    
    std::map<std::string, QMState> _hole_levels;
    std::map<std::string, QMState> _electron_levels;

        
};

}}
#endif	
