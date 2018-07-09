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

#ifndef _CALC_COUPLING_EXCL_H
#define	_CALC_COUPLING_EXCL_H

#include <votca/tools/property.h>

#include <votca/xtp/parallelxjobcalc.h>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <votca/xtp/xmapper.h>
#include <votca/xtp/xjob.h>

namespace votca { namespace xtp {
    
/**
* \brief Evaluates Transition Charge distributions classically
*
* Evaluates the electrostatic classical coupling between molecules in 
* their excited states.
* 

* 
* Callname: iexcitoncl
*/

class IEXCITON : public xtp::ParallelXJobCalc< vector<xtp::Job*>, xtp::Job*, xtp::Job::JobResult >
{
public:

    IEXCITON() {};
   ~IEXCITON() {};
   
    void    Initialize(tools::Property *options );
    
    string  Identify() { return "iexcitoncl"; }
    
    xtp::Job::JobResult EvalJob(xtp::Topology *top, xtp::Job *job, xtp::QMThread *Thread);

    void WriteJobFile(xtp::Topology *top);
    void ReadJobFile(xtp::Topology *top);

    

private:

    
    double                              _cutoff;
    double                              _epsilon;
    xtp::XMpsMap                        _mps_mapper;
    bool                           _induce;
    int                           _statenumber;
    string                         _emp_file;
    string                         _xml_file;

        

        
    
    void PreProcess(xtp::Topology *top);
    void CustomizeLogger(xtp::QMThread *thread);
    double EvaluatePair(xtp::Topology *top,xtp::PolarSeg* Seg1,xtp::PolarSeg* Seg2, xtp::Logger* pLog);
 
    
 
        
};

}}
#endif	/* _CALC_INTEGRALS_DFT_H */
