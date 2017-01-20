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

#ifndef _CALC_COUPLING_EXCL_H
#define	_CALC_COUPLING_EXCL_H

#include <votca/tools/property.h>

#include <votca/xtp/parallelxjobcalc.h>
#include <boost/numeric/ublas/io.hpp>
#include <sys/stat.h>
#include <boost/filesystem.hpp>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xjob.h>

namespace votca { namespace xtp {
    
    namespace CTP = votca::ctp;
/**
* \brief Evaluates Transition Charge distributions classically
*
* Evaluates the electrostatic classical coupling between molecules in 
* their excited states.
* 

* 
* Callname: iexcitoncl
*/

class IEXCITON : public ParallelXJobCalc< vector<CTP::Job*>, CTP::Job*, CTP::Job::JobResult >
{
public:

    IEXCITON() {};
   ~IEXCITON() {};
   
    void    Initialize(votca::tools::Property *options );
    
    string  Identify() { return "iexcitoncl"; }
    
    CTP::Job::JobResult EvalJob(CTP::Topology *top, CTP::Job *job, CTP::QMThread *Thread);

    void WriteJobFile(CTP::Topology *top);
    void ReadJobFile(CTP::Topology *top);

    

private:

    
    
   
    CTP::XMpsMap                        _mps_mapper;
    bool                           _induce;
    int                           _statenumber;
    string                         _emp_file;
    string                         _xml_file;

        

        
    
    void PreProcess(CTP::Topology *top);
    void CustomizeLogger(CTP::QMThread *thread);
    double EvaluatePair(CTP::Topology *top, CTP::PolarSeg* Seg1, CTP::PolarSeg* Seg2, CTP::Logger* pLog);
 
    
 
        
};

}}
#endif	/* _CALC_INTEGRALS_DFT_H */
