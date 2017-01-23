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

#ifndef __PARALLELXJOBCALC__H
#define __PARALLELXJOBCALC__H


#include <votca/ctp/jobcalculator.h>
#include <votca/ctp/qmthread.h>
#include <votca/tools/mutex.h>
#include <votca/ctp/job.h>
#include <votca/ctp/progressobserver.h>


// PATHWAYS TO A NEW THREADED CALCULATOR
// ... 1 Define 'JobContainer' (needs to define iterator), 'pJob' ( = *iterator)
// ... 2 Derive new calculator as ': public ParallelXJobCalc<JobContainer,pJob>'
// ... 3 Specialize XJOBS_FROM_TABLE< JobContainer, pJob> in xjob.cc
// ... 4 Register new calculator (see end of parallelxjobcalc.cc)

// REQUIRED MEMBERS FOR pJob
// pJob::JobResult (struct)


namespace votca { namespace xtp {

    namespace CTP = votca::ctp;
template<typename JobContainer, typename pJob, typename rJob> 
class ParallelXJobCalc : public CTP::JobCalculator
{

public:

    class JobOperator;
    
    ParallelXJobCalc() : _jobfile("__NOFILE__") {};
   ~ParallelXJobCalc() { ; };

    std::string       Identify() { return "Parallel XJob Calculator"; }

    bool         EvaluateFrame(CTP::Topology *top);
    virtual void LoadJobs() { ; }
    virtual void CustomizeLogger(CTP::QMThread* thread);
    virtual void PreProcess(CTP::Topology *top) { ; } 
    virtual rJob EvalJob(CTP::Topology *top, const pJob job, CTP::QMThread *thread) = 0;
    virtual void PostProcess(CTP::Topology *top) { ; }
    
    void         LockCout() { _coutMutex.Lock(); }
    void         UnlockCout() { _coutMutex.Unlock(); }
    void         LockLog() { _logMutex.Lock(); }
    void         UnlockLog() { _logMutex.Unlock(); }

    
    // ======================================== //
    // XJOB OPERATOR (THREAD)                   //
    // ======================================== //
    

    class JobOperator : public CTP::QMThread
    {
    public:

        JobOperator(int id,   CTP::Topology *top, ParallelXJobCalc<JobContainer,pJob,rJob> *master)
                      : _top(top),          _master(master) { _id = id; };
       ~JobOperator() {};

        void        InitData(CTP::Topology *top) { ; }
        void        Run(void);
        

    public:

        CTP::Topology         *_top;
        ParallelXJobCalc<JobContainer,pJob,rJob> *_master;
        pJob              _job;

    };

    
    
    


protected:

    JobContainer             _XJobs;
    Mutex                    _coutMutex;
    Mutex                    _logMutex;
    std::string                   _jobfile;
    int                      _subthreads;
    
    // ProgObserver< JobContainer, pJob > *_progObs;


};

}}

#endif
