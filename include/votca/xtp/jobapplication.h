/*
 *            Copyright 2009-2017 The VOTCA Development Team
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


#ifndef VOTCA_XTP_JOBAPPLICATION
#define	VOTCA_XTP_JOBAPPLICATION

#include <votca/xtp/xtpapplication.h>

#include <votca/ctp/progressobserver.h>
#include <votca/ctp/topology.h>

#include "statesaversqlite.h"
#include <votca/ctp/jobcalculator.h>


namespace votca { namespace xtp {



class JobApplication : public XtpApplication
{
public:
    JobApplication();
   ~JobApplication() {
       for (ctp::JobCalculator* calculator : _calculators) {
            delete calculator;
        } 
   };

   void Initialize();
   bool EvaluateOptions();
   void Run(void);

   virtual void BeginEvaluate(int nThreads, ctp::ProgObserver< std::vector<ctp::Job*>, ctp::Job*, ctp::Job::JobResult> *obs);
   virtual bool EvaluateFrame();
   virtual void EndEvaluate();
   void AddCalculator(ctp::JobCalculator *calculator);

protected:
    
    bool _generate_input, _run, _import;
    ctp::Topology           _top;
    std::list< ctp::JobCalculator* >   _calculators;

};

}}









#endif /* _QMApplication_H */














