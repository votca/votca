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


#ifndef _VOTCA_XTP_DMA_H
#define	_VOTCA_XTP_DMA_H



#include <votca/xtp/qmpackagefactory.h>

#include <votca/ctp/parallelxjobcalc.h>
#include <votca/ctp/segment.h>

#include <fstream>

#include <boost/filesystem.hpp>
#include <boost/format.hpp>

using boost::format;

namespace votca { namespace xtp {

/**
* \brief Distributed multipole analysis using Gaussian input
*
* Evaluates distributed multipoles
* Requires GAUSSIAN and GDMA
*
* Callname: dma
*/

class DMA : public ctp::ParallelXJobCalc< vector<ctp::Job*>, ctp::Job*, ctp::Job::JobResult >
{
public:

    DMA() {};
   ~DMA() {};

    string   Identify() { return "dma"; }
    void     Initialize(Property *options);
    void     WriteJobFile(ctp::Topology *top);
    
    ctp::Job::JobResult EvalJob(ctp::Topology *top, ctp::Job *job, ctp::QMThread *thread);

private:

    // what to do
    bool                _do_input;
    bool                _do_orbitals;
    bool                _do_dma;

    string _executable;
    string _jobFile;
    string _chkFile;
            
    string _package;
    Property _package_options;   
    
    // dma input file options
    string _density;
    int _limit;
    double _radius;
    int _switch;
    string _outFile;
};

}}

#endif	/* _VOTCA_XTP_DMA_H */
