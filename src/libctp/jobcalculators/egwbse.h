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

// UBLAS stops checking types and array bounds if this flag is defined
#define NDEBUG
#define BOOST_UBLAS_NDEBUG

#ifndef _CALC_GWBSE_TOOL_H
#define	_CALC_GWBSE_TOOL_H


#include <votca/ctp/gwbse.h> // including GWBSE functionality
#include <votca/ctp/segment.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/aobasis.h>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>

#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>

#include <fstream>
#include <sys/stat.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>
#include <votca/tools/linalg.h>

#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>


namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
/**
* \brief GWBSE implementation
*
* Requires a first-principles package, i.e. GAUSSIAN installation
*
* Callname: egwbse
*/

class EGWBSE : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    EGWBSE() {};
   ~EGWBSE() {};

    string  Identify() { return "egwbse"; }
    void    Initialize( Property *options);
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);
    
    string              _package;
    Property            _package_options; 
    Property            _gwbse_options; 
    
    void    CleanUp();

    
    // what to do
    bool                _do_dft_input;
    bool                _do_dft_run;
    bool                _do_dft_parse;
    bool                _do_gwbse;

    
    // all GWBSE functionality is in GWBSE object
    GWBSE _gwbse; 
    
    
    void ParseOptionsXML( Property *options);    

};


}}

#endif	/* _CALC_GWBSE_TOOL_H */
