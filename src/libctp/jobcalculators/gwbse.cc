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

// Overload of uBLAS prod function with MKL/GSL implementations
#include <votca/ctp/votca_ctp_config.h>

#include "gwbse.h"
#include <votca/ctp/mbpt.h>


#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>
// #include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

using boost::format;
using namespace boost::filesystem;

namespace votca {
    namespace ctp {
        namespace ub = boost::numeric::ublas;

        // +++++++++++++++++++++++++++++ //
        // GWBSE MEMBER FUNCTIONS         //
        // +++++++++++++++++++++++++++++ //

        void GWBSE::CleanUp() {

        }

        void GWBSE::Initialize(Property *options) {

            _maverick = (_nThreads == 1) ? true : false;

            // setting some defaults
            _mbpt.set_do_qp_diag( false );
            _mbpt.set_do_bse_singlets(false);
            _mbpt.set_do_bse_triplets(false);
            _mbpt.set_ranges("default");
            _mbpt.set_store_qp_pert(true);
            
            // _bse_nmax        = 100;
            
            string key = "options." + Identify() + ".job";
            _jobfile = options->get(key + ".file").as<string>();

            key = "options." + Identify();
            // getting level ranges 
            _mbpt.set_ranges ( options->get(key + ".ranges").as<string> () );
            // now check validity, and get rpa, qp, and bse level ranges accordingly
            if (_mbpt.get_ranges() == "factor") {
                // get factors
                _mbpt.set_rpamaxfactor ( options->get(key + ".rpamax").as<double> () );
                _mbpt.set_qpminfactor( options->get(key + ".qpmin").as<double> () );
                _mbpt.set_qpmaxfactor( options->get(key + ".qpmax").as<double> () );
                _mbpt.set_bseminfactor( options->get(key + ".bsemin").as<double> () );
                _mbpt.set_bsemaxfactor( options->get(key + ".bsemax").as<double> () );
            } else if (_mbpt.get_ranges() == "explicit") {
                //get explicit numbers
                _mbpt.set_rpamax( options->get(key + ".rpamax").as<unsigned int> () );
                _mbpt.set_qpmin( options->get(key + ".qpmin").as<unsigned int> () );
                _mbpt.set_qpmax( options->get(key + ".qpmax").as<unsigned int> () );
                _mbpt.set_bse_vmin( options->get(key + ".bsemin").as<unsigned int> () );
                _mbpt.set_bse_cmax( options->get(key + ".bsemax").as<unsigned int> () );
            } else if  ( _mbpt.get_ranges() == "" ){
                _mbpt.set_ranges("default");
            } else {
                cerr << "\nSpecified range option " << _mbpt.get_ranges() << " invalid. ";
                throw std::runtime_error("\nValid options are: default,factor,explicit");
            }
            
            _mbpt.set_bse_nmax( options->get(key + ".exctotal").as<int> () );
            
            
            _mbpt.set_gwbasis_name( options->get(key + ".gwbasis").as<string> () );
            _mbpt.set_dftbasis_name( options->get(key + ".dftbasis").as<string> () );
            _mbpt.set_shift( options->get(key + ".shift").as<double> () );


            // possible tasks
            // diagQP, singlets, triplets, all
            string _tasks_string = options->get(key+".tasks").as<string> ();
            if (_tasks_string.find("all") != std::string::npos) {
                _mbpt.set_do_qp_diag(true);
                _mbpt.set_do_bse_singlets(true);
                _mbpt.set_do_bse_triplets(true);
            }
            if (_tasks_string.find("qpdiag") != std::string::npos) _mbpt.set_do_qp_diag(true);
            if (_tasks_string.find("singlets") != std::string::npos) _mbpt.set_do_bse_singlets(true);
            if (_tasks_string.find("triplets") != std::string::npos) _mbpt.set_do_bse_triplets(true);
            
            // possible storage 
            // qpPert, qpdiag_energies, qp_diag_coefficients, bse_singlet_energies, bse_triplet_energies, bse_singlet_coefficients, bse_triplet_coefficients

            string _store_string = options->get(key+".store").as<string> ();
            if ((_store_string.find("all") != std::string::npos) ||(_store_string.find("") != std::string::npos))  {
                // store according to tasks choice
                if ( _mbpt.get_do_qp_diag() ) _mbpt.set_store_qp_diag(true);
                if ( _mbpt.get_do_bse_singlets() ) _mbpt.set_store_bse_singlets(true);
                if ( _mbpt.get_do_bse_triplets() ) _mbpt.set_store_bse_triplets(true);
            }
            if (_store_string.find("qpdiag") != std::string::npos) _mbpt.set_store_qp_diag(true);
            if (_store_string.find("singlets") != std::string::npos) _mbpt.set_store_bse_singlets(true);
            if (_store_string.find("triplets") != std::string::npos) _mbpt.set_store_bse_triplets(true);
    

        }

        Job::JobResult GWBSE::EvalJob(Topology *top, Job *job, QMThread *opThread) {

            cout << "Starting GW-BSE";
            Orbitals _orbitals;
            Job::JobResult jres = Job::JobResult();
            Property _job_input = job->getInput();
            list<Property*> lSegments = _job_input.Select("segment");

            vector < Segment* > segments;
            int segId = lSegments.front()->getAttribute<int>("id");
            string segType = lSegments.front()->getAttribute<string>("type");

            Segment *seg = top->getSegment(segId);
            assert(seg->Name() == segType);
            segments.push_back(seg);

            Logger* pLog = opThread->getLogger();
            LOG(logINFO, *pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush;

            // load the DFT data from serialized orbitals object
            string orb_file = (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
            string frame_dir = "frame_" + boost::lexical_cast<string>(top->getDatabaseId());
            string edft_work_dir = "OR_FILES";
            string DIR = edft_work_dir + "/molecules_gwbse/" + frame_dir;
            std::ifstream ifs((DIR + "/" + orb_file).c_str());
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loading DFT data from " << DIR << "/" << orb_file << flush;
            boost::archive::binary_iarchive ia(ifs);
            ia >> _orbitals;
            ifs.close();

            _mbpt.setLogger(pLog);
            bool _evaluate = _mbpt.Evaluate( &_orbitals );
            
            
            LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating site " << seg->getId() << flush; 
 
            Property _job_summary;
            Property *_output_summary = &_job_summary.add("output","");
            Property *_segment_summary = &_output_summary->add("segment","");
            string segName = seg->getName();
            segId = seg->getId();
            _segment_summary->setAttribute("id", segId);
            _segment_summary->setAttribute("type", segName);
            // output of the JOB 
            jres.setOutput( _job_summary );
            jres.setStatus(Job::COMPLETE);

            // dump the LOG
            cout << *pLog;
            return jres;
        }

    }
    
 
};
