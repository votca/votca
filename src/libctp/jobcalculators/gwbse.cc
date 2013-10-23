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


#include "gwbse.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <votca/ctp/aomatrix.h>
#include <votca/ctp/threecenters.h>
// #include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>

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

        int GWBSE::NumFuncShell(string shell_type) {
            int _nbf;
            if (shell_type == "S") {
                _nbf = 1;
            } else if (shell_type == "P") {
                _nbf = 3;
            } else if (shell_type == "D") {
                _nbf = 5;
            } else if (shell_type == "SP") {
                _nbf = 4;
            } else if (shell_type == "SPD") {
                _nbf = 9;
            }
            return _nbf;
        }

        int GWBSE::OffsetFuncShell(string shell_type) {
            int _nbf;
            if (shell_type == "S") {
                _nbf = 0;
            } else if (shell_type == "P") {
                _nbf = 1;
            } else if (shell_type == "D") {
                _nbf = 4;
            } else if (shell_type == "SP") {
                _nbf = 0;
            } else if (shell_type == "SPD") {
                _nbf = 0;
            }
            return _nbf;
        }

        void GWBSE::Initialize(Property *options) {

            _maverick = (_nThreads == 1) ? true : false;

            /* obsolete string key = "options." + Identify();
            _jobfile = options->get(key + ".job_file").as<string>(); */

            string key = "options." + Identify() + ".job";
            _jobfile = options->get(key + ".file").as<string>();

            /*_do_input = false;
            _do_run = false;
            _do_parse = false;
            _do_trim = false;
    
            // conversion to GW 
            _do_convert = false;
            _do_gwbse_input = false;
            _do_gwbse_run = false;
            _do_gwbse_parse = false;
        

            string _package_xml = options->get(key+".package").as<string> ();

    
            string _tasks_string = options->get(key+".tasks").as<string> ();
            if (_tasks_string.find("input") != std::string::npos) _do_input = true;
            if (_tasks_string.find("run") != std::string::npos) _do_run = true;
            if (_tasks_string.find("trim") != std::string::npos) _do_trim = true;
            if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;    
            // GW-BSE tasks
            if (_tasks_string.find("convert") != std::string::npos) _do_convert = true;   
            if (_tasks_string.find("gwbse_setup") != std::string::npos) _do_gwbse_input = true;
            if (_tasks_string.find("gwbse_exec") != std::string::npos) _do_gwbse_run = true;    
            if (_tasks_string.find("gwbse_read") != std::string::npos) _do_gwbse_parse = true;
    
            string _store_string = options->get(key+".store").as<string> ();
            if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
            if (_store_string.find("qppert") != std::string::npos) _store_qppert = true;
            if (_store_string.find("qpdiag") != std::string::npos) _store_qpdiag = true;
            if (_store_string.find("singlets") != std::string::npos) _store_singlets = true;
            if (_store_string.find("triplets") != std::string::npos) _store_triplets = true;
    
            load_property_from_xml( _package_options, _package_xml.c_str() );    
            key = "package";
            _package = _package_options.get(key+".name").as<string> ();


   
            // only required, if GWBSE is to be run
            if ( _do_gwbse_input || _do_gwbse_run || _do_gwbse_parse ){
                key = "options." + Identify();
                string _gwpackage_xml = options->get(key+".gwpackage").as<string> ();
                load_property_from_xml( _gwpackage_options, _gwpackage_xml.c_str() );  
                key = "package";
                _gwpackage = _gwpackage_options.get(key+".name").as<string> ();
            }
    
    
            // register all QM packages (Gaussian, turbomole, nwchem))
            QMPackageFactory::RegisterAll(); */
            cout << "I'm supposed to initialize GWBSE";

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

            // load the DFT data 
            string orb_file = (format("%1%_%2%%3%") % "molecule" % segId % ".orb").str();
            string frame_dir = "frame_" + boost::lexical_cast<string>(top->getDatabaseId());
            string edft_work_dir = "OR_FILES";
            string DIR = edft_work_dir + "/molecules_gwbse/" + frame_dir;
            std::ifstream ifs((DIR + "/" + orb_file).c_str());
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loading DFT data from " << DIR << "/" << orb_file << flush;
            boost::archive::binary_iarchive ia(ifs);
            ia >> _orbitals;
            ifs.close();
            string _dft_package = _orbitals.getQMpackage();
            LOG(logDEBUG, *pLog) << TimeStamp() << " DFT data was created by " << _dft_package << flush;

            // reorder DFT data, load DFT basis set
            BasisSet dftbs;
            string dftbasis_name("ubecppol");

            AOBasis dftbasis;

            dftbs.LoadBasisSet(dftbasis_name);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loaded DFT Basis Set " << dftbasis_name << flush;

            dftbasis.AOBasisFill(&dftbs, segments);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled DFT Basis of size " << dftbasis._AOBasisSize << flush;

            // do the reordering depending on the QM package used to obtain the DFT data
            ub::matrix<double> _dft_orbitals = *_orbitals.getOrbitals();
            if (_dft_package != "votca") {
                // get reordering vector _dft_package -> Votca 
                vector<int> neworder;
                dftbasis.getReorderVector(_dft_package, neworder);
                // and reorder rows of _orbitals->_mo_coefficients() accordingly
                AOBasis::ReorderMOs(_dft_orbitals, neworder);
                // NWChem inverted sign for xz d-orbital
                if (_dft_package == "nwchem") {
                    // get vector with multipliers, e.g. NWChem -> Votca (bloody sign for d_xz)
                    vector<int> multiplier;
                    dftbasis.getMultiplierVector(_dft_package, multiplier);
                    // and reorder rows of _orbitals->_mo_coefficients() accordingly
                    AOBasis::MultiplyMOs(_dft_orbitals, multiplier);
                }
            }
            LOG(logDEBUG, *pLog) << TimeStamp() << " Converted DFT orbital coefficient order " << flush;

            // setting up ao_overlap_matrix
            list<string> elements;
            BasisSet gwbs;
            string gwbasis_name("gwdefault");

            AOBasis gwbasis;
            bool PPM_symmetric = true; // only PPM supported


            gwbs.LoadBasisSet(gwbasis_name);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Loaded GW Basis Set " << gwbasis_name << flush;

            gwbasis.AOBasisFill(&gwbs, segments);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled GW Basis of size " << gwbasis._AOBasisSize << flush;

            // get overlap matrix as AOOverlap
            AOOverlap _gwoverlap;
            // initialize overlap matrix
            _gwoverlap.Initialize(gwbasis._AOBasisSize);
            // Fill overlap
            _gwoverlap.Fill(&gwbasis);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled GW Overlap matrix of dimension: " << _gwoverlap._aomatrix.size1() << flush;
            // _aooverlap.Print( "S" );



            // printing some debug info
            // _gwcoulomb.PrintIndexToFunction( &aobasis );

            // check eigenvalues of overlap matrix
            ub::vector<double> _eigenvalues;
            ub::matrix<double> _eigenvectors;
            _eigenvalues.resize(_gwoverlap._aomatrix.size1());
            _eigenvectors.resize(_gwoverlap._aomatrix.size1(), _gwoverlap._aomatrix.size1());
            linalg_eigenvalues(_gwoverlap._aomatrix, _eigenvalues, _eigenvectors);
            // cout << _eigenvalues << endl;
            sort(_eigenvalues.begin(), _eigenvalues.end());
            LOG(logDEBUG, *pLog) << TimeStamp() << " Smallest eigenvalue of GW Overlap matrix : " << _eigenvalues[0] << flush;



            // get Coulomb matrix as AOCoulomb
            AOCoulomb _gwcoulomb;
            // initialize Coulomb matrix
            _gwcoulomb.Initialize(gwbasis._AOBasisSize);
            // Fill Coulomb matrix
            _gwcoulomb.Fill(&gwbasis);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Filled GW Coulomb matrix of dimension: " << _gwcoulomb._aomatrix.size1() << flush;
            // _gwcoulomb.Print( "COU" );


            // PPM is symmetric, so we need to get the sqrt of the Coulomb matrix
            if (PPM_symmetric) {
                // _gwcoulomb.Symmetrize(  _gwoverlap , gwbasis  );
                ;
            }
            LOG(logDEBUG, *pLog) << TimeStamp() << " Prepared GW Coulomb matrix for symmetric PPM " << endl;
            // scissors-shift DFT energies

            // calculate 3-center integrals,  convoluted with DFT eigenvectors

            // --- prepare a vector (gwdacay) of matrices (orbitals, orbitals) as container => M_mn
            mmin = 1; // lowest index occ 
            mmax = 2 * _orbitals.getNumberOfElectrons();
            nmin = 1;
            nmax = _orbitals.getNumberOfLevels();
            maxf = gwbasis.getMaxFunctions(); // maximum number of functions per shell in basis set
            mtotal = mmax - mmin + 1;
            ntotal = nmax - nmin + 1;


            // prepare 3-center integral object
            TCMatrix _Mmn;
            _Mmn.Initialize(gwbasis._AOBasisSize, mmin, mmax, nmin, nmax);
            _Mmn.Fill(gwbasis, dftbasis, _dft_orbitals);
            LOG(logDEBUG, *pLog) << TimeStamp() << " Calculated Mmn_beta (3-center-overlap x orbitals)  " << endl;





            return jres;
        }





    }
};
