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
#include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>

#include <boost/numeric/ublas/symmetric.hpp>
#include <votca/tools/linalg.h>

using boost::format;
using namespace boost::filesystem;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
// +++++++++++++++++++++++++++++ //
// GWBSE MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //
/*
    void GWBSE::FillOverlap(ub::matrix<double>* overlap, BasisSet* bs, vector<Segment* > segments){
        // supposed to fill atomic orbital overlap matrix
        cout << "\nYahoo, I'm filling an atomic orbital overlap matrix..." << endl;
        cout << "... well, maybe at some point :( " << endl;
        
        cout << "First, let's try and determine the size of the AO basis:" << endl;
        
        vector< Atom* > _atoms;
        vector< Atom* > ::iterator ait;
        vector< Segment* >::iterator sit;

        unsigned _AObasis_size = 0;
        
        // loop over segments
        for (sit = segments.begin() ; sit != segments.end(); ++sit) {
        
            _atoms = (*sit)-> Atoms();
            // loop over atoms in segment
            for (ait = _atoms.begin(); ait < _atoms.end(); ++ait) {
                // get coordinates of this atom
                vec     pos = (*ait)->getQMPos();
                // get element type of the atom
                string  name = (*ait)->getElement();
                
                // get the basis set entry for this element
                Element* element = bs->getElement(name);
                // and loop over all shells
                for (Element::ShellIterator its = element->firstShell(); its != element->lastShell(); its++) {
                        
                    Shell* shell = (*its);
                    // we don't like contracted basis sets yet
                    if ( shell->getSize() > 1 ) {
                        cerr << "No contracted basis sets!" << flush;
                    } else {
                        _AObasis_size += NumFuncShell( shell->getType() );
                    }
                        
                    // shell type, number primitives, scale factor
                    //_com_file << shell->getType() << " " << shell->getSize() << " " << shell->getScale() << endl;
                    //shell->getType() 
                    /* for (Shell::GaussianIterator itg = shell->firstGaussian(); itg != shell->lastGaussian(); itg++) {
                        GaussianPrimitive* gaussian = *itg;
                        //_com_file << gaussian->decay << " " << gaussian->contraction << endl;
                        _el_file << FortranFormat( gaussian->decay )<< " " << FortranFormat( gaussian->contraction ) << endl; 
                    }
                }
            }
        }
        cout << "Atomic orbitals basis set size: " << _AObasis_size << endl;
    }
    */
    
    
   
 
    void GWBSE::CleanUp() {

}
    
    
    int GWBSE::NumFuncShell( string shell_type ) {
    int _nbf;
    if ( shell_type == "S" ){
        _nbf = 1;
    } else if ( shell_type == "P" ){
        _nbf = 3;
    } else if ( shell_type == "D" ){
        _nbf = 5;
    } else if ( shell_type == "SP" ) {
        _nbf = 4;
    } else if ( shell_type == "SPD" ) {
        _nbf = 9;
    }
    return _nbf;
}

    int GWBSE::OffsetFuncShell( string shell_type ) {
    int _nbf;
    if ( shell_type == "S" ){
        _nbf = 0;
    } else if ( shell_type == "P" ){
        _nbf = 1;
    } else if ( shell_type == "D" ){
        _nbf = 4;
    } else if ( shell_type == "SP" ) {
        _nbf = 0;
    } else if ( shell_type == "SPD" ) {
        _nbf = 0;
    }
    return _nbf;
}

    
    
void GWBSE::Initialize(Property *options) {

    _maverick = (_nThreads == 1) ? true : false;
    
    /* obsolete string key = "options." + Identify();
    _jobfile = options->get(key + ".job_file").as<string>(); */
    
    string key = "options." + Identify() +".job";
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

    /*string output;
    
    FILE *out;
    Orbitals _orbitals;
    Job::JobResult jres = Job::JobResult();



    // log, com, and orbital files will be stored in ORB_FILES/package_name/frame_x/mol_ID/
    // extracted information will be stored in  ORB_FILES/molecules/frame_x/molecule_ID.orb
    
    string edft_work_dir = "OR_FILES";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    string ID   = boost::lexical_cast<string>( seg->getId() );

   
  */
    
    cout << "Starting GW-BSE"   ;      
    Orbitals _orbitals;
    Job::JobResult jres = Job::JobResult();
    Property _job_input = job->getInput();  
    list<Property*> lSegments = _job_input.Select( "segment" );  
    
    vector < Segment* > segments;    
    int segId = lSegments.front()->getAttribute<int>( "id" );
    string segType = lSegments.front()->getAttribute<string>( "type" );
    
    Segment *seg = top->getSegment( segId );
    assert( seg->Name() == segType );
    segments.push_back( seg );

    Logger* pLog = opThread->getLogger();
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush; 
    
    

    /* obsolete Job::JobResult jres = Job::JobResult();

    vector < Segment* > segments;    
    Segment *seg = top->getSegment( boost::lexical_cast<int>( job->getTag() ));
    segments.push_back( seg );

    Logger* pLog = opThread->getLogger();
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating site " << seg->getId() << flush; 
     */ 
       
    // test setting up ao_overlap_matrix
    list<string> elements;
    BasisSet bs;
    string basis_name("gwdefault");
    
    AOBasis gwbasis;
    bool PPM_symmetric = true; // only PPM supported
    
    
    bs.LoadBasisSet( basis_name );
    LOG(logDEBUG,*pLog) << TimeStamp() << " Loaded Basis Set " << basis_name << flush;
    
    gwbasis.AOBasisFill( &bs, segments );
    LOG(logDEBUG,*pLog) << TimeStamp() << " Filled GW Basis of size " << gwbasis._AOBasisSize << flush;
    
    
    
    // test output
    /*for (vector< AOShell* >::iterator it = aobasis.firstShell(); it != aobasis.lastShell() ; it++ ) {
        AOShell* _shell = aobasis.getShell( it );
        cout << _shell->getType() << ":" << _shell->getScale() << ":" << _shell->getNumFunc() << ":"  << _shell->getStartIndex() << ":" << _shell->getPos() << endl;
        for (AOShell::GaussianIterator itg = _shell->firstGaussian(); itg != _shell->lastGaussian(); itg++) {
            cout << "\t" << (*itg)->decay << ":" << (*itg)->contraction << endl;
        }
    }*/

    // get overlap matrix as AOOverlap
    AOOverlap _gwoverlap;
    // initialize overlap matrix
    _gwoverlap.Initialize(gwbasis._AOBasisSize);
    // Fill overlap
    _gwoverlap.Fill( &gwbasis );
    LOG(logDEBUG,*pLog) << TimeStamp() << " Filled  GW Overlap matrix of dimension: " << _gwoverlap._aomatrix.size1() <<  flush;
    // _aooverlap.Print( "S" );
 
    // get Coulomb matrix as AOCoulomb
    AOCoulomb _gwcoulomb;
    // initialize Coulomb matrix
    _gwcoulomb.Initialize(gwbasis._AOBasisSize);
    // Fill Coulomb matrix
    _gwcoulomb.Fill( &gwbasis );
    LOG(logDEBUG,*pLog) << TimeStamp() << " Filled GW Coulomb matrix of dimension: " << _gwcoulomb._aomatrix.size1() <<  flush;
    // _gwcoulomb.Print( "COU" );
    
    // printing some debug info
    // _gwcoulomb.PrintIndexToFunction( &aobasis );
    
    // check eigenvalues of overlap matrix
    ub::vector<double>                  _eigenvalues;
    ub::matrix<double>                  _eigenvectors;
    _eigenvalues.resize( _gwoverlap._aomatrix.size1() );
    _eigenvectors.resize( _gwoverlap._aomatrix.size1(), _gwoverlap._aomatrix.size1() ); 
    linalg_eigenvalues(_gwoverlap._aomatrix , _eigenvalues, _eigenvectors);
    // cout << _eigenvalues << endl;
    sort (_eigenvalues.begin(), _eigenvalues.end() );   
    LOG(logDEBUG,*pLog) << TimeStamp() << " Smallest eigenvalue of GW Overlap matrix : " << _eigenvalues[0] <<  flush;
 
    
    // PPM is symmetric, so we need to get the sqrt of the Coulomb matrix
    if ( PPM_symmetric ){
        
        _gwcoulomb.Symmetrize(  _gwoverlap , gwbasis  );
        
        
 
    }
    return jres;
}

    
    
    
    
}};
  