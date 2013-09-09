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


#include "idft.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <votca/ctp/eigenvalues.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>

using boost::format;
using namespace boost::filesystem;

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize(ctp::Topology *top, tools::Property* options ) {
    
    _energy_difference = 0.0;
    
    _do_input = false;
    _do_run = false;
    _do_parse = false;
    _do_project = false;
    _do_trim = false;
    
    _store_orbitals = false;
    _store_overlap = false;
    _store_integrals = false;
    
    ParseOptionsXML( options  );
    
    // register all QM packages (Gaussian, turbomole, etc))
    QMPackageFactory::RegisterAll();

}

    
void IDFT::ParseOptionsXML( tools::Property *opt ) {
   
    // Orbitals are in fort.7 file; number of electrons in .log file
    
    string key = "options." + Identify();
    _energy_difference = opt->get( key + ".degeneracy" ).as< double > ();
    
    _jobfile = opt->get(key + ".job_file").as<string>();

    string _tasks_string = opt->get(key+".tasks").as<string> ();
    if (_tasks_string.find("input") != std::string::npos) _do_input = true;
    if (_tasks_string.find("run") != std::string::npos) _do_run = true;
    if (_tasks_string.find("parse") != std::string::npos) _do_parse = true;
    if (_tasks_string.find("project") != std::string::npos) _do_project = true;
    if (_tasks_string.find("trim") != std::string::npos) _do_trim = true;

    string _store_string = opt->get(key+".store").as<string> ();
    if (_store_string.find("orbitals") != std::string::npos) _store_orbitals = true;
    if (_store_string.find("overlap") != std::string::npos) _store_overlap = true;
    if (_store_string.find("integrals") != std::string::npos) _store_integrals = true;
    
    _max_occupied_levels = opt->get(key+".levels").as<int> ();
    _max_unoccupied_levels = _max_occupied_levels;

    _trim_factor = opt->get(key+".trim").as<int> ();
    
    string _package_xml = opt->get(key+".package").as<string> ();
    //cout << endl << "... ... Parsing " << _package_xml << endl ;

    load_property_from_xml( _package_options, _package_xml.c_str() );    
    
     key = "package";
    _package = _package_options.get(key+".name").as<string> ();

}


double inv_sqrt(double x) { return 1./sqrt(x); }

/*
 * Calculates S^{-1/2}
 */
void IDFT::SQRTOverlap(ub::symmetric_matrix<double> &S, ub::matrix<double> &S2 ) {
       
    double (*_inv_sqrt)(double);
    _inv_sqrt = &inv_sqrt;

    ub::vector<double>                  _eigenvalues;
    ub::matrix<double>                  _eigenvectors;

    int _size = S.size1(); 

    //cout << "... ... Calculating SQRT of the " << _size << "x" << _size  << " overlap matrix" << endl;

    _eigenvalues.resize( _size );
    _eigenvectors.resize( _size, _size ); 
    
    
//  test case  

/*  test of the eigendecomposition code
    int _basis_size = 3;
    ub::symmetric_matrix<double> _overlap(_basis_size);
    _overlap.resize( _basis_size ); 
    _eigenvalues.resize( _basis_size );
    _eigenvectors.resize( _basis_size, _basis_size ); 
    
    //eigenvalues 3, 6, 9
    //eigenvectors (1,2,2), (-2,-1,2), (2,-2,1)
   
    _overlap(0,0) = 7;   
    _overlap(1,0) =-2;  _overlap(1,1) = 6;  
    _overlap(2,0) = 0;  _overlap(2,1) =-2; _overlap(2,2) = 5;

    EigenvaluesSymmetric(_overlap, _eigenvalues, _eigenvectors);
    cout << "....eigenvalue problem solved " << endl;
    
    cout << "eigenvalues" << _eigenvalues << endl;
    cout << "eigenvectors" << _eigenvectors << endl;
    
    ub::diagonal_matrix<double> _diag( _eigenvalues.size(), _eigenvalues.data() );
    ub::matrix<double> _left = ub::prod( _eigenvectors, _diag );
    cout <<  ub::prod( _left, ub::trans( _eigenvectors ) );
    
    exit(0);
*/    
    /* for the test case above S2 has the following form 
    * [[0.3937418627,0.07087375404,0.0209304492],
    *  [0.07087375404,0.4501091889,0.0918042032],
    *  [0.0209304492,0.0918042032,0.4750808413]]
    */
    
    
    EigenvaluesSymmetric(S, _eigenvalues, _eigenvectors);
    //cout << "... ... Eigenvalue problem solved " << endl;
    
    //cout << "eigenvalues" << _eigenvalues << endl;
    //cout << _eigenvectors << endl;     
    
    // compute inverse sqrt of all eigenvalues
    std::transform(_eigenvalues.begin(), _eigenvalues.end(), _eigenvalues.begin(),  _inv_sqrt );

    // form a diagonal matrix S^{-1/2}
    ub::diagonal_matrix<double> _diagS2( _eigenvalues.size(), _eigenvalues.data() ); 

    // multiply from the left on the U
    ub::matrix<double> _temp = ub::prod( _eigenvectors, _diagS2 );
    
    // multiply from the right on the transpose U
    S2 = ub::prod( _temp, ub::trans( _eigenvectors ) );
    //cout << "... ... Projection matrix constructed  " << endl;
       


    // cleanup
    _diagS2.clear();
    _temp.clear();

    //cout << "S2: " << S2 << endl;
    //cout << "Overlap: " << _overlap << endl;
    
    //cout << "... ... Done with the sqrt of the overlap matrix" << endl;
    
    
 }

bool IDFT::CalculateIntegrals(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
    Orbitals* _orbitalsAB, ub::matrix<double>* _JAB, QMThread *opThread) {
          
    /* test case
    ub::matrix<double> _monomersAB (4, 5);
    ub::zero_matrix<double> _AB (4, 5);

    _monomersAB = _AB;
    
    std::cout << _monomersAB << std::endl;
    
    ub::matrix<double> C(2, 2);
    C(0,0) = 3; C(0,1) = 3;
    C(1,0) = 3; C(1,1) = 3;
    
    ub::matrix<double> B(2, 2);
    B(0,0) = 5; B(0,1) = 5;
    B(1,0) = 5; B(1,1) = 5;
    
    ub::project(_monomersAB, ub::range (2, 4), ub::range (3, 5)) = C;
    ub::project(_monomersAB, ub::range (0, 2), ub::range (0, 2)) = B;

    std::cout << _monomersAB << std::endl;
    */
    
    Logger* _pLog = opThread->getLogger();
    LOG(logDEBUG,*_pLog) << "Calculating electronic couplings" << flush;
        
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    if ( ( _basisA == 0 ) || ( _basisB == 0 ) ) {
        LOG(logERROR,*_pLog) << "No basis set size information is stored in monomers" << flush;
        return false;
    }
    
    LOG(logDEBUG,*_pLog) << "Basis [molA:molB] " << _basisA << ":" << _basisB << flush;
    
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    if ( ( _levelsA == 0 ) || (_levelsB == 0) ) {
        LOG(logERROR,*_pLog) << "No information about number of occupied/unoccupied levels is stored" << flush;
        return false;
    } 
     
    
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
        
    //cout << zeroB << endl;
    //cout << zeroA << endl;
    
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    
    // AxB = | A 0 |  //
    //       | 0 B |  //  
    LOG(logDEBUG,*_pLog) << "Constructing direct product AxB" << flush;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals();    
    //cout << "_psi_AxB: " << _psi_AxB << endl;
    
    // Fock matrix of a dimer   
    LOG(logDEBUG,*_pLog) << "Constructing the dimer Fock matrix" << flush;    
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB->getNumberOfLevels(), (*_orbitalsAB->getEnergies()).data() ); 

    // psi_AxB * S_AB * psi_AB
    LOG(logDEBUG,*_pLog) << "Projecting the dimer onto monomer orbitals" << flush;    
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB->getOverlap(), ub::trans( *_orbitalsAB->getOrbitals() ) );          
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );
    _psi_AB.clear();
   
    /*
    for (int i = 0; i < _psi_AxB_dimer_basis.size1(); i++ ) {
        for (int j = 0; j < _psi_AxB_dimer_basis.size2(); j++ ) {
            cout << i << " " << j << " " << _psi_AxB_dimer_basis.at_element(i, j) << endl;
            
        }
    }
    exit(0);
     */
    
    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    LOG(logDEBUG,*_pLog) << "Projecting the Fock matrix onto the dimer basis" << flush;    
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ;
    ub::matrix<double> JAB_dimer = ub::prod( _psi_AxB_dimer_basis, _temp);
    
    /* DEBUG 
    int levelA = _orbitalsA->getNumberOfElectrons() ;
    int levelB = _orbitalsB->getNumberOfElectrons() ;

    cout << "... ... Coupling before processing " 
            << JAB_dimer.at_element( levelA - 1  , levelB -1 + _levelsA ) * _conv_Hrt_eV << " "
            << JAB_dimer.at_element( levelA , levelB + _levelsA ) * _conv_Hrt_eV << "\n";
    */
    
    _temp.clear(); _fock_AB.clear();
    
    // S = psi_AxB_dimer_basis * psi_AxB_dimer_basis^T
    ub::symmetric_matrix<double> _S_AxB = ub::prod( _psi_AxB_dimer_basis, ub::trans( _psi_AxB_dimer_basis ));
    //cout << "SAxB: " << _S_AxB << endl;

    /* test of an assignment 
    ub::matrix<double> C(2,2);
    C(0,0) = 1; C(0,1) = 2;
    C(1,0) = 2; C(1,1) = 3;
    
    ub::symmetric_matrix<double> B = C;
    cout << C << endl; 
    cout << B << endl; 
    */
       
    ub::matrix<double> _S_AxB_2(_S_AxB.size1(), _S_AxB.size1() );
    
    /* test of the SQRT routine
    ub::symmetric_matrix<double> _test(3,3);
    ub::matrix<double> _test2(3,3); 
    _test(0,0) = 7;   
    _test(1,0) =-2;  _test(1,1) = 6;  
    _test(2,0) = 0;  _test(2,1) =-2; _test(2,2) = 5; 
    SQRTOverlap(_test, _test2 );
    cout << _test2;
    exit(0);
    */
     ub::trans( _S_AxB );
     LOG(logDEBUG,*_pLog) << "Calculating square root of the overlap matrix" << flush;    
     SQRTOverlap( _S_AxB , _S_AxB_2 );        
     _S_AxB.clear(); 
     
    ///if ( tools::globals::verbose ) *opThread << "... ... Calculating the effective overlap\n" ;
    //stringstream test ;
    //test << "BLA" << "BA";
    
    ub::matrix<double> JAB_temp = prod( JAB_dimer, _S_AxB_2 );
        
    (*_JAB) = ub::prod( _S_AxB_2, JAB_temp );
    
    // cleanup
    JAB_dimer.clear(); JAB_temp.clear(); _S_AxB_2.clear();
    
    //cout << JAB << endl;
    
    //cout << _S_AxB << endl;
    //_has_integrals = true;
    //cout << JAB_dimer.at_element( HOMO_A , HOMO_B + _levelsA ) * conv_Hrt_eV << endl; 
    //cout << JAB_dimer.at_element(_levelsA + HOMO_B, HOMO_A ) * conv_Hrt_eV << endl;

    LOG(logDEBUG,*_pLog) << "Done with electronic couplings" << flush;
    return true;   

}

double IDFT::getCouplingElement( int levelA, int levelB,  Orbitals* _orbitalsA,
    Orbitals* _orbitalsB, ub::matrix<double>* _JAB  ) {
       
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();    
    
    if ( _energy_difference != 0 ) {
        std::vector<int> list_levelsA = *_orbitalsA->getDegeneracy( levelA, _energy_difference );
        std::vector<int> list_levelsB = *_orbitalsA->getDegeneracy( levelB, _energy_difference );
        
        double _JAB_sq = 0; double _JAB_one_level;
        
        for (std::vector<int>::iterator iA = list_levelsA.begin()++; iA != list_levelsA.end(); iA++) {
                for (std::vector<int>::iterator iB = list_levelsB.begin()++; iB != list_levelsB.end(); iB++) { 
                    //cout << *iA << ':' << *iB << endl;
                    _JAB_one_level = _JAB->at_element( *iA - 1  , *iB -1 + _levelsA );
                    _JAB_sq +=  _JAB_one_level*_JAB_one_level ;
                }
        }
        
        return sqrt(_JAB_sq / ( list_levelsA.size() * list_levelsB.size() ) ) * _conv_Hrt_eV ;
        
    } else {
        return _JAB->at_element( levelA - 1  , levelB -1 + _levelsA ) * _conv_Hrt_eV;
    }
    // the  matrix should be symmetric, could also return this element
    // _JAB.at_element( _levelsA + levelB - 1  , levelA - 1 );
}
 
Job::JobResult IDFT::EvalJob(Topology *top, Job *job, QMThread *opThread) {

    bool _run_status = false;
    bool _parse_log_status = false ;
    bool _parse_orbitals_status = false;
    bool _calculate_integrals = false;
    stringstream sout;
    string output;
    
     // report back to the progress observer
    Job::JobResult jres = Job::JobResult();
    
    // get the logger from the thread
    Logger* pLog = opThread->getLogger();   
    
    // get the information about the job executed by the thread
    int _job_ID = job->getId();
    string _job_tag = job->getTag();
    
    // job tag is ID_A:ID_B
    Tokenizer _tok ( _job_tag, ":" ); 
    vector<string> _mol_ids;
    _tok.ToVector( _mol_ids );
    
    int ID_A   = boost::lexical_cast<int>( _mol_ids.front() );
    int ID_B   = boost::lexical_cast<int>( _mol_ids.back() );
    
    Segment *seg_A = top->getSegment( ID_A );
    Segment *seg_B = top->getSegment( ID_B );

    vector < Segment* > segments;
    segments.push_back( seg_A );
    segments.push_back( seg_B );
    
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());     

    /* if we really want that this pair exists in the neighbor list
    QMNBList &nblist = top->NBList(); 
    QMPair *qmpair = nblist.FindPair( seg_A, seg_B );
    if ( qmpair == NULL ) {
        output += (format("Pair %1%:%2% does not exist") % sID_A % sID_B ).str() ;
        LOG(logERROR,*pLog) << "Non-existing pair " << sID_A << ":" << sID_B << flush;
        cout << *pLog;
        jres.setOutput( output ); 
        jres.setStatus(Job::FAILED);
        return jres;
    } else {
        FILE *out;
        out = fopen((GAUSS_DIR + "/" + XYZ_FILE).c_str(),"w");
        qmpair->WriteXYZ(out);
        fclose(out);
    }
    */
    
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating pair "  
            << _job_ID << " ["  << ID_A << ":" << ID_B << "] out of " << 
           (top->NBList()).size()  << flush; 

    // set the folders 
    string idft_work_dir = "OR_FILES";
    string edft_work_dir = "OR_FILES";
    string _pair_dir = ( format("%1%%2%%3%%4%%5%") % "pair" % "_" % ID_A % "_" % ID_B ).str();
 
    string _qmpackage_work_dir = idft_work_dir + "/" + _package + "/" + frame_dir + "/" + _pair_dir ;
    //cout << endl << _qmpackage_work_dir << endl;
    
    // get the corresponding object from the QMPackageFactory
    QMPackage *_qmpackage =  QMPackages().Create( _package );

    // set a log file for the package
    _qmpackage->setLog( pLog );       

    // set the run dir 
    _qmpackage->setRunDir( _qmpackage_work_dir );
        
    // get the package options
    _qmpackage->Initialize( &_package_options );
    
    // if asked, prepare the input files
    if (_do_input) {
        

        boost::filesystem::create_directories( _qmpackage_work_dir );           
        // in case we do not want to do an SCF loop for a dimer
        if ( _qmpackage->GuessRequested() ) {
            LOG(logDEBUG, *pLog) << "Preparing the guess" << flush;
            
            Orbitals _orbitalsA;
            Orbitals _orbitalsB;
            Orbitals _orbitalsAB;
            
            // load the corresponding monomer orbitals and prepare the dimer guess 
            string orb_file_A = (format("%1%_%2%%3%") % "molecule" % ID_A % ".orb").str();
            string DIR_A  = edft_work_dir + "/" + "molecules/" + frame_dir;
            std::ifstream ifs_A( (DIR_A +"/" + orb_file_A).c_str() );
            LOG(logDEBUG,*pLog) << "Loading orbitals from " << orb_file_A << flush;   
            boost::archive::binary_iarchive ia_A( ifs_A );
            ia_A >> _orbitalsA;
            ifs_A.close();

            string orb_file_B = (format("%1%_%2%%3%") % "molecule" % ID_B % ".orb").str();
            LOG(logDEBUG,*pLog) << "Loading orbitals from " << orb_file_B << flush;    
            string DIR_B  = edft_work_dir + "/" + "molecules/" + frame_dir;
            std::ifstream ifs_B( (DIR_B +"/" + orb_file_B).c_str() );
            boost::archive::binary_iarchive ia_B( ifs_B );
            ia_B >> _orbitalsB;
            ifs_B.close();
                            
            PrepareGuess(&_orbitalsA, &_orbitalsB, &_orbitalsAB, opThread);

            boost::filesystem::create_directories( _qmpackage_work_dir );
            _qmpackage->WriteInputFile(segments, &_orbitalsAB);

        } else {
            _qmpackage->WriteInputFile(segments);
        }

    } // end of the input
    
    // run the executable
    if ( _do_run ) {
            _run_status = _qmpackage->Run( );
            if ( !_run_status ) {
                    output += "run failed; " ;
                    LOG(logERROR,*pLog) << _qmpackage->getPackageName() << " run failed" << flush;
                    cout << *pLog;
                    jres.setOutput( output ); 
                    jres.setStatus(Job::FAILED);
                    delete _qmpackage;
                    return jres;
            } 
    } // end of the run
 
    Orbitals _orbitalsAB; // This will be later used to write orbitals of the dimer to a file
   // parse the log/orbitals files
    if ( _do_parse ) {
             _parse_log_status = _qmpackage->ParseLogFile( &_orbitalsAB );

            if ( !_parse_log_status ) {
                    output += "log incomplete; ";
                    LOG(logERROR,*pLog) << "LOG parsing failed" << flush;
                    cout << *pLog;
                    jres.setOutput( output ); 
                    jres.setStatus(Job::FAILED);
                    delete _qmpackage;
                    return jres;
            } 
            
            _parse_orbitals_status = _qmpackage->ParseOrbitalsFile( &_orbitalsAB );

            if ( !_parse_orbitals_status ) {
                    output += "fort7 failed; " ;
                    LOG(logERROR,*pLog) << "GAUSSIAN orbitals (fort.7) parsing failed" << flush;
                    cout << *pLog;
                    jres.setOutput( output ); 
                    jres.setStatus(Job::FAILED);
                    delete _qmpackage;
                    return jres;
            } 
    } // end of the parse orbitals/log
        

    string _orbitals_storage_dir = idft_work_dir + "/" + "pairs/" + frame_dir;
   // orbital file used to archive parsed data
    string _pair_file = ( format("%1%%2%%3%%4%%5%") % "pair_" % ID_A % "_" % ID_B % ".orb" ).str();
   ub::matrix<double> _JAB;
   
   // trim virtual orbitals if too many are given
   if ( _do_trim ) {

       if ( !_do_parse ) { // orbitals must be loaded from a file
           LOG(logDEBUG,*pLog) << "Loading orbitals from " << _pair_file << flush;    
           std::ifstream ifs( (_orbitals_storage_dir + "/" + _pair_file).c_str() );
           boost::archive::binary_iarchive ia( ifs );
           ia >> _orbitalsAB;
           ifs.close();
       }     
       
       LOG(logDEBUG,*pLog) << "Trimming dimer virtual orbitals from " 
               << _orbitalsAB.getNumberOfLevels() - _orbitalsAB.getNumberOfElectrons() << " to " 
               << _orbitalsAB.getNumberOfElectrons()*(_trim_factor-1) << flush;   
       
       _orbitalsAB.Trim(_trim_factor);
   }
   
   if ( _do_project ) {
       
       if ( !_do_parse ) { // orbitals must be loaded from a file
           LOG(logDEBUG,*pLog) << "Loading orbitals from " << _pair_file << flush;    
           std::ifstream ifs( (_orbitals_storage_dir + "/" + _pair_file).c_str() );
           boost::archive::binary_iarchive ia( ifs );
           ia >> _orbitalsAB;
           ifs.close();
       }
       
       Orbitals _orbitalsA;
       Orbitals _orbitalsB;
 
        // load the corresponding monomer orbitals and prepare the dimer guess 
        string orb_file_A = (format("%1%_%2%%3%") % "molecule" % ID_A % ".orb").str();
        string DIR_A  = edft_work_dir + "/" + "molecules/" + frame_dir;
        std::ifstream ifs_A( (DIR_A +"/" + orb_file_A).c_str() );
        LOG(logDEBUG,*pLog) << "Loading orbitals from " << orb_file_A << flush;   
        boost::archive::binary_iarchive ia_A( ifs_A );
        ia_A >> _orbitalsA;
        ifs_A.close();

        string orb_file_B = (format("%1%_%2%%3%") % "molecule" % ID_B % ".orb").str();
        LOG(logDEBUG,*pLog) << "Loading orbitals from " << orb_file_B << flush;    
        string DIR_B  = edft_work_dir + "/" + "molecules/" + frame_dir;
        std::ifstream ifs_B( (DIR_B +"/" + orb_file_B).c_str() );
        boost::archive::binary_iarchive ia_B( ifs_B );
        ia_B >> _orbitalsB;
        ifs_B.close();
     
        if ( _do_trim ) {
             LOG(logDEBUG,*pLog) << "Trimming molecule A virtual orbitals from " 
                    << _orbitalsA.getNumberOfLevels() - _orbitalsA.getNumberOfElectrons() << " to " 
                    << _orbitalsA.getNumberOfElectrons()*(_trim_factor-1) << flush;  
            
            _orbitalsA.Trim(_trim_factor);
            
            LOG(logDEBUG,*pLog) << "Trimming molecule B virtual orbitals from " 
                    << _orbitalsB.getNumberOfLevels() - _orbitalsB.getNumberOfElectrons() << " to " 
                    << _orbitalsB.getNumberOfElectrons()*(_trim_factor-1) << flush;              
            _orbitalsB.Trim(_trim_factor);
        }
     
       _calculate_integrals = CalculateIntegrals( &_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB, opThread );

        if ( !_calculate_integrals ) {
                output += "integrals failed; " ;
                LOG(logERROR,*pLog) << "Calculating integrals failed" << flush;
                cout << *pLog;
                jres.setOutput( output ); 
                jres.setStatus(Job::FAILED);
                return jres;
        } 
    
        int HOMO_A = _orbitalsA.getNumberOfElectrons() ;
        int HOMO_B = _orbitalsB.getNumberOfElectrons() ;
    
        int LUMO_A = HOMO_A + 1;
        int LUMO_B = HOMO_B + 1;
    
        double J_h = getCouplingElement( HOMO_A , HOMO_B, &_orbitalsA, &_orbitalsB, &_JAB );
        double J_e = getCouplingElement( LUMO_A , LUMO_B, &_orbitalsA, &_orbitalsB, &_JAB );
    
        LOG(logINFO,*pLog) << "Couplings h/e " << ID_A << ":" << ID_B << " " << J_h  << ":" << J_e  << flush; 
    
        //qmpair->setJeff2( J_h,  1 );
        //qmpair->setJeff2( J_e, -1 );
       
        // Output the thread run summary and clean the Logger
        LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating pair " << ID_A << ":" << ID_B << flush; 
        cout << *pLog;

       // save orbitals 
       boost::filesystem::create_directories(_orbitals_storage_dir);  

       LOG(logDEBUG,*pLog) << "Saving orbitals to " << _pair_file << flush;
       std::ofstream ofs( (_orbitals_storage_dir + "/" + _pair_file).c_str() );
       boost::archive::binary_oarchive oa( ofs );

       if ( !( _store_orbitals && _do_parse && _parse_orbitals_status) )   _store_orbitals = false;
       if ( !( _store_overlap && _do_parse && _parse_log_status) )    _store_overlap = false;
       if ( !( _store_integrals && _do_project && _calculate_integrals) )  { _store_integrals = false; _orbitalsAB.setIntegrals( &_JAB ) ; }

       _orbitalsAB.setStorage( _store_orbitals, _store_overlap, _store_integrals );
       oa << _orbitalsAB;
       ofs.close();
   
    // save project summary    
    /* <pair idA="" idB="" typeA="" typeB="">
     *          <overlap orbA="" orbB="" enA="" enB="" ></overlap>
     * </pair>
     * 
     */
    Property _job_summary;
        Property *_pair_summary = &_job_summary.add("pair","");
         string nameA = seg_A->getName();
         string nameB = seg_B->getName();
        _pair_summary->setAttribute("idA", ID_A);
        _pair_summary->setAttribute("idB", ID_B);
        _pair_summary->setAttribute("homoA", HOMO_A);
        _pair_summary->setAttribute("homoB", HOMO_B);
        _pair_summary->setAttribute("typeA", nameA);
        _pair_summary->setAttribute("typeB", nameB);
        for (int levelA = HOMO_A - _max_occupied_levels +1; levelA <= LUMO_A + _max_unoccupied_levels; ++levelA ) {
                for (int levelB = HOMO_B - _max_occupied_levels + 1; levelB <= HOMO_B + _max_unoccupied_levels ; ++levelB ) {        
                        Property *_overlap_summary = &_pair_summary->add("overlap",""); 
                        double JAB = getCouplingElement( levelA , levelB, &_orbitalsA, &_orbitalsB, &_JAB );
                        double energyA = _orbitalsA.getEnergy( levelA );
                        double energyB = _orbitalsB.getEnergy( levelB );
                        _overlap_summary->setAttribute("orbA", levelA);
                        _overlap_summary->setAttribute("orbB", levelB);
                        _overlap_summary->setAttribute("jAB", JAB);
                        _overlap_summary->setAttribute("eA", energyA);
                        _overlap_summary->setAttribute("eB", energyB);
                }
        }
    
        sout <<  setlevel(1) << _job_summary;
   } // end of the projection loop

   // cleanup whatever is not needed
   _qmpackage->CleanUp();
   delete _qmpackage;
   
    jres.setOutput( sout.str() );   
    jres.setStatus(Job::COMPLETE);
    
    return jres;
}


void IDFT::PrepareGuess( Orbitals* _orbitalsA, Orbitals* _orbitalsB,
    Orbitals* _orbitalsAB, QMThread *opThread ) {
    
    LOG(logDEBUG,*opThread->getLogger())  << "Constructing the guess for the dimer orbitals" << flush;   
   
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
       
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    int _electronsA = _orbitalsA->getNumberOfElectrons();
    int _electronsB = _orbitalsB->getNumberOfElectrons();
    
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
    
    ub::matrix<double>* _mo_coefficients = _orbitalsAB->getOrbitals();    
    //cout << "MO coefficients " << *_mo_coefficients << endl;
    
    // AxB = | A 0 |  //   A = [EA, EB]  //
    //       | 0 B |  //                 //
    _mo_coefficients->resize( _levelsA + _levelsB, _basisA + _basisB  );
    _orbitalsAB->setBasisSetSize( _basisA + _basisB );
    _orbitalsAB->setNumberOfLevels( _electronsA - _electronsB , 
                                    _levelsA + _levelsB - _electronsA - _electronsB );
    _orbitalsAB->setNumberOfElectrons( _electronsA + _electronsB );
    
    ub::project( *_mo_coefficients, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( *_mo_coefficients, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( *_mo_coefficients, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    ub::project( *_mo_coefficients, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals();   

    //cout << "MO coefficients " << *_mo_coefficients << endl;
    
    ub::vector<double>* _energies = _orbitalsAB->getEnergies();
    _energies->resize( _levelsA + _levelsB );
     
    ub::project( *_energies, ub::range (0, _levelsA ) ) = *_orbitalsA->getEnergies();
    ub::project( *_energies, ub::range (_levelsA, _levelsA + _levelsB ) ) = *_orbitalsB->getEnergies();
    
    //cout << "MO energies " << *_energies << endl;
    
    ///"... ... Have now " >> _energies->size() >> " energies\n" >> *opThread;

}   

}};
  