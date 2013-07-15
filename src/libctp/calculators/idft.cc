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
#include <votca/ctp/eigenvalues.h>
#include <votca/ctp/logger.h>
#include <iostream>
#include <c++/4.7/bits/stl_vector.h>

namespace votca { namespace ctp {
    namespace ub = boost::numeric::ublas;
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize(ctp::Topology *top, tools::Property* options ) {
    
    _energy_difference = 0.0;
    
    ParseOptionsXML( options  );

}

    
void IDFT::ParseOptionsXML( tools::Property *opt ) {
   
    // Orbitals are in fort.7 file; number of electrons in .log file
    
    string key = "options." + Identify();
    if ( opt->exists( key + ".degeneracy" ) ) {
        _energy_difference = opt->get( key + ".degeneracy" ).as< double > ();
    }
    else {
        cout << "... ... NOT treating degenerate orbitals\n" ;
    }    

    _jobfile = opt->get(key + ".control.job_file").as<string>();

    string _package_xml = opt->get(key+".package").as<string> ();
    //cout << endl << "... ... Parsing " << _package_xml << endl ;

    load_property_from_xml( _package_options, _package_xml.c_str() );    
    
     key = "package";
    _package = _package_options.get(key+".name").as<string> ();

    
    /* --- ORBITALS.XML Structure ---
     * <options>
     *   <idft>
     *     <orbitals_A>fort.7</orbitals_A>
     *     <orbitals_B>fort.7</orbitals_B>
     *     <orbitals_AB>fort.7</orbitals_AB>
     *     <overlap_AB>dimer.log</overlap_AB>
     *   </idft>
     * </options>
     */

}

/*
class _inv_sqrt {
public:
  double operator()(double x) { return 1./x; }
};
*/

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

void IDFT::CalculateIntegrals(Orbitals* _orbitalsA, Orbitals* _orbitalsB, 
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
    
    ///if ( tools::globals::verbose ) *opThread << "\n... ... Calculating electronic couplings \n" ;
    
    // constructing the direct product orbA x orbB
    int _basisA = _orbitalsA->getBasisSetSize();
    int _basisB = _orbitalsB->getBasisSetSize();
    
    
    //cout << "... ... Basis [molA:molB] " << _basisA << ":" << _basisB << endl;
    
    int _levelsA = _orbitalsA->getNumberOfLevels();
    int _levelsB = _orbitalsB->getNumberOfLevels();
    
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
        
    //cout << zeroB << endl;
    //cout << zeroA << endl;
    
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    
    // AxB = | A 0 |  //
    //       | 0 B |  //      
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals();    
    //cout << "_psi_AxB: " << _psi_AxB << endl;
    
    // Fock matrix of a dimer   
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB->getNumberOfLevels(), (*_orbitalsAB->getEnergies()).data() ); 

    // psi_AxB * S_AB * psi_AB
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
    ///if ( tools::globals::verbose ) *opThread << "... ... Done calculating electronic couplings\n";
       
    //cout << JAB_dimer.at_element( HOMO_A , HOMO_B + _levelsA ) * conv_Hrt_eV << endl; 
    //cout << JAB_dimer.at_element(_levelsA + HOMO_B, HOMO_A ) * conv_Hrt_eV << endl;

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

    bool _run_status;
    bool _parse_log_status;
    bool _parse_orbitals_status;
    
    QMPair *qmpair = NULL;
    Logger* pLog = opThread->getLogger();

    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    
    
    int _job_ID = job->getId();
    string sID = boost::lexical_cast<string>( _job_ID ); 
    
    string _job_tag = job->getTag();
    
    Tokenizer _tok ( _job_tag, ":" ); 
    vector<string> _mol_ids;
    _tok.ToVector( _mol_ids );
    
    string sID_A = _mol_ids[0];
    string sID_B = _mol_ids[1];

    int ID_A   = boost::lexical_cast<int>( sID_A );
    int ID_B   = boost::lexical_cast<int>( sID_B );
    
    Segment *seg_A = top->getSegment( ID_A );
    Segment *seg_B = top->getSegment( ID_B );
    qmpair = nblist.FindPair( seg_A, seg_B );

    
    LOG(logDEBUG,*pLog) << TimeStamp() << " Evaluating pair "  
            << sID << " ["  << ID_A << ":" << ID_B << "] out of " << 
           (top->NBList()).size()  << flush; 
    
    FILE *out;
    vector < Segment* > segments;
    segments.push_back( qmpair->Seg1() );
    segments.push_back( qmpair->Seg2() );
    
    ub::matrix<double> _JAB;
    
    Orbitals _orbitalsA;
    Orbitals _orbitalsB;
    Orbitals _orbitalsAB;

    string idft_work_dir = "OR_FILES";
    string edft_work_dir = "OR_FILES";
    string gaussian_work_dir = "gaussian";
    string orbitals_storage_dir = "pairs";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    
    string GAUSS_DIR  = idft_work_dir + "/" + gaussian_work_dir + "/" + frame_dir + "/" + "pair_"  + sID_A + "_" + sID_B;
    string ORBIT_DIR = idft_work_dir + "/" + orbitals_storage_dir + "/" + frame_dir;
    
    // orbital file used to archive parsed data
    string ORB_FILE = "pair_" + sID_A + "_" + sID_B + ".orb";
    
    boost::filesystem::create_directories(GAUSS_DIR);     
    boost::filesystem::create_directories(ORBIT_DIR);        
    
    string fileName = "pair_"  + sID_A + "_" + sID_B ;
    //string DIR  = _outParent + "/" + "pair_"  + ID_A + "_" + ID_B;
    string XYZ_FILE = fileName + ".xyz";
    //string ORB_FILE = fileName + ".orb";
    string COM_FILE = fileName + ".com";
    string LOG_FILE = fileName + ".log"; 
    string SHL_FILE = fileName + ".sh";
    string GAUSSIAN_ORB_FILE = "fort.7" ;
        
    out = fopen((GAUSS_DIR + "/" + XYZ_FILE).c_str(),"w");
    qmpair->WriteXYZ(out);
    fclose(out);
 
    string ORB_FILE_A = "molecule_" + sID_A + ".orb";
    string ORB_FILE_B = "molecule_" + sID_B + ".orb";
    
    //string ORB_FILE_A = "monomer.orb";
    //string ORB_FILE_B = "monomer.orb";
   
    string DIR_A  = edft_work_dir + "/" + "molecules/" + frame_dir;
    //cout << "... ... " + DIR_A +"/" + ORB_FILE_A + "\n";
    std::ifstream ifs_A( (DIR_A +"/" + ORB_FILE_A).c_str() );
    boost::archive::binary_iarchive ia_A( ifs_A );
    ia_A >> _orbitalsA;
    ifs_A.close();
    //LOG(logDEBUG,*pLog) << "Basis size A " << _orbitalsA.getBasisSetSize() << flush;
    //LOG(logDEBUG,*pLog) << "Number of Levels A " << _orbitalsA.getNumberOfLevels() << flush; 
    
    
    string DIR_B  = edft_work_dir + "/" + "molecules/" + frame_dir;
    //cout << "... ... " << DIR_B +"/" + ORB_FILE_B << "\n";
    std::ifstream ifs_B( (DIR_B +"/" + ORB_FILE_B).c_str() );
    boost::archive::binary_iarchive ia_B( ifs_B );
    ia_B >> _orbitalsB;
    ifs_B.close();
    //LOG(logDEBUG,*pLog) << "Basis size B " << _orbitalsB.getBasisSetSize() << flush;
    //LOG(logDEBUG,*pLog) << "Number of Levels B " << _orbitalsB.getNumberOfLevels() << flush; 
      
   if ( _package == "gaussian" ) { 
        
        Gaussian _gaussian( &_package_options );
        
        _gaussian.setLog( pLog );       
        _gaussian.setRunDir( GAUSS_DIR );
        _gaussian.setInputFile( COM_FILE );

        //cout << GAUSS_DIR << endl;
        //cout << COM_FILE << endl;
        
        // provide a separate scratch dir for every thread
        if ( ( _gaussian.getScratchDir() ).size() != 0 ) {
          _gaussian.setShellFile( SHL_FILE );
           string SCR_DIR  = _gaussian.getScratchDir() + "/pair_" + sID;
          _gaussian.setScratchDir( SCR_DIR );
          _gaussian.WriteShellScript ();
        } 
        
        // in case we do not want to do an SCF loop for a dimer
        if ( _gaussian.GuessRequested() ) {
            
            LOG(logDEBUG,*pLog) << "Preparing the guess" << flush;
            
            PrepareGuess(&_orbitalsA, &_orbitalsB, &_orbitalsAB, opThread);
            
            //cout << *_orbitalsAB.getEnergies()  << endl;
            //cout << "Guess address " << &_orbitalsAB << endl;
            //Orbitals *_pOrbitalsAB = &_orbitalsAB;
            //cout << "Guess address " << _pOrbitalsAB << endl;
           
            //Orbitals* test;
            //test = NULL;
            //cout << "NULL address " << test << endl;          
            
            _gaussian.WriteInputFile( segments, &_orbitalsAB );
        } else {
            _gaussian.WriteInputFile( segments );
        }
        
        // Run the executable
        _run_status = _gaussian.Run( );

        // Collect information     
        _gaussian.setLogFile( GAUSS_DIR + "/" + LOG_FILE );
        _parse_log_status = _gaussian.ParseLogFile( &_orbitalsAB );
        
        _gaussian.setOrbitalsFile( GAUSS_DIR + "/" + GAUSSIAN_ORB_FILE );
        _parse_orbitals_status = _gaussian.ParseOrbitalsFile( &_orbitalsAB );
 
        // save orbitals 
        std::ofstream ofs( (ORBIT_DIR + "/" + ORB_FILE).c_str() );
        boost::archive::binary_oarchive oa( ofs );
        oa << _orbitalsAB;
        ofs.close();
        
        _gaussian.CleanUp();
        
   }      

    CalculateIntegrals( &_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB, opThread );
     
    int HOMO_A = _orbitalsA.getNumberOfElectrons() ;
    int HOMO_B = _orbitalsB.getNumberOfElectrons() ;
    
    int LUMO_A = _orbitalsA.getNumberOfElectrons() + 1;
    int LUMO_B = _orbitalsB.getNumberOfElectrons() + 1;
    
    LOG(logINFO,*pLog) << "Coupling h/e " << ID_A << ":" << ID_B << " " 
         << getCouplingElement( HOMO_A , HOMO_B, &_orbitalsA, &_orbitalsB, &_JAB ) << " "
         << getCouplingElement( LUMO_A , LUMO_B, &_orbitalsA, &_orbitalsB, &_JAB ) << flush; 
    
    qmpair->setJeff2( getCouplingElement( HOMO_A , HOMO_B, &_orbitalsA, &_orbitalsB, &_JAB ),  1 );
    qmpair->setJeff2( getCouplingElement( LUMO_A , LUMO_B, &_orbitalsA, &_orbitalsB, &_JAB ), -1 );
    
    // Output the thread run summary and clean the thread
    cout << *pLog;
 
     // GENERATE OUTPUT AND FORWARD TO PROGRESS OBSERVER (RETURN)
    Job::JobResult jres = Job::JobResult();
    string output = "GAUSSIAN: ";
    jres.setStatus(Job::COMPLETE);
    
    if ( !_run_status ) {
        output += "run failed; " ;
        LOG(logERROR,*pLog) << "GAUSSAIN run failed" << flush;
        jres.setStatus(Job::FAILED);
    } else {
        output += "run completed; " ;
    }
    
    if ( !_parse_log_status ) {
        output += "log incomplete; ";
        LOG(logERROR,*pLog) << "GAUSSIAN log incomplete" << flush;
        jres.setStatus(Job::FAILED);
    } else {
        output += "log parsed; " ;
    }

    if ( !_parse_orbitals_status ) {
        output += "fort7 failed; " ;
        LOG(logERROR,*pLog) << "GAUSSIAN orbitals (fort.7) not parsed" << flush;
    } else {
        output += "orbitals parsed; " ;
    }

    if ( _run_status && _parse_log_status &&  _parse_orbitals_status ) {
        output += "Coupling h/e " + sID_A + ":" + sID_B;
    }
    
    LOG(logINFO,*pLog) << TimeStamp() << " Finished evaluating pair " << sID_A << ":" << sID_B << flush; 
        
    // output of the JOB 
    jres.setOutput( output );   
    
    return Job::JobResult();
}

void IDFT::PrepareGuess( Orbitals* _orbitalsA, Orbitals* _orbitalsB,
    Orbitals* _orbitalsAB, QMThread *opThread ) {
    
    //LOG(logDEBUG,*opThread->getLogger())  << "... ... Constructing the guess for the dimer orbitals\n" ;   
   
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
