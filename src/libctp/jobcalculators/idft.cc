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


//#ifdef MKLROOT
//#include <mkl_boost_ublas_matrix_prod.hpp>
//#endif
         
#include <votca/ctp/logger.h>
#include <votca/ctp/qmpackagefactory.h>
#include <votca/ctp/overlap.h>

using boost::format;
using namespace boost::filesystem;
using namespace votca::tools;

namespace ub = boost::numeric::ublas;
    
namespace votca { namespace ctp {
    
// +++++++++++++++++++++++++++++ //
// IDFT MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IDFT::Initialize(votca::tools::Property* options ) {
    
    _energy_difference = 0.0;
    
    _do_input = false;
    _do_run = false;
    _do_parse = false;
    _do_project = false;
    _do_trim = false;
    
    _store_orbitals = false;
    _store_overlap = false;
    _store_integrals = false;
    
    cout << _options;
    
    ParseOptionsXML( options  );
    
    // register all QM packages (Gaussian, turbomole, etc))
    QMPackageFactory::RegisterAll();

}

void IDFT::ParseOptionsXML( votca::tools::Property *opt ) {
   
    // Orbitals are in fort.7 file; number of electrons in .log file
    
    string key = "options." + Identify();
    _energy_difference = opt->get( key + ".degeneracy" ).as< double > ();
    
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
    
    key = "options." + Identify() +".job";
    _jobfile = opt->get(key + ".file").as<string>();    

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
    Property _job_input = job->getInput();  
    list<Property*> lSegments = _job_input.Select( "segment" );    
    int ID_A   = lSegments.front()->getAttribute<int>( "id" );
    string type_A = lSegments.front()->getAttribute<string>( "type" );
    int ID_B   = lSegments.back()->getAttribute<int>( "id" );
    string type_B = lSegments.back()->getAttribute<string>( "type" );
    
    Segment *seg_A = top->getSegment( ID_A );   
    assert( seg_A->getName() == type_A );
    
    Segment *seg_B = top->getSegment( ID_B );
    assert( seg_B->getName() == type_B );
    
    vector < Segment* > segments;
    segments.push_back( seg_A );
    segments.push_back( seg_B );
    
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());     

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
            LOG(logDEBUG,*pLog) << "Loading " << orb_file_A << flush;   
            boost::archive::binary_iarchive ia_A( ifs_A );
            ia_A >> _orbitalsA;
            ifs_A.close();

            string orb_file_B = (format("%1%_%2%%3%") % "molecule" % ID_B % ".orb").str();
            LOG(logDEBUG,*pLog) << "Loading " << orb_file_B << flush;    
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
    
    
    // This will be later used to write orbitals of the dimer to a file 
    // SOMETHING TO CLEANUP
    Orbitals _orbitalsAB; 
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
   
   /* trimming DIMER orbitals is not giving accurate results 
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
       
       //_orbitalsAB.Trim(_trim_factor);
   } */
   
   Property _job_summary;
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
             LOG(logDEBUG,*pLog) << "Trimming virtual orbitals A:" 
                    << _orbitalsA.getNumberOfLevels() - _orbitalsA.getNumberOfElectrons() << "->" 
                    << _orbitalsA.getNumberOfElectrons()*(_trim_factor-1);  
            
            _orbitalsA.Trim(_trim_factor);
            
            LOG(logDEBUG,*pLog) << " B:" 
                    << _orbitalsB.getNumberOfLevels() - _orbitalsB.getNumberOfElectrons() << "->" 
                    << _orbitalsB.getNumberOfElectrons()*(_trim_factor-1) << flush;              
            _orbitalsB.Trim(_trim_factor);
        }
     
        Overlap _overlap; 
        _overlap.setLogger(pLog);
                
       _calculate_integrals = _overlap.CalculateIntegrals( &_orbitalsA, &_orbitalsB, &_orbitalsAB, &_JAB );

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
    
        double J_h = _overlap.getCouplingElement( HOMO_A , HOMO_B, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference );
        double J_e = _overlap.getCouplingElement( LUMO_A , LUMO_B, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference );
    
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
       if ( !( _store_integrals && _do_project && _calculate_integrals) )  {
           _store_integrals = false; 
       } else {
           _orbitalsAB.setIntegrals( &_JAB );
       }

       _orbitalsAB.setStorage( _store_orbitals, _store_overlap, _store_integrals );

       
       oa << _orbitalsAB;
       ofs.close();
   
    // save project summary    
    /* <pair idA="" idB="" typeA="" typeB="">
     *          <overlap orbA="" orbB="" enA="" enB="" ></overlap>
     * </pair>
     * 
     */

        Property *_job_output = &_job_summary.add("output","");
        Property *_pair_summary = &_job_output->add("pair","");
         string nameA = seg_A->getName();
         string nameB = seg_B->getName();
        _pair_summary->setAttribute("idA", ID_A);
        _pair_summary->setAttribute("idB", ID_B);
        _pair_summary->setAttribute("homoA", HOMO_A);
        _pair_summary->setAttribute("homoB", HOMO_B);
        _pair_summary->setAttribute("typeA", nameA);
        _pair_summary->setAttribute("typeB", nameB);
        for (int levelA = HOMO_A - _max_occupied_levels +1; levelA <= LUMO_A + _max_unoccupied_levels - 1; ++levelA ) {
                for (int levelB = HOMO_B - _max_occupied_levels + 1; levelB <= LUMO_B + _max_unoccupied_levels -1 ; ++levelB ) {        
                        Property *_overlap_summary = &_pair_summary->add("overlap",""); 
                        double JAB = _overlap.getCouplingElement( levelA , levelB, &_orbitalsA, &_orbitalsB, &_JAB, _energy_difference );
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
   
    jres.setOutput( _job_summary );   
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

void IDFT::WriteJobFile(Topology *top) {

    cout << endl << "... ... Writing job file " << flush;
    ofstream ofs;
    ofs.open(_jobfile.c_str(), ofstream::out);
    if (!ofs.is_open()) throw runtime_error("\nERROR: bad file handle: " + _jobfile);

 
    QMNBList::iterator pit;
    QMNBList &nblist = top->NBList();    

    int jobCount = 0;
    if (nblist.size() == 0) {
        cout << endl << "... ... No pairs in neighbor list, skip." << flush;
        return;
    }    

    ofs << "<jobs>" << endl;    
    string tag = "";
    
    for (pit = nblist.begin(); pit != nblist.end(); ++pit) {

        int id1 = (*pit)->Seg1()->getId();
        string name1 = (*pit)->Seg1()->getName();
        int id2 = (*pit)->Seg2()->getId();
        string name2 = (*pit)->Seg2()->getName();   

        int id = ++jobCount;

        Property Input;
        Property *pInput = &Input.add("input","");
        Property *pSegment =  &pInput->add("segment" , (format("%1$s") % id1).str() );
        pSegment->setAttribute<string>("type", name1 );
        pSegment->setAttribute<int>("id", id1 );

        pSegment =  &pInput->add("segment" , (format("%1$s") % id2).str() );
        pSegment->setAttribute<string>("type", name2 );
        pSegment->setAttribute<int>("id", id2 );
        
        Job job(id, tag, Input, Job::AVAILABLE );
        job.ToStream(ofs,"xml");
        
    }

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
    cout << endl << "... ... In total " << jobCount << " jobs" << flush;
    
}


void IDFT::Import( Topology *top ) 
{
    Property xml;

    QMNBList &nblist = top->NBList();   
    int _number_of_pairs = nblist.size();
    int _current_pairs = 0;
    int _incomplete_jobs = 0;
    
    Logger _log;
    _log.setReportLevel(logINFO);
    

    string _idft_jobs_file;
    load_property_from_xml(xml, _idft_jobs_file);
    
    list<Property*> jobProps = xml.Select("jobs.job");
    list<Property*> ::iterator it;

    for (it = jobProps.begin(); it != jobProps.end(); ++it) {
 
        if ( (*it)->exists("output") && (*it)->exists("output.pair") ) {
            
            //cout << **it;
            
            Property poutput = (*it)->get("output.pair");
            
            int homoA = poutput.getAttribute<int>("homoA");
            int homoB = poutput.getAttribute<int>("homoB");
            
            int idA = poutput.getAttribute<int>("idA");
            int idB = poutput.getAttribute<int>("idB");
                       
            string typeA = poutput.getAttribute<string>("typeA");
            string typeB = poutput.getAttribute<string>("typeB");
            //cout << idA << ":" << idB << "\n"; 
            Segment *segA = top->getSegment(idA);
            Segment *segB = top->getSegment(idB);

            QMPair *qmp = nblist.FindPair(segA,segB);
            
            if (qmp == NULL) { // there is no pair in the neighbor list with this name
                //LOG(logINFO, _log) << "No pair " <<  idA << ":" << idB << " found in the neighbor list. Ignoring" << flush; 
            }   else {
                
                _current_pairs++;
                
                list<Property*> pOverlap = poutput.Select("overlap");
                list<Property*> ::iterator itOverlap;

                    
                    for (itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {

                        double energyA = (*itOverlap)->getAttribute<double>("eA");
                        double energyB = (*itOverlap)->getAttribute<double>("eB");
                        double overlapAB = (*itOverlap)->getAttribute<double>("jAB");
                        int orbA = (*itOverlap)->getAttribute<double>("orbA");
                        int orbB = (*itOverlap)->getAttribute<double>("orbB");

                        if ( orbA == homoA && orbB == homoB ) {
                                qmp->setJeff2(overlapAB*overlapAB, 1);
                                qmp->setIsPathCarrier(true, 1);
                        }

                        if ( orbA == homoA+1 && orbB == homoB+1 ) {
                                qmp->setJeff2(overlapAB*overlapAB, -1);
                                qmp->setIsPathCarrier(true, -1);
                        }
                    }    
            }
            
        } else {
            _incomplete_jobs++;
            LOG(logINFO, _log) << "Job " << (*it)->get( "id" ).as<string>() << " is " << (*it)->get( "status" ).as<string>() << endl;
        }
    }
    
    LOG(logINFO, _log) << "Pairs [total:updated] " <<  _number_of_pairs << ":" << _current_pairs << " Incomplete jobs: " << _incomplete_jobs << flush; 
    cout << _log;
}

/* SUPEREXCHANGE 

void IImport::FromIDFTWithSuperExchange(Topology *top, string &_idft_jobs_file) {

    Property xml;

    vector<Property*> records;
            
    QMNBList &nblist = top->NBList();
    int _number_of_pairs = nblist.size();
    int _current_pairs = 0;
    int _incomplete_jobs = 0;
    
    Logger _log;
    _log.setReportLevel(logINFO);
    
    //generate lists of bridges for superexchange pairs
    nblist.GenerateSuperExchange();

    // load the QC results in a vector indexed by the pair ID
    load_property_from_xml(xml, _idft_jobs_file);
    list<Property*> jobProps = xml.Select("jobs.job");
    
    records.resize( jobProps.size() + 1  );
    
    for (list<Property*> ::iterator  it = jobProps.begin(); it != jobProps.end(); ++it) {
 
        if ( (*it)->exists("output") && (*it)->exists("output.pair") ) {
            
            Property poutput = (*it)->get("output.pair");
            
            int idA = poutput.getAttribute<int>("idA");
            int idB = poutput.getAttribute<int>("idB");
                       
            Segment *segA = top->getSegment(idA);
            Segment *segB = top->getSegment(idB);

            QMPair *qmp = nblist.FindPair(segA,segB);
            
            if (qmp == NULL) { // there is no pair in the neighbor list with this name
                ;//LOG(logINFO, _log) << "No pair " <<  idA << ":" << idB << " found in the neighbor list. Ignoring" << flush; 
            }   else {
                LOG(logINFO, _log) << "Store in record: " <<  idA << ":" << idB << flush; 
                records[qmp->getId()] = & ((*it)->get("output.pair"));
            }
        }
    } // finished loading from the file


    // loop over all pairs in the neighborlist
    for (QMNBList::iterator ipair = top->NBList().begin(); ipair != top->NBList().end(); ++ipair) {
        
        QMPair *pair = *ipair;
        Segment* segmentA = pair->Seg1PbCopy();
        Segment* segmentB = pair->Seg2PbCopy();
        
        double Jeff2_homo = 0;
        double Jeff2_lumo = 0;
        
        cout << "Processing pair " << segmentA->getId() << ":" << segmentB->getId() << endl;
        
        QMPair::PairType _ptype = pair->getType();
        Property* pair_property = records[ pair->getId() ];
 
        int homoA = pair_property->getAttribute<int>("homoA");
        int homoB = pair_property->getAttribute<int>("homoB");
       
        // If a pair is of a direct type 
        if ( _ptype == QMPair::Hopping ||  _ptype == QMPair::SuperExchangeAndHopping ) {
            cout << "Pair is hopping" << endl;
            list<Property*> pOverlap = pair_property->Select("overlap");
 
            for (list<Property*> ::iterator itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {

                double overlapAB = (*itOverlap)->getAttribute<double>("jAB");
                int orbA = (*itOverlap)->getAttribute<double>("orbA");
                int orbB = (*itOverlap)->getAttribute<double>("orbB");

                if ( orbA == homoA && orbB == homoB ) {
                    Jeff2_homo += overlapAB*overlapAB;
                }

                if ( orbA == homoA+1 && orbB == homoB+1 ) {
                    Jeff2_lumo += overlapAB*overlapAB;
                }
            }    
            
        }
        
        // if pair has bridges only
        if ( _ptype == QMPair::SuperExchange  ||  _ptype == QMPair::SuperExchangeAndHopping ) {
            
            list<Property*> pOverlap = pair_property->Select("overlap");
            
            // this is to select HOMO_A and HOMO_B 
            double overlapAB;
            int orbA;
            int orbB;
            double energyA;
            double energyB;
            
            for (list<Property*> ::iterator itOverlap = pOverlap.begin(); itOverlap != pOverlap.end(); ++itOverlap) {
                if ( orbA == homoA && orbB == homoB ) {  
                    overlapAB = (*itOverlap)->getAttribute<double>("jAB");
                    orbA = (*itOverlap)->getAttribute<double>("orbA");
                    orbB = (*itOverlap)->getAttribute<double>("orbB");
                    energyA = (*itOverlap)->getAttribute<double>("eA");
                    energyB = (*itOverlap)->getAttribute<double>("eB");
                }
            }
            
            
            
            // loop over the bridging segments
            for ( vector< Segment* >::iterator itBridge = pair->getBridgingSegments().begin() ; itBridge != pair->getBridgingSegments().end(); itBridge++ ) {

                Segment* Bridge = *itBridge;
                int IDBridge = Bridge->getId();

                // pairs from the bridge to the donor and acceptor
                QMPair* Bridge_A = nblist.FindPair( segmentA, Bridge );
                QMPair* Bridge_B = nblist.FindPair( segmentB, Bridge );

                Property* pBridge_A = records[ Bridge_A->getId() ];
                Property* pBridge_B = records[ Bridge_B->getId() ];

                list<Property*> pOverlapA = pBridge_A->Select("overlap");
                list<Property*> pOverlapB = pBridge_B->Select("overlap");

                // IDs of the Donor and Acceptor
                int IdA = segmentA->getId();
                int IdB = segmentB->getId();

                // IDs stored in the file
                int id1A = pBridge_A->getAttribute<int>("idA");
                int id2A = pBridge_A->getAttribute<int>("idB");

                int id1B = pBridge_B->getAttribute<int>("idA");
                int id2B = pBridge_B->getAttribute<int>("idB");

                // suffix for the donor and acceptor
                string suffixA = ( id1A == IDBridge ) ? "B" : "A"; // use "A" as a bridge 
                string suffixB = ( id1B == IDBridge ) ? "B" : "A"; // use "A" as a bridge 
                string suffixBridgeA = ( id1A == IDBridge ) ? "A" : "B";
                string suffixBridgeB = ( id1B == IDBridge ) ? "A" : "B";
                
                int homoBridgeA = pBridge_A->getAttribute<int>("orb" + suffixBridgeA );
                int homoBridgeB = pBridge_B->getAttribute<int>("orb" + suffixBridgeB );
                assert( homoBridgeA == homoBridgeB );
                int homoBridge = homoBridgeA;
               
                for (list<Property*> ::iterator itOverlapA = pOverlapA.begin(); itOverlapA != pOverlapA.end(); ++itOverlapA) {
                for (list<Property*> ::iterator itOverlapB = pOverlapB.begin(); itOverlapB != pOverlapB.end(); ++itOverlapB) {
                    
                    int orbDonor = (*itOverlapA)->getAttribute<int>( "orb" + suffixA );
                    int orbAcceptor = (*itOverlapB)->getAttribute<int>( "orb" + suffixB );
                    int orbBridgeA  = (*itOverlapA)->getAttribute<int>( "orb" + suffixBridgeA );
                    int orbBridgeB = (*itOverlapB)->getAttribute<int>( "orb" + suffixBridgeB );
                    
                    if (  orbDonor == homoA && orbAcceptor == homoB && orbBridgeA == orbBridgeB && orbBridgeA <= homoBridge) {
                        
                        double jDB = (*itOverlapA)->getAttribute<double>( "jAB" );
                        double jBA = (*itOverlapB)->getAttribute<double>( "jAB" );
                        double eA  = (*itOverlapA)->getAttribute<double>( "e" + suffixA );
                        double eB  = (*itOverlapB)->getAttribute<double>( "e" + suffixB );
                        
                        double eBridgeA  = (*itOverlapA)->getAttribute<double>( "e" + suffixBridgeA );
                        double eBridgeB  = (*itOverlapB)->getAttribute<double>( "e" + suffixBridgeB );
                        
                        //assert( eBridgeA - eBridgeB < 1e-50 );
                     
                        cout << homoA << " " << homoB << " " << (*itOverlapA)->getAttribute<int>( "orb" + suffixBridgeA )
                             << " JDB " << jDB 
                             << " JBA " << jBA << endl;
                        
                        // This in principle violates detailed balance. Any ideas?
                        Jeff2_homo += 0.5 * (jDB*jBA / (eA - eBridgeA) + jDB*jBA / (eB - eBridgeB));
                        
                                
                    }

                    if (  orbDonor == homoA+1 && orbAcceptor == homoB+1 && orbBridgeA == orbBridgeB && orbBridgeA > homoBridge) {
                        
                        double jDB = (*itOverlapA)->getAttribute<double>( "jAB" );
                        double jBA = (*itOverlapB)->getAttribute<double>( "jAB" );
                        double eA  = (*itOverlapA)->getAttribute<double>( "e" + suffixA );
                        double eB  = (*itOverlapB)->getAttribute<double>( "e" + suffixB );
                        
                        double eBridgeA  = (*itOverlapA)->getAttribute<double>( "e" + suffixBridgeA );
                        double eBridgeB  = (*itOverlapB)->getAttribute<double>( "e" + suffixBridgeB );
                        
                         // This in principle violates detailed balance. Any ideas?
                        Jeff2_lumo += 0.5 * (jDB*jBA / (eA - eBridgeA) + jDB*jBA / (eB - eBridgeB));
                        //jDB*jBA / (eB - eBridgeB);
                                
                    }
                    
                     
                }}
            } // end over bridges 
            
            
            
        } // end of if superexchange
         
        pair->setJeff2(Jeff2_homo, 1);
        pair->setIsPathCarrier(true, 1);
        
        pair->setJeff2(Jeff2_lumo, -1);
        pair->setIsPathCarrier(true, -1);
       
        break;
    }
                    
    LOG(logINFO, _log) << "Pairs [total:updated] " <<  _number_of_pairs << ":" << _current_pairs << " Incomplete jobs: " << _incomplete_jobs << flush; 
    cout << _log;
}


*/



}};
  

    
    /* OLD COMPREHENCIVE BUT SLOW WAY OF PROJECTING 
    // AxB = | A 0 |  //
    //       | 0 B |  //  
    ub::zero_matrix<double> zeroB( _levelsA, _basisB ) ;
    ub::zero_matrix<double> zeroA( _levelsB, _basisA ) ;
        
    //cout << zeroB << endl;
    //cout << zeroA << endl;
    
    ub::matrix<double> _psi_AxB ( _levelsA + _levelsB, _basisA + _basisB  );
    

     LOG(logDEBUG,*_pLog) << "Constructing direct product AxB [" 
            << _psi_AxB.size1() << "x" 
            << _psi_AxB.size2() << "]"<< flush;    
    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( _basisA, _basisA +_basisB ) ) = zeroB;
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( 0, _basisA ) ) = zeroA;    
    ub::project( _psi_AxB, ub::range (0, _levelsA ), ub::range ( 0, _basisA ) ) = *_orbitalsA->getOrbitals();
    ub::project( _psi_AxB, ub::range (_levelsA, _levelsA + _levelsB ), ub::range ( _basisA, _basisA + _basisB ) ) = *_orbitalsB->getOrbitals();    
    //cout << "_psi_AxB: " << _psi_AxB << endl;
    
    // psi_AxB * S_AB * psi_AB
    LOG(logDEBUG,*_pLog) << TimeStamp() << " Projecting the dimer onto monomer orbitals" << flush;    
    ub::matrix<double> _psi_AB = ub::prod( *_orbitalsAB->getOverlap(), ub::trans( *_orbitalsAB->getOrbitals() ) );          
    ub::matrix<double> _psi_AxB_dimer_basis = ub::prod( _psi_AxB, _psi_AB );
    LOG(logDEBUG,*_pLog) << TimeStamp() << " Done" << flush;    

    LOG(logDEBUG,*_pLog) << TimeStamp() << " Projecting the dimer onto monomer orbitals" << flush;   
    ub::matrix<double> _psi_AB = ub::prod( _psi_AxB, *_orbitalsAB->getOverlap() );          
    LOG(logDEBUG,*_pLog) << TimeStamp() << " Done" << flush;    

    // J = psi_AxB_dimer_basis * FAB * psi_AxB_dimer_basis^T
    LOG(logDEBUG,*_pLog) << TimeStamp() << " Projecting the Fock matrix onto the dimer basis" << flush;    
    ub::matrix<double> _temp = ub::prod( _fock_AB, ub::trans( _psi_AxB_dimer_basis ) ) ;
    ub::matrix<double> JAB_dimer = ub::prod( _psi_AxB_dimer_basis, _temp);
    LOG(logDEBUG,*_pLog) << TimeStamp() << " DOne" << flush;    

    // Fock matrix of a dimer   
    LOG(logDEBUG,*_pLog) << "Dimer Fock matrix [" 
            << _orbitalsAB->getNumberOfLevels() << "x" 
            << _orbitalsAB->getNumberOfLevels() << "]" << flush;    
    ub::diagonal_matrix<double> _fock_AB( _orbitalsAB->getNumberOfLevels(), (*_orbitalsAB->getEnergies()).data() ); 

    LOG(logDEBUG,*_pLog) << TimeStamp() << " Projecting the Fock matrix onto the dimer basis" << flush;    
    ub::matrix<double> _temp = ub::prod( _psi_AxB_dimer_basis, _fock_AB ) ;
    ub::matrix<double> JAB_dimer = ub::prod( _temp, ub::trans( _psi_AxB_dimer_basis ));
    LOG(logDEBUG,*_pLog) << TimeStamp() << " Done" << flush;    
    
    _temp.clear(); _fock_AB.clear();
     */
