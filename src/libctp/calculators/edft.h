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


#ifndef _CALC_DFT_ENERGIES_H
#define	_CALC_DFT_ENERGIES_H

#include <votca/ctp/segment.h>
#include <votca/ctp/gaussian.h>
#include <votca/ctp/orbitals.h>
#include <votca/ctp/parallelxjobcalc.h>
#include <unistd.h>
#include <fstream>
#include <sys/stat.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/filesystem.hpp>

namespace votca { namespace ctp {

/**
* \brief Site energies and orbitals for QM pairs
*
* QM orbitals and energies for all molecules
* Requires first-principles package, i.e. GAUSSIAN installation
*
* Callname: idft
*/

class EDFT : public ParallelXJobCalc< vector<Job*>, Job*, Job::JobResult >
{
public:

    EDFT() {};
   ~EDFT() {};

    string  Identify() { return "edft"; }
    void    Initialize(Topology *top, Property *options);
    void    ParseOrbitalsXML(Topology *top, Property *options);
    Job::JobResult EvalJob(Topology *top, Job *job, QMThread *thread);

    void    CleanUp();


private:

    //bool   _maverick;
    string _outParent;
    string _outMonDir;
    
    string _package;
    Property _package_options;   

};

void EDFT::CleanUp() {

}

void EDFT::Initialize(Topology *top, Property *options) {

    cout << endl << "... ... Initialize with " << _nThreads << " threads.";
    _maverick = (_nThreads == 1) ? true : false;
    
    string key = "options." + Identify();
    _jobfile = options->get(key + ".control.job_file").as<string>();
    /* ---- OPTIONS.XML Structure -----
     *
     * <edft>
     *      <control>
     *           <job_file>edft.jobs<job_file>
     *      </control>
     *      <package>gaussian.xml</package>
     *
     * </edft>
     *
     */

    this->ParseOrbitalsXML(top, options);

}


void EDFT::ParseOrbitalsXML(Topology *top, Property *opt) {

    string key = "options.edft";
    string _package_xml = opt->get(key+".package").as<string> ();
    //cout << endl << "... ... Parsing " << _package_xml << endl ;

    load_property_from_xml( _package_options, _package_xml.c_str() );    
    
    key = "package";
    _package = _package_options.get(key+".name").as<string> ();
    
    //cout << endl << "... ... Using " << _package << " package" << endl ;    

}


Job::JobResult EDFT::EvalJob(Topology *top, Job *job, QMThread *opThread) {

    bool _run_status;
    bool _parse_log_status;
    bool _parse_orbitals_status;
    
    Segment *seg = top->getSegment(job->getId());
    
    Logger* pLog = opThread->getLogger();
    pLog->setReportLevel(logDEBUG);
    pLog->setMultithreading( _maverick );
    
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating site " << seg->getId() << endl; 
    
    FILE *out;
    Orbitals _orbitals;
    vector < Segment* > segments;    
    
    segments.push_back( seg );

    // log, com, fort 7 files will be stored in ORB_FILES/gaussian/frame_x/mol_ID/
    // extracted information will be stored in  ORB_FILES/molecules/frame_x/molecule_ID.orb
    
    string edft_work_dir = "OR_FILES";
    string gaussian_work_dir = "gaussian";
    string orbitals_storage_dir = "molecules";
    string frame_dir =  "frame_" + boost::lexical_cast<string>(top->getDatabaseId());      
    
    string ID   = boost::lexical_cast<string>( seg->getId() );
    
    string DIR  = edft_work_dir + "/" + gaussian_work_dir + "/" + frame_dir + "/mol_" + ID;
    string ORB_DIR = edft_work_dir + "/" + orbitals_storage_dir + "/" + frame_dir;

    // GAUSSIAN filenames
    string fileName = "monomer";
    string XYZ_FILE = fileName + ".xyz";
    string COM_FILE = fileName + ".com";
    string LOG_FILE = fileName + ".log"; 
    string SHL_FILE = fileName + ".sh";
    string GAUSSIAN_ORB_FILE = "fort.7" ;

    // orbital file used to archive parsed data
    string ORB_FILE = "molecule_" + ID + ".orb";

    
    boost::filesystem::create_directories(DIR);     
    boost::filesystem::create_directories(ORB_DIR);        
    
    //out = fopen((DIR + "/" + XYZ_FILE).c_str(),"w");
    //seg->WriteXYZ(out);
    //fclose(out);
 
   if ( _package == "gaussian" ) { 
        
        Gaussian _gaussian( &_package_options );
        _gaussian.setLog( pLog );       
        _gaussian.setRunDir( DIR );
        _gaussian.setInputFile( COM_FILE );

        // provide a separate scratch dir for every thread
        if ( ( _gaussian.getScratchDir() ).size() != 0 ) {
          _gaussian.setShellFile( SHL_FILE );
           string SCR_DIR  = _gaussian.getScratchDir() + "/mol_" + ID;
          _gaussian.setScratchDir( SCR_DIR );
          _gaussian.WriteShellScript ();
        } 
        
        _gaussian.WriteInputFile( segments );
        
        // Run the executable
        _run_status = _gaussian.Run( );
        
        // Collect information     
        _gaussian.setLogFile( DIR + "/" + LOG_FILE );
        _parse_log_status = _gaussian.ParseLogFile( &_orbitals );
        
        _gaussian.setOrbitalsFile( DIR + "/" + GAUSSIAN_ORB_FILE );
        _parse_orbitals_status = _gaussian.ParseOrbitalsFile( &_orbitals );
 
        // save orbitals
        std::ofstream ofs( (ORB_DIR + "/" + ORB_FILE).c_str() );
        boost::archive::binary_oarchive oa( ofs );
        oa << _orbitals;
        ofs.close();
        
        /*
        std::vector< QMAtom* >* _atoms = _orbitals.getAtoms();
        int _nat = _atoms->size() ;
        LOG(logDEBUG,*pLog) << "Serializing " << _nat << " atoms" << endl; 
        */
       
        //Orbitals _orbitals_new;
        //std::ifstream ifs( (DIR +"/" + ORB_FILE).c_str() );
        //boost::archive::binary_iarchive ia( ifs );
        //ia >> _orbitals_new;
        //ifs.close();

        //_atoms = _orbitals_new.getAtoms();
        //cout << "Deserializing " << _atoms->size() << " atoms" << endl;
        //cout << (_atoms->front())->type << " ";
        //cout << (_atoms->front())->x << " ";
        //cout << (_atoms->front())->y << " ";
        //cout << (_atoms->front())->z << " ";
        //cout << (_atoms->front())->charge << endl;
        //exit(0);
        
        //cout << "ELECTRONS " << _orbitals_new.getNumberOfElectrons() << endl ;
        //cout << "BASIS SIZE " << _orbitals_new.getBasisSetSize() << endl ;
 
        _gaussian.CleanUp();
        
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
        
    // output of the JOB 
    jres.setOutput( output );

    // dump the LOG
    cout << *pLog;
    
    return jres;

   }    
   
}

}}

#endif	/* _CALC_DFT_ENERGIES_H */
