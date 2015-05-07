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


#include "iexcitoncl.h"

#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <votca/tools/propertyiomanipulator.h>
#include <votca/ctp/logger.h>
#include <votca/ctp/xmapper.h>
#include <votca/ctp/xinteractor.h>

using boost::format;
using namespace boost::filesystem;
using namespace votca::tools;

namespace ub = boost::numeric::ublas;
    
namespace votca { namespace ctp {
    
// +++++++++++++++++++++++++++++ //
// IEXCITON MEMBER FUNCTIONS         //
// +++++++++++++++++++++++++++++ //

void IEXCITON::Initialize(votca::tools::Property* opt ) {
    

    _options = opt;
    
    cout << endl
         << "... ... Initialized with " << _nThreads << " threads. "
         << flush;

    _maverick = (_nThreads == 1) ? true : false;


    _induce= false;
    _singlet=true;
    string key = "options." + Identify();
    if ( opt->exists(key+".job_file")) {
        _jobfile = opt->get(key+".job_file").as<string>();
        }
    else {
            throw std::runtime_error("Job-file not set. Abort.");
    }
    if ( opt->exists(key+".mapping")) {
        _xml_file = opt->get(key+".mapping").as<string>();
        }
    else {
            throw std::runtime_error("Mapping-file not set. Abort.");
        }
    if ( opt->exists(key+".emp_file")) {
            _emp_file   = opt->get(key+".emp_file").as<string>();
        }
    else {
            throw std::runtime_error("Emp-file not set. Abort.");
        }
    if ( opt->exists(key+".spin")) {
            string type=opt->get(key+".spin").as<string>();
            if (type=="singlet") _singlet=true;
            else if (type=="triplet") _singlet=false;
            else std::runtime_error("Spin state not known. Abort.");
        }
    else {
            throw std::runtime_error("Spin not set. Abort.");
        }
    if ( opt->exists(key+".induce")) {
            _induce   = opt->get(key+".induce").as<bool>();
        }     
    
    cout << "done"<< endl;
}



void IEXCITON::PreProcess(Topology *top) {

    // INITIALIZE MPS-MAPPER (=> POLAR TOP PREP)
    cout << endl << "... ... Initialize MPS-mapper: " << flush;
    _mps_mapper.GenerateMap(_xml_file, _emp_file, top);
}


void IEXCITON::CustomizeLogger(QMThread *thread) {
    
    // CONFIGURE LOGGER
    Logger* log = thread->getLogger();
    log->setReportLevel(logDEBUG);
    log->setMultithreading(_maverick);

    log->setPreface(logINFO,    (format("\nT%1$02d INF ...") % thread->getId()).str());
    log->setPreface(logERROR,   (format("\nT%1$02d ERR ...") % thread->getId()).str());
    log->setPreface(logWARNING, (format("\nT%1$02d WAR ...") % thread->getId()).str());
    log->setPreface(logDEBUG,   (format("\nT%1$02d DBG ...") % thread->getId()).str());        
}



Job::JobResult IEXCITON::EvalJob(Topology *top, Job *job, QMThread *opThread) {
    
     // report back to the progress observer
    Job::JobResult jres = Job::JobResult();
    
    // get the logger from the thread
    Logger* pLog = opThread->getLogger();   
    
    // get the information about the job executed by the thread
    int _job_ID = job->getId();
    Property _job_input = job->getInput();  
    list<Property*> segment_list = _job_input.Select( "segment" );    
    int ID_A   = segment_list.front()->getAttribute<int>( "id" );
    string type_A = segment_list.front()->getAttribute<string>( "type" );
    string mps_fileA = segment_list.front()->getAttribute<string>( "mps_file" );
    int ID_B   = segment_list.back()->getAttribute<int>( "id" );
    string type_B = segment_list.back()->getAttribute<string>( "type" );
    string mps_fileB = segment_list.back()->getAttribute<string>( "mps_file" );

    
  
    Segment *seg_A = top->getSegment( ID_A );   
    assert( seg_A->getName() == type_A );
    
    Segment *seg_B = top->getSegment( ID_B );
    assert( seg_B->getName() == type_B );
    
    LOG(logINFO,*pLog) << TimeStamp() << " Evaluating pair "  
            << _job_ID << " ["  << ID_A << ":" << ID_B << "]" << flush; 
    
   vector<APolarSite*> seg_A_raw=_mps_mapper.GetOrCreateRawSites(mps_fileA);
   vector<APolarSite*> seg_B_raw=_mps_mapper.GetOrCreateRawSites(mps_fileB);
   
   PolarSeg* seg_A_polar=_mps_mapper.MapPolSitesToSeg(seg_A_raw,seg_A);
   PolarSeg* seg_B_polar=_mps_mapper.MapPolSitesToSeg(seg_B_raw,seg_B);
   
   PolarTop cleanup;
   vector<PolarSeg*> cleanup2;
   cleanup2.push_back(seg_A_polar);
   cleanup2.push_back(seg_B_polar);
   cleanup.setMM1(cleanup2);
  
   
   double JAB=EvaluatePair(top,seg_A_polar,seg_B_polar);
   
    
    
    
    
   
   Property _job_summary;
 

    Property *_job_output = &_job_summary.add("output","");
    Property *_pair_summary = &_job_output->add("pair","");
    _pair_summary->setAttribute("jAB", JAB);
     


    votca::tools::PropertyIOManipulator iomXML(votca::tools::PropertyIOManipulator::XML, 1, "");
    cout <<  iomXML << _job_summary;
    // end of the projection loop

 
   
    jres.setOutput( _job_summary );   
    jres.setStatus(Job::COMPLETE);
    
    return jres;
}

double IEXCITON::EvaluatePair(Topology *top,PolarSeg* Seg1,PolarSeg* Seg2){
    double int2eV = 1/(4*M_PI*8.854187817e-12) * 1.602176487e-19 / 1.000e-9;
    XInteractor actor;
    actor.ResetEnergy();
    Seg1->CalcPos();
    Seg2->CalcPos();
    vec s=top->PbShortestConnect(Seg1->getPos(),Seg2->getPos())-Seg1->getPos()+Seg2->getPos();
    PolarSeg::iterator pit1;
    PolarSeg::iterator pit2;
    double E=0.0;
    for (pit1=Seg1->begin();pit1<Seg1->end();++pit1){
        for (pit2=Seg2->begin();pit2<Seg2->end();++pit2){
            actor.BiasIndu(*(*pit1), *(*pit2),s);
            E += actor.E_f(*(*pit1), *(*pit2));               
    }
    }

    
  return E*int2eV;  
}


void IEXCITON::WriteJobFile(Topology *top) {

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
        if ((*pit)->getType()==3){
            int id1 = (*pit)->Seg1()->getId();
            string name1 = (*pit)->Seg1()->getName();
            int id2 = (*pit)->Seg2()->getId();
            string name2 = (*pit)->Seg2()->getName();   

            int id = ++jobCount;
            string mps_file1="MP_FILES/"+name1;
            string mps_file2="MP_FILES/"+name2;
            if (_singlet){
                mps_file1=mps_file1+"_n2s.mps";
                mps_file2=mps_file2+"_n2s.mps";
            }
            else{
                mps_file1=mps_file1+"_n2t.mps";
                mps_file2=mps_file2+"_n2t.mps";
            }                
            Property Input;
            Property *pInput = &Input.add("input","");
            Property *pSegment =  &pInput->add("segment" , boost::lexical_cast<string>(id1) );
            pSegment->setAttribute<string>("type", name1 );
            pSegment->setAttribute<int>("id", id1 );
            pSegment->setAttribute<string>("mps_file",mps_file1);
            pSegment =  &pInput->add("segment" , boost::lexical_cast<string>(id2) );
            pSegment->setAttribute<string>("type", name2 );
            pSegment->setAttribute<int>("id", id2 );
            pSegment->setAttribute<string>("mps_file",mps_file2);

            Job job(id, tag, Input, Job::AVAILABLE );
            job.ToStream(ofs,"xml");
        }
    }

    // CLOSE STREAM
    ofs << "</jobs>" << endl;    
    ofs.close();
    
    cout << endl << "... ... In total " << jobCount << " jobs" << flush;
    
}





}};
