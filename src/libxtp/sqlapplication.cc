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


#include <votca/xtp/sqlapplication.h>
#include <votca/xtp/calculatorfactory.h>
#include <votca/xtp/version.h>
#include <boost/format.hpp>

#include <votca/ctp/calculatorfactory.h>


namespace votca { namespace xtp {

    namespace CTP = votca::ctp;
    
XSqlApplication::XSqlApplication() {
    XCalculatorfactory::RegisterAll();
    CTP::Calculatorfactory::RegisterAll();
}


void XSqlApplication::Initialize(void) {
    XtpApplication::Initialize();

    XCalculatorfactory::RegisterAll();
    CTP::Calculatorfactory::RegisterAll();

    namespace propt = boost::program_options;

    AddProgramOptions() ("file,f", propt::value<string>(),
        "  sqlight state file, *.sql");
    AddProgramOptions() ("first-frame,i", propt::value<int>()->default_value(1),
        "  start from this frame");
    AddProgramOptions() ("nframes,n", propt::value<int>()->default_value(1),
        "  number of frames to process");
    AddProgramOptions() ("nthreads,t", propt::value<int>()->default_value(1),
        "  number of threads to create");
    AddProgramOptions() ("save,s", propt::value<int>()->default_value(1),
        "  whether or not to save changes to state file");
}


bool XSqlApplication::EvaluateOptions(void) {
    CheckRequired("file", "Please provide the state file");
    return true;
}


void XSqlApplication::Run() {

    // load_property_from_xml(_options, _op_vm["options"].as<string>());

    // EVALUATE OPTIONS
    int nThreads = OptionsMap()["nthreads"].as<int>();
    int nframes = OptionsMap()["nframes"].as<int>();
    int fframe = OptionsMap()["first-frame"].as<int>();
    if (fframe-- == 0) throw runtime_error("ERROR: First frame is 0, counting "
                                           "in VOTCA::XTP starts from 1.");
    int  save = OptionsMap()["save"].as<int>();

    // STATESAVER & PROGRESS OBSERVER
    string statefile = OptionsMap()["file"].as<string>();
    StateSaverSQLite statsav;
    statsav.Open(_top, statefile);
    
    // INITIALIZE & RUN CALCULATORS
    cout << "Initializing calculators " << endl;
    BeginEvaluate(nThreads);

    int frameId = -1;
    int framesDone = 0;
    while (statsav.NextFrame() && framesDone < nframes) {
        frameId += 1;
        if (frameId < fframe) continue;
        cout << "Evaluating frame " << _top.getDatabaseId() << endl;
        EvaluateFrame();
        if (save == 1) { statsav.WriteFrame(); }
        else { cout << "Changes have not been written to state file." << endl; }
        framesDone += 1;
    }
    
    if (framesDone == 0)
        cout << "Input requires first frame = " << fframe+1 << ", # frames = " 
             << nframes << " => No frames processed.";
    
    statsav.Close();
    EndEvaluate();

}




void XSqlApplication::AddCalculator(XQMCalculator* calculator) {
    _calculators.push_back(calculator);
}

void XSqlApplication::AddCalculator(CTP::QMCalculator* calculator) {
    _ctp_calculators.push_back(calculator);
}

void XSqlApplication::BeginEvaluate(int nThreads = 1) {
    list< XQMCalculator* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        cout << "... " << (*it)->Identify() << " ";
        (*it)->setnThreads(nThreads);
        (*it)->Initialize(&_options); 
        cout << endl;
    }
    

    list< CTP::QMCalculator* > ::iterator cit;
    for (cit = _ctp_calculators.begin(); cit != _ctp_calculators.end(); cit++) {
        cout << "... " << (*cit)->Identify() << " ";
        (*cit)->setnThreads(nThreads);
        (*cit)->Initialize(&_options); 
        cout << endl;
    }
    
}

bool XSqlApplication::EvaluateFrame() {
    list< XQMCalculator* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        cout << "... " << (*it)->Identify() << " " << flush;
        (*it)->EvaluateFrame(&_top);
        cout << endl;
    }
    list< CTP::QMCalculator* > ::iterator cit;
    for (cit = _ctp_calculators.begin(); cit != _ctp_calculators.end(); cit++) {
        cout << "... " << (*cit)->Identify() << " " << flush;
        (*cit)->EvaluateFrame(&_top);
        cout << endl;
    }   
    
    
    return true;
}

void XSqlApplication::EndEvaluate() {
    list< XQMCalculator* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        (*it)->EndEvaluate(&_top);
    }
    
    list< CTP::QMCalculator* > ::iterator cit;
    for (cit = _ctp_calculators.begin(); cit != _ctp_calculators.end(); cit++) {
        (*cit)->EndEvaluate(&_top);
    }
    
    
}

}}
