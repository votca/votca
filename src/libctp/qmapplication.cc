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


#include <votca/ctp/qmapplication.h>
#include <votca/ctp/version.h>
#include <boost/format.hpp>

namespace votca { namespace ctp {

QMApplication::QMApplication() {
    Calculatorfactory::RegisterAll();
}


void QMApplication::Initialize(void) {
    votca::csg::TrajectoryWriter::RegisterPlugins();
    votca::csg::TrajectoryReader::RegisterPlugins();

    Calculatorfactory::RegisterAll();

    namespace propt = boost::program_options;

    AddProgramOptions() ("options,o", propt::value<string>(),
        "  calculator options");
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
    AddProgramOptions() ("restart,r", propt::value<string>()->default_value(""),
        "  parallelized job execution: job restart pattern; "
            "e.g. '-r host(mach1:1234,mach2:5678) stat(FAILED)'");
    AddProgramOptions() ("cache,c", propt::value<int>()->default_value(8),
        "  parallelized job execution: job cache size");
}


bool QMApplication::EvaluateOptions(void) {
    CheckRequired("options", "Please provide an xml file with calculator options");
    CheckRequired("file", "Please provide the state file");
    return true;
}


void QMApplication::Run() {

    load_property_from_xml(_options, _op_vm["options"].as<string>());

    // EVALUATE OPTIONS
    int nThreads = OptionsMap()["nthreads"].as<int>();
    int nframes = OptionsMap()["nframes"].as<int>();
    int fframe = OptionsMap()["first-frame"].as<int>();
    if (fframe-- == 0) throw runtime_error("ERROR: First frame is 0, counting "
                                           "in VOTCA::CTP starts from 1.");
    int  save = OptionsMap()["save"].as<int>();

    // STATESAVER & PROGRESS OBSERVER
    string statefile = OptionsMap()["file"].as<string>();
    StateSaverSQLite statsav;
    statsav.Open(_top, statefile);    

    ProgObserver< vector<Job*>, Job*, Job::JobResult > progObs
        = ProgObserver< vector<Job*>, Job*, Job::JobResult >();
    progObs.InitCmdLineOpts(OptionsMap());
    
    // INITIALIZE & RUN CALCULATORS
    cout << "Initializing calculators " << endl;
    BeginEvaluate(nThreads, &progObs);

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




void QMApplication::AddCalculator(QMCalculator* calculator) {
    _calculators.push_back(calculator);
}


void QMApplication::BeginEvaluate(int nThreads = 1, 
        ProgObserver< vector<Job*>, Job*, Job::JobResult > *obs = NULL) {
    list< QMCalculator* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        cout << "... " << (*it)->Identify() << " ";
        (*it)->setnThreads(nThreads);
        (*it)->setProgObserver(obs);
        (*it)->Initialize(&_top, &_options);        
        cout << endl;
    }
}

bool QMApplication::EvaluateFrame() {
    list< QMCalculator* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        cout << "... " << (*it)->Identify() << " " << flush;
        (*it)->EvaluateFrame(&_top);
        cout << endl;
    }
}

void QMApplication::EndEvaluate() {
    list< QMCalculator* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        (*it)->EndEvaluate(&_top);
    }
}

void QMApplication::ShowHelpText(std::ostream &out) {
    string name = ProgramName();
    if (VersionString() != "") name = name + ", version " + VersionString();
    votca::ctp::HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}

void QMApplication::PrintDescription(std::ostream &out, string name,  HelpOutputType _help_output_type) {
    
    // loading documentation from the xml file in VOTCASHARE
    char *votca_share = getenv("VOTCASHARE");
    if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/") + name + string(".xml");

    boost::format _format("%|3t|%1% %|20t|%2% \n");
    try {
        Property options;
        load_property_from_xml(options, xmlFile);

        switch (_help_output_type) {

            case _helpShort:
                _format % name % options.get("options." + name).getAttribute<string>("help");
                out << _format;
                break;
                
            case _helpLong:
                out << HLP << setlevel(2) << options;
        }

    } catch (std::exception &error) {
        out << _format % name % "Undocumented";
    }    
}

}}
