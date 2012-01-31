#include <votca/ctp/qmapplication2.h>
#include <votca/ctp/version.h>


namespace votca { namespace ctp {

QMApplication2::QMApplication2() {
    CalculatorFactory2::RegisterAll();
}


void QMApplication2::Initialize(void) {
    votca::csg::TrajectoryWriter::RegisterPlugins();
    votca::csg::TrajectoryReader::RegisterPlugins();

    CalculatorFactory2::RegisterAll();

    namespace propt = boost::program_options;
    AddProgramOptions() ("cg", propt::value<string>(),
                         " segment & fragment definitions");
    AddProgramOptions() ("options,o", propt::value<string>(),
                         " program and calculator options");
    AddProgramOptions() ("file,f", propt::value<string>(),
                         " SQLite state file");
    AddProgramOptions() ("first-frame,i", propt::value<int>()->default_value(1),
                         " start from this frame");
    AddProgramOptions() ("nframes,n", propt::value<int>()->default_value(-1),
                         " number of frames to process");    
}


bool QMApplication2::EvaluateOptions(void) {
    CheckRequired("cg", "Please provide segment & fragment def. xml file");
    CheckRequired("options", "Please provide calculator options xml file");
    CheckRequired("file", "Please provide SQLite database file");
    return true;
}


void QMApplication2::Run() {

    load_property_from_xml(_options, _op_vm["options"].as<string>());

    int nframes = OptionsMap()["nframes"].as<int>();
    int fframe = OptionsMap()["first-frame"].as<int>();
    if (fframe-- == 0) throw runtime_error("ERROR: First frame is 0, counting "
                                           "in VOTCA::CTP starts from 1.");

    BeginEvaluate();



    string statefile = OptionsMap()["file"].as<string>();
    StateSaverSQLite2 statsav;
    statsav.Open(_top, statefile);
    //if (statsav.FramesInDatabase() != 1) {
    //    throw runtime_error("ERROR: Database contains more than one frame.");
    //}

    while (statsav.NextFrame()) {
        cout << "Evaluating frame " << _top.getDatabaseId() << endl;
        cout << "... ";
        EvaluateFrame();
        statsav.WriteFrame();
    }
    statsav.Close();
    EndEvaluate();

}




void QMApplication2::AddCalculator(QMCalculator2* calculator) {
    _calculators.push_back(calculator);
}


void QMApplication2::BeginEvaluate() {
    list< QMCalculator2* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        (*it)->Initialize(&_top, &_options);
    }
}

bool QMApplication2::EvaluateFrame() {
    list< QMCalculator2* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        cout << (*it)->Identify() << " ";
        (*it)->EvaluateFrame(&_top);
    }
}

void QMApplication2::EndEvaluate() {
    list< QMCalculator2* > ::iterator it;
    for (it = _calculators.begin(); it != _calculators.end(); it++) {
        (*it)->EndEvaluate(&_top);
    }
}

void QMApplication2::ShowHelpText(std::ostream &out) {
    string name = ProgramName();
    if (VersionString() != "") name = name + ", version " + VersionString();
    votca::ctp::HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}



}}
