/* 
 * File:   ctp_map.cc
 */

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <votca/ctp/qmtopology.h>
#include "md2qm_observer.h"
#include <votca/csg/csgapplication.h>

using namespace std;
using namespace votca::csg;

class CtpMapApp
    : public CsgApplication
{
    string ProgramName() { return "ctp_map"; }
    void HelpText(ostream &out) {
        out << "Converts atomistic topology and trajectory into coarse-grained ones";
    }

    void Initialize();
    bool EvaluateOptions();
    void BeginEvaluate(Topology *top, Topology *top_ref);

    bool DoTrajectory() {return true;}
    bool DoMapping() {return true;}


protected:
    // we have one observer
    MD2QMObserver _observer;
    Property _options;
    QMTopology _qmtopol;
};

namespace po = boost::program_options;

void CtpMapApp::Initialize()
{
    CsgApplication::Initialize();
    AddProgramOptions("Mapping options")
            ("listcharges,l", po::value<string>(), "  conjugated segment definitions")
            ("db", po::value<string>()->default_value("state.db"), " state file");
}

bool CtpMapApp::EvaluateOptions()
{
    CsgApplication::EvaluateOptions();
    CheckRequired("listcharges");

    _observer.setOut(OptionsMap()["out"].as<string>());
    _observer.Initialize(_qmtopol, _options);

    // add our observer that it gets called to analyze frames
    AddObserver(dynamic_cast<CGObserver*>(&_observer));
    return true;
}

void CtpMapApp::BeginEvaluate(Topology *top, Topology *top_ref)
{
    _qmtopol.LoadListCharges(OptionsMap()["listcharges"].as<string>());
    CsgApplication::BeginEvaluate(top, top_ref);
}

int main(int argc, char** argv)
{
    CtpMapApp app;
    return app.Exec(argc, argv);
}
