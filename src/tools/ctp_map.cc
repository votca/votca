/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <fstream>
#include <stdexcept>
#include <votca/ctp/qmtopology.h>
#include "md2qm_observer.h"
#include <votca/csg/csgapplication.h>

using namespace std;
using namespace votca::csg;
using namespace votca::ctp;

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
            ("segments,s", po::value<string>(), "  conjugated segment definitions")
            ("file,f", po::value<string>(), " sqlite state file");
}

bool CtpMapApp::EvaluateOptions()
{
    CsgApplication::EvaluateOptions();
    CheckRequired("segments");

    _observer.setOut(OptionsMap()["file"].as<string>());
    _observer.Initialize(_qmtopol, _options);

    // add our observer that it gets called to analyze frames
    AddObserver(dynamic_cast<CGObserver*>(&_observer));
    return true;
}

void CtpMapApp::BeginEvaluate(Topology *top, Topology *top_ref)
{
    _qmtopol.LoadListCharges(OptionsMap()["segments"].as<string>());
    CsgApplication::BeginEvaluate(top, top_ref);
}

int main(int argc, char** argv)
{
    CtpMapApp app;
    return app.Exec(argc, argv);
}
