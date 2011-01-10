/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/program_options.hpp>
#include <cgengine.h>
#include "version.h"
#include "imc.h"

using namespace std;
using namespace votca::csg;
using namespace votca::tools;


#include <stdlib.h>
#include <csgapplication.h>

//using namespace votca::tools;
using namespace std;
using namespace votca::csg;

class CsgStatApp
    : public CsgApplication
{
public:
    CsgStatApp() : _write_every(0) {}
    
    string ProgramName() { return "csg_stat"; }
    void HelpText(ostream &out);

    bool DoTrajectory() {return true;}
    bool DoMapping() {return true;}
    bool DoMappingDefault(void) { return false; }
    bool DoThreaded() {return true; }
    bool SynchronizeThreads() {return true;}
    void Initialize();
    bool EvaluateOptions();

    void BeginEvaluate(Topology *top, Topology *top_ref);
    void EndEvaluate();

    CsgApplication::Worker *ForkWorker() {
        return _imc.ForkWorker();
    }

    void MergeWorker(CsgApplication::Worker *worker) {
        _imc.MergeWorker(worker);
    }

public:
    Imc _imc;
    int _write_every;
};

void CsgStatApp::HelpText(ostream &out)
{
    out << "Calculate all distribuions (bonded + non-bonded) specified in options file.\n"
            "Optionally calculates update matrix for invere Monte Carlo. This program\n"
            "is called inside the inverse scripts. Unlike csg_boltzmann, big systems\n"
            "can be treated as well as non-bonded interactions evaluated.";
}

void CsgStatApp::Initialize()
{
    CsgApplication::Initialize();
    AddProgramOptions("Specific options")
            ("options", boost::program_options::value<string>(), "  options file for coarse graining")
            ("do-imc", "  write out inverse monte carlo data")
            ("write-every", boost::program_options::value<int>(&_write_every), "  write afer every block of this length, " \
                "if --blocking   is set, the averages are cleared after every write")
            ("do-blocks", "  write output for blocking analysis");
}

bool CsgStatApp::EvaluateOptions()
{
    CsgApplication::EvaluateOptions();
    CheckRequired("options");
    CheckRequired("trj", "no trajectory file specified");
    
    _imc.LoadOptions(OptionsMap()["options"].as<string>());

    _imc.WriteEvery(_write_every);
    if(OptionsMap().count("do-blocks"))
        _imc.DoBlocks(true);
    if(OptionsMap().count("do-imc"))
    _imc.DoImc(true);

    _imc.Initialize();
    return true;
}

void CsgStatApp::BeginEvaluate(Topology *top, Topology *top_ref)
{
    _imc.BeginEvaluate(top, top_ref);
}

void CsgStatApp::EndEvaluate()
{
    _imc.EndEvaluate();
}

int main(int argc, char** argv)
{
    CsgStatApp app;
    app.Exec(argc, argv);
}

