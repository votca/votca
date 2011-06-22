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

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <votca/csg/csgapplication.h>
#include <votca/csg/trajectorywriter.h>

using namespace std;
using namespace votca::csg;

class CsgMapApp
    : public CsgApplication
{
public:
    string ProgramName() { return "csg_map"; }
    void HelpText(ostream &out) {
        out << "Map a reference trajectory to a coarse-grained trajectory.\n"
            "This program can be used to map a whole trajectory or to\n"
            "create an initial configuration for a coarse-grained run only.";
    }

    bool DoTrajectory() { return true;}
    bool DoMapping() { return true;}

    void Initialize() {
        CsgApplication::Initialize();
        AddProgramOptions()
            ("out", boost::program_options::value<string>(),
                "  output file for coarse-grained trajectory");
    }
    
    bool EvaluateOptions() {
        CsgApplication::EvaluateOptions();
        CheckRequired("trj", "no trajectory file specified");
        CheckRequired("out", "need to specify output trajectory");
        return true;
    }

    void BeginEvaluate(Topology *top, Topology *top_ref);

    void EvalConfiguration(Topology *top, Topology *top_ref) {
        _writer->Write(top);
    }
    void EndEvaluate() {
        _writer->Close();
        delete _writer;
    }

protected:
    TrajectoryWriter *_writer;
};

void CsgMapApp::BeginEvaluate(Topology *top, Topology *top_atom)
{
    string out = OptionsMap()["out"].as<string>();
    cout << "writing coarse-grained trajectory to " << out << endl;
    _writer = TrjWriterFactory().Create(out);
    if(_writer == NULL)
        throw runtime_error("output format not supported: " + out);
        
    _writer->Open(out);
};

int main(int argc, char **argv)
{
    CsgMapApp app;
    return app.Exec(argc, argv);
}

