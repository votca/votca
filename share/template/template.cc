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

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <version.h>
#include <stdexcept>

using namespace std;

void help_text(void)
{
    votca::csg::HelpTextHeader("template");
    cout << "Template for VOTCA application\n\n";
}
using namespace std;

class CGAnalyzer
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {
    };
    void EndCG() {
    };
    
    void EvalConfiguration(Topology *top, Topology *top_atom = 0) {
    }
    
protected:
};


int main(int argc, char** argv)
{    
    // we have one observer
    CGAnalyzer no;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    try {
        cg_engine.Initialize();

        // add observer that it gets called to analyze frames
        cg_engine.AddObserver((CGObserver*)&no);
    
        // lets read in some program options
        namespace po = boost::program_options;
        
        // Add a user option
        cg_engine.AddProgramOptions()
            ("myoption", po::value<string>(), "  Example for a new option");
    
        cg_engine.ParseCommandLine(argc, argv);

        // some shortcuts
        po::variables_map &vm
            = cg_engine.OptionsMap();
        po::options_description &desc
            = cg_engine.OptionsDesc();

        // does the user want help?
        if (vm.count("help")) {
            cout << "csg_nemat, lib version " << LIB_VERSION_STR << "\n\n";                
            cout << desc << endl;
            return 0;
        }

        // or asks for the program version?
        if (vm.count("myoption"))
            cout << "myoption = " << vm["myoption"].as<string>() << endl;

        // start coarse graining
        cg_engine.Run();
    }
    // did an error occour?
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
    }
    return 0;
}

