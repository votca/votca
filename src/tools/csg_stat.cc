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

void help_text(void)
{
    votca::csg::HelpTextHeader("csg_stat");
    cout << "Calculate all distribuions (bonded + non-bonded) specified in options file.\n"
            "Optionally calculates update matrix for invere Monte Carlo. This program\n"
            "is called inside the inverse scripts. Unlike csg_boltzmann, big systems\n"
            "can be treated as well as non-bonded interactions evaluated.\n\n";
}

int main(int argc, char** argv)
{    
    int write_every=0;
    // we have one observer
    Imc imc;        
    // The CGEngine does the work
    CGEngine cg_engine;
    namespace po=boost::program_options;
    
    try {

        // add our observer that it gets called to analyze frames
        cg_engine.AddObserver((CGObserver*)&imc);
    
        // let cg_engine add some program options
        cg_engine.Initialize();
   
        cg_engine.AddProgramOptions()
            ("options", po::value<string>(), "  options file for coarse graining")
            ("do-imc", "  write out inverse monte carlo data")
            ("write-every", po::value<int>(&write_every), "  write afer every block of this length, " \
                "if --blocking is set, the averages are cleared after every write")
            ("do-blocks", "  write output for blocking analysis");
     
        cg_engine.ParseCommandLine(argc, argv);

        po::variables_map &vm
            = cg_engine.OptionsMap();
    
        // does the user want help?
        if (vm.count("help")) {
            help_text();
            cout << cg_engine.OptionsDesc() << endl;
            return 0;
        }
    
        if(!vm.count("options")) {
            cout << "need to specify options file\n";
            cout << cg_engine.OptionsDesc() << endl;
            return -1;
        }

        imc.WriteEvery(write_every);
        if(vm.count("do-blocks"))
            imc.DoBlocks(true);
        if(vm.count("do-imc"))
            imc.DoImc(true);
        
        imc.LoadOptions(vm["options"].as<string>());
            
        // try to run the cg process, go through the frames, etc...
        cg_engine.Run();
    }
    // did an error occour?
    catch(exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
        return -1;
    }
    return 0;
}


