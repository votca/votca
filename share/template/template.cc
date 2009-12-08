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
// 
// File:   template.cc
// Author: ruehle
//
// Created on June 8, 2008, 10:41 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <libversion.h>

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
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&no);


    // initialize the readers/writers,
    // this will be combined in an initialize function later
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    
    // lets read in some program options
    namespace po = boost::program_options;
        
    
    // Declare the supported options.
    po::options_description desc("Allowed options");    
    
    // let cg_engine add some program options
    cg_engine.AddProgramOptions(desc);
    
    // now read in the command line
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);    
        po::notify(vm);
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }
    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    
    // try to run the cg process, go through the frames, etc...
    try {
        cg_engine.Run(desc, vm);
    }
    // did an error occour?
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
    }
    return 0;
}

