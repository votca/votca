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
/* 
 * File:   csg_fmatch_main.cc
 * Author: lukyanov
 *
 * Created on June 10, 2009, 5:00 PM
 */

#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include "csg_fmatch.h"
#include "version.h"

void help_text(void)
{
    votca::csg::HelpTextHeader("csg_fmatch");
    cout << "Perform force matching (also call multiscale coarse-graining)\n\n";
}

/*
 * 
 */
int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGForceMatching fmatch;        
    // The CGEngine does the work
    CGEngine cg_engine;

    namespace po = boost::program_options;

    try {
        cg_engine.Initialize();
        // add our observer that it gets called to analyze frames
        cg_engine.AddObserver((CGObserver*)&fmatch);
    
        cg_engine.AddProgramOptions()
            ("options", po::value<string>(), "  options file for coarse graining");
    
        cg_engine.ParseCommandLine(argc, argv);

        po::variables_map &vm
            = cg_engine.OptionsMap();
        po::options_description &desc
            = cg_engine.OptionsDesc();

        // does the user want help?
        if (vm.count("help")) {
            help_text();
            cout << desc << endl;
            return 0;
        }
    
        if(!vm.count("options")) {
            cout << "need to specify options file\n";
            cout << desc << endl;
            return -1;
        }
    
        fmatch.LoadOptions(vm["options"].as<string>());
    
        cg_engine.Run();
    }
    // did an error occour?
    catch(exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
        return -1;
    }
    return 0;
}
