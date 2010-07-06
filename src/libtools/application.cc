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
#include "application.h"
#include "version.h"

namespace votca { namespace tools {

Application::Application()
{
}

Application::~Application()
{
}

void Application::HelpText(std::ostream &out)
{
    out << "no help available";
}

string Application::VersionString()
{
    return "";
}

void Application::ShowHelpText(std::ostream &out)
{
    out << "\t------ VOTCA ( http://www.votca.org ) ------\n"
        << ProgramName();
    if(VersionString() != "")
        out << ", version " << VersionString();
    out << endl
        << "votca_tools, version " << ToolsVersionStr()
        << "\n\n";
}

int Application::Exec(int argc, char **argv)
{
    try {
        AddProgramOptions()("help,h", "  produce this help message");
        Initialize(); /// initialize program-specific parameters

        ParseCommandLine(argc, argv); /// initialize general parameters & read input file

        if (_op_vm.count("help")) {
            ShowHelpText(cout);
            return 0;
        }

        EvaluateOptions();
        Run();
    }
    catch(std::exception &error) {
         cerr << "an error occured:\n" << error.what() << endl;
         return -1;
    }
    return 0;
}

void Application::ParseCommandLine(int argc, char **argv)
{
    namespace po = boost::program_options;

    /// parse the command line
    try {
        po::store(po::parse_command_line(argc, argv, _op_desc), _op_vm);
        po::notify(_op_vm);
    }
    catch(boost::program_options::error err) {
        throw runtime_error(string("error parsing command line: ") + err.what());
    }
}

}}

