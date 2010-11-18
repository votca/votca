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
    : _op_desc("Allowed options")
{
}

Application::~Application()
{
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

    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}

int Application::Exec(int argc, char **argv)
{
    try {
        AddProgramOptions()("help,h", "  produce this help message");
        Initialize(); // initialize program-specific parameters

        ParseCommandLine(argc, argv); // initialize general parameters & read input file

        if (_op_vm.count("help")) {
            ShowHelpText(cout);
            return 0;
        }

        if(!EvaluateOptions()) {
            ShowHelpText(cout);
            return -1;
        }
        
        Run();
    }
    catch(std::exception &error) {
         cerr << "an error occurred:\n" << error.what() << endl;
         return -1;
    }
    return 0;
}

boost::program_options::options_description_easy_init
    Application::AddProgramOptions(const string &group)
{
    // if no group is given, add it to standard options
    if(group == "")
        return _op_desc.add_options();
    
    // does group already exist, if yes, add it there
    std::map<string, boost::program_options::options_description>::iterator iter;
    iter = _op_groups.find(group);
    if(iter!=_op_groups.end())
        return iter->second.add_options();

    // no group with given name was found -> create group
    _op_groups.insert(make_pair(group, boost::program_options::options_description(group)));

    return _op_groups[group].add_options();
}


void Application::ParseCommandLine(int argc, char **argv)
{
    namespace po = boost::program_options;

    std::map<string, boost::program_options::options_description>::iterator iter;

    // add all cathegories to list of available options
    for(iter=_op_groups.begin(); iter!=_op_groups.end(); ++iter)
        _op_desc.add(iter->second);
    
    // parse the command line
    try {
        po::store(po::parse_command_line(argc, argv, _op_desc), _op_vm);
        po::notify(_op_vm);
    }
    catch(boost::program_options::error err) {
        throw runtime_error(string("error parsing command line: ") + err.what());
    }
}

void Application::CheckRequired(const string &option_name, const string &error_msg)
{
    if(!_op_vm.count(option_name)) {
        ShowHelpText(cout);
        throw std::runtime_error("missing argument " + option_name + "\n" + error_msg);
    }
}


}}

