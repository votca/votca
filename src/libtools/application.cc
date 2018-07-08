/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#include <votca/tools/application.h>
#include <votca/tools/version.h>
#include <votca/tools/globals.h>
#include <votca/tools/propertyiomanipulator.h>

#include <boost/format.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/filesystem.hpp>

namespace votca { namespace tools {
    using namespace std;
Application::Application()
    : _op_desc("Allowed options"), _continue_execution(true)
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
    out << "==================================================\n";
    out << "========   VOTCA (http://www.votca.org)   ========\n";
    out << "==================================================\n\n";

    out << "please submit bugs to bugs@votca.org\n\n";
    out << ProgramName();
    if(VersionString() != "")
        out << ", version " << VersionString();
    out << endl
        << "votca_tools, version " << ToolsVersionStr()
        << "\n\n";

    HelpText(out);
    
    // remove Hidden group from the option list and print
    out << "\n\n" << VisibleOptions() << endl;
    
    //out << "\n\n" << OptionsDesc() << endl;
}
     
void Application::ShowManPage(std::ostream &out) {
    
        out << boost::format(globals::man::header) %  ProgramName() % VersionString();        
        out << boost::format(globals::man::name) % ProgramName() % globals::url;
        out << boost::format(globals::man::synopsis) % ProgramName();
        std::stringstream ss;
        HelpText(ss);
        out << boost::format(globals::man::description) % ss.str();
        out << boost::format(globals::man::options);

        typedef std::vector<boost::shared_ptr<boost::program_options::option_description> >::const_iterator OptionsIterator;
        OptionsIterator it = _op_desc.options().begin(), it_end = _op_desc.options().end();
        
        while(it < it_end) {
            string format_name = (*it)->format_name() + " " + (*it)->format_parameter();
            boost::replace_all(format_name, "-", "\\-");
            out << boost::format(globals::man::option) % format_name % (*it)->description();
            ++it;           
        }      

        out << boost::format(globals::man::authors) % globals::email;
        out << boost::format(globals::man::copyright) % globals::url;
    
}

void Application::ShowTEXPage(std::ostream &out) {
    string program_name = ProgramName();
    boost::replace_all(program_name, "_", "\\_");
        out << boost::format(globals::tex::section) %  program_name;        
        out << boost::format(globals::tex::label) % ProgramName();
        std::stringstream ss, os;
        HelpText(ss);
        out << boost::format(globals::tex::description) % ss.str();

       typedef std::vector<boost::shared_ptr<boost::program_options::option_description> >::const_iterator OptionsIterator;
       OptionsIterator it = _op_desc.options().begin(), it_end = _op_desc.options().end();
       while(it < it_end) {
            string format_name = (*it)->format_name() + " " + (*it)->format_parameter();
            boost::replace_all(format_name, "-", "{-}");
            os << boost::format(globals::tex::option) % format_name % (*it)->description();
            ++it;           
        }      
        out << boost::format(globals::tex::options) % os.str();
}

int Application::Exec(int argc, char **argv)
{
    try {
        //_continue_execution = true;
	AddProgramOptions()("help,h", "  display this help and exit");
	AddProgramOptions()("verbose,v", "  be loud and noisy");
	AddProgramOptions("Hidden")("man", "  output man-formatted manual pages");
	AddProgramOptions("Hidden")("tex", "  output tex-formatted manual pages");
	
	Initialize(); // initialize program-specific parameters

        ParseCommandLine(argc, argv); // initialize general parameters & read input file

        if (_op_vm.count("verbose")) {
	  globals::verbose = true;
        }
        
        if (_op_vm.count("man")) {
            ShowManPage(cout);
            return 0;
        }
        
        if (_op_vm.count("tex")) {
            ShowTEXPage(cout);
            return 0;
        }
        
        if (_op_vm.count("help")) {
            ShowHelpText(cout);
            return 0;
        }

        if(!EvaluateOptions()) {
            ShowHelpText(cout);
            return -1;
        }

        if(_continue_execution)
            Run();
	else cout << "nothing to be done - stopping here\n";
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
    
    // default options should be added to visible (the rest is handled via a map))
    _visible_options.add(_op_desc);
    
    // add all categories to list of available options
    for(iter=_op_groups.begin(); iter!=_op_groups.end(); ++iter) {       
        _op_desc.add(iter->second);
        if ( iter->first != "Hidden" ) _visible_options.add(iter->second);
    }
    
    // parse the command line
    try {
        po::store(po::parse_command_line(argc, argv, _op_desc), _op_vm);
        po::notify(_op_vm);
    }
    catch(boost::program_options::error& err) {
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

void Application::PrintDescription(std::ostream &out, const string &calculator_name, const string help_path, HelpType help_type )
{
    boost::format _format("%|3t|%1% %|20t|%2% \n");
    string help_string;
    boost::filesystem::path arg_path;  
    Property options;
    string xmlFile = (arg_path / string(getenv("VOTCASHARE")) / help_path / (boost::format("%1%.%2%") % calculator_name % "xml").str()).string().c_str();
    
    // loading the documentation xml file from VOTCASHARE
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");

    try {

                load_property_from_xml(options, xmlFile);
                Property &calculator_options = options.get("options." + calculator_name);
                Property::AttributeIterator atr_it = calculator_options.findAttribute("help");

                if (atr_it != calculator_options.lastAttribute()) {
                    help_string = (*atr_it).second;
                } else {
                    if (tools::globals::verbose) out << _format % calculator_name % "Undocumented";
                    return;
                }

                switch (help_type) {
                    default:
                        break;
                    case HelpShort: // short description of the calculator
                        out << _format % calculator_name % help_string;
                        break;
                    case HelpLong:
                        votca::tools::PropertyIOManipulator iom(votca::tools::PropertyIOManipulator::HLP, 2, "");
                        out << iom << options;
                        break;
                }

            } catch (std::exception &error) {
                if (tools::globals::verbose) out << _format % calculator_name % "Undocumented";
            }

}

}}
