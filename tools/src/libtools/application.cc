/*
 * Copyright 2009-2024 The VOTCA Development Team (http://www.votca.org)
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

// Standard includes
#include <iostream>

// Third party includes
#include <boost/algorithm/string/replace.hpp>

// Local VOTCA includes
#include "votca/tools/application.h"
#include "votca/tools/globals.h"
#include "votca/tools/propertyiomanipulator.h"
#include "votca/tools/version.h"

// Local private VOTCA includes
#include "votca_tools_config.h"

namespace votca {
namespace tools {
using namespace std;
Application::Application() : op_desc_("Allowed options") {}

Application::~Application() = default;

void Application::ShowHelpText(std::ostream &out) {
  out << "==================================================\n";
  out << "========   VOTCA (http://www.votca.org)   ========\n";
  out << "==================================================\n\n";

  out << "please sread and cite: "PROJECT_CITATION "\n";
  out << "and submit bugs to " PROJECT_BUGREPORT "\n\n";
  out << ProgramName();
  if (VersionString() != "") {
    out << ", version " << VersionString();
  }
  out << endl << "votca_tools, version " << ToolsVersionStr() << "\n\n";

  HelpText(out);

  // remove Hidden group from the option list and print
  out << "\n\n" << VisibleOptions() << endl;
}

int Application::Exec(int argc, char **argv) {
  try {
    AddProgramOptions()("help,h", "  display this help and exit");
    AddProgramOptions()("verbose", "  be loud and noisy");
    AddProgramOptions()("verbose1", "  be very loud and noisy");
    AddProgramOptions()("verbose2,v", "  be extremly loud and noisy");

    Initialize();  // initialize program-specific parameters

    ParseCommandLine(argc,
                     argv);  // initialize general parameters & read input file

    Log::current_level = Log::error;
    if (op_vm_.count("verbose")) {
      Log::current_level = Log::warning;
    }
    if (op_vm_.count("verbose1")) {
      Log::current_level = Log::info;
    }

    if (op_vm_.count("verbose2")) {
      Log::current_level = Log::debug;
    }

    if (op_vm_.count("help")) {
      ShowHelpText(cout);
      return 0;
    }

    if (!EvaluateOptions()) {
      ShowHelpText(cout);
      return -1;
    }

    if (continue_execution_) {
      Run();
    } else {
      cout << "Done - stopping here\n";
    }
  } catch (std::exception &error) {
    cerr << "an error occurred:\n" << error.what() << endl;
    return -1;
  }
  return 0;
}

boost::program_options::options_description_easy_init
    Application::AddProgramOptions(const string &group) {
  // if no group is given, add it to standard options
  if (group == "") {
    return op_desc_.add_options();
  }

  // does group already exist, if yes, add it there
  std::map<string, boost::program_options::options_description>::iterator iter =
      op_groups_.find(group);
  if (iter != op_groups_.end()) {
    return iter->second.add_options();
  }

  // no group with given name was found -> create group
  op_groups_.insert(
      make_pair(group, boost::program_options::options_description(group)));

  return op_groups_[group].add_options();
}

void Application::ParseCommandLine(int argc, char **argv) {
  namespace po = boost::program_options;

  // default options should be added to visible (the rest is handled via a map))
  visible_options_.add(op_desc_);

  // add all categories to list of available options
  for (const auto &pair : op_groups_) {
    op_desc_.add(pair.second);
    if (pair.first != "Hidden") {
      visible_options_.add(pair.second);
    }
  }

  // parse the command line
  try {
    po::store(po::parse_command_line(argc, argv, op_desc_), op_vm_);
    po::notify(op_vm_);
  } catch (boost::program_options::error &err) {
    throw runtime_error(string("error parsing command line: ") + err.what());
  }
}

void Application::CheckRequired(const string &option_name,
                                const string &error_msg) {
  if (!op_vm_.count(option_name)) {
    ShowHelpText(cout);
    throw std::runtime_error("missing argument " + option_name + "\n" +
                             error_msg);
  }
}

}  // namespace tools
}  // namespace votca
