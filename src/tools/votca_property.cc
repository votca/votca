/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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
#include <boost/format.hpp>
#include <boost/program_options.hpp>

// Local VOTCA includes
#include "votca/tools/application.h"
#include "votca/tools/globals.h"
#include "votca/tools/property.h"
#include "votca/tools/propertyiomanipulator.h"
#include "votca/tools/version.h"

using namespace std;
using namespace votca::tools;
namespace po = boost::program_options;

class VotcaProperty final : public Application {

 public:
  VotcaProperty() = default;
  ~VotcaProperty() = default;

  string ProgramName() { return "votca_property"; }

  void HelpText(ostream& out) { out << "Helper for parsing XML files"; }

  void Initialize() {

    format = "XML";
    level = 1;

    AddProgramOptions()("file", po::value<string>(), "xml file to parse")(
        "format", po::value<string>(),
        "output format [XML TXT]")("level", po::value<votca::Index>(),
                                   "output from this level ");
  };

  bool EvaluateOptions() {
    CheckRequired("file", "Missing XML file");
    return true;
  };

  void Run() {

    file = op_vm_["file"].as<string>();

    if (op_vm_.count("format")) {
      format = op_vm_["format"].as<string>();
    }
    if (op_vm_.count("level")) {
      level = op_vm_["level"].as<votca::Index>();
    }

    try {

      Property p;
      map<string, PropertyIOManipulator*> mformat_;
      mformat_["XML"] = &XML;
      mformat_["TXT"] = &TXT;
      mformat_["HLP"] = &HLP;
      p.LoadFromXML(file);

      if (mformat_.find(format) != mformat_.end()) {
        PropertyIOManipulator* piom = mformat_.find(format)->second;
        piom->setLevel(level);
        piom->setIndentation("");
        piom->setColorScheme<csRGB>();
        cout << *piom << p;
      } else {
        cout << "format " << format << " not supported \n";
      }

    } catch (std::exception& error) {
      cerr << "an error occurred:\n" << error.what() << endl;
    }
  };

 private:
  string file;
  string format;
  votca::Index level;
};

int main(int argc, char** argv) {
  VotcaProperty vp;
  return vp.Exec(argc, argv);
}
