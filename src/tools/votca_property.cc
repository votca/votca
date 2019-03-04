/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#include <boost/algorithm/string/replace.hpp>
#include <boost/format.hpp>
#include <boost/program_options.hpp>
#include <iostream>

#include <list>
#include <votca/tools/application.h>
#include <votca/tools/globals.h>
#include <votca/tools/property.h>
#include <votca/tools/propertyiomanipulator.h>
#include <votca/tools/version.h>

using namespace std;
using namespace votca::tools;
namespace po = boost::program_options;

class VotcaProperty : public Application {

 public:
  VotcaProperty();
  ~VotcaProperty();

  string ProgramName() { return "votca_property"; }

  void HelpText(ostream& out) { out << "Helper for parsing XML files"; }

  void Initialize() {

    format = "XML";
    level = 1;

    AddProgramOptions()("file", po::value<string>(), "xml file to parse")(
        "format", po::value<string>(), "output format [XML TXT TEX]")(
        "level", po::value<int>(), "output from this level ");
  };

  bool EvaluateOptions() {
    CheckRequired("file", "Missing XML file");
    return true;
  };

  void Run() {

    file = _op_vm["file"].as<string>();

    if (_op_vm.count("format")) format = _op_vm["format"].as<string>();
    if (_op_vm.count("level")) level = _op_vm["level"].as<int>();

    try {

      Property p;

      map<string, PropertyIOManipulator*> _mformat;
      map<string, PropertyIOManipulator*>::iterator it;

      _mformat["XML"] = &XML;
      _mformat["TXT"] = &TXT;
      _mformat["TEX"] = &TEX;
      _mformat["HLP"] = &HLP;

      load_property_from_xml(p, file);

      it = _mformat.find(format);
      if (it != _mformat.end()) {
        PropertyIOManipulator* piom = _mformat.find(format)->second;
        piom->setLevel(level);
        piom->setIndentation("");
        piom->setColorScheme<csRGB>();
        cout << *piom << p;
      } else {
        cout << "format " << format << " not supported \n";
      }

      // PropertyIOManipulator XML(PropertyIOManipulator::XML,0,"---");
      // cout << XML << p;
      // cout << TXT << p;
      // cout << T2T << p;
      // cout << LOG << p;
      // cout << TEX << p;

    } catch (std::exception& error) {
      cerr << "an error occurred:\n" << error.what() << endl;
    }
  };

 private:
  string file;
  string format;
  int level;
};

VotcaProperty::VotcaProperty(void) {}

VotcaProperty::~VotcaProperty(void) {}

int main(int argc, char** argv) {
  VotcaProperty vp;
  return vp.Exec(argc, argv);
}
