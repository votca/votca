/* 
 *            Copyright 2009-2017 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/xtp/sqlapplication.h>
#include <votca/xtp/extractorfactory.h>
#include <votca/ctp/extractorfactory.h>


using namespace std;
using namespace votca;


class XtpDump : public xtp::SqlApplication
{
public:

    string  ProgramName() { return "xtp_dump"; }    

    void    HelpText(ostream &out) { out <<"Extracts information from the state file"<< endl; }
    void    HelpText() { };

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    
    //void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void XtpDump::Initialize() {
    xtp::ExtractorFactory::RegisterAll();
    ctp::ExtractorFactory::RegisterAll();
    xtp::SqlApplication::Initialize();
   

    AddProgramOptions("Extractors") ("extract,e", propt::value<string>(),
                      "List of extractors separated by ',' or ' '");
    AddProgramOptions("Extractors") ("list,l",
                      "Lists all available extractors");
    AddProgramOptions("Extractors") ("description,d", propt::value<string>(),
                      "Short description of an extractor");
    return;
}

bool XtpDump::EvaluateOptions() {

  if (OptionsMap().count("list")) {
    cout << "Available XTP extractors: \n";
    for (const auto& extract:xtp::Extractors().getObjects()) {
      PrintDescription(std::cout, extract.first, "xtp/xml", Application::HelpShort);
    }
    cout << "Available (wrapped) CTP calculators: \n";
    for (const auto& extract:ctp::Extractors().getObjects()) {
      bool printctp = true;
      std::string ctpcalc = extract.first.c_str();
      for (const auto& xtpextract:xtp::Extractors().getObjects()) {
        if (ctpcalc.compare(xtpextract.first.c_str()) == 0) {
          printctp = false;
          break;
        }
      }

      if (printctp){
        PrintDescription(std::cout, extract.first, "ctp/xml", Application::HelpShort);
      }
    }
    StopExecution();
  
    return true;
  }


  if (OptionsMap().count("description")) {
    CheckRequired("description", "no extractor is given");
    tools::Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
    // loop over the names in the description string
    for (const string& n:tok) {
      // loop over calculators
      bool printerror = true;
      for (const auto& extract:xtp::Extractors().getObjects()) {
        if (n.compare(extract.first.c_str()) == 0) {
          PrintDescription(std::cout,extract.first, "xtp/xml", Application::HelpLong);
          printerror = false;
          break;
        }
      }
      for (const auto& extract:ctp::Extractors().getObjects()) {
        if (n.compare(extract.first.c_str()) == 0) {
          bool printctp = true;
          std::string ctpcalc = extract.first.c_str();
          for (const auto& xtpextract:xtp::Extractors().getObjects()) {
            if (ctpcalc.compare(xtpextract.first.c_str()) == 0) {
              printctp = false;
              break;
            }
            if (printctp) {
              PrintDescription(std::cout, extract.first, "ctp/xml", Application::HelpLong);
              printerror = false;
              break;
            }
          }
        }
        }
        if (printerror) cout << "Extractor " << n << " does not exist\n";
      }
      StopExecution();
      return true;
    }

    xtp::SqlApplication::EvaluateOptions();
    CheckRequired("extract", "Nothing to do here: Abort.");

    tools::Tokenizer calcs(OptionsMap()["extract"].as<string>(), " ,\n\t");
     for (const string& n:calcs) {
      
      bool _found_calc = false;
      for (const auto& extract:xtp::Extractors().getObjects()) {
        if (n.compare(extract.first.c_str()) == 0) {
          cout << " This is a XTP app" << endl;
          xtp::SqlApplication::AddCalculator(xtp::Extractors().Create(n.c_str()));
          _found_calc = true;
        }
      }
      if (!_found_calc) {
        for (const auto& extract:ctp::Extractors().getObjects()) {
          if (n.compare(extract.first.c_str()) == 0) {
            _found_calc = true;
            cout << " This is a CTP app" << endl;
            xtp::SqlApplication::AddCalculator(ctp::Extractors().Create(n.c_str()));
          }
      }
    }
    if (!_found_calc) {
      cout << "Extractor " << n << " does not exist\n";
      StopExecution();
    }

  }
  return true;
  }

int main(int argc, char** argv) {
    
    XtpDump xtpdump;
    return xtpdump.Exec(argc, argv);

}
