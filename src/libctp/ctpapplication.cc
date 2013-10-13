/*
 *            Copyright 2009-2012 The VOTCA Development Team
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


#include <votca/ctp/ctpapplication.h>
#include <votca/ctp/version.h>
#include <votca/tools/propertyiomanipulator.h>
#include <boost/format.hpp>

namespace votca { namespace ctp {

CtpApplication::CtpApplication() {
    ;
}

/**
 * \brief Adds program options to the executable
 * 
 * Every executable requires option file for calculators it is running
 * It is thus a part of the base CtpApplication class 
 * 
 */
void CtpApplication::Initialize(void) {

     AddProgramOptions() ("options,o", boost::program_options::value<string>(),
        "  calculator options");
}


bool CtpApplication::EvaluateOptions(void) {
    return true;
}


void CtpApplication::ShowHelpText(std::ostream &out) {
    string name = ProgramName();
    if (VersionString() != "") name = name + ", version " + VersionString();
    votca::ctp::HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}


void CtpApplication::PrintDescription(std::ostream &out, string name,  HelpOutputType _help_output_type) {
    
    // loading documentation from the xml file in VOTCASHARE
    char *votca_share = getenv("VOTCASHARE");
    if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/") + name + string(".xml");

    boost::format _format("%|3t|%1% %|20t|%2% \n");
    try {
        votca::tools::Property options;
        load_property_from_xml(options, xmlFile);

        switch (_help_output_type) {

            case _helpShort:
                _format % name % options.get("options." + name).getAttribute<string>("help");
                out << _format;
                break;
                
            case _helpLong:
                votca::tools::PropertyIOManipulator iom(votca::tools::PropertyIOManipulator::HLP, 2, "");
                out << iom << options;
        }

    } catch (std::exception &error) {
        out << _format % name % "Undocumented";
    }    
}


}}
