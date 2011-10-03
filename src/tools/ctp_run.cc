/*
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <stdlib.h>
#include <votca/ctp/qmapplication.h>
#include <votca/ctp/calculatorfactory.h>
#include <string>

using namespace votca::ctp;

/*
 *
 */
class QMAppRun : public QMApplication
{
public:
    void HelpText() {}

    string ProgramName() { return "ctp_run"; }

    void HelpText(std::ostream &out) {
        out << "Runs specified calculators." << endl;
    }

    void Initialize() {
        QMApplication::Initialize();
        AddProgramOptions("Calculators")
            ("exec", boost::program_options::value<string>(), "list of calculators separated by commas or spaces")
	    ("list", "lists all available calculators")
            ("desc", boost::program_options::value<string>(), "detailed description of a calculator");
    }

    //TODO: Support for XML-File based options
    bool EvaluateOptions() {
        if(OptionsMap().count("list")) {
            cout << "Available calculators: \n";
            for(CalculatorFactory::assoc_map::const_iterator iter=Calculators().getObjects().begin();
                    iter != Calculators().getObjects().end(); ++iter) {
                cout << Description( (iter->first).c_str(), _short );
            }
            StopExecution();
            return true;
        }


         if(OptionsMap().count("desc")) {
            CheckRequired("desc", "no calculator is given");
 	    Tokenizer tok(OptionsMap()["desc"].as<string>(), " ,\n\t");
            // loop over the names in the description string
            for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n) {
                // loop over calculators
                bool printerror = true;
                for(CalculatorFactory::assoc_map::const_iterator iter=Calculators().getObjects().begin(); 
                        iter != Calculators().getObjects().end(); ++iter) {

                    if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                         cout << Description( (iter->first).c_str(), _long );
                        printerror = false;
                        //break;
                    }
                 }
                 if ( printerror ) cout << "Calculator " << *n << " does not exist\n";
            }
            StopExecution();
            return true;
         }

        QMApplication::EvaluateOptions();
        CheckRequired("exec", "no calculator is given");
        
        Tokenizer tok(OptionsMap()["exec"].as<string>(), " ,\n\t");
        for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n)
            AddCalculator(Calculators().Create((*n).c_str()));
        return true;
    }

    
    const char *Description(const char *name, const bool length) {
        // loading the documentation xml file from VOTCASHARE
        char *votca_share = getenv("VOTCASHARE");
        if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
        string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/")+name+string(".xml");
        try {
            Property options;
            load_property_from_xml(options, xmlFile);

           string description = string(" ");
           if ( length ) { // short description of the calculator
               
                 description += string("  ") + _fwstring(string(name),18);
                 description += options.get(name+string(".description")).as<string>();

            } else { // long description of the calculator
                description += _fwstring(string(name),18);
                description += options.get(name+string(".description")).as<string>()+string("\n");
 
                list<Property *> items = options.Select(name+string(".item"));

                for(list<Property*>::iterator iter = items.begin(); iter!=items.end(); ++iter) {
                    //cout << "Long description" << endl;
                    Property *pname=&( (*iter)->get( string("name") ) );
                    Property *pdesc=&( (*iter)->get( string("description") ) );
                    //Property *pdflt=&( (*iter)->get( string("default") ) );

                    description += string("  -") + _fwstring(pname->value(), 12) +
                                   string(" [") +
                                   //(*iter)->path() + string(".") +
                                   //(*iter)->name() + string(".") +
                                    string("] ");
                    description += pdesc->value() +string("\n");
                 }
            }
            description += string("\n");
            return description.c_str();

        } catch(std::exception &error) {
            return (string("XML file or description tag missing: ")+xmlFile).c_str();
        }
    }

private:
    static const bool _short = true;
    static const bool _long = false;

    string _fwstring(string original, size_t charCount ) {
        original.resize( charCount, ' ' );
        return original;
    }

    
};

int main(int argc, char** argv) {
    QMAppRun qmapprun;
    return qmapprun.Exec(argc, argv);
}
