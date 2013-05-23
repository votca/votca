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
#include <votca/kmc/kmccalculatorfactory.h>
#include <votca/kmc/kmcapplication.h>
#include <string>
#include <cstdio>
#include <iostream>

namespace votca { namespace kmc {

KMCApplication::KMCApplication()
{
    KMCCalculatorFactory::RegisterAll();
}

KMCApplication::~KMCApplication()
{}

void KMCApplication::Initialize() 
{
    Application::Initialize();
    KMCCalculatorFactory::RegisterAll();

    AddProgramOptions()
            ("options,o", boost::program_options::value<string>(), "  program and calculator options")
            ("file,f", boost::program_options::value<string>(), "  sqlite state file")
            ("textfile,t", boost::program_options::value<string>(), "  output text file (otherwise: screen output)");
    
    AddProgramOptions("Calculators")
            ("execute,e", boost::program_options::value<string>(), "list of calculators separated by commas or spaces")
            ("list,l", "lists all available calculators")
            ("description,d", boost::program_options::value<string>(), "detailed description of a calculator");
}

void KMCApplication::ShowHelpText(std::ostream &out)
{
    string name =  ProgramName();
    if(VersionString() != "")
         name = name + ", version " + VersionString();

    votca::kmc::HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}


void KMCApplication::PrintDescription(const char *name, const bool length)
{
        // loading the documentation xml file from VOTCASHARE
        char *votca_share = getenv("VOTCASHARE");
        if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
        string xmlFile = string(getenv("VOTCASHARE")) + string("/kmc/xml/")+name+string(".xml");
        try {
            Property options;
            load_property_from_xml(options, xmlFile);

           if ( length ) { // short description of the calculator
               
                 cout << string("  ") << _fwstring(string(name),14);
                 cout << options.begin()->getAttribute<string>("help");

            } else { // long description of the calculator
                cout << HLP << options;
             }
            cout << endl;
        } catch(std::exception &error) {
            cout << string("XML file or description tag missing: ") << xmlFile << endl;
        }
}

// check if required options are provided
bool KMCApplication::EvaluateOptions() {

           if(OptionsMap().count("list")) {
            cout << "Available calculators: \n";
            for(KMCCalculatorFactory::assoc_map::const_iterator iter=Calculators().getObjects().begin();
                    iter != Calculators().getObjects().end(); ++iter) {
                PrintDescription( (iter->first).c_str(), _short );
            }
            StopExecution();
            return true;
        }


         if(OptionsMap().count("description")) {
            CheckRequired("description", "no calculator is given");
 	    Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
            // loop over the names in the description string
            for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n) {
                // loop over calculators
                bool printerror = true;
                for(KMCCalculatorFactory::assoc_map::const_iterator iter=Calculators().getObjects().begin(); 
                        iter != Calculators().getObjects().end(); ++iter) {

                    if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                         PrintDescription( (iter->first).c_str(), _long );
                        printerror = false;
                        break;
                    }
                 }
                 if ( printerror ) cout << "Calculator " << *n << " does not exist\n";
            }
            StopExecution();
            return true;
         }

        Application::EvaluateOptions();

        CheckRequired("execute", "no calculator is given");
        Tokenizer tok(OptionsMap()["execute"].as<string>(), " ,\n\t");
        for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n)
            AddCalculator(Calculators().Create((*n).c_str()));    
        
        CheckRequired("options", "please provide an xml file with program options");
        CheckRequired("file", "no database file specified");
        
        string _outputfile = "";
        //if(OptionsMap()["textfile"] == NULL)
        //{cout << "Mist!";}
        if(OptionsMap().count("textfile")) 
        {
            _outputfile = OptionsMap()["textfile"].as<string > ();
        }
        if(_outputfile != "")
        {
            //char char_outputfile[1024] = {_outputfile}; // max. 1024 characters for filename
            //char_outputfile = _outputfile;
            cout << "Output into file: " << _outputfile.c_str() << "." << endl;
            freopen(_outputfile.c_str(),"w",stdout);   
            //cout << "hier ist output" << endl;
        }
        else
        {
            cout << " Output to screen." << endl;
        }
        
        _filename = OptionsMap()["file"].as<string > ();
        cout << " Database file: " << _filename << endl;        
        return true;
}
        

void KMCApplication::AddCalculator(KMCCalculator* calculator)
{
        _calculators.push_back(calculator);
}


void KMCApplication::Run()
{
    load_property_from_xml(_options, _op_vm["options"].as<string>());
    BeginEvaluate();
    EvaluateFrame();
    EndEvaluate();
}

void KMCApplication::BeginEvaluate(){
    list<KMCCalculator *>::iterator iter;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->Initialize(_filename.c_str(), &_options, _outputfile.c_str());
    }
}

bool KMCApplication::EvaluateFrame(){
    list<KMCCalculator *>::iterator iter;
    int i=0;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->EvaluateFrame();
    }
}

void KMCApplication::EndEvaluate(){
    list<KMCCalculator *>::iterator iter;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->EndEvaluate();
    }
}

}}
