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

namespace votca { namespace kmc {

KMCApplication::KMCApplication()
{
    KMCCalculatorFactory::RegisterAll();
}

KMCApplication::~KMCApplication()
{}

void KMCApplication::Initialize(void)
{
    //Application::Initialize();
    KMCCalculatorFactory::RegisterAll();
    AddProgramOptions()
            ("options,o", boost::program_options::value<string>(), "  program and calculator options")
            ("file,f", boost::program_options::value<string>(), "  sqlite state file");

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

    ///votca::ctp::HelpTextHeader(name);
    HelpText(out);
    out << "\n\n" << OptionsDesc() << endl;
}

/*
void KMCApplication::PrintDescription(const char *name, const bool length)
{
}
*/

// check if required options are provided
bool KMCApplication::EvaluateOptions() {

    if(OptionsMap().count("list")) {
        cout << "Available calculators: \n";
        for(KMCCalculatorFactory::assoc_map::const_iterator iter=Calculators().getObjects().begin();
            iter != Calculators().getObjects().end(); ++iter) {
            //PrintDescription( (iter->first).c_str(), _short );
        }
        StopExecution();
        return true;
    }

    CheckRequired("options", "please provide an xml file with the program options");
    CheckRequired("file", "no database file specified");

    return true;
}

void KMCApplication::AddCalculator(KMCCalculator* calculator){
        _calculators.push_back(calculator);
}

void KMCApplication::Run()
{
}

void KMCApplication::BeginEvaluate(){
    list<KMCCalculator *>::iterator iter;
    for (iter = _calculators.begin(); iter != _calculators.end(); ++iter){
        (*iter)->Initialize(&_options);
    }
}

bool KMCApplication::EvaluateFrame(){
    list<KMCCalculator *>::iterator iter;
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
