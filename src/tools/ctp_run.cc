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
        out << "Runs a specified calculator(s)." << endl;
    }

    void Initialize() {
        QMApplication::Initialize();
        AddProgramOptions("calculator execution")
            ("exec", boost::program_options::value<string>(), "list of calculators separated by commas or spaces")
            ("list", "show available calculators");
    }

    //TODO: Support for XML-File based options
    bool EvaluateOptions() {
        if(OptionsMap().count("list")) {
            cout << "Available calculators:\n\n";
            for(CalculatorFactory::assoc_map::const_iterator iter=Calculators().getObjects().begin();
                    iter != Calculators().getObjects().end(); ++iter) {
                QMCalculator *tmp = iter->second();                
                cout << " * " << iter->first << ", " << tmp->Description() << endl;
                delete tmp;
            }
            cout << "\n";
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
};

int main(int argc, char** argv) {
    QMAppRun qmapprun;
    return qmapprun.Exec(argc, argv);
}
