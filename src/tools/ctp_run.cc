/* 
 * File:   ctp_run.cc
 * Author: schrader
 *
 * Created on July 12, 2010, 3:11 PM
 */

#include <stdlib.h>
#include <votca/ctp/qmapplication.h>
#include <votca/ctp/calculatorfactory.h>
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
