/* 
 * File:   ctp_run.cc
 * Author: schrader
 *
 * Created on July 12, 2010, 3:11 PM
 */

#include <stdlib.h>
#include "qmapplication.h"
#include <calculatorfactory.h>
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
            ("exec", boost::program_options::value<string>(), "list of calculators separated by commas or spaces");
    }

    //TODO: Support for XML-File based options
    bool EvaluateOptions() {
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