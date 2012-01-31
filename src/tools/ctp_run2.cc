#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/ctp/qmapplication2.h>
#include <votca/ctp/calculatorfactory2.h>



using namespace std;
using namespace votca::ctp;


class CtpRun : public QMApplication2
{
public:

    string  ProgramName() { return "ctp_run"; }    

    void    HelpText(ostream &out) { out <<"Runs CTP calculators"<< endl; }
    void    HelpText() { };
    void    PrintDescription();

    void    Initialize();
    bool    EvaluateOptions();

};

namespace propt = boost::program_options;

void CtpRun::Initialize() {

    QMApplication2::Initialize();

    AddProgramOptions("Calculators") ("execute,e", propt::value<string>(),
                      "List of calculators separated by ',' or ' '");
    AddProgramOptions("Calculators") ("list,l",
                      "Lists all available calculators");
    AddProgramOptions("Calculators") ("description,d", propt::value<string>(),
                      "Short description of a calculator");
}

bool CtpRun::EvaluateOptions() {

    if (OptionsMap().count("list")) {
        cout << "Sorry... Not implemented." << endl;
        Application::StopExecution();
    }
    if (OptionsMap().count("description")) {
        cout << "Sorry... Note implemented." << endl;
        Application::StopExecution();
    }

    QMApplication2::EvaluateOptions();
    CheckRequired("execute", "Nothing to do here: Abort.");

    Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = calcs.begin(); it != calcs.end(); it++) {
        QMApplication2::AddCalculator(Calculators().Create((*it).c_str()));
    }
    return 1;
}

void CtpRun::PrintDescription() {
    cout << "Sorry... Not implemented." << endl;
}



int main(int argc, char** argv) {

    std::cout << "CTP_RUN new version..." << std::endl;
    CtpRun ctprun;
    return ctprun.Exec(argc, argv);

}
