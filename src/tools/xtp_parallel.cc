#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/xtp/jobapplication.h>
#include <votca/xtp/jobcalculatorfactory.h>


using namespace std;
using namespace votca::xtp;


class XtpParallel : public JobApplication
{
public:

    string  ProgramName() { return "xtp_parallel"; }    

    void    HelpText(ostream &out) { out <<"Runs job-based heavy-duty calculators"<< endl; }
    void    HelpText() { };

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    
    //void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void XtpParallel::Initialize() {

    JobApplication::Initialize();

    AddProgramOptions("Calculators") ("execute,e", propt::value<string>(),
                      "List of calculators separated by ',' or ' '");
    AddProgramOptions("Calculators") ("list,l",
                      "Lists all available calculators");
    AddProgramOptions("Calculators") ("description,d", propt::value<string>(),
                      "Short description of a calculator");
}

bool XtpParallel::EvaluateOptions() {

    if (OptionsMap().count("list")) {
            cout << "Available calculators: \n";
            for(JobCalculatorfactory::assoc_map::const_iterator iter=
                    JobCalculators().getObjects().begin();
                    iter != JobCalculators().getObjects().end(); ++iter) {
                    PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpShort );
            }
            StopExecution();
            return true;
    }
 
    
    if (OptionsMap().count("description")) {
            CheckRequired("description", "no calculator is given");
 	    Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
            // loop over the names in the description string
            for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n) {
                // loop over calculators
                bool printerror = true;
                for(JobCalculatorfactory::assoc_map::const_iterator iter=JobCalculators().getObjects().begin(); 
                        iter != JobCalculators().getObjects().end(); ++iter) {

                    if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                        PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpLong ); 
                        printerror = false;
                        break;
                    }
                 }
                 if ( printerror ) cout << "Calculator " << *n << " does not exist\n";
            }
            StopExecution();
            return true;     
    }

    JobApplication::EvaluateOptions();
    CheckRequired("execute", "Nothing to do here: Abort.");

    Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = calcs.begin(); it != calcs.end(); it++) {
        JobApplication::AddCalculator(JobCalculators().Create((*it).c_str()));
    }
    return true;
}

int main(int argc, char** argv) {
    
    XtpParallel xtprun;
    return xtprun.Exec(argc, argv);

}
