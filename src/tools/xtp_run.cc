#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/xtp/sqlapplication.h>
#include <votca/xtp/calculatorfactory.h>


using namespace std;
using namespace votca::xtp;


class XtpRun : public SqlApplication
{
public:

    string  ProgramName() { return "xtp_run"; }    

    void    HelpText(ostream &out) { out <<"Runs charge transport calculators"<< endl; }
    void    HelpText() { };

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    
    //void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void XtpRun::Initialize() {

    SqlApplication::Initialize();

    AddProgramOptions("Calculators") ("execute,e", propt::value<string>(),
                      "List of calculators separated by ',' or ' '");
    AddProgramOptions("Calculators") ("list,l",
                      "Lists all available calculators");
    AddProgramOptions("Calculators") ("description,d", propt::value<string>(),
                      "Short description of a calculator");
}

bool XtpRun::EvaluateOptions() {

    string helpdir = "xtp/xml";
    
    if (OptionsMap().count("list")) {
            cout << "Available calculators: \n";
            for(Calculatorfactory::assoc_map::const_iterator iter=
                    Calculators().getObjects().begin();
                    iter != Calculators().getObjects().end(); ++iter) {
                PrintDescription( std::cout, (iter->first), helpdir, Application::HelpShort );
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
                for(Calculatorfactory::assoc_map::const_iterator iter=Calculators().getObjects().begin(); 
                        iter != Calculators().getObjects().end(); ++iter) {

                    if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                        PrintDescription( std::cout, (iter->first), helpdir, Application::HelpLong );
                        printerror = false;
                        break;
                    }
                 }
                 if ( printerror ) cout << "Calculator " << *n << " does not exist\n";
            }
            StopExecution();
            return true;     
    }

    SqlApplication::EvaluateOptions();
    CheckRequired("options", "Please provide an xml file with calculator options");
    CheckRequired("execute", "Nothing to do here: Abort.");

    Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = calcs.begin(); it != calcs.end(); it++) {
        SqlApplication::AddCalculator(Calculators().Create((*it).c_str()));
    }
    
    load_property_from_xml(_options, _op_vm["options"].as<string>());
    return true;
}

int main(int argc, char** argv) {
    
    XtpRun xtprun;
    return xtprun.Exec(argc, argv);

}
