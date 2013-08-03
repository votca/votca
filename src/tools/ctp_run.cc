#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/ctp/qmapplication.h>
#include <votca/ctp/calculatorfactory.h>
#include <boost/format.hpp>


using namespace std;
using namespace votca::ctp;
using boost::format;


class CtpRun : public QMApplication
{
public:

    string  ProgramName() { return "ctp_run"; }    

    void    HelpText(ostream &out) { out <<"Runs charge transport calculators"<< endl; }
    void    HelpText() { };

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    
    enum HelpOutputType { _helpShort, _helpLong };
    void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void CtpRun::Initialize() {

    QMApplication::Initialize();

    AddProgramOptions("Calculators") ("execute,e", propt::value<string>(),
                      "List of calculators separated by ',' or ' '");
    AddProgramOptions("Calculators") ("list,l",
                      "Lists all available calculators");
    AddProgramOptions("Calculators") ("description,d", propt::value<string>(),
                      "Short description of a calculator");
}

bool CtpRun::EvaluateOptions() {

    if (OptionsMap().count("list")) {
            cout << "Available calculators: \n";
            for(Calculatorfactory::assoc_map::const_iterator iter=
                    Calculators().getObjects().begin();
                    iter != Calculators().getObjects().end(); ++iter) {
                PrintDescription( (iter->first), _helpShort );
            }
            StopExecution();
            return true;
       //Application::StopExecution();
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
                        PrintDescription( (iter->first), _helpLong );
                        printerror = false;
                        break;
                    }
                 }
                 if ( printerror ) cout << "Calculator " << *n << " does not exist\n";
            }
            StopExecution();
            return true;     
        //cout << "Sorry... Note implemented." << endl;
        //Application::StopExecution();
    }

    QMApplication::EvaluateOptions();
    CheckRequired("execute", "Nothing to do here: Abort.");

    Tokenizer calcs(OptionsMap()["execute"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = calcs.begin(); it != calcs.end(); it++) {
        QMApplication::AddCalculator(Calculators().Create((*it).c_str()));
    }
    return 1;
}

void CtpRun::PrintDescription(string name, HelpOutputType _help_output_type) {
    // loading the documentation xml file from VOTCASHARE
    char *votca_share = getenv("VOTCASHARE");
    if (votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/") + name + string(".xml");
    string _format("%|3t|%1% %|20t|%2% \n");
    try {
        Property options;
        load_property_from_xml(options, xmlFile);

        string _help;
        switch (_help_output_type) {

            case _helpShort:
                _help = options.get("options." + name).getAttribute<string>("help");
                cout << format(_format) % name % _help;
                break;

            case _helpLong:
                cout << HLP << setlevel(2) << options;
        }

    } catch (std::exception &error) {
        cout << format(_format) % name % "Undocumented";
    }
}



int main(int argc, char** argv) {
    
    CtpRun ctprun;
    return ctprun.Exec(argc, argv);

}
