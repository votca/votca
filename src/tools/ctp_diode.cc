#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/ctp/qmapplication.h>
#include <votca/ctp/calculatorfactory.h>



using namespace std;
using namespace votca::ctp;


class CtpDiode : public QMApplication
{
public:

    string  ProgramName() { return "ctp_diode"; }    

    void    HelpText(ostream &out) { out <<"Runs charge transport calculators"<< endl; }
    void    HelpText() { };
    void    PrintDescription(const char *name, const bool length);

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    static const bool _short = true;
    static const bool _long = false;

    string _fwstring(string original, size_t charCount ) {
        original.resize( charCount, ' ' );
        return original;
    }


};

namespace propt = boost::program_options;

void CtpDiode::Initialize() {

    QMApplication::Initialize();

    AddProgramOptions("Calculators") ("execute,e", propt::value<string>(),
                      "List of calculators separated by ',' or ' '");
    AddProgramOptions("Calculators") ("list,l",
                      "Lists all available calculators");
    AddProgramOptions("Calculators") ("description,d", propt::value<string>(),
                      "Short description of a calculator");
}

bool CtpDiode::EvaluateOptions() {

    if (OptionsMap().count("list")) {
            cout << "Available calculators: \n";
            for(Calculatorfactory::assoc_map::const_iterator iter=
                    Calculators().getObjects().begin();
                    iter != Calculators().getObjects().end(); ++iter) {
                PrintDescription( (iter->first).c_str(), _short );
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
                         PrintDescription( (iter->first).c_str(), _long );
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

void CtpDiode::PrintDescription(const char *name, const bool length) {
        // loading the documentation xml file from VOTCASHARE
        char *votca_share = getenv("VOTCASHARE");
        if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
        string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/")+name+string(".xml");
        try {
            Property options;
            load_property_from_xml(options, xmlFile);

           if ( length ) { // short description of the calculator
                 cout << string("  ") << _fwstring(string(name),14);
                 cout << options.get("options."+string(name)).getAttribute<string>("help");
            } else { // long description of the calculator
                cout << HLP << setlevel(2) << options;
            }
            cout << endl;
        } catch(std::exception &error) {
            // cout << string("XML file or description tag missing: ") << xmlFile;
            cout << string("  ") << _fwstring(string(name),14);
            cout << "Undocumented" << endl;
            
        }
}



int main(int argc, char** argv) {
    
    CtpDiode ctprun;
    return ctprun.Exec(argc, argv);

}
