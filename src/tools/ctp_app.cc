#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/tools/application.h>
#include <votca/ctp/toolfactory.h>



using namespace std;
using namespace votca::ctp;


class CtpApp : public Application
{
public:

    string  ProgramName() { return "ctp_app"; }    

    void    HelpText(ostream &out) { out <<"Runs charge transport tools"<< endl; }
    void    HelpText() { };

    void    AddTool(QMTool *tool) { _tools.push_back(tool); }
    void    Initialize();
    bool    EvaluateOptions();
    
private:
    static const bool _short = true;
    static const bool _long = false;
    
    Property          _options;
    list< QMTool* >   _tools;

    string _fwstring(string original, size_t charCount ) {
        original.resize( charCount, ' ' );
        return original;
    }


};

namespace propt = boost::program_options;

void CtpApp::Initialize() {

    Application::Initialize();

    AddProgramOptions("Tools") ("execute,e", propt::value<string>(),
                      "List of tools separated by ',' or ' '");
    AddProgramOptions("Tools") ("list,l",
                      "Lists all available tools");
    AddProgramOptions("Tools") ("description,d", propt::value<string>(),
                      "Short description of a tool");
}

bool CtpApp::EvaluateOptions() {

    if (OptionsMap().count("list")) {
            cout << "Available tools: \n";
            for(QMToolFactory::assoc_map::const_iterator iter=
                    QMTools().getObjects().begin();
                    iter != QMTools().getObjects().end(); ++iter) {
                ; //PrintDescription( (iter->first).c_str(), _short );
            }
            StopExecution();
            return true;
    }
 
    
    if (OptionsMap().count("description")) {
            CheckRequired("description", "no tool is given");
 	    Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
            // loop over the names in the description string
            for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n) {
                // loop over tools
                bool printerror = true;
                for(QMToolFactory::assoc_map::const_iterator iter=QMTools().getObjects().begin(); 
                        iter != QMTools().getObjects().end(); ++iter) {

                    if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                        ; // PrintDescription( (iter->first).c_str(), _long );
                        printerror = false;
                        break;
                    }
                 }
                 if ( printerror ) cout << "Tool " << *n << " does not exist\n";
            }
            StopExecution();
            return true;
    }

    Application::EvaluateOptions();
    CheckRequired("execute", "Nothing to do here: Abort.");

    Tokenizer tools(OptionsMap()["execute"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = tools.begin(); it != tools.end(); it++) {
        this->AddTool(QMTools().Create((*it).c_str()));
    }
    return 1;
}


int main(int argc, char** argv) {
    
    CtpApp ctpapp;
    return ctpapp.Exec(argc, argv);

}
