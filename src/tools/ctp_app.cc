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
    
    CtpApp() {  }

    string  ProgramName() { return "ctp_app"; }    

    void    HelpText(ostream &out) { out <<"Runs charge transport tools"<< endl; }
    void    HelpText() { };
    void    PrintDescription(const char *name, const bool length);

    void    AddTool(QMTool *tool) { _tools.push_back(tool); }
    void    Initialize();
    bool    EvaluateOptions();
    void    Run(void);
    
    void BeginEvaluate(int nThreads);
    bool Evaluate();
    void EndEvaluate();
    
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
    QMToolFactory::RegisterAll();    

    // Tools-related
    AddProgramOptions("Tools") ("execute,e", propt::value<string>(),
                      "List of tools separated by ',' or ' '");
    AddProgramOptions("Tools") ("list,l",
                      "Lists all available tools");
    AddProgramOptions("Tools") ("description,d", propt::value<string>(),
                      "Short description of a tool");
    // Options-related
    AddProgramOptions() ("options,o", propt::value<string>(),
                         "  tool options");
    AddProgramOptions() ("nthreads,t", propt::value<int>()->default_value(1),
                         "  number of threads to create");
}

bool CtpApp::EvaluateOptions() {

    if (OptionsMap().count("list")) {
        cout << "Available tools: \n";
        for(QMToolFactory::assoc_map::const_iterator iter=
            QMTools().getObjects().begin();
            iter != QMTools().getObjects().end(); ++iter) {
            PrintDescription( (iter->first).c_str(), _short );
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
                    PrintDescription( (iter->first).c_str(), _long );
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
    CheckRequired("options", "Please provide an xml file with tool options");

    Tokenizer tools(OptionsMap()["execute"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = tools.begin(); it != tools.end(); it++) {
        this->AddTool(QMTools().Create((*it).c_str()));
    }
    return 1;
}


void CtpApp::PrintDescription(const char *name, const bool length) {
    // loading the documentation xml file from VOTCASHARE
    char *votca_share = getenv("VOTCASHARE");
    if(votca_share == NULL) throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
    string xmlFile = string(getenv("VOTCASHARE")) + string("/ctp/xml/")+name+string(".xml");
    try {
        Property options;
        load_property_from_xml(options, xmlFile);

       if ( length ) { // short description of the tool
            cout << string("  ") << _fwstring(string(name),14);
            cout << options.get(name+string(".description")).as<string>();
       } 
       else { // long description of the tool
            cout << " " << _fwstring(string(name),18);
            cout << options.get(name+string(".description")).as<string>() << endl;

            list<Property *> items = options.Select(name+string(".item"));

            for(list<Property*>::iterator iter = items.begin(); iter!=items.end(); ++iter) {
                //cout << "Long description" << endl;
                Property *pname=&( (*iter)->get( string("name") ) );
                Property *pdesc=&( (*iter)->get( string("description") ) );
                //Property *pdflt=&( (*iter)->get( string("default") ) );
                if ( ! (pname->value()).empty() ) {
                    string out_name = "  <" + pname->value() + ">";
                    cout << _fwstring(out_name, 20);
                    //cout << string("  <") << _fwstring(pname->value(), 20) << string(">");
                    cout << pdesc->value() << endl;
                }
             }
        }
        cout << endl;
    } catch(std::exception &error) {
        // cout << string("XML file or description tag missing: ") << xmlFile;
        cout << string("  ") << _fwstring(string(name),14);
        cout << "Undocumented" << endl;            
    }
}


void CtpApp::Run() {

    load_property_from_xml(_options, _op_vm["options"].as<string>());

    int nThreads = OptionsMap()["nthreads"].as<int>();
    
    cout << "Initializing tools " << endl;
    BeginEvaluate(nThreads);

    cout << "Evaluating tools " << endl;
    Evaluate();

    EndEvaluate();
}


void CtpApp::BeginEvaluate(int nThreads = 1) {
    list< QMTool* > ::iterator it;
    for (it = _tools.begin(); it != _tools.end(); it++) {
        cout << "... " << (*it)->Identify() << " ";
        (*it)->setnThreads(nThreads);
        (*it)->Initialize(&_options);        
        cout << endl;
    }
}

bool CtpApp::Evaluate() {
    list< QMTool* > ::iterator it;
    for (it = _tools.begin(); it != _tools.end(); it++) {
        cout << "... " << (*it)->Identify() << " " << flush;
        (*it)->Evaluate();
        cout << endl;
    }
}

void CtpApp::EndEvaluate() {
    list< QMTool* > ::iterator it;
    for (it = _tools.begin(); it != _tools.end(); it++) {
        (*it)->EndEvaluate();
    }
}


int main(int argc, char** argv) {
    
    CtpApp ctpapp;
    return ctpapp.Exec(argc, argv);

}
