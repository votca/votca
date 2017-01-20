/* 
 *            Copyright 2009-2016 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/tools/property.h>
#include <votca/xtp/xtpapplication.h>
#include <votca/xtp/toolfactory.h>
#include <votca/ctp/toolfactory.h>



using namespace std;
using namespace votca::ctp;
using namespace votca::tools;

//class XtpTools : public votca::ctp::XtpApplication
class XtpTools : public votca::xtp::XtpApplication
{
public:
    
    XtpTools() { votca::ctp::QMToolFactory::RegisterAll(); }

    string  ProgramName() { return "xtp_tools"; }    

    void    HelpText(ostream &out) { out <<"Runs excitation/charge transport tools"<< endl; }

    void    AddTool(votca::ctp::QMTool *tool) { _tools.push_back(tool); }
    //void    AddTool(QMTool *tool) { _ctp_tools.push_back(tool); }


    void    Initialize();
    bool    EvaluateOptions();
    void    Run(void);
    
    void BeginEvaluate(int nThreads);
    bool Evaluate();
    void EndEvaluate();
    
private:
    
    votca::tools::Property _options;
    list< votca::ctp::QMTool* >   _tools;
    //list< votca::ctp::QMTool* >   _ctp_tools;
};



void XtpTools::Initialize() {
    
    votca::xtp::XQMToolFactory::RegisterAll(); 
    votca::ctp::QMToolFactory::RegisterAll();

    namespace propt = boost::program_options;    
    // Tools-related
    AddProgramOptions("Tools") ("execute,e", propt::value<string>(),
                      "List of tools separated by ',' or ' '");
    AddProgramOptions("Tools") ("list,l",
                      "Lists all available tools");
    AddProgramOptions("Tools") ("description,d", propt::value<string>(),
                      "Short description of a tool");
    // Options-related
    AddProgramOptions() ("nthreads,t", propt::value<int>()->default_value(1),
                         "  number of threads to create");
    AddProgramOptions() ("options,o", propt::value<string>(),
                         "  calculator options");
}


bool XtpTools::EvaluateOptions() {

    namespace XTP=votca::xtp;
    
    if (OptionsMap().count("list")) {
        cout << "Available XTP tools: \n";
        for(XTP::XQMToolFactory::assoc_map::const_iterator iter=
            XTP::XQMTools().getObjects().begin();
            iter != XTP::XQMTools().getObjects().end(); ++iter) {
            PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpShort );
        }
        cout << "Available (wrapped) CTP tools: \n";
        
        // also include the CTP Tools
        for(CTP::QMToolFactory::assoc_map::const_iterator iter=
            QMTools().getObjects().begin();
            iter != QMTools().getObjects().end(); ++iter) {
            PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpShort );
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
            for(XTP::XQMToolFactory::assoc_map::const_iterator iter=XTP::XQMTools().getObjects().begin(); 
                iter != XTP::XQMTools().getObjects().end(); ++iter) {

                if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                    PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpLong );
                    printerror = false;
                    break;
                }
             }
            
            // also check CTP tools
            for(votca::ctp::QMToolFactory::assoc_map::const_iterator iter=votca::ctp::QMTools().getObjects().begin(); 
                iter != votca::ctp::QMTools().getObjects().end(); ++iter) {

                if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                    PrintDescription(std::cout, iter->first, "ctp/xml", Application::HelpLong );
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
       

        // check if XTP or CTP tool
        for(XTP::XQMToolFactory::assoc_map::const_iterator iter=XTP::XQMTools().getObjects().begin(); 
                iter != XTP::XQMTools().getObjects().end(); ++iter) {

                if ( (*it).compare( (iter->first).c_str() ) == 0 ) {
                    cout << "Registered XTP " << (*it).c_str() << endl;
                    this->AddTool(XTP::XQMTools().Create((*it).c_str()));
                    //PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpLong );
                    //printerror = false;
                    break;
                }
            }


        for(votca::ctp::QMToolFactory::assoc_map::const_iterator iter=votca::ctp::QMTools().getObjects().begin(); 
                iter != votca::ctp::QMTools().getObjects().end(); ++iter) {

                if ( (*it).compare( (iter->first).c_str() ) == 0 ) {
                    cout << "Registered CTP " << (*it).c_str() << endl;
                    this->AddTool(votca::ctp::QMTools().Create((*it).c_str()));
                    //PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpLong );
                    //printerror = false;
                    break;
                }
            }
        
    }

    
    return 1;
}


void XtpTools::Run() {

    string optionsFile = _op_vm["options"].as<string>();    
    load_property_from_xml(_options, optionsFile);   
    
    int nThreads = OptionsMap()["nthreads"].as<int>();
    
    cout << "Initializing tools " << endl;
    BeginEvaluate(nThreads);

    cout << "Evaluating tools " << endl;
    Evaluate();

    EndEvaluate();
}


void XtpTools::BeginEvaluate(int nThreads = 1) {
    list< votca::ctp::QMTool* > ::iterator it;
    for (it = _tools.begin(); it != _tools.end(); it++) {
        cout << "... " << (*it)->Identify() << " " << flush;
        (*it)->setnThreads(nThreads);
        (*it)->Initialize(&_options);        
        cout << endl;
    }
    
    // CTP tools 
    /* list< votca::ctp::QMTool* > ::iterator cit;
    for (cit = _ctp_tools.begin(); cit != _ctp_tools.end(); cit++) {
        cout << "... " << (*cit)->Identify() << " " << flush;
        (*cit)->setnThreads(nThreads);
        (*cit)->Initialize(&_options);        
        cout << endl;
    }  */  
}

bool XtpTools::Evaluate() {
    list< votca::ctp::QMTool* > ::iterator it;
    for (it = _tools.begin(); it != _tools.end(); it++) {
        cout << "... " << (*it)->Identify() << " " << flush;
        (*it)->Evaluate();
        cout << endl;
    }
    
    return true;
}

void XtpTools::EndEvaluate() {
    list< votca::ctp::QMTool* > ::iterator it;
    for (it = _tools.begin(); it != _tools.end(); it++) {
        (*it)->EndEvaluate();
    }
}


int main(int argc, char** argv) {
    
    XtpTools xtpapp;
    return xtpapp.Exec(argc, argv);

}
