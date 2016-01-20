#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/xtp/sqlapplication.h>
#include <votca/xtp/extractorfactory.h>


using namespace std;
using namespace votca::xtp;


class XtpDump : public SqlApplication
{
public:

    string  ProgramName() { return "xtp_dump"; }    

    void    HelpText(ostream &out) { out <<"Extracts information from the state file"<< endl; }
    void    HelpText() { };

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    
    //void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void XtpDump::Initialize() {
    ExtractorFactory::RegisterAll();
    SqlApplication::Initialize();

    AddProgramOptions("Extractors") ("extract,e", propt::value<string>(),
                      "List of extractors separated by ',' or ' '");
    AddProgramOptions("Extractors") ("list,l",
                      "Lists all available extractors");
    AddProgramOptions("Extractors") ("description,d", propt::value<string>(),
                      "Short description of an extractor");
}

bool XtpDump::EvaluateOptions() {

    if (OptionsMap().count("list")) {
            cout << "Available extractors: \n";
            for(ExtractorFactory::assoc_map::const_iterator iter=
                    Extractors().getObjects().begin();
                    iter != Extractors().getObjects().end(); ++iter) {
                    PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpShort );
            }
            StopExecution();
            return true;
    }
 
    
    if (OptionsMap().count("description")) {
            CheckRequired("description", "no extractor is given");
 	    Tokenizer tok(OptionsMap()["description"].as<string>(), " ,\n\t");
            // loop over the names in the description string
            for (Tokenizer::iterator n = tok.begin(); n != tok.end(); ++n) {
                // loop over calculators
                bool printerror = true;
                for(ExtractorFactory::assoc_map::const_iterator iter=Extractors().getObjects().begin(); 
                        iter != Extractors().getObjects().end(); ++iter) {

                    if ( (*n).compare( (iter->first).c_str() ) == 0 ) {
                        PrintDescription(std::cout, iter->first, "xtp/xml", Application::HelpLong );
                        printerror = false;
                        break;
                    }
                 }
                 if ( printerror ) cout << "Extractor " << *n << " does not exist\n";
            }
            StopExecution();
            return true;
    }

    SqlApplication::EvaluateOptions();
    CheckRequired("extract", "Nothing to do here: Abort.");

    Tokenizer calcs(OptionsMap()["extract"].as<string>(), " ,\n\t");
    Tokenizer::iterator it;
    for (it = calcs.begin(); it != calcs.end(); it++) {
        SqlApplication::AddCalculator(Extractors().Create((*it).c_str()));
    }
    return true;
}

int main(int argc, char** argv) {
    
    XtpDump xtpdump;
    return xtpdump.Exec(argc, argv);

}
