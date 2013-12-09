#include <stdlib.h>
#include <string>
#include <iostream>
#include <votca/ctp/sqlapplication.h>
#include <votca/ctp/extractorfactory.h>


using namespace std;
using namespace votca::ctp;


class CtpDump : public SqlApplication
{
public:

    string  ProgramName() { return "ctp_dump"; }    

    void    HelpText(ostream &out) { out <<"Extracts information from the state file"<< endl; }
    void    HelpText() { };

    void    Initialize();
    bool    EvaluateOptions();
    
private:
    
    //void    PrintDescription(string name, HelpOutputType _help_output_type);

};

namespace propt = boost::program_options;

void CtpDump::Initialize() {
    ExtractorFactory::RegisterAll();
    SqlApplication::Initialize();

    AddProgramOptions("Extractors") ("extract,e", propt::value<string>(),
                      "List of extractors separated by ',' or ' '");
    AddProgramOptions("Extractors") ("list,l",
                      "Lists all available extractors");
    AddProgramOptions("Extractors") ("description,d", propt::value<string>(),
                      "Short description of an extractor");
}

bool CtpDump::EvaluateOptions() {

    if (OptionsMap().count("list")) {
            cout << "Available extractors: \n";
            for(ExtractorFactory::assoc_map::const_iterator iter=
                    Extractors().getObjects().begin();
                    iter != Extractors().getObjects().end(); ++iter) {
                PrintDescription( std::cout, (iter->first), _helpShort );
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
                        PrintDescription( std::cout, (iter->first), _helpLong );
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
    
    CtpDump ctpdump;
    return ctpdump.Exec(argc, argv);

}
