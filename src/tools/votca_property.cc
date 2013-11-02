/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string/replace.hpp>

#include <votca/tools/version.h>
#include <votca/tools/property.h>
#include <votca/tools/propertyiomanipulator.h>
#include <list>

using namespace std;
using namespace votca::tools;

int main(int argc, char** argv)
{      
    string program_name = "votca_property";
    string help_text = "Helper for parsing xml files";
    string file;
    string format;
    int level;

    // read in the program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Program options");    
    desc.add_options()
        ("help", "produce this help message")
        ("man", "man pages")
        ("file", po::value<string>(&file), "xml file to parse")
        ("format", po::value<string>(&format), "output format [XML TXT TEX]")
        ("level", po::value<int>(&level), "output from this level ");

    // now read in the command line
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);    
        po::notify(vm);
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }
    // does the user want help?
    if (vm.count("help")) {
        cout << help_text << "\n\n";
        cout << desc << endl;
        return 0;
    }

    // does the user want man pages?
    if (vm.count("man")) {
        string option_format(".TP\n\\fB%1%\\fR\n%2%\n");
        string header_format(".TH \"%1%\" 1 \"\" \"Version: %2%\"\n\n");
        string name_format(".SH NAME\n\n.P\n%1% \\- Part of the VOTCA package\n"
                           "\n.P\nPlease visit the program site at __%2%__\n\n"
        );
        string footer_format ("\n.SH AUTHORS\n"
                              "\n.P\nWritten and maintained by the VOTCA Development Team <%1%>\n"
                              "\n.SH COPYRIGHT\n\n.P\n\n"
                              "Copyright 2009\\-2011 The VOTCA Development Team (%2%)\n"
                              "\n.P\nLicensed under the Apache License, Version 2.0 (the \"License\") "
                              "you may not use this file except in compliance with the License. "
                              "You may obtain a copy of the License at"
                              "\n.P\nhttp://www.apache.org/licenses/LICENSE\\-2.0\n"
                              "\n.P\nUnless required by applicable law or agreed to in writing, software "
                              "distributed under the License is distributed on an \"AS IS\" BASIS, "
                              "WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. "
                              "See the License for the specific language governing permissions and "
                              "limitations under the License.");
        
        
        string synopsis_format ("\n.SH SYNOPSIS\n\n.P\n\\fB%1%\\fR [\\fIOPTION\\fR] [\\fIPARAMETERS\\fR]\n");
        string helptext_format ("\n.SH DESCRIPTION\n\n.P\n%1%\n\n.SH OPTIONS\n\n");
        string url = "http://www.votca.org"; 
        string email = "devs@votca.org";
        
        std::cout << boost::format(header_format) %  program_name % ToolsVersionStr();        
        std::cout << boost::format(name_format) % program_name % url;
        std::cout << boost::format(synopsis_format) % program_name;       
        std::cout << boost::format(helptext_format) % help_text;

        typedef std::vector<boost::shared_ptr<boost::program_options::option_description> >::const_iterator OptionsIterator;
        OptionsIterator it = desc.options().begin(), it_end = desc.options().end();
        while(it < it_end) {
            
            string format_name = (*it)->format_name() + " " + (*it)->format_parameter();
            boost::replace_all(format_name, "-", "\\-");
            std::cout << boost::format(option_format) % format_name % (*it)->description();
            ++it;
            
            //std::cout << " longname: " << (*it)->long_name() << endl;
            //std::cout << " description: " << (*it)->description() << endl;
            //std::cout << " format: " << (*it)->format_parameter() << endl;
            //std::cout << " formatname: " << (*it)->format_name() << endl << endl;
            
        }      

        std::cout << boost::format(footer_format) % email % url;

        return 0;
    }    // file specified

    if (!vm.count("file")) {
        cout << "please specify file\n";                
        cout << desc << endl;
        return -1;
    }
    // format specified
    if (!vm.count("format")) {
        cout << "format not specified, using XML\n";  
        format = "XML";
    } 
    // format specified
    if (!vm.count("level")) {
        level = 1;
    } 
    
    try {

    Property p;
    
    map<string, PropertyIOManipulator* > _mformat;
    map<string, PropertyIOManipulator* >::iterator it;

    _mformat["XML"] = &XML;
    _mformat["TXT"] = &TXT;
    _mformat["TEX"] = &TEX;
    _mformat["HLP"] = &HLP;
    
    load_property_from_xml(p, file);
    
    it = _mformat.find( format );
    if ( it != _mformat.end() ) {
        PropertyIOManipulator *piom = _mformat.find( format )->second;
        piom->setLevel(level);
        piom->setIndentation("");
        piom->setColorScheme<csRGB>();
        cout << *piom  << p ;
    } else {
        cout << "format " << format << " not supported \n";
        cout << desc << endl;
        return -1;
    }

    //PropertyIOManipulator XML(PropertyIOManipulator::XML,0,"---"); 
    //cout << XML << p;
    //cout << TXT << p;
    //cout << T2T << p;
    //cout << LOG << p;
    //cout << TEX << p;
            
    } catch(std::exception &error) {
        cerr << "an error occurred:\n" << error.what() << endl;
        return -1;
    }
    return 0;
}

