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
#include <votca/tools/property.h>
#include <list>

using namespace std;
using namespace votca::tools;

void help_text()
{
    cout << "Helper program to parse xml file.\n\n";
}

int main(int argc, char** argv)
{      
    string file;
    bool short_output = false;
    bool with_path = false;

    // lets read in some program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Allowed options");    
    desc.add_options()
        ("help", "produce this help message")
        ("file", po::value<string>(&file), "xml file to parse");

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
        help_text();
        cout << desc << endl;
        return 0;
    }
    // file specified
    if (!vm.count("file")) {
        cout << "please specify file\n";                
        cout << desc << endl;
        return -1;
    }

    try {

    Property p;
    load_property_from_xml(p, file);
    
    cout << "loaded  "  << file << endl;

    list<Property *> sel = p.Select("*");

    Tokenizer tok("cg_molecule",".");
    Property pn;
    for( Tokenizer::iterator n = tok.begin(); n!=tok.end(); ++n) {
         pn = pn.get(*n);
    }


    cout << "info: " << p.name() <<":"<< p.path() <<":"<< p.value() << ":" << p.size() << endl;

        for(list<Property*>::iterator iter = sel.begin(); iter!=sel.end(); ++iter) {
        }
    } catch(std::exception &error) {
        cerr << "an error occurred:\n" << error.what() << endl;
        return -1;
    }
    return 0;
}

