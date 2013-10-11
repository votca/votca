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
#include <votca/tools/propertyformat.h>
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
    string format;
    int level;

    // read in the program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Program options");    
    desc.add_options()
        ("help", "produce this help message")
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
    
    map<string, PropertyFormat* > _mformat;
    map<string, PropertyFormat* >::iterator it;
     
    _mformat["XML"] = &XML;
    _mformat["TXT"] = &TXT;
    _mformat["LOG"] = &LOG;
    _mformat["T2T"] = &T2T;
    _mformat["TEX"] = &TEX;
    _mformat["HLP"] = &HLP;
    
    load_property_from_xml(p, file);
    
    it = _mformat.find( format );
    if ( it != _mformat.end() ) {
        cout << *(_mformat.find( format )->second) << setlevel(level) << p ;
    } else {
        cout << "format " << format << " not supported \n";
        cout << desc << endl;
        return -1;
    }
    
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

