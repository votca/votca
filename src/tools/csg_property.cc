// 
// File:   property.cc
// Author: ruehle
//
// Created on March 3, 2008, 15:13 PM
//

#include <iostream>
#include <boost/program_options.hpp>
#include <tools/property.h>

using namespace std;

int main(int argc, char** argv)
{      
    string filter, file;
    bool short_output = false;
    bool with_path = false;

    // lets read in some program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Allowed options");    
    desc.add_options()
        ("help", "produce this help message")
        ("values", po::value<string>(&filter)->default_value(""),
            "list option values that match given criteria")
        ("file", po::value<string>(&file), "xml file to parse")
        ("short", "short version of output")
        ("with-path", "include path of node in output");

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
        cout << "csg_property\n\n";                
        cout << desc << endl;
        return 0;
    }
    // file specified
    if (!vm.count("file")) {
        cout << "please specify file\n";                
        cout << desc << endl;
        return -1;
    }
    if(vm.count("short"))
        short_output = true;
    if(vm.count("with-path"))
        with_path = true;
    
    Property p;
    load_property_from_xml(p, file);
    
    list<Property *> sel = p.Select(filter);
    for(list<Property*>::iterator iter = sel.begin();
        iter!=sel.end(); ++iter) {
        if(!short_output && with_path)
            cout << (*iter)->path() << ".";
        if(!short_output)
            cout << (*iter)->name() << " = ";
        if(!(*iter)->HasChilds())
            cout << (*iter)->value();
        cout << endl;
    }
    return 0;
}

