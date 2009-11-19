// 
// File:   property.cc
// Author: ruehle
//
// Created on March 3, 2008, 15:13 PM
//

#include <iostream>
#include <boost/program_options.hpp>
#include <tools/property.h>
#include <tools/tokenizer.h>

using namespace std;


void help_text()
{
    cout << "csg_property\n\n";
    cout << "Helper program called by inverse scripts to parse xml file.\n\n";
}

int main(int argc, char** argv)
{      
    string filter, file, path, print;
    bool short_output = false;
    bool with_path = false;

    // lets read in some program options
    namespace po = boost::program_options;

    // Declare the supported options.
    po::options_description desc("Allowed options");    
    desc.add_options()
        ("help", "produce this help message")
        //("values", po::value<string>(&filter)->default_value(""),
        //    "list option values that match given criteria")
        ("path", po::value<string>(&path)->default_value(""),
            "list option values that match given criteria")
        ("filter", po::value<string>(&filter)->default_value(""),
            "list option values that match given criteria")
        ("print", po::value<string>(&print)->default_value(". "),
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
    if(vm.count("short"))
        short_output = true;
    if(vm.count("with-path"))
        with_path = true;
    
    Property p;
    load_property_from_xml(p, file);
    
    list<Property *> sel = p.Select(path);
    for(list<Property*>::iterator iter = sel.begin();
        iter!=sel.end(); ++iter) {
        if(filter!="") {
            Tokenizer tokenizer(filter, "=");
            Tokenizer::iterator tok;
            tok = tokenizer.begin();
            if(tok == tokenizer.end())
                throw std::invalid_argument("error, specified invalid filgter");
           
            string field = *tok;
            ++tok;
            if(tok == tokenizer.end()) 
                throw std::invalid_argument("error, specified invalid filgter");
            
            string value = *tok;
            if(!wildcmp(value.c_str(), (*iter)->get(field).value().c_str()))
                continue;
        }
        
        Property *p=&((*iter)->get(print));
       
        if(!short_output && with_path)
            cout << p->path() << ".";
        if(!short_output)
            cout << p->name() << " = ";
        if(!p->HasChilds())
            cout << p->value();
        cout << endl;
    }
    return 0;
}

