// 
// File:   csg_update.cc
// Author: ruehle
//
// Created on June 11, 2008, 10:31 AM
//


/// this is the most dirtyest program, clean it up, don't copy anything from here!!!

#define kbT  8.3109*300*0.01

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <table.h>

using namespace std;

void ReadMatrix(ifstream &in, vector<double> &matrix)
{
    int N;
    in >> N;
    
    cout << "entries: " << N << endl;
    matrix.resize(N*(N+1)/2);
    for(int i=0; i<N*(N+1)/2; ++i)
        in >> matrix[i];        
}

void InitialGuess(Table &out, Table &in)
{
    out.resize(in.size());
    
    for(int i=0; i<in.size(); ++i) {
        out.x(i)=in.x(i);
        if(in.y(i)>0)
            out.y(i) = - kbT*log(in.y(i));
        else
            out.y(i) = 10000;
    }
}

int main(int argc, char** argv)
{    
    string method, type;
    ifstream in;
    
    // lets read in some program options
    namespace po = boost::program_options;
    // Declare the supported options.
    po::options_description desc("Allowed options");    
         
    desc.add_options()
    ("help", "produce this help message")
    //("version", "show version info")
    ("in", boost::program_options::value<string>(), "file containing number distribution")
    ("pin", boost::program_options::value<string>(), "file containing number distribution")
    ("target", boost::program_options::value<string>(), "file containing target")
    ("method", boost::program_options::value<string>()->default_value("ibm"), "either ibm or imc")            
    ("type", boost::program_options::value<string>()->default_value("nb"), "nb, bond, ang, dih")
    ;
    
    // now read in the command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_update \n\n";                
        cout << desc << endl;
        return 0;
    }
    
    if (!vm.count("target")) {
        cout << desc << endl;
        cout << "give file with target distribution";
        return 1;
    }

    if (!vm.count("target")) {
        cout << desc << endl;
        cout << "give target distribution";
        return 1;
    }
    
    Table dist;     
    vector<double> matrix;     
    Table target_dist;     
    Table pot_in;
    Table pot_out;
        
    // Read in the target distribution
    target_dist.Load(vm["target"].as<string>());
    
    // if no input potential is given, make initial guess
    if (!vm.count("pin")) {
        InitialGuess(pot_out, target_dist);
    }
    // otherwise correct the potential
    else {
        // read in the current potential
        pot_in.Load(vm["pin"].as<string>());
        
        // read the current distribution
        in.open(vm["in"].as<string>().c_str());
        if(!in) {
            cerr << "cannot open in file";
            return 1;
        }

        in >> dist;
        
        if(method == "ibm") {
        //    UpdateIBM();
        }
        else if(method == "imc") {
         //   UpdateIMC();
            ReadMatrix(in, matrix);
        }
        else
            cerr << "unknown update method\n";
        in.close();            
    }

    cout << pot_out << endl;
    return 0;
}

