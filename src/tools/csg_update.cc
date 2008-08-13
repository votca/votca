// 
// File:   csg_update.cc
// Author: ruehle
//
// Created on June 11, 2008, 10:31 AM
//


/// this is the most dirtyest program, clean it up, don't copy anything from here!!!

#define kB  8.3109*0.01

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>
#include <tools/table.h>


using namespace std;

namespace po = boost::program_options;

int DoIBM(const string &in, const string &out, const string &target, const string &current, double T, double scale);

int main(int argc, char** argv)
{
    string method, type;
    string action;
    string in, out;
    string target, current;
    double T, scale;
    
    // lets read in some program options
    // Declare the supported options.
    po::options_description desc("Allowed options");    

    desc.add_options()
    ("help", "produce this help message")
    //("version", "show version info")
    ("in", boost::program_options::value<string>(&in), "file containing current potential")
    ("out", boost::program_options::value<string>(&out)->default_value("out.dat"), "file to write the new potential")    
    ("cur", boost::program_options::value<string>(&current), "file containing current rdf")
    ("target", boost::program_options::value<string>(&target), "file containing target rdf")
    ("action", boost::program_options::value<string>(&action)->default_value("ibm"), "ibm (to do: smooth, imc,...)")
    ("T", boost::program_options::value<double>(&T)->default_value(300.), "temperature");
    ("scale", boost::program_options::value<double>(&scale)->default_value(1.), "correction scaling");

//    ("type", boost::program_options::value<string>()->default_value("nb"), "nb, bond, ang, dih")
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

    if(action == "ibm") {
        return DoIBM(in, out, target, current, T, scale);
    }
    else
        cerr << "unknown action\n";

    return -1;
}


// do inverse boltzmann
int DoIBM(const string &in, const string &out, const string &target_dist, const string &current, double T, double scale)
{
    Table target;     
    Table pout;
    // 
    if(target_dist=="") {
        cerr << "error, not target given for iterative inverse boltzmann"; 
        return -1;
    }

    target.Load(target_dist);

    // if no input potential is given, do initial guess
    if(in=="") {
        pout.resize(target.size());
        pout.x() = target.x();
	pout.flags() = ub::scalar_vector<unsigned short>(target.size(), 0);

        for(int i=0; i<pout.size(); ++i) {            
            if(target.y(i) == 0) {
                pout.y(i) = 0; 
                pout.flags(i) |= TBL_INVALID;
            }
            else
                pout.y(i) = -kB*T*log(target.y(i));
        }
    }
    // otherwise do ibm update
    else {
        Table pin;
                
        
        if(current == "") {
            cerr << "error, give current distribution";
            return -1;
	}
        
        // read in the current potential
        pin.Load(in);        
        // read the current distribution
        Table cur;     
        cur.Load(current);

        pout = pin;
        
        for(int i=0; i<pout.size(); ++i) {            
            if(target.y(i) == 0 || cur.y(i) == 0) {
                pout.y(i) = 0; 
                pout.flags(i) |= TBL_INVALID;
            }
            else
                pout.y(i) += -kB*T*log(cur.y(i) / target.y(i));
        }                
    }
    
    pout.Save(out);
    // should not get here
    return 0;
}


