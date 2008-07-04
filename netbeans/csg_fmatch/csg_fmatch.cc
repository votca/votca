// 
// File:   csg_nemat.cc
// Author: ruehle
//
// Created on March 6, 2008, 4:35 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include "cgengine.h"
#include "libversion.h"
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

namespace ub = boost::numeric::ublas;
using namespace std;

class CGForceMatching
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom)
    {
        cout << "hey, someone wants to coarse grain\n";
    }
    
    void EndCG() {
        cout << "write out results\n";
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
        cout << "yea, new configuration!\n";
        InteractionContainer &ic = conf->getTopology()->getBondedInteractions();
        InteractionContainer::iterator ia;
        
        for(ia=ic.begin(); ia != ic.end(); ++ia) {            
            const string &name = (*ia)->getName();             
            cout << name << ": " << ((*ia)->EvaluateVar(*conf)) << endl;
            //break;
        }
    }
    
protected:
    // _A*_x = _b
    
  ub::mapped_matrix<double> _A;
  ub::vector<double> _b; // F_ref
  ub::vector<double> _x; // 
    // _A(i, j) = 10;
    // _A.resize(n, m, true);
};


int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGForceMatching fmatch;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&fmatch);


    // initialize the readers/writers,
    // this will be combined in an initialize function later
    TrajectoryWriter::RegisterPlugins();
    TrajectoryReader::RegisterPlugins();
    TopologyReader::RegisterPlugins();

    
    // lets read in some program options
    namespace po = boost::program_options;
        
    
    // Declare the supported options.
    po::options_description desc("Allowed options");    
    
    // let cg_engine add some program options
    cg_engine.AddProgramOptions(desc);
    
    // now read in the command line
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    // does the user want help?
    if (vm.count("help")) {
        cout << "csg_fmatch, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_fmatch, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    
    // try to run the cg process, go through the frames, etc...
    try {
        cg_engine.Run(desc, vm);
    }
    // did an error occour?
    catch(string error) {
        cerr << "An error occoured!" << endl << error << endl;
    }
    return 0;
}

