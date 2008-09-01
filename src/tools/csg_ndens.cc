/// \addtogroup csg
///@{
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
#include <cgengine.h>
#include <libversion.h>
#include <numberdist.h>
#include <imccoefficients.h>
// #include <modules/cg/bondedstatistics.h>

using namespace std;

class CGDistNb
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {
        _ndist.setCutoff(2);
        _ndist.setN(200);
        _imccoeff.setN(200);
    };
    void EndCG() {
//        cout << _ndist << endl;
        cout <<  endl << endl;
        _imccoeff.Average();
        _imccoeff.OutputDist(cout);
        cout << endl << endl;
        _imccoeff.OutputCross(cout);
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
        _ndist.clear();
        _ndist.Process(*conf);
        _imccoeff.Process(_ndist.getDist());
    }
    
protected:
    NumberDist _ndist;
    IMCCoefficients _imccoeff;
};

/*class CGDistBonded
    : public CGBondedStatistics
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {
        _imccoeff.setN(200);
        Histogram::options_t op;
        op.setN(200);
        _hist.
    };
    void EndCG() {
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
        _ndist.clear();
        _ndist.Process(*conf);
        _imccoeff.Process(_ndist.getDist());
    }
    
protected:
    IMCCoefficients _imccoeff;
};*/

int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGDistNb no;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&no);


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
        cout << "csg_nemat, lib version " << LIB_VERSION_STR << "\n\n";                
        cout << desc << endl;
        return 0;
    }
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR  << "\n";                        
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

/// @}
