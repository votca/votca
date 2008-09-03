// 
// File:   template.cc
// Author: ruehle
//
// Created on June 8, 2008, 10:41 PM
//

#include <math.h>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <cgengine.h>
#include <libversion.h>
#include <math.h>

using namespace std;

double thres;

class CGConj
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {
         for(int i=0;i<10;++i)
             _dist[i]=0;
    };
    void EndCG() {
        for(int i=0;i<10;++i) {
           cout << i+1 << " " << _dist[i] << endl;
        }
    };
    
    void EvalConfiguration(Configuration *conf, Configuration *conf_atom = 0) {
        Topology *top = conf->getTopology();
        for(MoleculeContainer::iterator iter=top->getMolecules().begin();
            iter!=top->getMolecules().end(); ++iter) {
            MoleculeInfo *mi;
            mi = *iter;
            vec no = conf->getU(mi->getBeadId(0));
            int c=0;
            for(int i=1; i<(*iter)->BeadCount(); ++i) {
                vec n = conf->getU(mi->getBeadId(i));
                ++c;
                // cout << acos(fabs(no*n)) << endl ;
                if(fabs(no*n) < thres) {
                    _dist[c-1]++;
                    c=0;
                }
                no = n; 
            }
            if(c!=0)
                _dist[c]++;
	}
    }
    
protected:
    int _dist[10];
};


int main(int argc, char** argv)
{    
    // we have one observer
    CGConj no;        
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
    
    desc.add_options()
        ("thres", boost::program_options::value<double>()->default_value(0.7), "conjugation threshold");
 
    // let cg_engine add some program options
    cg_engine.AddProgramOptions(desc);
    
    // now read in the command line
    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);    
        po::notify(vm);
        thres = vm["thres"].as<double>();
    }
    catch(po::error err) {
        cout << "error parsing command line: " << err.what() << endl;
        return -1;
    }
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

