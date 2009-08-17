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
         for(int i=0;i<10;++i)
             for(int j=0;j<1000;++j)
             _lifetime[i][j]=0;
    };
    void EndCG() {
        for(int i=0;i<10;++i) {
           cout << i+1 << " " << _dist[i] << endl;
        }
        BeginCycle();EndCycle();
        for(int j=0;j<1000;++j) {
            cout << j << " ";
            for(int i=0;i<10;++i)
                cout << _lifetime[i][j] << " ";
            cout << endl;
        }
             
    };
    
    void EvalConfiguration(Topology *conf, Topology *conf_atom = 0) {
        Topology *top = conf;
        BeginCycle();
        for(MoleculeContainer::iterator iter=top->Molecules().begin();
            iter!=top->Molecules().end(); ++iter) {
            Molecule *mi = *iter;
            if(mi->BeadCount() <= 1) continue;
            vector<int> beads;
            vec no = mi->getBead(0)->getU();
            beads.push_back(mi->getBead(0)->getId());

            for(int i=1; i<mi->BeadCount(); ++i) {
                vec n = mi->getBead(i)->getU();
                if(fabs(no*n) < thres) {
                    UpdateCrgUnit(beads);
                    _dist[beads.size()-1]++;
                    beads.clear();
                }
                beads.push_back(mi->getBead(i)->getId());
                no = n; 
            }
           _dist[beads.size()-1]++;
            UpdateCrgUnit(beads);
	    }
        EndCycle();
    }
    
protected:
    int _dist[10];
    int _lifetime[10][1000];
    
    struct crgunit_t {
        vector<int> _beads;
        crgunit_t(vector<int> beads)
            : _beads(beads), _alive(true), _lifetime(0) {}
        crgunit_t() {}
        
        bool operator==(const vector<int> &c2) {
            if(_beads.size()!=c2.size()) return false;
            for(int i=0; i<_beads.size(); ++i)
                if(_beads[i] != c2[i]) return false;
            return true;
        } 
        
        int _lifetime;
        bool _alive;
    };
    
    list<crgunit_t> _crgunits;

    void UpdateCrgUnit(vector<int> beads) {
        list<crgunit_t>::iterator iter
            = find(_crgunits.begin(), _crgunits.end(), beads);
        if(iter == _crgunits.end())
            _crgunits.push_back(crgunit_t(beads));
        else {
            (*iter)._alive = true;
            (*iter)._lifetime++;
        }                
    }
    
    void BeginCycle() {
        for(list<crgunit_t>::iterator iter = _crgunits.begin();
            iter != _crgunits.end(); ++iter)
            (*iter)._alive = false;
    }
    
    void EndCycle() {
        for(list<crgunit_t>::iterator iter = _crgunits.begin();
            iter != _crgunits.end();) {           
            if(!(*iter)._alive) {
                if((*iter)._lifetime < 1000)
                    _lifetime[(*iter)._beads.size()-1][(*iter)._lifetime]+=(*iter)._lifetime;
                else
                    _lifetime[(*iter)._beads.size()-1][999]+=(*iter)._lifetime;
                iter =  _crgunits.erase(iter);
            } else ++iter;
        }
    }
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

