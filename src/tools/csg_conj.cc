/* 
 * Copyright 2009 The VOTCA Development Team (http://www.votca.org)
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
#include <math.h>
#include "version.h"

#define TIME_ARRAY_SIZE 500
#define TIME_GROUP  5
#define TIME_RANGE  TIME_GROUP*TIME_ARRAY_SIZE

#define TIME_MULTIPLIER 2

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
             for(int j=0;j<TIME_ARRAY_SIZE;++j)
             _lifetime[i][j]=0;
    };
    void EndCG() {
        for(int i=0;i<10;++i) {
           cout << i+1 << " " << _dist[i] << endl;
        }
        BeginCycle();EndCycle();
        for(int j=0;j<TIME_ARRAY_SIZE;++j) {
            cout << j*TIME_MULTIPLIER*TIME_GROUP << " ";
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
    int _lifetime[10][TIME_ARRAY_SIZE];
    
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
                if((*iter)._lifetime < TIME_RANGE)
                    _lifetime[(*iter)._beads.size()-1][(*iter)._lifetime/TIME_GROUP]+=(*iter)._lifetime;
                else
                    _lifetime[(*iter)._beads.size()-1][TIME_ARRAY_SIZE-1]+=(*iter)._lifetime;
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
    
    namespace po = boost::program_options;
    try {
        cg_engine.Initialize();
        // add our observer that it gets called to analyze frames
        cg_engine.AddObserver((CGObserver*)&no);
    
        // lets read in some program options
        

        cg_engine.AddProgramOptions()
            ("thres", boost::program_options::value<double>(&thres)->default_value(0.7), "conjugation threshold");

        cg_engine.ParseCommandLine(argc, argv); 
    // let cg_engine add some program options
        po::variables_map &vm
            = cg_engine.OptionsMap();
        po::options_description &desc
            = cg_engine.OptionsDesc();
 
        // does the user want help?
        if (vm.count("help")) {
            votca::csg::HelpTextHeader("csg_nemat");
            cout << desc << endl;
            return 0;
        }
        cout << "using thresold: " << thres << endl;
        cg_engine.Run();
    }
    // did an error occour?
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
    }
    return 0;
}

