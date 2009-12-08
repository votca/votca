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
#include "version.h"
#include "nematicorder.h"

using namespace std;

void help_text()
{
    votca::csg::HelpTextHeader("csg_nemat");
    cout << "!! EXPERIMENTAL !! Calculate nematic order parameter.\n"
            " Needs non-spherical beads in mapping.\n\n";           
}

class CGNematicOrder
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {
        if(_filename != "") {
            _file.open(_filename.c_str());
            if(!_file)
                throw runtime_error("cannot open " + _filename + " for output");
        }
        _u.zero();
        _v.zero();
        _w.zero();
        _n = 0;
    }
    
    void EndCG() {
        _file.close();
        _u*=1./(double)_n;
        _v*=1./(double)_n;
        _w*=1./(double)_n;
        if(_bU) {
            cout << "nematic order of u: "<< endl        
            << _u.eigenvalues[0] << ": " << _u.eigenvecs[0].normalize() << "\n"
            << _u.eigenvalues[1] << ": " << _u.eigenvecs[1].normalize() << "\n"
            << _u.eigenvalues[2] << ": " << _u.eigenvecs[2].normalize() << "\n";
        }
        
        if(_bV) {
            cout << "nematic order of v: "<< endl        
            << _v.eigenvalues[0] << ": " << _v.eigenvecs[0] << "\n"
            << _v.eigenvalues[1] << ": " << _v.eigenvecs[1] << "\n"
            << _v.eigenvalues[2] << ": " << _v.eigenvecs[2] << "\n";
        }
        if(_bW) {
            cout << "nematic order of w: "<< endl        
            << _w.eigenvalues[0] << ": " << _w.eigenvecs[0] << "\n"
            << _w.eigenvalues[1] << ": " << _w.eigenvecs[1] << "\n"
            << _w.eigenvalues[2] << ": " << _w.eigenvecs[2] << "\n";
        }
    };
    
    void EvalConfiguration(Topology *conf, Topology*conf_atom = 0) {
        _bU = true; //conf->HasU();
        _bV = true; //conf->HasV();
        _bW = true; //conf->HasW();

        nemat.Process(*conf, _filter);
        if(_file) {
            _file << conf->getTime() << " ";        
            if(_bU)        
                _file << nemat.NematicU().eigenvalues[2] << " ";
            if(_bV)        
                _file << nemat.NematicV().eigenvalues[2] << " ";
            if(_bW)        
                _file << nemat.NematicW().eigenvalues[2] << " ";
            _file << endl;
        } 
        
        _u += nemat.NematicU();
        _v += nemat.NematicV();
        _w += nemat.NematicW();
        _n++;        
    }
    
    void setOut(string filename) { _filename = filename; }
    
    void setFilter(const string &filter) { _filter = filter; }
protected:
    NematicOrder nemat;
    bool _bU, _bV, _bW;
    ofstream _file;
    string _filename;
    int _n;
    matrix::eigensystem_t _u, _v, _w;
    string _filter;
};


int main(int argc, char** argv)
{    
    // we have one observer, this analyzes neamtic order
    CGNematicOrder no;        
    // The CGEngine does the work
    CGEngine cg_engine;
    string filter;
    
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
    
    desc.add_options()
        ("filter", boost::program_options::value<string>(&filter)->default_value("*"), "filter molecule names");
        ("out", boost::program_options::value<string>(), "output nematic order prameter into file");
    
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
    // or asks for the program version?
    if (vm.count("version")) {
        cout << "csg_nemat, lib version " << LIB_VERSION_STR  << "\n";                        
        return 0;
    }
    if(vm.count("out"))
        no.setOut(vm["out"].as<string>());    
    no.setFilter(filter);
    // try to run the cg process, go through the frames, etc...
    try {
        cg_engine.Run(desc, vm);
    }
    // did an error occour?
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
    }
    return 0;
}

