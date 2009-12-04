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
#include "version.h"

using namespace std;

void help_text()
{
    votca::csg::HelpTextHeader("csg_map");
    cout << "Map a reference trajectory to a coarse-grained trajectory\n"
            "This program can be used to map a whole trajectory or only\n"
            "create an initial configuration for a coarse-grained run.\n\n";     
}

class CGTrjWrite
    : public CGObserver
{
public:
    void BeginCG(Topology *top, Topology *top_atom) {
        cout << "writing coarse-grained trajectory to " << _out << endl;

        _writer = TrjWriterFactory().Create(_out);
        if(_writer == NULL)
            throw runtime_error("output format not supported: " + _out);
        
        _writer->Open(_out);
    };

    void EndCG() {
        _writer->Close();
        delete _writer;
    };
    
    void EvalConfiguration(Topology *top, Topology *top_atom = 0) {
        _writer->Write(top);
    }

    void setOut(const string &out) { _out = out; }
    
protected:
    string _out;
    TrajectoryWriter *_writer;
};

// lets read in some program options
namespace po = boost::program_options;

int main(int argc, char** argv)
{    
    // we have one observer
    CGTrjWrite writer;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // file to write cg trajectory to
    string out;
    
    try {
        // add our observer that it gets called to analyze frames
        cg_engine.Initialize();

        cg_engine.AddObserver((CGObserver*)&writer);


        cg_engine.AddProgramOptions()
            ("out", po::value<string>(&out), "  output file for coarse-grained trajectory");

        cg_engine.ParseCommandLine(argc, argv);
        
        po::variables_map &vm
            = cg_engine.OptionsMap();
        po::options_description &desc
            = cg_engine.OptionsDesc();

        // does the user want help?
        if (vm.count("help")) {
            help_text();
            cout << desc << endl;
            return 0;
        }

        if (!vm.count("out")) {
            cerr << "need to specify output trajectory\n";
            return -1;
        }
        writer.setOut(out);

        cg_engine.Run();
    }
    // did an error occour?
    catch(std::exception &error) {
        cerr << "An error occoured!" << endl << error.what() << endl;
    }
    return 0;
}

