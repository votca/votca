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

using namespace std;

void help_text()
{
    cout << "csg_map \n\n";                
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


int main(int argc, char** argv)
{    
    // we have one observer
    CGTrjWrite writer;        
    // The CGEngine does the work
    CGEngine cg_engine;
    
    // file to write cg trajectory to
    string out;
    
    // add our observer that it gets called to analyze frames
    cg_engine.AddObserver((CGObserver*)&writer);


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
      ("out", po::value<string>(&out), "output file for coarse-grained trajectory");
    
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

    if (!vm.count("out")) {
        cerr << "need to specify output trajectory\n";
        return -1;
    }
    writer.setOut(out);

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

