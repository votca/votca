#include <iostream>
#include <fstream>
#include <stdexcept>

#include <votca/ctp/qmtopology.h>
#include <votca/csg/csgapplication.h>
#include "votca/tools/application.h"
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/topologyreader.h>
#include "Md2QmEngine.h"

using namespace std;

namespace CSG = votca::csg;
namespace CTP = votca::ctp;


class CtpMapExp : public Application
{

public:
    string ProgramName()  { return "ctp_map_exp"; }
    void   HelpText(ostream &out) {out << "Generates QM|MD topology" << endl;}

    void Initialize();
    bool EvaluateOptions();
    void Run();

    void BeginEvaluate() { ; }
    bool DoTrajectory() { return 0; }
    bool DoMapping() { return 0; }


protected:
    Property            _options;
    CSG::Topology       _mdtopol;
    CTP::Topology       _qmtopol;

    Md2QmEngine         _md2qm;

};

namespace propt = boost::program_options;

void CtpMapExp::Initialize() {

    CSG::TrajectoryWriter::RegisterPlugins();
    CSG::TrajectoryReader::RegisterPlugins();
    CSG::TopologyReader::RegisterPlugins();

    AddProgramOptions() ("top", propt::value<string> (),
                         "  Atomistic topology file ");
    AddProgramOptions() ("trj", propt::value<string> (),
                         "  Atomistic trajetory file ");
    AddProgramOptions() ("cg",  propt::value<string> (),
                         "  Coarse-Graining definitions ");
    AddProgramOptions() ("file,f", propt::value<string> (),
                         "  SQLite state file ");
}

bool CtpMapExp::EvaluateOptions() {

    CheckRequired("top", "Missing topology file");
    CheckRequired("cg", "Missing CG definition file");
    CheckRequired("file");

    return 1;
}

void CtpMapExp::Run() {

    // ++++++++++++++++++++++++++++ //
    // Create MD topology from file //
    // ++++++++++++++++++++++++++++ //

    // Create topology reader
    string topfile = _op_vm["top"].as<string> ();
    CSG::TopologyReader *topread;
    topread = CSG::TopReaderFactory().Create(topfile);

    if (topread == NULL) {
        throw runtime_error( string("Input format not supported: ")
                           + _op_vm["top"].as<string> () );
    }

    topread->ReadTopology(topfile, this->_mdtopol);
    cout << "MD Topology from " << topfile << ": Found "
         << _mdtopol.BeadCount() << " atoms in "
         << _mdtopol.MoleculeCount() << " molecules. "
         << endl;

    // ++++++++++++++++++++++++++++++ //
    // Create MD trajectory from file //
    // ++++++++++++++++++++++++++++++ //

    // Create trajectory reader and initialize
    string trjfile =  _op_vm["trj"].as<string> ();
    CSG::TrajectoryReader *trjread;
    trjread = CSG::TrjReaderFactory().Create(trjfile);

    if (trjread == NULL) {
        throw runtime_error( string("Input format not supproated: ")
                           + _op_vm["trj"].as<string> () );
    }
    trjread->Open(trjfile);
    trjread->FirstFrame(this->_mdtopol);

    int    firstFrame = 1;
    int    frameNo    = 0;
    bool   beginAt    = 0;
    double startTime  = _mdtopol.getTime();

    if (_op_vm.count("nframes")) {
        frameNo = _op_vm["nframes"].as<int> ();
    }
    if (_op_vm.count("first-frame")) {
        firstFrame = _op_vm["first-frame"].as<int> ();
    }    
    if (_op_vm.count("begin")) {
        beginAt = true;
        startTime = _op_vm["begin"].as<double> ();
    }

    // Extract first frame specified
    bool hasFrame;

    for (hasFrame = true; hasFrame == true;
         hasFrame = trjread->NextFrame(this->_mdtopol)) {
         if (  ((_mdtopol.getTime() < startTime) && beginAt )
               || firstFrame > 1 ) {
             firstFrame--;
             continue;
         }
         break;
    }
    if ( ! hasFrame) {
        trjread->Close();
        delete trjread;

        throw runtime_error("Time or frame number exceeds trajectory length");
    }
    /*
    CSG::MoleculeContainer::iterator mit;
    for (mit = _mdtopol.Molecules().begin();
         mit != _mdtopol.Molecules().end();
         mit++ ) {

        CSG::Molecule *mol = *mit;
        for (int i = 0; i < mol->BeadCount(); i++) {
            CSG::Bead* btm = mol->getBead(i);
            cout << "ATOM getId  " << btm->getId();
            cout << " | getName  " << btm->getName();
            cout << " | getMol   " << btm->getMolecule()->getId();
            cout << " | MolName  " << btm->getMolecule()->getName();
            cout << " | getResnr " << btm->getResnr();
            cout << " | getType  " << btm->getType()->getName();
            cout << " | getPos   " << btm->getPos();
            cout << endl;
        }    
    }  */
    
    // ++++++++++++++++++++++++++ //
    // Convert MD- to QM Topology //
    // ++++++++++++++++++++++++++ //
    
    string cgfile = _op_vm["cg"].as<string> ();
    _md2qm.Initialize(cgfile);
    _md2qm.PrintInfo();
    _md2qm.Md2Qm(&_mdtopol, &_qmtopol);

    

    
}


int main(int argc, char** argv)
{
    CtpMapExp ctpmap;
    return ctpmap.Exec(argc, argv);
}
