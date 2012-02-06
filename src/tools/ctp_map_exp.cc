#include <iostream>
#include <fstream>
#include <stdexcept>

#include "votca/tools/application.h"
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/topologyreader.h>
#include <votca/ctp/statesaversqlite2.h>
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
    void Save(string mode);

    void BeginEvaluate() { ; }
    bool DoTrajectory() { return 0; }
    bool DoMapping() { return 0; }


protected:
    Property               _options;
    CSG::Topology          _mdtopol;
    CTP::Topology          _qmtopol;

    Md2QmEngine            _md2qm;
    CTP::StateSaverSQLite2 _statsav;
    string                 _outdb;

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
    AddProgramOptions() ("check", propt::value<string> (),
                         "  Extra info + PDB output ");
}

bool CtpMapExp::EvaluateOptions() {

    CheckRequired("top", "Missing topology file");
    CheckRequired("cg", "Missing CG definition file");
    CheckRequired("file");

    return 1;
}

void CtpMapExp::Run() {

    // +++++++++++++++++++++++++++++++++++++ //
    // Initialize MD2QM Engine and SQLite Db //
    // +++++++++++++++++++++++++++++++++++++ //

    _outdb = _op_vm["file"].as<string> ();
    _statsav.Open(_qmtopol, _outdb);
    _statsav.FramesInDatabase();
    _statsav.Close();

    string cgfile = _op_vm["cg"].as<string> ();
    _md2qm.Initialize(cgfile);


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
    //cout << "MD Topology from " << topfile << ": Found "
    //     << _mdtopol.BeadCount() << " atoms in "
    //     << _mdtopol.MoleculeCount() << " molecules. "
    //     << endl;



    // ++++++++++++++++++++++++++++++ //
    // Create MD trajectory from file //
    // ++++++++++++++++++++++++++++++ //

    // Create trajectory reader and initialize
    string trjfile =  _op_vm["trj"].as<string> ();
    CSG::TrajectoryReader *trjread;
    trjread = CSG::TrjReaderFactory().Create(trjfile);

    if (trjread == NULL) {
        throw runtime_error( string("Input format not supported: ")
                           + _op_vm["trj"].as<string> () );
    }
    trjread->Open(trjfile);
    trjread->FirstFrame(this->_mdtopol);

    int    firstFrame = 1;
    int    nFrames    = 1;
    bool   beginAt    = 0;
    double startTime  = _mdtopol.getTime();

    if (_op_vm.count("nframes")) {
        nFrames = _op_vm["nframes"].as<int> ();
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
    
    // +++++++++++++++++++++++++ //
    // Convert MD to QM Topology //
    // +++++++++++++++++++++++++ //

    for (int saved = 0; hasFrame && saved < nFrames;
         hasFrame = trjread->NextFrame(this->_mdtopol), saved++) {

        _md2qm.Md2Qm(&_mdtopol, &_qmtopol);

        //if (_op_vm.count("check")) {
        //    string pdbfile = _op_vm["check"].as<string> ();
        //    _md2qm.CheckProduct(&_qmtopol, "md_" + pdbfile);
        //}

    // +++++++++++++++++++++++++ //
    // Save to SQLite State File //
    // +++++++++++++++++++++++++ //

        this->Save("");
    }

    // trjread->Close();
    // delete trjread;

}

void CtpMapExp::Save(string mode) {    
    
    _statsav.Open(_qmtopol, _outdb);

    _statsav.WriteFrame();

    if (_op_vm.count("check")) {
        CTP::Topology *TopSQL = NULL;
        TopSQL = _statsav.getTopology();
        cout << endl << "Checking topology read from SQL file." << endl;
        string pdbfile = _op_vm["check"].as<string> ();
        _md2qm.CheckProduct(TopSQL, pdbfile);
    }

    _statsav.Close();

}


int main(int argc, char** argv)
{
    CtpMapExp ctpmap;
    return ctpmap.Exec(argc, argv);
}
