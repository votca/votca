#ifndef VOTCA_XTP_EQUILIBRATE_H
#define VOTCA_XTP_EQUILIBRATE_H

#include <votca/xtp/qmtool.h>
#include <votca/tools/tokenizer.h>
#include <math.h>

#include <votca/csg/topology.h>
#include <votca/csg/trajectorywriter.h>
#include <votca/csg/trajectoryreader.h>
#include <votca/csg/topologyreader.h>


namespace votca { namespace xtp {
    
namespace CSG = votca::csg;


class OrientedMolecule
{
public:
    OrientedMolecule() {};
   ~OrientedMolecule() {};
private:
        
};




class Equilibrate : public QMTool
{
public:

    Equilibrate() { };
   ~Equilibrate() { };

    string Identify() { return "equilibrate"; }

    void Initialize(Property *opt);
    bool Evaluate();

private:
    
    string _topfile;
    string _trjfile;
    
    CSG::Topology _mdtopol;
    CSG::TopologyReader *_topread;
    CSG::TrajectoryReader *_trjread;

};



void Equilibrate::Initialize(Property *opt) {
    
    string key = "options.equilibrate";
    _topfile = opt->get(key+".topfile").as<string>();
    _trjfile = opt->get(key+".trjfile").as<string>();
    
    CSG::TrajectoryWriter::RegisterPlugins();
    CSG::TrajectoryReader::RegisterPlugins();
    CSG::TopologyReader::RegisterPlugins();
    
    cout << endl;
    
    // READ TOPOLOGY
    _topread = CSG::TopReaderFactory().Create(_topfile);

    if (_topread == NULL) {
        throw runtime_error( string("Input format not supported: ")
                           + _topfile );
    }

    _topread->ReadTopology(_topfile, this->_mdtopol);
    if (votca::tools::globals::verbose) {
        cout << "Read MD topology from " << _topfile << ": Found "
             << _mdtopol.BeadCount() << " atoms in "
             << _mdtopol.MoleculeCount() << " molecules. "
             << endl;
    }
    
    // READ CONFIGURATION
    _trjread = CSG::TrjReaderFactory().Create(_trjfile);

    if (_trjread == NULL) {
        throw runtime_error( string("Input format not supported: ")
                           + _trjfile );
    }
    _trjread->Open(_trjfile);
    _trjread->FirstFrame(this->_mdtopol);
    
    return;
}

bool Equilibrate::Evaluate() {
    
    cout << endl << "Equilibrate says 'Hello!'" << flush;
    
    for (int i=0; i<_mdtopol.MoleculeCount(); ++i) {
        
        CSG::Molecule *mol = _mdtopol.getMolecule(i);
        cout << endl << "Molecule " << mol->getName() << " " << mol->getId() << flush;
        
        for (int j=0; j<mol->BeadCount(); ++j) {
            
            CSG::Bead *atm = mol->getBead(j);
            cout << endl << "    Atom " << atm->getName() << " " << atm->getId() << flush;
            
        }
        
    }
    
    return true;
}

}}

#endif // VOTCA_XTP_EQUILIBRATE_H