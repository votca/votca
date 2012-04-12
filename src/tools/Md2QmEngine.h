#ifndef _MD2QMENGINE_H
#define _MD2QMENGINE_H

#include <votca/ctp/calculatorfactory.h>
#include <votca/tools/property.h>
#include <votca/ctp/statesaversqlite.h>
#include <votca/moo/units.h>
#include <votca/ctp/topology.h>
#include <votca/csg/topology.h>


namespace CTP = votca::ctp;
namespace CSG = votca::csg;
namespace TOOLS = votca::tools;

class Md2QmEngine
{
public:

    Md2QmEngine() { };
   ~Md2QmEngine();

    void Initialize(const string &xmlfile);
    void PrintInfo();

    // Converts atomistic to QM topology
    void            Md2Qm(CSG::Topology *mdtop, CTP::Topology *qmtop);
    // Creates an QM molecule container based on MD molecule and the xml map 
    CTP::Molecule  *MoleculeFactory(CSG::Molecule *molMDTemplate);
    // Partitions the QM molecule on segments and fragments
    CTP::Molecule  *ExportMolecule(CTP::Molecule *molQM, CTP::Topology *qmtop);
    // outputs the mapping, and other useful info
    void CheckProduct(CTP::Topology *outtop, const string &pdbfile);


private:

    Property _typology;
    
    // Type vectors
    vector < CTP::Molecule* >    _molecule_types;
    vector < CTP::Segment* >     _segment_types;
    vector < CTP::SegmentType* > _qmUnits;
    vector < CTP::Fragment* >    _fragment_types;
    vector < CTP::Atom* >        _atom_types;    

    // MD <-> QM Maps
    map < string, map < int, map < string, CTP::Atom* > > >
                                            _map_mol_resNr_atm_atmType;
    map < string, CTP::Molecule* >          _map_MoleculeName_MoleculeType;
    map < string, string >                  _map_MoleculeMDName_MoleculeName;
    map < int, CTP::Segment* >              _map_id_segment;


    // Type Creators
    CTP::Molecule    *AddMoleculeType(int molecule_id, Property *property);
    CTP::Segment     *AddSegmentType(int segment_id, Property *property);
    CTP::SegmentType *AddQMUnit(int unit_id, Property *property);
    CTP::Fragment    *AddFragmentType(int fragment_id, Property *property);
    CTP::Atom        *AddAtomType(CTP::Molecule *owner,
                                  string residue_name,    int residue_number,
                                  string md_atom_name,    int md_atom_id,
                                  bool hasQMPart,         int qm_atom_id,
                                  vec qmpos,              string element,
                                  double weight);
    
    

    const string   &getMoleculeName(const string &mdname);
    CTP::Molecule  *getMoleculeType(const string &name);
    CTP::Atom      *getAtomType(const string &molMdName,
                                int resNr, const string &mdAtomName);
    void            getIntCoords(string &file,
                                 map<int, pair<string,vec> > &intCoords);


};

#endif /* _MD2QMENGINE_H */
