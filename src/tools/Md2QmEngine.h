#ifndef _MD2QMENGINE_H
#define _MD2QMENGINE_H

#include <votca/ctp/qmtopology.h>
#include <votca/tools/property.h>
#include <votca/ctp/statesaversqlite.h>
#include <votca/moo/units.h>
// #include <votca/ctp/atom.h>
// #include <votca/ctp/fragment.h>
// #include <votca/ctp/segment.h>
// #include <votca/ctp/segmenttype.h>

using namespace votca::ctp;
using namespace votca::tools;


namespace CTP = votca::ctp;


class Md2QmEngine
{
public:

    Md2QmEngine() { };
   ~Md2QmEngine() { };

    /**
     * Load properties from --cg xmlfile
     * @param xmlfile
     */
    void Initialize(const string &xmlfile);

    /**
     * Map mdtop to qmtop using typology
     * @param mdtop
     * @param qmtop
     */
    void Md2Qm(Topology *mdtop, Topology *qmtop);

private:

    Property _typology;
    /*
    // Type vectors
    vector < Molecule* > _molecule_types;
    vector < Segment* >  _segment_types;
    vector < Fragment* > _fragment_types;
    vector < Atom* >     _atom_types;

    // MD <-> QM Maps
    map < string, Molecule* > _map_MoleculeName_MoleculeType;
    map < string, string >    _map_MoleculeMDName_MoleculeName;
    map < int, Segment* >     _map_id_segment;

    // Type Creators
    Molecule  *AddMoleculeType(int molecule_id, Property *property);
    Segment   *AddSegmentType(int segment_id, Property *property);
    Fragment  *AddFragmentType(int fragment_id, Property *property);
    Atom      *AddAtomType(Molecule *owner, int atom_id, string atom_name,
                           int residue_number, double weight);

    Molecule *getMoleculeType(string name);
    */
};











#endif /* _MD2QMENGINE_H */
