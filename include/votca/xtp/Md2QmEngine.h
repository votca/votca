/* 
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef _MD2QMENGINE_H
#define _MD2QMENGINE_H

#include <votca/xtp/calculatorfactory.h>
#include <votca/tools/property.h>
#include <votca/xtp/statesaversqlite.h>
#include <votca/xtp/topology.h>
#include <votca/csg/topology.h>

namespace XTP = votca::xtp;
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
    void            Md2Qm(CSG::Topology *mdtop, XTP::Topology *qmtop);
    // Creates an QM molecule container based on MD molecule and the xml map 
    XTP::Molecule  *MoleculeFactory(CSG::Molecule *molMDTemplate);
    // Partitions the QM molecule on segments and fragments
    XTP::Molecule  *ExportMolecule(XTP::Molecule *molQM, XTP::Topology *qmtop);

private:

    Property _typology;
    
    // Type vectors
    vector < XTP::Molecule* >    _molecule_types;
    vector < XTP::Segment* >     _segment_types;
    vector < XTP::SegmentType* > _qmUnits;
    vector < XTP::Fragment* >    _fragment_types;
    vector < XTP::Atom* >        _atom_types;    

    // MD <-> QM Maps
    map < string, map < int, map < string, XTP::Atom* > > >
                                            _map_mol_resNr_atm_atmType;
    map < string, XTP::Molecule* >          _map_MoleculeName_MoleculeType;
    map < string, string >                  _map_MoleculeMDName_MoleculeName;
    map < int, XTP::Segment* >              _map_id_segment;


    // Type Creators
    XTP::Molecule    *AddMoleculeType(int molecule_id, Property *property);
    XTP::Segment     *AddSegmentType(int segment_id, Property *property);
    XTP::SegmentType *AddQMUnit(int unit_id, Property *property);
    XTP::Fragment    *AddFragmentType(int fragment_id, Property *property);
    XTP::Atom        *AddAtomType(XTP::Molecule *owner,
                                  string residue_name,    int residue_number,
                                  string md_atom_name,    int md_atom_id,
                                  bool hasQMPart,         int qm_atom_id,
                                  vec qmpos,              string element,
                                  double weight);
    
    

    const string   &getMoleculeName(const string &mdname);
    XTP::Molecule  *getMoleculeType(const string &name);
    XTP::Atom      *getAtomType(const string &molMdName,
                                int resNr, const string &mdAtomName);
    void            ReadXYZFile(string &file,
                                 map<int, pair<string,vec> > &intCoords);


};

#endif /* _MD2QMENGINE_H */
