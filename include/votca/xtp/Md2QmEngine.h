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

#ifndef VOTCA_XTP_MD2QMENGINE_H
#define VOTCA_XTP_MD2QMENGINE_H

#include <votca/xtp/calculatorfactory.h>
#include <votca/tools/property.h>
#include <votca/xtp/statesaversqlite.h>
#include <votca/xtp/topology.h>
#include <votca/csg/topology.h>

namespace XTP = votca::xtp;
namespace CSG = votca::csg;
namespace TOOLS = votca::tools;

namespace votca{
namespace xtp{


class Md2QmEngine
{
public:

    Md2QmEngine() { };
   ~Md2QmEngine();

    void Initialize(const std::string &xmlfile);
    void PrintInfo();

    // Converts atomistic to QM topology
    void            Md2Qm(csg::Topology *mdtop, Topology *qmtop);
    // Creates an QM molecule container based on MD molecule and the xml std::map 
    Molecule  *MoleculeFactory(csg::Molecule *molMDTemplate);
    // Partitions the QM molecule on segments and fragments
    Molecule  *ExportMolecule(Molecule *molQM, Topology *qmtop);

private:

    tools::Property _typology;
    
    // Type std::vectors
    std::vector < Molecule* >    _molecule_types;
    std::vector < Segment* >     _segment_types;
    std::vector < SegmentType* > _qmUnits;
    std::vector < Fragment* >    _fragment_types;
    std::vector < Atom* >        _atom_types;    

    // MD <-> QM Maps
    std::map < std::string, std::map < int, std::map < std::string, Atom* > > >
                                            _map_mol_resNr_atm_atmType;
    std::map < std::string, Molecule* >          _map_MoleculeName_MoleculeType;
    std::map < std::string, std::string >                  _map_MoleculeMDName_MoleculeName;
    std::map < int, Segment* >              _map_id_segment;


    // Type Creators
    Molecule    *AddMoleculeType(int molecule_id, tools::Property *property);
    Segment     *AddSegmentType(int segment_id, tools::Property *property);
    SegmentType *AddQMUnit(int unit_id, tools::Property *property);
    Fragment    *AddFragmentType(int fragment_id, tools::Property *property);
    Atom        *AddAtomType(Molecule *owner,
                                  std::string residue_name,    int residue_number,
                                  std::string md_atom_name,    int md_atom_id,
                                  bool hasQMPart,         int qm_atom_id,
                                  tools::vec qmpos,              std::string element,
                                  double weight);
    
    

    const std::string   &getMoleculeName(const std::string &mdname);
    Molecule  *getMoleculeType(const std::string &name);
    Atom      *getAtomType(const std::string &molMdName,
                                int resNr, const std::string &mdAtomName);
    void            ReadXYZFile(std::string &file,
                                 std::map<int, std::pair<std::string,tools::vec> > &intCoords);


};

}
}


#endif // VOTCA_XTP_MD2QMENGINE_H
