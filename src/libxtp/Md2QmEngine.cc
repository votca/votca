/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/csg/boundarycondition.h>
#include <votca/tools/globals.h>
#include <votca/xtp/Md2QmEngine.h>

#include <votca/xtp/atom.h>
#include <votca/xtp/fragment.h>
#include <votca/xtp/molecule.h>
#include <votca/xtp/segment.h>
#include <votca/xtp/segmenttype.h>
#include <votca/xtp/topology.h>

namespace votca {
namespace xtp {

/**
 * Clears all engine template ('type') containers.
 */
Md2QmEngine::~Md2QmEngine() {
  // Clean up list of molecule types
  for (Molecule *mol : _molecule_types) {
    delete mol;
  }
  _molecule_types.clear();

  // clean up maps
  _map_MoleculeName_MoleculeType.clear();
  _map_MoleculeMDName_MoleculeName.clear();
}

/**
 * Reads in mapping .xml file to generate molecule, segment, fragment and
 * atom templates ('types') that contain all necessary information to
 * start mapping, except coordinates.
 * @param xmlfile
 */
void Md2QmEngine::Initialize(const std::string &xmlfile) {

  tools::Property typology;
  load_property_from_xml(typology, xmlfile.c_str());

  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
  // XML to Types: Molecules => Segments => Fragments => Atoms //
  // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

  std::string key = "topology.molecules.molecule";
  std::list<tools::Property *> molecules = typology.Select(key);
  int molecule_id = 1;
  int qmunit_id = 1;  // Counts segment types

  for (tools::Property *mol : molecules) {

    Molecule *molecule = AddMoleculeType(molecule_id++, mol);
    std::string molMdName = mol->get("mdname").as<std::string>();

    // +++++++++++++ //
    // Load segments //
    // +++++++++++++ //

    key = "segments.segment";
    std::list<tools::Property *> segments = mol->Select(key);

    int segment_id = 1;
    int md_atom_id = 1;  // <- atom id count with respect to molecule

    for (tools::Property *segprop : segments) {

      // Create new segment + associated type (QM Unit)
      Segment *segment = AddSegmentType(segment_id++, segprop);
      SegmentType *qmUnit = AddQMUnit(qmunit_id, segprop);
      ++qmunit_id;

      segment->setType(qmUnit);
      molecule->AddSegment(segment);

      // Load internal (i.e. QM-) coord.s and MOO-related properties
      std::string qmcoordsFile = "";
      std::map<int, std::pair<std::string, tools::vec> > intCoords;
      if (segprop->exists("qmcoords")) {
        qmcoordsFile = segprop->get("qmcoords").as<std::string>();
        //  QM ID    Element   Position
        this->ReadXYZFile(qmcoordsFile, intCoords);
      }

      // ++++++++++++++ //
      // Load fragments //
      // ++++++++++++++ //

      std::map<std::string, bool> fragname_isTaken;

      key = "fragments.fragment";
      std::list<tools::Property *> fragments = segprop->Select(key);
      int fragment_id = 1;

      for (tools::Property *fragmentprop : fragments) {

        // Create new fragment
        Fragment *fragment = AddFragmentType(fragment_id++, fragmentprop);
        segment->AddFragment(fragment);

        // Verify that this fragment name is not taken already
        try {
          // This should throw
          bool taken = fragname_isTaken.at(fragment->getName());
          std::cout << "ERROR Fragment name '" << fragment->getName()
                    << "' in segment '" << segment->getName()
                    << "' occurs more than once." << std::endl;
          if (taken) throw std::runtime_error("(see above, naming collision)");
        } catch (const std::exception &out_of_range) {
          fragname_isTaken[fragment->getName()] = true;
        }

        // Load local-frame definition
        std::vector<int> trihedron =
            fragmentprop->get("localframe").as<std::vector<int> >();
        while (trihedron.size() < 3) {
          trihedron.push_back(-1);
        }
        fragment->setTrihedron(trihedron);

        // ++++++++++ //
        // Load atoms //
        // ++++++++++ //

        std::string mdatoms = fragmentprop->get("mdatoms").as<std::string>();
        std::string qmatoms = fragmentprop->get("qmatoms").as<std::string>();
        std::string weights = fragmentprop->get("weights").as<std::string>();

        tools::Tokenizer tok_md_atoms(mdatoms, " \t\n");
        tools::Tokenizer tok_qm_atoms(qmatoms, " \t\n");
        tools::Tokenizer tok_weights(weights, " \t\n");

        std::vector<std::string> md_atoms_info;
        std::vector<std::string> qm_atoms_info;
        std::vector<std::string> atom_weights;
        ;

        tok_md_atoms.ToVector(md_atoms_info);
        tok_qm_atoms.ToVector(qm_atoms_info);
        tok_weights.ToVector(atom_weights);

        if ((md_atoms_info.size() != qm_atoms_info.size()) ||
            (md_atoms_info.size() != atom_weights.size())) {
          std::cout << "ERROR: "
                    << "Could not allocate MD atoms to QM atoms or weights"
                    << " in fragment " << fragment->getName() << " in segment "
                    << segment->getName() << " in molecule " << molMdName
                    << " due to inconsistent number of columns"
                    << " (MD: " << md_atoms_info.size() << ","
                    << " QM: " << qm_atoms_info.size() << ")"
                    << " Weights: " << atom_weights.size() << ")." << std::endl;
          std::cout << "NOTE: "
                    << "To define an MD atom without QM counterpart, insert "
                    << "a single ':' in the associated QM-atoms column and "
                    << "specify a mapping weight of 0." << std::endl;
          throw std::runtime_error("Inconsistency in mapping file.");
        }

        std::vector<std::string>::iterator it_md_atom_name;
        std::vector<std::string>::iterator it_qm_atom_name;
        std::vector<std::string>::iterator it_atom_weight;

        for (it_md_atom_name = md_atoms_info.begin(),
            it_qm_atom_name = qm_atoms_info.begin(),
            it_atom_weight = atom_weights.begin();
             it_md_atom_name != md_atoms_info.end();
             ++it_md_atom_name, ++it_qm_atom_name, ++it_atom_weight) {

          // ++++++++++++++++++ //
          // Create single atom //
          // ++++++++++++++++++ //

          tools::Tokenizer tok_md((*it_md_atom_name), ":");
          tools::Tokenizer tok_qm((*it_qm_atom_name), ":");

          std::vector<std::string> md_atom_specs;
          std::vector<std::string> qm_atom_specs;

          tok_md.ToVector(md_atom_specs);
          tok_qm.ToVector(qm_atom_specs);

          // MD atom information
          int residue_number = boost::lexical_cast<int>(md_atom_specs[0]);
          std::string residue_name = md_atom_specs[1];
          std::string md_atom_name = md_atom_specs[2];

          // QM atom information
          bool hasQMPart = false;
          int qm_atom_id = -1;
          tools::vec qmPos = tools::vec(0, 0, 0);
          // The atomic element is first taken as first character of
          // the MD name; for atoms with QM counterpart, the element
          // is read from the coordinates file
          std::string element = md_atom_name.substr(0, 1);
          // Check whether MD atom has QM counterpart
          if (qm_atom_specs.size() == 2) {
            hasQMPart = true;
            qm_atom_id = boost::lexical_cast<int>(qm_atom_specs[0]);

            // Look up qm coordinates in table created from xyz file
            try {
              qmPos = intCoords.at(qm_atom_id).second;
              element = intCoords.at(qm_atom_id).first;
              // Check whether elements of md and qm match
              if (intCoords.at(qm_atom_id).first.substr(0, 1) !=
                  md_atom_name.substr(0, 1)) {
                std::cout << "WARNING: Atom " << md_atom_name << " in mol. "
                          << molMdName << " appears to have element type "
                          << md_atom_name.substr(0, 1)
                          << ", but QM partner (ID " << qm_atom_id
                          << ") has element type "
                          << intCoords.at(qm_atom_id).first << std::endl;
              }
            } catch (const std::exception &out_of_range) {
              ;  // No QM coordinates specified
            }
          }

          // Mapping weight
          double weight = stod(*it_atom_weight);
          if (!hasQMPart && weight != 0) {
            std::cout << "ERROR: "
                      << "Atom " << md_atom_name << " in residue "
                      << residue_name << " in molecule " << molMdName
                      << " has no QM counterpart despite non-zero weight. "
                      << std::endl;
            throw std::runtime_error("Error in mapping file");
          }

          // Create atom
          Atom *atom = AddAtomType(molecule, residue_name, residue_number,
                                   md_atom_name, md_atom_id++, hasQMPart,
                                   qm_atom_id, qmPos, element, weight);

          try {
            this->_map_mol_resNr_atm_atmType.at(molMdName)
                .at(residue_number)
                .at(md_atom_name);
            // If this succeeded, atom has already been defined:
            std::cout << "ERROR: "
                      << "Ambiguous atom definition in molecule " << molMdName
                      << ": "
                      << "Atom " << md_atom_name << " in residue "
                      << residue_number << " exists more than once. "
                      << std::endl;
            throw std::runtime_error("Ambiguity in atom definition.");
          } catch (const std::exception &out_of_range) {
            this->_map_mol_resNr_atm_atmType[molMdName][residue_number]
                                            [md_atom_name] = atom;
          }

          fragment->AddAtom(atom);
          segment->AddAtom(atom);
          molecule->AddAtom(atom);

        } /* exit loop over atoms */

      } /* exit loop over fragments */

    } /* exit loop over segments */

  } /* exit loop over molecules */
}

/**
 * Using the templates for molecules, segments, ... generated in
 * ::Initialize, Md2Qm maps the MD topology to the initially empty
 * MD/QM Topology. Calls ::MoleculeFactory and ::ExportMolecule.
 * @param mdtop
 * @param qmtop
 */
void Md2QmEngine::Md2Qm(csg::Topology *mdtop, Topology *qmtop) {

  qmtop->CleanUp();

  // Create periodic box
  qmtop->setBox(mdtop->getBox());

  // Set trajectory meta data
  qmtop->setStep(mdtop->getStep());
  qmtop->setTime(mdtop->getTime());
  qmtop->setCanRigidify(true);

  // Add types (=> Segment types / QM units)
  for (SegmentType *type : _qmUnits) {

    std::string name = type->getName();
    std::string basis = type->getBasisName();
    std::string orbitals = type->getOrbitalsFile();
    std::string qmcoords = type->getQMCoordsFile();
    bool canRigidify = type->canRigidify();
    SegmentType *segType = qmtop->AddSegmentType(name);
    segType->setBasisName(basis);
    segType->setOrbitalsFile(orbitals);
    segType->setQMCoordsFile(qmcoords);
    segType->setCanRigidify(canRigidify);
    if (!canRigidify) {
      qmtop->setCanRigidify(false);
    }
  }

  // Populate topology in a trickle-down manner
  // (i.e. molecules => ... ... => atoms)
  for (csg::Molecule *molMD : mdtop->Molecules()) {

    // MD molecule + name
    std::string nameMolMD = molMD->getName();

    // Find QM counterpart
    Molecule *molQM = this->MoleculeFactory(molMD);
    std::string nameMolQM = molQM->getName();

    // Generate and export
    // Molecule *product =
    (void)this->ExportMolecule(molQM, qmtop);
  }
}

/**
 * Takes MD molecule, finds QM counterpart, fills out atom positions.
 * @param molMDTemplate
 * @return molQMTemplate
 */
Molecule *Md2QmEngine::MoleculeFactory(csg::Molecule *molMDTemplate) {

  std::string nameMolQM = this->getMoleculeName(molMDTemplate->getName());
  Molecule *molQMTemplate = this->getMoleculeType(nameMolQM);

  int resnrOffset = molMDTemplate->getBead(0)->getResnr();

  for (int i = 0; i < molMDTemplate->BeadCount(); i++) {
    csg::Bead *atomMD = molMDTemplate->getBead(i);
    try {
      Atom *counterpart = this->getAtomType(
          molMDTemplate->getName(), atomMD->getResnr() - resnrOffset + 1,
          atomMD->getName());
      counterpart->setPos(atomMD->getPos());
    } catch (const std::exception &out_of_range) {
      std::cout << "WARNING: No mapping instruction found for atom "
                << atomMD->getName() << " in residue number "
                << atomMD->getResnr() << " in molecule "
                << molMDTemplate->getName() << ". Skipping..." << std::endl;
    }
  }
  return molQMTemplate;
}

/**
 * Takes QM molecule template as reference, and step by step adds 'new'
 * copies of all segments, fragments, atoms within that molecule to
 * the target topology.
 * @param refMol
 * @param qmtop
 * @return
 */
Molecule *Md2QmEngine::ExportMolecule(Molecule *refMol, Topology *qmtop) {

  Molecule *newMol = qmtop->AddMolecule(refMol->getName());

  // Note: The topology is responsible for allocating IDs to
  //       molecules, segments, ...

  for (Segment *refSeg : refMol->Segments()) {

    Segment *newSeg = qmtop->AddSegment(refSeg->getName());
    newSeg->setType(qmtop->getSegmentType(refSeg->getType()->getId()));

    for (Fragment *refFrag : refSeg->Fragments()) {

      Fragment *newFrag = qmtop->AddFragment(refFrag->getName());
      newFrag->setTrihedron(refFrag->getTrihedron());

      for (Atom *refAtom : refFrag->Atoms()) {

        Atom *newAtom = qmtop->AddAtom(refAtom->getName());

        newAtom->setWeight(refAtom->getWeight());
        newAtom->setResnr(refAtom->getResnr());
        newAtom->setResname(refAtom->getResname());
        newAtom->setElement(refAtom->getElement());
        if (refAtom->HasQMPart()) {
          newAtom->setQMPart(refAtom->getQMId(), refAtom->getQMPos());
        }
        newAtom->setPos(refAtom->getPos());

        newFrag->AddAtom(newAtom);
        newSeg->AddAtom(newAtom);
        newMol->AddAtom(newAtom);
      } /* exit loop over template atoms */
      newFrag->calcPos();
      newSeg->AddFragment(newFrag);
      newMol->AddFragment(newFrag);
    } /* exit loop over template fragments */
    newSeg->calcPos();
    newMol->AddSegment(newSeg);
  } /* exit loop over template molecules */

  return newMol;
}

Atom *Md2QmEngine::AddAtomType(Molecule *owner, std::string residue_name,
                               int residue_number, std::string md_atom_name,
                               int md_atom_id, bool hasQMPart, int qm_atom_id,
                               tools::vec qmPos, std::string element,
                               double weight) {
  Atom *atom =
      new Atom(owner, residue_name, residue_number, md_atom_name, md_atom_id,
               hasQMPart, qm_atom_id, qmPos, element, weight);
  _atom_types.push_back(atom);
  return atom;
}

Fragment *Md2QmEngine::AddFragmentType(int fragment_id,
                                       tools::Property *property) {
  std::string fragment_name = property->get("name").as<std::string>();
  Fragment *fragment = new Fragment(fragment_id, fragment_name);
  _fragment_types.push_back(fragment);
  return fragment;
}

Segment *Md2QmEngine::AddSegmentType(int segment_id,
                                     tools::Property *property) {
  std::string segment_name = property->get("name").as<std::string>();
  Segment *segment = new Segment(segment_id, segment_name);
  _segment_types.push_back(segment);
  return segment;
}

Molecule *Md2QmEngine::AddMoleculeType(int molecule_id,
                                       tools::Property *property) {
  std::string molecule_name = property->get("name").as<std::string>();
  std::string molecule_mdname = property->get("mdname").as<std::string>();

  Molecule *molecule = new Molecule(molecule_id, molecule_name);
  _molecule_types.push_back(molecule);
  // TODO check if the name is already there
  _map_MoleculeMDName_MoleculeName[molecule_mdname] = molecule_name;
  _map_MoleculeName_MoleculeType[molecule_name] = molecule;
  return molecule;
}

SegmentType *Md2QmEngine::AddQMUnit(int qmunit_id, tools::Property *property) {
  std::string qmCoordsFile = "nofile";
  std::string orbitalsFile = "nofile";
  std::string basisSetName = "noname";
  std::vector<int> torbNrs;

  bool canRigidify = false;

  std::string qmunit_name = property->get("name").as<std::string>();

  if (property->exists("qmcoords")) {
    qmCoordsFile = property->get("qmcoords").as<std::string>();
    canRigidify = true;
  }
  if (property->exists("orbitals")) {
    orbitalsFile = property->get("orbitals").as<std::string>();
  }
  if (property->exists("basisset")) {
    basisSetName = property->get("basisset").as<std::string>();
  }

  SegmentType *qmUnit =
      new SegmentType(qmunit_id, qmunit_name, basisSetName, orbitalsFile,
                      qmCoordsFile, canRigidify);
  _qmUnits.push_back(qmUnit);
  return qmUnit;
}

Molecule *Md2QmEngine::getMoleculeType(const std::string &name) {
  try {
    return _map_MoleculeName_MoleculeType.at(name);
  } catch (const std::exception &out_of_range) {
    std::cout << "WARNING: Molecule '" << name
              << "' not included in mapping definition. Skipping... ";
    std::cout << std::endl;
    return NULL;
  }
}

const std::string &Md2QmEngine::getMoleculeName(const std::string &mdname) {
  try {
    return _map_MoleculeMDName_MoleculeName.at(mdname);
  } catch (const std::exception &out_of_range) {
    return mdname;
  }
}

Atom *Md2QmEngine::getAtomType(const std::string &molMdName, int resNr,
                               const std::string &mdAtomName) {
  return this->_map_mol_resNr_atm_atmType.at(molMdName).at(resNr).at(
      mdAtomName);
}
// TODO move to filereader
void Md2QmEngine::ReadXYZFile(
    std::string &file,
    std::map<int, std::pair<std::string, tools::vec> > &intCoords) {

  std::string line;
  std::ifstream intt;
  intt.open(file.c_str());
  if (!intt)
    throw std::runtime_error(std::string("Error reading coordinates from: ") +
                             file);
  std::getline(intt, line);
  tools::Tokenizer tok1(line, " \t");
  std::vector<std::string> line1;
  tok1.ToVector(line1);
  if (line1.size() != 1) {
    throw std::runtime_error(
        "First line of xyz file should contain number of atoms, nothing else.");
  }
  std::getline(intt, line);  // Comment line

  int atomCount = 0;
  if (intt.is_open()) {
    while (intt.good()) {
      std::getline(intt, line);
      std::vector<std::string> split;
      tools::Tokenizer toker(line, " \t");
      toker.ToVector(split);
      if (split.size() < 4) {
        continue;
      }
      // Interesting information written here: e.g. 'C 0.000 0.000 0.000'
      atomCount++;
      std::string element = split[0];
      double x = stod(split[1]) / 10.;  //°A to NM
      double y = stod(split[2]) / 10.;
      double z = stod(split[3]) / 10.;
      tools::vec qmPos = tools::vec(x, y, z);
      std::pair<std::string, tools::vec> qmTypePos(element, qmPos);
      intCoords[atomCount] = qmTypePos;
    }
  } else {
    throw std::runtime_error("No such file: '" + file + "'.");
  }
}

void Md2QmEngine::PrintInfo() {

  std::cout << "Summary ~~~~~"
               "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << std::endl;
  std::cout << "Created " << _molecule_types.size() << " molecule type(s): ";
  for (Molecule *mol : _molecule_types) {
    std::cout << "[ " << mol->getName() << " ]";
  }
  std::cout << std::endl
            << "with a total of " << _segment_types.size() << " segment(s), "
            << _fragment_types.size() << " fragments, " << _atom_types.size()
            << " atoms. \n"
            << std::endl;

  for (const std::pair<std::string, std::string> &mssit :
       _map_MoleculeMDName_MoleculeName) {
    std::cout << "MD [ " << mssit.first << " ] mapped to "
              << "QM [ " << mssit.second << " ] \n"
              << std::endl;
  }

  for (Molecule *mol : _molecule_types) {
    std::cout << "[ " << mol->getName() << " ]" << std::endl;

    for (Segment *seg : mol->Segments()) {
      std::cout << " - Segment [ " << seg->getName() << " ]"
                << " ID " << seg->getId() << std::endl;

      for (Fragment *frag : seg->Fragments()) {
        std::cout << "   - Fragment [ " << frag->getName() << " ]"
                  << " ID " << frag->getId() << ": " << frag->Atoms().size()
                  << " atoms " << std::endl;
      }
    }
  }

  if (!votca::tools::globals::verbose) {
    return;
  }

  std::cout << std::endl
            << "Mapping table"
               " ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            << std::endl;
  for (auto &map0 : this->_map_mol_resNr_atm_atmType) {
    for (auto &map1 : map0.second) {
      for (auto &map2 : map1.second) {

        printf(
            "MD Molecule %4s | Residue %2d | Atom %3s "
            "| ID %3d => QM ID %3d \n",
            map2.first.c_str(), map2.second->getResnr(), map2.first.c_str(),
            map2.second->getId(), map2.second->getQMId());
      }
    }
  }
}

}  // namespace xtp
}  // namespace votca
