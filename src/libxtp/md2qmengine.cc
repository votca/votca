/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/md2qmengine.h>

namespace votca {
namespace xtp {

void Md2QmEngine::CheckMappingFile(tools::Property& topology_map) const {
  std::string molkey = "topology.molecules.molecule";
  std::vector<tools::Property*> molecules = topology_map.Select(molkey);
  if (SameValueForMultipleEntries<std::string>(molecules, "mdname")) {
    throw std::runtime_error("Multiple molecules have same mdname");
  }
  std::string segkey = "segments.segment";
  std::vector<tools::Property*> segments_all;
  for (tools::Property* mol : molecules) {
    std::vector<tools::Property*> segments = mol->Select(segkey);
    if (SameValueForMultipleEntries<std::string>(segments, "name")) {
      throw std::runtime_error("Multiple segments in molecule:" +
                               mol->get("mdname").as<std::string>() +
                               " have same name");
    }
    segments_all.insert(segments_all.end(), segments.begin(), segments.end());
    for (tools::Property* seg : segments) {
      std::string fragkey = "fragments.fragment";
      std::vector<tools::Property*> fragments = seg->Select(fragkey);
      if (SameValueForMultipleEntries<std::string>(fragments, "name")) {
        throw std::runtime_error(
            "Multiple fragments have same name in molecule " +
            mol->get("mdname").as<std::string>() + " segment " +
            seg->get("name").as<std::string>());
      }

      std::vector<std::string> atomnames_seg;
      for (tools::Property* frag : fragments) {
        std::string mdatoms = frag->get("mdatoms").as<std::string>();
        tools::Tokenizer tok_md_atoms(mdatoms, " \t\n");
        std::vector<std::string> atomnames = tok_md_atoms.ToVector();
        atomnames_seg.insert(atomnames_seg.end(), atomnames.begin(),
                             atomnames.end());
      }
      std::sort(atomnames_seg.begin(), atomnames_seg.end());  // O(N log N)
      if (adjacent_find(atomnames_seg.begin(), atomnames_seg.end()) !=
          atomnames_seg.end()) {
        throw std::runtime_error(
            "Multiple mdatoms have same identifier in molecule " +
            mol->get("mdname").as<std::string>() + " segment " +
            seg->get("name").as<std::string>());
      }
    }
  }
  if (SameValueForMultipleEntries<std::string>(segments_all, "name")) {
    throw std::runtime_error("Multiple segments have same name");
  }
}

long Md2QmEngine::DetermineAtomNumOffset(
    const csg::Molecule* mol, const std::vector<Index>& atom_ids_map) const {
  std::vector<Index> IDs;
  IDs.reserve(mol->BeadCount());
  for (const csg::Bead* bead : mol->Beads()) {
    IDs.push_back(bead->getId());
  }
  std::sort(IDs.begin(), IDs.end());
  Index offset = IDs[0] - atom_ids_map[0];
  for (Index i = 1; i < Index(IDs.size()); i++) {
    if (IDs[i] - atom_ids_map[i] != offset) {
      throw std::runtime_error(
          "AtomIds offset could not be determined, either our MD trajectory or "
          "your mapping file have wrong Atom ids");
    }
  }
  return offset;
}

bool Md2QmEngine::CheckMolWhole(const Topology& top, const Segment& seg) const {
  Eigen::Vector3d CoM = seg.getPos();
  bool whole = true;
  for (const Atom& a : seg) {
    Eigen::Vector3d r = a.getPos() - CoM;
    Eigen::Vector3d r_pbc = top.PbShortestConnect(CoM, a.getPos());
    Eigen::Vector3d shift = r_pbc - r;
    if (shift.norm() > 1e-9) {
      whole = false;
      break;
    }
  }
  return whole;
}

void Md2QmEngine::MakeSegmentsWholePBC(Topology& top) const {
  for (Segment& seg : top.Segments()) {
    seg.calcPos();
    while (!CheckMolWhole(top, seg)) {
      Eigen::Vector3d CoM = seg.getPos();
      for (Atom& a : seg) {
        Eigen::Vector3d r = a.getPos() - CoM;
        Eigen::Vector3d r_pbc = top.PbShortestConnect(CoM, a.getPos());
        Eigen::Vector3d shift = r_pbc - r;
        if (shift.norm() > 1e-9) {
          a.Translate(shift);
        }
      }
      seg.calcPos();
    }
  }
}

Topology Md2QmEngine::map(const csg::Topology& top) const {

  tools::Property topology_map;
  topology_map.LoadFromXML(_mapfile);
  CheckMappingFile(topology_map);
  Topology xtptop;
  xtptop.setStep(top.getStep());
  xtptop.setTime(top.getTime());
  xtptop.setBox(top.getBox() * tools::conv::nm2bohr, top.getBoxType());

  // which segmentname does an atom belong to molname atomid
  std::map<std::string, std::map<long, std::string> > MolToSegMap;

  // which atomids belong to molname
  std::map<std::string, std::vector<Index> > MolToAtomIds;

  // names of segments in one molecule;
  std::map<std::string, std::vector<std::string> > SegsinMol;

  std::string molkey = "topology.molecules.molecule";
  std::vector<tools::Property*> molecules = topology_map.Select(molkey);
  std::string segkey = "segments.segment";
  for (tools::Property* mol : molecules) {
    std::string molname = mol->get("mdname").as<std::string>();
    std::vector<tools::Property*> segments = mol->Select(segkey);
    std::vector<std::string> segnames;
    std::vector<Index> atomids;
    for (tools::Property* seg : segments) {
      std::string segname = seg->get("name").as<std::string>();
      segnames.push_back(segname);
      std::string fragkey = "fragments.fragment";
      std::vector<tools::Property*> fragments = seg->Select(fragkey);
      for (tools::Property* frag : fragments) {
        std::string mdatoms = frag->get("mdatoms").as<std::string>();
        tools::Tokenizer tok_md_atoms(mdatoms, " \t\n");
        std::vector<std::string> atomnames;
        tok_md_atoms.ToVector(atomnames);
        for (const std::string& atomname : atomnames) {
          tools::Tokenizer tok_atom_name(atomname, ":");
          std::vector<std::string> entries = tok_atom_name.ToVector();
          if (entries.size() != 3) {
            throw std::runtime_error("Atom entry " + atomname +
                                     " is not well formatted");
          }
          // format should be RESNUM:ATOMNAME:ATOMID we do not care about the
          // first two
          Index atomid = 0;
          try {
            atomid = std::stoi(entries[2]);
          } catch (std::invalid_argument& e) {
            throw std::runtime_error("Atom entry " + atomname +
                                     " is not well formatted");
          }
          atomids.push_back(atomid);
          MolToSegMap[molname][atomid] = segname;
        }
      }
    }
    std::sort(atomids.begin(), atomids.end());
    MolToAtomIds[molname] = atomids;
    SegsinMol[molname] = segnames;
  }

  for (const csg::Molecule* mol : top.Molecules()) {
    const std::vector<std::string> segnames = SegsinMol[mol->getName()];

    std::map<std::string, Segment*> segments;  // we first add them to topology
                                               // and then modify them via
                                               // pointers;
    for (const std::string& segname : segnames) {
      segments[segname] = &xtptop.AddSegment(segname);
      segments[segname]->AddMoleculeId(mol->getId());
    }

    Index IdOffset = DetermineAtomNumOffset(mol, MolToAtomIds[mol->getName()]);
    for (const csg::Bead* bead : mol->Beads()) {
      Segment* seg =
          segments[MolToSegMap[mol->getName()][bead->getId() - IdOffset]];

      Atom atom(bead->getResnr(), bead->getName(), bead->getId(),
                bead->getPos() * tools::conv::nm2bohr, bead->getType());
      seg->push_back(atom);
    }
  }

  MakeSegmentsWholePBC(xtptop);

  return xtptop;
}

template <class T>
bool Md2QmEngine::SameValueForMultipleEntries(
    const std::vector<tools::Property*>& props, std::string valuetag) const {
  std::vector<T> entries;
  for (tools::Property* prop : props) {
    entries.push_back(prop->get(valuetag).as<T>());
  }
  std::sort(entries.begin(), entries.end());  // O(N log N)
  return adjacent_find(entries.begin(), entries.end()) != entries.end();
}

}  // namespace xtp
}  // namespace votca
