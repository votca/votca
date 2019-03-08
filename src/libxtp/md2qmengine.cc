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
        std::vector<std::string> atomnames;
        tok_md_atoms.ToVector(atomnames);
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

int Md2QmEngine::DetermineResNumOffset(const csg::Molecule* mol,
                                       const std::vector<int>& resnums_map) {
  std::vector<int> resnums;
  for (const csg::Bead* bead : mol->Beads()) {
    resnums.push_back(bead->getResnr());
  }
  std::sort(resnums.begin(), resnums.end());
  int offset = resnums[0] - resnums_map[0];
  for (unsigned i = 0; i < resnums.size(); i++) {
    if (resnums[i] - resnums_map[i] != offset) {
      throw std::runtime_error(
          "Residue offset could not be determined, either our MD trajectory or "
          "your mapping file have wrong Residue ids");
    }
  }
  return offset;
}

Topology Md2QmEngine::map(const csg::Topology& top) {

  tools::Property topology_map;
  tools::load_property_from_xml(topology_map, _mapfile);
  CheckMappingFile(topology_map);

  Topology xtptop;
  xtptop.setStep(top.getStep());
  xtptop.setTime(top.getTime());
  xtptop.setBox(top.getBox(), top.getBoxType());

  // which segment does an atom belong to molname resnum name segid
  std::map<std::string, std::map<int, std::map<std::string, int> > >
      MolToSegMap;

  // residue IDs of Molecules
  std::map<std::string, std::vector<int> > MolToResNum;
  // names of segments in one molecule;
  std::map<std::string, std::vector<std::string> > SegsinMol;

  std::string molkey = "topology.molecules.molecule";
  std::vector<tools::Property*> molecules = topology_map.Select(molkey);
  std::string segkey = "segments.segment";
  for (tools::Property* mol : molecules) {
    std::string molname = mol->get("mdname").as<std::string>();
    std::vector<tools::Property*> segments = mol->Select(segkey);
    int segid = 0;
    std::vector<std::string> segnames;
    std::vector<int> resnums;
    for (tools::Property* seg : segments) {
      segnames.push_back(seg->get("name").as<std::string>());
      std::string fragkey = "fragments.fragment";
      std::vector<tools::Property*> fragments = seg->Select(fragkey);
      for (tools::Property* frag : fragments) {
        std::string mdatoms = frag->get("mdatoms").as<std::string>();
        tools::Tokenizer tok_md_atoms(mdatoms, " \t\n");
        std::vector<std::string> atomnames;
        tok_md_atoms.ToVector(atomnames);
        for (const std::string& atomname : atomnames) {
          tools::Tokenizer tok_atom_name(atomname, ":");
          std::vector<std::string> entries;
          tok_atom_name.ToVector(entries);
          if (entries.size() != 3) {
            throw std::runtime_error("Atom entry " + atomname +
                                     " is not well formatted");
          }
          int resnr = std::stoi(entries[0]);
          resnums.push_back(resnr);
          std::string resname = entries[1];
          std::string name = entries[2];
          MolToSegMap[molname][resnr][name] = segid;
        }
      }
      segid++;
    }
    std::sort(resnums.begin(), resnums.end());
    MolToResNum[molname] = resnums;
    SegsinMol[molname] = segnames;
  }

  int atomid = 0;
  for (const csg::Molecule* mol : top.Molecules()) {
    const std::vector<int> resnums_map = MolToResNum[mol->getName()];
    const std::vector<std::string> segnames = SegsinMol[mol->getName()];

    std::vector<Segment*> segments;  // we first add them to topology and then
                                     // modify them via pointers;
    for (const std::string& segname : segnames) {
      segments.push_back(&xtptop.AddSegment(segname));
    }

    // we have to figure out how the Residue numbers change from molecule to
    // molecule to get the correct mapping information this does not require for
    // the atoms to be sorted according to resnum
    int ResNumOffset = DetermineResNumOffset(mol, resnums_map);
    for (const csg::Bead* bead : mol->Beads()) {
      Segment* seg =
          segments[MolToSegMap[mol->getName()][bead->getResnr() - ResNumOffset]
                              [bead->getName()]];
      Atom atom(bead->getResnr(), bead->getName(), atomid,
                bead->getPos() * tools::conv::nm2bohr);
      seg->push_back(atom);
    }
  }

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
