/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
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

// Local VOTCA includes
#include "votca/xtp/md2qmengine.h"

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
        std::vector<std::string> atomnames =
            frag->get("mdatoms").as<std::vector<std::string>>();
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

Index Md2QmEngine::DetermineAtomNumOffset(
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
  topology_map.LoadFromXML(mapfile_);
  CheckMappingFile(topology_map);
  Topology xtptop;
  xtptop.setStep(top.getStep());
  xtptop.setTime(top.getTime());
  xtptop.setBox(top.getBox() * tools::conv::nm2bohr, top.getBoxType());

  // which segmentname does an atom belong to molname atomid
  std::map<std::string, std::map<Index, std::string>> MolToSegMap;

  // which atomids belong to molname
  std::map<std::string, std::vector<Index>> MolToAtomIds;

  // names of segments in one molecule;
  std::map<std::string, std::vector<std::string>> SegsinMol;

  std::string molkey = "topology.molecules.molecule";
  std::vector<tools::Property*> molecules = topology_map.Select(molkey);
  std::string segkey = "segments.segment";

  for (tools::Property* mol : molecules) {
    // get the name of this molecule
    std::string molname = mol->get("mdname").as<std::string>();
    // get all segment-mapping info
    std::vector<tools::Property*> segments = mol->Select(segkey);
    std::vector<std::string> segnames;
    std::vector<Index> atomids;
    // now go through all the defined segments
    for (tools::Property* seg : segments) {
      // get the name of this segment and add to segnames vector
      std::string segname = seg->get("name").as<std::string>();
      segnames.push_back(segname);
      std::string fragkey = "fragments.fragment";
      // get all fragement mapping info
      std::vector<tools::Property*> fragments = seg->Select(fragkey);
      // go over all fragments in this segement
      for (tools::Property* frag : fragments) {
        // get all mdatom names from this fragment
        std::vector<std::string> atomnames =
            frag->get("mdatoms").as<std::vector<std::string>>();
        // go over all atoms
        for (const std::string& atomname : atomnames) {
          // split atom entry at :
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
          if (votca::Log::verbose()) {
            std::cout << "... ... processing mapping information for atom "
                      << atomname << " with ID " << atomid << std::endl;
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

  // Build a direct, one-time "bead ID -> directly-bonded partner bead
  // IDs" lookup, from the MD-level topology's own real, actual bond
  // connectivity (csg::Topology::BondedInteractions(), confirmed
  // directly, by reading csg's own interaction.h/topology.h, to
  // already exist and be available exactly here -- this is a purely
  // *geometric* bond list, from the original MD topology itself, NOT
  // yet aware of segment/fragment membership at all). Only 2-bead
  // ("B"/IBond-style) interactions are relevant here -- angles (3-bead)
  // and dihedrals (4-bead) do not represent a direct bond between two
  // atoms at all, so BeadCount() != 2 entries are skipped.
  std::map<Index, std::vector<Index>> bead_bonded_partners;
  for (csg::Interaction* interaction : top.BondedInteractions()) {
    if (interaction->BeadCount() != 2) {
      continue;
    }
    Index id1 = interaction->getBeadId(0);
    Index id2 = interaction->getBeadId(1);
    bead_bonded_partners[id1].push_back(id2);
    bead_bonded_partners[id2].push_back(id1);
  }

  // go through all molecules in MD topology
  for (const csg::Molecule& mol : top.Molecules()) {

    // lookup all segment *names* in this molecule
    const std::vector<std::string> segnames = SegsinMol[mol.getName()];
    std::vector<Segment>& topology_segments = xtptop.Segments();
    Index IdOffset = DetermineAtomNumOffset(&mol, MolToAtomIds[mol.getName()]);

    if (votca::Log::verbose()) {
      std::cout << "... Mapping molecule " << mol.getId() << ", name "
                << mol.getName() << ", # of segments " << segnames.size()
                << ", atomID offset " << IdOffset << std::endl;
    }

    for (const std::string& segname : segnames) {

      Index segid = topology_segments.size();
      // construct a segment
      Segment this_segment = Segment(segname, segid);
      this_segment.AddMoleculeId(mol.getId());

      // create atomlist
      for (const csg::Bead* bead : mol.Beads()) {
        // check if it belongs to this segment, and add it
        if (segname == MolToSegMap[mol.getName()][bead->getId() - IdOffset]) {
          Atom atom(bead->getResnr(), bead->getName(), bead->getId(),
                    bead->getPos() * tools::conv::nm2bohr, bead->getType());

          // Check each of this bead's own, real, direct bonded
          // partners (from the lookup built above): if a partner's
          // own segment assignment differs from this atom's own
          // segname (or the partner has no segment assignment at
          // all -- e.g. an unmapped, non-charge-transport-relevant
          // spectator atom), the bond crosses the segment boundary --
          // record the direction toward it directly on this atom, for
          // later use (H-saturation) once this direction has also
          // been carried through SegmentMapper's own, later,
          // rigid-body MD->QM-template transform. Only the first such
          // partner found is recorded (an atom with more than one
          // external bond is rare, and not handled specially here).
          auto it = bead_bonded_partners.find(bead->getId());
          if (it != bead_bonded_partners.end()) {
            for (Index partner_id : it->second) {
              const csg::Bead* partner_bead = top.getBead(partner_id);
              // MoleculeByIndex() is not const-qualified (confirmed
              // directly, from a real compile error), and this
              // function only ever sees top as const -- so the
              // partner's own parent molecule is found directly here
              // instead, by matching getMoleculeId() against each
              // molecule's own getId() (the same, canonical way this
              // function itself already identifies molecules, per its
              // own, pre-existing this_segment.AddMoleculeId(mol.getId())
              // call above), via the const-compatible Molecules()
              // overload -- this also sidesteps a second, separate
              // uncertainty MoleculeByIndex() would have carried:
              // whether its own index parameter expects a molecule ID
              // or an array position, which are not necessarily the
              // same thing.
              const csg::Molecule* partner_mol = nullptr;
              for (const csg::Molecule& candidate : top.Molecules()) {
                if (candidate.getId() == partner_bead->getMoleculeId()) {
                  partner_mol = &candidate;
                  break;
                }
              }
              if (partner_mol == nullptr) {
                continue;
              }
              Index partner_offset =
                  DetermineAtomNumOffset(partner_mol,
                                         MolToAtomIds[partner_mol->getName()]);
              std::string partner_segname =
                  MolToSegMap[partner_mol->getName()][partner_bead->getId() -
                                                       partner_offset];
              if (partner_segname != segname) {
                Eigen::Vector3d direction =
                    (partner_bead->getPos() - bead->getPos()) *
                    tools::conv::nm2bohr;
                atom.setExternalBondDirection(direction);
                break;
              }
            }
          }

          this_segment.push_back(atom);
        }
      }
      // add segment to topology
      topology_segments.push_back(this_segment);
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
