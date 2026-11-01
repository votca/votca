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
      if (votca::Log::verbose()) {
        std::cout << "... ... processing mapping information for segment "
                  << segname << std::endl;
      }
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
    if (votca::Log::verbose()) {
      std::cout << "... molecule " << molname << " has " << segnames.size()
                << " segment(s)" << std::endl;
    }
    SegsinMol[molname] = segnames;
  }

  // Lookup bead ID -> IDs of directly bonded beads, built from the bonded
  // interactions of the MD topology (purely geometric, independent of the
  // segment/fragment mapping). Only 2-bead interactions (bonds) are used;
  // angles and dihedrals are skipped.
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

  // Always printed (not only in verbose mode): without bonds no external
  // bonds can be detected below, so the H-saturation of cut segment
  // boundaries in IPodCoupling::EvalJob (see transport_theory.rst,
  // "H-Saturation of Cut Segment Boundaries") silently leaves dangling
  // valences. This happened e.g. with GMXTopologyReader before it read bonds.
  if (bead_bonded_partners.empty()) {
    std::cout
        << "\nWARNING: the MD topology contains no bonds (no 2-bead bonded "
           "interactions found); the topology reader may only provide atom "
           "positions. External bonds at cut segment boundaries cannot be "
           "detected, so automatic H-saturation (e.g. in ipodcoupling) will "
           "do nothing. Use a topology file/reader that provides bonds, or "
           "xtp_map --guess-bonds."
        << std::endl;
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

          // A bonded partner in a different segment (or in no segment,
          // e.g. an unmapped atom) marks a bond crossing the segment
          // boundary. Its direction is stored on the atom (in MD
          // orientation; SegmentMapper rotates it with the QM template)
          // for later H-saturation. Only the first such partner is stored.
          auto it = bead_bonded_partners.find(bead->getId());
          if (it != bead_bonded_partners.end()) {
            for (Index partner_id : it->second) {
              const csg::Bead* partner_bead = top.getBead(partner_id);
              // MoleculeByIndex() is not const, so the partner's molecule
              // is looked up by matching getMoleculeId() against
              // Molecule::getId(), the same ID used for AddMoleculeId above.
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
              Index partner_offset = DetermineAtomNumOffset(
                  partner_mol, MolToAtomIds[partner_mol->getName()]);
              std::string partner_segname =
                  MolToSegMap[partner_mol->getName()]
                             [partner_bead->getId() - partner_offset];
              if (partner_segname != segname) {
                Eigen::Vector3d direction =
                    (partner_bead->getPos() - bead->getPos()) *
                    tools::conv::nm2bohr;
                atom.setExternalBondDirection(direction, partner_bead->getId());
                break;
              }
            }
          }

          // Store all bonded partners (not only external ones); the
          // fragment relaxation in FragmentSaturator needs the full internal
          // connectivity (see Atom::getBondedPartnerIds). These are MD-level
          // IDs, translated to QM-level IDs in SegmentMapper.
          if (it != bead_bonded_partners.end()) {
            for (Index partner_id : it->second) {
              atom.AddBondedPartner(partner_id);
            }
          }

          this_segment.push_back(atom);
        }
      }
      // add segment to topology
      topology_segments.push_back(this_segment);
    }
  }

  // Second pass: resolve the MD atom ID of each external-bond partner
  // (getExternalBondPartnerAtomId(), transient) to the ID of the segment it
  // belongs to (getExternalBondPartnerSegmentId(), persisted). This cannot
  // be done in the loop above, since later segments do not exist yet there.
  // The segment ID identifies which segment a bond crosses into, e.g. to
  // decide whether the bond is satisfied within an assembled supermolecule.
  std::map<Index, Index> md_atom_id_to_segment_id;
  for (const Segment& seg : xtptop.Segments()) {
    for (const Atom& segatom : seg) {
      md_atom_id_to_segment_id[segatom.getId()] = seg.getId();
    }
  }
  for (Segment& seg : xtptop.Segments()) {
    for (Atom& segatom : seg) {
      if (!segatom.hasExternalBond()) {
        continue;
      }
      auto it =
          md_atom_id_to_segment_id.find(segatom.getExternalBondPartnerAtomId());
      if (it != md_atom_id_to_segment_id.end()) {
        segatom.setExternalBondPartnerSegmentId(it->second);
      }
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
