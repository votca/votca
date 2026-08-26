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
  std::cout << "First loop done..." << std::endl;
  std::sort(IDs.begin(), IDs.end());

      std::cout << IDs[0] << std::endl;
    std::cout << atom_ids_map[0] << std::endl;

  Index offset = IDs[0] - atom_ids_map[0];
  std::cout << offset << std::endl;
  for (Index i = 1; i < Index(IDs.size()); i++) {
    if (IDs[i] - atom_ids_map[i] != offset) {
      throw std::runtime_error(
          "AtomIds offset could not be determined, either our MD trajectory or "
          "your mapping file have wrong Atom ids");
    }
  }
  std::cout << offset << std::endl; 
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
    std::cout << "adding " << segnames[0] << std::endl;
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

  // Real, direct, always-visible (not gated behind -v/verbose at all)
  // warning, worked through directly with the user: if the underlying
  // MD topology genuinely has no real bond connectivity at all (e.g.
  // a topology reader that only ever provides atom positions, no real
  // bond data at all -- this exact situation is exactly what this
  // whole session's own real, direct debugging arc started from,
  // csg::GMXTopologyReader itself never having read any real bonds at
  // all, before that specific fix), external-bond detection itself
  // (below) can never find anything at all, meaning
  // IPodCoupling::EvalJob's own H-saturation at cut segment boundaries
  // (see transport_theory.rst's own "H-Saturation of Cut Segment
  // Boundaries" section) can never actually fire at all either --
  // silently producing dangling valences at every cut segment
  // boundary instead, with no other, direct signal of this at all
  // until (if ever) a user separately, manually notices something is
  // wrong much further downstream. Warning here instead, directly, at
  // the earliest point this is actually knowable at all (right after
  // bead_bonded_partners is built, genuinely reflecting the real,
  // complete state of top.BondedInteractions() itself), gives a real,
  // direct, upfront signal instead.
  if (bead_bonded_partners.empty()) {
    std::cout
        << "\nWARNING: the MD topology being mapped contains no real bond "
           "connectivity at all (no bonded interactions were found within "
           "it) -- this topology reader may only provide atom positions, "
           "not real, actual bond data. Automatic H-saturation of cut "
           "segment boundaries (used by e.g. the ipodcoupling calculator) "
           "will not be able to detect any external bonds at all, and will "
           "silently do nothing at all, rather than saturating anything -- "
           "check that the real, actual topology file/reader used here "
           "genuinely provides real bond data, not just atom positions."
        << std::endl;
  }

  // go through all molecules in MD topology
  for (const csg::Molecule& mol : top.Molecules()) {

    // lookup all segment *names* in this molecule
    const std::vector<std::string> segnames = SegsinMol[mol.getName()];
    std::cout << segnames.size() << std::endl;
    for (auto segment : segnames) {
      std::cout << segment << std::endl;
    }

    std::vector<Segment>& topology_segments = xtptop.Segments();
    //std::cout << mol.getName() << " " << MolToAtomIds[mol.getName()] << std::endl;
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

          // Records ALL of this bead's own, real, direct bonded
          // partners (not just the ones crossing a segment boundary,
          // unlike the external-bond-direction detection above, which
          // deliberately stops at the first one found) -- needed for
          // FragmentSaturator's own, planned OpenBabel-based
          // relaxation step, which needs a fragment's own, full,
          // internal connectivity to set up its own force field
          // correctly at all (see Atom::getBondedPartnerIds's own
          // header comment for why). These are RAW, MD-level partner
          // IDs -- not yet translated into QM-level IDs here, matching
          // the same "raw here, translated later, in SegmentMapper"
          // split already used for the external-bond direction.
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

  // Second pass: resolves each atom's own external-bond partner
  // (recorded above only as a raw MD-level ATOM id,
  // getExternalBondPartnerAtomId() -- deliberately transient, see its
  // own header comment) to the actual SEGMENT that partner atom
  // belongs to. Cannot be done inline, in the loop above -- at the
  // point a given atom's own external bond is first detected, later
  // segments (in molecule/segname iteration order) do not exist yet
  // at all, so there is no way yet to know which segment a partner
  // atom that happens to belong to one of them will end up in.
  //
  // First builds a direct "MD-level atom id -> segment id" lookup,
  // spanning every segment in the whole, now-complete xtptop, then
  // uses it to resolve every atom's own, already-recorded
  // getExternalBondPartnerAtomId() into the real, actual, persisted
  // getExternalBondPartnerSegmentId() -- needed downstream (a planned
  // linking-segment graph, and the decision of whether a given
  // external bond is already satisfied within an assembled
  // supermolecule, worked through directly with the user before this
  // was implemented) to identify not just THAT a given bond crosses a
  // segment boundary, but SPECIFICALLY WHICH segment it crosses into.
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
