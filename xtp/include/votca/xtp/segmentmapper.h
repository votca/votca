/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#pragma once
#ifndef VOTCA_XTP_SEGMENTMAPPER_H
#define VOTCA_XTP_SEGMENTMAPPER_H

// Standard includes
#include <type_traits>

// VOTCA includes
#include <votca/csg/pdbwriter.h>
#include <votca/tools/property.h>

// Local VOTCA includes
#include "classicalsegment.h"
#include "logger.h"
#include "qmmolecule.h"
#include "segid.h"
#include "topology.h"

namespace votca {
namespace xtp {
template <class AtomContainer>
class SegmentMapper {
 public:
  SegmentMapper(Logger& log);

  void LoadMappingFile(const std::string& mapfile);

  AtomContainer map(const Segment& seg, const SegId& segid) const;

  AtomContainer map(const Segment& seg, QMState state) const;

  AtomContainer map(const Segment& seg, const std::string& coordfilename) const;

 private:
  using mapAtom = typename AtomContainer::Atom_Type;

  using atom_id = std::pair<Index, std::string>;

  struct FragInfo {
    std::vector<double> weights;
    std::vector<atom_id> mapatom_ids;
    std::vector<atom_id> mdatom_ids;
    std::vector<Index> map_local_frame;
  };

  struct Seginfo {
    std::pair<Index, Index> minmax;
    std::vector<Index> mdatoms;
    std::vector<FragInfo> fragments;
    bool map2md;
    std::string segname;
    std::vector<double> weights;
    std::vector<atom_id> mapatoms;
    std::map<std::string, std::string> coordfiles;
    // chargeresponse_<state>: .crm files of the charge-response model
    // (polar segments only), next to multipoles_<state>
    std::map<std::string, std::string> crfiles;
  };

  // Attaches the charge-response matrix of the state, if the mapping file
  // has one, in the atom order of the mapped segment. No-op except for
  // PolarSegment (specialized in segmentmapper.cc).
  void AttachChargeResponse(AtomContainer&, const Seginfo&,
                            const QMState&) const {}
  std::map<std::string, std::string> mapatom_xml_;
  std::map<std::string, Seginfo> segment_info_;

  Index FindVectorIndexFromAtomId(
      Index atomid, const std::vector<mapAtom*>& fragment_mapatoms) const;

  void ParseFragment(Seginfo& seginfo, const tools::Property& frag);

  template <typename T>
  Eigen::Vector3d CalcWeightedPos(const std::vector<double>& weights,
                                  const T& atoms) const;

  void PlaceMapAtomonMD(const std::vector<mapAtom*>& fragment_mapatoms,
                        const std::vector<const Atom*>& fragment_mdatoms) const;

  void MapMapAtomonMD(const FragInfo& frag,
                      const std::vector<mapAtom*>& fragment_mapatoms,
                      const std::vector<const Atom*>& fragment_mdatoms) const;

  // Transfers the external-bond direction of the MD atom, rotated by rot
  // (no translation, a direction is translation invariant), and the ID of
  // the segment the bond crosses into. No-op by default; only QMAtom has
  // this interface (see the QMMolecule specialization below).
  void TransferExternalBondDirection(mapAtom*, const Atom*,
                                     const Eigen::Matrix3d&) const {}

  // Translates the MD-level bonded-partner IDs of md_atom into mapped IDs
  // via the segment-wide lookup MD atom ID -> mapped atom, and adds them to
  // map_atom. Partners outside the segment are skipped. No-op by default;
  // specialized for QMMolecule below.
  void TransferBondedPartners(mapAtom*, const Atom*,
                              const std::map<Index, mapAtom*>&) const {}

  // Rewrites the bonded-partner IDs of atom from template atom IDs to the
  // IDs of the returned container (old_to_new). No-op by default;
  // specialized for QMMolecule below.
  void RemapBondedPartners(mapAtom&, const std::map<Index, Index>&) const {}

  Logger& log_;
  std::pair<Index, Index> CalcAtomIdRange(const Segment& seg) const;
  std::pair<Index, Index> CalcAtomIdRange(const std::vector<Index>& seg) const;

  atom_id StringToMapIndex(const std::string& map_string) const;

  atom_id StringToMDIndex(const std::string& md_string) const;

  Index getRank(const mapAtom& atom) const { return atom.getRank(); }

  std::vector<double> getWeights(const tools::Property& frag) const;

  std::string getFrame(const tools::Property& frag) const {
    if (frag.exists(mapatom_xml_.at("frame"))) {
      return frag.get(mapatom_xml_.at("frame")).template as<std::string>();
    }
    return frag.get("localframe").template as<std::string>();
  }

  void FillMap() {
    mapatom_xml_["tag"] = "MP";
    mapatom_xml_["name"] = "MPole";
    mapatom_xml_["atoms"] = "mpoles";
    mapatom_xml_["coords"] = "multipoles";
    mapatom_xml_["weights"] = "mp_weights";
    mapatom_xml_["frame"] = "mp_localframe";
  }
};

template <>
inline void SegmentMapper<QMMolecule>::FillMap() {
  mapatom_xml_["tag"] = "QM";
  mapatom_xml_["name"] = "QMAtom";
  mapatom_xml_["atoms"] = "qmatoms";
  mapatom_xml_["coords"] = "qmcoords";
  mapatom_xml_["weights"] = "qm_weights";
  mapatom_xml_["frame"] = "qm_localframe";
}

template <>
inline Index SegmentMapper<QMMolecule>::getRank(const QMAtom&) const {
  return 0;
}

template <>
inline void SegmentMapper<QMMolecule>::TransferExternalBondDirection(
    QMAtom* map_atom, const Atom* md_atom, const Eigen::Matrix3d& rot) const {
  if (md_atom->hasExternalBond()) {
    map_atom->setExternalBondDirection(rot *
                                       md_atom->getExternalBondDirection());
    map_atom->setExternalBondPartnerSegmentId(
        md_atom->getExternalBondPartnerSegmentId());
  }
}

template <>
inline void SegmentMapper<QMMolecule>::TransferBondedPartners(
    QMAtom* map_atom, const Atom* md_atom,
    const std::map<Index, QMAtom*>& md_id_to_map_atom) const {
  for (Index i = 0; i < Atom::kMaxBondedPartners; i++) {
    Index partner_md_id = md_atom->getBondedPartnerIds()[i];
    if (partner_md_id == -1) {
      continue;
    }
    auto it = md_id_to_map_atom.find(partner_md_id);
    if (it != md_id_to_map_atom.end()) {
      map_atom->AddBondedPartner(it->second->getId());
    }
  }
}

template <>
inline void SegmentMapper<QMMolecule>::RemapBondedPartners(
    QMAtom& atom, const std::map<Index, Index>& old_to_new) const {
  Index ids[QMAtom::kMaxBondedPartners];
  Index n = 0;
  for (Index i = 0; i < QMAtom::kMaxBondedPartners; i++) {
    auto it = old_to_new.find(atom.getBondedPartnerIds()[i]);
    if (it != old_to_new.end()) {
      ids[n++] = it->second;
    }
  }
  for (Index i = n; i < QMAtom::kMaxBondedPartners; i++) {
    ids[i] = -1;
  }
  atom.setBondedPartnerIds(ids);
}

template <>
void SegmentMapper<PolarSegment>::AttachChargeResponse(
    PolarSegment& result, const Seginfo& seginfo, const QMState& state) const;

using QMMapper = SegmentMapper<QMMolecule>;
using StaticMapper = SegmentMapper<StaticSegment>;
using PolarMapper = SegmentMapper<PolarSegment>;
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENTMAPPER_H
