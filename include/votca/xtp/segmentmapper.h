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
  };
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

using QMMapper = SegmentMapper<QMMolecule>;
using StaticMapper = SegmentMapper<StaticSegment>;
using PolarMapper = SegmentMapper<PolarSegment>;
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SEGMENTMAPPER_H
