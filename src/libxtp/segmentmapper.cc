/*
 *            Copyright 2016 The MUSCET Development Team
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

#include <regex>
#include <votca/xtp/segmentmapper.h>

namespace votca {
namespace xtp {

template <class AtomContainer>
SegmentMapper<AtomContainer>::SegmentMapper(Logger& log) : _log(log) {
  FillMap();
};

template <class AtomContainer>
AtomContainer SegmentMapper<AtomContainer>::map(const Segment& seg,
                                                const SegId& segid) const {
  if (segid.hasFile()) {
    std::string filename = segid.FileName();
    filename = std::regex_replace(filename, std::regex("\\$SEGID"),
                                  std::to_string(seg.getId()));
    filename =
        std::regex_replace(filename, std::regex("\\$SEGNAME"), seg.getName());
    return map(seg, filename);
  } else {
    QMState state = segid.getQMState();
    return map(seg, state);
  }
}

template <class AtomContainer>
template <typename T>
Eigen::Vector3d SegmentMapper<AtomContainer>::CalcWeightedPos(
    const std::vector<double>& weights, const T& atoms) const {
  Eigen::Vector3d map_pos = Eigen::Vector3d::Zero();
  double map_weight = 0.0;
  for (unsigned i = 0; i < atoms.size(); i++) {
    map_pos += atoms[i]->getPos() * weights[i];
    map_weight += weights[i];
  }
  return map_pos = map_pos / map_weight;
}

template <class AtomContainer>
void SegmentMapper<AtomContainer>::ParseFragment(Seginfo& seginfo,
                                                 const tools::Property& frag) {
  tools::Tokenizer tok_map_atoms(
      frag.get(_mapatom_xml["atoms"]).as<std::string>(), " \t\n");
  std::vector<std::string> map_atoms = tok_map_atoms.ToVector();
  tools::Tokenizer tok_md_atoms(frag.get("mdatoms").as<std::string>(), " \t\n");
  std::vector<std::string> md_atoms = tok_md_atoms.ToVector();

  if (md_atoms.size() != map_atoms.size()) {
    throw std::runtime_error(
        "Mapping for segment " + seginfo.segname + " fragment " +
        frag.get("name").as<std::string>() +
        " does not have same numbers of md and " + _mapatom_xml["atoms"] +
        "."
        "If you want to leave a qmatom out, place a ':' instead");
  }

  tools::Tokenizer tok_weights(getWeights(frag), " \t\n");
  std::vector<double> weights;
  tok_weights.ConvertToVector(weights);

  if (md_atoms.size() != weights.size()) {
    throw std::runtime_error("Mapping for segment " + seginfo.segname +
                             " fragment " + frag.get("name").as<std::string>() +
                             " does not have same numbers of md and weights. "
                             "If you want to leave a " +
                             _mapatom_xml["name"] +
                             " out, place a '0' instead");
  }

  FragInfo mapfragment;
  std::vector<int> mapatom_ids;

  for (unsigned i = 0; i < map_atoms.size(); i++) {
    const std::string& map_string = map_atoms[i];
    const std::string& md_string = md_atoms[i];
    const double& weight = weights[i];
    MD_atom_id md_result = StringToMDIndex(md_string);
    seginfo.mdatoms.push_back(md_result);

    if (map_string == ":") {
      continue;
    }
    mapsite_id map_result = StringToMapIndex(map_string);
    mapatom_ids.push_back(map_result.first);
    seginfo.mapatoms.push_back(map_result);
    if (Atom::GetElementFromMDName(md_result.second) != map_result.second) {
      XTP_LOG_SAVE(logINFO, _log)
          << "WARNING: mdatom'" << md_result.second << "' and "
          << _mapatom_xml["name"] << " '" << map_result.second
          << "' do not have same element" << std::flush;
    }
    mapfragment._mapatom_ids.push_back(map_result);
    mapfragment._weights.push_back(weight);
    mapfragment._mdatom_ids.push_back(md_result);
  }

  tools::Tokenizer tok_frame(getFrame(frag), " \t\n");
  std::vector<int> frame;
  tok_frame.ConvertToVector(frame);
  if (frame.size() > 3) {
    throw std::runtime_error(
        "Local frame for segment " + seginfo.segname + " fragment " +
        frag.get("name").as<std::string>() +
        " has more than 3 atoms, please specify only up to three atoms");
  } else if (frame.empty()) {
    throw std::runtime_error("No local frame for segment " + seginfo.segname +
                             " fragment " + frag.get("name").as<std::string>() +
                             " specified");
  }
  for (int atomid : frame) {
    if (std::find(mapatom_ids.begin(), mapatom_ids.end(), atomid) ==
        mapatom_ids.end()) {
      throw std::runtime_error("Atom " + std::to_string(atomid) +
                               " in local frame cannot be found in " +
                               _mapatom_xml["atoms"] + ".");
    }
  }

  mapfragment._map_local_frame = frame;
  seginfo.fragments.push_back(mapfragment);
}

template <class AtomContainer>
void SegmentMapper<AtomContainer>::LoadMappingFile(const std::string& mapfile) {

  tools::Property topology_map;
  tools::load_property_from_xml(topology_map, mapfile);

  std::string molkey = "topology.molecules.molecule";
  std::vector<tools::Property*> molecules = topology_map.Select(molkey);
  std::string segkey = "segments.segment";
  for (tools::Property* mol : molecules) {
    std::vector<tools::Property*> segments = mol->Select(segkey);
    for (tools::Property* seg : segments) {
      Seginfo seginfo;

      std::string coordfile_key = _mapatom_xml["coords"] + "_*";
      std::vector<tools::Property*> files = seg->Select(coordfile_key);
      for (tools::Property* file : files) {
        seginfo.coordfiles[file->name()] = file->as<std::string>();
      }
      std::string segname = seg->get("name").as<std::string>();
      seginfo.segname = segname;
      std::string fragkey = "fragments.fragment";

      seginfo.map2md = seg->ifExistsReturnElseThrowRuntimeError<bool>("map2md");

      std::vector<tools::Property*> fragments = seg->Select(fragkey);
      for (tools::Property* frag : fragments) {
        ParseFragment(seginfo, *frag);
      }

      int map_atom_min_id =
          std::min_element(seginfo.mapatoms.begin(), seginfo.mapatoms.end(),
                           [](const mapsite_id& a, const mapsite_id& b) {
                             return a.first < b.first;
                           })
              ->first;
      if (map_atom_min_id != 0) {
        throw std::runtime_error(_mapatom_xml["atoms"] + " for segment " +
                                 seginfo.segname +
                                 " do not start at zero index. Each segment "
                                 "should have its own coordinate file");
      }

      seginfo.minmax = CalcResidueRange(seginfo.mdatoms);
      _segment_info[segname] = seginfo;
    }
  }
}

template <class AtomContainer>
std::pair<int, std::string> SegmentMapper<AtomContainer>::StringToMapIndex(
    const std::string& map_string) const {
  tools::Tokenizer tok(map_string, ":");
  std::vector<std::string> result = tok.ToVector();
  if (result.size() != 2) {
    throw std::runtime_error("Entry " + map_string +
                             " is not properly formatted.");
  }
  return std::pair<int, std::string>(std::stoi(result[0]), result[1]);
}
template <class AtomContainer>
MD_atom_id SegmentMapper<AtomContainer>::StringToMDIndex(
    const std::string& md_string) const {
  tools::Tokenizer tok(md_string, ":");
  std::vector<std::string> result = tok.ToVector();
  if (result.size() != 3) {
    throw std::runtime_error("Entry " + md_string +
                             " is not properly formatted.");
  }
  return MD_atom_id(std::stoi(result[0]), result[2]);
}

template <class AtomContainer>
std::pair<int, int> SegmentMapper<AtomContainer>::CalcResidueRange(
    const std::vector<MD_atom_id>& seg) const {
  int max_res_id =
      std::min_element(seg.begin(), seg.end(),
                       [](const MD_atom_id& a, const MD_atom_id& b) {
                         return a.first < b.first;
                       })
          ->first;
  int min_res_id =
      std::min_element(seg.begin(), seg.end(),
                       [](const MD_atom_id& a, const MD_atom_id& b) {
                         return a.first < b.first;
                       })
          ->first;
  return std::pair<int, int>(min_res_id, max_res_id);
}
template <class AtomContainer>
std::pair<int, int> SegmentMapper<AtomContainer>::CalcResidueRange(
    const Segment& seg) const {
  int max_res_id = std::min_element(seg.begin(), seg.end(),
                                    [](const Atom& a, const Atom& b) {
                                      return a.getResnr() < b.getResnr();
                                    })
                       ->getResnr();

  int min_res_id = std::min_element(seg.begin(), seg.end(),
                                    [](const Atom& a, const Atom& b) {
                                      return a.getResnr() < b.getResnr();
                                    })
                       ->getResnr();
  return std::pair<int, int>(min_res_id, max_res_id);
}

template <class AtomContainer>
void SegmentMapper<AtomContainer>::PlaceMapAtomonMD(
    const std::vector<mapAtom*>& fragment_mapatoms,
    const std::vector<const Atom*>& fragment_mdatoms) const {
  for (unsigned i = 0; i < fragment_mapatoms.size(); i++) {
    const Atom* a = fragment_mdatoms[i];
    mapAtom* b = fragment_mapatoms[i];
    b->setPos(a->getPos());
  }
}

template <class AtomContainer>
int SegmentMapper<AtomContainer>::FindVectorIndexFromAtomId(
    int atomid, const std::vector<mapAtom*>& fragment_mapatoms) const {
  unsigned i = 0;
  for (; i < fragment_mapatoms.size(); i++) {
    if (fragment_mapatoms[i]->getId() == atomid) {
      break;
    }
  }
  return i;
}
template <class AtomContainer>
void SegmentMapper<AtomContainer>::MapMapAtomonMD(
    const FragInfo& frag, const std::vector<mapAtom*>& fragment_mapatoms,
    const std::vector<const Atom*>& fragment_mdatoms) const {
  std::vector<Eigen::Vector3d> local_map_frame;
  std::vector<Eigen::Vector3d> local_md_frame;
  for (int id : frag._map_local_frame) {
    int i = FindVectorIndexFromAtomId(id, fragment_mapatoms);
    local_map_frame.push_back(fragment_mapatoms[i]->getPos());
    local_md_frame.push_back(fragment_mdatoms[i]->getPos());
  }

  int symmetry = frag._map_local_frame.size();
  Eigen::Vector3d map_com = CalcWeightedPos(frag._weights, fragment_mapatoms);
  Eigen::Vector3d md_com = CalcWeightedPos(frag._weights, fragment_mdatoms);

  Eigen::Vector3d shift_map2md = md_com - map_com;

  Eigen::Matrix3d rot_map = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d rot_md = Eigen::Matrix3d::Identity();

  // building local frame
  if (symmetry > 1) {
    // middle atom is the origin
    Eigen::Vector3d x_map = local_map_frame[0] - local_map_frame[1];
    Eigen::Vector3d x_md = local_md_frame[0] - local_md_frame[1];
    Eigen::Vector3d y_map;
    Eigen::Vector3d y_md;
    Eigen::Vector3d z_map;
    Eigen::Vector3d z_md;
    if (symmetry == 3) {
      y_map = local_map_frame[2] - local_map_frame[1];
      y_md = local_md_frame[2] - local_md_frame[1];
      z_map = x_map.cross(y_map);
      z_md = x_md.cross(y_md);
      y_map = z_map.cross(x_map);
      y_md = z_md.cross(x_md);

    } else {
      Eigen::Vector3d unit = Eigen::Vector3d::UnitX();
      if (std::abs(unit.dot(x_map) / x_map.norm()) < 1e-6) {
        unit = Eigen::Vector3d::UnitY();
      }
      y_map = (x_map.cross(unit));
      if (std::abs(unit.dot(x_md) / x_md.norm()) < 1e-6) {
        unit = Eigen::Vector3d::UnitX();
      }
      y_md = (x_md.cross(unit));
      z_map = x_map.cross(y_map);
      z_md = x_md.cross(y_md);
    }

    if (x_map.squaredNorm() < 1e-18 || y_map.squaredNorm() < 1e-18 ||
        z_map.squaredNorm() < 1e-18) {
      throw std::runtime_error(
          _mapatom_xml.at("tag") +
          " basis vectors are very small, choose different local basis");
    }
    rot_map.col(0) = x_map.normalized();
    rot_map.col(1) = y_map.normalized();
    rot_map.col(2) = z_map.normalized();

    if (x_md.squaredNorm() < 1e-18 || y_md.squaredNorm() < 1e-18 ||
        z_md.squaredNorm() < 1e-18) {
      throw std::runtime_error(
          "MD basis vectors are very small, choose different local basis");
    }
    rot_md.col(0) = x_md.normalized();
    rot_md.col(1) = y_md.normalized();
    rot_md.col(2) = z_md.normalized();
  }
  Eigen::Matrix3d rotateMAP2MD = rot_md * rot_map.transpose();
  for (mapAtom* atom : fragment_mapatoms) {
    if (getRank(*atom) > 0 && symmetry < 3) {
      throw std::runtime_error(
          "Local frame has less than 3 atoms, thus higher rank multipoles "
          "cannot be mapped.");
    }
    atom->Translate(shift_map2md);
    atom->Rotate(rotateMAP2MD, md_com);
  }
}
template <class AtomContainer>
AtomContainer SegmentMapper<AtomContainer>::map(const Segment& seg,
                                                QMState state) const {
  Seginfo seginfo = _segment_info.at(seg.getName());
  std::string coordsfiletag =
      _mapatom_xml.at("coords") + "_" + state.ToString();
  if (seginfo.coordfiles.count(coordsfiletag) == 0) {
    throw std::runtime_error(
        "Could not find a coordinate file for segment/state: " + coordsfiletag);
  }
  std::string coordsfilename = seginfo.coordfiles.at(coordsfiletag);
  return map(seg, coordsfilename);
}

template <class AtomContainer>
AtomContainer SegmentMapper<AtomContainer>::map(
    const Segment& seg, const std::string& coordfilename) const {
  Seginfo seginfo = _segment_info.at(seg.getName());

  if (int(seginfo.mdatoms.size()) != seg.size()) {
    throw std::runtime_error(
        "Segment '" + seg.getName() +
        "' does not contain the same number of atoms as mapping file: " +
        std::to_string(seginfo.mdatoms.size()) + " vs. " +
        std::to_string(seg.size()));
  }

  std::pair<int, int> minmax_map = seginfo.minmax;
  std::pair<int, int> minmax = CalcResidueRange(seg);

  if ((minmax_map.first - minmax_map.second) !=
      (minmax.first - minmax.second)) {
    throw std::runtime_error("Residue range for segment " + seg.getName() +
                             ":" + std::to_string(seg.getId()) +
                             " and the mapping do not agree: Segment[" +
                             std::to_string(minmax.first) + "," +
                             std::to_string(minmax.second) + "] Map[" +
                             std::to_string(minmax_map.first) + "," +
                             std::to_string(minmax_map.second) + "]");
  }

  int residueoffset = minmax.first - minmax_map.first;

  AtomContainer Result(seg.getName(), seg.getId());
  Result.LoadFromFile(coordfilename);

  if (int(seginfo.mapatoms.size()) != Result.size()) {
    throw std::runtime_error(
        _mapatom_xml.at("tag") + "Segment '" + seg.getName() +
        "' does not contain the same number of atoms as mapping file: " +
        std::to_string(seginfo.mapatoms.size()) + " vs. " +
        std::to_string(Result.size()));
  }

  for (FragInfo& frag : seginfo.fragments) {
    for (MD_atom_id& id : frag._mdatom_ids) {
      id.first += residueoffset;
    }

    std::vector<mapAtom*> fragment_mapatoms;
    for (const mapsite_id& id : frag._mapatom_ids) {
      if (id.second != Result[id.first].getElement()) {
        throw std::runtime_error("Element of mapping atom " +
                                 std::to_string(id.first) + ":" + id.second +
                                 " does not agree with Element of parsed Atom" +
                                 Result[id.first].getElement());
      }
      fragment_mapatoms.push_back(&Result[id.first]);
    }
    std::vector<const Atom*> fragment_mdatoms;
    for (const MD_atom_id& id : frag._mdatom_ids) {
      const Atom* atom = seg.getAtom(id);
      if (atom == nullptr) {
        throw std::runtime_error("Could not find an atom with name " +
                                 std::to_string(id.first) + ":" + id.second +
                                 " in segment " + seg.getName());
      }
      fragment_mdatoms.push_back(atom);
    }

    if (seginfo.map2md) {
      PlaceMapAtomonMD(fragment_mapatoms, fragment_mdatoms);
    } else {
      MapMapAtomonMD(frag, fragment_mapatoms, fragment_mdatoms);
    }
  }

  Result.getPos();
  return Result;
}

template class SegmentMapper<QMMolecule>;
template class SegmentMapper<StaticSegment>;
template class SegmentMapper<PolarSegment>;

}  // namespace xtp
}  // namespace votca
