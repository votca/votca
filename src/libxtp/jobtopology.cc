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

#include <votca/xtp/checkpoint.h>
#include <votca/xtp/jobtopology.h>
#include <votca/xtp/mmregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/segmentmapper.h>

#include <numeric>

namespace votca {
namespace xtp {

void JobTopology::SortRegionsDefbyId(
    std::vector<const tools::Property*>& regions_def) const {
  std::sort(regions_def.begin(), regions_def.end(),
            [](const tools::Property* A, const tools::Property* B) {
              return (A->get("id").as<int>()) < (B->get("id").as<int>());
            });
}

void JobTopology::BuildRegions(const Topology& top,
                               const tools::Property& options) {

  std::vector<const tools::Property*> regions_def =
      options.Select("regions.region");

  CheckEnumerationOfRegions(regions_def);
  SortRegionsDefbyId(regions_def);

  std::vector<std::vector<int> > region_seg_ids =
      PartitionRegions(regions_def, top);

  // around this point the whole jobtopology will be centered
  CreateRegions(options, top, region_seg_ids);
  return;
}

template <class T>
T JobTopology::GetInputFromXMLorJob(const tools::Property* region_def,
                                    std::string keyword) const {
  T result;
  int id = region_def->ifExistsReturnElseThrowRuntimeError<int>("id");
  if (region_def->get(keyword).as<std::string>() != "jobfile") {
    result = region_def->get(keyword).as<T>();
  } else {
    const tools::Property& jobprop = _job.getInput();
    std::vector<const tools::Property*> regions =
        jobprop.Select("regions.region");
    for (const tools::Property* reg : regions) {
      int idregion_job = reg->ifExistsReturnElseThrowRuntimeError<int>("id");
      if (idregion_job == id) {
        result = reg->ifExistsReturnElseThrowRuntimeError<T>(keyword);
        break;
      }
    }
  }
  return result;
}

template <class T>
void JobTopology::ShiftPBC(const Topology& top, const Eigen::Vector3d& center,
                           T& mol) const {
  Eigen::Vector3d r_pbc = top.PbShortestConnect(center, mol.getPos());
  Eigen::Vector3d r = center - mol.getPos();
  Eigen::Vector3d shift = r_pbc - r;
  mol.Translate(shift);
}

void JobTopology::CreateRegions(
    const tools::Property& options, const Topology& top,
    const std::vector<std::vector<int> >& region_seg_ids) {
  std::string mapfile =
      options.ifExistsReturnElseThrowRuntimeError<std::string>("mapfile");
  std::vector<const tools::Property*> regions_def =
      options.Select("regions.region");
  // around this point the whole jobtopology will be centered for removing pbc
  Eigen::Vector3d center = top.getSegment(region_seg_ids[0][0]).getPos();

  for (const tools::Property* region_def : regions_def) {
    int id = region_def->ifExistsReturnElseThrowRuntimeError<int>("id");
    const std::vector<int>& seg_ids = region_seg_ids[id];
    std::string type =
        region_def->ifExistsReturnElseThrowRuntimeError<std::string>("type");
    std::string qmstate_string = "n";
    if (region_def->exists("geometry")) {
      qmstate_string =
          GetInputFromXMLorJob<std::string>(region_def, "geometry");
    }
    QMState geometry = QMState(qmstate_string);
    std::unique_ptr<Region> region;
    if (type == "gwbse" || type == "dft") {
      std::unique_ptr<QMRegion> qmregion =
          std::unique_ptr<QMRegion>(new QMRegion());
      QMMapper qmmapper(_log);
      qmmapper.LoadMappingFile(mapfile);
      for (int seg_index : seg_ids) {
        const Segment& segment = top.getSegment(seg_index);
        QMMolecule mol = qmmapper.map(segment, geometry);
        ShiftPBC(top, center, mol);
        qmregion->push_back(mol);
      }
      region = std::move(qmregion);
    } else if (type == "polar") {
      std::unique_ptr<PolarRegion> polarregion =
          std::unique_ptr<PolarRegion>(new PolarRegion());
      PolarMapper polmap(_log);
      polmap.LoadMappingFile(mapfile);
      for (int seg_index : seg_ids) {
        const Segment& segment = top.getSegment(seg_index);
        PolarSegment mol = polmap.map(segment, geometry);
        ShiftPBC(top, center, mol);
        polarregion->push_back(mol);
      }
      region = std::move(polarregion);
    } else if (type == "static") {
      std::unique_ptr<StaticRegion> staticregion =
          std::unique_ptr<StaticRegion>(new StaticRegion());
      StaticMapper staticmap(_log);
      staticmap.LoadMappingFile(mapfile);
      for (int seg_index : seg_ids) {
        const Segment& segment = top.getSegment(seg_index);
        StaticSegment mol = staticmap.map(segment, geometry);
        ShiftPBC(top, center, mol);
        staticregion->push_back(mol);
      }
      region = std::move(staticregion);

    } else {
      throw std::runtime_error("Region type not known!");
    }
    _regions.push_back(std::move(region));
  }
}

void JobTopology::WriteToPdb(std::string filename) const {

  csg::PDBWriter writer;
  writer.Open(filename, false);
  writer.WriteHeader("Job:" + std::to_string(this->_job.getId()));
  for (const std::unique_ptr<Region>& reg : _regions) {
    reg->WritePDB(writer);
  }
  writer.Close();
}

std::vector<std::vector<int> > JobTopology::PartitionRegions(
    const std::vector<const tools::Property*>& regions_def,
    const Topology& top) const {

  std::vector<std::vector<int> > segids_per_region;
  std::vector<bool> processed_segments =
      std::vector<bool>(false, top.Segments().size());

  for (const tools::Property* region_def : regions_def) {
    std::vector<int> seg_ids;
    if (!region_def->exists("segments") && region_def->exists("cutoff")) {
      throw std::runtime_error(
          "Region definition needs either segments or a cutoff to find "
          "segments");
    }

    if (region_def->exists("segments")) {
      seg_ids = GetInputFromXMLorJob<std::vector<int> >(region_def, "segments");

      for (int seg_id : seg_ids) {
        if (seg_id > int(top.Segments().size() - 1)) {
          throw std::runtime_error("Segment id is not in topology");
        }
        processed_segments[seg_id] = true;
      }
    }
    if (region_def->exists("cutoff")) {
      double cutoff = region_def->ifExistsReturnElseThrowRuntimeError<double>(
          "cutoff.radius");
      std::vector<int> center = seg_ids;
      if (region_def->exists("cutoff.region")) {
        int id = region_def->get("cutoff.region").as<int>();
        center = segids_per_region[id];
      }
      for (int id : center) {
        const Segment& center_seg = top.getSegment(id);
        for (const Segment& otherseg : top.Segments()) {
          if (center_seg.getId() == otherseg.getId() ||
              processed_segments[otherseg.getId()]) {
            continue;
          }
          if (top.GetShortestDist(center_seg, otherseg) < cutoff) {
            seg_ids.push_back(otherseg.getId());
            processed_segments[otherseg.getId()] = true;
          }
        }
      }
    }
    segids_per_region.push_back(seg_ids);
  }
  return segids_per_region;
}

void JobTopology::CheckEnumerationOfRegions(
    const std::vector<const tools::Property*>& regions_def) const {
  std::vector<int> reg_ids;
  for (const tools::Property* region_def : regions_def) {
    reg_ids.push_back(
        region_def->ifExistsReturnElseThrowRuntimeError<int>("id"));
  }

  std::vector<int> v(reg_ids.size());
  std::iota(v.begin(), v.end(), 0);
  if (std::is_permutation(reg_ids.begin(), reg_ids.end(), v.begin())) {
    throw std::runtime_error(
        "Region id definitions are not clear. You must start at id 0 and then "
        "ascending order. i.e. 0 1 2 3.");
  }
}

void JobTopology::WriteToHdf5(std::string filename) const {
  CheckpointFile cpf(filename, CheckpointAccessLevel::CREATE);
  CheckpointWriter a = cpf.getWriter();
  a(_job.getId(), "jobid");

  for (const auto& region : _regions) {
    CheckpointWriter w =
        cpf.getWriter("region_" + std::to_string(region->getId()));
    region->WriteToCpt(w);
  }
}

}  // namespace xtp
}  // namespace votca
