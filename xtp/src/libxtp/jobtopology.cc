/*
 *            Copyright 2009-2021 The VOTCA Development Team
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

// Standard includes
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>

// Local VOTCA includes
#include "votca/tools/property.h"
#include "votca/tools/version.h"
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/jobtopology.h"
#include "votca/xtp/polarregion.h"
#include "votca/xtp/qmregion.h"
#include "votca/xtp/segmentmapper.h"
#include "votca/xtp/staticregion.h"

namespace votca {
namespace xtp {

std::vector<const tools::Property*> JobTopology::SortRegionsDefbyId(
    const tools::Property& regions_def) const {
  std::vector<const tools::Property*> view;
  for (const auto& child : regions_def) {
    view.push_back(&child);
  }
  std::sort(view.begin(), view.end(),
            [](const tools::Property* A, const tools::Property* B) {
              return (A->get("id").as<Index>()) < (B->get("id").as<Index>());
            });
  return view;
}

std::vector<std::string> ModifyRegionOptionsFromJobFileRegion(
    tools::Property& opt, const tools::Property& job_opt) {
  std::vector<std::string> modified_options;
  for (auto& child : opt) {
    if (child.as<std::string>() == "jobfile") {
      if (job_opt.exists(child.name())) {
        const tools::Property& job_opt_child = job_opt.get(child.name());
        modified_options.push_back(child.path() + "." + child.name());
        if (job_opt_child.HasChildren()) {
          child.value() = "";
          for (const auto& job_childchild : job_opt_child) {
            child.add(job_childchild);
          }
        } else {
          child.value() = job_opt_child.value();
        }
      } else {
        throw std::runtime_error("Requested to replace:" + child.path() + "." +
                                 child.name() +
                                 " but no option in jobfile found");
      }
    } else if (job_opt.exists(child.name())) {
      std::vector<std::string> child_modified_options =
          ModifyRegionOptionsFromJobFileRegion(child,
                                               job_opt.get(child.name()));
      modified_options.insert(modified_options.end(),
                              child_modified_options.begin(),
                              child_modified_options.end());
    }
  }
  return modified_options;
}

void JobTopology::ModifyOptionsByJobFile(tools::Property& regions_def) const {

  const tools::Property& prop = job_.getInput();
  std::vector<const tools::Property*> regions_def_job =
      prop.Select("regions.*region");

  for (const auto& job_region : regions_def_job) {

    Index job_region_id = job_region->get("id").as<Index>();

    auto found_prop =
        std::find_if(regions_def.begin(), regions_def.end(),
                     [job_region_id](const tools::Property& p) {
                       return p.get("id").as<Index>() == job_region_id;
                     });
    if (found_prop == regions_def.end()) {
      throw std::runtime_error(
          "Your jobfile options want to modify a region with id " +
          std::to_string(job_region_id) +
          " it does not exist in your region specifications.");
    }

    tools::Property& region = *found_prop;
    if (region.name() != job_region->name()) {
      throw std::runtime_error("Types of region with id" +
                               std::to_string(job_region_id) +
                               " do not agree jobfile:" + job_region->name() +
                               " options:" + region.name());
    }

    std::vector<std::string> changed_options =
        ModifyRegionOptionsFromJobFileRegion(region, *job_region);

    if (!changed_options.empty()) {
      XTP_LOG(Log::info, log_) << " Region " << job_region_id
                               << " is modified by jobfile" << std::flush;
      XTP_LOG(Log::info, log_)
          << " Replacing the following paths with jobfile entries"
          << std::flush;
      for (const std::string& path : changed_options) {
        XTP_LOG(Log::info, log_) << " - " << path << std::flush;
      }
    }
  }
}

void JobTopology::BuildRegions(
    const Topology& top, std::pair<std::string, tools::Property> options) {

  CheckEnumerationOfRegions(options.second);
  ModifyOptionsByJobFile(options.second);

  std::vector<std::vector<SegId>> region_seg_ids =
      PartitionRegions(options.second, top);

  // // around this point the whole jobtopology will be centered
  CreateRegions(options, top, region_seg_ids);
  XTP_LOG(Log::error, log_) << " Regions created" << std::flush;
  for (const auto& region : regions_) {
    XTP_LOG(Log::error, log_) << *region << std::flush;
  }

  return;
}

template <class T>
void JobTopology::ShiftPBC(const Topology& top, const Eigen::Vector3d& center,
                           T& mol) const {
  Eigen::Vector3d r_pbc = top.PbShortestConnect(center, mol.getPos());
  Eigen::Vector3d r = mol.getPos() - center;
  Eigen::Vector3d shift = r_pbc - r;
  if (shift.norm() > 1e-9) {
    mol.Translate(shift);
  }
}

void JobTopology::CreateRegions(
    const std::pair<std::string, tools::Property>& options, const Topology& top,
    const std::vector<std::vector<SegId>>& region_seg_ids) {
  std::string mapfile = options.first;
  // around this point the whole jobtopology will be centered for removing pbc
  Eigen::Vector3d center = top.getSegment(region_seg_ids[0][0].Id()).getPos();

  for (const tools::Property& region_def : options.second) {
    Index id = region_def.get("id").as<Index>();
    const std::vector<SegId>& seg_ids = region_seg_ids[id];
    std::string type = region_def.name();
    std::unique_ptr<Region> region;
    QMRegion QMdummy(0, log_, "");
    StaticRegion Staticdummy(0, log_);
    PolarRegion Polardummy(0, log_);
    if (type == QMdummy.identify()) {
      std::unique_ptr<QMRegion> qmregion =
          std::make_unique<QMRegion>(id, log_, workdir_);
      QMMapper qmmapper(log_);
      qmmapper.LoadMappingFile(mapfile);
      for (const SegId& seg_index : seg_ids) {
        const Segment& segment = top.getSegment(seg_index.Id());
        QMMolecule mol = qmmapper.map(segment, seg_index);
        mol.setType("qm" + std::to_string(id));
        ShiftPBC(top, center, mol);
        qmregion->push_back(mol);
      }
      region = std::move(qmregion);
    } else if (type == Polardummy.identify()) {
      std::unique_ptr<PolarRegion> polarregion =
          std::make_unique<PolarRegion>(id, log_);
      PolarMapper polmap(log_);
      polmap.LoadMappingFile(mapfile);
      for (const SegId& seg_index : seg_ids) {
        const Segment& segment = top.getSegment(seg_index.Id());

        PolarSegment mol = polmap.map(segment, seg_index);

        ShiftPBC(top, center, mol);
        mol.setType("mm" + std::to_string(id));
        polarregion->push_back(mol);
      }
      region = std::move(polarregion);
    } else if (type == Staticdummy.identify()) {
      std::unique_ptr<StaticRegion> staticregion =
          std::make_unique<StaticRegion>(id, log_);
      StaticMapper staticmap(log_);
      staticmap.LoadMappingFile(mapfile);
      for (const SegId& seg_index : seg_ids) {
        const Segment& segment = top.getSegment(seg_index.Id());
        StaticSegment mol = staticmap.map(segment, seg_index);
        mol.setType("mm" + std::to_string(id));
        ShiftPBC(top, center, mol);
        staticregion->push_back(mol);
      }
      region = std::move(staticregion);

    } else {
      throw std::runtime_error("Region type not known!");
    }
    region->Initialize(region_def);
    regions_.push_back(std::move(region));
  }
}

void JobTopology::WriteToPdb(std::string filename) const {

  csg::PDBWriter writer;
  writer.Open(filename, false);
  writer.WriteHeader("Job:" + std::to_string(this->job_.getId()));
  for (const std::unique_ptr<Region>& reg : regions_) {
    reg->WritePDB(writer);
  }
  writer.Close();
}

std::vector<std::vector<SegId>> JobTopology::PartitionRegions(
    const tools::Property& regions_def, const Topology& top) const {

  std::vector<const tools::Property*> sorted_regions =
      SortRegionsDefbyId(regions_def);

  std::vector<Index> explicitly_named_segs_per_region;
  std::vector<std::vector<SegId>> segids_per_region;
  std::vector<bool> processed_segments =
      std::vector<bool>(top.Segments().size(), false);
  for (const tools::Property* region_def : sorted_regions) {

    if (!region_def->exists("segments") && !region_def->exists("cutoff")) {
      throw std::runtime_error(
          "Region definition needs either segments or a cutoff to find "
          "segments");
    }
    std::vector<SegId> seg_ids;
    if (region_def->exists("segments")) {
      seg_ids = region_def->get("segments").as<std::vector<SegId>>();
      for (const SegId& seg_id : seg_ids) {
        if (seg_id.Id() > Index(top.Segments().size() - 1)) {
          throw std::runtime_error("Segment id is not in topology");
        }
        processed_segments[seg_id.Id()] = true;
      }
    }
    explicitly_named_segs_per_region.push_back(Index(seg_ids.size()));

    if (region_def->exists("cutoff")) {
      double cutoff =
          tools::conv::nm2bohr * region_def->get("cutoff.radius").as<double>();

      std::string seg_geometry =
          region_def->get("cutoff.geometry").as<std::string>();
      double min = top.getBox().diagonal().minCoeff();
      if (cutoff > 0.5 * min) {
        throw std::runtime_error(
            (boost::format("Cutoff is larger than half the box size. Maximum "
                           "allowed cutoff is %1$1.1f") %
             (tools::conv::bohr2nm * 0.5 * min))
                .str());
      }
      std::vector<SegId> center = seg_ids;

      Index id = region_def->get("cutoff.region").as<Index>();
      bool only_explicit = region_def->get("cutoff.explicit_segs").as<bool>();
      if (id < Index(segids_per_region.size())) {
        center = segids_per_region[id];
        if (only_explicit) {
          Index no_of_segs = explicitly_named_segs_per_region[id];
          if (no_of_segs == 0) {
            throw std::runtime_error(
                "Region with id '" + std::to_string(id) +
                "' does not have explicitly named segments");
          }
          center.resize(no_of_segs, SegId(0, "n"));
          // need the second argument because resize can also increase
          // capacity of vector and then needs a constructor,
          // here we shrink, so should not happen
        }
      } else if (id != region_def->get("id").as<Index>()) {
        throw std::runtime_error("Region with id '" + std::to_string(id) +
                                 "' used for cutoff does not exist");
      }

      if (center.empty()) {
        throw std::runtime_error(
            "Region needs either a segment or another region to which to "
            "apply "
            "the cutoff");
      }
      for (const SegId& segid : center) {
        const Segment& center_seg = top.getSegment(segid.Id());
        for (const Segment& otherseg : top.Segments()) {
          if (center_seg.getId() == otherseg.getId() ||
              processed_segments[otherseg.getId()]) {
            continue;
          }
          if (top.PbShortestConnect(center_seg.getPos(), otherseg.getPos())
                  .norm() < cutoff) {
            seg_ids.push_back(SegId(otherseg.getId(), seg_geometry));
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
    const tools::Property& regions_def) const {
  std::vector<Index> reg_ids;
  for (const tools::Property& region_def : regions_def) {
    reg_ids.push_back(region_def.get("id").as<Index>());
  }

  std::vector<Index> v(reg_ids.size());
  std::iota(v.begin(), v.end(), 0);
  if (!std::is_permutation(reg_ids.begin(), reg_ids.end(), v.begin())) {
    throw std::runtime_error(
        "Region id definitions are not clear. You must start at id 0 and "
        "then "
        "ascending order. i.e. 0 1 2 3.");
  }
}

void JobTopology::WriteToHdf5(std::string filename) const {
  CheckpointFile cpf(filename, CheckpointAccessLevel::CREATE);
  CheckpointWriter a = cpf.getWriter();
  a(job_.getId(), "jobid");
  a(votca::tools::ToolsVersionStr(), "XTPVersion");
  a(jobtopology_version(), "version");
  for (const auto& region : regions_) {
    CheckpointWriter w =
        cpf.getWriter("region_" + std::to_string(region->getId()));
    region->WriteToCpt(w);
  }
}

void JobTopology::ReadFromHdf5(std::string filename) {
  CheckpointFile cpf(filename, CheckpointAccessLevel::READ);
  CheckpointReader a = cpf.getReader();
  Index id = -1;
  a(id, "jobid");
  if (id != job_.getId()) {
    throw std::runtime_error(
        "Jobid from checkpoint file does not agree with jobid" +
        std::to_string(id) + ":" + std::to_string(job_.getId()));
  }
  regions_.clear();
  Index no_regions = a.getNumDataSets();
  regions_.reserve(no_regions);

  QMRegion QMdummy(0, log_, "");
  StaticRegion Staticdummy(0, log_);
  PolarRegion Polardummy(0, log_);

  for (Index i = 0; i < no_regions; i++) {
    CheckpointReader r = cpf.getReader("region_" + std::to_string(i));
    std::string type;
    r(type, "type");

    std::unique_ptr<Region> region = nullptr;
    if (type == QMdummy.identify()) {
      region = std::make_unique<QMRegion>(i, log_, workdir_);
    } else if (type == Polardummy.identify()) {
      region = std::make_unique<PolarRegion>(i, log_);
    } else if (type == Staticdummy.identify()) {
      region = std::make_unique<StaticRegion>(i, log_);
    } else {
      throw std::runtime_error("Region type:" + type + " not known!");
    }
    region->ReadFromCpt(r);
    if (id != region->getId()) {
      throw std::runtime_error("read in the wrong region, ids are mismatched.");
    }
    regions_.push_back(std::move(region));
  }
}

}  // namespace xtp
}  // namespace votca
