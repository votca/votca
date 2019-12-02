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

#include <numeric>
#include <votca/xtp/checkpoint.h>
#include <votca/xtp/jobtopology.h>
#include <votca/xtp/polarregion.h>
#include <votca/xtp/qmregion.h>
#include <votca/xtp/segmentmapper.h>
#include <votca/xtp/staticregion.h>

namespace votca {
namespace xtp {

void JobTopology::SortRegionsDefbyId(
    std::vector<tools::Property*>& regions_def) const {
  std::sort(regions_def.begin(), regions_def.end(),
            [](const tools::Property* A, const tools::Property* B) {
              return (A->get("id").as<Index>()) < (B->get("id").as<Index>());
            });
}

void JobTopology::ModifyOptionsByJobFile(
    std::vector<tools::Property*>& regions_def) const {

  const tools::Property& jobinput = _job.getInput();
  std::vector<const tools::Property*> regions_def_job =
      jobinput.Select("regions.region");

  std::string tag = "jobfile";
  for (tools::Property* prop : regions_def) {
    Index id = prop->get("id").as<Index>();
    std::vector<std::string> paths = FindReplacePathsInOptions(*prop, tag);
    if (!paths.empty()) {
      XTP_LOG(Log::info, _log) << " Region " << std::to_string(id)
                               << " is modified by jobfile" << std::flush;
      XTP_LOG(Log::info, _log)
          << " Replacing the following paths with jobfile entries"
          << std::flush;
      for (const std::string& path : paths) {
        XTP_LOG(Log::info, _log) << " - " << path << std::flush;
      }

      bool found_region_in_jobfile = false;
      const tools::Property* job_prop = nullptr;
      for (const tools::Property* prop_job : regions_def_job) {
        Index id2 = prop_job->get("id").as<Index>();
        if (id2 == id) {
          job_prop = prop_job;
          found_region_in_jobfile = true;
        }
      }
      if (!found_region_in_jobfile) {
        throw std::runtime_error("Region " + std::to_string(id) +
                                 " was not found in jobfile.");
      }
      UpdateFromJobfile(*prop, *job_prop, paths);
    }
  }
}

void JobTopology::BuildRegions(const Topology& top, tools::Property options) {

  std::vector<tools::Property*> regions_def = options.Select("region");
  CheckEnumerationOfRegions(regions_def);
  SortRegionsDefbyId(regions_def);
  ModifyOptionsByJobFile(regions_def);

  std::vector<std::vector<SegId> > region_seg_ids =
      PartitionRegions(regions_def, top);

  // around this point the whole jobtopology will be centered
  CreateRegions(options, top, region_seg_ids);
  XTP_LOG(Log::error, _log) << " Regions created" << std::flush;
  for (const auto& region : _regions) {
    XTP_LOG(Log::error, _log) << *region << std::flush;
  }

  return;
}

std::vector<std::string> JobTopology::FindReplacePathsInOptions(
    const tools::Property& options, std::string tag) const {
  std::vector<std::string> result;
  std::string options_path = "options.qmmm.regions.region";
  for (const tools::Property& sub : options) {
    if (sub.HasChildren()) {
      std::vector<std::string> subresult = FindReplacePathsInOptions(sub, tag);
      result.insert(result.end(), subresult.begin(), subresult.end());
    } else if (sub.value() == tag) {
      std::string path = sub.path() + "." + sub.name();
      std::size_t pos = path.find(options_path);
      if (pos != std::string::npos) {
        path.replace(pos, options_path.size(), "");
      }
      result.push_back(path);
    }
  }
  return result;
}

void JobTopology::UpdateFromJobfile(
    tools::Property& options, const tools::Property& job_opt,
    const std::vector<std::string>& paths) const {
  for (const std::string& path : paths) {
    if (job_opt.exists(path)) {
      options.set(path, job_opt.get(path).value());
    } else {
      throw std::runtime_error("Jobfile does not contain options for " + path);
    }
  }
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
    const tools::Property& options, const Topology& top,
    const std::vector<std::vector<SegId> >& region_seg_ids) {
  std::string mapfile =
      options.ifExistsReturnElseThrowRuntimeError<std::string>("mapfile");
  std::vector<const tools::Property*> regions_def = options.Select("region");
  // around this point the whole jobtopology will be centered for removing pbc
  Eigen::Vector3d center = top.getSegment(region_seg_ids[0][0].Id()).getPos();

  for (const tools::Property* region_def : regions_def) {
    Index id = region_def->ifExistsReturnElseThrowRuntimeError<Index>("id");
    const std::vector<SegId>& seg_ids = region_seg_ids[id];
    std::string type =
        region_def->ifExistsReturnElseThrowRuntimeError<std::string>("type");
    std::unique_ptr<Region> region;
    QMRegion QMdummy(0, _log, "");
    StaticRegion Staticdummy(0, _log);
    PolarRegion Polardummy(0, _log);

    if (type == QMdummy.identify()) {
      std::unique_ptr<QMRegion> qmregion =
          std::unique_ptr<QMRegion>(new QMRegion(id, _log, _workdir));
      QMMapper qmmapper(_log);
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
          std::unique_ptr<PolarRegion>(new PolarRegion(id, _log));
      PolarMapper polmap(_log);
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
          std::unique_ptr<StaticRegion>(new StaticRegion(id, _log));
      StaticMapper staticmap(_log);
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
    region->Initialize(*region_def);
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

std::vector<std::vector<SegId> > JobTopology::PartitionRegions(
    const std::vector<tools::Property*>& regions_def,
    const Topology& top) const {

  std::vector<Index> explicitly_named_segs_per_region;
  std::vector<std::vector<SegId> > segids_per_region;
  std::vector<bool> processed_segments =
      std::vector<bool>(top.Segments().size(), false);
  for (const tools::Property* region_def : regions_def) {

    if (!region_def->exists("segments") && !region_def->exists("cutoff")) {
      throw std::runtime_error(
          "Region definition needs either segments or a cutoff to find "
          "segments");
    }
    std::vector<SegId> seg_ids;
    if (region_def->exists("segments")) {
      std::string seg_ids_string =
          region_def->get("segments").as<std::string>();
      tools::Tokenizer tok(seg_ids_string, " \n\t");
      for (const std::string& seg_id_string : tok.ToVector()) {
        seg_ids.push_back(SegId(seg_id_string));
      }
      for (const SegId& seg_id : seg_ids) {
        if (seg_id.Id() > Index(top.Segments().size() - 1)) {
          throw std::runtime_error("Segment id is not in topology");
        }
        processed_segments[seg_id.Id()] = true;
      }
    }
    explicitly_named_segs_per_region.push_back(int(seg_ids.size()));

    if (region_def->exists("cutoff")) {
      double cutoff = tools::conv::nm2bohr *
                      region_def->ifExistsReturnElseThrowRuntimeError<double>(
                          "cutoff.radius");

      std::string seg_geometry =
          region_def->ifExistsReturnElseReturnDefault<std::string>(
              "cutoff.geometry", "n");
      double min = top.getBox().diagonal().minCoeff();
      if (cutoff > 0.5 * min) {
        throw std::runtime_error(
            (boost::format("Cutoff is larger than half the box size. Maximum "
                           "allowed cutoff is %1$1.1f") %
             (tools::conv::bohr2nm * 0.5 * min))
                .str());
      }
      std::vector<SegId> center = seg_ids;
      if (region_def->exists("cutoff.region")) {
        Index id = region_def->get("cutoff.region").as<Index>();
        bool only_explicit = region_def->ifExistsReturnElseReturnDefault<bool>(
            "cutoff.relative_to_explicit_segs", false);
        if (id < Index(segids_per_region.size())) {
          center = segids_per_region[id];
          if (only_explicit) {
            Index no_of_segs = explicitly_named_segs_per_region[id];
            if (no_of_segs == 0) {
              throw std::runtime_error("Region with id '" + std::to_string(id) +
                                       "' does not have");
            }
            center.resize(no_of_segs, SegId(0, "n"));
            // need the second argument because resize can also increase
            // capacity of vector and then needs a constructor,
            // here we shrink, so should not happen
          }
        } else {
          throw std::runtime_error("Region with id '" + std::to_string(id) +
                                   "' used for cutoff does not exist");
        }
      }
      if (center.empty()) {
        throw std::runtime_error(
            "Region needs either a segment or another region to which to apply "
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
    const std::vector<tools::Property*>& regions_def) const {
  std::vector<Index> reg_ids;
  for (const tools::Property* region_def : regions_def) {
    reg_ids.push_back(
        region_def->ifExistsReturnElseThrowRuntimeError<Index>("id"));
  }

  std::vector<Index> v(reg_ids.size());
  std::iota(v.begin(), v.end(), 0);
  if (!std::is_permutation(reg_ids.begin(), reg_ids.end(), v.begin())) {
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
