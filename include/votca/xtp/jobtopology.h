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

#pragma once
#ifndef VOTCA_XTP_JOBTOPOLOGY_H
#define VOTCA_XTP_JOBTOPOLOGY_H

#include <votca/xtp/job.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/region.h>
#include <votca/xtp/topology.h>
/**
 * \brief Class to set up the topology, e.g division of molecules into different
 * regions for a specific job.
 *
 * How energies are evaluated depends critically on the id of the region.
 * a) Lower ids means being evaluated later. So expensive regions should have
 * low ids, as they can then already incorporate partially converged results
 * from higher id regions b) The energy of a region, includes all interactions
 * with regions of higher ids. i.e. E(region0)=E0+E01+E02, whereas
 * E(region1)=E1+E12 and not E10 The reason is that DFT codes only return E0+
 * all interaction energies and not the individual terms, this has to be changed
 * once we want to have multiple QM regions.
 */

namespace votca {
namespace xtp {
class SegId;
class JobTopology {
 public:
  JobTopology(Job& job, Logger& log, std::string workdir)
      : _job(job), _log(log), _workdir(workdir){};
  void BuildRegions(const Topology& top, tools::Property options);

  void WriteToHdf5(std::string filename) const;

  void WriteToPdb(std::string filename) const;

  std::vector<std::unique_ptr<Region> >::iterator begin() {
    return _regions.begin();
  }
  std::vector<std::unique_ptr<Region> >::iterator end() {
    return _regions.end();
  }

  const std::vector<std::unique_ptr<Region> >& Regions() const {
    return _regions;
  }

  std::vector<std::unique_ptr<Region> >& Regions() { return _regions; }

  int size() const { return _regions.size(); }

  std::vector<std::unique_ptr<Region> >::const_iterator begin() const {
    return _regions.begin();
  }
  std::vector<std::unique_ptr<Region> >::const_iterator end() const {
    return _regions.end();
  }

 private:
  std::vector<std::vector<SegId> > PartitionRegions(
      const std::vector<tools::Property*>& regions_def,
      const Topology& top) const;

  void CreateRegions(const tools::Property& options, const Topology& top,
                     const std::vector<std::vector<SegId> >& region_seg_ids);

  void UpdateFromJobfile(tools::Property& options,
                         const tools::Property& job_opt,
                         const std::vector<std::string>& paths) const;
  std::vector<std::string> FindReplacePathsInOptions(
      const tools::Property& options, std::string tag) const;
  void ModifyOptionsByJobFile(std::vector<tools::Property*>& regions_def) const;
  void UpdateFromJobfile(tools::Property& options,
                         const tools::Property& job_opt,
                         const std::string& tag) const;

  template <class T>
  void ShiftPBC(const Topology& top, const Eigen::Vector3d& center,
                T& mol) const;

  void CheckEnumerationOfRegions(
      const std::vector<tools::Property*>& regions_def) const;
  void SortRegionsDefbyId(std::vector<tools::Property*>& regions_def) const;

  Job& _job;
  Logger& _log;
  std::vector<std::unique_ptr<Region> > _regions;
  std::string _workdir = "";
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOBTOPOLOGY_H
