/*
 *            Copyright 2009-2020 The VOTCA Development Team
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

// Local VOTCA includes
#include "job.h"
#include "logger.h"
#include "region.h"
#include "topology.h"

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
      : job_(job), log_(log), workdir_(workdir){};
  void BuildRegions(const Topology& top,
                    std::pair<std::string, tools::Property> options);

  void WriteToHdf5(std::string filename) const;

  void ReadFromHdf5(std::string filename);

  void WriteToPdb(std::string filename) const;

  std::vector<std::unique_ptr<Region> >::iterator begin() {
    return regions_.begin();
  }
  std::vector<std::unique_ptr<Region> >::iterator end() {
    return regions_.end();
  }

  const std::vector<std::unique_ptr<Region> >& Regions() const {
    return regions_;
  }

  std::vector<std::unique_ptr<Region> >& Regions() { return regions_; }

  Index size() const { return Index(regions_.size()); }

  std::vector<std::unique_ptr<Region> >::const_iterator begin() const {
    return regions_.begin();
  }
  std::vector<std::unique_ptr<Region> >::const_iterator end() const {
    return regions_.end();
  }

 private:
  std::vector<std::vector<SegId> > PartitionRegions(
      const tools::Property& regions_def, const Topology& top) const;

  void CreateRegions(const std::pair<std::string, tools::Property>& options,
                     const Topology& top,
                     const std::vector<std::vector<SegId> >& region_seg_ids);

  void ModifyOptionsByJobFile(tools::Property& regions_def) const;

  template <class T>
  void ShiftPBC(const Topology& top, const Eigen::Vector3d& center,
                T& mol) const;

  void CheckEnumerationOfRegions(const tools::Property& regions_def) const;
  std::vector<const tools::Property*> SortRegionsDefbyId(
      const tools::Property& regions_def) const;

  Job& job_;
  Logger& log_;
  std::vector<std::unique_ptr<Region> > regions_;
  std::string workdir_ = "";

  static constexpr int jobtopology_version() { return 1; }
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOBTOPOLOGY_H
