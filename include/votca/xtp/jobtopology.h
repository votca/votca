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

#ifndef VOTCA_XTP_JOBTOPOLOGY_H
#define VOTCA_XTP_JOBTOPOLOGY_H
#include <votca/tools/elements.h>
#include <votca/xtp/checkpoint.h>
#include <votca/xtp/region.h>
#include <votca/xtp/topology.h>

/**
 * \brief Class to set up the toplogy, e.g division of molecules into different
 * regions for a specific job.
 */

namespace votca {
namespace xtp {

class JobTopology {
 public:
  void BuildRegions(Topology& top);
  void WriteToHdf5(std::string filename) const;

 private:
  std::vector<std::unique_ptr<Region> > _regions;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOBTOPOLOGY_H
