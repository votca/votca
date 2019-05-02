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

namespace votca {
namespace xtp {

void JobTopology::BuildRegions(const Topology& top, const Job& job,
                               const tools::Property& options) {
  std::string mapfile =
      options.ifExistsReturnElseThrowRuntimeError<std::string>("mapfile");
  std::vector<const tools::Property*> regions_def = options.Select("region");

  std::vector<int> reg_id;
  for (const tools::Property* region_def : regions_def) {
    reg_id.push_back(
        region_def->ifExistsReturnElseThrowRuntimeError<int>("id"));
  }
  CheckEnumerationOfRegions(reg_id);
}

void JobTopology::CheckEnumerationOfRegions(const std::vector<int>& reg_ids) {
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
  a(_job_id, "jobid");

  for (const auto& region : _regions) {
    CheckpointWriter w =
        cpf.getWriter("region_" + std::to_string(region->getId()));
    region->WriteToCpt(w);
  }
}

}  // namespace xtp
}  // namespace votca
