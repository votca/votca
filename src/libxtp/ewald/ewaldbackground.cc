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

#include "ewaldbackground.h"
#include "kspace.h"
#include "rspace.h"
#include "votca/xtp/tictoc.h"
#include <fstream>
#include <iostream>
#include <vector>
#include "votca/xtp/polarregion.h"

namespace votca {
namespace xtp {

void EwaldBackground::ImportBackgroundFromHdf5(std::string filename) {
  CheckpointFile cpf(filename, CheckpointAccessLevel::READ);
  CheckpointReader r = cpf.getReader("polar_background");
  for (Index i = 0; i < r.getNumDataSets(); ++i) {
    CheckpointReader rr = r.openChild("background_" + std::to_string(i));
    background.push_back(EwdSegment(rr));
  }
}

void EwaldBackground::ApplyBackgroundFields(
    std::vector<std::unique_ptr<votca::xtp::Region>>& regions,
    const std::vector<std::vector<SegId>>& region_seg_ids) {
  // Sanity check
  if (region_seg_ids.size() > 2) {
    throw std::runtime_error(
        "You requested a calculation with more than two inner regions. This is "
        "not allowed, the ewald background can only be used with a polar inner "
        "region or with a qm region inside a polar region.");
  }
  // Because we implement it step wise
  if (region_seg_ids.size() > 1) {
    throw std::runtime_error(
        "The qm region inside a polar region inside a ewald background is not "
        "yet implemented.");
  }

  // Apply background fields to sites in the polarization cloud
  if (region_seg_ids.size() == 1) {  // i.e. only a polariation cloud
    PolarRegion* pCloud = dynamic_cast<PolarRegion*>(regions[0].get());
    std::vector<SegId> pCloud_original_ids = region_seg_ids[0];
    for (Index i = 0; i < pCloud->size(); i++){
      // (*pCloud)[i] will be the ith segment in pCloud
      bgFieldAtSegment((*pCloud)[i], pCloud_original_ids);
    }
  }

}

void EwaldBackground::bgFieldAtSegment(PolarSegment& seg, std::vector<SegId> pCloud_indices){
std::cout << " My ID is " << seg.getId() << std::endl;
}

}  // namespace xtp
}  // namespace votca