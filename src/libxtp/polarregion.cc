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

#include <votca/xtp/polarregion.h>

namespace votca {
namespace xtp {

void PolarRegion::Initialize(const tools::Property& prop) { return; }

bool PolarRegion::Converged() const { return false; }

void PolarRegion::Evaluate() { return; }

void PolarRegion::ApplyInfluenceOfOtherRegions(
    const std::vector<std::unique_ptr<Region> >& regions) {
  return;
}

}  // namespace xtp
}  // namespace votca
