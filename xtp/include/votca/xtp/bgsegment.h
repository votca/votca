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
#ifndef VOTCA_XTP_BGSEGMENT_H
#define VOTCA_XTP_BGSEGMENT_H

#include <vector>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"

// Private VOTCA includes
#include "votca/xtp/bgsite.h"

namespace votca {
namespace xtp {

class BGSegment {
 public:
  BGSegment(const PolarSegment& pol) {
    for (const PolarSite& psite : pol) {
      BGSite esite(psite);
      sites_.push_back(esite);
    }
    id_ = pol.getId();
    position_ = pol.getPos();
  }

  ~BGSegment() = default;

 private:
  Index id_;
  std::vector<BGSite> sites_;
  Eigen::Vector3d position_;
};
}  // namespace xtp
}  // namespace votca
#endif