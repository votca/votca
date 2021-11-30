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
      BGSite bgsite(psite);
      sites_.push_back(bgsite);
    }
    id_ = pol.getId();
    position_ = pol.getPos();
  }

  Index getId() const { return id_;}
  const Eigen::Vector3d& getPos() const { return position_;}

  void calcPos() {
    tools::Elements element;
    Eigen::Vector3d pos = Eigen::Vector3d::Zero();
    double totalmass = 0.0;
    for (const auto& site : sites_) {
      double mass = element.getMass(site.getElement());
      totalmass += mass;
      pos += mass * site.getPos();
    }
    position_ = pos / totalmass;
  }

  ~BGSegment() = default;

  const BGSite& at(Index index) const { return sites_.at(index); }
  BGSite& at(Index index) { return sites_.at(index); }

  const BGSite& operator[](Index index) const { return sites_[index]; }
  BGSite& operator[](Index index) { return sites_[index]; }

  typename std::vector<BGSite>::iterator begin() { return sites_.begin(); }
  typename std::vector<BGSite>::iterator end() { return sites_.end(); }

  typename std::vector<BGSite>::const_iterator begin() const {
    return sites_.begin();
  }
  typename std::vector<BGSite>::const_iterator end() const {
    return sites_.end();
  }

 private:
  Index id_;
  std::vector<BGSite> sites_;
  Eigen::Vector3d position_;
};
}  // namespace xtp
}  // namespace votca
#endif