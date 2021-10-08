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
#ifndef VOTCA_XTP_EWALDSEGMENT_H
#define VOTCA_XTP_EWALDSEGMENT_H

#include <vector>

// Local VOTCA includes
#include "votca/xtp/classicalsegment.h"

// Private VOTCA includes
#include "votca/xtp/ewd_site.h"

namespace votca {
namespace xtp {

class EwdSegment {
 public:
  EwdSegment(const PolarSegment& pol) {
    for (const PolarSite& psite : pol) {
      EwdSite esite(psite);
      sites_.push_back(esite);
    }
    id_ = pol.getId();
    position_ = pol.getPos();
  }

  EwdSegment(CheckpointReader&r, Index id){
    CptTable table = r.openTable<EwdSite>("background_sites");
    sites_.clear();
    sites_.reserve(table.numRows());
    std::vector<typename EwdSite::data> dataVec(table.numRows());
    table.read(dataVec);
    for (std::size_t i = 0; i < table.numRows(); ++i) {
      sites_.push_back(EwdSite(dataVec[i]));
    }
    id_ = id;
    r(position_, "pos");
  };

  ~EwdSegment() = default;

  bool operator==(const EwdSegment& other){
    if (other.id_ != this->id_){
      return false;
    }
    if (other.position_ != this->position_){
      return false;
    }
    if (other.size() != this->size()){
      return false;
    } else {
      for(Index i = 0; i < other.size(); ++i){
        if(this->sites_[i] != other.sites_[i]){
          return false;
        }
      }
    }
    return true;
  }

  bool operator!=(const EwdSegment& other){
    return !operator==(other);
  }

  const Eigen::Vector3d& getPos() const { return position_; }

  Index getId() const {return id_;}

  const EwdSite& at(Index index) const { return sites_.at(index); }
  EwdSite& at(Index index) { return sites_.at(index); }

  const EwdSite& operator[](Index index) const { return sites_[index]; }
  EwdSite& operator[](Index index) { return sites_[index]; }

  typename std::vector<EwdSite>::iterator begin() { return sites_.begin(); }
  typename std::vector<EwdSite>::iterator end() { return sites_.end(); }

  typename std::vector<EwdSite>::const_iterator begin() const {
    return sites_.begin();
  }
  typename std::vector<EwdSite>::const_iterator end() const {
    return sites_.end();
  }

  Index size() const { return sites_.size();}

  void WriteToCpt(CheckpointWriter& w) {
    w(position_, "pos");
    CptTable table = w.openTable<EwdSite>("background_sites", sites_.size());
    std::vector<EwdSite::data> dataVec(sites_.size());
    for (std::size_t i = 0; i < sites_.size(); ++i) {
      sites_[i].WriteData(dataVec[i]);
    }
    table.write(dataVec);
  }


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

 private:
  Index id_;
  std::vector<EwdSite> sites_;
  Eigen::Vector3d position_;
};
}  // namespace xtp
}  // namespace votca
#endif