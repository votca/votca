/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_GLINK_H
#define VOTCA_XTP_GLINK_H

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

class GNode;
class GLink {

 public:
  GLink(GNode* dest, double rate, const Eigen::Vector3d& dr)
      : destination(dest), rate_(rate), dr_(dr){};

  GLink(double rate) : rate_(rate), decayevent_(true){};

  double getValue() const { return rate_; }
  double getRate() const { return rate_; }
  GNode* getDestination() const {
    assert(!decayevent_ && "Decay event has no destination");
    return destination;
  }
  const Eigen::Vector3d& getDeltaR() const { return dr_; }
  bool isDecayEvent() const { return decayevent_; }

 private:
  GNode* destination = nullptr;
  double rate_ = 0.0;
  Eigen::Vector3d dr_ = Eigen::Vector3d::Zero();
  bool decayevent_ = false;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GLINK_H
