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
 * author: Kordt
 */

#pragma once
#ifndef VOTCA_XTP_CHARGECARRIER_H
#define VOTCA_XTP_CHARGECARRIER_H

// Local VOTCA includes
#include "glink.h"
#include "gnode.h"

namespace votca {
namespace xtp {

class Chargecarrier {
 public:
  Chargecarrier(Index id)
      : id_(id),
        lifetime(0.0),
        steps(0),
        dr_travelled_(Eigen::Vector3d::Zero()),
        node(nullptr) {};
  bool hasNode() { return (node != nullptr); }
  void updateLifetime(double dt) { lifetime += dt; }
  void updateOccupationtime(double dt) { node->UpdateOccupationTime(dt); }
  void updateSteps(Index t) { steps += t; }
  void resetCarrier() {
    lifetime = 0;
    steps = 0;
    dr_travelled_ = Eigen::Vector3d::Zero();
  }
  double getLifetime() const { return lifetime; }
  Index getSteps() const { return steps; }
  Index getCurrentNodeId() const { return node->getId(); }
  double getCurrentEnergy() const { return node->getSitenergy(); }
  const Eigen::Vector3d& getCurrentPosition() const { return node->getPos(); }
  double getCurrentEscapeRate() const { return node->getEscapeRate(); }
  GNode& getCurrentNode() const { return *node; }

  void ReleaseNode() { node->setOccupation(false); }

  void settoNote(GNode* newnode) {
    node = newnode;
    node->setOccupation(true);
  }

  void jumpAccordingEvent(const GLink& event) {
    ReleaseNode();
    settoNote(event.getDestination());
    dr_travelled_ += event.getDeltaR();
  }

  const Eigen::Vector3d& get_dRtravelled() const { return dr_travelled_; }

  Index getId() const { return id_; }
  void setId(Index id) { id_ = id; }

 private:
  Index id_;
  double lifetime;
  Index steps;
  Eigen::Vector3d dr_travelled_;
  GNode* node;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CHARGECARRIER_H
