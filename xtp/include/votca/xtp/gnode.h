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
#ifndef VOTCA_XTP_GNODE_H
#define VOTCA_XTP_GNODE_H

// Local VOTCA includes
#include "glink.h"
#include "huffmantree.h"
#include "qmpair.h"
#include "segment.h"

namespace votca {
namespace xtp {

class GNode {
 public:
  GNode(const Segment& seg, QMStateType carriertype, bool injectable)
      : id_(seg.getId()),
        siteenergy_(seg.getSiteEnergy(carriertype)),
        position_(seg.getPos()),
        injectable_(injectable) {};

  bool isOccupied() const { return occupied_; }
  void setOccupation(bool occupied) { occupied_ = occupied; }
  bool isInjectable() const { return injectable_; }
  bool canDecay() const { return hasdecay_; }
  const Eigen::Vector3d& getPos() const { return position_; }
  Index getId() const { return id_; }
  void UpdateOccupationTime(double deltat) { occupationtime_ += deltat; }

  const std::vector<GLink>& Events() const { return events_; }
  double OccupationTime() const { return occupationtime_; }

  double getEscapeRate() const { return escape_rate_; }
  void InitEscapeRate();
  void AddDecayEvent(double decayrate);
  void AddEventfromQmPair(const QMPair& pair, std::vector<GNode>& nodes,
                          double rate);
  double getSitenergy() const { return siteenergy_; }

  GLink* findHoppingDestination(double p) const;
  void MakeHuffTree();
  void AddEvent(GNode* seg2, const Eigen::Vector3d& dr, double rate);

 private:
  Index id_ = 0;
  bool occupied_ = false;
  double occupationtime_ = 0.0;
  double escape_rate_ = 0.0;
  bool hasdecay_ = false;
  double siteenergy_;
  Eigen::Vector3d position_;
  bool injectable_ = true;
  std::vector<GLink> events_;

  huffmanTree<GLink> hTree;

  void organizeProbabilities(Index id, double add);
  void moveProbabilities(Index id);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GNODE_H
