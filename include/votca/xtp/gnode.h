/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_XTP_GNODE_H
#define VOTCA_XTP_GNODE_H

#include <vector>
#include <votca/xtp/glink.h>
#include <votca/xtp/huffmantree.h>
#include <votca/xtp/qmpair.h>
#include <votca/xtp/segment.h>

namespace votca {
namespace xtp {

class GNode {
 public:
  GNode(const Segment& seg, QMStateType carriertype, bool injectable)
      : _id(seg.getId()),
        _siteenergy(seg.getSiteEnergy(carriertype)),
        _position(seg.getPos()),
        _injectable(injectable){};

  bool isOccupied() const { return _occupied; }
  void setOccupation(bool occupied) { _occupied = occupied; }
  bool isInjectable() const { return _injectable; }
  bool canDecay() const { return _hasdecay; }
  const Eigen::Vector3d& getPos() const { return _position; }
  int getId() const { return _id; }
  void UpdateOccupationTime(double deltat) { _occupationtime += deltat; }

  const std::vector<GLink>& Events() const { return _events; }
  double OccupationTime() const { return _occupationtime; }
  void AddEvent(GNode* seg2, const Eigen::Vector3d& dr, double rate);
  double getEscapeRate() const { return _escape_rate; }
  void InitEscapeRate();
  void AddDecayEvent(double decayrate);
  void AddEventfromQmPair(const QMPair& pair, std::vector<GNode>& nodes,
                          double rate);
  double getSitenergy() const { return _siteenergy; }

  GLink* findHoppingDestination(double p) const;
  void MakeHuffTree();

 private:
  int _id = 0;
  bool _occupied = false;
  bool _injectable = true;
  double _occupationtime = 0.0;
  double _escape_rate = 0.0;
  bool _hasdecay = false;
  double _siteenergy;
  Eigen::Vector3d _position;
  std::vector<GLink> _events;

  huffmanTree<GLink> hTree;
  void organizeProbabilities(int id, double add);
  void moveProbabilities(int id);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_GNODE_H
