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
#include <queue>
#include <votca/xtp/gnode.h>

using namespace std;

namespace votca {
namespace xtp {
void GNode::AddDecayEvent(double decayrate) {
  _events.push_back(GLink(decayrate));
  _hasdecay = true;
}

void GNode::AddEvent(GNode* seg2, const Eigen::Vector3d& dr, double rate) {
  _events.push_back(GLink(seg2, rate, dr));
}

void GNode::InitEscapeRate() {
  _escape_rate = 0.0;
  for (const auto& event : _events) {
    _escape_rate += event.getRate();
  }
}

GLink* GNode::findHoppingDestination(double p) const {
  return hTree.findHoppingDestination(p);
}

void GNode::MakeHuffTree() {
  hTree.setEvents(&_events);
  hTree.makeTree();
}

void GNode::AddEventfromQmPair(const QMPair& pair, std::vector<GNode>& nodes,
                               double rate) {
  Index destination = 0;
  Eigen::Vector3d dr = Eigen::Vector3d::Zero();
  if (_id == pair.Seg1()->getId()) {
    destination = pair.Seg2()->getId();
    dr = pair.R();
  } else {
    destination = pair.Seg1()->getId();
    dr = -pair.R();
  }

  AddEvent(&nodes[destination], dr, rate);

  return;
}

}  // namespace xtp
}  // namespace votca
