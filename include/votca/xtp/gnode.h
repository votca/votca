/*
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_KMC_GNODE_H
#define _VOTCA_KMC_GNODE_H

#include <vector>
#include <votca/ctp/qmpair.h>
#include <votca/ctp/segment.h>
#include <votca/tools/vec.h>
#include <votca/xtp/glink.h>
#include <votca/xtp/huffmantree.h>

using namespace std;

namespace votca {
namespace xtp {

class GNode {
 public:
  GNode()
      : occupied(false),
        occupationtime(0.0),
        escape_rate(0.0),
        hasdecay(false){};

  ~GNode(){};

  int id;
  bool occupied;
  bool injectable;
  double occupationtime;
  double escape_rate;
  bool hasdecay;
  tools::vec position;
  std::vector<GLink> events;
  // stuff for Coulomb interaction:
  double siteenergy;
  double reorg_intorig;  // UnCnN
  double reorg_intdest;  // UcNcC
  void AddEvent(int seg2, double rate12, tools::vec dr, double Jeff2,
                double reorg_out);
  const double& getEscapeRate() { return escape_rate; }
  void InitEscapeRate();
  void AddDecayEvent(double decayrate);
  void ReadfromSegment(ctp::Segment* seg, int carriertype);
  void AddEventfromQmPair(ctp::QMPair* pair, int carriertype);

  GLink* findHoppingDestination(double p);
  void MakeHuffTree();

 private:
  huffmanTree<GLink> hTree;
  void organizeProbabilities(int id, double add);
  void moveProbabilities(int id);
};

}  // namespace xtp
}  // namespace votca

#endif /* _VOTCA_KMC_GNODE_H */
