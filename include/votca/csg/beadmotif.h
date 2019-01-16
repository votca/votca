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

#ifndef VOTCA_CSG_BEADMOTIF_H
#define VOTCA_CSG_BEADMOTIF_H

#include <list>
#include <map>
#include <memory>
#include <set>
#include <string>

#include <votca/csg/basebead.h>
#include <votca/tools/reducedgraph.h>

#include "beadstructure.h"

namespace votca {
namespace csg {

/**
 * \brief Designed determine what kind of structure a beadstructure has
 *
 * Wants a beadstructure is created there are several different forms it can
 * have, this class helps to classify and break structures up into the
 * appropriate sub class. The possible classes include:
 *
 * 1. Single bead
 * 2. Line
 * 3. Loop
 * 4. Fused Ring
 * 5. Single Structure
 * 6. Multiple Structures
 * 6. Undefined
 *
 * The Single, line, loop and Fused Ring types are all elementary types that
 * represent a fundamental structural unit.
 *
 * Single Structure represents a type that has not been broken up into its
 * fundamental components but consists of a single interconnected structure.
 *
 * Multiple Structures means that the structure consists of multiple
 * structures
 *
 * Undefined means the structure has not yet been categorized.
 *
 * Though the Beadmotif inherits from the Beadstructure its methods are kept
 * private the reason is that the type might be incorrect if an  AddBead
 * from the original class is called.
 **/

class BeadMotif : private BeadStructure {
 public:
  BeadMotif()
      : BeadStructure(),
        type_(motif_type::undefined),
        junctionsUpToDate_(false){};
  ~BeadMotif(){};

  enum motif_type {
    empty,
    single_bead,
    line,
    loop,
    fused_ring,
    single_structure,
    multiple_structures,
    undefined
  };

  motif_type getType();

  BaseBead *getBead(int id);
  void ConnectBeads(int bead1_id, int bead2_id);

  std::vector<BaseBead *> getNeighBeads(int index);

  void AddBead(BaseBead *bead);

  int BeadCount();

  bool isStructureEquivalent(BeadMotif &beadmotif);

 private:
  void InitializeGraph_();
  motif_type type_;
  bool junctionsUpToDate_;
  bool type_up_to_date_ = false;
  std::vector<int> junctions_;
  std::unique_ptr<ReducedGraph> reduced_graph_;
  bool junctionExist_();
  bool isSingle_();
  bool isLoop_();
  /**
   * \brief Calculates the type of the motif
   **/
  void CalculateType_();

  // One has to explore the whole tree to from each of the junctions to
  // determine if the model is a fused ring or not. For speed it might
  // make since to reduce the graph first to junctions of 3 or more.
  //
  // if There is not way back to the junction than you have something
  // like this:
  //
  // c1 - c2 - c5 - c6
  // |    |    |    |
  // c3 - c4   c7 - c8
  //
  // Say you start at c2 and head down tree c5 if you never find a way back
  // you can split it
  //
  // If you have something like this
  //
  // c1 - c2 - c3
  // |  /   \  |
  //  c4     c5
  //
  //  Then you do not have a fused ring, must be represented as a joint
  //  and two lines. Exploring a tree will only lead to one way back.
  //
  //         c6
  //        /  |
  // c1 - c2 - c3
  // |  /   \  |
  //  c4     c5
  //
  //  Still acts like a joint, For it not to be a joint exploring a single
  //  branch originating from the junction should lead to exploration of
  //  all the edges.
  bool isFusedRing_();
  bool isLine_();
};
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMOTIF_H
