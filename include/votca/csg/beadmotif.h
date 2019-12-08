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

#pragma once
#ifndef VOTCA_CSG_BEADMOTIF_H
#define VOTCA_CSG_BEADMOTIF_H

#include "basebead.h"
#include "beadstructure.h"
#include <votca/tools/reducedgraph.h>

namespace TOOLS = votca::tools;

namespace votca {
namespace csg {

/**
 * \brief Designed determine what kind of structure a beadstructure has
 *
 * Wants a beadstructure is created there are several different forms it can
 * have, this class helps to classify and break structures up into the
 * appropriate sub class. The possible classes include:
 *
 * Some of the structures are considered fundamental and cannot be broken up any
 * more than they are. These elements are listed as simple
 *
 * 1. Single bead          // Simple
 * 2. Line                 // Simple
 * 3. Loop                 // Simple
 * 4. Fused Ring           // Simple
 * 5. Single Structure     // Complex
 * 6. Multiple Structures  // Complex
 * 6. Undefined            // Undefined
 *
 * Examples of each type are shown below connections are indicated with lines
 *
 * Single Bead:
 *
 * H1
 *
 * Line:
 *
 * C1 - C2 - C3                   H1 - H2                   H1  H2
 *                                                           \  /
 *                                                            O1
 * Loop:
 *
 * C1 - C2                        C5 - C6
 * |    |                        /       \
 * C4 - C3                      C7       C8
 *                                \      /
 *                                C9 - C10
 *
 * Fused Ring:
 *
 *    C5 - C6
 *   /       \
 *  C7       C8 - C11
 *    \      /      \
 *    C9 - C10      C12
 *           \      /
 *           C13 - C14
 *
 * Single Structures:
 *
 *      H1                         H5      H6
 *      |                           \     /
 * H2 - C1 - H3                     C1 - C2
 *      |                           /      \
 *      H4                    H7 - C3      C4 - H8
 *                                   \     /
 *                                   C5 - C6
 *                                   /      \
 *                                  H9      H10
 *
 * Multiple Structures:
 *
 * Any combination of more than one structure
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

class BeadMotif : public BeadStructure<BaseBead> {
 public:
  enum MotifType {
    empty,
    single_bead,
    line,
    loop,
    fused_ring,
    single_structure,
    multiple_structures,
    undefined
  };

  BeadMotif() = default;

  BeadMotif(const BeadStructure &structure) : BeadStructure(structure){};

  /// Gets the motif type, calculates it first if it is not yet known
  MotifType getType();

  /**
   * \brief Determines if the motif type is a simple type
   *
   * If the motif is of type single_bead, line, loop or fused ring it will
   * return true.
   *
   * @return true if the bead motif type is simple
   **/
  bool isMotifSimple();

  /**
   * \brief Adds a new bead to the motif
   *
   * This method calls the beastructure AddBead method but it also switches an
   * attribute indicating that the beadtype is now out of date.
   *
   * @param[in] basebead pointer
   **/
  void AddBead(BaseBead *basebead);

  /**
   * \brief Adds a new connection to the motif
   *
   * Also switches an internal attribute to indicate the beadtype is no longer
   * up to date.
   *
   * \param[in] - id of the first and second beads that are connected
   **/
  void ConnectBeads(Index bead1_id, Index bead2_id);

 private:
  MotifType type_ = MotifType::undefined;
  bool junctionsUpToDate_ = false;
  bool type_up_to_date_ = false;

  std::vector<Index> junctions_;
  TOOLS::ReducedGraph reduced_graph_;

  void InitializeGraph_();
  bool junctionExist_();
  void CalculateType_();
  bool isSingle_();
  bool isLine_();
  bool isLoop_();
  bool isFusedRing_();
};
}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMOTIF_H
