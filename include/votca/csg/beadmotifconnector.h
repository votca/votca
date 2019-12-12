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
#ifndef VOTCA_CSG_BEADMOTIFCONNECTOR_H
#define VOTCA_CSG_BEADMOTIFCONNECTOR_H

#include <unordered_set>
#include <vector>

#include <boost/bimap.hpp>
#include <boost/bimap/multiset_of.hpp>
#include <boost/bimap/set_of.hpp>

#include <votca/tools/edge.h>

namespace votca {
namespace csg {

/**
 * @brief Maps a reduced edge to all the edges that make up the reduced edge
 *
 * Reduced edge
 *
 * 0 - 4
 *
 * Will point to explicit edges:
 *
 * 0 - 1 - 2 - 3 - 4
 *
 * Each of these explicit edges will point back to the reduced edge
 *
 * Explicit Edges       Reduced Edge
 *
 *      0 - 1      ->       0 - 4
 *      1 - 2      ->       0 - 4
 *      2 - 3      ->       0 - 4
 *      3 - 4      ->       0 - 4
 */
typedef boost::bimap<boost::bimaps::multiset_of<tools::Edge>,
                     boost::bimaps::set_of<tools::Edge>>
    reduced_edge_to_edges_map;

/**
 * \brief Simple class for storing the connections between motifs and the
 * underlying beads that are part of the connection.
 *
 * A graph like this
 *
 * 1 - 2 - 3 - 4
 *         |   |
 *         5 - 6
 *
 * Will be broken up into two motifs a line motif and a loop motif the
 * BeadMotifConnector tracks the connections between the now independ motifs
 *
 * Motif 0            Motif 1
 *
 * 1 - 2              3 - 4
 *                    |   |
 *                    5 - 6
 *
 * The motif connection is stored as edge:     0 - 1
 * The corresbonding bead connection as edge:  2 - 3
 *
 **/
class BeadMotifConnector {
 public:
  void AddMotifAndBeadEdge(tools::Edge motif_edge, tools::Edge bead_edge);
  /// Returns the bead edges connecting the motifs specified by motif_edge
  std::vector<tools::Edge> getBeadEdges(const tools::Edge& motif_edge) const;
  /// Returns all the bead edges connecting the motifs
  std::vector<tools::Edge> getBeadEdges() const;

  /// Returns the motifs involved between two beads given by bead_edge
  tools::Edge getMotifEdge(const tools::Edge bead_edge) const;
  /// Returns all the motif edges
  std::unordered_set<tools::Edge> getMotifEdges() const;

 private:
  reduced_edge_to_edges_map;
  motif_and_bead_edges_;
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_BEADMOTIFCONNECTOR_H
