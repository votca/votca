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
/ */

// Standard includes
#include <cassert>
#include <cstddef>
#include <utility>

// VOTCA includes
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/reducedgraph.h>

// Local VOTCA includes
#include "../../include/votca/csg/beadmotifconnector.h"
#include "votca/csg/beadmotifalgorithms.h"

namespace votca {
namespace csg {
class BeadMotif;
}  // namespace csg
namespace tools {
class Graph;
}  // namespace tools
}  // namespace votca

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

/******************************************************************************
 * Locally Global Variables
 ******************************************************************************/

typedef std::pair<Index, BeadMotif> IdMotif;

const Index unassigned_id = -1;

static Index motif_index = 0;

/******************************************************************************
 * Internal Private Functions
 ******************************************************************************/

/**
 * \brief Determines which edges are connected to the junction and stored them
 * in the remove_edges map.
 **/
void removeAllEdgesConnectedToVertex_(Index junction, Graph& full_graph,
                                      unordered_map<Edge, bool>& remove_edges) {

  vector<Edge> neigh_edges = full_graph.getNeighEdges(junction);
  for (Edge& edge : neigh_edges) {
    remove_edges[edge] = true;
  }
}
/**
 * The purpose of this function is to determine which edges connected to the
 * starting vertex in the full graph were explored for a branch. We are
 * starting with a contracted edge from the reduced graph. So it must first
 * be expanded and then it can be explored.
 **/
Index determineEdgesOfBranchesWithSingleConnectionToJunction_(
    Index junction, ReducedGraph& reduced_graph,
    unordered_map<Edge, bool>& remove_edges,
    unordered_map<Index, vector<Edge>>
        contracted_branch_edges_connected_to_junction) {

  Index number_branches_with_more_than_one_bond = 0;
  // Cycle the branches
  for (pair<Index, vector<Edge>> branch_and_contracted_edges :
       contracted_branch_edges_connected_to_junction) {
    bool edge_is_single_connection_to_junction = false;
    if (branch_and_contracted_edges.second.size() == 1) {
      // Make sure the branch is not a loop as well
      Edge contracted_edge = branch_and_contracted_edges.second.at(0);
      if (contracted_edge.loop() == false) {
        vector<vector<Edge>> expanded_edges =
            reduced_graph.expandEdge(contracted_edge);
        // Determine which of these edges are conneted to the vertex
        if (expanded_edges.size() == 1) {
          for (Edge& edge : expanded_edges.at(0)) {
            if (edge.contains(junction)) {
              remove_edges[edge] = true;
              edge_is_single_connection_to_junction = true;
              break;
            }
          }
        }
      }
    }
    if (edge_is_single_connection_to_junction == false) {
      ++number_branches_with_more_than_one_bond;
    }
  }

  return number_branches_with_more_than_one_bond;
}

/**
 * \brief Calculates which edges to remove from a graph, in order to split it up
 * correctly into simple motif types.
 *
 * This will depend on the number of junctions and the number of edges of each
 * branch that is connected to each of the junctions.
 **/
void calculateEdgesToRemove_(Graph& full_graph,
                             unordered_map<Edge, bool>& remove_edges) {

  ReducedGraph reduced_graph = reduceGraph(full_graph);

  vector<Index> junctions = reduced_graph.getJunctions();

  for (Index& junction : junctions) {
    vector<Edge> neighboring_edges = reduced_graph.getNeighEdges(junction);
    unordered_map<Edge, bool> contracted_edge_explored;
    for (Edge& edge : neighboring_edges) {
      contracted_edge_explored[edge] = false;
    }

    Index branch_index = 0;
    unordered_map<Index, vector<Edge>>
        contracted_branch_edges_connected_to_junction;
    for (Edge& branch_starting_edge : neighboring_edges) {
      if (contracted_edge_explored[branch_starting_edge] == false) {
        // The starting edge could be redundant but I guess it doesn't matter
        // redundant edges would mean they are part of the same branch anyway
        set<Edge> branch_edges =
            exploreBranch(reduced_graph, junction, branch_starting_edge);

        vector<Edge> starting_edges_explored_by_branch;
        for (Edge& branch_starting_edge2 : neighboring_edges) {
          if (branch_edges.count(branch_starting_edge2)) {
            contracted_edge_explored[branch_starting_edge2] = true;
            starting_edges_explored_by_branch.push_back(branch_starting_edge2);
          }
        }
        contracted_branch_edges_connected_to_junction[branch_index] =
            starting_edges_explored_by_branch;
        ++branch_index;
      }
    }  // for(Edge & brach_starting_edge : neighboring_edges)

    // If there is more than one branch
    if (branch_index > 1) {
      // We now know how many edges are connected to the starting node from each
      // branch We can use the full graph to decouple the edges we can begin by
      // snipping branches with only a single edge connecting them to starting
      // nodes.
      Index number_branches_with_more_than_one_bond =
          determineEdgesOfBranchesWithSingleConnectionToJunction_(
              junction, reduced_graph, remove_edges,
              contracted_branch_edges_connected_to_junction);

      // If more than one branch claims the junction, and each branch has more
      // than one connetion to the junction neigther branch can claim the
      // junction, in this case all connections to the junction will be cut.
      if (number_branches_with_more_than_one_bond > 1) {
        removeAllEdgesConnectedToVertex_(junction, full_graph, remove_edges);
      }
    }
  }  // for(Index & junction : junctions)
}

/******************************************************************************
 * Internal Private Class
 ******************************************************************************/

/**
 * \brief This is an internal class meant to deconstruct bead motifs into their
 * simple forms
 *
 * It is only meant to be used within the source file. It is simply more
 * convenient to write a class with shared internal variables then a series of
 * functions. However, everything this class can do is accessed via the
 * graph motif algorithm breakIntoSimpleMotifs. Do not pull this class out of
 * the file.
 **/
class MotifDeconstructor_ {
 public:
  /// Adds a motif with its corresponding id
  void AddMotif(IdMotif id_and_motif);

  /**
   * \brief Takes a bead motif of type single structure and deconstructs it
   * into smaller components.
   *
   * E.g.
   *
   * Scenario I
   *
   * 1 - 2
   * |   |
   * 4 - 3 - 5
   *     |   |
   *     6 _ 7
   *
   * Scenario II
   *
   * 1 - 2  8
   * |   | / \
   * 4 - 3 - 5
   *     |   |
   *     6 _ 7
   *
   * Scenario III
   *
   * 1   2 - 5           4 - 1     2 - 5
   * | \ | /
   * 4 - 3        =>            3
   *     | \
   *     6 -7                6 - 7
   *
   * Scenario IV
   *
   * 1 - 2 - 3  4 - 5       1 - 2 - 3    4 - 5
   *          \ /
   *           6 - 7    =>         6 - 7
   *           |   |               |   |
   *           8 _ 9               8 - 9
   *
   * In Scenario I and III all branches equally share the vertex 3 so it is
   * treated as a pivot vertex that is not owned by any of the motifs think
   * of AlQ3
   *
   * In Scenario II there are three edges connected to vertex 3 from one
   * branch and only two from the other we will also consider it as a pivot
   * vertex
   *
   * In scenario IV vertex 6 will be assigned to the branch that has two
   * connecting edges.
   *
   * @param[in,out] - BeadMotifConnector this keeps track of all the
   * connections between motifs. Both the bead ids and the motif ids are
   * tracked and stored
   **/
  void deconstructComplexSingleStructures(BeadMotifConnector& connector);

  /// Counts the number of complex motifs that have not yet been broken down
  /// into simple motifs
  size_t CountComplexMotifs() const;

  /// returns whether all the connections between the simple motifs have
  /// been correctly handled.
  bool ConnectionsCompleted() const { return bead_edges_removed_.size() == 0; }

  /// Returns a map of all the simple motifs with their ids
  unordered_map<Index, BeadMotif> getSimpleMotifs();

 private:
  list<IdMotif> motifs_simple_;
  list<IdMotif> motifs_complex_single_structure_;
  list<IdMotif> motifs_multiple_structures_;
  list<Edge> bead_edges_removed_;

  // Returns the ids of the simple motifs that were assigned ids
  void sortMotifsAndAssignIdsToSimpleMotifs_(list<BeadMotif>& bead_motifs);
  void determineMotifConnections_(BeadMotifConnector& connector);
};

unordered_map<Index, BeadMotif> MotifDeconstructor_::getSimpleMotifs() {
  unordered_map<Index, BeadMotif> simple_motifs;
  for (IdMotif& id_and_motif : motifs_simple_) {
    simple_motifs.insert(id_and_motif);
  }
  return simple_motifs;
}

size_t MotifDeconstructor_::CountComplexMotifs() const {
  return motifs_complex_single_structure_.size() +
         motifs_multiple_structures_.size();
}

void MotifDeconstructor_::AddMotif(IdMotif id_and_motif) {
  if (id_and_motif.second.getType() ==
      BeadMotif::MotifType::multiple_structures) {
    motifs_multiple_structures_.push_back(id_and_motif);
  } else if (id_and_motif.second.getType() ==
             BeadMotif::MotifType::single_structure) {
    motifs_complex_single_structure_.push_back(id_and_motif);
  } else if (id_and_motif.second.isMotifSimple()) {
    motifs_simple_.push_back(id_and_motif);
  } else {
    assert(false && "Error the motif type must be a simple or complex type.");
  }
}

void MotifDeconstructor_::deconstructComplexSingleStructures(
    BeadMotifConnector& connector) {
  list<BeadMotif> split_motifs;
  for (IdMotif& id_and_bead_motif : motifs_complex_single_structure_) {
    Graph full_graph = id_and_bead_motif.second.getGraph();
    vector<Edge> all_edges = full_graph.getEdges();
    unordered_map<Edge, bool> remove_edges;
    for (Edge& edge : all_edges) {
      remove_edges[edge] = false;
    }

    calculateEdgesToRemove_(full_graph, remove_edges);

    vector<Index> all_vertices = full_graph.getVertices();

    std::vector<Edge> new_beadstructure_edges;
    for (pair<const Edge, bool>& edge_and_remove : remove_edges) {
      if (edge_and_remove.second == false) {

        new_beadstructure_edges.push_back(edge_and_remove.first);
      } else {
        bead_edges_removed_.push_back(edge_and_remove.first);
      }
    }
    BeadStructure new_beadstructure = id_and_bead_motif.second.getSubStructure(
        all_vertices, new_beadstructure_edges);

    list<BeadMotif> split_motifs_temp =
        breakIntoMotifs<list<BeadMotif>>(new_beadstructure);

    split_motifs.splice(split_motifs.end(), split_motifs_temp);
  }

  motifs_complex_single_structure_.clear();
  sortMotifsAndAssignIdsToSimpleMotifs_(split_motifs);
  determineMotifConnections_(connector);
}

void MotifDeconstructor_::determineMotifConnections_(
    BeadMotifConnector& connector) {

  // Cycle the edges
  list<Edge>::iterator edge_iterator = bead_edges_removed_.begin();
  while (edge_iterator != bead_edges_removed_.end()) {
    assert(edge_iterator->loop() == false);
    Index bead_id1 = edge_iterator->getEndPoint1();
    Index bead_id2 = edge_iterator->getEndPoint2();
    // Cycle the motifs
    list<IdMotif>::iterator motif_bead1_iterator;
    for (motif_bead1_iterator = motifs_simple_.begin();
         motif_bead1_iterator != motifs_simple_.end(); ++motif_bead1_iterator) {

      if (motif_bead1_iterator->second.BeadExist(bead_id1)) {
        break;
      }
    }

    list<IdMotif>::iterator motif_bead2_iterator;
    for (motif_bead2_iterator = motifs_simple_.begin();
         motif_bead2_iterator != motifs_simple_.end(); ++motif_bead2_iterator) {

      if (motif_bead2_iterator->second.BeadExist(bead_id2)) {
        break;
      }
    }

    // We remove the edge from the list and add it as a connection
    if (motif_bead2_iterator != motifs_simple_.end() &&
        motif_bead1_iterator != motifs_simple_.end()) {
      Edge motif_edge(motif_bead1_iterator->first, motif_bead2_iterator->first);
      connector.AddMotifAndBeadEdge(motif_edge, *edge_iterator);
      edge_iterator = bead_edges_removed_.erase(edge_iterator);
    } else {
      ++edge_iterator;
    }
  }  // while( edge_iterator != bead_edges.end() )
}

void MotifDeconstructor_::sortMotifsAndAssignIdsToSimpleMotifs_(
    list<BeadMotif>& bead_motifs) {

  list<BeadMotif>::iterator bead_motif_iter = bead_motifs.begin();
  while (bead_motif_iter != bead_motifs.end()) {

    Index new_motif_id = unassigned_id;
    if (bead_motif_iter->isMotifSimple()) {
      new_motif_id = motif_index;
      ++motif_index;
    }
    pair<Index, BeadMotif> id_and_motif(new_motif_id,
                                        std::move(*bead_motif_iter));
    AddMotif(id_and_motif);
    bead_motif_iter = bead_motifs.erase(bead_motif_iter);
  }
}

/******************************************************************************
 * Public Functions
 ******************************************************************************/

pair<unordered_map<Index, BeadMotif>, BeadMotifConnector> breakIntoSimpleMotifs(
    BeadMotif bead_motif) {

  MotifDeconstructor_ motif_manipulator;
  motif_manipulator.AddMotif(IdMotif(unassigned_id, bead_motif));

  BeadMotifConnector connector;

  do {
    motif_manipulator.deconstructComplexSingleStructures(connector);
  } while (motif_manipulator.CountComplexMotifs());

  assert(motif_manipulator.ConnectionsCompleted());
  unordered_map<Index, BeadMotif> motifs_simple_map =
      motif_manipulator.getSimpleMotifs();
  auto simple_motifs_and_connector =
      pair<unordered_map<Index, BeadMotif>, BeadMotifConnector>(
          motifs_simple_map, connector);
  return simple_motifs_and_connector;
}

}  // namespace csg
}  // namespace votca
