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
#include "../../include/votca/csg/beadmotifalgorithms.h"
#include "../../include/votca/csg/beadstructurealgorithms.h"
#include <list>
#include <unordered_map>
#include <votca/csg/beadmotif.h>
#include <votca/tools/graphalgorithm.h>

using namespace votca::tools;
using namespace std;

namespace votca {
namespace csg {

/******************************************************************************
 * Declartion of Internal Private Functions
 ******************************************************************************/

/**
 * \brief This function will take a list of bead_motifs and split them up
 * depending on whether they are of type multiple_structures or not.
 *
 * Elements are moved out of the bead_motifs they are not destroyed.
 *
 * @param[in,out] - list of bead motifs
 * @return - list of bead motifs of type multiple structures.
 **/
list<BeadMotif> moveMultipleStructureMotifsToSeparateList_(
    list<BeadMotif>& bead_motifs);

/**
 * \brief Separates the bead motifs in a list into two separate list the first
 * containing simple motifs and the scond containing complex motifs
 *
 * @param[in,out] - bead motifs of all types
 * @return - pair of lists of bead motifs the first contains simple and second
 * complex bead motifs
 **/
pair<list<BeadMotif>, list<BeadMotif>> separateComplexAndSimpleMotifs_(
    list<BeadMotif>& bead_motifs);

/**
 * The purpose of this function is to determine which edges connected to the
 * starting vertex in the full graph were explored for a branch. We are
 * starting with a contracted edge from the reduced graph. So it must first
 * be expanded and then it can be explored.
 **/
int determineEdgesOfBranchesWithSingleConnectionToJunction_(
    int junction, ReducedGraph& reduced_graph,
    unordered_map<Edge, bool>& remove_edges,
    unordered_map<int, vector<Edge>>
        contracted_branch_edges_connected_to_junction);

/**
 * \brief Determines which edges are connected to the junction and stored them
 * in the rmoeve_edges map.
 **/
void removeAllEdgesConnectedToVertex_(int junction, Graph& full_graph,
                                      unordered_map<Edge, bool>& remove_edges);

/**
 * \brief takes a list of bead motifs and places them in a map while assigning
 * an id to them.
 **/
unordered_map<int, BeadMotif> convertToUnorderedMapAndAssignId_(
    list<BeadMotif> motifs);

/**
 * \brief Calculate which edges to remove from a graph to split it up correctly
 *
 * This will depend on the number of junctions and the number of edges of each
 * branch that is connected to each of the junctions.
 **/
void calculateEdgesToRemove_(Graph& full_graph,
                             unordered_map<Edge, bool>& remove_edges);

/**
 * \brief Takes a bead motif of type single structure and splits it into
 * smaller components
 *
 * This internal function can be called only when we are dealng with a bead
 * motif that is known be a single network and it is not a simple component
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
 **/
list<BeadMotif> splitConnectedBeadMotif_(BeadMotif bead_motif);

list<BeadMotif> breakMultipleMotifs_(list<BeadMotif> motifs);
/******************************************************************************
 * Internal Private Functions
 ******************************************************************************/

list<BeadMotif> moveMultipleStructureMotifsToSeparateList_(
    list<BeadMotif>& bead_motifs) {

  list<BeadMotif> bead_motif_multiple_structures;
  list<BeadMotif>::iterator bead_motif_iter = bead_motifs.begin();

  while (bead_motif_iter != bead_motifs.end()) {
    if (bead_motif_iter->getType() ==
        BeadMotif::MotifType::multiple_structures) {
      auto temp_iter = bead_motif_iter;
      ++temp_iter;
      bead_motif_multiple_structures.splice(
          bead_motif_multiple_structures.end(), bead_motifs, bead_motif_iter);
      bead_motif_iter = temp_iter;
    } else {
      ++bead_motif_iter;
    }
  }
  return bead_motif_multiple_structures;
}

pair<list<BeadMotif>, list<BeadMotif>> separateComplexAndSimpleMotifs_(
    list<BeadMotif>& bead_motifs) {
  list<BeadMotif> complex_motifs;
  list<BeadMotif> simple_motifs;

  list<BeadMotif>::iterator motif_iterator = bead_motifs.begin();
  while (motif_iterator != bead_motifs.end()) {
    if (motif_iterator->isMotifSimple()) {
      auto temp_iter = motif_iterator;
      ++temp_iter;
      simple_motifs.splice(simple_motifs.end(), bead_motifs, motif_iterator);
      motif_iterator = temp_iter;
    } else {
      assert(motif_iterator->getType() != BeadMotif::MotifType::undefined);
      auto temp_iter = motif_iterator;
      ++temp_iter;
      complex_motifs.splice(complex_motifs.end(), bead_motifs, motif_iterator);
      motif_iterator = temp_iter;
    }
  }
  return pair<list<BeadMotif>, list<BeadMotif>>{simple_motifs, complex_motifs};
}

int determineEdgesOfBranchesWithSingleConnectionToJunction_(
    int junction, ReducedGraph& reduced_graph,
    unordered_map<Edge, bool>& remove_edges,
    unordered_map<int, vector<Edge>>
        contracted_branch_edges_connected_to_junction) {

  int number_branches_with_more_than_one_bond = 0;
  // Cycle the branches
  for (pair<int, vector<Edge>> branch_and_contracted_edges :
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

void removeAllEdgesConnectedToVertex_(int junction, Graph& full_graph,
                                      unordered_map<Edge, bool>& remove_edges) {

  vector<Edge> neigh_edges = full_graph.getNeighEdges(junction);
  for (Edge& edge : neigh_edges) {
    remove_edges[edge] = true;
  }
  return;
}

unordered_map<int, BeadMotif> convertToUnorderedMapAndAssignId_(
    list<BeadMotif> motifs) {
  int motif_id = 0;
  unordered_map<int, BeadMotif> motifs2;
  for (BeadMotif& motif : motifs) {
    motifs2[motif_id] = motif;
  }
  return motifs2;
}

void calculateEdgesToRemove_(Graph& full_graph,
                             unordered_map<Edge, bool>& remove_edges) {

  ReducedGraph reduced_graph = reduceGraph(full_graph);

  vector<int> junctions = reduced_graph.getJunctions();
  for (int& junction : junctions) {

    vector<Edge> neighboring_edges = reduced_graph.getNeighEdges(junction);
    unordered_map<Edge, bool> contracted_edge_explored;
    for (Edge& edge : neighboring_edges) {
      contracted_edge_explored[edge] = false;
    }

    int branch_index = 0;
    unordered_map<int, vector<Edge>>
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

      int number_branches_with_more_than_one_bond =
          determineEdgesOfBranchesWithSingleConnectionToJunction_(
              junction, reduced_graph, remove_edges,
              contracted_branch_edges_connected_to_junction);

      // If more than one branch claims the junction and they both have more
      // than one connection the junction neigther branch can claim the junction
      // thus all connections to the junction will be cut.
      if (number_branches_with_more_than_one_bond > 1) {
        removeAllEdgesConnectedToVertex_(junction, full_graph, remove_edges);
      }
    }
  }  // for(int & junction : junctions)
}

list<BeadMotif> splitConnectedBeadMotif_(BeadMotif bead_motif) {

  assert(bead_motif.getType() == BeadMotif::MotifType::single_structure);

  Graph full_graph = bead_motif.getGraph();
  vector<Edge> all_edges = full_graph.getEdges();
  unordered_map<Edge, bool> remove_edges;
  for (Edge& edge : all_edges) {
    remove_edges[edge] = false;
  }

  calculateEdgesToRemove_(full_graph, remove_edges);

  vector<int> all_vertices = full_graph.getVertices();

  BeadStructure new_beadstructure;
  for (int& vertex : all_vertices) {
    new_beadstructure.AddBead(bead_motif.getBead(vertex));
  }

  for (pair<const Edge, bool>& edge_and_remove : remove_edges) {
    if (edge_and_remove.second == false) {
      Edge edge = edge_and_remove.first;
      new_beadstructure.ConnectBeads(edge.getEndPoint1(), edge.getEndPoint2());
    }
  }

  list<BeadMotif> split_motifs =
      breakIntoMotifs<list<BeadMotif>>(new_beadstructure);
  return split_motifs;
}

list<BeadMotif> breakMultipleMotifs_(list<BeadMotif> motifs) {
  list<BeadMotif> all_motifs;
  for (BeadMotif& motif : motifs) {
    list<BeadMotif> motifs_temp = breakIntoMotifs<list<BeadMotif>>(motif);
    all_motifs.splice(all_motifs.end(), motifs_temp);
  }
  return all_motifs;
}
/******************************************************************************
 * Public Functions
 ******************************************************************************/

/**
 * \brief This function will take a beadmotif and break it into its simple
 *motifs
 *
 * It will return the simple motif and the edges connecting each motif with
 *other motifs as well as the edges describing which vertices are connecting the
 *motifs
 **/
unordered_map<int, BeadMotif> breakIntoSimpleMotifs(BeadMotif bead_motif) {

  list<BeadMotif> motifs_simple;
  list<BeadMotif> motifs_complex;
  list<BeadMotif> motifs;
  motifs.push_back(bead_motif);

  do {  // while(motifs_complex.size()!=0)

    for (BeadMotif& complex_motif : motifs_complex) {
      list<BeadMotif> motifs_temp = splitConnectedBeadMotif_(complex_motif);
      motifs.splice(motifs.end(), motifs_temp);
    }

    {  // Sorting block
      list<BeadMotif> motifs_multiple_structure =
          moveMultipleStructureMotifsToSeparateList_(motifs);
      list<BeadMotif> motifs_temp =
          breakMultipleMotifs_(motifs_multiple_structure);
      motifs.splice(motifs.end(), motifs_temp);
      auto simple_and_complex_motifs = separateComplexAndSimpleMotifs_(motifs);
      motifs_simple = simple_and_complex_motifs.first;
      motifs_complex = simple_and_complex_motifs.second;
      assert(motifs.size() == 0);
    }  // End of Sorting block

  } while (motifs_complex.size() != 0);

  return convertToUnorderedMapAndAssignId_(motifs_simple);
}

}  // namespace csg
}  // namespace votca
