/*
 *            Copyright 2009-2019 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writingraph, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include <list>

#include "../../include/votca/tools/graph.h"
#include "../../include/votca/tools/graph_bf_visitor.h"
#include "../../include/votca/tools/graph_df_visitor.h"
#include "../../include/votca/tools/graphalgorithm.h"
#include "../../include/votca/tools/graphvisitor.h"

using namespace std;

namespace votca {
namespace tools {

/**********************
 * Internal Functions *
 **********************/

/********************
 * Public Functions *
 ********************/
bool singleNetwork(Graph& graph, GraphVisitor& graph_visitor) {
  exploreGraph(graph, graph_visitor);
  return graph_visitor.getExploredVertices().size() ==
             graph.getVertices().size() &&
         graph.getIsolatedNodes().size() == 0;
}

std::set<Edge> exploreBranch(Graph g, Index starting_vertex, const Edge& edge) {
  // Check if the starting vertex is in the graph
  if (!g.vertexExist(starting_vertex)) {
    throw invalid_argument(
        "Cannot explore a branch of the graph when the "
        "exploration is started from a vertex that does not exist within the "
        "graph.");
  }
  if (!g.edgeExist(edge)) {
    throw invalid_argument(
        "Edge does not exist in the graph so exploration of"
        " the graph branch cannot continue");
  }
  if (!edge.contains(starting_vertex)) {
    throw invalid_argument(
        "The edge determining which branch to explore does "
        "not contain the starting vertex.");
  }

  set<Edge> branch_edges;
  if (edge.getEndPoint1() == edge.getEndPoint2()) {
    branch_edges.insert(edge);
    return branch_edges;
  }
  Graph_BF_Visitor gv_breadth_first;
  GraphNode gn;
  pair<Index, GraphNode> pr_gn(starting_vertex, gn);
  gv_breadth_first.exploreNode(pr_gn, g);
  gv_breadth_first.setStartingVertex(edge.getOtherEndPoint(starting_vertex));
  gv_breadth_first.initialize(g);

  branch_edges.insert(edge);
  while (!gv_breadth_first.queEmpty()) {
    Edge ed = gv_breadth_first.nextEdge(g);
    branch_edges.insert(ed);
    gv_breadth_first.exec(g, ed);
  }

  vector<Edge> neigh_edges = g.getNeighEdges(starting_vertex);
  for (Edge& ed : neigh_edges) {
    Index neigh_vertex = ed.getOtherEndPoint(starting_vertex);
    if (neigh_vertex != starting_vertex) {
      if (gv_breadth_first.vertexExplored(neigh_vertex)) {
        branch_edges.insert(ed);
      }
    }
  }

  return branch_edges;
}

ReducedGraph reduceGraph(Graph graph) {

  /****************************
   * Internal Function Class
   ****************************/

  /**
   * \brief This class is to help keep track of which vertices have and have not
   * been explored.
   *
   * It is only used within the function and hence the encapsulation, it is not
   * a public class and is not meant to be.
   **/
  class ExplorationRecord {
   private:
    unordered_map<Index, std::pair<bool, Index>> vertex_explored_;
    size_t unexplored_vertex_count_;

   public:
    explicit ExplorationRecord(
        const unordered_map<Index, std::pair<bool, Index>>& vertex_explored)
        : vertex_explored_(vertex_explored),
          unexplored_vertex_count_(vertex_explored.size()){};

    void explore(Index vertex) {
      vertex_explored_[vertex].first = true;
      --unexplored_vertex_count_;
    }

    bool unexploredVerticesExist() { return unexplored_vertex_count_ > 0; }

    Index getUnexploredVertex() {

      vector<Index> remaining_unexplored;
      for (const pair<Index, pair<bool, Index>>& vertex_record :
           vertex_explored_) {
        bool vertex_explored = vertex_record.second.first;
        if (!vertex_explored) {
          Index degree = vertex_record.second.second;
          if (degree > 2) {
            return vertex_record.first;
          }
          remaining_unexplored.push_back(vertex_record.first);
        }
      }

      // Search tips next
      for (const Index& vertex : remaining_unexplored) {
        if (vertex_explored_[vertex].second == 1) {
          return vertex;
        }
      }

      // Finally if there are no tips or junctions left we will return a vertex
      // of degree 2 if one exists
      for (const Index& vertex : remaining_unexplored) {
        if (!vertex_explored_[vertex].first) {
          return vertex;
        }
      }

      throw runtime_error(
          "You cannot get an unexplored vertex as they have all"
          " been explored.");
    }
  };  // Class ExplorationRecord

  unordered_map<Index, pair<bool, Index>> unexplored_vertices;
  vector<Index> vertices = graph.getVertices();
  for (const Index& vertex : vertices) {
    unexplored_vertices[vertex] =
        pair<bool, Index>(false, graph.getDegree(vertex));
  }

  ExplorationRecord exploration_record(unexplored_vertices);

  vector<vector<Index>> chains;

  while (exploration_record.unexploredVerticesExist()) {
    Graph_DF_Visitor graph_visitor;
    Index starting_vertex = exploration_record.getUnexploredVertex();
    exploration_record.explore(starting_vertex);
    graph_visitor.setStartingVertex(starting_vertex);
    graph_visitor.initialize(graph);

    vector<Index> chain{starting_vertex};
    Index old_vertex = starting_vertex;
    bool new_chain = false;
    while (!graph_visitor.queEmpty()) {
      Edge edge = graph_visitor.nextEdge(graph);
      vector<Index> unexplored_vertex = graph_visitor.getUnexploredVertex(edge);

      if (new_chain) {
        if (unexplored_vertex.size() == 0) {
          old_vertex = edge.getEndPoint1();
          chain.push_back(old_vertex);
          new_chain = false;
        } else {
          old_vertex = edge.getOtherEndPoint(unexplored_vertex.at(0));
          chain.push_back(old_vertex);
          new_chain = false;
        }
      }
      Index new_vertex = edge.getOtherEndPoint(old_vertex);
      if (unexplored_vertex.size() == 0) {
        chain.push_back(new_vertex);
        chains.push_back(chain);
        chain.clear();
        new_chain = true;
      } else if (graph.getDegree(new_vertex) == 1) {
        chain.push_back(new_vertex);
        chains.push_back(chain);
        chain.clear();
        exploration_record.explore(new_vertex);
        new_chain = true;
      } else if (graph.getDegree(new_vertex) > 2) {
        chain.push_back(new_vertex);
        chains.push_back(chain);
        chain.clear();
        exploration_record.explore(new_vertex);
        new_chain = true;
      } else if (unexplored_vertex.size() == 1) {
        chain.push_back(new_vertex);
        old_vertex = new_vertex;
        exploration_record.explore(new_vertex);
      }

      graph_visitor.exec(graph, edge);
    }
  }
  vector<ReducedEdge> reduced_edges;
  for (vector<Index> chain : chains) {
    ReducedEdge reduced_edge(chain);
    reduced_edges.push_back(reduced_edge);
  }

  ReducedGraph reduced_g(reduced_edges);
  auto nodes_graph = graph.getNodes();
  reduced_g.clearNodes();
  reduced_g.copyNodes(graph);
  auto nodes_reduced_g = reduced_g.getNodes();
  return reduced_g;
}

vector<Graph> decoupleIsolatedSubGraphs(Graph graph) {

  const std::vector<Index>& vertices = graph.getVertices();
  // bool vector to see if vertex is already part of graph
  std::vector<bool> vertex_analysed = std::vector<bool>(vertices.size(), false);

  std::vector<Graph> subGraphs;
  Index i = 0;
  while (i < Index(vertices.size())) {
    if (vertex_analysed[i]) {
      i++;
      continue;
    }
    Graph_BF_Visitor graph_visitor_breadth_first;
    graph_visitor_breadth_first.setStartingVertex(vertices[i]);
    exploreGraph(graph, graph_visitor_breadth_first);
    set<Index> sub_graph_explored_vertices =
        graph_visitor_breadth_first.getExploredVertices();

    for (Index vertex : sub_graph_explored_vertices) {
      for (Index j = 0; j < Index(vertices.size()); j++) {
        if (vertex_analysed[j]) {
          continue;
        }
        if (vertex == vertices[j]) {
          vertex_analysed[j] = true;
          break;
        }
      }
    }

    set<Edge> sub_graph_edges;
    unordered_map<Index, GraphNode> sub_graph_nodes;
    for (Index vertex : sub_graph_explored_vertices) {

      for (const Edge& sub_graph_edge : graph.getNeighEdges(vertex)) {
        sub_graph_edges.insert(sub_graph_edge);
      }
      sub_graph_nodes[vertex] = graph.getNode(vertex);
    }
    subGraphs.push_back(
        Graph(std::vector<Edge>(sub_graph_edges.begin(), sub_graph_edges.end()),
              sub_graph_nodes));
  }
  return subGraphs;
}

void exploreGraph(Graph& graph, GraphVisitor& graph_visitor) {

  if (!graph.vertexExist(graph_visitor.getStartingVertex())) {
    string err = "Cannot explore graph starting at vertex " +
                 to_string(graph_visitor.getStartingVertex()) +
                 " because it does not exist in the "
                 "graph. You can change the starting vertex by using the "
                 "setStartingVertex method of the visitor instance.";
    throw invalid_argument(err);
  }
  // Create a list of all vertices and determine if they have all been
  // explored when the visitor que is empty
  graph_visitor.initialize(graph);
  while (!graph_visitor.queEmpty()) {
    Edge edge = graph_visitor.nextEdge(graph);
    graph_visitor.exec(graph, edge);
  }
}
}  // namespace tools
}  // namespace votca
