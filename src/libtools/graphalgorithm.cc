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

#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
#include <votca/tools/graph_df_visitor.h>
#include <votca/tools/graphalgorithm.h>
#include <votca/tools/graphvisitor.h>

using namespace std;

namespace votca {
namespace tools {

/**********************
 * Internal Functions *
 **********************/

list<int> vectorToList_(vector<int> values) {
  list<int> values_list;
  copy(values.begin(), values.end(), back_inserter(values_list));
  return values_list;
}

vector<Edge> edgeSetToVector_(set<Edge> edges) {
  vector<Edge> edges_vector;
  copy(edges.begin(), edges.end(), back_inserter(edges_vector));
  return edges_vector;
}

/********************
 * Public Functions *
 ********************/

bool singleNetwork(const Graph& graph, GraphVisitor& graph_visitor) {
  exploreGraph(graph, graph_visitor);
  return graph_visitor.getExploredVertices().size() ==
             graph.getVertices().size() &&
         graph.getIsolatedNodes().size() == 0;
}

std::set<Edge> exploreBranch(const Graph& g, const int starting_vertex,
                             const Edge& edge) {
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
  pair<int, GraphNode> pr_gn(starting_vertex, gn);
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
    int neigh_vertex = ed.getOtherEndPoint(starting_vertex);
    if (neigh_vertex != starting_vertex) {
      if (gv_breadth_first.vertexExplored(neigh_vertex)) {
        branch_edges.insert(ed);
      }
    }
  }

  return branch_edges;
}

ReducedGraph reduceGraph(const Graph& graph) {

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
    unordered_map<int, std::pair<bool, int>> vertex_explored_;
    size_t unexplored_vertex_count_;

   public:
    ExplorationRecord(unordered_map<int, std::pair<bool, int>> vertex_explored)
        : vertex_explored_(vertex_explored),
          unexplored_vertex_count_(vertex_explored.size()){};

    void explore(int vertex) {
      vertex_explored_[vertex].first = true;
      --unexplored_vertex_count_;
    }

    bool unexploredVerticesExist() { return unexplored_vertex_count_ > 0; }

    int getUnexploredVertex() {

      vector<int> remaining_unexplored;
      for (const pair<int, pair<bool, int>>& vertex_record : vertex_explored_) {
        bool vertex_explored = vertex_record.second.first;
        if (!vertex_explored) {
          int degree = vertex_record.second.second;
          if (degree > 2) {
            return vertex_record.first;
          }
          remaining_unexplored.push_back(vertex_record.first);
        }
      }

      // Search tips next
      for (const int& vertex : remaining_unexplored) {
        if (vertex_explored_[vertex].second == 1) {
          return vertex;
        }
      }

      // Finally if there are no tips or junctions left we will return a vertex
      // of degree 2 if one exists
      for (const int& vertex : remaining_unexplored) {
        if (!vertex_explored_[vertex].first) return vertex;
      }

      throw runtime_error(
          "You cannot get an unexplored vertex as they have all"
          " been explored.");
    }
  };  // Class ExplorationRecord

  unordered_map<int, pair<bool, int>> unexplored_vertices;
  vector<int> vertices = graph.getVertices();
  for (const int& vertex : vertices) {
    unexplored_vertices[vertex] =
        pair<bool, int>(false, graph.getDegree(vertex));
  }

  ExplorationRecord exploration_record(unexplored_vertices);

  vector<vector<int>> chains;

  while (exploration_record.unexploredVerticesExist()) {
    Graph_DF_Visitor graph_visitor;
    int starting_vertex = exploration_record.getUnexploredVertex();
    exploration_record.explore(starting_vertex);
    graph_visitor.setStartingVertex(starting_vertex);
    graph_visitor.initialize(graph);

    vector<int> chain{starting_vertex};
    int old_vertex = starting_vertex;
    bool new_chain = false;
    while (!graph_visitor.queEmpty()) {
      Edge edge = graph_visitor.nextEdge(graph);
      vector<int> unexplored_vertex = graph_visitor.getUnexploredVertex(edge);

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
      int new_vertex = edge.getOtherEndPoint(old_vertex);
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
  for (vector<int> chain : chains) {
    ReducedEdge reduced_edge(chain);
    reduced_edges.push_back(reduced_edge);
  }

  ReducedGraph reduced_g(reduced_edges);
  auto nodes_graph = graph.getNodes();
  reduced_g.copyNodes(graph);
  auto nodes_reduced_g = reduced_g.getNodes();
  return reduced_g;
}

vector<Graph> decoupleIsolatedSubGraphs(const Graph& graph) {

  list<int> vertices_list = vectorToList_(graph.getVertices());
  vector<Graph> subGraphs;
  {
    Graph_BF_Visitor graph_visitor_breadth_first;
    graph_visitor_breadth_first.setStartingVertex(*vertices_list.begin());
    if (singleNetwork(graph, graph_visitor_breadth_first)) {
      subGraphs.push_back(graph);
      return subGraphs;
    }
  }

  list<int>::iterator vertices_iterator = vertices_list.begin();
  unordered_map<int, GraphNode> sub_graph_nodes;
  while (vertices_iterator != vertices_list.end()) {
    Graph_BF_Visitor graph_visitor_breadth_first;
    graph_visitor_breadth_first.setStartingVertex(*vertices_iterator);
    exploreGraph(graph, graph_visitor_breadth_first);

    ++vertices_iterator;
    set<int> sub_graph_explored_vertices =
        graph_visitor_breadth_first.getExploredVertices();
    set<Edge> sub_graph_edges;

    set<int>::iterator sub_graph_vertex_it =
        sub_graph_explored_vertices.begin();
    while (sub_graph_vertex_it != sub_graph_explored_vertices.end()) {
      if (*vertices_iterator == *sub_graph_vertex_it) {
        ++vertices_iterator;
      }
      vertices_list.remove(*sub_graph_vertex_it);

      vector<Edge> sub_graph_neigh_edges =
          graph.getNeighEdges(*sub_graph_vertex_it);
      for (const Edge sub_graph_edge : sub_graph_neigh_edges) {
        sub_graph_edges.insert(sub_graph_edge);
      }

      sub_graph_nodes[*sub_graph_vertex_it] =
          graph.getNode(*sub_graph_vertex_it);
      ++sub_graph_vertex_it;
    }

    vector<Edge> sub_graph_vector_edges = edgeSetToVector_(sub_graph_edges);

    subGraphs.push_back(Graph(sub_graph_vector_edges, sub_graph_nodes));
  }
  return subGraphs;
}

void exploreGraph(const Graph& graph, GraphVisitor& graph_visitor) {

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
