/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#include <list>

#include <votca/tools/graph.h>
#include <votca/tools/graph_bf_visitor.h>
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

bool singleNetwork(Graph g, GraphVisitor& gv) {
  exploreGraph(g, gv);
  return gv.getExploredVertices().size() == g.getVertices().size() &&
         g.getIsolatedNodes().size() == 0;
}

vector<shared_ptr<Graph>> decoupleIsolatedSubGraphs(Graph g) {

  vector<shared_ptr<Graph>> subGraphs;
  {
    Graph_BF_Visitor gv_breadth_first;
    if (singleNetwork(g, gv_breadth_first)) {
      subGraphs.push_back(make_shared<Graph>(g));
      return subGraphs;
    }
  }

  auto vertices_list = vectorToList_(g.getVertices());

  auto v_it = vertices_list.begin();
  unordered_map<int, GraphNode> sub_graph_nodes;
  while (v_it != vertices_list.end()) {
    Graph_BF_Visitor gv_breadth_first;
    gv_breadth_first.setStartingVertex(*v_it);
    exploreGraph(g, gv_breadth_first);

    ++v_it;
    auto sub_graph_explored_vertices = gv_breadth_first.getExploredVertices();
    set<Edge> sub_graph_edges;
    for (auto sub_graph_vertex : sub_graph_explored_vertices) {
      if (*v_it == sub_graph_vertex) {
        ++v_it;
      }
      vertices_list.remove(sub_graph_vertex);

      auto sub_graph_neigh_edges = g.getNeighEdges(sub_graph_vertex);
      for (auto sub_graph_edge : sub_graph_neigh_edges) {
        sub_graph_edges.insert(sub_graph_edge);
      }

      sub_graph_nodes[sub_graph_vertex] = g.getNode(sub_graph_vertex);
    }

    auto sub_graph_vector_edges = edgeSetToVector_(sub_graph_edges);

    subGraphs.push_back(
        make_shared<Graph>(Graph(sub_graph_vector_edges, sub_graph_nodes)));
  }
  return subGraphs;
}

void exploreGraph(Graph& g, GraphVisitor& gv) {

  if (!g.vertexExist(gv.getStartingVertex())) {
    string err = "Cannot explore graph starting at vertex " +
                 to_string(gv.getStartingVertex()) +
                 " because it does not exist in the "
                 "graph. You can change the starting vertex by using the "
                 "setStartingVertex method of the visitor instance.";
    throw invalid_argument(err);
  }
  // Create a list of all vertices and determine if they have all been
  // explored when the visitor que is empty
  gv.initialize(g);
  while (!gv.queEmpty()) {
    auto ed = gv.nextEdge(g);
    gv.exec(g, ed);
  }
}
}
}
