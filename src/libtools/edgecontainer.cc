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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <algorithm>
#include <cassert>
#include <exception>
#include <set>
#include <vector>
#include <votca/tools/edge.h>
#include <votca/tools/edgecontainer.h>

namespace votca {
namespace tools {

using namespace std;

EdgeContainer::EdgeContainer(Edge edge) { addEdge(edge); }

EdgeContainer::EdgeContainer(vector<Edge> edges) {
  for (Edge& edge : edges) {
    addEdge(edge);
  }
}

long int EdgeContainer::getMaxDegree(void) const {
  long int max = 0;
  for (const auto& vertex_and_neigh_and_count : adj_list_) {
    long int degree = getDegree(vertex_and_neigh_and_count.first);
    if (degree > max) {
      max = degree;
    }
  }
  return max;
}

long int EdgeContainer::getDegree(const long int vertex) const {
  if (!adj_list_.count(vertex)) {
    throw invalid_argument("vertex is not defined");
  }
  long int degree_count = 0;
  if (adj_list_.at(vertex).size() == 0) {
    return degree_count;
  }
  for (const pair<long int, long int>& neighbor_and_count :
       adj_list_.at(vertex)) {
    if (neighbor_and_count.first == vertex) {
      degree_count += neighbor_and_count.second * 2;
    } else {
      degree_count += neighbor_and_count.second;
    }
  }
  return degree_count;
}

vector<long int> EdgeContainer::getVerticesDegree(long int degree) const {
  vector<long int> vertices;
  for (const auto& vertex_and_neigh_and_count : adj_list_) {
    long int degree_count = getDegree(vertex_and_neigh_and_count.first);
    if (degree_count == degree) {
      vertices.push_back(vertex_and_neigh_and_count.first);
    }
  }
  return vertices;
}

bool EdgeContainer::vertexExistWithDegree(long int degree) const {
  for (const auto& vertex_and_neigh_and_count : adj_list_) {
    long int degree_count = getDegree(vertex_and_neigh_and_count.first);
    if (degree_count == degree) {
      return true;
    }
  }
  return false;
}

bool EdgeContainer::edgeExist(const Edge& edge) const {
  if (adj_list_.count(edge.getEndPoint1())) {
    if (adj_list_.at(edge.getEndPoint1()).count(edge.getEndPoint2())) {
      return true;
    }
  }
  if (adj_list_.count(edge.getEndPoint2())) {
    if (adj_list_.at(edge.getEndPoint2()).count(edge.getEndPoint1())) {
      return true;
    }
  }
  return false;
}

bool EdgeContainer::vertexExist(const long int vertex) const {
  return adj_list_.count(vertex);
}

void EdgeContainer::addEdge(Edge edge) {

  long int point1 = edge.getEndPoint1();
  long int point2 = edge.getEndPoint2();
  if (adj_list_[point1].count(point2)) {
    ++adj_list_[point1][point2];
  } else {
    adj_list_[point1][point2] = 1;
  }

  // Do not add the same edge a second time if the points are the same
  if (point1 != point2) {
    if (adj_list_[point2].count(point1)) {
      ++adj_list_[point2][point1];
    } else {
      adj_list_[point2][point1] = 1;
    }
  }
  return;
}

void EdgeContainer::addVertex(long int vertex) {
  assert(adj_list_.count(vertex) == 0 && "Cannot add vertex already exists");
  unordered_map<long int, long int> empty_temp;
  adj_list_[vertex] = empty_temp;
}

vector<long int> EdgeContainer::getVertices() const {
  vector<long int> vertices;
  for (const pair<const long int, unordered_map<int long, int long>>&
           vertex_and_neigh_and_count : adj_list_) {
    vertices.push_back(vertex_and_neigh_and_count.first);
  }
  return vertices;
}

vector<long int> EdgeContainer::getNeighVertices(long int vertex) const {
  vector<long int> neigh_verts;
  if (adj_list_.count(vertex)) {
    for (const pair<long int, long int>& neigh_and_count :
         adj_list_.at(vertex)) {
      neigh_verts.push_back(neigh_and_count.first);
    }
  }
  return neigh_verts;
}

vector<Edge> EdgeContainer::getNeighEdges(long int vertex) const {
  vector<Edge> neigh_edges;
  if (adj_list_.count(vertex)) {
    for (const pair<long int, long int>& neigh_and_count :
         adj_list_.at(vertex)) {
      for (long int count = 0;
           count < adj_list_.at(vertex).at(neigh_and_count.first); ++count) {
        neigh_edges.push_back(Edge(vertex, neigh_and_count.first));
      }
    }
  }
  return neigh_edges;
}

vector<Edge> EdgeContainer::getEdges() const {
  unordered_map<Edge, long int> extra_edge_count;
  for (const auto& vertex_and_neigh_and_count : adj_list_) {
    for (const pair<long int, long int>& neigh_and_count :
         vertex_and_neigh_and_count.second) {
      extra_edge_count[Edge(vertex_and_neigh_and_count.first,
                            neigh_and_count.first)] = neigh_and_count.second;
    }
  }
  vector<Edge> edges;
  for (pair<const Edge, long int>& edge_count : extra_edge_count) {
    for (long int count = 0; count < edge_count.second; ++count) {
      edges.push_back(edge_count.first);
    }
  }
  return edges;
}

ostream& operator<<(ostream& os, const EdgeContainer edgecontainer) {
  auto edges = edgecontainer.getEdges();
  for (auto edge : edges) {
    os << edge << endl;
  }
  return os;
}

}  // namespace tools
}  // namespace votca
