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

#include <vector>
#include <votca/tools/edge.h>
#include <votca/tools/edgecontainer.h>
#include <set>
#include <algorithm>

using namespace votca::tools;
using namespace std;

EdgeContainer::EdgeContainer(Edge ed) {
  if (ed.getV1() != ed.getV2()) {
    vertices_.push_back(ed.getV1());
    vertices_.push_back(ed.getV2());
  }
  vertices_.push_back(ed.getV1());

  edges_.push_back(ed);
}

EdgeContainer::EdgeContainer(vector<Edge> eds) {
  set<int> f_verts;
  set<Edge> f_eds;
  for (auto ed : eds) {
    f_verts.insert(ed.getV1());
    f_verts.insert(ed.getV2());
    f_eds.insert(ed);
  }

  for (auto vert : f_verts) vertices_.push_back(vert);
  for (auto ed : f_eds) edges_.push_back(ed);
}

bool EdgeContainer::edgeExist(Edge ed) {
  return (find(edges_.begin(), edges_.end(), ed) != edges_.end());
}

bool EdgeContainer::vertexExist(int vert) {
  return (find(vertices_.begin(), vertices_.end(), vert) != vertices_.end());
}

void EdgeContainer::addEdge(Edge ed) {
  if (!edgeExist(ed)) edges_.push_back(ed);
  if (!vertexExist(ed.getV1())) vertices_.push_back(ed.getV1());
  if (!vertexExist(ed.getV2())) vertices_.push_back(ed.getV2());
  return;
}

vector<int> EdgeContainer::getVertices() { return vertices_; }

vector<int> EdgeContainer::getNeighVertices(int vert) {
  vector<int> neigh_verts;
  for (auto ed : edges_) {
    if (ed.contains(vert)) {
      neigh_verts.push_back(ed.getOtherV(vert));
    }
  }
  return neigh_verts;
}

vector<Edge> EdgeContainer::getNeighEdges(int vert) {
  vector<Edge> neigh_edges;
  for (auto ed : edges_) {
    if (ed.contains(vert)) neigh_edges.push_back(ed);
  }
  return neigh_edges;
}

vector<Edge> EdgeContainer::getEdges() { return edges_; }
