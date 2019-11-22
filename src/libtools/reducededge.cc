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
#include <votca/tools/reducededge.h>

namespace votca {
namespace tools {

using namespace std;

/**
 * Given a chain of vertices where the end points are the same this method
 * will rotate the chain so that the vertex with the smallest value is at the
 * end points e.g.
 *
 * 3 - 4 - 5 - 9 - 2 - 8 - 3
 *
 * The method will return the chain as
 *
 * 2 - 8 - 3 - 4 - 5 - 9 - 2
 *
 * Notice that the order is preserved, because it is assumed that the chain
 * makes a loop.
 **/
void moveVertexWithSmallestValueToEnds_(vector<Index>& vertices) {
  vertices.pop_back();
  vector<Index>::iterator starting_iterator = next(vertices.begin());
  Index min_vertex_index =
      min_element(starting_iterator, vertices.end()) - vertices.begin();
  if (vertices.front() > vertices.at(min_vertex_index)) {
    rotate(vertices.begin(), vertices.begin() + min_vertex_index,
           vertices.end());
  }
  vertices.push_back(vertices.front());
}

/**
 * This function will determine if it is necessary to reverse the order of
 * a chain of vertices. E.g.
 *
 * 1 - 2 - 4 - 5
 *
 * This chain will not need to be reversed because the end points have end
 * point 1 smaller than 5
 *
 * This chain will need to be reversed
 *
 * 4 - 1 - 8 - 2
 *
 * Because 2 is smaller than 4
 *
 * If the chain has the same end points the next two vertices will be compared
 *
 * 1 - 2 - 3 - 4 - 1
 *     ^       ^
 *     |       |
 *
 * The 2 is smaller than the 4 so there is no need to reverse the order.
 *
 * 4 - 9 - 8 - 6 - 4
 *
 * In this case th 6 is smaller than the 9 so the order will be reversed.
 *
 **/
bool verticesShouldBeReversed_(vector<Index>& vertices) {
  size_t length = vertices.size();
  size_t max_index = length / 2;
  for (size_t index = 0; index < max_index; ++index) {
    if (vertices.at(index) < vertices.at(length - 1 - index)) {
      return false;
    } else if (vertices.at(index) > vertices.at(length - 1 - index)) {
      return true;
    }
  }
  return false;
}

ReducedEdge::ReducedEdge(vector<Index> vertices) {
  assert(vertices.size() >= 2 &&
         "Edge vertices must consist of at least two vertices.");
  // Smallest value always placed at the front of the vector

  // If we have a loop will arange the edge so that the smallest valued vertex
  // sits at both ends
  if (vertices.back() == vertices.front() && vertices.size() > 2) {
    moveVertexWithSmallestValueToEnds_(vertices);
  }

  if (verticesShouldBeReversed_(vertices)) {
    reverse(vertices.begin(), vertices.end());
  }
  vertices_ = vertices;
}

vector<Edge> ReducedEdge::expand() const {
  vector<Edge> edges;
  for (size_t index = 0; index < (vertices_.size() - 1); ++index) {
    Edge ed(vertices_.at(index), vertices_.at(index + 1));
    edges.push_back(ed);
  }
  return edges;
}

bool ReducedEdge::vertexExistInChain(const int& vertex) const {
  vector<Index>::const_iterator vertex_iterator =
      find(vertices_.begin(), vertices_.end(), vertex);
  return vertex_iterator != vertices_.end();
}

bool ReducedEdge::operator==(const ReducedEdge& edge) const {
  if (edge.vertices_.size() != vertices_.size()) {
    return false;
  }
  for (size_t index = 0; index < vertices_.size(); ++index) {
    if (vertices_.at(index) != edge.vertices_.at(index)) {
      return false;
    }
  }
  return true;
}

bool ReducedEdge::operator!=(const ReducedEdge& edge) const {
  return !(*this == edge);
}

bool ReducedEdge::operator<(const ReducedEdge& edge) const {
  if (this->vertices_.front() < edge.vertices_.front()) {
    return true;
  }
  if (this->vertices_.front() > edge.vertices_.front()) {
    return false;
  }
  if (this->vertices_.back() < edge.vertices_.back()) {
    return true;
  }
  if (vertices_.size() < edge.vertices_.size()) {
    return true;
  }
  for (size_t index = 0; index < vertices_.size(); ++index) {
    if (vertices_.at(index) > edge.vertices_.at(index)) {
      return false;
    }
  }
  if (*this == edge) {
    return false;
  }
  return true;
}

bool ReducedEdge::operator<=(const ReducedEdge& edge) const {
  return (*this < edge || *this == edge);
}

bool ReducedEdge::operator>(const ReducedEdge& edge) const {
  return !(*this <= edge);
}

bool ReducedEdge::operator>=(const ReducedEdge& edge) const {
  return !(*this < edge);
}

ostream& operator<<(ostream& os, const ReducedEdge& edge) {
  os << "Vertices" << endl;
  for (auto vertex : edge.vertices_) {
    os << vertex << " ";
  }
  os << endl;
  return os;
}
}  // namespace tools
}  // namespace votca
