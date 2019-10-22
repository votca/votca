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
#include <iostream>
#include <vector>
#include <votca/tools/edge.h>

namespace votca {
namespace tools {

using namespace std;

Edge::Edge(long int ID1, long int ID2) {
  vertices_ = vector<long int>{min({ID1, ID2}), max({ID1, ID2})};
}

int Edge::getOtherEndPoint(int ver) const {
  if (ver == vertices_.front()) {
    return vertices_.back();
  } else {
    return vertices_.front();
  }
}

bool Edge::contains(int ID) const {
  return (vertices_.front() == ID || vertices_.back() == ID);
}

bool Edge::operator==(const Edge& ed) const {
  if (this->vertices_.front() == ed.vertices_.front() &&
      this->vertices_.back() == ed.vertices_.back()) {
    return true;
  }
  if (this->vertices_.back() == ed.vertices_.front() &&
      this->vertices_.front() == ed.vertices_.back()) {
    return true;
  }
  return false;
}

bool Edge::operator!=(const Edge& ed) const { return !(*this == ed); }

bool Edge::operator<(const Edge& ed) const {
  if (this->vertices_.front() < ed.vertices_.front()) {
    return true;
  }
  if (this->vertices_.front() > ed.vertices_.front()) {
    return false;
  }
  if (this->vertices_.back() < ed.vertices_.back()) {
    return true;
  }
  return false;
}

bool Edge::operator<=(const Edge& ed) const {
  return (*this < ed || *this == ed);
}

bool Edge::operator>(const Edge& ed) const { return !(*this <= ed); }

bool Edge::operator>=(const Edge& ed) const { return !(*this < ed); }

ostream& operator<<(ostream& os, const Edge& ed) {
  os << "Vertices" << endl;
  os << ed.vertices_.front() << " " << ed.vertices_.back() << endl;
  return os;
}
}  // namespace tools
}  // namespace votca
