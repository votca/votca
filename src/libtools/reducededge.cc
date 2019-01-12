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

#include <algorithm>
#include <votca/tools/reducededge.h>

namespace votca {
namespace tools {

using namespace std;

ReducedEdge::ReducedEdge(vector<int> vertices) {
  assert(vertices.size()>=2 && "Edge vertices must consist of at least two vertices.");
  // Smallest value always placed at the front of the vector

  // If we have a loop will arange the edge so that the smallest valued vertex 
  // sits at both ends
  if(vertices.back()==vertices.front()){
    int min_vertex = vertices.at(0);
    size_t min_index =0;
    for(size_t index=1;index<(vertices.size()-1);++index){
      if(vertices.at(index)<min_vertex){
        min_index = index;
        min_vertex = vertices.at(index);
      }
    }
    vector<int> temp_vertices;
    for(size_t count=0;count<vertices.size();++count){
      size_t index = (min_index+count) % (vertices.size()-1);
      temp_vertices.push_back(vertices.at(index));
    }
    vertices = temp_vertices;
  }
  bool reverse = false;
  size_t length = vertices.size();
  size_t max_index = length/2;
  for(auto index = 0; index< max_index;++index){
    if(vertices.at(index)<vertices.at(length-1-index)){
      break;
    }
    if(vertices.at(index)>vertices.at(length-1-index)){
      reverse = true;
      break;
    }
  }

  if(!reverse){


    vertices_=vertices;
  }else{

    for(auto it=vertices.rbegin();it!=vertices.rend();++it){
      vertices_.push_back(*it);
    }
  }


}

vector<Edge> ReducedEdge::expand() const {
  vector<Edge> edges;
  for(size_t index = 0; index<(vertices_.size()-1);++index){
    Edge ed(vertices_.at(index),vertices_.at(index+1));
    edges.push_back(ed);
  }
  return edges;
}

bool ReducedEdge::vertexExistInChain(int vertex){
  auto it = find(vertices_.begin(),vertices_.end(),vertex);
  return it!=vertices_.end();
}

bool ReducedEdge::operator==(const ReducedEdge ed) const {
  if(ed.vertices_.size()!=vertices_.size()) return false;
  for(auto index=0;index<vertices_.size();++index){
    if(vertices_.at(index)!=ed.vertices_.at(index)) return false;
  }
  return true;
}

bool ReducedEdge::operator!=(const ReducedEdge ed) const { return !(*this == ed); }

bool ReducedEdge::operator<(const ReducedEdge ed) const {
  if (this->vertices_.front() < ed.vertices_.front()) return true;
  if (this->vertices_.front() > ed.vertices_.front()) return false;
  if (this->vertices_.back() < ed.vertices_.back()) return true;

  if(vertices_.size()<ed.vertices_.size()) return true;

  for(auto index=0;index<vertices_.size();++index){
    if(vertices_.at(index)>ed.vertices_.at(index)) return false;
  }
  if(*this==ed) return false;

  return true;
}

bool ReducedEdge::operator<=(const ReducedEdge ed) const {
  return (*this < ed || *this == ed);
}

bool ReducedEdge::operator>(const ReducedEdge ed) const { return !(*this <= ed); }

bool ReducedEdge::operator>=(const ReducedEdge ed) const { return !(*this < ed); }

ostream& operator<<(ostream& os, const ReducedEdge ed) {
  os << "Vertices" << endl;
  for(auto vertex : ed.vertices_){
    os << vertex << " ";
  }
  os << endl;
  return os;
}
}
}
