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

#include <set>
#include <vector>
#include <exception>
#include <votca/tools/edge.h>
#include <votca/tools/edgecontainer.h>
#include <algorithm>

namespace votca {
  namespace tools {

using namespace std;

EdgeContainer::EdgeContainer(Edge ed) { addEdge(ed); }

EdgeContainer::EdgeContainer(vector<Edge> eds) {
  for (auto ed : eds) {
    addEdge(ed);
  }
}

int EdgeContainer::getMaxDegree(void) const{
  int max = 0;
  for(auto const& it : adj_list_) {
    int degree = getDegree(it.first);
    if(degree>max) {
      max = degree;
    }
  }
  return max;
}

int EdgeContainer::getDegree(const int vert) const{
  if(!adj_list_.count(vert)) throw invalid_argument("vertex is not defined");
  int degree_count=0;
  for(auto neighbor_and_count : adj_list_[vert]){
    degree_count+=neighbor_and_count.second;
  }
  return degree_count;
}

vector<int> EdgeContainer::getVerticesDegree(int degree) const{
  vector<int> verts;
  for(auto v_list : adj_list_){
    int degree_count= getDegree(v_list.first);
    if(degree_count==degree){
      verts.push_back(v_list.first);
    }
  }
  return verts;
}

bool EdgeContainer::edgeExist(Edge ed) {
  return adj_list_[ed.getEndPoint1()].count(ed.getEndPoint2()) || 
    adj_list_[ed.getEndPoint2()].count(ed.getEndPoint1()); 
}

bool EdgeContainer::vertexExist(int vert) { return adj_list_.count(vert); }

void EdgeContainer::addEdge(Edge ed) {
  int point1 = ed.getEndPoint1();
  int point2 = ed.getEndPoint2();
  if(adj_list_[point1].count(point2)){
    ++adj_list_[point1][point2];
  }else{
    adj_list_[point1][point2]=1;
  }
  if(adj_list_[point2].count(point1)){
    ++adj_list_[point2][point1];
  }else{
    adj_list_[point2][point1]=1;
  }
  return;
}

vector<int> EdgeContainer::getVertices() {
  vector<int> vertices;
  for (auto const& it : adj_list_) vertices.push_back(it.first);
  return vertices;
}

vector<int> EdgeContainer::getNeighVertices(int vert) {
  vector<int> neigh_verts;
  for (auto const& neigh_vert : adj_list_[vert]) {
    neigh_verts.push_back(neigh_vert.first);
  }
  return neigh_verts;
}

vector<Edge> EdgeContainer::getNeighEdges(int vert) {
  vector<Edge> neigh_edges;
  for (auto const& neigh_vert : adj_list_[vert]) {
    for(int count=0;count<adj_list_[vert][neigh_vert.first];++count){
      neigh_edges.push_back(Edge(vert, neigh_vert.first));
    }
  }
  return neigh_edges;
}

vector<Edge> EdgeContainer::getEdges() const {
  unordered_map<Edge,int> extra_edge_count;
  for (auto const& it : adj_list_) {
    for (auto const& vert : it.second) {
      extra_edge_count[Edge(it.first,vert.first)]=vert.second;
    }
  }
  vector<Edge> vec_edgs;
  for(auto edge_count : extra_edge_count){
    for(int count = 0; count < edge_count.second;++count){
      vec_edgs.push_back(edge_count.first);
    }
  }
  return vec_edgs;
}

ostream& operator<<(ostream& os, const EdgeContainer edgecontainer){
  auto edges = edgecontainer.getEdges();
  for(auto edge : edges){
    os << edge << endl;
  }
  return os;
}

}
}
