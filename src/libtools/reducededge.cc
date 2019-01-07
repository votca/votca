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

#include <votca/tools/reducededge.h>

namespace votca {
namespace tools {

using namespace std;

ReducedEdge::ReducedEdge(vector<int> chain) {
  if(chain.size()<2){
    throw invalid_argument("Edge chain must consist of at least two vertices.");
  }
  // Smallest value always placed at the front of the vectolr

  bool reverse = false;
  size_t length = chain.size();
  size_t max_index = length/2;
  for(auto index = 0; index< max_index;++index){
    if(chain.at(index)<chain.at(length-1-index)){
      break;
    }
    if(chain.at(index)>chain.at(length-1-index)){
      reverse = true;
      break;
    }
  }
  if(!reverse){
    chain_=chain;
  }else{
    for(auto it=chain.rbegin();it!=chain.rend();++it){
      chain_.push_back(*it);
    }
  }
}

bool ReducedEdge::operator==(const ReducedEdge ed) const {
  if(ed.vertices_.size()!=vertices_.size()) return false;
  for(auto index=0;index<verties_.size();++index){
    if(vertices_.at(index)!=ed.verties_.at(index)) return false;
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
    if(verties_.at(index)>ed.vertices_.at(index)) return false;
  }

  return true;
}

bool ReducedEdge::operator<=(const ReducedEdge ed) const {
  return (*this < ed || *this == ed);
}

bool ReducedEdge::operator>(const ReducedEdge ed) const { return !(*this <= ed); }

bool ReducedEdge::operator>=(const ReducedEdge ed) const { return !(*this < ed); }

ostream& operator<<(ostream& os, const ReducedEdge ed) {
  os << "Vertices" << endl;
  for(auto vertex : vertices_){
    os << vertex << " ";
  }
  os << endl;
  return os;
}
}
}
