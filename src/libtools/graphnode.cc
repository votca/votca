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

#include <iostream>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>
#include <votca/tools/graphnode.h>

namespace votca {
namespace tools {

using namespace std;

GraphNode::GraphNode(size_t id,
                     const unordered_map<string,int> int_vals,
                     const unordered_map<string,double> double_vals,
                     const unordered_map<string,string> str_vals){
    
  vector<string> keys;
  for(auto it : int_vals) keys.push_back(it->first);
  for(auto it : double_vals) keys.push_back(it->first);
  for(auto it : str_vals) keys.push_back(it->first);
  sort(keys.begin(),keys.end());
  for(auto key : keys) type_.append(key);
  setId(id);
}

bool GraphNode::operator!=(const GraphNode gn) const { 
  if(type 

}

bool GraphNode::operator==(const GraphNode gn) const {
  
  return false;
}


bool Edge::operator<(const Edge ed) const {
  if (this->vertices.first < ed.vertices.first) return true;
  if (this->vertices.first > ed.vertices.first) return false;
  if (this->vertices.second < ed.vertices.second) return true;
  return false;
}

ostream& operator<<(ostream& os, const Edge ed) {
  os << "Vertices" << endl;
  os << ed.vertices.first << " " << ed.vertices.second << endl;
  return os;
}
}
}
