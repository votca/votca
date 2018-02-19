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

GraphNode::GraphNode(const unordered_map<string,int> int_vals,
                     const unordered_map<string,double> double_vals,
                     const unordered_map<string,string> str_vals){

  vector<string> keys;
  // Grab integer keys sort alphabetically and store in type
  for(auto it : int_vals) keys.push_back(it->first);
  sort(keys.begin(),keys.end());
  for(auto key : keys ){
    type_.append(key); 
    type_.append(lexical_cast<string>(int_vals));
  }
  key.clear();

  // Grab double keys sort alphabetically and store in type
  for(auto it : double_vals) keys.push_back(it->first);
   sort(keys.begin(),keys.end());
  for(auto key : keys ){
    type_.append(key); 
    type_.append(lexical_cast<string>(int_vals));
  }
  key.clear();

  for(auto it : str_vals) keys.push_back(it->first);
  sort(keys.begin(),keys.end());
  for(auto key : keys ){
    type_.append(key); 
    type_.append(lexical_cast<string>(int_vals));
  }
}

bool GraphNode::operator!=(const GraphNode gn) const { 
  return (type_.compare(gn.type_)!=0);
}

bool GraphNode::operator==(const GraphNode gn) const {
  return !((*this)!=gn);
}

}
}
