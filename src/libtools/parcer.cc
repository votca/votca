/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <iostream>
#include <string>
#include <vector>

namespace votca {
namespace tools {

using namespace std;

// wrapper written with Denis to do MOO only on a file which give the unitary
// transformation + displacement one per line
void parce_string(string line, string delims, vector<string>* result) {
  string::size_type begIdx, endIdx;

  begIdx = line.find_first_not_of(delims);
  while (begIdx != string::npos) {
    endIdx = line.find_first_of(delims, begIdx);
    if (endIdx == string::npos) {
      endIdx = line.length();
    }
    result->push_back(line.substr(begIdx, endIdx - begIdx));
    if (endIdx == line.length()) {
      break;
      cout << "I am still here";
    }
    begIdx = line.find_first_not_of(delims, endIdx);
  }
}

}  // namespace tools
}  // namespace votca