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

#include <votca/tools/spline.h>

namespace votca {
namespace tools {

using namespace std;

long Spline::GenerateGrid(double min, double max, double h) {
  int vec_size = (int)((max - min) / h + 1.00000001);
  _r.resize(vec_size);
  int i;

  double r_init;

  for (r_init = min, i = 0; i < vec_size - 1; r_init += h) {
    _r[i++] = r_init;
  }
  _r[i] = max;
  _f.resize(_r.size());
  _f2.resize(_r.size());
  return _r.size();
}

}  // namespace tools
}  // namespace votca
