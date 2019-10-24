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

#include <votca/csg/potentialfunctions/potentialfunctionlj126.h>

using namespace std;

namespace votca {
namespace csg {
PotentialFunctionLJ126::PotentialFunctionLJ126(const string &name_,
                                               const double min_,
                                               const double max_)
    : PotentialFunction(name_, 2, min_, max_) {}

double PotentialFunctionLJ126::CalculateF(const double r) const {

  if (r >= _min && r <= _cut_off) {
    return _lam(0) / pow(r, 12) - _lam(1) / pow(r, 6);
  } else {
    return 0.0;
  }
}

// calculate first derivative w.r.t. ith parameter
double PotentialFunctionLJ126::CalculateDF(const long int i,
                                           const double r) const {

  if (r >= _min && r <= _cut_off) {

    switch (i) {
      case 0:
        return 1.0 / pow(r, 12);
      case 1:
        return -1.0 / pow(r, 6);
    }
  }
  return 0.0;
}

// calculate second derivative w.r.t. ith and jth parameters
double PotentialFunctionLJ126::CalculateD2F(const long int i, const long int j,
                                            const double r) const {

  return 0.0;
}
}  // namespace csg
}  // namespace votca
