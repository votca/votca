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

#ifndef _VOTCA_TOOLS_RANDOM2_H_
#define _VOTCA_TOOLS_RANDOM2_H_

#include <random>
namespace votca {
namespace tools {

class Random2 {
 public:
  void init(int seed) { _mt = std::mt19937(seed); }
  // draws a random double from [0,1)
  double rand_uniform() { return _distribution(_mt); }
  // sets maxint for a uniform integer distribution [0,maxint]
  void setMaxInt(int maxint) {
    _int_distribution = std::uniform_int_distribution<int>{0, maxint};
  }
  // draws from a uniform integer distribution [0,maxint]
  int rand_uniform_int() { return _int_distribution(_mt); }

 private:
  std::mt19937 _mt;
  std::uniform_real_distribution<double> _distribution{0.0, 1.0};
  std::uniform_int_distribution<int> _int_distribution;
};

}  // namespace tools
}  // namespace votca

#endif /* _RANMARS2_H_ */
