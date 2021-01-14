/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_TOOLS_RANDOM_H
#define VOTCA_TOOLS_RANDOM_H

// Standard includes
#include <random>

// Local VOTCA includes
#include "types.h"

namespace votca {
namespace tools {

class Random {
 public:
  void init(Index seed) {
    if (seed < 0) {
      throw std::runtime_error("seed integer must be positive.");
    }
    _mt = std::mt19937(unsigned(seed));
  }
  // draws a random double from [0,1)
  double rand_uniform() { return _distribution(_mt); }
  // sets maxint for a uniform integer distribution [0,maxint]
  void setMaxInt(Index maxint) {
    _int_distribution = std::uniform_int_distribution<Index>{0, maxint};
  }
  // draws from a uniform integer distribution [0,maxint]
  Index rand_uniform_int() { return _int_distribution(_mt); }

 private:
  std::mt19937 _mt;
  std::uniform_real_distribution<double> _distribution{0.0, 1.0};
  std::uniform_int_distribution<Index> _int_distribution;
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_RANDOM_H
