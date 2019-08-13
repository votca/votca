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

/*************************************************
     MARSAGLIA pseudo random number generator
     See: G. Marsaglia and A. Zaman. Toward a universal random number generator,
          Statistics & Probability Letters, 9(1):35–39, 1990.
     This function returns a double precision floating point number
     uniformly distributed in the range [0,1)
*************************************************/

#pragma once
#ifndef _VOTCA_TOOLS_RANDOM2_H_
#define _VOTCA_TOOLS_RANDOM2_H_

#include <vector>

namespace votca {
namespace tools {

/**
  \brief MARSAGLIA pseudo random number generator

  This class generates double precision floating point numbers
  uniformly distributed in the range [0,1)
*/

class Random2 {
 public:
  Random2(){};
  ~Random2(){};

  void init(int nA1, int nA2, int nA3, int nB1);

  double rand_uniform(void);
  int rand_uniform_int(int max_int);
  double rand_gaussian(double sigma);

 private:
  static const int MARS_FIELD_SIZE = 98;
  std::vector<double> MARSarray;
  double MARSc, MARScd, MARScm;
  int MARSi, MARSj;
};

}  // namespace tools
}  // namespace votca

#endif /* _RANMARS2_H_ */
