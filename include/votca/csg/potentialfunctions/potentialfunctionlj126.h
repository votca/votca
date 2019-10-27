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

#ifndef POTENTIALFUNCTIONLJ126_H
#define POTENTIALFUNCTIONLJ126_H

#include "potentialfunction.h"
#include <math.h>
#include <votca/tools/table.h>

namespace votca {
namespace csg {
// LJ 12-6 potential class
// with c12,c6 parameters
class PotentialFunctionLJ126 : public PotentialFunction {
 public:
  PotentialFunctionLJ126(const std::string &name, double min = 0.0,
                         double max = 10.0);
  ~PotentialFunctionLJ126() override = default;
  // calculate function value for given r
  double CalculateF(double r) const override;
  // calculate first derivative w.r.t. ith parameter
  double CalculateDF(long i, double r) const override;
  // calculate second derivative w.r.t. ith parameter
  double CalculateD2F(long i, long j, double r) const override;
};
}  // namespace csg
}  // namespace votca
#endif
