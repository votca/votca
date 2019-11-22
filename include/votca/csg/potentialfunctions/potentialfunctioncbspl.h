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

#ifndef POTENTIALFUNCTIONCBSPL_H
#define POTENTIALFUNCTIONCBSPL_H

#include "potentialfunction.h"
#include <math.h>
#include <votca/tools/table.h>

namespace votca {
namespace csg {

class PotentialFunctionCBSPL : public PotentialFunction {
 public:
  PotentialFunctionCBSPL(const std::string &name, const Index nlam,
                         const double min = 0.0, const double max = 10.0);
  ~PotentialFunctionCBSPL() override = default;
  // calculate function value for given r
  double CalculateF(const double r) const override;
  // calculate first derivative w.r.t. ith parameter
  double CalculateDF(const Index i, const double r) const override;
  // calculate second derivative w.r.t. ith parameter
  double CalculateD2F(const Index i, const Index j,
                      const double r) const override;

  Index getOptParamSize() const override;

  void setParam(std::string filename) override;

  void SaveParam(const std::string &filename) override;

  void SavePotTab(const std::string &filename, const double step) override;

  void SavePotTab(const std::string &filename, const double step,
                  const double rmin, const double rcut) override;
  void setOptParam(const Index i, const double val) override;

  double getOptParam(const Index i) const override;

  void extrapolExclParam();

 protected:
  // exclude these many first coefficients from optimization
  // since the region relevant to these coefficients is not sampled
  // the value of _nexcl is determined from rmin
  Index _nexcl;
  // fix these many coeff near the cut-off to zero to ensure
  // zero potential and force values near cut-off
  Index _ncutcoeff;

  Index _nbreak;
  double _dr;
  Eigen::VectorXd _rbreak;

  Eigen::Matrix4d _M;
};
}  // namespace csg
}  // namespace votca
#endif
