/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
#define	POTENTIALFUNCTIONCBSPL_H

#include <votca/tools/table.h>
#include <math.h>
#include "potentialfunction.h"

using namespace std;
using namespace votca::tools;

class PotentialFunctionCBSPL : public PotentialFunction {
 public:
  PotentialFunctionCBSPL(const string& name_,const int nlam_,
                         const double min_=0.0, const double max_=10.0);
  ~PotentialFunctionCBSPL(){}
  // calculate function value for given r
  double CalculateF (const double r) const;
  // calculate first derivative w.r.t. ith parameter
  double CalculateDF(const int i, const double r) const;
  // calculate second derivative w.r.t. ith parameter
  double CalculateD2F(const int i, const int j, const double r) const;

  int getOptParamSize() const ;

  void setParam(string filename);

  void SaveParam(const string& filename);

  void SavePotTab(const string& filename, const double step);

  void SavePotTab(const string& filename, const double step, const double rmin, const double rcut);
  void setOptParam(const int i, const double val);

  double getOptParam(const int i) const;

  void extrapolExclParam();

 protected:

  // exclude these many first coefficients from optimization
  // since the region relevant to these coefficients is not sampled
  // the value of _nexcl is determined from rmin
  int _nexcl;
  // fix these many coeff near the cut-off to zero to ensure
  // zero potential and force values near cut-off
  int _ncutcoeff;

  int _nbreak;
  double _dr;
  Eigen::VectorXd _rbreak;

  Eigen::MatrixXd _M;


};

#endif
