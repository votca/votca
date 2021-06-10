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

#ifndef POTENTIALFUNCTION_H
#define POTENTIALFUNCTION_H

#include <Eigen/Dense>
#include <string>
#include <votca/tools/types.h>

namespace votca {
namespace csg {
class PotentialFunction {
 public:
  virtual ~PotentialFunction() = default;
  // read parameters from the input file
  virtual void setParam(std::string filename);
  // save parameters to the file
  virtual void SaveParam(const std::string &filename);
  // write potential table
  virtual void SavePotTab(const std::string &filename, double step);
  // write potential table for specified interval
  virtual void SavePotTab(const std::string &filename, double step, double rmin,
                          double rcut);
  // set all parameters
  void setParam(const Eigen::VectorXd &param) { lam_ = param; }
  // set ith parameter
  void setParam(Index i, double val) { lam_(i) = val; }
  // set ith parameter among those to be optimized
  virtual void setOptParam(Index i, double val) { setParam(i, val); }
  // set minimum r value to avoid large values
  void setMinDist(double min) { min_ = min; }
  // set cut-off value
  void setCutOffDist(double cutoff) { cut_off_ = cutoff; }
  // calculate function
  virtual double CalculateF(double r) const = 0;
  // calculate first derivative w.r.t. ith parameter
  virtual double CalculateDF(Index i, double r) const = 0;
  // calculate second derivative w.r.t. ith parameter
  virtual double CalculateD2F(Index i, Index j, double r) const = 0;
  // return parameter
  Eigen::VectorXd &Params() { return lam_; }
  // return ith parameter
  double getParam(Index i) const { return lam_(i); }
  // return ith parameter among those to be optimized
  virtual double getOptParam(Index i) const { return getParam(i); }
  // return size of parameters
  Index getParamSize() const { return lam_.size(); }
  // return size of parameters to be optimized
  virtual Index getOptParamSize() const { return getParamSize(); }
  // return cut-off value
  double getCutOff() const { return cut_off_; }
  double getMinDist() const { return min_; }

 protected:
  PotentialFunction(const std::string &name, Index nlam, double min,
                    double max);

  std::string name_;
  Eigen::VectorXd lam_;
  double cut_off_;
  double min_;
};
}  // namespace csg
}  // namespace votca
#endif
