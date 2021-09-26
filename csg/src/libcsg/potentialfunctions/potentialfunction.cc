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

#include "../../../include/votca/csg/potentialfunctions/potentialfunction.h"
#include <votca/tools/table.h>

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace csg {

PotentialFunction::PotentialFunction(const string &name, Index nlam, double min,
                                     double max) {

  name_ = name;
  lam_ = Eigen::VectorXd::Zero(nlam);
  min_ = min;
  cut_off_ = max;
}

void PotentialFunction::setParam(string filename) {

  Table param;
  param.Load(filename);

  if (param.size() != lam_.size()) {

    throw std::runtime_error("In potential " + name_ +
                             ": parameters size mismatch!\n"
                             "Check input parameter file \"" +
                             filename + "\" \nThere should be " +
                             boost::lexical_cast<string>(lam_.size()) +
                             " parameters");
  } else {
    for (Index i = 0; i < lam_.size(); i++) {
      lam_(i) = param.y(i);
    }
  }
}

void PotentialFunction::SaveParam(const string &filename) {

  Table param;
  param.SetHasYErr(false);
  param.resize(lam_.size());

  for (Index i = 0; i < lam_.size(); i++) {
    param.set(i, double(i), lam_(i), 'i');
  }

  param.Save(filename);
}

void PotentialFunction::SavePotTab(const string &filename, double step) {
  Index ngrid = (Index)((cut_off_ - min_) / step + 1.00000001);
  Table pot_tab;
  pot_tab.SetHasYErr(false);
  pot_tab.resize(ngrid);
  double r_init;
  Index i;

  for (r_init = min_, i = 0; i < ngrid - 1; r_init += step) {
    pot_tab.set(i++, r_init, CalculateF(r_init), 'i');
  }

  pot_tab.set(i, cut_off_, CalculateF(cut_off_), 'i');
  pot_tab.Save(filename);
}

void PotentialFunction::SavePotTab(const string &filename, double step,
                                   double rmin, double rcut) {
  Index ngrid = (Index)((rcut - rmin) / step + 1.00000001);
  Table pot_tab;
  pot_tab.SetHasYErr(false);
  pot_tab.resize(ngrid);
  double r_init;
  Index i;
  char flag = 'i';

  for (r_init = rmin, i = 0; i < ngrid - 1; r_init += step) {
    pot_tab.set(i++, r_init, CalculateF(r_init), flag);
  }

  pot_tab.set(i, rcut, CalculateF(rcut), flag);
  pot_tab.Save(filename);
}
}  // namespace csg
}  // namespace votca
