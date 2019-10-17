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

#include <votca/csg/potentialfunctions/potentialfunction.h>
#include <votca/tools/table.h>

#include <boost/lexical_cast.hpp>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace csg {

PotentialFunction::PotentialFunction(const string &name_, const int nlam_,
                                     const double min_, const double max_) {

  _name = name_;
  _lam = Eigen::VectorXd::Zero(nlam_);
  _min = min_;
  _cut_off = max_;
}

void PotentialFunction::setParam(string filename) {

  Table param;
  param.Load(filename);

  if (param.size() != _lam.size()) {

    throw std::runtime_error("In potential " + _name +
                             ": parameters size mismatch!\n"
                             "Check input parameter file \"" +
                             filename + "\" \nThere should be " +
                             boost::lexical_cast<string>(_lam.size()) +
                             " parameters");
  } else {
    for (unsigned int i = 0; i < _lam.size(); i++) {
      _lam(i) = param.y(i);
    }
  }
}

void PotentialFunction::SaveParam(const string &filename) {

  Table param;
  param.SetHasYErr(false);
  param.resize(_lam.size());

  for (unsigned int i = 0; i < _lam.size(); i++) {
    param.set(i, i, _lam(i), 'i');
  }

  param.Save(filename);
}

void PotentialFunction::SavePotTab(const string &filename, const double step) {
  int ngrid = (int)((_cut_off - _min) / step + 1.00000001);
  Table pot_tab;
  pot_tab.SetHasYErr(false);
  pot_tab.resize(ngrid);
  double r_init;
  int i;

  for (r_init = _min, i = 0; i < ngrid - 1; r_init += step) {
    pot_tab.set(i++, r_init, CalculateF(r_init), 'i');
  }

  pot_tab.set(i, _cut_off, CalculateF(_cut_off), 'i');
  pot_tab.Save(filename);
}

void PotentialFunction::SavePotTab(const string &filename, const double step,
                                   const double rmin, const double rcut) {
  int ngrid = (int)((rcut - rmin) / step + 1.00000001);
  Table pot_tab;
  pot_tab.SetHasYErr(false);
  pot_tab.resize(ngrid);
  double r_init;
  int i;
  char flag = 'i';

  for (r_init = rmin, i = 0; i < ngrid - 1; r_init += step) {
    pot_tab.set(i++, r_init, CalculateF(r_init), flag);
  }

  pot_tab.set(i, rcut, CalculateF(rcut), flag);
  pot_tab.Save(filename);
}
}  // namespace csg
}  // namespace votca
