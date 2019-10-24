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

#include <boost/lexical_cast.hpp>
#include <iostream>
#include <votca/csg/potentialfunctions/potentialfunctioncbspl.h>
#include <votca/tools/table.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace csg {

PotentialFunctionCBSPL::PotentialFunctionCBSPL(const string &name_,
                                               const int nlam_,
                                               const double min_,
                                               const double max_)
    : PotentialFunction(name_, nlam_, min_, max_) {

  /* Here nlam_ is the total number of coeff values that are to be optimized
   * To ensure that potential and force go to zero smoothly near cut-off,
   * as suggested in Ref. PCCP, 11, 1901, 2009, coeff values leading up to
   * cut-off and beyond take a value of zero.
   *
   * Since region less than rmin is not sampled sufficiently for stability
   * first _nexcl coefficients are not optimized instead their values are
   * extrapolated from first statistically significant knot values near rmin
   */

  long int nknots;

  nknots = _lam.size();

  _nbreak = nknots - 2;

  _dr = (_cut_off) / (double(_nbreak - 1));

  // break point locations
  // since ncoeff = nbreak +2 , r values for last two coefficients are also
  // computed
  _rbreak = Eigen::VectorXd::Zero(nknots);

  for (int i = 0; i < nknots; i++) {
    _rbreak(i) = i * _dr;
  }

  // exclude knots corresponding to r <= _min
  _nexcl = min((long int)(_min / _dr), _nbreak - 2) + 1;

  // account for finite numerical division of _min/_dr
  // e.g. 0.24/0.02 may result in 11.99999999999999
  if (_rbreak(_nexcl) == _min) {
    _nexcl++;
  }

  // fixing last 4 knots to zeros is reasonable
  _ncutcoeff = 4;

  // check if we have enough parameters to optimize
  if ((int(_lam.size()) - _nexcl - _ncutcoeff) < 1) {
    throw std::runtime_error(
        "In potential " + _name +
        ": no parameters to optimize!\n"
        "All the knot values fall in the range of either excluded (due to high "
        "repulsive region) or cut-off region.\n"
        "This issue can be resolved by one or combination of following steps:\n"
        "1. Make sure you are using large-enough cut-off for this CG "
        "potential.\n"
        "2. Make sure the CG-MD runs are sufficiently long and CG-MD RDF are "
        "statistically reliable.\n"
        "3. Use more knot values.\n");
  }

  _M = Eigen::MatrixXd::Zero(4, 4);
  _M(0, 0) = 1.0;
  _M(0, 1) = 4.0;
  _M(0, 2) = 1.0;
  _M(0, 3) = 0.0;
  _M(1, 0) = -3.0;
  _M(1, 1) = 0.0;
  _M(1, 2) = 3.0;
  _M(1, 3) = 0.0;
  _M(2, 0) = 3.0;
  _M(2, 1) = -6.0;
  _M(2, 2) = 3.0;
  _M(2, 3) = 0.0;
  _M(3, 0) = -1.0;
  _M(3, 1) = 3.0;
  _M(3, 2) = -3.0;
  _M(3, 3) = 1.0;
  _M /= 6.0;
}

long int PotentialFunctionCBSPL::getOptParamSize() const {

  return _lam.size() - _nexcl - _ncutcoeff;
}

void PotentialFunctionCBSPL::setParam(string filename) {

  Table param;
  param.Load(filename);
  _lam.setZero();

  if (param.size() != _lam.size()) {

    throw std::runtime_error("In potential " + _name +
                             ": parameters size mismatch!\n"
                             "Check input parameter file \"" +
                             filename + "\" \nThere should be " +
                             boost::lexical_cast<string>(_lam.size()) +
                             " parameters");
  } else {

    // force last _ncutcoeff to zero
    for (unsigned int i = 0; i < _lam.size() - _ncutcoeff; i++) {
      _lam(i) = param.y(i);
    }
  }
}

void PotentialFunctionCBSPL::SaveParam(const string &filename) {

  extrapolExclParam();

  Table param;
  param.SetHasYErr(false);
  param.resize(_lam.size());

  // write extrapolated knots with flag 'o'
  // points close to rmin can also be stastically not reliable
  // so flag 3 more points next to rmin as 'o'
  for (int i = 0; i < _nexcl + 3; i++) {
    param.set(i, _rbreak(i), _lam(i), 'o');
  }

  for (long int i = _nexcl + 3; i < _lam.size(); i++) {
    param.set(i, _rbreak(i), _lam(i), 'i');
  }

  param.Save(filename);
}

void PotentialFunctionCBSPL::SavePotTab(const string &filename,
                                        const double step, const double rmin,
                                        const double rcut) {
  extrapolExclParam();
  PotentialFunction::SavePotTab(filename, step, rmin, rcut);
}

void PotentialFunctionCBSPL::SavePotTab(const string &filename,
                                        const double step) {
  extrapolExclParam();
  PotentialFunction::SavePotTab(filename, step);
}

void PotentialFunctionCBSPL::extrapolExclParam() {

  double u0 = _lam(_nexcl);
  double m = (_lam(_nexcl + 1) - _lam(_nexcl)) /
             (_rbreak(_nexcl + 1) - _rbreak(_nexcl));
  double r0 = _rbreak(_nexcl);

  /* If the slope m is positive then the potential core
   * will be attractive. So, artificially forcing core to be
   * repulsive by setting m = -m
   */
  if (m > 0) {
    cout << _name << " potential's extrapolated core is attractive!" << endl;
    cout << "Artifically enforcing repulsive core.\n" << endl;
    m *= -1.0;
  }
  // using linear extrapolation
  // u(r) = ar + b
  // a = m
  // b = - m*r0 + u0
  // m = (u1-u0)/(r1-r0)

  double a = m;
  double b = -1.0 * m * r0 + u0;
  for (int i = 0; i < _nexcl; i++) {
    _lam(i) = a * _rbreak(i) + b;
  }
}

void PotentialFunctionCBSPL::setOptParam(const long int i, const double val) {

  _lam(i + _nexcl) = val;
}

double PotentialFunctionCBSPL::getOptParam(const long int i) const {

  return _lam(i + _nexcl);
}

double PotentialFunctionCBSPL::CalculateF(const double r) const {

  if (r <= _cut_off) {

    double u = 0.0;
    long int indx = min((long int)(r / _dr), _nbreak - 2);
    double rk = (double)indx * _dr;
    double t = (r - rk) / _dr;

    Eigen::VectorXd R = Eigen::VectorXd::Zero(4);
    R(0) = 1.0;
    R(1) = t;
    R(2) = t * t;
    R(3) = t * t * t;
    Eigen::VectorXd B = _lam.segment(indx, 4);
    u += ((R.transpose() * _M) * B).value();
    return u;

  } else {
    return 0.0;
  }
}

// calculate first derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateDF(const long int i,
                                           const double r) const {

  // since first _nexcl parameters are not optimized for stability reasons

  if (r <= _cut_off) {

    long int i_opt = i + _nexcl;
    long int indx;
    double rk;

    indx = min((long int)(r / _dr), _nbreak - 2);
    rk = (double)indx * _dr;

    if (i_opt >= indx && i_opt <= indx + 3) {

      Eigen::VectorXd R = Eigen::VectorXd::Zero(4);

      double t = (r - rk) / _dr;

      R(0) = 1.0;
      R(1) = t;
      R(2) = t * t;
      R(3) = t * t * t;

      Eigen::VectorXd RM = R.transpose() * _M;

      return RM(i_opt - indx);

    } else {
      return 0.0;
    }

  } else {
    return 0.0;
  }
}

// calculate second derivative w.r.t. ith parameter
double PotentialFunctionCBSPL::CalculateD2F(const long int i, const long int j,
                                            const double r) const {

  return 0.0;
}

}  // namespace csg
}  // namespace votca
