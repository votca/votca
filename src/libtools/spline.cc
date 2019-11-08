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

#include <votca/tools/spline.h>

namespace votca {
namespace tools {

using namespace std;

long int Spline::GenerateGrid(double min, double max, double h) {
  int vec_size = (int)((max - min) / h + 1.00000001);
  _r.resize(vec_size);
  int i;

  double r_init;

  for (r_init = min, i = 0; i < vec_size - 1; r_init += h) {
    _r[i++] = r_init;
  }
  _r[i] = max;
  return _r.size();
}

Eigen::VectorXd Spline::Calculate(const Eigen::VectorXd &x) {
  Eigen::VectorXd y(x.size());
  for (long int i = 0; i < x.size(); ++i) {
    y(i) = Calculate(x(i));
  }
  return y;
}

Eigen::VectorXd Spline::CalculateDerivative(const Eigen::VectorXd &x) {
  Eigen::VectorXd y(x.size());
  for (long int i = 0; i < x.size(); ++i) {
    y(i) = CalculateDerivative(x(i));
  }
  return y;
}

void Spline::Print(std::ostream &out, double interval) {
  for (double x = _r[0]; x < _r[_r.size() - 1]; x += interval) {
    out << x << " " << Calculate(x) << "\n";
  }
}

long int Spline::getInterval(double r) {
  if (r < _r[0]) {
    return 0;
  }
  if (r > _r[_r.size() - 2]) {
    return _r.size() - 2;
  }
  long int i;
  for (i = 0; i < _r.size(); ++i) {
    if (_r[i] > r) {
      break;
    }
  }
  return i - 1;
}

double Spline::getGridPoint(int i) {
  if (i >= _r.size()) {
    return 0;
  }
  return _r[i];
}

}  // namespace tools
}  // namespace votca
