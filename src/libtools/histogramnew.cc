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

#include <algorithm>
#include <votca/tools/histogramnew.h>

using namespace std;
namespace votca {
namespace tools {

void HistogramNew::Initialize_() {
  if (_periodic) {
    _step = (_max - _min) / double(_nbins);
  } else {
    _step = (_max - _min) / (double(_nbins) - 1.0);
  }

  if (_nbins == 1) {
    _step = 1;
  }
  _data.resize(_nbins);
  double v = _min;
  for (Index i = 0; i < _nbins; v += _step, ++i) {
    _data.x(i) = v;
  }
  _data.y() = Eigen::VectorXd::Zero(_nbins);
  _data.yerr() = Eigen::VectorXd::Zero(_nbins);
  _data.flags() = std::vector<char>(_nbins, 'i');
}

void HistogramNew::Initialize(double min, double max, Index nbins) {
  _min = min;
  _max = max;
  _nbins = nbins;
  Initialize_();
}

void HistogramNew::Process(const double &v, double scale) {

  Index i = (Index)floor((v - _min) / _step + 0.5);
  if (i < 0 || i >= _nbins) {
    if (_periodic) {
      if (i < 0) {
        i = _nbins - ((-i) % _nbins);
      } else {
        i = i % _nbins;
      }
    } else {
      return;
    }
  }
  _data.y(i) += scale;
}

double HistogramNew::getMinBinVal() const { return _data.getMinY(); }

double HistogramNew::getMaxBinVal() const { return _data.getMaxY(); }

pair<double, double> HistogramNew::getInterval(Index bin) const {
  pair<double, double> bounds;
  double value = static_cast<double>(bin);
  bounds.first = value * _step + _min;
  bounds.second += _step;
  return bounds;
}

void HistogramNew::Normalize() {
  double area = 0;
  area = _data.y().cwiseAbs().sum() * _step;
  double scale = 1. / area;
  _data.y() *= scale;
}

void HistogramNew::Clear() {
  _data.y() = Eigen::VectorXd::Zero(_nbins);
  _data.yerr() = Eigen::VectorXd::Zero(_nbins);
}

}  // namespace tools
}  // namespace votca
