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

#include <votca/tools/histogramnew.h>
#include <algorithm>

using namespace std;
namespace votca {
  namespace tools {

    HistogramNew::HistogramNew() {
      _min = _max = _step = _step_p = 0;
      _weight = 1.;
      _periodic = false;
    }

    HistogramNew::HistogramNew(const HistogramNew &hist)
    : _min(hist._min),
    _max(hist._max),
    _step(hist._step),
    _step_p(hist._step_p),
    _weight(hist._weight),
    _periodic(hist._periodic) {
    }

    void HistogramNew::Initialize_(double min, double max) {
      _step = (_max - _min) / (_nbins - 1.0);
      _data.resize(_nbins);
      for (double v = _min, i = 0; i < _nbins; v += _step, ++i) {
        _data.x(i) = v;
      }
      _data.y() = Eigen::VectorXd::Zero(_nbins);
      _data.yerr() = Eigen::VectorXd::Zero(_nbins);
      _data.flags() = std::vector<char>(_nbins, 'i');
    }

    void HistogramNew::InitializeP_(double min, double max) {
      _step_p = (_max - _min) / (_nbins);
      _data_p.resize(_nbins);
      for (double v_p = _min, i = 0; i < _nbins; v_p += _step_p, ++i) {
        _data_p.x(i) = v_p;
      }
      _data_p.y() = Eigen::VectorXd::Zero(_nbins);
      _data_p.yerr() = Eigen::VectorXd::Zero(_nbins);
      _data_p.flags() = std::vector<char>(_nbins, 'i');
    }

    void HistogramNew::Initialize(double min, double max, int nbins) {
      _min = min;
      _max = max;
      _weight = 1.;
      _nbins = nbins;
      Initialize_(min, max);
      InitializeP_(min, max);
    }

    void HistogramNew::Process_(const double &v, double scale) {
      int i = (int) floor((v - _min) / _step + 0.5);
      if (i < 0 || i >= _nbins) {
        return;
      }
      _data.y(i) += _weight * scale;
    }

    void HistogramNew::ProcessP_(const double &v, double scale) {
      int i = (int) floor((v - _min) / _step_p + 0.5);
      if (i < 0 || i >= _nbins) {
        if (i < 0)
          i = _nbins - ((-i) % _nbins);
        else
          i = i % _nbins;
      }
      _data_p.y(i) += _weight * scale;
    }

    void HistogramNew::Process(const double &v, double scale) {
      Process_(v, scale);
      ProcessP_(v, scale);
    }

    double HistogramNew::getMinBinVal(void) const {
      return _periodic ? _data_p.getMinY() : _data.getMinY();
    }

    double HistogramNew::getMaxBinVal(void) const {
      return _periodic ? _data_p.getMaxY() : _data.getMaxY();
    }

    pair<double, double> HistogramNew::getInterval(int bin) const {
      pair<double, double> bounds;
      double value = static_cast<double> (bin);
      if (_periodic) {
        bounds.first = value * _step_p + _min;
        bounds.second += _step_p;
      } else {
        bounds.first = value * _step + _min;
        bounds.second += _step;
      }
      return bounds;
    }

    void HistogramNew::Normalize_(const double step, Table &data) {
      double area = 0;
      area = data.y().cwiseAbs().sum() * step;
      double scale = 1. / area;
      data.y() *= scale;
    }

    void HistogramNew::Normalize() {
      Normalize_(_step, _data);
      Normalize_(_step_p, _data_p);
    }

    void HistogramNew::Clear() {
      _weight = 1.;
      _data.y() = Eigen::VectorXd::Zero(_nbins);
      _data.yerr() = Eigen::VectorXd::Zero(_nbins);
      _data_p.y() = Eigen::VectorXd::Zero(_nbins);
      _data_p.yerr() = Eigen::VectorXd::Zero(_nbins);
    }
  }
}
