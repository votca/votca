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

#include <limits>
#include <math.h>
#include <numeric>
#include <votca/tools/histogram.h>

namespace votca {
namespace tools {

Histogram::Histogram() : _min(0), _max(0) {}

Histogram::Histogram(options_t op) : _min(0), _max(0), _options(op) {}

Histogram::~Histogram() {}

void Histogram::ProcessData(DataCollection<double>::selection *data) {
  DataCollection<double>::selection::iterator array;
  DataCollection<double>::array::iterator iter;
  int ii;
  int ndata = 0;

  _pdf.assign(_options._n, 0);

  if (_options._auto_interval) {
    _min = std::numeric_limits<double>::max();
    _max = std::numeric_limits<double>::min();
    _options._extend_interval = true;
  } else {
    _min = _options._min;
    _max = _options._max;
  }

  for (array = data->begin(); array != data->end(); ++array) {
    ndata += (*array)->size();
    if (_options._extend_interval || _options._auto_interval) {
      for (iter = (*array)->begin(); iter != (*array)->end(); ++iter) {
        _min = std::min(*iter, _min);
        _max = std::max(*iter, _max);
      }
    }
  }

  // make that the highes value fits into interval
  // if(_options._auto_interval || _max!=_options._max)
  //    _max = _max + 0.5*(_max - _min)/(double)(_options._n);

  _interval = (_max - _min) / (double)(_options._n - 1);

  double v = 1.;
  for (array = data->begin(); array != data->end(); ++array) {
    for (iter = (*array)->begin(); iter != (*array)->end(); ++iter) {
      ii = (int)floor((*iter - _min) / _interval + 0.5);  // the interval should
                                                          // be centered around
                                                          // the sampling point
      if (ii < 0 || ii >= _options._n) {
        if (_options._periodic) {
          while (ii < 0) ii += _options._n;
          ii = ii % _options._n;
        } else {
          continue;
        }  // cout << "[histogram.cc]: out of bounds" << endl; continue;}
      }
      _pdf[ii] += v;
    }
  }

  // cout << _pdf.size() << " " << _options._periodic << endl;
  if (_options._scale == "bond") {
    for (size_t i = 0; i < _pdf.size(); ++i) {
      double r = _min + _interval * (double)i;
      if (abs(r) < 1e-10) {
        r = _min + _interval * (double)(i + 1);
        _pdf[i] = _pdf[i + 1];
      }
      _pdf[i] /= (r * r);
    }
  } else if (_options._scale == "angle") {
    for (size_t i = 0; i < _pdf.size(); ++i) {
      double alpha = _min + _interval * (double)i;
      double sa = sin(alpha);
      if (abs(sa) < 1e-5) {
        if (i < _pdf.size() - 1) {
          alpha = _min + _interval * (double)(i + 1);
          _pdf[i] = _pdf[i + 1] / sin(alpha);
        } else {
          _pdf[i] = _pdf[i - 1];
        }

      } else
        _pdf[i] /= sa;
    }
  }

  if (_options._periodic) {
    _pdf[0] = (_pdf[0] + _pdf[_options._n - 1]);
    _pdf[_options._n - 1] = _pdf[0];
  }
  if (_options._normalize) Normalize();
}

void Histogram::Normalize(void) {
  double norm = 1. / (_interval * accumulate(_pdf.begin(), _pdf.end(), 0.0));
  std::transform(_pdf.begin(), _pdf.end(), _pdf.begin(),
                 std::bind2nd(std::multiplies<double>(), norm));
}

}  // namespace tools
}  // namespace votca
