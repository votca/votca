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

#pragma once
#ifndef _histogram_H
#define _histogram_H

#include "datacollection.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace votca {
namespace tools {

using namespace std;

/**
    \brief class to generate histograms

    This class produces a histogram out of a vector of values. This class
    is obsolate and only used in csg_boltzman. Use HistogramNew instead!
*/
class Histogram {
 public:
  struct options_t;

  /// constructor
  Histogram(options_t op);
  Histogram();
  /// destructor
  ~Histogram();

  /**
      process data and generate histogram
   */
  void ProcessData(DataCollection<double>::selection *data);

  /// returns the minimum value
  double getMin() const { return _min; }
  /// return the maximum value
  double getMax() const { return _max; }
  /// return the number of grid points
  double getN() const { return _options._n; }
  vector<double> &getPdf() { return _pdf; }
  double getInterval() const { return _interval; }

  void Normalize(void);

  struct options_t {
    int _n;
    bool _auto_interval;
    bool _extend_interval;
    double _min, _max;
    bool _periodic;
    bool _normalize;
    string _scale;

    options_t() {
      _n = 101;
      _auto_interval = true;
      _extend_interval = false;
      _min = 0.;
      _max = 1.;
      _periodic = false;
      _normalize = true;
      _scale = "no";
    }
  };

 private:
  vector<double> _pdf;
  double _min, _max;
  double _interval;

  options_t _options;
};

inline ostream &operator<<(ostream &out, Histogram &h) {
  for (int i = 0; i < h.getN(); i++) {
    out << h.getMin() + h.getInterval() * ((double)i + 0.0) << " "
        << h.getPdf()[i] << endl;
  }
  return out;
}

}  // namespace tools
}  // namespace votca

#endif /* _histogram_H */
