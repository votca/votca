/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef VOTCA_TOOLS_HISTOGRAM_H
#define VOTCA_TOOLS_HISTOGRAM_H

// Standard includes
#include <cmath>
#include <limits>
#include <vector>

// Local VOTCA includes
#include "datacollection.h"

namespace votca {
namespace tools {

/**
    \brief class to generate histograms

    This class produces a histogram out of a vector of values. This class
    is obsolate and only used in csg_boltzman. Use HistogramNew instead!
*/
class Histogram {
 public:
  struct options_t;

  /// constructor
  Histogram(const options_t &op);
  Histogram();
  /// destructor
  ~Histogram();

  /**
      process data and generate histogram
   */
  void ProcessData(DataCollection<double>::selection *data);

  /// returns the minimum value
  double getMin() const { return min_; }
  /// return the maximum value
  double getMax() const { return max_; }
  /// return the number of grid points
  Index getN() const { return options_.n_; }
  std::vector<double> &getPdf() { return pdf_; }
  double getInterval() const { return interval_; }

  void Normalize(void);

  struct options_t {
    Index n_ = 101;
    bool auto_interval_ = true;
    bool extend_interval_ = false;
    double min_ = 0;
    double max_ = 1.;
    bool periodic_ = false;
    bool normalize_ = true;
    std::string scale_ = "no";
  };

 private:
  std::vector<double> pdf_;
  double min_ = 0;
  double max_ = 0;
  double interval_;

  options_t options_;
};

inline std::ostream &operator<<(std::ostream &out, Histogram &h) {
  for (Index i = 0; i < h.getN(); i++) {
    out << h.getMin() + h.getInterval() * ((double)i + 0.0) << " "
        << h.getPdf()[i] << std::endl;
  }
  return out;
}

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_HISTOGRAM_H
