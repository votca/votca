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

// Standard includes
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>

// Local VOTCA includes
#include "votca/tools/histogram.h"

namespace votca {
namespace tools {

Histogram::Histogram() = default;

Histogram::Histogram(const options_t& op) : options_(op) {}

Histogram::~Histogram() = default;

void Histogram::ProcessData(DataCollection<double>::selection* data) {

  pdf_.assign(options_.n_, 0);

  if (options_.auto_interval_) {
    min_ = std::numeric_limits<double>::max();
    max_ = std::numeric_limits<double>::min();
    options_.extend_interval_ = true;
  } else {
    min_ = options_.min_;
    max_ = options_.max_;
  }
  for (auto& array : *data) {
    if (options_.extend_interval_ || options_.auto_interval_) {
      for (auto& value : *array) {
        min_ = std::min(value, min_);
        max_ = std::max(value, max_);
      }
    }
  }

  interval_ = (max_ - min_) / (double)(options_.n_ - 1);

  double v = 1.;
  for (auto& array : *data) {
    for (auto& value : *array) {
      Index ii = (Index)floor((value - min_) / interval_ +
                              0.5);  // the interval should
                                     // be centered around
                                     // the sampling point
      if (ii < 0 || ii >= options_.n_) {
        if (options_.periodic_) {
          while (ii < 0) {
            ii += options_.n_;
          }
          ii = ii % options_.n_;
        } else {
          continue;
        }
      }
      pdf_[ii] += v;
    }
  }

  if (options_.scale_ == "bond") {
    for (size_t i = 0; i < pdf_.size(); ++i) {
      double r = min_ + interval_ * (double)i;
      if (std::fabs(r) < 1e-10) {
        r = min_ + interval_ * (double)(i + 1);
        pdf_[i] = pdf_[i + 1];
      }
      pdf_[i] /= (r * r);
    }
  } else if (options_.scale_ == "angle") {
    for (size_t i = 0; i < pdf_.size(); ++i) {
      double alpha = min_ + interval_ * (double)i;
      double sa = sin(alpha);
      if (std::fabs(sa) < 1e-5) {
        if (i < pdf_.size() - 1) {
          alpha = min_ + interval_ * (double)(i + 1);
          pdf_[i] = pdf_[i + 1] / sin(alpha);
        } else {
          pdf_[i] = pdf_[i - 1];
        }

      } else {
        pdf_[i] /= sa;
      }
    }
  }

  if (options_.periodic_) {
    pdf_[0] = (pdf_[0] + pdf_[options_.n_ - 1]);
    pdf_[options_.n_ - 1] = pdf_[0];
  }
  if (options_.normalize_) {
    Normalize();
  }
}

void Histogram::Normalize() {
  double norm = 1. / (interval_ * accumulate(pdf_.begin(), pdf_.end(), 0.0));
  std::transform(pdf_.begin(), pdf_.end(), pdf_.begin(),
                 [norm](double a) { return a * norm; });
}

}  // namespace tools
}  // namespace votca
