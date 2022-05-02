/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_HIST_H
#define VOTCA_XTP_HIST_H

// Local VOTCA includes
#include "checkpoint.h"
#include "eigen.h"
#include "energy_terms.h"

/**
 * \brief Small container for convergence checking, stores current and last
 * value and gives the diff, or just the current value in the first iteration
 *
 *
 */

namespace votca {
namespace xtp {
template <class T>
class hist {
 public:
  T getDiff() const {
    if (filled_ > 1) {
      return metric_ - metric_old_;
    } else if (filled_ == 1) {
      return metric_;
    } else {
      throw std::runtime_error("hist is not filled yet");
    }
  }

  const T& back() const { return metric_; }
  void push_back(const T& metric) {
    metric_old_ = std::move(metric_);
    metric_ = metric;
    filled_++;
  }
  void push_back(T&& metric) {
    metric_old_ = std::move(metric_);
    metric_ = std::move(metric);
    filled_++;
  }

  bool filled() const { return filled_ > 0; }

  void WriteToCpt(CheckpointWriter& w) const {
    w(filled_, "filled");
    if (filled_ > 0) {
      WriteMetric(metric_, "metric", w);
    }
    if (filled_ > 1) {
      WriteMetric(metric_old_, "metric_old", w);
    }
  }

  void ReadFromCpt(CheckpointReader& r) {
    r(filled_, "filled");
    if (filled_ > 0) {
      ReadMetric(metric_, "metric", r);
    }
    if (filled_ > 1) {
      ReadMetric(metric_old_, "metric_old", r);
    }
  }

 private:
  void ReadMetric(T&, std::string tag, CheckpointReader& r);
  void WriteMetric(const T&, std::string tag, CheckpointWriter& w) const;
  Index filled_ = 0;
  T metric_;
  T metric_old_;
};

template <class T>
inline void hist<T>::ReadMetric(T& metric, std::string tag,
                                CheckpointReader& r) {
  r(metric, tag);
}
template <class T>
inline void hist<T>::WriteMetric(const T& metric, std::string tag,
                                 CheckpointWriter& w) const {
  w(metric, tag);
}

template <>
inline void hist<Energy_terms>::ReadMetric(Energy_terms& metric,
                                           std::string tag,
                                           CheckpointReader& r) {
  r(metric.data(), tag);
}
template <>
inline void hist<Energy_terms>::WriteMetric(const Energy_terms& metric,
                                            std::string tag,
                                            CheckpointWriter& w) const {
  w(metric.data(), tag);
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_HIST_H
