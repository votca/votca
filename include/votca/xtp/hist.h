/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <votca/xtp/checkpoint.h>
#include <votca/xtp/eigen.h>
#include <votca/xtp/energy_terms.h>

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
    if (_filled > 1) {
      return _metric - _metric_old;
    } else if (_filled == 1) {
      return _metric;
    } else {
      throw std::runtime_error("hist is not filled yet");
    }
  }

  const T& back() const { return _metric; }
  void push_back(const T& metric) {
    _metric_old = std::move(_metric);
    _metric = metric;
    _filled++;
  }
  void push_back(T&& metric) {
    _metric_old = std::move(_metric);
    _metric = std::move(metric);
    _filled++;
  }

  bool filled() const { return _filled > 0; }

  void WriteToCpt(CheckpointWriter& w) const {
    w(_filled, "filled");
    WriteMetric(_metric, "metric", w);
    WriteMetric(_metric_old, "metric_old", w);
  }

  void ReadFromCpt(CheckpointReader& r) {
    r(_filled, "filled");
    ReadMetric(_metric, "metric", r);
    ReadMetric(_metric_old, "metric_old", r);
  }

 private:
  void ReadMetric(T&, std::string tag, CheckpointReader& r);
  void WriteMetric(const T&, std::string tag, CheckpointWriter& w) const;
  int _filled = 0;
  T _metric;
  T _metric_old;
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

#endif  // VOTCA_XTP_JOB_H
