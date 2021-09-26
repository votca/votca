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

#ifndef VOTCA_TOOLS_AVERAGE_H
#define VOTCA_TOOLS_AVERAGE_H

// Standard includes
#include <cmath>

namespace votca {
namespace tools {

template <typename T>
class Average {
 public:
  void Process(const T &value);
  void Clear();
  template <typename iterator_type>
  void ProcessRange(const iterator_type &begin, const iterator_type &end);

  T CalcDev() const;
  T CalcSig2() const;
  const T &getAvg() const;
  const T getM2() const;
  size_t getN() const;

 private:
  size_t n_ = 0;
  T av_ = 0;  // average
  T m2_ = 0;  // second moment
};

template <typename T>
inline void Average<T>::Process(const T &value) {
  av_ = av_ * (double)n_ / (double)(n_ + 1) + value / (double)(n_ + 1);
  n_++;
  m2_ += value * value;
}

template <typename T>
inline void Average<T>::Clear() {
  av_ = 0;
  n_ = 0;
  m2_ = 0;
}

template <typename T>
template <typename iterator_type>
void Average<T>::ProcessRange(const iterator_type &begin,
                              const iterator_type &end) {
  for (iterator_type iter = begin; iter != end; ++iter) {
    Process(*iter);
  }
}

template <typename T>
T Average<T>::CalcDev() const {
  return std::sqrt((m2_ - n_ * av_ * av_) / (n_ - 1));
}

template <typename T>
T Average<T>::CalcSig2() const {
  double dev = 0.0;
  dev = m2_ / n_ - av_ * av_;
  return dev;
}

template <typename T>
const T &Average<T>::getAvg() const {
  return av_;
}

template <typename T>
const T Average<T>::getM2() const {
  double m2 = 0.0;
  m2 = m2_ / n_;
  return m2;
}

template <typename T>
size_t Average<T>::getN() const {
  return n_;
}

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_AVERAGE_H
