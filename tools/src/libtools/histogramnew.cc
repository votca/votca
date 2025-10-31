/*
 * Copyright 2009-2025 The VOTCA Development Team (http://www.votca.org)
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
#include <algorithm>

// Local VOTCA includes
#include "votca/tools/histogram.h"

using namespace std;
namespace votca {
namespace tools {

void Histogram::Initialize_() {
  if (periodic_) {
    step_ = (max_ - min_) / double(nbins_);
  } else {
    step_ = (max_ - min_) / (double(nbins_) - 1.0);
  }

  if (nbins_ == 1) {
    step_ = 1;
  }
  data_.resize(nbins_);
  double v = min_;
  for (Index i = 0; i < nbins_; v += step_, ++i) {
    data_.x(i) = v;
  }
  data_.y() = Eigen::VectorXd::Zero(nbins_);
  data_.yerr() = Eigen::VectorXd::Zero(nbins_);
  data_.flags() = std::vector<char>(nbins_, 'i');
}

void Histogram::Initialize(double min, double max, Index nbins) {
  min_ = min;
  max_ = max;
  nbins_ = nbins;
  Initialize_();
}

void Histogram::Process(const double &v, double scale) {

  Index i = (Index)floor((v - min_) / step_ + 0.5);
  if (i < 0 || i >= nbins_) {
    if (periodic_) {
      if (i < 0) {
        i = nbins_ - ((-i) % nbins_);
      } else {
        i = i % nbins_;
      }
    } else {
      return;
    }
  }
  data_.y(i) += scale;
}

double Histogram::getMinBinVal() const { return data_.getMinY(); }

double Histogram::getMaxBinVal() const { return data_.getMaxY(); }

pair<double, double> Histogram::getInterval(Index bin) const {
  pair<double, double> bounds;
  double value = static_cast<double>(bin);
  bounds.first = value * step_ + min_;
  bounds.second += step_;
  return bounds;
}

void Histogram::Normalize() {
  double area = 0;
  area = data_.y().cwiseAbs().sum() * step_;
  double scale = 1. / area;
  data_.y() *= scale;
}

void Histogram::Clear() {
  data_.y() = Eigen::VectorXd::Zero(nbins_);
  data_.yerr() = Eigen::VectorXd::Zero(nbins_);
}

}  // namespace tools
}  // namespace votca
