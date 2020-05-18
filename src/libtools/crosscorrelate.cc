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

#include "../../include/votca/tools/crosscorrelate.h"

namespace votca {
namespace tools {

void CrossCorrelate::AutoCorrelate(DataCollection<double>::selection& data) {
  Index N = data[0].size();
  Eigen::Map<Eigen::VectorXd> input(data[0].data(), N);
  Eigen::FFT<double> fft;
  Eigen::VectorXcd frequency = fft.fwd(input);
  Eigen::VectorXcd magnitude = frequency.cwiseAbs2();

  _corrfunc.resize(N);
  Eigen::Map<Eigen::VectorXd> corr_map(_corrfunc.data(), N);
  corr_map = fft.inv(magnitude);
  double d = corr_map(0);
  corr_map.array() /= d;
}

}  // namespace tools
}  // namespace votca
