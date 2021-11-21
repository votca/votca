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
#ifndef VOTCA_XTP_PPM_H
#define VOTCA_XTP_PPM_H

// Local VOTCA includes
#include "eigen.h"
#include "rpa.h"

namespace votca {
namespace xtp {

class PPM {
 public:
  PPM() : screening_r(0.0), screening_i(0.5){};

  // This sets the screening frequencies for real and imaginary part in hartree

  void PPM_construct_parameters(const RPA& rpa);

  const Eigen::VectorXd& getPpm_weight() const { return ppm_weight_; }

  const Eigen::VectorXd& getPpm_freq() const { return ppm_freq_; }

  const Eigen::MatrixXd& getPpm_phi() const { return ppm_phi_; }

  void FreeMatrix() { ppm_phi_.resize(0, 0); }

 private:
  double screening_r;
  double screening_i;

  // PPM related variables and functions
  Eigen::MatrixXd ppm_phi_;
  Eigen::VectorXd ppm_freq_;
  Eigen::VectorXd ppm_weight_;
};
}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_PPM_H
