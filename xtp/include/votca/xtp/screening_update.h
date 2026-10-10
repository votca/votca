/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_SCREENING_UPDATE_H
#define VOTCA_XTP_SCREENING_UPDATE_H

// Standard includes
#include <stdexcept>
#include <string>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief When evGW rebuilds the screening W.
 *
 * evGW updates the energies in G and in W. Sigma_c reads the energies in G
 * when it is evaluated, but W only changes when the screening is rebuilt
 * (RPA, PPM parameters or RPA pair eigenproblem), which is most of the cost
 * of an iteration.
 *
 *  - every:    rebuild W in every iteration (the classic scheme);
 *  - adaptive: iterate G with W fixed and rebuild W once the last change of
 *              the energies is small against how far they have moved since
 *              W was built (below ratio times that drift), or after
 *              max_inner iterations. Converged when the energies change by
 *              less than the limit AND differ from those W was built with by
 *              less than the limit, i.e. at the same evGW fixed point.
 */
class ScreeningUpdate {
 public:
  enum class Mode { Every, Adaptive };
  enum class Next { Continue, Rebuild, Converged };

  static Mode Parse(const std::string& name) {
    if (name == "every") {
      return Mode::Every;
    } else if (name == "adaptive") {
      return Mode::Adaptive;
    }
    throw std::runtime_error("Unknown gw.screening_update '" + name +
                             "'; choices are every, adaptive");
  }
  static std::string Name(Mode mode) {
    return (mode == Mode::Adaptive) ? "adaptive" : "every";
  }

  void Configure(Mode mode, double ratio, Index max_inner) {
    mode_ = mode;
    ratio_ = ratio;
    max_inner_ = max_inner;
  }
  Mode mode() const { return mode_; }

  /// W has just been built from these energies
  void Built(const Eigen::VectorXd& energies) {
    built_from_ = energies;
    inner_ = 0;
    ++builds_;
  }
  Index builds() const { return builds_; }
  /// QP iterations since W was last built (before the current one)
  Index inner() const { return inner_; }

  /// After a QP iteration: old and new energies (the input of this and of
  /// the next iteration). Returns what to do next and the drift since W was
  /// built (largest change of an energy).
  Next Decide(const Eigen::VectorXd& old_energies,
              const Eigen::VectorXd& new_energies, double limit,
              double* drift_out = nullptr) {
    ++inner_;
    const double step = (new_energies - old_energies).cwiseAbs().maxCoeff();
    const double drift = (new_energies - built_from_).cwiseAbs().maxCoeff();
    if (drift_out != nullptr) {
      *drift_out = drift;
    }
    if (mode_ == Mode::Every) {
      return (step < limit) ? Next::Converged : Next::Rebuild;
    }
    if (step < limit) {
      return (drift < limit) ? Next::Converged : Next::Rebuild;
    }
    if (step < ratio_ * drift || inner_ >= max_inner_) {
      return Next::Rebuild;
    }
    return Next::Continue;
  }

 private:
  Mode mode_ = Mode::Every;
  double ratio_ = 0.25;
  Index max_inner_ = 10;
  Eigen::VectorXd built_from_;
  Index inner_ = 0;
  Index builds_ = 0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_SCREENING_UPDATE_H
