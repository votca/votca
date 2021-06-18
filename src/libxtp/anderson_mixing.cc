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

// Third party includes
#include <boost/format.hpp>

// Local VOTCA includes
#include "votca/xtp/anderson_mixing.h"

namespace votca {
namespace xtp {

void Anderson::Configure(const Index order, const double alpha) {
  order_ = order + 1;
  alpha_ = alpha;
}

void Anderson::UpdateOutput(const Eigen::VectorXd &newOutput) {

  // Check if max mixing history is reached and adding new step to history
  Index size = output_.size();
  if (size > order_ - 1) {
    output_.erase(output_.begin());
  }
  output_.push_back(newOutput);
}

void Anderson::UpdateInput(const Eigen::VectorXd &newInput) {
  Index size = output_.size();
  if (size > order_ - 1) {
    input_.erase(input_.begin());
  }
  input_.push_back(newInput);
}

const Eigen::VectorXd Anderson::MixHistory() {

  const Index iteration = output_.size();
  const Index used_history = iteration - 1;
  Eigen::VectorXd OutMixed = output_.back();
  Eigen::VectorXd InMixed = input_.back();

  if (iteration > 1 && order_ > 1) {

    Eigen::VectorXd DeltaN = OutMixed - InMixed;

    // Building Linear System for Coefficients
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(used_history, used_history);
    Eigen::VectorXd c = Eigen::VectorXd::Zero(used_history);

    for (Index m = 1; m < iteration; m++) {

      c(m - 1) = (DeltaN - output_[used_history - m] + input_[used_history - m])
                     .dot(DeltaN);

      for (Index j = 1; j < iteration; j++) {
        A(m - 1, j - 1) =
            (DeltaN - output_[used_history - m] + input_[used_history - m])
                .dot((DeltaN - output_[used_history - j] +
                      input_[used_history - j]));
      }
    }
    // Solving the System to obtain coefficients
    Eigen::VectorXd coefficients = A.fullPivHouseholderQr().solve(c);

    // Mixing the Potentials
    for (Index n = 1; n < iteration; n++) {

      OutMixed += coefficients(n - 1) *
                  (output_[used_history - n] - output_[used_history]);
      InMixed += coefficients(n - 1) *
                 (input_[used_history - n] - input_[used_history]);
    }
  }

  // Returning the linear Mix of Input and Output
  return alpha_ * OutMixed + (1 - alpha_) * InMixed;
}
}  // namespace xtp
}  // namespace votca
