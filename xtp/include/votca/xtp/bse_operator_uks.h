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

#include <Eigen/Core>
#include <votca/xtp/tcmatrix.h>
#include <votca/xtp/tools/eigensystem.h>

namespace votca {
namespace xtp {

class BSEOperatorUKS {
public:
  using Matrix = Eigen::MatrixXd;
  using Vector = Eigen::VectorXd;

  BSEOperatorUKS(const TCMatrix_gwbse& Mmn,
                 const Matrix& Hqp,
                 const Vector& omega,
                 const Matrix& V,
                 Index homo,
                 Index lumo,
                 Index vmin,
                 Index cmax);

  // matrix-vector product
  void Apply(const Vector& x, Vector& y) const;

  Index size() const { return size_; }

private:
  // spin-resolved inputs
  const TCMatrix_gwbse& Mmn_;
  const Matrix& Hqp_;

  // screening (shared)
  const Vector& omega_;
  const Matrix& V_;

  // orbital ranges
  Index homo_;
  Index lumo_;
  Index vmin_;
  Index cmax_;

  Index size_;

  // helpers
  void ApplyHqp(const Vector& x, Vector& y) const;
  void ApplyExchange(const Vector& x, Vector& y) const;
  void ApplyDirect(const Vector& x, Vector& y) const;
};

}  // namespace xtp
}  // namespace votca