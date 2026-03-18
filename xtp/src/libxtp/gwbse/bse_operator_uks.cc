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

 #include <votca/xtp/bse_operator_uks.h>

namespace votca {
namespace xtp {

BSEOperatorUKS::BSEOperatorUKS(const TCMatrix_gwbse& Mmn,
                               const Matrix& Hqp,
                               const Vector& omega,
                               const Matrix& V,
                               Index homo,
                               Index lumo,
                               Index vmin,
                               Index cmax)
    : Mmn_(Mmn),
      Hqp_(Hqp),
      omega_(omega),
      V_(V),
      homo_(homo),
      lumo_(lumo),
      vmin_(vmin),
      cmax_(cmax) {

  size_ = (homo_ - vmin_ + 1) * (cmax_ - lumo_ + 1);
}

void BSEOperatorUKS::Apply(const Vector& x, Vector& y) const {
  y.setZero(size_);

  ApplyHqp(x, y);
  ApplyExchange(x, y);
  ApplyDirect(x, y);
}


