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

// Local VOTCA includes
#include "votca/xtp/hirshfeldpartition.h"
#include "votca/xtp/aoshell.h"

namespace votca {
namespace xtp {

double HirshfeldPartition::EvaluateAtomicDensity(
    const AOBasis& atom_basis, const Eigen::MatrixXd& reference_density,
    const Eigen::Vector3d& point) {
  Eigen::VectorXd ao_values = Eigen::VectorXd::Zero(atom_basis.AOBasisSize());
  for (const AOShell& shell : atom_basis) {
    AOShell::AOValues vals = shell.EvalAOspace(point);
    ao_values.segment(shell.getStartIndex(), shell.getNumFunc()) =
        vals.values;
  }
  return ao_values.dot(reference_density * ao_values);
}

}  // namespace xtp
}  // namespace votca
