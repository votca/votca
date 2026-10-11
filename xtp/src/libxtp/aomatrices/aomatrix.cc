/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
 * distributed under the License is distributed on an "A_ol I_ol" BA_olI_ol,
 * WITHOUT WARRANTIE_ol OR CONDITION_ol OF ANY KIND, either express or implied.
 *  olee_ the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Standard includes
#include <vector>
// Local VOTCA includes
#include "votca/xtp/aomatrix.h"
#include "votca/xtp/screening_kernels.h"

namespace votca {
namespace xtp {

Eigen::MatrixXd AOOverlap::Pseudo_InvSqrt(double etol) {
  const SymmetricEigenSystem es = SymmetricEigen(aomatrix_);
  smallestEigenvalue = es.values(0);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.values.size());
  removedfunctions = 0;
  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.values(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.values(i));
    }
  }

  return es.vectors * diagonal.asDiagonal() * es.vectors.transpose();
}

Eigen::MatrixXd AOOverlap::Sqrt() {
  const SymmetricEigenSystem es = SymmetricEigen(aomatrix_);
  smallestEigenvalue = es.values(0);
  return es.vectors * es.values.cwiseSqrt().asDiagonal() *
         es.vectors.transpose();
}

// This converts V into ((S-1/2 V S-1/2)-1/2 S-1/2)T, which is needed to
// construct 4c integrals,
Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt_GWBSE(const AOOverlap& auxoverlap,
                                                double etol) {

  const SymmetricEigenSystem eo = SymmetricEigen(auxoverlap.Matrix());
  removedfunctions = 0;
  Eigen::VectorXd diagonal_overlap = Eigen::VectorXd::Zero(eo.values.size());
  for (Index i = 0; i < diagonal_overlap.size(); ++i) {
    if (eo.values(i) < etol) {
      removedfunctions++;
    } else {
      diagonal_overlap(i) = 1.0 / std::sqrt(eo.values(i));
    }
  }
  Eigen::MatrixXd Ssqrt =
      eo.vectors * diagonal_overlap.asDiagonal() * eo.vectors.transpose();

  Eigen::MatrixXd ortho = Ssqrt * aomatrix_ * Ssqrt;
  const SymmetricEigenSystem es = SymmetricEigen(ortho);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.values.size());

  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.values(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.values(i));
    }
  }

  Eigen::MatrixXd Vm1 =
      es.vectors * diagonal.asDiagonal() * es.vectors.transpose();
  Eigen::MatrixXd result = (Vm1 * Ssqrt).transpose();
  return result;
}

Eigen::MatrixXd AOCoulomb::Pseudo_InvSqrt(double etol) {
  const SymmetricEigenSystem es = SymmetricEigen(aomatrix_);
  Eigen::VectorXd diagonal = Eigen::VectorXd::Zero(es.values.size());
  removedfunctions = 0;
  for (Index i = 0; i < diagonal.size(); ++i) {
    if (es.values(i) < etol) {
      removedfunctions++;
    } else {
      diagonal(i) = 1.0 / std::sqrt(es.values(i));
    }
  }

  return es.vectors * diagonal.asDiagonal() * es.vectors.transpose();
}

}  // namespace xtp
}  // namespace votca
