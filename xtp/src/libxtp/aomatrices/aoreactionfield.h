/*
 *            Copyright 2009-2024 The VOTCA Development Team
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
#ifndef VOTCA_XTP_AOREACTIONFIELD_H
#define VOTCA_XTP_AOREACTIONFIELD_H

// Local VOTCA includes
#include "aobasis.h"
#include "classicalsegment.h"
#include "eeinteractor.h"
#include "eigen.h"
#include "qmmolecule.h"
#include "vxc_grid.h"

namespace votca {
namespace xtp {

/**
 * \brief Computes the MM reaction-field correction to the 2-centre
 *        Coulomb matrix in the auxiliary basis,
 *
 *          (mu|nu)_reac = sum_alpha  Delta_d_mu(alpha) · grad_R <xi_mu|1/|r-R||xi_nu>|_{R_alpha}
 *
 *        where Delta_d_mu is the Thole-model induced dipole at MM site alpha
 *        driven by the electrostatic field of aux-basis function xi_mu.
 *
 *        The field of xi_mu at each MM site is evaluated numerically on a
 *        Vxc_Grid (step 1).  Once the induced dipoles are known the integral
 *        against xi_nu is evaluated analytically via AOMultipole (step 2).
 *        The resulting matrix is symmetric (enforced by explicit
 *        symmetrisation after assembly).
 *
 *        Usage:
 *          AOReactionField rf;
 *          rf.Fill(auxbasis, qmatoms, polar_segments, "xfine");
 *          // then add rf.Matrix() to the AOCoulomb matrix before
 *          // calling Pseudo_InvSqrt_GWBSE.
 */
class AOReactionField {
 public:
  AOReactionField() = default;

  /**
   * Compute the (Naux x Naux) reaction-field correction matrix.
   *
   * @param auxbasis      The auxiliary basis set (RI basis).
   * @param qmatoms       The QM molecule (needed to centre the numerical
   *                      grid for step 1).
   * @param polar_segs    The MM polar segments containing Thole
   *                      polarisabilities.  Their V_ and V_noE_ fields are
   *                      used as scratch and are reset on entry and exit.
   * @param grid_quality  Vxc_Grid quality string, e.g. "xfine".
   *                      "xfine" is recommended because aux-basis functions
   *                      can be more contracted than the DFT basis.
   */
  void Fill(const AOBasis& auxbasis, const QMMolecule& qmatoms,
            std::vector<PolarSegment>& polar_segs,
            const std::string& grid_quality = "xfine");

  Index Dimension() const { return matrix_.rows(); }

  const Eigen::MatrixXd& Matrix() const { return matrix_; }

 private:
  Eigen::MatrixXd matrix_;

  // Evaluate the electrostatic field of a single aux-basis function xi_mu
  // at all MM polar-site positions, using the supplied pre-built grid.
  // The grid must have been built from auxbasis / qmatoms before this call.
  // Returns a (3*N_sites) vector: [Ex_0, Ey_0, Ez_0, Ex_1, ...].
  Eigen::VectorXd CalcAuxFunctionField(
      const AOShell& shell, Index func_in_shell, const Vxc_Grid& grid,
      const std::vector<PolarSegment>& polar_segs) const;

  // Given induced dipoles stored in polar_segs (Induced_Dipole()),
  // compute the full (Naux x Naux) potential matrix analytically via
  // AOMultipole and return column mu_idx of that matrix.
  Eigen::VectorXd CalcDipolePotentialColumn(
      Index mu_idx, const AOBasis& auxbasis,
      const std::vector<PolarSegment>& polar_segs) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOREACTIONFIELD_H
