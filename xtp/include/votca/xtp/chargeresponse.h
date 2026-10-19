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
#ifndef VOTCA_XTP_CHARGERESPONSE_H
#define VOTCA_XTP_CHARGERESPONSE_H

// Standard includes
#include <string>
#include <vector>

// Local VOTCA includes
#include "checkpoint.h"
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Charge-response matrix Pi of one molecule type.
 *
 * In the charge-response (CR) environment model the induced atomic charges of
 * a molecule are dq = -Pi phi, with phi the external electrostatic potential
 * at its atoms. Pi (N x N, atomic units e^2/Hartree) is symmetric, positive
 * semidefinite and conserves charge (every row sums to zero). It already
 * contains the molecule's own Coulomb screening, so induced charges of the
 * same molecule do not interact. What charges cannot carry is in residual
 * induced dipoles, the ordinary polarizabilities of the molecule's .mps file.
 *
 * Pi is invariant under rotations, so it is stored once per molecule type in
 * the atom order of the type's .mps file, in a .crm file:
 *
 *   ! comment lines start with '!'
 *   N
 *   El x y z        (N lines, Angstrom, the fit geometry; informational)
 *   Pi_11 ... Pi_1N (N lines of N numbers, e^2/Hartree)
 *
 * written by xtp_tools -e crfit (ChargeResponseFit) or the PySCF script
 * xtp_cr_fit: the response to Hirshfeld-localized potentials, measured
 * with Hirshfeld charges.
 */
class ChargeResponse {
 public:
  ChargeResponse() = default;
  ChargeResponse(std::vector<std::string> elements, Eigen::MatrixXd pi);

  /// Reads a .crm file and checks it (see Check). Positions are in bohr.
  void LoadFromFile(const std::string& filename);

  bool empty() const { return pi_.size() == 0; }
  Index size() const { return Index(pi_.rows()); }
  const Eigen::MatrixXd& Pi() const { return pi_; }
  const std::vector<std::string>& Elements() const { return elements_; }

  /// The response of the atoms order[0], order[1], ... of this one, in that
  /// order. order must be a permutation of 0..N-1: a subset of atoms would
  /// not conserve charge.
  ChargeResponse Permuted(const std::vector<Index>& order) const;

  void WriteToCpt(CheckpointWriter& w) const;
  void ReadFromCpt(CheckpointReader& r);

 private:
  // Throws unless pi_ is square, matches elements_, symmetric, charge
  // conserving and positive semidefinite to a relative tolerance; then
  // removes the rounding of the file (symmetrizes and projects onto the
  // charge-conserving vectors).
  void Check(const std::string& source);

  std::vector<std::string> elements_;
  Eigen::MatrixXd pi_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CHARGERESPONSE_H
