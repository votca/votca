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
#ifndef VOTCA_XTP_DFTENGINE_H
#define VOTCA_XTP_DFTENGINE_H

// VOTCA includes
#include <votca/tools/property.h>

// Local VOTCA includes
#include "ERIs.h"
#include "convergenceacc.h"
#include "ecpaobasis.h"
#include "logger.h"
#include "staticsite.h"
#include "vxc_grid.h"
#include "vxc_potential.h"

namespace votca {
namespace xtp {
class Orbitals;

/**
 * \brief Electronic ground-state via Density-Functional Theory
 *
 * Evaluates electronic ground state in molecular systems based on
 * density functional theory with Gaussian Orbitals.
 *
 */

class DFTEngine {
 public:
  void Initialize(tools::Property& options);

  void setLogger(Logger* pLog) { pLog_ = pLog; }

  void setExternalcharges(
      std::vector<std::unique_ptr<StaticSite> >* externalsites) {
    externalsites_ = externalsites;
    addexternalsites_ = true;
  }

  bool Evaluate(Orbitals& orb);

  std::string getDFTBasisName() const { return dftbasis_name_; };

 private:
  void Prepare(QMMolecule& mol);

  Vxc_Potential<Vxc_Grid> SetupVxc(const QMMolecule& mol);

  Eigen::MatrixXd OrthogonalizeGuess(const Eigen::MatrixXd& GuessMOs) const;
  void PrintMOs(const Eigen::VectorXd& MOEnergies, Log::Level level);
  void CalcElDipole(const Orbitals& orb) const;

  std::array<Eigen::MatrixXd, 2> CalcERIs_EXX(const Eigen::MatrixXd& MOCoeff,
                                              const Eigen::MatrixXd& Dmat,
                                              double error) const;

  Eigen::MatrixXd CalcERIs(const Eigen::MatrixXd& Dmat, double error) const;

  void ConfigOrbfile(Orbitals& orb);
  void SetupInvariantMatrices();

  Mat_p_Energy SetupH0(const QMMolecule& mol) const;
  Mat_p_Energy IntegrateExternalMultipoles(
      const QMMolecule& mol,
      const std::vector<std::unique_ptr<StaticSite> >& multipoles) const;
  Mat_p_Energy IntegrateExternalDensity(const QMMolecule& mol,
                                        const Orbitals& extdensity) const;

  Eigen::MatrixXd IntegrateExternalField(const QMMolecule& mol) const;

  tools::EigenSystem IndependentElectronGuess(const Mat_p_Energy& H0) const;
  tools::EigenSystem ModelPotentialGuess(
      const Mat_p_Energy& H0, const QMMolecule& mol,
      const Vxc_Potential<Vxc_Grid>& vxcpotential) const;

  Eigen::MatrixXd AtomicGuess(const QMMolecule& mol) const;

  Eigen::MatrixXd RunAtomicDFT_unrestricted(const QMAtom& uniqueAtom) const;

  double NuclearRepulsion(const QMMolecule& mol) const;
  double ExternalRepulsion(
      const QMMolecule& mol,
      const std::vector<std::unique_ptr<StaticSite> >& multipoles) const;
  Eigen::MatrixXd SphericalAverageShells(const Eigen::MatrixXd& dmat,
                                         const AOBasis& dftbasis) const;
  Logger* pLog_;

  // basis sets
  std::string auxbasis_name_;
  std::string dftbasis_name_;
  std::string ecp_name_;
  AOBasis dftbasis_;
  AOBasis auxbasis_;
  ECPAOBasis ecp_;

  bool with_ecp_;

  Index fock_matrix_reset_;
  // Pre-screening
  double screening_eps_;

  // numerical integration Vxc
  std::string grid_name_;

  // AO Matrices
  AOOverlap dftAOoverlap_;

  bool with_guess_;
  std::string initial_guess_;

  // Convergence
  Index numofelectrons_ = 0;
  Index max_iter_;
  ConvergenceAcc::options conv_opt_;
  // DIIS variables
  ConvergenceAcc conv_accelerator_;
  // Electron repulsion integrals
  ERIs ERIs_;

  // external charges
  std::vector<std::unique_ptr<StaticSite> >* externalsites_;
  bool addexternalsites_ = false;

  // exchange and correlation
  double ScaHFX_;
  std::string xc_functional_name_;

  bool integrate_ext_density_ = false;
  // integrate external density
  std::string orbfilename_;
  std::string gridquality_;
  std::string state_;

  Eigen::Vector3d extfield_ = Eigen::Vector3d::Zero();
  bool integrate_ext_field_ = false;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTENGINE_H
