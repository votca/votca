/*
 *            Copyright 2009-2019 The VOTCA Development Team
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

#include <boost/filesystem.hpp>
#include <votca/xtp/ERIs.h>
#include <votca/xtp/classicalsegment.h>
#include <votca/xtp/convergenceacc.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/topology.h>

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
  DFTEngine(Orbitals& orbitals) : _orbitals(orbitals){};

  void Initialize(tools::Property& options);

  void setLogger(Logger* pLog) { _pLog = pLog; }

  void setExternalcharges(
      std::vector<std::unique_ptr<StaticSite> >* externalsites) {
    _externalsites = externalsites;
    _addexternalsites = true;
  }

  bool Evaluate();

  void Prepare();

  std::string getDFTBasisName() const { return _dftbasis_name; };

 private:
  Eigen::MatrixXd OrthogonalizeGuess(const Eigen::MatrixXd& GuessMOs) const;
  void PrintMOs(const Eigen::VectorXd& MOEnergies);
  void CalcElDipole() const;
  Mat_p_Energy CalculateERIs(const Eigen::MatrixXd& DMAT) const;
  Mat_p_Energy CalcEXXs(const Eigen::MatrixXd& MOs,
                        const Eigen::MatrixXd& DMAT) const;
  void ConfigOrbfile();
  void SetupInvariantMatrices();
  Eigen::MatrixXd AtomicGuess() const;
  std::string ReturnSmallGrid(const std::string& largegrid);

  Mat_p_Energy IntegrateExternalDensity(const Orbitals& extdensity) const;

  Eigen::MatrixXd RunAtomicDFT_unrestricted(const QMAtom& uniqueAtom) const;

  double NuclearRepulsion() const;
  double ExternalRepulsion() const;
  Eigen::MatrixXd SphericalAverageShells(const Eigen::MatrixXd& dmat,
                                         const AOBasis& dftbasis) const;
  Logger* _pLog;

  int _openmp_threads;

  // atoms
  Orbitals& _orbitals;

  // basis sets
  std::string _auxbasis_name;
  std::string _dftbasis_name;
  std::string _ecp_name;
  AOBasis _dftbasis;
  AOBasis _auxbasis;
  AOBasis _ecp;

  bool _with_ecp;
  bool _with_RI;

  std::string _four_center_method;  // direct | cache

  // Pre-screening
  bool _with_screening;
  double _screening_eps;

  // numerical integration Vxc
  std::string _grid_name;
  std::string _grid_name_small;
  bool _use_small_grid;
  NumericalIntegration _gridIntegration;
  NumericalIntegration _gridIntegration_small;

  // AO Matrices
  AOOverlap _dftAOoverlap;
  AOKinetic _dftAOkinetic;
  AOESP _dftAOESP;
  AOECP _dftAOECP;
  AODipole_Potential _dftAODipole_Potential;
  AOQuadrupole_Potential _dftAOQuadrupole_Potential;
  AOPlanewave _dftAOplanewave;

  bool _with_guess;
  std::string _initial_guess;

  // Convergence
  int _numofelectrons = 0;
  int _max_iter = 100;
  ConvergenceAcc::options _conv_opt;
  // DIIS variables
  ConvergenceAcc _conv_accelerator;
  // Electron repulsion integrals
  ERIs _ERIs;

  // external charges
  std::vector<std::unique_ptr<StaticSite> >* _externalsites;
  bool _addexternalsites = false;

  // exchange and correlation
  double _ScaHFX;
  std::string _xc_functional_name;

  bool _integrate_ext_density = false;
  // integrate external density
  std::string _orbfilename;
  std::string _gridquality;
  std::string _state;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_DFTENGINE_H
