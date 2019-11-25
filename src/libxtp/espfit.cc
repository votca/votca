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

#include <votca/tools/constants.h>
#include <votca/xtp/aomatrix.h>
#include <votca/xtp/espfit.h>
#include <votca/xtp/grid.h>
#include <votca/xtp/numerical_integrations.h>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
using std::flush;

StaticSegment Espfit::Fit2Density(const Orbitals& orbitals,
                                  const QMState& state, std::string gridsize) {

  const Eigen::MatrixXd dmat = orbitals.DensityMatrixFull(state);
  // setting up grid
  Grid grid;
  grid.setupCHELPGGrid(orbitals.QMAtoms());
  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Done setting up CHELPG grid with " << grid.size()
      << " points " << flush;

  // Calculating nuclear potential at gridpoints
  AOBasis basis = orbitals.SetupDftBasis();
  AOOverlap overlap;
  overlap.Fill(basis);
  double N_comp = dmat.cwiseProduct(overlap.Matrix()).sum();

  NumericalIntegration numway;

  numway.GridSetup(gridsize, orbitals.QMAtoms(), basis);
  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Setup " << gridsize << " Numerical Grid with "
      << numway.getGridSize() << " gridpoints." << flush;
  double N = numway.IntegrateDensity(dmat);
  XTP_LOG_SAVE(Log::error, _log)
      << TimeStamp()
      << " Calculated Densities at Numerical Grid, Number of electrons is " << N
      << flush;

  if (std::abs(N - N_comp) > 0.001) {
    XTP_LOG_SAVE(Log::error, _log) << "=======================" << flush;
    XTP_LOG_SAVE(Log::error, _log)
        << "WARNING: Calculated Densities at Numerical Grid, Number of "
           "electrons "
        << N << " is far away from the the real value " << N_comp
        << ", you should increase the accuracy of the integration grid."
        << flush;
    N = N_comp;
    XTP_LOG_SAVE(Log::error, _log)
        << "WARNING: Electronnumber set to " << N << flush;
    XTP_LOG_SAVE(Log::error, _log) << "=======================" << flush;
  }

  XTP_LOG_SAVE(Log::error, _log)
      << TimeStamp() << " Calculating ESP at CHELPG grid points" << flush;
#pragma omp parallel for
  for (Index i = 0; i < grid.size(); i++) {
    grid.getGridValues()(i) =
        numway.IntegratePotential(grid.getGridPositions()[i]);
  }

  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Electron contribution calculated" << flush;
  double netcharge = 0.0;
  if (!state.isTransition()) {
    EvalNuclearPotential(orbitals.QMAtoms(), grid);
    Index Znuc = 0;
    for (const QMAtom& atom : orbitals.QMAtoms()) {
      Znuc += atom.getNuccharge();
    }
    netcharge = double(Znuc) - N;
  }
  netcharge = std::round(netcharge);
  XTP_LOG_SAVE(Log::error, _log)
      << TimeStamp() << " Netcharge constrained to " << netcharge << flush;
  return FitPartialCharges(orbitals, grid, netcharge);
  ;
}

void Espfit::EvalNuclearPotential(const QMMolecule& atoms, Grid& grid) {

  const std::vector<Eigen::Vector3d>& gridpoints = grid.getGridPositions();
  Eigen::VectorXd& gridvalues = grid.getGridValues();
  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Calculating ESP of nuclei at CHELPG grid points"
      << flush;

  for (Index i = 0; i < Index(gridpoints.size()); i++) {
    for (Index j = 0; j < atoms.size(); j++) {
      const Eigen::Vector3d& posatom = atoms[j].getPos();
      Index Znuc = atoms[j].getNuccharge();
      double dist_j = (gridpoints[i] - posatom).norm();
      gridvalues(i) += double(Znuc) / dist_j;
    }
  }
  return;
}

StaticSegment Espfit::FitPartialCharges(const Orbitals& orbitals,
                                        const Grid& grid, double netcharge) {
  const QMMolecule& atomlist = orbitals.QMAtoms();
  const Index NoOfConstraints =
      1 + Index(_regionconstraint.size()) + Index(_pairconstraint.size());
  const Index matrixSize = atomlist.size() + NoOfConstraints;
  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Setting up Matrices for fitting of size "
      << matrixSize << " x " << matrixSize << flush;

  const std::vector<Eigen::Vector3d>& gridpoints = grid.getGridPositions();
  const Eigen::VectorXd& potential = grid.getGridValues();
  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Using " << atomlist.size() << " Fittingcenters and "
      << gridpoints.size() << " Gridpoints." << flush;

  Eigen::MatrixXd Amat = Eigen::MatrixXd::Zero(matrixSize, matrixSize);
  Eigen::VectorXd Bvec = Eigen::VectorXd::Zero(matrixSize);
// setting up _Amat
#pragma omp parallel for
  for (Index i = 0; i < atomlist.size(); i++) {
    for (Index j = i; j < atomlist.size(); j++) {
      for (const auto& gridpoint : gridpoints) {
        double dist_i = (atomlist[i].getPos() - gridpoint).norm();
        double dist_j = (atomlist[j].getPos() - gridpoint).norm();

        Amat(i, j) += 1.0 / dist_i / dist_j;
      }
      Amat(j, i) = Amat(i, j);
    }
  }

  // setting up Bvec
#pragma omp parallel for
  for (Index i = 0; i < atomlist.size(); i++) {
    for (Index j = 0; j < Index(gridpoints.size()); j++) {
      double dist_i = (atomlist[i].getPos() - gridpoints[j]).norm();
      Bvec(i) += potential(j) / dist_i;
    }
  }
  // Total charge constraint
  for (Index i = 0; i < atomlist.size() + 1; i++) {
    Amat(i, atomlist.size()) = 1.0;
    Amat(atomlist.size(), i) = 1.0;
  }
  Amat(atomlist.size(), atomlist.size()) = 0.0;
  Bvec(atomlist.size()) = netcharge;  // netcharge!!!!

  // Pairconstraint
  for (Index i = 0; i < Index(_pairconstraint.size()); i++) {
    const std::pair<Index, Index>& pair = _pairconstraint[i];
    Amat(pair.first, atomlist.size() + 1 + i) = 1.0;
    Amat(atomlist.size() + 1 + i, pair.first) = 1.0;
    Amat(pair.second, atomlist.size() + 1 + i) = -1.0;
    Amat(atomlist.size() + 1 + i, pair.second) = -1.0;
  }

  // Regionconstraint
  for (Index i = 0; i < Index(_regionconstraint.size()); i++) {
    const QMFragment<double>& reg = _regionconstraint[i];
    for (Index index : reg) {
      Amat(index, atomlist.size() + i + 1 + _pairconstraint.size()) = 1.0;
      Amat(atomlist.size() + i + 1 + _pairconstraint.size(), index) = 1.0;
    }
    Bvec(atomlist.size() + i + 1 + _pairconstraint.size()) = reg.value();
  }

  XTP_LOG_SAVE(Log::info, _log)
      << TimeStamp() << " Solving linear Equation " << flush;
  Eigen::VectorXd charges;
  if (_do_svd) {
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    svd.setThreshold(_conditionnumber);
    svd.compute(Amat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    charges = svd.solve(Bvec);
    XTP_LOG_SAVE(Log::info, _log) << TimeStamp() << " SVD Done. " << flush;
    if ((Bvec.size() - svd.nonzeroSingularValues()) != 0) {
      XTP_LOG_SAVE(Log::error, _log)
          << TimeStamp() << Bvec.size() - svd.nonzeroSingularValues()
          << " Sites could not be fitted and are set to zero." << flush;
    }
  } else {
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QR(Amat);
    if (QR.info() != Eigen::ComputationInfo::Success) {
      throw std::runtime_error(
          "Espfit: Solving the constrained equation failed. Maybe try SVD.");
    }
    charges = QR.solve(Bvec);
    XTP_LOG_SAVE(Log::info, _log)
        << TimeStamp() << " Solved linear least square fit ." << flush;
  }
  // remove constraints from charges
  charges.conservativeResize(atomlist.size());
  StaticSegment seg =
      StaticSegment(orbitals.QMAtoms().getType(), orbitals.QMAtoms().getId());

  XTP_LOG_SAVE(Log::error, _log)
      << " Sum of fitted charges: " << charges.sum() << flush;
  for (Index i = 0; i < atomlist.size(); i++) {
    seg.push_back(StaticSite(atomlist[i], charges(i)));
  }
  // get RMSE
  double rmse = 0.0;
  double totalPotSq = 0.0;
  for (Index k = 0; k < Index(gridpoints.size()); k++) {
    double temp = 0.0;
    for (const StaticSite& atom : seg) {
      double dist = (gridpoints[k] - atom.getPos()).norm();
      temp += atom.getCharge() / dist;
    }
    rmse += (potential(k) - temp) * (potential(k) - temp);
    totalPotSq += potential(k) * potential(k);
  }
  XTP_LOG_SAVE(Log::error, _log)
      << " RMSE of fit:  " << std::sqrt(rmse / double(gridpoints.size()))
      << flush;
  XTP_LOG_SAVE(Log::error, _log)
      << " RRMSE of fit: " << std::sqrt(rmse / totalPotSq) << flush;

  return seg;
}

}  // namespace xtp
}  // namespace votca
