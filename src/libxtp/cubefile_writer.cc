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

#include <votca/xtp/amplitude_integration.h>
#include <votca/xtp/cubefile_writer.h>
#include <votca/xtp/density_integration.h>
#include <votca/xtp/regular_grid.h>
namespace votca {
namespace xtp {

std::vector<std::vector<double>> CubeFile_Writer::CalculateValues(
    const Orbitals& orb, QMState state, bool dostateonly,
    const Regular_Grid& grid) const {

  bool do_amplitude = (state.Type().isSingleParticleState());
  if (do_amplitude) {
    Eigen::VectorXd amplitude;
    if (state.Type() == QMStateType::DQPstate) {
      Index amplitudeindex = state.StateIdx() - orb.getGWAmin();
      amplitude = orb.CalculateQParticleAORepresentation().col(amplitudeindex);
    } else {
      amplitude = orb.MOs().eigenvectors().col(state.StateIdx());
    }
    AmplitudeIntegration<Regular_Grid> ampl(grid);
    return ampl.IntegrateAmplitude(amplitude);
  } else {
    Eigen::MatrixXd mat;
    if (state.Type().isExciton() && dostateonly) {
      std::array<Eigen::MatrixXd, 2> DMAT =
          orb.DensityMatrixExcitedState(state);
      mat = DMAT[1] - DMAT[0];
    } else {
      mat = orb.DensityMatrixFull(state);
    }
    DensityIntegration<Regular_Grid> density(grid);
    density.IntegrateDensity(mat);
    return density.getDensities();
  }
}

std::vector<double> FlattenValues(
    const std::vector<std::vector<double>>& gridvalues) {
  Index size = 0;
  for (const auto& box : gridvalues) {
    size += Index(box.size());
  }
  std::vector<double> result;
  result.reserve(size);
  for (const auto& box : gridvalues) {
    result.insert(result.end(), box.begin(), box.end());
  }
  return result;
}

void CubeFile_Writer::WriteFile(const std::string& filename,
                                const Orbitals& orb, QMState state,
                                bool dostateonly) const {

  Regular_Grid grid;
  Eigen::Array3d padding = Eigen::Array3d::Ones() * _padding;
  AOBasis basis = orb.SetupDftBasis();
  XTP_LOG(Log::info, _log) << " Loaded DFT Basis Set " << orb.getDFTbasisName()
                           << std::flush;
  grid.GridSetup(_steps, padding, orb.QMAtoms(), basis);
  XTP_LOG(Log::info, _log) << " Calculating Gridvalues " << std::flush;
  auto temp = CalculateValues(orb, state, dostateonly, grid);
  std::vector<double> grid_values = FlattenValues(temp);
  XTP_LOG(Log::info, _log) << " Calculated Gridvalues " << std::flush;
  bool do_amplitude = (state.Type().isSingleParticleState());
  std::ofstream out(filename);
  if (!out.is_open()) {
    throw std::runtime_error("Bad file handle: " + filename);
  }

  // write cube header
  if (state.isTransition()) {
    out << boost::format("Transition state: %1$s \n") % state.ToString();
  } else if (do_amplitude) {
    out << boost::format("%1$s with energy %2$f eV \n") % state.ToString() %
               (orb.getExcitedStateEnergy(state) * tools::conv::hrt2ev);
  } else {
    if (dostateonly) {
      out << boost::format(
                 "Difference electron density of excited state %1$s \n") %
                 state.ToString();
    } else {
      out << boost::format("Total electron density of %1$s state\n") %
                 state.ToString();
    }
  }
  Eigen::Vector3d start = grid.getStartingPoint();
  out << "Created by VOTCA-XTP \n";
  if (do_amplitude) {
    out << boost::format("-%1$lu %2$f %3$f %4$f \n") % orb.QMAtoms().size() %
               start.x() % start.y() % start.z();
  } else {
    out << boost::format("%1$lu %2$f %3$f %4$f \n") % orb.QMAtoms().size() %
               start.x() % start.y() % start.z();
  }
  Eigen::Array<Index, 3, 1> steps = grid.getSteps();
  Eigen::Array3d stepsizes = grid.getStepSizes();
  out << boost::format("%1$d %2$f 0.0 0.0 \n") % steps.x() % stepsizes.x();
  out << boost::format("%1$d 0.0 %2$f 0.0 \n") % steps.y() % stepsizes.y();
  out << boost::format("%1$d 0.0 0.0 %2$f \n") % steps.z() % stepsizes.z();
  tools::Elements elements;
  for (const QMAtom& atom : orb.QMAtoms()) {
    double x = atom.getPos().x();
    double y = atom.getPos().y();
    double z = atom.getPos().z();
    std::string element = atom.getElement();
    Index atnum = elements.getEleNum(element);
    Index crg = atom.getNuccharge();
    out << boost::format("%1$d %2$d %3$f %4$f %5$f\n") % atnum % crg % x % y %
               z;
  }

  if (do_amplitude) {
    out << boost::format("  1 %1$d \n") % (state.StateIdx() + 1);
  }

  Eigen::TensorMap<Eigen::Tensor<double, 3>> gridvalues(
      grid_values.data(), steps.z(), steps.y(), steps.x());
  for (Index ix = 0; ix < steps.x(); ix++) {
    for (Index iy = 0; iy < steps.y(); iy++) {
      Index Nrecord = 0;
      for (Index iz = 0; iz < steps.z(); iz++) {
        Nrecord++;
        if (Nrecord == 6 || iz == (steps.z() - 1)) {
          out << boost::format("%1$E \n") % gridvalues(iz, iy, ix);
          Nrecord = 0;
        } else {
          out << boost::format("%1$E ") % gridvalues(iz, iy, ix);
        }
      }
    }
  }
  out.close();
}

}  // namespace xtp
}  // namespace votca
