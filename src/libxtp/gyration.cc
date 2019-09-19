/*
 *            Copyright 2016 The MUSCET Development Team
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

#include <boost/format.hpp>
#include <votca/tools/elements.h>
#include <votca/xtp/gyration.h>

using namespace std;
using namespace votca::tools;

namespace votca {
namespace xtp {

void Density2Gyration::Initialize(tools::Property& options) {
  string key = Identify();

  std::string statestring =
      options.ifExistsReturnElseThrowRuntimeError<string>(key + ".state");
  _state.FromString(statestring);
  _dostateonly = options.ifExistsReturnElseReturnDefault<bool>(
      key + ".difference_to_groundstate", false);
  _gridsize = options.ifExistsReturnElseReturnDefault<string>(key + ".gridsize",
                                                              "medium");
}

void Density2Gyration::AnalyzeDensity(const Orbitals& orbitals) {
  XTP_LOG_SAVE(logDEBUG, _log) << "===== Running on " << OPENMP::getMaxThreads()
                               << " threads ===== " << flush;

  const QMMolecule& Atomlist = orbitals.QMAtoms();
  Eigen::MatrixXd DMAT_tot;
  BasisSet bs;
  bs.Load(orbitals.getDFTbasisName());
  AOBasis basis;
  basis.Fill(bs, Atomlist);
  AnalyzeGeometry(Atomlist);
  std::vector<Eigen::MatrixXd> DMAT;

  // setup numerical integration grid
  NumericalIntegration numway;
  numway.GridSetup(_gridsize, Atomlist, basis);

  if (!_dostateonly) {
    Eigen::MatrixXd DMATGS = orbitals.DensityMatrixFull(_state);
    Gyrationtensor gyro = numway.IntegrateGyrationTensor(DMAT_tot);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.computeDirect(gyro.gyration);
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " Converting to Eigenframe " << flush;
    XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << " Reporting " << flush;
    ReportAnalysis(_state.ToLongString(), gyro, es);

  } else {
    // hole density first
    std::array<Eigen::MatrixXd, 2> DMAT =
        orbitals.DensityMatrixExcitedState(_state);
    Gyrationtensor gyro_hole = numway.IntegrateGyrationTensor(DMAT[0]);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es_h;
    es_h.computeDirect(gyro_hole.gyration);
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " Converting to Eigenframe " << flush;
    XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << " Reporting " << flush;
    ReportAnalysis("hole", gyro_hole, es_h);

    // electron density
    Gyrationtensor gyro_electron = numway.IntegrateGyrationTensor(DMAT[1]);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es_e;
    es_e.computeDirect(gyro_electron.gyration);
    XTP_LOG_SAVE(logDEBUG, _log)
        << TimeStamp() << " Converting to Eigenframe " << flush;
    XTP_LOG_SAVE(logDEBUG, _log) << TimeStamp() << " Reporting " << flush;
    ReportAnalysis("electron", gyro_electron, es_e);
  }
  return;
}

void Density2Gyration::AnalyzeGeometry(const QMMolecule& atoms) {

  tools::Elements elements;
  double mass = 0.0;
  Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
  Eigen::Matrix3d gyration = Eigen::Matrix3d::Zero();
  for (const QMAtom& atom : atoms) {
    double m = elements.getMass(atom.getElement());
    const Eigen::Vector3d& pos = atom.getPos();
    mass += m;
    centroid += m * pos;
    gyration += m * pos * pos.transpose();
  }
  centroid /= mass;
  gyration /= mass;
  gyration -= centroid * centroid.transpose();
  Gyrationtensor gyro;
  gyro.mass = mass;
  gyro.centroid = centroid;
  gyro.gyration = gyration;
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
  es.computeDirect(gyro.gyration);
  ReportAnalysis("geometry", gyro, es);
}

void Density2Gyration::ReportAnalysis(
    string label, const Gyrationtensor& gyro,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& es) {

  XTP_LOG_SAVE(logINFO, _log)
      << "---------------- " << label << " ----------------" << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Norm               = %1$9.4f ") % (gyro.mass))
      << flush;

  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Centroid x         = %1$9.4f Ang") %
          (gyro.centroid.x() * tools::conv::bohr2ang))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Centroid y         = %1$9.4f Ang") %
          (gyro.centroid.y() * tools::conv::bohr2ang))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Centroid y         = %1$9.4f Ang") %
          (gyro.centroid.z() * tools::conv::bohr2ang))
      << flush;

  double RA2 = tools::conv::bohr2ang * tools::conv::bohr2ang;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor xx = %1$9.4f Ang^2") %
          (gyro.gyration(0, 0) * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor xy = %1$9.4f Ang^2") %
          (gyro.gyration(0, 1) * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor xz = %1$9.4f Ang^2") %
          (gyro.gyration(0, 2) * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor yy = %1$9.4f Ang^2") %
          (gyro.gyration(1, 1) * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor yz = %1$9.4f Ang^2") %
          (gyro.gyration(1, 2) * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor zz = %1$9.4f Ang^2") %
          (gyro.gyration(2, 2) * RA2))
      << flush;

  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor D1 = %1$9.4f Ang^2") %
          (es.eigenvalues()[0] * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor D2 = %1$9.4f Ang^2") %
          (es.eigenvalues()[1] * RA2))
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Gyration Tensor D3 = %1$9.4f Ang^2") %
          (es.eigenvalues()[2] * RA2))
      << flush;

  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Radius of Gyration = %1$9.4f Ang") %
          (std::sqrt(es.eigenvalues().sum()) * tools::conv::bohr2ang))
      << flush;

  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 1 1 = %1$9.4f ") %
          es.eigenvectors().col(0).x())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 1 2 = %1$9.4f ") %
          es.eigenvectors().col(0).y())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 1 3 = %1$9.4f ") %
          es.eigenvectors().col(0).z())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 2 1 = %1$9.4f ") %
          es.eigenvectors().col(1).x())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 2 2 = %1$9.4f ") %
          es.eigenvectors().col(1).y())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 2 3 = %1$9.4f ") %
          es.eigenvectors().col(1).z())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 3 1 = %1$9.4f ") %
          es.eigenvectors().col(2).x())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 3 2 = %1$9.4f ") %
          es.eigenvectors().col(2).y())
      << flush;
  XTP_LOG_SAVE(logINFO, _log)
      << (boost::format("  Tensor EF Axis 3 3 = %1$9.4f ") %
          es.eigenvectors().col(2).z())
      << flush;
  return;
}

}  // namespace xtp
}  // namespace votca
