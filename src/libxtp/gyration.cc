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

// Third party includes
#include <boost/format.hpp>

// VOTCA includes
#include <votca/tools/elements.h>

// Local VOTCA includes
#include "votca/xtp/gyration.h"
#include "votca/xtp/vxc_grid.h"

using namespace votca::tools;

namespace votca {
namespace xtp {

void Density2Gyration::Initialize(tools::Property& options) {
  std::string key = Identify();
  state_ = options.get(key + ".state").as<QMState>();
  dostateonly_ = options.ifExistsReturnElseReturnDefault<bool>(
      key + ".difference_to_groundstate", false);
  gridsize_ = options.ifExistsReturnElseReturnDefault<std::string>(
      key + ".gridsize", "medium");
}

void Density2Gyration::AnalyzeDensity(const Orbitals& orbitals) {
  XTP_LOG(Log::error, log_) << "===== Running on " << OPENMP::getMaxThreads()
                            << " threads ===== " << std::flush;

  const QMMolecule& Atomlist = orbitals.QMAtoms();
  BasisSet bs;
  bs.Load(orbitals.getDFTbasisName());
  AOBasis basis;
  basis.Fill(bs, Atomlist);
  AnalyzeGeometry(Atomlist);

  // setup numerical integration grid
  Vxc_Grid grid;
  grid.GridSetup(gridsize_, Atomlist, basis);
  DensityIntegration<Vxc_Grid> numway(grid);

  if (!dostateonly_) {
    Eigen::MatrixXd DMAT_tot = orbitals.DensityMatrixFull(state_);
    Gyrationtensor gyro = numway.IntegrateGyrationTensor(DMAT_tot);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
    es.computeDirect(gyro.gyration);
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Converting to Eigenframe " << std::flush;
    XTP_LOG(Log::error, log_) << TimeStamp() << " Reporting " << std::flush;
    ReportAnalysis(state_.ToLongString(), gyro, es);

  } else {
    // hole density first
    std::array<Eigen::MatrixXd, 2> DMAT =
        orbitals.DensityMatrixExcitedState(state_);
    Gyrationtensor gyro_hole = numway.IntegrateGyrationTensor(DMAT[0]);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es_h;
    es_h.computeDirect(gyro_hole.gyration);
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Converting to Eigenframe " << std::flush;
    XTP_LOG(Log::error, log_) << TimeStamp() << " Reporting " << std::flush;
    ReportAnalysis("hole", gyro_hole, es_h);

    // electron density
    Gyrationtensor gyro_electron = numway.IntegrateGyrationTensor(DMAT[1]);
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es_e;
    es_e.computeDirect(gyro_electron.gyration);
    XTP_LOG(Log::error, log_)
        << TimeStamp() << " Converting to Eigenframe " << std::flush;
    XTP_LOG(Log::error, log_) << TimeStamp() << " Reporting " << std::flush;
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
    std::string label, const Gyrationtensor& gyro,
    const Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d>& es) {

  XTP_LOG(Log::error, log_)
      << "---------------- " << label << " ----------------" << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Norm               = %1$9.4f ") % (gyro.mass))
      << std::flush;

  XTP_LOG(Log::error, log_)
      << (boost::format("  Centroid x         = %1$9.4f Ang") %
          (gyro.centroid.x() * tools::conv::bohr2ang))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Centroid y         = %1$9.4f Ang") %
          (gyro.centroid.y() * tools::conv::bohr2ang))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Centroid y         = %1$9.4f Ang") %
          (gyro.centroid.z() * tools::conv::bohr2ang))
      << std::flush;

  double RA2 = tools::conv::bohr2ang * tools::conv::bohr2ang;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor xx = %1$9.4f Ang^2") %
          (gyro.gyration(0, 0) * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor xy = %1$9.4f Ang^2") %
          (gyro.gyration(0, 1) * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor xz = %1$9.4f Ang^2") %
          (gyro.gyration(0, 2) * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor yy = %1$9.4f Ang^2") %
          (gyro.gyration(1, 1) * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor yz = %1$9.4f Ang^2") %
          (gyro.gyration(1, 2) * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor zz = %1$9.4f Ang^2") %
          (gyro.gyration(2, 2) * RA2))
      << std::flush;

  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor D1 = %1$9.4f Ang^2") %
          (es.eigenvalues()[0] * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor D2 = %1$9.4f Ang^2") %
          (es.eigenvalues()[1] * RA2))
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Gyration Tensor D3 = %1$9.4f Ang^2") %
          (es.eigenvalues()[2] * RA2))
      << std::flush;

  XTP_LOG(Log::error, log_)
      << (boost::format("  Radius of Gyration = %1$9.4f Ang") %
          (std::sqrt(es.eigenvalues().sum()) * tools::conv::bohr2ang))
      << std::flush;

  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 1 1 = %1$9.4f ") %
          es.eigenvectors().col(0).x())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 1 2 = %1$9.4f ") %
          es.eigenvectors().col(0).y())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 1 3 = %1$9.4f ") %
          es.eigenvectors().col(0).z())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 2 1 = %1$9.4f ") %
          es.eigenvectors().col(1).x())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 2 2 = %1$9.4f ") %
          es.eigenvectors().col(1).y())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 2 3 = %1$9.4f ") %
          es.eigenvectors().col(1).z())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 3 1 = %1$9.4f ") %
          es.eigenvectors().col(2).x())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 3 2 = %1$9.4f ") %
          es.eigenvectors().col(2).y())
      << std::flush;
  XTP_LOG(Log::error, log_)
      << (boost::format("  Tensor EF Axis 3 3 = %1$9.4f ") %
          es.eigenvectors().col(2).z())
      << std::flush;
}

}  // namespace xtp
}  // namespace votca
