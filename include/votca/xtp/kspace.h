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
#ifndef VOTCA_XTP_KSPACE_H
#define VOTCA_XTP_KSPACE_H

#include "votca/xtp/segid.h"
#include <boost/math/constants/constants.hpp>
#include <complex>
#include <iomanip>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/ewaldoptions.h"
#include "votca/xtp/ewd_segment.h"
#include "votca/xtp/unitcell.h"

namespace votca {
namespace xtp {

class KVector;

class KSpace {
 public:
  KSpace(Logger& log) : _log(log) { ; }
  ~KSpace() = default;

  void Initialize(const EwaldOptions& options, const UnitCell& unitcell,
                  std::vector<EwdSegment>& segments) {
    this->options = options;
    this->_unit_cell = unitcell;
    a1 = options.alpha;
    a2 = a1 * a1;
    a3 = a1 * a2;
    a4 = a2 * a2;
    a5 = a4 * a1;
    fourPiVolume =
        4.0 * boost::math::constants::pi<double>() / _unit_cell.getVolume();
    cutoff = options.k_cutoff;
    cutoff2 = cutoff * cutoff;

    _unit_cell = unitcell;

    // compute max k-space vectors
    const Eigen::Matrix3d& inverseCellMatrix = _unit_cell.getInverseMatrix();
    for (Index i = 0; i < 3; ++i) {
      max_K[i] = static_cast<Index>(
          std::ceil(cutoff / inverseCellMatrix.col(i).norm()));
    }

    XTP_LOG(Log::error, _log)
        << "************* KSPACE: PARAMETERS *************" << std::endl;
    XTP_LOG(Log::error, _log)
        << "kspace cutoff: " << cutoff << "a.u. (" << (1 / 0.05291) * cutoff
        << " nm-1)" << std::endl;
    switch (options.shape) {
      case Shape::sphere:
        XTP_LOG(Log::error, _log) << "shape: sphere" << std::endl;
        break;
      case Shape::cube:
        XTP_LOG(Log::error, _log) << "shape: cube" << std::endl;
        break;
      case Shape::xyslab:
        XTP_LOG(Log::error, _log) << "shape: xyslab" << std::endl;
        break;
    }

    XTP_LOG(Log::error, _log)
        << "Max K copies: [" << max_K[0] << ", " << max_K[1] << ", " << max_K[2]
        << "]" << std::endl
        << std::endl;

    systemSize = 0;
    for (const auto& seg : _ewaldSegments) {
      segmentOffSet.push_back(systemSize);
      systemSize += 3 * seg.size();
    }

    _ewaldSegments = segments;

    systemSize = 0;
    for (const auto& seg : _ewaldSegments) {
      segmentOffSet.push_back(systemSize);
      systemSize += 3 * seg.size();
    }
    // precompute the K-Vectors
    computeKVectors();
  }

  void computeStaticField();
  void applyStaticShapeField();
  Eigen::Vector3d computeShapeField();
  void computeIntraMolecularCorrection();

  // Total Field Calculations
  void computeTotalField(PolarSegment& seg);
  void applyTotalShapeField(PolarSegment& seg);
  void applySICorrection(PolarSegment& seg);
  void applyAPeriodicCorrection(PolarSegment& seg,
                                std::vector<SegId> pCloud_indices);

  void addInducedDipoleInteractionTo(Eigen::MatrixXd& result);
  void addShapeCorrectionTo(Eigen::MatrixXd& result);
  void addSICorrectionTo(Eigen::MatrixXd& result);

  double backgroundInteractionEnergy() { return 0.0; }

  double selfInteractionEnergy() { return 0.0; }

  double aPeriodicCorrectionEnergy() { return 0.0; }

  double shapeCorrectionEnergy() { return 0.0; }

 private:
  void computeTholeVariables(const Eigen::Matrix3d& pol1,
                             const Eigen::Matrix3d& pol2);
  std::complex<double> computeSk(const Eigen::Vector3d& kvector) const;
  double computeAk(const Eigen::Vector3d& kvector) const;
  void computeKVectors();
  Eigen::VectorXcd getSkInteractionVector(const Eigen::Vector3d& kvector);

  void computeScreenedInteraction();
  void computeDistanceVariables(Eigen::Vector3d distVec);

  Eigen::Vector3d staticFieldAtBy(const EwdSite& site, const EwdSite& nbSite);
  Eigen::Vector3d totalFieldAtBy(const EwdSite& site, const EwdSite& nbSite);

  Eigen::Matrix3d inducedDipoleInteractionAtBy(EwdSite& site,
                                               const EwdSite& nbSite);
  Eigen::Vector3d kspaceCorrectionFieldAtBy(EwdSite& site,
                                            const EwdSite& nbSite);

  std::vector<Index> segmentOffSet;

  Index systemSize;

  double a1, a2, a3, a4, a5;  // alpha (splitting param) and its powers
  double l3, l5, l7, l9;
  double thole, thole2, thole3, thole_u3;
  UnitCell _unit_cell;
  std::vector<EwdSegment> _ewaldSegments;
  std::vector<KVector> _kvector_list;
  double fourPiVolume;
  double cutoff, cutoff2;
  Eigen::Vector3d dr = Eigen::Vector3d::Zero();
  const std::complex<double> ii =
      std::complex<double>(0.0, 1.0);  // imaginary i
  std::array<Index, 3> max_K;
  // rRns = reciprocal R, of order n, screened with erfc
  double rR1s, rR3s, rR5s, rR7s;
  double R1, R2;    // distance and powers
  double rR1, rR2;  // reciprocal (i.e. 1.0/ ...) distance and powers
  double pi = boost::math::constants::pi<double>();
  double rSqrtPi = 1.0 / std::sqrt(pi);
  Logger& _log;
  EwaldOptions options;
};

/**
 * \brief Class that contains everything related to a single k-vector:
 *  - k-vector
 *  - A(k) value
 *  - S(k) value (the structure factor)
 * Implements comparison operators for std::sort, such that vectors can be
 * sorted based on importance, current ordering is the distance from central
 * cell.
 */
class KVector {
 public:
  KVector(const Eigen::Vector3d& kvector) : _kvector(kvector){};
  ~KVector() = default;

  const Eigen::Vector3d& getVector() const { return _kvector; }
  double getAk() const { return _Ak; }
  std::complex<double> getSk() const { return _Sk; }

  void setAk(double Ak) { _Ak = Ak; }
  void setSk(std::complex<double> Sk) { _Sk = Sk; }

  friend std::ostream& operator<<(std::ostream& out, const KVector& kvector) {
    out << std::scientific << std::setprecision(5)
        << "vec: " << kvector.getVector().transpose()
        << " Ak: " << kvector.getAk() << " Sk: " << kvector.getSk();
    return out;
  }

  bool operator<(const KVector& other) {
    if (this->_kvector.norm() < other.getVector().norm()) {
      return true;
    }
    return false;
  }

  bool operator>(const KVector& other) {
    if (this->_kvector.norm() > other.getVector().norm()) {
      return true;
    }
    return false;
  }

  bool operator==(const KVector& other) {
    if (this->_kvector.norm() == other.getVector().norm()) {
      return true;
    }
    return false;
  }

 private:
  Eigen::Vector3d _kvector;
  double _Ak = 0;
  std::complex<double> _Sk;
};

}  // namespace xtp
}  // namespace votca

#endif