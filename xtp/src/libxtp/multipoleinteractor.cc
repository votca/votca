

#include "votca/xtp/multipoleinteractor.h"

namespace votca {
namespace xtp {

inline double dyadic3d(const Eigen::Matrix3d& mat1, const Eigen::Matrix3d& mat2) {
  return (mat1.array() * mat2.array()).sum();
}

std::array<double, 5> MultipoleInteractor::orientationDependence(
    Multipole mp1, Multipole mp2, Eigen::Vector3d dr) {
  std::array<double, 5> results {0, 0, 0, 0, 0};
  Index maxRank = std::max(mp1.rank, mp2.rank);
  // always add the partical charge contributions
  results[0] = (mp1.charge * mp2.charge);
  if (maxRank > 0) {  // add all contributions containing dipoles
    results[1] = mp1.dipole.dot(mp2.dipole) +
                 (mp1.charge * mp2.dipole - mp2.charge * mp1.dipole).dot(dr);
    results[2] = -mp1.dipole.dot(dr) * mp2.dipole.dot(dr);
    if (maxRank > 1) {  // we add all contributions with quadrupoles
      Eigen::Matrix3d dr_outer = dr * dr.transpose();
      Eigen::Matrix3d chargeQuad =
          (mp1.charge * mp2.quadrupole + mp2.charge * mp1.quadrupole);
      Eigen::Matrix3d dip1R = mp1.dipole * dr.transpose();
      Eigen::Matrix3d dip2R = mp2.dipole * dr.transpose();
      // += because we did the first part for maxRank = 0
      results[2] += 2 * dyadic3d(mp1.quadrupole, mp2.quadrupole) +
                    dyadic3d(chargeQuad, dr_outer) -
                    2 * dyadic3d(dip2R, mp1.quadrupole) +
                    2 * dyadic3d(dip1R, mp2.quadrupole);
      results[3] = -(mp1.dipole.dot(dr))*dyadic(mp2.quadrupole, dr_outer) +
                    (mp2.dipole.dot(dr))*dyadic(mp1.quadrupole, dr_outer) -
                    4 * dr.transpose() * mp1.quadrupole * mp2.quadrupole * dr;
      results[4] = dyadic3d(mp1.quadrupole, dr_outer) *
                    dyadic3d(mp2.quadrupole, dr_outer);
    }
  }
  return results
}

}  // namespace xtp
}  // namespace votca