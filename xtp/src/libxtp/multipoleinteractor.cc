

#include "votca/xtp/multipoleinteractor.h"

namespace votca {
namespace xtp {

inline double dyadic3d(const Eigen::Matrix3d& mat1,
                       const Eigen::Matrix3d& mat2) {
  return (mat1.array() * mat2.array()).sum();
}

MultipoleInteractor::MultipoleInteractor(double alpha, double thole_damping) {
  a1 = alpha;
  a2 = alpha * alpha;
  a3 = a2 * a1;
  a4 = a2 * a2;
  a5 = a4 * a1;
  thole = thole_damping;
  thole2 = thole * thole;
  thole3 = thole * thole2;
}

std::array<double, 5> MultipoleInteractor::orientationDependence(
    const Multipole& mp1, const Multipole& mp2,
    const Eigen::Vector3d& dr) const {
  std::array<double, 5> results{0, 0, 0, 0, 0};
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
      results[3] = -(mp1.dipole.dot(dr)) * dyadic3d(mp2.quadrupole, dr_outer) +
                   (mp2.dipole.dot(dr)) * dyadic3d(mp1.quadrupole, dr_outer) -
                   4 * dr.transpose() * mp1.quadrupole * mp2.quadrupole * dr;
      results[4] = dyadic3d(mp1.quadrupole, dr_outer) *
                   dyadic3d(mp2.quadrupole, dr_outer);
    }
  }
  return results;
}

std::array<Eigen::Vector3d, 3> MultipoleInteractor::orientationDependence(
    const Multipole& mp1, const Eigen::Vector3d& dr) const {
  std::array<Eigen::Vector3d, 3> results;
  results[0] = -mp1.charge * dr;
  if (mp1.rank > 0) {
    results[0] += mp1.dipole;
    results[1] = -(mp1.dipole.dot(dr)) * dr;
    if (mp1.rank > 1) {
      results[1] += 2 * mp1.quadrupole * dr;
      Eigen::Matrix3d dr_outer = dr * dr.transpose();
      results[2] = -dyadic3d(mp1.quadrupole, dr_outer) * dr;
    }
  }
  return results;
}

std::array<double, 4> MultipoleInteractor::erfcDistanceDependence(
    const Eigen::Vector3d& dr) const {
  std::array<double, 4> results{0, 0, 0, 0};

  double R = dr.norm();
  double R2 = R * R;
  double rR1 = 1.0 / R;
  double rR2 = rR1 * rR1;

  double rSqrtPiExp = rSqrtPi * std::exp(-a2 * R2);

  results[0] = std::erfc(a1 * R) * rR1;
  results[1] = rR2 * (results[0] + 2.0 * a1 * rSqrtPiExp);
  results[2] = rR2 * (3.0 * results[1] + 4.0 * a3 * rSqrtPiExp);
  results[3] = rR2 * (5.0 * results[2] + 8.0 * a5 * rSqrtPiExp);

  return results;
}

std::array<double, 4> MultipoleInteractor::erfDistanceDependence(
    const Eigen::Vector3d& dr) const {
  std::array<double, 4> results{0, 0, 0, 0};

  double R = dr.norm();
  double R2 = R * R;
  double rR1 = 1.0 / R;
  double rR2 = rR1 * rR1;

  double rSqrtPiExp = rSqrtPi * std::exp(-a2 * R2);

  results[0] = std::erf(a1 * R) * rR1;
  results[1] = rR2 * (results[0] - 2.0 * a1 * rSqrtPiExp);
  results[2] = rR2 * (3.0 * results[1] - 4.0 * a3 * rSqrtPiExp);
  results[3] = rR2 * (5.0 * results[2] - 8.0 * a5 * rSqrtPiExp);

  return results;
}

std::array<double, 4> MultipoleInteractor::tholeDamping(
    const Eigen::Matrix3d& pol1, const Eigen::Matrix3d& pol2,
    const Eigen::Vector3d& dr) const {
  double R = dr.norm();
  double R2 = R * R;

  std::array<double, 4> results{1.0, 1.0, 1.0, 1.0};

  double thole_u3 =
      (R * R2) / std::sqrt((1.0 / 3.0) * (pol1.array() * pol2.array()).sum());

  if (thole * thole_u3 < 40) {
    double thole_exp = std::exp(-thole * thole_u3);
    double thole_u6 = thole_u3 * thole_u3;
    results[0] = 1 - thole_exp;
    results[1] = 1 - (1 + thole * thole_u3) * thole_exp;
    results[2] = 1 - (1 + thole * thole_u3 + (3. / 5.) * thole2 * thole_u6);
    results[3] = 1 - (1 + thole * thole_u3 + (18. / 35.) * thole2 * thole_u6 +
                      (9. / 35.) * thole3 * thole_u6 * thole_u3) *
                         thole_exp;
  }  // else: the exponent is close to zero and hence the expression is close to
     // 1 and we return an array of ones.
  return results;
}

}  // namespace xtp

}  // namespace votca