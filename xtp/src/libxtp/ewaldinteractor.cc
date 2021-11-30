

#include "votca/xtp/ewaldinteractor.h"

namespace votca {
namespace xtp {

Eigen::Vector3d EwaldInteractor::r_staticFieldAtBy(
    const BGSite& site, const BGSite& nbSite, const Eigen::Vector3d& shift) {
  Eigen::Vector3d dr = (nbSite.getPos() - site.getPos()) + shift;
  return Eigen::Vector3d::Zero();

}

}  // namespace xtp

}  // namespace votca