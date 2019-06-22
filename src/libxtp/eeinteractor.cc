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
#include "votca/xtp/eeinteractor.h"
#include <iostream>

namespace votca {
namespace xtp {

Eigen::Matrix4d eeInteractor::FillInteraction_noQuadrupoles(
    const StaticSite& site1, const StaticSite& site2) const {
  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();
  int rankA = site1.getRank();
  int rankB = site2.getRank();
  Eigen::Matrix4d interaction = Eigen::Matrix4d::Zero();
  const Eigen::Vector3d r_AB =
      posB - posA;               // Vector of the distance between polar sites
  const double R = r_AB.norm();  // Norm of distance vector
  const Eigen::Vector3d pos_a = r_AB / R;
  const double fac0 = 1 / R;

  // Charge-Charge Interaction
  interaction(0, 0) = fac0;  // T_00,00
  const double fac1 = std::pow(fac0, 2);
  const double fac2 = std::pow(fac0, 3);

  // Dipole-Charge Interaction
  interaction.block<3, 1>(1, 0) = fac1 * pos_a;  // T_1alpha,00 (alpha=x,y,z)
  // Charge-Dipole Interaction
  interaction.block<1, 3>(0, 1) =
      -fac1 * pos_a.transpose();  // T_00,1alpha (alpha=x,y,z)

  if (rankA > 0 && rankB > 0) {
    // Dipole-Dipole Interaction
    interaction.block<3, 3>(1, 1) =
        fac2 *
        (-3 * pos_a * pos_a.transpose() +
         Eigen::Matrix3d::Identity());  // T_1alpha,1beta (alpha,beta=x,y,z)
  }
  return interaction;
}

Matrix9d eeInteractor::FillInteraction(const StaticSite& site1,
                                       const StaticSite& site2) const {

  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();
  int rankA = site1.getRank();
  int rankB = site2.getRank();

  Matrix9d interaction = Matrix9d::Zero();
  const Eigen::Vector3d r_AB =
      posB - posA;               // Vector of the distance between polar sites
  const double R = r_AB.norm();  // Norm of distance vector
  const Eigen::Vector3d pos_a =
      r_AB /
      R;  // unit vector on the sites reciprocal direction; This points toward A
  const Eigen::Vector3d pos_b = -pos_a;  // unit vector on the sites reciprocal
                                         // direction; This points toward B

  const double fac0 = 1 / R;

  interaction.block<4, 4>(0, 0) = FillInteraction_noQuadrupoles(site1, site2);
  const double sqr3 = std::sqrt(3);
  const Eigen::Matrix3d AxA = pos_a * pos_a.transpose();
  const Eigen::Matrix3d& BxB = AxA;
  const double fac2 = std::pow(fac0, 3);

  // Quadrupole-Charge interaction
  Eigen::Matrix<double, 1, 5> block;
  block(0) = fac2 * 0.5 * (3 * AxA(2, 2) - 1);             // T20,00
  block(1) = fac2 * sqr3 * AxA(0, 2);                      // 21c,00
  block(2) = fac2 * sqr3 * AxA(1, 2);                      // T21s,000
  block(3) = fac2 * 0.5 * sqr3 * (AxA(0, 0) - AxA(1, 1));  // T22c,00
  block(4) = fac2 * sqr3 * AxA(0, 1);                      // T22s,00
  if (rankA > 1) {
    interaction.block<5, 1>(4, 0) = block.transpose();
  }
  if (rankB > 1) {
    interaction.block<1, 5>(0, 4) = block;
  }

  if (rankA > 0 && rankB > 0) {
    const Eigen::Matrix3d AxB = -AxA;
    const Eigen::Matrix3d c = Eigen::Matrix3d::Identity();
    const double fac3 = std::pow(fac0, 4);
    if (rankA > 1 || rankB > 1) {

      Eigen::Matrix<double, 3, 5> block;
      // Quadrupole-Dipole Interaction
      block.col(0) = 0.5 * fac3 *
                     (15 * AxA(2, 2) * pos_b + 6 * pos_a.z() * c.col(2) -
                      3 * pos_b);  // T20-1beta (beta=x,y,z)
      block.col(1) = fac3 * sqr3 *
                     (pos_a.x() * c.col(2) + c.col(0) * pos_a.z() +
                      5 * AxA(0, 2) * pos_b);  // T21c-1beta (beta=x,y,z)
      block.col(2) = fac3 * sqr3 *
                     (pos_a.y() * c.col(2) + c.col(1) * pos_a.z() +
                      5 * AxA(1, 2) * pos_b);  // T21s-1beta (beta=x,y,z)
      block.col(3) =
          fac3 * 0.5 * sqr3 *
          (5 * (AxA(0, 0) - AxA(1, 1)) * pos_b + 2 * pos_a.x() * c.col(0) -
           2 * pos_a.y() * c.col(1));  // T22c-1beta (beta=x,y,z)
      block.col(4) = fac3 * sqr3 *
                     (5 * AxA(0, 1) * pos_b + pos_a.x() * c.col(1) +
                      pos_a.y() * c.col(0));  // T22s-1beta (beta=x,y,z)

      if (rankA > 1) {
        interaction.block<5, 3>(4, 1) = block.transpose();
      }
      if (rankB > 1) {
        interaction.block<3, 5>(1, 4) = -block;
      }
    }

    if (rankA > 1 && rankB > 1) {
      const double fac4 = std::pow(fac0, 5);
      // Quadrupole-Quadrupole Interaction
      Eigen::Matrix<double, 5, 5> block;
      block(0, 0) =
          fac4 * (3. / 4.) *
          (35 * AxA(2, 2) * BxB(2, 2) - 5 * AxA(2, 2) - 5 * BxB(2, 2) +
           20 * AxB(2, 2) * c(2, 2) + 2 * c(2, 2) * c(2, 2) + 1);  // T20,20
      block(1, 0) = 0.5 * fac4 * sqr3 *
                    (35 * AxA(2, 2) * BxB(0, 2) - 5 * BxB(0, 2) +
                     10 * AxB(2, 0) * c(2, 2) + 10 * AxB(2, 2) * c(2, 1) +
                     2 * c(2, 0) * c(2, 2));  // T20,21c
      block(2, 0) = 0.5 * fac4 * sqr3 *
                    (35 * AxA(2, 2) * BxB(1, 2) - 5 * BxB(1, 2) +
                     10 * AxB(2, 1) * c(2, 2) + 10 * AxB(2, 2) * c(2, 1) +
                     2 * c(2, 1) * c(2, 2));  // T20,21s
      block(3, 0) =
          0.25 * fac4 * sqr3 *
          (35 * AxB(2, 0) - 35 * AxB(2, 1) - 5 * BxB(0, 0) + 5 * BxB(1, 1) +
           20 * AxB(2, 0) * c(2, 0) - 20 * AxB(2, 1) * c(2, 1) +
           2 * c(2, 0) * c(2, 0) - 2 * c(2, 1) * c(2, 1));  // T20,22c
      block(4, 0) = 0.5 * fac4 * sqr3 *
                    (35 * AxA(2, 2) * BxB(0, 1) - 5 * BxB(0, 1) +
                     10 * AxB(2, 0) * c(2, 1) + 10 * AxB(2, 1) * c(2, 0) +
                     2 * c(2, 0) * c(2, 1));  // T20,22s
      block(1, 1) = fac4 * (35 * AxA(0, 2) * BxB(0, 2) +
                            5 * AxB(0, 0) * c(2, 2) + 5 * AxB(0, 2) * c(2, 0) +
                            5 * AxB(2, 0) * c(0, 2) + 5 * AxB(2, 2) * c(0, 0) +
                            c(0, 0) * c(2, 2) + c(0, 2) * c(2, 0));  // T21c,21c
      block(2, 1) = fac4 * (35 * AxA(0, 2) * BxB(1, 2) +
                            5 * AxB(0, 1) * c(2, 2) + 5 * AxB(0, 2) * c(2, 1) +
                            5 * AxB(2, 1) * c(0, 2) + 5 * AxB(2, 2) * c(0, 1) +
                            c(0, 1) * c(2, 2) + c(0, 2) * c(2, 1));  // T21c,21s
      block(3, 1) =
          0.5 * fac4 *
          (35 * AxA(0, 2) * BxB(0, 0) - 35 * AxA(0, 2) * BxB(1, 1) +
           10 * AxB(0, 0) * c(2, 0) - 10 * AxB(0, 1) * c(2, 1) +
           10 * AxB(0, 0) * c(0, 0) - 10 * AxB(2, 1) * c(0, 1) +
           2 * c(0, 0) * c(2, 0) - 2 * c(0, 1) * c(2, 1));  // T21c,22c
      block(4, 1) = fac4 * (35 * AxA(0, 2) * BxB(0, 1) +
                            5 * AxB(0, 0) * c(2, 1) + 5 * AxB(0, 1) * c(2, 0) +
                            5 * AxB(2, 0) * c(0, 1) + 5 * AxB(2, 1) * c(0, 0) +
                            c(0, 0) * c(2, 1) + c(0, 1) * c(2, 0));  // T21c,22s
      block(2, 2) = fac4 * (35 * AxA(1, 2) * BxB(1, 2) +
                            5 * AxB(1, 1) * c(2, 2) + 5 * AxB(1, 2) * c(2, 1) +
                            5 * AxB(2, 1) * c(1, 2) + 5 * AxB(2, 2) * c(1, 1) +
                            c(1, 1) * c(2, 2) + c(1, 2) * c(2, 1));  // T21s,21s
      block(3, 2) =
          0.5 * fac4 *
          (35 * AxA(1, 2) * BxB(0, 0) - 35 * AxA(1, 2) * BxB(1, 1) +
           10 * AxB(1, 2) * c(2, 0) - 10 * AxB(1, 1) * c(2, 1) +
           10 * AxB(2, 0) * c(1, 0) - 10 * AxB(2, 1) * c(1, 1) +
           2 * c(1, 0) * c(2, 0) - 2 * c(1, 1) * c(2, 1));  // T21s,22c
      block(4, 2) = fac4 * (35 * AxA(1, 2) * BxB(0, 1) +
                            5 * AxB(1, 0) * c(2, 1) + 5 * AxB(1, 1) * c(2, 1) +
                            5 * AxB(2, 0) * c(1, 1) + 5 * AxB(2, 1) * c(1, 2) +
                            c(1, 0) * c(2, 1) + c(1, 1) * c(2, 0));  // T21s,22s
      block(3, 3) =
          0.25 * fac4 *
          (35 * AxA(0, 0) * BxB(0, 0) - 35 * AxA(0, 0) * BxB(1, 1) -
           35 * AxA(1, 1) * BxB(0, 0) + 35 * AxA(1, 1) * BxB(1, 1) +
           20 * AxB(0, 0) * c(0, 0) - 20 * AxB(0, 1) * c(0, 1) -
           20 * AxB(1, 0) * c(1, 0) + 20 * AxB(0, 0) * c(1, 1) +
           2 * c(0, 0) * c(0, 0) - 2 * c(0, 1) * c(0, 1) -
           2 * c(1, 0) * c(1, 0) + 2 * c(1, 1) * c(1, 1));  // T22c,22c
      block(4, 3) =
          0.5 * fac4 *
          (35 * BxB(0, 1) * AxA(0, 0) - 35 * BxB(1, 2) * AxA(1, 1) +
           10 * AxB(0, 0) * c(0, 1) + 10 * AxB(0, 1) * c(0, 0) -
           10 * AxB(1, 0) * c(1, 1) - 10 * AxB(1, 1) * c(1, 2) +
           2 * c(0, 0) * c(0, 1) - 2 * c(1, 0) * c(1, 1));  // T22c,22s
      block(4, 4) = 0.5 * fac4 *
                    (35 * AxA(0, 1) * BxB(0, 1) + 5 * AxB(0, 0) * c(1, 1) +
                     5 * AxB(0, 1) * c(1, 0) + 5 * AxB(1, 0) * c(0, 1) +
                     5 * AxB(1, 1) * c(0, 0) + c(0, 0) * c(1, 1) +
                     c(0, 1) * c(1, 0));  // T22s,22s

      block.triangularView<Eigen::StrictlyUpper>() =
          block.triangularView<Eigen::StrictlyLower>().transpose();
      interaction.block<5, 5>(4, 4) = block;
    }
  }
  return interaction;  // in units of 4piepsilon0
}

Eigen::Matrix3d eeInteractor::FillTholeInteraction_diponly(
    const PolarSite& site1, const PolarSite& site2) const {

  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();
  const Eigen::Vector3d r_AB =
      posB - posA;               // Vector of the distance between polar sites
  const double R = r_AB.norm();  // Norm of distance vector
  const Eigen::Vector3d pos_a =
      r_AB /
      R;  // unit vector on the sites reciprocal direction; This points toward A
  const double fac2 = std::pow(R, -3);
  const double au3 =
      _expdamping /
      (fac2 * std::sqrt(site1.getEigenDamp() *
                        site2.getEigenDamp()));  // dimensionless eigendamp is
  double lambda3 = fac2;
  double lambda5 = fac2;
  if (au3 < 40) {
    const double exp_ua = std::exp(-au3);
    lambda3 *= 1 - exp_ua;
    lambda5 *= 1 - (1 + au3) * exp_ua;
  }
  Eigen::Matrix3d result = -3 * lambda5 * pos_a * pos_a.transpose();
  result.diagonal().array() += lambda3;
  return result;  // T_1alpha,1beta (alpha,beta=x,y,z)
}

Eigen::Matrix4d eeInteractor::FillTholeInteraction_noQuadrupoles(
    const PolarSite& site1, const PolarSite& site2) const {
  Eigen::Matrix4d interaction = Eigen::Matrix4d::Zero();
  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();
  const Eigen::Vector3d r_AB =
      posB - posA;               // Vector of the distance between polar sites
  const double R = r_AB.norm();  // Norm of distance vector
  const Eigen::Vector3d pos_a =
      r_AB /
      R;  // unit vector on the sites reciprocal direction; This points toward A

  const double fac0 = 1 / R;
  const double fac1 = std::pow(fac0, 2);
  const double fac2 = std::pow(fac0, 3);
  const double au3 =
      _expdamping /
      (fac2 * std::sqrt(site1.getEigenDamp() *
                        site2.getEigenDamp()));  // dimensionless eigendamp is
  double lambda3 = 1.0;

  if (au3 < 40) {
    const double exp_ua = std::exp(-au3);
    lambda3 = 1 - exp_ua;
  }
  // Dipole-Charge Interaction
  interaction.block<3, 1>(1, 0) =
      lambda3 * fac1 * pos_a;  // T_1alpha,00 (alpha=x,y,z)
  // Charge-Dipole Interaction
  interaction.block<1, 3>(0, 1) =
      -lambda3 * fac1 * pos_a;  // T_00,1alpha (alpha=x,y,z)
  interaction.block<3, 3>(1, 1) = FillTholeInteraction_diponly(site1, site2);
  return interaction;
}

Matrix9d eeInteractor::FillTholeInteraction(const PolarSite& site1,
                                            const PolarSite& site2) const {

  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();
  int rankA = site1.getRank();
  int rankB = site2.getRank();

  Matrix9d interaction = Matrix9d::Zero();
  const Eigen::Vector3d r_AB =
      posB - posA;               // Vector of the distance between polar sites
  const double R = r_AB.norm();  // Norm of distance vector
  const Eigen::Vector3d pos_a =
      r_AB /
      R;  // unit vector on the sites reciprocal direction; This points toward A
  const Eigen::Vector3d pos_b = -pos_a;  // unit vector on the sites reciprocal
                                         // direction; This points toward B

  const double sqr3 = std::sqrt(3);
  const double fac0 = 1 / R;
  const double fac2 = std::pow(fac0, 3);
  const double au3 =
      _expdamping /
      (fac2 * std::sqrt(site1.getEigenDamp() *
                        site2.getEigenDamp()));  // dimensionless eigendamp is
  double lambda5 = 1.0;
  double lambda7 = 1.0;
  if (au3 < 40) {
    const double exp_ua = std::exp(-au3);
    lambda5 = 1 - (1 + au3) * exp_ua;
    lambda7 = 1 - (1 + au3 + 0.6 * au3 * au3) * exp_ua;
  }
  interaction.block<4, 4>(0, 0) =
      FillTholeInteraction_noQuadrupoles(site1, site2);
  const double fac3 = std::pow(fac0, 4);
  const Eigen::Matrix3d AxA = pos_a * pos_a.transpose();
  const Eigen::Matrix3d c = Eigen::Matrix3d::Identity();
  Eigen::Matrix<double, 3, 5> block;
  // Quadrupole-Dipole Interaction
  block.col(0) = 0.5 * fac3 *
                 (lambda7 * 15 * AxA(2, 2) * pos_b +
                  lambda5 * (6 * pos_a.z() * c.col(2) -
                             3 * pos_b));  // T20-1beta (beta=x,y,z)
  block.col(1) = fac3 * sqr3 *
                 (lambda5 * (pos_a.x() * c.col(2) + c.col(0) * pos_a.z()) +
                  lambda7 * 5 * AxA(0, 2) * pos_b);  // T21c-1beta (beta=x,y,z)
  block.col(2) = fac3 * sqr3 *
                 (lambda5 * (pos_a.y() * c.col(2) + c.col(1) * pos_a.z()) +
                  lambda7 * 5 * AxA(1, 2) * pos_b);  // T21s-1beta (beta=x,y,z)
  block.col(3) =
      fac3 * 0.5 * sqr3 *
      (lambda7 * 5 * (AxA(0, 0) - AxA(1, 1)) * pos_b +
       lambda5 * (2 * pos_a.x() * c.col(0) -
                  2 * pos_a.y() * c.col(1)));  // T22c-1beta (beta=x,y,z)
  block.col(4) = fac3 * sqr3 *
                 (lambda7 * 5 * AxA(0, 1) * pos_b +
                  lambda5 * (pos_a.x() * c.col(1) +
                             pos_a.y() * c.col(0)));  // T22s-1beta (beta=x,y,z)

  if (rankA > 1) {
    interaction.block<5, 3>(4, 1) = block.transpose();
  }
  if (rankB > 1) {
    interaction.block<3, 5>(1, 4) = -block;
  }
  return interaction;
}

double eeInteractor::InteractStatic_site(StaticSite& site1,
                                         StaticSite& site2) const {

  if (site1.getRank() < 2 && site2.getRank() < 2) {
    const Eigen::Matrix4d Tab = FillInteraction_noQuadrupoles(site1, site2);
    auto V1 = Tab * (site2.Q().segment<4>(0));
    const Eigen::Vector4d V2 = Tab.transpose() * (site1.Q().segment<4>(0));
    site1.V().segment<4>(0) += V1;
    site2.V().segment<4>(0) += V2;
    return (site2.Q().segment<4>(0)).dot(V2);
  } else {
    const Matrix9d Tab = FillInteraction(site1, site2);
    auto V1 = Tab * site2.Q();
    const Vector9d V2 = Tab.transpose() * site1.Q();
    site1.V() += V1;
    site2.V() += V2;
    return site2.Q().dot(V2);
  }
}

double eeInteractor::InteractStatic_site(const StaticSite& site1,
                                         StaticSite& site2) const {
  if (site1.getRank() < 2 && site2.getRank() < 2) {
    const Eigen::Matrix4d Tab = FillInteraction_noQuadrupoles(site1, site2);
    const Eigen::Vector4d V2 = Tab.transpose() * (site1.Q().segment<4>(0));
    site2.V().segment<4>(0) += V2;
    return (site2.Q().segment<4>(0)).dot(V2);
  } else {
    const Matrix9d Tab = FillInteraction(site1, site2);
    const Vector9d V2 = Tab.transpose() * site1.Q();
    site2.V() += V2;
    return site2.Q().dot(V2);
  }
}

double eeInteractor::InteractPolar_site(const PolarSite& site1,
                                        PolarSite& site2) const {

  const Matrix9d tTab =
      FillTholeInteraction(site1, site2);  // \tilde{T}^(ab)_tu
  const Vector9d V2 =
      tTab.block<3, 9>(1, 0).transpose() * site1.Induced_Dipole();
  auto V1 = tTab.block<9, 3>(0, 1) * site2.Induced_Dipole();
  site2.V() += V2;
  // indu -indu
  double e = site1.Induced_Dipole().transpose() * tTab.block<3, 3>(1, 1) *
             site2.Induced_Dipole();
  // indu-stat
  e += site2.Q().dot(V2);
  // stat-indu
  e += site1.Q().dot(V1);
  return e;
}

double eeInteractor::InteractPolar_site(PolarSite& site1,
                                        PolarSite& site2) const {
  const Matrix9d tTab =
      FillTholeInteraction(site1, site2);  // \tilde{T}^(ab)_tu
  const Vector9d V2 =
      site1.Induced_Dipole().transpose() * tTab.block<3, 9>(1, 0);
  const Vector9d V1 = tTab.block<9, 3>(0, 1) * site2.Induced_Dipole();

  site1.V() += V1;
  site2.V() += V2;
  double e = site1.Induced_Dipole().transpose() * tTab.block<3, 3>(1, 1) *
             site2.Induced_Dipole();
  // indu-stat
  e += site2.Q().dot(V2);
  // stat-indu
  e += site1.Q().dot(V1);
  return e;
}

template <class T1, class T2>
double eeInteractor::InteractStatic(T1& seg1, T2& seg2) const {
  assert(&seg1 != reinterpret_cast<T1*>(&seg2)  &&
         "InteractStatic(seg1,seg2) needs two distinct objects");
  double energy = 0.0;
  for (auto& site1 : seg1) {
    for (auto& site2 : seg2) {
      energy += InteractStatic_site(site1, site2);
    }
  }
  return energy;
}

template double eeInteractor::InteractStatic(const StaticSegment& seg1,
                                             StaticSegment& seg2) const;

template double eeInteractor::InteractStatic(StaticSegment& seg1,
                                             StaticSegment& seg2) const;
template double eeInteractor::InteractStatic(const StaticSegment& seg1,
                                             PolarSegment& seg2) const;
template double eeInteractor::InteractStatic(const PolarSegment& seg1,
                                             PolarSegment& seg2) const;
template double eeInteractor::InteractStatic(PolarSegment& seg1,
                                             PolarSegment& seg2) const;

template <class T>
double eeInteractor::InteractStatic_IntraSegment(T& seg) const {
  double energy = 0.0;
  for (int i = 1; i < seg.size(); i++) {
    for (int j = 0; j < i; j++) {
      energy += InteractStatic_site(seg[i], seg[j]);
    }
  }
  return energy;
}

template double eeInteractor::InteractStatic_IntraSegment(
    StaticSegment& seg) const;
template double eeInteractor::InteractStatic_IntraSegment(
    PolarSegment& seg) const;

double eeInteractor::InteractPolar_IntraSegment(PolarSegment& seg) const {
  double energy = 0.0;
  for (int i = 1; i < seg.size(); i++) {
    for (int j = 0; j < i; j++) {
      energy += InteractPolar_site(seg[i], seg[j]);
    }
  }
  return energy;
}

double eeInteractor::InteractPolar(PolarSegment& seg1,
                                   PolarSegment& seg2) const {
  assert(&seg1 != &seg2 &&
         "InteractPolar(seg1,seg2) needs two distinct objects");
  double energy = 0.0;
  for (PolarSite& site1 : seg1) {
    for (PolarSite& site2 : seg2) {
      energy += InteractPolar_site(site1, site2);
    }
  }
  return energy;
}

double eeInteractor::InteractPolar(const PolarSegment& seg1,
                                   PolarSegment& seg2) const {
  assert(&seg1 != &seg2 &&
         "InteractPolar(seg1,seg2) needs two distinct objects");
  double energy = 0.0;
  for (const PolarSite& site1 : seg1) {
    for (PolarSite& site2 : seg2) {
      energy += InteractPolar_site(site1, site2);
    }
  }
  return energy;
}

Eigen::VectorXd eeInteractor::Cholesky_IntraSegment(
    const PolarSegment& seg) const {
  int size = 3 * seg.size();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);
  for (int i = 1; i < seg.size(); i++) {
    for (int j = 0; j < i; j++) {
      A.block<3, 3>(3 * i, 3 * j) =
          FillTholeInteraction_diponly(seg[i], seg[j]);
    }
  }

  for (int i = 0; i < seg.size(); i++) {
    A.block<3, 3>(3 * i, 3 * i) = seg[i].getPInv();
  }
  Eigen::VectorXd b = Eigen::VectorXd(size);
  for (int i = 0; i < seg.size(); i++) {
    b.segment<3>(3 * i) = -seg[i].V().segment<3>(1);
  }

  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> lltOfA(A);
  return lltOfA.solve(b);
}

}  // namespace xtp
}  // namespace votca
