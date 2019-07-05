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

template <int N, int M>
Eigen::Matrix<double, N, M> eeInteractor::FillInteraction(
    const StaticSite& site1, const StaticSite& site2) const {

  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();

  Eigen::Vector3d a =
      posB - posA;  // Vector of the distance between polar sites
  const double R = a.norm();
  const double fac1 = 1.0 / R;
  a *= fac1;  // unit vector pointing from A to B
  Eigen::Matrix<double, N, M> interaction;
  interaction(0, 0) = fac1;
  if (N > 1 || M > 1) {
    const double fac2 = std::pow(fac1, 2);

    // Dipole-Charge Interaction
    if (N > 1) {
      interaction.block(1, 0, 3, 1) = fac2 * a;  // T_1alpha,00 (alpha=x,y,z)
    }
    // Charge-Dipole Interaction
    if (M > 1) {
      interaction.block(0, 1, 1, 3) =
          -fac2 * a.transpose();  // T_00,1alpha (alpha=x,y,z)
    }

    const double fac3 = std::pow(fac1, 3);
    if (N > 1 && M > 1) {
      // Dipole-Dipole Interaction
      Eigen::Matrix3d dd = -3 * fac3 * a * a.transpose();
      dd.diagonal().array() += fac3;
      interaction.block(1, 1, 3, 3) = dd;
      // T_1alpha,1beta // (alpha,beta=x,y,z)
    }

    const double sqr3 = std::sqrt(3);
    if (N > 4 || M > 4) {
      const AxA r(a);
      {
        // Quadrupole-Charge interaction
        Eigen::Matrix<double, 1, 5> Qq;
        Qq(0) = fac3 * (1.5 * r.zz() - 0.5);  // T20,00
        Qq(1) = r.xz();                       // T21c,00
        Qq(2) = r.yz();                       // T21s,000
        Qq(3) = 0.5 * (r.xx() - r.yy());      // T22c,00
        Qq(4) = r.xy();                       // T22s,00
        Qq.tail<4>() *= (fac3 * sqr3);
        if (N > 4) {
          interaction.block(4, 0, 5, 1) = Qq.transpose();
        }
        if (M > 4) {
          interaction.block(0, 4, 1, 5) = Qq;
        }
      }
      if (N > 1 && M > 1) {

        Eigen::Matrix<double, 3, 5> dQ;
        const double fac4 = std::pow(fac1, 4);
        // Quadrupole-Dipole Interaction
        const Eigen::Vector3d afac = a * fac4;

        dQ.col(0) = (1.5 - 7.5 * r.zz()) * afac;
        dQ.col(0).z() += 3 * afac.z();  // T20-1beta (beta=x,y,z)

        dQ.col(1) = -5 * r.xz() * afac;
        dQ.col(1).z() += afac.x();
        dQ.col(1).x() += afac.z();  // T21c-1beta (beta=x,y,z)

        dQ.col(2) = -5 * r.yz() * afac;
        dQ.col(2).z() += afac.y();
        dQ.col(2).y() += afac.z();

        dQ.col(3) = -2.5 * (r.xx() - r.yy()) * afac;
        dQ.col(3).x() += afac.x();
        dQ.col(3).y() -= afac.y();  // T22c-1beta (beta=x,y,z)

        dQ.col(4) = -5 * r.xy() * a;
        dQ.col(4).y() += afac.x();
        dQ.col(4).x() += afac.y();  // T22s-1beta (beta=x,y,z)

        dQ.rightCols<4>() *= sqr3;

        if (N > 4) {
          interaction.block(4, 1, 5, 3) = dQ.transpose();
        }
        if (M > 4) {
          interaction.block(1, 4, 3, 5) = -dQ;
        }
      }

      if (N > 4 && M > 4) {

        // Quadrupole-Quadrupole Interaction
        Eigen::Matrix<double, 5, 5> QQ;
        QQ(0, 0) = 0.75 * ((35 * r.zz() - 30) * r.zz() + 3);  // T20,20
        double temp = 0.5 * sqr3 * (35 * r.zz() - 15);
        QQ(1, 0) = temp * r.xz();  // T20,21c
        QQ(2, 0) = temp * r.yz();  // T20,21s

        temp = 5 * (7 * r.zz() - 1);
        QQ(3, 0) = sqr3 * 0.25 * temp * (r.xx() - r.yy());  // T20,22c
        QQ(4, 0) = sqr3 * 0.5 * temp * r.xy();              // T20,22s
        QQ(1, 1) =
            35 * r.zz() * r.xx() - 5 * (r.xx() + r.zz()) + 1;     // T21c,21c
        QQ(2, 1) = r.xy() * temp;                                 // T21c,21s
        QQ(3, 1) = 0.5 * r.xz() * (35 * (r.xx() - r.yy()) - 10);  // T21c,22c
        QQ(4, 1) = r.yz() * temp;                                 // T21c,22s
        QQ(2, 2) =
            5 * (7 * r.yy() * r.zz() - (r.yy() + r.zz())) + 1;    // T21s,21s
        QQ(3, 2) = 0.5 * r.yz() * (35 * (r.xx() - r.yy()) + 10);  // T21s,22c
        QQ(4, 2) = r.xz() * temp;                                 // T21s,22s
        QQ(3, 3) = 8.75 * std::pow(r.xx() - r.yy(), 2) - 5 * (r.xx() + r.yy()) +
                   1;                                  // T22c,22c
        QQ(4, 3) = 17.5 * r.xy() * (r.xx() - r.yy());  // T22c,22s
        QQ(4, 4) =
            5 * (7 * r.xx() * r.yy() - (r.xx() + r.yy())) + 1;  // T22s,22s
        const double fac5 = std::pow(fac1, 5);
        QQ.triangularView<Eigen::Lower>() *= fac5;
        QQ.triangularView<Eigen::StrictlyUpper>() =
            QQ.triangularView<Eigen::StrictlyLower>().transpose();

        interaction.block(4, 4, 5, 5) = QQ;
      }
    }
  }
  return interaction;  // in units of 4piepsilon0
}

Eigen::Matrix3d eeInteractor::FillTholeInteraction(
    const PolarSite& site1, const PolarSite& site2) const {

  const Eigen::Vector3d& posB = site2.getPos();
  const Eigen::Vector3d& posA = site1.getPos();
  Eigen::Vector3d a =
      posB - posA;            // Vector of the distance between polar sites
  const double R = a.norm();  // Norm of distance vector
  const double fac1 = 1 / R;
  a *= fac1;  // unit vector pointing from A to B

  double lambda3 = std::pow(fac1, 3);
  double lambda5 = lambda3;
  const double au3 =
      _expdamping /
      (lambda3 * std::sqrt(site1.getEigenDamp() *
                           site2.getEigenDamp()));  // au3 is dimensionless
  if (au3 < 40) {
    const double exp_ua = std::exp(-au3);
    lambda3 *= (1 - exp_ua);
    lambda5 *= (1 - (1 + au3) * exp_ua);
  }
  Eigen::Matrix3d result = -3 * lambda5 * a * a.transpose();
  result.diagonal().array() += lambda3;
  return result;  // T_1alpha,1beta (alpha,beta=x,y,z)
}

template <enum Estatic CE>
double eeInteractor::ApplyStaticField_site(const StaticSite& site1,
                                           PolarSite& site2) const {
  Eigen::Vector3d V = Eigen::Vector3d::Zero();
  double e = 0.0;
  if (site1.getRank() == 0) {
    if (site2.getRank() < 2) {
      const Eigen::RowVector4d Tab0 = FillInteraction<1, 4>(site1, site2);
      V = site1.Q()[0] * Tab0.tail<3>();
      e = site1.Q()[0] * Tab0.dot(site2.Q().head<4>());
    } else {
      const Eigen::Matrix<double, 1, 9> Tab0 =
          FillInteraction<1, 9>(site1, site2);
      V = site1.Q()[0] * Tab0.segment<3>(1);
      e = site1.Q()[0] * Tab0.dot(site2.Q());
    }
  } else if (site1.getRank() == 1) {
    if (site2.getRank() < 2) {
      const Eigen::Matrix4d Tab1 = FillInteraction<4, 4>(site1, site2);
      V = Tab1.block<4, 3>(0, 1).transpose() * (site1.Q().head<4>());
      e = site1.Q().head<4>().transpose() * Tab1 * site2.Q().head<4>();
    } else {
      const Eigen::Matrix<double, 4, 9> Tab1 =
          FillInteraction<4, 9>(site1, site2);
      V = Tab1.block<4, 3>(0, 1).transpose() * (site1.Q().head<4>());
      e = site1.Q().head<4>().transpose() * Tab1 * site2.Q();
    }
  } else if (site1.getRank() == 2) {
    if (site2.getRank() < 2) {
      const Eigen::Matrix<double, 9, 4> Tab2 =
          FillInteraction<9, 4>(site1, site2);
      V = Tab2.block<9, 3>(0, 1).transpose() * (site1.Q());
      e = site1.Q().transpose() * Tab2 * site2.Q().head<4>();
    } else {
      const Eigen::Matrix<double, 9, 9> Tab2 =
          FillInteraction<9, 9>(site1, site2);
      V = Tab2.block<9, 3>(0, 1).transpose() * (site1.Q());
      e = site1.Q().transpose() * Tab2 * site2.Q();
    }
  } else {
    throw std::runtime_error("Ranks higher 2 not implemented");
  }

  if (CE == Estatic::noE_V) {
    site2.V_noE() += V;
  } else {
    site2.V() += V;
  }
  return e;
}

template <enum Estatic CE>
double eeInteractor::ApplyInducedField_site(const PolarSite& site1,
                                            PolarSite& site2) const {
  const Eigen::Matrix3d tTab = FillTholeInteraction(site1, site2);
  const Eigen::Vector3d V = tTab.transpose() * site1.Induced_Dipole();
  double e = 0.0;
  if (CE == Estatic::noE_V) {
    site2.V_noE() += V;
  } else {
    site2.V() += V;
    e = CalcPolar_stat_Energy_site(site1, site2);
  }
  return e;
}

double eeInteractor::CalcStaticEnergy_site(const StaticSite& site1,
                                           const StaticSite& site2) const {
  // you can optimise these for all possible configs;
  if (site1.getRank() < 1 && site2.getRank() < 1) {
    return site1.getCharge() * FillInteraction<1, 1>(site1, site2)(0, 0) *
           site2.getCharge();
  } else if (site1.getRank() == 0 && site2.getRank() == 2) {
    const Eigen::Matrix<double, 1, 9> Tab = FillInteraction<1, 9>(site1, site2);
    return site1.Q()(0) * Tab.dot(site2.Q());
  } else if (site1.getRank() == 2 && site2.getRank() == 0) {
    const Eigen::Matrix<double, 9, 1> Tab = FillInteraction<9, 1>(site1, site2);
    return site1.Q().dot(Tab) * site2.Q()(0);
  } else if (site1.getRank() < 2 && site2.getRank() < 2) {
    const Eigen::Matrix<double, 4, 4> Tab = FillInteraction<4, 4>(site1, site2);
    return site1.Q().head<4>().transpose() * Tab * site2.Q().head<4>();
  } else {
    const Eigen::Matrix<double, 9, 9> Tab = FillInteraction<9, 9>(site1, site2);
    return site1.Q().transpose() * Tab * site2.Q();
  }
}

double eeInteractor::CalcPolar_stat_Energy_site(const PolarSite& site1,
                                                const StaticSite& site2) const {
  if (site2.getRank() == 0) {
    const Eigen::Vector4d Tab0 = FillInteraction<4, 1>(site1, site2);
    return site1.Induced_Dipole().dot(Tab0.tail<3>()) * site2.Q()[0];
  } else if (site2.getRank() == 1) {
    const Eigen::Matrix4d Tab1 = FillInteraction<4, 4>(site1, site2);
    return site1.Induced_Dipole().transpose() * Tab1.block<3, 4>(1, 0) *
           (site2.Q().head<4>());
  } else if (site2.getRank() == 2) {
    const Eigen::Matrix<double, 4, 9> Tab2 =
        FillInteraction<4, 9>(site1, site2);
    return site1.Induced_Dipole().transpose() * Tab2.block<3, 9>(1, 0) *
           (site2.Q());
  } else {
    throw std::runtime_error("Ranks higher 2 not implemented");
  }
}

eeInteractor::E_terms eeInteractor::CalcPolarEnergy_site(
    const PolarSite& site1, const StaticSite& site2) const {
  eeInteractor::E_terms val;
  val.E_indu_stat() = CalcPolar_stat_Energy_site(site1, site2);
  return val;
}

eeInteractor::E_terms eeInteractor::CalcPolarEnergy_site(
    const PolarSite& site1, const PolarSite& site2) const {
  // contributions are stat-induced, induced-stat and induced induced
  eeInteractor::E_terms val;
  val.E_indu_indu() = site1.Induced_Dipole().transpose() *
                      FillTholeInteraction(site1, site2) *
                      site2.Induced_Dipole();
  val.E_indu_stat() = CalcPolar_stat_Energy_site(site1, site2);
  val.E_indu_stat() += CalcPolar_stat_Energy_site(site2, site1);
  return val;
}

template <class T, enum Estatic CE>
double eeInteractor::ApplyStaticField(const T& segment1,
                                      PolarSegment& segment2) const {
  double e = 0.0;
  for (PolarSite& s2 : segment2) {
    for (const auto& s1 : segment1) {
      e += ApplyStaticField_site<CE>(s1, s2);
    }
  }
  return e;
}

template double eeInteractor::ApplyStaticField<StaticSegment, Estatic::V>(
    const StaticSegment& seg1, PolarSegment& seg2) const;
template double eeInteractor::ApplyStaticField<StaticSegment, Estatic::noE_V>(
    const StaticSegment& seg1, PolarSegment& seg2) const;
template double eeInteractor::ApplyStaticField<PolarSegment, Estatic::V>(
    const PolarSegment& seg1, PolarSegment& seg2) const;
template double eeInteractor::ApplyStaticField<PolarSegment, Estatic::noE_V>(
    const PolarSegment& seg1, PolarSegment& seg2) const;

template <enum Estatic CE>
double eeInteractor::ApplyInducedField(const PolarSegment& segment1,
                                       PolarSegment& segment2) const {
  double e = 0;
  for (PolarSite& s2 : segment2) {
    for (const PolarSite& s1 : segment1) {
      e += ApplyInducedField_site<CE>(s1, s2);
    }
  }
  return e;
}

template double eeInteractor::ApplyInducedField<Estatic::V>(
    const PolarSegment& seg1, PolarSegment& seg2) const;
template double eeInteractor::ApplyInducedField<Estatic::noE_V>(
    const PolarSegment& seg1, PolarSegment& seg2) const;

template <class S1, class S2>
double eeInteractor::CalcStaticEnergy(const S1& segment1,
                                      const S2& segment2) const {
  double e = 0;
  for (const auto& s1 : segment2) {
    for (const auto& s2 : segment1) {
      e += CalcStaticEnergy_site(s2, s1);
    }
  }
  return e;
}

template double eeInteractor::CalcStaticEnergy(const StaticSegment& seg1,
                                               const PolarSegment& seg2) const;
template double eeInteractor::CalcStaticEnergy(const StaticSegment& seg1,
                                               const StaticSegment& seg2) const;
template double eeInteractor::CalcStaticEnergy(const PolarSegment& seg1,
                                               const PolarSegment& seg2) const;
template double eeInteractor::CalcStaticEnergy(const PolarSegment& seg1,
                                               const StaticSegment& seg2) const;
template <class S>
double eeInteractor::CalcStaticEnergy_IntraSegment(const S& seg) const {
  double e = 0;
  for (int i = 0; i < seg.size(); i++) {
    for (int j = 0; j < i; j++) {
      e += CalcStaticEnergy_site(seg[i], seg[j]);
    }
  }
  return e;
}
template double eeInteractor::CalcStaticEnergy_IntraSegment(
    const PolarSegment& seg1) const;
template double eeInteractor::CalcStaticEnergy_IntraSegment(
    const StaticSegment& seg2) const;

template <class S1, class S2>
eeInteractor::E_terms eeInteractor::CalcPolarEnergy(const S1& segment1,
                                                    const S2& segment2) const {
  eeInteractor::E_terms e;
  for (const auto& s1 : segment2) {
    for (const auto& s2 : segment1) {
      e += CalcPolarEnergy_site(s2, s1);
    }
  }
  return e;
}

template eeInteractor::E_terms eeInteractor::CalcPolarEnergy(
    const PolarSegment& seg1, const PolarSegment& seg2) const;
template eeInteractor::E_terms eeInteractor::CalcPolarEnergy(
    const PolarSegment& seg1, const StaticSegment& seg2) const;

double eeInteractor::CalcPolarEnergy_IntraSegment(
    const PolarSegment& seg) const {
  double e = 0;
  for (int i = 0; i < seg.size(); i++) {
    for (int j = 0; j < i; j++) {
      e += seg[i].Induced_Dipole().transpose() *
           FillTholeInteraction(seg[i], seg[j]) * seg[j].Induced_Dipole();
    }
  }
  return e;
}

Eigen::VectorXd eeInteractor::Cholesky_IntraSegment(
    const PolarSegment& seg) const {
  int size = 3 * seg.size();

  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(size, size);
  for (int i = 1; i < seg.size(); i++) {
    for (int j = 0; j < i; j++) {
      A.block<3, 3>(3 * i, 3 * j) = FillTholeInteraction(seg[i], seg[j]);
    }
  }

  for (int i = 0; i < seg.size(); i++) {
    A.block<3, 3>(3 * i, 3 * i) = seg[i].getPInv();
  }
  Eigen::VectorXd b = Eigen::VectorXd(size);
  for (int i = 0; i < seg.size(); i++) {
    const Eigen::Vector3d V = seg[i].V() + seg[i].V_noE();
    b.segment<3>(3 * i) = -V;
  }

  Eigen::LLT<Eigen::Ref<Eigen::MatrixXd>> lltOfA(A);
  return lltOfA.solve(b);
}

}  // namespace xtp
}  // namespace votca
