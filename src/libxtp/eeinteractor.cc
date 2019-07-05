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

  Eigen::Matrix<double, N, M> interaction = Eigen::Matrix<double, N, M>::Zero();
  const Eigen::Vector3d r_AB =
      posB - posA;               // Vector of the distance between polar sites
  const double R = r_AB.norm();  // Norm of distance vector
  const Eigen::Vector3d pos_a = r_AB / R;  // unit vector pointing from A to B

  const double fac1 = 1.0 / R;
  interaction(0, 0) = fac1;
  if (N > 1 || M > 1) {
    const double fac2 = std::pow(fac1, 2);

    // Dipole-Charge Interaction
    if (N > 1) {
      interaction.block(1, 0, 3, 1) =
          fac2 * pos_a;  // T_1alpha,00 (alpha=x,y,z)
    }
    // Charge-Dipole Interaction
    if (M > 1) {
      interaction.block(0, 1, 1, 3) =
          -fac2 * pos_a.transpose();  // T_00,1alpha (alpha=x,y,z)
    }

    const double fac3 = std::pow(fac1, 3);
    if (N > 1 && M > 1) {
      // Dipole-Dipole Interaction
      Eigen::Matrix3d dd = -3 * fac3 * pos_a * pos_a.transpose();
      dd.diagonal().array() += fac3;
      interaction.block(1, 1, 3, 3) = dd;
      // T_1alpha,1beta // (alpha,beta=x,y,z)
    }

    const double sqr3 = std::sqrt(3);
    if (N > 4 || M > 4) {
      const AxA aa(pos_a);
      // Quadrupole-Charge interaction
      Eigen::Matrix<double, 1, 5> Qq;
      Qq(0) = fac3 * 0.5 * (3 * aa.zz() - 1);           // T20,00
      Qq(1) = fac3 * sqr3 * aa.xz();                    // T21c,00
      Qq(2) = fac3 * sqr3 * aa.yz();                    // T21s,000
      Qq(3) = fac3 * 0.5 * sqr3 * (aa.xx() - aa.yy());  // T22c,00
      Qq(4) = fac3 * sqr3 * aa.xy();                    // T22s,00
      if (N > 4) {
        interaction.block(4, 0, 5, 1) = Qq.transpose();
      }
      if (M > 4) {
        interaction.block(0, 4, 1, 5) = Qq;
      }

      if (N > 1 && M > 1) {

        const double fac4 = std::pow(fac1, 4);

        Eigen::Matrix<double, 3, 5> dQ;
        // Quadrupole-Dipole Interaction
        dQ.col(0) = -0.5 * fac4 * (15 * aa.zz() - 3) * pos_a;
        dQ.col(0).z() += 3 * fac4 * pos_a.z();  // T20-1beta (beta=x,y,z)

        double faccol = fac4 * sqr3;
        dQ.col(1) = -faccol * 5 * aa.xz() * pos_a;
        dQ.col(1).z() += faccol * pos_a.x();
        dQ.col(1).x() += faccol * pos_a.z();  // T21c-1beta (beta=x,y,z)

        dQ.col(2) = -faccol * 5 * aa.yz() * pos_a;
        dQ.col(2).z() += faccol * pos_a.y();
        dQ.col(2).y() += faccol * pos_a.z();

        dQ.col(3) = -faccol * 2.5 * (aa.xx() - aa.yy()) * pos_a;
        dQ.col(3).x() += faccol * pos_a.x();
        dQ.col(3).y() -= faccol * pos_a.y();  // T22c-1beta (beta=x,y,z)

        dQ.col(4) = -faccol * 5 * aa.xy() * pos_a;
        dQ.col(4).y() += faccol * pos_a.x();
        dQ.col(4).x() += faccol * pos_a.y();  // T22s-1beta (beta=x,y,z)

        if (N > 4) {
          interaction.block(4, 1, 5, 3) = dQ.transpose();
        }
        if (M > 4) {
          interaction.block(1, 4, 3, 5) = -dQ;
        }
      }

      if (N > 4 && M > 4) {
        const double fac5 = std::pow(fac1, 5);
        // Quadrupole-Quadrupole Interaction
        Eigen::Matrix<double, 5, 5> QQ;
        QQ(0, 0) =
            fac5 * (3. / 4.) * ((35 * aa.zz() - 30) * aa.zz() + 3);  // T20,20
        QQ(1, 0) =
            0.5 * fac5 * sqr3 * aa.xz() * (35 * aa.zz() - 15);  // T20,21c
        QQ(2, 0) =
            0.5 * fac5 * sqr3 * aa.yz() * (35 * aa.zz() - 15);  // T20,21s
        QQ(3, 0) = 0.25 * fac5 * sqr3 * 5 *
                   (7 * (aa.yz() - aa.xz()) - aa.xx() + aa.yy());  // T20,22c
        QQ(4, 0) =
            0.5 * fac5 * sqr3 * 5 * aa.xy() * (7 * aa.zz() - 1);  // T20,22s
        QQ(1, 1) = fac5 * (35 * aa.xz() * aa.xz() - 5 * (aa.xx() + aa.zz()) +
                           1);                                    // T21c,21c
        QQ(2, 1) = fac5 * 5 * (7 * aa.xz() * aa.yz() - aa.xy());  // T21c,21s
        QQ(3, 1) =
            0.5 * fac5 *
            (35 * aa.xz() * (aa.xx() - aa.yy()) - 10 * aa.xx());  // T21c,22c
        QQ(4, 1) = fac5 * 7 * (5 * aa.xz() * aa.xy() - aa.yz());  // T21c,22s
        QQ(2, 2) = fac5 * (35 * aa.yz() * aa.yz() - 5 * (aa.yy() + aa.zz()) +
                           1);  // T21s,21s
        QQ(3, 2) =
            0.5 * fac5 * 35 * aa.yz() * (aa.xx() - aa.yy() + 10);  // T21s,22c
        QQ(4, 2) = fac5 * 5 * (7 * aa.yz() * aa.xy() - aa.xz());   // T21s,22s
        QQ(3, 3) = 0.25 * fac5 *
                   (35 * std::pow(aa.xx() - aa.yy(), 2) - 40 * aa.xx() +
                    4);  // T22c,22c
        QQ(4, 3) = 0.5 * fac5 * 35 *
                   (aa.xy() * aa.xx() - aa.yz() * aa.yy());  // T22c,22s
        QQ(4, 4) =
            0.5 * fac5 *
            (35 * aa.xy() * aa.xy() - 5 * (aa.xx() + aa.yy()) + 1);  // T22s,22s

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
