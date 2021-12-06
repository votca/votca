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
#ifndef VOTCA_XTP_INTERACTIONTENSORS_H
#define VOTCA_XTP_INTERACTIONTENSORS_H

#include <optional>
#include <votca/tools/eigen.h>
#include <votca/tools/types.h>

#include <boost/math/constants/constants.hpp>

namespace votca {
namespace xtp {

enum class Screening { none, erfc, erf };

template <Screening screen, Index maxRank>
class InteractionTensor {
 public:
  InteractionTensor(double alpha = 1.0, double thole = 1.0)
      : rank3_(3, 3, 3), rank4_(3, 3, 3, 3), thole_(thole) {
    static_assert(maxRank >= 0 && maxRank < 5);
    a1 = alpha;
    a2 = a1 * a1;
    a3 = a1 * a2;
    if constexpr (maxRank > 2) {
      a5 = a3 * a2;
      a7 = a5 * a2;
    }

    thole2_ = thole_ * thole_;
    thole3_ = thole_ * thole2_;
  }

  void applyTholeDamping(const Eigen::Matrix3d& pol1,
                         const Eigen::Matrix3d& pol2);

  /**
   * \brief this function computes the interaction tensors.
   *
   * NB If you need thole damping, apply it before calculating the tensors.
   */
  void computeTensors(const Eigen::Vector3d& dr) {
    dr_ = dr;
    rank1_ = Eigen::Vector3d::Zero();
    rank2_ = Eigen::Matrix3d::Zero();
    if constexpr (maxRank > 2) {
      rank3_.setZero();
      rank4_.setZero();
    }
    // Perform the actual tensor calculations
    if constexpr (screen == Screening::none) {
      noneTensors();
    } else if constexpr (screen == Screening::erfc) {
      erfcTensors();
    } else if constexpr (screen == Screening::erf) {
      erfTensors();
    }
  }

  double rank0() { return rank0_; }

  const Eigen::Vector3d& rank1() {
    static_assert(maxRank >= 1);
    return rank1_;
  }

  const Eigen::Matrix3d& rank2() {
    static_assert(maxRank >= 2);
    return rank2_;
  }

  const Eigen::Tensor<double, 3>& rank3() {
    static_assert(maxRank >= 3);
    return rank3_;
  }

  const Eigen::Tensor<double, 4>& rank4() {
    static_assert(maxRank >= 4);
    return rank4_;
  }

  Eigen::Vector3d interact(double charge) const;
  Eigen::Vector3d interact(const Eigen::Vector3d& dipole) const;
  Eigen::Vector3d interact(const Eigen::Matrix3d& quadrupole) const;

  double interact(double charge1, double charge2) const;
  double interact(double charge, const Eigen::Vector3d& dipole) const;
  double interact(const Eigen::Vector3d& dipole1,
                  const Eigen::Vector3d& dipole2) const;
  double interact(double charge, const Eigen::Matrix3d& quadrupole) const;
  double interact(const Eigen::Vector3d& dipole,
                  const Eigen::Matrix3d& quadrupole) const;
  double interact(const Eigen::Matrix3d& quadrupole1,
                  const Eigen::Matrix3d& quadrupole2) const;

 private:
  Eigen::Vector3d dr_ = Eigen::Vector3d::Zero();
  double rank0_ = 0.0;
  Eigen::Vector3d rank1_ = Eigen::Vector3d::Zero();
  Eigen::Matrix3d rank2_ = Eigen::Matrix3d::Zero();
  Eigen::Tensor<double, 3> rank3_;
  Eigen::Tensor<double, 4> rank4_;
  double l3 = 1.0;
  double l5 = 1.0;
  double l7 = 1.0;
  double l9 = 1.0;
  double a1 = 1.0;
  double a2 = 1.0;
  double a3 = 1.0;
  double a5 = 1.0;
  double a7 = 1.0;
  void noneTensors();
  void erfcTensors();
  void erfTensors();

  double thole_ = 0.0;
  double thole2_ = 0.0;
  double thole3_ = 0.0;

  double pi = boost::math::constants::pi<double>();
  double rSqrtPi = 1.0 / std::sqrt(pi);
};

template <Screening screen, Index maxRank>
void InteractionTensor<screen, maxRank>::applyTholeDamping(
    const Eigen::Matrix3d& pol1, const Eigen::Matrix3d& pol2) {
  double R = dr_.norm();
  double R2 = R * R;
  double thole_u3 =
      (R * R2) / std::sqrt((1.0 / 3.0) * (pol1.array() * pol2.array()).sum());

  if (thole_ * thole_u3 < 40) {
    double thole_exp = std::exp(-thole_ * thole_u3);
    double thole_u6 = thole_u3 * thole_u3;
    l3 = 1 - thole_exp;
    l5 = 1 - (1 + thole_ * thole_u3) * thole_exp;
    l7 = 1 - (1 + thole_ * thole_u3 + (3. / 5.) * thole2_ * thole_u6);
    l9 = 1 - (1 + thole_ * thole_u3 + (18. / 35.) * thole2_ * thole_u6 +
              (9. / 35.) * thole3_ * thole_u6 * thole_u3) *
                 thole_exp;
  }  // else: the exponent is close to zero and hence the expression is close to
     // 1 and we use the default 1.0 values.
}

template <Screening screen, Index maxRank>
Eigen::Vector3d InteractionTensor<screen, maxRank>::interact(
    double charge) const {
  static_assert(maxRank > 0);
  return charge * rank1_;
}

template <Screening screen, Index maxRank>
Eigen::Vector3d InteractionTensor<screen, maxRank>::interact(
    const Eigen::Vector3d& dipole) const {
  static_assert(maxRank > 1);
  return rank2_ * dipole;
}

template <Screening screen, Index maxRank>
Eigen::Vector3d InteractionTensor<screen, maxRank>::interact(
    const Eigen::Matrix3d& quadrupole) const {
  static_assert(maxRank > 2);
  Eigen::Vector3d result = Eigen::Vector3d::Zero();
  for (Index i = 0; i < 3; ++i) {
    for (Index j = 0; j < 3; ++j) {
      for (Index k = 0; k < 3; ++k) {
        result[k] += quadrupole(j, i) * rank3_(k, j, i);
      }
    }
  }
  return result;
}

template <Screening screen, Index maxRank>
double InteractionTensor<screen, maxRank>::interact(double charge1,
                                                    double charge2) const {
  return charge1 * rank0_ * charge2;
}

template <Screening screen, Index maxRank>
double InteractionTensor<screen, maxRank>::interact(
    double charge, const Eigen::Vector3d& dipole) const {
  static_assert(maxRank > 0);
  return charge * rank1_.dot(dipole);
}

template <Screening screen, Index maxRank>
double InteractionTensor<screen, maxRank>::interact(
    const Eigen::Vector3d& dipole1, const Eigen::Vector3d& dipole2) const {
  static_assert(maxRank > 1);
  return dipole1.transpose() * rank2_ * dipole2;
}

template <Screening screen, Index maxRank>
double InteractionTensor<screen, maxRank>::interact(
    double charge, const Eigen::Matrix3d& quadrupole) const {
  static_assert(maxRank > 1);
  // (mat1.array() * mat2.array()).sum() is the dyadic product
  return charge * (rank2_.array() * quadrupole.array()).sum();
}

template <Screening screen, Index maxRank>
double InteractionTensor<screen, maxRank>::interact(
    const Eigen::Vector3d& dipole, const Eigen::Matrix3d& quadrupole) const {
  static_assert(maxRank > 2);
  double result = 0.0;
  for (Index i = 0; i < 3; ++i) {
    for (Index j = 0; j < 3; ++j) {
      for (Index k = 0; k < 3; ++k) {
        result += dipole[k] * rank3_(k, j, i) * quadrupole(j, i);
      }
    }
  }
  return result;
}

template <Screening screen, Index maxRank>
double InteractionTensor<screen, maxRank>::interact(
    const Eigen::Matrix3d& quadrupole1,
    const Eigen::Matrix3d& quadrupole2) const {
  static_assert(maxRank > 3);
  double result = 0.0;
  for (Index i = 0; i < 3; ++i) {
    for (Index j = 0; j < 3; ++j) {
      for (Index k = 0; k < 3; ++k) {
        for (Index l = 0; l < 3; ++l) {
          result += quadrupole1(l, k) * rank4_(l, k, j, i) * quadrupole2(j, i);
        }
      }
    }
  }
  return result;
}

template <Screening screen, Index maxRank>
void InteractionTensor<screen, maxRank>::erfcTensors() {
  double R = dr_.norm();
  double rR = 1.0 / R;
  double rSqrtPiExp = rSqrtPi * std::exp(-a2 * R * R);
  double rR1s = std::erfc(a1 * R) * rR;
  rank0_ = rR1s;
  if constexpr (maxRank > 0) {
    double rR2 = rR * rR;
    double rR3s = rR2 * (rR1s + 2.0 * a1 * rSqrtPiExp);
    rank1_ = -l3 * rR3s * dr_;
    if constexpr (maxRank > 1) {
      double rR5s = rR2 * (rR3s + (4.0 / 3.0) * a3 * rSqrtPiExp);
      rank2_ = 3 * dr_ * dr_.transpose() * rR5s * l5 -
               Eigen::Matrix3d::Identity() * rR3s * l3;
      if constexpr (maxRank > 2) {
        double rR7s = rR2 * (rR5s + (8.0 / (5.0 * 3.0)) * a5 * rSqrtPiExp);
        for (Index i = 0; i < 3; ++i) {
          for (Index j = 0; j < 3; ++j) {
            for (Index k = 0; k < 3; ++k) {
              rank3_(k, j, i) =
                  -15 * dr_[k] * dr_[j] * dr_[i] * rR7s * l7 +
                  3 *
                      (dr_[k] * (double)(j == i) + dr_[j] * (double)(k == i) +
                       dr_[i] * (double)(k == j)) *
                      rR5s * l5;
            }
          }
        }
        if constexpr (maxRank > 3) {
          double rR9s =
              rR2 * (rR7s + (16.0 / (3.0 * 5.0 * 7.0)) * a7 * rSqrtPiExp);
          for (Index i = 0; i < 3; ++i) {
            for (Index j = 0; j < 3; ++j) {
              for (Index k = 0; k < 3; ++k) {
                for (Index l = 0; l < 3; ++l) {
                  double part1 = dr_[l] * dr_[k] * dr_[j] * dr_[i];
                  double part2 = dr_[l] * dr_[k] * (double)(j == i) +
                                 dr_[l] * dr_[j] * (double)(k == i) +
                                 dr_[l] * dr_[i] * (double)(k == j) +
                                 dr_[k] * dr_[j] * (double)(l == i) +
                                 dr_[k] * dr_[i] * (double)(l == j) +
                                 dr_[j] * dr_[i] * (double)(l == k);
                  double part3 =
                      (double)((l == k) * (j == i) + (l == j) * (k == i) +
                               (l == i) * (k == j));
                  rank4_(l, k, j, i) = 105 * part1 * rR9s * l9 -
                                       15 * part2 * rR7s * l7 +
                                       3 * part3 * rR5s * l5;
                }
              }
            }
          }
        }
      }
    }
  }
}

template <Screening screen, Index maxRank>
void InteractionTensor<screen, maxRank>::erfTensors() {
  double R = dr_.norm();
  double rR = 1.0 / R;
  double rSqrtPiExp = rSqrtPi * std::exp(-a2 * R * R);
  double rR1s = std::erf(a1 * R) * rR;
  rank0_ = rR1s;
  if constexpr (maxRank > 0) {
    double rR2 = rR * rR;
    double rR3s = rR2 * (rR1s - 2.0 * a1 * rSqrtPiExp);
    rank1_ = -rR3s * l3 * dr_;
    if constexpr (maxRank > 1) {
      double rR5s = rR2 * (rR3s - (4.0 / 3.0) * a3 * rSqrtPiExp);
      rank2_ = 3 * dr_ * dr_.transpose() * rR5s * l5 -
               Eigen::Matrix3d::Identity() * rR3s * l3;
      if constexpr (maxRank > 2) {
        double rR7s = rR2 * (rR5s - (8.0 / (5.0 * 3.0)) * a5 * rSqrtPiExp);
        for (Index i = 0; i < 3; ++i) {
          for (Index j = 0; j < 3; ++j) {
            for (Index k = 0; k < 3; ++k) {
              rank3_(k, j, i) =
                  -15 * dr_[k] * dr_[j] * dr_[i] * rR7s * l7 +
                  3 *
                      (dr_[k] * (double)(j == i) + dr_[j] * (double)(k == i) +
                       dr_[i] * (double)(k == j)) *
                      rR5s * l5;
            }
          }
        }
        if constexpr (maxRank > 3) {
          double rR9s =
              rR2 * (rR7s - (16.0 / (3.0 * 5.0 * 7.0)) * a7 * rSqrtPiExp);
          for (Index i = 0; i < 3; ++i) {
            for (Index j = 0; j < 3; ++j) {
              for (Index k = 0; k < 3; ++k) {
                for (Index l = 0; l < 3; ++l) {
                  double part1 = dr_[l] * dr_[k] * dr_[j] * dr_[i];
                  double part2 = dr_[l] * dr_[k] * (double)(j == i) +
                                 dr_[l] * dr_[j] * (double)(k == i) +
                                 dr_[l] * dr_[i] * (double)(k == j) +
                                 dr_[k] * dr_[j] * (double)(l == i) +
                                 dr_[k] * dr_[i] * (double)(l == j) +
                                 dr_[j] * dr_[i] * (double)(l == k);
                  double part3 =
                      (double)((l == k) * (j == i) + (l == j) * (k == i) +
                               (l == i) * (k == j));
                  rank4_(l, k, j, i) = 105 * part1 * rR9s * l9 -
                                       15 * part2 * rR7s * l7 +
                                       3 * part3 * rR5s * l5;
                }
              }
            }
          }
        }
      }
    }
  }
}

template <Screening screen, Index maxRank>
void InteractionTensor<screen, maxRank>::noneTensors() {
  double R = dr_.norm();
  double rR = 1.0 / R;
  rank0_ = rR;
  if constexpr (maxRank > 0) {
    double rR2 = rR * rR;
    double rR3 = rR * rR2;
    rank1_ = -l3 * rR3 * dr_;
    if constexpr (maxRank > 1) {
      double rR5 = rR3 * rR2;
      Eigen::Matrix3d dr_outer = dr_ * dr_.transpose();
      rank2_ = 3 * dr_outer * l5 * rR5 - Eigen::Matrix3d::Identity() * l3 * rR3;
      if constexpr (maxRank > 2) {
        double rR7 = rR5 * rR2;
        for (Index i = 0; i < 3; ++i) {
          for (Index j = 0; j < 3; ++j) {
            for (Index k = 0; k < 3; ++k) {
              rank3_(k, j, i) =
                  -15 * dr_[k] * dr_[j] * dr_[i] * l7 * rR7 +
                  3 *
                      (dr_[k] * (double)(j == i) + dr_[j] * (double)(k == i) +
                       dr_[i] * (double)(k == j)) *
                      l5 * rR5;
            }
          }
        }
        if constexpr (maxRank > 3) {
          double rR9 = rR7 * rR2;
          for (Index i = 0; i < 3; ++i) {
            for (Index j = 0; j < 3; ++j) {
              for (Index k = 0; k < 3; ++k) {
                for (Index l = 0; l < 3; ++l) {
                  double part1 = dr_[l] * dr_[k] * dr_[j] * dr_[i];

                  double part2 = dr_[l] * dr_[k] * (double)(j == i) +
                                 dr_[l] * dr_[j] * (double)(k == i) +
                                 dr_[l] * dr_[i] * (double)(k == j) +
                                 dr_[k] * dr_[j] * (double)(l == i) +
                                 dr_[k] * dr_[i] * (double)(l == j) +
                                 dr_[j] * dr_[i] * (double)(l == k);
                  double part3 =
                      (double)((l == k) * (j == i) + (l == j) * (k == i) +
                               (l == i) * (k == j));
                  rank4_(l, k, j, i) = 105 * part1 * rR9 * l9 -
                                       15 * part2 * rR7 * l7 +
                                       3 * part3 * rR5 * l5;
                }
              }
            }
          }
        }
      }
    }
  }
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_JOB_H
