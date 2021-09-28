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
#ifndef VOTCA_XTP_AOTRANSFORM_H
#define VOTCA_XTP_AOTRANSFORM_H

// Local VOTCA includes
#include "aoshell.h"
#include "eigen.h"

namespace votca {
namespace xtp {
// clang-format off
//clang format puts one entry on each line
namespace Cart{
enum Cart {
  s,  // s
  x,  y,  z,  // p
  xx,  xy,  xz,  yy,  yz,  zz,  // d
  xxx,  xxy,  xxz,  xyy,  xyz,  xzz,  yyy,  yyz,  yzz,  zzz,  // f
  
  xxxx,  xxxy,  xxxz,  xxyy,  xxyz,  xxzz,  xyyy,  xyyz,  xyzz,  xzzz,  yyyy,
  yyyz,  yyzz,  yzzz,  zzzz,  // g
  
  xxxxx,  xxxxy,  xxxxz,  xxxyy,  xxxyz,  xxxzz,  xxyyy,  xxyyz,  xxyzz,
  xxzzz,  xyyyy,  xyyyz,  xyyzz,  xyzzz,  xzzzz,  yyyyy,  yyyyz,  yyyzz,
  yyzzz,  yzzzz,  zzzzz,  // h

  xxxxxx,  xxxxxy,  xxxxxz,  xxxxyy,  xxxxyz,  xxxxzz,  xxxyyy,  xxxyyz,
  xxxyzz,  xxxzzz,  xxyyyy,  xxyyyz,  xxyyzz,  xxyzzz,  xxzzzz,  xyyyyy,
  xyyyyz,  xyyyzz,  xyyzzz,  xyzzzz,  xzzzzz,  yyyyyy,  yyyyyz,  yyyyzz,
  yyyzzz,  yyzzzz,  yzzzzz,  zzzzzz, //i
  
  xxxxxxx,  xxxxxxy,  xxxxxxz,  xxxxxyy,  xxxxxyz,  xxxxxzz,  xxxxyyy,
  xxxxyyz,  xxxxyzz,  xxxxzzz,  xxxyyyy,  xxxyyyz,  xxxyyzz,  xxxyzzz,
  xxxzzzz,  xxyyyyy,  xxyyyyz,  xxyyyzz,  xxyyzzz,  xxyzzzz,  xxzzzzz,
  xyyyyyy,  xyyyyyz,  xyyyyzz,  xyyyzzz,  xyyzzzz,  xyzzzzz,  xzzzzzz,
  yyyyyyy,  yyyyyyz,  yyyyyzz,  yyyyzzz,  yyyzzzz,  yyzzzzz,  yzzzzzz,
  zzzzzzz, //j

  xxxxxxxx,  xxxxxxxy,  xxxxxxxz,  xxxxxxyy,  xxxxxxyz,  xxxxxxzz,  xxxxxyyy,
  xxxxxyyz,  xxxxxyzz,  xxxxxzzz,  xxxxyyyy,  xxxxyyyz,  xxxxyyzz,  xxxxyzzz,
  xxxxzzzz,  xxxyyyyy,  xxxyyyyz,  xxxyyyzz,  xxxyyzzz,  xxxyzzzz,  xxxzzzzz,
  xxyyyyyy,  xxyyyyyz,  xxyyyyzz,  xxyyyzzz,  xxyyzzzz,  xxyzzzzz,  xxzzzzzz,
  xyyyyyyy,  xyyyyyyz,  xyyyyyzz,  xyyyyzzz,  xyyyzzzz,  xyyzzzzz,  xyzzzzzz,
  xzzzzzzz,  yyyyyyyy,  yyyyyyyz,  yyyyyyzz,  yyyyyzzz,  yyyyzzzz,  yyyzzzzz,
  yyzzzzzz,  yzzzzzzz,  zzzzzzzz, //k
};
}

// clang-format on

/* contains cartesian to spherical conversion
 */

class AOTransform {

 public:
  static std::array<int, 9> n_orbitals();
  static std::array<int, 165> nx();
  static std::array<int, 165> ny();
  static std::array<int, 165> nz();
  static std::array<int, 165> i_less_x();
  static std::array<int, 165> i_less_y();
  static std::array<int, 165> i_less_z();
  static std::array<int, 120> i_more_x();
  static std::array<int, 120> i_more_y();
  static std::array<int, 120> i_more_z();

  static Index getBlockSize(Index lmax);
  static Eigen::VectorXd XIntegrate(Index size, double U);

  static double getNorm(L l, const AOGaussianPrimitive& gaussian);

  /// transforms a cartesian shell to a spherical cartesian shell
  template <typename Matrix>
  static Matrix tform(L l_row, L l_col, const Matrix& cartesian);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_AOTRANSFORM_H
