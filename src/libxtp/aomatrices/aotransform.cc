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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include <votca/xtp/aotransform.h>
namespace votca {
namespace xtp {

Index AOTransform::getCartesianSize(Index l) { return (l + 1) * (l + 2) / 2; }
Index AOTransform::getSphericalSize(Index l) { return 2 * l + 1; }

Eigen::MatrixXd AOTransform::getPrimitiveShellTrafo(Index l, double decay,
                                                    double contraction) {
  switch (l) {
    case 0: {
      return contraction * Eigen::MatrixXd::Ones(1, 1);
    }
    case 1: {
      Eigen::MatrixXd trafo = Eigen::MatrixXd::Zero(3, 3);
      const double factor = 2. * sqrt(decay) * contraction;
      trafo(P::z, 0) = factor;  // Y 1,0
      trafo(P::y, 1) = factor;  // Y 1,-1
      trafo(P::x, 2) = factor;  // Y 1,1
      return trafo;
    }
    case 2: {
      Eigen::MatrixXd trafo = Eigen::MatrixXd::Zero(6, 5);
      const double factor = 2. * decay * contraction;
      const double factor_1 = factor / sqrt(3.);
      trafo(D::xx, 0) = -factor_1;
      trafo(D::yy, 0) = -factor_1;  //   Y 2,0
      trafo(D::zz, 0) = 2. * factor_1;

      trafo(D::yz, 1) = 2. * factor;  //           Y 2,-1

      trafo(D::xz, 2) = 2. * factor;  //           Y 2,1

      trafo(D::xy, 3) = 2. * factor;  //          Y 2,-2

      trafo(D::xx, 4) = factor;  //   Y 2,2
      trafo(D::yy, 4) = -factor;
      return trafo;
    }
    case 3: {
      Eigen::MatrixXd trafo = Eigen::MatrixXd::Zero(10, 7);
      const double factor = 2. * pow(decay, 1.5) * contraction;
      const double factor_1 = factor * 2. / sqrt(15.);
      const double factor_2 = factor * sqrt(2.) / sqrt(5.);
      const double factor_3 = factor * sqrt(2.) / sqrt(3.);

      trafo(F::xxz, 0) = -3. * factor_1;
      trafo(F::yyz, 0) = -3. * factor_1;  //        Y 3,0
      trafo(F::zzz, 0) = 2. * factor_1;

      trafo(F::xxy, 1) = -factor_2;
      trafo(F::yyy, 1) = -factor_2;  //    Y 3,-1
      trafo(F::yzz, 1) = 4. * factor_2;

      trafo(F::xxx, 2) = -factor_2;
      trafo(F::xyy, 2) = -factor_2;  //    Y 3,1
      trafo(F::xzz, 2) = 4. * factor_2;

      trafo(F::xyz, 3) = 4. * factor;  //      Y 3,-2

      trafo(F::xxz, 4) = 2. * factor;
      trafo(F::yyz, 4) = -2. * factor;  //   Y 3,2

      trafo(F::xxy, 5) = 3. * factor_3;
      trafo(F::yyy, 5) = -factor_3;  //    Y 3,-3

      trafo(F::xxx, 6) = factor_3;
      trafo(F::xyy, 6) = -3. * factor_3;  //    Y 3,3
      return trafo;
    }
    case 4: {
      Eigen::MatrixXd trafo = Eigen::MatrixXd::Zero(15, 9);
      const double factor = 2. / sqrt(3.) * decay * decay * contraction;
      const double factor_1 = factor / sqrt(35.);
      const double factor_2 = factor * 4. / sqrt(14.);
      const double factor_3 = factor * 2. / sqrt(7.);
      const double factor_4 = factor * 2. * sqrt(2.);

      trafo(0, 0) = 3. * factor_1;  /// Y 4,0
      trafo(G::xxyy, 0) = 6. * factor_1;
      trafo(G::xxzz, 0) = -24. * factor_1;
      trafo(G::yyyy, 0) = 3. * factor_1;
      trafo(G::yyzz, 0) = -24. * factor_1;
      trafo(G::zzzz, 0) = 8. * factor_1;

      trafo(G::xxyz, 1) = -3. * factor_2;  /// Y 4,-1
      trafo(G::yyyz, 1) = -3. * factor_2;
      trafo(G::yzzz, 1) = 4. * factor_2;

      trafo(G::xxxz, 2) = -3. * factor_2;  /// Y 4,1
      trafo(G::xyyz, 2) = -3. * factor_2;
      trafo(G::xzzz, 2) = 4. * factor_2;

      trafo(G::xxxy, 3) = -2. * factor_3;  /// Y 4,-2
      trafo(G::xyyy, 3) = -2. * factor_3;
      trafo(G::xyzz, 3) = 12. * factor_3;

      trafo(G::xxxx, 4) = -factor_3;  /// Y 4,2
      trafo(G::xxzz, 4) = 6. * factor_3;
      trafo(G::yyyy, 4) = factor_3;
      trafo(G::yyzz, 4) = -6. * factor_3;

      trafo(G::xxyz, 5) = 3. * factor_4;  /// Y 4,-3
      trafo(G::yyyz, 5) = -factor_4;

      trafo(G::xxxz, 6) = factor_4;  /// Y 4,3
      trafo(G::xyyz, 6) = -3. * factor_4;

      trafo(G::xxxy, 7) = 4. * factor;  /// Y 4,-4
      trafo(G::xyyy, 7) = -4. * factor;

      trafo(G::xxxx, 8) = factor;  /// Y 4,4
      trafo(G::xxyy, 8) = -6. * factor;
      trafo(G::yyyy, 8) = factor;
      return trafo;
    }
    case 5: {
      Eigen::MatrixXd trafo = Eigen::MatrixXd::Zero(21, 11);
      const double factor = (2. / 3.) * std::pow(decay, 2.5) * contraction;
      const double factor_1 = factor * 2. / sqrt(105.);
      const double factor_2 = factor * 2. / sqrt(7.);
      const double factor_3 = factor * sqrt(6.) / 3.;
      const double factor_4 = factor * 2. * sqrt(3.);
      const double factor_5 = factor * .2 * sqrt(30.);

      trafo(H::xxxxz, 0) = 15. * factor_1;  /// Y 5,0
      trafo(H::xxyyz, 0) = 30. * factor_1;
      trafo(H::xxzzz, 0) = -40. * factor_1;
      trafo(H::yyyyz, 0) = 15. * factor_1;
      trafo(H::yyzzz, 0) = -40. * factor_1;
      trafo(H::zzzzz, 0) = 8. * factor_1;

      trafo(H::xxxxy, 1) = factor_2;  /// Y 5,-1
      trafo(H::xxyyy, 1) = 2. * factor_2;
      trafo(H::xxyzz, 1) = -12. * factor_2;
      trafo(H::yyyyy, 1) = factor_2;
      trafo(H::yyyzz, 1) = -12. * factor_2;
      trafo(H::yzzzz, 1) = 8. * factor_2;

      trafo(H::xxxxx, 2) = factor_2;  /// Y 5,1
      trafo(H::xxxyy, 2) = 2. * factor_2;
      trafo(H::xxxzz, 2) = -12. * factor_2;
      trafo(H::xyyyy, 2) = factor_2;
      trafo(H::xyyzz, 2) = -12. * factor_2;
      trafo(H::xzzzz, 2) = 8. * factor_2;

      trafo(H::xxxyz, 3) = -8. * factor;  /// Y 5,-2
      trafo(H::xyyyz, 3) = -8. * factor;
      trafo(H::xyzzz, 3) = 16. * factor;

      trafo(H::xxxxz, 4) = -4. * factor;  /// Y 5,2
      trafo(H::xxzzz, 4) = 8. * factor;
      trafo(H::yyyyz, 4) = 4. * factor;
      trafo(H::yyzzz, 4) = -8. * factor;

      trafo(H::xxxxy, 5) = -3. * factor_3;  /// Y 5,-3
      trafo(H::xxyyy, 5) = -2. * factor_3;
      trafo(H::xxyzz, 5) = 24. * factor_3;
      trafo(H::yyyyy, 5) = factor_3;
      trafo(H::yyyzz, 5) = -8. * factor_3;

      trafo(H::xxxxx, 6) = -factor_3;  /// Y 5,3
      trafo(H::xxxyy, 6) = 2. * factor_3;
      trafo(H::xxxzz, 6) = 8. * factor_3;
      trafo(H::xyyyy, 6) = 3. * factor_3;
      trafo(H::xyyzz, 6) = -24. * factor_3;

      trafo(H::xxxyz, 7) = 4. * factor_4;  /// Y 5,-4
      trafo(H::xyyyz, 7) = -4. * factor_4;

      trafo(H::xxxxz, 8) = factor_4;  /// Y 5,4
      trafo(H::xxyyz, 8) = -6. * factor_4;
      trafo(H::yyyyz, 8) = factor_4;

      trafo(H::xxxxy, 9) = 5. * factor_5;  /// Y 5,-5
      trafo(H::xxyyy, 9) = -10. * factor_5;
      trafo(H::yyyyy, 9) = factor_5;

      trafo(H::xxxxx, 10) = factor_5;  /// Y 5,5
      trafo(H::xxxyy, 10) = -10. * factor_5;
      trafo(H::xyyyy, 10) = 5. * factor_5;
      return trafo;
    }
    case 6: {
      Eigen::MatrixXd trafo = Eigen::MatrixXd::Zero(28, 13);
      const double factor = (2. / 3.) * decay * decay * decay * contraction;
      const double factor_1 = factor * 2. / sqrt(1155.);
      const double factor_2 = factor * 4. / sqrt(55.);
      const double factor_3 = factor * sqrt(22.) / 11.;
      const double factor_4 = factor * 2. * sqrt(165.) / 55.;
      const double factor_5 = factor * .4 * sqrt(30.);
      const double factor_6 = factor * .2 * sqrt(10.);

      trafo(I::xxxxxx, 0) = -5. * factor_1;  /// Y 6,0
      trafo(I::xxxxyy, 0) = -15. * factor_1;
      trafo(I::xxxxzz, 0) = 90. * factor_1;
      trafo(I::xxyyyy, 0) = -15. * factor_1;
      trafo(I::xxyyzz, 0) = 180. * factor_1;
      trafo(I::xxzzzz, 0) = -120. * factor_1;
      trafo(I::yyyyyy, 0) = -5. * factor_1;
      trafo(I::yyyyzz, 0) = 90. * factor_1;
      trafo(I::yyzzzz, 0) = -120. * factor_1;
      trafo(I::zzzzzz, 0) = 16. * factor_1;

      trafo(I::xxxxyz, 1) = 5. * factor_2;  /// Y 6,-1
      trafo(I::xxyyyz, 1) = 10. * factor_2;
      trafo(I::xxyzzz, 1) = -20. * factor_2;
      trafo(I::yyyyyz, 1) = 5. * factor_2;
      trafo(I::yyyzzz, 1) = -20. * factor_2;
      trafo(I::yzzzzz, 1) = 8. * factor_2;

      trafo(I::xxxxxz, 2) = 5. * factor_2;  /// Y 6,1
      trafo(I::xxxyyz, 2) = 10. * factor_2;
      trafo(I::xxxzzz, 2) = -20. * factor_2;
      trafo(I::xyyyyz, 2) = 5. * factor_2;
      trafo(I::xyyzzz, 2) = -20. * factor_2;
      trafo(I::xzzzzz, 2) = 8. * factor_2;

      trafo(I::xxxxxy, 3) = 2. * factor_3;  /// Y 6,-2
      trafo(I::xxxyyy, 3) = 4. * factor_3;
      trafo(I::xxxyzz, 3) = -32. * factor_3;
      trafo(I::xyyyyy, 3) = 2. * factor_3;
      trafo(I::xyyyzz, 3) = -32. * factor_3;
      trafo(I::xyzzzz, 3) = 32. * factor_3;

      trafo(I::xxxxxy, 4) = factor_3;  /// Y 6,2
      trafo(I::xxxxyy, 4) = factor_3;
      trafo(I::xxxxzz, 4) = -16. * factor_3;
      trafo(I::xxyyyy, 4) = -factor_3;
      trafo(I::xxzzzz, 4) = 16. * factor_3;
      trafo(I::yyyyyy, 4) = -factor_3;
      trafo(I::yyyyzz, 4) = 16. * factor_3;
      trafo(I::yyzzzz, 4) = -16. * factor_3;

      trafo(I::xxxxyz, 5) = -18. * factor_3;  /// Y 6,-3
      trafo(I::xxyyyz, 5) = -12. * factor_3;
      trafo(I::xxyzzz, 5) = 48. * factor_3;
      trafo(I::yyyyyz, 5) = 6. * factor_3;
      trafo(I::yyyzzz, 5) = -16. * factor_3;

      trafo(I::xxxxxz, 6) = -6. * factor_3;  /// Y 6,3
      trafo(I::xxxyyz, 6) = 12. * factor_3;
      trafo(I::xxxzzz, 6) = 16. * factor_3;
      trafo(I::xyyyyz, 6) = 18. * factor_3;
      trafo(I::xyyzzz, 6) = -48. * factor_3;

      trafo(I::xxxxxy, 7) = -4. * factor_4;  /// Y 6,-4
      trafo(I::xxxyzz, 7) = 40. * factor_4;
      trafo(I::xyyyyy, 7) = 4. * factor_4;
      trafo(I::xyyyzz, 7) = -40. * factor_4;

      trafo(I::xxxxxx, 8) = -factor_4;  /// Y 6,4
      trafo(I::xxxxyy, 8) = 5. * factor_4;
      trafo(I::xxxxzz, 8) = 10. * factor_4;
      trafo(I::xxyyyy, 8) = 5. * factor_4;
      trafo(I::xxyyzz, 8) = -60. * factor_4;
      trafo(I::yyyyyy, 8) = -factor_4;
      trafo(I::yyyyzz, 8) = 10. * factor_4;

      trafo(I::xxxxyz, 9) = 5. * factor_5;  /// Y 6,-5
      trafo(I::xxyyyz, 9) = -10. * factor_5;
      trafo(I::yyyyyz, 9) = factor_5;

      trafo(I::xxxxxz, 10) = factor_5;  /// Y 6,5
      trafo(I::xxxyyz, 10) = -10. * factor_5;
      trafo(I::xyyyyz, 10) = 5. * factor_5;

      trafo(I::xxxxxy, 11) = 6. * factor_6;  /// Y 6,-6
      trafo(I::xxxyyy, 11) = -20. * factor_6;
      trafo(I::xyyyyy, 11) = 6. * factor_6;

      trafo(I::xxxxxx, 12) = factor_6;  /// Y 6,6
      trafo(I::xxxxyy, 12) = -15. * factor_6;
      trafo(I::xxyyyy, 12) = 15. * factor_6;
      trafo(I::yyyyyy, 12) = -factor_6;
      return trafo;
    }
    default:
      throw std::runtime_error("No transforms for shells higher than i (l=6)");
  }

  return Eigen::MatrixXd(0, 0);
}

Eigen::MatrixXd AOTransform::getTrafo(const AOGaussianPrimitive& gaussian) {

  const AOShell& shell = gaussian.getShell();
  if (!shell.isCombined()) {
    return AOTransform::getPrimitiveShellTrafo(
        shell.getLmax(), gaussian.getDecay(),
        gaussian.getContraction()(shell.getLmax()));
  }

  Eigen::MatrixXd trafo =
      Eigen::MatrixXd::Zero(shell.getCartesianNumFunc(), shell.getNumFunc());

  Index rowoffset = 0;
  Index coloffset = 0;
  for (Index l = shell.getLmin(); l <= shell.getLmax(); l++) {
    Index cart_size = AOTransform::getCartesianSize(l);
    Index spherical_size = AOTransform::getSphericalSize(l);
    trafo.block(rowoffset, coloffset, cart_size, spherical_size) =
        AOTransform::getPrimitiveShellTrafo(l, gaussian.getDecay(),
                                            gaussian.getContraction()(l));
    rowoffset += cart_size;
    coloffset += spherical_size;
  }
  return trafo;
}

Eigen::VectorXd AOTransform::XIntegrate(Index size, double U) {
  Eigen::VectorXd FmU = Eigen::VectorXd::Zero(size);
  const Index mm = size - 1;
  const double pi = boost::math::constants::pi<double>();
  if (mm < 0) {
    throw std::runtime_error("mm is: " + std::to_string(mm) +
                             " This should not have happened!");
  }

  if (U < 0.0) {
    throw std::runtime_error("U is: " + std::to_string(U) +
                             " This should not have happened!");
  }

  if (U >= 10.0) {
    // forward iteration
    FmU[0] = 0.50 * std::sqrt(pi / U) * std::erf(std::sqrt(U));

    const double expU = std::exp(-U);
    for (Index m = 1; m < FmU.size(); m++) {
      FmU[m] =
          (2.0 * double(m) - 1) * FmU[m - 1] / (2.0 * U) - expU / (2.0 * U);
    }
  }

  else if (U < 1e-10) {
    for (Index m = 0; m < FmU.size(); m++) {
      FmU[m] = 1.0 / (2.0 * double(m) + 1.0) - U / (2.0 * double(m) + 3.0);
    }
  }

  else if (U >= 1e-10 && U < 10.0) {
    // backward iteration
    double fm = 0.0;
    const double expU = std::exp(-U);
    for (Index m = 60; m >= mm; m--) {
      fm = (2.0 * U) / (2.0 * double(m) + 1.0) * (fm + expU / (2.0 * U));
    }
    FmU[mm] = fm;
    for (Index m = mm - 1; m >= 0; m--) {
      FmU[m] =
          (2.0 * U) / (2.0 * double(m) + 1.0) * (FmU[m + 1] + expU / (2.0 * U));
    }
  }

  return FmU;
}

Index AOTransform::getBlockSize(Index lmax) {
  // Each cartesian shells has (l+1)(l+2)/2 elements
  // Sum of all shells up to _lmax leads to blocksize=1+11/6 l+l^2+1/6 l^3
  Index blocksize = 6 + 11 * lmax + 6 * lmax * lmax + lmax * lmax * lmax;
  blocksize /= 6;
  return blocksize;
}

// blockSize till l=8
std::array<int, 9> AOTransform::n_orbitals() {
  return {1, 4, 10, 20, 35, 56, 84, 120, 165};
}

std::array<int, 165> AOTransform::nx() {
  return {0, 1, 0, 0, 2, 1, 1, 0, 0, 0, 3, 2, 2, 1, 1, 1, 0, 0, 0, 0, 4,
          3, 3, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0, 0, 5, 4, 4, 3, 3, 3, 2,
          2, 2, 2, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 6, 5, 5, 4, 4, 4, 3,
          3, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
          7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,
          1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 8, 7, 7, 6, 6, 6,
          5, 5, 5, 5, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2,
          2, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}

std::array<int, 165> AOTransform::ny() {
  return {0, 0, 1, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 0,
          1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3,
          2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0, 3,
          2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1, 0,
          0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0,
          6, 5, 4, 3, 2, 1, 0, 7, 6, 5, 4, 3, 2, 1, 0, 0, 1, 0, 2, 1, 0,
          3, 2, 1, 0, 4, 3, 2, 1, 0, 5, 4, 3, 2, 1, 0, 6, 5, 4, 3, 2, 1,
          0, 7, 6, 5, 4, 3, 2, 1, 0, 8, 7, 6, 5, 4, 3, 2, 1, 0};
}

std::array<int, 165> AOTransform::nz() {
  return {0, 0, 0, 1, 0, 0, 1, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0,
          0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 0, 1, 0, 1, 2, 0,
          1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 0, 1, 0, 1, 2, 0,
          1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 6,
          0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5,
          0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 5, 6, 7, 0, 0, 1, 0, 1, 2,
          0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5,
          6, 0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 5, 6, 7, 8};
}

std::array<int, 165> AOTransform::i_less_x() {
  return {0,   0,   0,   0,   1,   2,   3,   0,   0,   0,   4,   5,   6,   7,
          8,   9,   0,   0,   0,   0,   10,  11,  12,  13,  14,  15,  16,  17,
          18,  19,  0,   0,   0,   0,   0,   20,  21,  22,  23,  24,  25,  26,
          27,  28,  29,  30,  31,  32,  33,  34,  0,   0,   0,   0,   0,   0,
          35,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,
          49,  50,  51,  52,  53,  54,  55,  0,   0,   0,   0,   0,   0,   0,
          56,  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,
          70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,
          0,   0,   0,   0,   0,   0,   0,   0,   84,  85,  86,  87,  88,  89,
          90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103,
          104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
          118, 119, 0,   0,   0,   0,   0,   0,   0,   0,   0};
}
std::array<int, 165> AOTransform::i_less_y() {
  return {0,   0,  0,   0,   0,   1,   0,   2,   3,   0,   0,   4,   0,   5,
          6,   0,  7,   8,   9,   0,   0,   10,  0,   11,  12,  0,   13,  14,
          15,  0,  16,  17,  18,  19,  0,   0,   20,  0,   21,  22,  0,   23,
          24,  25, 0,   26,  27,  28,  29,  0,   30,  31,  32,  33,  34,  0,
          0,   35, 0,   36,  37,  0,   38,  39,  40,  0,   41,  42,  43,  44,
          0,   45, 46,  47,  48,  49,  0,   50,  51,  52,  53,  54,  55,  0,
          0,   56, 0,   57,  58,  0,   59,  60,  61,  0,   62,  63,  64,  65,
          0,   66, 67,  68,  69,  70,  0,   71,  72,  73,  74,  75,  76,  0,
          77,  78, 79,  80,  81,  82,  83,  0,   0,   84,  0,   85,  86,  0,
          87,  88, 89,  0,   90,  91,  92,  93,  0,   94,  95,  96,  97,  98,
          0,   99, 100, 101, 102, 103, 104, 0,   105, 106, 107, 108, 109, 110,
          111, 0,  112, 113, 114, 115, 116, 117, 118, 119, 0};
}

std::array<int, 165> AOTransform::i_less_z() {
  return {0,   0,   0,  0,   0,   0,   1,   0,   2,   3,   0,   0,   4,   0,
          5,   6,   0,  7,   8,   9,   0,   0,   10,  0,   11,  12,  0,   13,
          14,  15,  0,  16,  17,  18,  19,  0,   0,   20,  0,   21,  22,  0,
          23,  24,  25, 0,   26,  27,  28,  29,  0,   30,  31,  32,  33,  34,
          0,   0,   35, 0,   36,  37,  0,   38,  39,  40,  0,   41,  42,  43,
          44,  0,   45, 46,  47,  48,  49,  0,   50,  51,  52,  53,  54,  55,
          0,   0,   56, 0,   57,  58,  0,   59,  60,  61,  0,   62,  63,  64,
          65,  0,   66, 67,  68,  69,  70,  0,   71,  72,  73,  74,  75,  76,
          0,   77,  78, 79,  80,  81,  82,  83,  0,   0,   84,  0,   85,  86,
          0,   87,  88, 89,  0,   90,  91,  92,  93,  0,   94,  95,  96,  97,
          98,  0,   99, 100, 101, 102, 103, 104, 0,   105, 106, 107, 108, 109,
          110, 111, 0,  112, 113, 114, 115, 116, 117, 118, 119};
}

std::array<int, 120> AOTransform::i_more_x() {
  return {1,   4,   5,   6,   10,  11,  12,  13,  14,  15,  20,  21,  22,  23,
          24,  25,  26,  27,  28,  29,  35,  36,  37,  38,  39,  40,  41,  42,
          43,  44,  45,  46,  47,  48,  49,  56,  57,  58,  59,  60,  61,  62,
          63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,
          84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,
          98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
          120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
          134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147,
          148, 149, 150, 151, 152, 153, 154, 155};
}

std::array<int, 120> AOTransform::i_more_y() {
  return {2,   5,   7,   8,   11,  13,  14,  16,  17,  18,  21,  23,  24,  26,
          27,  28,  30,  31,  32,  33,  36,  38,  39,  41,  42,  43,  45,  46,
          47,  48,  50,  51,  52,  53,  54,  57,  59,  60,  62,  63,  64,  66,
          67,  68,  69,  71,  72,  73,  74,  75,  77,  78,  79,  80,  81,  82,
          85,  87,  88,  90,  91,  92,  94,  95,  96,  97,  99,  100, 101, 102,
          103, 105, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
          121, 123, 124, 126, 127, 128, 130, 131, 132, 133, 135, 136, 137, 138,
          139, 141, 142, 143, 144, 145, 146, 148, 149, 150, 151, 152, 153, 154,
          156, 157, 158, 159, 160, 161, 162, 163};
}

std::array<int, 120> AOTransform::i_more_z() {
  return {3,   6,   8,   9,   12,  14,  15,  17,  18,  19,  22,  24,  25,  27,
          28,  29,  31,  32,  33,  34,  37,  39,  40,  42,  43,  44,  46,  47,
          48,  49,  51,  52,  53,  54,  55,  58,  60,  61,  63,  64,  65,  67,
          68,  69,  70,  72,  73,  74,  75,  76,  78,  79,  80,  81,  82,  83,
          86,  88,  89,  91,  92,  93,  95,  96,  97,  98,  100, 101, 102, 103,
          104, 106, 107, 108, 109, 110, 111, 113, 114, 115, 116, 117, 118, 119,
          122, 124, 125, 127, 128, 129, 131, 132, 133, 134, 136, 137, 138, 139,
          140, 142, 143, 144, 145, 146, 147, 149, 150, 151, 152, 153, 154, 155,
          157, 158, 159, 160, 161, 162, 163, 164};
}

}  // namespace xtp
}  // namespace votca
