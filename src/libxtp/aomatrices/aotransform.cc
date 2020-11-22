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
 *Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

// Local VOTCA includes
#include "votca/xtp/aotransform.h"

namespace votca {
namespace xtp {

Eigen::VectorXd AOTransform::getNorm(L l, double contraction, double decay) {
  switch (l) {
    case L::S: {
      return contraction * Eigen::VectorXd::Ones(1);
    }
    case L::P: {
      return 2. * std::sqrt(decay) * contraction * Eigen::VectorXd::Ones(3);
    }
    case L::D: {
      Eigen::VectorXd norm = Eigen::VectorXd::Zero(6);
      const double factor = 2. * decay * contraction;
      const double factor_1 = factor / sqrt(3.);
      norm(D::xx) = factor_1;
      norm(D::xy) = factor;
      norm(D::xz) = factor;
      norm(D::yy) = factor_1;
      norm(D::yz) = factor;
      norm(D::zz) = factor_1;
    }
    case L::F: {
      Eigen::VectorXd norm = Eigen::VectorXd::Zero(10);
      const double factor = 2. * pow(decay, 1.5) * contraction;
      const double factor_1 = factor * 2. / sqrt(15.);
      const double factor_2 = factor * sqrt(2.) / sqrt(5.);
      const double factor_3 = factor * sqrt(2.) / sqrt(3.);

      norm(F::xxx) = factor_3;
      norm(F::xxy) = -factor_2;
      norm(F::xxz) = -3. * factor_1;
      norm(F::xyz) = 4. * factor;  //      Y 3,-2

      norm(F::yyy) = factor_3;
      norm(F::yyz) = -3. * factor_1;  //        Y 3,0
      norm(F::yzz) = 4. * factor_2;

      norm(F::zzz) = 2. * factor_1;
    }
  }
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
        FmU[m] = (2.0 * U) / (2.0 * double(m) + 1.0) *
                 (FmU[m + 1] + expU / (2.0 * U));
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
    return {0,   0,   0,   0,   1,   2,   3,   0,   0,   0,   4,   5,   6,
            7,   8,   9,   0,   0,   0,   0,   10,  11,  12,  13,  14,  15,
            16,  17,  18,  19,  0,   0,   0,   0,   0,   20,  21,  22,  23,
            24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  0,   0,
            0,   0,   0,   0,   35,  36,  37,  38,  39,  40,  41,  42,  43,
            44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  0,
            0,   0,   0,   0,   0,   0,   56,  57,  58,  59,  60,  61,  62,
            63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,
            76,  77,  78,  79,  80,  81,  82,  83,  0,   0,   0,   0,   0,
            0,   0,   0,   84,  85,  86,  87,  88,  89,  90,  91,  92,  93,
            94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106,
            107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
            0,   0,   0,   0,   0,   0,   0,   0,   0};
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
    return {1,   4,   5,   6,   10,  11,  12,  13,  14,  15,  20,  21,
            22,  23,  24,  25,  26,  27,  28,  29,  35,  36,  37,  38,
            39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  56,
            57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,
            69,  70,  71,  72,  73,  74,  75,  76,  84,  85,  86,  87,
            88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,
            100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111,
            120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131,
            132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143,
            144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155};
  }

  std::array<int, 120> AOTransform::i_more_y() {
    return {2,   5,   7,   8,   11,  13,  14,  16,  17,  18,  21,  23,
            24,  26,  27,  28,  30,  31,  32,  33,  36,  38,  39,  41,
            42,  43,  45,  46,  47,  48,  50,  51,  52,  53,  54,  57,
            59,  60,  62,  63,  64,  66,  67,  68,  69,  71,  72,  73,
            74,  75,  77,  78,  79,  80,  81,  82,  85,  87,  88,  90,
            91,  92,  94,  95,  96,  97,  99,  100, 101, 102, 103, 105,
            106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118,
            121, 123, 124, 126, 127, 128, 130, 131, 132, 133, 135, 136,
            137, 138, 139, 141, 142, 143, 144, 145, 146, 148, 149, 150,
            151, 152, 153, 154, 156, 157, 158, 159, 160, 161, 162, 163};
  }

  std::array<int, 120> AOTransform::i_more_z() {
    return {3,   6,   8,   9,   12,  14,  15,  17,  18,  19,  22,  24,
            25,  27,  28,  29,  31,  32,  33,  34,  37,  39,  40,  42,
            43,  44,  46,  47,  48,  49,  51,  52,  53,  54,  55,  58,
            60,  61,  63,  64,  65,  67,  68,  69,  70,  72,  73,  74,
            75,  76,  78,  79,  80,  81,  82,  83,  86,  88,  89,  91,
            92,  93,  95,  96,  97,  98,  100, 101, 102, 103, 104, 106,
            107, 108, 109, 110, 111, 113, 114, 115, 116, 117, 118, 119,
            122, 124, 125, 127, 128, 129, 131, 132, 133, 134, 136, 137,
            138, 139, 140, 142, 143, 144, 145, 146, 147, 149, 150, 151,
            152, 153, 154, 155, 157, 158, 159, 160, 161, 162, 163, 164};
  }

}  // namespace xtp
}  // namespace xtp
