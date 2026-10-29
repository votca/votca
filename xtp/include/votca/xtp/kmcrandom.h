/*
 *            Copyright 2009-2026 The VOTCA Development Team
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
#ifndef VOTCA_XTP_KMCRANDOM_H
#define VOTCA_XTP_KMCRANDOM_H

// Standard includes
#include <cmath>
#include <cstdint>
#include <random>

// Local VOTCA includes
#include "eigen.h"

namespace votca {
namespace xtp {

/**
 * \brief Random numbers for kinetic Monte Carlo: one independent stream per
 * (seed, stream) pair. std::seed_seq and mt19937_64 are fully specified by
 * the standard and the uniform numbers are built from the raw bits, so a
 * stream is identical with every standard library (unlike
 * std::uniform_real_distribution / std::normal_distribution).
 */
class KMCRandom {
 public:
  KMCRandom(std::uint64_t seed, std::uint64_t stream);
  /// uniform in (0,1]
  double Uniform() {
    return (double(engine_() >> 11) + 1.0) * (1.0 / 9007199254740992.0);
  }
  /// Exp(rate), i.e. -ln(u)/rate
  double Exponential(double rate) { return -std::log(Uniform()) / rate; }
  /// uniform integer in [0, n)
  Index Integer(Index n);
  /// standard normal (Box-Muller)
  double Normal();

 private:
  std::mt19937_64 engine_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCRANDOM_H
