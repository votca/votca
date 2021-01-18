/*
 *            Copyright 2009-2021 The VOTCA Development Team
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
#ifndef VOTCA_XTP_GAUSSIANWRITER_H
#define VOTCA_XTP_GAUSSIANWRITER_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>

namespace votca {
namespace xtp {
class GaussianWriter {
 public:
  GaussianWriter(Logger& log) : _log(log){};

  ~GaussianWriter() = default;

  void WriteFile(const std::string& filename, const Orbitals& orbitals) const;

 private:
  Logger& _log;
  Index toGaussianL(L l) const;
  std::string reorderedMOCoefficients(const Orbitals& orbitals) const;
  std::string densityMatrixToString(const Orbitals& orbitals) const;
};

}  // namespace xtp
}  // namespace votca

#endif