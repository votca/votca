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
#ifndef VOTCA_XTP_CUBEFILE_WRITER_H
#define VOTCA_XTP_CUBEFILE_WRITER_H

// Local VOTCA includes
#include "logger.h"
#include "orbitals.h"
#include "regular_grid.h"

/**
 * \brief Writes an orbital file to a .cube file
 */

namespace votca {
namespace xtp {
class CubeFile_Writer {

 public:
  CubeFile_Writer(Eigen::Array<Index, 3, 1> steps, double padding, Logger& log)
      : steps_(steps), padding_(padding), log_(log) {};

  void WriteFile(const std::string& filename, const Orbitals& orb,
                 QMState state, bool dostateonly) const;

 private:
  std::vector<std::vector<double> > CalculateValues(
      const Orbitals& orb, QMState state, bool dostateonly,
      const Regular_Grid& grid) const;

  Eigen::Array<Index, 3, 1> steps_;
  double padding_;
  Logger& log_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CUBEFILE_WRITER_H
