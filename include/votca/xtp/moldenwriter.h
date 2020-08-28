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
#ifndef VOTCA_XTP_MOLDEN_WRITER_H
#define VOTCA_XTP_MOLDEN_WRITER_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
class MoldenWriter {
 public:
  void WriteFile(const std::string& filename, const Orbitals& orbitals);

  MoldenWriter(Logger& log) : _log(log){};

  ~MoldenWriter() = default;

 private:
  // clang-format off
  Logger& _log;
  std::array<Index,25> _multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,-1,-1, //f 
            1,1,1,1,1,-1,-1,-1,-1 //g
            };
  std::array<Index, 25> _reorderList={
            0, //s
            1,-1,0, //p
            0,1,-1,2,-2, //d
            0,1,-1,2,-2,3,-3, //f 
            0,1,-1,2,-2,3,-3,4,-4 //g
            };
  // clang-format on
  AOBasis _basis;
  BasisSet _bs;

  void writeAtoms(const Orbitals& orbitals, std::ofstream& outFile) const;
  void writeMOs(const Orbitals& orbitals, std::ofstream& outFile) const;
  void writeBasisSet(const Orbitals& orbitals, std::ofstream& outFile) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOLDEN_WRITER_H