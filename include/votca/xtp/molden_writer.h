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

/**
 * \brief Writes an orbital file to a .cube file
 */

namespace votca {
namespace xtp {
class Molden_Writer {
  Molden_Writer(Logger& log) : _log(log){};

 public:
  void WriteFile(const std::string& filename, const Orbitals& orbitals) const;

 private:  // clang-format off
  Logger& _log;
  std::array<Index,25> _multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,-1,-1, //f 
            1,1,1,1,1,-1,-1,-1,-1 //g
            };

  OrbTranspositions _transpositions { 
    std::vector<std::array<Index, 2>> {}, //s
    std::vector<std::array<Index, 2>> {   //p
      {0, 2}
    }, 
    std::vector<std::array<Index, 2>> {   //d
      {1, 2},
      {3, 4}
      }, 
    std::vector<std::array<Index, 2>> {   //f
      {1, 2},  
      {3, 4},
      {5, 6}
    }, 
    std::vector<std::array<Index, 2>> {   //g
      {1, 2},
      {3, 4},
      {5, 6},
      {7, 8}
    }
  };
  // clang-format on
  std::string _basisset_name;
  std::string _aux_basisset_name;
  AOBasis _basis;
  BasisSet _bs;
  Logger& _log;

  void writeAtoms(const Orbitals& orbitals, std::ofstream& outFile);
  void writeMOs(const Orbitals& orbitals, std::ofstream& outFile);
  void writeBasisSet(const Orbitals& orbitals, std::ofstream& outFile);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOLDEN_WRITER_H