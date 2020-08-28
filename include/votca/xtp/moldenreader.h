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
#ifndef VOTCA_XTP_MOLDEN_READER_H
#define VOTCA_XTP_MOLDEN_READER_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>
#include <votca/xtp/orbitals.h>

namespace votca {
namespace xtp {
class MoldenReader {
 public:
  MoldenReader(Logger& log) : _log(log){};

  ~MoldenReader() = default;

  void setBasissetInfo(std::string basisset_name,
                       std::string aux_basisset_name) {
    _basisset_name = basisset_name;
    _aux_basisset_name = aux_basisset_name;
  }

  void parseMoldenFile(const std::string& filename, Orbitals& orbitals);

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
  std::string _basisset_name;
  std::string _aux_basisset_name;
  AOBasis _basis;

  std::string readAtoms(QMMolecule& mol, const std::string& units,
                        std::ifstream& input_file) const;
  std::string readMOs(Orbitals& orbitals, std::ifstream& input_file) const;
  void addBasissetInfo(Orbitals& orbitals);
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOLDEN_READER_H