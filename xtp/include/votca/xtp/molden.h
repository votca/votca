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
#ifndef VOTCA_XTP_MOLDEN_H
#define VOTCA_XTP_MOLDEN_H

// Local VOTCA includes
#include "votca/xtp/logger.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmtool.h"
#include <fstream>

namespace votca {
namespace xtp {
class Molden {
 public:
  Molden(Logger& log) : log_(log) {};

  ~Molden() = default;

  void WriteFile(const std::string& filename, const Orbitals& orbitals) const;

  void setBasissetInfo(const std::string& basisset_name,
                       const std::string& aux_basisset_name = "") {
    basisset_name_ = basisset_name;
    aux_basisset_name_ = aux_basisset_name;
  }

  void parseMoldenFile(const std::string& filename, Orbitals& orbitals) const;

 private:
  // clang-format off
  Logger&  log_;

  std::array<Index,49>  multipliers_={{
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            -1,1,1,1,1,1,-1, //f 
            -1,-1,1,1,1,1,1,-1,-1, //g
            -1,-1,-1,1,1,1,1,1,-1,-1,-1, //h
            -1,-1,-1,-1,1,1,1,1,1,-1,-1,-1,-1 //i
            }};
  std::array<Index, 49>  reorderList_={{
            0, //s
            1,-1,0, //p
            0,1,-1,2,-2, //d
            0,1,-1,2,-2,3,-3, //f 
            0,1,-1,2,-2,3,-3,4,-4, //g
            0,1,-1,2,-2,3,-3,4,-4,5,-5, //h
            0,1,-1,2,-2,3,-3,4,-4,5,-5,6,-6 //i
            }};
  // clang-format on
  std::string basisset_name_ = "";
  std::string aux_basisset_name_ = "";

  void writeAtoms(const Orbitals& orbitals, std::ofstream& outFile) const;
  void writeMOs(const Orbitals& orbitals, std::ofstream& outFile) const;
  void writeBasisSet(const Orbitals& orbitals, std::ofstream& outFile) const;

  std::string readAtoms(QMMolecule& mol, const std::string& units,
                        std::ifstream& input_file) const;
  std::string readMOs(Orbitals& orbitals, std::ifstream& input_file) const;
  void addBasissetInfo(Orbitals& orbitals) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_MOLDEN_H
