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
#ifndef VOTCA_XTP_ORCA_H
#define VOTCA_XTP_ORCA_H

// Local VOTCA includes
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/qmpackage.h"

namespace votca {
namespace xtp {
/**
    \brief Wrapper for the ORCA program

    The ORCA class executes the ORCA package
    and extracts information from its log and io files

*/
class Orbitals;
class Orca : public QMPackage {
 public:
  std::string getPackageName() const override { return "orca"; }

  void Initialize(const tools::Property& options) override;

  bool WriteInputFile(const Orbitals& orbitals) override;

  bool WriteShellScript();

  bool Run() override;

  void CleanUp() override;

  bool CheckLogFile();

  bool ParseLogFile(Orbitals& orbitals) override;

  bool ParseMOsFile(Orbitals& orbitals) override;

  StaticSegment GetCharges() const override;

  Eigen::Matrix3d GetPolarizability() const override;

 protected:
  const std::array<Index, 25>& ShellMulitplier() const final {
    return _multipliers;
  }
  const std::array<Index, 25>& ShellReorder() const final {
    return _reorderList;
  }

 private:
  // clang-format off
  std::array<Index,25> _multipliers={
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            -1,1,1,1,1,1,-1, //f 
            -1,-1,1,1,1,1,1,-1,-1 //g
            };
  std::array<Index, 25> _reorderList={
            0, //s
            0, 1,-1, //p
            0,1,-1,2,-2, //d
            0,1,-1,2,-2,3,-3, //f 
            0,1,-1,2,-2,3,-3,4,-4 //g
            };

  // clang-format on
  std::string indent(const double& number);

  void WriteBasisset(const QMMolecule& qmatoms, std::string& bs_name,
                     std::string& el_file_name);
  void WriteCoordinates(std::ofstream& inp_file, const QMMolecule&);
  void WriteECP(std::ofstream& inp_file, const QMMolecule&);
  void WriteBackgroundCharges();

  void WriteChargeOption() override;
  template <class T>
  void GetCoordinates(T& mol, std::string& line,
                      std::ifstream& input_file) const;
  std::string WriteMethod() const;
  std::string CreateInputSection(const std::string& key) const;
  bool KeywordIsSingleLine(const std::string& key) const;
  std::string GetOrcaFunctionalName() const;

  std::unordered_map<std::string, std::string> _convergence_map{
      {"low", "LooseSCF"},
      {"normal", "StrongSCF"},
      {"tight", "TightSCF"},
      {"verytight", "VeryTightSCF"},
      {"none", ""}};
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_ORCA_H
