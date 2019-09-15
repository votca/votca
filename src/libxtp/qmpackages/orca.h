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
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef __VOTCA_XTP_ORCA_H
#define __VOTCA_XTP_ORCA_H

#include <votca/xtp/qmpackage.h>

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
  std::string getPackageName() const { return "orca"; }

  void Initialize(tools::Property& options);

  bool WriteInputFile(const Orbitals& orbitals);

  bool WriteShellScript();

  bool Run();

  void CleanUp();

  bool CheckLogFile();

  bool ParseLogFile(Orbitals& orbitals);

  bool ParseMOsFile(Orbitals& orbitals);

  StaticSegment GetCharges() const;

  Eigen::Matrix3d GetPolarizability() const;

 private:
  std::string indent(const double& number);
  std::string getLName(int lnum);

  void WriteBasisset(const QMMolecule& qmatoms, std::string& _bs_name,
                     std::string& el_file_name);
  void WriteCoordinates(std::ofstream& com_file, const QMMolecule&);
  void WriteECP(std::ofstream& com_file, const QMMolecule&);
  void WriteBackgroundCharges();

  void WriteChargeOption();
  template <class T>
  void GetCoordinates(T& mol, std::string& line,
                      std::ifstream& input_file) const;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_XTP_ORCA_H */
