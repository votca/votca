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
#ifndef VOTCA_XTP_XTPDFT_H
#define VOTCA_XTP_XTPDFT_H

// Standard includes
#include <string>

// Local VOTCA includes
#include "votca/xtp/dftengine.h"
#include "votca/xtp/orbitals.h"
#include "votca/xtp/orbreorder.h"
#include "votca/xtp/polarsite.h"
#include "votca/xtp/qmpackage.h"

namespace votca {
namespace xtp {

/**
    \brief Wrapper for the internal XTP DFT engine


 */

class XTPDFT : public QMPackage {
 public:
  std::string getPackageName() const final { return "xtp"; }

  void Initialize(const tools::Property& options) final;

  bool WriteInputFile(const Orbitals& orbitals) final;

  bool Run() final;

  void CleanUp() final;

  bool CheckLogFile();

  bool ParseLogFile(Orbitals& orbitals) final;

  bool ParseMOsFile(Orbitals& orbitals) final;

  StaticSegment GetCharges() const final {
    throw std::runtime_error(
        "If you want partial charges just run the 'partialcharges' calculator");
  }

  Eigen::Matrix3d GetPolarizability() const final {
    throw std::runtime_error(
        "GetPolarizability() is not implemented for xtpdft");
  }

 protected:
  const std::array<Index, 49>& ShellMulitplier() const final {
    return _multipliers;
  }
  const std::array<Index, 49>& ShellReorder() const final {
    return _reorderList;
  }

 private:
  // clang-format off
  std::array<Index,49> _multipliers={{
            1, //s
            1,1,1, //p
            1,1,1,1,1, //d
            1,1,1,1,1,1,1, //f 
            1,1,1,1,1,1,1,1,1, //g
            1,1,1,1,1,1,1,1,1,1,1, //h
            1,1,1,1,1,1,1,1,1,1,1,1,1 //i
  }};
  std::array<Index,49> _reorderList={{
            0, //s
            0,0,0, //p
            0,0,0,0,0, //d
            0,0,0,0,0,0,0, //f 
            0,0,0,0,0,0,0,0,0, //g
            0,0,0,0,0,0,0,0,0,0,0, //h
            0,0,0,0,0,0,0,0,0,0,0,0,0, //i
            }};
  // clang-format on

  void WriteChargeOption() final { return; }
  tools::Property _xtpdft_options;

  Orbitals _orbitals;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_XTPDFT_H
