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
#ifndef VOTCA_XTP_EINTERNAL_H
#define VOTCA_XTP_EINTERNAL_H

// Local VOTCA includes
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/qmstate.h"

namespace votca {
namespace xtp {

class EInternal final : public QMCalculator {
 public:
  EInternal() = default;
  ~EInternal() = default;

  std::string Identify() const { return "einternal"; }

  bool WriteToStateFile() const { return true; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Evaluate(Topology &top);

 private:
  void ParseEnergies();

  std::map<std::string, QMStateCarrierStorage<double> > seg_U_xX_nN_;
  std::map<std::string, QMStateCarrierStorage<double> > seg_U_nX_nN_;
  std::map<std::string, QMStateCarrierStorage<double> > seg_U_xN_xX_;

  std::map<std::string, QMStateCarrierStorage<bool> > seg_has_state_;

  std::map<std::string, bool> has_seg_;

  std::string energiesXML_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EINTERNAL_H
