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
#ifndef _VOTCA_XTP_EINTERNAL_H
#define _VOTCA_XTP_EINTERNAL_H

#include "votca/xtp/qmstate.h"
#include <votca/xtp/qmcalculator.h>

namespace votca {
namespace xtp {

class EInternal : public QMCalculator {
 public:
  EInternal(){};
  ~EInternal(){};

  std::string Identify() { return "einternal"; }
  bool WriteToStateFile() const { return true; }
  void Initialize(tools::Property &options);
  void ParseEnergies();
  bool EvaluateFrame(Topology &top);

 private:
  std::map<std::string, QMStateCarrierStorage<double> > _seg_U_xX_nN;
  std::map<std::string, QMStateCarrierStorage<double> > _seg_U_nX_nN;
  std::map<std::string, QMStateCarrierStorage<double> > _seg_U_xN_xX;

  std::map<std::string, QMStateCarrierStorage<bool> > _seg_has_state;

  std::map<std::string, bool> _has_seg;

  std::string _energiesXML;
};

}  // namespace xtp
}  // namespace votca

#endif  //_VOTCA_XTP_EINTERNAL_H
