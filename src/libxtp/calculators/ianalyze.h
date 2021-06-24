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
#ifndef VOTCA_XTP_IANALYZE_H
#define VOTCA_XTP_IANALYZE_H

// Local VOTCA includes
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/qmpair.h"
#include "votca/xtp/qmstate.h"

namespace votca {
namespace xtp {

class IAnalyze final : public QMCalculator {
 public:
  std::string Identify() const { return "ianalyze"; }
  bool WriteToStateFile() const { return false; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Evaluate(Topology &top);

 private:
  void IHist(Topology &top, QMStateType state);
  void IRdependence(Topology &top, QMStateType state);

  double resolution_logJ2_;
  std::vector<QMStateType> states_;
  double resolution_spatial_;
  std::vector<QMPair::PairType> pairtype_;
  bool do_pairtype_ = false;
  bool do_IRdependence_ = false;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_IANALYZE_H
