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
#ifndef VOTCA_XTP_IANALYZE_H
#define VOTCA_XTP_IANALYZE_H
#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/qmpair.h>
#include <votca/xtp/qmstate.h>

namespace votca {
namespace xtp {

class IAnalyze : public QMCalculator {
 public:
  std::string Identify() { return "ianalyze"; }

  void Initialize(tools::Property &options);
  bool EvaluateFrame(Topology &top);
  void IHist(Topology &top, QMStateType state);
  void IRdependence(Topology &top, QMStateType state);

 private:
  double _resolution_logJ2;
  std::vector<QMStateType> _states;
  double _resolution_space;
  std::vector<QMPair::PairType> _pairtype;
  bool _do_pairtype;
  bool _do_IRdependence;
};

}  // namespace xtp
}  // namespace votca

#endif
