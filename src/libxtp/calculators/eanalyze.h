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
#ifndef _VOTCA_XTP_EANALYZE_H
#define _VOTCA_XTP_EANALYZE_H

#include <fstream>
#include <math.h>
#include <numeric>
#include <votca/tools/histogramnew.h>
#include <votca/tools/tokenizer.h>
#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/qmstate.h>
#include <votca/xtp/topology.h>

namespace votca {
namespace xtp {

class EAnalyze : public QMCalculator {
 public:
  EAnalyze() = default;
  ;

  ~EAnalyze() override = default;
  ;
  bool WriteToStateFile() const override { return false; }
  std::string Identify() override { return "eanalyze"; }

  void Initialize(tools::Property &opt) override;
  bool EvaluateFrame(Topology &top) override;

 private:
  void SiteHist(QMStateType state) const;
  void PairHist(const Topology &top, QMStateType state) const;
  void SiteCorr(const Topology &top, QMStateType state) const;

  double _resolution_pairs;
  double _resolution_sites;
  double _resolution_space;
  bool _atomdistances = false;

  std::vector<QMStateType> _states;

  bool _skip_corr = false;
  bool _skip_sites = false;
  bool _skip_pairs = false;

  bool _doenergy_landscape;
  int _first_seg;
  int _last_seg;

  std::string _seg_pattern;
  std::vector<Segment *> _seg_shortlist;
};

}  // namespace xtp
}  // namespace votca

#endif  // _VOTCA_XTP_EANALYZE_H
