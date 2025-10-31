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
#ifndef VOTCA_XTP_EANALYZE_H
#define VOTCA_XTP_EANALYZE_H

// Standard includes
#include <cmath>
#include <fstream>
#include <numeric>

// VOTCA includes
#include <votca/tools/histogram.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/qmstate.h"
#include "votca/xtp/topology.h"

namespace votca {
namespace xtp {

class EAnalyze final : public QMCalculator {
 public:
  EAnalyze() = default;

  ~EAnalyze() = default;
  bool WriteToStateFile() const { return false; }
  std::string Identify() const { return "eanalyze"; }

 protected:
  void ParseOptions(const tools::Property &user_options);
  bool Evaluate(Topology &top);

 private:
  void SiteHist(QMStateType state) const;
  void PairHist(const Topology &top, QMStateType state) const;
  void SiteCorr(const Topology &top, QMStateType state) const;

  double resolution_pairs_;
  double resolution_sites_;
  double resolution_spatial_;
  bool atomdistances_ = false;

  std::vector<QMStateType> states_;

  bool doenergy_landscape_;
  Index first_seg_;
  Index last_seg_;

  std::string seg_pattern_;
  std::vector<Segment *> seg_shortlist_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_EANALYZE_H
