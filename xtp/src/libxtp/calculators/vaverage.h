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
#ifndef VOTCA_XTP_VAVERAGE_H
#define VOTCA_XTP_VAVERAGE_H

// Standard includes
#include <cstdio>
#include <string>
#include <vector>

// Local VOTCA includes
#include "votca/xtp/eigen.h"
#include "votca/xtp/logger.h"
#include "votca/xtp/qmcalculator.h"
#include "votca/xtp/rate_engine.h"

namespace votca {
namespace xtp {

class VAverage final : public QMCalculator {
 public:
  VAverage() = default;

  ~VAverage() = default;

  std::string Identify() const { return "vaverage"; }
  bool WriteToStateFile() const { return false; }

 protected:
  void ParseOptions(const tools::Property& user_options);
  bool Evaluate(Topology& top);

 private:
  Logger log_;
  std::string ratefile_;
  std::string occfile_;
  std::string outputfile_;

  std::vector<double> ReadOccfile(std::string filename) const;
  std::vector<Rate_Engine::PairRates> ReadRatefile(std::string filename) const;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_VAVERAGE_H
