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
#ifndef _VOTCA_XTP_VAVERAGE_H
#define _VOTCA_XTP_VAVERAGE_H

#include <stdio.h>

#include <votca/xtp/eigen.h>
#include <votca/xtp/logger.h>
#include <votca/xtp/qmcalculator.h>
#include <votca/xtp/rate_engine.h>

namespace votca {
namespace xtp {

class VAverage : public QMCalculator {
 public:
  VAverage(){};

  ~VAverage() override{};

  std::string Identify() override { return "vaverage"; }
  bool WriteToStateFile() const override { return false; }
  void Initialize(tools::Property& options) override;
  bool EvaluateFrame(Topology& top) override;

 private:
  Logger _log;
  std::string _ratefile;
  std::string _occfile;
  std::string _outputfile;

  std::vector<double> ReadOccfile(std::string filename) const;
  std::vector<Rate_Engine::PairRates> ReadRatefile(std::string filename) const;
};

}  // namespace xtp
}  // namespace votca

#endif
