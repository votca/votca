/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#pragma once
#ifndef VOTCA_XTP_KMCCALCULATOR_H
#define VOTCA_XTP_KMCCALCULATOR_H

// VOTCA includes
#include <votca/tools/globals.h>
#include <votca/tools/random.h>
#include <votca/tools/tokenizer.h>

// Local VOTCA includes
#include "chargecarrier.h"
#include "gnode.h"
#include "logger.h"
#include "qmcalculator.h"
#include "qmstate.h"

namespace votca {
namespace xtp {
class QMNBList;
class KMCCalculator : public QMCalculator {
 public:
  ~KMCCalculator() override = default;

  void ParseOptions(const tools::Property& options) final {
    ParseCommonOptions(options);
    ParseSpecificOptions(options);
  }

 protected:
  virtual void ParseSpecificOptions(const tools::Property& options) = 0;

  QMStateType _carriertype;

  void LoadGraph(Topology& top);
  virtual void RunVSSM() = 0;

  void ParseCommonOptions(const tools::Property& options);

  double Promotetime(double cumulated_rate);
  void ResetForbiddenlist(std::vector<GNode*>& forbiddenid) const;
  void AddtoForbiddenlist(GNode& node, std::vector<GNode*>& forbiddenid) const;
  bool CheckForbidden(const GNode& node,
                      const std::vector<GNode*>& forbiddenlist) const;
  bool CheckSurrounded(const GNode& node,
                       const std::vector<GNode*>& forbiddendests) const;
  const GLink& ChooseHoppingDest(const GNode& node);
  Chargecarrier* ChooseAffectedCarrier(double cumulated_rate);

  void WriteOccupationtoFile(double simtime, std::string filename) const;
  void WriteRatestoFile(std::string filename, const QMNBList& nblist) const;

  void RandomlyCreateCharges();
  void RandomlyAssignCarriertoSite(Chargecarrier& Charge);
  std::vector<GNode> _nodes;
  std::vector<Chargecarrier> _carriers;

  tools::Random _RandomVariable;
  std::string _injection_name;
  std::string _injectionmethod;
  Index _seed;
  Index _numberofcarriers;
  Eigen::Vector3d _field = Eigen::Vector3d::Zero();
  double _maxrealtime;
  std::string _trajectoryfile;
  std::string _ratefile;
  std::string _occfile;

  Logger _log;

  double _temperature;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCCALCULATOR_H
