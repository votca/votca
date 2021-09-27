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
#ifndef VOTCA_XTP_KMCLIFETIME_H
#define VOTCA_XTP_KMCLIFETIME_H

// Local VOTCA includes
#include "votca/xtp/kmccalculator.h"

namespace votca {
namespace xtp {

class KMCLifetime final : public KMCCalculator {
 public:
  KMCLifetime() = default;
  ~KMCLifetime() = default;
  bool WriteToStateFile() const { return false; }
  std::string Identify() const { return "kmclifetime"; }

 protected:
  void ParseSpecificOptions(const tools::Property& user_options);
  bool Evaluate(Topology& top);

 private:
  void WriteDecayProbability(std::string filename);

  void RunVSSM();
  void WriteToTraj(std::fstream& traj, unsigned long insertioncount,
                   double simtime, const Chargecarrier& affectedcarrier) const;

  void ReadLifetimeFile(std::string filename);
  std::string probfile_;
  bool do_carrierenergy_;
  std::string energy_outputfile_;
  double alpha_;
  unsigned long outputsteps_;
  unsigned long insertions_;
  std::string lifetimefile_;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCLIFETIME_H
