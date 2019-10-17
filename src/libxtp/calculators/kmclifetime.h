/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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
#ifndef __VOTCA_KMC_LIFETIME_H
#define __VOTCA_KMC_LIFETIME_H

#include <votca/xtp/kmccalculator.h>
using namespace std;

namespace votca {
namespace xtp {

class KMCLifetime : public KMCCalculator {
 public:
  KMCLifetime() = default;
  ;
  ~KMCLifetime() override = default;
  ;
  bool WriteToStateFile() const override { return false; }
  std::string Identify() override { return "kmclifetime"; }
  void Initialize(tools::Property& options) override;
  bool EvaluateFrame(Topology& top) override;

 private:
  void WriteDecayProbability(std::string filename);

  void RunVSSM() override;
  void WriteToTraj(fstream& traj, unsigned insertioncount, double simtime,
                   const Chargecarrier& carrier) const;

  void ReadLifetimeFile(std::string filename);
  std::string _probfile;
  bool _do_carrierenergy;
  std::string _energy_outputfile;
  double _alpha;
  unsigned _outputsteps;
  unsigned int _insertions;
  std::string _lifetimefile;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_KMC_MULTIPLE_H */
