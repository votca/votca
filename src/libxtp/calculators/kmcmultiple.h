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
#ifndef VOTCA_XTP_KMCMULTIPLE_H
#define VOTCA_XTP_KMCMULTIPLE_H

// Standard includes
#include <fstream>

// Local VOTCA includes
#include "votca/xtp/kmccalculator.h"

namespace votca {
namespace xtp {

class KMCMultiple final : public KMCCalculator {
 public:
  KMCMultiple() = default;
  ~KMCMultiple() = default;
  bool WriteToStateFile() const { return false; }
  std::string Identify() const { return "kmcmultiple"; }

 protected:
  void ParseSpecificOptions(const tools::Property& user_options);
  bool Evaluate(Topology& top);

 private:
  void RunVSSM();
  void PrintChargeVelocity(double simtime);

  void PrintDiagDandMu(const Eigen::Matrix3d& avgdiffusiontensor,
                       double simtime, unsigned long step);

  void WriteToEnergyFile(std::fstream& tfile, double simtime,
                         unsigned long step) const;

  void WriteToTrajectory(std::fstream& traj,
                         std::vector<Eigen::Vector3d>& startposition,
                         double simtime, unsigned long step) const;

  void PrintDiffandMu(const Eigen::Matrix3d& avgdiffusiontensor, double simtime,
                      unsigned long step);

  double runtime_;
  double outputtime_;
  std::string timefile_ = "";
  Index intermediateoutput_frequency_ = 10000;
  unsigned long diffusionresolution_ = 1000;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_KMCMULTIPLE_H
