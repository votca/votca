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
 * author: Kordt
 */

#ifndef __VOTCA_KMC_LIFETIME_H
#define __VOTCA_KMC_LIFETIME_H

#include <votca/xtp/kmccalculator.h>
using namespace std;

namespace votca {
namespace xtp {

class KMCLifetime : public KMCCalculator {
 public:
  KMCLifetime(){};
  ~KMCLifetime(){};
  std::string Identify() { return "kmclifetime"; }
  void Initialize(tools::Property &options);
  bool EvaluateFrame(Topology &top);

 private:
  void WriteDecayProbability(std::string filename);

  void RunVSSM(Topology &top);

  void ReadLifetimeFile(std::string filename);
  std::string _probfile;
  bool _do_carrierenergy;
  std::string _energy_outputfile;
  double _alpha;
  unsigned _outputsteps;
  unsigned int _insertions;
  std::string _lifetimefile;
  double _maxrealtime;
  std::string _trajectoryfile;
  std::string _outputfile;
  std::string _filename;
  std::string _occfile;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_KMC_MULTIPLE_H */
