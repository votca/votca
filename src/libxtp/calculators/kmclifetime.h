/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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
  ~KMCLifetime() {
    for (auto &node : _nodes) {
      delete node;
    }
    for (auto &carrier : _carriers) {
      delete carrier;
    }
  };
  std::string Identify() { return "kmclifetime"; }
  void Initialize(tools::Property *options);
  bool EvaluateFrame(ctp::Topology *top);

 private:
  void WriteDecayProbability(string filename);

  void RunVSSM(ctp::Topology *top);

  void ReadLifetimeFile(string filename);

  // tools::vec _field;
  string _probfile;
  bool _do_carrierenergy;
  string _energy_outputfile;
  double _alpha;
  unsigned _outputsteps;
  unsigned int _insertions;
  std::string _lifetimefile;
  double _maxrealtime;
  string _trajectoryfile;
  string _outputfile;
  string _filename;
};

}  // namespace xtp
}  // namespace votca

#endif /* __VOTCA_KMC_MULTIPLE_H */
