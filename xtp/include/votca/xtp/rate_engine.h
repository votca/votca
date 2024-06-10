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
#ifndef VOTCA_XTP_RATE_ENGINE_H
#define VOTCA_XTP_RATE_ENGINE_H

// Local VOTCA includes
#include "eigen.h"
#include "qmpair.h"
#include "qmstate.h"

namespace votca {
namespace xtp {

class Rate_Engine {

 public:
  struct PairRates {
    double rate12 = 0.0;
    double rate21 = 0.0;
  };

  Rate_Engine(double temperature, const Eigen::Vector3d& field)
      : temperature_(temperature), field_(field) {};

  PairRates Rate(const QMPair& pair, QMStateType carriertype) const;

  friend std::ostream& operator<<(std::ostream& out,
                                  const Rate_Engine& rate_engine);

 private:
  double Marcusrate(double Jeff2, double deltaG, double reorg) const;
  std::string ratetype_ = "marcus";
  double temperature_ = 0.0;                         // units:Hartree
  Eigen::Vector3d field_ = Eigen::Vector3d::Zero();  // units:Hartree/bohr
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_RATE_ENGINE_H
