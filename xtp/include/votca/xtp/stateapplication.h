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
#ifndef VOTCA_XTP_STATEAPPLICATION_H
#define VOTCA_XTP_STATEAPPLICATION_H

// Local VOTCA includes
#include "qmcalculator.h"
#include "topology.h"
#include "xtpapplication.h"

namespace votca {
namespace xtp {

class StateApplication : public XtpApplication {
 public:
  StateApplication();

  ~StateApplication() override = default;

  void Initialize() override;
  bool EvaluateOptions() override;
  void Run() override;

  void BeginEvaluate(Index nThreads);
  bool EvaluateFrame(Topology& top);

  void SetCalculator(std::unique_ptr<QMCalculator>&& calculator);

 protected:
  std::unique_ptr<QMCalculator> _calculator;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATEAPPLICATION_H
