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
  StateApplication() = default;

  ~StateApplication() override = default;

 protected:

 virtual bool savetoStateFile()const =0;

  virtual void ConfigCalculator()=0;
  virtual bool EvaluateFrame(Topology& top)=0;

  void EvaluateSpecificOptions() final;
  virtual void CheckOptions() =0;
  void execute() final;
  void AddCommandLineOptions() final;
  virtual void AddCommandLineOpt()=0;

};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_STATEAPPLICATION_H
