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
#ifndef VOTCA_XTP_QMCALCULATOR_H
#define VOTCA_XTP_QMCALCULATOR_H

// VOTCA includes
#include <votca/tools/calculator.h>

namespace votca {
namespace xtp {

class Topology;

class QMCalculator : public tools::Calculator {
 public:
  QMCalculator() = default;
  ~QMCalculator() override = default;

  std::string Identify() const override = 0;

  std::string Package() const final { return "xtp"; }

  virtual bool WriteToStateFile() const = 0;

  bool EvaluateFrame(Topology& top);

  void Initialize(const tools::Property& opt) final;

 protected:
  virtual void ParseOptions(const tools::Property& opt) = 0;
  virtual bool Evaluate(Topology& top) = 0;
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_QMCALCULATOR_H
