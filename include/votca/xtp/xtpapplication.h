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

#include "votca/tools/propertyiomanipulator.h"
#ifndef VOTCA_XTP_XTPAPPLICATION_H
#define VOTCA_XTP_XTPAPPLICATION_H

// Third party includes
#include <algorithm>
#include <votca/tools/application.h>
#include <votca/tools/property.h>

namespace votca {
namespace xtp {

class XtpApplication : public votca::tools::Application {
 public:
  XtpApplication();
  ~XtpApplication() override;

  void Initialize() final;
  bool EvaluateOptions() final;
  void Run() final;
  void ShowHelpText(std::ostream &out) final;


 protected:
  virtual void execute() = 0;
  virtual std::string CalculatorType() const = 0;

  virtual void EvaluateSpecificOptions() = 0;

  virtual std::vector<std::string> CalculatorNames() const = 0;

  virtual void CreateCalculator(const std::string &name) = 0;

  virtual void AddCommandLineOptions() = 0;

  votca::tools::Property options_;

 private:
  void PrintShortHelp(std::ostream &out,
                      const std::string &calculator_name) const;

  void PrintLongHelp(std::ostream &out, const std::string &calculator_name,
                     tools::PropertyIOManipulator::Type format) const;

  bool CalcExists(const std::string &name) const {
    std::vector<std::string> names = CalculatorNames();
    return (std::find(names.begin(), names.end(), name) != names.end());
  }
};

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_XTPAPPLICATION_H
