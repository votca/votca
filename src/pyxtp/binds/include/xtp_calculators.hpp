/*
 * Copyright 2009-2021 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#ifndef XTP_CALCULATORS_H_
#define XTP_CALCULATORS_H_

#include "votca/tools/property.h"
#include "votca/xtp/calculatorfactory.h"
#include "votca/xtp/qmcalculator.h"
#include <memory>
#include <string>

using namespace votca;

namespace pyxtp {

int call_calculator(const std::string& name, int nThreads,
                    std::string xml_file);

class XTPCalculators {
 public:
  void Initialize(const std::string& name, int nThreads,
                  votca::tools::Property prop);

 private:
  std::unique_ptr<xtp::QMCalculator> _calculator;
};

}  // namespace pyxtp

#endif  // XTP_CALCULATORS_H_
