/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#ifndef PYXTP_H_
#define PYXTP_H_

#include "factory.hpp"
#include "parent.hpp"
#include "votca/xtp/calculatorfactory.h"
#include "votca/xtp/qmcalculator.h"
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace votca;

namespace pyxtp {

void printer(const std::vector<std::string>& keys);

int call_calculator(const std::string& name, int nThreads);

class PyXTP {
 public:
  void Initialize(const std::string& name, int nThreads);

 private:
  std::unique_ptr<xtp::QMCalculator> _calculator;
};

}  // namespace pyxtp

#endif  // PYXTP_H_
