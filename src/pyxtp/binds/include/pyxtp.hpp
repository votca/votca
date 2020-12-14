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

#if !defined(PYXTP_H_)
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

void printer(const std::vector<std::string>& keys) {
  for (const auto& key : keys) {
    std::cout << "key: " << key << "\n";
  }
  std::cout << "size: " << keys.size() << "\n";
}

class PyXTP {
 public:
  void Initialize(const std::string& name, int nThreads) {
    std::cout << "Votca Factory:\n";
    xtp::Calculatorfactory::RegisterAll();
    std::vector<std::string> keys = xtp::Calculators().getKeys();
    printer(keys);
    std::cout << "Mock Factory:\n";
    pyxtp::Factory<Parent>::RegisterAll();
    keys = pyxtp::Calculators<Parent>().getKeys();
    printer(keys);
  }

 private:
  std::unique_ptr<xtp::QMCalculator> _calculator;
};

}  // namespace pyxtp 

#endif  // PYXTP_H_
