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
#ifndef VOTCA_XTP_CALCULATORFACTORY_H
#define VOTCA_XTP_CALCULATORFACTORY_H

// Standard includes
#include <map>

// VOTCA includes
#include <votca/tools/objectfactory.h>

// Local VOTCA includes
#include "qmcalculator.h"

namespace votca {
namespace xtp {

class Calculatorfactory
    : public tools::ObjectFactory<std::string, QMCalculator> {
 private:
  Calculatorfactory() = default;

 public:
  static void RegisterAll(void);

  friend Calculatorfactory &Calculators();
};

inline Calculatorfactory &Calculators() {
  static Calculatorfactory instance;
  return instance;
}

}  // namespace xtp
}  // namespace votca

#endif  // VOTCA_XTP_CALCULATORFACTORY_H
