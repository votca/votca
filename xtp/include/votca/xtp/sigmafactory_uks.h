/*
 *            Copyright 2009-2023 The VOTCA Development Team
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

#ifndef VOTCA_XTP_SIGMAFACTORY_UKS_H
#define VOTCA_XTP_SIGMAFACTORY_UKS_H

#include <memory>
#include <string>

#include <votca/tools/objectfactory.h>

#include "votca/xtp/rpa_uks.h"
#include "votca/xtp/sigma_base_uks.h"
#include "votca/xtp/threecenter.h"

namespace votca {
namespace xtp {

class SigmaFactory_UKS
    : public tools::ObjectFactory<std::string, Sigma_base_UKS,
                                  TCMatrix_gwbse_spin&, RPA_UKS&,
                                  TCMatrix::SpinChannel> {
 private:
  void RegisterAll();

 public:
  SigmaFactory_UKS() { RegisterAll(); }
};

}  // namespace xtp
}  // namespace votca

#endif