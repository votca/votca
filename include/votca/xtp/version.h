/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#ifndef VOTCA_MD2QM_VERSION_H
#define VOTCA_MD2QM_VERSION_H

#include <string>
/**
 * \namespace votca::xtp
 * \brief Charge transport classes
 *
 * Classes used for charge and exciton transport simulations
 */
namespace votca {
namespace xtp {
const std::string &XtpVersionStr();
void HelpTextHeader(const std::string &tool_name);
}  // namespace xtp
}  // namespace votca

#endif  //  VOTCA_MD2QM_VERSION_H
