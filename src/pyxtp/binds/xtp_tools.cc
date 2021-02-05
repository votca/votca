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

#include "xtp_bind_tools.h"
#include <iostream>

using namespace votca;

namespace pyxtp {

void call_tool(const std::string& name, Index nThreads, std::string xmlfile) {
  votca::tools::Property prop;
  prop.LoadFromXML(xmlfile);
  // Call calculator
  pyxtp::XTPTools tool;
  tool.Initialize(name, nThreads, prop);
}

void XTPTools::Initialize(const std::string& name, Index nThreads,
                          votca::tools::Property prop) {
  xtp::QMToolFactory factory;
  _tool = factory.Create(name);
  _tool->setnThreads(nThreads);
  _tool->Initialize(prop);
  _tool->Evaluate();
}

}  // namespace pyxtp