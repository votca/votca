/*
 * Copyright 2009-2023 The VOTCA Development Team (http://www.votca.org)
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
#include <votca/tools/optionshandler.h>

using namespace votca;

namespace pyxtp {

void call_tool(const std::string& name, Index nThreads, std::string xmlfile) {
  votca::tools::Property prop;
  prop.LoadFromXML(xmlfile);

  votca::tools::OptionsHandler handler(tools::GetVotcaShare() + "/xtp/xml/");
  votca::tools::Property options =
      handler.ProcessUserInput(prop, name).get("options." + name);

  // Call calculator
  pyxtp::XTPTools tool;
  tool.Initialize(name, nThreads, options);
}

void XTPTools::Initialize(const std::string& name, Index nThreads,
                          votca::tools::Property prop) {
  tool_ = xtp::QMToolFactory().Create(name);
  tool_->setnThreads(nThreads);
  tool_->Initialize(prop);
  tool_->Evaluate();
}

}  // namespace pyxtp