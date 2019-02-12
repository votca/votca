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

#ifndef VOTCA_CSG_ANALYSISTOOL_H
#define VOTCA_CSG_ANALYSISTOOL_H

#include "bondedstatistics.h"
#include <map>
#include <string>
#include <votca/csg/cgengine.h>

/**
    \brief base class for all analasys tools

    This is the base class for all analasys tool.
    \todo do option functions!!!
*/
namespace votca {
namespace csg {

class AnalysisTool {
public:
  AnalysisTool() {}
  virtual ~AnalysisTool() {}

  virtual void Register(std::map<std::string, AnalysisTool *> &lib) {}
  virtual void Command(BondedStatistics &bs, std::string cmd,
                       std::vector<std::string> &args){};
  virtual void Help(std::string cmd, std::vector<std::string> &args){};
};

} // namespace csg
} // namespace votca
#endif // VOTCA_CSG_ANALYSISTOOL_H
