/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

// Local VOTCA includes
#include "votca/tools/globals.h"
#include <stdexcept>

namespace votca {

Log::Level Log::current_level = Log::info;
namespace tools {

std::string GetVotcaShare() {
  char *votca_share = getenv("VOTCASHARE");
  if (votca_share == nullptr) {
    throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
  }
  return std::string(votca_share);
}

bool VotcaShareSet() { return (getenv("VOTCASHARE") != nullptr); }

std::string globals::url = "http://www.votca.org";
std::string globals::email = "devs@votca.org";

}  // namespace tools
}  // namespace votca
