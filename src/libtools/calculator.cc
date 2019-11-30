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

#include <votca/tools/calculator.h>

namespace votca {
namespace tools {

void Calculator::UpdateWithDefaults(Property &options, std::string package) {

  // copy options from the object supplied by the Application
  std::string id = Identify();
  Property options_id = options.get("options." + id);

  // add default values if specified in VOTCASHARE
  char *votca_share = getenv("VOTCASHARE");
  if (votca_share == nullptr) {
    throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
  }
  // load the xml description of the calculator (with defaults and test values)
  std::string xmlFile = std::string(getenv("VOTCASHARE")) + std::string("/") +
                        package + std::string("/xml/") + id +
                        std::string(".xml");

  Property defaults_all;
  defaults_all.LoadFromXML(xmlFile);
  Property defaults = defaults_all.get("options." + id);

  // if a value not given or a tag not present, provide default values
  AddDefaults(options_id, defaults);

  // output calculator options
  std::string indent("          ");
  Index level = 1;
  votca::tools::PropertyIOManipulator IndentedText(PropertyIOManipulator::TXT,
                                                   level, indent);
  if (Log::verbose()) {
    std::cout << "\n... ... options\n"
              << IndentedText << options_id << "... ... options\n"
              << std::flush;
  }
}

void Calculator::AddDefaults(Property &p, const Property &defaults) {

  for (const Property &prop : defaults) {
    std::string name = prop.path() + "." + prop.name();

    Property rootp = *p.begin();
    if (prop.hasAttribute("default")) {
      if (rootp.exists(name)) {
        if (rootp.HasChildren()) {
          rootp.value() = prop.value();
        }
      } else {
        rootp.add(prop.name(), prop.value());
      }
    }
    AddDefaults(p, prop);
  }
}

}  // namespace tools
}  // namespace votca
