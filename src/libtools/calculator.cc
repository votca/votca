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

#include "../../include/votca/tools/calculator.h"

namespace votca {
namespace tools {

std::string Calculator::GetVotcaShare() {
  char *votca_share = getenv("VOTCASHARE");
  if (votca_share == nullptr) {
    throw std::runtime_error("VOTCASHARE not set, cannot open help files.");
  }
  return std::string(votca_share);
}

Property Calculator::LoadDefaults(const std::string package) {

  std::string calculator_name = Identify();
  // add default values if specified in VOTCASHARE
  std::string votca_share = Calculator::GetVotcaShare();

  // load the xml description of the calculator (with defaults and test values)
  std::string xmlFile = votca_share + std::string("/") + package +
                        std::string("/xml/") + calculator_name +
                        std::string(".xml");

  Property defaults_all;
  defaults_all.LoadFromXML(xmlFile);
  return defaults_all.get("options." + calculator_name);
}

void Calculator::UpdateWithUserOptions(Property &default_options,
                                       const Property &user_options) {

  // copy options from the object supplied by the Application
  std::string calculator_name = Identify();
  Property options_id = user_options.get("options." + calculator_name);

  // if a value is given override default values
  OverwriteDefaultsWithUserInput(options_id, default_options);
}

void Calculator::OverwriteDefaultsWithUserInput(const Property &p,
                                                Property &defaults) {

  // Go through everything that is defined in user option
  for (const Property &prop : p) {
    if (prop.HasChildren()) {
      if (defaults.exists(prop.name())) {
        OverwriteDefaultsWithUserInput(prop, defaults.get(prop.name()));
      } else {
        Property &new_prop = defaults.add(prop.name(), "");
        new_prop = prop;
      }
    } else if (defaults.exists(prop.name())) {
      defaults.set(prop.name(), prop.value());
    } else {
      defaults.add(prop.name(), prop.value());
    }
  }
}

std::vector<std::string> Calculator::GetPropertyChoices(const Property &p) {
  if (p.hasAttribute("choices")) {
    std::string att = p.getAttribute<std::string>("choices");
    Tokenizer tok{att, " ,"};
    return tok.ToVector();
  } else {
    return {""};
  }
}

void Calculator::RecursivelyCheckOptions(const Property &p) {
  for (const Property &prop : p) {
    if (prop.HasChildren()) {
      RecursivelyCheckOptions(prop);
    } else {
      std::vector<std::string> choices = Calculator::GetPropertyChoices(prop);
      const std::string &head = choices.front();
      if (head != "") {
        try {
          Calculator::CheckOption(prop, choices);
        } catch (const std::runtime_error &e) {
          std::ostringstream oss;
          oss << "\nInput option: " << prop.name() << "\n";
          oss << "Must be one of: ";
          for (const std::string c : choices) {
            oss << c << "\n";
          }
          std::cout << oss.str();
          throw e;
        }
      }
    }
  }
}

void Calculator::CheckOption(const Property &prop,
                             const std::vector<std::string> &choices) {
  const std::string &head = choices.front();
  if (head == "float") {
    prop.as<double>();
  } else if (head == "float+") {
    if (prop.as<double>() < 0.0) {
      throw std::runtime_error("");
    }
  } else if (head == "int") {
    prop.as<Index>();
  } else if (head == "int+") {
    if (prop.as<Index>() < 0) {
      throw std::runtime_error("");
    }
  } else {
    std::string value = prop.as<std::string>();
    auto it = std::find(std::cbegin(choices), std::cend(choices), value);
    if (it == std::cend(choices)) {
      throw std::runtime_error("");
    }
  }
}

}  // namespace tools
}  // namespace votca
