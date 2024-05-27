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
#include "votca/tools/optionshandler.h"
#include "votca/tools/propertyiomanipulator.h"
#include "votca/tools/tokenizer.h"
#include <algorithm>
#include <stdexcept>
#include <string>

namespace votca {
namespace tools {

template <typename T>
static bool IsValidCast(const tools::Property &prop) {
  try {
    prop.as<T>();
    return true;
  } catch (const std::runtime_error &e) {
    return false;
  }
}

void OptionsHandler::ResolveLinks(Property &prop) const {

  if (prop.hasAttribute("link")) {
    Tokenizer tok(prop.getAttribute<std::string>("link"), " ,");
    for (std::string path : tok) {
      std::string relative_path = "subpackages/" + path;
      std::string file_path = defaults_path_ + relative_path;
      tools::Property package;
      package.LoadFromXML(file_path);
      const tools::Property &options = *(package.begin());
      for (Property::const_AttributeIterator attr = options.firstAttribute();
           attr != options.lastAttribute(); ++attr) {
        if (!prop.hasAttribute(attr->first)) {
          prop.setAttribute(attr->first, attr->second);
        }
      }

      for (const auto &child : options) {
        prop.add(child);
      }
    }
  }

  for (tools::Property &child : prop) {
    ResolveLinks(child);
  }
}

Property OptionsHandler::ProcessUserInput(const Property &user_input,
                                          const std::string &calcname) const {
  Property print = LoadDefaults(calcname);

  CheckUserInput(user_input, print);

  OverwriteDefaultsWithUserInput(user_input.get("options"),
                                 print.get("options"));
  RemoveOptional(print);
  CheckRequired(print);
  InjectDefaultsAsValues(print);
  RecursivelyCheckOptions(print);
  return print;
}

void OptionsHandler::CheckUserInput(const Property &user_input,
                                    const Property &defaults) const {
  if (user_input.hasAttribute("unchecked")) {
    return;
  } else {
    for (const auto &child : user_input) {
      if (defaults.exists(child.name())) {
        CheckUserInput(child, defaults.get(child.name()));
      } else {
        throw std::runtime_error("Votca has no option:" + child.path() + "." +
                                 child.name());
      }
    }
  }
}

void OptionsHandler::CheckRequired(const Property &options) const {
  for (const auto &child : options) {
    CheckRequired(child);
  }
  if (options.hasAttribute("default") &&
      options.getAttribute<std::string>("default") == "REQUIRED" &&
      !options.hasAttribute("injected")) {
    throw std::runtime_error("Please specify an input for:" + options.path() +
                             "." + options.name());
  }
}

void OptionsHandler::RemoveOptional(Property &options) const {
  options.deleteChildren([](const Property &p) {
    return p.hasAttribute("default") &&
           p.getAttribute<std::string>("default") == "OPTIONAL" &&
           !p.hasAttribute("injected");
  });
  for (auto &child : options) {
    RemoveOptional(child);
  }
}

Property OptionsHandler::CalculatorOptions(const std::string &calcname) const {
  Property print = LoadDefaults(calcname);
  InjectDefaultsAsValues(print);
  CleanAttributes(print, {"link"});
  return print;
}

void OptionsHandler::CleanAttributes(
    Property &options, const std::vector<std::string> &attributes) const {
  for (const auto &attribute : attributes) {
    if (options.hasAttribute(attribute)) {
      options.deleteAttribute(attribute);
    }
  }
  for (auto &child : options) {
    CleanAttributes(child, attributes);
  }
}

// load the xml description of the calculator (with defaults and test values)
Property OptionsHandler::LoadDefaults(const std::string &calculatorname) const {
  Property defaults_all;
  std::string defaults_file_path =
      defaults_path_ + "/" + calculatorname + ".xml";
  defaults_all.LoadFromXML(defaults_file_path);
  ResolveLinks(defaults_all);
  return defaults_all;
}

void OptionsHandler::InjectDefaultsAsValues(Property &defaults) const {
  for (Property &prop : defaults) {
    if (prop.HasChildren()) {
      InjectDefaultsAsValues(prop);
    } else if (prop.hasAttribute("default") && !prop.hasAttribute("injected")) {
      std::string value = prop.getAttribute<std::string>("default");
      if (std::none_of(
              reserved_keywords_.begin(), reserved_keywords_.end(),
              [value](const std::string &keyword) { return value == keyword; }))
        prop.value() = value;
      ;
    }
  }
}

void OptionsHandler::OverwriteDefaultsWithUserInput(const Property &user_input,
                                                    Property &defaults) const {
  // There are 3 distinct cases
  // a) normal option that can be copied over
  // b) a list attribute is discovered
  // c) an unchecked attriute is found,then the whole set is simply copied
  // over from the user_options, should be done afterwards, as an unchecked
  // session is not checked and simply copied over

  if (!defaults.hasAttribute("list")) {
    defaults.value() = user_input.value();
    defaults.setAttribute("injected", "true");
    for (auto &child : defaults) {
      if (user_input.exists(child.name())) {
        OverwriteDefaultsWithUserInput(user_input.get(child.name()), child);
      }
    }
  } else {
    defaults.setAttribute("injected", "true");
    std::map<std::string, Index> tags;
    for (const auto &child : defaults) {
      tags[child.name()]++;
    }

    for (const auto &tag : tags) {
      if (tag.second > 1) {
        throw std::runtime_error(
            "Developers: Each distinct tag in list should only appear once.");
      }
      std::vector<const tools::Property *> inputs =
          user_input.Select(tag.first);

      // if the input has no elements with that tag, remove all children from
      // default with that tag as well
      if (inputs.empty()) {
        defaults.deleteChildren(
            [&](const Property &prop) { return prop.name() == tag.first; });
      } else {
        // copy the element from defaults if the user_input has the element more
        // than once
        Property copy = defaults.get(tag.first);
        for (Index i = 1; i < Index(inputs.size()); i++) {
          defaults.add(copy);
        }
        std::vector<Property *> newdefault_elements =
            defaults.Select(tag.first);
        // for each element in both lists do the overwrite again
        for (Index i = 0; i < Index(inputs.size()); i++) {
          OverwriteDefaultsWithUserInput(*inputs[i], *newdefault_elements[i]);
        }
      }
    }
  }

  if (defaults.hasAttribute("unchecked")) {
    for (const auto &child : user_input) {
      defaults.add(child);
    }
  }
}

std::vector<std::string> OptionsHandler::GetPropertyChoices(const Property &p) {
  if (p.hasAttribute("choices")) {
    std::string att = p.getAttribute<std::string>("choices");
    std::size_t start_bracket = att.find('[');
    if (start_bracket != std::string::npos) {
      std::size_t end_bracket = att.find(']');
      att = att.substr(start_bracket + 1, end_bracket - start_bracket - 1);
    }
    return Tokenizer{att, " ,"}.ToVector();
  } else {
    return {};
  }
}

void OptionsHandler::RecursivelyCheckOptions(const Property &p) const {

  for (const Property &prop : p) {
    if (prop.HasChildren()) {
      RecursivelyCheckOptions(prop);
    } else {
      std::vector<std::string> choices = GetPropertyChoices(prop);
      if (choices.empty()) {
        continue;
      }
      if (!IsValidOption(prop, choices)) {
        std::ostringstream oss;
        oss << "\nThe input value for \"" << prop.name() << "\"";
        if (choices.size() == 1) {
          oss << " should be a \"" << choices.front() << "\"";
        } else {
          oss << " should be one of the following values: ";
          for (const std::string &c : choices) {
            oss << "\"" << c << "\"" << " ";
          }
        }
        oss << " But \"" << prop.value()
            << "\" cannot be converted into one.\n";
        throw std::runtime_error(oss.str());
      }
    }
  }
}

bool OptionsHandler::IsValidOption(
    const Property &prop, const std::vector<std::string> &choices) const {

  const std::string &head = choices.front();
  std::ostringstream oss;
  bool is_valid = true;
  if (head == "bool") {
    is_valid = IsValidCast<bool>(prop);
  } else if (head == "float") {
    is_valid = IsValidCast<double>(prop);
  } else if (head == "float+") {
    is_valid = IsValidCast<double>(prop) && (prop.as<double>() >= 0.0);
  } else if (head == "int") {
    is_valid = IsValidCast<Index>(prop);
  } else if (head == "int+") {
    is_valid = IsValidCast<Index>(prop) && (prop.as<Index>() >= 0);
  } else {
    std::string value = prop.as<std::string>();
    std::string att = prop.getAttribute<std::string>("choices");
    std::size_t start_bracket = att.find('[');
    if (start_bracket == std::string::npos) {
      // There is a single choice out of multiple default valid choices
      auto it = std::find(std::cbegin(choices), std::cend(choices), value);
      is_valid = (it != std::cend(choices));
    } else {
      // there are multiple valid choices
      Tokenizer tok{value, " ,"};
      for (const std::string &word : tok) {
        auto it = std::find(std::cbegin(choices), std::cend(choices), word);
        if (it == std::cend(choices)) {
          is_valid = false;
          break;
        }
      }
    }
  }
  if (std::find(additional_choices_.begin(), additional_choices_.end(),
                prop.as<std::string>()) != additional_choices_.end()) {
    is_valid = true;
  }
  return is_valid;
}

}  // namespace tools
}  // namespace votca
