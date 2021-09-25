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

#ifndef VOTCA_TOOLS_PROPERTYIOMANIPULATOR_H
#define VOTCA_TOOLS_PROPERTYIOMANIPULATOR_H

// Local VOTCA includes
#include "colors.h"
#include "property.h"

namespace votca {
namespace tools {

/**
 * \brief Manipulates the format state of the output stream
 *
 * Changes the state of the output stream. Property class reads this state
 * and formats its output according to this state (XML, TXT, etc)
 */
class PropertyIOManipulator {

 public:
  enum Type { XML, HLP, TXT };

  explicit PropertyIOManipulator(Type type = XML, Index level = 0,
                                 std::string indentation = "",
                                 ColorSchemeBase *color_scheme = nullptr)
      : type_(type),
        level_(level),
        indentation_(indentation),
        color_scheme_(color_scheme) {
    ;
  }

  ~PropertyIOManipulator() { delete color_scheme_; }
  friend std::ostream &operator<<(std::ostream &os,
                                  PropertyIOManipulator &piom) {
    os.pword(int(Property::getIOindex())) = &piom;
    return os;
  }

  const Type &getType() { return type_; }
  void setType(Type type) { type_ = type; }
  const Index &getLevel() { return level_; }
  void setLevel(Index level) { level_ = level; }
  const std::string &getIndentation() { return indentation_; }
  void setIndentation(std::string indentation) { indentation_ = indentation; }
  const ColorSchemeBase *getColorScheme() {
    if (!color_scheme_) {
      return &DEFAULT_COLORS;
    }
    return color_scheme_;
  }

  template <typename T>
  const ColorSchemeBase *setColorScheme() {
    if (color_scheme_) {
      delete color_scheme_;
    }
    color_scheme_ = new Color<T>();
    return color_scheme_;
  }

 private:
  Type type_;
  Index level_;
  std::string indentation_;
  ColorSchemeBase *color_scheme_;
};

extern PropertyIOManipulator XML;
extern PropertyIOManipulator HLP;
extern PropertyIOManipulator TXT;

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_PROPERTYIOMANIPULATOR_H
