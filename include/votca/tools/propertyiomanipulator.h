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

#ifndef _VOTCA_TOOLS_PROPERTY_IO_MANIPULATOR_H
#define _VOTCA_TOOLS_PROPERTY_IO_MANIPULATOR_H

#include <iostream>

#include "colors.h"
#include "property.h"

namespace votca {
namespace tools {

/**
 * \brief Manipulates the format state of the output stream
 *
 * Changes the state of the output stream. Property class reads this state
 * and formats its output according to this state (XML, TXT, T2T, etc)
 */
class PropertyIOManipulator {

 public:
  enum Type { XML, HLP, TEX, TXT };

  explicit PropertyIOManipulator(Type type = XML, int level = 0,
                                 std::string indentation = "",
                                 ColorSchemeBase *color_scheme = NULL)
      : _type(type),
        _level(level),
        _indentation(indentation),
        _color_scheme(color_scheme) {
    ;
  }

  ~PropertyIOManipulator() { delete _color_scheme; }
  friend std::ostream &operator<<(std::ostream &os,
                                  PropertyIOManipulator &piom) {
    os.pword(Property::getIOindex()) = &piom;
    return os;
  }

  const Type &getType() { return _type; }
  void setType(Type type) { _type = type; }
  const int &getLevel() { return _level; }
  void setLevel(int level) { _level = level; }
  const std::string &getIndentation() { return _indentation; }
  void setIndentation(std::string indentation) { _indentation = indentation; }
  const ColorSchemeBase *getColorScheme() {
    if (!_color_scheme) return &DEFAULT_COLORS;
    return _color_scheme;
  }

  template <typename T>
  const ColorSchemeBase *setColorScheme() {
    if (_color_scheme) delete _color_scheme;
    _color_scheme = new Color<T>();
    return _color_scheme;
  }

 private:
  Type _type;
  int _level;
  std::string _indentation;
  ColorSchemeBase *_color_scheme;
};

extern PropertyIOManipulator XML;
extern PropertyIOManipulator TXT;
extern PropertyIOManipulator TEX;
extern PropertyIOManipulator HLP;

}  // namespace tools
}  // namespace votca

#endif /* _VOTCA_TOOLS_PROPERTY_FORMAT_H */
