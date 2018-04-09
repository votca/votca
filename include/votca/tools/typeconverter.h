/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __TOOLS_TYPECONVERTER__H
#define __TOOLS_TYPECONVERTER__H

#include <string>
#include <unordered_map>
#include <vector>

namespace votca {
namespace tools {

/**
    \brief Provides a means to reduce cross dependency of header files
 *
 * This templace class allows the conversion of one type of object to another
 * this is done by breaking the said object up into the standard types.
 *
 * For the type converter to work each object it is converting between must
 * have defined at least one write function and or one read function.
 *
 * This tempalte object only allows converting between doubles,strings,ints.
 */
template <class T, class G>
class TypeConverter {
 private:
  std::unordered_map<std::string, std::vector<int>>         integers;
  std::unordered_map<std::string, std::vector<double>>      doubles;
  std::unordered_map<std::string, std::vector<std::string>> strings;

  /// Generic functions that must be defined in the T class
  std::unordered_map<std::string, std::vector<int>> (*writeInt_)(T t);
  std::unordered_map<std::string, std::vector<double>> (*writeDouble_)(T t);
  std::unordered_map<std::string, std::vector<std::string>> (*writeStr_)(T t);

  /// Generic function pointers that must be defined in the G class
  void (*readData_)(std::unordered_map<std::string, std::vector<int>>,
                    std::unordered_map<std::string, std::vector<double>>,
                    std::unordered_map<std::string, std::vector<std::string>>,
                    G& g);

 public:
  /// Constructor
  TypeConverter() {
    writeInt_    = &(T::writeInt);
    writeDouble_ = &T::writeDouble;
    writeStr_    = &T::writeStr;
    readData_    = &G::readData;
  }

  /// Imports data from type T
  void importData(T temp) {
    this->integers = writeInt_(temp);
    this->doubles  = writeDouble_(temp);
    this->strings  = writeStr_(temp);
  }

  /// Exports the data from type T and places it in type G
  void exportData(G& g) {
    return readData_(this->integers, this->doubles, this->strings, g);
  }
};
}
}
#endif  // # __TOOLS_TYPECONVERTER__H_
