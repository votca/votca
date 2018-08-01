/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */
#ifndef __VOTCA_TOOLS_IDENTITY_H
#define __VOTCA_TOOLS_IDENTITY_H
#include <exception>
namespace votca {
namespace tools {

/**
    \brief Information about Identity

    The identity object is meant to provide functionality for storing the id of
    an object it primariy meant to be used in child classes and provides a more
    safety than other imlementations.
*/
template<typename T>
class Identity {
 private:
  T id_;
  bool id_set_;

 public:
  /// Constructor
  Identity() : id_set_(false) {}
  /// Constructor that takes initial id
  Identity(const T &id) : id_(id), id_set_(true) {};
  /// Gets the id returns error of the id has not been set
  const T getId() {
    return (id_set_ ? id_ : throw std::runtime_error("ID not set"));
  }
  /// Set the id
  void setId(T id) { _id_set_ = true; id_ = id; }
};
}
}

#endif
