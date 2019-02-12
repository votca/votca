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

#ifndef _VOTCA_TOOLS_PROPERTY_H
#define _VOTCA_TOOLS_PROPERTY_H

#include "lexical_cast.h"
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <iostream>
#include <list>
#include <map>
#include <stdexcept>
#include <stdlib.h>
#include <string>

#include "tokenizer.h"
#include "vec.h"

namespace votca {
namespace tools {

/**
 * \brief class to manage program options with xml serialization functionality
 *
 * This class stores tags and content in a hierarchical (tree) structure similar
 * to the one used in the XML format. The structure can be either filled
 * manually or read in from an XML file using load_property_from_xml. The
 * supported XML constructs are TAGS, ATTRIBUTES, and CONTENT: <tag
 * attribute_name="attribute_value"> content
 * </tag>
 * The property object can be output to an ostream using format modifiers:
 * cout << XML << property;
 * Supported formats are XML, TXT, TEX, HLP
 */
class Property {

  /// \brief outputs the property to the ostream
  friend std::ostream &operator<<(std::ostream &out, Property &p);

 public:
  Property() : _path("") {}

  Property(const std::string &name, const std::string &value,
           const std::string &path)
      : _name(name), _value(value), _path(path) {}

  /**
   * \brief add a new property to structure
   * @param key identifier
   * @param value value
   * @return reference to the created Property object
   */
  Property &add(const std::string &key, const std::string &value);
  /**
   * \brief set value of existing property
   * @param key identifier
   * @param value value
   * @return reference to the created Property object
   */
  Property &set(const std::string &key, const std::string &value);

  /**
   * \brief get existing property
   * @param key identifier
   * @return Reference to property object
   *
   * This function tries to find a property specified by key separated
   * by "." to step down hierarchy. If the property is not
   * found a runtime_exception is thrown.
   */
  Property &get(const std::string &key);

  /**
   * \brief check whether property exists
   * @param key identifier
   * @return true or false
   */
  bool exists(const std::string &key);

  /**
   * \brief select property based on a filter
   * @param filter
   * @return list of pointers to property objects
   *
   * returns a list of properties that match the key criteria including
   * wildcards "*" and "?". Example: "base.item*.value"
   */
  std::list<Property *> Select(const std::string &filter);

  /**
   * \brief reference to value of property
   * @return std::string content
   */
  std::string &value() { return _value; }
  /**
   * \brief name of property
   * @return name
   */
  std::string name() { return _name; }
  /**
   * \brief full path of property (including parents)
   * @return path
   *
   * e.g. cg.inverse.value
   */
  std::string path() { return _path; }
  /**
   * \brief return value as type
   *
   * returns the value after type conversion, e.g.
   * p.as<int>() returns an integer
   */
  template <typename T>
  T as() const;

  template <typename T>
  T ifExistsReturnElseReturnDefault(const std::string &key, T defaultvalue);

  template <typename T>
  T ifExistsReturnElseThrowRuntimeError(const std::string &key);

  template <typename T>
  T ifExistsAndinListReturnElseThrowRuntimeError(
      const std::string &key, std::vector<T> possibleReturns);

  /**
   * \brief does the property has childs?
   * \return true or false
   */
  bool HasChilds() { return !_map.empty(); }

  /// iterator to iterate over properties
  typedef std::list<Property>::iterator iterator;
  /// \brief iterator to first child property
  iterator begin() { return _properties.begin(); }
  /// \brief end iterator for child properties
  iterator end() { return _properties.end(); }
  /// \brief number of child properties
  std::list<Property>::size_type size() { return _properties.size(); }

  // throw error and comment (with filename+code line)
  void throwRuntimeError(std::string message);
  /**
   * \brief return attribute as type
   *
   * returns an attribute after type conversion, e.g.
   * p.getAttribute<int>() returns an integer
   */
  template <typename T>
  T getAttribute(const std::string &attribute);
  /**
   * \brief set an attribute
   */
  template <typename T>
  void setAttribute(const std::string &attribute, const T &value);
  /**
   * \brief return true if a node has attributes
   */
  bool hasAttributes() { return _attributes.size() > 0; }
  /**
   * \brief return true if an attribute exists
   */
  bool hasAttribute(const std::string &attribute);
  /** for iterator-based access of Attributes */
  typedef std::map<std::string, std::string>::iterator AttributeIterator;
  /**
   * \brief returns an iterator to an attribute
   */
  AttributeIterator findAttribute(const std::string &attribute) {
    return _attributes.find(attribute);
  }
  /**
   * \brief returns an iterator to the first attribute
   */
  AttributeIterator firstAttribute() { return _attributes.begin(); }
  /**
   * \brief returns an iterator to the last attribute
   */
  AttributeIterator lastAttribute() { return _attributes.end(); }
  /**
   * \brief return attribute as type
   *
   * returns an attribute after type conversion, e.g.
   * p.getAttribute<int>() returns an integer
   */
  template <typename T>
  T getAttribute(AttributeIterator it);

  static int getIOindex() { return IOindex; };

 private:
  std::map<std::string, Property *> _map;
  std::map<std::string, std::string> _attributes;
  std::list<Property> _properties;

  std::string _name;
  std::string _value;
  std::string _path;

  static const int IOindex;
};

inline Property &Property::set(const std::string &key,
                               const std::string &value) {
  Property &p = get(key);
  p.value() = value;
  return p;
}

inline Property &Property::add(const std::string &key,
                               const std::string &value) {
  std::string path = _path;
  if (path != "") path = path + ".";
  _properties.push_back(Property(key, value, path + _name));
  _map[key] = &(_properties.back());
  return _properties.back();
}

inline bool Property::exists(const std::string &key) {
  try {
    get(key);
  } catch (std::exception &err) {
    return false;
  }
  return true;
}

bool load_property_from_xml(Property &p, std::string file);

// TO DO: write a better function for this!!!!
template <>
inline bool Property::as<bool>() const {
  if (_value == "true" || _value == "TRUE" || _value == "1")
    return true;
  else
    return false;
}

template <typename T>
inline T Property::as() const {
  return lexical_cast<T>(_value, "wrong type in " + _path + "." + _name + "\n");
}

template <typename T>
inline T Property::ifExistsReturnElseReturnDefault(const std::string &key,
                                                   T defaultvalue) {
  T result;
  if (this->exists(key)) {
    result = this->get(key).as<T>();
  } else {
    result = defaultvalue;
  }
  return result;
}

template <typename T>
inline T Property::ifExistsReturnElseThrowRuntimeError(const std::string &key) {
  T result;
  if (this->exists(key)) {
    result = this->get(key).as<T>();
  } else {
    throw runtime_error((boost::format("Error: %s is not found") % key).str());
  }
  return result;
}

template <typename T>
inline T Property::ifExistsAndinListReturnElseThrowRuntimeError(
    const std::string &key, std::vector<T> possibleReturns) {
  T result;
  result = ifExistsReturnElseThrowRuntimeError<T>(key);
  if (std::find(possibleReturns.begin(), possibleReturns.end(), result) ==
      possibleReturns.end()) {
    cerr << "Allowed options are: ";
    for (unsigned i = 0; i < possibleReturns.size(); ++i) {
      cerr << possibleReturns[i] << " ";
    }
    cerr << endl;
    throw runtime_error(
        (boost::format("Error: %s is not allowed") % key).str());
  }
  return result;
}

template <>
inline std::string Property::as<std::string>() const {
  std::string tmp(_value);
  boost::trim(tmp);
  return tmp;
}

template <>
inline vec Property::as<vec>() const {
  std::vector<double> tmp;
  Tokenizer tok(as<std::string>(), " ,");
  tok.ConvertToVector<double>(tmp);
  if (tmp.size() != 3)
    throw std::runtime_error("Vector has " +
                             boost::lexical_cast<std::string>(tmp.size()) +
                             " instead of three entries");
  return vec(tmp[0], tmp[1], tmp[2]);
}

template <>
inline Eigen::VectorXd Property::as<Eigen::VectorXd>() const {
  std::vector<double> tmp;
  Tokenizer tok(as<std::string>(), " ,");
  tok.ConvertToVector<double>(tmp);
  Eigen::VectorXd result;
  result.resize(tmp.size());
  for (int i = 0; i < result.size(); i++) {
    result(i) = tmp[i];
  }
  return result;
}

template <>
inline Eigen::Vector3d Property::as<Eigen::Vector3d>() const {
  std::vector<double> tmp;
  Tokenizer tok(as<std::string>(), " ,");
  tok.ConvertToVector<double>(tmp);
  Eigen::Vector3d result;
  if (int(tmp.size()) != result.size()) {
    throw std::runtime_error("Vector has " +
                             boost::lexical_cast<std::string>(tmp.size()) +
                             " instead of three entries");
  }
  result << tmp[0], tmp[1], tmp[2];
  return result;
}

template <>
inline std::vector<unsigned int> Property::as<std::vector<unsigned int> >()
    const {
  std::vector<unsigned int> tmp;
  Tokenizer tok(as<std::string>(), " ,");
  tok.ConvertToVector<unsigned int>(tmp);
  return tmp;
}

template <>
inline std::vector<int> Property::as<std::vector<int> >() const {
  std::vector<int> tmp;
  Tokenizer tok(as<std::string>(), " ,\n\t");
  tok.ConvertToVector<int>(tmp);
  return tmp;
}

template <>
inline std::vector<double> Property::as<std::vector<double> >() const {
  std::vector<double> tmp;
  Tokenizer tok(as<std::string>(), " ,\n\t");
  tok.ConvertToVector<double>(tmp);
  return tmp;
}

inline bool Property::hasAttribute(const std::string &attribute) {
  std::map<std::string, std::string>::iterator it;
  it = _attributes.find(attribute);
  if (it == _attributes.end()) return false;
  return true;
}

template <typename T>
inline T Property::getAttribute(
    std::map<std::string, std::string>::iterator it) {
  if (it != _attributes.end()) {
    return lexical_cast<T>((*it).second);
  } else {
    throw std::runtime_error("attribute " + (*it).first + " not found\n");
  }
}

template <typename T>
inline T Property::getAttribute(const std::string &attribute) {
  std::map<std::string, std::string>::iterator it;

  it = _attributes.find(attribute);

  if (it != _attributes.end()) {
    return lexical_cast<T>(_attributes[attribute],
                           "wrong type in attribute " + attribute +
                               " of element " + _path + "." + _name + "\n");
  } else {
    throw std::runtime_error("attribute " + attribute + " not found\n");
  }
}
template <typename T>
inline void Property::setAttribute(const std::string &attribute,
                                   const T &value) {
  _attributes[attribute] =
      lexical_cast<std::string>(value, "wrong type to set attribute");
}

inline void throwRuntimeError(std::string message) {}

}  // namespace tools
}  // namespace votca

#endif /* _VOTCA_TOOLS_PROPERTY_H */
