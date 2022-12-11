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

#ifndef VOTCA_TOOLS_PROPERTY_H
#define VOTCA_TOOLS_PROPERTY_H

// Standard includes
#include <cstdlib>
#include <exception>
#include <iostream>
#include <list>
#include <map>
#include <stdexcept>
#include <string>

// Third party includes
#include <boost/algorithm/string/trim.hpp>
#include <boost/format.hpp>
#include <utility>

// Local VOTCA includes
#include "eigen.h"
#include "lexical_cast.h"
#include "tokenizer.h"
#include "types.h"

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
 * Supported formats are XML, TXT, HLP
 */
class Property {

  /// \brief outputs the property to the ostream
  friend std::ostream &operator<<(std::ostream &out, const Property &p);

 public:
  Property() = default;
  Property(const std::string &name, const std::string &value,
           const std::string &path)
      : name_(name), value_(value), path_(path) {}

  /**
   * \brief add a new property to structure
   * @param key identifier
   * @param value value
   * @return reference to the created Property object
   */
  Property &add(const std::string &key, const std::string &value);

  /**
   * \brief add a copy of an existing property to a property
   * @param other other property
   */
  Property &add(const Property &other);

  /**
   * \brief add a new property tree to structure
   * @param key identifier
   * @param value value
   * @return reference to the created Property object
   * This function adds a property specified by key separated
   * by "." to step down hierarchy. If the intermediate properties are not
   * found they are also created.
   */
  Property &addTree(const std::string &key, const std::string &value);
  Property &addTree(const std::vector<std::string> &key,
                    const std::string &value);

  /**
   * \brief set value of existing property
   * @param key identifier
   * @param value value
   * @return reference to the created Property object
   * If more than property with this name exists, return the last added one.
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
   * If more than property with this name exists, return the last added one.
   */
  Property &get(const std::string &key);
  const Property &get(const std::string &key) const;

  /**
   * \brief adds new or gets existing property
   * @param key identifier
   * @return Reference to property object
   *
   * This function tries to find a property specified by key separated
   * by "." to step down hierarchy. If the property is not
   * found a property with that name is added and returned.
   * If more than property with this name exists, return the last added one.
   */
  Property &getOradd(const std::string &key);

  /**
   * \brief check whether property exists
   * @param key identifier
   * @return true or false
   */
  bool exists(const std::string &key) const;

  template <typename T>
  T ifExistsReturnElseReturnDefault(const std::string &key,
                                    T defaultvalue) const;

  /**
   * \brief select property based on a filter
   * @param filter
   * @return list of pointers to property objects
   *
   * returns a list of properties that match the key criteria including
   * wildcard "*". Example: "base.item*.value"
   */
  std::vector<Property *> Select(const std::string &filter);
  std::vector<const Property *> Select(const std::string &filter) const;

  /**
   * \brief reference to value of property
   * @return std::string content
   */
  std::string &value() { return value_; }
  const std::string &value() const { return value_; }
  /**
   * \brief name of property
   * @return name
   */
  std::string &name() { return name_; }
  const std::string &name() const { return name_; }
  /**
   * \brief full path of property (including parents)
   * @return path
   *
   * e.g. cg.inverse.value
   */
  std::string &path() { return path_; }
  const std::string &path() const { return path_; }
  /**
   * \brief return value as type
   *
   * returns the value after type conversion, e.g.
   * p.as<Index>() returns an integer
   */
  template <typename T>
  T as() const;

  /**
   * \brief does the property have children?
   * \return true or false
   */
  bool HasChildren() const { return !map_.empty(); }

  /// iterator to iterate over properties
  using iterator = std::vector<Property>::iterator;
  using const_iterator = std::vector<Property>::const_iterator;
  /// \brief iterator to first child property
  iterator begin() { return properties_.begin(); }
  const_iterator begin() const { return properties_.begin(); }
  /// \brief end iterator for child properties
  iterator end() { return properties_.end(); }
  const_iterator end() const { return properties_.end(); }
  /// \brief number of child properties
  Index size() const { return Index(properties_.size()); }

  /**
   * \brief deletes all children that fulfill a condition
   * @param condition unary function which takes a const reference to a property
   * and returns a bool
   *
   */
  template <class cond>
  void deleteChildren(cond condition);

  /**
   * \brief return attribute as type
   *
   * returns an attribute after type conversion, e.g.
   * p.getAttribute<Index>() returns an integer
   */
  template <typename T>
  T getAttribute(const std::string &attribute) const;
  /**
   * \brief set an attribute
   */
  template <typename T>
  void setAttribute(const std::string &attribute, const T &value);

  /**
   * \brief deletes an attribute
   */
  void deleteAttribute(const std::string &attribute);
  /**
   * \brief return true if a node has attributes
   */
  bool hasAttributes() const { return attributes_.size() > 0; }
  /**
   * \brief return true if an attribute exists
   */
  bool hasAttribute(const std::string &attribute) const;
  /** for iterator-based access of Attributes */
  typedef std::map<std::string, std::string>::iterator AttributeIterator;
  typedef std::map<std::string, std::string>::const_iterator
      const_AttributeIterator;
  /**
   * \brief returns an iterator to an attribute
   */
  AttributeIterator findAttribute(const std::string &attribute) {
    return attributes_.find(attribute);
  }
  const_AttributeIterator findAttribute(const std::string &attribute) const {
    return attributes_.find(attribute);
  }
  /**
   * \brief returns an iterator to the first attribute
   */
  AttributeIterator firstAttribute() { return attributes_.begin(); }
  const_AttributeIterator firstAttribute() const { return attributes_.begin(); }
  /**
   * \brief returns an iterator to the last attribute
   */
  AttributeIterator lastAttribute() { return attributes_.end(); }
  const_AttributeIterator lastAttribute() const { return attributes_.end(); }
  /**
   * \brief return attribute as type
   *
   * returns an attribute after type conversion, e.g.
   * p.getAttribute<Index>() returns an integer
   */
  template <typename T>
  T getAttribute(AttributeIterator it);

  template <typename T>
  T getAttribute(const_AttributeIterator it) const;

  void LoadFromXML(std::string filename);

  static Index getIOindex() { return IOindex; };

 private:
  std::map<std::string, std::vector<Index>> map_;
  std::map<std::string, std::string> attributes_;
  std::vector<Property> properties_;

  std::string name_ = "";
  std::string value_ = "";
  std::string path_ = "";

  static const Index IOindex;
};

template <typename T>
inline T Property::as() const {
  std::string trimmed = value_;
  boost::trim(trimmed);
  try {
    return convertFromString<T>(trimmed);
  } catch (std::runtime_error &e) {
    throw std::runtime_error("Property with name '" + name() + "' in path '" +
                             path() + "' and value :" + e.what());
  }
}

template <typename T>
inline T Property::getAttribute(
    std::map<std::string, std::string>::const_iterator it) const {
  if (it != attributes_.end()) {
    return convertFromString<T>(it->second);
  } else {
    std::stringstream s;
    s << *this << std::endl;
    throw std::runtime_error(s.str() + "attribute " + it->first +
                             " not found\n");
  }
}

bool load_property_from_xml(Property &p, std::string file);

// TO DO: write a better function for this!!!!
template <>
inline bool Property::as<bool>() const {
  if (value_ == "true" || value_ == "TRUE" || value_ == "1")
    return true;
  else
    return false;
}

template <typename T>
inline T Property::getAttribute(const std::string &attribute) const {
  std::map<std::string, std::string>::const_iterator it =
      attributes_.find(attribute);
  return getAttribute<T>(it);
}
template <typename T>
inline void Property::setAttribute(const std::string &attribute,
                                   const T &value) {
  attributes_[attribute] =
      lexical_cast<std::string>(value, "wrong type to set attribute");
}

template <typename T>
inline T Property::ifExistsReturnElseReturnDefault(const std::string &key,
                                                   T defaultvalue) const {
  T result;
  if (this->exists(key)) {
    result = this->get(key).as<T>();
  } else {
    result = defaultvalue;
  }
  return result;
}

template <>
inline std::string Property::ifExistsReturnElseReturnDefault(
    const std::string &key, std::string defaultvalue) const {
  std::string result;
  if (this->exists(key)) {
    result = this->get(key).as<std::string>();
    if (result.empty()) {
      result = defaultvalue;
    }
  } else {
    result = defaultvalue;
  }
  return result;
}

template <class cond>
void Property::deleteChildren(cond condition) {

  properties_.erase(
      std::remove_if(properties_.begin(), properties_.end(), condition),
      properties_.end());

  // rebuild map_
  map_.clear();
  Index index = 0;
  for (const auto &prop : properties_) {
    map_[prop.name()].push_back(index);
    index++;
  }
  return;
}

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_PROPERTY_H
