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

#ifndef VOTCA_TOOLS_OBJECTFACTORY
#define VOTCA_TOOLS_OBJECTFACTORY

#include <boost/lexical_cast.hpp>
#include <iostream>
#include <list>
#include <map>
#include <stdexcept>

namespace votca {
namespace tools {

/**
    \brief template class for object factory

    This class is a template for an object factory. The factory creates an
   instance of an derived class, be giving a key (e.g. a string) for the object
   which identifies it uniquely. This allows the implementation of new features
   (e.g. new file formats, new mapping algorithms) without touching or
   recompiling existing bits of code.

    If you don't understand this, read the book by Alexandresku (Modern C++
   design) everything explained there in detail!
*/
template <typename key_t, typename T>
class ObjectFactory {
 private:
  typedef T *(*creator_t)();

 public:
  typedef T abstract_type;
  typedef std::map<key_t, creator_t> assoc_map;

  ObjectFactory() {}
  ~ObjectFactory(){};

  /**
   * \brief register an object
   * \param key identifier
   * \param creator create policy
   *
   * This function is called to register an object in the factory. After an
   * object is registered, an instance of it can be created by calling Create
   * specifying the corresponding key.
   */
  void Register(const key_t &key, creator_t creator);

  template <typename obj_t>
  void Register(const key_t &key);

  /**
     Create an instance of the object identified by key.
  */
  T *Create(const key_t &key);
  bool IsRegistered(const key_t &_id) const;

  static ObjectFactory<key_t, T> &Instance() {
    static ObjectFactory<key_t, T> _this;
    return _this;
  }

  const assoc_map &getObjects() { return _objects; }

 private:
  assoc_map _objects;
};

template <class parent, class T>
parent *create_policy_new() {
  return new T();
}

template <typename key_t, typename T>
inline void ObjectFactory<key_t, T>::Register(const key_t &key,
                                              creator_t creator) {
  (void)_objects.insert(typename assoc_map::value_type(key, creator)).second;
}

template <typename key_t, typename T>
template <typename obj_t>
inline void ObjectFactory<key_t, T>::Register(const key_t &key) {
  Register(key, create_policy_new<abstract_type, obj_t>);
}

template <typename key_t, typename T>
inline T *ObjectFactory<key_t, T>::Create(const key_t &key) {
  typename assoc_map::const_iterator it(_objects.find(key));
  if (it != _objects.end())
    return (it->second)();
  else
    throw std::runtime_error(
        "factory key " + boost::lexical_cast<std::string>(key) + " not found.");
}

template <typename key_t, typename T>
inline bool ObjectFactory<key_t, T>::IsRegistered(const key_t &_id) const {
  return (_objects.find(_id) != _objects.end());
}

template <typename object_type>
class ObjectFactoryRegister {
 public:
  template <typename factory_type, typename key_type>
  ObjectFactoryRegister(factory_type &factory, key_type &key) {
    factory.Register(
        key,
        &create_policy_new<typename factory_type::abstract_type, object_type>);
  }
};

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_OBJECTFACTORY
