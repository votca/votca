
#if !defined(FACTORY_H_)
#define FACTORY_H_

#include "example.hpp"
#include "factory.hpp"
#include "parent.hpp"
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace pyxtp {

template <typename T>
class Factory {
 public:
  Factory() = default;
  template <typename U>
  friend Factory<U> &Calculators();

 private:
  using creator_t = std::unique_ptr<T> (*)();
  std::map<std::string, creator_t> _objects;

 public:
  static void RegisterAll();

  /**
   * \brief register an object
   * \param key identifier
   * \param creator create policy
   *
   * This function is called to register an object in the factory.
   */
  void Register(const std::string &key, creator_t creator);

  template <typename obj_t>
  void Register(const std::string &key);

  std::vector<std::string> getKeys() const {
    std::vector<std::string> key;
    key.reserve(_objects.size());
    for (const auto &pair : _objects) {
      key.push_back(pair.first);
    }
    return key;
  }
};

template <typename T>
inline Factory<T> &Calculators() {
  static Factory<T> instance;
  return instance;
}

template <typename parent, typename T>
std::unique_ptr<parent> create_policy_new() {
  return std::make_unique<T>();
}

template <typename T>
inline void Factory<T>::Register(const std::string &key, creator_t creator) {
  (void)_objects
      .insert(
          typename std::map<std::string, creator_t>::value_type(key, creator))
      .second;
}

template <typename T>
template <typename obj_t>
inline void Factory<T>::Register(const std::string &key) {
  Register(key, create_policy_new<T, obj_t>);
}

template <>
void Factory<Parent>::RegisterAll() {
  Calculators<Parent>().Register<Example>("example");
}

}  // namespace pyxtp
#endif  // FACTORY_H_
