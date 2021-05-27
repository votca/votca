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

#ifndef VOTCA_TOOLS_TOKENIZER_H
#define VOTCA_TOOLS_TOKENIZER_H

// Standard includes
#include <string>
#include <vector>

// Third party includes
#include "eigen.h"
#include "lexical_cast.h"
#include <boost/tokenizer.hpp>
#include <type_traits>
namespace votca {
namespace tools {

namespace internal {

template <typename T>
struct type {};

template <typename T,
          typename std::enable_if_t<std::is_arithmetic<T>::value &&
                                        !std::is_same<T, bool>::value,
                                    bool> = true>
inline T convert_impl(const std::string &s, type<T>) {
  return lexical_cast<T>(s, "Cannot create arithmetic type from" + s);
}

template <typename T,
          typename std::enable_if_t<
              std::is_constructible<T, std::string>::value, bool> = true>
inline T convert_impl(const std::string &s, type<T>) {
  return T(s);
}

inline bool convert_impl(const std::string &s, type<bool>) {
  if (s == "true" || s == "TRUE" || s == "1") {
    return true;

  } else if (s == "false" || s == "FALSE" || s == "0") {
    return false;
  } else {
    throw std::runtime_error("'" + s + "' cannot be converted to bool.");
  }
}

}  // namespace internal

/**
 * \brief break string into words
 *
 * This class wraps boost::tokenizer to break a string into words. A list of
 * delimeters can be freely choosen.
 */
class Tokenizer {
 public:
  using iterator = boost::tokenizer<boost::char_separator<char>>::iterator;

  /**
   * \brief startup tokenization
   * @param str string to break up
   * @param separators list of separators
   *
   * After initialization,the words can be accessed using the iterator
   * interface or directly transferred to a vector ToVector of ConvertToVector.
   */

  Tokenizer(const std::string &str, const char *separators) : _str(str) {
    boost::char_separator<char> sep(separators);
    tok_ = std::make_unique<boost::tokenizer<boost::char_separator<char>>>(_str,
                                                                           sep);
  }
  Tokenizer(const std::string &str, const std::string &separators)
      : Tokenizer(str, separators.c_str()){};

  /**
   * \brief iterator to first element
   * @return begin iterator
   */
  iterator begin() { return tok_->begin(); }
  /**
   * \brief end iterator
   * @return end iterator
   */
  iterator end() { return tok_->end(); }

  /**
   * \brief store all words in a vector of type T, does type conversion.
   * @return storage vector
   */
  template <class T = std::string>
  std::vector<T> ToVector() {
    std::vector<T> result;
    for (auto &seg : *this) {
      result.push_back(internal::convert_impl(seg, internal::type<T>{}));
    }
    return result;
  }

 private:
  std::unique_ptr<boost::tokenizer<boost::char_separator<char>>> tok_;
  std::string _str;
};

// Matches a string against a wildcard string such as &quot;*.*&quot; or
// &quot;bl?h.*&quot; etc. This is good for file globbing or to match hostmasks.
int wildcmp(const char *wild, const char *string);
int wildcmp(const std::string &wild, const std::string &string);

namespace internal {

template <class T>
inline std::vector<T> convert_impl(const std::string &s, type<std::vector<T>>) {
  return Tokenizer(s, " ,\n\t").ToVector<T>();
}

template <class T>
inline Eigen::Matrix<T, Eigen::Dynamic, 1> convert_impl(
    const std::string &s, type<Eigen::Matrix<T, Eigen::Dynamic, 1>>) {
  std::vector<T> tmp = convert_impl(s, type<std::vector<T>>{});
  return Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1>>(tmp.data(),
                                                         tmp.size());
}

template <class T>
inline Eigen::Matrix<T, 3, 1> convert_impl(const std::string &s,
                                           type<Eigen::Matrix<T, 3, 1>>) {
  std::vector<T> tmp = convert_impl(s, type<std::vector<T>>{});
  if (tmp.size() != 3) {
    throw std::runtime_error("Vector has " + std::to_string(tmp.size()) +
                             " instead of 3 entries");
  }
  return {tmp[0], tmp[1], tmp[2]};
}
}  // namespace internal

template <class T>
T convertFromString(const std::string &s) {
  return internal::convert_impl(s, internal::type<T>{});
}

}  // namespace tools
}  // namespace votca

#endif  // VOTCA_TOOLS_TOKENIZER_H
