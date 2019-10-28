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

#ifndef _VOTCA_TOOLS_TOKENIZER_H
#define _VOTCA_TOOLS_TOKENIZER_H

#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <string>
#include <vector>

namespace votca {
namespace tools {

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
    tok = std::make_unique<boost::tokenizer<boost::char_separator<char>>>(_str,
                                                                          sep);
  }
  Tokenizer(const std::string &str, const std::string &separators)
      : Tokenizer(str, separators.c_str()){};

  /**
   * \brief iterator to first element
   * @return begin iterator
   */
  iterator begin() { return tok->begin(); }
  /**
   * \brief end iterator
   * @return end iterator
   */
  iterator end() { return tok->end(); }

  /**
   * \brief store all words in a vector of strings.
   * @param v storage vector
   *
   * This class appends all words to a vector of strings.
   */
  void ToVector(std::vector<std::string> &v) {
    for (iterator iter = begin(); iter != end(); ++iter) {
      v.push_back(*iter);
    }
  }

  std::vector<std::string> ToVector() {
    std::vector<std::string> result;
    for (iterator iter = begin(); iter != end(); ++iter) {
      result.push_back(*iter);
    }
    return result;
  }

  /**
   * \brief store all words in a vector with type conversion.
   * @param v storage vector
   *
   * This class appends all words to a vector of arbitrary type (e.g. double)
   * and also does type conversion.
   */
  template <typename T>
  void ConvertToVector(std::vector<T> &v) {
    std::vector<std::string> tmp;
    ToVector(tmp);
    v.resize(tmp.size());
    typename std::vector<T>::iterator viter = v.begin();
    typename std::vector<std::string>::iterator iter;
    for (iter = tmp.begin(); iter != tmp.end(); ++iter, ++viter) {
      *viter = boost::lexical_cast<T, std::string>(*iter);
    }
  }

 private:
  std::unique_ptr<boost::tokenizer<boost::char_separator<char>>> tok;
  std::string _str;
};

// Matches a string against a wildcard string such as &quot;*.*&quot; or
// &quot;bl?h.*&quot; etc. This is good for file globbing or to match hostmasks.
int wildcmp(const char *wild, const char *string);
int wildcmp(const std::string &wild, const std::string &string);

}  // namespace tools
}  // namespace votca

#endif /* _tools_H */
