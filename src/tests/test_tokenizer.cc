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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE tokenizer_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include <votca/tools/tokenizer.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(tokenizer_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  string str1 = "blah,ya";
  string separators = ",";
  Tokenizer tok(str1, separators.c_str());
  Tokenizer tok2(str1, separators);
}

BOOST_AUTO_TEST_CASE(tokenizer_test) {
  string str1 = "blah,ya";
  string separators = ",";
  Tokenizer tok(str1, separators);
  std::vector<std::string> result = tok.ToVector();
  BOOST_CHECK_EQUAL(result[0], "blah");
  BOOST_CHECK_EQUAL(result[1], "ya");

  string str2 = "";
  Tokenizer tok2(str2, separators);
  BOOST_CHECK_EQUAL(tok2.begin() == tok2.end(), true);
  std::vector<std::string> result2 = tok2.ToVector();
  BOOST_CHECK_EQUAL(result2.size(), 0);

  string str3 = "hello";
  Tokenizer tok3(str3, separators);
  BOOST_CHECK_EQUAL(tok3.begin() == tok3.end(), false);
  std::vector<std::string> result3 = tok3.ToVector();
  BOOST_CHECK_EQUAL(result3[0], "hello");
}

BOOST_AUTO_TEST_CASE(wildcmp_test) {
  string wildcard = "";
  string potential_match = "";

  auto out = wildcmp(wildcard.c_str(), potential_match.c_str());
  BOOST_CHECK_EQUAL(out, 1);

  string wildcard2 = "file";
  out = wildcmp(wildcard2.c_str(), potential_match.c_str());
  BOOST_CHECK_EQUAL(out, 0);

  string potential_match2 = "file2";
  out = wildcmp(wildcard.c_str(), potential_match2.c_str());
  BOOST_CHECK_EQUAL(out, 0);
  out = wildcmp(wildcard2.c_str(), potential_match2.c_str());
  BOOST_CHECK_EQUAL(out, 0);

  string wildcard3 = "file*";
  out = wildcmp(wildcard3.c_str(), potential_match2.c_str());
  BOOST_CHECK_EQUAL(out, 1);

  string wildcard4 = "file*.txt";
  out = wildcmp(wildcard4.c_str(), potential_match2.c_str());
  BOOST_CHECK_EQUAL(out, 0);

  string potential_match3 = "file1.txt";
  out = wildcmp(wildcard4.c_str(), potential_match3.c_str());
  BOOST_CHECK_EQUAL(out, 1);
}

BOOST_AUTO_TEST_CASE(wildcmp_test2) {
  string wildcard = "";
  string potential_match = "";

  auto out = wildcmp(wildcard, potential_match);
  BOOST_CHECK_EQUAL(out, 1);

  string wildcard2 = "file";
  out = wildcmp(wildcard2, potential_match);
  BOOST_CHECK_EQUAL(out, 0);

  string potential_match2 = "file2";
  out = wildcmp(wildcard, potential_match2);
  BOOST_CHECK_EQUAL(out, 0);
  out = wildcmp(wildcard2, potential_match2);
  BOOST_CHECK_EQUAL(out, 0);

  string wildcard3 = "file*";
  out = wildcmp(wildcard3, potential_match2);
  BOOST_CHECK_EQUAL(out, 1);

  string wildcard4 = "file*.txt";
  out = wildcmp(wildcard4, potential_match2);
  BOOST_CHECK_EQUAL(out, 0);

  string potential_match3 = "file1.txt";
  out = wildcmp(wildcard4, potential_match3);
  BOOST_CHECK_EQUAL(out, 1);
}

BOOST_AUTO_TEST_SUITE_END()
