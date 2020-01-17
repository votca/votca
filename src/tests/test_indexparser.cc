/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE index_parser_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/IndexParser.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(index_parser_test)

BOOST_AUTO_TEST_CASE(readstring_test) {
  IndexParser parser;

  std::string test = "1 5 3...11 15...17";
  std::vector<votca::Index> result = parser.CreateIndexVector(test);

  std::vector<votca::Index> ref = {1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 15, 16, 17};

  BOOST_CHECK_EQUAL(result.size(), ref.size());
  for (votca::Index i = 0; i < votca::Index(ref.size()); i++) {
    BOOST_CHECK_EQUAL(result[i], ref[i]);
  }
}

BOOST_AUTO_TEST_CASE(generatestring_test) {

  std::vector<votca::Index> input = {1, 3,  4,  5,  6,  7, 8,
                                     9, 10, 11, 15, 16, 17};
  IndexParser parser;
  std::string result = parser.CreateIndexString(input);
  std::string ref = "1 3...11 15...17";
  BOOST_CHECK_EQUAL(result, ref);
}

BOOST_AUTO_TEST_SUITE_END()
