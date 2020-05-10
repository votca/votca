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

#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE graphnode_test
#include "../../include/votca/tools/graphnode.h"
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <iostream>
using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(graphnode_test)

BOOST_AUTO_TEST_CASE(constructors_test) { GraphNode gn; }

BOOST_AUTO_TEST_CASE(accessors_test) {
  GraphNode gn;

  BOOST_CHECK_EQUAL(gn == gn, true);
  BOOST_CHECK_EQUAL(gn.getContentLabel().get(), "");

  unordered_map<string, votca::Index> int_vals = {{"Num", 134}};
  unordered_map<string, double> double_vals = {{"Height", 159.32}};
  unordered_map<string, string> str_vals = {{"Name", "George"}};
  GraphNode gn2(int_vals, double_vals, str_vals);
  GraphNode gn3(int_vals, double_vals, str_vals);
  BOOST_CHECK_EQUAL(gn != gn2, true);
  BOOST_CHECK_EQUAL(gn2 == gn2, true);
  BOOST_CHECK_EQUAL(gn3 == gn2, true);
}

BOOST_AUTO_TEST_CASE(add_test) {
  unordered_map<string, votca::Index> int_vals = {{"Num", 134}};
  unordered_map<string, double> double_vals = {{"Height", 159.32}};
  unordered_map<string, string> str_vals = {{"Name", "George"}};
  GraphNode gn2(int_vals, double_vals, str_vals);

  string str{"Height=159.32,Name=George,Num=134;"};
  BOOST_CHECK_EQUAL(gn2.getContentLabel().get(), str);

  unordered_map<string, votca::Index> int_vals2 = {{"Second", 2}, {"First", 1}};
  /// Remove all the indices from the graphnode
  gn2.removeAll<votca::Index>();
  gn2.add(int_vals2);
  str = "First=1,Height=159.32,Name=George,Second=2;";
  BOOST_CHECK_EQUAL(gn2.getContentLabel().get(), str);

  unordered_map<string, double> double_vals2 = {{"Height", 159.32},
                                                {"Weight", 101.43}};
  gn2.remove(string("Height"));
  gn2.add(double_vals2);
  str = "First=1,Height=159.32,Name=George,Second=2,Weight=101.43;";
  BOOST_CHECK_EQUAL(gn2.getContentLabel().get(), str);

  gn2.reset(string("Name"), string("George"));
  gn2.add(string("Address"), string("Koogler St"));
  str =
      "Address=Koogler "
      "St,First=1,Height=159.32,Name=George,Second=2,Weight=101.43;";
  BOOST_CHECK_EQUAL(gn2.getContentLabel().get(), str);
}

BOOST_AUTO_TEST_CASE(comparisontest) {
  unordered_map<string, votca::Index> int_vals1 = {{"a", 134}};
  unordered_map<string, votca::Index> int_vals2 = {{"b", 134}};

  unordered_map<string, double> double_vals;
  unordered_map<string, string> str_vals;

  GraphNode gn1(int_vals1, double_vals, str_vals);
  GraphNode gn2(int_vals2, double_vals, str_vals);

  BOOST_CHECK_EQUAL(cmpNode(gn1, gn2), true);
  BOOST_CHECK_EQUAL(cmpNode(gn2, gn1), false);

  vector<GraphNode> vec_gn = {gn1, gn2};
  sort(vec_gn.begin(), vec_gn.end(), cmpNode);

  string str1{"a=134;"};
  string str2{"b=134;"};

  BOOST_CHECK_EQUAL(vec_gn.at(0).getContentLabel().get(), str1);
  BOOST_CHECK_EQUAL(vec_gn.at(1).getContentLabel().get(), str2);

  vector<GraphNode> vec_gn2 = {gn2, gn1};
  sort(vec_gn2.begin(), vec_gn2.end(), cmpNode);

  BOOST_CHECK_EQUAL(vec_gn2.at(0).getContentLabel().get(), str1);
  BOOST_CHECK_EQUAL(vec_gn2.at(1).getContentLabel().get(), str2);
}

BOOST_AUTO_TEST_SUITE_END()
