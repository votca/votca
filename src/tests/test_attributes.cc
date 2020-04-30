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

#define BOOST_TEST_MODULE attributes_test
#include "votca/tools/attributes.h"
#include "votca/tools/types.h"
#include <boost/any.hpp>
#include <boost/test/unit_test.hpp>
#include <cmath>
#include <exception>
#include <iostream>
#include <typeinfo>
using namespace std;
using namespace votca::tools;
using namespace votca;

BOOST_AUTO_TEST_SUITE(attributes_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Attributes attr; }

BOOST_AUTO_TEST_CASE(accessors_test) {
  Attributes attr;

  BOOST_CHECK_EQUAL(attr == attr, true);
  BOOST_CHECK_EQUAL(attr.getContentLabel(), "");

  unordered_map<string, boost::any> values;
  int num = 134;
  double height = 159.32;
  string name = "George";
  values["Num"] = num;
  values["Height"] = height;
  values["Name"] = name;

  Attributes attr2(values);
  Attributes attr3(values);
  BOOST_CHECK_EQUAL(attr != attr2, true);
  BOOST_CHECK_EQUAL(attr2 == attr2, true);
  BOOST_CHECK_EQUAL(attr3 == attr2, true);
}

BOOST_AUTO_TEST_CASE(setters_test) {

  unordered_map<string, boost::any> values;
  int num = 134;
  double height = 159.32;
  string name = "George";
  values["Num"] = num;
  values["Height"] = height;
  values["Name"] = name;

  Attributes attr(values);

  string str{"Height=159.32,Name=George,Num=134;"};
  BOOST_CHECK_EQUAL(attr.getContentLabel(), str);

  unordered_map<string, boost::any> values2;
  Index ind2 = 2;
  Index ind1 = 1;
  values2["Second"] = ind2;
  values2["First"] = ind1;

  attr.set(values2);
  str = "First=1,Second=2;";
  BOOST_CHECK_EQUAL(attr.getContentLabel(), str);

  unordered_map<string, boost::any> values3;
  double weight = 101.43;
  values3["Height"] = height;
  values3["Weight"] = weight;

  attr.set(values3);
  str = "Height=159.32,Weight=101.43;";
  BOOST_CHECK_EQUAL(attr.getContentLabel(), str);

  unordered_map<string, boost::any> values4;
  string street = "Koogler St";
  values4["Name"] = name;
  values4["Address"] = street;

  attr.set(values4);
  str = "Address=Koogler St,Name=George;";
  BOOST_CHECK_EQUAL(attr.getContentLabel(), str);
}

BOOST_AUTO_TEST_CASE(comparisontest) {
  int val1 = 134;
  int val2 = 134;
  unordered_map<string, boost::any> values;
  values["a"] = val1;
  unordered_map<string, boost::any> values2;
  values2["b"] = val2;

  Attributes attr1(values);
  Attributes attr2(values2);

  BOOST_CHECK_EQUAL(cmpAttributes(attr1, attr2), true);
  BOOST_CHECK_EQUAL(cmpAttributes(attr2, attr1), false);

  vector<Attributes> vec_attr = {attr1, attr2};
  sort(vec_attr.begin(), vec_attr.end(), cmpAttributes);

  string str1{"a=134;"};
  string str2{"b=134;"};

  BOOST_CHECK_EQUAL(vec_attr.at(0).getContentLabel(), str1);
  BOOST_CHECK_EQUAL(vec_attr.at(1).getContentLabel(), str2);

  vector<Attributes> vec_attr2 = {attr2, attr1};
  sort(vec_attr2.begin(), vec_attr2.end(), cmpAttributes);

  BOOST_CHECK_EQUAL(vec_attr2.at(0).getContentLabel(), str1);
  BOOST_CHECK_EQUAL(vec_attr2.at(1).getContentLabel(), str2);
}

BOOST_AUTO_TEST_CASE(addtest) {

  Attributes attr1;
  Attributes attr2;
  int val1 = 134;
  int val2 = 134;
  std::string key1 = "a";
  std::string key2 = "b";

  attr1.add(key1, val1);
  attr2.add(key2, val2);

  BOOST_CHECK_EQUAL(cmpAttributes(attr1, attr2), true);
  BOOST_CHECK_EQUAL(cmpAttributes(attr2, attr1), false);

  vector<Attributes> vec_attr = {attr1, attr2};
  sort(vec_attr.begin(), vec_attr.end(), cmpAttributes);

  string str1{"a=134;"};
  string str2{"b=134;"};

  BOOST_CHECK_EQUAL(vec_attr.at(0).getContentLabel(), str1);
  BOOST_CHECK_EQUAL(vec_attr.at(1).getContentLabel(), str2);

  vector<Attributes> vec_attr2 = {attr2, attr1};
  sort(vec_attr2.begin(), vec_attr2.end(), cmpAttributes);

  BOOST_CHECK_EQUAL(vec_attr2.at(0).getContentLabel(), str1);
  BOOST_CHECK_EQUAL(vec_attr2.at(1).getContentLabel(), str2);
}

BOOST_AUTO_TEST_SUITE_END()
