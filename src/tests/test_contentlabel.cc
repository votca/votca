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

#define BOOST_TEST_MODULE contentlabel_test

// Standard includes
#include <iostream>
#include <string>
#include <unordered_map>

// Third party includes
#include <boost/any.hpp>
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/contentlabel.h"

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(contentlabel_test)

BOOST_AUTO_TEST_CASE(contentlabel_constructor) { ContentLabel label(); }

BOOST_AUTO_TEST_CASE(contentlabel_initializer) {

  unordered_map<string, boost::any> contents;
  contents["Name"] = string("Joe");
  contents["Age"] = int(32);
  contents["Height"] = double(1.54);
  ContentLabel label(contents);

  BOOST_CHECK(label.isEmpty() == false);
  label.clear();
  BOOST_CHECK(label.isEmpty() == true);
}

BOOST_AUTO_TEST_CASE(contentlabel_get) {

  unordered_map<string, boost::any> contents;
  contents["Name"] = string("Joe");
  contents["Age"] = int(32);
  contents["Height"] = double(1.54);
  ContentLabel label(contents);

  string full_str_label = "Age=32,Height=1.54,Name=Joe;";
  string brief_str_label = "32,1.54,Joe;";
  BOOST_CHECK(full_str_label.compare(label.get()) == 0);
  BOOST_CHECK(brief_str_label.compare(label.get(LabelType::concise)) == 0);
}

BOOST_AUTO_TEST_CASE(contentlabel_operators) {

  unordered_map<string, boost::any> contents;
  contents["Name"] = string("Joe");
  contents["Age"] = int(32);
  contents["Height"] = double(1.54);
  ContentLabel label(contents);

  ContentLabel label2 = label;

  BOOST_CHECK(label == label2);

  unordered_map<string, boost::any> contents2;
  contents2["Name"] = string("Randy");
  contents2["Age"] = int(21);
  contents2["Height"] = double(1.64);
  ContentLabel label3(contents2);

  // Because sorted total string length first
  BOOST_CHECK(label < label3);
  BOOST_CHECK(label <= label3);

  BOOST_CHECK(!(label > label3));
  BOOST_CHECK(!(label >= label3));

  BOOST_CHECK(label3 > label);
  BOOST_CHECK(label3 >= label);

  BOOST_CHECK(!(label3 < label));
  BOOST_CHECK(!(label3 <= label));

  BOOST_CHECK(label3 != label);
}

BOOST_AUTO_TEST_CASE(contentlabel_append) {

  unordered_map<string, boost::any> contents;
  contents["Name"] = string("Joe");
  contents["Age"] = int(32);
  contents["Height"] = double(1.54);
  ContentLabel label(contents);

  unordered_map<string, boost::any> contents2;
  contents2["Name"] = string("Randy");
  contents2["Age"] = int(21);
  contents2["Height"] = double(1.64);

  ContentLabel label2(contents2);

  label.append(label2);

  string full_str_label =
      "Age=32,Height=1.54,Name=Joe;Age=21,Height=1.64,Name=Randy;";
  string brief_str_label = "32,1.54,Joe;21,1.64,Randy;";
  BOOST_CHECK(full_str_label.compare(label.get()) == 0);
  BOOST_CHECK(brief_str_label.compare(label.get(LabelType::concise)) == 0);
}

BOOST_AUTO_TEST_SUITE_END()
