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

  /// Ensuring that the clear and isEmpty methods work
  BOOST_CHECK(label.isEmpty() == false);
  label.clear();
  BOOST_CHECK(label.isEmpty() == true);
}

BOOST_AUTO_TEST_CASE(contentlabel_get) {

  /*
   * Labels are built by first sorting the key value pairs in order
   * of thier keys.
   *
   * The key value pairs are separated by commas, and closure of the 
   * content group is indicated with a semicolon. 
   *
   * The brief label will omit the keys and only display the values
   */
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

  /*
   * This test demonstrates what happens when two content labels are
   * appended. 
   *
   * The first label has the form
   *
   * Age=32,Height=1.54,Name=Joe;
   *
   * And the second
   *
   * Age=21,Height=1.64,Name=Randy;
   * 
   * Appending them is a simple matter of appending thier strings
   */ 
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
 
  size_t num_equals = 6; 
  std::cout << "Length of label " << label.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label.length()-num_equals) << std::endl;
  BOOST_CHECK(label.getCharLen() == full_str_label.length()-num_equals);

  string brief_str_label = "32,1.54,Joe;21,1.64,Randy;";
  BOOST_CHECK(full_str_label.compare(label.get()) == 0);
  BOOST_CHECK(brief_str_label.compare(label.get(LabelType::concise)) == 0);
}

BOOST_AUTO_TEST_CASE(contentlabel_make_branch) {

  /*
   * This test demonstrates what happens when a content label is turned into
   * a branch 
   *
   * Age=32,Height=1.54,Name=Joe;Age=21,Height=1.64,Name=Randy;
   * 
   * Should become
   *
   * {{Age=32,Height=1.54,Name=Joe}{Age=21,Height=1.64,Name=Randy}}
   */ 
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

  // Before making it a branch lets check the different labels
  string starting_node_str_label = "Age=32,Height=1.54,Name=Joe;";
  BOOST_CHECK(starting_node_str_label.compare(label.get(LabelType::verbose_start_node)) == 0);
  string terminating_node_str_label = "Age=21,Height=1.64,Name=Randy;";
  BOOST_CHECK(terminating_node_str_label.compare(label.get(LabelType::verbose_terminating_node)) == 0);
  string middle_nodes = "";
  BOOST_CHECK(middle_nodes.compare(label.get(LabelType::verbose_without_end_nodes)) == 0);

  label.makeBranchLabel();

  string full_str_label =
      "(Age=32,Height=1.54,Name=Joe)(Age=21,Height=1.64,Name=Randy)";

  size_t num_equals = 6; 
  BOOST_CHECK(label.getCharLen() == full_str_label.length()-num_equals);

  std::cout << "Length of label " << label.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label.length()-num_equals) << std::endl;
  string brief_str_label = "(32,1.54,Joe)(21,1.64,Randy)";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << label.get() << std::endl;
  std::cout << full_str_label << std::endl;
  BOOST_CHECK(full_str_label.compare(label.get()) == 0);
  BOOST_CHECK(brief_str_label.compare(label.get(LabelType::concise)) == 0);
  BOOST_CHECK(label.isBranch());

  unordered_map<string, boost::any> contents3;
  contents3["Name"] = string("Sarah");
  contents3["Age"] = int(30);
  contents3["Height"] = double(1.34);

  ContentLabel label3(contents3);

  ContentLabel labelAll(contents);
  labelAll.append(label2);
  labelAll.append(label3);

  labelAll.makeBranchLabel();
  string full_str_label2 =
      "(Age=32,Height=1.54,Name=Joe){Age=21,Height=1.64,Name=Randy}(Age=30,Height=1.34,Name=Sarah)";
  
  size_t num_equals2 = 9; 
  BOOST_CHECK(labelAll.getCharLen() == full_str_label2.length()-num_equals2);

  std::cout << "Length of label " << labelAll.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label2.length()-num_equals) << std::endl;
  string brief_str_label2 = "(32,1.54,Joe){21,1.64,Randy}(30,1.34,Sarah)";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << labelAll.get() << std::endl;
  std::cout << full_str_label2 << std::endl;
  BOOST_CHECK(full_str_label2.compare(labelAll.get()) == 0);
  BOOST_CHECK(brief_str_label2.compare(labelAll.get(LabelType::concise)) == 0);
  BOOST_CHECK(labelAll.isBranch());

}

BOOST_AUTO_TEST_CASE(contentlabel_make_furcation_stem) {

  /*
   * This test demonstrates what happens when a content label is turned into
   * a branch 
   *
   * Age=32,Height=1.54,Name=Joe;Age=21,Height=1.64,Name=Randy;
   * 
   * Should become
   *
   * {{Age=32,Height=1.54,Name=Joe}{Age=21,Height=1.64,Name=Randy}}
   */ 
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

  // Before making it a branch lets check the different labels
  string starting_node_str_label = "Age=32,Height=1.54,Name=Joe;";
  BOOST_CHECK(starting_node_str_label.compare(label.get(LabelType::verbose_start_node)) == 0);
  string terminating_node_str_label = "Age=21,Height=1.64,Name=Randy;";
  BOOST_CHECK(terminating_node_str_label.compare(label.get(LabelType::verbose_terminating_node)) == 0);
  string middle_nodes = "";
  BOOST_CHECK(middle_nodes.compare(label.get(LabelType::verbose_without_end_nodes)) == 0);

  label.makeBranchLabel();
  label.makeFurcationStem();

  string full_str_label =
      "(Age=32,Height=1.54,Name=Joe)(Age=21,Height=1.64,Name=Randy)<";

  size_t num_equals = 6; 
  BOOST_CHECK(label.getCharLen() == full_str_label.length()-num_equals);

  std::cout << "Length of label " << label.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label.length()-num_equals) << std::endl;
  string brief_str_label = "(32,1.54,Joe)(21,1.64,Randy)<";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << label.get() << std::endl;
  std::cout << full_str_label << std::endl;
  BOOST_CHECK(full_str_label.compare(label.get()) == 0);
  BOOST_CHECK(brief_str_label.compare(label.get(LabelType::concise)) == 0);
  BOOST_CHECK(label.isBranch() == false);

  unordered_map<string, boost::any> contents3;
  contents3["Name"] = string("Sarah");
  contents3["Age"] = int(30);
  contents3["Height"] = double(1.34);

  ContentLabel label3(contents3);

  ContentLabel labelAll(contents);
  labelAll.append(label2);
  labelAll.append(label3);

  string full_str_label2 =
      "(Age=32,Height=1.54,Name=Joe){Age=21,Height=1.64,Name=Randy}(Age=30,Height=1.34,Name=Sarah)<";
 
  labelAll.makeBranchLabel(); 
  labelAll.makeFurcationStem();
  size_t num_equals2 = 9; 
  BOOST_CHECK(labelAll.getCharLen() == full_str_label2.length()-num_equals2);

  std::cout << "Length of label " << labelAll.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label2.length()-num_equals) << std::endl;
  string brief_str_label2 = "(32,1.54,Joe){21,1.64,Randy}(30,1.34,Sarah)<";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << labelAll.get() << std::endl;
  std::cout << full_str_label2 << std::endl;
  BOOST_CHECK(full_str_label2.compare(labelAll.get()) == 0);
  BOOST_CHECK(brief_str_label2.compare(labelAll.get(LabelType::concise)) == 0);
  BOOST_CHECK(labelAll.isBranch() == false);
}

BOOST_AUTO_TEST_CASE(contentlabel_make_furcation_branch) {

  /*
   * This test demonstrates what happens when a content label is turned into
   * a branch 
   *
   * Age=32,Height=1.54,Name=Joe;Age=21,Height=1.64,Name=Randy;
   * 
   * Should become
   *
   * {{Age=32,Height=1.54,Name=Joe}{Age=21,Height=1.64,Name=Randy}}
   */ 
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

  // Before making it a branch lets check the different labels
  string starting_node_str_label = "Age=32,Height=1.54,Name=Joe;";
  BOOST_CHECK(starting_node_str_label.compare(label.get(LabelType::verbose_start_node)) == 0);
  string terminating_node_str_label = "Age=21,Height=1.64,Name=Randy;";
  BOOST_CHECK(terminating_node_str_label.compare(label.get(LabelType::verbose_terminating_node)) == 0);
  string middle_nodes = "";
  BOOST_CHECK(middle_nodes.compare(label.get(LabelType::verbose_without_end_nodes)) == 0);

  label.makeBranchLabel();
  label.makeFurcationBranch();

  string full_str_label =
      "[(Age=21,Height=1.64,Name=Randy)]";

  size_t num_equals = 3; 
  BOOST_CHECK(label.getCharLen() == full_str_label.length()-num_equals);

  std::cout << "Length of label " << label.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label.length()-num_equals) << std::endl;
  string brief_str_label = "[(21,1.64,Randy)]";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << label.get() << std::endl;
  std::cout << full_str_label << std::endl;
  BOOST_CHECK(full_str_label.compare(label.get()) == 0);
  BOOST_CHECK(brief_str_label.compare(label.get(LabelType::concise)) == 0);

  unordered_map<string, boost::any> contents3;
  contents3["Name"] = string("Sarah");
  contents3["Age"] = int(30);
  contents3["Height"] = double(1.34);

  ContentLabel label3(contents3);

  ContentLabel labelAll(contents);
  labelAll.append(label2);
  labelAll.append(label3);

  labelAll.makeBranchLabel();
  labelAll.makeFurcationBranch();
  string full_str_label2 =
      "[{Age=21,Height=1.64,Name=Randy}(Age=30,Height=1.34,Name=Sarah)]";
  
  size_t num_equals2 = 6; 
  BOOST_CHECK(labelAll.getCharLen() == full_str_label2.length()-num_equals2);

  std::cout << "Length of label " << labelAll.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label2.length()-num_equals) << std::endl;
  string brief_str_label2 = "[{21,1.64,Randy}(30,1.34,Sarah)]";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << labelAll.get() << std::endl;
  std::cout << full_str_label2 << std::endl;
  BOOST_CHECK(full_str_label2.compare(labelAll.get()) == 0);
  BOOST_CHECK(brief_str_label2.compare(labelAll.get(LabelType::concise)) == 0);
}

BOOST_AUTO_TEST_CASE(contentlabel_make_tree) {

  /*
   * This test demonstrates what happens when a content label is turned into
   * a branch 
   *
   * Age=32,Height=1.54,Name=Joe;Age=21,Height=1.64,Name=Randy;
   * 
   * Should become
   *
   * {{Age=32,Height=1.54,Name=Joe}{Age=21,Height=1.64,Name=Randy}}
   */ 
  unordered_map<string, boost::any> contents;
  contents["Name"] = string("Joe");
  contents["Age"] = int(32);
  contents["Height"] = double(1.54);

  unordered_map<string, boost::any> contents2;
  contents2["Name"] = string("Randy");
  contents2["Age"] = int(21);
  contents2["Height"] = double(1.64);

  unordered_map<string, boost::any> contents3;
  contents3["Name"] = string("Sally");
  contents3["Age"] = int(20);
  contents3["Height"] = double(1.14);

  ContentLabel stem(contents);
  ContentLabel cl2(contents2);
  ContentLabel cl3(contents3);
  stem.append(cl2);
  stem.append(cl3);

  stem.makeBranchLabel();
  stem.makeFurcationStem();

  unordered_map<string, boost::any> contents4;
  contents4["Name"] = string("Felix");
  contents4["Age"] = int(38);
  contents4["Height"] = double(1.44);

  unordered_map<string, boost::any> contents5;
  contents5["Name"] = string("Malcom");
  contents5["Age"] = int(56);
  contents5["Height"] = double(1.62);

  ContentLabel furcation(contents3);
  ContentLabel cl4(contents4);
  ContentLabel cl5(contents5);
  furcation.append(cl4); 
  furcation.append(cl5); 

  furcation.makeBranchLabel();
  furcation.makeFurcationBranch();

  stem.append(furcation);

  string full_str_label =
  "(Age=32,Height=1.54,Name=Joe){Age=21,Height=1.64,Name=Randy}(Age=20,Height=1.14,Name=Sally)<[{Age=38,Height=1.44,Name=Felix}(Age=56,Height=1.62,Name=Malcom)]";

  std::cout << "Expected label" << std::endl;
  std::cout << full_str_label << std::endl;
  std::cout << "Actual label" << std::endl;
  std::cout << stem.get() << std::endl; 
  int num_equals = 15; 
  std::cout << "Length of label " << stem.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label.length()-num_equals) << std::endl;
  BOOST_CHECK(stem.getCharLen() == full_str_label.length()-num_equals);

  string brief_str_label = 
    "(32,1.54,Joe){21,1.64,Randy}(20,1.14,Sally)<[{38,1.44,Felix}(56,1.62,Malcom)]";

  std::cout << full_str_label << std::endl;
  BOOST_CHECK(full_str_label.compare(stem.get()) == 0);

  std::cout << "Expected label" << std::endl;
  std::cout << brief_str_label << std::endl;
  std::cout << "Actual label" << std::endl;
  std::cout << stem.get(LabelType::concise) << std::endl; 

  BOOST_CHECK(brief_str_label.compare(stem.get(LabelType::concise)) == 0);

  // Lets append another furcation branch

  unordered_map<string, boost::any> contents6;
  contents6["Name"] = string("Mitchell");
  contents6["Age"] = int(72);
  contents6["Height"] = double(1.34);

  unordered_map<string, boost::any> contents7;
  contents7["Name"] = string("Tina");
  contents7["Age"] = int(18);
  contents7["Height"] = double(1.1);

  ContentLabel furcation2(contents3);
  ContentLabel cl6(contents6);
  ContentLabel cl7(contents7);
  furcation2.append(cl6); 
  furcation2.append(cl7); 

  furcation2.makeBranchLabel();
  furcation2.makeFurcationBranch();

  stem.append(furcation2);

  full_str_label =
  "(Age=32,Height=1.54,Name=Joe){Age=21,Height=1.64,Name=Randy}(Age=20,Height=1.14,Name=Sally)<[{Age=38,Height=1.44,Name=Felix}(Age=56,Height=1.62,Name=Malcom),{Age=72,Height=1.34,Name=Mitchell}(Age=18,Height=1.1,Name=Tina)]";

  std::cout << "Expected label" << std::endl;
  std::cout << full_str_label << std::endl;
  std::cout << "Actual label" << std::endl;
  std::cout << stem.get() << std::endl; 

  num_equals = 21;
  BOOST_CHECK(stem.getCharLen() == full_str_label.length()-num_equals);

  std::cout << "Label after making into a branch" << std::endl;
  std::cout << stem.get() << std::endl;
  std::cout << full_str_label << std::endl;
  BOOST_CHECK(full_str_label.compare(stem.get()) == 0);

  brief_str_label = 
    "(32,1.54,Joe){21,1.64,Randy}(20,1.14,Sally)<[{38,1.44,Felix}(56,1.62,Malcom),{72,1.34,Mitchell}(18,1.1,Tina)]";

  std::cout << "Expected label" << std::endl;
  std::cout << brief_str_label << std::endl;
  std::cout << "Actual label" << std::endl;
  std::cout << stem.get(LabelType::concise) << std::endl; 
  BOOST_CHECK(brief_str_label.compare(stem.get(LabelType::concise)) == 0);


/*
  unordered_map<string, boost::any> contents3;
  contents3["Name"] = string("Sarah");
  contents3["Age"] = int(30);
  contents3["Height"] = double(1.34);

  ContentLabel label3(contents3);

  ContentLabel labelAll(contents);
  labelAll.append(label2);
  labelAll.append(label3);

  labelAll.makeBranchLabel();
  labelAll.makeFurcationBranch();
  string full_str_label2 =
      "[{Age=21,Height=1.64,Name=Randy}(Age=30,Height=1.34,Name=Sarah)]";
  
  size_t num_equals2 = 6; 
  BOOST_CHECK(labelAll.getCharLen() == full_str_label2.length()-num_equals2);

  std::cout << "Length of label " << labelAll.getCharLen() << std::endl;
  std::cout << "Expected size " << (full_str_label2.length()-num_equals) << std::endl;
  string brief_str_label2 = "[{21,1.64,Randy}(30,1.34,Sarah)]";
  std::cout << "Label after making into a branch" << std::endl;
  std::cout << labelAll.get() << std::endl;
  std::cout << full_str_label2 << std::endl;
  BOOST_CHECK(full_str_label2.compare(labelAll.get()) == 0);
  BOOST_CHECK(brief_str_label2.compare(labelAll.get(LabelType::concise)) == 0);*/
}

BOOST_AUTO_TEST_SUITE_END()
