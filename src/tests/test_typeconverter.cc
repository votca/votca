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
#define BOOST_TEST_MODULE typeconverter_test

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <votca/tools/typeconverter.h>

using namespace std;
using namespace votca::tools;

// Create test classes
class Person {
 private:
  int phone_num_;
  string name_;
  string address_;
  int age_;
  string favorite_color_;
  double height_;

 public:
  Person(){};
  Person(int phone_num, string name, string address, int age,
         string favorite_color, double height)
      : phone_num_(phone_num),
        name_(name),
        address_(address),
        age_(age),
        favorite_color_(favorite_color),
        height_(height){};

  /////////////////////////////////////////////////////////////////////////////
  // Required to use the converter template                                  //
  /////////////////////////////////////////////////////////////////////////////
  static unordered_map<string, vector<string>> writeStr(Person p) {
    unordered_map<string, vector<string>> str_map;
    str_map["Person"].push_back(p.name_);
    str_map["Person"].push_back(p.address_);
    str_map["Person"].push_back(p.favorite_color_);
    return str_map;
  }
  static unordered_map<string, vector<int>> writeInt(Person p) {
    unordered_map<string, vector<int>> int_map;
    int_map["Person"].push_back(p.phone_num_);
    int_map["Person"].push_back(p.age_);
    return int_map;
  }
  static unordered_map<string, vector<double>> writeDouble(Person p) {
    unordered_map<string, vector<double>> double_map;
    double_map["Person"].push_back(p.height_);
    return double_map;
  }
  /////////////////////////////////////////////////////////////////////////////
};

class ContactInfo {
 private:
  int phone_num_;
  string name_;
  string address_;

 public:
  ContactInfo(){};
  ~ContactInfo(){};
  ContactInfo(int phone_num, string name, string address)
      : phone_num_(phone_num), name_(name), address_(address){};
  int getPhoneNum(void) { return phone_num_; }
  string getName(void) { return name_; }
  string getAddress(void) { return address_; }
  /////////////////////////////////////////////////////////////////////////////
  // Required to use the converter template                                  //
  /////////////////////////////////////////////////////////////////////////////
  static void readData(unordered_map<string, vector<int>> int_map,
                       unordered_map<string, vector<double>> double_map,
                       unordered_map<string, vector<string>> str_map,
                       ContactInfo& cont_info) {

    cont_info.phone_num_ = int_map["Person"].at(0);
    cont_info.name_ = str_map["Person"].at(0);
    cont_info.address_ = str_map["Person"].at(1);
  }
  /////////////////////////////////////////////////////////////////////////////
};

BOOST_AUTO_TEST_SUITE(typeconverter_test)

BOOST_AUTO_TEST_CASE(constructors_test) {
  TypeConverter<Person, ContactInfo> converter;
}

BOOST_AUTO_TEST_CASE(import_test) {
  int num = 87094;
  string name = "John Doe";
  double height = 167.8;
  string address = "12 Koogler St, Dayton, OH 32345";
  int age = 24;
  string fav_col = "mauve";
  Person John(num, name, address, age, fav_col, height);

  TypeConverter<Person, ContactInfo> converter;
  converter.importData(John);
}

BOOST_AUTO_TEST_CASE(export_test) {
  int num = 87094;
  string name = "John Doe";
  double height = 167.8;
  string address = "12 Koogler St, Dayton, OH 32345";
  int age = 24;
  string fav_col = "mauve";
  Person John(num, name, address, age, fav_col, height);

  TypeConverter<Person, ContactInfo> converter;

  converter.importData(John);
  ContactInfo cont_info;
  converter.exportData(cont_info);

  BOOST_CHECK_EQUAL(cont_info.getPhoneNum(), num);
  BOOST_CHECK_EQUAL(cont_info.getName(), name);
  BOOST_CHECK_EQUAL(cont_info.getAddress(), address);
}

BOOST_AUTO_TEST_SUITE_END()
