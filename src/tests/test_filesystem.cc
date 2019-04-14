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

#include <boost/test/unit_test.hpp>
#include <fstream>
#include <votca/tools/filesystem.h>
using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(filesystem_test)

BOOST_AUTO_TEST_CASE(FileExists) {

  std::string filename = "blablablablablabla.txt";
  BOOST_CHECK_EQUAL(filesystem::FileExists(filename), false);
  std::string filename2 = "test_exist.txt";
  std::ofstream outfile(filename2);

  outfile << "my text here!" << std::endl;

  outfile.close();
  BOOST_CHECK_EQUAL(filesystem::FileExists(filename2), true);
}

BOOST_AUTO_TEST_CASE(GetFileExtension) {

  std::string filename = "bla.txt";
  BOOST_CHECK_EQUAL(filesystem::GetFileExtension(filename), "txt");

  std::string filename2 = "blubb";
  BOOST_CHECK_EQUAL(filesystem::GetFileExtension(filename2), "");

  std::string filename3 = "a.b.c.d.gro";
  BOOST_CHECK_EQUAL(filesystem::GetFileExtension(filename3), "gro");
}

BOOST_AUTO_TEST_CASE(GetFileBase) {

  std::string filename = "bla.txt";
  BOOST_CHECK_EQUAL(filesystem::GetFileBase(filename), "bla");

  std::string filename2 = "blubb";
  BOOST_CHECK_EQUAL(filesystem::GetFileBase(filename2), "blubb");

  std::string filename3 = "a.b.c.d.gro";
  BOOST_CHECK_EQUAL(filesystem::GetFileBase(filename3), "a.b.c.d");
}

BOOST_AUTO_TEST_SUITE_END()
