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

#include <boost/test/tools/old/interface.hpp>
#include <stdexcept>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE name_test

// Standard includes
#include <exception>
#include <string>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/tools/objectfactory.h"

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(objectfactory_test)

class base {

 public:
  base() = default;
  virtual ~base() = default;

  virtual std::string identify() = 0;
};

class A : public base {
 public:
  std::string identify() override { return "A"; }
};

class B : public base {
 public:
  std::string identify() override { return "B"; }
};

BOOST_AUTO_TEST_CASE(construction_test) {

  ObjectFactory<std::string, base> factory;

  factory.Register<A>("A");
  factory.Register<B>("B");

  std::unique_ptr<base> A = factory.Create("A");
  BOOST_REQUIRE_EQUAL(A->identify(), "A");

  BOOST_REQUIRE(factory.IsRegistered("A"));

  std::unique_ptr<base> B = factory.Create("B");
  BOOST_REQUIRE_EQUAL(B->identify(), "B");
  std::vector<std::string> keys = factory.getKeys();
  bool A_present = (std::find(keys.begin(), keys.end(), "A") != keys.end());
  bool B_present = (std::find(keys.begin(), keys.end(), "B") != keys.end());

  BOOST_REQUIRE(A_present);
  BOOST_REQUIRE(B_present);
  BOOST_REQUIRE_THROW(factory.Create("C"), std::runtime_error);
  BOOST_REQUIRE_EQUAL(factory.IsRegistered("D"), false);
}

BOOST_AUTO_TEST_SUITE_END()
