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
#define BOOST_TEST_MODULE structureparameters_test

#include "../../include/votca/tools/structureparameters.h"
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(structureparameters_test)

BOOST_AUTO_TEST_CASE(constructors_test) { StructureParameters parameters; }

BOOST_AUTO_TEST_CASE(set_test) {
  StructureParameters parameters;

  double mass = 2.0;
  parameters.set(StructureParameter::CSG_Mass, mass);
  string element_type = "C";
  parameters.set(StructureParameter::Element, element_type);
}

BOOST_AUTO_TEST_CASE(get_test) {
  StructureParameters parameters;

  double mass = 2.0;
  parameters.set(StructureParameter::CSG_Mass, mass);
  double mass_check = parameters.get<double>(StructureParameter::CSG_Mass);
  BOOST_CHECK_EQUAL(mass, mass_check);

  string element_type = "C";
  parameters.set(StructureParameter::Element, element_type);
  string element_type_check =
      parameters.get<string>(StructureParameter::Element);
  BOOST_CHECK_EQUAL(element_type, element_type_check);

  int molecule_id = 201;
  parameters.set(StructureParameter::MoleculeId, molecule_id);
  int molecule_id_check = parameters.get<int>(StructureParameter::MoleculeId);
  BOOST_CHECK_EQUAL(molecule_id, molecule_id_check);
}

BOOST_AUTO_TEST_SUITE_END()
