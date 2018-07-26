/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE datacollection_test
#include <boost/test/unit_test.hpp>
#include <iostream>
#include <vector>

#include <votca/tools/datacollection.h>

using namespace std;
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(datacollection_test)

BOOST_AUTO_TEST_CASE(constructors_test) { DataCollection<double> datacollection; }

BOOST_AUTO_TEST_CASE(test) {

	DataCollection<int> datacollection;

	string name_tag1 = "x positions";
	string name_tag2 = "y positions";
	string name_tag3 = "z positions";

	auto xpositions = datacollection.CreateArray(name_tag1);
	auto ypositions = datacollection.CreateArray(name_tag2);
	auto zpositions = datacollection.CreateArray(name_tag3);

	bool x_compare_name = !name_tag1.compare(xpositions->getName()); 
	bool y_compare_name = !name_tag2.compare(ypositions->getName()); 
	bool z_compare_name = !name_tag3.compare(zpositions->getName()); 
	
	BOOST_CHECK(x_compare_name);
	BOOST_CHECK(y_compare_name);
	BOOST_CHECK(z_compare_name);

	xpositions->push_back(1);
	xpositions->push_back(2);
	xpositions->push_back(3);
	
	ypositions->push_back(-4);
	ypositions->push_back(4);
	ypositions->push_back(4);

	zpositions->push_back(5);
	zpositions->push_back(-5);
	zpositions->push_back(5);

	x_compare_name = !name_tag1.compare(datacollection.Data().at(0)->getName());
	y_compare_name = !name_tag2.compare(datacollection.Data().at(1)->getName());
	z_compare_name = !name_tag3.compare(datacollection.Data().at(2)->getName());

	BOOST_CHECK(x_compare_name);
	BOOST_CHECK(y_compare_name);
	BOOST_CHECK(z_compare_name);

	auto xPosArray = datacollection.ArrayByName(name_tag1);

	BOOST_CHECK_EQUAL(xPosArray->at(0),1);
	BOOST_CHECK_EQUAL(xPosArray->at(1),2);
	BOOST_CHECK_EQUAL(xPosArray->at(2),3);

	auto yPosArray = datacollection.ArrayByName(name_tag2);

	BOOST_CHECK_EQUAL(yPosArray->at(0),-4);
	BOOST_CHECK_EQUAL(yPosArray->at(1),4);
	BOOST_CHECK_EQUAL(yPosArray->at(2),4);

	auto zPosArray = datacollection.ArrayByName(name_tag3);

	BOOST_CHECK_EQUAL(zPosArray->at(0),5);
	BOOST_CHECK_EQUAL(zPosArray->at(1),-5);
	BOOST_CHECK_EQUAL(zPosArray->at(2),5);

}

BOOST_AUTO_TEST_SUITE_END()
