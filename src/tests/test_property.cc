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

#define BOOST_TEST_MODULE property_test
#include <boost/test/unit_test.hpp>
#include <votca/tools/property.h>
#include <fstream>
#include <sstream>
using namespace votca::tools;

BOOST_AUTO_TEST_SUITE(property_test)

BOOST_AUTO_TEST_CASE(eigen_test) {
  Property prop;
  prop.add("vec", "1 2 3 4 5 6 7 8");
  Eigen::VectorXd result = prop.get("vec").as<Eigen::VectorXd>();

  Eigen::VectorXd ref;
  ref.resize(8);
  ref << 1, 2, 3, 4, 5, 6, 7, 8;
  BOOST_CHECK_EQUAL(ref.isApprox(result, 0.0001), true);

  BOOST_CHECK_THROW(prop.get("vec").as<Eigen::Vector3d>(), std::runtime_error);

  Property prop2;
  prop.add("vec", "1 2 3");
  Eigen::Vector3d result2 = prop.get("vec").as<Eigen::Vector3d>();
  Eigen::Vector3d ref2{1, 2, 3};
  BOOST_CHECK_EQUAL(ref2.isApprox(result2, 0.0001), true);
}


BOOST_AUTO_TEST_CASE(readin){
  std::ofstream xmlfile("notnormalized.xml");
  xmlfile << "<basis name=\"def2-TZVP\">" << std::endl;
  xmlfile << "  <element name=\"Al\">" << std::endl;
  xmlfile << "    <shell scale=\"1.0\" type=\"D\">" << std::endl;
  xmlfile << "      <constant decay=\"1.570000e+00\">" << std::endl;
  xmlfile << "        <contractions factor=\"2.000000e-01\" type=\"D\"/>"
            << std::endl;
  xmlfile << "      </constant>" << std::endl;
  xmlfile << "      <constant decay=\"3.330000e-01\">" << std::endl;
  xmlfile << "        <contractions factor=\"1.000000e+00\" type=\"D\"/>"
            << std::endl;
  xmlfile << "      </constant>" << std::endl;
  xmlfile << "    </shell> " << std::endl;
  xmlfile << "  </element> " << std::endl;
  xmlfile << "</basis> " << std::endl;
  xmlfile.close();
  
  Property prop;
  prop.LoadFromXML("notnormalized.xml");
  
  Property& e=prop.get("basis.element");
  BOOST_CHECK_EQUAL(e.getAttribute<std::string>("name"),"Al");
  
  BOOST_REQUIRE_THROW(prop.get("basis.blabla"),std::runtime_error);
    
  Property& c=prop.get("basis.element.shell.constant.contractions");
  BOOST_CHECK_EQUAL(c.getAttribute<std::string>("type"),"D");
  
  BOOST_REQUIRE_THROW(c.getAttribute<std::string>("bla"),std::runtime_error);
    
}



BOOST_AUTO_TEST_CASE(printtostream){
    
   Property prop;
   Property& p=prop.add("hello","");
   Property& pp=p.add("hi","5");
   Property& pd=p.add("ho","bumm");
   
   std::stringstream dudu;
   dudu<<prop;
   std::string printout=dudu.str();
   
   std::string ref="<>\n"
"	<hello>\n"
"		<hi>5</hi>\n"
"		<ho>bumm</ho>\n"
"	</hello>\n"
"</>\n";
   
   BOOST_CHECK_EQUAL(printout,ref);
   
}

BOOST_AUTO_TEST_SUITE_END()
