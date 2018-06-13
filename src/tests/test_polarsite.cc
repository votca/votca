/*
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
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

#define BOOST_TEST_MODULE polararsite_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/polarsite.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(polararsite_test)

BOOST_AUTO_TEST_CASE(constructors_test) { PolarSite ps(1, "ps1"); }

BOOST_AUTO_TEST_CASE(getters_test) {
  PolarSite ps(1,"ps2");
  BOOST_CHECK_EQUAL(ps.getId(),1);
  BOOST_CHECK_EQUAL(ps.getName(),"ps2");
}

BOOST_AUTO_TEST_CASE(getters_test) {
  PolarSite ps(1,"ps2");
  BOOST_CHECK_EQUAL(ps.getId(),1);
  BOOST_CHECK_EQUAL(ps.getName(),"ps2");
}

BOOST_AUTO_TEST_CASE(multipole_test) {
  PolarSite ps(1,"ps2");
  Eigen::VectorXd multipoles=Eigen::VectorXd::Zero(9);
  multipoles<<1,2,3,4,8,7,2,3.3,-0.5;
  ps.setMultipoles(multipoles);
  bool check_mpoles=multipoles.isApprox(ps.getMultipoles(),0.0001);
   BOOST_CHECK_EQUAL(check_mpoles,true);
   
   bool check_rank=(ps.getRank()==2);
   BOOST_CHECK_EQUAL(check_rank,true);
  
  
}


BOOST_AUTO_TEST_SUITE_END()
