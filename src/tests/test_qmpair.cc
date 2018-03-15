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

#define BOOST_TEST_MODULE qmpair_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/qmpair.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(qmpair_test)

BOOST_AUTO_TEST_CASE(constructors_test) { QMPair qm_p; }

BOOST_AUTO_TEST_CASE(getters_test) { 
  QMPair qm_p; 
  BOOST_CHECK_EQUAL(qm_p.getId(),-1);
  auto r = qm_p.R();
  BOOST_CHECK_EQUAL(r.x(),0.0);
  BOOST_CHECK_EQUAL(r.y(),0.0);
  BOOST_CHECK_EQUAL(r.z(),0.0);
  BOOST_CHECK_EQUAL(qm_p.hasGhost(),false);

  int state = -1;
  BOOST_CHECK_EQUAL(qm_p.getLambda0(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate21(state),0.0);
  BOOST_CHEKC_EQUAL(qm_p.getJeff2(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getdE12(state),0.0);
  
  state = 1;
  BOOST_CHECK_EQUAL(qm_p.getLambda0(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate21(state),0.0);
  BOOST_CHEKC_EQUAL(qm_p.getJeff2(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getdE12(state),0.0);

  int state = 2;
  BOOST_CHECK_EQUAL(qm_p.getLambda0(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate21(state),0.0);
  BOOST_CHEKC_EQUAL(qm_p.getJeff2(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getdE12(state),0.0);
  
  int state = 3;
  BOOST_CHECK_EQUAL(qm_p.getLambda0(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg12_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getReorg21_x(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate12(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getRate21(state),0.0);
  BOOST_CHEKC_EQUAL(qm_p.getJeff2(state),0.0);
  BOOST_CHECK_EQUAL(qm_p.getdE12(state),0.0);
  
}

BOOST_AUTO_TEST_SUITE_END()
