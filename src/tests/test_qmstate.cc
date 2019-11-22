/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE qmstate_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/qmstate.h>

#include <fstream>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(qmstate_test)

BOOST_AUTO_TEST_CASE(QMStatetype_test) {
  QMStateType type;
  type.FromString("S");

  BOOST_CHECK_EQUAL(type.Type() == QMStateType::Singlet, true);
  BOOST_CHECK_EQUAL(type.ToString(), "s");
  BOOST_CHECK_EQUAL(type.ToLongString(), "singlet");

  QMStateType type2;
  type2.FromString("Singlet");
  BOOST_CHECK_EQUAL(type2.Type() == QMStateType::Singlet, true);

  BOOST_CHECK_EQUAL(type == type2, true);

  QMStateType type3;
  type3.FromString("N");
  BOOST_CHECK_EQUAL(type3 == type2, false);
  BOOST_CHECK_EQUAL(type3 == QMStateType::Gstate, true);

  QMStateType type4 = QMStateType(QMStateType::KSstate);
  BOOST_CHECK_EQUAL(type4 == QMStateType::KSstate, true);
}

BOOST_AUTO_TEST_CASE(QMState_test) {
  QMState state;
  state.FromString("S1");
  BOOST_CHECK_EQUAL(state.Type() == QMStateType::Singlet, true);
  BOOST_CHECK_EQUAL(state.StateIdx(), 0);

  QMState state2;
  state.FromString("N");
  BOOST_CHECK_EQUAL(state.StateIdx(), -1);

  QMState state3;
  state3.FromString("n2s12");
  BOOST_CHECK_EQUAL(state3.isTransition(), true);

  BOOST_CHECK_EQUAL(state3.Type() == QMStateType::Singlet, true);
  BOOST_CHECK_EQUAL(state3.StateIdx(), 11);

  QMState state4;
  state4.FromString("n2s16");
  QMState state5;
  state5.FromString("groundstate to singlet 16");

  BOOST_CHECK_EQUAL(state4 == state5, true);
  std::string result = state4.ToLongString();
  BOOST_CHECK_EQUAL(result, "Groundstate to singlet 16");

  QMState state6 = QMState(QMStateType::Singlet, 15, true);
  BOOST_CHECK_EQUAL(state6 == state5, true);

  QMState hole;
  hole.FromString("h1");
  BOOST_CHECK_EQUAL(hole.StateIdx(), 0);
  QMState hole2;
  hole2.FromString("h");
  BOOST_CHECK_EQUAL(hole2.StateIdx(), 0);
  QMState electron;
  electron.FromString("e1");
  BOOST_CHECK_EQUAL(electron.StateIdx(), 0);
  QMState electron2;
  electron2.FromString("e");
  BOOST_CHECK_EQUAL(electron2.StateIdx(), 0);
}

BOOST_AUTO_TEST_SUITE_END()
