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

#define BOOST_TEST_MODULE qmnblist_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/qmnblist.h>

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(qmnblist_test)

BOOST_AUTO_TEST_CASE(constructors_test) { QMNBList qmnb; }

BOOST_AUTO_TEST_CASE(add_pair) {

  Segment seg1("one", 1);
  Atom atm1(1, "C", Eigen::Vector3d::Ones());
  seg1.push_back(atm1);

  Segment seg2("two", 2);
  Atom atm2(2, "C", -Eigen::Vector3d::Ones());
  seg2.push_back(atm2);

  Segment seg3("three", 3);
  Atom atm3(3, "Ca", -Eigen::Vector3d::UnitX());
  seg3.push_back(atm3);

  QMNBList qmnb;
  BOOST_CHECK_EQUAL(qmnb.empty(), true);
  qmnb.Add(seg1, seg3, 0.5 * Eigen::Vector3d::Ones());
  qmnb.Add(seg1, seg2, -0.5 * Eigen::Vector3d::Ones());
  BOOST_CHECK_EQUAL(qmnb.size(), 2);
  BOOST_CHECK_EQUAL(qmnb.empty(), false);

  QMPair* p0 = qmnb.FindPair(&seg1, &seg2);

  QMPair* p1 = qmnb.FindPair(&seg1, &seg3);
  BOOST_CHECK_EQUAL(p0->R().isApprox(-0.5 * Eigen::Vector3d::Ones(), 1e-5),
                    true);
  BOOST_CHECK_EQUAL(p0->getId(), 1);
  BOOST_CHECK_EQUAL(p1->getId(), 0);

  // sort qmpairs by seg1id and then by seg2id then reindex the pair id
  // according to that.
  qmnb.sortAndReindex([](QMPair* a, QMPair* b) {
    if (a->Seg1()->getId() != b->Seg1()->getId()) {
      return a->Seg1()->getId() < b->Seg1()->getId();
    }
    return a->Seg2()->getId() < b->Seg2()->getId();
  });

  BOOST_CHECK_EQUAL(p0->getId(), 0);
  BOOST_CHECK_EQUAL(p1->getId(), 1);

  votca::csg::PairList<const Segment*, QMPair>::partners* Partners =
      qmnb.FindPartners(&seg1);

  QMPair* p_ref0 = Partners->at(&seg2);
  QMPair* p_ref1 = Partners->at(&seg3);
  BOOST_CHECK_EQUAL(p_ref0, p0);
  BOOST_CHECK_EQUAL(p_ref1, p1);
}

BOOST_AUTO_TEST_SUITE_END()
