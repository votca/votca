/*
 * Copyright 2009-2020 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE topology_test

// Third party includes
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/csg/boundarycondition.h>

// Local VOTCA includes
#include "votca/xtp/atom.h"
#include "votca/xtp/qmpair.h"
#include "votca/xtp/segment.h"
#include "votca/xtp/topology.h"

using namespace votca::tools;

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(topology_test)

BOOST_AUTO_TEST_CASE(constructors_test) { Topology top; }

BOOST_AUTO_TEST_CASE(box_test) {

  Eigen::Matrix3d box = 2 * Eigen::Matrix3d::Identity();
  Topology top;
  top.setBox(box);

  auto vol = top.BoxVolume();
  BOOST_CHECK_CLOSE(vol, 8, 0.0001);
  auto box2 = top.getBox();

  BOOST_CHECK_EQUAL(box.isApprox(box2, 1e-6), true);
}

BOOST_AUTO_TEST_CASE(simple_test) {

  Topology top;
  top.setStep(1);
  BOOST_CHECK_EQUAL(top.getStep(), 1);
  top.setTime(1.21);
  BOOST_CHECK_CLOSE(top.getTime(), 1.21, 0.0001);
}

BOOST_AUTO_TEST_CASE(copy_keeps_pair_data) {
  Topology top;
  top.setBox(Eigen::Matrix3d::Identity() * 10.0);
  for (Index i = 0; i < 3; ++i) {
    Segment& seg = top.AddSegment("s");
    seg.push_back(Atom(0, "C", Eigen::Vector3d(double(i), 0.0, 0.0)));
  }
  QMPair& a = top.NBList().Add(top.Segments()[0], top.Segments()[1],
                               Eigen::Vector3d(1.0, 0.0, 0.0));
  a.setJeff2(1e-4, QMStateType::Hole);
  a.setJeff2(2e-4, QMStateType::Electron);
  a.setJeff(3e-3, QMStateType::Singlet);
  a.setLambdaO(0.05, QMStateType::Hole);
  QMPair& b = top.NBList().Add(top.Segments()[1], top.Segments()[2],
                               Eigen::Vector3d(1.0, 0.0, 0.0));
  b.setType(QMPair::Excitoncl);
  b.setJeff2(5e-5, QMStateType::Triplet);

  auto check = [&top](const Topology& copy) {
    BOOST_REQUIRE_EQUAL(copy.NBList().size(), top.NBList().size());
    for (Index p = 0; p < top.NBList().size(); ++p) {
      const QMPair& orig = *top.NBList()[p];
      const QMPair& c = *copy.NBList()[p];
      BOOST_CHECK_EQUAL(c.getId(), orig.getId());
      BOOST_CHECK_EQUAL(c.getType(), orig.getType());
      BOOST_CHECK(c.R() == orig.R());
      // the copy points to its own segments
      BOOST_CHECK(c.Seg1() == &copy.Segments()[orig.Seg1()->getId()]);
      BOOST_CHECK(c.Seg2() == &copy.Segments()[orig.Seg2()->getId()]);
      for (QMStateType s :
           {QMStateType(QMStateType::Electron), QMStateType(QMStateType::Hole),
            QMStateType(QMStateType::Singlet),
            QMStateType(QMStateType::Triplet)}) {
        BOOST_CHECK_EQUAL(c.getJeff2(s), orig.getJeff2(s));
        BOOST_CHECK_EQUAL(c.getJeff(s), orig.getJeff(s));
        BOOST_CHECK_EQUAL(c.getLambdaO(s), orig.getLambdaO(s));
      }
    }
  };
  Topology copied(top);
  check(copied);
  Topology assigned;
  assigned.AddSegment("x");
  assigned = top;
  check(assigned);
}

BOOST_AUTO_TEST_SUITE_END()
