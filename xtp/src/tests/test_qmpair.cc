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

#define BOOST_TEST_MODULE qmpair_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include <votca/xtp/qmpair.h>

using namespace votca::tools;

using namespace votca::xtp;

BOOST_AUTO_TEST_SUITE(qmpair_test)

BOOST_AUTO_TEST_CASE(getters_test) {

  Segment seg("one", 1);
  Atom atm1(1, "C", Eigen::Vector3d::Ones());
  seg.push_back(atm1);

  Segment seg2("two", 2);
  Atom atm2(2, "C", -Eigen::Vector3d::Ones());
  seg2.push_back(atm2);

  QMPair pair(0, &seg, &seg2, 0.5 * Eigen::Vector3d::Ones());

  pair.setType(QMPair::PairType::Hopping);

  BOOST_CHECK_EQUAL(pair.getType(), QMPair::Hopping);
}

BOOST_AUTO_TEST_CASE(state_storage_test) {
  Segment seg1("one", 0);
  seg1.push_back(Atom(1, "C", Eigen::Vector3d::Zero()));

  Segment seg2("two", 1);
  seg2.push_back(Atom(2, "C", Eigen::Vector3d::Ones()));

  QMPair pair(7, &seg1, &seg2, Eigen::Vector3d(1.0, 2.0, 3.0));

  pair.setJeff(+0.11, QMStateType::Electron);
  pair.setJeff(-0.22, QMStateType::Hole);
  pair.setJeff(+0.33, QMStateType::Singlet);
  pair.setJeff(-0.44, QMStateType::Triplet);

  pair.setJeff2(1.1, QMStateType::Electron);
  pair.setJeff2(2.2, QMStateType::Hole);
  pair.setJeff2(3.3, QMStateType::Singlet);
  pair.setJeff2(4.4, QMStateType::Triplet);

  pair.setLambdaO(5.1, QMStateType::Electron);
  pair.setLambdaO(5.2, QMStateType::Hole);
  pair.setLambdaO(5.3, QMStateType::Singlet);
  pair.setLambdaO(5.4, QMStateType::Triplet);

  BOOST_CHECK_CLOSE(pair.getJeff(QMStateType::Electron), 0.11, 1e-12);
  BOOST_CHECK_CLOSE(pair.getJeff(QMStateType::Hole), -0.22, 1e-12);
  BOOST_CHECK_CLOSE(pair.getJeff(QMStateType::Singlet), 0.33, 1e-12);
  BOOST_CHECK_CLOSE(pair.getJeff(QMStateType::Triplet), -0.44, 1e-12);

  BOOST_CHECK_CLOSE(pair.getJeff2(QMStateType::Electron), 1.1, 1e-12);
  BOOST_CHECK_CLOSE(pair.getJeff2(QMStateType::Hole), 2.2, 1e-12);
  BOOST_CHECK_CLOSE(pair.getJeff2(QMStateType::Singlet), 3.3, 1e-12);
  BOOST_CHECK_CLOSE(pair.getJeff2(QMStateType::Triplet), 4.4, 1e-12);

  BOOST_CHECK_CLOSE(pair.getLambdaO(QMStateType::Electron), 5.1, 1e-12);
  BOOST_CHECK_CLOSE(pair.getLambdaO(QMStateType::Hole), 5.2, 1e-12);
  BOOST_CHECK_CLOSE(pair.getLambdaO(QMStateType::Singlet), 5.3, 1e-12);
  BOOST_CHECK_CLOSE(pair.getLambdaO(QMStateType::Triplet), 5.4, 1e-12);
}

BOOST_AUTO_TEST_CASE(identity_and_geometry_test) {
  Segment seg1("one", 3);
  seg1.push_back(Atom(1, "C", Eigen::Vector3d::Zero()));

  Segment seg2("two", 4);
  seg2.push_back(Atom(2, "C", Eigen::Vector3d(2.0, 0.0, 0.0)));

  Eigen::Vector3d dr(1.0, 2.0, 2.0);
  QMPair pair(17, &seg1, &seg2, dr);

  BOOST_CHECK_EQUAL(pair.getId(), 17);
  BOOST_CHECK(pair.R().isApprox(dr));
  BOOST_CHECK_CLOSE(pair.Dist(), 3.0, 1e-12);

  BOOST_CHECK_EQUAL(pair.Seg1()->getId(), 3);
  BOOST_CHECK_EQUAL(pair.Seg2()->getId(), 4);
  BOOST_CHECK_EQUAL(pair.first()->getId(), 3);
  BOOST_CHECK_EQUAL(pair.second()->getId(), 4);

  pair.setType(QMPair::Excitoncl);
  BOOST_CHECK_EQUAL(pair.getType(), QMPair::Excitoncl);
}

BOOST_AUTO_TEST_CASE(derived_energy_quantities_test) {
  Segment seg1("one", 0);
  seg1.push_back(Atom(1, "C", Eigen::Vector3d::Zero()));
  seg1.setU_nX_nN(0.10, QMStateType::Hole);
  seg1.setU_xN_xX(0.20, QMStateType::Hole);
  seg1.setU_xX_nN(1.00, QMStateType::Hole);
  seg1.setEMpoles(QMStateType::Hole, 2.00);

  Segment seg2("two", 1);
  seg2.push_back(Atom(2, "C", Eigen::Vector3d::Ones()));
  seg2.setU_nX_nN(0.30, QMStateType::Hole);
  seg2.setU_xN_xX(0.40, QMStateType::Hole);
  seg2.setU_xX_nN(1.50, QMStateType::Hole);
  seg2.setEMpoles(QMStateType::Hole, 3.00);

  QMPair pair(0, &seg1, &seg2, Eigen::Vector3d::Zero());

  BOOST_CHECK_CLOSE(pair.getReorg12(QMStateType::Hole), 0.10 + 0.40, 1e-12);
  BOOST_CHECK_CLOSE(pair.getReorg21(QMStateType::Hole), 0.20 + 0.30, 1e-12);

  // site energy = EMpoles + U_xX_nN
  BOOST_CHECK_CLOSE(pair.getdE12(QMStateType::Hole),
                    (2.00 + 1.00) - (3.00 + 1.50), 1e-12);
}

BOOST_AUTO_TEST_CASE(write_read_data_roundtrip_test) {
  Segment seg1("one", 0);
  seg1.push_back(Atom(1, "C", Eigen::Vector3d::Zero()));

  Segment seg2("two", 1);
  seg2.push_back(Atom(2, "C", Eigen::Vector3d::Ones()));

  QMPair pair(42, &seg1, &seg2, Eigen::Vector3d(0.5, -1.5, 2.5));
  pair.setType(QMPair::Excitoncl);

  pair.setLambdaO(1.1, QMStateType::Electron);
  pair.setLambdaO(1.2, QMStateType::Hole);
  pair.setLambdaO(1.3, QMStateType::Singlet);
  pair.setLambdaO(1.4, QMStateType::Triplet);

  pair.setJeff(2.1, QMStateType::Electron);
  pair.setJeff(-2.2, QMStateType::Hole);
  pair.setJeff(2.3, QMStateType::Singlet);
  pair.setJeff(-2.4, QMStateType::Triplet);

  pair.setJeff2(3.1, QMStateType::Electron);
  pair.setJeff2(3.2, QMStateType::Hole);
  pair.setJeff2(3.3, QMStateType::Singlet);
  pair.setJeff2(3.4, QMStateType::Triplet);

  QMPair::data d{};
  pair.WriteData(d);

  std::vector<Segment> segments;
  segments.push_back(seg1);
  segments.push_back(seg2);

  QMPair restored(d, segments);

  BOOST_CHECK_EQUAL(restored.getId(), 42);
  BOOST_CHECK(restored.R().isApprox(Eigen::Vector3d(0.5, -1.5, 2.5)));
  BOOST_CHECK_EQUAL(restored.getType(), QMPair::Excitoncl);

  BOOST_CHECK_EQUAL(restored.Seg1()->getId(), 0);
  BOOST_CHECK_EQUAL(restored.Seg2()->getId(), 1);

  BOOST_CHECK_CLOSE(restored.getLambdaO(QMStateType::Electron), 1.1, 1e-12);
  BOOST_CHECK_CLOSE(restored.getLambdaO(QMStateType::Hole), 1.2, 1e-12);
  BOOST_CHECK_CLOSE(restored.getLambdaO(QMStateType::Singlet), 1.3, 1e-12);
  BOOST_CHECK_CLOSE(restored.getLambdaO(QMStateType::Triplet), 1.4, 1e-12);

  BOOST_CHECK_CLOSE(restored.getJeff(QMStateType::Electron), 2.1, 1e-12);
  BOOST_CHECK_CLOSE(restored.getJeff(QMStateType::Hole), -2.2, 1e-12);
  BOOST_CHECK_CLOSE(restored.getJeff(QMStateType::Singlet), 2.3, 1e-12);
  BOOST_CHECK_CLOSE(restored.getJeff(QMStateType::Triplet), -2.4, 1e-12);

  BOOST_CHECK_CLOSE(restored.getJeff2(QMStateType::Electron), 3.1, 1e-12);
  BOOST_CHECK_CLOSE(restored.getJeff2(QMStateType::Hole), 3.2, 1e-12);
  BOOST_CHECK_CLOSE(restored.getJeff2(QMStateType::Singlet), 3.3, 1e-12);
  BOOST_CHECK_CLOSE(restored.getJeff2(QMStateType::Triplet), 3.4, 1e-12);
}

BOOST_AUTO_TEST_CASE(pair_type_conversion_test) {
  BOOST_CHECK_EQUAL(QMPair::get_name(QMPair::Hopping), "Hopping");
  BOOST_CHECK_EQUAL(QMPair::get_name(QMPair::Excitoncl), "Excitoncl");

  BOOST_CHECK_EQUAL(QMPair::get_Enum("Hopping"), QMPair::Hopping);
  BOOST_CHECK_EQUAL(QMPair::get_Enum("Excitoncl"), QMPair::Excitoncl);

  BOOST_CHECK_THROW(QMPair::get_Enum("NotAType"), std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
