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
#include <libint2/initialize.h>
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE multipole_interactions_test

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/bgsite.h"
#include "votca/xtp/dipoledipoleinteraction.h"
#include "votca/xtp/eeinteractor.h"
#include "votca/xtp/mpinteractions.h"

#include "votca/xtp/eigen.h"

using namespace votca::xtp;
using namespace votca;

BOOST_AUTO_TEST_SUITE(multipole_interactions_test)

BOOST_AUTO_TEST_CASE(field_calculations) {
  Vector9d mpoles1;
  mpoles1.Random();
  Vector9d mpoles2;
  mpoles2.Random();

  StaticSegment seg1("one", 1);
  PolarSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setMultipole(mpoles1, 2);
  PolarSite two(2, "H");
  two.setPos(Eigen::Vector3d(2, 1, 2.1));
  two.setMultipole(mpoles2, 2);

  seg1.push_back(one);
  seg2.push_back(two);

  BGSite bg_one(one);
  BGSite bg_two(two);

  // original spherical implementation for field calculations
  eeInteractor interactor;
  interactor.ApplyStaticField<StaticSegment, Estatic::V>(seg1, seg2);
  // cartesian implementation for the the field calculation
  MPField<Screening::none, 2> mpfield;
  Eigen::Vector3d cartField = mpfield.fieldAtBy(bg_two, bg_one);

  // should be more or less equal
  BOOST_CHECK((cartField - seg2[0].V()).norm() < 1e-8);
}

BOOST_AUTO_TEST_CASE(energy_calculations) {
  Vector9d mpoles1;
  Vector9d mpoles2;
  mpoles2.Random();
  
  mpoles1.Random();
  StaticSegment seg1("one", 1);
  PolarSegment seg2("two", 2);
  StaticSite one(1, "H");
  one.setPos(Eigen::Vector3d::Zero());
  one.setMultipole(mpoles1, 2);
  PolarSite two(2, "H");
  two.setPos(Eigen::Vector3d(2, 1, 2.1));
  two.setMultipole(mpoles2, 2);

  seg1.push_back(one);
  seg2.push_back(two);

  BGSite bg_one(one);
  BGSite bg_two(two);

  // original spherical implementation for energy evaluation
  eeInteractor interactor;
  double sphericalEnergy =
      interactor.CalcStaticEnergy<StaticSegment, PolarSegment>(seg1, seg2);
  // cartesian implementation for the energy evaluation
  MPEnergy<Screening::none, 2> mpenergy;
  double cartesianEnergy = mpenergy.energy(bg_two, bg_one);

  std::cout << sphericalEnergy << std::endl;
  std::cout << cartesianEnergy << std::endl;

  // should be more or less equal
  BOOST_CHECK_CLOSE(cartesianEnergy, sphericalEnergy, 1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
