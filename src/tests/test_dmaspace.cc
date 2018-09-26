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

#define BOOST_TEST_MODULE dmaspace_test
#include <vector>

#include <boost/test/unit_test.hpp>
#include <votca/xtp/dmaspace.h>
#include <votca/tools/vec.h>

using namespace std;
using namespace votca::tools;
using namespace votca::DMA;

BOOST_AUTO_TEST_SUITE(dmaspace_test)

BOOST_AUTO_TEST_CASE(constructors_test) { 
  vec r(1.0,2.0,3.0);
  RegularSphericalHarmonics RegSph(r);
  BOOST_CHECK_EQUAL(static_cast<int>(RegSph.R00()),1);
  BOOST_CHECK_EQUAL(static_cast<int>(RegSph.R10()),3);
  BOOST_CHECK_EQUAL(static_cast<int>(RegSph.R11c()),1);
  BOOST_CHECK_EQUAL(static_cast<int>(RegSph.R11s()),2);

  vector<double> r2 = {1.0, 2.0, 3.0,4.0,5.0,6.0,7.0,8.0,9.0};
  ComplexSphericalMoments CompSph(r2);
  
}

BOOST_AUTO_TEST_SUITE_END()
