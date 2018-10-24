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

#define BOOST_TEST_MODULE qmmolecule_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/qmmolecule.h>


using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(qmmolecule_test)

BOOST_AUTO_TEST_CASE(constructors_test) { QMMolecule seg("seg1",1); }

BOOST_AUTO_TEST_CASE(load_xyz_test) {
  ofstream xyzfile("molecule.xyz");
  xyzfile << " 5" << endl;
  xyzfile << " methane" << endl;
  xyzfile << " C            .000000     .000000     .000000" << endl;
  xyzfile << " H            .629118     .629118     .629118" << endl;
  xyzfile << " H           -.629118    -.629118     .629118" << endl;
  xyzfile << " H            .629118    -.629118    -.629118" << endl;
  xyzfile << " H           -.629118     .629118    -.629118" << endl;
  xyzfile.close();
  
  QMMolecule seg( "seg1",1);
  seg.LoadFromXYZ("molecule.xyz");
  
  
 
}



BOOST_AUTO_TEST_SUITE_END()
