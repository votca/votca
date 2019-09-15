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

#define BOOST_TEST_MODULE polarsegment_test
#include <boost/test/unit_test.hpp>
#include <votca/xtp/classicalsegment.h>

using namespace votca::xtp;
using namespace std;

BOOST_AUTO_TEST_SUITE(polarsegment_test)

BOOST_AUTO_TEST_CASE(constructors_test) { PolarSegment("seg1", 0); }

BOOST_AUTO_TEST_CASE(load_mps) {

  PolarSegment seg = PolarSegment("seg1", 0);
  std::ofstream mpsfile("polarsite.mps");
  mpsfile << "! Two Sites" << endl;
  mpsfile << "! N=2 " << endl;
  mpsfile << "Units angstrom" << endl;
  mpsfile << "C +0 0 3 Rank 2" << endl;
  mpsfile << "+1" << endl;
  mpsfile << "10 0 0" << endl;
  mpsfile << "100 0 0 0 0" << endl;
  mpsfile
      << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 +1.9445387"
      << endl;
  mpsfile << "H +0 0 1 Rank 0" << endl;
  mpsfile << "-1" << endl;
  mpsfile << "P +1.0000000" << endl;

  seg.LoadFromFile("polarsite.mps");
  Eigen::Vector3d ref_pos =
      Eigen::Vector3d(0, 0, 2.845154333 * votca::tools::conv::ang2bohr);
  bool is_equal = seg.getPos().isApprox(ref_pos, 0.0001);
  BOOST_CHECK_EQUAL(is_equal, true);
  if (!is_equal) {
    std::cout << "result" << std::endl;
    std::cout << seg.getPos() << std::endl;
    std::cout << "reference" << std::endl;
    std::cout << ref_pos << std::endl;
  }

  BOOST_CHECK_EQUAL(seg.size(), 2);
  BOOST_CHECK_EQUAL(seg[0].getRank(), 2);
  BOOST_CHECK_EQUAL(seg[0].getElement(), "C");
  BOOST_CHECK_EQUAL(seg[1].getRank(), 0);
  BOOST_CHECK_EQUAL(seg[1].getElement(), "H");

  Eigen::VectorXd mul_ref = Eigen::VectorXd::Zero(9);
  mul_ref << 1, 0, 0, 10, 100, 0, 0, 0, 0;
  bool multipoles_equal = mul_ref.isApprox(seg[0].Q(), 1e-5);
  if (!multipoles_equal) {
    std::cout << "result" << std::endl;
    std::cout << seg[0].Q() << std::endl;
    std::cout << "reference" << std::endl;
    std::cout << mul_ref << std::endl;
  }
  BOOST_CHECK_EQUAL(multipoles_equal, true);
  Eigen::VectorXd mul_ref2 = Eigen::VectorXd::Zero(9);
  mul_ref2 << -1, 0, 0, 0, 0, 0, 0, 0, 0;
  bool multipoles_equal2 = mul_ref2.isApprox(seg[1].Q(), 1e-5);
  if (!multipoles_equal2) {
    std::cout << "result" << std::endl;
    std::cout << seg[1].Q() << std::endl;
    std::cout << "reference" << std::endl;
    std::cout << mul_ref2 << std::endl;
  }
  BOOST_CHECK_EQUAL(multipoles_equal2, true);

  std::string ref_string =
      "  C +0.0000000 +0.0000000 +3.0000000 Rank 2\n"
      "    +1.0000000\n"
      "    +10.0000000 +0.0000000 +0.0000000\n"
      "    +100.0000000 +0.0000000 +0.0000000 +0.0000000 +0.0000000\n"
      "     P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 "
      "+1.9445387\n";
  bool string_equal = (ref_string == seg[0].WriteMpsLine("angstrom"));
  BOOST_CHECK_EQUAL(string_equal, true);
  if (!string_equal) {
    std::string result = seg[0].WriteMpsLine("angstrom");
    std::cout << "result" << std::endl;
    std::cout << result << std::endl;
    std::cout << "reference" << std::endl;
    std::cout << ref_string << std::endl;
  }
  std::string ref_string2 =
      "  H +0.0000000 +0.0000000 +1.0000000 Rank 0\n"
      "    -1.0000000\n"
      "     P +1.0000000 +0.0000000 +0.0000000 +1.0000000 +0.0000000 "
      "+1.0000000\n";

  bool string_equal2 = (ref_string2 == seg[1].WriteMpsLine("angstrom"));
  BOOST_CHECK_EQUAL(string_equal2, true);
  if (!string_equal2) {
    std::string result = seg[1].WriteMpsLine("angstrom");
    std::cout << "result" << std::endl;
    std::cout << result << std::endl;
    std::cout << "reference" << std::endl;
    std::cout << ref_string2 << std::endl;
  }
}

BOOST_AUTO_TEST_CASE(add_atom_test) {
  PolarSegment seg = PolarSegment("seg1", 0);

  bool check_pos = seg.getPos().isApprox(Eigen::Vector3d::Zero());
  BOOST_CHECK_EQUAL(check_pos, true);
  Eigen::Vector3d pos = Eigen::Vector3d::Ones();
  PolarSite site = PolarSite(0, "C", pos);
  Eigen::VectorXd poles = Vector9d::Ones(9);
  int rank = 2;
  site.setMultipole(poles, rank);
  seg.push_back(site);

  bool check_pos2 = seg.getPos().isApprox(Eigen::Vector3d::Ones());
  BOOST_CHECK_EQUAL(check_pos2, true);
}

BOOST_AUTO_TEST_CASE(readwritehdf) {
  PolarSegment seg = PolarSegment("seg1", 0);
  std::ofstream mpsfile("polarsite.mps");
  mpsfile << "! Two Sites" << endl;
  mpsfile << "! N=2 " << endl;
  mpsfile << "Units angstrom" << endl;
  mpsfile << "C +0 0 3 Rank 2" << endl;
  mpsfile << "+1" << endl;
  mpsfile << "10 0 0" << endl;
  mpsfile << "100 0 0 0 0" << endl;
  mpsfile
      << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 +1.9445387"
      << endl;
  mpsfile << "H +0 0 1 Rank 0" << endl;
  mpsfile << "-1" << endl;
  mpsfile << "P +1.0000000" << endl;

  seg.LoadFromFile("polarsite.mps");

  CheckpointFile ff("polarsegment_test.hdf5");
  CheckpointWriter ww = ff.getWriter();
  seg.WriteToCpt(ww);

  CheckpointReader rr = ff.getReader();
  PolarSegment seg2(rr);
  BOOST_CHECK_EQUAL(seg.size(), seg2.size());
  for (int i = 0; i < seg.size(); i++) {
    BOOST_CHECK_EQUAL(seg2[i].Q().isApprox(seg[i].Q(), 1e-5), true);
    BOOST_CHECK_EQUAL(seg2[i].getPos().isApprox(seg[i].getPos(), 1e-5), true);
    BOOST_CHECK_EQUAL(seg2[i].getElement(), seg[i].getElement());
    BOOST_CHECK_EQUAL(
        seg2[i].getPolarisation().isApprox(seg[i].getPolarisation(), 1e-5),
        true);
  }
}

BOOST_AUTO_TEST_CASE(readwritemps) {
  PolarSegment seg = PolarSegment("seg1", 0);
  std::ofstream mpsfile("polarsite.mps");
  mpsfile << "! Two Sites" << endl;
  mpsfile << "! N=2 " << endl;
  mpsfile << "Units angstrom" << endl;
  mpsfile << "C +0 0 3 Rank 2" << endl;
  mpsfile << "+1" << endl;
  mpsfile << "10 0 0" << endl;
  mpsfile << "100 0 0 0 0" << endl;
  mpsfile
      << "P +1.9445387 +0.0000000 +0.0000000 +1.9445387 +0.0000000 +1.9445387"
      << endl;
  mpsfile << "H +0 0 1 Rank 0" << endl;
  mpsfile << "-1" << endl;
  mpsfile << "P +1.0000000" << endl;

  seg.LoadFromFile("polarsite.mps");
  seg.WriteMPS("test.mps", "test");
  PolarSegment seg2("seg2", 1);
  seg2.LoadFromFile("test.mps");
  BOOST_CHECK_EQUAL(seg.size(), seg2.size());
  for (int i = 0; i < seg.size(); i++) {
    BOOST_CHECK_EQUAL(seg2[i].Q().isApprox(seg[i].Q(), 1e-5), true);
    BOOST_CHECK_EQUAL(seg2[i].getPos().isApprox(seg[i].getPos(), 1e-5), true);
    BOOST_CHECK_EQUAL(seg2[i].getElement(), seg[i].getElement());
    BOOST_CHECK_EQUAL(
        seg2[i].getPolarisation().isApprox(seg[i].getPolarisation(), 1e-5),
        true);
  }
}

BOOST_AUTO_TEST_SUITE_END()
