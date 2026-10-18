/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
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

#define BOOST_TEST_MODULE chargeresponse_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/chargeresponse.h"
#include "votca/xtp/checkpoint.h"
#include "votca/xtp/classicalsegment.h"
#include "votca/xtp/segmentmapper.h"

using namespace votca::xtp;
using namespace votca;

namespace {
// A charge-conserving, positive semidefinite 5x5 matrix for methane
// (C, H, H, H, H) with distinguishable entries, so that a wrong
// permutation shows: Pi = L^T diag(w) L with L the bond incidence
// matrix of the four C-H bonds.
Eigen::MatrixXd MethanePi() {
  Eigen::MatrixXd L = Eigen::MatrixXd::Zero(4, 5);
  for (Index b = 0; b < 4; ++b) {
    L(b, 0) = 1.0;
    L(b, b + 1) = -1.0;
  }
  Eigen::Vector4d w(0.5, 0.7, 1.1, 1.3);
  return L.transpose() * w.asDiagonal() * L;
}

void WriteCrm(const std::string& filename, const Eigen::MatrixXd& pi,
              const std::vector<std::string>& el) {
  std::ofstream f(filename);
  f << "! test charge-response file\n" << pi.rows() << "\n";
  for (Index i = 0; i < pi.rows(); ++i) {
    f << el[std::size_t(i)] << " 0.0 0.0 " << double(i) << "\n";
  }
  f.precision(17);
  for (Index i = 0; i < pi.rows(); ++i) {
    for (Index j = 0; j < pi.cols(); ++j) {
      f << pi(i, j) << " ";
    }
    f << "\n";
  }
}

const std::vector<std::string> methane = {"C", "H", "H", "H", "H"};
}  // namespace

BOOST_AUTO_TEST_SUITE(chargeresponse_test)

BOOST_AUTO_TEST_CASE(load_and_check) {
  WriteCrm("methane.crm", MethanePi(), methane);
  ChargeResponse cr;
  cr.LoadFromFile("methane.crm");
  BOOST_CHECK_EQUAL(cr.size(), 5);
  BOOST_CHECK(cr.Elements() == methane);
  BOOST_CHECK(cr.Pi().isApprox(MethanePi(), 1e-12));
  BOOST_CHECK_SMALL(cr.Pi().rowwise().sum().cwiseAbs().maxCoeff(), 1e-14);
}

BOOST_AUTO_TEST_CASE(rejects_bad_matrices) {
  ChargeResponse cr;
  // not charge conserving
  Eigen::MatrixXd pi = MethanePi();
  pi(0, 0) += 0.1;
  WriteCrm("bad1.crm", pi, methane);
  BOOST_CHECK_THROW(cr.LoadFromFile("bad1.crm"), std::runtime_error);
  // not symmetric
  pi = MethanePi();
  pi(0, 1) += 0.1;
  pi(0, 0) -= 0.1;
  WriteCrm("bad2.crm", pi, methane);
  BOOST_CHECK_THROW(cr.LoadFromFile("bad2.crm"), std::runtime_error);
  // not positive semidefinite
  WriteCrm("bad3.crm", -MethanePi(), methane);
  BOOST_CHECK_THROW(cr.LoadFromFile("bad3.crm"), std::runtime_error);
  // truncated
  std::ofstream f("bad4.crm");
  f << "2\nC 0 0 0\nH 0 0 1\n1 -1\n";
  f.close();
  BOOST_CHECK_THROW(cr.LoadFromFile("bad4.crm"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(permutation) {
  ChargeResponse cr(methane, MethanePi());
  std::vector<Index> order = {2, 0, 4, 1, 3};
  ChargeResponse p = cr.Permuted(order);
  for (Index i = 0; i < 5; ++i) {
    BOOST_CHECK_EQUAL(p.Elements()[std::size_t(i)],
                      methane[std::size_t(order[std::size_t(i)])]);
    for (Index j = 0; j < 5; ++j) {
      BOOST_CHECK_EQUAL(p.Pi()(i, j),
                        cr.Pi()(order[std::size_t(i)], order[std::size_t(j)]));
    }
  }
  BOOST_CHECK_THROW(cr.Permuted({0, 1, 2, 3}), std::runtime_error);
  BOOST_CHECK_THROW(cr.Permuted({0, 1, 2, 3, 3}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(mapping_and_checkpoint) {
  // template files in the atom order C H H H H; the mapping file lists the
  // multipole sites in another order (C, 2, 1, 4, 3), so the mapped segment
  // and its Pi must follow that order
  std::ofstream mps("methane_n.mps");
  mps << "Units angstrom\n";
  const double xyz[5][3] = {{0.1, 0, 0},
                            {0.729118, 0.629118, 0.629118},
                            {-0.529118, -0.629118, 0.629118},
                            {0.729118, -0.629118, -0.629118},
                            {-0.529118, 0.629118, -0.629118}};
  for (Index i = 0; i < 5; ++i) {
    mps << methane[std::size_t(i)] << " " << xyz[i][0] << " " << xyz[i][1]
        << " " << xyz[i][2] << " Rank 0\n"
        << (i == 0 ? "-0.4" : "0.1") << "\nP 1.0\n";
  }
  mps.close();
  WriteCrm("methane_n.crm", MethanePi(), methane);
  std::ofstream map("methane_cr_map.xml");
  map << "<topology><molecules><molecule><name>Methane</name>"
         "<mdname>Methane</mdname><segments><segment><name>Methane</name>"
         "<multipoles_n>methane_n.mps</multipoles_n>"
         "<chargeresponse_n>methane_n.crm</chargeresponse_n>"
         "<map2md>0</map2md><fragments><fragment><name>all</name>"
         "<mdatoms>1:CB:0 1:HB1:1 1:HB2:2 1:HB3:3 1:HB4:4</mdatoms>"
         "<mpoles>0:C 2:H 1:H 4:H 3:H</mpoles>"
         "<localframe>0 1 2</localframe>"
         "</fragment></fragments></segment></segments></molecule>"
         "</molecules></topology>\n";
  map.close();

  Segment seg("Methane", 1);
  const char* names[5] = {"CB", "HB1", "HB2", "HB3", "HB4"};
  for (Index i = 0; i < 5; ++i) {
    Eigen::Vector3d pos(xyz[i][0], xyz[i][1], xyz[i][2]);
    seg.push_back(Atom(1, names[i], i, pos * tools::conv::ang2bohr,
                       methane[std::size_t(i)]));
  }
  Logger log;
  PolarMapper mapper(log);
  mapper.LoadMappingFile("methane_cr_map.xml");
  PolarSegment mol = mapper.map(seg, QMState("n"));
  BOOST_REQUIRE(mol.hasChargeResponse());
  const std::vector<Index> order = {0, 2, 1, 4, 3};
  const Eigen::MatrixXd ref = MethanePi();
  for (Index i = 0; i < 5; ++i) {
    for (Index j = 0; j < 5; ++j) {
      BOOST_CHECK_CLOSE(
          mol.getChargeResponse().Pi()(i, j) + 10.0,
          ref(order[std::size_t(i)], order[std::size_t(j)]) + 10.0, 1e-10);
    }
  }

  // round trip through a checkpoint
  {
    CheckpointFile ff("chargeresponse.hdf5");
    CheckpointWriter ww = ff.getWriter();
    mol.WriteToCpt(ww);
  }
  CheckpointFile ff("chargeresponse.hdf5", CheckpointAccessLevel::READ);
  CheckpointReader rr = ff.getReader();
  PolarSegment mol2(rr);
  BOOST_REQUIRE(mol2.hasChargeResponse());
  BOOST_CHECK(mol2.getChargeResponse().Pi().isApprox(
      mol.getChargeResponse().Pi(), 1e-14));
  BOOST_CHECK(mol2.getChargeResponse().Elements() ==
              mol.getChargeResponse().Elements());

  // a Thole segment stays one through the checkpoint
  PolarSegment thole("thole", 2);
  thole.LoadFromFile("methane_n.mps");
  {
    CheckpointFile f2("thole.hdf5");
    CheckpointWriter w2 = f2.getWriter();
    thole.WriteToCpt(w2);
  }
  CheckpointFile f2("thole.hdf5", CheckpointAccessLevel::READ);
  CheckpointReader r2 = f2.getReader();
  PolarSegment thole2(r2);
  BOOST_CHECK(!thole2.hasChargeResponse());

  // element order must agree
  ChargeResponse cr(methane, MethanePi());
  PolarSegment wrong("wrong", 3);
  wrong.LoadFromFile("methane_n.mps");
  BOOST_CHECK_THROW(wrong.setChargeResponse(cr.Permuted({1, 0, 2, 3, 4})),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
