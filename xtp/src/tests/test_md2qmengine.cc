/*
 * Copyright 2009-2026 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE md2qmengine_test

// Standard includes
#include <cstdio>
#include <cstdlib>
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// VOTCA includes
#include <votca/csg/interaction.h>
#include <votca/csg/topology.h>

// Local VOTCA includes
#include "votca/xtp/md2qmengine.h"

using namespace votca::xtp;
using namespace votca;
BOOST_AUTO_TEST_SUITE(md2qmengine_test)

BOOST_AUTO_TEST_CASE(external_bond_direction_detection_test) {
  // Directly exercises Md2QmEngine::map()'s own computation of the
  // RAW external-bond direction from real, MD-level bond connectivity
  // -- the one piece of the whole external-bond-direction pipeline
  // (Atom's own storage, this computation, QMAtom's own parallel
  // storage, SegmentMapper's own transform-and-transfer) that had no
  // direct test coverage at all until now -- test_segmentmapper.cc's
  // own external_bond_direction_transfer_test only ever exercises the
  // transform-and-transfer half, given an already-set input
  // direction, never this detection half at all.
  //
  // Per the user's own, direct suggestion: based on the bithiophene
  // geometry already used earlier this session, with each ring
  // mapped to its own, SEPARATE segment (matching a real, plausible
  // "polymer terminal unit" scenario) -- confirmed directly, this
  // session, that this specific, segment-level (not fragment-level)
  // granularity is exactly what this feature is designed to detect:
  // a real bond crossing an actual *segment* boundary, not merely an
  // internal *fragment* boundary within one, single segment (the
  // "3 fragments, 1 segment" case, matching test_segmentmapper.cc's
  // own ch4.xml setup, never triggers this feature at all, confirmed
  // directly, by inspecting Md2QmEngine::map()'s own code, since it
  // never tracks fragment membership when building its own output at
  // all, only segment membership).
  //
  // Single, real, synthetic csg::Topology: one molecule ("Bithiophene"),
  // 16 beads, matching this session's own, already-established
  // bithiophene.xyz atom order and geometry exactly (positions
  // converted from the original Angstrom to nm here, since
  // Md2QmEngine::map() itself directly converts bead positions from
  // nm to bohr via tools::conv::nm2bohr) -- with real bonded
  // interactions (csg::IBond) matching that same, already-confirmed
  // connectivity, including the key, single C-C bond connecting the
  // two rings (atom 1 <-> atom 9). The mapping file splits this one
  // molecule into two, separate segments, "ThiopheneLeft" (atoms
  // 0-7) and "ThiopheneRight" (atoms 8-15) -- so the 1<->9 bond is
  // exactly the one, real bond that should be detected as crossing a
  // segment boundary; every other bond in the structure stays within
  // one segment or the other.
  csg::Topology top;
  csg::Molecule* mol = top.CreateMolecule("Bithiophene");

  struct BeadSpec {
    std::string element;
    Eigen::Vector3d pos_angstrom;
  };
  // Same order/positions as this session's own bithiophene.xyz.
  std::vector<BeadSpec> beads = {
      {"S", {0.000000, 0.000000, 0.000000}},
      {"C", {1.714000, 0.000000, 0.000000}},
      {"C", {2.216107, 1.274672, 0.000000}},
      {"C", {1.192486, 2.263171, 0.000000}},
      {"C", {-0.063886, 1.716885, 0.000000}},
      {"H", {3.274448, 1.504460, 0.000000}},
      {"H", {1.385199, 3.328887, 0.000000}},
      {"H", {-0.980470, 2.293742, 0.000000}},
      {"S", {4.249695, -1.206821, 0.000000}},
      {"C", {2.535695, -1.206821, 0.000000}},
      {"C", {2.033589, -2.481493, 0.000000}},
      {"C", {3.057209, -3.469992, 0.000000}},
      {"C", {4.313581, -2.923706, 0.000000}},
      {"H", {0.975247, -2.711282, 0.000000}},
      {"H", {2.864496, -4.535708, 0.000000}},
      {"H", {5.230165, -3.500563, 0.000000}},
  };

  for (Index i = 0; i < Index(beads.size()); i++) {
    const BeadSpec& b = beads[size_t(i)];
    csg::Bead* bead =
        top.CreateBead(csg::Bead::spherical, b.element, b.element, 1, 1.0, 0.0);
    // Angstrom -> nm (Md2QmEngine::map() itself converts nm -> bohr).
    bead->setPos(b.pos_angstrom * 0.1);
    mol->AddBead(bead, b.element);
  }

  std::vector<std::pair<Index, Index>> bonds = {
      // Ring A (atoms 0-7)
      {0, 1},
      {1, 2},
      {2, 3},
      {3, 4},
      {4, 0},
      {2, 5},
      {3, 6},
      {4, 7},
      // Ring B (atoms 8-15)
      {8, 9},
      {9, 10},
      {10, 11},
      {11, 12},
      {12, 8},
      {10, 13},
      {11, 14},
      {12, 15},
      // The one, real, inter-ring bond -- the only one crossing a
      // segment boundary in this test's own mapping below.
      {1, 9},
  };
  for (const auto& b : bonds) {
    // Real, direct fix for a real, genuine, repeated instance of the
    // exact same setGroup() bug already fixed twice earlier this same
    // session, in gmxtopologyreader.cc and xtp_map.cc's own
    // GuessBonds -- confirmed directly, this third time, from the
    // user's own real, direct CI failure report. A freshly-
    // constructed IBond's own group_ starts out empty by default;
    // getGroup() itself directly asserts this is non-empty.
    csg::Interaction* ic = new csg::IBond(b.first, b.second);
    ic->setGroup("BONDS");
    top.AddBondedInteraction(ic);
  }

  // Mapping file: splits the single "Bithiophene" molecule into two,
  // separate segments, one per ring -- deliberately omits
  // qmcoords/map2md entirely, confirmed directly, by reading
  // Md2QmEngine::map()/CheckMappingFile()'s own code first, that
  // neither is ever accessed there at all (only SegmentMapper::map(),
  // not exercised by this test, uses qmcoords).
  std::string mapfile = "md2qmengine_bithiophene_test.xml";
  std::ofstream f(mapfile);
  f << R"(<topology>
  <molecules>
    <molecule>
      <name>Bithiophene</name>
      <mdname>Bithiophene</mdname>
      <segments>
        <segment>
          <name>ThiopheneLeft</name>
          <fragments>
            <fragment>
              <name>ring_a</name>
              <mdatoms>1:S:0 1:C:1 1:C:2 1:C:3 1:C:4 1:H:5 1:H:6 1:H:7</mdatoms>
              <qmatoms>0:S 1:C 2:C 3:C 4:C 5:H 6:H 7:H</qmatoms>
              <weights>32 12 12 12 12 1 1 1</weights>
              <localframe>0 1 4</localframe>
            </fragment>
          </fragments>
        </segment>
        <segment>
          <name>ThiopheneRight</name>
          <fragments>
            <fragment>
              <name>ring_b</name>
              <mdatoms>1:S:8 1:C:9 1:C:10 1:C:11 1:C:12 1:H:13 1:H:14 1:H:15</mdatoms>
              <qmatoms>0:S 1:C 2:C 3:C 4:C 5:H 6:H 7:H</qmatoms>
              <weights>32 12 12 12 12 1 1 1</weights>
              <localframe>0 1 4</localframe>
            </fragment>
          </fragments>
        </segment>
      </segments>
    </molecule>
  </molecules>
</topology>
)";
  f.close();

  Md2QmEngine engine(mapfile);
  Topology xtptop = engine.map(top);

  BOOST_CHECK_EQUAL(xtptop.Segments().size(), 2);

  const Segment* left = nullptr;
  const Segment* right = nullptr;
  for (const Segment& seg : xtptop.Segments()) {
    if (seg.getType() == "ThiopheneLeft") {
      left = &seg;
    } else if (seg.getType() == "ThiopheneRight") {
      right = &seg;
    }
  }
  BOOST_REQUIRE(left != nullptr);
  BOOST_REQUIRE(right != nullptr);
  BOOST_CHECK_EQUAL(left->size(), 8);
  BOOST_CHECK_EQUAL(right->size(), 8);

  // Atom 1 (ring-junction C, on "ThiopheneLeft") must have its own
  // external bond correctly detected -- direction (raw, MD-level,
  // untransformed -- this test does not go through
  // SegmentMapper::map() at all) pointing toward atom 9's own real
  // position.
  const Atom* left_junction = nullptr;
  for (const Atom& a : *left) {
    if (a.getId() == 1) {
      left_junction = &a;
    }
  }
  BOOST_REQUIRE(left_junction != nullptr);
  BOOST_CHECK_EQUAL(left_junction->hasExternalBond(), true);

  Eigen::Vector3d expected_left_dir =
      (beads[9].pos_angstrom - beads[1].pos_angstrom).normalized();
  Eigen::Vector3d actual_left_dir = left_junction->getExternalBondDirection();
  BOOST_CHECK_EQUAL(actual_left_dir.isApprox(expected_left_dir, 1e-5), true);
  if (!actual_left_dir.isApprox(expected_left_dir, 1e-5)) {
    std::cout << "left junction direction" << std::endl;
    std::cout << actual_left_dir << std::endl;
    std::cout << "expected" << std::endl;
    std::cout << expected_left_dir << std::endl;
  }

  // Not just THAT an external bond exists, but WHICH segment it
  // crosses into -- must be the real "ThiopheneRight" segment's own,
  // actual, real id (right->getId()), not merely non-negative/set.
  BOOST_CHECK_EQUAL(left_junction->getExternalBondPartnerSegmentId(),
                    right->getId());

  // Symmetric check on atom 9 (ring-junction C, on "ThiopheneRight"),
  // pointing back toward atom 1.
  const Atom* right_junction = nullptr;
  for (const Atom& a : *right) {
    if (a.getId() == 9) {
      right_junction = &a;
    }
  }
  BOOST_REQUIRE(right_junction != nullptr);
  BOOST_CHECK_EQUAL(right_junction->hasExternalBond(), true);

  Eigen::Vector3d expected_right_dir =
      (beads[1].pos_angstrom - beads[9].pos_angstrom).normalized();
  Eigen::Vector3d actual_right_dir = right_junction->getExternalBondDirection();
  BOOST_CHECK_EQUAL(actual_right_dir.isApprox(expected_right_dir, 1e-5), true);
  if (!actual_right_dir.isApprox(expected_right_dir, 1e-5)) {
    std::cout << "right junction direction" << std::endl;
    std::cout << actual_right_dir << std::endl;
    std::cout << "expected" << std::endl;
    std::cout << expected_right_dir << std::endl;
  }

  // Symmetric check: must point back to the real "ThiopheneLeft"
  // segment's own, actual id.
  BOOST_CHECK_EQUAL(right_junction->getExternalBondPartnerSegmentId(),
                    left->getId());

  // Every other atom, in either segment, has no real external bond at
  // all (every other bond in the structure stays within its own,
  // single segment) -- must not report one.
  for (const Segment* seg : {left, right}) {
    for (const Atom& a : *seg) {
      if (a.getId() == 1 || a.getId() == 9) {
        continue;
      }
      BOOST_CHECK_EQUAL(a.hasExternalBond(), false);
    }
  }
}
BOOST_AUTO_TEST_SUITE_END()
