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
#include "xtp_libint2.h"
#define BOOST_TEST_MAIN

#define BOOST_TEST_MODULE activeregion_test

// Standard includes
#include <fstream>

// Third party includes
#include <boost/test/unit_test.hpp>

// Local VOTCA includes
#include "votca/xtp/activeregion.h"
#include "votca/xtp/orbitals.h"

using namespace votca::xtp;
using votca::Index;

BOOST_AUTO_TEST_SUITE(activeregion_test)

BOOST_AUTO_TEST_CASE(connectivity_per_molecule) {
  // two waters 5 A apart
  std::ofstream xyz("dimer.xyz");
  xyz << "6\n\n"
      << "O 0.000 0.000 0.118\nH 0.000 0.755 -0.471\nH 0.000 -0.755 -0.471\n"
      << "O 5.000 0.000 0.118\nH 5.000 0.755 -0.471\nH 5.000 -0.755 -0.471\n";
  xyz.close();
  QMMolecule mol("dimer", 0);
  mol.LoadFromFile("dimer.xyz");
  const auto graph = ActiveRegion::BondGraph(mol);
  BOOST_CHECK_EQUAL(graph[0].size(), 2);
  BOOST_CHECK_EQUAL(graph[1].size(), 1);
  // whole first water plus part of the second: one piece per molecule
  BOOST_CHECK_NO_THROW(ActiveRegion::CheckConnected(mol, {0, 1, 2, 3, 4}));
  BOOST_CHECK_NO_THROW(ActiveRegion::CheckConnected(mol, {0, 1, 2, 4}));
  // the two hydrogens of one water without their oxygen: two pieces
  BOOST_CHECK_THROW(ActiveRegion::CheckConnected(mol, {0, 4, 5}),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(grow_from_frontier_orbitals) {
  libint2::initialize();
  Orbitals orb;
  orb.ReadFromCpt(std::string(XTP_TEST_DATA_FOLDER) +
                  "/activedftengine/molecule.orb");
  ActiveRegion::Options opt;

  // one hydrogen of water: the HOMO (oxygen lone pair) is not on it
  const ActiveRegion::Report before = ActiveRegion::Evaluate(orb, {1}, opt);
  BOOST_CHECK(before.added.empty());
  BOOST_REQUIRE_EQUAL(before.orbitals.size(), 1);
  BOOST_CHECK_EQUAL(before.orbitals[0], orb.getNumberOfAlphaElectrons() - 1);
  BOOST_CHECK_LT(before.coverage_input(0), 0.1);

  // growing adds the oxygen and with it the second hydrogen
  opt.grow = true;
  const ActiveRegion::Report grown = ActiveRegion::Evaluate(orb, {1}, opt);
  BOOST_CHECK_EQUAL(grown.atoms.size(), 3);
  BOOST_CHECK_EQUAL(grown.added.size(), 2);
  BOOST_CHECK_CLOSE(grown.coverage_final(0), 1.0, 1e-6);
  BOOST_TEST_MESSAGE(grown.Format());

  // the oxygen alone between the two active hydrogens would be cut out
  // twice: refused as disconnected input, reported as enclosed otherwise
  opt.grow = false;
  const ActiveRegion::Report one = ActiveRegion::Evaluate(orb, {0, 1}, opt);
  BOOST_CHECK(one.enclosed.empty());
  opt.grow = true;

  // the oxygen without its hydrogens: two X-H bonds cut, reported; grow
  // adds the hydrogens
  opt.grow = false;
  const ActiveRegion::Report oxygen = ActiveRegion::Evaluate(orb, {0}, opt);
  BOOST_CHECK_EQUAL(oxygen.cut_hydrogens.size(), 2);
  opt.grow = true;

  // a disconnected input is refused
  BOOST_CHECK_THROW(ActiveRegion::Evaluate(orb, {1, 2}, opt),
                    std::runtime_error);
  libint2::finalize();
}

BOOST_AUTO_TEST_SUITE_END()
