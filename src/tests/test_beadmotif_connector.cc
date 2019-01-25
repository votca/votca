/*
 * Copyright 2009-2019 The VOTCA Development Team (http://www.votca.org)
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

#define BOOST_TEST_MODULE beadmotifconnector_test
#include <boost/test/unit_test.hpp>

#include <votca/csg/beadmotifconnector.h>

using namespace std;
using namespace votca::csg;
using namespace votca::tools;

// used for rounding doubles so we can compare them
double round_(double v, int p) {
  v *= pow(10, p);
  v = round(v);
  v /= pow(10, p);
  return v;
}

BOOST_AUTO_TEST_SUITE(beadmotifconnector_test)

BOOST_AUTO_TEST_CASE(test_beadmotifconnector_constructor) {
  BeadMotifConnector bead_connector;
}

BOOST_AUTO_TEST_CASE(test_beadmotifconnector_addmotifandbeadedge) {
  BeadMotifConnector bead_connector;
  // Connects two beads
  Edge ed_bead(1, 2);
  // Connects motifs
  Edge ed_motif(4, 5);
  bead_connector.AddMotifAndBeadEdge(ed_motif, ed_bead);
}

BOOST_AUTO_TEST_CASE(test_beadmotifconnector_getedges) {
  BeadMotifConnector bead_connector;
  // Connects two beads
  Edge ed_bead(1, 2);
  // Connects motifs
  Edge ed_motif(4, 5);
  bead_connector.AddMotifAndBeadEdge(ed_motif, ed_bead);

  vector<Edge> bead_edges = bead_connector.getBeadEdges(ed_motif);
  BOOST_CHECK_EQUAL(bead_edges.size(), 1);
  BOOST_CHECK_EQUAL(bead_edges.at(0), ed_bead);

  Edge motif_edge = bead_connector.getMotifEdge(ed_bead);
  BOOST_CHECK_EQUAL(motif_edge, ed_motif);

  Edge ed_bead2(7, 8);
  bead_connector.AddMotifAndBeadEdge(ed_motif, ed_bead2);

  bead_edges = bead_connector.getBeadEdges(ed_motif);
  BOOST_CHECK_EQUAL(bead_edges.size(), 2);

  Edge ed_motif2(4, 7);
  Edge ed_bead3(9, 11);
  bead_connector.AddMotifAndBeadEdge(ed_motif2, ed_bead3);

  unordered_set<Edge> motif_edges = bead_connector.getMotifEdges();
  BOOST_CHECK_EQUAL(motif_edges.size(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
