/*
 *            Copyright 2009-2020 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writingraph, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

//#ifndef VOTCA_TOOLS_BRANCHCANONIZER_H
//#define VOTCA_TOOLS_BRANCHCANONIZER_H
//#pragma once
//
//#include <string>
//#include <vector>
//
//#include "votca/tools/graph.h"
//
// namespace votca {
//  namespace tools {
//
//
//    /**
//     * @brief The purpose of this class is to organize the branches, and nodes
//     * of a graph into unique sequence
//     */
//    class BranchCanonizer {
//
//      private:
//
//        Graph starting_graph_;
//        /// The nodes in the chosen graph will be different than from the
//        starting
//        /// graph
//        Graph chosen_graph_;
//        BranchTree chosen_tree_;
//        std::vector<Index> canonized_vertex_sequence_;
//      public:
//
//        BranchCanonizer(Graph graph) : graph_(graph) {};
//        void canonize();
//
//        // Get full string id of the canonized graph, should be possible from
//        // the string id to build the graph
//        std::string getStringId();
//        // Get a hash of the unique string id, should make comparisons faster
//        // the hash should be a much smaller value with the caveat that the
//        // structure cannot be determined from it alone
//        std::string getHash();
//        // Get the canonized sequence of the vertices, they will appear in an
//        // order specific to the structure, the sequence should be invariante
//        to
//        // symmetry that exists in the graph
//        std::vector<Index> getSequence();
//    };
//  }
//}
//#endif // VOTCA_TOOLS_BRANCHCANONIZER_H
