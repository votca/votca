/*
 *            Copyright 2009-2018 The VOTCA Development Team
 *                       (http://www.votca.org)
 *
 *      Licensed under the Apache License, Version 2.0 (the "License")
 *
 * You may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *              http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#ifndef __VOTCA_TOOLS_GRAPH_VISITOR_H
#define __VOTCA_TOOLS_GRAPH_VISITOR_H

#include <vector>
#include <set>
#include <votca/tools/edge.h>
#include <votca/tools/graphnode.h>

/**
 * \brief a graph visitor outline for creating graph visitor objects
 */
namespace votca {
namespace tools {

class Graph;

class GraphVisitor{
  protected:
    std::set<int> explored_;
    int startingVertex_;
    // Determine which vertices have been unexplored 
    std::vector<int> getUnexploredVertex_(Edge ed);
    // What is done to an individual graph node as it is explored
    virtual void addEdges_(Graph& g, int vertex);
    virtual Edge getEdge_(Graph g);
    // Edge(0,0) is a dummy value
    virtual void exploreNode_(std::pair<int,GraphNode&> p_gn,Graph g,Edge ed = DUMMY_EDGE);
  public:

    GraphVisitor() : startingVertex_(0) {};
    
    // Determine if the exploration is complete
    virtual bool queEmpty();
    // Which node the exploration begins at. 
    void startingVertex(Graph& g, int vertex=0);
    // What the visitor does to each node as it is visited
    void exec(Graph& g, Edge ed);    
    // The next node to be explored
    Edge nextEdge(Graph g); 
    // Get the set of all the vertices that have been explored
    std::set<int> getExploredVertices();
};

}}

#endif // __VOTCA_TOOLS_GRAPH_VISITOR_H
