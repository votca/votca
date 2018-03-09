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

/**
 * \brief a graph visitor outline for creating graph visitor objects
 */
namespace votca {
namespace tools {

class Grapth;
class Edge;

class GraphVisitor{
  private:
    // First element is the distance for breadth first search
    // Second element is a queue to pop the edges out as we go along
    std::map<int,queue<Edge>> edge_que_;
    // Keep track of whether the vertex has been explored or not
    std::vector<bool> explored_;
 
  public:

    // What the visitor does to each node as it is visited
    virtual void exec(Graph g, Edge ed){};    
    // The next node to be explored
    virtual Edge nextEdge(Graph g) {}; 
};

}}

#endif // __VOTCA_TOOLS_GRAPH_VISITOR_H
