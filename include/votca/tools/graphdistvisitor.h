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

#ifndef __VOTCA_TOOLS_GRAPH_DIST_VISITOR_H
#define __VOTCA_TOOLS_GRAPH_DIST_VISITOR_H

#include <iostream>
#include <votca/tools/graph_bf_visitor.h>
#include <deque>
#include <queue>
/**
 * \brief A graph visitor that will calculate the distance of each node
 * it is built on top of the graph basic visitor which explores in a breadth
 * first manner. The graph nodes themselves should carry the distance
 * that they are from the starting node on completion of the exploration. 
 *
 */
namespace votca {
namespace tools {

class Graph;
class Edge;
class GraphNode;
class Graph_BF_Visitor;

class GraphDistVisitor : public Graph_BF_Visitor {
  protected:

  public:
    void exploreNode_(std::pair<int,GraphNode&> p_gn, Graph& g, Edge ed = DUMMY_EDGE);    
    GraphDistVisitor(){};
};

}} 
#endif // __VOTCA_TOOLS_GRAPH_DIST_VISITOR_H
