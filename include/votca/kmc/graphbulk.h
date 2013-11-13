/*
 * Copyright 2009-2013 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_GRAPHBULK_H_
#define __VOTCA_KMC_GRAPHBULK_H_

#include <vector>
#include <votca/kmc/graphsql.h>
//#include <votca/kmc/graphcubic.h>

namespace votca { namespace kmc {

template<class TGraph, class TNode, class TLink>    
class GraphBulk : public TGraph {

public:
    
    ///add backwards versions of the hopping links to the links vector
    void LinkBackwards(){};
    ///sort all links per node, ordered by node ID (automatic check for duplicates) and declare the link information to all seperate nodes
    void LinkSort(){};
    ///calculate the maximum of all degrees in the graph
    int MaxPairDegree(){};
    ///prepare links vector for binary search tree algorithm by adding per node zero rate links to this node till the node-based link vectors have a size equal to max_pair_degree
    void LinkBinary(){};
    /*typename std::vector<TLink*>::iterator it;
    for (it = graph->_links.begin(); it != graph->_links.end(); it++ ) {;}*/

    ///sort links in links vector and split forwards and backwards links in seperate links
   // void LinkSort();
    
    //int test;
    
};


}}



#endif

