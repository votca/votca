/* 
 * Copyright 2009-2016 The VOTCA Development Team (http://www.votca.org)
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

#ifndef __VOTCA_KMC_BSUMTREE_H_
#define __VOTCA_KMC_BSUMTREE_H_
#include <vector>


//nrelements is number of leaves
//treesize is number of nodes

namespace votca { namespace xtp {
  


class Bsumtree {
    
public:
 
    /// initialize binary sumtree objects
    void initialize(int nrelements);
    
    /// set rate of element i
    void setrate(int i, double value);
    
    /// obtain rate of element i
    double getrate(int i);
    
    /// calculate all partial sums and obtain total sum of all elements in the tree
    double compute_sum();
    
    /// search element with searchkey
    long search(double searchkey);
    
    /// resize binary sum trees
    void resize(int newsize);
    
    /// get nr of elements
    long getnrrates();
 
    double partsum(int i);
    
private:
    
    bool dirty(int i);

    std::vector<bool> dirty_array; // Are the subtrees dirty?
    std::vector<double> element_array; // The elements (summands)
    std::vector<double> partsum_array; // Array of partial sums
    int treesize;
    int nrelements;
    
};


  
}}

#endif
