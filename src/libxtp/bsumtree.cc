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

#include "votca/xtp/bsumtree.h"
#include <cmath>
#include <cstdlib>

namespace votca {
    namespace xtp {
        
void Bsumtree::initialize(int nr_elements) { // Must be called before use

    // treesize is the smallest power of two above nrelements minus 1
    nrelements = nr_elements;
    treesize = std::pow(2,std::ceil(log((double) nrelements)/std::log((double) 2)))-1; // number of nodes
  
    // clear the arrays
    dirty_array.clear();
    element_array.clear();
    partsum_array.clear();

    // Initialize arrays
    for (int i=0;i<treesize;i++) {
        dirty_array.push_back(false);
        partsum_array.push_back(0.0);
    }
    for (int i=0;i<nrelements;i++) {
        element_array.push_back(0.0);
    }
}

void Bsumtree::setrate(int i, double value) { // 0 <= i < nrelements
    element_array[i] = value;
    int j = i+treesize;
    j = std::div(j-1, 2).quot; // Parent node
    while (!dirty_array[j]) { // Mark this node and all parents dirty if not already
        dirty_array[j] = true;
        if (j != 0) { // Make sure we stop at the root node
            j = std::div(j-1, 2).quot; // Parent node
        }
    }
}

double Bsumtree::getrate(int i) {
    return element_array[i];
}
  
double Bsumtree::compute_sum() { // Returns total sum of all elements
    // recursively recompute all dirty nodes
    int i = 0; // Start at root node
    while (dirty(i)) {
        if (dirty(2*i + 1)) { // Is left subtree dirty ?
            i = 2*i + 1;
        }
        else {
            if (dirty(2*i + 2)) { // Is right subtree dirty?
                i = 2*i + 2;
            }
            else { // Both subtrees are clean, update this node
                partsum_array[i] = partsum(2*i+1) + partsum(2*i+2);
                dirty_array[i] = false;
                if (i != 0) { // Make sure we stop at the root node
                    i = std::div(i-1, 2).quot; // Parent node
                }
            }
        }
    }
    return partsum_array[0];
}

// Search returns index to element i: sum(0..i) <= searchkey < sum(0..i+1),
// where the sum is taken over the succesive elements.
long Bsumtree::search(double searchkey) { // Returns index to element
    int maxindex = treesize + nrelements;

    int i = 0; // value must be located in subtree denoted by index i
    while (2*i+2<=maxindex) {
        if (searchkey <= partsum(2*i+1)) { // value is located in left subtree
            i = 2*i+1;
        }
        else { // value is located in right subtree
            searchkey -= partsum(2*i+1); // values are relative
            i = 2*i+2;
        }
    }    

    i -= treesize;

    return i;
}


void Bsumtree::resize(int newsize) { // Resize arrays. Expensive, so use with care!
    /*
     *  When newsize >= oldsize: all elements are copied, new elements are 0.
     *  When newsize < oldsize: excess elements are thrown away.
     */
    std::vector<double> temp_element_array; // Temporary storage
    temp_element_array.clear();

    for (int i=0;i<nrelements && i<newsize;i++) {
        temp_element_array.push_back(element_array[i]);
    }
  
    int oldsize = nrelements;
    initialize(newsize);
  
    for (int i=0;i<oldsize;i++) {
        setrate(i, temp_element_array[i]);
    }
}

long Bsumtree::getnrrates() {
    return nrelements;
}

bool Bsumtree::dirty(int i) {
    if (i<treesize) {
        return dirty_array[i];
    }
    else {
        return false; // Nodes outside the partial rate sum tree are always clean
    }
}

double Bsumtree::partsum(int i) {
    if (i < treesize) {
        return partsum_array[i];
    }
    else {
        if (i<treesize + nrelements) {
            return element_array[i-treesize];
        }
        else {
            return 0.0; // Non-existing nodes have partial rate sum equal to 0
        }
    }
}


        
    }
}
