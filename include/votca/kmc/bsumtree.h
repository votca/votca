/* 
 * Copyright 2009-2011 The VOTCA Development Team (http://www.votca.org)
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

#include <valarray>
//nrelements is number of leaves
//treesize is number of nodes

namespace votca { namespace kmc {
  
using namespace std;

class Bsumtree {
public:
  void initialize(unsigned long nrelements);
  void setrate(unsigned long i, double value);
  double getrate(long i);
  double compute_sum();
  long search(double searchkey);
  void resize(unsigned long newsize);
  long getnrrates();
  
private:
  bool dirty(unsigned long i);
  double partsum(unsigned long i);
  valarray<bool> dirty_array; // Are the subtrees dirty?
  valarray<double> element_array; // The elements (summands)
  valarray<double> partsum_array; // Array of partial sums
  unsigned long treesize;
  unsigned long nrelements;
};

void Bsumtree::initialize(unsigned long nrelements) { // Must be called before use
  // treesize is the smallest power of two above nrelements minus 1
  treesize = (unsigned long) pow(2,ceil(log((double) nrelements)/log((double) 2)))-1; // number of nodes
    
  dirty_array.resize(treesize);
  if (element_array.size()<nrelements) { element_array.resize(nrelements); }
  partsum_array.resize(treesize);
    
  // Initialize arrays
  for (unsigned long i=0;i<treesize;i++) {
    dirty_array[i] = false;
    partsum_array[i] = 0.0;
  }
  for (unsigned long i=0;i<nrelements;i++) {
    element_array[i] = 0.0;
  }
}

void Bsumtree::setrate(unsigned long i, double value) { // 0 <= i < nrelements
  element_array[i] = value;
  long j = i+treesize;
  j = div(j-1,(long) 2).quot; // Parent node
  while (!dirty_array[j]) { // Mark this node and all parents dirty if not already
    dirty_array[j] = true;
    if (j != 0) { // Make sure we stop at the root node
      j = div(j-1,(long) 2).quot; // Parent node
    }
  }
}

double Bsumtree::getrate(long i) {
  return element_array[i];
}
  
double Bsumtree::compute_sum() { // Returns total sum of all elements
  // recursively recompute all dirty nodes
  long i = 0; // Start at root node
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
          i = div(i-1,(long) 2).quot; // Parent node
        }
      }
    }
  }
  return partsum_array[0];
}

// Search returns index to element i: sum(0..i) <= searchkey < sum(0..i+1),
// where the sum is taken over the succesive elements.
long Bsumtree::search(double searchkey) { // Returns index to element
  long maxindex = treesize + nrelements;
  long i = 0; // value must be located in subtree denoted by index i
  while (2*i+2<maxindex) {
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

void Bsumtree::resize(unsigned long newsize) { // Resize arrays. Expensive, so use with care!
  /*
   *  When newsize >= oldsize: all elements are copied, new elements are 0.
   *  When newsize < oldsize: excess elements are thrown away.
   */
  valarray<double> temp_element_array(double(0),newsize); // Temporary storage
  for (unsigned long i=0;i<nrelements && i<newsize;i++) {
    temp_element_array[i] = element_array[i];
  }
  initialize(newsize);
  for (unsigned long i=0;i<newsize;i++) {
    setrate(i, temp_element_array[i]);
  }
}

long Bsumtree::getnrrates() {
  return nrelements;
}

bool Bsumtree::dirty(unsigned long i) {
  if (i<treesize) {
    return dirty_array[i];
  }
  else {
    return false; // Nodes outside the partial rate sum tree are always clean
  }
}

double Bsumtree::partsum(unsigned long i) {
  if (i < treesize) {
    return partsum_array[i];
  }
  else {
    if (i<treesize + nrelements) {
      return element_array[i-treesize];
    }
    else {
      return 0.0; // Non-existent nodes have partial rate sum equal to 0
    }
  }
}
  
}}

#endif
