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

#ifndef _VOTCA_KMC_LINK_H
#define	_VOTCA_KMC_LINK_H

//#include <votca/xtp/node.h>
 using namespace std;
namespace votca { namespace xtp {
   
class Node;

/**
 * \brief A link between two nodes
 * 
 * Container for pair properties: rates, couplings, separations 
 * 
 */
class Link
{

public:
    Link( long id, Node* node1, Node* node2, votca::tools::vec r12) {
        _id  = id;
        _node1 = node1;
        _node2 = node2;
        _r12 = r12;
    } 
    virtual ~Link(){};
    
    /// link ID
    const long &id() const { return _id; }
    
    const long &reverse_id() const {return _reverse_id; }
    
    /// r2 - r1
    const votca::tools::vec &r12() const { return _r12; }    
    
    /// node 1
    Node* &node1() { return _node1; }
    
    /// node 2
    Node* &node2() { return _node2; }
    
    /// set initial and final nodes of link
    void SetNodes(Node* node1, Node* node2) { _node1 = node1; _node2 = node2; }
    
    /// (re)set ID
    void SetID(long id) { _id = id; } 
    void SetReverseID(long reverse_id) { _reverse_id = reverse_id; }
    
    // print Link info
    virtual void Print(std::ostream &out) {
        out << _id << " " << _r12 << endl ;
    }
    
    void setcount(int type) {_count = type;}
    void inccount() {_count++;}
    void incval(int val) {_count += val;}
    void deccount() {_count--;}
    void decval(int val) {_count -= val;}
    const int &count() const { return _count; }
    
private:
    
    long _id;
    long _reverse_id;
    
    Node *_node1;
    Node *_node2;
    
    int _count;
    
    /// r2 - r1
    votca::tools::vec _r12;       
};

}} 

#endif // _VOTCA_KMC_LINK_H
