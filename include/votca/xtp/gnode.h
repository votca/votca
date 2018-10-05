/*
 * Copyright 2009-2017 The VOTCA Development Team (http://www.votca.org)
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

#ifndef _VOTCA_KMC_GNODE_H
#define _VOTCA_KMC_GNODE_H

#include <votca/tools/vec.h>
#include <votca/xtp/glink.h>
#include <votca/ctp/segment.h>
#include <votca/ctp/qmpair.h>
#include <vector>


using namespace std;


namespace votca { namespace xtp {


//hnode as huffman-node
struct hnode{
    //id's of the child nodes
    //if node is a leaf, both are null
    //hnode * small;
    //hnode * big;

    //own index/id
    int id;
    //id's (indexes) of the children
    int leftId;
    int rightId;
    //probability of this node being chosen
    double prob;
    //pointer to the gnode, if hnode is a leaf
    //GLink * edge; not needed, as the index in the edge list is the same as the hnodes index
    //in the htree vector
};



class GNode
{
    private:
        void organizeProbabilities(int id, double add);
        void moveProbabilities(int id);

    public:
        GNode():occupied(false),occupationtime(0.0),escape_rate(0.0),hasdecay(false){};
        
        ~GNode(){};

        int id;
        bool occupied;
        bool injectable;
        double occupationtime;
        double escape_rate;
        bool hasdecay;
        tools::vec position;
        std::vector<GLink> events;
        // stuff for Coulomb interaction:
        double siteenergy;
        double reorg_intorig; // UnCnN
        double reorg_intdest; // UcNcC
        void AddEvent(int seg2, double rate12, tools::vec dr, double Jeff2, double reorg_out);
        const double &getEscapeRate(){return escape_rate;}
        void InitEscapeRate();
        void AddDecayEvent(double decayrate);
        void ReadfromSegment(ctp::Segment* seg, int carriertype);
        void AddEventfromQmPair(ctp::QMPair* pair,int carriertype);
        
        GLink* findHoppingDestination(double p);
        hnode * root;
        void MakeHuffTree();
        std::vector <hnode> htree;
};





}}

#endif  /* _VOTCA_KMC_GNODE_H */

