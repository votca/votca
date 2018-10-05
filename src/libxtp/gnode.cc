/* 
 * Copyright 2009-2018 The VOTCA Development Team (http://www.votca.org)
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

#include "votca/xtp/gnode.h"
#include <boost/format.hpp>
#include <vector>
#include "votca/xtp/glink.h"
#include <queue>

using namespace std;




namespace votca {
    namespace xtp {

        struct Comparer
{
    bool operator() (hnode * a, hnode * b)
    {
        return (a->prob>b->prob);
    }
};


void GNode::AddDecayEvent(double decayrate)
{
    GLink newEvent;
    newEvent.destination = -1;
    newEvent.rate = decayrate;
    newEvent.initialrate = decayrate;
    newEvent.dr = votca::tools::vec(0.0);
    newEvent.Jeff2 = 0.0;
    newEvent.decayevent=true;
    newEvent.reorg_out = 0.0;
    this->events.push_back(newEvent);
    hasdecay=true;
}

void GNode::AddEvent(int seg2, double rate12, tools::vec dr, double Jeff2, double reorg_out)
{
    GLink newEvent;
    newEvent.destination = seg2;
    newEvent.rate = rate12;
    newEvent.initialrate = rate12;
    newEvent.dr = dr;
    newEvent.Jeff2 = Jeff2;
    newEvent.decayevent=false;
    newEvent.reorg_out = reorg_out;
    this->events.push_back(newEvent);
}


void GNode::InitEscapeRate()
{
    double newEscapeRate = 0.0;
    for(unsigned int i=0; i<this->events.size();i++)
    {
        newEscapeRate += this->events[i].rate;
    }
    this->escape_rate = newEscapeRate;
    // cout << "Escape rate for segment " << this->id << " was set to " << newEscapeRate << endl;
}

GLink* GNode::findHoppingDestination(double p){
    hnode * node=root;
    while (node->leftId!=-1){
        if (p>node->prob) node=&(htree[node->leftId]);
        else node=&(htree[node->rightId]);
    }
    return &events[node->id];
}

void GNode::MakeHuffTree(){
    //queue of the nodes, sorted by probability
    priority_queue<hnode *,vector<hnode *>,Comparer> queue;
    htree=vector<hnode>(events.size()*2-1);
    //first, make a leaf for every GLink
    int index=0;
    for (GLink L:events){
        htree[index].prob=L.rate/escape_rate;
        htree[index].leftId=-1;
        htree[index].rightId=-1;
        htree[index].id=index;
        queue.push(&(htree[index]));
        index++;
    }

     //now connect the hnodes, making a new one for every connection:
     //always take the two nodes with the smallest probability and "combine" them, repeat, until just one node (the root) is left.
     hnode * h1;
     hnode * h2;
     while (queue.size()>1){
         h1=queue.top();
         queue.pop();
         h2=queue.top();
         queue.pop();
         htree[index].prob=h1->prob+h2->prob;
         htree[index].leftId=h1->id;
         htree[index].rightId=h2->id;
         htree[index].id=index;
         queue.push(&(htree[index]));
         index++;
     }
     //save the root
     root=&(htree[index-1]);
     //reorganize the probabilities: in every node, add the probability of one subtree ("small")
     //to all nodes of the other subtree.
     organizeProbabilities(index-1,0);
    moveProbabilities(index-1);
}

void GNode::organizeProbabilities(int id,double add){

    //adds "add" to the probability, then calls itself recursively.
    //this calculates the probabilities needed to traverse the tree quickly
    htree[id].prob+=add;
    //if leftId=-1 (=> node is leaf), returns
    if (htree[id].leftId==-1) return;

    organizeProbabilities(htree[id].leftId,add+htree[htree[id].rightId].prob);
    organizeProbabilities(htree[id].rightId,add);
}

void GNode::moveProbabilities(int id){
    //moves the probabilities "one up" so that htree[id].prob can be checked instead of htree[htree[id].rightId].prob
    if (htree[id].rightId!=-1){
        htree[id].prob=htree[htree[id].rightId].prob;
        moveProbabilities(htree[id].rightId);
        moveProbabilities(htree[id].leftId);
    }
    else htree[id].prob=-1;
}


 void GNode::ReadfromSegment(ctp::Segment* seg,int carriertype){
     
     position=seg->getPos();
     id=seg->getId()-1;
     siteenergy=seg->getSiteEnergy(carriertype);
     
     if (carriertype<2){
         reorg_intorig=seg->getU_nC_nN(carriertype);
         reorg_intdest=seg->getU_cN_cC(carriertype);
     }
     else{
         reorg_intorig=seg->getU_nX_nN(carriertype);
         reorg_intdest=seg->getU_xN_xX(carriertype);
     }
     
    return; 
 }
 
 
 void GNode::AddEventfromQmPair(ctp::QMPair* pair,int carriertype){
     double Jeff2=pair->getJeff2(carriertype);
     if(pair->getType()==ctp::QMPair::PairType::Excitoncl && carriertype!=2){
         return;
     }
     int destination=0;
     double rate12=0.0;
     tools::vec dr=tools::vec(0.0);
     if(id==pair->Seg1()->getId()-1){
        destination=pair->Seg2()->getId()-1;
        rate12=pair->getRate12(carriertype);
        dr=pair->getR();
     }
     else{
         destination=pair->Seg1()->getId()-1;
         rate12=pair->getRate21(carriertype);
         dr=-pair->getR();
     }
    
    double reorg_out=pair->getLambdaO(carriertype);
    AddEvent(destination,rate12,dr,Jeff2,reorg_out);
     
    return; 
 }
        
    }
}
