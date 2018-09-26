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

using namespace std;

namespace votca {
    namespace xtp {
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

 void GNode::ReadfromSegment(Segment* seg,int carriertype){
     
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
 
 
 void GNode::AddEventfromQmPair(QMPair* pair,int carriertype){
     double Jeff2=pair->getJeff2(carriertype);
     if(pair->getType()==QMPair::PairType::Excitoncl && carriertype!=2){
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
