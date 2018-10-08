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

#include <votca/xtp/kmccalculator.h>
#include <votca/xtp/gnode.h>
#include <votca/tools/constants.h>
#include <boost/format.hpp>
#include <votca/ctp/topology.h>
#include <locale>

using namespace std;

namespace votca {
    namespace xtp {
        
        KMCCalculator::KMCCalculator(){};

    void KMCCalculator::LoadGraph(ctp::Topology *top) {

        std::vector< ctp::Segment* >& seg = top->Segments();
        
        if(seg.size()<1){
          throw std::runtime_error("Your sql file contains no segments!");
        }
        
        for (unsigned i = 0; i < seg.size(); i++) {
            GNode *newNode = new GNode();
            newNode->ReadfromSegment(seg[i], _carriertype);
            if (tools::wildcmp(_injection_name.c_str(), seg[i]->getName().c_str())) {
                newNode->injectable = true;
            } else {
                newNode->injectable = false;
            }
            _nodes.push_back(newNode);
        }

        ctp::QMNBList &nblist = top->NBList();
        if(nblist.size()<1){
          throw std::runtime_error("Your sql file contains no pairs!");    
        }
        
        
        for (ctp::QMNBList::iterator it = nblist.begin(); it < nblist.end(); ++it) {
            _nodes[(*it)->Seg1()->getId()-1]->AddEventfromQmPair(*it, _carriertype);
            _nodes[(*it)->Seg2()->getId()-1]->AddEventfromQmPair(*it, _carriertype);
        }
        
        unsigned events=0;
        unsigned max=std::numeric_limits<unsigned>::min();
        unsigned min=std::numeric_limits<unsigned>::max();
        minlength=std::numeric_limits<double>::max();
        double maxlength=0;
        for(const auto& node:_nodes){
            
            unsigned size=node->events.size();
            for( const auto& event:node->events){
                if(event.decayevent){continue;}
                double dist=abs(event.dr);
                if(dist>maxlength){
                    maxlength=dist;
                } else if(dist<minlength){
                    minlength=dist;        
                }
            }
            
            events+=size;
            if(size==0){
                cout<<"Node "<<node->id<<" has 0 jumps"<<endl;
            }
            else if(size<min){
                min=size;
            }
            else if(size>max){
                max=size;
            }
        }
        double avg=double(events)/double(_nodes.size());
        double deviation=0.0;
        for(const auto& node:_nodes){
            double size=node->events.size();
            deviation+=(size-avg)*(size-avg);
        }
        deviation=std::sqrt(deviation/double(_nodes.size()));
        
        cout<<"Nblist has "<<nblist.size()<<" pairs. Nodes contain "<<events<<" jump events"<<endl;
        cout<<"with avg="<<avg<<" std="<<deviation<<" max="<<max<<" min="<<min<<endl;
        cout<<"Minimum jumpdistance ="<<minlength<<" nm Maximum distance ="<<maxlength<<" nm"<<endl;
        cout<<"Grouping into "<<lengthdistribution<<" boxes"<<endl;
        lengthresolution=(1.00001*maxlength-minlength)/double(lengthdistribution);
        cout<<"Resolution is "<<lengthresolution<<" nm"<<endl;   
       
        _jumplengthdistro=std::vector<long unsigned>(lengthdistribution,0);
        _jumplengthdistro_weighted=std::vector<double>(lengthdistribution,0);

       
        cout << "spatial density: " << _numberofcharges / top->BoxVolume() << " nm^-3" << endl;

        for (unsigned int i = 0; i < _nodes.size(); i++) {
            _nodes[i]->InitEscapeRate();
        }
            
        return;
    }
    

        void KMCCalculator::ResetForbiddenlist(std::vector<int> &forbiddenid) {
            forbiddenid.clear();
            return;
        }

        void KMCCalculator::AddtoForbiddenlist(int id, std::vector<int> &forbiddenid) {
            forbiddenid.push_back(id);
            return;
        }

        bool KMCCalculator::CheckForbidden(int id,const std::vector<int> &forbiddenlist) {
            bool forbidden = false;
            for (unsigned int i = 0; i < forbiddenlist.size(); i++) {
                if (id == forbiddenlist[i]) {
                    forbidden = true;
                    break;
                }
            }
            return forbidden;
        }

        bool KMCCalculator::CheckSurrounded(GNode* node,const std::vector<int> & forbiddendests) {
            bool surrounded = true;
            for (unsigned  i = 0; i < node->events.size(); i++) {
                bool thisevent_possible = true;
                for (unsigned int j = 0; j < forbiddendests.size(); j++) {
                    if (node->events[i].destination == forbiddendests[j]) {
                        thisevent_possible = false;
                        break;
                    }
                }
                if (thisevent_possible == true) {
                    surrounded = false;
                    break;
                }
            }
            return surrounded;
        }

        
        
        
        std::string KMCCalculator::CarrierInttoLongString(int carriertype){
            std::string name="";
            if (carriertype==-1){
                name="electron";
            }
            else if(carriertype==1){
                name="hole";
            }
            else if(carriertype==2){
                name="singlet";
            }
            else if(carriertype==3){
                name="triplet";
            }
            else{
                throw runtime_error((boost::format("Carriertype %i not known") % carriertype).str());
            }
            return name;
        }
        
        std::string KMCCalculator::CarrierInttoShortString(int carriertype){
            std::string name="";
            if (carriertype==-1){
                name="e";
            }
            else if(carriertype==1){
                name="h";
            }
            else if(carriertype==2){
                name="s";
            }
            else if(carriertype==3){
                name="t";
            }
            else{
                throw runtime_error((boost::format("Carriertype %i not known") % carriertype).str());
            }
            return name;
        }
        
         int KMCCalculator::StringtoCarriertype(std::string name){
             char firstcharacter=std::tolower(name.at(0), std::locale());
             int carriertype=0;
            if (firstcharacter=='e'){
                carriertype=-1;
            }
            else if(firstcharacter=='h'){
                carriertype=1;
            }
            else if(firstcharacter=='s'){
                carriertype=2;
            }
            else if(firstcharacter=='t'){
                carriertype=3;
            }
            else{
                throw runtime_error((boost::format("Carriername %s not known") % name).str());
            }
            return carriertype;
        }
         
         
         void KMCCalculator::RandomlyCreateCharges(){
         
        
        cout << "looking for injectable nodes..." << endl;
        for (unsigned int i = 0; i < _numberofcharges; i++) {
            Chargecarrier *newCharge = new Chargecarrier;
            newCharge->id = i;
            RandomlyAssignCarriertoSite(newCharge);
            
            cout << "starting position for charge " << i + 1 << ": segment " << newCharge->getCurrentNodeId()+1 << endl;
            _carriers.push_back(newCharge);
        }
        return;
         }
         
         void KMCCalculator::RandomlyAssignCarriertoSite(Chargecarrier* Charge){
            int nodeId_guess=-1;
            do{
            nodeId_guess=_RandomVariable.rand_uniform_int(_nodes.size());   
            }
            while (_nodes[nodeId_guess]->occupied || _nodes[nodeId_guess]->injectable==false ); // maybe already occupied? or maybe not injectable?
            if (Charge->hasNode()){
                Charge->jumpfromCurrentNodetoNode(_nodes[nodeId_guess]);
            }
            else{
            Charge->settoNote(_nodes[nodeId_guess]);
            }
             return;
         }
        
        void KMCCalculator::InitialRates() {
            
            cout << endl << "Calculating initial Marcus rates." << endl;
            cout << "    Temperature T = " << _temperature << " K." << endl;
           
            cout << "    carriertype: " << CarrierInttoLongString(_carriertype) << endl;
            unsigned numberofsites = _nodes.size();
            cout << "    Rates for " << numberofsites << " sites are computed." << endl;
            double charge=0.0;
            if (_carriertype == -1)
            {
                charge = -1.0;
            }
            else if (_carriertype == 1)
            {
                charge = 1.0;
            }
            cout<<"electric field ="<<_field<<" V/nm"<<endl;
            
            double maxreldiff = 0;
            double maxrate=0;
            double minrate=std::numeric_limits<double>::max();
            int totalnumberofrates = 0;
            for (unsigned int i = 0; i < numberofsites; i++) {
                unsigned numberofneighbours = _nodes[i]->events.size();
                for (unsigned int j = 0; j < numberofneighbours; j++) {
                    if(_nodes[i]->events[j].decayevent){
                        //if event is a decay event there is no point in calculating its rate, because it already has that from the reading in.
                        continue;
                    }

                    double destindex = _nodes[i]->events[j].destination;
                    double reorg = _nodes[i]->reorg_intorig + _nodes[destindex]->reorg_intdest + _nodes[i]->events[j].reorg_out;
                     if(std::abs(reorg)<1e-12){
                        throw std::runtime_error("Reorganisation energy for a pair is extremly close to zero,\n"
                                " you probably forgot to import reorganisation energies into your sql file.");
                    }
                    double dG_Field =0.0;
                    if(charge!=0.0){
                        dG_Field=charge * (_nodes[i]->events[j].dr*_field);
                    }
                    double dG_Site = _nodes[destindex]->siteenergy - _nodes[i]->siteenergy;
                    double dG=dG_Site-dG_Field;
                    double J2 = _nodes[i]->events[j].Jeff2;

                    double rate = 2 * tools::conv::Pi / tools::conv::hbar * J2 / sqrt(4 * tools::conv::Pi * reorg * tools::conv::kB * _temperature) 
                    * exp(-(dG + reorg)*(dG + reorg) / (4 * reorg * tools::conv::kB * _temperature));

                    // calculate relative difference compared to values in the table
                    double reldiff = (_nodes[i]->events[j].rate - rate) / _nodes[i]->events[j].rate;
                    if (reldiff > maxreldiff) {
                        maxreldiff = reldiff;
                    }
                    reldiff = (_nodes[i]->events[j].rate - rate) / rate;
                    if (reldiff > maxreldiff) {
                        maxreldiff = reldiff;
                    }

                    // set rates to calculated values
                    _nodes[i]->events[j].rate = rate;
                    _nodes[i]->events[j].initialrate = rate;
                    
                    if(rate>maxrate){
                        maxrate=rate;
                    }
                    else if(rate<minrate){
                        minrate=rate;
                    }

                    totalnumberofrates++;
                }

                

            }
            
            // Initialise escape rates
                for (auto* node:_nodes) {
                    node->InitEscapeRate();
                }
            
            cout << "    " << totalnumberofrates << " rates have been calculated." << endl;
            cout<< " Largest rate="<<maxrate<<" 1/s  Smallest rate="<<minrate<<" 1/s"<<endl;
            if (maxreldiff < 0.01) {
                cout << "    Good agreement with rates in the state file. Maximal relative difference: " << maxreldiff * 100 << " %" << endl;
            } else {
                cout << "    WARNING: Rates differ from those in the state file up to " << maxreldiff * 100 << " %." << " If the rates in the state file are calculated for a different temperature/field or if they are not Marcus rates, this is fine. Otherwise something might be wrong here." << endl;
            }
            
            return;
        }
        
        
        
        
        double KMCCalculator::Promotetime(double cumulated_rate){
            double dt = 0;
                double rand_u = 1 - _RandomVariable.rand_uniform();
                while (rand_u == 0) {
                    cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                    rand_u = 1 - _RandomVariable.rand_uniform();
                }
                dt = -1 / cumulated_rate * log(rand_u);
            return dt;
        }
        
        
        GLink* KMCCalculator::ChooseHoppingDest(GNode* node){
            double u = 1 - _RandomVariable.rand_uniform();
            
            for (unsigned int j = 0; j < node->events.size(); j++) {
                u -= node->events[j].rate / node->getEscapeRate();
                if (u <= 0 || j==node->events.size()-1) {    
                    return &(node->events[j]);
                }
            }
            throw runtime_error("Choose Hopping Destination, somehow no event was found");
            return NULL;
        }
        
        Chargecarrier* KMCCalculator::ChooseAffectedCarrier(double cumulated_rate){
            if(_carriers.size()==1){
                return _carriers[0];
            }
            Chargecarrier* carrier=NULL;
            double u = 1 - _RandomVariable.rand_uniform();
            for (unsigned int i = 0; i < _numberofcharges; i++) {
                u -= _carriers[i]->getCurrentEscapeRate() / cumulated_rate;

                if (u <= 0 || i==_numberofcharges-1) {

                    carrier = _carriers[i];
                    break;}  

            }
            return carrier;
        }
        
        void KMCCalculator::AddtoJumplengthdistro(const GLink* event,double dt){
            if(dolengthdistributon){
            double dist=abs(event->dr)-minlength;
            int index=int(dist/lengthresolution);
           
            _jumplengthdistro[index]++;
            _jumplengthdistro_weighted[index]+=dt;
           
            
            }
            return; 
        }
        
        void KMCCalculator::PrintJumplengthdistro(){
            if(dolengthdistributon){
            long unsigned noofjumps=0;
            
            double weightedintegral=0;
            for(unsigned i=0;i<_jumplengthdistro.size();++i){
                noofjumps+=_jumplengthdistro[i];   
                weightedintegral+=_jumplengthdistro_weighted[i];
            }
            double noofjumps_double=double(noofjumps);
            cout<<"Total number of jumps: "<<noofjumps<<endl;
            cout<<" distance[nm] |   # of jumps   | # of jumps [%] | .w. by dist [nm] | w. by timestep [%]"<<endl;
            cout<<"------------------------------------------------------------------------------------"<<endl;
            for(unsigned i=0;i<_jumplengthdistro.size();++i){
                double dist=lengthresolution*(i+0.5)+minlength;
                double percent=_jumplengthdistro[i]/noofjumps_double;
                double rtimespercent=percent*dist;
                cout<<(boost::format("    %4.3f    | %15d |    %04.2f    |     %4.3e     |     %04.2f")
                            % (dist) % (_jumplengthdistro[i]) % (percent*100) %(rtimespercent) % (_jumplengthdistro_weighted[i]/weightedintegral*100)).str()<<endl;                
            }
            cout<<"------------------------------------------------------------------------------------"<<endl;
           
            }
            return;
        }
 
        
    }
}
