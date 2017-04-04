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

#include "calculators/kmclifetime.h"
#include <votca/xtp/gnode.h>
#include <votca/tools/property.h>
#include <votca/tools/constants.h>
#include <boost/format.hpp>
#include <votca/ctp/topology.h>
#include <locale>


using namespace std;

namespace votca {
    namespace xtp {

        
        
    void KMCLifetime::Initialize( tools::Property *options) {
            
    std::string key = "options." + Identify();

    _insertions=options->ifExistsReturnElseThrowRuntimeError<unsigned int>(key+".numberofinsertions");
    _seed=options->ifExistsReturnElseThrowRuntimeError<int>(key+"seed");
    _numberofcharges=options->ifExistsReturnElseThrowRuntimeError<int>(key+".numberofcharges");
    _injection_name=options->ifExistsReturnElseThrowRuntimeError<int>(key+".injectionpattern");
    _lifetimefile=options->ifExistsReturnElseThrowRuntimeError<int>(key+".lifetimefile");

    _maxrealtime=options->ifExistsReturnElseReturnDefault<double>(key+".maxrealtime",1E10);
     _trajectoryfile=options->ifExistsReturnElseReturnDefault<std::string>(key+".trajectoryfile","trajectory.csv");
    _temperature=options->ifExistsReturnElseReturnDefault<double>(key+".temperature",300);
    _rates=options->ifExistsReturnElseReturnDefault<std::string>(key+".rates","statefile");

     std::string subkey=key+".carrierenergy";
    if (options->exists(subkey)) {
        _do_carrierenergy=options->ifExistsReturnElseReturnDefault<bool>(subkey+".run",false);
        _energy_outputfile = options->ifExistsReturnElseReturnDefault<std::string>(subkey+".outputfile","energy.tab");
         _alpha = options->ifExistsReturnElseReturnDefault<double>(subkey+".alpha",0.3);
         _outputsteps=options->ifExistsReturnElseReturnDefault<double>(subkey+".outputsteps",100);

    } else {
        _do_carrierenergy=false;
    }

    _carriertype = 2;
    cout << "carrier type: singlets" << endl;



    if (_rates != "statefile" && _rates != "calculate") {
        cout << "WARNING in kmclifetime: Invalid option rates. Valid options are 'statefile' or 'calculate'. Setting it to 'statefile'." << endl;
        _rates = "statefile";
    }
    
    return;
    }


    void KMCLifetime::LoadGraph(ctp::Topology *top) {

        vector< ctp::Segment* >& seg = top->Segments();

        for (unsigned i = 0; i < seg.size(); i++) {
            GNode *newNode = new GNode();
            _nodes.push_back(newNode);
            _nodes[i]->ReadfromSegment(seg[i], _carriertype);
            if (tools::wildcmp(_injection_name.c_str(), seg[i]->getName().c_str())) {
                _nodes[i]->injectable = 1;
            } else {
                _nodes[i]->injectable = 0;
            }
        }

        ctp::QMNBList &nblist = top->NBList();

        for (ctp::QMNBList::iterator it = nblist.begin(); it < nblist.end(); ++it) {
            _nodes[(*it)->Seg1()->getId()]->AddEventfromQmPair(*it, _carriertype);
        }


        if (votca::tools::globals::verbose) {
            cout << "pairs: " << nblist.size() / 2 << endl;
        }
        cout << "spatial density: " << _numberofcharges / top->BoxVolume() << " nm^-3" << endl;

        for (unsigned int i = 0; i < _nodes.size(); i++) {
            _nodes[i]->InitEscapeRate();
        }
            
        return;
    }
    

        void KMCLifetime::ResetForbiddenlist(std::vector<int> &forbiddenid) {
            forbiddenid.clear();
            return;
        }

        void KMCLifetime::AddtoForbiddenlist(int id, std::vector<int> &forbiddenid) {
            forbiddenid.push_back(id);
            return;
        }

        bool KMCLifetime::CheckForbidden(int id,const std::vector<int> &forbiddenlist) {
            // cout << "forbidden list has " << forbiddenlist.size() << " entries" << endl;
            bool forbidden = false;
            for (unsigned int i = 0; i < forbiddenlist.size(); i++) {
                if (id == forbiddenlist[i]) {
                    forbidden = true;
                    //cout << "ID " << id << " has been found as element " << i << " (" << forbiddenlist[i]<< ") in the forbidden list." << endl;
                    break;
                }
            }
            return forbidden;
        }

        bool KMCLifetime::CheckSurrounded(GNode* node,const std::vector<int> & forbiddendests) {
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

        
        void KMCLifetime::ReadLifetimeFile(std::string filename){
            tools::Property xml;
            load_property_from_xml(xml, filename);
            list<tools::Property*> jobProps = xml.Select("lifetimes.site");
            if (jobProps.size()!=_nodes.size()){
                throw  runtime_error((boost::format("The number of sites in the sqlfile: %i does not match the number in the lifetimefile: %i")
                        % _nodes.size() % jobProps.size()).str()); 
            }
            
            for (list<tools::Property*> ::iterator  it = jobProps.begin(); it != jobProps.end(); ++it) {
                int site_id =(*it)->getAttribute<int>("id");
                double lifetime=boost::lexical_cast<double>((*it)->value());
                bool check=false;
                for (unsigned i=0;i<_nodes.size();i++){
                    if (_nodes[i]->id==site_id-1 && !(_nodes[i]->hasdecay)){
                        _nodes[i]->AddDecayEvent(1.0/lifetime);
                        check=true;
                        break;
                    }
                    else if(_nodes[i]->id==site_id && _nodes[i]->hasdecay){
                        throw runtime_error((boost::format("Node %i appears twice in your list") %site_id).str());
                    } 
                    
                }
                if (!check){
                throw runtime_error((boost::format("Site from file with id: %i not found in sql") %site_id).str());
                }
            }
            
            return;
        }

        
        
        std::string KMCLifetime::CarrierInttoLongString(int carriertype){
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
        
        std::string KMCLifetime::CarrierInttoShortString(int carriertype){
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
        
         int KMCLifetime::StringtoCarriertype(std::string name){
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
        
        void KMCLifetime::InitialRates() {
            
            cout << endl << "Calculating initial Marcus rates." << endl;
            cout << "    Temperature T = " << _temperature << " K." << endl;
           
            cout << "    carriertype: " << CarrierInttoLongString(_carriertype) << endl;
            unsigned numberofsites = _nodes.size();
            cout << "    Rates for " << numberofsites << " sites are computed." << endl;
            double maxreldiff = 0;
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

                    double dG = _nodes[destindex]->siteenergy - _nodes[i]->siteenergy;

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

                    totalnumberofrates++;
                }

                // Initialise escape rates
                for (unsigned int i = 0; i < _nodes.size(); i++) {
                    _nodes[i]->InitEscapeRate();
                }

            }
            
            cout << "    " << totalnumberofrates << " rates have been calculated." << endl;
            if (maxreldiff < 0.01) {
                cout << "    Good agreement with rates in the state file. Maximal relative difference: " << maxreldiff * 100 << " %" << endl;
            } else {
                cout << "    WARNING: Rates differ from those in the state file up to " << maxreldiff * 100 << " %." << " If the rates in the state file are calculated for a different temperature/field or if they are not Marcus rates, this is fine. Otherwise something might be wrong here." << endl;
            }
            
            return;
        }
        
        
        
        double KMCLifetime::Promotetime(double cumulated_rate){
            double dt = 0;
                double rand_u = 1 - _RandomVariable->rand_uniform();
                while (rand_u == 0) {
                    cout << "WARNING: encountered 0 as a random variable! New try." << endl;
                    rand_u = 1 - _RandomVariable->rand_uniform();
                }
                dt = -1 / cumulated_rate * log(rand_u);
            return dt;
        }
        
        
        GLink* KMCLifetime::ChooseHoppingDest(GNode* node){
            double u = 1 - _RandomVariable->rand_uniform();
            
            for (unsigned int j = 0; j < node->events.size(); j++) {
                u -= node->events[j].rate / node->getEscapeRate();
                if (u <= 0 || j==node->events.size()-1) {                   
                    return &(node->events[j]);
                }
            }
            throw runtime_error("Choose Hopping Destination, somehow no event was found");
            return NULL;
        }
        

        void  KMCLifetime::RunVSSM(ctp::Topology *top) {
            
            int realtime_start = time(NULL);
            cout << endl << "Algorithm: VSSM for Multiple Charges with finite Lifetime" << endl;
            cout << "number of charges: " << _numberofcharges << endl;
            cout << "number of nodes: " << _nodes.size() << endl;

            if (_numberofcharges > _nodes.size()) {
                throw runtime_error("ERROR in kmclifetime: specified number of charges is greater than the number of nodes. This conflicts with single occupation.");
            }

            fstream traj;
            fstream energyfile;

            cout << "Writing trajectory to " <<  _trajectoryfile << "." << endl; 
            traj.open ( _trajectoryfile.c_str(), fstream::out);
            traj << "#Simtime [s]\t Insertion\t Carrier ID\t Lifetime[s]\tSteps\t Last Segment\t x_travelled[nm]\t y_travelled[nm]\t z_travelled[nm]"<<endl;
            
            if(_do_carrierenergy){

                cout << "Tracking the energy of one charge carrier and exponential average with alpha=" << _alpha << " to "<<_energy_outputfile << endl;
                energyfile.open(_energy_outputfile.c_str(),fstream::out);
                energyfile << "Simtime [s]\tSteps\tCarrier ID\tEnergy_a="<<_alpha<<"[eV]"<<endl;
            }
            
            // Injection
            cout << endl << "injection method: " << _injectionmethod << endl;
            
            std::vector< Chargecarrier* > carrier;
            std::vector<tools::vec> startposition(_numberofcharges, tools::vec(0, 0, 0));
            cout << "looking for injectable nodes..." << endl;
            for (unsigned int i = 0; i < _numberofcharges; i++) {
                Chargecarrier *newCharge = new Chargecarrier;
                newCharge->id = i;
                do {
                newCharge->node = _nodes[_RandomVariable->rand_uniform_int(_nodes.size())];     
                }
                while (newCharge->node->occupied == 1 || newCharge->node->injectable != 1 ); // maybe already occupied? or maybe not injectable?
                    
                
                // cout << "selected segment " << newCharge->node->id+1 << " which has energy " << newCharge->node->siteenergy << " within the interval [" << energypercarrier-0*deltaE << ", " << energypercarrier+2*deltaE << "]" << endl;
                newCharge->node->occupied = 1;
                cout << "starting position for charge " << i + 1 << ": segment " << newCharge->node->id + 1 << endl;
                carrier.push_back(newCharge);
            }

            unsigned insertioncount = 0;
            unsigned long step=0;
            double simtime=0.0;

            std::vector<int> forbiddennodes;
            std::vector<int> forbiddendests;
            
            time_t now = time(0);
            tm* localtm = localtime(&now);
            cout << "Run started at " << asctime(localtm) << endl;
 
            double avlifetime=0.0;
            double meanfreepath=0.0;
            tools::vec difflength=tools::vec(0,0,0);
            unsigned steps=0;
            double avgenergy=carrier[0]->node->siteenergy;
            int     carrieridold=carrier[0]->id;
            
            while (insertioncount < _insertions) {
                if ((time(NULL) - realtime_start) > _maxrealtime * 60. * 60.) {
                    cout << endl << "Real time limit of " << _maxrealtime << " hours (" << int(_maxrealtime * 60 * 60 + 0.5) << " seconds) has been reached. Stopping here." << endl << endl;
                    break;
                }

                step += 1;
                double cumulated_rate = 0;

                for (unsigned int i = 0; i < carrier.size(); i++) {
                    cumulated_rate += carrier[i]->node->getEscapeRate();
                }
                if (cumulated_rate == 0) { // this should not happen: no possible jumps defined for a node
                    throw runtime_error("ERROR in kmclifetime: Incorrect rates in the database file. All the escape rates for the current setting are 0.");
                }
                // go forward in time
                double dt=Promotetime(cumulated_rate);
                
                if(_do_carrierenergy){
                    bool print=false;
                    if (carrier[0]->id>carrieridold){
                        avgenergy=carrier[0]->node->siteenergy;
                        print=true;
                        carrieridold=carrier[0]->id;
                    }
                    else if(step%_outputsteps==0){
                        avgenergy=_alpha*carrier[0]->node->siteenergy+(1-_alpha)*avgenergy;
                        print=true;
                    }
                    if(print){
                        energyfile << simtime<<"\t"<<steps <<"\t"<<carrier[0]->id<<"\t"<<avgenergy<<endl;                  
                    }
                }
                
                simtime += dt;
                steps++;
                for (unsigned int i = 0; i < carrier.size(); i++) {
                    carrier[i]->updateLifetime(dt);
                    carrier[i]->updateSteps(1);
                    carrier[i]->node->occupationtime += dt;
                }
                
                ResetForbiddenlist(forbiddennodes);
                bool secondlevel=true;
                while (secondlevel){

                    // determine which carrier will escape
                    GNode* oldnode=NULL;
                    GNode* newnode=NULL;
                    Chargecarrier* affectedcarrier=NULL;

                    double u = 1 - _RandomVariable->rand_uniform();
                    for (unsigned int i = 0; i < _numberofcharges; i++) {
                        u -= carrier[i]->node->getEscapeRate() / cumulated_rate;
                       
                        if (u <= 0 || i==_numberofcharges-1) {
                            oldnode = carrier[i]->node;
                            affectedcarrier = carrier[i];
                            break;}  
                       
                    }
                    tools::vec *dr= NULL;
                   
                    if (CheckForbidden(oldnode->id, forbiddennodes)) {
                        continue;
                    }

                    // determine where it will jump to
                    ResetForbiddenlist(forbiddendests);
                    
                    while (true) {
                        // LEVEL 2
 
                        newnode = NULL;
                        GLink* event=ChooseHoppingDest(oldnode);
                       
                        if (event->decayevent){
                            oldnode->occupied = 0;
                            avlifetime+=affectedcarrier->getLifetime();
                            meanfreepath+=tools::abs(affectedcarrier->dr_travelled);
                            difflength+=tools::elementwiseproduct(affectedcarrier->dr_travelled,affectedcarrier->dr_travelled);
                            traj << simtime<<"\t"<<insertioncount<< "\t"<< affectedcarrier->id<<"\t"<< affectedcarrier->getLifetime()<<"\t"<<affectedcarrier->getSteps()<<"\t"<< (oldnode->id)+1<<"\t"<<affectedcarrier->dr_travelled.getX()<<"\t"<<affectedcarrier->dr_travelled.getY()<<"\t"<<affectedcarrier->dr_travelled.getZ()<<endl;
                            if(_insertions<1500 ||insertioncount% (_insertions/1000)==0 || insertioncount<0.001*_insertions){
                            std::cout << "\rInsertion " << insertioncount<<" of "<<_insertions;
                            std::cout << std::flush;
                            }
                            do {
                                affectedcarrier->node = _nodes[_RandomVariable->rand_uniform_int(_nodes.size())];     
                            }
                            while (affectedcarrier->node->occupied == 1 || affectedcarrier->node->injectable != 1 );
                            affectedcarrier->node->occupied=1;
                            affectedcarrier->resetCarrier();
                            insertioncount++;
                            affectedcarrier->id=_numberofcharges-1+insertioncount;
                            secondlevel=false;
                            break;
                                }
                        else{
                        newnode = _nodes[event->destination];
                        dr = &(event->dr);
                        }
                    
                        // check after the event if this was allowed
                        if (CheckForbidden(newnode->id, forbiddendests)) {
                            continue;
                        }

                        // if the new segment is unoccupied: jump; if not: add to forbidden list and choose new hopping destination
                        if (newnode->occupied == 1) {
                            if (CheckSurrounded(oldnode, forbiddendests)) {     
                                AddtoForbiddenlist(oldnode->id, forbiddennodes);
                                break; // select new escape node (ends level 2 but without setting level1step to 1)
                            }
                            AddtoForbiddenlist(newnode->id, forbiddendests);
                            continue; // select new destination
                        } else {
                            newnode->occupied = 1;
                            oldnode->occupied = 0;
                            affectedcarrier->node = newnode;
                            affectedcarrier->dr_travelled += *(dr);
                            secondlevel=false;

                            break; // this ends LEVEL 2 , so that the time is updated and the next MC step started
                        }

                       
                        // END LEVEL 2
                    }
                    // END LEVEL 1
                }
            }
            
            
            cout<<endl;
            cout << "Total runtime:\t\t\t\t\t"<< simtime << " s"<< endl;
            cout << "Total KMC steps:\t\t\t\t"<< steps << endl;
            cout << "Average lifetime:\t\t\t\t"<<avlifetime/_insertions<< " s"<<endl;
            cout << "Mean freepath\t l=<|r_x-r_o|> :\t\t"<<(meanfreepath/_insertions)<< " nm"<<endl;
            cout << "Average diffusionlength\t d=sqrt(<(r_x-r_o)^2>)\t"<<sqrt(abs(difflength)/_insertions)<< " nm"<<endl;
            cout<<endl;
            
            vector< ctp::Segment* >& seg = top->Segments();

            for (unsigned i = 0; i < seg.size(); i++) {
                double occupationprobability=_nodes[i]->occupationtime / simtime;
                seg[i]->setOcc(occupationprobability,_carriertype);
            }
    
            return;
        }



    bool KMCLifetime::EvaluateFrame(ctp::Topology *top) {
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "      KMCLIFETIME started" << std::endl;
        std::cout << "-----------------------------------" << std::endl << std::endl;

        // Initialise random number generator
        if (votca::tools::globals::verbose) {
            cout << endl << "Initialising random number generator" << endl;
        }
        std::srand(_seed); // srand expects any integer in order to initialise the random number generator
        _RandomVariable = new tools::Random2();
        _RandomVariable->init(rand(), rand(), rand(), rand());

        LoadGraph(top);
        ReadLifetimeFile(_lifetimefile);
        InitialRates();
        RunVSSM(top);

        time_t now = time(0);
        tm* localtm = localtime(&now);
        std::cout << "      KMCLIFETIME finished at:" <<asctime(localtm) <<  std::endl;

        return true;
    }
    
    }
}
